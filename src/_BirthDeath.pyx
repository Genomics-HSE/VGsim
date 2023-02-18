# cython: language_level=3
# cython: initializedcheck = False
# distutils: language = c++

cimport cython

from numpy cimport npy_int64

from libc.math cimport log, floor, abs
from libcpp.vector cimport vector
from mc_lib.rndm cimport RndmWrapper

from prettytable import PrettyTable
import numpy as np
cimport numpy as np
import sys
import os
import time
import tskit

from numpy.random.c_distributions cimport random_poisson, random_hypergeometric

include "fast_choose.pxi"
include "models.pxi"
include "events.pxi"

#pi - population ID, pn - popoulation number, spi - source population ID, tpi - target population ID
#hi - haplotype ID, hn - haplotype number, nhi - new haplotype number
#si - susceptibility ID, sn - susceptibility number, ssi - source susceptibility ID, tsi - target susceptibility ID
cdef class BirthDeathModel:
    cdef:
        RndmWrapper seed

        bint first_simulation, sampling_probability, memory_optimization
        Py_ssize_t user_seed, sites, hapNum, currentHapNum, maxHapNum, addMemoryNum, popNum, susNum, bCounter, dCounter, sCounter, \
        mCounter, iCounter, swapLockdown, migPlus, migNonPlus, globalInfectious, countsPerStep, good_attempt, genome_length
        double currentTime, totalRate, totalMigrationRate, rn, tau_l, recombination, SuperSpreadRate

        Events events
        multiEvents multievents
        Mutations mut
        Migrations mig
        Lockdowns loc
        Recombination rec

        npy_int64[::1] suscType, sizes, totalSusceptible, totalInfectious, lockdownON, hapToNum, numToHap, sitesPosition
        # npy_int64[::1] tree, tree_pop
        npy_int64[:,::1] susceptible, infectious, initial_susceptible, initial_infectious
        npy_int64[:,::1] tree, tree_pop

        double[::1] bRate, dRate, sRate, tmRate, maxEffectiveBirthMigration, suscepCumulTransition, immunePopRate, infectPopRate, \
        popRate, migPopRate, actualSizes, contactDensity, contactDensityBeforeLockdown, contactDensityAfterLockdown, startLD, endLD, \
        samplingMultiplier, times, birthInf
        double[:,::1] mRate, susceptibility, tEventHapPopRate, suscepTransition, immuneSourcePopRate, hapPopRate, migrationRates, \
        effectiveMigration
        double[:,:,::1] hapMutType, eventHapPopRate, susceptHapPopRate

        double[:,:,:,::1] PropensitiesMigr, PropensitiesMutatations
        double[:,:,::1] PropensitiesSuscep, PropensitiesTransmission
        double[:,::1] PropensitiesRecovery, PropensitiesSampling

        npy_int64[:,:,:,::1] eventsMigr, eventsMutatations
        npy_int64[:,:,::1] eventsSuscep, eventsTransmission
        npy_int64[:,::1] eventsRecovery, eventsSampling

        double[:,::1] infectiousAuxTau, susceptibleAuxTau
        npy_int64[:,::1] infectiousDelta, susceptibleDelta

        vector[double] super_spread_rate
        vector[Py_ssize_t] super_spread_left, super_spread_right, super_spread_pop


    def __init__(self, number_of_sites, populations_number, number_of_susceptible_groups, seed, sampling_probability, \
        memory_optimization, genome_length, recombination_probability):
        self.check_amount(seed, 'seed', zero=False)
        self.user_seed = seed
        self.seed = RndmWrapper(seed=(self.user_seed, 0))

        self.first_simulation = False
        if sampling_probability != True and sampling_probability != False:
            raise ValueError('Incorrect value of sampling probability. Value of sampling probability should be True or False.')
        self.sampling_probability = sampling_probability
        if memory_optimization != True and memory_optimization != False:
            raise ValueError('Incorrect value of memory optimization. Value of memory optimization should be True or False.')
        self.memory_optimization = memory_optimization

        self.check_amount(number_of_sites, 'number of sites', zero=False)
        self.sites = number_of_sites
        self.hapNum = 4**self.sites
        self.check_amount(number_of_susceptible_groups, 'number of susceptible groups')
        self.susNum = number_of_susceptible_groups
        self.check_amount(populations_number, 'populations number')
        self.popNum = populations_number

        self.check_value(recombination_probability, 'recombination probability', edge=1)
        self.recombination = recombination_probability
        self.check_amount(genome_length, 'genome length')
        self.genome_length = genome_length
        self.sitesPosition = np.zeros(self.sites, dtype=np.int64)
        if self.sites > self.genome_length:
            raise ValueError('Incorrect value of number of sites or genome length. Genome length should be more or equal number of sites.')
        if self.sites > 1:
            for s in range(self.sites):
                self.sitesPosition[s] = int(s * self.genome_length / (self.sites - 1))
            self.birthInf = np.zeros(self.hapNum, dtype=float)


        if self.memory_optimization:
            if self.sites > 2:
                self.maxHapNum = 4**(self.sites-2)
                self.addMemoryNum = 4**(self.sites-2)
            else:
                self.maxHapNum = 4
                self.addMemoryNum = 4
        else:
            self.maxHapNum = self.hapNum
            self.addMemoryNum = 0

        self.SuperSpreadRate = 0

        self.hapToNum = np.zeros(self.hapNum, dtype=np.int64) # from haplotype to program number
        self.numToHap = np.zeros(self.maxHapNum, dtype=np.int64) # from program number to haplotype

        if self.memory_optimization:
            self.currentHapNum = 0
        else:
            self.currentHapNum = self.hapNum
            for hn in range(self.hapNum):
                self.hapToNum[hn] = hn
                self.numToHap[hn] = hn

        self.bCounter = 0
        self.dCounter = 0
        self.sCounter = 0
        self.mCounter = 0
        self.iCounter = 0
        self.swapLockdown = 0
        self.migPlus = 0
        self.migNonPlus = 0
        self.globalInfectious = 0

        self.currentTime = 0.0
        self.tau_l=0.01
        self.totalRate = 0.0
        self.totalMigrationRate = 0.0

        self.events = Events()
        self.multievents = multiEvents()
        self.mut = Mutations()
        self.mig = Migrations()
        self.loc = Lockdowns()
        self.rec = Recombination()

        # memory_optimization
        self.infectious = np.zeros((self.popNum, self.maxHapNum), dtype=np.int64)
        self.tmRate = np.zeros(self.maxHapNum, dtype=float)
        self.tEventHapPopRate = np.zeros((self.popNum, self.maxHapNum), dtype=float)
        self.hapPopRate = np.zeros((self.popNum, self.maxHapNum), dtype=float)
        self.eventHapPopRate = np.zeros((self.popNum, self.maxHapNum, 4), dtype=float)
        self.susceptHapPopRate = np.zeros((self.popNum, self.maxHapNum, self.susNum), dtype=float)

        self.suscType = np.zeros(self.hapNum, dtype=np.int64)
        self.bRate = np.zeros(self.hapNum, dtype=float)
        self.dRate = np.zeros(self.hapNum, dtype=float)
        self.sRate = np.zeros(self.hapNum, dtype=float)
        self.mRate = np.zeros((self.hapNum, self.sites), dtype=float)
        self.susceptibility = np.zeros((self.hapNum, self.susNum), dtype=float)
        self.hapMutType = np.ones((self.hapNum, self.sites, 3), dtype=float)

        for hn in range(self.hapNum):
            self.bRate[hn] = 2.0
            self.dRate[hn] = 1.0
            self.sRate[hn] = 0.01
            for s in range(self.sites):
                self.mRate[hn, s] = 0.01
            self.susceptibility[hn, 0] = 1.0

        self.sizes = np.zeros(self.popNum, dtype=np.int64)
        self.totalSusceptible = np.zeros(self.popNum, dtype=np.int64)
        self.totalInfectious = np.zeros(self.popNum, dtype=np.int64)
        self.lockdownON = np.zeros(self.popNum, dtype=np.int64)

        self.susceptible = np.zeros((self.popNum, self.susNum), dtype=np.int64)
        self.initial_susceptible = np.zeros((self.popNum, self.susNum), dtype=np.int64)
        self.initial_infectious = np.zeros((self.popNum, self.hapNum), dtype=np.int64)

        self.maxEffectiveBirthMigration = np.zeros(self.popNum, dtype=float)
        self.suscepCumulTransition = np.zeros(self.susNum, dtype=float)
        self.infectPopRate = np.zeros(self.popNum, dtype=float)
        self.immunePopRate = np.zeros(self.popNum, dtype=float)
        self.popRate = np.zeros(self.popNum, dtype=float)
        self.migPopRate = np.zeros(self.popNum, dtype=float)
        self.actualSizes = np.zeros(self.popNum, dtype=float)
        self.contactDensity = np.ones(self.popNum, dtype=float)
        self.contactDensityBeforeLockdown = np.ones(self.popNum, dtype=float)
        self.contactDensityAfterLockdown = np.zeros(self.popNum, dtype=float)
        self.startLD = np.ones(self.popNum, dtype=float)
        self.endLD = np.ones(self.popNum, dtype=float)
        self.samplingMultiplier = np.ones(self.popNum, dtype=float)

        self.suscepTransition = np.zeros( (self.susNum, self.susNum), dtype=float)
        self.immuneSourcePopRate = np.zeros((self.popNum, self.susNum), dtype=float)
        self.migrationRates = np.zeros((self.popNum, self.popNum), dtype=float)
        self.effectiveMigration = np.zeros((self.popNum, self.popNum), dtype=float)

        for pn in range(self.popNum):
            self.sizes[pn] = 1000000
            self.totalSusceptible[pn] = 1000000
            self.susceptible[pn, 0] = 1000000

        # self.tree = np.zeros(1, dtype=np.int64)
        # self.tree_pop = np.zeros(1, dtype=np.int64)
        self.tree = np.zeros((1, 1), dtype=np.int64)
        self.tree_pop = np.zeros((1, 1), dtype=np.int64)

        #Init propensities
        self.PropensitiesMigr = np.zeros((self.popNum, self.popNum, self.susNum, self.hapNum), dtype=float)
        self.PropensitiesSuscep = np.zeros((self.popNum, self.susNum, self.susNum), dtype=float)
        self.PropensitiesRecovery = np.zeros((self.popNum, self.hapNum), dtype=float)
        self.PropensitiesSampling = np.zeros((self.popNum, self.hapNum), dtype=float)
        self.PropensitiesMutatations = np.zeros((self.popNum, self.hapNum, self.sites, 3), dtype=float)
        self.PropensitiesTransmission = np.zeros((self.popNum, self.hapNum, self.susNum), dtype=float)

        #Sampled number of events per step
        self.eventsMigr = np.zeros((self.popNum, self.popNum, self.susNum, self.hapNum), dtype=np.int64)
        self.eventsSuscep = np.zeros((self.popNum, self.susNum, self.susNum), dtype=np.int64)
        self.eventsRecovery = np.zeros((self.popNum, self.hapNum), dtype=np.int64)
        self.eventsSampling = np.zeros((self.popNum, self.hapNum), dtype=np.int64)
        self.eventsMutatations = np.zeros((self.popNum, self.hapNum, self.sites, 3), dtype=np.int64)
        self.eventsTransmission = np.zeros((self.popNum, self.hapNum, self.susNum), dtype=np.int64)

        self.infectiousAuxTau = np.zeros((self.popNum, self.hapNum), dtype=float)
        self.susceptibleAuxTau = np.zeros((self.popNum, self.susNum), dtype=float)

        self.infectiousDelta = np.zeros((self.popNum, self.hapNum), dtype=np.int64)
        self.susceptibleDelta = np.zeros((self.popNum, self.susNum), dtype=np.int64)


    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void FirstInfection(self):
        if self.globalInfectious == 0:
            for sn in range(self.susNum):
                if self.susceptible[0, sn] != 0:
                    if self.memory_optimization:
                        self.AddHaplotype(0, 0)
                    self.NewInfections(0, sn, 0)
                    return

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef inline void NewInfections(self, Py_ssize_t pi, Py_ssize_t si, Py_ssize_t hi, Py_ssize_t num=1):
        self.susceptible[pi, si] -= num
        self.totalSusceptible[pi] -= num
        self.infectious[pi, hi] += num
        self.totalInfectious[pi] += num
        self.globalInfectious += num

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef inline void NewRecoveries(self, Py_ssize_t pi, Py_ssize_t si, Py_ssize_t hi, Py_ssize_t num=1):
        self.susceptible[pi, si] += num
        self.totalSusceptible[pi] += num
        self.infectious[pi, hi] -= num
        self.totalInfectious[pi] -= num
        self.globalInfectious -= num

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void AddMemory(self):
        if self.addMemoryNum + self.maxHapNum > self.hapNum:
            self.addMemoryNum = self.hapNum - self.maxHapNum
        self.numToHap = np.concatenate((self.numToHap, np.zeros(self.addMemoryNum, dtype=np.int64)))
        self.infectious = np.concatenate((self.infectious, np.zeros((self.popNum, self.addMemoryNum), dtype=np.int64)), axis=1)
        self.tmRate = np.concatenate((self.tmRate, np.zeros(self.addMemoryNum, dtype=float)))
        self.tEventHapPopRate = np.concatenate((self.tEventHapPopRate, np.zeros((self.popNum, self.addMemoryNum), dtype=float)), axis=1)
        self.hapPopRate = np.concatenate((self.hapPopRate, np.zeros((self.popNum, self.addMemoryNum), dtype=float)), axis=1)
        self.eventHapPopRate = np.concatenate((self.eventHapPopRate, np.zeros((self.popNum, self.addMemoryNum, 4), dtype=float)), axis=1)
        self.susceptHapPopRate = np.concatenate((self.susceptHapPopRate, np.zeros((self.popNum, self.addMemoryNum, self.susNum), dtype=float)), axis=1)
        self.maxHapNum += self.addMemoryNum

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void UpdateAllRates(self):
        cdef:
            double maxEffectiveBirth
            double[::1] maxEffectiveMigration

        for sn1 in range(self.susNum):
            self.suscepCumulTransition[sn1] = 0
            for sn2 in range(self.susNum):
                self.suscepCumulTransition[sn1] += self.suscepTransition[sn1, sn2]

        for pn1 in range(self.popNum):
            self.migrationRates[pn1, pn1] = 1.0
            self.actualSizes[pn1] = 0.0
            for pn2 in range(self.popNum):
                if pn1 == pn2:
                    continue
                self.migrationRates[pn1, pn1] -= self.migrationRates[pn1, pn2]
                self.actualSizes[pn1] += self.migrationRates[pn2, pn1] * self.sizes[pn2]
            self.actualSizes[pn1] += self.migrationRates[pn1, pn1] * self.sizes[pn1]

        self.totalRate = 0.0
        for pn in range(self.popNum):
            self.infectPopRate[pn] = 0
            self.immunePopRate[pn] = 0
            self.popRate[pn] = 0.
        for pn in range(self.popNum):
            for hn in range(self.currentHapNum): # hn - program number
                self.tmRate[hn] = 0
                for s in range(self.sites):
                    self.tmRate[hn] += self.mRate[self.numToHap[hn], s]

                self.eventHapPopRate[pn, hn, 0] = self.BirthRate(pn, hn)
                self.eventHapPopRate[pn, hn, 1] = self.dRate[self.numToHap[hn]]
                self.eventHapPopRate[pn, hn, 2] = self.sRate[self.numToHap[hn]] * self.samplingMultiplier[pn]
                self.eventHapPopRate[pn, hn, 3] = self.tmRate[hn]
                self.tEventHapPopRate[pn, hn] = 0
                for i in range(4):
                    self.tEventHapPopRate[pn, hn] += self.eventHapPopRate[pn, hn, i]
                self.hapPopRate[pn, hn] = self.tEventHapPopRate[pn, hn] * self.infectious[pn, hn]
                self.infectPopRate[pn] += self.hapPopRate[pn, hn]
            for sn in range(self.susNum):
                self.immuneSourcePopRate[pn, sn] = self.suscepCumulTransition[sn] * self.susceptible[pn, sn]
                self.immunePopRate[pn] += self.immuneSourcePopRate[pn, sn]
            self.popRate[pn] = self.infectPopRate[pn] + self.immunePopRate[pn]
            self.totalRate += self.popRate[pn]

        maxEffectiveMigration = np.zeros(self.popNum, dtype=float)
        for pn1 in range(self.popNum):
            for pn2 in range(self.popNum):
                if pn1 == pn2:
                    continue
                self.effectiveMigration[pn1, pn2] = 0.0
                for pn3 in range(self.popNum):
                    self.effectiveMigration[pn1, pn2] += self.migrationRates[pn1, pn3] * self.migrationRates[pn2, pn3] * \
                    self.contactDensity[pn3] / self.actualSizes[pn3]
                if self.effectiveMigration[pn1, pn2] > maxEffectiveMigration[pn2]:
                # if pn1 != pn2 and self.effectiveMigration[pn1, pn2] > maxEffectiveMigration[pn2]:
                    maxEffectiveMigration[pn2] = self.effectiveMigration[pn1, pn2]

        maxEffectiveBirth = 0.0
        for hn in range(self.currentHapNum):
            for sn in range(self.susNum):
                if self.bRate[self.numToHap[hn]] * self.susceptibility[self.numToHap[hn], sn] > maxEffectiveBirth:
                    maxEffectiveBirth = self.bRate[self.numToHap[hn]] * self.susceptibility[self.numToHap[hn], sn]

        self.totalMigrationRate = 0.0
        for pn in range(self.popNum):
            self.maxEffectiveBirthMigration[pn] = maxEffectiveMigration[pn] * maxEffectiveBirth
            self.migPopRate[pn] = self.maxEffectiveBirthMigration[pn] * self.totalSusceptible[pn] * (self.globalInfectious - \
                self.totalInfectious[pn])
            self.totalMigrationRate += self.migPopRate[pn]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void AddHaplotype(self, Py_ssize_t nhi, Py_ssize_t pi): # nhi - haplotype
        cdef:
            bint check = False
            Py_ssize_t mem_h, mem_inf

        if self.currentHapNum == self.maxHapNum and self.maxHapNum < self.hapNum:
            self.AddMemory()

        for i in range(1, self.currentHapNum+1):
            if check:
                self.hapToNum[mem_h] += 1
                self.infectious[pi, i], mem_inf = mem_inf, self.infectious[pi, i]
                self.numToHap[i], mem_h = mem_h, self.numToHap[i]
            if (self.numToHap[i] > nhi or self.numToHap[i] == 0) and self.numToHap[i-1] < nhi:
                mem_h = self.numToHap[i]
                check = True
                self.numToHap[i] = nhi
                mem_inf = self.infectious[pi, i]
                self.infectious[pi, i] = 0
                self.hapToNum[nhi] = i

        self.currentHapNum += 1
        self.UpdateAllRates()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef inline double BirthRate(self, Py_ssize_t pi, Py_ssize_t hi): # hi - program number
        cdef double ps = 0.0 

        for sn in range(self.susNum):
            self.susceptHapPopRate[pi, hi, sn] = self.susceptible[pi, sn] * self.susceptibility[self.numToHap[hi], sn]
            for pn in range(self.popNum):
                ps += self.susceptHapPopRate[pi, hi, sn] * self.migrationRates[pi, pn] * self.migrationRates[pi, pn] * \
                self.contactDensity[pn] / self.actualSizes[pn]

        return self.bRate[self.numToHap[hi]] * ps

        # for sn in range(self.susNum):
        #     self.susceptHapPopRate[pi, hi, sn] = self.susceptible[pi, sn] * self.susceptibility[self.numToHap[hi], sn]
        #     for pn in range(self.popNum):
        #         ps += self.susceptHapPopRate[pi, hi, sn] * self.effectiveMigration[pi, pn]

        # return self.bRate[self.numToHap[hi], cn] * ps


    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef void SimulatePopulation(self, Py_ssize_t iterations, Py_ssize_t sample_size, float time, Py_ssize_t attempts):
        cdef Py_ssize_t pi 

        self.PrepareParameters(iterations)
        self.CheckSizes()

        for i in range(attempts):
            self.seed = RndmWrapper(seed=(self.user_seed, i))
            if self.totalRate+self.totalMigrationRate != 0.0 and self.globalInfectious != 0:
                gs_index = 0
                gs_time = self.general_sampling_times[gs_index]
                while (self.events.ptr<self.events.size and (sample_size==-1 or self.sCounter<=sample_size) and (time==-1 or self.currentTime<time)):
                    self.SampleTime()
                    pi = self.GenerateEvent()
                    # self.Debug()
                    if self.totalRate == 0.0 or self.globalInfectious == 0:
                        break
                    self.CheckLockdown(pi)
                    if gs_index < self.sampling_event_number and self.currentTime > gs_time:
                        gs_index += 1
                        gs_time = self.general_sampling_times[gs_index]
                        self.general_sampling()







            if self.events.ptr <= 100 and iterations > 100:
                self.Restart()
            else:
                self.good_attempt = i+1
                break

        if self.totalRate == 0.0 or self.globalInfectious == 0:
            print('Simulation finished because no infections individuals remain!')
        if self.events.ptr>=self.events.size:
            print("Achieved maximal number of iterations.")
        if self.sCounter>sample_size and sample_size!=-1:
            print("Achieved sample size.")
        if self.currentTime>time and time != -1:
            print("Achieved internal time limit.")
        if self.sCounter <= 1:
            print('\033[41m{}\033[0m'.format('WARNING!'), 'Simulated less 2 samples, so genealogy will not work!')

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void PrepareParameters(self, Py_ssize_t iterations):
        self.events.CreateEvents(iterations)
        if self.first_simulation == False:
            self.FirstInfection()
            for pn in range(self.popNum):
                for sn in range(self.susNum):
                    self.initial_susceptible[pn, sn] = self.susceptible[pn, sn]
                for hn in range(self.hapNum):
                    self.initial_infectious[pn, hn] = self.infectious[pn, hn]
            self.first_simulation = True
        for pn in range(self.popNum):
            self.CheckLockdown(pn)
        self.UpdateAllRates()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void CheckSizes(self):
        check = False
        list_pop = []
        print('Actual sizes: ', end='')
        for pn in range(self.popNum):
            print(self.actualSizes[pn], end=' ')
            if abs(self.actualSizes[pn]/self.sizes[pn]-1) >= 0.1:
                check = True
                list_pop.append(str(pn))
        print()
        if check:
            print('\033[41m{}\033[0m'.format('WARNING!'), 'Actual population size in deme: ', end='')
            print(", ".join(list_pop) )
            print("\tis more than 10% different from the population size. The migration probabilities might be unrealistically high.")
            print("\tWe recommend to check your model with print_populations() method before proceding to simulation.")
            print("\tCheck the documentation file:https://vg-sim.readthedocs.io/en/latest/Migration.html for more details.")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef inline void SampleTime(self):
        cdef double tau = - log(self.seed.uniform()) / (self.totalRate + self.totalMigrationRate)
        self.currentTime += tau

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef Py_ssize_t GenerateEvent(self):
        cdef:
            Py_ssize_t pi, hi, ei
            double choose

        self.rn = self.seed.uniform()
        # choose = self.rn * (self.totalRate + self.totalMigrationRate)
        choose = self.rn * (self.totalRate + self.totalMigrationRate + self.SuperSpreadRate)
        if self.totalRate > choose:
            self.rn = choose / self.totalRate
            pi, self.rn = fastChoose(self.popRate, self.totalRate, self.rn)
            choose = self.rn * self.popRate[pi]
            if self.immunePopRate[pi] > choose:
                self.rn = choose / self.immunePopRate[pi]
                self.ImmunityTransition(pi)
            else:
                self.rn = (choose - self.immunePopRate[pi]) / self.infectPopRate[pi]
                hi, self.rn = fastChoose(self.hapPopRate[pi], self.infectPopRate[pi], self.rn) # hi - program number
                ei, self.rn = fastChoose(self.eventHapPopRate[pi, hi], self.tEventHapPopRate[pi, hi], self.rn)
                if ei == BIRTH:
                    self.Birth(pi, hi)
                elif ei == DEATH:
                    self.Death(pi, hi)
                elif ei == SAMPLING:
                    self.Sampling(pi, hi)
                else:
                    self.Mutation(pi, hi)
        elif self.totalMigrationRate + self.totalRate > choose:
            self.rn = (choose - self.totalRate) / self.totalMigrationRate
            pi = self.GenerateMigration()
        else:
            self.rn = (choose - self.totalRate - self.totalMigrationRate) / self.SuperSpreadRate
            pi = self.SuperSpreadEvent()
        return pi

    # @cython.boundscheck(False)
    # @cython.wraparound(False)
    # cdef Py_ssize_t GenerateEvent(self):
    #     cdef:
    #         Py_ssize_t pi, hi, ci, ei

    #     self.rn = self.seed.uniform()
    #     pi, self.rn = fastChoose(self.popRate, self.totalRate, self.rn)
    #     ei, self.rn = fastChoose(self.eventPopRate[pi], self.popRate[pi], self.rn)
    #     if ei == 1:
    #         hi, self.rn = fastChoose(self.infectHapPopRate[pi], self.infectPopRate[pi], self.rn) # hi - program number
    #         ci, self.rn = fastChoose(self.infectConHapPopRate[pi, hi], self.infectHapPopRate[pi, hi], self.rn)
    #         ei, self.rn = fastChoose(self.eventConHapPopRate[pi, hi, ci], self.cumulEventConHapPopRate[pi, hi, ci], self.rn)
    #         if ei == BIRTH:
    #             self.Birth(pi, hi, ci)
    #         elif ei == DEATH:
    #             self.Death(pi, hi, ci)
    #         elif ei == SAMPLING:
    #             self.Sampling(pi, hi, ci)
    #         else ei == MUTATION:
    #             self.Mutation(pi, hi, ci)
    #         else:
    #             self.ConditionTransition(pi)
    #     elif ei == 2:
    #         self.GenerateMigration(pi)
    #     elif ei == 3:
    #         self.ImmunityTransition(pi)
    #     else:
    #         self.SuperSpread(pi)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void UpdateRates(self, Py_ssize_t pi, bint infect, bint immune, bint migration):
        cdef double tmp
        if infect:
            self.infectPopRate[pi] = 0.0
            for hn in range(self.currentHapNum):
                self.eventHapPopRate[pi, hn, 0] = self.BirthRate(pi, hn)
                tmp = (self.eventHapPopRate[pi, hn, 0] +
                       self.eventHapPopRate[pi, hn, 1] +
                       self.eventHapPopRate[pi, hn, 2] +
                       self.eventHapPopRate[pi, hn, 3] )
                self.tEventHapPopRate[pi, hn] = tmp
                self.hapPopRate[pi, hn] = self.tEventHapPopRate[pi, hn] * self.infectious[pi, hn]
                self.infectPopRate[pi] += self.hapPopRate[pi, hn]

        if immune:
            self.immunePopRate[pi] = 0
            for sn in range(self.susNum):
                self.immunePopRate[pi] += self.immuneSourcePopRate[pi, sn]

        if infect or immune:
            self.popRate[pi] = self.infectPopRate[pi] + self.immunePopRate[pi]
            self.totalRate = 0.0
            for pn in range(self.popNum):
                self.totalRate += self.popRate[pn]

        if migration:
            self.totalMigrationRate = 0.0
            for pn in range(self.popNum):
                self.migPopRate[pn] = self.maxEffectiveBirthMigration[pn] * self.totalSusceptible[pn] * \
                (self.globalInfectious - self.totalInfectious[pn])
                self.totalMigrationRate += self.migPopRate[pn]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void ImmunityTransition(self, Py_ssize_t pi):
        cdef:
            Py_ssize_t ssi, tsi

        ssi, self.rn = fastChoose(self.immuneSourcePopRate[pi], self.immunePopRate[pi], self.rn)
        tsi, self.rn = fastChoose(self.suscepTransition[ssi], self.suscepCumulTransition[ssi], self.rn)

        self.susceptible[pi, ssi] -= 1
        self.susceptible[pi, tsi] += 1
        self.immuneSourcePopRate[pi, ssi] = self.susceptible[pi, ssi]*self.suscepCumulTransition[ssi]
        self.immuneSourcePopRate[pi, tsi] = self.susceptible[pi, tsi]*self.suscepCumulTransition[tsi]
        self.UpdateRates(pi, False, True, False)

        self.iCounter += 1
        self.events.AddEvent(self.currentTime, SUSCCHANGE, ssi, pi, tsi, 0)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void Birth(self, Py_ssize_t pi, Py_ssize_t hi): # hi - program number
        cdef double ws = 0.0, hs = 0.0

        for sn in range(self.susNum):
            ws += self.susceptHapPopRate[pi, hi, sn]
        si, self.rn = fastChoose(self.susceptHapPopRate[pi, hi], ws, self.rn)

        if self.rn < self.recombination and self.totalInfectious[pi] > 1:
            self.rn = self.rn /  self.recombination

            self.infectious[pi, hi] -= 1
            for hn in range(self.hapNum):
                self.birthInf[hn] = self.eventHapPopRate[pi, hn, 0]*self.infectious[pi, hn]
                hs += self.birthInf[hn]
            self.infectious[pi, hi] += 1
            hi2, self.rn = fastChoose(self.birthInf, hs, self.rn)

            posRecomb = int(self.genome_length * self.rn)
            nhi = 0
            for s in range(self.sites):
                if self.sitesPosition[s] < posRecomb:
                    nhi += 4**(self.sites-s-1)*int(floor(hi/4**(self.sites-s-1)) % 4)
                else:
                    nhi += 4**(self.sites-s-1)*int(floor(hi2/4**(self.sites-s-1)) % 4)

            self.rec.AddRecombination_forward(self.events.ptr, hi, hi2, posRecomb, nhi)

            self.NewInfections(pi, si, nhi)
            self.events.AddEvent(self.currentTime, BIRTH, self.numToHap[hi], pi, si, self.numToHap[hi2])
        else:
            self.NewInfections(pi, si, hi)
            self.events.AddEvent(self.currentTime, BIRTH, self.numToHap[hi], pi, si, 0)
            # self.events.AddEvent(self.currentTime, BIRTH, self.numToHap[hi], pi, si, self.hapNum)

        self.immuneSourcePopRate[pi, si] = self.suscepCumulTransition[si]*self.susceptible[pi, si]
        self.UpdateRates(pi, True, True, True)

        self.bCounter += 1

    def print_recomb(self, left, right):
        for i in range(left, right):
            print('hi(', self.calculate_string(self.rec.his[i]),  ') = ', self.rec.his[i], \
                ', hi2(', self.calculate_string(self.rec.hi2s[i]),  ') = ', self.rec.hi2s[i], \
                ', nhi(', self.calculate_string(self.rec.nhis[i]),  ') = ', self.rec.nhis[i], \
                ', pos = ', self.rec.posRecombs[i], sep='')

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void Death(self, Py_ssize_t pi, Py_ssize_t hi, bint add_event = True): # hi - program number
        self.NewRecoveries(pi, self.suscType[self.numToHap[hi]], hi)
        self.immuneSourcePopRate[pi, self.suscType[self.numToHap[hi]]] = self.susceptible[pi, self.suscType[self.numToHap[hi]]] * \
        self.suscepCumulTransition[self.suscType[self.numToHap[hi]]]
        self.UpdateRates(pi, True, True, True)

        if add_event:
            self.dCounter += 1
            self.events.AddEvent(self.currentTime, DEATH, self.numToHap[hi], pi, self.suscType[self.numToHap[hi]], 0)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void Sampling(self, Py_ssize_t pi, Py_ssize_t hi): # hi - program number
        self.Death(pi, hi, False)

        self.sCounter += 1
        self.events.AddEvent(self.currentTime, SAMPLING, self.numToHap[hi], pi, self.suscType[self.numToHap[hi]], 0)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void Mutation(self, Py_ssize_t pi, Py_ssize_t hi): # hi - program number
        cdef:
            bint check = True
            Py_ssize_t ohi, mi, digit4, AS, DS, nhi

        ohi = self.numToHap[hi] # ohi - haplotype
        mi, self.rn = fastChoose(self.mRate[ohi], self.tmRate[hi], self.rn)
        DS, self.rn = fastChoose(self.hapMutType[ohi, mi], self.hapMutType[ohi, mi, 0] \
            + self.hapMutType[ohi, mi, 1] + self.hapMutType[ohi, mi, 2], self.rn)
        nhi = self.Mutate(ohi, mi, DS)

        if self.memory_optimization:
            #TODO binary search
            for hn in range(self.currentHapNum+1): # hn - program number
                if self.numToHap[hn] == nhi:
                    check = False
                    break
            if check:
                self.AddHaplotype(nhi, pi)
                if nhi < ohi:
                    hi += 1

        self.infectious[pi, self.hapToNum[nhi]] += 1
        self.infectious[pi, hi] -= 1
        self.UpdateRates(pi, True, False, False)

        self.mCounter += 1
        self.events.AddEvent(self.currentTime, MUTATION, ohi, pi, nhi, 0)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef Py_ssize_t SuperSpreadEvent(self):
        cdef:
            bint is_any_inf
            Py_ssize_t ei, pn, pi, sii, sn, hn
            npy_int64 i
            npy_int64[::1] sus_inf, SuperSpreadGroup
            double cumulSusHap
            double[::1] SusHap

        ei, self.rn = fastChoose_vec(self.super_spread_rate, sum(self.super_spread_rate), self.rn)
        pn = self.super_spread_left[ei] + int(self.rn * (self.super_spread_right[ei] - self.super_spread_left[ei]))
        pi = self.super_spread_pop[ei]
        sus_inf = np.zeros(self.susNum + self.hapNum, dtype=np.int64)
        SuperSpreadGroup = np.zeros(self.susNum + self.hapNum, dtype=np.int64)
        for sn in range(self.susNum):
            sus_inf[sn] = self.susceptible[pi, sn]
        for hn in range(self.hapNum):
            sus_inf[hn + self.susNum] = self.infectious[pi, hn]
        for i in range(pn):
            rn = self.seed.uniform()
            sii, rn = fastChoose(sus_inf, self.sizes[pi] - i, rn)
            sus_inf[sii] -= 1
            SuperSpreadGroup[sii] += 1
        is_any_inf = False
        for i in range(self.hapNum):
            if SuperSpreadGroup[i] != 0:
                is_any_inf = True
                break
        if is_any_inf:
            for sn in range(self.susNum):
                SusHap = np.zeros(self.hapNum, dtype=float)
                cumulSusHap = 0
                for hn in range(self.hapNum):
                    SusHap[hn] = self.susceptibility[hn, sn] * self.bRate[hn] * SuperSpreadGroup[hn]
                    cumulSusHap += SusHap[hn]
                for _ in range(self.SuperSpreadGroup[sn]):
                    rn = self.seed.uniform()
                    hi, rn = fastChoose(SusHap, cumulSusHap, rn)
                    self.NewInfections(pi, sn, hi)
                    self.events.AddEvent(self.currentTime, BIRTH, self.numToHap[hi], pi, sn, 0)
        return pi

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef Py_ssize_t GenerateMigration(self):
        cdef:
            Py_ssize_t tpi, spi, hi, si
            double p_accept

        tpi, self.rn = fastChoose(self.migPopRate, self.totalMigrationRate, self.rn)
        spi, self.rn = fastChoose_skip(self.totalInfectious, self.globalInfectious-self.totalInfectious[tpi], self.rn, skip=tpi)
        hi, self.rn = fastChoose(self.infectious[spi], self.totalInfectious[spi], self.rn) # hi - program number
        si, self.rn = fastChoose(self.susceptible[tpi], self.totalSusceptible[tpi], self.rn)

        p_accept = self.effectiveMigration[spi, tpi] * self.bRate[self.numToHap[hi]] * self.susceptibility[self.numToHap[hi], si] / \
        self.maxEffectiveBirthMigration[tpi]
        if self.rn < p_accept:
            self.NewInfections(tpi, si, hi)
            self.UpdateRates(tpi, True, True, True)

            self.migPlus += 1
            self.events.AddEvent(self.currentTime, MIGRATION, self.numToHap[hi], spi, si, tpi)
        else:
            self.migNonPlus += 1

        return tpi

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void CheckLockdown(self, Py_ssize_t pi):
        if self.totalInfectious[pi] > self.startLD[pi]*self.sizes[pi] and self.lockdownON[pi] == 0:
            self.contactDensity[pi] = self.contactDensityAfterLockdown[pi]
            self.swapLockdown += 1
            self.lockdownON[pi] = 1
            self.UpdateAllRates()
            self.loc.AddLockdown(True, pi, self.currentTime)
        if self.totalInfectious[pi] < self.endLD[pi]*self.sizes[pi] and self.lockdownON[pi] == 1:
            self.contactDensity[pi] = self.contactDensityBeforeLockdown[pi]
            self.swapLockdown += 1
            self.lockdownON[pi] = 0
            self.UpdateAllRates()
            self.loc.AddLockdown(False, pi, self.currentTime)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void Restart(self):
        self.events.ptr = 0
        self.multievents.ptr = 0
        self.bCounter = 0
        self.dCounter = 0
        self.sCounter = 0
        self.mCounter = 0
        self.iCounter = 0
        self.migPlus = 0
        self.migNonPlus = 0
        self.currentTime = 0.0
        self.globalInfectious = 0
        for pn in range(self.popNum):
            self.totalSusceptible[pn] = 0
            self.totalInfectious[pn] = 0
            for sn in range(self.susNum):
                self.susceptible[pn, sn] = self.initial_susceptible[pn, sn]
                self.totalSusceptible[pn] += self.initial_susceptible[pn, sn]
            for hn in range(self.hapNum):
                self.infectious[pn, hn] = self.initial_infectious[pn, hn]
                self.totalInfectious[pn] += self.initial_infectious[pn, hn]
                self.globalInfectious += self.initial_infectious[pn, hn]
        for pn in range(self.popNum):
            self.CheckLockdown(pn)
        self.UpdateAllRates()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef GetGenealogy(self, seed):
        cdef:
            Py_ssize_t ptrTreeAndTime, n1, n2, id1, id2, id3, lbs, lbs_e, ns, nt, idt, ids, lbss, propNum
            double p, p1
            vector[vector[vector[Py_ssize_t]]] liveBranchesS, newLineages
            vector[vector[Py_ssize_t]] vecint2
            vector[Py_ssize_t] vecint1
            npy_int64[:,::1] liveBranches

            double e_time, me_time
            Py_ssize_t e_type_, e_population, e_haplotype, e_newHaplotype, e_newPopulation
            Py_ssize_t me_num, me_type_, me_population, me_haplotype, me_newHaplotype, me_newPopulation
            Py_ssize_t mt_ev_num, mt_ev_num2

        propNum = self.popNum * ((self.popNum - 1) * self.hapNum * self.susNum + self.susNum * (self.susNum - 1) + self.hapNum * \
            (2 + self.sites * 3 + self.susNum))

        if self.sCounter < 2:
            print("Less than two cases were sampled...")
            print("_________________________________")
            sys.exit(0)
        else:
            if seed != None:
                self.seed = RndmWrapper(seed=(seed, 0))

            ptrTreeAndTime = 0
            # self.tree = np.zeros(2 * self.sCounter - 1, dtype=np.int64)
            # self.tree_pop = np.zeros(2 * self.sCounter - 1, dtype=np.int64)
            self.tree = np.zeros((2 * self.sCounter - 1, 2), dtype=np.int64)
            self.tree_pop = np.zeros((2 * self.sCounter - 1, 2), dtype=np.int64)
            self.times = np.zeros(2 * self.sCounter - 1, dtype=float)

            for i in range( self.popNum ):
                liveBranchesS.push_back(vecint2)
                for _ in range( self.hapNum ):
                    liveBranchesS[i].push_back(vecint1)

            for i in range( self.popNum ):
                newLineages.push_back(vecint2)
                for _ in range( self.hapNum ):
                    newLineages[i].push_back(vecint1)

            for i in range( self.popNum ):
                for j in range( self.hapNum ):
                    self.infectiousDelta[i, j] = 0

            for e_id in range(self.events.ptr-1, -1, -1):
                # this event
                e_time = self.events.times[e_id]
                e_type_ = self.events.types[e_id]
                e_haplotype = self.events.haplotypes[e_id]
                e_population = self.events.populations[e_id]
                e_newHaplotype = self.events.newHaplotypes[e_id]
                e_newPopulation = self.events.newPopulations[e_id]
                if e_type_ == BIRTH:
                    lbs = liveBranchesS[e_population][e_haplotype].size()
                    lbs_e = self.infectious[e_population, self.hapToNum[e_haplotype]]
                    p = float(lbs)*(float(lbs)-1.0)/ float(lbs_e) / (float(lbs_e) - 1.0)
                    if self.seed.uniform() < p:
                        n1 = int(floor( lbs*self.seed.uniform() ))
                        n2 = int(floor( (lbs-1)*self.seed.uniform() ))
                        if n2 >= n1:
                            n2 += 1
                        id1 = liveBranchesS[e_population][e_haplotype][n1]
                        id2 = liveBranchesS[e_population][e_haplotype][n2]
                        id3 = ptrTreeAndTime
                        liveBranchesS[e_population][e_haplotype][n1] = id3
                        liveBranchesS[e_population][e_haplotype][n2] = liveBranchesS[e_population][e_haplotype][lbs-1]
                        liveBranchesS[e_population][e_haplotype].pop_back()
                        # self.tree[id1] = id3
                        # self.tree[id2] = id3
                        # self.tree[ptrTreeAndTime] = -1
                        # self.tree_pop[ptrTreeAndTime] = e_population
                        self.tree[id1, 0] = id3
                        self.tree[id2, 0] = id3
                        self.tree[ptrTreeAndTime, 0] = -1
                        self.tree_pop[ptrTreeAndTime, 0] = e_population
                        self.times[ptrTreeAndTime] = e_time
                        ptrTreeAndTime += 1
                    self.infectious[e_population, self.hapToNum[e_haplotype]] -= 1

                    #TODO
                    # if e_newPopulation != self.hapNum:
                    #     lbs1 = liveBranchesS[e_population][e_newPopulation].size()
                    #     lbs1_e = self.infectious[e_population, self.hapToNum[e_newPopulation]]
                    #     p = float(lbs1) * (float(lbs1) - 1.0)/ float(lbs1_e) / (float(lbs1_e) - 1.0)
                    #     if self.seed.uniform() < p:
                    #         n1 = int(floor( lbs1 * self.seed.uniform() ))
                    #         n2 = int(floor( (lbs1 - 1) * self.seed.uniform() ))
                    #         if n2 >= n1:
                    #             n2 += 1
                    #         id1 = liveBranchesS[e_population][e_newPopulation][n1]
                    #         id2 = liveBranchesS[e_population][e_newPopulation][n2]
                    #         id3 = ptrTreeAndTime
                    #         liveBranchesS[e_population][e_newPopulation][n1] = id3
                    #         liveBranchesS[e_population][e_newPopulation][n2] = liveBranchesS[e_population][e_newPopulation][lbs-1]
                    #         liveBranchesS[e_population][e_newPopulation].pop_back()
                    #         # self.tree[id1] = id3
                    #         # self.tree[id2] = id3
                    #         # self.tree[ptrTreeAndTime] = -1
                    #         # self.tree_pop[ptrTreeAndTime] = e_population
                    #         self.tree[id1, 1] = id3
                    #         self.tree[id2, 1] = id3
                    #         self.tree[ptrTreeAndTime, 1] = -1
                    #         self.tree_pop[ptrTreeAndTime, 1] = e_population
                    #     # self.rec.AddRecombination(posRecomb)
                    #     self.infectious[e_population, self.hapToNum[e_newPopulation]] -= 1

                elif e_type_ == DEATH:
                    self.infectious[e_population, self.hapToNum[e_haplotype]] += 1
                elif e_type_ == SAMPLING:
                    self.infectious[e_population, self.hapToNum[e_haplotype]] += 1
                    liveBranchesS[e_population][e_haplotype].push_back( ptrTreeAndTime )
                    # self.tree[ptrTreeAndTime] = -1
                    # self.tree_pop[ptrTreeAndTime] = e_population
                    self.tree[ptrTreeAndTime, 0] = -1
                    self.tree_pop[ptrTreeAndTime, 0] = e_population
                    self.times[ptrTreeAndTime] = e_time
                    ptrTreeAndTime += 1
                elif e_type_ == MUTATION:
                    lbs = liveBranchesS[e_population][e_newHaplotype].size()
                    p = float(lbs)/self.infectious[e_population, e_newHaplotype]
                    if self.seed.uniform() < p:
                        n1 = int(floor( lbs*self.seed.uniform() ))
                        id1 = liveBranchesS[e_population][e_newHaplotype][n1]
                        liveBranchesS[e_population][e_newHaplotype][n1] = liveBranchesS[e_population][e_newHaplotype][lbs-1]
                        liveBranchesS[e_population][e_newHaplotype].pop_back()
                        liveBranchesS[e_population][e_haplotype].push_back(id1)
                        self.mut.AddMutation(id1, e_haplotype, e_newHaplotype, e_time)
                    self.infectious[e_population, e_newHaplotype] -= 1
                    self.infectious[e_population, e_haplotype] += 1
                elif e_type_ == SUSCCHANGE:
                    pass
                elif e_type_ == MIGRATION:
                    lbs = liveBranchesS[e_newPopulation][e_haplotype].size()
                    p = float(lbs)/self.infectious[e_newPopulation, e_haplotype]
                    if self.seed.uniform() < p:
                        nt = int(floor( lbs*self.seed.uniform() ))
                        lbss = liveBranchesS[e_population][e_haplotype].size()
                        p1 = float(lbss)/self.infectious[e_population, e_haplotype]
                        if self.seed.uniform() < p1:
                            ns = int(floor( lbss*self.seed.uniform() ))
                            idt = liveBranchesS[e_newPopulation][e_haplotype][nt]
                            ids = liveBranchesS[e_population][e_haplotype][ns]
                            id3 = ptrTreeAndTime
                            liveBranchesS[e_population][e_haplotype][ns] = id3
                            liveBranchesS[e_newPopulation][e_haplotype][nt] = liveBranchesS[e_newPopulation][e_haplotype][lbs-1]
                            liveBranchesS[e_newPopulation][e_haplotype].pop_back()
                            # self.tree[idt] = id3
                            # self.tree[ids] = id3
                            # self.tree[ptrTreeAndTime] = -1
                            # self.tree_pop[ptrTreeAndTime] = e_population
                            self.tree[idt, 0] = id3
                            self.tree[ids, 0] = id3
                            self.tree[ptrTreeAndTime, 0] = -1
                            self.tree_pop[ptrTreeAndTime, 0] = e_population
                            self.times[ptrTreeAndTime] = e_time
                            ptrTreeAndTime += 1
                            self.mig.AddMigration(idt, e_time, e_population, e_newPopulation)
                        else:
                            liveBranchesS[e_population][e_haplotype].push_back(liveBranchesS[e_newPopulation][e_haplotype][nt])
                            liveBranchesS[e_newPopulation][e_haplotype][nt] = liveBranchesS[e_newPopulation][e_haplotype][lbs-1]
                            liveBranchesS[e_newPopulation][e_haplotype].pop_back()
                    self.infectious[e_newPopulation, e_haplotype] -= 1
                elif e_type_ == MULTITYPE:
                    for me_id in range(e_haplotype, e_population):
                        me_num = self.multievents.num[me_id]
                        me_time = self.multievents.times[me_id]
                        me_type_ = self.multievents.types[me_id]
                        me_haplotype = self.multievents.haplotypes[me_id]
                        me_population = self.multievents.populations[me_id]
                        me_newHaplotype = self.multievents.newHaplotypes[me_id]
                        me_newPopulation = self.multievents.newPopulations[me_id]

                        if me_type_ == BIRTH:
                            lbs = liveBranchesS[me_population][me_haplotype].size()
                            lbs_e = self.infectious[me_population, me_haplotype]
                            if me_num == 0 or lbs == 0:
                                mt_ev_num = 0
                            else:
                                mt_ev_num = random_hypergeometric(self.seed.rng, int(lbs*(lbs-1.0)/2.0), int(lbs_e*(lbs_e-1)/2-lbs*(lbs-1)/2), me_num)
                            #print("me_num=", me_num, "  mt_ev_num", mt_ev_num)
                            for i in range(mt_ev_num):
                                n1 = int(floor( lbs*self.seed.uniform() ))
                                n2 = int(floor( (lbs-1)*self.seed.uniform() ))
                                if n2 >= n1:
                                    n2 += 1
                                id1 = liveBranchesS[me_population][me_haplotype][n1]
                                id2 = liveBranchesS[me_population][me_haplotype][n2]
                                id3 = ptrTreeAndTime
                                newLineages[me_population][me_haplotype].push_back(id3)
                                if n1 == lbs-1:#TODO: need to check if we really need these if and elif. Most likely - no!
                                    liveBranchesS[me_population][me_haplotype].pop_back()
                                    liveBranchesS[me_population][me_haplotype][n2] = liveBranchesS[me_population][me_haplotype][lbs-2]
                                    liveBranchesS[me_population][me_haplotype].pop_back()
                                elif n2 == lbs-1:
                                    liveBranchesS[me_population][me_haplotype].pop_back()
                                    liveBranchesS[me_population][me_haplotype][n1] = liveBranchesS[me_population][me_haplotype][lbs-2]
                                    liveBranchesS[me_population][me_haplotype].pop_back()
                                else:
                                    liveBranchesS[me_population][me_haplotype][n1] = liveBranchesS[me_population][me_haplotype][lbs-1]
                                    liveBranchesS[me_population][me_haplotype].pop_back()
                                    liveBranchesS[me_population][me_haplotype][n2] = liveBranchesS[me_population][me_haplotype][lbs-2]
                                    liveBranchesS[me_population][me_haplotype].pop_back()
                                # self.tree[id1] = id3
                                # self.tree[id2] = id3
                                # self.tree[ptrTreeAndTime] = -1
                                # self.tree_pop[ptrTreeAndTime] = me_population
                                self.tree[id1, 0] = id3
                                self.tree[id2, 0] = id3
                                self.tree[ptrTreeAndTime, 0] = -1
                                self.tree_pop[ptrTreeAndTime, 0] = me_population
                                self.times[ptrTreeAndTime] = me_time
                                ptrTreeAndTime += 1
                                lbs -= 2
                            self.infectiousDelta[me_population, me_haplotype] -= me_num
                        elif me_type_ == DEATH:
                            self.infectiousDelta[me_population, me_haplotype] += me_num
                        elif me_type_ == SAMPLING:
                            self.infectiousDelta[me_population, me_haplotype] += me_num
                            for i in range(me_num):
                                #liveBranchesS[me_population][me_haplotype].push_back( ptrTreeAndTime )
                                newLineages[me_population][me_haplotype].push_back( ptrTreeAndTime )
                                # self.tree[ptrTreeAndTime] = -1
                                # self.tree_pop[ptrTreeAndTime] = me_population
                                self.tree[ptrTreeAndTime, 0] = -1
                                self.tree_pop[ptrTreeAndTime, 0] = me_population
                                self.times[ptrTreeAndTime] = me_time
                                ptrTreeAndTime += 1
                        elif me_type_ == MUTATION:
                            lbs = liveBranchesS[me_population][me_newHaplotype].size()
                            if me_num == 0 or lbs == 0:
                                mt_ev_num = 0
                            else:
                                mt_ev_num = random_hypergeometric(self.seed.rng, lbs, self.infectious[me_population, me_newHaplotype]-lbs, me_num)
                            for i in range(mt_ev_num):
                                n1 = int(floor( lbs*self.seed.uniform() ))
                                id1 = liveBranchesS[me_population][me_newHaplotype][n1]
                                liveBranchesS[me_population][me_newHaplotype][n1] = liveBranchesS[me_population][me_newHaplotype][lbs-1]
                                liveBranchesS[me_population][me_newHaplotype].pop_back()
                                #liveBranchesS[me_population][me_haplotype].push_back(id1)
                                newLineages[me_population][me_haplotype].push_back( id1 )
                                self.mut.AddMutation(id1, me_haplotype, me_newHaplotype, me_time)
                                lbs-=1
                            self.infectiousDelta[me_population, me_newHaplotype] -= me_num
                            self.infectiousDelta[me_population, me_haplotype] += me_num
                        elif me_type_ == SUSCCHANGE:
                            pass
                        elif me_type_ == MIGRATION:
                            lbs = liveBranchesS[me_newPopulation][me_haplotype].size()
                            if me_num == 0 or lbs == 0:
                                mt_ev_num = 0
                            else:
                                mt_ev_num = random_hypergeometric(self.seed.rng, lbs, self.infectious[me_newPopulation, me_haplotype]-lbs, me_num)
                            #for i in range(mt_ev_num):
                                lbss = liveBranchesS[me_population][me_haplotype].size()
                                if mt_ev_num == 0 or lbss == 0:
                                    mt_ev_num2 = 0
                                else:
                                    mt_ev_num2 = random_hypergeometric(self.seed.rng, lbss, self.infectious[me_population, me_haplotype]-lbss, mt_ev_num)
                                #print("me_num: ", me_num, "mt_ev_num: ", mt_ev_num, "mt_ev_num2: ", mt_ev_num2)
                                for i in range(mt_ev_num2):
                                    nt = int(floor( lbs*self.seed.uniform() ))
                                    ns = int(floor( lbss*self.seed.uniform() ))
                                    idt = liveBranchesS[me_newPopulation][me_haplotype][nt]
                                    ids = liveBranchesS[me_population][me_haplotype][ns]
                                    id3 = ptrTreeAndTime
                                    liveBranchesS[me_population][me_haplotype][ns] = liveBranchesS[me_population][me_haplotype][  liveBranchesS[me_population][me_haplotype].size()-1 ]
                                    liveBranchesS[me_population][me_haplotype].pop_back()
                                    liveBranchesS[me_newPopulation][me_haplotype][nt] = liveBranchesS[me_newPopulation][me_haplotype][lbs-1]
                                    liveBranchesS[me_newPopulation][me_haplotype].pop_back()
                                    newLineages[me_population][me_haplotype].push_back( id3 )
                                    # self.tree[idt] = id3
                                    # self.tree[ids] = id3
                                    # self.tree[ptrTreeAndTime] = -1
                                    # self.tree_pop[ptrTreeAndTime] = me_population
                                    self.tree[idt, 0] = id3
                                    self.tree[ids, 0] = id3
                                    self.tree[ptrTreeAndTime, 0] = -1
                                    self.tree_pop[ptrTreeAndTime, 0] = me_population
                                    self.times[ptrTreeAndTime] = me_time
                                    ptrTreeAndTime += 1
                                    self.mig.AddMigration(idt, me_time, me_population, me_newPopulation)
                                    lbss -= 1
                                    lbs -= 1
                                for i in range(mt_ev_num-mt_ev_num2):
                                    nt = int(floor( lbs*self.seed.uniform() ))
                                    newLineages[me_population][me_haplotype].push_back(liveBranchesS[me_newPopulation][me_haplotype][nt])
                                    liveBranchesS[me_newPopulation][me_haplotype][nt] = liveBranchesS[me_newPopulation][me_haplotype][lbs-1]
                                    liveBranchesS[me_newPopulation][me_haplotype].pop_back()
                                    lbs -= 1
                            self.infectiousDelta[me_newPopulation, me_haplotype] -= me_num
                        else:
                            print("Unknown event type: ", me_type_)
                            print("_________________________________")
                            sys.exit(0)
                        for pi in range(self.popNum):
                            for hi in range(self.hapNum):
                                self.infectious[pi, hi] += self.infectiousDelta[pi, hi]
                                self.infectiousDelta[pi, hi] = 0
                                while newLineages[pi][hi].size() > 0:
                                    liveBranchesS[pi][hi].push_back( newLineages[pi][hi][newLineages[pi][hi].size()-1] )
                                    newLineages[pi][hi].pop_back()
                                #for i in range(newLineages[pi][hi].size()-1, -1, -1):
                                #    liveBranchesS[pi][hi].push_back(newLineages[pi][hi][i])
                                #    newLineages[pi][hi].pop_back()
                else:
                    print("Unknown event type: ", e_type_)
                    print("_________________________________")
                    sys.exit(0)
            for i in range(self.sCounter * 2 - 2):
                # if self.tree_pop[self.tree[i]] != self.tree_pop[i]:
                #     self.mig.AddMigration(i, self.times[i], self.tree_pop[self.tree[i]], self.tree_pop[i])

                if self.tree_pop[self.tree[i, 0], 0] != self.tree_pop[i, 0]:
                    self.mig.AddMigration(i, self.times[i], self.tree_pop[self.tree[i, 0], 0], self.tree_pop[i, 0])

            #for i in range(self.sCounter * 2 - 1):
            #    print(self.tree[i], end=" ")
            #print("")
            #deg = [0 for i in range(self.sCounter * 2 - 1)]
            #for i in range(self.sCounter * 2 - 1):
            #    deg[self.tree[i]] += 1
            #for i in range(self.sCounter * 2 - 1):
            #    print(deg[i], end=" ")

            #self.CheckTree()


    def print_basic_parameters(self):
        print("*****************")
        print("***Basic rates***")
        print("*****************")
        table = PrettyTable()

        field = ["H", "TR", "RR", "SR", "ST"]
        for s in range(self.sites):
            field.append("M" + str(s))
            field.append("MW" + str(s))
        table.field_names = field
        for hn in range(self.currentHapNum):
            list = ["\n" + self.calculate_string(hn), "\n" + str(self.bRate[hn]), "\n" + str(self.dRate[hn]), "\n" + \
            str(self.sRate[hn]), "\n" + str(self.suscType[hn])]
            for s in range(self.sites):
                list.append("\n" + str(self.mRate[hn, s]))
                list.append(self.create_mutations(hn, s))
            table.add_row(list)

        print(table)
        print("Legend:")
        print("H - haplotype")
        print("TR - transmission rate")
        print("RR - recovery rate")
        print("SR - sampling rate")
        print("ST - susceptibility type")
        for s in range(self.sites):
            print("M" + str(s) + " - " + str(s) + " mutation rate")
            print("MW" + str(s) + " - " + str(s) + " mutation weights")
        print()

    def GetCurrentIndividuals(self):
        current_susceptible = [[0 for _ in range(self.susNum)] for _ in range(self.popNum)]
        current_infectious = [[0 for _ in range(self.hapNum)] for _ in range(self.popNum)]

        for pi in range(self.popNum):
            for si in range(self.susNum):
                current_susceptible[pi][si] = self.initial_susceptible[pi, si]
            for hi in range(self.hapNum):
                current_infectious[pi][hi] = self.initial_infectious[pi, hi]
        for i in range(self.events.ptr):
            if self.events.types[i] == BIRTH:
                current_susceptible[self.events.populations[i]][self.events.newHaplotypes[i]] -= 1
                current_infectious[self.events.populations[i]][self.events.haplotypes[i]] += 1
            elif self.events.types[i] == DEATH:
                current_susceptible[self.events.populations[i]][self.events.newHaplotypes[i]] += 1
                current_infectious[self.events.populations[i]][self.events.haplotypes[i]] -= 1
            elif self.events.types[i] == SAMPLING:
                current_susceptible[self.events.populations[i]][self.events.newHaplotypes[i]] += 1
                current_infectious[self.events.populations[i]][self.events.haplotypes[i]] -= 1
            elif self.events.types[i] == MUTATION:
                current_infectious[self.events.populations[i]][self.events.haplotypes[i]] -= 1
                current_infectious[self.events.populations[i]][self.events.newHaplotypes[i]] += 1
            elif self.events.types[i] == SUSCCHANGE:
                current_susceptible[self.events.populations[i]][self.events.haplotypes[i]] -= 1
                current_susceptible[self.events.populations[i]][self.events.newHaplotypes[i]] += 1
            else:
                current_susceptible[self.events.newPopulations[i]][self.events.newHaplotypes[i]] -= 1
                current_infectious[self.events.newPopulations[i]][self.events.haplotypes[i]] += 1

        return [current_susceptible, current_infectious]

    def print_populations(self):
        current_susceptible, current_infectious = self.GetCurrentIndividuals()

        print("*****************")
        print("***Populations***")
        print("*****************")
        table_populations = PrettyTable()

        table_populations.field_names = ["ID", "Size", 'Actual size', "CD",'CDBLC', "CDALD", "SLD", "ELD", "SM"]
        for pn in range(self.popNum):
            table_populations.add_row([pn, self.sizes[pn], self.actualSizes[pn], self.contactDensity[pn], \
                self.contactDensityBeforeLockdown[pn], self.contactDensityAfterLockdown[pn], self.startLD[pn], self.endLD[pn], \
                self.samplingMultiplier[pn]])

        print(table_populations)
        print("Legend:")
        print("ID - number of population")
        print("Size - size of population")
        print("Actual size - actual size of population")
        print("CD - contact density")
        print("CDBLD - contact density without lockdown")
        print("CDALD - contact density at lockdown")
        print("SLD - start of lockdown")
        print("ELD - end of lockdown")
        print("SM - sampling multiplier")
        print()

        print("*****************")
        print("***Susceptible***")
        print("*****************")
        table_susceptible = PrettyTable()

        field = ["ST\\ID"]
        for pn in range(self.popNum):
            field.append(pn)
        for sn in range(self.susNum):
            row = [sn]
            for pn in range(self.popNum):
                row.append(current_susceptible[pn][sn])
            table_susceptible.add_row(row)
        table_susceptible.field_names = field

        print(table_susceptible)
        print("Legend:")
        print("ID - ID population")
        print("ST - susceptibility type")
        print()

        print("****************")
        print("***Infectious***")
        print("****************")
        table_infectious = PrettyTable()

        field = ["H\\ID"]
        for pn in range(self.popNum):
            field.append(pn)
        for hn in range(self.hapNum):
            row = [self.calculate_string(hn)]
            for pn in range(self.popNum):
                row.append(current_infectious[pn][hn])
            table_infectious.add_row(row)
        table_infectious.field_names = field

        print(table_infectious)
        print("Legend:")
        print("ID - ID population")
        print("H - haplotype")
        print()

        print("**********************")
        print("***Migration matrix***")
        print("**********************")
        table_migration = PrettyTable()

        field = ["S\\T"]
        for pn1 in range(self.popNum):
            field.append(pn1)
            row = [pn1]
            for pn2 in range(self.popNum):
                row.append(self.migrationRates[pn1, pn2])
            table_migration.add_row(row)
        table_migration.field_names = field

        print(table_migration)
        print("Legend:")
        print("S - ID source population")
        print("T - ID target population")
        print()

    def print_immunity_model(self):
        print("********************")
        print("***Immunity model***")
        print("********************")
        table_immunity = PrettyTable()

        field = ["H\\ST"]
        for sn in range(self.susNum):
            field.append("S" + str(sn))
        table_immunity.field_names = field
        for hn in range(self.currentHapNum):
            row = [self.calculate_string(hn)]
            for sn in range(self.susNum):
                row.append(self.susceptibility[hn, sn])
            table_immunity.add_row(row)

        print(table_immunity)
        print("Legend:")
        print("H - haplotype")
        print("ST - susceptibility type")
        print()

        print("*******************************")
        print("***Immunity transition rates***")
        print("*******************************")
        table_immunity_transition = PrettyTable()

        field = ["ID"]
        for sn1 in range(self.susNum):
            field.append(sn1)
            row = [sn1]
            for sn2 in range(self.susNum):
                row.append(self.suscepTransition[sn1, sn2])
            table_immunity_transition.add_row(row)
        table_immunity_transition.field_names = field

        print(table_immunity_transition)
        print("Legend:")
        print("ID - ID susceptibility type")
        print()

    def print_mutations(self):
        print('nodeId\tDS\tAS\tsite\ttime')
        for i in range(self.mut.nodeId.size()):
            print(self.mut.get_mutation(i))

    def print_migrations(self):
        print('nodeId\ttime\tsource population\ttarget population')
        for i in range(self.mig.nodeId.size()):
            print(self.mig.get_migration(i))


    def create_list_for_cycles(self, index, edge):
        if isinstance(index, str):
            haplotypes = [index]
            for s in range(self.sites):
                for i in range(len(haplotypes)):
                    haplotype_old = haplotypes[i]
                    if haplotype_old[s] == "*":
                        index = haplotype_old.replace("*", "A", 1)
                        haplotypes.append(index)
                        index = haplotype_old.replace("*", "T", 1)
                        haplotypes.append(index)
                        index = haplotype_old.replace("*", "C", 1)
                        haplotypes.append(index)
                        index = haplotype_old.replace("*", "G", 1)
                        haplotypes.append(index)
            for i in range(len(haplotypes)-1, -1, -1):
                if haplotypes[i].count("*") != 0:
                    haplotypes.remove(haplotypes[i])
            for i in range(len(haplotypes)):
                haplotypes[i] = self.calculate_haplotype(haplotypes[i])

            return haplotypes
        elif isinstance(index, int):
            return [index]
        else:
            return range(edge)


    def create_mutations(self, haplotype, site):
        hap = self.calculate_string(haplotype)
        haplotypes = [hap[:site] + "A" + hap[site+1:], hap[:site] + "T" + hap[site+1:], hap[:site] + "C" + hap[site+1:], \
        hap[:site] + "G" + hap[site+1:]]
        for i in range(4):
            if haplotypes[i] == hap:
                haplotypes.remove(haplotypes[i])
                haplotypes.append(hap)
                break
        color_hap=[]
        for hapl in haplotypes:
            a = ""
            for s in range(self.sites):
                if s == site:
                    a = a + "\033[31m{}\033[0m" .format(hapl[s:s+1])
                else:
                    a = a + hapl[s:s+1]
            color_hap.append(a)
        string = color_hap[3] + "->" + color_hap[0] + ": " + str(self.hapMutType[haplotype, site, 0]) + "\n" + color_hap[3] + \
        "->" + color_hap[1] + ": " + str(self.hapMutType[haplotype, site, 1]) + "\n" + color_hap[3] + "->" + color_hap[2] + ": " + \
        str(self.hapMutType[haplotype, site, 2]) + "\n"
        return string

    def calculate_string(self, hapNum):
        letters = ["A", "T", "C", "G"]
        string = ""
        for s in range(self.sites):
            string = letters[hapNum%4] + string
            hapNum = hapNum // 4
        return string

    def calculate_haplotype(self, string):
        string = string[::-1]
        haplotype = 0
        for s in range(self.sites):
            if string[s]=="T":
                haplotype += (4**s)
            elif string[s]=="C":
                haplotype += 2*(4**s)
            elif string[s]=="G":
                haplotype += 3*(4**s)
        return haplotype

    def calculate_allele(self, haplotype, site):
        for _ in range(self.sites-site):
            allele = haplotype % 4
            haplotype = haplotype // 4
        return allele

    def general_sampling(self):
        for pi in range(self.popNum):
            sampling_size = int(self.sizes[pi] * self.sampling_proportion)
            for i in range(sampling_size):
                self.rn = self.seed.uniform()
                hi, self.rn = fastChoose(self.hapPopRate[pi], self.infectPopRate[pi], self.rn)
                self.Sampling(pi, hi)

    @property
    def seed(self):
        return self.user_seed

    @property
    def sampling_probability(self):
        return self.sampling_probability

    @property
    def memory_optimization(self):
        return self.memory_optimization

    @property
    def number_of_sites(self):
        return self.sites

    @property
    def haplotype_number(self):
        return self.hapNum

    @property
    def populations_number(self):
        return self.popNum

    @property
    def number_of_susceptible_groups(self):
        return self.susNum


    def check_amount(self, amount, smth, zero=True):
        if isinstance(amount, int) == False:
            raise TypeError('Incorrect type of ' + smth + '. Type should be int.')
        elif amount <= 0 and zero:
            raise ValueError('Incorrect value of ' + smth + '. Value should be more 0.')
        elif amount < 0 and zero == False:
            raise ValueError('Incorrect value of ' + smth + '. Value should be more or equal 0.')

    def check_value(self, value, smth, edge=None, none=False):
        if none:
            if isinstance(value, (int, float)) == False and value != None:
                raise TypeError('Incorrect type of ' + smth + '. Type should be int or float or None.')
        else:
            if isinstance(value, (int, float)) == False:
                raise TypeError('Incorrect type of ' + smth + '. Type should be int or float.')
        if edge == None:
            if isinstance(value, (int, float)):
                if value < 0:
                    raise ValueError('Incorrect value of ' + smth + '. Value should be more or equal 0.')
        else:
            if isinstance(value, (int, float)):
                if value < 0 or value > edge:
                    raise ValueError('Incorrect value of ' + smth + '. Value should be more or equal 0 and equal or less ' + str(edge) + '.')

    def check_index(self, index, edge, smth, hap=False, none=True):
        if none == False and index == None:
            raise TypeError('Incorrect type of ' + smth + '. Type should be int.')
        elif isinstance(index, int):
            if index < 0 or index >= edge:
                raise IndexError('There are no such ' + smth + '!')
        elif isinstance(index, str) and hap:
            if index.count("A") + index.count("T") + index.count("C") + index.count("G") + index.count("*") \
            != self.sites:
                raise ValueError('Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"*\" and length of haplotype should be equal number of mutations sites.')
        elif index != None:
            if hap:
                raise TypeError('Incorrect type of haplotype. Type should be int or str or None.')
            else:
                raise TypeError('Incorrect type of ' + smth + '. Type should be int or None.')

    def check_list(self, data, smth, length):
        if isinstance(data, list):
            if len(data) != length:
                raise ValueError('Incorrect length of ' + smth + '. Length should be equal ' + str(length) + '.')
        else:
            raise TypeError('Incorrect type of ' + smth + '. Type should be list.')

    def check_amount_sus(self, amount, source_type, target_type, population):
        if self.initial_susceptible[population, source_type] - amount < 0:
            raise ValueError('Number of susceptible minus amount should be more or equal 0.')
        if self.initial_susceptible[population, target_type] + amount > self.sizes[population]:
            raise ValueError('Number of susceptible plus amount should be equal or less population size.')

    def check_amount_inf(self, amount, source_type, target_haplotype, population):
        if self.initial_susceptible[population, source_type] - amount < 0:
            raise ValueError('Number of susceptible minus amount should be more or equal 0.')
        if self.initial_infectious[population, target_haplotype] + amount > self.sizes[population]:
            raise ValueError('Number of infectious plus amount should be equal or less population size.')

    def check_mig_rate(self):
        for pn1 in range(self.popNum):
            summa = 0
            self.migrationRates[pn1, pn1] = 1.0
            for pn2 in range(self.popNum):
                if pn1 != pn2:
                    summa += self.migrationRates[pn1, pn2]
                    self.migrationRates[pn1, pn1] -= self.migrationRates[pn1, pn2]
            if summa > 1:
                raise ValueError('Incorrect the sum of migration probabilities. The sum of migration probabilities from each population should be equal or less 1.')
        for pn in range(self.popNum):
            if self.migrationRates[pn, pn] <= 1e-15:
                raise ValueError('Incorrect value of migration probability. Value of migration probability from source population to target population should be more 0.')
            
 

    @property
    def initial_haplotype(self):
        return self.maxHapNum

    def set_initial_haplotype(self, amount):
        if self.memory_optimization == False:
            raise ValueError('Incorrect value of memory optimization. Value should be equal \'True\' for work this function.')
        self.check_amount(amount, 'amount of initial haplotype')

        if amount >= self.hapNum:
            self.maxHapNum = self.hapNum
        else:
            self.maxHapNum = amount

    @property
    def step_haplotype(self):
        return self.addMemoryNum

    def set_step_haplotype(self, amount):
        if self.memory_optimization == False:
            raise ValueError('Incorrect value of memory optimization. Value should be equal \'True\' for work this function.')
        self.check_amount(amount, 'amount of step haplotype')

        self.addMemoryNum = amount

    @property
    def genome_length(self):
        return self.genome_length

    def set_genome_length(self, genome_length):
        self.check_amount(genome_length, 'genome length')
        if self.sites > genome_length:
            raise ValueError('Incorrect value of number of sites or genome length. Genome length should be more or equal number of sites.')

        self.genome_length = genome_length
        for s in range(self.sites):
            self.sitesPosition[s] = int(s * self.genome_length / (self.sites - 1))

    @property
    def coinfection_parameters(self):
        return self.recombination

    def set_coinfection_parameters(self, recombination):
        self.check_value(recombination, 'recombination probability', edge=1)

        self.recombination = recombination

    @property
    def super_spread_rate(self):
        return [np.asarray(<double [:self.super_spread_rate.size()]>self.super_spread_rate.data()),
        np.asarray(<Py_ssize_t [:self.super_spread_left.size()]>self.super_spread_left.data()),
        np.asarray(<Py_ssize_t [:self.super_spread_right.size()]>self.super_spread_right.data()),
        np.asarray(<Py_ssize_t [:self.super_spread_pop.size()]>self.super_spread_pop.data())]

    def set_super_spread_rate(self, rate,left,right, population):
        self.check_value(rate, 'super spread rate')
        self.check_index(population, self.popNum, 'population')
        self.check_amount(left,'left point')
        self.check_amount(right, 'right point')
        if right < left:
            raise ValueError('Incorrect value of max value. Value should be more then min value.')
        populations = self.create_list_for_cycles(population, self.popNum)
        for pn in populations:
            if right > self.sizes[pn]:
                raise ValueError('Incorrect value of max value. Value should be less then population size value.')

        for pn in populations:
            self.SuperSpreadRate += rate
            self.super_spread_rate.push_back(rate)
            self.super_spread_left.push_back(left)
            self.super_spread_right.push_back(right)
            self.super_spread_pop.push_back(pn)

    @property
    def transmission_rate(self):
        return self.bRate

    def set_transmission_rate(self, rate, haplotype):
        self.check_value(rate, 'transmission rate')
        self.check_index(haplotype, self.hapNum, 'haplotype', True)

        haplotypes = self.create_list_for_cycles(haplotype, self.hapNum)
        for hn in haplotypes:
            self.bRate[hn] = rate

    # def set_transmission_rate(self, rate, haplotype, condition):
    #     self.check_value(rate, 'transmission rate')
    #     haplotypes = self.check_index_hap(haplotype)
    #     conditions = self.check_index(condition, self.conNum, 'number of states')

    #     for hn in haplotypes:
    #         for cn in conditions:
    #             self.bRate[hn, cn] = rate

    @property
    def recovery_rate(self):
        return self.dRate

    def set_recovery_rate(self, rate, haplotype):
        self.check_value(rate, 'recovery rate')
        self.check_index(haplotype, self.hapNum, 'haplotype', True)

        haplotypes = self.create_list_for_cycles(haplotype, self.hapNum)
        for hn in haplotypes:
            self.dRate[hn] = rate

    @property
    def sampling_rate(self):
        return self.sRate

    def set_sampling_rate(self, rate, haplotype):
        self.check_index(haplotype, self.hapNum, 'haplotype', True)

        if self.sampling_probability == True:
            self.check_value(rate, 'sampling probability', edge=1)
            
            haplotypes = self.create_list_for_cycles(haplotype, self.hapNum)
            for hn in haplotypes:
                deathRate = self.dRate[hn] + self.sRate[hn]
                self.dRate[hn] = (1 - rate) * deathRate
                self.sRate[hn] = rate * deathRate

        elif self.sampling_probability == False:
            self.check_value(rate, 'sampling rate')

            haplotypes = self.create_list_for_cycles(haplotype, self.hapNum)
            for hn in haplotypes:
                self.sRate[hn] = rate

    @property
    def mutation_rate(self):
        return self.mRate

    def set_mutation_rate(self, rate, haplotype, mutation):
        self.check_value(rate, 'mutation rate')
        self.check_index(haplotype, self.hapNum, 'haplotype', True)
        self.check_index(mutation, self.sites, 'mutation site')
        
        haplotypes = self.create_list_for_cycles(haplotype, self.hapNum)
        sites = self.create_list_for_cycles(mutation, self.sites)
        for hn in haplotypes:
            for s in sites:
                self.mRate[hn, s] = rate  

    @property
    def mutation_probabilities(self):
        return self.hapMutType 

    def set_mutation_probabilities(self, probabilities, haplotype, mutation):
        self.check_list(probabilities, 'probabilities list', 4)
        if isinstance(probabilities, list):
            for i in range(4):
                self.check_value(probabilities[i], 'mutation probabilities')
        self.check_index(haplotype, self.hapNum, 'haplotype', True)
        self.check_index(mutation, self.sites, 'mutation site')

        haplotypes = self.create_list_for_cycles(haplotype, self.hapNum)
        sites = self.create_list_for_cycles(mutation, self.sites)
        for hn in haplotypes:
            for s in sites:
                probabilities_allele = list(probabilities)
                del probabilities_allele[self.calculate_allele(hn, s)]
                if sum(probabilities_allele) == 0:
                    raise ValueError('Incorrect probabilities list. The sum of three elements without mutation allele should be more 0.')
                self.hapMutType[hn, s, 0] = probabilities_allele[0]
                self.hapMutType[hn, s, 1] = probabilities_allele[1]
                self.hapMutType[hn, s, 2] = probabilities_allele[2]

    @property
    def mutation_position(self):
        return self.sitesPosition

    def set_mutation_position(self, mutation, position):
        self.check_index(mutation, self.sites, 'number of site', none=False)
        self.check_index(position, self.genome_length, 'mutation position', none=False)
        for s in range(self.sites):
            if self.sitesPosition[s] == position and s != mutation:
                raise IndexError('Incorrect value of position. Two mutations can\'t have the same position.')

        self.sitesPosition[mutation] = position


    @property
    def susceptibility_type(self):
        return self.suscType

    def set_susceptibility_type(self, susceptibility_type, haplotype):
        if isinstance(susceptibility_type, int) == False:
            raise TypeError('Incorrect type of susceptibility type. Type should be int.')
        elif susceptibility_type < 0 or susceptibility_type >= self.susNum:
            raise IndexError('There are no such susceptibility type!')
        self.check_index(haplotype, self.hapNum, 'haplotype', True)

        haplotypes = self.create_list_for_cycles(haplotype, self.hapNum)
        for hn in haplotypes:
            self.suscType[hn] = susceptibility_type

    @property
    def susceptibility(self):
        return self.susceptibility

    def set_susceptibility(self, rate, haplotype, susceptibility_type):
        self.check_value(rate, 'susceptibility rate')
        self.check_index(haplotype, self.hapNum, 'haplotype', True)
        self.check_index(susceptibility_type, self.susNum, 'susceptibility type')

        haplotypes = self.create_list_for_cycles(haplotype, self.hapNum)
        sus_types = self.create_list_for_cycles(susceptibility_type, self.susNum)
        for hn in haplotypes:
            for sn in sus_types:
                self.susceptibility[hn, sn] = rate

    @property
    def immunity_transition(self):
        return self.suscepTransition

    def set_immunity_transition(self, rate, source, target):
        self.check_value(rate, 'immunity transition rate')
        self.check_index(source, self.susNum, 'susceptibility type')
        self.check_index(target, self.susNum, 'susceptibility type')

        sus_types_1 = self.create_list_for_cycles(source, self.susNum)
        sus_types_2 = self.create_list_for_cycles(target, self.susNum)
        for sn1 in sus_types_1:
            for sn2 in sus_types_2:
                if sn1 != sn2:
                    self.suscepTransition[sn1, sn2] = rate

    @property
    def population_size(self):
        return self.sizes

    def set_population_size(self, amount, population):
        if self.first_simulation == True:
            raise ValueError('Changing population size is available only before first simulation!')
        self.check_amount(amount, 'population size')
        self.check_index(population, self.popNum, 'population')

        populations = self.create_list_for_cycles(population, self.popNum)
        for pn in populations:
            self.sizes[pn] = amount
            self.initial_susceptible[pn, 0] = amount
            for sn in range(1, self.susNum):
                self.initial_susceptible[pn, 0] = 0

    @property
    def susceptible(self):
        return self.initial_susceptible

    def set_susceptible(self, amount, source_type, target_type, population):
        if self.first_simulation:
            raise ValueError('This function is available only before first simulation!')
        self.check_amount(amount)
        self.check_index(source_type, self.susNum, 'susceptibility type')
        self.check_index(target_type, self.susNum, 'susceptibility type')
        if source_type == target_type:
            raise ValueError('Source and target susceptibility type shouldn\'t be equal!')
        self.check_index(population, self.popNum, 'population')

        populations = self.create_list_for_cycles(population, self.popNum)
        for pn in populations:
            self.check_amount_sus(amount, source_type, target_type, pn)

            self.initial_susceptible[pn, source_type] -= amount
            self.initial_susceptible[pn, target_type] += amount

    @property
    def infectious(self):
        return self.initial_infectious

    def set_infectious(self, amount, source_type, target_haplotype, population):
        if self.first_simulation:
            raise ValueError('This function is available only before first simulation!')
        self.check_amount(amount)
        self.check_index(source_type, self.susNum, 'susceptibility type')
        self.check_index(target_haplotype, self.hapNum, 'haplotype')
        self.check_index(population, self.popNum, 'population')

        populations = self.create_list_for_cycles(population, self.popNum)
        for pn in populations:
            self.check_amount_inf(amount, source_type, target_haplotype, pn)

            self.initial_susceptible[pn, source_type] -= amount
            self.initial_infectious[pn, target_haplotype] += amount

    @property
    def contact_density(self):
        return self.contactDensity

    def set_contact_density(self, value, population):
        self.check_value(value, 'contact density')
        self.check_index(population, self.popNum, 'population')

        populations = self.create_list_for_cycles(population, self.popNum)
        for pn in populations:
            self.contactDensity[pn] = value
            self.contactDensityBeforeLockdown[pn] = value

    @property
    def npi(self):
        return [self.contactDensityAfterLockdown, self.startLD, self.endLD]

    def set_npi(self, parameters, population):
        self.check_list(parameters, 'npi parameters', 3)
        if isinstance(parameters, list):
            self.check_value(parameters[0], 'first npi parameter')
            self.check_value(parameters[1], 'second npi parameter', edge=1)
            self.check_value(parameters[2], 'third npi parameter', edge=1)
        self.check_index(population, self.popNum, 'population')

        populations = self.create_list_for_cycles(population, self.popNum)
        for pn in populations:
            self.contactDensityAfterLockdown[pn] = parameters[0]
            self.startLD[pn] = parameters[1]
            self.endLD[pn] = parameters[2]

    @property
    def sampling_multiplier(self):
        return self.samplingMultiplier

    def set_sampling_multiplier(self, multiplier, population):
        self.check_value(multiplier, 'sampling multiplier')
        self.check_index(population, self.popNum, 'population')

        populations = self.create_list_for_cycles(population, self.popNum)
        for pn in populations:
            self.samplingMultiplier[pn] = multiplier

    @property
    def migration_probability(self):
        return self.migrationRates

    def set_migration_probability(self, probability, source, target):
        self.check_value(probability, 'migration probability', edge=1)
        self.check_index(source, self.popNum, 'population')
        self.check_index(target, self.popNum, 'population')

        populations_1 = self.create_list_for_cycles(source, self.popNum)
        populations_2 = self.create_list_for_cycles(target, self.popNum)
        for pn1 in populations_1:
            for pn2 in populations_2:
                if pn1 != pn2:
                    self.migrationRates[pn1, pn2] = probability

        self.check_mig_rate()

    def set_total_migration_probability(self, total_probability):
        self.check_value(total_probability, 'total migration probability', edge=1)

        source_rate = 1.0 - total_probability
        target_rate = total_probability/(self.popNum-1)
        for pn1 in range(self.popNum):
            for pn2 in range(self.popNum):
                if pn1 != pn2:
                    self.migrationRates[pn1, pn2] = target_rate
                else:
                    self.migrationRates[pn1, pn2] = source_rate

        self.check_mig_rate()


    def set_chain_events(self, name_file):
        if isinstance(name_file, str) == False:
            pass
            #TODO
        tokens = np.load(name_file + '.npy')
        #TODO Check sizes and resize if needed
        self.size = tokens.shape[1]
        self.ptr = tokens.shape[1]

        self.events.times = tokens[0]
        self.events.types = tokens[1]
        self.events.haplotypes = tokens[2]
        self.events.populations = tokens[3]
        self.events.newHaplotypes = tokens[4]
        self.events.newPopulations = tokens[5]

    def set_settings(self, file_template):
        pass
        pass

    def set_general_sampling(self, sampling_proportion, sampling_events_number, sampling_times):
        self.check_value(sampling_proportion, "general sampling proportion")
        self.check_index(sampling_proportion, 1, "general sampling proportion")
        self.check_value(sampling_events_number, "general sampling events number")
        for sampling_time in sampling_times:
            self.check_value(sampling_time, "general sampling time")
        self.sampling_proportion = sampling_proportion
        self.sampling_events_number = sampling_events_number
        self.sampling_times = sampling_times


    def output_tree_mutations(self):
        if self.tree.shape[0] == 1:
            print('Genealogy was not simulated. Use VGsim.genealogy() method to simulate it.')
            sys.exit(1)
        mut = [[], [], [], [], []]
        for i in range(self.mut.nodeId.size()):
            mut[0].append(self.mut.nodeId[i])
            mut[1].append(self.mut.AS[i])
            mut[2].append(self.mut.site[i])
            mut[3].append(self.mut.DS[i])
            mut[4].append(self.mut.time[i])
        times_dict = {self.events.times[i]: i for i in range(len(self.events.times))}
        populations = {}
        for time in self.times:
            populations[time] = self.events.populations[times_dict[time]]

        return self.tree, self.times, mut, populations

    def export_migrations(self, name_file, file_path):
        if self.tree.shape[0] == 1:
            print('Genealogy was not simulated. Use VGsim.genealogy() method to simulate it.')
            sys.exit(1)
        if file_path != None:
            f_mig = open(file_path + '/' + name_file + '.tsv', 'w')
        else:
            f_mig = open(name_file + '.tsv', 'w')
        f_mig.write("Node\tTime\tOld_population\tNew_population\n")
        for i in range(self.mig.nodeId.size()):
            f_mig.write(str(self.mig.nodeId[i]) + '\t' + str(self.mig.time[i]) + '\t' + str(self.mig.oldPop[i]) + '\t' + str(self.mig.newPop[i]) + "\n")
        f_mig.close()

    def output_sample_data(self):
        time, pop, hap = [], [], []
        for i in range(self.events.ptr):
            if self.events.types[i] == SAMPLING:
                time.append(self.events.times[i])
                pop.append(self.events.populations[i])
                hap.append(self.events.haplotypes[i])
        return time, pop, hap

    def output_epidemiology_timelines(self, step_num, output_file):
        time_points = [i*self.currentTime/step_num for i in range(step_num+1)]
        suscepDate = np.zeros((self.popNum, self.susceptible_num), dtype=np.int64)
        hapDate = np.zeros((self.popNum, self.hapNum), dtype=np.int64)
        for i in range(self.popNum):
            for sn in range(self.susceptible_num):
                if sn == 0:
                    suscepDate[i, sn] = self.sizes[i]
                else:
                    suscepDate[i, sn] = 0
            for hn in range(self.hapNum):
                hapDate[i, hn] = 0
        hapDate[0, 0] += 1
        suscepDate[0, 0] -= 1
        if output_file == True:
            if not os.path.isdir("logs"):
                os.mkdir("logs")
            logDynamics = []
            for i in range(self.popNum):
                logDynamics.append(open('logs/PID' + str(i) + '.log', 'w'))
                logDynamics[i].write("time")
                for sn in range(self.susceptible_num):
                    logDynamics[i].write(" S" + str(sn))
                for hn in range(self.hapNum):
                    logDynamics[i].write(" H" + str(hn))
                logDynamics[i].write("\n")
        else:
            log = dict()
            log["time"] = list()
            for i in range(self.popNum):
                log["P" + str(i)] = dict()
                for j  in range(self.susceptible_num):
                    log["P" + str(i)]["S" + str(j)] = list()
                for j  in range(self.hapNum):
                    log["P" + str(i)]["H" + str(j)] = list()

        point = 0
        for j in range(self.events.ptr):
            if self.events.types[j] == BIRTH:
                hapDate[self.events.populations[j], self.events.haplotypes[j]] += 1
                suscepDate[self.events.populations[j], self.events.newHaplotypes[j]] -= 1
            elif self.events.types[j] == DEATH:
                hapDate[self.events.populations[j], self.events.haplotypes[j]] -= 1
                suscepDate[self.events.populations[j], self.events.newHaplotypes[j]] += 1
            elif self.events.types[j] == SAMPLING:
                hapDate[self.events.populations[j], self.events.haplotypes[j]] -= 1
                suscepDate[self.events.populations[j], self.events.newHaplotypes[j]] += 1
            elif self.events.types[j] == MUTATION:
                hapDate[self.events.populations[j], self.events.haplotypes[j]] -= 1
                hapDate[self.events.populations[j], self.events.newHaplotypes[j]] += 1
            elif self.events.types[j] == SUSCCHANGE:
                suscepDate[self.events.populations[j], self.events.haplotypes[j]] -= 1
                suscepDate[self.events.populations[j], self.events.newHaplotypes[j]] += 1
            elif self.events.types[j] == MIGRATION:
                suscepDate[self.events.newPopulations[j], self.events.newHaplotypes[j]] -= 1
                hapDate[self.events.newPopulations[j], self.events.haplotypes[j]] += 1
            elif self.events.types[j] == MULTITYPE:
                pass
                #TODO
            if time_points[point] <= self.events.times[j]:
                if output_file == True:
                    for i in range(self.popNum):
                        logDynamics[i].write(str(time_points[point]) + " ")
                        for k in range(self.susceptible_num):
                            logDynamics[i].write(str(suscepDate[i, k]) + " ")
                        for k in range(self.hapNum):
                            logDynamics[i].write(str(hapDate[i, k]) + " ")
                        logDynamics[i].write("\n")
                    point += 1
                else:
                    log["time"].append(time_points[point])
                    for i in range(self.popNum):
                        for j  in range(self.susceptible_num):
                            log["P" + str(i)]["S" + str(j)].append(suscepDate[i, j])
                        for j  in range(self.hapNum):
                            log["P" + str(i)]["H" + str(j)].append(hapDate[i, j])
                    point += 1

        if output_file == True:
            for i in range(self.popNum-1, -1, -1):
                logDynamics[i].close()
        else:
            return log

    def export_chain_events(self, name_file):
        np.save(name_file, [self.events.times, self.events.types, self.events.haplotypes, self.events.populations, \
            self.events.newHaplotypes, self.events.newPopulations])

    def export_settings(self, file_template):
        if isinstance(file_template, str) == False:
            pass
            #TODO
        if not os.path.isdir(file_template):
            os.mkdir(file_template)
        os.chdir(file_template)
        comand = 'Command line command: '
        with open(file_template + ".rt", "w") as file:
            file.write("#Rates_format_version 0.0.1\nH B D S")
            for s in range(self.sites):
                file.write(" M" + str(s))
            file.write("\n")
            for hn in range(self.hapNum):
                file.write(self.calculate_string(hn) + " " + str(self.bRate[hn]) + " " + str(self.dRate[hn]) + " " + str(self.sRate[hn]) + ' ')
                for s in range(self.sites):
                    file.write(str(self.mRate[hn, s]) + "," + str(self.hapMutType[hn, s, 0]) + "," + str(self.hapMutType[hn, s, 1]) + "," + str(self.hapMutType[hn, s, 2]) + " ")
                file.write("\n")
            comand += (file_template + '/' + file_template + '.rt ')

        with open(file_template + ".pp", "w") as file:
            file.write("#Population_format_version 0.0.1\nid size contactDensity conDenAfterLD startLD endLD samplingMulriplier\n")
            for pn in range(self.popNum):
                file.write(str(pn) + " " + str(self.sizes[pn]) + " " + str(self.contactDensity[pn]) + " " + str(self.contactDensityAfterLockdown[pn]) + "," + str(self.startLD[pn]) + "," + str(self.endLD[pn]) + " " + str(self.samplingMultiplier[pn]) + "\n")
            comand += ('-pm ' + file_template + '/' + file_template + '.pp ')

        with open(file_template + ".mg", "w") as file:
            file.write("#Migration_format_version 0.0.1\n")
            for pn1 in range(self.popNum):
                for pn2 in range(self.popNum):
                    file.write(str(self.migrationRates[pn1, pn2]) + " ")
                file.write("\n")
            comand += (file_template + '/' + file_template + '.mg ')

        with open(file_template + ".su", "w") as file:
            file.write("#Susceptibility_format_version 0.0.1\nH T")
            for sn in range(self.susNum):
                file.write(" S" + str(sn))
            file.write("\n")
            for hn in range(self.hapNum):
                file.write(self.calculate_string(hn) + " " + str(self.suscType[hn]))
                for sn in range(self.susNum):
                    file.write(" " + str(self.susceptibility[hn, sn]))
                file.write("\n")
            comand += ('-su ' + file_template + '/' + file_template + '.su ')

        with open(file_template + ".st", "w") as file:
            file.write("#Susceptibility_format_version 0.0.1\n")
            for sn1 in range(self.susNum):
                for sn2 in range(self.susNum):
                    file.write(str(self.suscepTransition[sn1, sn2]) + " ")
                file.write("\n")
            comand += ('-st ' + file_template + '/' + file_template + '.st ')
        print(comand)
        os.chdir('../')

    def export_ts(self):
        tc = tskit.TableCollection()
        tc.sequence_length = self.genome_length

        for i in range(self.mig.nodeId.size()):
            mig_nodeId, mig_time, mig_oldPop, mig_newPop = self.mig.get_migration(i)
            tc.migrations.add_row(0.0, 1.0, mig_nodeId, mig_oldPop, mig_newPop, self.times[0] - mig_time)

        for i in range(self.popNum):
            tc.populations.add_row(None)

        for i in range(2 * self.sCounter - 2):
            tc.edges.add_row(0.0, self.genome_length, self.tree[i, 0], i)

        child_or_parent = [1 for _ in range(2 * self.sCounter - 1)]
        for i in range(2 * self.sCounter - 2):
            child_or_parent[self.tree[i, 0]] = 0

        for i in range(2 * self.sCounter - 1):
            tc.nodes.add_row(child_or_parent[i], self.times[0] - self.times[i], self.tree_pop[i, 0])

        for i in range(self.sites):
            if self.sitesPosition[i] == 0:
                tc.sites.add_row(self.sitesPosition[i] + 1, 'A')
            elif self.sitesPosition[i] == self.genome_length:
                tc.sites.add_row(self.sitesPosition[i] - 1, 'A')
            else:
                tc.sites.add_row(self.sitesPosition[i], 'A')

        allele = ['A', 'T', 'C', 'G']
        for i in range(self.mut.nodeId.size()):
            mut_nodeId, mut_DS, mut_AS, mut_site, mut_time = self.mut.get_mutation(i)
            tc.mutations.add_row(site=mut_site, node=mut_nodeId, derived_state=allele[mut_DS], time=self.times[0] - mut_time)

        tc.sort()
        return tc.tree_sequence()

    def get_tree(self):
        if self.tree.shape[0] == 1:
            print('Genealogy was not simulated. Use VGsim.genealogy() method to simulate it.')
            sys.exit(1)
        return self.tree, self.times

    def get_tree_mutations(self):
        pass

    def get_tree_migrations(self):
        pass

    def get_tree_npis(self):
        pass

    def get_tree_recombinations(self):
        pass

    def get_data_infectious(self, pop, hap, step_num):
        time_points = [i*self.currentTime/step_num for i in range(step_num+1)]
        Data = np.zeros(step_num+1)
        Sample = np.zeros(step_num+1)
        Data[0] = self.initial_infectious[pop, hap]

        point = 0
        for i in range(self.events.ptr):
            while point != step_num and time_points[point] < self.events.times[i]:
                Data[point+1] = Data[point]
                Sample[point+1] = Sample[point]
                point += 1

            if self.events.types[i] == BIRTH and self.events.populations[i] == pop and self.events.haplotypes[i] == hap:
                Data[point] += 1
            elif self.events.types[i] == DEATH or self.events.types[i] == SAMPLING or self.events.types[i] == MUTATION and self.events.populations[i] == pop and self.events.haplotypes[i] == hap:
                Data[point] -= 1
                if self.events.types[i] == SAMPLING:
                    Sample[point] += 1
            elif self.events.types[i] == MUTATION and self.events.newHaplotypes[i] == hap and self.events.populations[i] == pop:
                Data[point] += 1
            elif self.events.types[i] == MIGRATION and self.events.newPopulations[i] == pop and self.events.haplotypes[i] == hap:
                Data[point] += 1
            elif self.events.types[i] == MULTITYPE:
                for j in range(self.events.haplotypes[i], self.events.populations[i]):
                    if self.multievents.types[j] == BIRTH and self.multievents.haplotypes[j] == hap and self.multievents.populations[j] == pop:
                        Data[point] += self.multievents.num[j]
                    elif self.multievents.types[j] == DEATH or self.multievents.types[j] == SAMPLING or self.multievents.types[j] == MUTATION and self.multievents.haplotypes[j] == hap and self.multievents.populations[j] == pop:
                        Data[point] -= self.multievents.num[j]
                        if self.multievents.types[j] == SAMPLING:
                            Sample[point] += self.multievents.num[j]
                    elif self.multievents.types[j] == MUTATION and self.multievents.newHaplotypes[j] == hap and self.multievents.populations[j] == pop:
                        Data[point] += self.multievents.num[j]
                    elif self.multievents.types[j] == MIGRATION and self.multievents.newPopulations[j] == pop and self.multievents.haplotypes[j] == hap:
                        Data[point] += self.multievents.num[j]

        Lockdowns = []
        for i in range(self.loc.times.size()):
            if self.loc.populationsId[i] == pop:
                Lockdowns.append([self.loc.states[i], self.loc.times[i]])

        return Data, Sample, time_points, Lockdowns

    def get_data_susceptible(self, pop, sus, step_num):
        time_points = [i*self.currentTime/step_num for i in range(step_num+1)]
        Data = np.zeros(step_num+1)
        Data[0] = self.initial_susceptible[pop, sus]

        point = 0
        for i in range(self.events.ptr):
            while point != step_num and time_points[point] < self.events.times[i]:
                Data[point+1] = Data[point]
                point += 1

            if self.events.types[i] == BIRTH and self.events.populations[i] == pop and self.events.newHaplotypes[i] == sus:
                Data[point] -= 1
            elif (self.events.types[i] == DEATH or self.events.types[i] == SAMPLING or self.events.types[i] == SUSCCHANGE) and self.events.populations[i] == pop and self.events.newHaplotypes[i] == sus:
                Data[point] += 1
            elif self.events.types[i] == SUSCCHANGE and self.events.haplotypes[i] == sus and self.events.populations[i] == pop:
                Data[point] -= 1
            elif self.events.types[i] == MIGRATION and self.events.newPopulations[i] == pop and self.events.newHaplotypes[i] == sus:
                Data[point] -= 1
            elif self.events.types[i] == MULTITYPE:
                for j in range(self.events.haplotypes[i], self.events.populations[i]):
                    if self.multievents.types[j] == BIRTH and self.multievents.newHaplotypes[j] == sus and self.multievents.populations[j] == pop:
                        Data[point] -= self.multievents.num[j]
                    elif (self.multievents.types[j] == DEATH or self.multievents.types[j] == SAMPLING or self.multievents.types[j] == SUSCCHANGE) and self.multievents.newHaplotypes[j] == sus and self.multievents.populations[j] == pop:
                        Data[point] += self.multievents.num[j]
                    elif self.multievents.types[j] == SUSCCHANGE and self.multievents.haplotypes[j] == sus and self.multievents.populations[j] == pop:
                        Data[point] -= self.multievents.num[j]
                    elif self.multievents.types[j] == MIGRATION and self.multievents.newPopulations[j] == pop and self.multievents.haplotypes[j] == sus:
                        Data[point] -= self.multievents.num[j]

        Lockdowns = []
        for i in range(self.loc.times.size()):
            if self.loc.populationsId[i] == pop:
                Lockdowns.append([self.loc.states[i], self.loc.times[i]])

        return Data, time_points, Lockdowns


    def Stats(self, time_simulation):
        print("Number of samples:", self.sCounter)
        print("Total number of iterations:", self.events.ptr)
        print('Success number:', self.good_attempt)
        print("Epidemic time:", self.currentTime)
        print('Simulation time:', time_simulation)
        print('Number of infections:', self.bCounter)
        print('Number of recoveries:', self.dCounter)
        if self.sites >= 1:
            print('Number of mutations:', self.mCounter)
        if self.popNum >= 2:
            print('Number of accepted migrations:', self.migPlus)
            print('Number of rejected migrations:', self.migNonPlus)
        check = False
        for sn in range(self.susNum):
            if self.suscepCumulTransition[sn] != 0.0:
                check = True

        if check:
            print('Number of immunity transitions:', self.iCounter)
        print('----------------------------------')

    def Debug(self):
        print("Parameters")
        print('first_simulation(mutable): ', self.first_simulation)
        print("sampling_probability(const): ", self.sampling_probability)
        print("memory_optimization(const): ", self.memory_optimization)
        print()
        print("sites(const): ", self.sites)
        print("hapNum(const): ", self.hapNum)
        print("currentHapNum(mutable): ", self.currentHapNum)
        print("maxHapNum(mutable): ", self.maxHapNum)
        print("popNum(const): ", self.popNum)
        print("susNum(const): ", self.susNum)
        print("bCounter(mutable): ", self.bCounter)
        print("dCounter(mutable): ", self.dCounter)
        print("sCounter(mutable): ", self.sCounter)
        print("mCounter(mutable): ", self.mCounter)
        print("iCounter(mutable):", self.iCounter)
        print("swapLockdown(mutable): ", self.swapLockdown)
        print("migPlus(mutable): ", self.migPlus)
        print("migNonPlus(mutable): ", self.migNonPlus)
        print("globalInfectious(mutable): ", self.globalInfectious)
        print()
        print("currentTime(mutable): ", self.currentTime)
        print("totalRate(mutable): ", self.totalRate)
        print("totalMigrationTate(mutable): ", self.totalMigrationRate)
        print("rn(mutable): ", self.rn)
        print()
        print("suscType(const): ", sep=" ", end="")
        for hn in range(self.hapNum):
            print(self.suscType[hn], end=" ")
        print()
        print()
        print("hapToNum(mutable): ", sep="", end="")
        for hn in range(self.hapNum):
            print(self.hapToNum[hn], end=" ")
        print()
        print("numToHap(mutable): ", sep="", end="")
        for hn in range(self.maxHapNum):
            print(self.numToHap[hn], end=" ")
        print()
        print("sizes(const): ", end="")
        for pn in range(self.popNum):
            print(self.sizes[pn], end=" ")
        print()
        print("totalSusceptible(mutable): ", end="")
        for pn in range(self.popNum):
            print(self.totalSusceptible[pn], end=" ")
        print()
        print("totalInfectious(mutable): ", end="")
        for pn in range(self.popNum):
            print(self.totalInfectious[pn], end=" ")
        print()
        print("lockdownON(mutable): ", end="")
        for pn in range(self.popNum):
            print(self.lockdownON[pn], end=" ")
        print()
        print()
        print("susceptible(mutable)----")
        for pn in range(self.popNum):
            for sn in range(self.susNum):
                print(self.susceptible[pn, sn], end=" ")
            print()
        print()
        print("infectious(mutable)----")
        for pn in range(self.popNum):
            for hn in range(self.currentHapNum):
                print(self.infectious[pn, hn], end=" ")
            print()
        print()
        print()
        print("Birth rate(const): ", sep="", end="")
        for hn in range(self.hapNum):
            print(self.bRate[hn], end=" ")
        print()
        print("Death rate(const): ", sep="", end="")
        for hn in range(self.hapNum):
            print(self.dRate[hn], end=" ")
        print()
        print("Sampling rate(const): ", sep="", end="")
        for hn in range(self.hapNum):
            print(self.sRate[hn], end=" ")
        print()
        print("tmRate(const): ", sep="", end="")
        for hn in range(self.currentHapNum):
            print(self.tmRate[hn], end=" ")
        print()
        print("maxEffectiveBirthMigration(const): ", sep="", end="")
        for pn in range(self.popNum):
            print(self.maxEffectiveBirthMigration[pn], end=" ")
        print()
        print("suscepCumulTransition(const): ", sep="", end="")
        for sn in range(self.susNum):
            print(self.suscepCumulTransition[sn], end=" ")
        print()
        print("immunePopRate(mutable): ", sep="", end="")
        for pn in range(self.popNum):
            print(self.immunePopRate[pn], end=" ")
        print()
        print("infectPopRate(mutable): ", sep="", end="")
        for pn in range(self.popNum):
            print(self.infectPopRate[pn], end=" ")
        print()
        print("popRate(mutable): ", sep="", end="")
        for pn in range(self.popNum):
            print(self.popRate[pn], end=" ")
        print()
        print("migPopRate(mutable): ", sep="", end="")
        for pn in range(self.popNum):
            print(self.migPopRate[pn], end=" ")
        print()
        print("actualSizes(const): ", end=" ")
        for pn in range(self.popNum):
            print(self.actualSizes[pn], end=" ")
        print()
        print("contactDensity(const): ", end=" ")
        for pn in range(self.popNum):
            print(self.contactDensity[pn], end=" ")
        print()
        print("contactDensityBeforeLockdown(const): ", end=" ")
        for pn in range(self.popNum):
            print(self.contactDensityBeforeLockdown[pn], end=" ")
        print()
        print("contactDensityBeforeLockdown(const): ", end=" ")
        for pn in range(self.popNum):
            print(self.contactDensityAfterLockdown[pn], end=" ")
        print()
        print("startLD(const): ", end=" ")
        for pn in range(self.popNum):
            print(self.startLD[pn], end=" ")
        print()
        print("endLD(const): ", end=" ")
        for pn in range(self.popNum):
            print(self.endLD[pn], end=" ")
        print()
        print("samplingMultiplier(const): ", end=" ")
        for pn in range(self.popNum):
            print(self.samplingMultiplier[pn], end=" ")
        print()
        print()
        print("mRate(const)----")
        for hn in range(self.hapNum):
            for s in range(self.sites):
                print(self.mRate[hn, s], end=" ")
            print()
        print()
        print("susceptibility(const)----")
        for hn in range(self.currentHapNum):
            for sn in range(self.susNum):
                print(self.susceptibility[hn, sn], end=" ")
            print()
        print()
        print("tEventHapPopRate(mutable)----")
        for pn in range(self.popNum):
            for hn in range(self.currentHapNum):
                print(self.tEventHapPopRate[pn, hn], end=" ")
            print()
        print()
        print("suscepTransition(const)----")
        for sn1 in range(self.susNum):
            for sn2 in range(self.susNum):
                print(self.suscepTransition[sn1, sn2], end=" ")
            print()
        print()
        print("immuneSourcePopRate(mutable)----")
        for pn in range(self.popNum):
            for sn in range(self.susNum):
                print(self.immuneSourcePopRate[pn, sn], end=" ")
            print()
        print()
        print("hapPopRate(mutable)----")
        for pn in range(self.popNum):
            for hn in range(self.currentHapNum):
                print(self.hapPopRate[pn, hn], end=" ")
            print()
        print()
        print("migrationRates(const)----")
        for pn1 in range(self.popNum):
            for pn2 in range(self.popNum):
                print(self.migrationRates[pn1, pn2], end=" ")
            print()
        print()
        print("effectiveMigration(const)----")
        for pn1 in range(self.popNum):
            for pn2 in range(self.popNum):
                print(self.effectiveMigration[pn1, pn2], end=" ")
            print()
        print()
        print()
        print("hapMutType(const)----")
        for hn in range(self.currentHapNum):
            for s in range(self.sites):
                for i in range(3):
                    print(self.hapMutType[hn, s, i], end=" ")
                print()
            print()
        print()
        print("eventHapPopRate(mutable)----")
        for pn in range(self.popNum):
            for hn in range(self.currentHapNum):
                for i in range(4):
                    print(self.eventHapPopRate[pn, hn, i], end=" ")
                print()
            print()
        print()
        print("susceptHapPopRate(mutable)----")
        for pn in range(self.popNum):
            for hn in range(self.currentHapNum):
                for sn in range(self.susNum):
                    print(self.susceptHapPopRate[pn, hn, sn], end=" ")
                print()
            print()
        print()


    def get_proportion(self):
        return self.migNonPlus / (self.events.ptr-1)

    #############################
    ### TAU LEAPING ALGORITHM ###
    #############################

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef void SimulatePopulation_tau(self, Py_ssize_t iterations, Py_ssize_t sample_size, float time, Py_ssize_t attempts):
        cdef:
            Py_ssize_t i, pn, propNum
            bint success

        self.PrepareParameters(iterations)
        self.CheckSizes()

        propNum = self.popNum*((self.popNum-1)*self.hapNum*self.susNum+self.susNum*(self.susNum-1)+self.hapNum*(2+self.sites*3+self.susNum))
        if self.globalInfectious == 0:
            self.FirstInfection()
        #TODO less memory for multievents
        self.multievents.CreateEvents(iterations*propNum)
        self.events.CreateEvents(iterations)
        self.UpdateAllRates()

        for i in range(attempts):
            self.seed = RndmWrapper(seed=(self.user_seed, i))
            if self.totalRate+self.totalMigrationRate != 0.0 and self.globalInfectious != 0:
                while ((self.events.ptr<self.events.size) and (sample_size==-1 or self.sCounter<sample_size) and (time==-1 or self.currentTime<time)):
                    self.Propensities()
                    self.ChooseTau()
                    success = False
                    while True:
                        success = self.GenerateEvents_tau()
                        if success:
                            break
                        else:
                            self.tau_l /= 2
                    self.currentTime += self.tau_l

                    self.UpdateCompartmentCounts_tau()
                    self.events.AddEvent(self.currentTime, MULTITYPE, self.multievents.ptr-propNum, self.multievents.ptr, 0, 0)
                    if self.globalInfectious == 0:
                        break
                    for pn in range(self.popNum):
                        self.CheckLockdown(pn)

            if self.events.ptr <= 100 and iterations > 100:
                self.Restart()
            else:
                self.good_attempt = i+1
                break

        if self.totalRate == 0.0 or self.globalInfectious == 0:
            print('Simulation finished because no infections individuals remain!')
        if self.events.ptr>=self.events.size:
            print("Achieved maximal number of iterations.")
        if self.sCounter>sample_size and sample_size!=-1:
            print("Achieved sample size.")
        if self.currentTime>time and time != -1:
            print("Achieved internal time limit.")
        if self.sCounter <= 1:
            print('\033[41m{}\033[0m'.format('WARNING!'), 'Simulated less 2 samples, so genealogy will not work!')

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void Propensities(self):
        cdef Py_ssize_t pn, spn, tpn, hn, sn, ssn, tsn, s, i

        for pn in range(self.popNum):
            for hn in range(self.hapNum):
                self.infectiousAuxTau[pn, hn] = 0.0
            for sn in range(self.susNum):
                self.susceptibleAuxTau[pn, sn] = 0.0

        for spn in range(self.popNum):
            for tpn in range(self.popNum):
                if spn == tpn:
                    continue
                for sn in range(self.susNum):
                    for hn in range(self.hapNum):
                        self.PropensitiesMigr[spn, tpn, sn, hn] = self.effectiveMigration[tpn, spn]*self.susceptible[tpn, sn]*\
                        self.infectious[spn, hn]*self.bRate[hn]*self.susceptibility[hn, sn]*self.migrationRates[spn, spn]

                        self.infectiousAuxTau[tpn, hn] += self.PropensitiesMigr[spn, tpn, sn, hn]
                        self.susceptibleAuxTau[tpn, sn] -= self.PropensitiesMigr[spn, tpn, sn, hn]

        for pn in range(self.popNum):
            #Susceptibility transition
            for ssn in range(self.susNum):
                for tsn in range(self.susNum):
                    if ssn == tsn:
                        continue
                    self.PropensitiesSuscep[pn, ssn, tsn] = self.suscepTransition[ssn, tsn]*self.susceptible[pn, ssn]

                    self.susceptibleAuxTau[pn, tsn] += self.PropensitiesSuscep[pn, ssn, tsn]
                    self.susceptibleAuxTau[pn, ssn] -= self.PropensitiesSuscep[pn, ssn, tsn]

            #Infectious-realted event
            for hn in range(self.hapNum):
                #Recovery
                self.PropensitiesRecovery[pn, hn] = self.dRate[hn]*self.infectious[pn, hn]

                self.susceptibleAuxTau[pn, self.suscType[hn]] += self.PropensitiesRecovery[pn, hn]
                self.infectiousAuxTau[pn, hn] -= self.PropensitiesRecovery[pn, hn]

                #Sampling
                self.PropensitiesSampling[pn, hn] = self.sRate[hn]*self.infectious[pn, hn]*self.samplingMultiplier[pn]

                self.susceptibleAuxTau[pn, self.suscType[hn]] += self.PropensitiesSampling[pn, hn]
                self.infectiousAuxTau[pn, hn] -= self.PropensitiesSampling[pn, hn]

                #Mutation
                for s in range(self.sites):
                    for i in range(3):
                        self.PropensitiesMutatations[pn, hn, s, i] = self.mRate[hn, s]*self.hapMutType[hn, s, i]/\
                        (self.hapMutType[hn, s, 0] + self.hapMutType[hn, s, 1] + self.hapMutType[hn, s, 2])*self.infectious[pn, hn]

                        self.infectiousAuxTau[pn, self.Mutate(hn, s, i)] += self.PropensitiesMutatations[pn, hn, s, i]
                        self.infectiousAuxTau[pn, hn] -= self.PropensitiesMutatations[pn, hn, s, i]
                
        #Transmission
        for tpn in range(self.popNum):
            for hn in range(self.hapNum):
                for sn in range(self.susNum):
                    self.PropensitiesTransmission[tpn, hn, sn] = 0.0
                    for spn in range(self.popNum):
                        self.PropensitiesTransmission[tpn, hn, sn] += self.bRate[hn]*self.susceptibility[hn, sn]*self.migrationRates[tpn, spn]*\
                        self.migrationRates[tpn, spn]*self.contactDensity[spn]*self.susceptible[tpn, sn]*self.infectious[tpn, hn]/\
                        self.actualSizes[spn]

                    self.infectiousAuxTau[tpn, hn] += self.PropensitiesTransmission[tpn, hn, sn]
                    self.susceptibleAuxTau[tpn, sn] -= self.PropensitiesTransmission[tpn, hn, sn]

    @cython.cdivision(True)
    cdef Py_ssize_t Mutate(self, Py_ssize_t hi, Py_ssize_t s, Py_ssize_t DS):
        cdef Py_ssize_t digit4, AS

        digit4 = 4**(self.sites-s-1)
        AS = int(floor(hi/digit4) % 4)
        if DS >= AS:
            DS += 1
        return hi + (DS-AS)*digit4

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void ChooseTau(self, float epsilon=0.03):
        cdef:
            Py_ssize_t pn, hn, sn
            double tmp

        self.tau_l = 1.0
        for pn in range(self.popNum):
            for hn in range(self.hapNum):
                if abs(self.infectiousAuxTau[pn, hn]) < 1e-8:
                    continue
                tmp = max(epsilon*self.infectious[pn,hn]/2.0,1.0)/abs(self.infectiousAuxTau[pn, hn])
                if tmp < self.tau_l:
                    self.tau_l = tmp
            for sn in range(self.susNum):
                if abs(self.susceptibleAuxTau[pn, sn]) < 1e-8:
                    continue
                tmp = max(epsilon*self.susceptible[pn,sn]/2.0,1.0)/abs(self.susceptibleAuxTau[pn, sn])
                if tmp < self.tau_l:
                    self.tau_l = tmp

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef bint GenerateEvents_tau(self):
        cdef Py_ssize_t event_num, pn, spn, tpn, hn, sn, ssn, tsn, s, i

        for pn in range(self.popNum):
            for sn in range(self.susNum):
                self.susceptibleDelta[pn, sn] = 0
            for hn in range(self.hapNum):
                self.infectiousDelta[pn, hn] = 0

        #Migrations
        for spn in range(self.popNum):
            for tpn in range(self.popNum):
                if spn == tpn:
                    continue
                for sn in range(self.susNum):
                    for hn in range(self.hapNum):
                        event_num = self.DrawEventsNum(self.PropensitiesMigr[spn, tpn, sn, hn]*self.tau_l)
                        self.eventsMigr[spn, tpn, sn, hn] = event_num

                        self.infectiousDelta[spn, hn] += event_num
                        self.susceptibleDelta[tpn, sn] -= event_num


        for pn in range(self.popNum):
            #Susceptibility transition
            for ssn in range(self.susNum):
                for tsn in range(self.susNum):
                    if ssn == tsn:
                        continue
                    event_num = self.DrawEventsNum(self.PropensitiesSuscep[pn, ssn, tsn]*self.tau_l)
                    self.eventsSuscep[pn, ssn, tsn] = event_num

                    self.susceptibleDelta[pn, tsn] += event_num
                    self.susceptibleDelta[pn, ssn] -= event_num

            #Infectious-realted event
            for hn in range(self.hapNum):
                #Recovery
                event_num = self.DrawEventsNum(self.PropensitiesRecovery[pn, hn]*self.tau_l)
                self.eventsRecovery[pn, hn] = event_num

                self.susceptibleDelta[pn, self.suscType[hn]] += event_num
                self.infectiousDelta[pn, hn] -= event_num

                #Sampling
                event_num = self.DrawEventsNum(self.PropensitiesSampling[pn, hn]*self.tau_l)
                self.eventsSampling[pn, hn] = event_num
                
                self.susceptibleDelta[pn, self.suscType[hn]] += event_num
                self.infectiousDelta[pn, hn] -= event_num

                #Mutation
                for s in range(self.sites):
                    for i in range(3):
                        event_num = self.DrawEventsNum(self.PropensitiesMutatations[pn, hn, s, i]*self.tau_l)
                        self.eventsMutatations[pn, hn, s, i] = event_num
                        
                        self.infectiousDelta[pn, self.Mutate(hn, s, i)] += event_num
                        self.infectiousDelta[pn, hn] -= event_num

                #Transmission
                for sn in range(self.susNum):
                    event_num = self.DrawEventsNum(self.PropensitiesTransmission[pn, hn, sn]*self.tau_l)
                    self.eventsTransmission[pn, hn, sn] = event_num
                    
                    self.infectiousDelta[pn, hn] += event_num
                    self.susceptibleDelta[pn, sn] -= event_num

        for pn in range(self.popNum):
            for sn in range(self.susNum):
                if self.susceptibleDelta[pn, sn] + self.susceptible[pn, sn] < 0 or self.susceptibleDelta[pn, sn] + self.susceptible[pn, sn] > self.sizes[pn]:
                    return False
            for hn in range(self.hapNum):
                if self.infectiousDelta[pn, hn] + self.infectious[pn, hn] < 0 or self.infectiousDelta[pn, hn] + self.infectious[pn, hn] > self.sizes[pn]:
                    return False
        return True

    cdef inline Py_ssize_t DrawEventsNum(self, double tau):
        return random_poisson(self.seed.rng, lam=tau)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void UpdateCompartmentCounts_tau(self):
        #TODO not write events with 0 num
        cdef Py_ssize_t event_num, pn, spn, tpn, hn, sn, ssn, tsn, s, i

        #Migrations
        for spn in range(self.popNum):
            for tpn in range(self.popNum):
                if spn == tpn:
                    continue
                for sn in range(self.susNum):
                    for hn in range(self.hapNum):
                        event_num = self.eventsMigr[spn, tpn, sn, hn]
                        self.NewInfections(tpn, sn, hn, event_num)
                        self.multievents.AddEvents(event_num, self.currentTime, MIGRATION, hn, spn, sn, tpn)
                        self.migPlus += event_num

        for pn in range(self.popNum):
            #Susceptibility transition
            for ssn in range(self.susNum):
                for tsn in range(self.susNum):
                    if ssn == tsn:
                        continue
                    event_num = self.eventsSuscep[pn, ssn, tsn]
                    self.susceptible[pn, tsn] += event_num
                    self.susceptible[pn, ssn] -= event_num
                    self.multievents.AddEvents(event_num, self.currentTime, SUSCCHANGE, ssn, pn, tsn, 0)
                    self.iCounter += event_num

            #Infectious-realted event
            for hn in range(self.hapNum):
                #Recovery
                event_num = self.eventsRecovery[pn, hn]
                self.NewRecoveries(pn, self.suscType[hn], hn, event_num)
                self.multievents.AddEvents(event_num, self.currentTime, DEATH, hn, pn, self.suscType[hn], 0)
                self.dCounter += event_num

                #Sampling
                event_num = self.eventsSampling[pn, hn]
                self.NewRecoveries(pn, self.suscType[hn], hn, event_num)
                self.multievents.AddEvents(event_num, self.currentTime, SAMPLING, hn, pn, self.suscType[hn], 0)
                self.sCounter += event_num

                #Mutation
                for s in range(self.sites):
                    for i in range(3):
                        nhn = self.Mutate(hn, s, i)
                        event_num = self.eventsMutatations[pn, hn, s, i]
                        self.infectious[pn, nhn] += event_num
                        self.infectious[pn, hn] -= event_num
                        self.multievents.AddEvents(event_num, self.currentTime, MUTATION, hn, pn, nhn, 0)
                        self.mCounter += event_num

                #Transmission
                for sn in range(self.susNum):
                    event_num = self.eventsTransmission[pn, hn, sn]
                    self.NewInfections(pn, sn, hn, event_num)
                    self.multievents.AddEvents(event_num, self.currentTime, BIRTH, hn, pn, sn, 0)
                    self.bCounter += event_num

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def Get_MultiEvents(self, id = None):
        if id is None:
            for i in range(self.multievents.ptr):
                ev = self.multievents.GetEvent(i)
                ev.PrintEvent()
        else:
            ev = self.multievents.GetEvent(id)
            ev.PrintEvent()


    def PrintCounters(self):
        print("Birth counter(mutable): ", self.bCounter)
        print("Death counter(mutable): ", self.dCounter)
        print("Sampling counter(mutable): ", self.sCounter)
        print("Mutation counter(mutable): ", self.mCounter)
        print("Immunity transition counter(mutable):", self.iCounter)
        print("Migration counter(mutable):", self.migPlus)

    def PrintPropensities(self):
        #self.FirstInfection()
        self.UpdateAllRates()
        self.Propensities()
        print("Migrations")
        for s in range(self.popNum):
            for r in range(self.popNum):
                if s == r:
                    continue
                for i in range(self.susNum):
                    for h in range(self.hapNum):
                        print (s, r, i, h, self.PropensitiesMigr[s, r, i, h])


        for s in range(self.popNum):
            print("Susceptibility transition")
            for i in range(self.susNum):
                for j in range(self.susNum):
                    if i == j:
                        continue
                    print(s, i, j, self.PropensitiesSuscep[s, i, j])

            #Infectious-realted event
            for h in range(self.hapNum):
                print("Recovery ", s, h, self.suscType[h], self.PropensitiesRecovery[s, h])
                print("Sampling ", s, h, self.suscType[h], self.PropensitiesSampling[s, h])

                for site in range(self.sites):
                    for i in range(3):
                        #ht = self.Mutate(h, site, i)
                        print("Mutation", s, h, site, i, self.PropensitiesMutatations[s, h, site, i])
                #Transmission
                for i in range(self.susNum):
                    print("Transmission", s, h, i, self.PropensitiesTransmission[s, h, i])
                    #print("migr=", self.migrationRates[s, s], "  prop", prop)
