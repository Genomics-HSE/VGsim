# cython: language_level=3
# cython: initializedcheck = False
# distutils: language = c++

cimport cython

from libc.math cimport log, floor, abs
from libcpp.vector cimport vector
from mc_lib.rndm cimport RndmWrapper

from prettytable import PrettyTable
import numpy as np
cimport numpy as np
import sys
import os
import time

from numpy.random.c_distributions cimport random_poisson, random_hypergeometric

include "fast_choose.pxi"
include "models.pxi"
include "events.pxi"

#pi - population ID, pn - popoulation number, spi - source population ID, tpi - target population ID
#hi - haplotype ID, hn - haplotype number, nhi - new haplotype number
#si - susceptibility ID, sn - susceptibility number, ssi - source susceptibility ID, tsi - target susceptibility ID
cdef class BirthDeathModel:
    cdef:
        RndmWrapper rndm

        bint first_simulation, sampling_probability, strong_migration
        Py_ssize_t sites, hapNum, popNum, susNum, bCounter, dCounter, sCounter, mCounter, iCounter, swapLockdown, migPlus, migNonPlus, migCounter, globalInfectious, countsPerStep
        double currentTime, seed, maxEffectiveBirth, totalRate, totalMigrationRate, totalLen, rn, tau_l

        Events events
        multiEvents multievents
        Mutations mut
        Migrations mig
        Lockdowns loc

        long[::1] suscType, sizes, totalSusceptible, totalInfectious, lockdownON, tree
        long[:,::1] susceptible, liveBranches

        double[::1] bRate, dRate, sRate, tmRate, maxEffectiveBirthMigration, maxSusceptibility, suscepCumulTransition, immunePopRate, infectPopRate, popRate, migPopRate, effectiveSizes, contactDensity, contactDensityBeforeLockdown, contactDensityAfterLockdown, startLD, endLD, samplingMultiplier, maxEffectiveMigration, times
        double[:,::1] mRate, totalHapMutType, tEventHapPopRate, susceptibility, suscepTransition, immuneSourcePopRate, hapPopRate, migrationRates, effectiveMigration
        double[:,:,::1] hapMutType, eventHapPopRate, susceptHapPopRate

        double[:,:,:,::1] PropensitiesMigr, PropensitiesMutatations
        double[:,:,::1] PropensitiesSuscep, PropensitiesTransmission
        double[:,::1] PropensitiesRecovery, PropensitiesSampling

        long[:,:,:,::1] eventsMigr, eventsMutatations
        long[:,:,::1] eventsSuscep, eventsTransmission
        long[:,::1] eventsRecovery, eventsSampling
        double[:,:,::1] infectiousAuxTau, susceptibleAuxTau
        long[:,::1] infectiousDelta, susceptibleDelta

    def __init__(self, sites_number, populations_number, susceptibility_types, seed, sampling_probability, strong_migration):
        self.rndm = RndmWrapper(seed=(seed, 0))

        self.first_simulation = False
        self.sampling_probability = sampling_probability
        self.strong_migration = strong_migration

        self.sites = sites_number
        self.hapNum = 4**self.sites
        self.susNum = susceptibility_types
        self.popNum = populations_number

        self.bCounter = 0
        self.dCounter = 0
        self.sCounter = 0
        self.mCounter = 0
        self.iCounter = 0
        self.swapLockdown = 0
        self.migPlus = 0
        self.migNonPlus = 0
        self.migCounter = 0
        self.globalInfectious = 0

        self.currentTime = 0.0
        self.tau_l=0.01
        self.seed = seed
        self.maxEffectiveBirth = 0.0
        self.totalRate = 0.0
        self.totalMigrationRate = 0.0

        self.events = Events()
        self.multievents = multiEvents()
        self.mut = Mutations()
        self.mig = Migrations()
        self.loc = Lockdowns()

        self.suscType = np.zeros(self.hapNum, dtype=np.int64)
        self.sizes = np.zeros(self.popNum, dtype=np.int64)
        self.totalSusceptible = np.zeros(self.popNum, dtype=np.int64)
        self.totalInfectious = np.zeros(self.popNum, dtype=np.int64)
        self.lockdownON = np.zeros(self.popNum, dtype=np.int64)

        self.susceptible = np.zeros((self.popNum, self.susNum), dtype=np.int64)
        self.liveBranches = np.zeros((self.popNum, self.hapNum), dtype=np.int64)

        self.bRate = np.zeros(self.hapNum, dtype=float)
        self.dRate = np.zeros(self.hapNum, dtype=float)
        self.sRate = np.zeros(self.hapNum, dtype=float)
        self.tmRate = np.zeros(self.hapNum, dtype=float)
        self.maxEffectiveBirthMigration = np.zeros(self.popNum, dtype=float)
        self.maxSusceptibility = np.zeros(self.hapNum, dtype=float)
        self.suscepCumulTransition = np.zeros(self.susNum, dtype=float)
        self.infectPopRate = np.zeros(self.popNum, dtype=float)
        self.immunePopRate = np.zeros(self.popNum, dtype=float)
        self.popRate = np.zeros(self.popNum, dtype=float)
        self.migPopRate = np.zeros(self.popNum, dtype=float)
        self.effectiveSizes = np.zeros(self.popNum, dtype=float)
        self.contactDensity = np.ones(self.popNum, dtype=float)
        self.contactDensityBeforeLockdown = np.ones(self.popNum, dtype=float)
        self.contactDensityAfterLockdown = np.zeros(self.popNum, dtype=float)
        self.startLD = np.ones(self.popNum, dtype=float)
        self.endLD = np.ones(self.popNum, dtype=float)
        self.samplingMultiplier = np.ones(self.popNum, dtype=float)
        self.maxEffectiveMigration = np.zeros(self.popNum, dtype=float)

        self.mRate = np.zeros((self.hapNum, self.sites), dtype=float)
        self.totalHapMutType = np.zeros((self.hapNum, self.sites), dtype=float)
        self.tEventHapPopRate = np.zeros((self.popNum, self.hapNum), dtype=float)
        self.susceptibility = np.zeros((self.hapNum, self.susNum), dtype=float)
        self.suscepTransition = np.zeros( (self.susNum, self.susNum), dtype=float)
        self.immuneSourcePopRate = np.zeros((self.popNum, self.susNum), dtype=float)
        self.hapPopRate = np.zeros((self.popNum, self.hapNum), dtype=float)
        self.migrationRates = np.zeros((self.popNum, self.popNum), dtype=float)
        self.effectiveMigration = np.zeros((self.popNum, self.popNum), dtype=float)

        self.hapMutType = np.ones((self.hapNum, self.sites, 3), dtype=float)
        self.eventHapPopRate = np.zeros((self.popNum, self.hapNum, 4), dtype=float)
        self.susceptHapPopRate = np.zeros((self.popNum, self.hapNum, self.susNum), dtype=float)

        for hn in range(self.hapNum):
            self.bRate[hn] = 2.0
            self.dRate[hn] = 1.0
            self.sRate[hn] = 0.01
            for s in range(self.sites):
                self.mRate[hn, s] = 0.01
                self.totalHapMutType[hn, s] = 3
            self.susceptibility[hn, 0] = 1.0
        for pn in range(self.popNum):
            self.sizes[pn] = 1000000
            self.totalSusceptible[pn] = 1000000
            self.susceptible[pn, 0] = 1000000

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

        self.infectiousAuxTau = np.zeros((self.popNum, self.hapNum, 2), dtype=float)
        self.susceptibleAuxTau = np.zeros((self.popNum, self.susNum, 2), dtype=float)

        self.infectiousDelta = np.zeros((self.popNum, self.hapNum), dtype=np.int64)
        self.susceptibleDelta = np.zeros((self.popNum, self.susNum), dtype=np.int64)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void FirstInfection(self):
        if self.globalInfectious == 0:
            for sn in range(self.susNum):
                if self.susceptible[0, sn] != 0:
                    self.NewInfection(0, sn, 0)
                    return


    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef inline void NewInfection(self, Py_ssize_t pi, Py_ssize_t si, Py_ssize_t hi):
        self.susceptible[pi, si] -= 1
        self.totalSusceptible[pi] -= 1
        self.liveBranches[pi, hi] += 1
        self.totalInfectious[pi] += 1
        self.globalInfectious += 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef inline void NewInfections(self, Py_ssize_t num, Py_ssize_t pi, Py_ssize_t si, Py_ssize_t hi):
        self.susceptible[pi, si] -= num
        self.totalSusceptible[pi] -= num
        self.liveBranches[pi, hi] += num
        self.totalInfectious[pi] += num
        self.globalInfectious += num

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef inline void NewRecovery(self, Py_ssize_t pi, Py_ssize_t si, Py_ssize_t hi):
        self.susceptible[pi, si] += 1
        self.totalSusceptible[pi] += 1
        self.liveBranches[pi, hi] -= 1
        self.totalInfectious[pi] -= 1
        self.globalInfectious -= 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef inline void NewRecoveries(self, Py_ssize_t num, Py_ssize_t pi, Py_ssize_t si, Py_ssize_t hi):
        self.susceptible[pi, si] += num
        self.totalSusceptible[pi] += num
        self.liveBranches[pi, hi] -= num
        self.totalInfectious[pi] -= num
        self.globalInfectious -= num

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void UpdateAllRates(self):
        cdef double summa
        for pn in range(self.popNum):
            self.maxEffectiveMigration[pn] = 0.0

        for pn1 in range(self.popNum):
            self.migrationRates[pn1, pn1] = 1.0
            self.effectiveSizes[pn1] = 0.0
            for pn2 in range(self.popNum):
                if pn1 == pn2:
                    continue
                self.migrationRates[pn1, pn1] -= self.migrationRates[pn1, pn2]
                self.effectiveSizes[pn1] += self.migrationRates[pn2, pn1]*self.sizes[pn2]
            self.effectiveSizes[pn1] += self.migrationRates[pn1, pn1]*self.sizes[pn1]

        for pn1 in range(self.popNum):
            for pn2 in range(self.popNum):
                if pn1==pn2:
                    continue
                self.effectiveMigration[pn1, pn2] = self.migrationRates[pn1, pn2]*self.contactDensity[pn2]/self.effectiveSizes[pn2]+self.migrationRates[pn2, pn1]*self.contactDensity[pn1]/self.effectiveSizes[pn1]
                if self.strong_migration:
                    self.effectiveMigration[pn1, pn2] *= self.migrationRates[pn2, pn2]
                if self.effectiveMigration[pn1, pn2] > self.maxEffectiveMigration[pn2]:
                    self.maxEffectiveMigration[pn2] = self.effectiveMigration[pn1, pn2]

        for sn1 in range(self.susNum):
            self.suscepCumulTransition[sn1] = 0
            for sn2 in range(self.susNum):
                self.suscepCumulTransition[sn1] += self.suscepTransition[sn1, sn2]

        self.totalRate = 0.0
        for pn in range(self.popNum):
            self.infectPopRate[pn] = 0
            self.immunePopRate[pn] = 0
            self.popRate[pn] = 0
        for pn in range(self.popNum):
            for hn in range(self.hapNum):
                self.tmRate[hn] = 0
                for s in range(self.sites):
                    self.tmRate[hn] += self.mRate[hn, s]

                self.eventHapPopRate[pn, hn, 0] = self.BirthRate(pn, hn)
                self.eventHapPopRate[pn, hn, 1] = self.dRate[hn]
                self.eventHapPopRate[pn, hn, 2] = self.sRate[hn] * self.samplingMultiplier[pn]
                self.eventHapPopRate[pn, hn, 3] = self.tmRate[hn]
                self.tEventHapPopRate[pn, hn] = 0
                for i in range(4):
                    self.tEventHapPopRate[pn, hn] += self.eventHapPopRate[pn, hn, i]
                self.hapPopRate[pn, hn] = self.tEventHapPopRate[pn, hn] * self.liveBranches[pn, hn]
                self.infectPopRate[pn] += self.hapPopRate[pn, hn]
            for sn in range(self.susNum):
                self.immuneSourcePopRate[pn, sn] = self.suscepCumulTransition[sn] * self.susceptible[pn, sn]
                self.immunePopRate[pn] += self.immuneSourcePopRate[pn, sn]
            self.popRate[pn] = self.infectPopRate[pn] + self.immunePopRate[pn]
            self.totalRate += self.popRate[pn]


        self.maxEffectiveBirth = 0.0
        for hn in range(self.hapNum):
            self.maxSusceptibility[hn] = 0.0
            for sn in range(self.susNum):
                if self.susceptibility[hn, sn] > self.maxSusceptibility[hn]:
                    self.maxSusceptibility[hn] = self.susceptibility[hn, sn]
            if self.maxEffectiveBirth < self.bRate[hn]*self.maxSusceptibility[hn]:
                self.maxEffectiveBirth = self.bRate[hn]*self.maxSusceptibility[hn]

        self.totalMigrationRate = 0.0
        for pn in range(self.popNum):
            self.maxEffectiveBirthMigration[pn] = self.maxEffectiveMigration[pn]*self.maxEffectiveBirth
            self.migPopRate[pn] = self.maxEffectiveBirthMigration[pn]*self.totalSusceptible[pn]*(self.globalInfectious-self.totalInfectious[pn])
            self.totalMigrationRate += self.migPopRate[pn]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef inline double BirthRate(self, Py_ssize_t pi, Py_ssize_t hi):
        cdef double ws = 0.0
        for sn in range(self.susNum):
            self.susceptHapPopRate[pi, hi, sn] = self.susceptible[pi, sn]*self.susceptibility[hi, sn]
            ws += self.susceptHapPopRate[pi, hi, sn]
        if self.strong_migration:
            return self.bRate[hi]*ws*self.migrationRates[pi, pi]*self.migrationRates[pi, pi]/self.effectiveSizes[pi]*self.contactDensity[pi]
        else:
            return self.bRate[hi]*ws/self.effectiveSizes[pi]*self.contactDensity[pi]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef void SimulatePopulation(self, Py_ssize_t iterations, Py_ssize_t sample_size, float time):
        cdef Py_ssize_t pi
        self.events.CreateEvents(iterations)
        if self.first_simulation == False:
            self.FirstInfection()
            self.first_simulation = True
        for pn in range(self.popNum):
            self.CheckLockdown(pn)
        self.UpdateAllRates()
        #self.totalLen = 0.0
        if self.totalRate != 0.0 and self.globalInfectious != 0:
            while (self.events.ptr<self.events.size and (sample_size==-1 or self.sCounter<=sample_size) and (time==-1 or self.currentTime<time)):
                self.SampleTime()
                pi = self.GenerateEvent()
                if self.totalRate == 0.0 or self.globalInfectious == 0:
                    print('#TODO')
                    break
                self.CheckLockdown(pi)

        if self.events.ptr>=self.events.size:
            print("Achieved maximal number of iterations.")
        if self.sCounter>sample_size and sample_size!=-1:
            print("Achieved sample size.")
        if self.currentTime>time and time != -1:
            print("Achieved internal time limit.")

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef inline void SampleTime(self):
        cdef double tau = - log(self.rndm.uniform()) / self.totalRate
        self.currentTime += tau
        #self.totalLen += tau*self.globalInfectious

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef Py_ssize_t GenerateEvent(self):
        cdef:
            Py_ssize_t pi, hi, ei
            double choose

        self.rn = self.rndm.uniform()
        choose = self.rn * (self.totalRate + self.totalMigrationRate)
        if self.totalRate > choose:
            self.rn = choose / self.totalRate
            pi, self.rn = fastChoose1(self.popRate, self.totalRate, self.rn)
            choose = self.rn * self.popRate[pi]
            if self.immunePopRate[pi] > choose:
                self.rn = choose / self.immunePopRate[pi]
                self.ImmunityTransition(pi)
            else:
                self.rn = (choose - self.immunePopRate[pi]) / self.infectPopRate[pi]
                hi, self.rn = fastChoose1(self.hapPopRate[pi], self.infectPopRate[pi], self.rn)
                ei, self.rn = fastChoose1(self.eventHapPopRate[pi, hi], self.tEventHapPopRate[pi, hi], self.rn)
                if ei == BIRTH:
                    self.Birth(pi, hi)
                elif ei == DEATH:
                    self.Death(pi, hi)
                elif ei == SAMPLING:
                    self.Sampling(pi, hi)
                else:
                    self.Mutation(pi, hi)
        else:
            self.rn = (choose - self.totalRate) / self.totalMigrationRate
            pi = self.GenerateMigration()
        return pi

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void UpdateRates(self, Py_ssize_t pi, bint infect, bint immune, bint migration):
        cdef double tmp
        if infect:
            self.infectPopRate[pi] = 0.0
            for hn in range(self.hapNum):
                self.eventHapPopRate[pi, hn, 0] = self.BirthRate(pi, hn)
                tmp = (self.eventHapPopRate[pi, hn, 0] +
                       self.eventHapPopRate[pi, hn, 1] +
                       self.eventHapPopRate[pi, hn, 2] +
                       self.eventHapPopRate[pi, hn, 3] )
                self.tEventHapPopRate[pi, hn] = tmp
                self.hapPopRate[pi, hn] = self.tEventHapPopRate[pi, hn] * self.liveBranches[pi, hn]
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
                self.migPopRate[pn] = self.maxEffectiveBirthMigration[pn]*self.totalSusceptible[pn]*(self.globalInfectious-self.totalInfectious[pn])
                self.totalMigrationRate += self.migPopRate[pn]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void ImmunityTransition(self, Py_ssize_t pi):
        cdef:
            Py_ssize_t si, ti
        si, self.rn = fastChoose1(self.immuneSourcePopRate[pi], self.immunePopRate[pi], self.rn)
        ti, self.rn = fastChoose1(self.suscepTransition[si], self.suscepCumulTransition[si], self.rn)

        self.susceptible[pi, si] -= 1
        self.susceptible[pi, ti] += 1
        self.immuneSourcePopRate[pi, si] = self.susceptible[pi, si]*self.suscepCumulTransition[si]
        self.immuneSourcePopRate[pi, ti] = self.susceptible[pi, ti]*self.suscepCumulTransition[ti]
        self.UpdateRates(pi, False, True, False)

        self.iCounter += 1
        self.events.AddEvent(self.currentTime, SUSCCHANGE, si, pi, ti, 0)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void Birth(self, Py_ssize_t pi, Py_ssize_t hi):
        cdef double ws = 0.0
        for sn in range(self.susNum):
            ws += self.susceptHapPopRate[pi, hi, sn]
        si, self.rn = fastChoose1(self.susceptHapPopRate[pi, hi], ws, self.rn)

        self.NewInfection(pi, si, hi)
        self.immuneSourcePopRate[pi, si] = self.suscepCumulTransition[si] * self.susceptible[pi, si]
        self.UpdateRates(pi, True, True, True)

        self.bCounter += 1
        self.events.AddEvent(self.currentTime, BIRTH, hi, pi, si, 0)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void Death(self, Py_ssize_t pi, Py_ssize_t hi, bint add_event = True):
        self.NewRecovery(pi, self.suscType[hi], hi)
        self.immuneSourcePopRate[pi, self.suscType[hi]] = self.susceptible[pi, self.suscType[hi]]*self.suscepCumulTransition[self.suscType[hi]]
        self.UpdateRates(pi, True, True, True)

        if add_event:
            self.dCounter += 1
            self.events.AddEvent(self.currentTime, DEATH, hi, pi, self.suscType[hi], 0)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void Sampling(self, Py_ssize_t pi, Py_ssize_t hi):
        self.Death(pi, hi, False)

        self.sCounter += 1
        self.events.AddEvent(self.currentTime, SAMPLING, hi, pi, self.suscType[hi], 0)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void Mutation(self, Py_ssize_t pi, Py_ssize_t hi):
        cdef:
            Py_ssize_t mi, digit4, AS, DS, nhi

        mi, self.rn = fastChoose1(self.mRate[hi], self.tmRate[hi], self.rn)
        digit4 = 4**mi
        AS = int(floor(hi/digit4) % 4)
        DS, self.rn = fastChoose1(self.hapMutType[hi, mi], self.totalHapMutType[hi, mi], self.rn)
        if DS >= AS:
            DS += 1
        nhi = hi + (DS-AS)*digit4

        self.liveBranches[pi, nhi] += 1
        self.liveBranches[pi, hi] -= 1
        self.UpdateRates(pi, True, False, False)

        self.mCounter += 1
        self.events.AddEvent(self.currentTime, MUTATION, hi, pi, nhi, 0)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef Py_ssize_t GenerateMigration(self):
        cdef:
            Py_ssize_t tpi, spi, hi, si
            double p_accept
        tpi, self.rn = fastChoose1(self.migPopRate, self.totalMigrationRate, self.rn)
        spi, self.rn = fastChoose2_skip(self.totalInfectious, self.globalInfectious-self.totalInfectious[tpi], self.rn, skip = tpi)
        hi, self.rn = fastChoose2(self.liveBranches[spi], self.totalInfectious[spi], self.rn)
        si, self.rn = fastChoose2(self.susceptible[tpi], self.totalSusceptible[tpi], self.rn)
        p_accept = self.effectiveMigration[spi, tpi]*self.bRate[hi]*self.susceptibility[hi, si]/self.maxEffectiveBirthMigration[tpi]
        if self.rn < p_accept:
            self.NewInfection(tpi, si, hi)
            self.UpdateRates(tpi, True, True, True)

            self.migPlus += 1
            self.migCounter += 1
            self.events.AddEvent(self.currentTime, MIGRATION, hi, spi, si, tpi)
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
    @cython.cdivision(True)
    cpdef GetGenealogy(self, seed):
        cdef:
            Py_ssize_t ptrTreeAndTime, n1, n2, id1, id2, id3, lbs, lbs_e, ns, nt, idt, ids, lbss, propNum
            double p, p1
            vector[vector[vector[Py_ssize_t]]] liveBranchesS, newLineages
            vector[vector[Py_ssize_t]] vecint2
            vector[Py_ssize_t] vecint1
            long[:,::1] liveBranches

            double e_time, me_time
            Py_ssize_t e_type_, e_population, e_haplotype, e_newHaplotype, e_newPopulation
            Py_ssize_t me_num, me_type_, me_population, me_haplotype, me_newHaplotype, me_newPopulation
            Py_ssize_t mt_ev_num, mt_ev_num2
        propNum = self.PropensitiesNumber()

        if self.sCounter < 2: #TODO if number of sampled leaves is 0 (probably 1 as well), then GetGenealogy seems to go to an infinite cycle
            print("Less than two cases were sampled...")
            print("_________________________________")
            sys.exit(0)

        if seed != None:
            self.rndm = RndmWrapper(seed=(seed, 0))

        ptrTreeAndTime = 0
        self.tree = np.zeros(2 * self.sCounter - 1, dtype=np.int64)
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
                lbs_e = self.liveBranches[e_population, e_haplotype]
                p = float(lbs)*(float(lbs)-1.0)/ float(lbs_e) / (float(lbs_e) - 1.0)
                if self.rndm.uniform() < p:
                    n1 = int(floor( lbs*self.rndm.uniform() ))
                    n2 = int(floor( (lbs-1)*self.rndm.uniform() ))
                    if n2 >= n1:
                        n2 += 1
                    id1 = liveBranchesS[e_population][e_haplotype][n1]
                    id2 = liveBranchesS[e_population][e_haplotype][n2]
                    id3 = ptrTreeAndTime
                    liveBranchesS[e_population][e_haplotype][n1] = id3
                    liveBranchesS[e_population][e_haplotype][n2] = liveBranchesS[e_population][e_haplotype][lbs-1]
                    liveBranchesS[e_population][e_haplotype].pop_back()
                    self.tree[id1] = id3
                    self.tree[id2] = id3
                    self.tree[ptrTreeAndTime] = -1
                    self.times[ptrTreeAndTime] = e_time
                    ptrTreeAndTime += 1
                self.liveBranches[e_population, e_haplotype] -= 1
            elif e_type_ == DEATH:
                self.liveBranches[e_population, e_haplotype] += 1
            elif e_type_ == SAMPLING:
                self.liveBranches[e_population, e_haplotype] += 1
                liveBranchesS[e_population][e_haplotype].push_back( ptrTreeAndTime )
                self.tree[ptrTreeAndTime] = -1
                self.times[ptrTreeAndTime] = e_time
                ptrTreeAndTime += 1
            elif e_type_ == MUTATION:
                lbs = liveBranchesS[e_population][e_newHaplotype].size()
                p = float(lbs)/self.liveBranches[e_population, e_newHaplotype]
                if self.rndm.uniform() < p:
                    n1 = int(floor( lbs*self.rndm.uniform() ))
                    id1 = liveBranchesS[e_population][e_newHaplotype][n1]
                    liveBranchesS[e_population][e_newHaplotype][n1] = liveBranchesS[e_population][e_newHaplotype][lbs-1]
                    liveBranchesS[e_population][e_newHaplotype].pop_back()
                    liveBranchesS[e_population][e_haplotype].push_back(id1)
                    self.mut.AddMutation(id1, e_haplotype, e_newHaplotype, e_time)
                self.liveBranches[e_population, e_newHaplotype] -= 1
                self.liveBranches[e_population, e_haplotype] += 1
            elif e_type_ == SUSCCHANGE:
                pass
            elif e_type_ == MIGRATION:
                lbs = liveBranchesS[e_newPopulation][e_haplotype].size()
                p = float(lbs)/self.liveBranches[e_newPopulation][e_haplotype]
                if self.rndm.uniform() < p:
                    nt = int(floor( lbs*self.rndm.uniform() ))
                    lbss = liveBranchesS[e_population][e_haplotype].size()
                    p1 = float(lbss)/self.liveBranches[e_population, e_haplotype]
                    if self.rndm.uniform() < p1:
                        ns = int(floor( lbss*self.rndm.uniform() ))
                        idt = liveBranchesS[e_newPopulation][e_haplotype][nt]
                        ids = liveBranchesS[e_population][e_haplotype][ns]
                        id3 = ptrTreeAndTime
                        liveBranchesS[e_population][e_haplotype][ns] = id3
                        liveBranchesS[e_newPopulation][e_haplotype][nt] = liveBranchesS[e_newPopulation][e_haplotype][lbs-1]
                        liveBranchesS[e_newPopulation][e_haplotype].pop_back()
                        self.tree[idt] = id3
                        self.tree[ids] = id3
                        self.tree[ptrTreeAndTime] = -1
                        self.times[ptrTreeAndTime] = e_time
                        ptrTreeAndTime += 1
                        self.mig.AddMigration(idt, e_time, e_population, e_newPopulation)
                    else:
                        liveBranchesS[e_population][e_haplotype].push_back(liveBranchesS[e_newPopulation][e_haplotype][nt])
                        liveBranchesS[e_newPopulation][e_haplotype][nt] = liveBranchesS[e_newPopulation][e_haplotype][lbs-1]
                        liveBranchesS[e_newPopulation][e_haplotype].pop_back()
                self.liveBranches[e_newPopulation, e_haplotype] -= 1
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
                        lbs_e = self.liveBranches[me_population, me_haplotype]
                        if me_num == 0 or lbs == 0:
                            mt_ev_num = 0
                        else:
                            mt_ev_num = random_hypergeometric(self.rndm.rng, int(lbs*(lbs-1.0)/2.0), int(lbs_e*(lbs_e-1)/2-lbs*(lbs-1)/2), me_num)
                        #print("me_num=", me_num, "  mt_ev_num", mt_ev_num)
                        for i in range(mt_ev_num):
                            n1 = int(floor( lbs*self.rndm.uniform() ))
                            n2 = int(floor( (lbs-1)*self.rndm.uniform() ))
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
                            self.tree[id1] = id3
                            self.tree[id2] = id3
                            self.tree[ptrTreeAndTime] = -1
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
                            self.tree[ptrTreeAndTime] = -1
                            self.times[ptrTreeAndTime] = me_time
                            ptrTreeAndTime += 1
                    elif me_type_ == MUTATION:
                        lbs = liveBranchesS[me_population][me_newHaplotype].size()
                        if me_num == 0 or lbs == 0:
                            mt_ev_num = 0
                        else:
                            mt_ev_num = random_hypergeometric(self.rndm.rng, lbs, self.liveBranches[me_population, me_newHaplotype]-lbs, me_num)
                        for i in range(mt_ev_num):
                            n1 = int(floor( lbs*self.rndm.uniform() ))
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
                            mt_ev_num = random_hypergeometric(self.rndm.rng, lbs, self.liveBranches[me_newPopulation, me_haplotype]-lbs, me_num)
                        #for i in range(mt_ev_num):
                            lbss = liveBranchesS[me_population][me_haplotype].size()
                            if mt_ev_num == 0 or lbss == 0:
                                mt_ev_num2 = 0
                            else:
                                mt_ev_num2 = random_hypergeometric(self.rndm.rng, lbss, self.liveBranches[me_population, me_haplotype]-lbss, mt_ev_num)
                            #print("me_num: ", me_num, "mt_ev_num: ", mt_ev_num, "mt_ev_num2: ", mt_ev_num2)
                            for i in range(mt_ev_num2):
                                nt = int(floor( lbs*self.rndm.uniform() ))
                                ns = int(floor( lbss*self.rndm.uniform() ))
                                idt = liveBranchesS[me_newPopulation][me_haplotype][nt]
                                ids = liveBranchesS[me_population][me_haplotype][ns]
                                id3 = ptrTreeAndTime
                                liveBranchesS[me_population][me_haplotype][ns] = liveBranchesS[me_population][me_haplotype][  liveBranchesS[me_population][me_haplotype].size()-1 ]
                                liveBranchesS[me_population][me_haplotype].pop_back()
                                liveBranchesS[me_newPopulation][me_haplotype][nt] = liveBranchesS[me_newPopulation][me_haplotype][lbs-1]
                                liveBranchesS[me_newPopulation][me_haplotype].pop_back()
                                newLineages[me_population][me_haplotype].push_back( id3 )
                                self.tree[idt] = id3
                                self.tree[ids] = id3
                                self.tree[ptrTreeAndTime] = -1
                                self.times[ptrTreeAndTime] = me_time
                                ptrTreeAndTime += 1
                                self.mig.AddMigration(idt, me_time, me_population, me_newPopulation)
                                lbss -= 1
                                lbs -= 1
                            for i in range(mt_ev_num-mt_ev_num2):
                                nt = int(floor( lbs*self.rndm.uniform() ))
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
                            self.liveBranches[pi, hi] += self.infectiousDelta[pi, hi]
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
        #for i in range(self.sCounter * 2 - 1):
        #    print(self.tree[i], end=" ")
        #print("")
        #deg = [0 for i in range(self.sCounter * 2 - 1)]
        #for i in range(self.sCounter * 2 - 1):
        #    deg[self.tree[i]] += 1
        #for i in range(self.sCounter * 2 - 1):
        #    print(deg[i], end=" ")

        #self.CheckTree()

    cdef void CheckTree(self):
        cdef Py_ssize_t counter
        counter = 0
        for i in range(self.sCounter * 2 - 1):
            if self.tree[i] == 0:
                print("Error 1")
                print("_________________________________")
                sys.exit(0)
            if self.tree[i] == 1:
                counter += 1
            if counter >= 2:
                print("Error 2")
                print("_________________________________")
                sys.exit(0)
            if self.tree[i] == i:
                print("Error 3")
                print("_________________________________")
                sys.exit(0)

                
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
        for hn in range(self.hapNum):
            list = ["\n" + self.calculate_string(hn), "\n" + str(self.bRate[hn]), "\n" + str(self.dRate[hn]), "\n" + str(self.sRate[hn]), "\n" + str(self.suscType[hn])]
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

    def print_populations(self):
        print("*****************")
        print("***Populations***")
        print("*****************")
        table_populations = PrettyTable()

        table_populations.field_names = ["ID", "Size", "CD",'CDBLC', "CDALD", "SLD", "ELD", "SM"]
        for pn in range(self.popNum):
            table_populations.add_row([pn, self.sizes[pn], self.contactDensity[pn], self.contactDensityBeforeLockdown[pn], self.contactDensityAfterLockdown[pn], self.startLD[pn], self.endLD[pn], self.samplingMultiplier[pn]])

        print(table_populations)
        print("Legend:")
        print("ID - number of population")
        print("Size - size of population")
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
                row.append(self.susceptible[pn, sn])
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
                row.append(self.liveBranches[pn, hn])
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
        for hn in range(self.hapNum):
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


    def Error(self, text):
        print(text)
        sys.exit(1)

    def create_list_haplotypes(self, haplotype):
        if haplotype.count("A") + haplotype.count("T") + haplotype.count("C") + haplotype.count("G") + haplotype.count("*") != self.sites:
            self.Error("Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"*\" and lenght of haplotype should be equal number of mutations sites.")

        haplotypes = [haplotype]
        for s in range(self.sites):
            for i in range(len(haplotypes)):
                haplotype_old = haplotypes[i]
                if haplotype_old[s] == "*":
                    haplotype = haplotype_old.replace("*", "A", 1)
                    haplotypes.append(haplotype)
                    haplotype = haplotype_old.replace("*", "T", 1)
                    haplotypes.append(haplotype)
                    haplotype = haplotype_old.replace("*", "C", 1)
                    haplotypes.append(haplotype)
                    haplotype = haplotype_old.replace("*", "G", 1)
                    haplotypes.append(haplotype)
        for i in range(len(haplotypes)-1, -1, -1):
            if haplotypes[i].count("*") != 0:
                haplotypes.remove(haplotypes[i])
        for i in range(len(haplotypes)):
            haplotypes[i] = self.calculate_haplotype(haplotypes[i])

        return haplotypes

    def create_mutations(self, haplotype, site):
        hap = self.calculate_string(haplotype)
        haplotypes = [hap[:site] + "A" + hap[site+1:], hap[:site] + "T" + hap[site+1:], hap[:site] + "C" + hap[site+1:], hap[:site] + "G" + hap[site+1:]]
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
        string = color_hap[3] + "->" + color_hap[0] + ": " + str(self.hapMutType[haplotype, site, 0]) + "\n" + color_hap[3] + "->" + color_hap[1] + ": " + str(self.hapMutType[haplotype, site, 1]) + "\n" + color_hap[3] + "->" + color_hap[2] + ": " + str(self.hapMutType[haplotype, site, 2]) + "\n"
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

    def get_sites(self):
        return self.sites

    def get_hapNum(self):
        return self.hapNum

    def get_popNum(self):
        return self.popNum

    def get_susNum(self):
        return self.susNum


    def set_transmission_rate(self, rate, haplotype):
        if isinstance(rate, (int, float)) == False:
            self.Error("Incorrect type of infection rate. Value should be int or float.")
        if rate<0:
            self.Error("Incorrect value of infection rate. Value should be more or equal 0.")

        if isinstance(haplotype, str):
            haplotypes = self.create_list_haplotypes(haplotype)
            for haplotype in haplotypes:
                self.bRate[haplotype] = rate
        elif isinstance(haplotype, int):
            if haplotype<0 or haplotype>=self.hapNum:
                self.Error("There are no such haplotype!")

            self.bRate[haplotype] = rate
        elif haplotype == None:
            for hn in range(self.hapNum):
                self.bRate[hn] = rate
        else:
            self.Error("Incorrect value of haplotype. Value should be string or int or None.")

    def set_recovery_rate(self, rate, haplotype):
        if isinstance(rate, (int, float)) == False:
            self.Error("Incorrect type of uninfection rate. Value should be int or float.")
        if rate<0:
            self.Error("Incorrect value of uninfection rate. Value should be more or equal 0.")

        if isinstance(haplotype, str):
            haplotypes = self.create_list_haplotypes(haplotype)
            for haplotype in haplotypes:
                self.dRate[haplotype] = rate
        elif isinstance(haplotype, int):
            if haplotype<0 or haplotype>=self.hapNum:
                self.Error("There are no such haplotype!")

            self.dRate[haplotype] = rate
        elif haplotype == None:
            for hn in range(self.hapNum):
                self.dRate[hn] = rate
        else:
            self.Error("Incorrect value of haplotype. Value should be string or int or None.")

    def set_sampling_rate(self, rate, haplotype):
        if self.sampling_probability == True:
            if isinstance(rate, (int, float)) == False:
                self.Error("Incorrect type of sampling probability. Value should be int or float.")
            if rate<0 or rate>1:
                self.Error("Incorrect value of sampling probability. Value should be more or equal 0 and less or equal 1.")

            if isinstance(haplotype, str):
                haplotypes = self.create_list_haplotypes(haplotype)
                for haplotype in haplotypes:
                    deathRate = self.dRate[haplotype]
                    self.dRate[haplotype] = (1-rate) * deathRate
                    self.sRate[haplotype] = rate * deathRate
            elif isinstance(haplotype, int):
                if haplotype<0 or haplotype>=self.hapNum:
                    self.Error("There are no such haplotype!")

                deathRate = self.dRate[haplotype]
                self.dRate[haplotype] = (1-rate) * deathRate
                self.sRate[haplotype] = rate * deathRate
            elif haplotype == None:
                for hn in range(self.hapNum):
                    deathRate = self.dRate[hn]
                    self.dRate[hn] = (1-rate) * deathRate
                    self.sRate[hn] = rate * deathRate
            else:
                self.Error("Incorrect value of haplotype. Value should be string or int or None.")

        elif self.sampling_probability == False:
            if isinstance(rate, (int, float)) == False:
                self.Error("Incorrect type of sampling rate. Value should be int or float.")
            if rate<0:
                self.Error("Incorrect value of sampling rate. Value should be more or equal 0.")

            if isinstance(haplotype, str):
                haplotypes = self.create_list_haplotypes(haplotype)
                for haplotype in haplotypes:
                    self.sRate[haplotype] = rate
            elif isinstance(haplotype, int):
                if haplotype<0 or haplotype>=self.hapNum:
                    self.Error("There are no such haplotype!")

                self.sRate[haplotype] = rate
            elif haplotype == None:
                for hn in range(self.hapNum):
                    self.sRate[hn] = rate
            else:
                self.Error("Incorrect value of haplotype. Value should be string or int or None.")
        else:
            self.Error("#TODO")

    def set_mutation_rate(self, rate, probabilities, haplotype, mutation):
        if isinstance(rate, (int, float)) and isinstance(probabilities, list) and isinstance(haplotype, str) and isinstance(mutation,int):#DONE
            if rate<0:
                self.Error("#TODO_mr1")
            if len(probabilities)!=4:
                self.Error("#TODO_mr2")
            for i in range(4):
                if isinstance(probabilities[i], (int, float)) == False:
                    self.Error("#TODO_mr3")
                if probabilities[i]<0:
                    self.Error("#TODO_mr4")
            haplotypes = self.create_list_haplotypes(haplotype)
            if mutation<0 or mutation>=self.sites:
                self.Error("There are no such mutation!")

            for haplotype in haplotypes:
                probabilities_allele = list(probabilities)
                del probabilities_allele[self.calculate_allele(haplotype, mutation)]
                if sum(probabilities_allele) == 0:
                    self.Error("#TODO_mr5")
                self.mRate[haplotype, mutation] = rate
                self.hapMutType[haplotype, mutation, 0] = probabilities_allele[0]
                self.hapMutType[haplotype, mutation, 1] = probabilities_allele[1]
                self.hapMutType[haplotype, mutation, 2] = probabilities_allele[2]
                self.totalHapMutType[haplotype, mutation] = sum(probabilities_allele)
        elif rate==None and isinstance(probabilities, list) and isinstance(haplotype, str) and isinstance(mutation,int):#DONE
            if len(probabilities)!=4:
                self.Error("#TODO_mr6")
            for i in range(4):
                if isinstance(probabilities[i], (int, float)) == False:
                    self.Error("#TODO_mr7")
                if probabilities[i]<0:
                    self.Error("#TODO_mr8")
            haplotypes = self.create_list_haplotypes(haplotype)
            if mutation<0 or mutation>=self.sites:
                self.Error("There are no such mutation!")

            for haplotype in haplotypes:
                probabilities_allele = list(probabilities)
                del probabilities_allele[self.calculate_allele(haplotype, mutation)]
                if sum(probabilities_allele) == 0:
                    self.Error("#TODO_mr9")
                self.hapMutType[haplotype, mutation, 0] = probabilities_allele[0]
                self.hapMutType[haplotype, mutation, 1] = probabilities_allele[1]
                self.hapMutType[haplotype, mutation, 2] = probabilities_allele[2]
                self.totalHapMutType[haplotype, mutation] = sum(probabilities_allele)
        elif isinstance(rate, (int, float)) and probabilities==None and isinstance(haplotype, str) and isinstance(mutation,int):#DONE
            if rate<0:
                self.Error("#TODO_mr10")
            haplotypes = self.create_list_haplotypes(haplotype)
            if mutation<0 or mutation>=self.sites:
                self.Error("There are no such mutation!")

            for haplotype in haplotypes:
                self.mRate[haplotype, mutation] = rate
        elif isinstance(rate, (int, float)) and isinstance(probabilities, list) and isinstance(haplotype, int) and isinstance(mutation,int):#DONE
            if rate<0:
                self.Error("#TODO_mr11")
            if len(probabilities)!=4:
                self.Error("#TODO_mr12")
            for i in range(4):
                if isinstance(probabilities[i], (int, float)) == False:
                    self.Error("#TODO_mr13")
                if probabilities[i]<0:
                    self.Error("#TODO_mr14")
            del probabilities[self.calculate_allele(haplotype, mutation)]
            if sum(probabilities) == 0:
                self.Error("#TODO_mr15")
            if haplotype<0 or haplotype>=self.hapNum:
                self.Error("There are no such haplotype!")
            if mutation<0 or mutation>=self.sites:
                self.Error("There are no such mutation!")

            self.mRate[haplotype, mutation] = rate
            self.hapMutType[haplotype, mutation, 0] = probabilities[0]
            self.hapMutType[haplotype, mutation, 1] = probabilities[1]
            self.hapMutType[haplotype, mutation, 2] = probabilities[2]
            self.totalHapMutType[haplotype, mutation] = sum(probabilities)
        elif rate==None and isinstance(probabilities, list) and isinstance(haplotype, int) and isinstance(mutation,int):#DONE
            if len(probabilities)!=4:
                self.Error("#TODO_mr16")
            for i in range(4):
                if isinstance(probabilities[i], (int, float)) == False:
                    self.Error("#TODO_mr17")
                if probabilities[i]<0:
                    self.Error("#TODO_mr18")
            del probabilities[self.calculate_allele(haplotype, mutation)]
            if sum(probabilities) == 0:
                self.Error("#TODO_mr19")
            if haplotype<0 or haplotype>=self.hapNum:
                self.Error("There are no such haplotype!")
            if mutation<0 or mutation>=self.sites:
                self.Error("There are no such mutation!")

            self.hapMutType[haplotype, mutation, 0] = probabilities[0]
            self.hapMutType[haplotype, mutation, 1] = probabilities[1]
            self.hapMutType[haplotype, mutation, 2] = probabilities[2]
            self.totalHapMutType[haplotype, mutation] = sum(probabilities)
        elif isinstance(rate, (int, float)) and probabilities==None and isinstance(haplotype, int) and isinstance(mutation,int):#DONE
            if rate<0:
                self.Error("#TODO_mr20")
            if haplotype<0 or haplotype>=self.hapNum:
                self.Error("There are no such haplotype!")
            if mutation<0 or mutation>=self.sites:
                self.Error("There are no such mutation!")

            self.mRate[haplotype, mutation] = rate
        elif isinstance(rate, (int, float)) and isinstance(probabilities, list) and haplotype==None and isinstance(mutation,int):#DONE
            if rate<0:
                self.Error("#TODO")
            if len(probabilities)!=4:
                self.Error("#TODO")
            for i in range(4):
                if isinstance(probabilities[i], (int, float)) == False:
                    self.Error("#TODO")
                if probabilities[i]<0:
                    self.Error("#TODO")
            if mutation<0 or mutation>=self.sites:
                self.Error("There are no such mutation!")

            for haplotype in range(self.hapNum):
                probabilities_allele = list(probabilities)
                del probabilities_allele[self.calculate_allele(haplotype, mutation)]
                if sum(probabilities_allele) == 0:
                    self.Error("#TODO")
                self.mRate[haplotype, mutation] = rate
                self.hapMutType[haplotype, mutation, 0] = probabilities_allele[0]
                self.hapMutType[haplotype, mutation, 1] = probabilities_allele[1]
                self.hapMutType[haplotype, mutation, 2] = probabilities_allele[2]
                self.totalHapMutType[haplotype, mutation] = sum(probabilities_allele)
        elif rate==None and isinstance(probabilities, list) and haplotype==None and isinstance(mutation,int):#DONE
            if len(probabilities)!=4:
                self.Error("#TODO")
            for i in range(4):
                if isinstance(probabilities[i], (int, float)) == False:
                    self.Error("#TODO")
                if probabilities[i]<0:
                    self.Error("#TODO")
            if mutation<0 or mutation>=self.sites:
                self.Error("There are no such mutation!")

            for haplotype in range(self.hapNum):
                probabilities_allele = list(probabilities)
                del probabilities_allele[self.calculate_allele(haplotype, mutation)]
                if sum(probabilities_allele) == 0:
                    self.Error("#TODO")
                self.hapMutType[haplotype, mutation, 0] = probabilities_allele[0]
                self.hapMutType[haplotype, mutation, 1] = probabilities_allele[1]
                self.hapMutType[haplotype, mutation, 2] = probabilities_allele[2]
                self.totalHapMutType[haplotype, mutation] = sum(probabilities_allele)
        elif isinstance(rate, (int, float)) and probabilities==None and haplotype==None and isinstance(mutation,int):#DONE
            if rate<0:
                self.Error("#TODO")
            if mutation<0 or mutation>=self.sites:
                self.Error("There are no such mutation!")

            for haplotype in range(self.hapNum):
                self.mRate[haplotype, mutation] = rate
        elif isinstance(rate, (int, float)) and isinstance(probabilities, list) and isinstance(haplotype, str) and mutation==None:
            if rate<0:
                self.Error("#TODO")
            if len(probabilities)!=4:
                self.Error("#TODO")
            for i in range(4):
                if isinstance(probabilities[i], (int, float)) == False:
                    self.Error("#TODO")
                if probabilities[i]<0:
                    self.Error("#TODO")
            haplotypes = self.create_list_haplotypes(haplotype)

            for haplotype in haplotypes:
                for s in range(self.sites):
                    probabilities_allele = list(probabilities)
                    del probabilities_allele[self.calculate_allele(haplotype, mutation)]
                    if sum(probabilities_allele) == 0:
                        self.Error("#TODO")

                    self.mRate[haplotype, s] = rate
                    self.hapMutType[haplotype, s, 0] = probabilities_allele[0]
                    self.hapMutType[haplotype, s, 1] = probabilities_allele[1]
                    self.hapMutType[haplotype, s, 2] = probabilities_allele[2]
                    self.totalHapMutType[haplotype, s] = sum(probabilities_allele)
        elif rate==None and isinstance(probabilities, list) and isinstance(haplotype, str) and mutation==None:#DONE
            if len(probabilities)!=4:
                self.Error("#TODO")
            for i in range(4):
                if isinstance(probabilities[i], (int, float)) == False:
                    self.Error("#TODO")
                if probabilities[i]<0:
                    self.Error("#TODO")
            haplotypes = self.create_list_haplotypes(haplotype)

            for haplotype in haplotypes:
                for s in range(self.sites):
                    probabilities_allele = list(probabilities)
                    del probabilities_allele[self.calculate_allele(haplotype, mutation)]
                    if sum(probabilities_allele) == 0:
                        self.Error("#TODO")
                    self.hapMutType[haplotype, s, 0] = probabilities_allele[0]
                    self.hapMutType[haplotype, s, 1] = probabilities_allele[1]
                    self.hapMutType[haplotype, s, 2] = probabilities_allele[2]
                    self.totalHapMutType[haplotype, s] = sum(probabilities_allele)
        elif isinstance(rate, (int, float)) and probabilities==None and isinstance(haplotype, str) and mutation==None:
            if rate<0:
                self.Error("#TODO")
            haplotypes = self.create_list_haplotypes(haplotype)

            for haplotype in haplotypes:
                for s in range(self.sites):
                    self.mRate[haplotype, s] = rate
        elif isinstance(rate, (int, float)) and isinstance(probabilities, list) and isinstance(haplotype, int) and mutation==None:#DONE
            if rate<0:
                self.Error("#TODO")
            if len(probabilities)!=4:
                self.Error("#TODO")
            for i in range(4):
                if isinstance(probabilities[i], (int, float)) == False:
                    self.Error("#TODO")
                if probabilities[i]<0:
                    self.Error("#TODO")
            del probabilities[self.calculate_allele(haplotype, mutation)]
            if sum(probabilities) == 0:
                self.Error("#TODO")
            if haplotype<0 or haplotype>=self.hapNum:
                self.Error("There are no such haplotype!")

            for s in range(self.sites):
                self.mRate[haplotype, s] = rate
                self.hapMutType[haplotype, s, 0] = probabilities[0]
                self.hapMutType[haplotype, s, 1] = probabilities[1]
                self.hapMutType[haplotype, s, 2] = probabilities[2]
                self.totalHapMutType[haplotype, s] = sum(probabilities)
        elif rate==None and isinstance(probabilities, list) and isinstance(haplotype, int) and mutation==None:#DONE
            if len(probabilities)!=4:
                self.Error("#TODO")
            for i in range(4):
                if isinstance(probabilities[i], (int, float)) == False:
                    self.Error("#TODO")
                if probabilities[i]<0:
                    self.Error("#TODO")
            del probabilities[self.calculate_allele(haplotype, mutation)]
            if sum(probabilities) == 0:
                self.Error("#TODO")
            if haplotype<0 or haplotype>=self.hapNum:
                self.Error("There are no such haplotype!")

            for s in range(self.sites):
                self.hapMutType[haplotype, s, 0] = probabilities[0]
                self.hapMutType[haplotype, s, 1] = probabilities[1]
                self.hapMutType[haplotype, s, 2] = probabilities[2]
                self.totalHapMutType[haplotype, s] = sum(probabilities)
        elif isinstance(rate, (int, float)) and probabilities==None and isinstance(haplotype, int) and mutation==None:#DONE
            if rate<0:
                self.Error("#TODO")
            if haplotype<0 or haplotype>=self.hapNum:
                self.Error("There are no such haplotype!")

            for s in range(self.sites):
                self.mRate[haplotype, s] = rate
        elif isinstance(rate, (int, float)) and isinstance(probabilities, list) and haplotype==None and mutation==None:#DONE
            if rate<0:
                self.Error("#TODO")
            if len(probabilities)!=4:
                self.Error("#TODO")
            for i in range(4):
                if isinstance(probabilities[i], (int, float)) == False:
                    self.Error("#TODO")
                if probabilities[i]<0:
                    self.Error("#TODO")

            for hn in range(self.hapNum):
                for s in range(self.sites):
                    probabilities_allele = list(probabilities)
                    del probabilities_allele[self.calculate_allele(hn, s)]
                    if sum(probabilities_allele) == 0:
                        self.Error("#TODO")
                    self.mRate[hn, s] = rate
                    self.hapMutType[hn, s, 0] = probabilities_allele[0]
                    self.hapMutType[hn, s, 1] = probabilities_allele[1]
                    self.hapMutType[hn, s, 2] = probabilities_allele[2]
                    self.totalHapMutType[hn, s] = sum(probabilities_allele)
        elif rate==None and isinstance(probabilities, list) and haplotype==None and mutation==None:#DONE
            if len(probabilities)!=4:
                self.Error("#TODO")
            for i in range(4):
                if isinstance(probabilities[i], (int, float)) == False:
                    self.Error("#TODO")
                if probabilities[i]<0:
                    self.Error("#TODO")

            for hn in range(self.hapNum):
                for s in range(self.sites):
                    probabilities_allele = list(probabilities)
                    del probabilities_allele[self.calculate_allele(hn, s)]
                    if sum(probabilities_allele) == 0:
                        self.Error("#TODO")
                    self.hapMutType[hn, s, 0] = probabilities_allele[0]
                    self.hapMutType[hn, s, 1] = probabilities_allele[1]
                    self.hapMutType[hn, s, 2] = probabilities_allele[2]
                    self.totalHapMutType[hn, s] = sum(probabilities_allele)
        elif isinstance(rate, (int, float)) and probabilities==None and haplotype==None and mutation==None:#DONE
            if rate<0:
                self.Error("#TODO")

            for hn in range(self.hapNum):
                for s in range(self.sites):
                    self.mRate[hn, s] = rate
        else:
            self.Error("#TODO")


    def set_susceptibility_type(self, susceptibility_type, haplotype):
        if isinstance(susceptibility_type, int) == False:
            self.Error("Incorrect value of susceptibility type. Value should be int.")
        if susceptibility_type<0 or susceptibility_type>=self.susNum:
            self.Error("There are no such susceptibility type!")

        if isinstance(haplotype, str):
            haplotypes = self.create_list_haplotypes(haplotype)
            for haplotype in haplotypes:
                self.suscType[haplotype] = susceptibility_type
        elif isinstance(haplotype, int):
            if haplotype<0 or haplotype>=self.hapNum:
                self.Error("There are no such haplotype!")

            self.suscType[haplotype] = susceptibility_type
        elif haplotype == None:
            for hn in range(self.hapNum):
                self.suscType[hn] = susceptibility_type
        else:
            self.Error("Incorrect value of haplotype. Value should be int or str or None.")

    def set_susceptibility(self, rate, haplotype, susceptibility_type):
        if isinstance(rate, (int, float)) == False:
            self.Error("Incorrect value of susceptibility rate. Value should be int or float.")
        if rate<0:
            self.Error("Incorrect value of susceptibility rate. Value should be more or equal 0.")

        if isinstance(haplotype, str) and isinstance(susceptibility_type, int):
            haplotypes = self.create_list_haplotypes(haplotype)
            if susceptibility_type<0 or susceptibility_type>=self.susNum:
                self.Error("There are no such susceptibility type!")

            for haplotype in haplotypes:
                self.susceptibility[haplotype, susceptibility_type] = rate
        elif isinstance(haplotype, int) and isinstance(susceptibility_type, int):
            if haplotype<0 or haplotype>=self.hapNum:
                self.Error("There are no such haplotype!")
            if susceptibility_type<0 or susceptibility_type>=self.susNum:
                self.Error("There are no such susceptibility type!")

            self.susceptibility[haplotype, susceptibility_type] = rate
        elif haplotype==None and isinstance(susceptibility_type, int):
            if susceptibility_type<0 or susceptibility_type>=self.susNum:
                self.Error("There are no such susceptibility type!")

            for hn in range(self.hapNum):
                self.susceptibility[hn, susceptibility_type] = rate
        elif isinstance(haplotype, str) and susceptibility_type==None:
            haplotypes = self.create_list_haplotypes(haplotype)

            for haplotype in haplotypes:
                for sn in range(self.susNum):
                    self.susceptibility[haplotype, sn] = rate
        elif isinstance(haplotype, int) and susceptibility_type==None:
            if haplotype<0 or haplotype>=self.hapNum:
                self.Error("There are no such haplotype!")

            for sn in range(self.susNum):
                self.susceptibility[haplotype, sn] = rate
        elif haplotype==None and susceptibility_type==None:
            for hn in range(self.hapNum):
                for sn in range(self.susNum):
                    self.susceptibility[hn, sn] = rate
        else:
            self.Error("Incorrect value of haplotype or susceptibility rate. Value should be int or None.")

    def set_immunity_transition(self, rate, source, target):
        if isinstance(rate, (int, float)) == False:
            self.Error("Incorrect value of rate. Value should be int or float.")
        if rate<0:
            self.Error("Incorrect value of rate. Value should be more or equal 0.")

        if isinstance(source, int) and isinstance(target, int):
            if source<0 or source>=self.susNum:
                self.Error("There are no such source susceptibility type!")
            if target<0 or target>=self.susNum:
                self.Error("There are no such target susceptibility type!")
            if source==target:
                self.Error("Source and target susceptibility type shouldn't be equal!")

            self.suscepTransition[source, target] = rate
        elif source==None and isinstance(target, int):
            if target<0 or target>=self.susNum:
                self.Error("There are no such target susceptibility type!")

            for sn1 in range(self.susNum):
                if sn1 != target:
                        self.suscepTransition[sn1, target] = rate
        elif isinstance(source, int) and target==None:
            if source<0 or source>=self.susNum:
                self.Error("There are no such source susceptibility type!")

            for sn2 in range(self.susNum):
                if source != sn2:
                    self.suscepTransition[source, sn2] = rate
        elif source==None and target==None:
            for sn1 in range(self.susNum):
                for sn2 in range(self.susNum):
                    if sn1 != sn2:
                        self.suscepTransition[sn1, sn2] = rate
        else:
            self.Error("Incorrect value of source or target susceptibility type. Value should be int or None.")


    def set_population_size(self, amount, population):
        if self.first_simulation == True:
            self.Error("Changing population size isn't available after first simulation!")
        if isinstance(amount, int) == False:
            self.Error("Incorrect value of amount. Value should be int.")
        if amount<0:
            self.Error("#TODO")
        if isinstance(population, int):
            if population<0 or population>=self.popNum:
                self.Error("There are no such population!")

            self.sizes[population] = amount
            self.totalSusceptible[population] = amount
            for sn in range(self.susNum):
                self.susceptible[population, 0] = amount
        elif population==None:
            for pn in range(self.popNum):
                self.sizes[pn] = amount
                self.totalSusceptible[pn] = amount
                for sn in range(self.susNum):
                    self.susceptible[pn, 0] = amount
        else:
            self.Error("Incorrect value of population. Value should be int or None.")

    def set_contact_density(self, value, population):
        if isinstance(value, (int, float)) == False:
            self.Error("Incorrect type of contact density. Value should be int or float.")
        if value<0:
            self.Error("Incorrect value of contact density. Value should be more or equal 0.")

        if isinstance(population, int):
            if population<0 or population>=self.popNum:
                self.Error("There are no such population!")

            self.contactDensity[population] = value
            self.contactDensityBeforeLockdown[population] = value
        elif population == None:
            for pn in range(self.popNum):
                self.contactDensity[pn] = value
                self.contactDensityBeforeLockdown[pn] = value
        else:
            self.Error("Incorrect value of population. Value should be int or None.")

    def set_lockdown(self, parameters, population):
        if isinstance(parameters, list) == False:
            self.Error("#TODO")
        if len(parameters) != 3:
            self.Error("#TODO")
        if parameters[0]<0:
            self.Error("#TODO")
        if parameters[1]<0:
            self.Error("#TODO")
        if parameters[2]<0:
            self.Error("#TODO")

        if isinstance(population, int):
            if population<0 or population>=self.popNum:
                self.Error("There are no such population!")

            self.contactDensityAfterLockdown[population] = parameters[0]
            self.startLD[population] = parameters[1]
            self.endLD[population] = parameters[2]
        elif population == None:
            for pn in range(self.popNum):
                self.contactDensityAfterLockdown[pn] = parameters[0]
                self.startLD[pn] = parameters[1]
                self.endLD[pn] = parameters[2]
        else:
            self.Error("#TODO")

    def set_sampling_multiplier(self, multiplier, population):
        if isinstance(multiplier, (int, float)) == False:
            self.Error("Incorrect type of multiplier. Value should be int or float.")
        if multiplier<0:
            self.Error("Incorrect value of multiplier. Value should be more or equal 0.")

        if isinstance(population, int):
            if population<0 or population>=self.popNum:
                self.Error("There are no such population!")

            self.samplingMultiplier[population] = multiplier
        elif population == None:
            for pn in range(self.popNum):
                self.samplingMultiplier[pn] = multiplier
        else:
            self.Error("Incorrect value of population. Value should be int or None.")

    def set_migration_probability(self, probability, total_probability, source, target):
        if isinstance(probability, float) == True:
            if probability<0 or probability>1:
                self.Error("Incorrect probability. Value should be between 0 and 1!")

            if isinstance(source, int) and isinstance(target, int):
                if source<0 or source>=self.popNum:
                    self.Error("There are no such population!")
                if target<0 or target>=self.popNum:
                    self.Error("There are no such population!")
                if source==target:
                    self.Error("Source and target population shouldn't be equal!")

                self.migrationRates[source, target] = probability
                summa = 0
                for pn in range(self.popNum):
                    summa += self.migrationRates[source, pn]
                if summa > 1:
                    self.Error("#TODO_mp1")
            elif source==None and isinstance(target, int):
                if target<0 or target>=self.popNum:
                    self.Error("There are no such population!")

                for pn1 in range(self.popNum):
                    if pn1 != target:
                        self.migrationRates[pn1, target] = probability
                    summa = 0
                    for pn2 in range(self.popNum):
                        summa += self.migrationRates[pn1, pn2]
                    if summa > 1:
                        self.Error("#TODO_mp2")
            elif isinstance(source, int) and target==None:
                if source<0 or source>=self.popNum:
                    self.Error("There are no such population!")

                for pn2 in range(self.popNum):
                    if source != pn2:
                        self.migrationRates[source, pn2] = probability
                summa = 0
                for pn in range(self.popNum):
                    summa += self.migrationRates[source, pn]
                if summa > 1:
                    self.Error("#TODO_mp3")
            elif source==None and target==None:
                for pn1 in range(self.popNum):
                    for pn2 in range(self.popNum):
                        if pn1 != pn2:
                            self.migrationRates[pn1, pn2] = probability
                    summa = 0
                    for pn2 in range(self.popNum):
                        summa += self.migrationRates[pn1, pn2]
                    if summa > 1:
                        self.Error("#TODO_mp4")
            else:
                self.Error("Incorrect value of population. Value should be int or None.")
        elif isinstance(total_probability, float) == True:
            if total_probability<0 or total_probability>1:
                self.Error("Incorrect total probability. Value should be between 0 and 1!")

            for pn1 in range(self.popNum):
                for pn2 in range(self.popNum):
                    if pn1 != pn2:
                        self.migrationRates[pn1, pn2] = total_probability/(self.popNum-1)
        else:
            self.Error("#TODO_mp5")


    def set_susceptible_individuals(self, amount, source_type, target_type, population):
        if self.first_simulation:
            self.Error('#TODO')
        if isinstance(amount, int) == False:
            self.Error("Incorrect value of amount. Value should be int.")
        if amount<0:
            self.Error("#TODO")
        if isinstance(source_type, int) == False:
            self.Error("Incorrect value of source susceptibility type. Value should be int.")
        if source_type<0 or source_type>=self.susNum:
            self.Error("There are no such susceptibility type!")
        if isinstance(target_type, int) == False:
            self.Error("Incorrect value of target susceptibility type. Value should be int.")
        if target_type<0 or target_type>=self.susNum:
            self.Error("There are no such susceptibility type!")
        if source_type==target_type:
            self.Error("Source and target susceptibility type shouldn't be equal!")

        if isinstance(population, int):
            if population<0 or population>=self.popNum:
                self.Error("There are no such population!")
            if self.susceptible[population, source_type] - amount < 0:
                self.Error('#TODO')
            if self.susceptible[population, target_type] + amount > self.sizes[population]:
                self.Error('#TODO')

            self.susceptible[population, source_type] -= amount
            self.susceptible[population, target_type] += amount
        elif population==None:
            for pn in range(self.popNum):
                if self.susceptible[pn, source_type] - amount < 0:
                    self.Error('#TODO')
                if self.susceptible[pn, target_type] + amount > self.sizes[pn]:
                    self.Error('#TODO')

                self.susceptible[pn, source_type] -= amount
                self.susceptible[pn, target_type] += amount
        else:
            self.Error("Incorrect value of population. Value should be int or None.")

    def set_infected_individuals(self, amount, source_haplotype, target_haplotype, population):
        if self.first_simulation:
            self.Error('#TODO')
        if isinstance(amount, int) == False:
            self.Error("Incorrect value of amount. Value should be int.")
        if amount<0:
            self.Error("#TODO")
        if isinstance(source_haplotype, int) == False:
            self.Error("Incorrect value of source haplotype. Value should be int.")
        if source_haplotype<0 or source_haplotype>=self.hapNum:
            self.Error("There are no such haplotype!")
        if isinstance(target_haplotype, int) == False:
            self.Error("Incorrect value of target haplotype. Value should be int.")
        if target_haplotype<0 or target_haplotype>=self.hapNum:
            self.Error("There are no such haplotype!")
        if source_haplotype==target_haplotype:
            self.Error("Source and target haplotype shouldn't be equal!")

        if isinstance(population, int):
            if population<0 or population>=self.popNum:
                self.Error("There are no such population!")
            if self.liveBranches[population, source_haplotype] - amount < 0:
                self.Error('#TODO')
            if self.liveBranches[population, target_haplotype] + amount > self.sizes[population]:
                self.Error('#TODO')

            self.liveBranches[population, source_haplotype] -= amount
            self.liveBranches[population, target_haplotype] += amount
        elif population==None:
            for pn in range(self.popNum):
                if self.liveBranches[pn, source_haplotype] - amount < 0:
                    self.Error('#TODO')
                if self.liveBranches[pn, target_haplotype] + amount > self.sizes[pn]:
                    self.Error('#TODO')

                self.liveBranches[pn, source_haplotype] -= amount
                self.liveBranches[pn, target_haplotype] += amount
        else:
            self.Error("Incorrect value of population. Value should be int or None.")

    def set_infection(self, amount, source_type, target_haplotype, population):
        if self.first_simulation:
            self.Error('#TODO')
        if isinstance(amount, int) == False:
            self.Error("Incorrect value of amount. Value should be int.")
        if amount<0:
            self.Error("#TODO")
        if isinstance(source_type, int) == False:
            self.Error("Incorrect value of source susceptibility type. Value should be int.")
        if source_type<0 or source_type>=self.susNum:
            self.Error("There are no such susceptibility type!")
        if isinstance(target_haplotype, int) == False:
            self.Error("Incorrect value of target haplotype. Value should be int.")
        if target_haplotype<0 or target_haplotype>=self.hapNum:
            self.Error("There are no such haplotype!")

        if isinstance(population, int):
            if population<0 or population>=self.popNum:
                self.Error("There are no such population!")
            if self.susceptible[population, source_type] - amount < 0:
                self.Error('#TODO')
            if self.liveBranches[population, target_haplotype] + amount > self.sizes[population]:
                self.Error('#TODO')

            self.susceptible[population, source_type] -= amount
            self.totalSusceptible[population] -= amount
            self.liveBranches[population, target_haplotype] += amount
            self.totalInfectious[population] += amount
            self.globalInfectious += amount
        elif population==None:
            for pn in range(self.popNum):
                if self.liveBranches[pn, source_type] - amount < 0:
                    self.Error('#TODO')
                if self.liveBranches[pn, target_haplotype] + amount > self.sizes[pn]:
                    self.Error('#TODO')

                self.susceptible[pn, source_type] -= amount
                self.totalSusceptible[pn] -= amount
                self.liveBranches[pn, target_haplotype] += amount
                self.totalInfectious[pn] += amount
                self.globalInfectious += amount
        else:
            self.Error("Incorrect value of population. Value should be int or None.")


    def set_chain_events(self, name_file):
        if isinstance(name_file, str) == False:
            self.Error('#TODO')
        tokens = np.load(name_file + '.npy')

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


    def output_tree_mutations(self):
        tree = []
        times = []
        for i in range(self.tree.shape[0]):
            tree.append(self.tree[i])
            times.append(self.times[i])
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

        return tree, times, mut, populations

    def output_migrations(self, name_file):
        #with open(name_file + '.mig', 'w') as file:
        file = open(name_file + '.mig', 'w')
        file.write("Node Time Old_population New_population\n")
        for i in range(self.mig.nodeId.size()):
            file.write(str(self.mig.nodeId[i]) + " " + str(self.mig.time[i]) + " " + str(self.mig.oldPop[i]) + " " + str(self.mig.newPop[i]) + "\n")
        file.close()

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

    def output_chain_events(self, name_file):
        np.save(name_file, [self.events.times, self.events.types, self.events.haplotypes, self.events.populations, self.events.newHaplotypes, self.events.newPopulations])

    def output_settings(self, file_template):
        if isinstance(file_template, str) == False:
            self.Error('#TODO')
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


    def get_data_infectious(self, pop, hap, step_num):
        time_points = [i*self.currentTime/step_num for i in range(step_num+1)]
        Date = np.zeros(step_num+1)
        Sample = np.zeros(step_num+1)

        point = 1
        for j in range(self.events.ptr):
            if point != step_num and time_points[point] < self.events.times[j]:
                Date[point+1] = Date[point]
                Sample[point+1] = Sample[point]
                point += 1
            if self.events.populations[j] == pop and self.events.haplotypes[j] == hap:
                if self.events.types[j] == BIRTH:
                    Date[point] += 1
                elif self.events.types[j] == DEATH:
                    Date[point] -= 1
                elif self.events.types[j] == SAMPLING:
                    Date[point] -= 1
                    Sample[point] += 1
                elif self.events.types[j] == MUTATION:
                    Date[point] -= 1
            elif self.events.types[j] == MUTATION and self.events.newHaplotypes[j] == hap and self.events.populations[j] == pop:
                Date[point] += 1
            elif self.events.types[j] == MIGRATION and self.events.newPopulations[j] == pop and self.events.haplotypes[j] == hap:
                Date[point] += 1

        ###
        # time_points = [i*self.currentTime/step_num for i in range(step_num+1)]
        # Date = np.zeros(step_num+1)
        # Sample = np.zeros(step_num+1)
        # Date[step_num-1] = self.liveBranches[pop, hap]

        # point = step_num-1
        # for j in range(self.events.ptr-1, -1, -1):
        #     if point != step_num and time_points[point] > self.events.times[j]:
        #         Date[point-1] = Date[point]
        #         point -= 1
        #     if self.events.populations[j] == pop and self.events.haplotypes[j] == hap:
        #         if self.events.types[j] == BIRTH:
        #             Date[point] -= 1
        #         elif self.events.types[j] == DEATH:
        #             Date[point] += 1
        #         elif self.events.types[j] == SAMPLING:
        #             Date[point] += 1
        #             Sample[point] -= 1
        #         elif self.events.types[j] == MUTATION:
        #             Date[point] += 1
        #     elif self.events.types[j] == MUTATION and self.events.newHaplotypes[j] == hap and self.events.populations[j] == pop:
        #         Date[point] -= 1
        #     elif self.events.types[j] == MIGRATION and self.events.newPopulations[j] == pop and self.events.haplotypes[j] == hap:
        #         Date[point] -= 1

        # for i in range(self.events.ptr):
        #     if point != step_num and time_points[point] < self.events.times[i]:
        #         Sample[point+1] = Sample[point]
        #         point += 1
        #     if self.events.populations[i] == pop and self.events.haplotypes[i] == hap and self.events.types[i] == SAMPLING:
        #         Sample[point] += 1
        ###

        Lockdowns = []
        for i in range(self.loc.times.size()):
            if self.loc.populationsId[i] == pop:
                Lockdowns.append([self.loc.states[i], self.loc.times[i]])

        return Date, Sample, time_points, Lockdowns

    def get_data_susceptible(self, pop, sus, step_num):
        time_points = [i*self.currentTime/step_num for i in range(step_num+1)]
        Date = np.zeros(step_num+1)
        if sus == 0:
            Date[0] = self.sizes[pop]
            Date[1] = self.sizes[pop]

        point = 1
        for j in range(self.events.ptr):
            if point != step_num and time_points[point] < self.events.times[j]:
                Date[point+1] = Date[point]
                point += 1
            if self.events.populations[j] == pop and self.events.newHaplotypes[j] == sus:
                if self.events.types[j] == BIRTH:
                    Date[point] -= 1
                elif self.events.types[j] == DEATH:
                    Date[point] += 1
                elif self.events.types[j] == SAMPLING:
                    Date[point] += 1
                elif self.events.types[j] == SUSCCHANGE:
                    Date[point] += 1
            elif self.events.types[j] == SUSCCHANGE and self.events.haplotypes[j] == sus and self.events.populations[j] == pop:
                Date[point] -= 1
            elif self.events.types[j] == MIGRATION and self.events.newPopulations[j] == pop and self.events.newHaplotypes[j] == sus:
                Date[point] -= 1

        Lockdowns = []
        for i in range(self.loc.times.size()):
            if self.loc.populationsId[i] == pop:
                Lockdowns.append([self.loc.states[i], self.loc.times[i]])

        return Date, time_points, Lockdowns


    def Stats(self):
        print("Seed: ", self.seed)
        print("Number of samples:", self.sCounter)
        print("Total number of iterations: ", self.events.ptr)
        print("Size events: ", self.events.size)
        print("Current time: ", self.currentTime)

    def Debug(self):
        print("Parameters")
        print("Mutations number(const): ", self.sites)
        print("Haplotypes number(const): ", self.hapNum)
        print("Populations number(const): ", self.popNum)
        print("Susceptible number(const): ", self.susNum)
        print("Samping probability(const): ", self.sampling_probability)
        print("Strong migration(const): ", self.strong_migration)

        print("Birth counter(mutable): ", self.bCounter)
        print("Death counter(mutable): ", self.dCounter)
        print("Sampling counter(mutable): ", self.sCounter)
        print("Mutation counter(mutable): ", self.mCounter)
        print("Immunity transition counter(mutable):", self.iCounter)
        print("swapLockdown(mutable): ", self.swapLockdown)
        print("Migration plus(mutable): ", self.migPlus)
        print("Migration non plus(mutable): ", self.migNonPlus)
        print("Migration counter(mutable): ", self.migCounter)
        print("globalInfectious(mutable): ", self.globalInfectious)

        print("Current time(mutable): ", self.currentTime)
        print("Seed(const): ", self.seed)
        print("Max effective birth(const): ", self.maxEffectiveBirth)
        print("Total rate(mutable): ", self.totalRate)
        print("Total migration rate(mutable): ", self.totalMigrationRate)
        print("Random number(mutable): ", self.rn)

        print("Susceptible type(const): ", sep=" ", end="")
        for hn in range(self.hapNum):
            print(self.suscType[hn], end=" ")
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
        print("Total mutation rate(const): ", sep="", end="")
        for hn in range(self.hapNum):
            print(self.tmRate[hn], end=" ")
        print()
        print("maxEffectiveBirthMigration(const): ", sep="", end="")
        for pn in range(self.popNum):
            print(self.maxEffectiveBirthMigration[pn], end=" ")
        print()
        print("Max susceptibility(const): ", end=" ")
        for hn in range(self.hapNum):
            print(self.maxSusceptibility[hn], end=" ")
        print()
        print("suscepCumulTransition(const): ", sep="", end="")
        for sn in range(self.susNum):
            print(self.suscepCumulTransition[sn], end=" ")
        print()
        print("infectPopRate(mutable): ", sep="", end="")
        for pn in range(self.popNum):
            print(self.infectPopRate[pn], end=" ")
        print()
        print("immunePopRate(mutable): ", sep="", end="")
        for pn in range(self.popNum):
            print(self.immunePopRate[pn], end=" ")
        print()
        print("Population rate(mutable): ", sep="", end="")
        for pn in range(self.popNum):
            print(self.popRate[pn], end=" ")
        print()
        print("Migration population rate(mutable): ", sep="", end="")
        for pn in range(self.popNum):
            print(self.migPopRate[pn], end=" ")
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
        print("effectiveSizes(const): ", end=" ")
        for pn in range(self.popNum):
            print(self.effectiveSizes[pn], end=" ")
        print()
        print("contact density(const): ", end=" ")
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
        print("max effective migration(const): ", end=" ")
        for pn in range(self.popNum):
            print(self.maxEffectiveMigration[pn], end=" ")
        print()

        print("Mutation rate(const)----")
        for hn in range(self.hapNum):
            for s in range(self.sites):
                print(self.mRate[hn, s], end=" ")
            print()
        print()
        print("Total haplotype mutation type(const)----")
        for hn in range(self.hapNum):
            for s in range(self.sites):
                print(self.totalHapMutType[hn, s], end=" ")
            print()
        print()
        print("Total event haplotype population rate(mutable)----")
        for pn in range(self.popNum):
            for hn in range(self.hapNum):
                print(self.tEventHapPopRate[pn, hn], end=" ")
            print()
        print()
        print("susceptible(mutable)----")
        for pn in range(self.popNum):
            for sn in range(self.susNum):
                print(self.susceptible[pn, sn], end=" ")
            print()
        print()
        print("liveBranches(mutable)----")
        for pn in range(self.popNum):
            for hn in range(self.hapNum):
                print(self.liveBranches[pn, hn], end=" ")
            print()
        print()
        print("Population model - migration rates(const)----")
        for pn1 in range(self.popNum):
            for pn2 in range(self.popNum):
                print(self.migrationRates[pn1, pn2], end=" ")
            print()
        print()
        print("Population model - effective migration(const)----")
        for pn1 in range(self.popNum):
            for pn2 in range(self.popNum):
                print(self.effectiveMigration[pn1, pn2], end=" ")
            print()
        print()

        print("Susceptibility(const)----")
        for hn in range(self.hapNum):
            for sn in range(self.susNum):
                print(self.susceptibility[hn, sn], end=" ")
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
        print("Haplotypes populations rates(mutable)----")
        for pn in range(self.popNum):
            for hn in range(self.hapNum):
                print(self.hapPopRate[pn, hn], end=" ")
            print()
        print()

        print("hapMutType(const)----")
        for hn in range(self.hapNum):
            for s in range(self.sites):
                for i in range(3):
                    print(self.hapMutType[hn, s, i], end=" ")
                print()
            print()
        print()
        print("Event haplotypes populations rate(mutable)----")
        for pn in range(self.popNum):
            for hn in range(self.hapNum):
                for i in range(4):
                    print(self.eventHapPopRate[pn, hn, i], end=" ")
                print()
            print()
        print()
        print("Susceptible haplotypes populations rate(mutable)----")
        for pn in range(self.popNum):
            for hn in range(self.hapNum):
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
    cdef double MigrationPropensity(self, Py_ssize_t s, Py_ssize_t r, Py_ssize_t i, Py_ssize_t h):#TODO: use EffevtiveMigration
        cdef double prop
        prop = self.migrationRates[r, s]*self.contactDensity[s]/self.effectiveSizes[s]
        prop += self.migrationRates[s, r]*self.contactDensity[r]/self.effectiveSizes[r]
        prop *= self.susceptible[r, i]*self.liveBranches[s, h]*self.bRate[h]*self.susceptibility[h, i]*self.migrationRates[s, s]
        return prop

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef Py_ssize_t PropensitiesNumber(self):
        cdef Py_ssize_t prop_num

        #Migrations
        prop_num = self.popNum*(self.popNum-1)*self.hapNum*self.susNum
        prop_num += self.popNum*self.susNum*(self.susNum-1)
        prop_num += self.popNum*self.hapNum*(1+1+self.sites*3+self.susNum)
        return prop_num

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void Propensities(self):
        cdef double prop

        #Migrations
        for s in range(self.popNum):
            for r in range(self.popNum):
                if s == r:
                    continue
                for i in range(self.susNum):
                    for h in range(self.hapNum):
                        self.PropensitiesMigr[s, r, i, h] = self.MigrationPropensity(s, r, i, h)


        for s in range(self.popNum):
            #Susceptibility transition
            for i in range(self.susNum):
                for j in range(self.susNum):
                    if i == j:
                        continue
                    self.PropensitiesSuscep[s, i, j] = self.suscepTransition[i, j]*self.susceptibility[h, i]

            #Infectious-realted event
            for h in range(self.hapNum):
                #Recovery
                self.PropensitiesRecovery[s, h] = self.dRate[h]*self.liveBranches[s, h]

                #Sampling
                self.PropensitiesSampling[s, h] = self.sRate[h]*self.liveBranches[s, h]

                #Mutation
                for site in range(self.sites):
                    for i in range(3):
                        #ht = self.Mutate(h, site, i)
                        self.PropensitiesMutatations[s, h, site, i] = self.mRate[h, site]*self.hapMutType[h, site, i]/sum(self.hapMutType[h, site])*self.liveBranches[s, h]#FIXME Normalize self.hapMutType !!!
                #Transmission
                for i in range(self.susNum):
                    prop = self.bRate[h]*self.contactDensity[s]
                    prop *= self.migrationRates[s, s]*self.migrationRates[s, s]
                    prop *= self.susceptible[s, i]*self.liveBranches[s, h]/self.effectiveSizes[s]
                    #print(self.effectiveSizes[s])
                    self.PropensitiesTransmission[s, h, i] = prop
                    #print("migr=", self.migrationRates[s, s], "  prop", prop)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef Py_ssize_t GenerateEvents_tau(self, float tau, Py_ssize_t debug=0):
        cdef Py_ssize_t event_num
        #self.infectiousDelta = np.zeros((self.popNum, self.hapNum), dtype=np.int64)
        #self.susceptibleDelta = np.zeros((self.popNum, self.susNum), dtype=np.int64)
        for s in range(self.popNum):
            for i in range(self.susNum):
                self.susceptibleDelta[s, i] = 0
            for h in range(self.hapNum):
                self.infectiousDelta[s, h] = 0
        #Migrations
        for s in range(self.popNum):
            for r in range(self.popNum):
                if s == r:
                    continue
                for i in range(self.susNum):
                    for h in range(self.hapNum):
                        event_num = self.DrawEventsNum(self.PropensitiesMigr[s, r, i, h], tau)
                        self.eventsMigr[s, r, i, h] = event_num
                        self.susceptibleDelta[r, i] -= event_num


        for s in range(self.popNum):
            #Susceptibility transition
            for i in range(self.susNum):
                for j in range(self.susNum):
                    if i == j:
                        continue
                    event_num = self.DrawEventsNum(self.PropensitiesSuscep[s, i, j], tau)
                    self.eventsSuscep[s, i, j] = event_num
                    self.susceptibleDelta[s, i] -= event_num
                    #self.susceptibleDelta[s, j] += event_num

            #Infectious-realted event
            for h in range(self.hapNum):
                #Recovery
                event_num = self.DrawEventsNum(self.PropensitiesRecovery[s, h], tau)
                self.eventsRecovery[s, h] = event_num
                self.infectiousDelta[s, h] -= event_num
                #self.susceptibleDelta[s, self.suscType[h]] += event_num

                #Sampling
                event_num = self.DrawEventsNum(self.PropensitiesSampling[s, h], tau)
                self.eventsSampling[s, h] = event_num
                self.infectiousDelta[s, h] -= event_num
                #self.susceptibleDelta[s, self.suscType[h]] += event_num

                #Mutation
                for site in range(self.sites):
                    for i in range(3):
                        event_num = self.DrawEventsNum(self.PropensitiesMutatations[s, h, site, i], tau)
                        self.eventsMutatations[s, h, site, i] = event_num
                        self.infectiousDelta[s, h] -= event_num
                        #if event_num != 0:
                        #    print(self.PropensitiesMutatations[s, h, site, i]*tau, "  |  ", event_num)
                        #self.infectiousDelta[s, self.Mutate(h, site, i)] += event_num
                #Transmission
                for i in range(self.susNum):
                    event_num = self.DrawEventsNum( self.PropensitiesTransmission[s, h, i], tau)
                    self.eventsTransmission[s, h, i] = event_num
                    self.susceptibleDelta[s, i] -= event_num
                    #self.infectiousDelta[s, h] += event_num
        self.countsPerStep = 0
        for s in range(self.popNum):
            for i in range(self.susNum):
                self.countsPerStep += self.susceptibleDelta[s, i]
                if self.susceptibleDelta[s, i] + self.susceptible[s, i] < 0:
                    return 0
            for h in range(self.hapNum):
                self.countsPerStep += self.infectiousDelta[s, h]
                if self.infectiousDelta[s, h] + self.liveBranches[s, h] < 0:
                    return 0
        return 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void UpdateCompartmentCounts_tau(self):
        cdef Py_ssize_t event_num

        #Migrations
        for s in range(self.popNum):
            for r in range(self.popNum):
                if s == r:
                    continue
                for i in range(self.susNum):
                    for h in range(self.hapNum):
                        event_num = self.eventsMigr[s, r, i, h]
                        self.NewInfections(event_num, r, i, h)
                        self.multievents.AddEvents(event_num, self.currentTime, MIGRATION, h, s, i, r)
                        self.migPlus += event_num

        for s in range(self.popNum):
            #Susceptibility transition
            for i in range(self.susNum):
                for j in range(self.susNum):
                    if i == j:
                        continue
                    event_num = self.eventsSuscep[s, i, j]
                    self.susceptible[s, i] -= event_num
                    self.susceptible[s, j] += event_num
                    self.multievents.AddEvents(event_num, self.currentTime, SUSCCHANGE, i, s, j, 0)
                    self.iCounter += event_num

            #Infectious-realted event
            for h in range(self.hapNum):
                #Recovery
                event_num = self.eventsRecovery[s, h]
                self.NewRecoveries(event_num, s, self.suscType[h], h)
                self.multievents.AddEvents(event_num, self.currentTime, DEATH, h, s, self.suscType[h], 0)
                self.dCounter += event_num

                #Sampling
                event_num = self.eventsSampling[s, h]
                self.NewRecoveries(event_num, s, self.suscType[h], h)
                self.multievents.AddEvents(event_num, self.currentTime, SAMPLING, h, s, self.suscType[h], 0)
                self.sCounter += event_num

                #Mutation
                for site in range(self.sites):
                    for i in range(3):
                        ht = self.Mutate(h, site, i)
                        event_num = self.eventsMutatations[s, h, site, i]
                        self.liveBranches[s, ht] += event_num
                        self.liveBranches[s, h] -= event_num
                        self.multievents.AddEvents(event_num, self.currentTime, MUTATION, h, s, ht, 0)
                        self.mCounter += event_num
                #Transmission
                for i in range(self.susNum):
                    event_num = self.eventsTransmission[s, h, i]
                    self.NewInfections(event_num, s, i, h)
                    self.multievents.AddEvents(event_num, self.currentTime, BIRTH, h, s, i, 0)
                    self.bCounter += event_num

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef inline Py_ssize_t DrawEventsNum(self, double prop, double tau):
        cdef Py_ssize_t n
        #n = np.random.poisson(prop*tau)
        n = random_poisson(self.rndm.rng, lam=prop*tau)
        #if not n ==0:
        #    print(n)
        return n

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void ChooseTau(self, double epsilon=0.03):
        cdef float tmp
        #self.infectiousAuxTau = np.zeros((self.popNum, self.hapNum, 2), dtype=float)
        #self.susceptibleAuxTau = np.zeros((self.popNum, self.susNum, 2), dtype=float)
        for s in range(self.popNum):
            for i in range(self.susNum):
                self.susceptibleAuxTau[s, i, 0] = 0.0
            for h in range(self.hapNum):
                self.infectiousAuxTau[s, h, 0] = 0.0
        #Migrations
        for s in range(self.popNum):
            for r in range(self.popNum):
                if s == r:
                    continue
                for i in range(self.susNum):
                    for h in range(self.hapNum):
                        self.infectiousAuxTau[r, h, 0] += self.PropensitiesMigr[s, r, i, h]
                        self.susceptibleAuxTau[r, i, 0] -= self.PropensitiesMigr[s, r, i, h]

        for s in range(self.popNum):
            #Susceptibility transition
            for i in range(self.susNum):
                for j in range(self.susNum):
                    if i == j:
                        continue
                    self.susceptibleAuxTau[s, i, 0] -= self.PropensitiesSuscep[s, i, j]
                    self.susceptibleAuxTau[s, j, 0] += self.PropensitiesSuscep[s, i, j]

            #Infectious-realted event
            for h in range(self.hapNum):
                #Recovery
                self.infectiousAuxTau[s, h, 0] -= self.PropensitiesRecovery[s, h]
                self.susceptibleAuxTau[r, self.suscType[h], 0] += self.PropensitiesRecovery[s, h]

                #Sampling
                self.infectiousAuxTau[s, h, 0] -= self.PropensitiesSampling[s, h]
                self.susceptibleAuxTau[r, self.suscType[h], 0] += self.PropensitiesSampling[s, h]
                #Mutation
                for site in range(self.sites):
                    for i in range(3):
                        ht = self.Mutate(h, site, i)
                        self.infectiousAuxTau[s, h, 0] -= self.PropensitiesMutatations[s, h, site, i]
                        self.infectiousAuxTau[s, ht, 0] += self.PropensitiesMutatations[s, h, site, i]
                #Transmission
                for i in range(self.susNum):
                    self.infectiousAuxTau[s, h, 0] += self.PropensitiesTransmission[s, h, i]
                    self.susceptibleAuxTau[s, i, 0] -= self.PropensitiesTransmission[s, h, i]
                    #print("migr=", self.migrationRates[s, s], "  prop", prop)
        self.tau_l = 1.0
        for s in range(self.popNum):
            for h in range(self.hapNum):
                if abs(self.infectiousAuxTau[s, h, 0]) < 1e-8:
                    continue
                tmp = max(epsilon*self.liveBranches[s,h]/2.0,1.0)/abs(self.infectiousAuxTau[s, h, 0])
                if tmp < self.tau_l:
                    self.tau_l = tmp
            for i in range(self.susNum):
                if abs(self.susceptibleAuxTau[s, i, 0]) < 1e-8:
                    continue
                tmp = max(epsilon*self.susceptible[s,i]/2.0,1.0)/abs(self.susceptibleAuxTau[s, i, 0])
                if tmp < self.tau_l:
                    self.tau_l = tmp

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef Py_ssize_t Mutate(self, Py_ssize_t h, Py_ssize_t site, Py_ssize_t DS):
        cdef:
            Py_ssize_t mi, digit4, AS
        digit4 = 4**site
        AS = int(floor(h/digit4) % 4)
        if DS >= AS:
            DS += 1
        return h + (DS-AS)*digit4

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef void SimulatePopulation_tau(self, Py_ssize_t iterations, Py_ssize_t sampleSize=-1, float time=-1):
        cdef Py_ssize_t pi, propNum, success
        propNum = self.PropensitiesNumber()
        if self.globalInfectious == 0:
            self.FirstInfection()
        self.multievents.CreateEvents(iterations*propNum)
        self.events.CreateEvents(iterations)
        self.UpdateAllRates()
        while (self.multievents.ptr<self.multievents.size and (sampleSize==-1 or self.sCounter<sampleSize) and (time==-1 or self.currentTime<time)):
            self.Propensities()
            self.ChooseTau()
            success = 0
            while success == 0:
                success = self.GenerateEvents_tau(self.tau_l)
                self.tau_l /= 2
            self.tau_l *= 2
            self.currentTime += self.tau_l

            self.UpdateCompartmentCounts_tau()
            self.events.AddEvent(self.currentTime, MULTITYPE, self.multievents.ptr-propNum, self.multievents.ptr, 0, 0)
            if self.globalInfectious == 0:
                break
            for s in range(self.popNum):
                self.CheckLockdown(s)
        print("Number of iterations: ", int(self.multievents.ptr/propNum))
        print("Simulation model time: ", self.currentTime)


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
                print("Recovery ", s, h, self.PropensitiesRecovery[s, h])
                print("Sampling ", s, h, self.PropensitiesSampling[s, h])

                for site in range(self.sites):
                    for i in range(3):
                        #ht = self.Mutate(h, site, i)
                        print("Mutation", s, h, site, i, self.PropensitiesMutatations[s, h, site, i])
                #Transmission
                for i in range(self.susNum):
                    print("Transmission", s, h, i, self.PropensitiesTransmission[s, h, i])
                    #print("migr=", self.migrationRates[s, s], "  prop", prop)

