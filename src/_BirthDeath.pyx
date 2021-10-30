# cython: language_level=3
# cython: initializedcheck = False
# distutils: language = c++

cimport cython

from libc.math cimport log, floor, abs
from libcpp.vector cimport vector
from mc_lib.rndm cimport RndmWrapper

from prettytable import PrettyTable
import numpy as np
import sys
import os

include "fast_choose.pxi"
include "models.pxi"


def print_err(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# use named constants for event types
DEF BIRTH = 0
DEF DEATH = 1
DEF SAMPLING = 2
DEF MUTATION = 3
DEF SUSCCHANGE = 4
DEF MIGRATION = 5


cdef class Event:
    cdef:
        Py_ssize_t type_, haplotype, population, newHaplotype, newPopulation
        double time

    def __init__(self, double time, Py_ssize_t type_, Py_ssize_t haplotype, Py_ssize_t population, Py_ssize_t newHaplotype, Py_ssize_t newPopulation):
        self.time = time
        self.type_ = type_
        self.haplotype = haplotype
        self.population = population
        self.newHaplotype = newHaplotype
        self.newPopulation = newPopulation


cdef class Events:
    cdef:
        Py_ssize_t size, ptr

        Py_ssize_t[::1] types, haplotypes, populations, newHaplotypes, newPopulations
        double[::1] times

    def __init__(self):
        self.size = 0
        self.ptr = 0#pointer to the first empty cell

        self.times = np.zeros(1, dtype=float)
        self.types = np.zeros(1, dtype=int)
        self.haplotypes = np.zeros(1, dtype=int)
        self.populations = np.zeros(1, dtype=int)
        self.newHaplotypes = np.zeros(1, dtype=int)
        self.newPopulations = np.zeros(1, dtype=int)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void AddEvent(self, double time_, Py_ssize_t type_, Py_ssize_t haplotype, Py_ssize_t population, Py_ssize_t newHaplotype, Py_ssize_t newPopulation):
        self.times[ self.ptr ] = time_
        self.types[ self.ptr ] = type_
        self.haplotypes[ self.ptr ] = haplotype
        self.populations[ self.ptr ] = population
        self.newHaplotypes[ self.ptr ] = newHaplotype
        self.newPopulations[ self.ptr ] = newPopulation
        self.ptr += 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef Event GetEvent(self, Py_ssize_t e_id):
        ev = Event( self.times[ e_id ], self.types[ e_id ], self.haplotypes[ e_id ], self.populations[ e_id ], self.newHaplotypes[ e_id ], self.newPopulations[ e_id ])
        return( ev )

    cdef void CreateEvents(self, Py_ssize_t iterations):
        self.size = iterations + self.ptr
        self.times = np.resize(self.times, self.size)
        self.types = np.resize(self.types, self.size)
        self.haplotypes = np.resize(self.haplotypes, self.size)
        self.populations = np.resize(self.populations, self.size)
        self.newHaplotypes = np.resize(self.newHaplotypes, self.size)
        self.newPopulations = np.resize(self.newPopulations, self.size)


#pi - population ID, pn - popoulation number, spi - source population ID, tpi - target population ID
#hi - haplotype ID, hn - haplotype number, nhi - new haplotype number
#si - susceptibility ID, sn - susceptibility number, ssi - source susceptibility ID, tsi - target susceptibility ID
cdef class BirthDeathModel:
    cdef:
        RndmWrapper rndm

        Py_ssize_t bCounter, dCounter, sCounter, mCounter, iCounter, popNum, sites, hapNum, susNum
        double currentTime, totalLen, rn, totalRate, maxEffectiveBirth, totalMigrationRate

        Events events
        PopulationModel pm
        Mutations mut

        long[::1] tree, suscType
        long[:,::1] liveBranches

        double[::1] bRate, dRate, sRate, tmRate, migPopRate, popRate, times, pm_maxEffectiveMigration, maxSusceptibility, immunePopRate, infectPopRate, sourceSuscepTransition, suscepCumulTransition
        double[:,::1] tEventHapPopRate, hapPopRate, mRate, susceptibility, totalHapMutType, suscepTransition, immuneSourcePopRate
        double[:,:,::1] eventHapPopRate, susceptHapPopRate, hapMutType

    def __init__(self, sites_number, population_sizes, susceptibility_types, seed):
        self.currentTime = 0.0
        self.bCounter = 0
        self.dCounter = 0
        self.sCounter = 0
        self.mCounter = 0
        self.iCounter = 0
        self.events = Events()
        self.mut = Mutations()

        self.sites = sites_number
        self.hapNum = 4**self.sites
        self.susNum = susceptibility_types
        self.popNum = len(population_sizes)

        self.pm = PopulationModel(population_sizes, self.susNum, self.hapNum)
        self.pm.NewInfection(0, 0, 0)
 
        self.susceptibility = np.zeros((self.hapNum, self.susNum), dtype=float)
        self.suscType = np.zeros(self.hapNum, dtype=np.int64)
        for hn in range(self.hapNum):
            self.susceptibility[hn][0] = 1.0
        self.maxSusceptibility = np.zeros(self.hapNum, dtype=float)
        self.susceptHapPopRate = np.zeros((self.popNum, self.hapNum, self.susNum), dtype=float)
        self.suscepTransition = np.zeros( (self.susNum, self.susNum), dtype=float)

        #Set rates
        self.bRate = np.zeros(self.hapNum, dtype=float)
        self.dRate = np.zeros(self.hapNum, dtype=float)
        self.sRate = np.zeros(self.hapNum, dtype=float)
        self.mRate = np.zeros((self.hapNum, self.sites), dtype=float)
        self.tmRate = np.zeros(self.hapNum, dtype=float)
        self.hapMutType = np.zeros((self.hapNum, self.sites, 3), dtype=float)
        self.totalHapMutType = np.zeros((self.hapNum, self.sites), dtype=float)
        self.eventHapPopRate = np.zeros((self.popNum, self.hapNum, 4), dtype=float)
        self.tEventHapPopRate = np.zeros((self.popNum, self.hapNum), dtype=float)
        self.hapPopRate = np.zeros((self.popNum, self.hapNum), dtype=float)
        self.infectPopRate = np.zeros(self.popNum, dtype=float)
        self.suscepCumulTransition = np.zeros(self.susNum, dtype=float)
        self.immuneSourcePopRate = np.zeros((self.popNum, self.susNum), dtype=float)
        self.immunePopRate = np.zeros(self.popNum, dtype=float)
        self.popRate = np.zeros(self.popNum, dtype=float)
        self.migPopRate = np.zeros(self.popNum, dtype=float)
        self.totalMigrationRate = 0.0
        self.totalRate = 0
        self.maxEffectiveBirth = 0.0

        for hn in range(self.hapNum):
            self.bRate[hn] = 2
            self.dRate[hn] = 1
            self.sRate[hn] = 0.01
            self.tmRate[hn] = 0.01 * self.sites
            for s in range(self.sites):
                self.mRate[hn, s] = 0.01
                self.hapMutType[hn, s, 0] = 1
                self.hapMutType[hn, s, 1] = 1
                self.hapMutType[hn, s, 2] = 1
                self.totalHapMutType[hn, s] = 3

        for hn in range(self.hapNum):
            self.tmRate[hn] = 0
            for s in range(self.sites):
                self.totalHapMutType[hn, s] = 0
                self.totalHapMutType[hn, s] += self.hapMutType[hn, s, 0]
                self.totalHapMutType[hn, s] += self.hapMutType[hn, s, 1]
                self.totalHapMutType[hn, s] += self.hapMutType[hn, s, 2]
                self.tmRate[hn] += self.mRate[hn, s]

        for pn in range(self.popNum):
            for hn in range(self.hapNum):
                self.eventHapPopRate[pn, hn, 0] = self.BirthRate(pn, hn)
                self.eventHapPopRate[pn, hn, 1] = self.dRate[hn]
                self.eventHapPopRate[pn, hn, 2] = self.sRate[hn]*self.pm.samplingMultiplier[pn]
                self.eventHapPopRate[pn, hn, 3] = self.tmRate[hn]
                for i in range(4):
                    self.tEventHapPopRate[pn, hn] += self.eventHapPopRate[pn, hn, i]
                self.hapPopRate[pn, hn] = self.tEventHapPopRate[pn, hn]*self.pm.liveBranches[pn, hn]
                self.infectPopRate[pn] += self.hapPopRate[pn, hn]

        for sn1 in range(self.susNum):
            for sn2 in range(self.susNum):
                self.suscepCumulTransition[sn1] += self.suscepTransition[sn1, sn2]

        for pn in range(self.popNum):
            self.popRate[pn] = 0
            self.immunePopRate[pn] = 0
            for sn in range(self.susNum):
                self.immuneSourcePopRate[pn, sn] = 0
                self.immuneSourcePopRate[pn, sn] += self.suscepCumulTransition[sn]*self.pm.susceptible[pn, sn]
                self.immunePopRate[pn] += self.immuneSourcePopRate[pn, sn]
            self.popRate[pn] = self.infectPopRate[pn] + self.immunePopRate[pn]
            self.totalRate += self.popRate[pn]

        for hn in range(self.hapNum):
            self.maxSusceptibility[hn] = 0.0
            for sn in range(self.susNum):
                if self.susceptibility[hn, sn] > self.maxSusceptibility[hn]:
                    self.maxSusceptibility[hn] = self.susceptibility[hn, sn]
        for hn in range(self.hapNum):
            if self.maxEffectiveBirth < self.bRate[hn]*self.maxSusceptibility[hn]:
                self.maxEffectiveBirth = self.bRate[hn]*self.maxSusceptibility[hn]
        
        for pn in range(self.popNum):
            self.migPopRate[pn] = self.pm.maxEffectiveMigration[pn]*self.maxEffectiveBirth*self.pm.totalSusceptible[pn]*(self.pm.globalInfectious-self.pm.totalInfectious[pn])
            self.totalMigrationRate += self.migPopRate[pn]

        #Set random generator
        self.rndm = RndmWrapper(seed=(seed, 0))

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef inline double BirthRate(self, Py_ssize_t pi, Py_ssize_t hi):
        cdef double ws = 0.0
        for sn in range(self.susNum):
            self.susceptHapPopRate[pi, hi, sn] = self.pm.susceptible[pi, sn]*self.susceptibility[hi, sn]
            ws += self.susceptHapPopRate[pi, hi, sn]

        return self.bRate[hi]*ws/self.pm.sizes[pi]*self.pm.contactDensity[pi]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef void SimulatePopulation(self, Py_ssize_t iterations, Py_ssize_t sampleSize, float time=-1):
        cdef Py_ssize_t pi
        self.events.CreateEvents(iterations)
        #self.totalLen = 0.0
        while (self.events.ptr<=self.events.size and self.sCounter<=sampleSize and (time==-1 or self.currentTime<time)):
            self.SampleTime()
            pi = self.GenerateEvent()
            if self.totalRate == 0.0 or self.pm.globalInfectious == 0:
                break
            self.CheckLockdown(pi)

    cdef void SetRates(self):
        self.pm.SetEffectiveMigration()

        for hn in range(self.hapNum):
            self.tmRate[hn] = 0
            for s in range(self.sites):
                self.tmRate[hn] += self.mRate[hn, s]

        for pn in range(self.popNum):
            self.infectPopRate[pn] = 0.0
            for hn in range(self.hapNum):
                self.eventHapPopRate[pn, hn, 0] = self.BirthRate(pn, hn)
                self.eventHapPopRate[pn, hn, 1] = self.dRate[hn]
                self.eventHapPopRate[pn, hn, 2] = self.pm.samplingMultiplier[pn] * self.sRate[hn]
                self.eventHapPopRate[pn, hn, 3] = self.tmRate[hn]
                self.tEventHapPopRate[pn, hn] = 0
                for i in range(4):
                    self.tEventHapPopRate[pn, hn] += self.eventHapPopRate[pn, hn, i]
                self.hapPopRate[pn, hn] = self.tEventHapPopRate[pn, hn]*self.pm.liveBranches[pn, hn]
                self.infectPopRate[pn] += self.hapPopRate[pn, hn]

        for sn1 in range(self.susNum):
            self.suscepCumulTransition[sn1] = 0
            for sn2 in range(self.susNum):
                self.suscepCumulTransition[sn1] += self.suscepTransition[sn1, sn2]
        
        for pn in range(self.popNum):
            self.popRate[pn] = 0.0
        self.totalRate = 0.0
        for pn in range(self.popNum):
            self.immunePopRate[pn] = 0.0
            for sn in range(self.susNum):
                self.immuneSourcePopRate[pn, sn] = self.suscepCumulTransition[sn]*self.pm.susceptible[pn, sn]
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

        for pn in range(self.popNum):
            self.migPopRate[pn] = 0.0
        self.totalMigrationRate = 0.0
        for pn in range(self.popNum):
            self.migPopRate[pn] = self.pm.maxEffectiveMigration[pn]*self.maxEffectiveBirth*self.pm.totalSusceptible[pn]*(self.pm.globalInfectious-self.pm.totalInfectious[pn])
            self.totalMigrationRate += self.migPopRate[pn]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef inline void SampleTime(self):
        cdef double tau = - log(self.rndm.uniform()) / self.totalRate
        self.currentTime += tau
        #self.totalLen += tau*self.pm.globalInfectious

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef Py_ssize_t GenerateEvent(self):
        cdef:
            Py_ssize_t pi, hi, ei

        self.rn = self.rndm.uniform()
        if self.totalRate > self.rn * (self.totalRate + self.totalMigrationRate):

            self.rn = self.rn * (self.totalRate + self.totalMigrationRate) / self.totalRate
            pi, self.rn = fastChoose1( self.popRate, self.totalRate, self.rn)
            if self.immunePopRate[pi] > self.rn * self.popRate[pi]:
                self.rn = self.rn * self.popRate[pi] / self.immunePopRate[pi]
                self.ImmunityTransition(pi)
            else:
                self.rn = (self.rn * self.popRate[pi] - self.immunePopRate[pi]) / self.infectPopRate[pi]
                hi, self.rn = fastChoose1( self.hapPopRate[pi], self.infectPopRate[pi], self.rn)
                ei, self.rn = fastChoose1( self.eventHapPopRate[pi, hi], self.tEventHapPopRate[pi, hi], self.rn)
                if ei == BIRTH:
                    self.Birth(pi, hi)
                elif ei == DEATH:
                    self.Death(pi, hi)
                elif ei == SAMPLING:
                    self.Sampling(pi, hi)
                else:
                    self.Mutation(pi, hi)
        else:
            self.rn = (self.rn * (self.totalRate + self.totalMigrationRate) - self.totalRate) / self.totalMigrationRate
            pi = self.GenerateMigration()
        return pi

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void ImmunityTransition(self, Py_ssize_t pi):
        cdef:
            Py_ssize_t si, ti
        si, self.rn = fastChoose1( self.immuneSourcePopRate[pi], self.immunePopRate[pi], self.rn)
        ti, self.rn = fastChoose1( self.suscepTransition[si], self.suscepCumulTransition[si], self.rn)

        self.pm.susceptible[pi, si] -= 1
        self.pm.susceptible[pi, ti] += 1

        self.immuneSourcePopRate[pi, si] = self.pm.susceptible[pi, si]*self.suscepCumulTransition[si]
        self.immuneSourcePopRate[pi, ti] = self.pm.susceptible[pi, ti]*self.suscepCumulTransition[ti]
        self.immunePopRate[pi] = 0.0
        for sn in range(self.susNum):
            self.immunePopRate[pi] += self.immuneSourcePopRate[pi, sn]

        self.iCounter += 1
        self.events.AddEvent(self.currentTime, SUSCCHANGE, si, pi, ti, 0)

        self.UpdateRates(pi)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void Birth(self, Py_ssize_t pi, Py_ssize_t hi):
        cdef double ws = 0.0
        for sn in range(self.susNum):
            ws += self.susceptHapPopRate[pi, hi, sn]
        si, self.rn = fastChoose1(self.susceptHapPopRate[pi, hi], ws, self.rn)

        self.pm.NewInfection(pi, si, hi)

        self.immuneSourcePopRate[pi, si] = self.pm.susceptible[pi, si]*self.suscepCumulTransition[si]
        self.immunePopRate[pi] = 0.0
        for sn in range(self.susNum):
            self.immunePopRate[pi] += self.immuneSourcePopRate[pi, sn]

        self.bCounter += 1
        self.events.AddEvent(self.currentTime, BIRTH, hi, pi, si, 0)

        self.UpdateRates(pi)
        self.MigrationRates()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void Death(self, Py_ssize_t pi, Py_ssize_t hi, bint add_event = True):
        self.pm.NewRecovery(pi, self.suscType[hi], hi)

        self.immuneSourcePopRate[pi, self.suscType[hi]] = self.pm.susceptible[pi, self.suscType[hi]]*self.suscepCumulTransition[self.suscType[hi]]
        self.immunePopRate[pi] = 0.0
        for sn in range(self.susNum):
            self.immunePopRate[pi] += self.immuneSourcePopRate[pi, sn]

        if add_event:
            self.dCounter += 1
            self.events.AddEvent(self.currentTime, DEATH, hi, pi, self.suscType[hi], 0)

        self.UpdateRates(pi)
        self.MigrationRates()

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

        mi, self.rn = fastChoose1( self.mRate[hi], self.tmRate[hi], self.rn)
        digit4 = 4**mi
        AS = int(floor(hi/digit4) % 4)
        DS, self.rn = fastChoose1(self.hapMutType[hi, mi], self.totalHapMutType[hi, mi], self.rn)
        if DS >= AS:
            DS += 1
        nhi = hi + (DS-AS)*digit4

        self.pm.liveBranches[pi, nhi] += 1
        self.pm.liveBranches[pi, hi] -= 1

        self.hapPopRate[pi, hi] = self.tEventHapPopRate[pi, hi]*self.pm.liveBranches[pi, hi]
        self.hapPopRate[pi, nhi] = self.tEventHapPopRate[pi, nhi]*self.pm.liveBranches[pi, nhi]
        self.infectPopRate[pi] = 0
        for hn in range(self.hapNum):
            self.infectPopRate[pi] += self.hapPopRate[pi, hn]
        self.popRate[pi] = self.infectPopRate[pi] + self.immunePopRate[pi]
        self.totalRate = 0
        for pn in range(self.popNum):
            self.totalRate += self.popRate[pn]

        self.mCounter += 1
        self.events.AddEvent(self.currentTime, MUTATION, hi, pi, nhi, 0)

        self.MigrationRates()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef Py_ssize_t GenerateMigration(self):
        cdef:
            Py_ssize_t tpi, spi, hi, si
        tpi, self.rn = fastChoose1( self.migPopRate, self.totalMigrationRate, self.rn)
        spi, self.rn = fastChoose2_skip( self.pm.totalInfectious, self.pm.globalInfectious-self.pm.totalInfectious[tpi], self.rn, skip = tpi)
        hi, self.rn = fastChoose2( self.pm.liveBranches[spi], self.pm.totalInfectious[spi], self.rn)
        si, self.rn = fastChoose2( self.pm.susceptible[tpi], self.pm.totalSusceptible[tpi], self.rn)
        if self.rn < self.pm.effectiveMigration[spi, tpi]*self.bRate[hi]*self.susceptibility[hi, si]/self.pm.maxEffectiveMigration[tpi]/self.maxEffectiveBirth:
            self.pm.NewInfection(tpi, si, hi)

            self.immuneSourcePopRate[tpi, si] = self.pm.susceptible[tpi, si]*self.suscepCumulTransition[si]
            self.immunePopRate[tpi] = 0.0
            for sn in range(self.susNum):
                self.immunePopRate[tpi] += self.immuneSourcePopRate[tpi, sn]

            self.pm.migPlus += 1
            self.pm.migCounter += 1
            self.events.AddEvent(self.currentTime, MIGRATION, hi, spi, si, tpi)

            self.UpdateRates(tpi)
            self.MigrationRates()
        else:
            self.pm.migNonPlus += 1
        return tpi

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void UpdateRates(self, Py_ssize_t pi):
        cdef double tmp

        for hn in range(self.hapNum):
            self.eventHapPopRate[pi, hn, 0] = self.BirthRate(pi, hn)
            tmp = (self.eventHapPopRate[pi, hn, 0] +
                   self.eventHapPopRate[pi, hn, 1] +
                   self.eventHapPopRate[pi, hn, 2] +
                   self.eventHapPopRate[pi, hn, 3] )
            self.tEventHapPopRate[pi, hn] = tmp
            self.hapPopRate[pi, hn] = self.tEventHapPopRate[pi, hn]*self.pm.liveBranches[pi, hn]

        self.infectPopRate[pi] = 0.0
        for hn in range(self.hapNum):
            self.infectPopRate[pi] += self.hapPopRate[pi, hn]
        self.popRate[pi] = self.infectPopRate[pi] + self.immunePopRate[pi]
        self.totalRate = 0.0
        for pn in range(self.popNum):
            self.totalRate += self.popRate[pn]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void CheckLockdown(self, Py_ssize_t pi):
        if self.pm.totalInfectious[pi] > self.pm.startLD[pi]*self.pm.sizes[pi] and self.pm.lockdownON[pi] == 0:
            self.UpdateContactDensity(pi, self.pm.contactDensityAfterLockdown[pi])
            self.pm.swapLockdown += 1
            self.pm.lockdownON[pi] = 1
        if self.pm.totalInfectious[pi] < self.pm.endLD[pi]*self.pm.sizes[pi] and self.pm.lockdownON[pi] == 1:
            self.UpdateContactDensity(pi, self.pm.contactDensityBeforeLockdown[pi])
            self.pm.swapLockdown += 1
            self.pm.lockdownON[pi] = 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void UpdateContactDensity(self, Py_ssize_t pi, float newCD):
        self.pm.contactDensity[pi] = newCD
        self.pm.SetEffectiveMigration()
        self.UpdateRates(pi)
        self.SetMaxBirth()
        self.MigrationRates()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void SetMaxBirth(self):
        for hn in range(self.hapNum):
            self.maxSusceptibility[hn] = 0.0
            for sn in range(self.susNum):
                if self.susceptibility[hn, sn] > self.maxSusceptibility[hn]:
                    self.maxSusceptibility[hn] = self.susceptibility[hn, sn]
        self.maxEffectiveBirth = 0.0
        for hn in range(self.hapNum):
            if self.maxEffectiveBirth < self.bRate[hn]*self.maxSusceptibility[hn]:
                self.maxEffectiveBirth = self.bRate[hn]*self.maxSusceptibility[hn]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef inline double MigrationRates(self):
        self.totalMigrationRate = 0.0
        for pn in range(self.popNum):
            self.migPopRate[pn] = self.pm.maxEffectiveMigration[pn]*self.maxEffectiveBirth*self.pm.totalSusceptible[pn]*(self.pm.globalInfectious-self.pm.totalInfectious[pn])
            self.totalMigrationRate += self.migPopRate[pn]

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

    def set_infectious_rate(self, rate, haplotype):
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

        self.SetRates()

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

        self.SetRates()

    def set_sampling_rate(self, rate, haplotype):
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

        self.SetRates()

    #TODO
    def set_mutation_rate(self, rate, probabilities, haplotype, mutation):
        if isinstance(rate, (int, float)) and isinstance(probabilities, list) and isinstance(haplotype, str) and isinstance(mutation,int):#DONE
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
            if mutation<0 or mutation>=self.sites:
                self.Error("There are no such mutation!")

            for haplotype in haplotypes:
                probabilities_allele = probabilities
                del probabilities_allele[self.calculate_allele(haplotype, mutation)]
                if sum(probabilities_allele) == 0:
                    self.Error("#TODO")
                self.mRate[haplotype, mutation] = rate
                self.hapMutType[haplotype, mutation, 0] = probabilities_allele[0]
                self.hapMutType[haplotype, mutation, 1] = probabilities_allele[1]
                self.hapMutType[haplotype, mutation, 2] = probabilities_allele[2]
                self.totalHapMutType[haplotype, mutation] = sum(probabilities_allele)
        elif rate==None and isinstance(probabilities, list) and isinstance(haplotype, str) and isinstance(mutation,int):#DONE
            if len(probabilities)!=4:
                self.Error("#TODO")
            for i in range(4):
                if isinstance(probabilities[i], (int, float)) == False:
                    self.Error("#TODO")
                if probabilities[i]<0:
                    self.Error("#TODO")
            haplotypes = self.create_list_haplotypes(haplotype)
            if mutation<0 or mutation>=self.sites:
                self.Error("There are no such mutation!")

            for haplotype in haplotypes:
                probabilities_allele = probabilities
                del probabilities_allele[self.calculate_allele(haplotype, mutation)]
                if sum(probabilities_allele) == 0:
                    self.Error("#TODO")
                self.hapMutType[haplotype, mutation, 0] = probabilities_allele[0]
                self.hapMutType[haplotype, mutation, 1] = probabilities_allele[1]
                self.hapMutType[haplotype, mutation, 2] = probabilities_allele[2]
                self.totalHapMutType[haplotype, mutation] = sum(probabilities_allele)
        elif isinstance(rate, (int, float)) and probabilities==None and isinstance(haplotype, str) and isinstance(mutation,int):#DONE
            if rate<0:
                self.Error("#TODO")
            haplotypes = self.create_list_haplotypes(haplotype)
            if mutation<0 or mutation>=self.sites:
                self.Error("There are no such mutation!")

            for haplotype in haplotypes:
                self.mRate[haplotype, mutation] = rate
        elif isinstance(rate, (int, float)) and isinstance(probabilities, list) and isinstance(haplotype, int) and isinstance(mutation,int):#DONE
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
            if mutation<0 or mutation>=self.sites:
                self.Error("There are no such mutation!")

            self.mRate[haplotype, mutation] = rate
            self.hapMutType[haplotype, mutation, 0] = probabilities[0]
            self.hapMutType[haplotype, mutation, 1] = probabilities[1]
            self.hapMutType[haplotype, mutation, 2] = probabilities[2]
            self.totalHapMutType[haplotype, mutation] = sum(probabilities)
        elif rate==None and isinstance(probabilities, list) and isinstance(haplotype, int) and isinstance(mutation,int):#DONE
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
            if mutation<0 or mutation>=self.sites:
                self.Error("There are no such mutation!")

            self.hapMutType[haplotype, mutation, 0] = probabilities[0]
            self.hapMutType[haplotype, mutation, 1] = probabilities[1]
            self.hapMutType[haplotype, mutation, 2] = probabilities[2]
            self.totalHapMutType[haplotype, mutation] = sum(probabilities)
        elif isinstance(rate, (int, float)) and probabilities==None and isinstance(haplotype, int) and isinstance(mutation,int):#DONE
            if rate<0:
                self.Error("#TODO")
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
                probabilities_allele = probabilities
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
                probabilities_allele = probabilities
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
                    probabilities_allele = probabilities
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
                    probabilities_allele = probabilities
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
            del probabilities[self.calculate_allele(haplotype, mutation)]
            if sum(probabilities) == 0:
                self.Error("#TODO")

            for hn in range(self.hapNum):
                for s in range(self.sites):
                    probabilities_allele = probabilities
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
                    probabilities_allele = probabilities
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

        self.SetRates()

    def set_contact_density(self, value, population):
        if isinstance(value, (int, float)) == False:
            self.Error("Incorrect type of contact density. Value should be int or float.")
        if value<0:
            self.Error("Incorrect value of contact density. Value should be more or equal 0.")

        if isinstance(population, int):
            if population<0 or population>=self.popNum:
                self.Error("There are no such population!")

            self.pm.contactDensity[population] = value
            self.pm.contactDensityBeforeLockdown[population] = value
        elif population == None:
            for pn in range(self.popNum):
                self.pm.contactDensity[pn] = value
                self.pm.contactDensityBeforeLockdown[pn] = value
        else:
            self.Error("Incorrect value of population. Value should be int or None.")

        self.SetRates()

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

            self.pm.contactDensityAfterLockdown[population] = parameters[0]
            self.pm.startLD[population] = parameters[1]
            self.pm.endLD[population] = parameters[2]
        elif population == None:
            for pn in range(self.popNum):
                self.pm.contactDensityAfterLockdown[pn] = parameters[0]
                self.pm.startLD[pn] = parameters[1]
                self.pm.endLD[pn] = parameters[2]
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

            self.pm.samplingMultiplier[population] = multiplier
        elif population == None:
            for pn in range(self.popNum):
                self.pm.samplingMultiplier[pn] = multiplier
        else:
            self.Error("Incorrect value of population. Value should be int or None.")

        self.SetRates()

    def set_migration_probability(self, rate, source, target):
        if isinstance(rate, (int, float)) == False:
            self.Error("Incorrect type of rate. Value should be int or float!")
        if rate<0 and rate>1:
            self.Error("Incorrect type of rate. Value should be between 0 and 1!")

        if isinstance(source, int) and isinstance(target, int):
            if source<0 or source>=self.popNum:
                self.Error("There are no such population!")
            if target<0 or target>=self.popNum:
                self.Error("There are no such population!")
            if source==target:
                self.Error("Source and target population type shouldn't be equal!")

            self.pm.migrationRates[source, target] = rate
            summa = 0
            for pn in range(self.popNum):
                summa += self.pm.migrationRates[source, pn]
            if summa > 1:
                self.Error("#TODO")
        elif source==None and isinstance(target, int):
            if target<0 or target>=self.popNum:
                self.Error("There are no such population!")

            for pn1 in range(self.popNum):
                if pn1 != target:
                    self.pm.migrationRates[pn1, target] = rate
                summa = 0
                for pn2 in range(self.popNum):
                    summa += self.pm.migrationRates[pn1, pn2]
                if summa > 1:
                    self.Error("#TODO")
        elif isinstance(source, int) and target==None:
            if source<0 or source>=self.popNum:
                self.Error("There are no such population!")

            for pn2 in range(self.popNum):
                if source != pn2:
                    self.pm.migrationRates[source, pn2] = rate
            summa = 0
            for pn in range(self.popNum):
                summa += self.pm.migrationRates[source, pn]
            if summa > 1:
                self.Error("#TODO")
        elif source==None and target==None:
            for pn1 in range(self.popNum):
                for pn2 in range(self.popNum):
                    if pn1 != pn2:
                        self.pm.migrationRates[pn1, pn2] = rate
                summa = 0
                for pn2 in range(self.popNum):
                    summa += self.pm.migrationRates[pn1, pn2]
                if summa > 1:
                    self.Error("#TODO")
        else:
            self.Error("Incorrect value of population. Value should be int or None.")

        self.pm.SetEffectiveMigration()
        self.MigrationRates()

    def set_susceptible(self, amount, source_type, target_type, population):
        if isinstance(amount, int) == False:
            self.Error("Incorrect value of amount. Value should be int.")
        if isinstance(source_type, int) == False:
            self.Error("Incorrect value of source susceptibility type. Value should be int.")
        if isinstance(target_type, int) == False:
            self.Error("Incorrect value of target susceptibility type. Value should be int.")
        if source_type<0 or source_type>=self.susNum:
            self.Error("There are no such susceptibility type!")
        if target_type<0 or target_type>=self.susNum:
            self.Error("There are no such susceptibility type!")
        if source_type==target_type:
            self.Error("Source and target susceptibility type shouldn't be equal!")

        if isinstance(population, int):
            if population<0 or population>=self.popNum:
                self.Error("There are no such population!")
            if amount<0 or amount>self.pm.susceptible[population, source_type]:
                self.Error("#TODO")

            self.pm.susceptible[population, source_type] -= amount
            self.pm.susceptible[population, target_type] += amount
        elif population==None:
            for pn in range(self.popNum):
                if amount<0 or amount>self.pm.susceptible[pn, source_type]:
                    self.Error("#TODO")
                self.pm.susceptible[pn, source_type] -= amount
                self.pm.susceptible[pn, target_type] += amount
        else:
            self.Error("Incorrect value of population. Value should be int or None.")

        self.SetRates()
       
    def set_immunity_type(self, susceptibility_type, haplotype):
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

        self.SetRates()

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

        self.SetRates()

    def print_basic_parameters(self):
        print("Basic rates")
        table = PrettyTable()

        field = ["H", "IR", "RR", "SR", "ST"]
        for s in range(self.sites):
            field.append("M" + str(s))
            field.append("MP" + str(s))
        table.field_names = field
        for hn in range(self.hapNum):
            list = [self.calculate_string(hn), self.bRate[hn], self.dRate[hn], self.sRate[hn], self.suscType[hn]]
            for s in range(self.sites):
                list.append(self.mRate[hn, s])
                list.append([self.hapMutType[hn, s, 0], self.hapMutType[hn, s, 1], self.hapMutType[hn, s, 2]])
            table.add_row(list)

        print(table)
        print("Legend:")
        print("H - haplotype")
        print("IR - infection rate")
        print("RR - recovery rate")
        print("SR - sampling rate")
        print("ST - susceptibility type")
        for s in range(self.sites):
            print("M" + str(s) + " - " + str(s) + " mutation rate")
            print("MP" + str(s) + " - " + str(s) + " mutation probabilities")
        print()

    def print_populations(self):
        print("Populations")
        table_populations = PrettyTable()

        table_populations.field_names = ["ID", "Size", "CD", "CDALD", "SLD", "ELD", "SM"]
        for pn in range(self.popNum):
            table_populations.add_row([pn, self.pm.sizes[pn], self.pm.contactDensity[pn], self.pm.contactDensityAfterLockdown[pn], self.pm.startLD[pn], self.pm.endLD[pn], self.pm.samplingMultiplier[pn]])
        
        print(table_populations)
        print("Legend:")
        print("ID - number of population")
        print("Size - size of population")
        print("CD - contact density")
        print("CDALD - contact density at lockdown")
        print("SLD - start of lockdown")
        print("ELD - end of lockdown")
        print("SM - sampling multiplier")
        print()

        print("Susceptible")
        table_susceptible = PrettyTable()

        field = ["ST\\ID"]
        for pn in range(self.popNum):
            field.append(pn)
        for sn in range(self.susNum):
            row = [sn]
            for pn in range(self.popNum):
                row.append(self.pm.susceptible[pn, sn])
            table_susceptible.add_row(row)
        table_susceptible.field_names = field

        print(table_susceptible)
        print("Legend:")
        print("ID - ID population")
        print("ST - susceptibility type")
        print()

        print("Migration matrix")
        table_migration = PrettyTable()

        field = ["S\\T"]
        for pn1 in range(self.popNum):
            field.append(pn1)
            row = [pn1]
            for pn2 in range(self.popNum):
                row.append(self.pm.migrationRates[pn1, pn2])
            table_migration.add_row(row)
        table_migration.field_names = field

        print(table_migration)
        print("Legend:")
        print("S - ID source population")
        print("T - ID target population")
        print()

    def print_immunity_model(self):
        print("Immunity model")
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

        print("Immunity transition rates")
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
        allele = 0
        for s in range(site):
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

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef GetGenealogy(self, seed):
        cdef:
            Py_ssize_t ptrTreeAndTime, n1, n2, id1, id2, id3, lbs, lbs_e, ns, nt, idt, ids, lbss
            double p
            vector[vector[vector[Py_ssize_t]]] liveBranchesS
            vector[vector[Py_ssize_t]] vecint2
            vector[Py_ssize_t] vecint1

            double e_time
            Py_ssize_t e_type_, e_population, e_haplotype, e_newHaplotype, e_newPopulation

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

        for e_id in range(self.events.ptr-1, -1, -1):
            # this event
            e_time = self.events.times[e_id]
            e_type_ = self.events.types[e_id]
            e_haplotype = self.events.haplotypes[e_id]
            e_population = self.events.populations[e_id]
            e_newHaplotype = self.events.newHaplotypes[e_id]
            e_newPopulation = self.events.newPopulations[e_id]
            if e_id == 7169:

                print("new event: ", e_id)
                print(e_type_)
                print(e_haplotype)
                print(e_population)
                print(e_newHaplotype)
                print(e_newPopulation)
                for pn in range(self.popNum):
                    for hn in range(self.hapNum):
                        print(self.pm.liveBranches[pn, hn], end=" ")
                    print()

                print()
                for pn in range(self.popNum):
                    for hn in range(self.hapNum):
                        print(len(liveBranchesS[pn][hn]), end=" ")
                    print()


            if e_type_ == BIRTH:
                lbs = liveBranchesS[e_population][e_haplotype].size()
                lbs_e = self.pm.liveBranches[e_population, e_haplotype]
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
                self.pm.liveBranches[e_population, e_haplotype] -= 1
            elif e_type_ == DEATH:
                self.pm.liveBranches[e_population, e_haplotype] += 1
            elif e_type_ == SAMPLING:
                self.pm.liveBranches[e_population, e_haplotype] += 1
                liveBranchesS[e_population][e_haplotype].push_back( ptrTreeAndTime )
                self.tree[ptrTreeAndTime] = -1
                self.times[ptrTreeAndTime] = e_time
                ptrTreeAndTime += 1
            elif e_type_ == MUTATION:
                lbs = liveBranchesS[e_population][e_newHaplotype].size()
                p = float(lbs)/self.pm.liveBranches[e_population, e_newHaplotype]
                if self.rndm.uniform() < p:
                    n1 = int(floor( lbs*self.rndm.uniform() ))
                    id1 = liveBranchesS[e_population][e_newHaplotype][n1]
                    liveBranchesS[e_population][e_newHaplotype][n1] = liveBranchesS[e_population][e_newHaplotype][lbs-1]
                    liveBranchesS[e_population][e_newHaplotype].pop_back()
                    liveBranchesS[e_population][e_haplotype].push_back(id1)
                    self.mut.AddMutation(id1, e_haplotype, e_newHaplotype, e_time)
                self.pm.liveBranches[e_population, e_newHaplotype] -= 1
                self.pm.liveBranches[e_population, e_haplotype] += 1
            elif e_type_ == SUSCCHANGE:
                pass
            elif e_type_ == MIGRATION:
                lbs = liveBranchesS[e_newPopulation][e_haplotype].size()
                p = float(lbs)/self.pm.liveBranches[e_newPopulation][e_haplotype]
                if self.rndm.uniform() < p:
                    nt = int(floor( lbs*self.rndm.uniform() ))
                    lbss = liveBranchesS[e_population][e_haplotype].size()
                    p1 = float(lbss)/self.pm.liveBranches[e_population, e_haplotype]
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
                        self.pm.mig.AddMigration(idt, e_time, e_population, e_newPopulation)
                    else:
                        liveBranchesS[e_population][e_haplotype].push_back(liveBranchesS[e_newPopulation][e_haplotype][nt])
                        liveBranchesS[e_newPopulation][e_haplotype][nt] = liveBranchesS[e_newPopulation][e_haplotype][lbs-1]
                        liveBranchesS[e_newPopulation][e_haplotype].pop_back()
                self.pm.liveBranches[e_newPopulation, e_haplotype] -= 1
            else:
                print("Unknown event type: ", e_type_)
                print("_________________________________")
                sys.exit(0)
        self.CheckTree()

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

    def LogDynamics(self, step_num, output_file):
        time_points = [i*self.currentTime/step_num for i in range(step_num+1)]
        suscepDate = np.zeros((self.popNum, self.susceptible_num), dtype=np.int64)
        hapDate = np.zeros((self.popNum, self.hapNum), dtype=np.int64)
        for i in range(self.popNum):
            for sn in range(self.susceptible_num):
                if sn == 0:
                    suscepDate[i, sn] = self.pm.sizes[i]
                else:
                    suscepDate[i, sn] = 0
            for hn in range(self.hapNum):
                hapDate[i, hn] = 0
        hapDate[0][0] += 1
        suscepDate[0][0] -= 1
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

    def get_data_infectious(self, pop, hap, step_num):
        time_points = [i*self.currentTime/step_num for i in range(step_num+1)]
        Date = np.zeros(step_num+1)
        Sample = np.zeros(step_num+1)

        point = 0
        for j in range(self.events.ptr):
            if time_points[point] < self.events.times[j]:
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
        return Date, Sample, time_points

    def get_data_susceptible(self, pop, sus, step_num):
        time_points = [i*self.currentTime/step_num for i in range(step_num+1)]
        Date = np.zeros(step_num+1)
        if sus == 0:
            Date[0] = self.pm.sizes[pop]

        point = 0
        for j in range(self.events.ptr):
            if time_points[point] < self.events.times[j]:
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

        return Date, time_points

    def sampleDate(self):
        time, pop, hap = [], [], []
        for i in range(self.events.ptr):
            if self.events.types[i] == SAMPLING:
                time.append(self.events.times[i])
                pop.append(self.events.populations[i])
                hap.append(self.events.haplotypes[i])
        return time, pop, hap

    def Stats(self):
        print("Number of samples:", self.sCounter)
        print("Total number of iterations: ", self.events.ptr-1)
        print("Size events: ", self.events.size)
        print("Current time: ", self.currentTime)

    def Debug(self):
        print("Parameters")
        print("swapLockdown: ", self.pm.swapLockdown)
        print("Migration plus: ", self.pm.migPlus)
        print("Migration non plus: ", self.pm.migNonPlus)
        print("Current time(mutable): ", self.currentTime)
        print("Random number(mutable): ", self.rn)
        print("Total rate(mutable): ", self.totalRate)
        print("Max effective birth(const): ", self.maxEffectiveBirth)
        print("Total migration rate(mutable): ", self.totalMigrationRate)
        print("Birth counter(mutable): ", self.bCounter)
        print("Death counter(mutable): ", self.dCounter)
        print("Sampling counter(mutable): ", self.sCounter)
        print("Migration counter(mutable): ", self.pm.migCounter)
        print("Mutation counter(mutable): ", self.mCounter)
        print("Immunity transition counter(mutable):", self.iCounter)
        print("Populations number(const): ", self.popNum)
        print("Mutations number(const): ", self.sites)
        print("Haplotypes number(const): ", self.hapNum)
        print("Susceptible number(const): ", self.susNum)
        print("Population model - globalInfectious(mutable): ", self.pm.globalInfectious)
        print("Susceptible type(): ", sep=" ", end="")
        for i in range(self.suscType.shape[0]):
            print(self.suscType[i], end=" ")
        print()
        print("Birth rate(const): ", sep="", end="")
        for i in range(self.hapNum):
            print(self.bRate[i], end=" ")
        print()
        print("Death rate(const): ", sep="", end="")
        for i in range(self.hapNum):
            print(self.dRate[i], end=" ")
        print()
        print("Sampling rate(const): ", sep="", end="")
        for i in range(self.hapNum):
            print(self.sRate[i], end=" ")
        print()
        print("Total mutation rate(const): ", sep="", end="")
        for i in range(self.hapNum):
            print(self.tmRate[i], end=" ")
        print()
        print("Migration population rate(mutable): ", sep="", end="")
        for i in range(self.popNum):
            print(self.migPopRate[i], end=" ")
        print()
        print("Population rate(mutable): ", sep="", end="")
        for i in range(self.popNum):
            print(self.popRate[i], end=" ")
        print()
        print("Population model - sizes(const): ", end="")
        for i in range(self.pm.sizes.shape[0]):
            print(self.pm.sizes[i], end=" ")
        print()
        print("Population model - totalSusceptible(mutable): ", end="")
        for i in range(self.pm.totalSusceptible.shape[0]):
            print(self.pm.totalSusceptible[i], end=" ")
        print()
        print("Population model - totalInfectious(mutable): ", end="")
        for i in range(self.pm.totalInfectious.shape[0]):
            print(self.pm.totalInfectious[i], end=" ")
        print()
        print("Population model - contac density(const): ", end=" ")
        for i in range(self.pm.sizes.shape[0]):
            print(self.pm.contactDensity[i], end=" ")
        print()
        print("Population model - max effective migration(const): ", end=" ")
        for i in range(self.pm.maxEffectiveMigration.shape[0]):
            print(self.pm.maxEffectiveMigration[i], end=" ")
        print()
        print("Population model - max susceptibility(const): ", end=" ")
        for i in range(self.maxSusceptibility.shape[0]):
            print(self.maxSusceptibility[i], end=" ")
        print()

        print("Population model - contactDensityAfterLockdown(const): ", end=" ")
        for i in range(self.pm.contactDensityAfterLockdown.shape[0]):
            print(self.pm.contactDensityAfterLockdown[i], end=" ")
        print()
        print("Population model - startLD(const): ", end=" ")
        for i in range(self.pm.startLD.shape[0]):
            print(self.pm.startLD[i], end=" ")
        print()
        print("Population model - endLD(const): ", end=" ")
        for i in range(self.pm.endLD.shape[0]):
            print(self.pm.endLD[i], end=" ")
        print()
        print("Population model - samplingMultiplier(const): ", end=" ")
        for i in range(self.pm.samplingMultiplier.shape[0]):
            print(self.pm.samplingMultiplier[i], end=" ")
        print()
        print("Population model - susceptible(mutable)----")
        for i in range(self.pm.sizes.shape[0]):
            for j in range(self.susNum):
                print(self.pm.susceptible[i, j], end=" ")
            print()
        print()
        print("Population model - migration rates(const)----")
        for i in range(self.pm.migrationRates.shape[0]):
            for j in range(self.pm.migrationRates.shape[1]):
                print(self.pm.migrationRates[i, j], end=" ")
            print()
        print()
        print("Population model - effective migration(const)----")
        for i in range(self.pm.effectiveMigration.shape[0]):
            for j in range(self.pm.effectiveMigration.shape[1]):
                print(self.pm.effectiveMigration[i, j], end=" ")
            print()
        print()
        print("Total event haplotype population rate(mutable)----")
        for i in range(self.popNum):
            for j in range(self.hapNum):
                print(self.tEventHapPopRate[i, j], end=" ")
            print()
        print()
        print("Haplotypes populations rates(mutable)----")
        for i in range(self.popNum):
            for j in range(self.hapNum):
                print(self.hapPopRate[i, j], end=" ")
            print()
        print()
        print("Mutation rate(const)----")
        for i in range(self.hapNum):
            for j in range(self.sites):
                print(self.mRate[i, j], end=" ")
            print()
        print()
        print("Susceptibility(const)----")
        for i in range(self.susceptibility.shape[0]):
            for j in range(self.susceptibility.shape[1]):
                print(self.susceptibility[i, j], end=" ")
            print()
        print()
        print("Event haplotypes populations rate(mutable)----")
        for i in range(self.popNum):
            for j in range(self.hapNum):
                for k in range(4):
                    print(self.eventHapPopRate[i, j, k], end=" ")
                print()
            print()
        print()
        print("Susceptible haplotypes populations rate(mutable)----")
        for i in range(self.popNum):
            for j in range(self.hapNum):
                for k in range(self.susNum):
                    print(self.susceptHapPopRate[i, j, k], end=" ")
                print()
            print()
        print()
        print("Probabilities of mutations(const)----")
        for i in range(self.hapMutType.shape[0]):
            for j in range(self.hapMutType.shape[1]):
                for k in range(self.hapMutType.shape[2]):
                    print(self.hapMutType[i, j, k], end=" ")
                print()
            print()
        print()

    def Output_tree_mutations(self):
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

    def writeMigrations(self, name_file):
        #with open(name_file + '.mig', 'w') as file:
        file = open(name_file + '.mig', 'w')
        file.write("Node Time Old_population New_population\n")
        for i in range(self.mig.nodeId.size()):
            file.write(str(self.mig.nodeId[i]) + " " + str(self.mig.time[i]) + " " + str(self.mig.oldPop[i]) + " " + str(self.mig.newPop[i]) + "\n")
        file.close()

    # def save_data(self):
    #     events, event_times, times, tree = [], [], [], []
    #     events, event_times = [], []
        
    #     for i in range(self.events.ptr):
    #         events.append(self.events.types[i])
    #         event_times.append(self.events.times[i])

    #     for i in range(2 * self.sCounter - 1):
    #         times.append(self.times[i])
    #         tree.append(self.tree[i])

    #     return events, event_times
    #     return events, event_times, times, tree
