# cython: language_level=3
# cython: initializedcheck = False
# distutils: language = c++

cimport cython

from libc.math cimport log, floor, abs
from libcpp.vector cimport vector
from mc_lib.rndm cimport RndmWrapper

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd

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



cdef class BirthDeathModel:
    cdef:
        RndmWrapper rndm

        Py_ssize_t bCounter, dCounter, sCounter, migCounter, mutCounter, popNum, dim, hapNum, susceptible_num, migPlus, migNonPlus, swapLockdown
        double currentTime, totalLen, rn, totalRate, maxEffectiveBirth, totalMigrationRate

        Events events
        PopulationModel pm
        Mutations mut
        Migrations mig

        long[::1] tree, suscType
        long[:,::1] liveBranches

        double[::1] bRate, dRate, sRate, tmRate, migPopRate, popRate, times, pm_maxEffectiveMigration, maxSusceptibility, elementsArr2, immunePopRate, infectPopRate, sourceSuscepTransition, suscepCumulTransition
        double[:,::1] pm_migrationRates, pm_effectiveMigration, tEventHapPopRate, hapPopRate, mRate, susceptibility, totalHapMutType, suscepTransition, immuneSourcePopRate
        double[:,:,::1] eventHapPopRate, susceptHapPopRate, hapMutType

    def __init__(self, bRate, dRate, sRate, mRate, populationModel=None, susceptible=None, suscepTransition=None, lockdownModel=None, samplingMultiplier=None, rndseed=1256, **kwargs):
        self.currentTime = 0.0
        self.sCounter = 0 #sample counter
        self.bCounter = 0
        self.dCounter = 0
        self.migCounter = 0
        self.mutCounter = 0
        self.events = Events()
        self.mut = Mutations()
        self.mig = Migrations()
        self.migPlus = 0
        self.migNonPlus = 0
        self.swapLockdown = 0

        if susceptible is None:
            self.susceptible_num = 2
        else:
            self.susceptible_num = len( susceptible[0][0] )

        #Set population model
        if populationModel is None:
            self.pm = PopulationModel( [ Population() ], self.susceptible_num)
            self.pm_migrationRates = np.asarray((0, 0), dtype=float)
        else:
            self.pm = PopulationModel( populationModel[0] , self.susceptible_num, lockdownModel, samplingMultiplier)
            self.pm_migrationRates = np.asarray(populationModel[1])
        self.popNum = self.pm.sizes.shape[0]
        self.pm_effectiveMigration = np.zeros((self.popNum, self.popNum), dtype=float)
        self.pm_maxEffectiveMigration = np.zeros(self.popNum, dtype=float)
        self.SetEffectiveMigration()

        #Initialise haplotypes
        if len(mRate) > 0:
            self.dim = len(mRate[0])
        else:
            self.dim = 0
        self.hapNum = int(4**self.dim)

        self.InitLiveBranches()

        self.elementsArr2 = np.zeros(2, dtype=float)

        if susceptible is None:
            self.susceptibility = np.asarray( [ [1.0, 0.0] for _ in range(self.hapNum) ] )
            self.suscType = np.ones(int(self.hapNum), dtype=np.int64)
        else:
            self.susceptibility = np.asarray( susceptible[0], dtype=float)
            self.suscType = np.asarray( susceptible[1], dtype=np.int64 )

        self.susceptHapPopRate = np.zeros((self.popNum, self.hapNum, self.susceptible_num), dtype=float)

        if suscepTransition is None:
            self.suscepTransition = np.zeros( (self.susceptible_num, self.susceptible_num), dtype=float)
        else:
            self.suscepTransition = np.asarray( suscepTransition )

        #Set rates
        self.SetRates(bRate, dRate, sRate, mRate)
        self.maxSusceptibility = np.zeros(self.hapNum, dtype=float)
        self.SetMaxBirth()
        self.migPopRate = np.zeros(len(self.pm_migrationRates), dtype=float)
        self.MigrationRates()

        #Set random generator
        self.rndm = RndmWrapper(seed=(rndseed, 0))

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void InitLiveBranches(self):
        self.liveBranches = np.zeros((self.popNum, self.hapNum), dtype=np.int64)
        self.events.AddEvent(0, 0, 0, 0, 0, 0)
        self.liveBranches[0, 0] += 2
        self.pm.NewInfection(0, 0)
        self.pm.NewInfection(0, 0)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void SetMaxBirth(self):
        for hn in range(self.hapNum):
            self.maxSusceptibility[hn] = 0.0
            for sType in range(self.susceptible_num):
                if self.susceptibility[hn, sType] > self.maxSusceptibility[hn]:
                    self.maxSusceptibility[hn] = self.susceptibility[hn, sType]
        self.maxEffectiveBirth = 0.0
        for hn in range(self.hapNum):
            if self.maxEffectiveBirth < self.bRate[hn]*self.maxSusceptibility[hn]:
                self.maxEffectiveBirth = self.bRate[hn]*self.maxSusceptibility[hn]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void SetEffectiveMigration(self):
        for pn in range(self.popNum):
            self.pm_maxEffectiveMigration[pn] = 0.0
        for pn1 in range(self.popNum):
            for pn2 in range(self.popNum):
                self.pm_effectiveMigration[pn1, pn2] = self.pm_migrationRates[pn1, pn2]*self.pm.contactDensity[pn2]/self.pm.sizes[pn2]+self.pm_migrationRates[pn2, pn1]*self.pm.contactDensity[pn1]/self.pm.sizes[pn1]
                if self.pm_effectiveMigration[pn1, pn2] > self.pm_maxEffectiveMigration[pn2]:
                    self.pm_maxEffectiveMigration[pn2] = self.pm_effectiveMigration[pn1, pn2]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void SetRates(self, bRate, dRate, sRate, mRate):
        self.bRate, self.dRate, self.sRate = np.asarray(bRate), np.asarray(dRate), np.asarray(sRate)

        self.mRate = np.zeros((len(mRate), len(mRate[0])), dtype=float)
        self.hapMutType = np.zeros((self.hapNum, len(mRate[0]), 3), dtype=float)
        self.totalHapMutType = np.zeros((self.hapNum, len(mRate[0])), dtype=float)
        for i in range(len(mRate)):
            for j in range(len(mRate[0])):
                self.totalHapMutType[i, j] = 0
                self.mRate[i, j] = mRate[i][j][0]
                self.hapMutType[i, j, 0] = mRate[i][j][1]
                self.totalHapMutType[i, j] += mRate[i][j][1]
                self.hapMutType[i, j, 1] = mRate[i][j][2]
                self.totalHapMutType[i, j] += mRate[i][j][2]
                self.hapMutType[i, j, 2] = mRate[i][j][3]
                self.totalHapMutType[i, j] += mRate[i][j][3]

        self.tmRate = np.zeros(len(mRate), dtype=float)
        for i in range(self.mRate.shape[0]):
            for j in range(self.mRate.shape[1]):
                self.tmRate[i] += self.mRate[i, j]

        self.eventHapPopRate = np.zeros((self.popNum, self.hapNum, 4), dtype=float)
        self.tEventHapPopRate = np.zeros((self.popNum, self.hapNum), dtype=float)
        for pn in range(self.popNum):
            for hn in range(self.hapNum):
                self.eventHapPopRate[pn, hn, 0] = self.BirthRate(pn, hn)
                self.eventHapPopRate[pn, hn, 1] = self.dRate[hn]
                self.eventHapPopRate[pn, hn, 2] = self.pm.samplingMultiplier[pn] * self.sRate[hn]
                self.eventHapPopRate[pn, hn, 3] = self.tmRate[hn]
                for i in range(4):
                    self.tEventHapPopRate[pn, hn] += self.eventHapPopRate[pn, hn, i]

        self.hapPopRate = np.zeros((self.popNum, self.hapNum), dtype=float)
        for pn in range(self.popNum):
            for hn in range(self.hapNum):
                self.HapPopRate(pn, hn)


        self.infectPopRate = np.zeros(self.popNum, dtype=float)
        for pn in range(self.popNum):
            for hn in range(self.hapNum):
                self.infectPopRate[pn] += self.hapPopRate[pn, hn]


        self.suscepCumulTransition = np.zeros(self.susceptible_num, dtype=float)
        for i in range(self.susceptible_num):
            for j in range(self.susceptible_num):
                self.suscepCumulTransition[i] += self.suscepTransition[i, j]

        self.immuneSourcePopRate = np.zeros((self.popNum, self.susceptible_num), dtype=float)
        self.immunePopRate = np.zeros(self.popNum, dtype=float)

        for pn in range(self.popNum):
            for j in range(self.susceptible_num):
                self.immuneSourcePopRate[pn, j] += self.suscepCumulTransition[j]*self.pm.susceptible[pn, j]
                self.immunePopRate[pn] += self.immuneSourcePopRate[pn, j]

        self.popRate = np.zeros(self.popNum, dtype=float)
        for pn in range(self.popNum):
            self.popRate[pn] = self.infectPopRate[pn] + self.immunePopRate[pn]

        self.totalRate = 0
        for pn in range(self.popNum):
            self.totalRate += self.popRate[pn]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef inline void HapPopRate(self, Py_ssize_t pi, Py_ssize_t hi):
        self.hapPopRate[pi, hi] = self.tEventHapPopRate[pi, hi]*self.liveBranches[pi, hi]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef inline double BirthRate(self, Py_ssize_t popId, Py_ssize_t haplotype):
        cdef double ws = 0.0
        for i in range(self.susceptible_num):
            self.susceptHapPopRate[popId, haplotype, i] = self.pm.susceptible[popId, i]*self.susceptibility[haplotype, i]
            ws += self.susceptHapPopRate[popId, haplotype, i]

        return self.bRate[haplotype]*ws/self.pm.sizes[popId]*self.pm.contactDensity[popId]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef inline double MigrationRates(self):
        self.totalMigrationRate = 0.0
        for p in range(self.popNum):
            self.migPopRate[p] = self.pm_maxEffectiveMigration[p]*self.maxEffectiveBirth*self.pm.totalSusceptible[p]*(self.pm.globalInfectious-self.pm.totalInfectious[p])
            self.totalMigrationRate += self.migPopRate[p]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef void SimulatePopulation(self, Py_ssize_t iterations, Py_ssize_t sampleSize, float time):
        cdef Py_ssize_t popId
        self.events.CreateEvents(iterations)
        #self.totalLen = 0.0
        print("Check!")
        while (self.events.ptr < self.events.size and self.sCounter < sampleSize and (time == -1 or self.currentTime < time)):
            self.SampleTime()
            popId = self.GenerateEvent()
            if self.totalRate == 0.0 or self.pm.globalInfectious == 0:
                break
            self.CheckLockdown(popId)
            
        print("Total number of iterations: ", self.events.ptr-1)
        print("Size events: ", self.events.size-1)
        if self.sCounter < 2: #TODO if number of sampled leaves is 0 (probably 1 as well), then GetGenealogy seems to go to an infinite cycle
            print("Less than two cases were sampled...")
            print("_________________________________")
            sys.exit(0)

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
            Py_ssize_t popId, haplotype, eventType, et

        self.rn = self.rndm.uniform()

        self.elementsArr2[0] = self.totalRate
        self.elementsArr2[1] = self.totalMigrationRate
        et, self.rn = fastChoose1( self.elementsArr2, self.totalRate+self.totalMigrationRate, self.rn)
        if et == 0:

        # if self.totalRate > self.rn * (self.totalRate + self.totalMigrationRate):
        #     self.rn = self.rn / self.totalRate

            popId, self.rn = fastChoose1( self.popRate, self.totalRate, self.rn)
            self.elementsArr2[0] = self.immunePopRate[popId]
            self.elementsArr2[1] = self.infectPopRate[popId]
            immune_vs_infect, self.rn = fastChoose1( self.elementsArr2, self.popRate[popId], self.rn)
            if immune_vs_infect == 0:
                self.ImmunityTransition(popId)
            else:
                haplotype, self.rn = fastChoose1( self.hapPopRate[popId], self.infectPopRate[popId], self.rn)
                eventType, self.rn = fastChoose1( self.eventHapPopRate[popId, haplotype], self.tEventHapPopRate[popId, haplotype], self.rn)
                if eventType == BIRTH:
                    self.Birth(popId, haplotype)
                elif eventType == DEATH:
                    self.Death(popId, haplotype)
                elif eventType == SAMPLING:
                    self.Sampling(popId, haplotype)
                else:
                    self.Mutation(popId, haplotype)
        else:
            popId = self.GenerateMigration()
        return popId

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void ImmunityTransition(self, Py_ssize_t popId):
        cdef:
            Py_ssize_t sourceImmune, targetImmune
        sourceImmune, self.rn = fastChoose1( self.immuneSourcePopRate[popId], self.immunePopRate[popId], self.rn)
        targetImmune, self.rn = fastChoose1( self.suscepTransition[sourceImmune], self.suscepCumulTransition[sourceImmune], self.rn)
        self.events.AddEvent(self.currentTime, SUSCCHANGE, sourceImmune, popId, targetImmune, 0)
        self.pm.susceptible[popId, sourceImmune] -= 1
        self.pm.susceptible[popId, targetImmune] += 1

        self.immuneSourcePopRate[popId, sourceImmune] = self.pm.susceptible[popId, sourceImmune]*self.suscepCumulTransition[sourceImmune]
        self.immuneSourcePopRate[popId, targetImmune] = self.pm.susceptible[popId, targetImmune]*self.suscepCumulTransition[targetImmune]
        self.immunePopRate[popId] = 0.0
        for j in range(self.susceptible_num):
            self.immunePopRate[popId] += self.immuneSourcePopRate[popId, j]
        self.UpdateRates(popId)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void Birth(self, Py_ssize_t popId, Py_ssize_t haplotype):
        self.liveBranches[popId, haplotype] += 1

        cdef double ws = 0.0
        for i in range(self.susceptible_num):
            ws += self.susceptHapPopRate[popId, haplotype, i]
        st, self.rn = fastChoose1(self.susceptHapPopRate[popId, haplotype], ws, self.rn)

        self.pm.NewInfection(popId, st)

        self.events.AddEvent(self.currentTime, BIRTH, haplotype, popId, 0, 0)
        self.immuneSourcePopRate[popId, st] = self.pm.susceptible[popId, st]*self.suscepCumulTransition[st]
        self.immunePopRate[popId] = 0.0
        for j in range(self.susceptible_num):
            self.immunePopRate[popId] += self.immuneSourcePopRate[popId, j]
        self.UpdateRates(popId)
        self.bCounter += 1
        self.MigrationRates()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void Death(self, Py_ssize_t popId, Py_ssize_t haplotype, bint add_event = True):
        self.liveBranches[popId, haplotype] -= 1
        self.pm.NewRecovery(popId, self.suscType[haplotype])

        self.immuneSourcePopRate[popId, self.suscType[haplotype]] = self.pm.susceptible[popId, self.suscType[haplotype]]*self.suscepCumulTransition[self.suscType[haplotype]]
        self.immunePopRate[popId] = 0.0
        for j in range(self.susceptible_num):
            self.immunePopRate[popId] += self.immuneSourcePopRate[popId, j]

        if add_event:
            self.dCounter += 1
            self.events.AddEvent(self.currentTime, DEATH, haplotype, popId, self.suscType[haplotype], 0)

        self.UpdateRates(popId)
        self.MigrationRates()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void Sampling(self, Py_ssize_t popId, Py_ssize_t haplotype):
        self.events.AddEvent(self.currentTime, SAMPLING, haplotype, popId, self.suscType[haplotype], 0)

        self.Death(popId, haplotype, False)
        self.sCounter += 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void Mutation(self, Py_ssize_t popId, Py_ssize_t haplotype):
        cdef:
            Py_ssize_t mutationType, digit4, AS, DS, newHaplotype

        mutationType, self.rn = fastChoose1( self.mRate[haplotype], self.tmRate[haplotype], self.rn)
        digit4 = 4**mutationType
        AS = int(floor(haplotype/digit4) % 4) #А тут хорошо всё?
        DS, self.rn = fastChoose1(self.hapMutType[haplotype, mutationType], self.totalHapMutType[haplotype, mutationType], self.rn)#TODO non-uniform rates???
        if DS >= AS:
            DS += 1
        newHaplotype = haplotype + (DS-AS)*digit4

        self.liveBranches[popId, newHaplotype] += 1
        self.liveBranches[popId, haplotype] -= 1

        self.events.AddEvent(self.currentTime, MUTATION, haplotype, popId, newHaplotype, 0)

        self.HapPopRate(popId, haplotype)
        self.HapPopRate(popId, newHaplotype)

        self.infectPopRate[popId] = 0
        for i in range(self.hapNum):
            self.infectPopRate[popId] += self.hapPopRate[popId, i]
        self.popRate[popId] = self.infectPopRate[popId] + self.immunePopRate[popId]

        self.totalRate = 0
        for i in range(self.popNum):
            self.totalRate += self.popRate[i]
        self.mutCounter += 1
        self.MigrationRates()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef Py_ssize_t GenerateMigration(self):
        cdef:
            Py_ssize_t targetPopId, sourcePopId, haplotype, suscType
            double p_accept
        targetPopId, self.rn = fastChoose1( self.migPopRate, self.totalMigrationRate, self.rn)
        sourcePopId, self.rn = fastChoose2_skip( self.pm.totalInfectious, self.pm.globalInfectious-self.pm.totalInfectious[targetPopId], self.rn, skip = targetPopId)
        haplotype, self.rn = fastChoose2( self.liveBranches[sourcePopId], self.pm.totalInfectious[sourcePopId], self.rn)
        suscType, self.rn = fastChoose2( self.pm.susceptible[targetPopId], self.pm.totalSusceptible[targetPopId], self.rn)
        p_accept = self.pm_effectiveMigration[sourcePopId, targetPopId]*self.bRate[haplotype]*self.susceptibility[haplotype, suscType]/self.pm_maxEffectiveMigration[targetPopId]/self.maxEffectiveBirth
        if self.rn < p_accept:
            self.liveBranches[targetPopId, haplotype] += 1
            self.pm.NewInfection(targetPopId, suscType)

            self.immuneSourcePopRate[targetPopId, suscType] = self.pm.susceptible[targetPopId, suscType]*self.suscepCumulTransition[suscType]
            self.immunePopRate[targetPopId] = 0.0
            for j in range(self.susceptible_num):
                self.immunePopRate[targetPopId] += self.immuneSourcePopRate[targetPopId, j]

            self.UpdateRates(targetPopId)
            self.MigrationRates()
            self.events.AddEvent(self.currentTime, MIGRATION, haplotype, sourcePopId, 0, targetPopId)
            self.migPlus += 1
            self.migCounter += 1
        else:
            self.migNonPlus += 1
        return targetPopId

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
            self.HapPopRate(pi, hn)

        self.infectPopRate[pi] = 0.0
        for hn in range(self.hapNum):
            self.infectPopRate[pi] += self.hapPopRate[pi, hn]

        self.popRate[pi] = self.infectPopRate[pi] + self.immunePopRate[pi]

        self.totalRate = 0.0
        for pn in range(self.popNum):
            self.totalRate += self.popRate[pn]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void UpdateContactDensity(self, Py_ssize_t pi, float newCD):
        self.pm.contactDensity[pi] = newCD
        self.SetEffectiveMigration()
        self.SetMaxBirth()
        self.MigrationRates()
        self.UpdateRates(pi)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void CheckLockdown(self, Py_ssize_t pi):
        if self.pm.totalInfectious[pi] > self.pm.startLD[pi] and self.pm.lockdownON[pi] == 0:
            self.UpdateContactDensity(pi, self.pm.contactDensityAfterLockdown[pi] )
            self.swapLockdown += 1
            self.pm.lockdownON[pi] = 1
        if self.pm.totalInfectious[pi] < self.pm.endLD[pi] and self.pm.lockdownON[pi] == 1:
            self.UpdateContactDensity(pi, self.pm.contactDensityBeforeLockdown[pi] )
            self.swapLockdown += 1
            self.pm.lockdownON[pi] = 0

    cpdef void UpdateMigration(self, float total_mig_rate):
        for pn1 in range(self.popNum):
            for pn2 in range(self.popNum):
                if pn1 != pn2:
                    self.pm_migrationRates[pn1, pn2] = total_mig_rate / (self.popNum-1)
        self.SetEffectiveMigration()
        self.SetMaxBirth()
        self.MigrationRates()
        for pn in range(self.popNum):
            self.UpdateRates(pn)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef GetGenealogy(self, rndseed):
        cdef:
            Py_ssize_t ptrTreeAndTime, n1, n2, id1, id2, id3, lbs, lbs_e, ns, nt, idt, ids, lbss
            double p
            vector[vector[vector[Py_ssize_t]]] liveBranchesS
            vector[vector[Py_ssize_t]] vecint2
            vector[Py_ssize_t] vecint1

            double e_time
            Py_ssize_t e_type_, e_population, e_haplotype, e_newHaplotype, e_newPopulation

        if rndseed != None:
            self.rndm = RndmWrapper(seed=(rndseed, 0))

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
            e_population = self.events.populations[e_id]
            e_haplotype = self.events.haplotypes[e_id]
            e_newHaplotype = self.events.newHaplotypes[e_id]
            e_newPopulation = self.events.newPopulations[e_id]
            if e_type_ == BIRTH:
                lbs = liveBranchesS[e_population][e_haplotype].size()
                lbs_e = self.liveBranches[e_population][e_haplotype]
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
                self.liveBranches[e_population][e_haplotype] -= 1
            elif e_type_ == DEATH:
                self.liveBranches[e_population][e_haplotype] += 1
            elif e_type_ == SAMPLING:
                self.liveBranches[e_population][e_haplotype] += 1
                liveBranchesS[e_population][e_haplotype].push_back( ptrTreeAndTime )
                self.tree[ptrTreeAndTime] = -1
                self.times[ptrTreeAndTime] = e_time
                ptrTreeAndTime += 1
            elif e_type_ == MUTATION:
                lbs = liveBranchesS[e_population][e_newHaplotype].size()
                p = float(lbs)/self.liveBranches[e_population][e_newHaplotype]
                if self.rndm.uniform() < p:
                    n1 = int(floor( lbs*self.rndm.uniform() ))
                    id1 = liveBranchesS[e_population][e_newHaplotype][n1]
                    liveBranchesS[e_population][e_newHaplotype][n1] = liveBranchesS[e_population][e_newHaplotype][lbs-1]
                    liveBranchesS[e_population][e_newHaplotype].pop_back()
                    liveBranchesS[e_population][e_haplotype].push_back(id1)
                    self.mut.AddMutation(id1, e_haplotype, e_newHaplotype, e_time)
                self.liveBranches[e_population][e_newHaplotype] -= 1
                self.liveBranches[e_population][e_haplotype] += 1
            elif e_type_ == SUSCCHANGE:
                pass
            elif e_type_ == MIGRATION:
                lbs = liveBranchesS[e_newPopulation][e_haplotype].size()
                p = float(lbs)/self.liveBranches[e_newPopulation][e_haplotype]
                if self.rndm.uniform() < p:
                    nt = int(floor( lbs*self.rndm.uniform() ))
                    lbss = liveBranchesS[e_population][e_haplotype].size()
                    p1 = float(lbss)/self.liveBranches[e_population][e_haplotype]
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
                self.liveBranches[e_newPopulation][e_haplotype] -= 1
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

        if output_file == True:
            for i in range(self.popNum-1, -1, -1):
                logDynamics[i].close()
        else: 
            return log

    def get_data(self, pop, hap, step_num):
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
                elif self.events.types[j] == MIGRATION:
                    Date[point] -= 1
                elif self.events.types[j] == MUTATION:
                    Date[point] -= 1
            elif self.events.types[j] == MUTATION and self.events.newHaplotypes[j] == hap and self.events.populations[j] == pop:
                Date[point] += 1
            elif self.events.types[j] == MIGRATION and self.events.newPopulations[j] == pop and self.events.haplotypes[j] == hap:
                Date[point] += 1
        return Date, Sample, time_points

    def sampleDate(self):
        time, pop, hap = [], [], []
        for i in range(self.events.ptr):
            if self.events.types[i] == SAMPLING:
                time.append(self.events.times[i])
                pop.append(self.events.populations[i])
                hap.append(self.events.haplotypes[i])
        return time, pop, hap

    def Report(self):
        print("Number of samples:", self.sCounter)

    def Debug(self):
        print("Parameters")
        print("swapLockdown: ", self.swapLockdown)
        print("Migration plus: ", self.migPlus)
        print("Migration non plus: ", self.migNonPlus)
        print("Current time(mutable): ", self.currentTime)
        print("Random number(mutable): ", self.rn)
        print("Total rate(mutable): ", self.totalRate)
        print("Max effective birth(const): ", self.maxEffectiveBirth)
        print("Total migration rate(mutable): ", self.totalMigrationRate)
        print("Birth counter(mutable): ", self.bCounter)
        print("Death counter(mutable): ", self.dCounter)
        print("Sampling counter(mutable): ", self.sCounter)
        print("Migration counter(mutable): ", self.migCounter)
        print("Mutation counter(mutable): ", self.mutCounter)
        print("Populations number(const): ", self.popNum)
        print("Mutations number(const): ", self.dim)
        print("Haplotypes number(const): ", self.hapNum)
        print("Susceptible number(const): ", self.susceptible_num)
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
        for i in range(self.pm_maxEffectiveMigration.shape[0]):
            print(self.pm_maxEffectiveMigration[i], end=" ")
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
            for j in range(self.susceptible_num):
                print(self.pm.susceptible[i, j], end=" ")
            print()
        print()
        print("Population model - migration rates(const)----")
        for i in range(self.pm_migrationRates.shape[0]):
            for j in range(self.pm_migrationRates.shape[1]):
                print(self.pm_migrationRates[i, j], end=" ")
            print()
        print()
        print("Population model - effective migration(const)----")
        for i in range(self.pm_effectiveMigration.shape[0]):
            for j in range(self.pm_effectiveMigration.shape[1]):
                print(self.pm_effectiveMigration[i, j], end=" ")
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
            for j in range(self.dim):
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
                for k in range(self.susceptible_num):
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
        file = open(name_file + '.mig', 'w')
        file.write("Node Time Old_population New_population\n")
        for i in range(self.mig.nodeId.size()):
            file.write(str(self.mig.nodeId[i]) + " " + str(self.mig.time[i]) + " " + str(self.mig.oldPop[i]) + " " + str(self.mig.newPop[i]) + "\n")
        file.close()