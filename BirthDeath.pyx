# cython: language_level=3
# distutils: language = c++

cimport cython

from libc.math cimport log
from libcpp.vector cimport vector
from mc_lib.rndm cimport RndmWrapper

import numpy as np
np.random.seed(1256)
import random
import sys
from math import floor

def print_err(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class Event:
    def __init__(self, time, type, population, haplotype, newHaplotype = None, newPopulation = None):
        self.time = time
        self.type = type
        self.population = population
        self.haplotype = haplotype
        self.newHaplotype = newHaplotype
        self.newPopulation = newPopulation

cdef class Mutation:
    cdef:
        int _nodeId, _AS, _DS
        double _time

    def __init__(self, int nodeId, double time, int AS, int DS):#AS = ancestral state, DS = derived state
        self._nodeId = nodeId
        self._time = time
        self._AS = AS
        self._DS = DS

    @property
    def nodeId(self):
        return self._nodeId

    @nodeId.setter
    def nodeId(self, nodeId):
        self._nodeId = nodeId

    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, time):
        self._time = time

    @property
    def AS(self):
        return self._AS

    @AS.setter
    def AS(self, AS):
        self._AS = AS

    @property
    def DS(self):
        return self._DS

    @DS.setter
    def DS(self, DS):
        self._DS = DS

cdef class Population:
    cdef:
        int _size, _susceptible, _infectious
        double _contactDensity

    def __init__(self, int size = 1000000, double contactDensity = 1.0):
        self._size = size
        self._susceptible = size
        self._infectious = 0
        self._contactDensity = contactDensity

    @property
    def size(self):
        return self._size

    @size.setter
    def size(self, size):
        self._size = size

    @property
    def susceptible(self):
        return self._susceptible

    @susceptible.setter
    def susceptible(self, susceptible):
        self._susceptible = susceptible

    @property
    def infectious(self):
        return self._infectious

    @infectious.setter
    def infectious(self, infectious):
        self._infectious = infectious

    @property
    def contactDensity(self):
        return self._contactDensity

    @contactDensity.setter
    def contactDensity(self, contactDensity):
        self._contactDensity = contactDensity

cdef class PopulationModel:
    cdef:
        list _populations, _migrationRates

    def __init__(self, list populations, list migrationRates):
        self._populations = populations #list/array with n populations список классов
        self._migrationRates = migrationRates #n * n-1 matrix with migration rates and zeroes on the diagonal список списков

    @property
    def populations(self):
        return self._populations

    @populations.setter
    def populations(self, populations):
        self._populations = populations

    @property
    def migrationRates(self):
        return self._migrationRates

    @migrationRates.setter
    def migrationRates(self, migrationRates):
        self._migrationRates = migrationRates

class NeutralMutations:
    def __init__(self):
        self.muRate = 0.0

    def muRate(self):
        return self.muRate

def fastChoose(w, tw):
    rn = tw*np.random.rand()
    i = 0
    total = w[0]
    while total < rn and i < len(w) - 1:
        i += 1
#        if i == len(w):
#            return -1
        total += w[i]
    if w[i] == 0.0:
        print_err("fastChoose() alert: 0-weight sampled")
    return i

cdef list fastChoose1(vector[double] w, double tw, double rn):
    rn = tw*rn
    i = 0
    total = w[0]
    while total < rn and i < w.size() - 1:
        i += 1
        total += w[i]
    if w[i] == 0.0:
        print_err("fastChoose() alert: 0-weight sampled")
    return [ i, ( rn-(total-w[i]) )/w[i] ]

def fastRandint(n):
    return floor( n*np.random.rand() )

cdef class BirthDeathModel:
    cdef:
        RndmWrapper rndm

        int sCounter, popNum, dim, hapNum
        double currentTime, totalRate, rn
        Py_ssize_t susceptible, lbCounter
        vector[double] popRate, tmRate, migPopRate, bRate, dRate, sRate, vecdob1
        vector[vector[double]] hapPopRate, birthHapPopRate, tEventHapPopRate, liveBranches, vecdob2
        vector[vector[vector[double]]] eventHapPopRate

        vector[int] tree, times, vecint1
        vector[vector[int]] vecint2
        vector[vector[vector[int]]] liveBranchesS


        object mRate
        object events, populationModel



    def __init__(self, bRate, dRate, sRate, mRate = [], **kwargs):
        self.currentTime = 0.0
        self.sCounter = 0 #sample counter
        self.events = []

        #Set population model
        if "populationModel" in kwargs:
            self.populationModel = kwargs["populationModel"]
        else:
            self.populationModel = PopulationModel([ Population() ], [[0.0]])
        self.susceptible = sum([pop.susceptible for pop in self.populationModel.populations])
        self.popNum = len( self.populationModel.populations )

        #Initialise haplotypes
        if len(mRate) > 0:
            self.dim = len(mRate[0])
        else:
            self.dim = 0
        self.hapNum = int(4**self.dim)
        self.InitLiveBranches()

        #Set rates
        self.SetRates(bRate, dRate, sRate, mRate)

        self.rndm = RndmWrapper(seed=(1256, 0))

    cdef void InitLiveBranches(self):
        for i in range( len( self.populationModel.populations) ):
            self.liveBranches.push_back(self.vecdob1)
            for j in range(self.hapNum):
                self.liveBranches[i].push_back(0)

        self.liveBranches[0][0] = 2
        self.lbCounter = 2 #live branch counter
        self.populationModel.populations[0].susceptible -= 2
        self.susceptible -= 2
        self.populationModel.populations[0].infectious = 2




    cdef void SetRates(self, object bRate, object dRate, object sRate, object mRate):
        self.bRate, self.dRate, self.sRate, self.mRate = bRate, dRate, sRate, mRate
        # for i in range(len(bRate)):
        #     self.bRate.push_back(bRate[i])
        #     self.dRate.push_back(dRate[i])
        #     self.sRate.push_back(sRate[i])
        # self.mRate = mRate

        #self.tmRate = [ sum(mr) for mr in self.mRate ]
        for i in range(len(self.mRate)):
            self.tmRate.push_back(0)
            for j in range(len(self.mRate[i])):
                self.tmRate[i] += self.mRate[i][j]

        #self.migPopRate = [sum(mr) for mr in self.populationModel.migrationRates]
        for i in range( len(self.populationModel.migrationRates) ):
            self.migPopRate.push_back(0)
            for j in range( len(self.populationModel.migrationRates[i]) ):
                self.migPopRate[i] += self.populationModel.migrationRates[i][j]

        #self.birthHapPopRate = [ [0.0]*self.hapNum for _ in range( self.popNum ) ]
        for i in range(self.popNum):
            self.birthHapPopRate.push_back(self.vecdob1)
            for _ in range(self.hapNum):
                self.birthHapPopRate[i].push_back(0)

        #self.eventHapPopRate = [ [ [0.0 for _ in range(5) ] for _ in range(self.hapNum) ] for _ in range(self.popNum) ]
        for i in range(self.popNum):
            self.eventHapPopRate.push_back(self.vecdob2)
            for j in range(self.hapNum):
                self.eventHapPopRate[i].push_back(self.vecdob1)
                for _ in range(5):
                    self.eventHapPopRate[i][j].push_back(0)

        #self.tEventHapPopRate = [ [ 0.0 for _ in range(self.hapNum) ] for _ in range(self.popNum) ]
        for i in range(self.popNum):
            self.tEventHapPopRate.push_back(self.vecdob1)
            for _ in range(self.hapNum):
                self.tEventHapPopRate[i].push_back(0)

        for pn in range(self.popNum):
            for hn in range(self.hapNum):
                self.birthHapPopRate[pn][hn] = self.BirthRate(pn, hn)
                self.eventHapPopRate[pn][hn] = [self.birthHapPopRate[pn][hn], self.dRate[hn], self.sRate[hn], self.migPopRate[pn], self.tmRate[hn] ]
                
                #self.tEventHapPopRate[pn][hn] = sum(self.eventHapPopRate[pn][hn])
                self.tEventHapPopRate[pn][hn] = 0
                for i in range(self.eventHapPopRate[pn][hn].size()):
                    self.tEventHapPopRate[pn][hn] += self.eventHapPopRate[pn][hn][i]

        
        #self.hapPopRate = [ [0.0 for hn in range(self.hapNum)] for pn in range(self.popNum) ]
        for i in range(self.popNum):
            self.hapPopRate.push_back(self.vecdob1)
            for j in range(self.hapNum):
                self.hapPopRate[i].push_back(0)

        for pn in range(self.popNum):
            for hn in range(self.hapNum):
                self.HapPopRate(pn, hn)

        #self.popRate = [ sum(self.hapPopRate[pn]) for pn in range(self.popNum) ]
        #self.totalRate = sum( self.popRate )
        self.totalRate = 0.0
        for i in range(self.popNum):
            self.popRate.push_back(0)
            for j in range(self.hapPopRate[i].size()):
                self.popRate[i] += self.hapPopRate[i][j]
                self.totalRate += self.hapPopRate[i][j]


    cdef void HapPopRate(self, Py_ssize_t popId, Py_ssize_t haplotype):
        self.hapPopRate[popId][haplotype] = self.tEventHapPopRate[popId][haplotype]*self.liveBranches[popId][haplotype]

    cdef double BirthRate(self, Py_ssize_t popId, Py_ssize_t haplotype):
        cdef object pop = self.populationModel.populations[popId]
        return self.bRate[haplotype]*pop.susceptible/pop.size*pop.contactDensity

    cdef void GenerateEvent(self):
        cdef:
            Py_ssize_t popId, haplotype, eventType

        self.rn = self.rndm.uniform()

        popId, self.rn = fastChoose1( self.popRate, self.totalRate, self.rn)
        haplotype, self.rn = fastChoose1( self.hapPopRate[popId], self.popRate[popId], self.rn)
        eventType, self.rn = fastChoose1( self.eventHapPopRate[popId][haplotype], self.tEventHapPopRate[popId][haplotype], self.rn)

        if eventType == 0:
            self.Birth(popId, haplotype)
        elif eventType == 1:
            self.Death(popId, haplotype)
        elif eventType == 2:
            self.Sampling(popId, haplotype)
        elif eventType == 3:
            self.Migration(popId, haplotype)
        else:
            self.Mutation(popId, haplotype)

    cdef void Migration(self, Py_ssize_t popId, Py_ssize_t haplotype):
        cdef Py_ssize_t targetPopId
        targetPopId, self.rn = fastChoose1(self.populationModel.migrationRates[popId], self.migPopRate[popId], self.rn)
        if targetPopId >= popId:
            targetPopId += 1

        self.liveBranches[targetPopId][haplotype] += 1
        self.liveBranches[popId][haplotype] -= 1

        self.populationModel.populations[popId].susceptible += 1
        self.populationModel.populations[popId].infectious -= 1

        self.populationModel.populations[targetPopId].susceptible -= 1
        self.populationModel.populations[targetPopId].infectious += 1

        event = Event(self.currentTime, 3, popId, haplotype, newPopulation = targetPopId)
        self.events.append(event)

        self.HapPopRate(popId, haplotype)
        self.HapPopRate(targetPopId, haplotype)

        #self.popRate[popId] = sum( self.hapPopRate[popId] )
        self.popRate[popId] = 0
        for i in range(self.hapPopRate[popId].size()):
            self.popRate[popId] += self.hapPopRate[popId][i]

        #self.popRate[targetPopId] = sum(self.hapPopRate[targetPopId])
        self.popRate[targetPopId] = 0
        for i in range(self.hapPopRate[targetPopId].size()):
            self.popRate[targetPopId] += self.hapPopRate[targetPopId][i]
        
        #self.totalRate = sum( self.popRate )
        self.totalRate = 0
        for i in range(self.popRate.size()):
            self.totalRate += self.popRate[i]



    cdef void Mutation(self, Py_ssize_t popId, Py_ssize_t haplotype):
        cdef:
            int mutationType, digit4, AS, DS, newHaplotype

        mutationType, self.rn = fastChoose1( self.mRate[haplotype], self.tmRate[haplotype], self.rn)
        digit4 = 4**mutationType
        AS = floor(haplotype/digit4)%4
        DS, self.rn = fastChoose1([1.0, 1.0, 1.0], 3.0, self.rn)#TODO non-uniform rates???
        if DS >= AS:
            DS += 1
        #self.mutations.append(Mutation(self.liveBranches[popId][haplotype][affectedBranch], self.currentTime, AS, DS))
        newHaplotype = haplotype + (DS-AS)*digit4

        self.liveBranches[popId][newHaplotype] += 1
        self.liveBranches[popId][haplotype] -= 1

        event = Event(self.currentTime, 4, popId, haplotype, newHaplotype = newHaplotype)
        self.events.append(event)

        self.HapPopRate(popId, haplotype)
        self.HapPopRate(popId, newHaplotype)
        
        #self.popRate[popId] = sum( self.hapPopRate[popId] )
        self.popRate[popId] = 0
        for i in range(self.hapPopRate[popId].size()):
            self.popRate[popId] += self.hapPopRate[popId][i]
        
        #self.totalRate = sum( self.popRate )
        self.totalRate = 0
        for i in range(self.popRate.size()):
            self.totalRate += self.popRate[i]

    #Сломалось время !! Починить
    cdef void SampleTime(self):
        cdef double tau = self.totalRate * log(1 + self.rndm.uniform())
        self.currentTime += tau

    cdef void Birth(self, Py_ssize_t popId, Py_ssize_t haplotype):
        self.liveBranches[popId][haplotype] += 1

        self.populationModel.populations[popId].susceptible -= 1
        self.populationModel.populations[popId].infectious += 1

        event = Event(self.currentTime, 0, popId, haplotype)
        self.events.append(event)

        for h in range(self.hapNum):
            self.birthHapPopRate[popId][h] = self.BirthRate(popId, h)
            self.eventHapPopRate[popId][h][0] = self.birthHapPopRate[popId][h]
            
            #self.tEventHapPopRate[popId][h] = sum(self.eventHapPopRate[popId][h])
            self.tEventHapPopRate[popId][h] = 0
            for i in range(self.eventHapPopRate[popId][h].size()):
                self.tEventHapPopRate[popId][h] += self.eventHapPopRate[popId][h][i]

            self.HapPopRate(popId, h)
        
        #self.popRate[popId] = sum( self.hapPopRate[popId] )
        self.popRate[popId] = 0
        for i in range(self.hapPopRate[popId].size()):
            self.popRate[popId] += self.hapPopRate[popId][i]
        
        #self.totalRate = sum( self.popRate )
        self.totalRate = 0
        for i in range(self.popRate.size()):
            self.totalRate += self.popRate[i]

        self.lbCounter += 1

    cdef void Death(self, Py_ssize_t popId, Py_ssize_t haplotype, add_event = True):
        self.liveBranches[popId][haplotype] -= 1

        if add_event:
            event = Event(self.currentTime, 1, popId, haplotype)
            self.events.append(event)

        self.populationModel.populations[popId].infectious -= 1

        self.HapPopRate(popId, haplotype)
        
        #self.popRate[popId] = sum(self.hapPopRate[popId])
        self.popRate[popId] = 0
        for i in range(self.hapPopRate[popId].size()):
            self.popRate[popId] += self.hapPopRate[popId][i]

        #self.totalRate = sum( self.popRate )
        self.totalRate = 0
        for i in range(self.popRate.size()):
            self.totalRate += self.popRate[i]

        self.lbCounter -= 1
        self.susceptible -= 1

    cdef void Sampling(self, Py_ssize_t popId, Py_ssize_t haplotype):
        event = Event(self.currentTime, 2, popId, haplotype)
        self.events.append(event)

        self.Death(popId, haplotype, False)
        self.sCounter += 1

#    def UpdateRate(self):
#        self.totalRate = self.B_rate[0] + self.D_rate[0] + self.S_rate[0] #TODO

    cpdef void SimulatePopulation(self, Py_ssize_t iterations):
        cdef:
            int max_time, sCounter
        max_time = 0
        sCounter = 0
        for j in range(0, iterations):
            #self.UpdateRate()
            self.SampleTime()
            self.GenerateEvent()
            if self.lbCounter == 0 or self.susceptible == 0:
                break
        if self.sCounter < 2: #TODO if number of sampled leaves is 0 (probably 1 as well), then GetGenealogy seems to go to an infinite cycle
            print("Bo-o-o-oring... Less than two cases were sampled.")
            sys.exit(1)

    cpdef void GetGenealogy(self):
        cdef:
            Py_ssize_t lbs, n1, n2, id1, id2, id3
            double p

        for i in range( len( self.populationModel.populations) ):
            self.liveBranchesS.push_back( self.vecint2 )
            for _ in range(self.hapNum):
                self.liveBranchesS[i].push_back( self.vecint1 )

        for event in reversed(self.events):
            if event.type == 0:
                lbs = self.liveBranchesS[event.population][event.haplotype].size()
                p = lbs*(lbs-1)/self.liveBranches[event.population][event.haplotype]/(self.liveBranches[event.population][event.haplotype]-1)
                if np.random.rand() < p:
                    n1 = floor( lbs*np.random.rand() )
                    n2 = floor( (lbs-1)*np.random.rand() )
                    if n2 >= n1:
                        n2 += 1
                    id1 = self.liveBranchesS[event.population][event.haplotype][n1]
                    id2 = self.liveBranchesS[event.population][event.haplotype][n2]
                    id3 = self.tree.size()
                    self.liveBranchesS[event.population][event.haplotype][n1] = id3
                    self.liveBranchesS[event.population][event.haplotype][n2] = self.liveBranchesS[event.population][event.haplotype][-1]
                    self.liveBranchesS[event.population][event.haplotype].pop_back()
                    self.tree[id1] = id3
                    self.tree[id2] = id3
                    self.tree.push_back(-1)
                    self.times.push_back( event.time )
                self.liveBranches[event.population][event.haplotype] -= 1
            elif event.type == 1:
                self.liveBranches[event.population][event.haplotype] += 1
            elif event.type == 2:
                self.liveBranches[event.population][event.haplotype] += 1
                self.liveBranchesS[event.population][event.haplotype].push_back( self.tree.size() )
                self.tree.push_back(-1)
                self.times.push_back( event.time )
            elif event.type == 3:
                lbs = self.liveBranchesS[event.newPopulation][event.haplotype].size()
                p = lbs/self.liveBranches[event.newPopulation][event.haplotype]
                if np.random.rand() < p:
                    n1 = floor( lbs*np.random.rand() )
                    id1 = self.liveBranchesS[event.newPopulation][event.haplotype][n1]
                    self.liveBranchesS[event.newPopulation][event.haplotype][n1] = self.liveBranchesS[event.newPopulation][event.haplotype][-1]
                    self.liveBranchesS[event.newPopulation][event.haplotype].pop_back()
                    self.liveBranchesS[event.population][event.haplotype].push_back(id1)
                self.liveBranches[event.newPopulation][event.haplotype] -= 1
                self.liveBranches[event.population][event.haplotype] += 1
            elif event.type == 4:
                lbs = self.liveBranchesS[event.population][event.newHaplotype].size()
                p = lbs/self.liveBranches[event.population][event.newHaplotype]
                if np.random.rand() < p:
                    n1 = floor( lbs*np.random.rand() )
                    id1 = self.liveBranchesS[event.population][event.newHaplotype][n1]
                    self.liveBranchesS[event.population][event.newHaplotype][n1] = self.liveBranchesS[event.population][event.newHaplotype][-1]
                    self.liveBranchesS[event.population][event.newHaplotype].pop_back()
                    self.liveBranchesS[event.population][event.haplotype].push_back(id1)
                self.liveBranches[event.population][event.newHaplotype] -= 1
                self.liveBranches[event.population][event.haplotype] += 1
            else:
                print("Unknown event type: ", event.type)
                sys.exit(1)

    def LogDynamics(self):
        lg = str(self.currentTime) + " " + str(self.susceptible)
        for pop in self.liveBranches:
            for hap in pop:
                lg += " " + str(len(hap))
        print_err(lg)

    def Report(self):
        print("Number of lineages of each hyplotype: ", end = "")
        for el in self.liveBranches:
            print(len(el), " ", end="")
        print("")
        print("Current time: ", self.currentTime)
        print("Number of sampled elements: ", self.sCounter)
