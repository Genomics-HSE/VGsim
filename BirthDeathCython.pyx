# cython: language_level=3
# distutils: language = c++

cimport cython
from math import floor
import numpy as np

from libc.math cimport log
from libcpp.vector cimport vector

from mc_lib.rndm cimport RndmWrapper


cdef class NodeS():
    cdef:
        int _state, _genealogyIndex;

    def __init__(self):
        self._state = 0
        self._genealogyIndex = -1

    @property
    def state(self):
        return self._state

    @state.setter
    def state(self, state):
        self._state = state

    @property
    def genealogyIndex(self):
        return self._genealogyIndex

    @genealogyIndex.setter
    def genealogyIndex(self, genealogyIndex):
        self._genealogyIndex = genealogyIndex

"""
class NodeS:
    def __init__(self):
        self.state = 0
        self.genealogyIndex = -1
"""

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

"""
class Mutation:
    def __init__(self, nodeId, time, AS, DS):#AS = ancestral state, DS = derived state
        self.nodeId = nodeId
        self.time = time
        self.AS = AS
        self.DS = DS
"""

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

"""
class Population:
    def __init__(self, size = 1000000, contactDensity = 1.0):
        self.size = size
        self.susceptible = size
        self.infectious = 0
        self.contactDensity = contactDensity
"""

cdef class PopulationModel: #Возможно перевести в структуру
    cdef:
        list _populations, _migrationRates

    def __init__(self, list populations, list migrationRates):
        self._populations = populations #list/array with n populations список классов
        self._migrationRates = migrationRates #n * n-1 matrix with migration rates and zeroes on the diagonal список списков
        #self.totalMigrationRate = [sum(mr) for mr in self.migrationRates]

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

"""
class PopulationModel:
    def __init__(self, populations, migrationRates):
        self.populations = populations #list/array with n populations
        self.migrationRates = migrationRates #n * n matrix with migration rates and zeroes on the diagonal
        #self.totalMigrationRate = [sum(mr) for mr in self.migrationRates]
"""

class NeutralMutations:
    def __init__(self):
        self.muRate = 0.0

    def muRate(self):
        return self.muRate

"""
def fastChoose(a, w, tw = None):
    if not tw:
        tw = sum(w)
    rn = tw*np.random.rand()
    i = 0
    total = w[0]
    while total < rn and i < len(a) - 1:
        i += 1
        total += w[i]
    if w[i] == 0.0:
        print("fastChoose() alert: 0-weight sampled")
    return a[i]

def fastRandint(n):
    return floor( n*np.random.rand() )
"""

cdef class BirthDeathModel:

    cdef:
        RndmWrapper rndm

        int sCounter, popNum, dim, hapNum
        double totalRate, currentTime
        Py_ssize_t susceptible, lbCounter
        #mutations, nodeSampling, migPopRate -- std::vector
        vector[int] Tree, genealogy, events
        vector[double] times, genealogyTimes, tPopRate, migPopRate
        vector[vector[vector[int]]] liveBranches
        vector[vector[int]] b
        vector[int] c

        vector[vector[double]] bPopHaplotypeRate, tPopHaplotypeRate
        vector[double] a

        #новый liveBranches (unsigned int) будет фиксированным и двумерным(полностью), но возможно будут редкие увеличения размера

        object populationModel, nodeSampling, mutations, bRate, dRate, sRate, mRate

        object debug

        #populationModel, bRate, dRate, sRate - одномерный неизменяемый
        #mRate, bPopHaplotypeRate, tPopHaplotypeRate - двумерный неизменяемый
        #liveBranches - трёхмерный изменяемый
    def __init__(self, bRate, dRate, sRate, mRate = [], **kwargs):
        #self.Tree = [-1, 0, 0]
        self.Tree.push_back(-1)
        self.Tree.push_back(0)
        self.Tree.push_back(0)
        if "populationModel" in kwargs:
            self.populationModel = kwargs["populationModel"]
        else:
            self.populationModel = PopulationModel([ Population() ], [[0.0]])
        self.susceptible = sum([pop.susceptible for pop in self.populationModel.populations])
        self.popNum = len( self.populationModel.populations )
        #self.genealogy = []
        self.mutations = []
        self.nodeSampling = [NodeS() for _ in range(3)]
        #self.times = [0]*3
        self.times.push_back(0)
        self.times.push_back(0)
        self.times.push_back(0)
        #self.genealogyTimes = []
        #self.liveBranches = [1,2]
        self.currentTime = 0.0
        self.SetRates(bRate, dRate, sRate, mRate)
        self.sCounter = 0 #sample counter
        self.debug = False
        #self.events = []
        if "debug" in kwargs:
            if kwargs["debug"]:
                self.debug = True
        #self.TypeOfMutations = [[0] * len(_M_rate[0])]
        self.rndm = RndmWrapper(seed=(1256, 0))

    cdef Py_ssize_t choice(self, Py_ssize_t len_a, vector[double] w, double tw):
        # fastCoose replacement.
        # XXX: can be *much* faster if w is double[::1] or some such
        cdef:
            double rn = tw*self.rndm.uniform()
            Py_ssize_t idx = 0
            double total = w[0]
        while total < rn and idx < len_a - 1:
            idx += 1
            total += w[idx]
        if w[idx] == 0.0:
            print("fastChoose() alert: 0-weight sampled")
        return idx
 
    cdef int fastRandint(self, int n):
        return floor( n * self.rndm.uniform() )

    def SetRates(self, bRate, dRate, sRate, mRate):
        if len(mRate) > 0:
            self.dim = len(mRate[0])
        else:
            self.dim = 0
        self.hapNum = int(4**self.dim)
        if len(bRate) != self.hapNum or len(dRate) != self.hapNum or len(sRate) != self.hapNum:
            print("BirthDeathModel.SetRates() fatal error: inconsistent dimension")
            quit()
        for el in mRate:
            if len(bRate) != self.hapNum:
                print("BirthDeathModel.SetRates() fatal error: inconsistent dimension")
                quit()
        self.InitLiveBranches()
        self.bRate = bRate
        self.dRate = dRate
        self.sRate = sRate
        self.mRate = mRate

        #self.migPopRate = [sum(mr) for mr in self.populationModel.migrationRates]
        for mr in self.populationModel.migrationRates:
            self.migPopRate.push_back(sum(mr))

        #self.bPopHaplotypeRate = [ [0]*self.hapNum for _ in range( self.popNum ) ]
        for i in range(self.popNum):
            self.bPopHaplotypeRate.push_back(self.a)
            for j in range(self.hapNum):
                self.bPopHaplotypeRate[i].push_back(0)

        self.bPopHaplotypeRate[0][0] = self.BirthRate(0, 0)

        #self.tPopHaplotypeRate = [ [0]*self.hapNum for p in self.populationModel.populations ]
        for i in range( len(self.populationModel.populations) ):
            self.tPopHaplotypeRate.push_back(self.a)
            for j in range(self.hapNum):
                self.tPopHaplotypeRate[i].push_back(0)

        #self.tPopRate = [0]*self.popNum
        for i in range(self.popNum):
            self.tPopRate.push_back(0)

        self.totalPopHaplotypeRate(0, 0)
        self.tPopRate[0] = self.tPopHaplotypeRate[0][0]
        self.totalRate = self.tPopRate[0]


    cdef void totalPopHaplotypeRate(self, Py_ssize_t popId, Py_ssize_t haplotype):
        cdef:
            double r
        r = self.bPopHaplotypeRate[popId][haplotype] + self.dRate[haplotype] + self.sRate[haplotype] + self.migPopRate[popId] + sum( self.mRate[haplotype] )
        self.tPopHaplotypeRate[popId][haplotype] = r*self.liveBranches[popId][haplotype].size()

    cdef double BirthRate(self, Py_ssize_t popId, Py_ssize_t haplotype):
        pop = self.populationModel.populations[popId]
        return self.bRate[haplotype]*pop.susceptible/pop.size*pop.contactDensity

    """
    cdef void InitLiveBranches(self):
        cdef: 
            Py_ssize_t i
        self.liveBranches = []
        for i in range( len( self.populationModel.populations) ):
            self.liveBranches.append( [[] for _ in range(self.hapNum)] )
        self.liveBranches[0][0] += [1,2]
        self.lbCounter = 2 #live branch counter
        self.populationModel.populations[0].infectious = 2
    """

    cdef void InitLiveBranches(self):
        cdef: 
            Py_ssize_t i

        for i in range( len( self.populationModel.populations) ):
            self.liveBranches.push_back(self.b)
            for j in range(self.hapNum):
                self.liveBranches[i].push_back(self.c)

        self.liveBranches[0][0].push_back(1)
        self.liveBranches[0][0].push_back(2)
        self.lbCounter = 2 #live branch counter
        self.populationModel.populations[0].infectious = 2

    cpdef void GenerateEvent(self):
        cdef:
            Py_ssize_t popId, haplotype, eventType, affectedBranch
            double tRate
            vector[double] prbs


        popId = self.choice( len(self.populationModel.populations), self.tPopRate, self.totalRate)
        haplotype = self.choice(self.hapNum, self.tPopHaplotypeRate[popId], self.tPopRate[popId])

        tRate =  self.tPopHaplotypeRate[popId][haplotype] / self.liveBranches[popId][haplotype].size()

        #prbs = [ self.bPopHaplotypeRate[popId][haplotype], self.dRate[haplotype], self.sRate[haplotype], self.migPopRate[popId] ] + self.mRate[haplotype]
        prbs.push_back(self.bPopHaplotypeRate[popId][haplotype])
        prbs.push_back(self.dRate[haplotype])
        prbs.push_back(self.sRate[haplotype])
        prbs.push_back(self.migPopRate[popId])
        for i in range(len(self.mRate[haplotype])):
            prbs.push_back(self.mRate[haplotype][i])

        eventType = self.choice( 4+self.dim, prbs , tRate)
        affectedBranch = self.fastRandint(self.liveBranches[popId][haplotype].size())

        if eventType == 0:
            self.Birth(popId, haplotype, affectedBranch)
        elif eventType == 1:
            self.Death(popId, haplotype, affectedBranch)
        elif eventType == 2:
            self.Sampling(popId, haplotype, affectedBranch)
        elif eventType == 3:
            self.Migration(popId, haplotype, affectedBranch)
        else:
            self.Mutation(popId, haplotype, affectedBranch, eventType - 4)
        self.events.push_back(eventType)#TODO - clean me!!!
        #print("ET=", eventType, "   ", self.liveBranches)

    cdef void Migration(self, Py_ssize_t popId, Py_ssize_t haplotype, Py_ssize_t affectedBranch):
        cdef:
            Py_ssize_t i
            int targetPopId
        targetPopId = self.choice(self.popNum-1, self.populationModel.migrationRates[popId], self.migPopRate[popId])
        if targetPopId >= popId:
            targetPopId += 1
        self.liveBranches[targetPopId][haplotype].push_back( self.liveBranches[popId][haplotype][affectedBranch] )

        self.liveBranches[popId][haplotype][affectedBranch] = self.liveBranches[popId][haplotype][-1]
        self.liveBranches[popId][haplotype].pop_back()

        self.totalPopHaplotypeRate(popId, haplotype)
        self.totalPopHaplotypeRate(targetPopId, haplotype)

        #self.tPopRate[popId] = sum(self.tPopHaplotypeRate[popId])
        self.tPopRate[popId] = 0
        for i in range(self.tPopHaplotypeRate[popId].size()):
            self.tPopRate[popId] += self.tPopHaplotypeRate[popId][i]

        #self.tPopRate[targetPopId] = sum(self.tPopHaplotypeRate[targetPopId])
        self.tPopRate[targetPopId] = 0
        for i in range(self.tPopHaplotypeRate[targetPopId].size()):
            self.tPopRate[targetPopId] += self.tPopHaplotypeRate[targetPopId][i]

        #self.totalRate = sum( self.tPopRate )
        self.totalRate = 0
        for i in range(self.tPopRate.size()):
            self.totalRate += self.tPopRate[i]

    cdef void Mutation(self, Py_ssize_t popId, Py_ssize_t haplotype, Py_ssize_t affectedBranch, Py_ssize_t mutationType):
        cdef: 
            int digit4, AS, DS
        digit4 = 4**mutationType
        AS = floor(haplotype/digit4)%4
        #DS = np.random.choice(range(3))#TODO non-uniform rates???
        DS = self.fastRandint(2) # 3 - 1 = 2
        if DS >= AS:
            DS += 1
        self.mutations.append(Mutation(self.liveBranches[popId][haplotype][affectedBranch], self.currentTime, AS, DS))
        newHaplotype = haplotype + (DS-AS)*digit4

        self.liveBranches[popId][newHaplotype].push_back( self.liveBranches[popId][haplotype][affectedBranch] )

        self.liveBranches[popId][haplotype][affectedBranch] = self.liveBranches[popId][haplotype][-1]
        self.liveBranches[popId][haplotype].pop_back()

        self.totalPopHaplotypeRate(popId, haplotype)
        self.totalPopHaplotypeRate(popId, newHaplotype)
        
        #self.tPopRate[popId] = sum(self.tPopHaplotypeRate[popId])
        self.tPopRate[popId] = 0
        for i in range(self.tPopHaplotypeRate[popId].size()):
            self.tPopRate[popId] += self.tPopHaplotypeRate[popId][i]

        #self.totalRate = sum( self.tPopRate )
        self.totalRate = 0
        for i in range(self.tPopRate.size()):
            self.totalRate += self.tPopRate[i]

    """
    def SampleTime(self):
        tau = np.random.exponential(1.0/self.totalRate)
        self.currentTime += tau
    """

    cdef void SampleTime(self):
        #cdef double tau = np.random.exponential(1.0/self.totalRate)
        cdef double tau = self.totalRate * log(1.0 - self.rndm.uniform())
        self.currentTime += tau

    cdef void Birth(self, Py_ssize_t popId, Py_ssize_t haplotype, Py_ssize_t affectedBranch):
        cdef:
            int parentId
        parentId = self.liveBranches[popId][haplotype][affectedBranch]
        self.times[parentId] = self.currentTime
        self.liveBranches[popId][haplotype].push_back(self.Tree.size())
        self.liveBranches[popId][haplotype][affectedBranch] = self.Tree.size() + 1
        for j in range(2):
            #self.Tree.append(parentId)
            self.Tree.push_back(parentId)
            self.times.push_back(0)
            self.nodeSampling.append( NodeS() )
        self.populationModel.populations[popId].susceptible -= 1
        self.populationModel.populations[popId].infectious += 1

        for h in range(self.hapNum):
            self.bPopHaplotypeRate[popId][h] = self.BirthRate(popId, h)
            self.totalPopHaplotypeRate(popId, h)
        
        #self.tPopRate[popId] = sum( self.tPopHaplotypeRate[popId] )
        self.tPopRate[popId] = 0
        for i in range(self.tPopHaplotypeRate[popId].size()):
            self.tPopRate[popId] += self.tPopHaplotypeRate[popId][i]
        
        #self.totalRate = sum( self.tPopRate )
        self.totalRate = 0
        for i in range(self.tPopRate.size()):
            self.totalRate += self.tPopRate[i]

        self.lbCounter += 1

    cdef void Death(self, Py_ssize_t popId, Py_ssize_t haplotype, Py_ssize_t affectedBranch):
        self.times[ self.liveBranches[popId][haplotype][affectedBranch] ] = self.currentTime
        self.liveBranches[popId][haplotype][affectedBranch] = self.liveBranches[popId][haplotype][-1]
        self.liveBranches[popId][haplotype].pop_back()

        self.populationModel.populations[popId].infectious -= 1

        self.totalPopHaplotypeRate(popId, haplotype)
        
        #self.tPopRate[popId] = sum(self.tPopHaplotypeRate[popId])
        self.tPopRate[popId] = 0
        for i in range(self.tPopHaplotypeRate[popId].size()):
            self.tPopRate[popId] += self.tPopHaplotypeRate[popId][i]

        #self.totalRate = sum( self.tPopRate )
        self.totalRate = 0
        for i in range(self.tPopRate.size()):
            self.totalRate += self.tPopRate[i]

        self.lbCounter -= 1
        self.susceptible -= 1

    cdef void Sampling(self, Py_ssize_t popId, Py_ssize_t haplotype, Py_ssize_t affectedBranch):
        self.nodeSampling[ self.liveBranches[popId][haplotype][affectedBranch] ].state = -1
        self.Death(popId, haplotype, affectedBranch)
        self.sCounter += 1

    cpdef void SimulatePopulation(self, Py_ssize_t iterations):
        cdef:
          double max_time = 0
          int sCounter = 0
          Py_ssize_t j

        for j in range(0, iterations):
            #self.UpdateRate()
            self.SampleTime()
            self.GenerateEvent()
            if self.lbCounter == 0 or self.susceptible == 0:
                break

        if self.debug:
            print(self.Tree)
            print(self.nodeSampling)
        if self.sCounter < 2: #TODO if number of sampled leaves is 0 (probably 1 as well), then GetGenealogy seems to go to an infinite cycle
            print("Bo-o-o-oring... Less than two cases were sampled.")
            quit()
        for popLB in self.liveBranches:
            for lb in popLB:
                for br in lb:
                    self.times[br] = self.currentTime

    cpdef void GetGenealogy(self):
        cdef: 
            Py_ssize_t i
            int parent
        for i in range(self.Tree.size() - 1, 0, -1):
            if self.nodeSampling[i].state == -1:
                parent = self.Tree[i]
                self.nodeSampling[i].genealogyIndex = self.genealogy.size()
                self.genealogy.push_back(-1)
                while parent != -1 and self.nodeSampling[parent].state != 2:
                    self.nodeSampling[parent].state += 1
                    if self.nodeSampling[parent].state == 2:
                        self.nodeSampling[parent].genealogyIndex = self.genealogy.size()
                        self.genealogy.push_back(-1)
                        break
                    parent = self.Tree[parent]
        if self.debug:
            print(self.nodeSampling)

        #self.genealogyTimes = [0]*self.genealogy.size()
        for i in range(self.genealogy.size()):
            self.genealogyTimes.push_back(0)
        for i in range(self.Tree.size() - 1, 0, -1):
            if self.nodeSampling[i].state == -1 or self.nodeSampling[i].state == 2:
                parent = self.Tree[i]
                while parent != -1 and self.nodeSampling[parent].state != 2:
                    parent = self.Tree[parent]
                self.genealogy[ self.nodeSampling[i].genealogyIndex ] = self.nodeSampling[parent].genealogyIndex
                self.genealogyTimes[ self.nodeSampling[i].genealogyIndex ] = self.times[i]

    """
    def LogDynamics(self):
        lg = str(self.currentTime) + " " + str(self.susceptible)
        for pop in self.liveBranches:
            for hap in pop:
                lg += " " + str(len(hap))
        print(lg)

    def Report(self):
        print("Number of lineages of each hyplotype: ")
        for el in self.liveBranches:
            print(len(el), " ")
        print("")
        print("Tree size: ", len(self.Tree))
        print("Number of sampled elements: ", self.sCounter)
    """    
