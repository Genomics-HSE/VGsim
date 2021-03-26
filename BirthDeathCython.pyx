# cython: language_level=3
# cython: initializedcheck = False
# distutils: language = c++

cimport cython

from libc.math cimport log, floor
from libcpp.vector cimport vector
from mc_lib.rndm cimport RndmWrapper

import numpy as np
np.random.seed(1256)
import sys


def print_err(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


# use named constants for event types
DEF BIRTH = 0
DEF DEATH = 1
DEF SAMPLING = 2
DEF MIGRATION = 3
DEF MUTATION = 4


cdef class Event:
    cdef:
        double time
        Py_ssize_t type_, population, haplotype, newHaplotype, newPopulation

    def __init__(self, double time, Py_ssize_t type_, Py_ssize_t population, Py_ssize_t haplotype, Py_ssize_t newHaplotype, Py_ssize_t newPopulation):
        self.time = time
        self.type_ = type_
        self.population = population
        self.haplotype = haplotype
        self.newHaplotype = newHaplotype
        self.newPopulation = newPopulation


cdef class Events:
    cdef:
        double[::1] times
        Py_ssize_t size, ptr
        Py_ssize_t[::1] types, populations, haplotypes, newHaplotypes, newPopulations

    def __init__(self, Py_ssize_t size_):
        self.size = size_
        self.ptr = 0#pointer to the first empty cell

        self.times = np.zeros(self.size, dtype=float)
        self.types = np.zeros(self.size, dtype=int)
        self.populations = np.zeros(self.size, dtype=int)
        self.haplotypes = np.zeros(self.size, dtype=int)
        self.newHaplotypes = np.zeros(self.size, dtype=int)
        self.newPopulations = np.zeros(self.size, dtype=int)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void AddEvent(self, double time_, Py_ssize_t type_, Py_ssize_t population, Py_ssize_t haplotype, Py_ssize_t newHaplotype, Py_ssize_t newPopulation):
        self.times[ self.ptr ] = time_
        self.types[ self.ptr ] = type_
        self.populations[ self.ptr ] = population
        self.haplotypes[ self.ptr ] = haplotype
        self.newHaplotypes[ self.ptr ] = newHaplotype
        self.newPopulations[ self.ptr ] = newPopulation
        self.ptr += 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef Event GetEvent(self, Py_ssize_t e_id):
        ev = Event( self.times[ e_id ], self.types[ e_id ], self.populations[ e_id ], self.haplotypes[ e_id ], self.newHaplotypes[ e_id ], self.newPopulations[ e_id ])
        return( ev )

class Mutation:
    def __init__(self, nodeId, time, AS, DS):#AS = ancestral state, DS = derived state
        self.nodeId = nodeId
        self.time = time
        self.AS = AS
        self.DS = DS

class Population:
    def __init__(self, size = 1000000, contactDensity = 1.0):
        self.size = size
        self.contactDensity = contactDensity

cdef class PopulationModel:
    cdef:
        int[::1] sizes
        int[:,::1] susceptible
        double[::1] contactDensity

    def __init__(self, populations, susceptible_num):
        sizePop = len(populations)

        self.sizes = np.zeros(sizePop, dtype=np.int32)
        for i in range(sizePop):
            self.sizes[i] = populations[i].size

        self.susceptible = np.zeros((sizePop, susceptible_num), dtype=np.int32)
        for i in range(sizePop):
            self.susceptible[i, 0] = populations[i].size

        self.contactDensity = np.zeros(sizePop, dtype=float)
        for i in range(sizePop):
            self.contactDensity[i] = populations[i].contactDensity

class NeutralMutations:
    def __init__(self):
        self.muRate = 0.0

    def muRate(self):
        return self.muRate

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline (Py_ssize_t, double) fastChoose1(double[::1] w, double tw, double rn):
    cdef:
        Py_ssize_t i
        double total

    rn = tw*rn
    i = 0
    total = w[0]
    while total < rn and i < w.shape[0] - 1:
        i += 1
        total += w[i]
    if w[i] == 0.0:
        print_err("fastChoose() alert: 0-weight sampled")
    return [ i, ( rn-(total-w[i]) )/w[i] ]


cdef class BirthDeathModel:
    cdef:
        RndmWrapper rndm

        double currentTime, rn, totalRate
        Py_ssize_t sCounter, popNum, dim, hapNum, lbCounter, susceptible_num
        Events events
        PopulationModel pm

        int[::1] tree, suscType
        int[:,::1] liveBranches

        double[::1] bRate, dRate, sRate, tmRate, migPopRate, popRate, elementsArr3, times
        double[:,::1] pm_migrationRates, birthHapPopRate, tEventHapPopRate, hapPopRate, mRate, susceptibility
        double[:,:,::1] eventHapPopRate, susceptHapPopRate

    def __init__(self, iterations, bRate, dRate, sRate, mRate = [], populationModel = None, susceptible = None, rndseed = 1256, **kwargs):
        self.currentTime = 0.0
        self.sCounter = 0 #sample counter
        self.events = Events(iterations+1)

        ##BEGIN:TOVADIM
        if susceptible == None:
            self.susceptible_num = 2
        else:
            self.susceptible_num = len( susceptible[0] )#CHECK

        #Set population model
        if populationModel == None:
            self.pm = PopulationModel( [ Population() ], self.susceptible_num)
            self.pm_migrationRates = np.asarray((0, 0), dtype=float)
        else:
            self.pm = PopulationModel( populationModel[0] , self.susceptible_num)
            self.pm_migrationRates = np.asarray(populationModel[1])
        ##END:TOVADIM

        self.popNum = self.pm.sizes.shape[0]

        #Initialise haplotypes
        if len(mRate) > 0:
            self.dim = len(mRate[0])
        else:
            self.dim = 0
        self.hapNum = int(4**self.dim)

        self.InitLiveBranches()

        self.elementsArr3 = np.ones(3)

        ##BEGIN:TOVADIM
        if susceptible == None:
            self.susceptibility = np.asarray( [ [1.0 for _ in range(self.hapNum)], [0.0 for _ in range(self.hapNum)] ] )
            self.suscType = np.asarray( [1 for _ in range(self.hapNum)], dtype=np.int32 )
        else:
            self.susceptibility = np.asarray( susceptible[0] )
            self.suscType = np.asarray( susceptible[1], dtype=np.int32 )
        ##END:TOVADIM

        self.susceptHapPopRate = np.zeros((self.popNum, self.hapNum, self.susceptible_num), dtype=float)
        for i in range(self.popNum):
            for j in range(self.hapNum):
                for k in range(self.susceptible_num):
                  print(self.pm.susceptible[i, k])
                  print(self.susceptibility[j, k])
                  print(self.pm.susceptible.shape)
                  print(self.susceptibility.shape)
                  print("Index: ", i, j, k)
                  self.susceptHapPopRate[i, j, k] = self.pm.susceptible[i, k]*self.susceptibility[j, k]

        #Set rates
        self.SetRates(bRate, dRate, sRate, mRate)

        #Set random generator
        self.rndm = RndmWrapper(seed=(rndseed, 0))

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void InitLiveBranches(self):

        # for i in range( self.popNum ):
        #     self.liveBranches.append( [0 for _ in range(self.hapNum)] )
        self.liveBranches = np.zeros((self.popNum, self.hapNum), dtype=np.int32)
        self.events.AddEvent(self.currentTime, 0, 0, 0, 0, 0)
        self.liveBranches[0, 0] += 2
        self.lbCounter = 2 #live branch counter
        self.pm.susceptible[0,0] -= 2

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void SetRates(self, bRate, dRate, sRate, mRate):
        self.bRate, self.dRate, self.sRate = np.asarray(bRate), np.asarray(dRate), np.asarray(sRate)


        #self.mRate = mRate
        self.mRate = np.zeros((len(mRate), len(mRate[0])), dtype=float)
        for i in range(len(mRate)):
            for j in range(len(mRate[0])):
                self.mRate[i, j] = mRate[i][j]


        #self.tmRate = [ sum(mr) for mr in self.mRate ]
        self.tmRate = np.zeros(len(mRate), dtype=float)
        for i in range(self.mRate.shape[0]):
            for j in range(self.mRate.shape[1]):
                self.tmRate[i] += self.mRate[i, j]


        #self.migPopRate = [sum(mr) for mr in self.pm_migrationRates]
        self.migPopRate = np.zeros(len(self.pm_migrationRates), dtype=float)
        for i in range(len(self.pm_migrationRates)):
            for j in range(len(self.pm_migrationRates[0])):
                self.migPopRate[i] += self.pm_migrationRates[i, j]



        #self.birthHapPopRate = [ [0.0]*self.hapNum for _ in range( self.popNum ) ]
        self.birthHapPopRate = np.zeros((self.popNum, self.hapNum), dtype=float)

        #self.eventHapPopRate = [ [ [0.0 for _ in range(5) ] for _ in range(self.hapNum) ] for _ in range(self.popNum) ]
        self.eventHapPopRate = np.zeros((self.popNum, self.hapNum, 5), dtype=float)

        #self.tEventHapPopRate = [ [ 0.0 for _ in range(self.hapNum) ] for _ in range(self.popNum) ]
        self.tEventHapPopRate = np.zeros((self.popNum, self.hapNum), dtype=float)

        for pn in range(self.popNum):
            for hn in range(self.hapNum):
                self.birthHapPopRate[pn, hn] = self.BirthRate(pn, hn)

                #self.eventHapPopRate[pn][hn] = [self.birthHapPopRate[pn][hn], self.dRate[hn], self.sRate[hn], self.migPopRate[pn], self.tmRate[hn] ]
                self.eventHapPopRate[pn, hn, 0] = self.birthHapPopRate[pn, hn]
                self.eventHapPopRate[pn, hn, 1] = self.dRate[hn]
                self.eventHapPopRate[pn, hn, 2] = self.sRate[hn]
                self.eventHapPopRate[pn, hn, 3] = self.migPopRate[pn]
                self.eventHapPopRate[pn, hn, 4] = self.tmRate[hn]

                #self.tEventHapPopRate[pn][hn] = sum(self.eventHapPopRate[pn][hn])
                for i in range(5):
                    self.tEventHapPopRate[pn, hn] += self.eventHapPopRate[pn, hn, i]

        #self.hapPopRate = [ [0.0 for hn in range(self.hapNum)] for pn in range(self.popNum) ]
        self.hapPopRate = np.zeros((self.popNum, self.hapNum), dtype=float)

        for pn in range(self.popNum):
            for hn in range(self.hapNum):
                self.HapPopRate(pn, hn)

        #self.popRate = [ sum(self.hapPopRate[pn]) for pn in range(self.popNum) ]
        self.popRate = np.zeros(self.popNum, dtype=float)
        for i in range(self.popNum):
            for j in range(self.hapNum):
                self.popRate[i] += self.hapPopRate[i, j]

        #self.totalRate = sum( self.popRate )
        self.totalRate = 0
        for i in range(self.popNum):
            self.totalRate += self.popRate[i]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef inline void HapPopRate(self, Py_ssize_t popId, Py_ssize_t haplotype):
        self.hapPopRate[popId, haplotype] = self.tEventHapPopRate[popId, haplotype]*self.liveBranches[popId, haplotype]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef inline double BirthRate(self, Py_ssize_t popId, Py_ssize_t haplotype):
        cdef double ws = 0.0
        for i in range(self.susceptible_num):
            self.susceptHapPopRate[popId, haplotype, i] = self.pm.susceptible[popId,i]*self.susceptibility[haplotype, i]
            ws += self.susceptHapPopRate[popId, haplotype, i]#TOVADIM
        return self.bRate[haplotype]*ws/self.pm.sizes[popId]*self.pm.contactDensity[popId]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void GenerateEvent(self, useNumpy = False):
        cdef:
            Py_ssize_t popId, haplotype, eventType

        #self.rn = np.random.rand()
        self.rn = self.rndm.uniform()

        popId, self.rn = fastChoose1( self.popRate, self.totalRate, self.rn)
        haplotype, self.rn = fastChoose1( self.hapPopRate[popId], self.popRate[popId], self.rn)
        eventType, self.rn = fastChoose1( self.eventHapPopRate[popId, haplotype], self.tEventHapPopRate[popId, haplotype], self.rn)

        if eventType == BIRTH:
            self.Birth(popId, haplotype)
        elif eventType == DEATH:
            self.Death(popId, haplotype)
        elif eventType == SAMPLING:
            self.Sampling(popId, haplotype)
        elif eventType == MIGRATION:
            self.Migration(popId, haplotype)
        else:
            self.Mutation(popId, haplotype)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void Migration(self, Py_ssize_t popId, Py_ssize_t haplotype):
        cdef Py_ssize_t targetPopId
        targetPopId, self.rn = fastChoose1(self.pm_migrationRates[popId], self.migPopRate[popId], self.rn)
        if targetPopId >= popId:
            targetPopId += 1

        self.liveBranches[targetPopId, haplotype] += 1

        cdef double ws = 0.0
        for i in range(self.susceptible_num):
            ws += self.susceptHapPopRate[targetPopId, haplotype, i]#TOVADIM
        st, self.rn = fastChoose1(self.susceptHapPopRate[targetPopId, haplotype], ws, self.rn)
        self.pm.susceptible[targetPopId, st] -= 1
        self.UpdateRates(targetPopId)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void Mutation(self, Py_ssize_t popId, Py_ssize_t haplotype):
        cdef:
            Py_ssize_t mutationType, digit4, AS, DS, newHaplotype

        mutationType, self.rn = fastChoose1( self.mRate[haplotype], self.tmRate[haplotype], self.rn)
        digit4 = 4**mutationType
        AS = int(floor(haplotype/digit4) % 4)
        DS, self.rn = fastChoose1(self.elementsArr3, 3.0, self.rn)#TODO non-uniform rates???
        if DS >= AS:
            DS += 1
        #self.mutations.append(Mutation(self.liveBranches[popId][haplotype][affectedBranch], self.currentTime, AS, DS))
        newHaplotype = haplotype + (DS-AS)*digit4

        self.liveBranches[popId, newHaplotype] += 1
        self.liveBranches[popId, haplotype] -= 1

        #event = Event(self.currentTime, 4, popId, haplotype, newHaplotype = newHaplotype)
        #self.events.append(event)
        self.events.AddEvent(self.currentTime, 4, popId, haplotype, newHaplotype, 0)

        self.HapPopRate(popId, haplotype)
        self.HapPopRate(popId, newHaplotype)

        #self.popRate[popId] = sum(self.hapPopRate[popId])
        self.popRate[popId] = 0
        for i in range(self.hapNum):
            self.popRate[popId] += self.hapPopRate[popId, i]

        #self.totalRate = sum( self.popRate )
        self.totalRate = 0
        for i in range(self.popNum):
            self.totalRate += self.popRate[i]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef inline void SampleTime(self):
        cdef double tau = - log(self.rndm.uniform()) / self.totalRate
        self.currentTime += tau

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void Birth(self, Py_ssize_t popId, Py_ssize_t haplotype):
        self.liveBranches[popId, haplotype] += 1

        cdef double ws = 0.0
        for i in range(self.susceptible_num):
            ws += self.susceptHapPopRate[popId, haplotype, i]#TOVADIM
        st, self.rn = fastChoose1(self.susceptHapPopRate[popId, haplotype], ws, self.rn)

        self.pm.susceptible[popId, st] -= 1

        #event = Event(self.currentTime, 0, popId, haplotype)
        #self.events.append(event)
        self.events.AddEvent(self.currentTime, 0, popId, haplotype, 0, 0)
        self.UpdateRates(popId)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void UpdateRates(self, Py_ssize_t popId):
        cdef double tmp

        for h in range(self.hapNum):
            self.birthHapPopRate[popId, h] = self.BirthRate(popId, h)
            self.eventHapPopRate[popId, h, 0] = self.birthHapPopRate[popId, h]
            tmp = (self.eventHapPopRate[popId, h, 0] +
                   self.eventHapPopRate[popId, h, 1] +
                   self.eventHapPopRate[popId, h, 2] +
                   self.eventHapPopRate[popId, h, 3] +
                   self.eventHapPopRate[popId, h, 4] )
            self.tEventHapPopRate[popId, h] = tmp
            self.HapPopRate(popId, h)

        self.popRate[popId] = 0
        for i in range(self.hapNum):
            self.popRate[popId] += self.hapPopRate[popId, i]

        self.totalRate = 0
        for i in range(self.popNum):
            self.totalRate += self.popRate[i]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void Death(self, Py_ssize_t popId, Py_ssize_t haplotype, bint add_event = True):
        self.liveBranches[popId, haplotype] -= 1
        self.pm.susceptible[popId, self.suscType[haplotype] ] += 1

        if add_event:
            #event = Event(self.currentTime, 1, popId, haplotype)
            #self.events.append(event)
            self.events.AddEvent(self.currentTime, 1, popId, haplotype, 0, 0)

        self.UpdateRates(popId)


    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void Sampling(self, Py_ssize_t popId, Py_ssize_t haplotype):
        #event = Event(self.currentTime, 2, popId, haplotype)
        #self.events.append(event)
        self.events.AddEvent(self.currentTime, 2, popId, haplotype, 0, 0)

        self.Death(popId, haplotype, False)
        self.sCounter += 1

#    def UpdateRate(self):
#        self.totalRate = self.B_rate[0] + self.D_rate[0] + self.S_rate[0] #TODO

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cpdef void SimulatePopulation(self, Py_ssize_t iterations):
        max_time = 0
        sCounter = 0
        for j in range(0, iterations):
            #self.UpdateRate()
            self.SampleTime()
            self.GenerateEvent()
            if self.totalRate == 0.0:
                break
        print("Total number of iterations: ", self.events.ptr)
        if self.sCounter < 2: #TODO if number of sampled leaves is 0 (probably 1 as well), then GetGenealogy seems to go to an infinite cycle
            print("Less than two cases were sampled...")
            sys.exit(1)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cpdef GetGenealogy(self):
        cdef:
            Py_ssize_t ptrTreeAndTime, n1, n2, id1, id2, id3, lbs, lbs_e
            double p
            vector[vector[vector[Py_ssize_t]]] liveBranchesS
            vector[vector[Py_ssize_t]] vecint2
            vector[Py_ssize_t] vecint1

            double e_time
            Py_ssize_t e_type_, e_population, e_haplotype, e_newHaplotype, e_newPopulation

        ptrTreeAndTime = 0
        self.tree = np.zeros(2 * self.sCounter - 1, dtype=np.int32)
        self.times = np.zeros(2 * self.sCounter - 1, dtype=float)

        # liveBranchesS = []
        # for i in range( self.popNum ):
        #     liveBranchesS.append( [[] for _ in range(self.hapNum)] )
        for i in range( self.popNum ):
            liveBranchesS.push_back(vecint2)
            for _ in range( self.hapNum ):
                liveBranchesS[i].push_back(vecint1)

        #for event in reversed(self.events):
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
                p = lbs*(lbs-1)/ lbs_e / (lbs_e - 1)

                # if np.random.rand() < p:
                #     n1 = int(floor( lbs*np.random.rand() ))
                #     n2 = int(floor( (lbs-1)*np.random.rand() ))
                if self.rndm.uniform() < p:
                    n1 = int(floor( lbs*self.rndm.uniform() ))
                    n2 = int(floor( (lbs-1)*self.rndm.uniform() ))

                    if n2 >= n1:
                        n2 += 1
                    id1 = liveBranchesS[e_population][e_haplotype][n1]
                    id2 = liveBranchesS[e_population][e_haplotype][n2]

                    #id3 = len( self.tree )
                    id3 = ptrTreeAndTime

                    liveBranchesS[e_population][e_haplotype][n1] = id3
                    liveBranchesS[e_population][e_haplotype][n2] = liveBranchesS[e_population][e_haplotype][lbs-1]
                    liveBranchesS[e_population][e_haplotype].pop_back()
                    self.tree[id1] = id3
                    self.tree[id2] = id3

                    #self.tree.append(-1)
                    #self.times.append( event.time )
                    self.tree[ptrTreeAndTime] = -1
                    self.times[ptrTreeAndTime] = e_time
                    ptrTreeAndTime += 1

                self.liveBranches[e_population][e_haplotype] -= 1
            elif e_type_ == DEATH:
                self.liveBranches[e_population][e_haplotype] += 1
            elif e_type_ == SAMPLING:
                self.liveBranches[e_population][e_haplotype] += 1
                liveBranchesS[e_population][e_haplotype].push_back( ptrTreeAndTime )

                # self.tree.append(-1)
                # self.times.append( event.time )
                self.tree[ptrTreeAndTime] = -1
                self.times[ptrTreeAndTime] = e_time
                ptrTreeAndTime += 1

            elif e_type_ == MIGRATION:
                lbs = liveBranchesS[e_newPopulation][e_haplotype].size()
                p = lbs/self.liveBranches[e_newPopulation][e_haplotype]

                # if np.random.rand() < p:
                #     n1 = int(floor( lbs*np.random.rand() ))
                if self.rndm.uniform() < p:
                    nt = int(floor( lbs*self.rndm.uniform() ))

                    lbss = liveBranchesS[e_population][e_haplotype].size()
                    ns = int(floor( lbss*self.rndm.uniform() ))

                    idt = liveBranchesS[e_newPopulation][e_haplotype][nt]
                    ids = liveBranchesS[e_population][e_haplotype][ns]

                    id3 = ptrTreeAndTime
                    liveBranchesS[e_population][e_haplotype][ns] = id3
                    liveBranchesS[e_newPopulation][e_haplotype][nt] = liveBranchesS[e_newPopulation][e_haplotype][lbs-1]
                    liveBranchesS[e_newPopulation][e_haplotype].pop_back()
                    self.tree[id1] = id3
                    self.tree[id2] = id3

                    self.tree[ptrTreeAndTime] = -1
                    self.times[ptrTreeAndTime] = e_time
                    ptrTreeAndTime += 1
                self.liveBranches[e_newPopulation][e_haplotype] -= 1
            elif e_type_ == MUTATION:
                lbs = liveBranchesS[e_population][e_newHaplotype].size()
                p = lbs/self.liveBranches[e_population][e_newHaplotype]

                # if np.random.rand() < p:
                #     n1 = int(floor( lbs*np.random.rand() ))
                if self.rndm.uniform() < p:
                    n1 = int(floor( lbs*self.rndm.uniform() ))

                    id1 = liveBranchesS[e_population][e_newHaplotype][n1]
                    liveBranchesS[e_population][e_newHaplotype][n1] = liveBranchesS[e_population][e_newHaplotype][lbs-1]
                    liveBranchesS[e_population][e_newHaplotype].pop_back()
                    liveBranchesS[e_population][e_haplotype].push_back(id1)
                self.liveBranches[e_population][e_newHaplotype] -= 1
                self.liveBranches[e_population][e_haplotype] += 1
            else:
                print("Unknown event type: ", e_type_)
                sys.exit(1)

    def LogDynamics(self, fn, step_num = 1000):
        count = 0
        for e_id in range(self.events.ptr-1, -1, -1):
            e_time = self.events.times[e_id]



    def Report(self):
        print("Number of lineages of each hyplotype: ", end = "")
        for el in self.liveBranches:
            print(len(el), " ", end="")
        print("")
        print("popNum: ", self.popNum)
        print("dim: ", self.dim)
        print("hapNum: ", self.hapNum)
        print("totalRate: ", self.totalRate)
        print("rn: ", self.rn)
        #print("susceptible: ", self.susceptible)
        print("lbCounter: ", self.lbCounter)
        print("Current time: ", self.currentTime)
        print("Tree size: ", len(self.tree))
        print("Number of sampled elements: ", self.sCounter)
        print("")
