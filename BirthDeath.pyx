# cython: language_level=3
# distutils: language = c++

cimport cython

from libcpp cimport bool
from libc.math cimport log, floor
from libcpp.vector cimport vector
from mc_lib.rndm cimport RndmWrapper

import numpy as np
np.random.seed(1256)
import sys
import math

def print_err(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

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

        #self.times = [0.0]*self.size
        self.times = np.zeros(self.size, dtype=float)

        #self.types = [0]*self.size
        self.types = np.zeros(self.size, dtype=int)

        #self.populations = [0]*self.size
        self.populations = np.zeros(self.size, dtype=int)

        #self.haplotypes = [0]*self.size
        self.haplotypes = np.zeros(self.size, dtype=int)

        #self.newHaplotypes = [0]*self.size
        self.newHaplotypes = np.zeros(self.size, dtype=int)

        #self.newPopulations = [0]*self.size
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
        self.susceptible = self.size
        self.infectious = 0
        self.contactDensity = contactDensity

cdef class PopulationModel:
    cdef:
        int[::1] sizes, susceptible, infectious
        double[::1] contactDensity

    def __init__(self, populations):
        sizePop = len(populations)

        #self.sizes = [pop.size for pop in populations]
        self.sizes = np.zeros(sizePop, dtype=np.int32)
        for i in range(sizePop):
            self.sizes[i] = populations[i].size

        #self.susceptible = [pop.susceptible for pop in populations]
        self.susceptible = np.zeros(sizePop, dtype=np.int32)
        for i in range(sizePop):
            self.susceptible[i] = populations[i].susceptible

        #self.infectious = [pop.infectious for pop in populations]
        self.infectious = np.zeros(sizePop, dtype=np.int32)
        for i in range(sizePop):
            self.infectious[i] = populations[i].infectious

        #self.contactDensity = [pop.contactDensity for pop in populations]
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
cdef (Py_ssize_t, double) fastChoose1(double[::1] w, double tw, double rn):
    cdef:
        Py_ssize_t i
        double total

    rn = tw*rn
    i = 0
    total = w[0]
    while total < rn and i < len(w) - 1:
        i += 1
        total += w[i]
    if w[i] == 0.0:
        print_err("fastChoose() alert: 0-weight sampled")
    return [ i, ( rn-(total-w[i]) )/w[i] ]

cdef class BirthDeathModel:
    cdef:
        RndmWrapper rndm

        double currentTime, rn, totalRate
        Py_ssize_t sCounter, susceptible, popNum, dim, hapNum, lbCounter
        Events events
        PopulationModel pm

        int[:,::1] liveBranches

        double[::1] bRate, dRate, sRate, tmRate, migPopRate, popRate, elementsArr3
        double[:,::1] pm_migrationRates, birthHapPopRate, tEventHapPopRate, hapPopRate, mRate
        double[:,:,::1] eventHapPopRate

        object tree, times

    def __init__(self, iterations, bRate, dRate, sRate, mRate = [], **kwargs):
        self.currentTime = 0.0
        self.sCounter = 0 #sample counter
        self.events = Events(iterations)

        #Set population model
        if "populationModel" in kwargs:
            self.pm = PopulationModel( kwargs["populationModel"][0] )
            self.pm_migrationRates = np.asarray(kwargs["populationModel"][1])
        else:
            #self.populationModel = PopulationModel([ Population() ], [[]])
            self.pm = PopulationModel( [ Population() ] )
            self.pm_migrationRates = np.asarray((0, 0), dtype=float)

        self.popNum = len( self.pm.sizes )

        #self.susceptible = sum( self.pm.susceptible )
        self.susceptible = 0
        for i in range(self.popNum):
            self.susceptible += self.pm.susceptible[i]

        #Initialise haplotypes
        if len(mRate) > 0:
            self.dim = len(mRate[0])
        else:
            self.dim = 0
        self.hapNum = int(4**self.dim)
        self.InitLiveBranches()


        self.elementsArr3 = np.ones(3)

        #Set rates
        self.SetRates(bRate, dRate, sRate, mRate)

        #Set random generator
        self.rndm = RndmWrapper(seed=(1256, 0))
        
        #Set "GetGenealogy"
        self.tree = []
        self.times = []

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void InitLiveBranches(self):
        
        # for i in range( self.popNum ):
        #     self.liveBranches.append( [0 for _ in range(self.hapNum)] )
        self.liveBranches = np.zeros((self.popNum, self.hapNum), dtype=np.int32)

        self.liveBranches[0][0] += 2
        self.lbCounter = 2 #live branch counter
        self.pm.susceptible[0] -= 2
        self.susceptible -= 2
        self.pm.infectious[0] = 2

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void SetRates(self, bRate, dRate, sRate, mRate):
        self.bRate, self.dRate, self.sRate = np.asarray(bRate), np.asarray(dRate), np.asarray(sRate)
        

        #self.mRate = mRate
        self.mRate = np.zeros((len(mRate), len(mRate[0])), dtype=float)
        for i in range(len(mRate)):
            for j in range(len(mRate[0])):
                self.mRate[i][j] = mRate[i][j]


        #self.tmRate = [ sum(mr) for mr in self.mRate ]
        self.tmRate = np.zeros(len(mRate), dtype=float)
        for i in range(len(self.mRate)):
            for j in range(len(self.mRate[0])):
                self.tmRate[i] += self.mRate[i][j]


        #self.migPopRate = [sum(mr) for mr in self.pm_migrationRates]
        self.migPopRate = np.zeros(len(self.pm_migrationRates), dtype=float)
        for i in range(len(self.pm_migrationRates)):
            for j in range(len(self.pm_migrationRates[0])):
                self.migPopRate[i] += self.pm_migrationRates[i][j]



        #self.birthHapPopRate = [ [0.0]*self.hapNum for _ in range( self.popNum ) ]
        self.birthHapPopRate = np.zeros((self.popNum, self.hapNum), dtype=float)

        #self.eventHapPopRate = [ [ [0.0 for _ in range(5) ] for _ in range(self.hapNum) ] for _ in range(self.popNum) ]
        self.eventHapPopRate = np.zeros((self.popNum, self.hapNum, 5), dtype=float)

        #self.tEventHapPopRate = [ [ 0.0 for _ in range(self.hapNum) ] for _ in range(self.popNum) ]
        self.tEventHapPopRate = np.zeros((self.popNum, self.hapNum), dtype=float)

        for pn in range(self.popNum):
            for hn in range(self.hapNum):
                self.birthHapPopRate[pn][hn] = self.BirthRate(pn, hn)

                #self.eventHapPopRate[pn][hn] = [self.birthHapPopRate[pn][hn], self.dRate[hn], self.sRate[hn], self.migPopRate[pn], self.tmRate[hn] ]
                self.eventHapPopRate[pn][hn][0] = self.birthHapPopRate[pn][hn]
                self.eventHapPopRate[pn][hn][1] = self.dRate[hn]
                self.eventHapPopRate[pn][hn][2] = self.sRate[hn]
                self.eventHapPopRate[pn][hn][3] = self.migPopRate[pn]
                self.eventHapPopRate[pn][hn][4] = self.tmRate[hn]

                #self.tEventHapPopRate[pn][hn] = sum(self.eventHapPopRate[pn][hn])
                for i in range(len(self.eventHapPopRate[pn][hn])):
                    self.tEventHapPopRate[pn][hn] += self.eventHapPopRate[pn][hn][i]

        #self.hapPopRate = [ [0.0 for hn in range(self.hapNum)] for pn in range(self.popNum) ]
        self.hapPopRate = np.zeros((self.popNum, self.hapNum), dtype=float)

        for pn in range(self.popNum):
            for hn in range(self.hapNum):
                self.HapPopRate(pn, hn)

        #self.popRate = [ sum(self.hapPopRate[pn]) for pn in range(self.popNum) ]
        self.popRate = np.zeros(self.popNum, dtype=float)
        for i in range(self.popNum):
            for j in range(len(self.hapPopRate[i])):
                self.popRate[i] += self.hapPopRate[i][j]

        #self.totalRate = sum( self.popRate )
        self.totalRate = 0
        for i in range(len(self.popRate)):
            self.totalRate += self.popRate[i]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void HapPopRate(self, Py_ssize_t popId, Py_ssize_t haplotype):
        self.hapPopRate[popId][haplotype] = self.tEventHapPopRate[popId][haplotype]*self.liveBranches[popId][haplotype]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef double BirthRate(self, Py_ssize_t popId, Py_ssize_t haplotype):
        return self.bRate[haplotype]*self.pm.susceptible[popId]/self.pm.sizes[popId]*self.pm.contactDensity[popId]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void GenerateEvent(self, useNumpy = False):
        cdef:
            Py_ssize_t popId, haplotype, eventType

        #self.rn = np.random.rand()
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

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void Migration(self, Py_ssize_t popId, Py_ssize_t haplotype):
        cdef Py_ssize_t targetPopId
        targetPopId, self.rn = fastChoose1(self.pm_migrationRates[popId], self.migPopRate[popId], self.rn)
        if targetPopId >= popId:
            targetPopId += 1

        self.liveBranches[targetPopId][haplotype] += 1
        self.liveBranches[popId][haplotype] -= 1

        self.pm.susceptible[popId] += 1
        self.pm.infectious[popId] -= 1

        self.pm.susceptible[targetPopId] -= 1
        self.pm.infectious[targetPopId] += 1

        #event = Event(self.currentTime, 3, popId, haplotype, newPopulation = targetPopId)
        #self.events.append(event)
        self.events.AddEvent(self.currentTime, 3, popId, haplotype, 0, targetPopId)

        self.HapPopRate(popId, haplotype)
        self.HapPopRate(targetPopId, haplotype)
        
        #self.popRate[popId] = sum(self.hapPopRate[popId])
        self.popRate[popId] = 0
        for i in range(len(self.hapPopRate[popId])):
            self.popRate[popId] += self.hapPopRate[popId][i]

        #self.popRate[targetPopId] = sum(self.hapPopRate[targetPopId])
        self.popRate[targetPopId] = 0
        for i in range(len(self.hapPopRate[targetPopId])):
            self.popRate[targetPopId] += self.hapPopRate[targetPopId][i]

        #self.totalRate = sum( self.popRate )
        self.totalRate = 0
        for i in range(len(self.popRate)):
            self.totalRate += self.popRate[i]

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

        self.liveBranches[popId][newHaplotype] += 1
        self.liveBranches[popId][haplotype] -= 1

        #event = Event(self.currentTime, 4, popId, haplotype, newHaplotype = newHaplotype)
        #self.events.append(event)
        self.events.AddEvent(self.currentTime, 4, popId, haplotype, newHaplotype, 0)

        self.HapPopRate(popId, haplotype)
        self.HapPopRate(popId, newHaplotype)

        #self.popRate[popId] = sum(self.hapPopRate[popId])
        self.popRate[popId] = 0
        for i in range(len(self.hapPopRate[popId])):
            self.popRate[popId] += self.hapPopRate[popId][i]

        #self.totalRate = sum( self.popRate )
        self.totalRate = 0
        for i in range(len(self.popRate)):
            self.totalRate += self.popRate[i]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void SampleTime(self):
        #cdef double tau = - log(np.random.rand()) / self.totalRate
        cdef double tau = - log(self.rndm.uniform()) / self.totalRate
        self.currentTime += tau

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void Birth(self, Py_ssize_t popId, Py_ssize_t haplotype):
        self.liveBranches[popId][haplotype] += 1

        self.pm.susceptible[popId] -= 1
        self.pm.infectious[popId] += 1

        #event = Event(self.currentTime, 0, popId, haplotype)
        #self.events.append(event)
        self.events.AddEvent(self.currentTime, 0, popId, haplotype, 0, 0)

        for h in range(self.hapNum):
            self.birthHapPopRate[popId][h] = self.BirthRate(popId, h)
            self.eventHapPopRate[popId][h][0] = self.birthHapPopRate[popId][h]

            #self.tEventHapPopRate[popId][h] = sum(self.eventHapPopRate[popId][h])
            self.tEventHapPopRate[popId][h] = 0
            for i in range(len(self.eventHapPopRate[popId][h])):
                self.tEventHapPopRate[popId][h] += self.eventHapPopRate[popId][h][i]

            self.HapPopRate(popId, h)
        
        #self.popRate[popId] = sum(self.hapPopRate[popId])
        self.popRate[popId] = 0
        for i in range(len(self.hapPopRate[popId])):
            self.popRate[popId] += self.hapPopRate[popId][i]

        #self.totalRate = sum( self.popRate )
        self.totalRate = 0
        for i in range(len(self.popRate)):
            self.totalRate += self.popRate[i]

        self.lbCounter += 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void Death(self, Py_ssize_t popId, Py_ssize_t haplotype, bool add_event = True):
        self.liveBranches[popId][haplotype] -= 1

        if add_event:
            #event = Event(self.currentTime, 1, popId, haplotype)
            #self.events.append(event)
            self.events.AddEvent(self.currentTime, 1, popId, haplotype, 0, 0)

        self.pm.infectious[popId] -= 1

        self.HapPopRate(popId, haplotype)
        
        #self.popRate[popId] = sum(self.hapPopRate[popId])
        self.popRate[popId] = 0
        for i in range(len(self.hapPopRate[popId])):
            self.popRate[popId] += self.hapPopRate[popId][i]

        #self.totalRate = sum( self.popRate )
        self.totalRate = 0
        for i in range(len(self.popRate)):
            self.totalRate += self.popRate[i]

        self.lbCounter -= 1
        self.susceptible -= 1

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
            if self.lbCounter == 0 or self.susceptible == 0:
                break
        if self.sCounter < 2: #TODO if number of sampled leaves is 0 (probably 1 as well), then GetGenealogy seems to go to an infinite cycle
            print("Bo-o-o-oring... Less than two cases were sampled.")
            sys.exit(1)

    def GetGenealogy(self):
        liveBranchesS = []
        for i in range( self.popNum ):
            liveBranchesS.append( [[] for _ in range(self.hapNum)] )

        #for event in reversed(self.events):
        for e_id in range(self.events.ptr-1, 0, -1):
            event = self.events.GetEvent(e_id)
            if event.type_ == 0:
                lbs = len( liveBranchesS[event.population][event.haplotype] )
                p = lbs*(lbs-1)/self.liveBranches[event.population][event.haplotype]/(self.liveBranches[event.population][event.haplotype]-1)
                if np.random.rand() < p:
                    n1 = math.floor( lbs*np.random.rand() )
                    n2 = math.floor( (lbs-1)*np.random.rand() )
                    if n2 >= n1:
                        n2 += 1
                    id1 = liveBranchesS[event.population][event.haplotype][n1]
                    id2 = liveBranchesS[event.population][event.haplotype][n2]
                    id3 = len( self.tree )
                    liveBranchesS[event.population][event.haplotype][n1] = id3
                    liveBranchesS[event.population][event.haplotype][n2] = liveBranchesS[event.population][event.haplotype][-1]
                    liveBranchesS[event.population][event.haplotype].pop()
                    self.tree[id1] = id3
                    self.tree[id2] = id3
                    self.tree.append(-1)
                    self.times.append( event.time )
                self.liveBranches[event.population][event.haplotype] -= 1
            elif event.type_ == 1:
                self.liveBranches[event.population][event.haplotype] += 1
            elif event.type_ == 2:
                self.liveBranches[event.population][event.haplotype] += 1
                liveBranchesS[event.population][event.haplotype].append( len(self.tree) )
                self.tree.append(-1)
                self.times.append( event.time )
            elif event.type_ == 3:
                lbs = len( liveBranchesS[event.newPopulation][event.haplotype] )
                p = lbs/self.liveBranches[event.newPopulation][event.haplotype]
                if np.random.rand() < p:
                    n1 = math.floor( lbs*np.random.rand() )
                    id1 = liveBranchesS[event.newPopulation][event.haplotype][n1]
                    liveBranchesS[event.newPopulation][event.haplotype][n1] = liveBranchesS[event.newPopulation][event.haplotype][-1]
                    liveBranchesS[event.newPopulation][event.haplotype].pop()
                    liveBranchesS[event.population][event.haplotype].append(id1)
                self.liveBranches[event.newPopulation][event.haplotype] -= 1
                self.liveBranches[event.population][event.haplotype] += 1
            elif event.type_ == 4:
                lbs = len( liveBranchesS[event.population][event.newHaplotype] )
                p = lbs/self.liveBranches[event.population][event.newHaplotype]
                if np.random.rand() < p:
                    n1 = math.floor( lbs*np.random.rand() )
                    id1 = liveBranchesS[event.population][event.newHaplotype][n1]
                    liveBranchesS[event.population][event.newHaplotype][n1] = liveBranchesS[event.population][event.newHaplotype][-1]
                    liveBranchesS[event.population][event.newHaplotype].pop()
                    liveBranchesS[event.population][event.haplotype].append(id1)
                self.liveBranches[event.population][event.newHaplotype] -= 1
                self.liveBranches[event.population][event.haplotype] += 1
            else:
                print("Unknown event type: ", event.type_)
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
        print("popNum: ", self.popNum)
        print("dim: ", self.dim)
        print("hapNum: ", self.hapNum)
        print("totalRate: ", self.totalRate)
        print("rn: ", self.rn)
        print("susceptible: ", self.susceptible)
        print("lbCounter: ", self.lbCounter)
        print("Current time: ", self.currentTime)
        print("Tree size: ", len(self.tree))
        print("Number of sampled elements: ", self.sCounter)
