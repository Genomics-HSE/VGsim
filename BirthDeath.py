import numpy as np
import random
import sys
from math import floor

def print_err(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class NodeS: #C-like structure
    def __init__(self):
        self.state = 0
        self.genealogyIndex = -1

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

class PopulationModel:
    def __init__(self, populations, migrationRates):
        self.populations = populations #list/array with n populations
        self.migrationRates = migrationRates #n * n matrix with migration rates and zeroes on the diagonal
        #self.totalMigrationRate = [sum(mr) for mr in self.migrationRates]

class NeutralMutations:
    def __init__(self):
        self.muRate = 0.0

    def muRate(self):
        return self.muRate

def fastChoose(a, w, tw = None):
    if not tw:
        tw = sum(w)
    rn = tw*np.random.rand()
    i = 0
    total = w[0]
    while total < rn and i < len(a) - 1:
        i += 1
#        if i == len(w):
#            return -1
        total += w[i]
    if w[i] == 0.0:
        print_err("fastChoose() alert: 0-weight sampled")
    return a[i]

def fastRandint(n):
    return floor( n*np.random.rand() )

class BirthDeathModel:
    def __init__(self, bRate, dRate, sRate, mRate = [], **kwargs):
        self.Tree = [-1, 0, 0]
        if "populationModel" in kwargs:
            self.populationModel = kwargs["populationModel"]
        else:
            self.populationModel = PopulationModel([ Population() ], [[0.0]])
        self.susceptible = sum([pop.susceptible for pop in self.populationModel.populations])
        self.popNum = len( self.populationModel.populations )
        self.genealogy = []
        self.mutations = []
        self.nodeSampling = [NodeS() for _ in range(3)]
        self.times = [0]*3
        self.genealogyTimes = []
        #self.liveBranches = [1,2]
        self.currentTime = 0.0
        self.SetRates(bRate, dRate, sRate, mRate)
        self.sCounter = 0 #sample counter
        self.debug = False
        self.events = []
        if "debug" in kwargs:
            if kwargs["debug"]:
                self.debug = True
        #self.TypeOfMutations = [[0] * len(_M_rate[0])]

    def SetRates(self, bRate, dRate, sRate, mRate):
        if len(mRate) > 0:
            self.dim = len(mRate[0])
        else:
            self.dim = 0
        self.hapNum = int(4**self.dim)
        if len(bRate) != self.hapNum or len(dRate) != self.hapNum or len(sRate) != self.hapNum:
            print_err("BirthDeathModel.SetRates() fatal error: inconsistent dimension")
            sys.exit(1)
        for el in mRate:
            if len(el) != self.dim:
                print_err("BirthDeathModel.SetRates() fatal error: inconsistent dimension")
                sys.exit(1)
        self.InitLiveBranches()
        self.bRate = bRate
        self.dRate = dRate
        self.sRate = sRate
        self.mRate = mRate
        self.migPopRate = [sum(mr) for mr in self.populationModel.migrationRates]

        self.bPopHaplotypeRate = [ [0]*self.hapNum for _ in range( self.popNum ) ]
        self.bPopHaplotypeRate[0][0] = self.BirthRate(0, 0)

        self.tPopHaplotypeRate = [ [0]*self.hapNum for p in self.populationModel.populations ]
        self.tPopRate = [0]*self.popNum
        self.totalPopHaplotypeRate(0, 0)
        self.tPopRate[0] = self.tPopHaplotypeRate[0][0]
        self.totalRate = self.tPopRate[0]


    def totalPopHaplotypeRate(self, popId, haplotype):
        r = self.bPopHaplotypeRate[popId][haplotype] + self.dRate[haplotype] + self.sRate[haplotype] + self.migPopRate[popId] + sum( self.mRate[haplotype] )
        self.tPopHaplotypeRate[popId][haplotype] = r*len(self.liveBranches[popId][haplotype])

    def BirthRate(self, popId, haplotype):
        pop = self.populationModel.populations[popId]
        return self.bRate[haplotype]*pop.susceptible/pop.size*pop.contactDensity

    def InitLiveBranches(self):
        self.liveBranches = []
        for i in range( len( self.populationModel.populations) ):
            self.liveBranches.append( [[] for _ in range(self.hapNum)] )
        self.liveBranches[0][0] += [1,2]
        self.lbCounter = 2 #live branch counter
        self.populationModel.populations[0].infectious = 2

    def GenerateEvent(self, useNumpy = False):
        if useNumpy:#TODO
            haplotype = np.random.choice( range(self.hapNum), p = [el/self.totalRate for el in self.tRate] )
            tRate =  self.tRate[haplotype] / len( self.liveBranches[haplotype] )
            prbs = [ self.bRate[haplotype]/tRate, self.dRate[haplotype]/tRate, self.sRate[haplotype]/tRate ] + [ el/tRate for el in self.mRate[haplotype] ]
            eventType = np.random.choice( range(3+self.dim), p = prbs )
            affectedBranch = np.random.randint(len(self.liveBranches[haplotype]))
        else:
            popId = fastChoose( range(len(self.populationModel.populations)), self.tPopRate, self.totalRate)
            haplotype = fastChoose(range(self.hapNum), self.tPopHaplotypeRate[popId], self.tPopRate[popId])

            tRate =  self.tPopHaplotypeRate[popId][haplotype] / len( self.liveBranches[popId][haplotype] )
            prbs = [ self.bPopHaplotypeRate[popId][haplotype], self.dRate[haplotype], self.sRate[haplotype], self.migPopRate[popId] ] + self.mRate[haplotype]
            eventType = fastChoose( range(4+self.dim), prbs , tRate)
            affectedBranch = fastRandint(len(self.liveBranches[popId][haplotype]))

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
        self.events.append(eventType)#TODO - clean me!!!
        #print("ET=", eventType, "   ", self.liveBranches)

    def Migration(self, popId, haplotype, affectedBranch):
        targetPopId = fastChoose(range(self.popNum-1), self.populationModel.migrationRates[popId], self.migPopRate[popId])
        if targetPopId >= popId:
            targetPopId += 1
        self.liveBranches[targetPopId][haplotype].append( self.liveBranches[popId][haplotype][affectedBranch] )

        self.liveBranches[popId][haplotype][affectedBranch] = self.liveBranches[popId][haplotype][-1]
        self.liveBranches[popId][haplotype].pop()

        self.totalPopHaplotypeRate(popId, haplotype)
        self.totalPopHaplotypeRate(targetPopId, haplotype)
        self.tPopRate[popId] = sum(self.tPopHaplotypeRate[popId])
        self.tPopRate[targetPopId] = sum(self.tPopHaplotypeRate[targetPopId])
        self.totalRate = sum( self.tPopRate )

    def Mutation(self, popId, haplotype, affectedBranch, mutationType):
        digit4 = 4**mutationType
        AS = floor(haplotype/digit4)%4
        DS = np.random.choice(range(3))#TODO non-uniform rates???
        if DS >= AS:
            DS += 1
        self.mutations.append(Mutation(self.liveBranches[popId][haplotype][affectedBranch], self.currentTime, AS, DS))
        newHaplotype = haplotype + (DS-AS)*digit4

        self.liveBranches[popId][newHaplotype].append( self.liveBranches[popId][haplotype][affectedBranch] )

        self.liveBranches[popId][haplotype][affectedBranch] = self.liveBranches[popId][haplotype][-1]
        self.liveBranches[popId][haplotype].pop()

        self.totalPopHaplotypeRate(popId, haplotype)
        self.totalPopHaplotypeRate(popId, newHaplotype)
        self.tPopRate[popId] = sum(self.tPopHaplotypeRate[popId])
        self.totalRate = sum( self.tPopRate )

    def SampleTime(self):
        tau = np.random.exponential(1.0/self.totalRate)
        self.currentTime += tau

    def Birth(self, popId, haplotype, affectedBranch):
        parentId = self.liveBranches[popId][haplotype][affectedBranch]
        self.times[parentId] = self.currentTime
        self.liveBranches[popId][haplotype].append(len(self.Tree))
        self.liveBranches[popId][haplotype][affectedBranch] = len(self.Tree) + 1
        for j in range(2):
            self.Tree.append(parentId)
            self.times.append(0)
            self.nodeSampling.append( NodeS() )
        self.populationModel.populations[popId].susceptible -= 1
        self.populationModel.populations[popId].infectious += 1

        for h in range(self.hapNum):
            self.bPopHaplotypeRate[popId][h] = self.BirthRate(popId, h)
            self.totalPopHaplotypeRate(popId, h)
        self.tPopRate[popId] = sum( self.tPopHaplotypeRate[popId] )
        self.totalRate = sum( self.tPopRate )

        self.lbCounter += 1

    def Death(self, popId, haplotype, affectedBranch):
        self.times[ self.liveBranches[popId][haplotype][affectedBranch] ] = self.currentTime
        self.liveBranches[popId][haplotype][affectedBranch] = self.liveBranches[popId][haplotype][-1]
        self.liveBranches[popId][haplotype].pop()

        self.populationModel.populations[popId].infectious -= 1

        self.totalPopHaplotypeRate(popId, haplotype)
        self.tPopRate[popId] = sum(self.tPopHaplotypeRate[popId])
        self.totalRate = sum( self.tPopRate )

        self.lbCounter -= 1
        self.susceptible -= 1

    def Sampling(self, popId, haplotype, affectedBranch):
        self.nodeSampling[ self.liveBranches[popId][haplotype][affectedBranch] ].state = -1
        self.Death(popId, haplotype, affectedBranch)
        self.sCounter += 1

#    def UpdateRate(self):
#        self.totalRate = self.B_rate[0] + self.D_rate[0] + self.S_rate[0] #TODO

    def SimulatePopulation(self, iterations):
        max_time = 0
        sCounter = 0
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
            sys.exit(1)
        for popLB in self.liveBranches:
            for lb in popLB:
                for br in lb:
                    self.times[br] = self.currentTime

    def GetGenealogy(self):
        for i in range(len(self.Tree) - 1, 0, -1):
            if self.nodeSampling[i].state == -1:
                parent = self.Tree[i]
                self.nodeSampling[i].genealogyIndex = len(self.genealogy)
                self.genealogy.append(-1)
                while parent != -1 and self.nodeSampling[parent].state != 2:
                    self.nodeSampling[parent].state += 1
                    if self.nodeSampling[parent].state == 2:
                        self.nodeSampling[parent].genealogyIndex = len(self.genealogy)
                        self.genealogy.append(-1)
                        break
                    parent = self.Tree[parent]
        if self.debug:
            print(self.nodeSampling)

        self.genealogyTimes = [0]*len(self.genealogy)
        for i in range(len(self.Tree) - 1, 0, -1):
            if self.nodeSampling[i].state == -1 or self.nodeSampling[i].state == 2:
                parent = self.Tree[i]
                while parent != -1 and self.nodeSampling[parent].state != 2:
                    parent = self.Tree[parent]
                self.genealogy[ self.nodeSampling[i].genealogyIndex ] = self.nodeSampling[parent].genealogyIndex
                self.genealogyTimes[ self.nodeSampling[i].genealogyIndex ] = self.times[i]

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
        print("Tree size: ", len(self.Tree))
        print("Number of sampled elements: ", self.sCounter)
