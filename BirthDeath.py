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
    while total < rn:
        i += 1
        total += w[i]
    return a[i]

class BirthDeathModel:
    def __init__(self, bRate, dRate, sRate, mRate = [], **kwargs):
        self.Tree = [-1, 0, 0]
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
            if len(bRate) != self.hapNum:
                print_err("BirthDeathModel.SetRates() fatal error: inconsistent dimension")
                sys.exit(1)
        self.bRate = bRate
        self.dRate = dRate
        self.sRate = sRate
        self.mRate = mRate
        self.liveBranches = [[] for _ in range(self.hapNum)]
        self.liveBranches[0] += [1,2]
        self.lbCounter = 2 #live branch counter
        self.tRate = [0]*self.hapNum

    def UpdateRate(self):
        self.totalRate = 0
        for i in range(self.hapNum):
            tRate = self.bRate[i] + self.dRate[i] + self.sRate[i]
            if self.dim > 0:
                for el in self.mRate[i]:
                    tRate += el
            self.tRate[i] = tRate * len( self.liveBranches[i] )
            self.totalRate += self.tRate[i]

    def GenerateEvent(self):
        #haplotype = np.random.choice( range(self.hapNum), p = [el/self.totalRate for el in self.tRate] )
        haplotype = fastChoose(range(self.hapNum), self.tRate, self.totalRate)
        tRate =  self.tRate[haplotype] / len( self.liveBranches[haplotype] )
        #prbs = [ self.bRate[haplotype]/tRate, self.dRate[haplotype]/tRate, self.sRate[haplotype]/tRate ] + [ el/tRate for el in self.mRate[haplotype] ]
        #eventType = np.random.choice( range(3+self.dim), p = prbs )
        prbs = [ self.bRate[haplotype], self.dRate[haplotype], self.sRate[haplotype] ] + self.mRate[haplotype]
        eventType = fastChoose( range(3+self.dim), prbs , tRate)
        affectedBranch = np.random.randint(len(self.liveBranches[haplotype]))
        #print("EVENT", eventType, haplotype, affectedBranch, len(self.liveBranches[haplotype]), sep = " ")
        if eventType == 0:
            self.Birth(haplotype, affectedBranch)
        elif eventType == 1:
            self.Death(haplotype, affectedBranch)
        elif eventType == 2:
            self.Sampling(haplotype, affectedBranch)
            self.sCounter += 1
        else:
            self.Mutation(haplotype, affectedBranch, eventType - 3)

    def Mutation(self, haplotype, affectedBranch, mutationType):
        digit4 = 4**mutationType
        AS = floor(haplotype/digit4)%4
        DS = np.random.choice(range(3))#TODO non-uniform rates???
        if DS >= AS:
            DS += 1
        self.mutations.append(Mutation(self.liveBranches[haplotype][affectedBranch], self.currentTime, AS, DS))
        newHaplotype = haplotype + (DS-AS)*digit4
        self.liveBranches[newHaplotype].append( self.liveBranches[haplotype][affectedBranch] )
        self.Death(haplotype, affectedBranch, False)

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

        for i in range(len(self.Tree) - 1, 0, -1):
            if self.nodeSampling[i].state == -1 or self.nodeSampling[i].state == 2:
                parent = self.Tree[i]
                while parent != -1 and self.nodeSampling[parent].state != 2:
                    parent = self.Tree[parent]
                self.genealogy[ self.nodeSampling[i].genealogyIndex ] = self.nodeSampling[parent].genealogyIndex
                self.genealogyTimes[ self.nodeSampling[i].genealogyIndex ] = self.times[i]

    def SampleTime(self):
        tau = np.random.exponential(self.totalRate)
        self.currentTime += tau

    def Birth(self, haplotype, affectedBranch):
        parentId = self.liveBranches[haplotype][affectedBranch]
        self.times[parentId] = self.currentTime
        self.liveBranches[haplotype].append(len(self.Tree))
        self.liveBranches[haplotype][affectedBranch] = len(self.Tree) + 1
        for j in range(2):
            self.Tree.append(parentId)
            self.times.append(0)
            self.nodeSampling.append( NodeS() )
        self.lbCounter += 1

    def Death(self, haplotype, affectedBranch, timeUpdate = True):
        if timeUpdate:
            self.times[ self.liveBranches[haplotype][affectedBranch] ] = self.currentTime
        self.liveBranches[haplotype][affectedBranch] = self.liveBranches[haplotype][-1]
        self.liveBranches[haplotype].pop()
        self.lbCounter -= 1

    def Sampling(self,haplotype, affectedBranch):
        self.nodeSampling[ self.liveBranches[haplotype][affectedBranch] ].state = -1
        self.Death(haplotype, affectedBranch)

#    def UpdateRate(self):
#        self.totalRate = self.B_rate[0] + self.D_rate[0] + self.S_rate[0] #TODO

    def SimulatePopulation(self, iterations):
        max_time = 0
        sCounter = 0
        for j in range(0, iterations):
            self.UpdateRate()
            self.SampleTime()
            self.GenerateEvent()
            if self.lbCounter == 0:
                break
        if self.debug:
            print(self.Tree)
            print(self.nodeSampling)
        if self.sCounter < 2: #TODO if number of sampled leaves is 0 (probably 1 as well), then GetGenealogy seems to go to an infinite cycle
            print("Bo-o-o-oring... Less than two cases were sampled.")
            sys.exit(1)
        for lb in self.liveBranches:
            for br in lb:
                self.times[br] = self.currentTime

    def Report(self):
        print("Number of lineages of each hyplotype: ", end = "")
        for el in self.liveBranches:
            print(len(el), " ", end="")
        print("")
        print("Tree size: ", len(self.Tree))
        print("Number of sampled elements: ", self.sCounter)
