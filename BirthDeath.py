import numpy as np
import random
import sys

def print_err(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class NodeS: #C-like structure
    def __init__(self):
        self.state = 0
        self.genealogyIndex = -1

class NeutralMutations:
    def __init__(self):
        self.muRate = 0.0

    def muRate(self):
        return self.muRate

class BirthDeathModel:
    def __init__(self, bRate, dRate, sRate, mRate = 0, **kwargs):

        self.Tree = [-1, 0, 0]
        self.genealogy = []
        self.nodeSampling = [NodeS() for _ in range(3)]
        self.times = [0]*3
        self.genealogyTimes = []
        #self.liveBranches = [1,2]
        self.totalTime = 0.0
        self.UpdateRate()
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
        self.tRate = [0]*self.hapNum

    def CalculateRate(self):
        self.totalRate = 0
        for i in range(hapNum):
            tRate = self.bRate[i] + self.dRate[i] + self.sRate[i]
            if self.dim > 0:
                for el in self.mRate[i]:
                    tRate += el
            self.tRate[i] = tRate * len( self.liveBranches[i] )
            self.totalRate += self.tRate[i]

    def GenerateEvent(self):
        haplotype = np.random.choice( range(self.hapNum), p = [el/self.totalRate for el in self.tRate] )
        tRate =  self.tRate[haplotype] / len( self.liveBranches[i] )
        prbs = [ self.bRate[haplotype]/tRate, self.dRate[haplotype]/tRate, self.sRate[haplotype]/tRate ] + [ el/tRate for el in self.mRate[haplotype] ]
        eventType = np.random.choice( range(self.hapNum), p = prbs )
        affectedBranch = np.random.randint(len(self.liveBranches[haplotype]))

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

    def SampleTime(self, affectedBranch):
        #Проверить параметр в expo
        tau = np.random.exponential(self.totalRate)
        self.totalTime += tau
        self.times[affectedBranch] = self.totalTime

    def Birth(self, affectedBranch):
        self.liveBranches.append(len(self.Tree))
        self.liveBranches[affectedBranch] = len(self.Tree) + 1
        for j in range(2):
            self.Tree.append(affectedBranch)
            self.times.append(0)
            self.nodeSampling.append( NodeS() )

    def Death(self, affectedBranch):
        self.liveBranches[affectedBranch] = self.liveBranches[-1]
        self.liveBranches.pop()

    def Sampling(self, affectedBranch):
        self.nodeSampling[affectedBranch].state = -1
        self.Death(affectedBranch)

    def UpdateRate(self):
        self.totalRate = self.B_rate[0] + self.D_rate[0] + self.S_rate[0] #TODO

    def SimulatePopulation(self, iterations):
        max_time = 0
        sCounter = 0
        for j in range(0, iterations):
            self.UpdateRate()
            eventType = random.random()
            affectedBranch = np.random.randint(len(self.liveBranches)) #random.choice(liveBranches)
            self.SampleTime(affectedBranch)
            if eventType < self.B_rate[0]/self.totalRate:
                self.Birth(affectedBranch)
            elif eventType < self.B_rate[0]/self.totalRate + self.D_rate[0]/self.totalRate:
                self.Death(affectedBranch)
            else:
                self.Sampling(affectedBranch)
                sCounter += 1
            if len(self.liveBranches) == 0:
                break
        if self.debug:
            print(self.Tree)
            print(self.nodeSampling)
        if sCounter < 2: #TODO if number of sampled leaves is 0 (probably 1 as well), then GetGenealogy seems to go to an infinite cycle
            sys.exit(1)
