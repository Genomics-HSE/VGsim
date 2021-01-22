import numpy
import random
import sys

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
    def __init__(self, _B_rate, _D_rate, _S_rate, _M_rate = 0, **kwargs):
        self.B_rate = _B_rate
        self.D_rate = _D_rate
        self.S_rate = _S_rate
        self.Tree = [-1, 0, 0]
        self.genealogy = []
        self.nodeSampling = [NodeS() for _ in range(3)]
        self.times = [0]*3
        self.genealogyTimes = []
        self.liveBranches = [1,2]
        self.totalTime = 0.0
        self.UpdateRate()
        self.debug = False
        if "debug" in kwargs:
            if kwargs["debug"]:
                self.debug = True
        #self.TypeOfMutations = [[0] * len(_M_rate[0])]

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
        tau = numpy.random.exponential(self.totalRate) ## А точно ли тут R_rate, если что упростить
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
            affectedBranch = numpy.random.randint(len(self.liveBranches)) #random.choice(liveBranches)
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
