import numpy
import random
import sys

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
        self.nodeSampling = [[0], [0], [0]]
        self.times = [0, 0, 0]
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
        for i in range(len(self.Tree) - 1, -1, -1):
            if self.nodeSampling[i][0] == -1:
                parent = self.Tree[i]
                while parent != -1:
                    if self.nodeSampling[parent][0] == 2:
                        break
                    elif self.nodeSampling[parent][0] == 1:
                        self.nodeSampling[parent] = [2]
                        break
                    elif self.nodeSampling[parent][0] == 0:
                        self.nodeSampling[parent] = [1]
                    parent = self.Tree[parent]
        if self.debug:
            print(self.nodeSampling)
        for i in range(len(self.Tree)):
            if self.nodeSampling[i][0] == 2 or self.nodeSampling[i][0] == -1:
                self.nodeSampling[i].append(-1)
                break
        if self.debug:
            print(self.nodeSampling)
        for i in range(len(self.Tree) - 1, 0, -1):
            if self.nodeSampling[i][0] == -1 or self.nodeSampling[i][0] == 2:
                child = i
                parent = self.Tree[child]
                while self.nodeSampling[parent][0] != 2:
                    parent = self.Tree[parent]
                self.nodeSampling[child].append(parent)
        if self.debug:
                print(self.nodeSampling)
        newPos = 0
        for i in range(len(self.Tree)):
            if len(self.nodeSampling[i]) == 2:
                self.nodeSampling[i].append(newPos)
                newPos += 1
        if self.debug:
            print(self.nodeSampling)
        self.genealogy = [-1] * newPos
        for i in range(len(self.Tree)):
            if len(self.nodeSampling[i]) == 3 and self.nodeSampling[i][1] != -1:
                self.genealogy[self.nodeSampling[i][2]] = self.nodeSampling[self.nodeSampling[i][1]][2]

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
            self.nodeSampling.append([0])

    def Death(self, affectedBranch):
        self.liveBranches[affectedBranch] = self.liveBranches[-1]
        self.liveBranches.pop()

    def Sampling(self, affectedBranch):
        self.nodeSampling[affectedBranch] = [-1]
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
            #self.SampleTime(affectedBranch)
            if eventType < self.B_rate[0]/self.totalRate:
                self.Birth(affectedBranch)
            elif eventType < self.B_rate[0]/self.totalRate + self.D_rate[0]/self.totalRate:
                self.Death(affectedBranch)
            else:
                self.Sampling(affectedBranch)
                sCounter += 1
            #max_time = self.update_time(liveBranches, ActiveBranch, max_time)
            if len(self.liveBranches) == 0:
                break
        if self.debug:
            print(self.Tree)
            print(self.nodeSampling)
        if sCounter < 2: #TODO if number of sampled leaves is 0 (probably 1 as well), then GetGenealogy seems to go to an infinite cycle
            sys.exit(1)
