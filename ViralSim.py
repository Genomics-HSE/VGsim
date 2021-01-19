import random
import sys
import time


class birth_death_model:
    def __init__(self, _B_rate, _D_rate, _S_rate, _M_rate = 0):
        self.B_rate = _B_rate
        self.D_rate = _D_rate
        self.S_rate = _S_rate
        self.Tree = [-1]
        self.newTree = []
        self.NodeSampling = [[0]]
        self.Time = [0]
        self.newTime = []
        #self.TypeOfMutations = [[0] * len(_M_rate[0])]
        
    def update_time(self, liveBranches, ActiveBranch, max_time):
        #Проверить параметр в expo
        tau = random.expovariate(1/len(liveBranches)) ## А точно ли тут R_rate, если что упростить
        self.Time[ActiveBranch] = max_time
        max_time += tau
        liveBranches[liveBranches.index(ActiveBranch)] = liveBranches[len(liveBranches) - 1]
        del liveBranches[-1]
        return max_time
        
#Может ли корень сразу получить мутацию или нет
#Возможны ли поэтапные мутации
#
    def create_new_tree(self):
        for i in range(len(self.Tree) - 1, -1, -1):
            if self.NodeSampling[i][0] == -1:
                parent = self.Tree[i]
                while parent != -1:
                    if self.NodeSampling[parent][0] == 2:
                        break
                    elif self.NodeSampling[parent][0] == 1:
                        self.NodeSampling[parent] = [2]
                        break
                    elif self.NodeSampling[parent][0] == 0:
                        self.NodeSampling[parent] = [1]
                    parent = self.Tree[parent]
        print(self.NodeSampling)
        for i in range(len(tree1.Tree)):
            if tree1.NodeSampling[i][0] == 2 or tree1.NodeSampling[i][0] == -1:
                tree1.NodeSampling[i].append(-1)
                break
        print(self.NodeSampling)
        for i in range(len(tree1.Tree) - 1, 0, -1):
            if tree1.NodeSampling[i][0] == -1 or tree1.NodeSampling[i][0] == 2:
                child = i
                parent = tree1.Tree[child]
                while tree1.NodeSampling[parent][0] != 2:
                    parent = tree1.Tree[parent]
                tree1.NodeSampling[child].append(parent)
        print(self.NodeSampling)
        newPos = 0
        for i in range(len(tree1.Tree)):
            if len(tree1.NodeSampling[i]) == 2:
                tree1.NodeSampling[i].append(newPos)
                newPos += 1
        print(self.NodeSampling)
        tree1.newTree = [-1] * newPos
        for i in range(len(tree1.Tree)):
            if len(tree1.NodeSampling[i]) == 3 and tree1.NodeSampling[i][1] != -1:
                tree1.newTree[tree1.NodeSampling[i][2]] = tree1.NodeSampling[tree1.NodeSampling[i][1]][2]
        
    def create_tree(self, iterations):
        liveBranches = [0]
        max_time = 0
        for j in range(0, iterations):
            if len(liveBranches) > 0:
                ProbabilityBDS = random.random()
                ActiveBranch = random.choice(liveBranches)
                if 0 <= ProbabilityBDS <= self.B_rate[0]:
                    # Потом делаем новые вершины
                    for j in range(0,2):
                        liveBranches.append(len(self.Tree))
                        self.Tree.append(ActiveBranch)
                        self.Time.append(0)
                        self.NodeSampling.append([0])
                elif self.B_rate[0] + self.D_rate[0] <= ProbabilityBDS:
                    self.NodeSampling[ActiveBranch] = [-1]
                #max_time = self.update_time(liveBranches, ActiveBranch, max_time)
            else:
                break
        self.create_new_tree()

def getData(B_rate_data, D_rate_data, S_rate_data):
    data = open('data.txt', 'r')
    Number_of_Strings = 0
    for line in data:
        Number_of_Strings += 1
        line = line.split()
        R = int(line[0]) + int(line[1]) + int(line[2])
        B_rate_data.append(int(line[0]) / R)
        D_rate_data.append(int(line[1]) / R)
        S_rate_data.append(int(line[2]) / R)
    data.close()
        

        
##Начало программы
B_rate_data = [25]
D_rate_data = [9]
S_rate_data = [1]
#getData(B_rate_data, D_rate_data, S_rate_data)
t1 = time.time()
tree1 = birth_death_model(B_rate_data, D_rate_data, S_rate_data)
tree1.create_tree(10000000)
t2 = time.time()
print(t2 - t1)
# print(tree1.Tree)
# print(tree1.newTree)
# print(tree1.NodeSampling)
# print(tree1.Time)
# print(tree1.newTime)
