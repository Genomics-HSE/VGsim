#!/usr/bin/env python3

import argparse
import sys
import time
from BirthDeath import BirthDeathModel

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

parser = argparse.ArgumentParser(description='Migration inference from PSMC.')

parser.add_argument('--iterations', '-it', nargs=1, type=int, default=1000,
                    help='number of iternations (default is 1000)')
parser.add_argument('--debug', action='store_true',
                    help='Debug mode, more input enabled')

clargs = parser.parse_args()

if isinstance(clargs.iterations, list):
    clargs.iterations = clargs.iterations[0]
if isinstance(clargs.debug, list):
    clargs.debug = clargs.debug[0]

##Начало программы
B_rate_data = [25]
D_rate_data = [9]
S_rate_data = [1]
#getData(B_rate_data, D_rate_data, S_rate_data)
tree1 = BirthDeathModel(B_rate_data, D_rate_data, S_rate_data, debug = clargs.debug)
t1 = time.time()
tree1.SimulatePopulation(clargs.iterations)
t2 = time.time()
tree1.GetGenealogy()
t3 = time.time()
print(t2 - t1)
print(t3 - t2)
# print(tree1.Tree)
# print(tree1.newTree)
# print(tree1.nodeSampling)
# print(tree1.times)
# print(tree1.newTimes)
