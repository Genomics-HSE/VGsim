#!/usr/bin/env python3

import argparse
import sys
import time
from BirthDeath import BirthDeathModel
from IO import ReadRates

parser = argparse.ArgumentParser(description='Migration inference from PSMC.')

parser.add_argument('frate',
                    help='file with rates')

parser.add_argument('--iterations', '-it', nargs=1, type=int, default=1000,
                    help='number of iternations (default is 1000)')
parser.add_argument('--populationSize', '-ps', nargs=1, type=int, default=1000000,
                    help='population size (number of individuals)')
parser.add_argument('--debug', action='store_true',
                    help='Debug mode, more input enabled')

clargs = parser.parse_args()

if isinstance(clargs.frate, list):
    clargs.frate = clargs.frate[0]
if isinstance(clargs.iterations, list):
    clargs.iterations = clargs.iterations[0]
if isinstance(clargs.populationSize, list):
    clargs.populationSize = clargs.populationSize[0]
if isinstance(clargs.debug, list):
    clargs.debug = clargs.debug[0]

bRate, dRate, sRate, mRate = ReadRates(clargs.frate)

B_rate_data = [25]
D_rate_data = [9]
S_rate_data = [1]

tree1 = BirthDeathModel(bRate, dRate, sRate, mRate, debug = clargs.debug, populationSize = clargs.populationSize)
t1 = time.time()
tree1.SimulatePopulation(clargs.iterations)
t2 = time.time()
tree1.GetGenealogy()
t3 = time.time()
tree1.Report()
print(t2 - t1)
print(t3 - t2)
# print(tree1.Tree)
# print(tree1.newTree)
# print(tree1.nodeSampling)
# print(tree1.times)
# print(tree1.newTimes)
