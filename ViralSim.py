#!/usr/bin/env python3

import argparse
import sys
import time
from BirthDeathCython import BirthDeathModel, PopulationModel, Population
from IO import ReadRates, ReadPopulations, ReadMigrationRates, ReadSusceptibility

parser = argparse.ArgumentParser(description='Migration inference from PSMC.')

parser.add_argument('frate',
                    help='file with rates')


parser.add_argument('--iterations', '-it', nargs=1, type=int, default=1000,
                    help='number of iterations (default is 1000)')
parser.add_argument('--populationModel', '-pm', nargs=2, default=None,
                    help='population model: a file with population sizes etc, and a file with migration rate matrix')
parser.add_argument('--susceptibility', '-su', nargs=1, default=None,
                    help='susceptibility file')
parser.add_argument('--seed', '-seed', nargs=1, type=int, default=None,
                    help='random seed')

clargs = parser.parse_args()

if isinstance(clargs.frate, list):
    clargs.frate = clargs.frate[0]
if isinstance(clargs.iterations, list):
    clargs.iterations = clargs.iterations[0]
if isinstance(clargs.susceptibility, list):
    clargs.susceptibility = clargs.susceptibility[0]
if isinstance(clargs.seed, list):
    clargs.seed = clargs.seed[0]

bRate, dRate, sRate, mRate = ReadRates(clargs.frate)

if clargs.populationModel == None:
    popModel = None
else:
    populations = ReadPopulations(clargs.populationModel[0])
    migrationRates = ReadMigrationRates(clargs.populationModel[1])
    popModel = [populations, migrationRates]

if clargs.susceptibility == None:
    susceptible = None
else:
    susceptible = ReadSusceptibility(clargs.susceptibility)

if clargs.seed == None:
    rndseed = int(time.time())
else: 
    rndseed = clargs.seed
print(rndseed)

simulation = BirthDeathModel(clargs.iterations, bRate, dRate, sRate, mRate, populationModel=popModel, susceptible=susceptible, rndseed=rndseed)
# simulation.Debug()
# t1 = time.time()
simulation.SimulatePopulation(clargs.iterations)
# simulation.Debug()
# t2 = time.time()
simulation.GetGenealogy()
simulation.Debug()
# t3 = time.time()
# simulation.Report()
# print(t2 - t1)
# print(t3 - t2)
print("_________________________________")
# print(tree1.Tree)
# print(tree1.newTree)
# print(tree1.nodeSampling)
# print(tree1.times)
# print(tree1.newTimes)
