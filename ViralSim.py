#!/usr/bin/env python3

import argparse
import sys
import time
from BirthDeath import BirthDeathModel, PopulationModel, Population
from IO import ReadRates, ReadPopulations, ReadMigrationRates

parser = argparse.ArgumentParser(description='Migration inference from PSMC.')

parser.add_argument('frate',
                    help='file with rates')


parser.add_argument('--iterations', '-it', nargs=1, type=int, default=1000,
                    help='number of iternations (default is 1000)')
parser.add_argument('--populationModel', '-pm', nargs=2, default=None,
                    help='number of iternations (default is 1000)')
parser.add_argument('--debug', action='store_true',
                    help='Debug mode, more input enabled')

clargs = parser.parse_args()

if isinstance(clargs.frate, list):
    clargs.frate = clargs.frate[0]
if isinstance(clargs.iterations, list):
    clargs.iterations = clargs.iterations[0]
if isinstance(clargs.debug, list):
    clargs.debug = clargs.debug[0]

bRate, dRate, sRate, mRate = ReadRates(clargs.frate)

if clargs.populationModel == None:
    populationModel = PopulationModel([Population()], [[]])
else:
    populations = ReadPopulations(clargs.populationModel[0])
    migrationRates = ReadMigrationRates(clargs.populationModel[1])
    populationModel = PopulationModel(populations, migrationRates)

simulation = BirthDeathModel(bRate, dRate, sRate, mRate, debug = clargs.debug, populationModel= populationModel)
t1 = time.time()
simulation.SimulatePopulation(clargs.iterations)
t2 = time.time()
simulation.GetGenealogy()
t3 = time.time()
simulation.Report()
print(t2 - t1)
print(t3 - t2)
# print(tree1.Tree)
# print(tree1.newTree)
# print(tree1.nodeSampling)
# print(tree1.times)
# print(tree1.newTimes)
