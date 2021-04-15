#!/usr/bin/env python3

import argparse
import sys
import time
from BirthDeathCython import BirthDeathModel, PopulationModel, Population, Lockdown
from IO import ReadRates, ReadPopulations, ReadMigrationRates, ReadSusceptibility, ReadSusceptibilityTransition, writeGenomeNewick, writeMutations
from random import randrange
import numpy as np

parser = argparse.ArgumentParser(description='Migration inference from PSMC.')

parser.add_argument('frate',
                    help='file with rates')


parser.add_argument('--iterations', '-it', nargs=1, type=int, default=1000,
                    help='number of iterations (default is 1000)')
parser.add_argument('--populationModel', '-pm', nargs=2, default=None,
                    help='population model: a file with population sizes etc, and a file with migration rate matrix')
parser.add_argument('--susceptibility', '-su', nargs=1, default=None,
                    help='susceptibility file')
parser.add_argument('--suscepTransition', '-st', nargs=1, default=None,
                    help='susceptibility transition file')

# parser.add_argument('--lockdownModel', '-ld', nargs=1, default=None,
#                     help='lockdown model: a file with parameters for lockdowns')
parser.add_argument('--seed', '-seed', nargs=1, type=int, default=None,
                    help='random seed')
parser.add_argument("--createNewick", '-nwk',
                    help="Create a newick file of tree *.nwk ",
                    action="store_true")
parser.add_argument("--writeMutations", '-tsv',
                    help="Create a mutation file *.tsv ",
                    action="store_true")

clargs = parser.parse_args()

if isinstance(clargs.frate, list):
    clargs.frate = clargs.frate[0]
if isinstance(clargs.iterations, list):
    clargs.iterations = clargs.iterations[0]
if isinstance(clargs.susceptibility, list):
    clargs.susceptibility = clargs.susceptibility[0]
if isinstance(clargs.suscepTransition, list):
    clargs.suscepTransition = clargs.suscepTransition[0]
if isinstance(clargs.seed, list):
    clargs.seed = clargs.seed[0]

bRate, dRate, sRate, mRate = ReadRates(clargs.frate)

if clargs.populationModel == None:
    popModel = None
    lockdownModel = None
else:
    populations, lockdownModel = ReadPopulations(clargs.populationModel[0])
    migrationRates = ReadMigrationRates(clargs.populationModel[1])
    popModel = [populations, migrationRates]

if clargs.susceptibility == None:
    susceptible = None
else:
    susceptible = ReadSusceptibility(clargs.susceptibility)

if clargs.suscepTransition == None:
    suscepTransition = None
else:
    suscepTransition = ReadSusceptibilityTransition(clargs.suscepTransition)

if clargs.seed == None:
    rndseed = randrange(sys.maxsize)
else:
    rndseed = clargs.seed
print("Seed: ", rndseed)

simulation = BirthDeathModel(clargs.iterations, bRate, dRate, sRate, mRate, populationModel=popModel, susceptible=susceptible, lockdownModel=lockdownModel, suscepTransition=suscepTransition, rndseed=rndseed)
# simulation.Debug()
t1 = time.time()
simulation.SimulatePopulation(clargs.iterations)
#simulation.Debug()
t2 = time.time()
simulation.GetGenealogy()
# simulation.Debug()
t3 = time.time()
simulation.Report()
print(t2 - t1)
print(t3 - t2)
print("_________________________________")

pruferSeq, times, mut = simulation.Output_tree_mutations()


if clargs.createNewick:
    writeGenomeNewick(pruferSeq, times)
if clargs.writeMutations:
    writeMutations(mut)
