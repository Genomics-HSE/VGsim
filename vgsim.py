# #!/usr/bin/env python3

# import argparse
# import sys
# import time
# from VGsim import BirthDeathModel, PopulationModel, Population, Lockdown
# from VGsim.IO import ReadRates, ReadPopulations, ReadMigrationRates, ReadSusceptibility, ReadSusceptibilityTransition, writeGenomeNewick, writeMutations
# from random import randrange
# import numpy as np

# parser = argparse.ArgumentParser(description='Migration inference from PSMC.')

# parser.add_argument('frate',
#                     help='file with rates')


# parser.add_argument('--iterations', '-it', nargs=1, type=int, default=1000,
#                     help='number of iterations (default is 1000)')
# parser.add_argument('--sampleSize', '-s', nargs=1, type=int, default=None,
#                     help='number of sample (default is None)')
# parser.add_argument('--time', '-t', nargs=1, type=float, default=None,
#                     help='time for stopping simulation (default is None)')
# parser.add_argument('--populationModel', '-pm', nargs=2, default=None,
#                     help='population model: a file with population sizes etc, and a file with migration rate matrix')
# parser.add_argument('--susceptibility', '-su', nargs=1, default=None,
#                     help='susceptibility file')
# parser.add_argument('--suscepTransition', '-st', nargs=1, default=None,
#                     help='susceptibility transition file')

# parser.add_argument('--seed', '-seed', nargs=1, type=int, default=None,
#                     help='random seed')
# parser.add_argument("--createNewick", '-nwk',
#                     help="Create a newick file of tree *.nwk ",
#                     action="store_true")
# parser.add_argument("--writeMutations", '-tsv',
#                     help="Create a mutation file *.tsv ",
#                     action="store_true")
# parser.add_argument("--writeMigrations",
#                     help="Create a migration file *.txt ",
#                     action="store_true")

# parser.add_argument("-citation", '-c', help="Information for citation.")

# clargs = parser.parse_args()

# if clargs.citation != None:
#     print("VGsim: scalable viral genealogy simulator for global pandemic")
#     print("Vladimir Shchur, Vadim Spirin, Victor Pokrovskii, Evgeni Burovski, Nicola De Maio, Russell Corbett-Detig")
#     print("medRxiv 2021.04.21.21255891; doi: https://doi.org/10.1101/2021.04.21.21255891")
#     sys.exit(0)

# if isinstance(clargs.frate, list):
#     clargs.frate = clargs.frate[0]
# if isinstance(clargs.iterations, list):
#     clargs.iterations = clargs.iterations[0]
# if isinstance(clargs.sampleSize, list):
#     clargs.sampleSize = clargs.sampleSize[0]
# if isinstance(clargs.susceptibility, list):
#     clargs.susceptibility = clargs.susceptibility[0]
# if isinstance(clargs.suscepTransition, list):
#     clargs.suscepTransition = clargs.suscepTransition[0]
# if isinstance(clargs.seed, list):
#     clargs.seed = clargs.seed[0]

# bRate, dRate, sRate, mRate = ReadRates(clargs.frate)

# if clargs.sampleSize == None:
#     clargs.sampleSize = clargs.iterations

# if clargs.populationModel == None:
#     popModel = None
#     lockdownModel = None
# else:
#     populations, lockdownModel, samplingMulti = ReadPopulations(clargs.populationModel[0])
#     migrationRates = ReadMigrationRates(clargs.populationModel[1])
#     popModel = [populations, migrationRates]

# if clargs.susceptibility == None:
#     susceptible = None
# else:
#     susceptible = ReadSusceptibility(clargs.susceptibility)

# if clargs.suscepTransition == None:
#     suscepTransition = None
# else:
#     suscepTransition = ReadSusceptibilityTransition(clargs.suscepTransition)

# if clargs.seed == None:
#     rndseed = randrange(sys.maxsize)
# else:
#     rndseed = clargs.seed
# print("Seed: ", rndseed)

# simulation = BirthDeathModel(bRate, dRate, sRate, mRate, populationModel=popModel, susceptible=susceptible, suscepTransition=suscepTransition, lockdownModel=lockdownModel, samplingMultiplier=samplingMulti, rndseed=rndseed)
# # simulation.Debug()
# # t1 = time.time()
# simulation.SimulatePopulation(clargs.iterations, clargs.sampleSize)
# # simulation.Debug()
# # t2 = time.time()
# simulation.GetGenealogy(rndseed)
# # simulation.Debug()
# # t3 = time.time()
# # simulation.Report()
# # print(t2 - t1)
# # print(t3 - t2)
# print("_________________________________")

# if clargs.createNewick or clargs.writeMutations:
#     pruferSeq, times, mut, populations = simulation.Output_tree_mutations()

# if clargs.createNewick:
#     writeGenomeNewick(pruferSeq, times, populations)
# if clargs.writeMutations:
#     writeMutations(mut, len(pruferSeq))
# if clargs.writeMigrations:
#     simulation.writeMigrations()

import VGsim
simulator = VGsim.Simulator(2, 3, 3, 1234)
simulator.set_transmission_rate(40)
simulator.set_recovery_rate(15)
simulator.set_sampling_rate(4)
simulator.set_migration_probability(0.01)
simulator.set_lockdown([0.5, 0.2, 0.05])
simulator.set_sampling_multiplier(1.8)
simulator.set_susceptibility_type(2)
simulator.set_susceptibility(0.5, susceptibility_type=1)
simulator.set_immunity_transition(0.00001)
simulator.simulate(10000)
simulator.genealogy(1234)
