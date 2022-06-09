#!/usr/bin/env python3

import argparse
import sys
import time
import VGsim
from random import randrange
import numpy as np
import math

parser = argparse.ArgumentParser(description='Migration inference from PSMC.')

parser.add_argument('--iterations', '-it', nargs=1, type=int, default=1000,
                    help='number of iterations (default is 1000)')
parser.add_argument('--sampleSize', '-s', nargs=1, type=int, default=None,
                    help='number of sample (default is None)')
parser.add_argument('--time', '-t', nargs=1, type=float, default=None,
                    help='time for stopping simulation (default is None)')
parser.add_argument('--max_haplotypes_number', '-mhn', nargs=1, type=int, default=None,
                    help='#TODO')
parser.add_argument('--seed', '-seed', nargs=1, type=float, default=None,
                    help='random seed')

parser.add_argument('--rates', '-rt', nargs=1, default=None,
                    help='rate: a file with rates for each haplotype')
parser.add_argument('--populationModel', '-pm', nargs=2, default=None,
                    help='population model: a file with population sizes etc, and a file with migration rate matrix')
parser.add_argument('--susceptibility', '-su', nargs=1, default=None,
                    help='susceptibility file')
parser.add_argument('--suscepTransition', '-st', nargs=1, default=None,
                    help='susceptibility transition file')

parser.add_argument('--sampling_probability', help="#TODO", action="store_true")

parser.add_argument("--createNewick", '-nwk', nargs=1, default=False, 
                    help="Create a newick file of tree *.nwk ")
parser.add_argument("--writeMutations", '-tsv', nargs=1, default=False, 
                    help="Create a mutation file *.tsv ")
parser.add_argument("--writeMigrations", nargs=1, default=False, 
                    help="Create a migration file *.tsv ")
parser.add_argument("--output_chain_events", nargs=1, default=False,
                    help="#TODO")

parser.add_argument("-citation", '-c', help="Information for citation.")

clargs = parser.parse_args()

if clargs.citation != None:
    print("VGsim: scalable viral genealogy simulator for global pandemic")
    print("Vladimir Shchur, Vadim Spirin, Victor Pokrovskii, Evgeni Burovski, Nicola De Maio, Russell Corbett-Detig")
    print("medRxiv 2021.04.21.21255891; doi: https://doi.org/10.1101/2021.04.21.21255891")
    sys.exit(0)

if isinstance(clargs.rates, list):
    clargs.rates = clargs.rates[0]
if isinstance(clargs.iterations, list):
    clargs.iterations = clargs.iterations[0]
if isinstance(clargs.sampleSize, list):
    clargs.sampleSize = clargs.sampleSize[0]
if isinstance(clargs.time, list):
	clargs.time = clargs.time[0]
if isinstance(clargs.susceptibility, list):
    clargs.susceptibility = clargs.susceptibility[0]
if isinstance(clargs.suscepTransition, list):
    clargs.suscepTransition = clargs.suscepTransition[0]
if isinstance(clargs.seed, list):
    clargs.seed = clargs.seed[0]
if isinstance(clargs.output_chain_events, list):
    clargs.output_chain_events = clargs.output_chain_events[0]

if clargs.rates == None:
    bRate, dRate, sRate, mRate = [2], [1], [0.1], [[]]
else:
    bRate, dRate, sRate, mRate = VGsim.IO.read_rates(clargs.rates)

if clargs.sampleSize == None:
    clargs.sampleSize = clargs.iterations
if clargs.time == None:
    clargs.time = -1

if clargs.populationModel == None:
    sizes, contactDensity, contactAfter, startLD, endLD, samplingMultiplier = [1000000], [1], [1], [1], [1], [1]
    migrationRates = [[0.0]]
else:
    sizes, contactDensity, contactAfter, startLD, endLD, samplingMultiplier = VGsim.IO.read_populations(clargs.populationModel[0])
    migrationRates = VGsim.IO.read_matrix(clargs.populationModel[1])

if clargs.susceptibility == None:
    susceptible = [[1.0] for _ in range(len(bRate))]
    susType = [0 for _ in range(len(bRate))]
else:
    susceptible, susType = VGsim.IO.read_susceptibility(clargs.susceptibility)

if clargs.suscepTransition == None:
    suscepTransition = [[0.0]]
else:
    suscepTransition = VGsim.IO.read_matrix(clargs.suscepTransition)

if clargs.seed == None:
    seed = randrange(sys.maxsize)
else:
	seed = clargs.seed

simulator = VGsim.Simulator(number_of_sites=int(math.log(len(bRate), 4)), populations_number=len(sizes), number_of_susceptible_groups=len(susceptible[0]), seed=int(seed), sampling_probability=clargs.sampling_probability)

for i in range(len(bRate)):
	simulator.set_transmission_rate(bRate[i], i)
	simulator.set_recovery_rate(dRate[i], i)
	simulator.set_sampling_rate(sRate[i], i)
	for j in range(len(mRate[0])):
		simulator.set_mutation_rate(mRate[i][j][0], [mRate[i][j][1], mRate[i][j][2], mRate[i][j][3], mRate[i][j][4]], i, j)

for i in range(len(sizes)):
    simulator.set_population_size(sizes[i], i)
    simulator.set_contact_density(contactDensity[i], i)
    simulator.set_npi([contactAfter[i], startLD[i], endLD[i]], i)
    simulator.set_sampling_multiplier(samplingMultiplier[i], i)
    for j in range(len(sizes)):
        if i != j:
            simulator.set_migration_probability(probability=migrationRates[i][j], source=i, target=j)
for i in range(len(susceptible)):
    for j in range(len(susceptible[i])):
        simulator.set_susceptibility(float(susceptible[i][j]), i, j)

for i in range(len(susType)):
    simulator.set_susceptibility_type(susType[i], i)

for i in range(len(suscepTransition)):
    for j in range(len(suscepTransition[i])):
        if i != j:
            simulator.set_immunity_transition(suscepTransition[i][j], i, j)

simulator.simulate(clargs.iterations, clargs.sampleSize, clargs.time)
simulator.genealogy(int(seed))

if clargs.createNewick:
    simulator.output_newick(clargs.createNewick)
if clargs.writeMutations:
    simulator.output_mutations(clargs.writeMutations)
if clargs.writeMigrations:
    simulator.output_migrations(clargs.writeMigrations)
if clargs.output_chain_events:
    simulator.output_chain_events(clargs.output_chain_events)