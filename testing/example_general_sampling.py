import VGsim
import sys
sys.path.clear()
sys.path.insert(0, '/home/lev/PycharmProjects')
print(sys.path)
#import os
#import os.path
#from src import _interface

#from src._interface import Simulator

number_of_sites = 2
populations_number = 3
number_of_susceptible_groups = 3
simulator = VGsim.Simulator(number_of_sites, populations_number, number_of_susceptible_groups, seed=1234)

#Set epidemiological parameters
simulator.set_transmission_rate(0.25)
simulator.set_transmission_rate(0.5, haplotype="GG")
simulator.set_recovery_rate(0.099)
simulator.set_sampling_rate(0.001)

#Set mutation rates
mutation_rate=0.00003
substitution_weights=[1,1,1,2]#ATCG
simulator.set_mutation_rate(mutation_rate, substitution_weights)
simulator.set_mutation_rate(3*mutation_rate, haplotype='G*', mutation=1)

#Set host immunity types triggered by infection
simulator.set_susceptibility_type(1)
simulator.set_susceptibility_type(2, haplotype='G*')
simulator.set_susceptibility(0.1, susceptibility_type=1)
simulator.set_susceptibility(0.5, susceptibility_type=1, haplotype='G*')
simulator.set_susceptibility(0.0, susceptibility_type=2)

#Set loss of immunity
simulator.set_immunity_transition(1/90, source=1, target=0)
simulator.set_immunity_transition(1/180, source=2, target=0)

#Set host population structure with three populations and migration
simulator.set_population_size(10000000, population=0)
simulator.set_population_size(5000000, population=1)
simulator.set_population_size(1000000, population=2)
simulator.set_migration_probability(10/365/2)

#Set specific sampling efforts in different populations
simulator.set_sampling_multiplier(3, population=1)
simulator.set_sampling_multiplier(0, population=2)

simulator.set_general_sampling(0.1, 3, [20, 10, 30])

#simulator.set_super_spread_rate(0.1, 10, 20, 0)

#Set non-pharmasutical interventions (same for each population)
simulator.set_npi([0.1, 0.01, 0.002])

#Run simulation with the exact algorithm
simulator.simulate(10000000, epidemic_time=110)

simulator.genealogy()
