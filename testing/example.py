import VGsim
import os

number_of_sites = 2
populations_number = 3
number_of_susceptible_groups = 3
simulator = VGsim.Simulator(number_of_sites, populations_number, number_of_susceptible_groups, seed=1234)

simulator.set_transmission_rate(0.25)
simulator.set_recovery_rate(0.099)
simulator.set_sampling_rate(0.001)
simulator.set_transmission_rate(0.5, haplotype="GG")
mutation_rate=0.000003
substitution_weights=[1,1,1,2]#ATCG
simulator.set_mutation_rate(mutation_rate, substitution_weights)
simulator.set_mutation_rate(3*mutation_rate, haplotype='G*', mutation=1)

simulator.set_susceptibility_type(1)
simulator.set_susceptibility_type(2, haplotype='G*')
simulator.set_susceptibility(0.1, susceptibility_type=1)
simulator.set_susceptibility(0.5, susceptibility_type=1, haplotype='G*')
simulator.set_immunity_transition(1/90, source=1, target=0)
simulator.set_susceptibility(0.0, susceptibility_type=2)
simulator.set_immunity_transition(1/180, source=2, target=0)

simulator.set_population_size(10000000, population=0)
simulator.set_population_size(5000000, population=1)
simulator.set_population_size(1000000, population=2)
simulator.set_sampling_multiplier(3, population=1)
simulator.set_sampling_multiplier(0, population=2)
simulator.set_npi([0.1, 0.01, 0.002])
simulator.set_migration_probability(10/365/2)

simulator.simulate(10000000, time=110)

simulator.set_immunity_transition(0.05, source=0, target=1)
simulator.set_immunity_transition(0.05, source=0, target=2)
simulator.set_contact_density(0.7, population=0)
simulator.set_contact_density(0.7, population=1)
simulator.set_migration_probability(2/365/2, source=0, target=2)
simulator.set_migration_probability(2/365/2, source=1, target=2)

simulator.simulate(1000, method='tau')

simulator.genealogy()

os.chdir('testing')
os.mkdir('output_example')
os.chdir('output_example')

simulator.output_newick()
simulator.output_mutations('mutations')
simulator.output_migrations('migrations')
