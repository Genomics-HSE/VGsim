import VGsim
import sys

def check(num):
    ref = np.load("reference_" + str(num) + ".npy")
	test = np.load("test_" + str(num) + ".npy")

	for i in range(ref.shape[1]):
		if ref[0, i] != test[0, i] or ref[1, i] != test[1, i] or ref[2, i] != test[2, i] or ref[3, i] != test[3, i] or ref[4, i] != test[4, i] or ref[5, i] != test[5, i]:
			print("Error!")
			print(num)
			print("--------------")
			print('ref-time: ' + str(ref[0, i]))
			print('ref-type: ' + str(ref[1, i]))
			print('ref-haplotype: ' + str(ref[2, i]))
			print('ref-population: ' + str(ref[3, i]))
			print('ref-newHaplotype: ' + str(ref[4, i]))
			print('ref-newPopulation: ' + str(ref[5, i]))
			print("--------------")
			print('test-time: ' + str(test[0, i]))
			print('test-type: ' + str(test[1, i]))
			print('test-haplotype: ' + str(test[2, i]))
			print('test-population: ' + str(test[3, i]))
			print('test-newHaplotype: ' + str(test[4, i]))
			print('test-newPopulation: ' + str(test[5, i]))
			print("--------------")
			sys.exit(1)


#First test
simulator = VGsim.Simulator(seed=1234)

name = 1
simulator.simulate(10000000)
simulator.output_chain_events('test_' + str(name))
check_data(name)


#Second test
number_of_sites = 2
simulator = VGsim.Simulator(number_of_sites, seed=1234)

simulator.set_transmission_rate(0.25)
simulator.set_recovery_rate(0.099)
simulator.set_sampling_rate(0.001)
simulator.set_transmission_rate(0.5, haplotype="GG")
mutation_rate=0.000003
substitution_weights=[1,1,1,2]#ATCG
simulator.set_mutation_rate(mutation_rate, substitution_weights)
simulator.set_mutation_rate(3*mutation_rate, haplotype="G*", mutation=1)

name = 2
simulator.simulate(10000000)
simulator.output_chain_events('test_' + str(name))
check_data(name)


#Third test
populations_number = 4
simulator = VGsim.Simulator(population_number=populations_number, seed=1234)

simulator.set_population_size(10000000, population=0)
simulator.set_population_size(5000000, population=1)
simulator.set_population_size(1000000, population=2)
simulator.set_sampling_multiplier(3, population=1)
simulator.set_sampling_multiplier(0, population=2)
simulator.set_lockdown([0.1, 0.01, 0.002])
simulator.set_migration_probability(10/365/2)
simulator.set_contact_density(0.7, population=0)
simulator.set_contact_density(0.7, population=1)
simulator.set_migration_probability(2/365/2, source=0, target=2)
simulator.set_migration_probability(2/365/2, source=1, target=2)

name = 3
simulator.simulate(10000000)
simulator.output_chain_events('test_' + str(name))
check_data(name)


#Fourth test
number_of_susceptible_groups = 3
simulator = VGsim.Simulator(number_of_susceptible_groups=number_of_susceptible_groups, seed=1234)

simulator.set_susceptibility_type(1)
simulator.set_susceptibility_type(2, haplotype="G*")
simulator.set_susceptibility(0.1, susceptibility_type=1)
simulator.set_susceptibility(0.5, susceptibility_type=1, haplotype="G*")
simulator.set_immunity_transition(1/90, source=1, target=0)
simulator.set_immunity_transition(1/180, source=2, target=0)
simulator.set_susceptibility(0.0, susceptibility_type=2)
simulator.set_immunity_transition(0.05, source=0, target=1)
simulator.set_immunity_transition(0.05, source=0, target=2)

name = 4
simulator.simulate(10000000)
simulator.output_chain_events('test_' + str(name))
check_data(name)


#Fifth test
number_of_sites = 2
populations_number = 3
simulator = VGsim.Simulator(number_of_sites=number_of_sites, populations_number=populations_number, seed=1234)

simulator.set_transmission_rate(0.25)
simulator.set_recovery_rate(0.099)
simulator.set_sampling_rate(0.001)
simulator.set_transmission_rate(0.5, haplotype="GG")
mutation_rate=0.000003
substitution_weights=[1,1,1,2]#ATCG
simulator.set_mutation_rate(mutation_rate, substitution_weights)
simulator.set_mutation_rate(3*mutation_rate, haplotype="G*", mutation=1)

simulator.set_population_size(10000000, population=0)
simulator.set_population_size(5000000, population=1)
simulator.set_population_size(1000000, population=2)
simulator.set_sampling_multiplier(3, population=1)
simulator.set_sampling_multiplier(0, population=2)
simulator.set_lockdown([0.1, 0.01, 0.002])
simulator.set_migration_probability(10/365/2)
simulator.set_contact_density(0.7, population=0)
simulator.set_contact_density(0.7, population=1)
simulator.set_migration_probability(2/365/2, source=0, target=2)
simulator.set_migration_probability(2/365/2, source=1, target=2)

name = 5
simulator.simulate(10000000)
simulator.output_chain_events('test_' + str(name))
check_data(name)

#Sixth test
number_of_sites = 2
number_of_susceptible_groups = 3
simulator = VGsim.Simulator(number_of_sites=number_of_sites, number_of_susceptible_groups=number_of_susceptible_groups, seed=1234)

simulator.set_transmission_rate(0.25)
simulator.set_recovery_rate(0.099)
simulator.set_sampling_rate(0.001)
simulator.set_transmission_rate(0.5, haplotype="GG")
mutation_rate=0.000003
substitution_weights=[1,1,1,2]#ATCG
simulator.set_mutation_rate(mutation_rate, substitution_weights)
simulator.set_mutation_rate(3*mutation_rate, haplotype="G*", mutation=1)

simulator.set_susceptibility_type(1)
simulator.set_susceptibility_type(2, haplotype="G*")
simulator.set_susceptibility(0.1, susceptibility_type=1)
simulator.set_susceptibility(0.5, susceptibility_type=1, haplotype="G*")
simulator.set_immunity_transition(1/90, source=1, target=0)
simulator.set_susceptibility(0.0, susceptibility_type=2)
simulator.set_immunity_transition(1/180, source=2, target=0)
simulator.set_immunity_transition(0.05, source=0, target=1)
simulator.set_immunity_transition(0.05, source=0, target=2)

name = 6
simulator.simulate(10000000)
simulator.output_chain_events('test_' + str(name))
check_data(name)

#Seventh test
populations_number = 3
number_of_susceptible_groups = 3
simulator = VGsim.Simulator(populations_number=populations_number, number_of_susceptible_groups=number_of_susceptible_groups, seed=1234)

simulator.set_population_size(10000000, population=0)
simulator.set_population_size(5000000, population=1)
simulator.set_population_size(1000000, population=2)
simulator.set_sampling_multiplier(3, population=1)
simulator.set_sampling_multiplier(0, population=2)
simulator.set_lockdown([0.1, 0.01, 0.002])
simulator.set_migration_probability(10/365/2)
simulator.set_contact_density(0.7, population=0)
simulator.set_contact_density(0.7, population=1)
simulator.set_migration_probability(2/365/2, source=0, target=2)
simulator.set_migration_probability(2/365/2, source=1, target=2)

simulator.set_susceptibility_type(1)
simulator.set_susceptibility_type(2, haplotype="G*")
simulator.set_susceptibility(0.1, susceptibility_type=1)
simulator.set_susceptibility(0.5, susceptibility_type=1, haplotype="G*")
simulator.set_immunity_transition(1/90, source=1, target=0)
simulator.set_susceptibility(0.0, susceptibility_type=2)
simulator.set_immunity_transition(1/180, source=2, target=0)
simulator.set_immunity_transition(0.05, source=0, target=1)
simulator.set_immunity_transition(0.05, source=0, target=2)

name = 7
simulator.simulate(10000000)
simulator.output_chain_events('test_' + str(name))
check_data(name)

#Eighth test
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
simulator.set_mutation_rate(3*mutation_rate, haplotype="G*", mutation=1)

simulator.set_population_size(10000000, population=0)
simulator.set_population_size(5000000, population=1)
simulator.set_population_size(1000000, population=2)
simulator.set_sampling_multiplier(3, population=1)
simulator.set_sampling_multiplier(0, population=2)
simulator.set_lockdown([0.1, 0.01, 0.002])
simulator.set_migration_probability(10/365/2)
simulator.set_contact_density(0.7, population=0)
simulator.set_contact_density(0.7, population=1)
simulator.set_migration_probability(2/365/2, source=0, target=2)
simulator.set_migration_probability(2/365/2, source=1, target=2)

simulator.set_susceptibility_type(1)
simulator.set_susceptibility_type(2, haplotype="G*")
simulator.set_susceptibility(0.1, susceptibility_type=1)
simulator.set_susceptibility(0.5, susceptibility_type=1, haplotype="G*")
simulator.set_immunity_transition(1/90, source=1, target=0)
simulator.set_susceptibility(0.0, susceptibility_type=2)
simulator.set_immunity_transition(1/180, source=2, target=0)
simulator.set_immunity_transition(0.05, source=0, target=1)
simulator.set_immunity_transition(0.05, source=0, target=2)

name = 8
simulator.simulate(10000000)
simulator.output_chain_events('test_' + str(name))

check_data(name)


#Ninth test
populations_number = 3
simulator = VGsim.Simulator(population_number=populations_number, seed=1234)

simulator.set_migration_probability(0.01, source=0, target=1)
simulator.set_migration_probability(0.02, source=2, target=1)

name = 9
simulator.simulate(10000000)
simulator.output_chain_events('test_' + str(name))
check_data(name)











