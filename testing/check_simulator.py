import VGsim
import sys
import os

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

	os.remove("test_" + str(num) + ".npy")

#1 model
simulator = VGsim.Simulator(seed=1234)

name = 1
simulator.set_transmission_rate(4.0)
simulator.set_recovery_rate(1.5)
simulator.set_sampling_rate(0.3)
simulator.simulate(100000)
simulator.output_chain_events('test_' + str(name))
# check_data(name)

#2 model
simulator = VGsim.Simulator(number_of_sites=1, seed=1234)

name = 2
simulator.set_transmission_rate(4, haplotype=3)
simulator.simulate(100000)
simulator.output_chain_events('test_' + str(name))
# check_data(name)

#3 model
simulator = VGsim.Simulator(number_of_susceptible_groups=2, seed=1234)

name = 3
simulator.set_susceptibility_type(1)
simulator.simulate(100000)
simulator.output_chain_events('test_' + str(name))
# check_data(name)

#4 model
simulator = VGsim.Simulator(number_of_sites=1, number_of_susceptible_groups=3, seed=1234)

name = 4
simulator.set_mutation_rate(0.01)
simulator.set_susceptibility_type(1, haplotype=0)
simulator.set_susceptibility_type(2, haplotype=1)
simulator.set_susceptibility_type(2, haplotype=2)
simulator.set_susceptibility_type(2, haplotype=3)
simulator.set_immunity_transition(0.01, source=0, target=1)
simulator.set_immunity_transition(0.01, source=1, target=2)
simulator.set_immunity_transition(0.02, source=2, target=1)
simulator.simulate(100000)
simulator.output_chain_events('test_' + str(name))
# check_data(name)

#5 model
simulator = VGsim.Simulator(populations_number=2, seed=1234)

name = 5
simulator.set_population_size(2000000)
simulator.set_contact_density(1.3, population=0)
simulator.set_contact_density(0.8, population=1)
simulator.set_migration_probability(0.01, source=0, target=1)
simulator.set_migration_probability(0.005, source=1, target=0)
simulator.simulate(100000)
simulator.output_chain_events('test_' + str(name))
# check_data(name)

#6 model
simulator = VGsim.Simulator(populations_number=3, seed=1234)

name = 6
simulator.set_migration_probability(0.01, source=0, target=1)
simulator.set_migration_probability(0.005, source=2, target=1)
simulator.simulate(100000)
simulator.output_chain_events('test_' + str(name))
# check_data(name)

#7 model
simulator = VGsim.Simulator(populations_number=3, seed=1234)

name = 7
simulator.set_migration_probability(0.01, source=0, target=1)
simulator.set_migration_probability(0.005, source=2, target=1)
simulator.set_sampling_multiplier(2.5, population=1)
simulator.set_sampling_multiplier(2, population=2)
simulator.set_npi([0.5, 30, 15], population=0)
simulator.simulate(100000)
simulator.output_chain_events('test_' + str(name))
# check_data(name)

#8 model
simulator = VGsim.Simulator(number_of_sites=2, seed=1234)

name = 8
simulator.set_mutation_rate(probabilities=[1, 0, 0, 1])
simulator.simulate(100000)
simulator.output_chain_events('test_' + str(name))
# check_data(name)

#9 model
simulator = VGsim.Simulator(number_of_sites=2, populations_number=3, number_of_susceptible_groups=3, seed=1234)

name = 9
simulator.set_transmission_rate(4.0, haplotype=12)
simulator.set_recovery_rate(1.5)
simulator.set_sampling_rate(0.3)
simulator.set_mutation_rate(0.01, probabilities=[1, 0, 0, 1])
simulator.set_migration_probability(0.01, source=0, target=1)
simulator.set_migration_probability(0.005, source=2, target=1)
simulator.set_sampling_multiplier(2.5, population=1)
simulator.set_sampling_multiplier(2, population=2)
simulator.set_npi([0.5, 30, 15], population=1)
simulator.set_susceptibility_type(1, haplotype=0)
simulator.set_susceptibility_type(2, haplotype=1)
simulator.set_susceptibility_type(2, haplotype=2)
simulator.set_susceptibility_type(2, haplotype=3)
simulator.set_immunity_transition(0.01, source=0, target=1)
simulator.set_immunity_transition(0.01, source=1, target=2)
simulator.set_immunity_transition(0.02, source=2, target=1)
simulator.simulate(100000)
simulator.output_chain_events('test_' + str(name))
# check_data(name)

