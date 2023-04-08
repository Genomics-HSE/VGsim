from ._BirthDeath import BirthDeathModel
from .IO import writeGenomeNewick, writeMutations
from random import randrange
import sys
import matplotlib.pyplot as plt
import numpy as np
import time

class Simulator:
	def __init__(self, number_of_sites=0, populations_number=1, number_of_susceptible_groups=1, seed=None, \
		sampling_probability=False, memory_optimization=False, genome_length=int(1e6), recombination_probability=0.0):
		self.fig = None
		if seed == None:
			seed = int(randrange(sys.maxsize))
		print('User seed:', seed)

		self.simulation = BirthDeathModel(number_of_sites=number_of_sites, populations_number=populations_number, \
			number_of_susceptible_groups=number_of_susceptible_groups, seed=seed, sampling_probability=sampling_probability, \
			memory_optimization=memory_optimization, genome_length=genome_length, recombination_probability=recombination_probability)


	def print_basic_parameters(self):
		self.simulation.print_basic_parameters()

	def print_populations(self):
		self.simulation.print_populations()

	def print_immunity_model(self):
		self.simulation.print_immunity_model()

	def print_mutations(self):
		self.simulation.print_mutations()

	def print_migrations(self):
		self.simulation.print_migrations()

	def print_all(self, basic_parameters=True, populations=True, immunity_model=True):
		if basic_parameters:
			self.simulation.print_basic_parameters()
		if populations:
			self.simulation.print_populations()
		if immunity_model:
			self.simulation.print_immunity_model()


	def get_indexes_from_haplotype(self, haplotype):
		return np.array(self.simulation.create_list_for_cycles(haplotype, self.simulation.haplotype_number))

	@property
	def seed(self):
		return self.simulation.seed

	@property
	def sampling_probability(self):
		return self.simulation.sampling_probability

	@property
	def memory_optimization(self):
		return self.simulation.memory_optimization

	@property
	def number_of_sites(self):
		return self.simulation.number_of_sites

	@property
	def haplotype_number(self):
		return self.simulation.haplotype_number

	@property
	def populations_number(self):
		return self.simulation.populations_number

	@property
	def number_of_susceptible_groups(self):
		return self.simulation.number_of_susceptible_groups

	@property
	def initial_haplotype(self):
		return self.simulation.initial_haplotype

	def set_initial_haplotype(self, amount):
		self.simulation.set_initial_haplotype(amount)

	@property
	def step_haplotype(self):
		return self.simulation.step_haplotype

	def set_step_haplotype(self, amount):
		self.simulation.set_step_haplotype(amount)

	@property
	def genome_length(self):
		return self.simulation.genome_length

	def set_genome_length(self, genome_length):
		self.simulation.set_genome_length(genome_length)

	@property
	def coinfection_parameters(self):
		return self.simulation.coinfection_parameters

	def set_coinfection_parameters(self, recombination):
		self.simulation.set_coinfection_parameters(recombination)

	@property
	def super_spread_rate(self):
		return self.simulation.super_spread_rate
		
	def set_super_spread_rate(self, rate, left, right, population=None):
		self.simulation.set_super_spread_rate(rate, left, right, population)

	def set_general_sampling(self, sampling_proportion, sampling_times):
		self.simulation.set_general_sampling(sampling_proportion, sampling_times)

	@property
	def transmission_rate(self):
		return self.simulation.transmission_rate

	def set_transmission_rate(self, rate, haplotype=None):
		self.simulation.set_transmission_rate(rate, haplotype)

	@property
	def recovery_rate(self):
		return self.simulation.recovery_rate

	def set_recovery_rate(self, rate, haplotype=None):
		self.simulation.set_recovery_rate(rate, haplotype)

	@property
	def sampling_rate(self):
		return self.simulation.sampling_rate

	def set_sampling_rate(self, rate, haplotype=None):
		self.simulation.set_sampling_rate(rate, haplotype)

	@property
	def mutation_rate(self):
		return self.simulation.mutation_rate
	
	def set_mutation_rate(self, rate, haplotype=None, mutation=None):
		self.simulation.set_mutation_rate(rate, haplotype, mutation)

	@property
	def mutation_probabilities(self):
		return self.simulation.mutation_probabilities

	def set_mutation_probabilities(self, probabilities, haplotype=None, mutation=None):
		self.simulation.set_mutation_probabilities(probabilities, haplotype, mutation)

	@property
	def mutation_position(self):
		return self.simulation.mutation_position

	def set_mutation_position(self, mutation, position):
		self.simulation.set_mutation_position(mutation, position)


	@property
	def susceptibility_type(self):
		return self.simulation.susceptibility_type

	def set_susceptibility_type(self, susceptibility_type, haplotype=None):
		self.simulation.set_susceptibility_type(susceptibility_type, haplotype)

	@property
	def susceptibility(self):
		return self.simulation.susceptibility

	def set_susceptibility(self, rate, haplotype=None, susceptibility_type=None):
		self.simulation.set_susceptibility(rate, haplotype, susceptibility_type)

	@property
	def immunity_transition(self):
		return self.simulation.immunity_transition

	def set_immunity_transition(self, rate, source=None, target=None):
		self.simulation.set_immunity_transition(rate, source, target)


	@property
	def population_size(self):
		return self.simulation.population_size

	def set_population_size(self, size, population=None):
		self.simulation.set_population_size(size, population)

	@property
	def contact_density(self):
		return self.simulation.contact_density

	def set_contact_density(self, value, population=None):
		self.simulation.set_contact_density(value, population)

	@property
	def npi(self):
		return self.simulation.npi

	def set_npi(self, parameters, population=None):
		self.simulation.set_npi(parameters, population)

	@property
	def sampling_multiplier(self):
		return self.simulation.sampling_multiplier

	def set_sampling_multiplier(self, multiplier, population=None):
		self.simulation.set_sampling_multiplier(multiplier, population)

	@property
	def migration_probability(self):
		return self.simulation.migration_probability

	def set_migration_probability(self, probability, source=None, target=None):
		self.simulation.set_migration_probability(probability, source, target)

	def set_total_migration_probability(self, total_probability):
		self.simulation.set_total_migration_probability(total_probability)

	@property
	def susceptible(self):
		return self.simulation.susceptible

	def set_susceptible(self, amount, source_type, target_type, population=None):
		self.simulation.set_susceptible(amount, source_type, target_type, population)

	@property
	def infectious(self):
		return self.simulation.infectious

	def set_infectious(self, amount, source_type, target_haplotype, population=None):
		self.simulation.set_infectious(amount, source_type, target_haplotype, population)


	def set_SIS_model(self, number_of_sites=0, populations_number=1, number_of_susceptible_groups=1, seed=None, \
		sampling_probability=False, memory_optimization=False, genome_length=int(1e6), recombination_probability=0.0):
		pass

	def set_SIR_model(self, number_of_sites=0, populations_number=1, number_of_susceptible_groups=1, seed=None, \
		sampling_probability=False, memory_optimization=False, genome_length=int(1e6), recombination_probability=0.0):
		pass


	def set_chain_events(self, file_name):
		self.simulation.set_chain_events(file_name)

	def set_settings(self, file_template):
		self.simulation.set_settings(file_template)

	def set_state(self, file_template):
		self.set_chain_events(file_template)
		self.set_settings(file_template)

	# def export_newick всё что в файлы это export
	def export_newick(self, file_template=None, file_path = None):
		pruferSeq, times, mut, populations = self.simulation.output_tree_mutations()
		writeGenomeNewick(pruferSeq, times, populations, file_template, file_path)

	def export_mutations(self, file_template=None, file_path = None):
		pruferSeq, times, mut, populations = self.simulation.output_tree_mutations()
		writeMutations(mut, len(pruferSeq), file_template, file_path)

	def export_migrations(self, file_template=None, file_path = None):
		self.simulation.export_migrations(file_template, file_path)

	def output_sample_data(self, output_print=False):
		time, pop, hap = self.simulation.output_sample_data()
		if output_print:
			return time, pop, hap
		else:
			print(time)
			print(pop)
			print(hap)

	def output_epidemiology_timelines(self, step=1000, output_file=False):
		if output_file:
			self.simulation.output_epidemiology_timelines(step, output_file)
		else:
			return self.simulation.output_epidemiology_timelines(step, output_file)

	def export_chain_events(self, file_name="chain_events"):
		self.simulation.export_chain_events(file_name)

	def export_settings(self, file_template="parameters"):
		self.simulation.export_settings(file_template)

	def export_state(self):
		self.output_chain_events()
		self.output_settings()

	def export_ts(self):
		return self.simulation.export_ts()

	def get_tree(self):
		return self.simulation.get_tree()

	def get_data_susceptible(self, population, susceptibility_type,
							 step_num):  # returns susceptible, time_points, lockdowns
		return self.simulation.get_data_susceptible(population, susceptibility_type, step_num)

	def get_data_infectious(self, population, haplotype,
							step_num):  # returns infections, sample, time_points, lockdowns
		return self.simulation.get_data_infectious(population, haplotype, step_num)

	def add_plot_infectious(self, population, haplotype, step_num=100, label_infectious=None, label_samples=None):
		if self.fig == None:
			self.fig, self.ax = plt.subplots(figsize=(8, 6))
			self.ax.set_ylabel('Number of samples')
			self.ax.set_xlabel('Time')
			self.ax_2 = self.ax.twinx()
			self.ax_2.set_ylabel('Number of individuals')

		if isinstance(haplotype, int) == True:
			self.plot_infectious(population, haplotype, step_num, label_infectious, label_samples)
		elif isinstance(haplotype, str) == True:
			haplotypes = self.simulation.create_list_haplotypes(haplotype)
			for hi in haplotypes:
				self.plot_infectious(population, hi, step_num, label_infectious, label_samples)
		else:
			print("#TODO")

	def plot_infectious(self, population, haplotype, step_num, label_infectious, label_samples):
		infections, sample, time_points, lockdowns = self.simulation.get_data_infectious(population, haplotype, step_num)
		if label_infectious == None:
			self.ax_2.plot(time_points, infections, label='Infectious pop:' + str(population) + ' hap:' + self.simulation.calculate_string(haplotype))
		elif isinstance(label_infectious, str) == True:
			self.ax_2.plot(time_points, infections, label=label_infectious)
		else:
			print("#TODO")

		if label_samples == None:
			self.ax.plot(time_points, sample, "--", label='Samples pop:' + str(population) + ' hap:' + self.simulation.calculate_string(haplotype))
		elif isinstance(label_samples, str) == True:
			self.ax.plot(time_points, sample, "--", label=label_samples)
		else:
			print("#TODO")

		if len(lockdowns) != 0:
			point = 0
			pointEnd = 0
			for ld in range(0, len(lockdowns), 2):
				if ld+1 == len(lockdowns):
					while time_points[point] < lockdowns[ld][1]:
						point += 1
					plt.fill_between(time_points[point:], infections[point:], alpha=0.2)
				else:
					while time_points[point] < lockdowns[ld][1]:
						point += 1
					while time_points[pointEnd] < lockdowns[ld+1][1]:
						pointEnd += 1
					if pointEnd == point:
						continue
					plt.fill_between(time_points[point:pointEnd+1], infections[point:pointEnd+1], alpha=0.2)

	def add_plot_susceptible(self, population, susceptibility_type, step_num=100, label_susceptible=None):
		if self.fig == None:
			self.fig, self.ax = plt.subplots(figsize=(8, 6))
			self.ax.set_ylabel('Number of samples')
			self.ax.set_xlabel('Time')
			self.ax_2 = self.ax.twinx()
			self.ax_2.set_ylabel('Number of individuals')

		susceptible, time_points, lockdowns = self.simulation.get_data_susceptible(population, susceptibility_type, step_num)
		if label_susceptible == None:
			self.ax_2.plot(time_points, susceptible, label='Susceptible pop:' + str(population) + ' sus:' + str(susceptibility_type))
		elif isinstance(label_susceptible, str) == True:
			self.ax_2.plot(time_points, susceptible, label=label_susceptible)
		else:
			print("#TODO")

		if len(lockdowns) != 0:
			point = 0
			pointEnd = 0
			for ld in range(0, len(lockdowns), 2):
				if ld+1 == len(lockdowns):
					while time_points[point] < lockdowns[ld][1]:
						point += 1
					plt.fill_between(time_points[point:], susceptible[point:], alpha=0.2)
				else:
					while time_points[point] < lockdowns[ld][1]:
						point += 1
					while time_points[pointEnd] < lockdowns[ld+1][1]:
						pointEnd += 1
					if pointEnd == point:
						continue
					plt.fill_between(time_points[point:pointEnd+1], susceptible[point:pointEnd+1], alpha=0.2)

	def add_legend(self):
		lines_1, labels_1 = self.ax.get_legend_handles_labels()
		lines_2, labels_2 = self.ax_2.get_legend_handles_labels()
		lines = lines_1 + lines_2
		labels = labels_1 + labels_2
		self.ax.legend(lines, labels, loc=0)

	def add_title(self, name="Plot"):
		self.ax.set_title(name)

	def plot(self, name_file=None):
		if name_file:
			plt.savefig(name_file)
		else:
			plt.show()
		self.fig = None


	def simulate(self, iterations=1000, sample_size=None, epidemic_time=-1, method='direct', attempts=200):
		if sample_size is None:
			sample_size = iterations
		if epidemic_time is None:
			epidemic_time = -1

		start_time = time.time()
		if method == 'direct':
			self.simulation.SimulatePopulation(iterations, sample_size, epidemic_time, attempts)
			self.simulation.Stats(time.time() - start_time)
		elif method == 'tau':
			self.simulation.SimulatePopulation_tau(iterations, sample_size, epidemic_time, attempts)
			self.simulation.Stats(time.time() - start_time)
		else:
			print("Unknown method. Choose between 'direct' and 'tau'.")

	def genealogy(self, seed=None):
		self.simulation.GetGenealogy(seed)

	def print_recomb(self, left, right):
		self.simulation.print_recomb(left, right)


	def citation(self):
		print("VGsim: scalable viral genealogy simulator for global pandemic")
		print("Vladimir Shchur, Vadim Spirin, Dmitry Sirotkin, EvgeniBurovski, Nicola De Maio, Russell Corbett-Detig")
		print("medRxiv 2021.04.21.21255891; doi: https://doi.org/10.1101/2021.04.21.21255891")

	def debug(self):
		self.simulation.Debug()

	def get_proportion(self):
		return self.simulation.get_proportion()

	def print_counters(self):
		self.simulation.PrintCounters()

	def print_propensities(self):
		self.simulation.PrintPropensities()
