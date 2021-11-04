from ._BirthDeath import BirthDeathModel
from .IO import writeGenomeNewick, writeMutations
from random import randrange
import sys
import matplotlib.pyplot as plt
import numpy as np

class Simulator:
	def __init__(self, number_of_sites=0, populations_number=1, number_of_susceptible_groups=2, seed=None, sampling_probability=False):
		if seed == None:
			seed = randrange(sys.maxsize)

		self.simulation = BirthDeathModel(number_of_sites, populations_number, number_of_susceptible_groups, seed, sampling_probability)


	def print_basic_parameters(self):
		self.simulation.print_basic_parameters()

	def print_populations(self):
		self.simulation.print_populations()

	def print_immunity_model(self):
		self.simulation.print_immunity_model()

	def print_all(self, basic_parameters=False, populations=False, immunity_model=False):
		if basic_parameters:
			self.simulation.print_basic_parameters()
		if populations:
			self.simulation.print_populations()
		if immunity_model:
			self.simulation.print_immunity_model()


	def set_transmission_rate(self, rate, haplotype=None):
		self.simulation.set_transmission_rate(rate, haplotype)

	def set_recovery_rate(self, rate, haplotype=None):
		self.simulation.set_recovery_rate(rate, haplotype)

	def set_sampling_rate(self, rate, haplotype=None):
		self.simulation.set_sampling_rate(rate, haplotype)

	def set_mutation_rate(self, rate=None, probabilities=None, haplotype=None, mutation=None):
		self.simulation.set_mutation_rate(rate, probabilities, haplotype, mutation)


	def set_immunity_type(self, susceptibility_type, haplotype=None):
		self.simulation.set_immunity_type(susceptibility_type, haplotype)

	def set_susceptibility(self, rate, haplotype=None, susceptibility_type=None):
		self.simulation.set_susceptibility(rate, haplotype, susceptibility_type)

	def set_immunity_transition(self, rate, source=None, target=None):
		self.simulation.set_immunity_transition(rate, source, target)


	def set_population_size(self, size, population=None):
		self.simulation.set_population_size(size, population)

	def set_contact_density(self, value, population=None):
		self.simulation.set_contact_density(value, population)

	def set_susceptible_individuals(self, amount, source_type, target_type, population=None):
		self.simulation.set_susceptible_individuals(amount, source_type, target_type, population)

	def set_lockdown(self, parameters, population=None):
		self.simulation.set_lockdown(parameters, population)

	def set_sampling_multiplier(self, multiplier, population=None):
		self.simulation.set_sampling_multiplier(multiplier, population)

	def set_migration_probability(self, probability=None, total_probability=None, source=None, target=None):
		self.simulation.set_migration_probability(probability, total_probability, source, target)


	def simulate(self, iterations=1000, sample_size=None, time=-1):
		if sample_size==None:
			sample_size = iterations
		self.simulation.SimulatePopulation(iterations, sample_size, time)
		self.simulation.Stats()

	def genealogy(self, seed=None):
		self.simulation.GetGenealogy(seed)


	def output_newick(self, name_file="newick_output"):
		pruferSeq, times, mut, populations = self.simulation.output_tree_mutations()
		writeGenomeNewick(pruferSeq, times, populations, name_file)

	def output_mutations(self, name_file="mutation_output"):
		pruferSeq, times, mut, populations = self.simulation.output_tree_mutations()
		writeMutations(mut, len(pruferSeq), name_file)

	def output_migrations(self, name_file="migrations"):
		self.simulation.writeMigrations(name_file)

	def sample_data(self, output_print=False):
		time, pop, hap = self.simulation.sample_data()
		if output_print:
			return time, pop, hap
		else:
			print(time)
			print(pop)
			print(hap)

	def epidemiology_timelines(self, step=1000, output_file=False):
		if output_file == True:
			self.simulation.LogDynamics(step, output_file)
		else:
			return self.simulation.LogDynamics(step, output_file)


	def add_plot_infectious(self, population, haplotype, step_num):
		if fig == None:
			fig, ax = plt.subplots(figsize=(8, 6))
			ax.set_ylabel('Number of people')
			ax.set_xlabel('Time')
			ax_2 = ax.twinx()
			ax_2.set_ylabel('Sampling')

		infections, sample, time_points = self.simulation.get_data_infectious(population, haplotype, step_num)
		ax.plot(time_points, infections, label='Infections')
		ax_2.plot(time_points, sample, label='Sampling')

	def add_plot_susceptible(self, population, susceptibility_type, step_num):
		if fig == None:
			fig, ax = plt.subplots(figsize=(8, 6))
			ax.set_ylabel('Number of people')
			ax.set_xlabel('Time')
			ax_2 = ax.twinx()
			ax_2.set_ylabel('Sampling')

		susceptible, time_points = self.simulation.get_data_susceptible(population, susceptibility_type, step_num)
		ax.plot(time_points, susceptible, label='Susceptible')

	def plot(self):
		plt.show()

	def add_title(self, name="Plot"):
		ax.set_title(name)

	def add_legend_infectious(self):
		lines_1, labels_1 = ax.get_legend_handles_labels()
		lines_2, labels_2 = ax_2.get_legend_handles_labels()
		lines = lines_1 + lines_2
		labels = labels_1 + labels_2
		ax.legend(lines, labels, loc=0)


	def citation(self):
		print("VGsim: scalable viral genealogy simulator for global pandemic")
		print("Vladimir Shchur, Vadim Spirin, Victor Pokrovskii, Evgeni Burovski, Nicola De Maio, Russell Corbett-Detig")
		print("medRxiv 2021.04.21.21255891; doi: https://doi.org/10.1101/2021.04.21.21255891")


	def debug(self):
		self.simulation.Debug()
