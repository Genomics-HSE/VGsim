from ._BirthDeath import BirthDeathModel
from .IO import writeGenomeNewick, writeMutations
from random import randrange
import sys
import matplotlib.pyplot as plt
import numpy as np

class Simulator:
	def __init__(self, number_of_sites=0, populations_number=1, number_of_susceptible_groups=2, seed=None, sampling_probability=False, strong_migration=False):
		self.fig = None
		self.first_sim = False
		if seed == None:
			seed = randrange(sys.maxsize)

		self.simulation = BirthDeathModel(number_of_sites, populations_number, number_of_susceptible_groups, seed, sampling_probability, strong_migration)


	def set_tau(self, tau):
		self.simulation.set_tau(tau)

	def print_basic_parameters(self):
		self.simulation.print_basic_parameters()

	def print_populations(self):
		self.simulation.print_populations()

	def print_immunity_model(self):
		self.simulation.print_immunity_model()

	def print_all(self, basic_parameters=True, populations=True, immunity_model=True):
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


	def set_susceptibility_type(self, susceptibility_type, haplotype=None):
		self.simulation.set_susceptibility_type(susceptibility_type, haplotype)

	def set_susceptibility(self, rate, haplotype=None, susceptibility_type=None):
		self.simulation.set_susceptibility(rate, haplotype, susceptibility_type)

	def set_immunity_transition(self, rate, source=None, target=None):
		self.simulation.set_immunity_transition(rate, source, target)


	def set_population_size(self, size, population=None):
		if self.first_sim == False:
			self.simulation.set_population_size(size, population)
		else:
			print("#TODO")

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


	def simulate(self, iterations=1000, sample_size=None, time=-1, method='direct'):
		self.first_sim = True
		if sample_size==None:
			sample_size = iterations
		if method == 'direct':
			self.simulation.SimulatePopulation(iterations, sample_size, time)
			self.simulation.Stats()
		elif method == 'tau':
			self.simulation.SimulatePopulation_tau(iterations)
		else:
			print("Unknown method. Choose between 'direct' and 'tau'.")

	def get_multi_events(self, id=None):
		self.simulation.Get_MultiEvents(id)

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

	def output_parameters(self, name_file="parameters"):
		self.simulation.output_parameters(name_file)

	def get_chain_events(self, name_file=None):
		self.simulation.get_chain_events(name_file)

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


	def add_plot_infectious(self, population, haplotype, step_num=100):
		if self.fig == None:
			self.fig, self.ax = plt.subplots(figsize=(8, 6))
			self.ax.set_ylabel('Number of people')
			self.ax.set_xlabel('Time')
			self.ax_2 = self.ax.twinx()
			self.ax_2.set_ylabel('Sampling')

		infections, sample, time_points = self.simulation.get_data_infectious(population, haplotype, step_num)
		self.ax.plot(time_points, infections, label='Infections-' + str(population) + '-' + str(haplotype))
		self.ax_2.plot(time_points, sample, "--", label='Sampling-' + str(population) + '-' + str(haplotype))

	def add_plot_susceptible(self, population, susceptibility_type, step_num=100):
		if self.fig == None:
			self.fig, self.ax = plt.subplots(figsize=(8, 6))
			self.ax.set_ylabel('Number of people')
			self.ax.set_xlabel('Time')
			self.ax_2 = self.ax.twinx()
			self.ax_2.set_ylabel('Sampling')

		susceptible, time_points = self.simulation.get_data_susceptible(population, susceptibility_type, step_num)
		self.ax.plot(time_points, susceptible, label='Susceptible-' + str(population) + '-' + str(susceptibility_type))

	def plot(self):
		plt.show()

	def add_title(self, name="Plot"):
		self.ax.set_title(name)

	def add_legend(self):
		lines_1, labels_1 = self.ax.get_legend_handles_labels()
		lines_2, labels_2 = self.ax_2.get_legend_handles_labels()
		lines = lines_1 + lines_2
		labels = labels_1 + labels_2
		self.ax.legend(lines, labels, loc=0)


	def citation(self):
		print("VGsim: scalable viral genealogy simulator for global pandemic")
		print("Vladimir Shchur, Vadim Spirin, Victor Pokrovskii, Evgeni Burovski, Nicola De Maio, Russell Corbett-Detig")
		print("medRxiv 2021.04.21.21255891; doi: https://doi.org/10.1101/2021.04.21.21255891")


	def debug(self):
		self.simulation.Debug()
