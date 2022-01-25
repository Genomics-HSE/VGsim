from ._BirthDeath import BirthDeathModel
from .IO import writeGenomeNewick, writeMutations
from random import randrange
import sys
import matplotlib.pyplot as plt
import numpy as np

class Simulator:
	def __init__(self, number_of_sites=0, populations_number=1, number_of_susceptible_groups=1, seed=None, sampling_probability=False):
		self.fig = None
		if seed == None:
			seed = randrange(sys.maxsize)

		self.simulation = BirthDeathModel(number_of_sites, populations_number, number_of_susceptible_groups, seed, sampling_probability)


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
		self.simulation.set_population_size(size, population)

	def set_contact_density(self, value, population=None):
		self.simulation.set_contact_density(value, population)

	def set_lockdown(self, parameters, population=None):
		self.simulation.set_lockdown(parameters, population)

	def set_sampling_multiplier(self, multiplier, population=None):
		self.simulation.set_sampling_multiplier(multiplier, population)

	def set_migration_probability(self, probability=None, total_probability=None, source=None, target=None):
		self.simulation.set_migration_probability(probability, total_probability, source, target)


	def set_susceptible_individuals(self, amount, source_type, target_type, population=None):
		self.simulation.set_susceptible_individuals(amount, source_type, target_type, population)

	def set_infected_individuals(self, amount, source_haplotype, target_haplotype, population=None):
		self.simulation.set_infected_individuals(amount, source_haplotype, target_haplotype, population)

	def set_infection(self, amount, source_type, target_haplotype, population=None):
		self.simulation.set_infection(amount, source_type, target_haplotype, population)

	def set_recovery(self, amount, source_haplotype, target_type, population=None):
		self.simulation.set_recovery(amount, source_haplotype, target_type, population)


	def set_chain_events(self, file_name):
		self.simulation.set_chain_events(file_name)

	def set_settings(self, file_template):
		self.simulation.set_settings(file_template)

	def set_state(self, file_template):
		self.set_chain_events(file_template)
		self.set_settings(file_template)



	def output_newick(self, file_template="newick_output"):
		pruferSeq, times, mut, populations = self.simulation.output_tree_mutations()
		writeGenomeNewick(pruferSeq, times, populations, file_template)

	def output_mutations(self, file_template="mutation_output"):
		pruferSeq, times, mut, populations = self.simulation.output_tree_mutations()
		writeMutations(mut, len(pruferSeq), file_template)

	def output_migrations(self, file_template="migrations"):
		self.simulation.output_migrations(file_template)

	def output_sample_data(self, output_print=False):
		time, pop, hap = self.simulation.output_sample_data()
		if output_print:
			return time, pop, hap
		else:
			print(time)
			print(pop)
			print(hap)

	def output_epidemiology_timelines(self, step=1000, output_file=False):
		if output_file == True:
			self.simulation.output_epidemiology_timelines(step, output_file)
		else:
			return self.simulation.output_epidemiology_timelines(step, output_file)

	def output_chain_events(self, file_name="chain_events"):
		self.simulation.output_chain_events(file_name)

	def output_settings(self, file_template="parameters"):
		self.simulation.output_settings(file_template)

	def output_state(self):
		self.output_chain_events()
		self.output_settings()


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
			self.ax_2.plot(time_points, infections, label='Individuals-' + str(population) + '-' + str(haplotype))
		elif isinstance(label_infectious, str) == True:
			self.ax_2.plot(time_points, infections, label=label_infectious)
		else:
			print("#TODO")

		if label_samples == None:
			self.ax.plot(time_points, sample, "--", label='Samples-' + str(population) + '-' + str(haplotype))
		elif isinstance(label_label_samples, str) == True:
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
			self.ax_2.plot(time_points, susceptible, label='Susceptible-' + str(population) + '-' + str(susceptibility_type))
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

	def plot(self):
		plt.show()
		self.fig = None


	def simulate(self, iterations=1000, sample_size=None, time=-1, method='direct'):
		self.first_sim = True
		if sample_size is None:
			sample_size = -1
		if time is None:
			time = -1
		if method == 'direct':
			self.simulation.SimulatePopulation(iterations, sample_size, time)
			self.simulation.Stats()
		elif method == 'tau':
			self.simulation.SimulatePopulation_tau(iterations, sample_size, time)
		else:
			print("Unknown method. Choose between 'direct' and 'tau'.")

	def genealogy(self, seed=None):
		self.simulation.GetGenealogy(seed)


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
