from ._BirthDeath import BirthDeathModel
from .IO import writeGenomeNewick, writeMutations
from random import randrange
import sys
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

class Population:
    def __init__(self, size=1000000, contactDensity=1.0):
        self.size = size
        self.contactDensity = contactDensity

class Lockdown:
    def __init__(self, conDenAfterLD=0.1, startLD=2, endLD=1):
        self.conDenAfterLD = conDenAfterLD
        self.startLD = startLD
        self.endLD = endLD

class Simulator:
	def __init__(self, sites_number=0, population_sizes=[1000000], susceptibility_types=2, seed=None):
		if seed == None:
			seed = randrange(sys.maxsize)
		self.simulation = BirthDeathModel(sites_number, population_sizes, susceptibility_types, seed)

	def simulate(self, iterations=1000, sampleSize=None, time=-1):
		if sampleSize==None:
			sampleSize = iterations
		self.simulation.SimulatePopulation(iterations, sampleSize, time)
		self.simulation.Stats()

	def set_infectious_rate(self, rate, haplotype=None):
		self.simulation.set_infectious_rate(rate, haplotype)

	def set_uninfectious_rate(self, rate, haplotype=None):
		self.simulation.set_uninfectious_rate(rate, haplotype)

	def set_sampling_rate(self, rate, haplotype=None):
		self.simulation.set_sampling_rate(rate, haplotype)

	def set_mutation_rate(self, rate=None, probabilities=None, haplotype=None, mutation=None):
		self.simulation.set_mutation_rate(rate, probabilities, haplotype, mutation)

	def set_contact_density(self, value, population=None):
		self.simulation.set_contact_density(value, population)

	def set_lockdown(self, parameters, population=None):
		self.simulation.set_lockdown(parameters, population)

	def set_sampling_multiplier(self, multiplier, population=None):
		self.simulation.set_sampling_multiplier(multiplier, population)

	def set_migration_rate(self, rate, from_population=None, to_population=None):
		self.simulation.set_migration_rate(rate, from_population, to_population)

	def set_susceptible(self, amount, from_type, to_type, population=None):
		self.simulation.set_susceptible(amount, from_type, to_type, population)

	def set_infections(self, amount, from_type, to_haplotype, population=None):
		self.simulation.set_infections(amount, from_type, to_haplotype, population)

	def set_infections_2(self, amount, from_haplotype, to_haplotype, population=None):
		self.simulation.set_infections(amount, from_haplotype, to_haplotype, population)

	def set_susceptibility_type(self, susceptibility_type, haplotype=None):
		self.simulation.set_susceptibility_type(susceptibility_type, haplotype)

	def set_susceptibility(self, rate, haplotype=None, susceptibility_type=None):
		self.simulation.set_susceptibility(rate, haplotype, susceptibility_type)

	def set_immunity_transition(self, rate, from_population=None, to_population=None):
		self.simulation.set_immunity_transition(rate, from_population, to_population)

	def print_basic_rates(self):
		print("Infectious_rate Uninfectious_rate Sampling_rate", end="")
		if self.mutations != None:
			for i in range(len(self.mutations[0])):
				print(" M" + str(i), end="")
			a = 3 + len(self.mutations[0])
		else: 
			a = 3
		print()
		for i in range(len(self.B_rate)):
			print(self.B_rate[i], end=" ")
			print(self.D_rate[i], end=" ")
			print(self.S_rate[i], end=" ")
			if self.mutations != None:
				for j in range(len(self.mutations[0])):
					print(self.mutations[i][j], end=" ")
			print()

	def print_populations(self):
		print("Id Size Contact_dencity ", end="")
		if self.lockdowns != None:
			print("Contact_dencity_after_lockdown Start_lockdown End_lockdown ", end="")
		if self.samplingMultiplier != None:
			print("Sampling_multiplier", end="")
		print()
		for i in range(len(self.populations)):
			print(i, self.populations[i].size, self.populations[i].contactDensity, end=" ")
			if self.lockdowns != None:
				print(self.lockdowns[i].conDenAfterLD, self.lockdowns[i].startLD, self.lockdowns[i].endLD, end=" ")
			if self.samplingMultiplier != None:
				print(self.samplingMultiplier[i], end="")
			print()
		print("Migration matrix")
		for i in range(len(self.migration)):
			for j in range(len(self.migration)):
				print(self.migration[i][j], end=" ")
			print()

	def print_immunity_model(self):
		if self.suscType != None:
			print("Immunity model", end="")
			for i in range(len(self.susceptible[0])):
				print(" S" + str(i), end="")
			print()
			for i in range(len(self.susceptible)):
				print(self.suscType[i], self.susceptible[i])

		if self.suscTrans != None:
			print("Immunity transition rates")
			for i in range(len(self.susceptible[0])):
				for j in range(len(self.susceptible[0])):
					print(self.suscTrans[i][j], end=" ")
				print()

	def citation(self):
		print("VGsim: scalable viral genealogy simulator for global pandemic")
		print("Vladimir Shchur, Vadim Spirin, Victor Pokrovskii, Evgeni Burovski, Nicola De Maio, Russell Corbett-Detig")
		print("medRxiv 2021.04.21.21255891; doi: https://doi.org/10.1101/2021.04.21.21255891")

	def genealogy(self, seed=None):
		self.simulation.GetGenealogy(seed)

	def debug(self):
		self.simulation.Debug()

	def epidemiology_timelines(self, step=1000, output_file=False):
		if output_file == True:
			self.simulation.LogDynamics(step, output_file)
		else:
			return self.simulation.LogDynamics(step, output_file)

	def plot_infectious(self, population=None, haplotype=None, step_num=100):
		if population == None and haplotype == None:
			for i in range(0, self.simulation.get_popNum()):
				for j in range(0, self.simulation.get_hapNum()):
					infections, sample, time_points = self.simulation.get_data_infectious(i, j, step_num)
					self.paint_infections(i, j, time_points, infections, sample)
		elif population == None:
			for i in range(0, self.simulation.get_popNum()):
				infections, sample, time_points = self.simulation.get_data_infectious(i, haplotype, step_num)
				self.paint_infections(i, haplotype, time_points, infections, sample)
		elif haplotype == None:
			for i in range(0, self.simulation.get_hapNum()):
				infections, sample, time_points = self.simulation.get_data_infectious(population, i, step_num)
				self.paint_infections(population, i, time_points, infections, sample)
		else:
			infections, sample, time_points = self.simulation.get_data_infectious(population, haplotype, step_num)
			self.paint_infections(population, haplotype, time_points, infections, sample)
		
	def paint_infections(self, pop, hap, time_points, infections, sample):
		# figure, axis_1 = plt.subplots(figsize=(8, 6))
		# axis_1.plot(time_points, infections, color='blue', label='Infections')
		# axis_1.set_ylabel('Infections')
		# axis_1.set_xlabel('Time')
		# axis_1.set_title('Population ' + str(pop) + ' and hapotype ' + str(hap))
		# axis_2 = axis_1.twinx()
		# axis_2.plot(time_points, sample, color='orange', label='Sampling')
		# axis_2.set_ylabel('Sampling')
		# lines_1, labels_1 = axis_1.get_legend_handles_labels()
		# lines_2, labels_2 = axis_2.get_legend_handles_labels()
		# lines = lines_1 + lines_2
		# labels = labels_1 + labels_2
		# axis_1.legend(lines, labels, loc=0)

		# plt.show()
		pass

	def plot_susceptible(self, population=None, susceptibility_type=None, step_num=100):
		if population == None and susceptibility_type == None:
			for i in range(0, self.simulation.get_popNum()):
				for j in range(0, self.simulation.get_susNum()):
					susceptible, time_points = self.simulation.get_data_susceptible(i, j, step_num)
					self.paint_susceptible(i, j, susceptible, time_points)
		elif population == None:
			for i in range(0, self.simulation.get_popNum()):
				susceptible, time_points = self.simulation.get_data_susceptible(i, susceptibility_type, step_num)
				self.paint_susceptible(i, susceptibility_type, susceptible, time_points)	
		elif susceptibility_type == None:
			for i in range(0, self.simulation.get_susNum()):
				susceptible, time_points = self.simulation.get_data_susceptible(population, i, step_num)
				self.paint_susceptible(population, i, susceptible, time_points)
		else:
			susceptible, time_points = self.simulation.get_data_susceptible(population, susceptibility_type, step_num)
			self.paint_susceptible(population, susceptibility_type, susceptible, time_points)
		
	def paint_susceptible(self, population, susceptibility_type, susceptible, time_points): 
		# figure, axis_1 = plt.subplots(figsize=(8, 6))
		# axis_1.plot(time_points, susceptible, color='blue', label='Susceptible')
		# axis_1.set_ylabel('Susceptible')
		# axis_1.set_xlabel('Time')
		# axis_1.set_title('Population ' + str(population) + ' and susceptibility type ' + str(susceptibility_type))
		# axis_1.legend()

		# plt.show()
		pass

	def output_newick(self, name_file="newick_output"):
		pruferSeq, times, mut, populations = self.simulation.Output_tree_mutations()
		writeGenomeNewick(pruferSeq, times, populations, name_file)

	def output_mutations(self, name_file="mutation_output"):
		pruferSeq, times, mut, populations = self.simulation.Output_tree_mutations()
		writeMutations(mut, len(pruferSeq), name_file)

	def output_migrations(self, name_file="migrations"):
		self.simulation.writeMigrations(name_file)

	def sample_list(self, output_print=False):
		time, pop, hap = self.simulation.sampleDate()
		if output_print:
			return time, pop, hap
		else:
			print(time)
			print(pop)
			print(hap)

	def check_migration(self):
		self.simulation.check_ratio()

	def print_times(self):
		self.simulation.times_print()


