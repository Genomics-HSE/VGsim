from ._BirthDeath import BirthDeathModel
from .IO import writeGenomeNewick, writeMutations
from random import randrange
import sys
import matplotlib.pyplot as plt
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
	def __init__(self, infectious_rate=25, uninfectious_rate=10, sampling_rate=1, sampling_probability=None, sites_number=0, mutation_rate=None, mutation_probabilities=None, populations_number=1, population_size=1000000, contact_density=1.0, total_migration_probability=0.0, lockdown=None, sampling_multiplier=None, immunity_type=None, susceptibility=None, total_immunity_transition=None):
		#Get rates
		if infectious_rate >= 0 and uninfectious_rate >= 0 and sampling_rate >= 0:
			self.B_rate = [float(infectious_rate) for _ in range(4**sites_number)]
			if sampling_probability != None:
				if 0 <= sampling_probability <= 1:
					self.D_rate = [float(uninfectious_rate) * (1 - sampling_probability) for _ in range(4**sites_number)]
					self.S_rate = [float(uninfectious_rate) * sampling_probability for _ in range(4**sites_number)]
				else: 
					self.Error("\"Sr\" should be more 0 and less 1!")
			else:
				self.D_rate = [float(uninfectious_rate) for _ in range(4**sites_number)]
				self.S_rate = [float(sampling_rate) for _ in range(4**sites_number)]
		else:
			self.Error("\"B\" or \"D\" or \"S\" should be more 0!")

		#Get mutations
		self.mutations = [[]]
		self.num_sites = sites_number
		if sites_number != 0:
			if mutation_rate != None and mutation_probabilities != None:
				if 0<=mutation_rate and len(mutation_probabilities)==3 and mutation_probabilities[0]>=0 and mutation_probabilities[1]>=0 and mutation_probabilities[2]>=0:
					self.mutations = [[[mutation_rate] + mutation_probabilities for _ in range(sites_number)] for _ in range(4**sites_number)]
				else:
					self.Error("\"mutation_rate\" should be more 0 and less 1, lenght of \"mutation_probabilities\" should be equal 3 and each element should be more or equal 0!")
			else:
				self.Error("Mutations model consist of \"sites_number\", \"mutation_rate\" and \"mutation_probabilities\"!")

		#Get populations
		self.num_pop = populations_number
		if population_size>=1 and contact_density>=0 and 0<=total_migration_probability<=1:
			self.populations = []
			self.migration = []
			for i in range(populations_number):
				self.migration.append([])
				self.populations.append(Population(population_size, contact_density))
				for j in range(populations_number):
					if i == j:
						self.migration[i].append(float(0))
					else:
						self.migration[i].append(total_migration_probability/(populations_number-1))
		else: 
			self.Error("\"population_size\" should be more 1, \"contact_dencity\" should be more 0 and \"total_migration_probability\" should be more or equal 0 and less or equal 1!")

		self.lockdowns = None
		#Get lockdowns
		if lockdown!=None:
			if lockdown[0]>=0 and 0<=lockdown[1]<=100 and 0<=lockdown[2]<=100 and len(lockdown)==3:
				self.lockdowns = []
				for _ in range(populations_number):
					self.lockdowns.append(Lockdown(lockdown[0], lockdown[1], lockdown[2]))
			else:
				self.Error("lenght of \"lockdown\" should be equal 3, first element should be more or equal 0 and second and third elements should be more or equal 0 and less or equal 100!")

		self.samplingMultiplier = None
		#Get sampling multipliers
		if sampling_multiplier != None:
			if sampling_multiplier>=0:
				self.samplingMultiplier = [sampling_multiplier for _ in range(populations_number)]
			else:
				self.Error("\"sampling_multiplier\" should be more or equal 0!")

		self.suscType = None
		self.susceptible = None
		self.susc = None
		#Get types of susceptible
		if immunity_type != None and susceptibility != None:
			for i in susceptibility:
				if 0<=i<=1:
					pass
				else:
					self.Error("Each element of \"susceptibility\" should be more or equal 0 and less or equal 1!")
			if 0<=immunity_type<len(susceptibility):
				self.suscType = [immunity_type for _ in range(4**sites_number)]
				self.susceptible = [susceptibility for _ in range(4**sites_number)]
			else:
				self.Error("\"immunity_type\" should be more or equal 0 and less lenght of \"susceptibility\"!")	
			self.susc = [self.suscType, self.susceptible]			

		self.suscTrans = None
		#Get matrix of transition of type of susceptible
		if total_immunity_transition != None:
			if 0<=total_immunity_transition<=1:
				self.suscTrans = []
				for i in range(len(susceptibility)):
					self.suscTrans.append([])
					for j in range(len(susceptibility)):
						if i == j:
							self.suscTrans[i].append(0)
						else:
							self.suscTrans[i].append(total_immunity_transition/(len(susceptibility)-1))
			else:
				self.Error("\"total_immunity_transition\" should be more or equal 0 and less or equal 1!")

		#Other date
		self.pruferSeq = None

	def Error(self, text):
		print("ERROR:", text)
		sys.exit(1)

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
		self.print_migration_matrix()

	def print_migration_matrix(self):
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

	def initialize(self, _seed=None):
		if _seed == None:
		    _seed = randrange(sys.maxsize)
		self.simulation = BirthDeathModel(self.B_rate, self.D_rate, self.S_rate, self.mutations, populationModel=[self.populations, self.migration], susceptible=self.susc, suscepTransition=self.suscTrans, lockdownModel=self.lockdowns, samplingMultiplier=self.samplingMultiplier, rndseed=int(_seed))

	def update_migration(self, total_migration_probability):
		self.simulation.UpdateMigration(float(total_migration_probability))

	def simulate(self, _iterations=1000, _sampleSize=None, _time=-1):
		if _sampleSize == None:
			_sampleSize = _iterations
		self.simulation.SimulatePopulation(_iterations, _sampleSize, _time)
		self.simulation.Report()

	def epidemiology_timeline(self, _seed=None):
		self.simulation.GetGenealogy(_seed)

	def debug(self):
		self.simulation.Debug()

	def log_dynamics(self, step=1000, output_file=False):
		if output_file == True:
			self.simulation.LogDynamics(step, output_file)
		else:
			return self.simulation.LogDynamics(step, output_file)

	def plot_infectious(self, step_num=100, population=None, haplotype=None):
		if population == None and haplotype == None:
			infections, sample, time_points = self.simulation.get_data_infectious(0, 0, step_num)
			for i in range(1, self.num_pop):
				for j in range(1, 4**self.num_sites):
					intections_2, sample_2, time_points_2 = self.simulation.get_data_infectious(i, j, step_num)
					for k in range(step_num + 1):
						intections[k] += intections_2[k]
						sample[k] += sample_2[k]

		elif population == None:
			infections, sample, time_points = self.simulation.get_data_infectious(population, 0, step_num)
			for i in range(1, self.num_pop):
				intections_2, sample_2, time_points_2 = self.simulation.get_data_infectious(i, haplotype, step_num)
				for j in range(step_num + 1):
					intections[j] += intections_2[j]
					sample[j] += sample_2[j]

		elif haplotype == None:
			infections, sample, time_points = self.simulation.get_data_infectious(population, 0, step_num)
			for i in range(1, 4**self.num_sites):
				intections_2, sample_2, time_points_2 = self.simulation.get_data_infectious(pop, i, step_num)
				for j in range(step_num + 1):
					intections[j] += intections_2[j]
					sample[j] += sample_2[j]

		else:
			infections, sample, time_points = self.simulation.get_data_infectious(population, haplotype, step_num)
		
		figure, axis_1 = plt.subplots(figsize=(8, 6))
		axis_1.plot(time_points, infections, color='blue', label='Infections')
		axis_1.set_ylabel('Infections')
		axis_1.set_xlabel('Time')
		axis_2 = axis_1.twinx()
		axis_2.plot(time_points, sample, color='orange', label='Sampling')
		axis_2.set_ylabel('Sampling')
		lines_1, labels_1 = axis_1.get_legend_handles_labels()
		lines_2, labels_2 = axis_2.get_legend_handles_labels()
		lines = lines_1 + lines_2
		labels = labels_1 + labels_2
		axis_1.legend(lines, labels, loc=0)

		plt.show()


	def plot_susceptible(self, step_num=100, population=None, susceptibility_type=None):
		if population == None and susceptibility_type == None:
			susceptible, time_points = self.simulation.get_data_susceptible(0, 0, step_num)
			for i in range(1, self.num_pop):
				for j in range(1, len(self.susceptible[0])):
					susceptible_2, time_points_2 = self.simulation.get_data_susceptible(i, j, step_num)
					for k in range(step_num + 1):
						susceptible[k] += susceptible_2[k]

		elif population == None:
			susceptible, time_points = self.simulation.get_data_susceptible(population, 0, step_num)
			for i in range(1, self.num_pop):
				susceptible_2, time_points_2 = self.simulation.get_data_susceptible(i, haplotype, step_num)
				for j in range(step_num + 1):
					susceptible[j] += susceptible_2[j]

		elif susceptibility_type == None:
			susceptible, time_points = self.simulation.get_data_susceptible(population, 0, step_num)
			for i in range(1, len(self.susceptible[0])):
				susceptible_2, time_points_2 = self.simulation.get_data_susceptible(pop, i, step_num)
				for j in range(step_num + 1):
					susceptible[j] += susceptible_2[j]

		else:
			susceptible, time_points = self.simulation.get_data_susceptible(population, susceptibility_type, step_num)
		
		figure, axis_1 = plt.subplots(figsize=(8, 6))
		axis_1.plot(time_points, susceptible, color='blue', label='Susceptible')
		axis_1.set_ylabel('Susceptible')
		axis_1.set_xlabel('Time')
		axis_1.legend()

		plt.show()

	def output_newick(self, name_file="newick_output"):
		if self.pruferSeq == None:
			pruferSeq, times, mut, populations = self.simulation.Output_tree_mutations()
		writeGenomeNewick(pruferSeq, times, populations, name_file)

	def output_mutations(self, name_file="mutation_output"):
		if self.pruferSeq == None:
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

	def set_infectious_rate(self, haplotype, rate):
		if rate < 0:
			print("Infectious rate is less than 0!")
			sys.exit(1)
		if haplotype < 0 and haplotype >= len(self.B_rate):
			print("There is not this haplotype!")
			sys.exit(1)
		self.B_rate[haplotype] = float(rate)

	def set_uninfectious_rate(self, haplotype, rate):
		if rate < 0:
			print("Uninfectious rate is less than 0!")
			sys.exit(1)
		if haplotype < 0 and haplotype >= len(self.D_rate):
			print("There is not this haplotype!")
			sys.exit(1)
		self.D_rate[haplotype] = float(rate)

	def set_sampling_rate(self, haplotype, rate):
		if rate < 0:
			print("Sampling rate is less than 0!")
			sys.exit(1)
		if haplotype < 0 and haplotype >= len(self.S_rate):
			print("There is not this haplotype!")
			sys.exit(1)
		self.S_rate[haplotype] = float(rate)

	def set_mutation_rate(self, haplotype, site, rate=None, probabilities=None): #Question
		if haplotype < 0 and haplotype >= len(self.mutations):
			print("There is not this haplotype!")
			sys.exit(1)
		if site < 0 and site >= len(self.mutations[0]):
			print("There is not this site")
			sys.exit(1)
		if rate != None:
			if rate < 0:
				print("Mutation rate is less than 0!")
				sys.exit(1)
			else:
				self.mutations[haplotype][site][0] = rate
		if probabilities != None:
			for i in range(3):
				if probabilities[i] < 0:
					print("Probability is less than 0!")
					sys.exit(1)
				else:
					self.mutations[haplotype][site][1+i] = probabilities[i]

	def set_migration_probability(self, source_population, target_population, probability):
		if probability <= 0 and probability > 1:
			print("Migration probability is less than 0 or more than 1!")
			sys.exit(1)
		if source_population < 0 and source_population >= len(self.migRate): #To ask how better
			print("There is not this source population!")
			print("No given source population!")
			sys.exit(1)
		if target_population < 0 and target_population >= len(self.migRate[0]):
			print("There is not this target population!")
			sys.exit(1)
		total_sum = 0.0
		for i in range(self.num_pop):
			total_sum += self.migration[i][target_population]
		if total_sum > 1:
			print("Total sum is more than 1!")
			sys.exit(1)
		self.migration[source_population][target_population] = float(probability)

	def set_start_lockdown(self, population, infectious_fraction=None, contact_density=None):
		if infectious_fraction != None:
			if infectious_fraction <= 0 and infectious_fraction > 1:
				print("Proportion of lockdown start can't be less than 0 or more than 1!")
				sys.exit(1)
			if population < 0 and population >= len(self.lockdowns):
				print("There is not this population!")
				sys.exit(1)
			self.lockdowns[population].startLD = float(infectious_fraction)
		if contact_density != None:
			if contact_density <= 0:
				print("Contact density after lockdown can't be less than 0!")
				sys.exit(1)
			if population < 0 and population >= len(self.lockdowns):
				print("There is not this population!")
				sys.exit(1)
			self.lockdowns[population].conDenAfterLD = float(contact_density)

	def set_end_lockdown(self, population, infectious_fraction):
		if infectious_fraction <= 0 and infectious_fraction > 1:
			print("Proportion of end lockdown can't be less than 0 or more than 1!")
			sys.exit(1)
		if population < 0 and population >= len(self.lockdowns):
			print("There is not this population!")
			sys.exit(1)
		self.lockdowns[population].endLD = float(infectious_fraction)

	def set_immunity_type(self, haplotype, immunity):
		if immunity < 0 and immunity >= len(self.susceptible[0]):
			print("Immunity type can't be less than 0 or more number of type susceptible!")
			sys.exit(1)
		if haplotype < 0 and haplotype >= len(self.susceptible):
			print("There is not this haplotype!")
			sys.exit(1)
		self.suscType[haplotype] = int(immunity)

	def set_susceptibility(self, haplotype, immunity, susceptibility):
		if susceptibility < 0:
			print("Susceptible rate is less than 0 or more than 1!")
			sys.exit(1)
		if haplotype < 0 and haplotype >= len(self.susceptible):
			print("There are no such haplotype!")
			sys.exit(1)
		if immunity < 0 and immunity >= len(self.susceptible[0]):
			print("There is not this immunity!")
			sys.exit(1)
		self.susceptible[haplotype][immunity] = susceptibility

	def set_immunity_transition(self, source_immunity, target_immunity, probability):
		if probability < 0 and probability > 1:
			print("Immunity transition rate is less than 0 or more than 1!")
			sys.exit(1)
		total = 0.0
		for i in range(len(self.susceptible)):
			total += self.suscTrans[i][target_immunity]
		if total > 1:
			print("Total sum is more 1!")
			sys.exit(1)
		if source_immunity < 0 and source_immunity >= len(self.susceptible[0]):
			print("There is not this source immunity!")
			sys.exit(1)
		if target_immunity < 0 and target_immunity >= len(self.susceptible[0]):
			print("There is not this target immunity!")
			sys.exit(1)
		self.suscTrans[source_immunity][target_immunity] = probability

	def check_migration(self):
		self.simulation.check_ratio()

	def print_times(self):
		self.simulation.times_print()


