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
	def __init__(self, infection=25, uninfection=10, Sr=1, Sp=None, num_sites=0, mut_rate=None, mut_target_rate=None, num_pop=1, size_pop=1000000, contact_density=1.0, total_mig_rate=0.0, lockdown=None, sampling_multiplier=None, susc_type=None, susceptible=None, susc_trans=None):
		#Get rates
		if infection > 0 and uninfection > 0 and Sr > 0:
			self.B_rate = [float(infection) for _ in range(4**num_sites)]
			if Sp != None:
				if 0 < Sp < 1:
					self.D_rate = [float(uninfection) * (1 - Sp) for _ in range(4**num_sites)]
					self.S_rate = [float(uninfection) * Sp for _ in range(4**num_sites)]
				else: 
					self.Error("\"Sr\" should be more 0 and less 1!")
			else:
				self.D_rate = [float(uninfection) for _ in range(4**num_sites)]
				self.S_rate = [float(Sr) for _ in range(4**num_sites)]
		else:
			self.Error("\"B\" or \"D\" or \"S\" should be more 0!")

		#Get mutations
		self.mutations = [[]]
		if num_sites != 0:
			if mut_rate != None and mut_target_rate != None:
				if 0<=mut_rate<=1 and len(mut_target_rate)==3 and mut_target_rate[0]>=0 and mut_target_rate[1]>=0 and mut_target_rate[2]>=0:
					self.mutations = [[[mut_rate] + mut_target_rate for _ in range(num_sites)] for _ in range(4**num_sites)]
				else:
					self.Error("\"mut_rate\" should be more 0 and less 1, lenght of \"mut_target_rate\" should be equal 3 and each element should be more or equal 0!")
			else:
				self.Error("Mutations model consist of \"num_sites\", \"mut_rate\" and \"mut_target_rate\"!")

		#Get populations
		if size_pop>=1 and contact_density>=0 and 0<=total_mig_rate<=1:
			self.populations = []
			self.migration = []
			for i in range(num_pop):
				self.migration.append([])
				self.populations.append(Population(size_pop, contact_density))
				for j in range(num_pop):
					if i == j:
						self.migration[i].append(float(0))
					else:
						self.migration[i].append(total_mig_rate/(num_pop-1))
		else: 
			self.Error("\"size_pop\" should be more 1, \"contact_dencity\" should be more 0 and \"total_mig_rate\" should be more or equal 0 and less or equal 1!")

		self.lockdowns = None
		#Get lockdowns
		if lockdown!=None:
			if lockdown[0]>=0 and 0<=lockdown[1]<=100 and 0<=lockdown[2]<=100 and len(lockdown)==3:
				self.lockdowns = []
				for _ in range(num_pop):
					self.lockdowns.append(Lockdown(lockdown[0], lockdown[1], lockdown[2]))
			else:
				self.Error("lenght of \"lockdown\" should be equal 3, first element should be more or equal 0 and second and third elements should be more or equal 0 and less or equal 100!")

		self.samplingMultiplier = None
		#Get sampling multipliers
		if sampling_multiplier != None:
			if sampling_multiplier>=0:
				self.samplingMultiplier = [sampling_multiplier for _ in range(num_pop)]
			else:
				self.Error("\"sampling_multiplier\" should be more or equal 0!")

		self.suscType = None
		self.susceptible = None
		self.susc = None
		#Get types of susceptible
		if susc_type != None and susceptible != None:
			for i in susceptible:
				if 0<=i<=1:
					pass
				else:
					self.Error("Each element of \"susceptible\" should be more or equal 0 and less or equal 1!")
			if 0<=susc_type<len(susceptible):
				self.suscType = [susc_type for _ in range(4**num_sites)]
				self.susceptible = [susceptible for _ in range(4**num_sites)]
			else:
				self.Error("\"susc_type\" should be more or equal 0 and less lenght of \"susceptible\"!")	
			self.susc = [self.suscType, self.susceptible]			

		self.suscTrans = None
		#Get matrix of transition of type of susceptible
		if susc_trans != None:
			if 0<=susc_trans<=1:
				self.suscTrans = []
				for i in range(len(susceptible)):
					self.suscTrans.append([])
					for j in range(len(susceptible)):
						if i == j:
							self.suscTrans[i].append(0)
						else:
							self.suscTrans[i].append(susc_trans/(len(susceptible)-1))
			else:
				self.Error("\"susc_trans\" should be more or equal 0 and less or equal 1!")

		#Other date
		self.pruferSeq = None

	def Error(self, text):
		print("ERROR:", text)
		sys.exit(1)

	def print_Rates(self):
		print("B D S", end="")
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

	def print_Pop(self):
		print("id size conDen ", end="")
		if self.lockdowns != None:
			print("conDenAfterLD startLD endLD ", end="")
		if self.samplingMultiplier != None:
			print("samplingMultiplier", end="")
		print()
		for i in range(len(self.populations)):
			print(i, self.populations[i].size, self.populations[i].contactDensity, end=" ")
			if self.lockdowns != None:
				print(self.lockdowns[i].conDenAfterLD, self.lockdowns[i].startLD, self.lockdowns[i].endLD, end=" ")
			if self.samplingMultiplier != None:
				print(self.samplingMultiplier[i], end="")
			print()

	def print_Mig(self):
		print("Migration rates")
		for i in range(len(self.migration)):
			for j in range(len(self.migration)):
				print(self.migration[i][j], end=" ")
			print()

	def print_Susc(self):
		if self.suscType != None:
			print("ST", end="")
			for i in range(len(self.susceptible[0])):
				print(" S" + str(i), end="")
			print()
			for i in range(len(self.susceptible)):
				print(self.suscType[i], self.susceptible[i])

	def print_SuscTrans(self):
		if self.suscTrans != None:
			print("SuscTrans")
			for i in range(len(self.susceptible[0])):
				for j in range(len(self.susceptible[0])):
					print(self.suscTrans[i][j], end=" ")
				print()

	def create_class(self, _seed=None):
		if _seed == None:
		    _seed = randrange(sys.maxsize)
		self.simulation = BirthDeathModel(self.B_rate, self.D_rate, self.S_rate, self.mutations, populationModel=[self.populations, self.migration], susceptible=self.susc, suscepTransition=self.suscTrans, lockdownModel=self.lockdowns, samplingMultiplier=self.samplingMultiplier, rndseed=int(_seed))

	def update_rates(self, total_mig_rate):
		self.simulation.UpdateMigration(float(total_mig_rate))

	def simulate(self, _iterations=1000, _sampleSize=None, _time=-1):
		if _sampleSize == None:
			_sampleSize = _iterations
		self.simulation.SimulatePopulation(_iterations, _sampleSize, _time)
		self.simulation.Report()

	def get_genealogy(self):
		self.simulation.GetGenealogy()

	def debug(self):
		self.simulation.Debug()

	def log_dynamics(self, step=1000):
		self.simulation.LogDynamics(step)

	def graph(self, pop, hap, step_num=10, option="True scale"):
		self.simulation.Graph(pop, hap,  step_num, option)

	def newick(self):
		if self.pruferSeq == None:
			pruferSeq, times, mut, populations = self.simulation.Output_tree_mutations()
		writeGenomeNewick(pruferSeq, times, populations)

	def mut(self):
		if self.pruferSeq == None:
			pruferSeq, times, mut, populations = self.simulation.Output_tree_mutations()
		writeMutations(mut, len(pruferSeq))

	def mig(self):
		self.simulation.writeMigrations()

	def sampleDate(self):
		time, pop, hap = self.simulation.sampleDate()
		print(time)
		print(pop)
		print(hap)

	def set_Infection(self, target, value):
		if value <= 0:
			print("Birth rate less than 0!")
			sys.exit(1)
		if target < 0 and target >= len(self.B_rate):
			print("TODO")
			sys.exit(1)
		self.B_rate[target] = float(value)

	def set_Uninfection(self, target, value):
		if value <= 0:
			print("Death rate less than 0!")
			sys.exit(1)
		if target < 0 and target >= len(self.D_rate):
			print("TODO")
			sys.exit(1)
		self.D_rate[target] = float(value)

	def set_S(self, target, value):
		if value <= 0:
			print("Sampling rate less than 0!")
			sys.exit(1)
		if target < 0 and target >= len(self.S_rate):
			print("TODO")
			sys.exit(1)
		self.S_rate[target] = float(value)

	def set_M(self, target_1, target_2, value):
		if value <= 0 and value > 1:
			print("Mutation rate less than 0 or more than 1!")
			sys.exit(1)
		if target_1 < 0 and target_1 >= len(self.mutations):
			print("TODO")
			sys.exit(1)
		if target_2 < 0 and target_2 >= len(self.mutations[0]):
			print("TODO")
			sys.exit(1)
		self.mutations[target_1][target_2][0] = float(value)

	def set_MutRate(self, target_1, target_2, target_3, value):
		if value < 0:
			print("Target mutation rate less than 0") #TODO
			sys.exit(1)
		if target_1 < 0 and target_1 >= len(self.mutations):
			print("TODO")
			sys.exit(1)
		if target_2 < 0 and target_2 >= len(self.mutations[0]):
			print("TODO")
			sys.exit(1)
		if target_3 < 1 and target_3 > 3:
			print("TODO")
			sys.exit(1)
		self.mutations[target_1][target_2][target_3] = float(value)

	def set_Migration(self, target_1, target_2, value):
		if value <= 0 and value > 1:
			print("Migration rate less than 0 or more than 1!")
			sys.exit(1)
		if target_1 < 0 and target_1 >= len(self.migRate):
			print("TODO")
			sys.exit(1)
		if target_2 < 0 and target_2 >= len(self.migRate[0]):
			print("TODO")
			sys.exit(1)
		self.migration[target_1][target_2] = float(value)

	def set_startLD(self, target, value):
		if value <= 0 and value > 100:
			print("Proportion of start LD less than 0 or more than 100!")
			sys.exit(1)
		if target < 0 and target >= len(self.lockdowns):
			print("TODO")
			sys.exit(1)
		self.lockdowns[target].startLD = float(value)

	def set_endLD(self, target, value):
		if value <= 0 and value > 100:
			print("Proportion of end LD less than 0 or more than 100!")
			sys.exit(1)
		if target < 0 and target >= len(self.lockdowns):
			print("TODO")
			sys.exit(1)
		self.lockdowns[target].endLD = float(value)

	def set_conDenAfterLD(self, target, value):
		if value <= 0:
			print("Contact density after LD less than 0!")
			sys.exit(1)
		if target < 0 and target >= len(self.lockdowns):
			print("TODO")
			sys.exit(1)
		self.lockdowns[target].conDenAfterLD = float(value)

	def set_suscType(self, target, value):
		if value < 0 and value >= len(self.susceptible[0]):
			print("Susc type less than 0 or more number of type susceptible!") #TODO
			sys.exit(1)
		if target < 0 and target >= len(self.susceptible):
			print("TODO")
			sys.exit(1)
		self.suscType[target] = int(value)

	def set_susceptible(self, target_1, target_2, value):
		if value < 0 and value > 1:
			print("Susceptible rate less than 0 or more than 1!") #TODO
			sys.exit(1)
		if target_1 < 0 and target_1 >= len(self.susceptible):
			print("TODO")
			sys.exit(1)
		if target_2 < 0 and target_2 >= len(self.susceptible[0]):
			print("TODO")
			sys.exit(1)
		self.susceptible[target_1][target_2] = value

	def set_suscTrans(self, target_1, target_2, value):
		if value < 0 and value > 1:
			print("Susc rate less than 0 or more than 1!") #TODO
			sys.exit(1)
		if target_1 < 0 and target_1 >= len(self.susceptible[0]):
			print("TODO")
			sys.exit(1)
		if target_2 < 0 and target_2 >= len(self.susceptible[0]):
			print("TODO")
			sys.exit(1)
		self.suscRate[target_1][target_2] = value