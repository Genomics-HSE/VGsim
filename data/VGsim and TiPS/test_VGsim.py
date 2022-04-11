import time
import VGsim

number_of_sites = 0
populations_number = 1
number_of_susceptible_groups = 2
simulator = VGsim.Simulator(number_of_sites, populations_number, number_of_susceptible_groups, seed=12312)

simulator.set_transmission_rate(0.25)
simulator.set_recovery_rate(0.099)
simulator.set_sampling_rate(0.001)
simulator.set_susceptibility_type(1, haplotype=None)

simulator.set_population_size(1000000000, population=0)
t1 = time.time()

simulator.simulate(100000000)

all_time_forward = 0
t2 = time.time()
all_time_forward += t2-t1
print("Forward time is", all_time_forward)

all_time_backward = 0
t1 = time.time()
simulator.genealogy()
t2 = time.time()
all_time_backward += t2 - t1
print("Backward time is", all_time_backward)