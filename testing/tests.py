import VGsim

# set_initial_haplotype
model = VGsim.Simulator(number_of_sites=3, memory_optimization=True)
model.set_initial_haplotype(0)
model.set_initial_haplotype(1)
model.set_initial_haplotype(64)
model.set_initial_haplotype(65)

try:
    model = VGsim.Simulator(number_of_sites=3, memory_optimization=True)
    model.set_initial_haplotype(1/2)
except TypeError:
    pass
    
try:
    model = VGsim.Simulator(number_of_sites=3, memory_optimization=True)
    model.set_initial_haplotype(-1)
except ValueError:
    pass
    

# set_step_haplotype
model = VGsim.Simulator(number_of_sites=3, memory_optimization=True)
model.set_step_haplotype(0)
model.set_step_haplotype(1)
model.set_step_haplotype(64)
model.set_step_haplotype(65)

try:
    model = VGsim.Simulator(number_of_sites=3, memory_optimization=True)
    model.set_step_haplotype(1/2)
except TypeError:
    pass
    
try:
    model = VGsim.Simulator(number_of_sites=3, memory_optimization=True)
    model.set_step_haplotype(-1)
except ValueError:
    pass


# set_transmission_rate
model = VGsim.Simulator(number_of_sites=2)
model.set_transmission_rate(0)
model.set_transmission_rate(0.002, 'A*')
model.set_transmission_rate(0.003, 'AT')
model.set_transmission_rate(0.004, 0)
model.set_transmission_rate(0.005, 15)

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_transmission_rate('str')
except TypeError:
    pass

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_transmission_rate(-1)
except ValueError:
    pass

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_transmission_rate(0.01, haplotype=[1, 2])
except TypeError:
    pass
    
try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_transmission_rate(0.01, haplotype=-1)
except IndexError:
    pass
    
try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_transmission_rate(0.01, haplotype=16)
except IndexError:
    pass
    
try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_transmission_rate(0.01, haplotype='str')
except ValueError:
    pass
    
try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_transmission_rate(0.01, haplotype='AAA')
except ValueError:
    pass


# set_recovery_rate
model = VGsim.Simulator(number_of_sites=2)
model.set_recovery_rate(0)
model.set_recovery_rate(0.002, 'A*')
model.set_recovery_rate(0.003, 'AT')
model.set_recovery_rate(0.004, 0)
model.set_recovery_rate(0.005, 15)

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_recovery_rate('str')
except TypeError:
    pass

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_recovery_rate(-1)
except ValueError:
    pass

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_recovery_rate(0.01, haplotype=[1, 2])
except TypeError:
    pass
    
try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_recovery_rate(0.01, haplotype=-1)
except IndexError:
    pass
    
try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_recovery_rate(0.01, haplotype=16)
except IndexError:
    pass
    
try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_recovery_rate(0.01, haplotype='str')
except ValueError:
    pass
    
try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_recovery_rate(0.01, haplotype='AAA')
except ValueError:
    pass


# set_sampling_rate
model = VGsim.Simulator(number_of_sites=2)
model.set_sampling_rate(0)
model.set_sampling_rate(0.002, 'A*')
model.set_sampling_rate(0.003, 'AT')
model.set_sampling_rate(0.004, 0)
model.set_sampling_rate(0.005, 15)

model = VGsim.Simulator(number_of_sites=2, sampling_probability=True)
model.set_sampling_rate(0)
model.set_sampling_rate(0.002, 'A*')
model.set_sampling_rate(0.003, 'AT')
model.set_sampling_rate(0.004, 0)
model.set_sampling_rate(0.005, 15)
model.set_sampling_rate(1)

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_sampling_rate('str')
except TypeError:
    pass

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_sampling_rate(-1)
except ValueError:
    pass

try:
    model = VGsim.Simulator(number_of_sites=2, sampling_probability=True)
    model.set_sampling_rate('str')
except TypeError:
    pass

try:
    model = VGsim.Simulator(number_of_sites=2, sampling_probability=True)
    model.set_sampling_rate(-1)
except ValueError:
    pass

try:
    model = VGsim.Simulator(number_of_sites=2, sampling_probability=True)
    model.set_sampling_rate(2)
except ValueError:
    pass

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_sampling_rate(0.01, haplotype=[1, 2])
except TypeError:
    pass
    
try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_sampling_rate(0.01, haplotype=-1)
except IndexError:
    pass
    
try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_sampling_rate(0.01, haplotype=16)
except IndexError:
    pass
    
try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_sampling_rate(0.01, haplotype='str')
except ValueError:
    pass
    
try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_sampling_rate(0.01, haplotype='AAA')
except ValueError:
    pass


model = VGsim.Simulator(number_of_sites=2)
model.set_mutation_rate(rate=0)
model.set_mutation_rate(rate=0.002, haplotype='A*')
model.set_mutation_rate(rate=0.003, haplotype='AT')
model.set_mutation_rate(rate=0.004, haplotype=0)
model.set_mutation_rate(rate=0.005, haplotype=15)
model.set_mutation_rate(rate=0.006, mutation=0)
model.set_mutation_rate(rate=0.007, mutation=1)
model.set_mutation_rate(probabilities=[5, 4, 8, 2], haplotype='A*')
model.set_mutation_rate(probabilities=[5, 4, 8, 2], haplotype='AT')
model.set_mutation_rate(probabilities=[5, 4, 8, 2], haplotype=0)
model.set_mutation_rate(probabilities=[5, 4, 8, 2], haplotype=15)
model.set_mutation_rate(probabilities=[5, 4, 8, 2], mutation=0)
model.set_mutation_rate(probabilities=[5, 4, 8, 2], mutation=1)

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_mutation_rate(rate='str')
except TypeError:
    pass

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_mutation_rate(rate=-1)
except ValueError:
    pass

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_mutation_rate(probabilities=[5, -1, 8, 2])
except ValueError:
    pass

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_mutation_rate(probabilities=[0, 0, 0, 0])
except ValueError:
    pass

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_mutation_rate(probabilities=2)
except TypeError:
    pass

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_mutation_rate(rate=0.002, haplotype='R*')
except ValueError:
    pass

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_mutation_rate(rate=0.002, haplotype=-1)
except IndexError:
    pass

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_mutation_rate(rate=0.002, haplotype=[1, 2])
except TypeError:
    pass

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_mutation_rate(rate=0.002, haplotype=16)
except IndexError:
    pass

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_mutation_rate(rate=0.002, mutation=-1)
except IndexError:
    pass

try:
    model = VGsim.Simulator(number_of_sites=2)
    model.set_mutation_rate(rate=0.002, mutation=2)
except IndexError:
    pass

