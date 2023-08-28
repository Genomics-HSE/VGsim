import VGsim
number_of_sites = 0
populations_number = 4
number_of_susceptible_groups = 2
simulator = VGsim.Simulator(number_of_sites, populations_number, number_of_susceptible_groups, seed=1234)

simulator.set_transmission_rate(0.25)
simulator.set_recovery_rate(0.099)
simulator.set_sampling_rate(0.001)

simulator.set_susceptibility_type(1)
simulator.set_susceptibility(0, susceptibility_type=1)

simulator.set_population_size(10000000, population=0)
simulator.set_population_size(10000000, population=1)
simulator.set_population_size(10000000, population=2)
simulator.set_population_size(10000000, population=3)


simulator.simulate(100000000, epidemic_time=180)

simulator.add_plot_susceptible(population=0, susceptibility_type=0, step_num=100)
simulator.add_plot_susceptible(population=0, susceptibility_type=1, step_num=100)
simulator.add_title(name="Susceptible group sizes")
simulator.add_legend()
simulator.plot()