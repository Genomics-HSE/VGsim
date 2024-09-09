import VGsim

simulator = VGsim.Simulator(2, 2, 3, 1234)
simulator.simulate(1000000, "direct", 100)