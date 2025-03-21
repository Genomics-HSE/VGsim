import VGsim

# pip install --editable .

simulator = VGsim.Simulator(1, 1, 1, 1234)
# simulator.Debug()
simulator.simulate(1000000)
# simulator.Debug()
simulator.genealogy()