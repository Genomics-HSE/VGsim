Examples
========

.. code-block:: python

    import VGsim
    
    simulator = VGsim.Simulator(2,[300000, 300000, 300000], 3, 1234)
    simulator.set_infectious_rate(40)
    simulator.set_uninfectious_rate(15)
    simulator.set_sampling_rate(4)
    simulator.set_mutation_rate(0.2, [4,2,2])
    simulator.set_migration_probability(0.01)
    simulator.set_lockdown([0.5, 0.2, 0.05])
    simulator.set_sampling_multiplier(1.8)
    simulator.set_immunity_type(2)
    simulator.set_susceptibility(0.5, susceptibility_type=1)
    simulator.set_immunity_transition(0.00001)
    simulator.simulate(10000)
    simulator.debug()
        

