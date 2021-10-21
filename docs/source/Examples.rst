Examples
========

Базовый пример
--------------

В данном примере создаётся самая базовая модель, которая имеет 1 популяцию, 1 гаплотип и 2 типа уязвимости.

.. code-block:: python

    import VGsim.interface

    simulator = VGsim.interface.Simulator()
    simulator.initialize()
    simulator.simulate()
    simulator.epidemiology_timeline()
    simulator.log_dynamics()
    simulator.plot_infectious()


simulator.initialize() - инициализирует модель с теми параметрами, которые были заданы для класса Simulator()

simulator.simulate() - симулирует последовательность событий

simulator.epidemiology_timeline() - воссоздаёт по цепи событий генеалогическое древо на основе сэмплированных образцов

simulator.log_dynamics() - записывает информацию о динамике развития симуляции по всем популяциям

simulator.plot_infectious() - отображает общую динамику заболевших по всем гаплотипам и всем популяциям

Более сложная модель
--------------------

В данном примере создаётся модель, которая имеет 1 мутацию, то есть 4 гаплотипа, 3 популяции и 2 типа уязвимости.

.. code-block:: python

    import VGsim.interface

    simulator = VGsim.interface.Simulator(infectious_rate=40, uninfectious_rate=15, sampling_rate=4,
     sites_number=1, mutation_rate=0.2, mutation_probabilities=[4,2,2], populations_number=3, 
     population_size=300000, contact_density=1.0, total_migration_probability=0.2)
    simulator.initialize(seed=1234)
    simulator.simulate(iterations=200000, time=0.3)
    simulator.epidemiology_timeline()
    simulator.log_dynamics()
    simulator.plot_infectious()

    

