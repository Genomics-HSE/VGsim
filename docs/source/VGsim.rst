VGsim
======

Our model of epidemiological spread is built as a compartmental model, and the random realisations of the corresponding stochastic process are drawn using the Gillespie algorithm. The different compartments in our model are defined based on several real-world complexities that interact with each other in a complex way and can affect epidemiological dynamics: population structure (which means defining separate host populations and assigning different frequencies to within-population and between-population contacts), separate infectious groups (which means that host individuals carrying different viral haplotypes are modeled differently, since some viral variants might be more transmissible than others), and different susceptible groups (different hosts having different types of immunity response to different haplotypes).


Interface
=========

VGsim.Simulator(sites_number=0, population_sizes=[1000000], susceptibility_types=2, seed=None)
----------------------------------------------------------------------------------------------

.. list-table::
   :widths: 15 15 70
   :header-rows: 1

   * - Basic parameters
     - Type
     - Meaning
   * - ``sites_number``
     - ``int``
     - Amount of parts of the genome where mutations can occur
   * - ``population_sizes``
     - ``list``
     - Population sizes
   * - ``susceptibility_types``
     - ``int``
     - Number of susceptibility types
   * - ``seed``
     - ``double``
     - #TODO


Show initial parameters
-----------------------

For example, this is the introductory paragraph

.. code-block:: python

  print_basic_parameters()
  print_populations()
  print_immunity_model()
  print_all(basic_parameters=False, populations=False, immunity_model=False)

Model simulation
================

simulate
--------


Parameters of model halt


.. list-table::
   :widths: 15 15 70
   :header-rows: 1

   * - Parameter
     - Type
     - Meaning
   * - ``iterations``
     - ``int``
     - Amount of iterations which will take place before the halt
   * - ``sampleSize``
     - ``int``
     - Halt the simulation when a particular number of participants was sampled
   * - ``time``
     - ``double``
     - time-based halt

genealogy
---------

.. list-table::
   :widths: 15 15 70
   :header-rows: 1

   * - Parameter
     - Type
     - Meaning
   * - ``step_num``
     - ``int``
     - A number of time intervals
   * - ``population``
     - ``int``
     - Population
   * - ``haplotype``
     - ``int``
     - Haplotype

plot
----
.. list-table::
   :widths: 15 25 70

   * - ``add_plot_infectious``
     - ``population, haplotype, step_num=1000``
     - #TODO
   * - ``add_plot_susceptible``
     - ``population, susceptibility_type, step_num=None``
     - #TODO
   * - ``plot``
     - ``None``
     - #TODO
   * - ``add_title``
     - ``name="Plot"``
     - #TODO
   * - ``add_legend_infectious``
     - ``None``
     - #TODO

System methods
--------------

.. code-block:: python

  output_newick(name_file="newick_output"): - record format of binary trees
  output_mutations(name_file="mutation_output"): - information about all mutations
  output_migrations(name_file="migrations"): - information about all migrations
  sampleDate(): - show all information about sampling, time, place, etc.
  epidemiology_timelines(self, step=1000, output_file=False): - records simulation state changes over some period of time. step - a number of parts epidemiology_timelines is split on.

Change values
-------------

.. list-table::
   :widths: 15 35 70
   :header-rows: 1

   * - Function
     - Parameters
     - Meaning
   * - ``set_infectious_rate``
     - ``rate, haplotype=None``
     - #TODO
   * - ``set_recovery_rate``
     - ``rate, haplotype=None``
     - #TODO
   * - ``set_sampling_rate``
     - ``rate, haplotype=None``
     - #TODO
   * - ``set_mutation_rate``
     - ``rate=None, probabilities=None, haplotype=None, mutation=None``
     - #TODO
   * - ``set_contact_density``
     - ``value, population=None``
     - #TODO
   * - ``set_lockdown``
     - ``parameters, population=None``
     - #TODO
   * - ``set_sampling_multiplier``
     - ``multiplier, population=None``
     - #TODO
   * - ``set_migration_probability``
     - ``rate, source=None, target=None``
     - #TODO
   * - ``set_susceptible``
     - ``amount, source_type, target_type, population=None``
     - #TODO
   * - ``set_immunity_type``
     - ``susceptibility_type, haplotype=None``
     - #TODO
   * - ``set_susceptibility``
     - ``rate, haplotype=None, susceptibility_type=None``
     - #TODO
   * - ``set_immunity_transition``
     - ``rate, source=None, target=None``
     - #TODO