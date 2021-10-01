VGsim
======

Our model of epidemiological spread is built as a compartmental model, and the random realisations of the corresponding stochastic process are drawn using the Gillespie algorithm. The different compartments in our model are defined based on several real-world complexities that interact with each other in a complex way and can affect epidemiological dynamics: population structure (which means defining separate host populations and assigning different frequencies to within-population and between-population contacts), separate infectious groups (which means that host individuals carrying different viral haplotypes are modeled differently, since some viral variants might be more transmissible than others), and different susceptible groups (different hosts having different types of immunity response to different haplotypes).

Interface
=========

init
----

.. list-table::
   :widths: 15 15 70
   :header-rows: 1

   * - Basic parameters
     - Type
     - Meaning
   * - ``infectious_rate``
     - ``double``
     - Frequency of incidence, how many people could become sick per period of time
   * - ``uninfectious_rate``
     - ``double``
     - A number of people recovered from the disease per period of time
   * - ``sampling_rate``
     - ``double``
     - Sampling frequency, getting a virus dna through testing people to track any virus mutations
   * - ``sampling_proportion``
     - ``double``
     - Second recovery model, percentage of people from a total sample of infected who recovered
   * - ``sites_number``
     - ``int``
     - Amount of parts of the genome where mutations can occur
   * - ``mutation_rate``
     - ``double``
     - Frequency rate of certain mutations
   * - ``mutation_probabilities``
     - ``[double, double, double]``
     - Frequency of transition from a particular nucleic acid to the other one  (e.g. Adenine to Thymine). Example [2,3,5]
   * - **Population parameters**
     -
     -
   * - ``populations_number``
     - ``int``
     - Population scope (e.g. five countries)
   * - ``population_size``
     - ``int``
     - Population size
   * - ``contact_density``
     - ``double``
     - Frequency of contacts between people in one population.
   * - ``total_migration_probability``
     - ``double``
     - Migration rate (e.g. probability of movement from one population to another one)
   * - ``lockdown``
     - ``[double, int, int]``
     - Set of 3 criterias (threshold of incidents which leads to lockdown, lifting of lockdown, updated contact_density) Example [2.0, 5, 2]
   * - ``sampling_multiplier``
     - ``double``
     - Virus sampling is multiplied by sampling rate to get the final amount for the population.
   * - **Susceptibility**
     -
     -
   * - ``immunity_type``
     - ``int``
     - Post recovery type of individual
   * - ``susceptibility``
     - ``list``
     - #TODO
   * - ``total_immunity_transition``
     - ``double``
     - Susceptibility matrix, the way susceptibility transitions




Show initial parameters
-----------------------

For example, this is the introductory paragraph
::
    print_basic_rates()
    print_populations()
    print_migration_matrix()
    print_immunity_model()


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
   * - ``_iterations``
     - ``int``
     - Amount of iterations which will take place before the halt
   * - ``_sampleSize``
     - ``int``
     - Halt the simulation when a particular number of participants was sampled
   * - ``_time``
     - ``double``
     - time-based halt
   * - ``_seed``
     - ``double``
     - parameter for recovering the simulation

plot
----

**Graph for particular population and particular haplotype**

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
     - haplotype


System methods
::
   output_newick(name_file="newick_output"): - record format of binary trees
   output_mutations(name_file="mutation_output"): - information about all mutations
   output_migrations(name_file="migrations")(): - information about all migrations
   sampleDate(): - show all information about sampling, time, place, etc.
   log_dynamics(self, step=1000, output_file=False): - records simulation state changes over some period of time. step - a number of parts log_dynamics is split on.

Change values
-------------

.. list-table::
   :widths: 15 25 70
   :header-rows: 1

   * - Function
     - Parameters
     - Meaning
   * - ``set_infectious_rate``
     - ``haplotype, rate``
     -
   * - ``set_uninfectious_rate``
     - ``haplotype, rate``
     -
   * - ``set_sampling_rate``
     - ``haplotype, rate``
     -
   * - ``set_mutation_rate``
     - ``haplotype, site, rate=None, probabilities=None``
     - frequency of transition between sets of mutations
   * - ``set_migration_probability``
     - ``source_population, target_population, probability``
     - displays shift from one migration to another one
   * - ``set_start_lockdown``
     - ``population, infectious_fraction=None, contact_density=None``
     - initial population
   * - ``set_end_lockdown``
     - ``population, infectious_fraction``
     - final population
   * - ``set_immunity_type``
     - ``haplotype, immunity``
     -  susceptibility type for a particular haplotype
   * - ``set_susceptibility``
     - ``haplotype, immunity, susceptibility``
     -  haplotype to be changed
   * - ``set_immunity_transition``
     - ``source_immunity, target_immunity, probability``
     -   matrix of susceptibility type