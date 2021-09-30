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

   * - Parameter
     - Type
     - Meaning
   * - ``infection``
     - ``int``
     - Frequency of incidence, how many people could become sick per period of time.
   * - ``uninfection``
     - ``int``
     - A number of people recovered from the disease per period of time
   * - ``Sr``
     - int
     - Sampling frequency, getting a virus dna through testing people to track any virus mutations
   * - ``Sp``
     - ``int``
     - second recovery model, percentage of people from a total sample of infected who recovered
   * - ``num_sites``
     - ``int``
     - Amount of parts of the genome where mutations can occur
   * - ``mut_rate``
     - ``int``
     - Frequency rate of certain mutations
   * - ``mut_target_rate``
     - ``int``
     - Frequency of transition from a particular nucleic acid to the other one  (e.g. Adenine to Thymine)
   * - **Population parameters**
     -
     -
   * - ``num_pop``
     - ``int``
     - Population scope (e.g. five countries)
   * - ``size_pop``
     - ``int``
     - Population size
   * - ``contact_density``
     - ``int``
     - Frequency of contacts between people in one population.
   * - ``contact_density``
     - ``int``
     - Frequency of contacts between people in one population.
   * - ``total_mig_rate``
     - ``int``
     - Migration rate (e.g. probability of movement from one population to another one)
   * - ``lockdown``
     - ``int``
     - Set of 3 criterias (threshold of incidents which leads to lockdown, lifting of lockdown, updated contact_density)
   * - ``sampling_multiplier``
     - ``int``
     - Virus sampling is multiplied by sampling rate to get the final amount for the population.
   * - **Susceptibility**
     -
     -
   * - ``susc_type``
     - ``int``
     - Post recovery type of individual
   * - ``susceptible``
     - ``int``
     - Amount of mutations per each haplotype
   * - ``susc_trans``
     - ``int``
     - Susceptibility matrix, the way susceptibility transitions




Show initial parameters
-----------------------

For example, this is the introductory paragraph
::
    print_Rates()
    print_Pop()
    print_Mig()
    print_Susc()
    print_SuscTrans()


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
     - int
     - time-based halt
   * - ``_seed``
     - ``int``
     - parameter for recovering the simulation

--

displays all simulator parameters
::
    debug(self)


log_dynamics(self, step=1000): - records simulation state changes over some period of time.
step - a number of parts log_dynamics is split on.

Graph
-----

**Graph for particular population and particular haplotype**

.. list-table::
   :widths: 15 15 70
   :header-rows: 1

   * - Parameter
     - Type
     - Meaning
   * - ``pop``
     - ``int``
     - Population
   * - ``hap``
     - ``int``
     - haplotype
   * - ``step_num``
     - int
     - A number of time intervals
   * - ``option``
     - ``int``
     - A graph type


System methods
::
   newick(): - record format of binary trees
   mut(): - information about all mutations
   mig(): - information about all migrations
   sampleDate(): - show all information about sampling, time, place, etc.

Change values
-------------

**target - current haplotype value, value - new haplotype value**

.. list-table::
   :widths: 15 70
   :header-rows: 1

   * - Parameter
     - Meaning
   * - ``set_Infection``
     -
   * - ``set_Uninfection``
     -
   * - ``set_S``
     -
   * - ``set_M``
     - frequency of particular mutation in particular haplotype
   * - ``set_MutRate``
     - frequency of transition between sets of mutations
   * - ``set_Migration``
     - displays shift from one migration to another one
   * - ``set_startLD``
     - initial population
   * - ``set_endLD``
     - final population
   * - ``set_conDenAfterLD``
     -
   * - ``set_suscType``
     -  susceptibility type for a particular haplotype
   * - ``set_susceptible``
     -  haplotype to be changed
   * - ``set_suscTrans``
     -   matrix of susceptibility type
