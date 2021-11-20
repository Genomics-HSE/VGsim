Command line
============

Stopping criterion
------------------

Exists two variants to stop simulation.
First `-it` or `--iterations` defined lenght of the chain of events. Second `-s` or `--sampleSize` additionally defined stopping criterion by the number of sample. If this parameter is not defined, then the default value is equal to the number of iterations.


Setting haplotype (strain) model
--------------------------------

You should provide a space-separated file with the per-haplotype rates (all the rates are assumed to be `float` non-negative numbers). Example:

.. code-block:: python

	#Rates_format_version 0.0.1
	H B D S M1 M2
	AA 25 9 1 0.1 0.1
	AT 25 9 1 0.1,3,4,1 0.1,0.4,0.6,0.0
	AC 25 9 1 0.1 0.1

Column `H` (optional) specifies the haplotype. In the example file there are two sites, so the total of 16 haplotypes should appear in the file in lexicographical order (we assume the following nucleotide alphabetic order: `ATCG`).

Column `B` (for Birth rate) is the base rate for an individual carrying the corresponding strain to transmit the infection to a new individual.

Column `D` (for Death rate) is the rate for an individual to get uninfectious.

Column `S` is the rate for an individual to be sampled. Sampling also means that the individual got isolated and treated, hence become uninfectious immediately.

`M1` and `M2` are two sites where non-neutral mutations can be ovserved. The user can specify mutation rates only (see haplotype `AA`), and the substitution rates are assumed to be unifrom. Also weights (probabilities, but not necessariliy normalised) can be specified for each substitution (see haplotype `AT` - in the first site susbstitions `A->T`, `A->C`, `A->G` occur with probabilities 3/8, 4/8 and 1/8 respectively, and for the second site substitutions `T->A`, `T->C`, `T->G` occur with probabilities 0.4, 0.6 and 0.0 respectively).

Setting population model
------------------------

In order to set population model use `--populationModel` or `-pm` flag followed by two files. The first file contains information abouch each population, and the second file contains the matrix of "visit" rates.

**Part 1 - setting populations**

Here is an example of the file with populations.

.. code-block:: python

	#Population_format_version 0.0.1
	id size contactDensity conDenAfterLD startLD endLD samplingMultiplier
	0 20000000 1.0 0.1,2,1 5 
	1 10000000 1.0 0.1,2,1 1

`id` is the population number (for convinience). `size` is the total number of individuals in the population. Contact density in the relative number of contacts per time unit. It can be used to model different social, cultural, economical and other aspects (e.g. population density in a city or a country side, holiday times etc.) Currently we provide only one out-of-the-box solution to change contact density during simulation. Those are optinal fields (included in the example) to tune lockdowns. `conDenAfterLD` is the contact density during lockdown, `startLD` and `endLD` are condition to impose and lift the lockdown. These two numbers are the percentage of individuals which are simultaneously infected. If the percentage of simultanelously infected individuals in a population becomes larger than `startLD`, lockdown is imposed. As soon as the percentage of simultaneously infected individuals drops below `endLD` the lockdown is lifted.

**Part 2 - setting migration matrix**

In the model we consider migration as some individuals spending some time in other population and then returning to their home population. This is modelled by the probabilities that an individual from their home population is found in a destination population. During such visits individuals are assumed to contact with individuals of another population according to contact density of the destination population (so lockdown in a destination country also changes the contact rate of a visitor). Example:

.. code-block:: python

	#Migration_format_version 0.0.1
	0.0 0.002 0.002 0.001
	0.0 0.0 0.002 0.0
	0.0 0.005 0.0 0.001
	0.0 0.002 0.001 0.0

Susceptibility types (immunity)
-------------------------------

Multiple susceptibility classes can be modelled by adding `--susceptibility` or, `-su` flag followed by the file with susceptibility constants (non-negative `float`). Example:

.. code-block:: python

	#Susceptibility_format_version 0.0.1
	H T S0 S1 S2
	AA 1 1.0 0.0 0.0
	AT 2 1.0 0.4 0.0
	AC 2 1.0 1.1 0.0

There are three susceptibility types `S0`, `S1` and `S2` in this example. All the individuals start in `S0`, and they have susceptibility of `1.0` to all strains. `T` is the type of susceptibility (immunity) caused by recovering from a particular strain. Inficting by haplotype `AA` leads to susceptibility type `S1`, which gives total resistance (susceptibility 0.0) to haplotype `AA`, partial resitance to `AT` (susceptibility 0.4) and increases susceptibility (1.1) to strain `AC`.

*NB* There is no "immunity memory" - the immunity does not depend on the whole illness history of an individual, but only on the **latest** infection.

**Susceptibility transition**

The user can specify the rates of direct transitions between susceptibility types. This can be used for example to model vaccination or immunity loss. Use `--suscepTransition` or, `-st` flag followed by the file with susceptibility transition rate matrix (non-negative `float` entrys). Example:

.. code-block:: python

	#Susceptibility_format_version 0.0.1
	0.0 0.0 0.0001
	0.001 0.0 0.0001
	0.0 0.0 0.0

Flags
-----

.. list-table::
   :widths: 15 15 70
   :header-rows: 1

   * - Flag
     - Parameter
     - Meaning
   * - -it or --iterations
     - int
     - number of iterations (default is 1000)
   * - -s or --sampleSize
     - int
     - number of sample (default is None)
   * - -t or --time
     - float
     - time for stopping simulation (default is None)
   * - -seed or --seed
     - float
     - random seed
   * - -pm or --populationModel
     - path to two file
     - population model: a file with population sizes etc, and a file with migration rate matrix
   * - -su or --susceptibility
     - path to file
     - susceptibility file
   * - -st or --suscepTransition
     - path to file
     - susceptibility transition file
   * - --sampling_probability
     - None
     - #TODO
   * - --strong_migration
     - None
     - #TODO
   * - -nwk or --createNewick
     - None
     - Create a newick file of tree *.nwk
   * - -tsv or --writeMutations
     - None
     - Create a mutation file *.tsv
   * - --writeMigrations
     - None
     - Create a migration file *.txt
   * - -c or -citation
     - None
     - Information for citation

Output
------

The final genealogical tree can be exported into Newick format with `--createNewick` or `-nwk`, and mutations can be exported in a tsv file with `--writeMutations` or `-tsv`. This format is compatible with phastSim (https://github.com/NicolaDM/phastSim) and Usher (REF).
