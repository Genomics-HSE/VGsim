Migration
---------

The migrations are modelled as between-population transmissions: an infected individual from the source population `s` contacts a susceptible individual from the recipient population `r`. So, the individuals should contact each other while being in the same population. That can happen in three ways: an infected individual travels to population `r`, susceptible individual travels to population `s`, or both individuals travel to a third population `q` simultaneously. These travels are not modelled directly, but instead there is a matrix of probabilities to find an individual from one population into another.

For example, if an individuals from population 0 spends on average 7 days per year in population 1 and 1 day per year in population 2, the corresponding matrix entries row will be 

.. code-block:: python

	Simulator.set_migration_probability(probability=7/365, source=0, target=1)
	Simulator.set_migration_probability(probability=1/365, source=0, target=2)

=========== =========== ===========
1-(7+1)/365 7/365       1/365
=========== =========== ===========

Assume that individuals from population 1 do not travel at all: 

= = =
0 1 0
= = =

And travelers from the last population 2 travel to population 0 and 1 symmetrically: 

.. code-block:: python

	Simualator.set_migration_probability(probability=1/365, source=2, target=0)

======= ======= =======
1/365   0       1-1/365
======= ======= =======

The resulting migration matrix is

===== ===== =====
0.978 0.019 0.003
0     1     0
0.003 0     0.997
===== ===== =====

Notice, that travels might considerably change the number of individuals which present in a deme at a given time. For example, assume that population sizes are `N0=1,000,000,000`, `N1=100,000,000` and `N2=1,000,000`, respectively. Then the expected number of individuals in each deme is

* D0=357/365*10\ :sup:`9`\+1/365*10\ :sup:`6`\=978,084,932 (98% of N0)
* D1=7/365*10\ :sup:`9`\+10\ :sup:`8`\=119,178,082 (119% of N1)
* D2=1/365*10\ :sup:`9`\+364/365*10\ :sup:`6`\=3,736,986 (373% of N2)

So, when the population sizes are considerably different, even small migration probabilities might lead to dramatic effects of the number of individuals in a deme. It might be unrealistic, and in most cases will probably be a sign to verify the simulation parameters.

One can check the number of individuals in each deme and their changes relative to corresponding population sizes by the following command print_populations() and pay attention to the `Size` and `Actual size` collums.