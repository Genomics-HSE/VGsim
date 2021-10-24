Tutorial
========

In this tutorial we will create a simulation scenario: we will choose parameters and set the model using the Python API.

Initialisation
--------------

We begin with choosing the number of mutable sites, where single nucleotide substitution can occur. We will have 2 such sites, leading to 16 different haplotypes might appear following mutations.

Next, we need to decide how many populations (demes) there are in the simulation, and choose the population sizes. Let’s create three populations with sizes 10000000, 5000000 and 1000000.

We will also create 3 types of susceptibility. Initially, all the individuals start with the absence of any specific immunity types.

The user can also specify a random seed.

.. code-block:: python

	import VGsim
	number_of_sites = 2
	population_sizes = [10000000, 5000000, 1000000]
	number_of_susceptible_groups = 3
	simulator = VGsim.Simulator(number_of_sites, population_sizes, number_of_susceptible_groups, seed=1234)

Choosing simulation parameters
------------------------------

**Transmission, recovery and sampling rates**

Assume that a patient becomes uninfectious within 10 days on average. Remember, that sampling an individual also means that this person becomes uninfectious. So, the sum of recovery and sampling rates is 1/10=0.1.

If 1 in 100 cases are sampled, then the recovery rate to sampling rate ratio is 99:1. Hence, recovery rate is 0.099 and sampling rate is 0.001.

Now, if in the absence of immune individuals in the population (e.g., in the beginning of the pandemic) an infectious person transmits the disease to 2.5 individuals on average, the birth rate is 2.5/10=0.25 (because if there are 2.5 transmissions in 10 days, hence there are 0.25 transmission in a single day).

.. code-block:: python

	simulator.set_infectious_rate(0.25)
	simulator.set_uninfectious_rate(0.099)
	simulator.set_sampling_rate(0.001)

Now let us assume that the original virus has haplotype AA. And the haplotype GG is two times more contagious. So, the transmission rate for that haplotype is 0.5.

.. code-block:: python

	simulator.set_infectious_rate(0.5, haplotype=”GG”)


**Mutation rates and substitution probabilities**

The GG haplotype can appear in the population only through mutations an both sites. For example, for SARS-CoV-2, there is roughly on average one new mutation per three transmissions (average time between transmissions are 1/0.25=4 days in our model) per approximately 30000bp genome. So, the per-site mutation rate is 1/(3*4*30000)≈0.000003. We also set the substitution probabilities. In fact, it might be more convenient sometimes to use unnormalised weights instead. Here we assume that the A->G substitution is twice more likely than A->T or A->C (the A->A substitution is not allowed, so the first entry of the array will be dropped for the corresponding haplotypes). Let’s also assume that having allele G at the first site increases the mutation rate at the second site by three times.

.. code-block:: python

	mutation_rate=0.000003
	substitution_weights=[1,1,1,2]#ATCG
	simulator.set_mutation_rate(mutation_rate, substitution_weights)
	simulator.set_mutation_rate(3*mutation_rate, XXX)

**Susceptibility types**

At the beginning of the epidemics all the individuals do not have any special immunity, and thus we assign them by default susceptibility type 0. Now let us assume that after recovery individuals get some immunity. If an individual was infected with haplotype carrying G at the first site (we will denote them by G*), the susceptibility type after recovery will be 2. For all other haplotypes it will be 1.

.. code-block:: python
	
	simulator.set_susceptibility_type(1)
	simulator.set_susceptibility_type(2, haplotype=”G*”)

Now for each susceptibility type we will specify how it changes the risk of infection. We will assume that susceptibility type 1 gives strong (but not perfect) protection against all haplotypes except G* by reducing the risk of a new infection by ten times. And it gives only some protection against haplotypes G*. The immunity is lost on average after 90 days, and the individual is hence moved to the susceptibility type 0.

.. code-block:: python

	simulator.set_susceptibility(0.1, susceptibility_type=1)
	simulator.set_susceptibility(0.5, susceptibility_type=1, haplotype=”G*”)
	simulator.set_immunity_transition(1/90, from_population=1, to_population=0)

Susceptibility type 2 will be a result of recovery from haplotype G* or vaccination. It gives complete protection from all the haplotypes, and it is lost after 180 days on average. At the beginning of the epidemics there is no vaccine, we will add it a bit later.

.. code-block:: python

	simulator.set_susceptibility(0.0, susceptibility_type=2)
	simulator.set_immunity_transition(1/180, from_population=2, to_population=0)

**Population model**

We have already set the population sizes. Now let us add some more heterogeneity. First of all assume that population 1 has a three times higher sampling rate than population 0, while population 2 does not sequence at all.

.. code-block:: python

	simulator.set_sampling_multiplier(3, population=1)
	simulator.set_sampling_multiplier(0, population=2)

Now, all the countries impose lockdowns when 20% of its population is infected simultaneously. The lockdowns are lifted if this number drops to 5%. The amount of contacts is ten times less during the lockdown.

.. code-block:: python

	simulator.set_lockdown([0.5, 0.2, 0.05])

**Migration**

There are 3 populations in our scenario. Assume that an average individual from one population spends on average 10 days per year outside of its population of origin. The chance to travel to any destination is the same. The cumulative migration probability is 10/365. The probability to travel to a particular destination is 10 (days)/365 (days)/2 (the number of possible destinations). So, we can set

.. code-block:: python

	simulator.set_migration_rate(10/365)