 .. py:currentmodule:: VGsim

API
===

.. autoclass:: Simulator

Epidemiological model
---------------------

.. autofunction:: VGsim.Simulator.set_transmission_rate
.. note::
	In many methods skipping an argument (or setting it to None) allows filling in multiple entries. Methods will apply the settings to all possible values of this argument.
.. autofunction:: VGsim.Simulator.set_recovery_rate
.. autofunction:: VGsim.Simulator.set_sampling_rate
.. autofunction:: VGsim.Simulator.set_mutation_rate
.. note::
	Haplotypes can be set either by their id (which can be calculated by assuming that the haplotype is a quaternary number with A=0, T=1, C=2, G=3), or directly by the string. E.g. haplotype “AGT” can be referred by the same string or by id=0*4\ :sup:`2`\+3*4\ :sup:`1`\+1*4\ :sup:`0`\=13. One can also use “*” to indicate “any nucleotide”. In this case the changes will be applied to all the haplotypes matching the pattern. E.g. “A*T” corresponds to four haplotypes “AAT”, “ATT”, “ACT”, “AGT”.

Host immunity
-------------

.. autofunction:: VGsim.Simulator.set_susceptibility_type
.. autofunction:: VGsim.Simulator.set_susceptibility
.. autofunction:: VGsim.Simulator.set_immunity_transition

Population model
----------------

.. autofunction:: VGsim.Simulator.set_population_size
.. autofunction:: VGsim.Simulator.set_contact_density
.. autofunction:: VGsim.Simulator.set_npi
.. autofunction:: VGsim.Simulator.set_sampling_multiplier
.. autofunction:: VGsim.Simulator.set_migration_probability
.. note::
	total_probability=float means that probability is the total probability to find an individual from population source outside of its population. All the entries (except for the diagonal) of the corresponding row of the migration probability matrix will be set to total_probability/(K-1), where K is the number of populations in the simulation.

.. Simulate
.. --------

.. .. autofunction:: functions.simulate
.. .. autofunction:: functions.genealogy

.. Output
.. ------

.. .. autofunction:: functions.output_newick
.. .. autofunction:: functions.output_mutations
.. .. autofunction:: functions.output_migrations
.. .. autofunction:: functions.sample_data
.. .. autofunction:: functions.epidemiology_timelines

.. Visualization
.. -------------

.. .. autofunction:: functions.add_plot_infectious
.. .. autofunction:: functions.add_plot_susceptible
.. .. autofunction:: functions.add_legend
.. .. autofunction:: functions.add_title
.. .. autofunction:: functions.plot

.. Parameters
.. ----------

.. .. autofunction:: functions.print_basic_parameters
.. .. autofunction:: functions.print_populations
.. .. autofunction:: functions.print_immunity_model
.. .. autofunction:: functions.print_all
