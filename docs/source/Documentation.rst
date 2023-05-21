API
===

.. autoclass:: VGsim.Simulator

Epidemiological model
---------------------

.. autofunction:: VGsim.Simulator.set_transmission_rate
.. note::
	In many methods skipping an argument (or setting it to None) allows filling in multiple entries. Methods will apply the settings to all possible values of this argument.
.. autofunction:: VGsim.Simulator.set_recovery_rate
.. autofunction:: VGsim.Simulator.set_sampling_rate
.. autofunction:: VGsim.Simulator.set_mutation_rate
.. autofunction:: VGsim.Simulator.set_mutation_probabilities
.. note::
	Haplotypes can be set either by their id (which can be calculated by assuming that the haplotype is a quaternary number with A=0, T=1, C=2, G=3), or directly by the string. E.g. haplotype “AGT” can be referred by the same string or by id=0*4\ :sup:`2`\+3*4\ :sup:`1`\+1*4\ :sup:`0`\=13. One can also use “*” to indicate “any nucleotide”. In this case the changes will be applied to all the haplotypes matching the pattern. E.g. “A*T” corresponds to four haplotypes “AAT”, “ATT”, “ACT”, “AGT”. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.

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

Simulate
--------

.. autofunction:: VGsim.Simulator.simulate
.. autofunction:: VGsim.Simulator.genealogy
.. autofunction:: VGsim.Simulator.set_chain_events

Output
------

.. autofunction:: VGsim.Simulator.export_newick
.. autofunction:: VGsim.Simulator.export_mutations
.. autofunction:: VGsim.Simulator.export_migrations
.. autofunction:: VGsim.Simulator.output_sample_data
.. autofunction:: VGsim.Simulator.output_epidemiology_timelines
.. autofunction:: VGsim.Simulator.export_chain_events
.. autofunction:: VGsim.Simulator.export_settings
.. autofunction:: VGsim.Simulator.export_state

Visualization
-------------

.. autofunction:: VGsim.Simulator.add_plot_infectious
.. autofunction:: VGsim.Simulator.add_plot_susceptible
.. autofunction:: VGsim.Simulator.add_legend
.. autofunction:: VGsim.Simulator.add_title
.. autofunction:: VGsim.Simulator.plot

Data retrieval
--------------

.. autofunction:: VGsim.Simulator.get_data_susceptible
.. autofunction:: VGsim.Simulator.get_data_infectious

Parameters
----------

.. autofunction:: VGsim.Simulator.print_basic_parameters
.. autofunction:: VGsim.Simulator.print_populations
.. autofunction:: VGsim.Simulator.print_immunity_model
.. autofunction:: VGsim.Simulator.print_all

Miscellaneous
-------------

.. autofunction:: VGsim.Simulator.citation