API
===

.. autoclass:: sim.Simulator
   :members:

   .. automethod:: __init__

Epidemiological model
---------------------

.. autofunction:: Simulator.set_transmission_rate
.. note::
	In many methods skipping an argument (or setting it to None) allows filling in multiple entries. Methods will apply the settings to all possible values of this argument.
.. autofunction:: Simulator.set_recovery_rate
.. autofunction:: Simulator.set_sampling_rate
.. autofunction:: Simulator.set_mutation_rate
.. note::
	Haplotypes can be set either by their id (which can be calculated by assuming that the haplotype is a quaternary number with A=0, T=1, C=2, G=3), or directly by the string. E.g. haplotype “AGT” can be referred by the same string or by id=0*4\ :sup:`2`\+3*4\ :sup:`1`\+1*4\ :sup:`0`\=13. One can also use “*” to indicate “any nucleotide”. In this case the changes will be applied to all the haplotypes matching the pattern. E.g. “A*T” corresponds to four haplotypes “AAT”, “ATT”, “ACT”, “AGT”. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.

Host immunity
-------------

.. autofunction:: Simulator.set_susceptibility_type
.. autofunction:: Simulator.set_susceptibility
.. autofunction:: Simulator.set_immunity_transition

Population model
----------------

.. autofunction:: Simulator.set_population_size
.. autofunction:: Simulator.set_contact_density
.. autofunction:: Simulator.set_npi
.. autofunction:: Simulator.set_sampling_multiplier
.. autofunction:: Simulator.set_migration_probability
.. note::
	total_probability=float means that probability is the total probability to find an individual from population source outside of its population. All the entries (except for the diagonal) of the corresponding row of the migration probability matrix will be set to total_probability/(K-1), where K is the number of populations in the simulation.

Simulate
--------

.. autofunction:: Simulator.simulate
.. autofunction:: Simulator.genealogy
.. autofunction:: Simulator.set_chain_events

Output
------

.. autofunction:: Simulator.output_newick
.. autofunction:: Simulator.output_mutations
.. autofunction:: Simulator.output_migrations
.. autofunction:: Simulator.output_sample_data
.. autofunction:: Simulator.output_epidemiology_timelines
.. autofunction:: Simulator.output_chain_events
.. autofunction:: Simulator.output_state

Visualization
-------------

.. autofunction:: Simulator.plot_infectious
.. autofunction:: Simulator.add_plot_infectious
.. autofunction:: Simulator.add_plot_susceptible
.. autofunction:: Simulator.add_legend
.. autofunction:: Simulator.add_title
.. autofunction:: Simulator.plot

Data retrieval
--------------

.. autofunction:: Simulator.get_data_susceptible
.. autofunction:: Simulator.get_data_infectious

Parameters
----------

.. autofunction:: Simulator.print_basic_parameters
.. autofunction:: Simulator.print_populations
.. autofunction:: Simulator.print_immunity_model
.. autofunction:: Simulator.print_all

Miscellaneous
-------------

.. autofunction:: Simulator.citation