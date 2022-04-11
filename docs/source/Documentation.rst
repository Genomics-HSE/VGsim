API
===

.. autoclass:: functions.Simulator
   :members:

   .. automethod:: __init__

Epidemiological model
---------------------

.. autofunction:: functions.set_transmission_rate
.. note::
    In many methods skipping an argument (or setting it to None) allows filling in multiple entries. Methods will apply the settings to all possible values of this argument.
.. autofunction:: functions.set_recovery_rate
.. autofunction:: functions.set_sampling_rate
.. autofunction:: functions.set_mutation_rate
.. note::
    Haplotypes can be set either by their id (which can be calculated by assuming that the haplotype is a quaternary number with A=0, T=1, C=2, G=3), or directly by the string. E.g. haplotype “AGT” can be referred by the same string or by id=0*4^2+3*4^1+1*4*0=13. One can also use “*” to indicate “any nucleotide”. In this case the changes will be applied to all the haplotypes matching the pattern. E.g. “A*T” corresponds to four haplotypes “AAT”, “ATT”, “ACT”, “AGT”. :ref:`Haplotypes` page explains it in detail.

Host immunity
-------------

.. autofunction:: functions.set_susceptibility_type
.. autofunction:: functions.set_susceptibility
.. autofunction:: functions.set_immunity_transition


Population model
----------------

.. autofunction:: functions.set_population_size
.. autofunction:: functions.set_contact_density
.. autofunction:: functions.set_lockdown
.. autofunction:: functions.set_sampling_multiplier
.. autofunction:: functions.set_migration_probability
.. autofunction:: functions.set_susceptible_individuals
.. autofunction:: functions.set_infected_individuals
.. autofunction:: functions.set_infection

Simulate
--------

.. autofunction:: functions.simulate
.. autofunction:: functions.genealogy
.. autofunction:: functions.set_chain_events

Output
------

.. autofunction:: functions.output_newick
.. autofunction:: functions.output_mutations
.. autofunction:: functions.output_migrations
.. autofunction:: functions.output_sample_data
.. autofunction:: functions.output_epidemiology_timelines
.. autofunction:: functions.output_chain_events
.. autofunction:: functions.output_state

Visualization
-------------

.. autofunction:: functions.plot_infectious
.. autofunction:: functions.add_plot_infectious
.. autofunction:: functions.add_plot_susceptible
.. autofunction:: functions.add_legend
.. autofunction:: functions.add_title
.. autofunction:: functions.plot

Data retrieval
--------------

.. autofunction:: functions.get_data_susceptible
.. autofunction:: functions.get_data_infectious

Parameters
----------

.. autofunction:: functions.print_basic_parameters
.. autofunction:: functions.print_populations
.. autofunction:: functions.print_immunity_model
.. autofunction:: functions.print_all

Miscellaneous
-------------

.. autofunction:: functions.citation
