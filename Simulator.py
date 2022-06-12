def print_basic_parameters():
	"""
	This methods prints the basic parameters of the epidemiological model.
	"""

def print_populations():
	"""
	This methods prints parameters of the population model.
	"""

def print_immunity_model():
	"""
	This methods prints the basic parameters of the immunity model.
	"""

def print_all(basic_parameters=True, populations=True, immunity_model=True):
	"""
	This methods prints all the parameters of the simulation.

	:param basic_parameters: whether to output the basic parameters of the epidemiological model.
	:type basic_parameters: bool

	:param populations: whether to output parameters of the population model.
	:type populations: bool

	:param immunity_model: whether to output parameters of the immunity model.
	:type immunity_model: bool
	"""

def set_transmission_rate(rate, haplotype=None):
	"""
	.. _Wikipedia: https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

	Transmission rate is the the expected number of new infections from a single infected individual per time unit in the beginning of the epidemics when all but one of the hosts are susceptible. See `Wikipedia`_ - that is parameter beta in SIR model.

	:param rate: transmission rate value.
	:type rate: float

	:param haplotype: haplotypes for which the new value is being set. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
	:type haplotype: int or string or None
	"""


def set_recovery_rate(rate, haplotype=None):
	"""
	.. _Wikipedia: https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

	Recovery rate is the inverse of the expected recovery time - non-negative float value. See `Wikipedia`_ - that is parameter gamma in SIR model.

	:param rate: recovery rate value.
	:type rate: float

	:param haplotype: haplotypes for which the new value is being set. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
	:type haplotype: int or string or None
	"""

def set_sampling_rate(rate, haplotype=None):
	"""
	Sampling rate: the rate at which infected individuals are sampled. The genealogy will be generated only for these samples. Alternatively, if *sampling_probability* is **True**, it will specify the fraction of recovered individuals which are sampled.

	:param rate: sampling rate value.
	:type rate: float

	:param haplotype: haplotypes for which the new value is being set. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
	:type haplotype: int or string or None
	"""

def set_mutation_rate(rate=None, probabilities=None, haplotype=None, mutation=None):
	"""
	This method allows setting the mutation rate and the weights for each single nucleotide substitution (given the mutation happened).

	:param rate: mutation rate value (None would not change the old value).
	:type rate: float or None

	:param probabilities: weights of each single nucleotide substitution given mutation occured (None would not change the old values). `See for example <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
	:type probabilities: list of four non-negative integers or None

	:param haplotype: the ancestral haplotype for the mutation. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
	:type haplotype: int or string or None

	:param mutation: the id of position on the haplotype where mutation arises.
	:type mutation: int or None
	"""

def set_susceptibility_type(susceptibility_type, haplotype=None):
	"""
	The type of immunity (or the susceptibility group) which an individual gets after being infected with a pathogen of a haplotype.

	:param susceptibility_type: immunity group id.
	:type susceptibility_type: int

	:param haplotype: haplotypes for which the new value is being set. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
	:type haplotype: int or string or None
	"""

def set_susceptibility(rate, haplotype=None, susceptibility_type=None):
	"""
	Susceptibility is a multiplicative modifier of the transmission rate based on the susceptible host immunity with a certain *susceptibility_type*. It can decrease or increase the transmission rate of a particular haplotype to the individuals of said *susceptibility_type*.

	:param rate: susceptibility value.
	:type rate: float

	:param haplotype: haplotypes for which the new value is being set. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
	:type haplotype: int or string or None

	:param susceptibility_type: immunity group id for which the new susceptibility value is being set.
	:type susceptibility_type: int or None
	"""

def set_immunity_transition(rate, source=None, target=None):
	"""
	The change of immunity without infection (e.g. due to vaccination or immunity loss with time).

	:param rate: transition rate value.
	:type rate: float

	:param source: source immunity group id (None means that the new value will be set to all immunity groups as source).
	:type source: int or None

	:param target: target immunity group id (None means that the new value will be set to all immunity groups as source).
	:type target: int or None
	"""


def set_population_size(size, population=None):
	"""
	Set the number of individuals in the population.

	:param size: total number of individuals in the population.
	:type size: int

	:param population: population for which the new population size is being set (in case of None the value will be updated for all populations).
	:type population: int or None
	"""

def set_contact_density(value, population=None):
	"""
	The relative number of contacts per time units.

	:param value: contact density value.
	:type value: float

	:param population: population for which the new contact density size is being set (in case of None the value will be updated for all populations).
	:type population: int or None
	"""

def set_npi(parameters, population=None):
	"""
	Setting conditions when lockdown in a population is imposed and lifted with the certain contact density during the lockdown.

	:param parameters: list with three elements: contact density value during lockdown, fraction of infectious population when the lockdown is set, fraction of infectious population when lockdown is lifted.
	:type parameters: list of length 3 with float values

	:param population: population for which the lockdown parameters are being set (in case of None the value will be updated for all populations).
	:type population: int or None
	"""

def set_sampling_multiplier(multiplier, population=None):
	"""
	The relative sampling in the population (multiplicative modifier). Sampling rate of each haplotype is modified by this factor.

	:param value: sampling multiplier value.
	:type value: float

	:param population: population for which the new sampling multiplier is being set (in case of None the value will be updated for all populations).
	:type population: int or None
	"""

def set_migration_probability(probability=None, total_probability=None, source=None, target=None):
	"""
	The probability that an individual from the population source is travelling to the population target.

	:param probability: probability of migration value.
	:type probability: float

	:param total_probability: if True, all the entries (except for the diagonal element) of row corresponding to source population is filled with probabilitie/(K-1). That represents individual travelling out of its population to a random destination.
	:type total_probability: bool

	:param source: source population with infectious individual (None means that the new value will be set to all populations as source).
	:type source: int or None

	:param target: target population with infectious individual (None means that the new value will be set to all populations as source).
	:type target: int or None
	"""

def set_susceptible_individuals(amount, source_type, target_type, population=None):
	"""
	The method allows moving susceptible individuals between susceptibility groups.
	
	:param amount: the number of individuals to be moved.
	:type amount: int

	:param source_type: source susceptibility group from which individuals are moved.
	:type source_type: int

	:param target_type: target susceptibility group into which individuals are moved.
	:type target_type: int
	
	:param population: population for which the new value is being set (in case of None the value will be updated for all populations).
	:type population: int or None
	"""

def set_infected_individuals(amount, source_haplotype, target_haplotype, population=None):
	"""
	The method allows moving infected individuals between susceptibility groups.

	:param amount: the number of individuals to be moved.
	:type amount: int

	:param source_haplotype: source haplotype from which individuals are moved.
	:type source_haplotype: int

	:param target_haplotype: target haplotype into which individuals are moved.
	:type target_haplotype: int

	:param population: population for which the new value is being set (in case of None the value will be updated for all populations).
	:type population: int or None
	"""

def set_infection(amount, source_type, target_haplotype, population=None):
	"""
	The method allows moving individuals from susceptibility group to target haplotype.

	:param amount: the number of individuals to be moved.
	:type amount: int

	:param source_type: source susceptibility group from which individuals are moved.
	:type source_type: int

	:param target_haplotype: target haplotype into which individuals are moved.
	:type target_haplotype: int

	:param population: population for which the new value is being set (in case of None the value will be updated for all populations).
	:type population: int or None
	"""

def set_chain_events(file_name):
	"""
	Allows to import chain of events directly from the .npy file

	:param file_name: a name of the file (without the .npy extension)
	:type file_name: str
	"""

def set_settings(file_template):
	"""
	Currently is a placeholder and does nothing. Do not include in actual docs until it does something.
	"""

def set_state(set_state):
	"""
	Calls set_chain_events and set_setting. Now identical to set_events. Do not include in docs until it starts to differ.
	"""


def output_newick(file_template="newick_output", file_path = None):
	"""
	Outputs the result of simulation in the newick format.

    :param file_template: template for the file name
    :type file_template: str

    :param file_path: path to output file
    :type file_path: str
	"""

def output_mutations(file_template="mutation_output", file_path = None):
	"""
	Outputs the information about all mutations

    :param file_template: template for the file name
    :type file_template: str

    :param file_path: path to output file
    :type file_path: str
	"""

def output_migrations(file_template="migrations", file_path = None):
	"""
	Outputs the information about all migrations

	:param file_template: template for the file name
	:type file_template: str

	:param file_path: path to output file
	:type file_path: str
	"""

def output_sample_data(output_print=False):
	"""
	Show all information about sampling, time, place, etc.

	:param output_print: if *Thue*, it outputs data about each sample to the console. Data includes time of sampling, population and haplotype of sampled individual.
	:type output_print: bool

	:return: 3 lists with times, population ids and haplotypes of the sampled individuals
	"""


def output_epidemiology_timelines(step=1000, output_file=False):
	"""
	Records simulation state changes over some period of time.

	:param step: a number of steps which split the epidemiology timelines.
	:type step: str

	:param output_file: if **True**, the timeline is written in the output file
	:type output_file: bool
	"""

def output_chain_events(file_name="chain_events"):
	"""
	Outputs event chain.

	:param file_name: file name for the output
	:param file_name: str
	"""

def output_settings(self, file_template="parameters"):
	"""
	Outputs settings.

	:param file_template: file name for the output
	:param file_template: str
	"""

def output_state():
	"""
	Outputs event chain and settings.
	"""


def get_data_susceptible(population, susceptibility_type, step_num):
	"""
	Returns the list with the amount of sucseptible individuals over some period of time.

	:param population: population to retrieve data from
	:type population: int

	:param susceptibility_type: susceptibility type to retrieve data from
	:type susceptibility_type: int

	:param step_num: number of steps which split the epidemiology timeline.
	:type step_num: int

	:return: 3 lists - amounts of susceptible individuals on the time points, time points themselves at which we retrieve amounts of individuals, and lockdowns data for the current epidemiology
	"""

def get_data_infectious(self, population, haplotype, step_num):
	"""
	Returns the list with the amount of infectious individuals over some period of time.

	:param population: population to retrieve data from
	:type population: int

	:param haplotype: haplotype to retrieve data from
	:type haplotype: int

	:param step_num: number of steps which split the epidemiology timeline
	:type step_num: int

	:return: 3 lists - amounts of infectious individuals on the time points, time points themselves at which we retrieve amounts of individuals, and lockdowns data for the current epidemiology
	"""


def plot_infectious(population, haplotype, step_num, label_infectious, label_samples):
	"""
	Creates plots for infectious and sampled individuals.

	:param population: population id for the plot
	:type population: int

	:param haplotype: haplotype for the plot
	:type haplotype: int or str

	:param step_num: number of steps which split the epidemiology timeline
	:type step_num: int

	:param label_infectious: the label for the infectious plot in the plot legend
	:type label_infectious: str

	:param label_samples: the label for the samples plot in the plot legend
	:type label_samples: str
	"""

def add_plot_infectious(population, haplotype, step_num=100, label_infectious=None, label_samples=None):
	"""
	Add to plot the trajectories of the change of the number of infectious individuals over time.

	:param population: population id for the plot
	:type population: int

	:param haplotype: haplotype for the plot
	:type haplotype: int or str

	:param step_num: number of steps which split the epidemiology timeline
	:type step_num: int

	:param label_infectious: the label for the infectious plot
	:type label_infectious: str or None

	:param label_samples: the label for the samples plot
	:type label_samples: str or None
	"""

def add_plot_susceptible(population, susceptibility_type, step_num=100, label=None):
	"""
	Add to plot the trajectories of the change of the number of susceptible individuals over time.

	:param population: population id for the plot
	:type population: int

	:param haplotype: haplotype for the plot
	:type haplotype: int or str

	:param step_num: number of steps which split the epidemiology timeline
	:type step_num: int

	:param label: the label for the susceptible plot
	:type label: str
	"""


def add_legend():
	"""
	Add legend to the plot.
	"""

def add_title(name="Plot"):
	"""
	Add plot title.

	:param name: the title for the plot
	:type name: str
	"""

def plot():
	"""
	Show plot with all the trajectories added by other methods.
	"""

def simulate(iterations=1000, sample_size=None, epidemic_time=None, method='direct'):
	"""
	This methods performs the simulation. The simulation interrupts when either one of the following conditions is satisfied: the number of iterations equals *iterations*, the number of collected samples equals *sample_size*, ot the total virtual time of the epidemic exceeds *time*
	It can be called multiple times, changes of most parameters are allowed between simulation runs.

    :param iterations: maximal number of iterations.
    :type iterations: int

    :param sample_size: desired sample size.
    :type sample_size: int or None

    :param epidemic_time: virtual (model) time to be simulated.
    :type epidemic_time: float or None

    :param method: 'direct' for the exact algorithm, 'tau' for tau-leaping approximated algorithm.
    :type method: string
	"""


def genealogy(seed=None):
	"""
	Generates a genealogy based on the chain of events generated during all the instances of simulate() method.

    :param seed: seed value (None for random seed).
    :type seed: int or None
	"""

def citation():
	"""
	Prints a citation of the paper in the console
	"""

def debug():
	"""
	Returns all of the available information about simulation. Mainly used for the debug purposes. Not documented.
	"""

def get_proportion():
	"""
	Used for paper preparation. Not documented.
	"""

def print_counters():
	"""
	Used for paper preparation. Not documented.
	"""

def print_propensities():
	"""
	Used for paper preparation. Not documented.
	"""
