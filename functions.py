class Simulator():
	"""
	This is the most important class which creates and proceeds the simulation.
	"""

	def __init__(self, number_of_sites=0, populations_number=1, number_of_susceptible_groups=1, seed=None, sampling_probability=False, memory_optimization=None):
		"""

		:param number_of_sites: the number of mutable sites
		:type number_of_sites: int

		:param populations_number: the number of populations (demes)
		:type populations_number: int

		:param number_of_susceptible_groups: the number of susceptible groups (groups with different immunity response)
		:type number_of_susceptible_groups: int

		:param seed: seed to generate simulation from. If **None**, then chosen at random
		:type seed: float or None

		:param sampling_probability: whether we set sampling probability as a share of recovered individuals (**True** value) or we will set it explicitly. Default is **False**.
		:type sampling_probability: bool or None

		:param memory_optimization: if True, then memory optimization is conducted (useful for large number of possible haplotypes)
		:type memory_optimization: bool or None
		"""

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

def print_all(self, basic_parameters=True, populations=True, immunity_model=True):
	"""
	This methods prints all the parameters of the simulation.

	:param basic_parameters: whether to output the basic parameters of the epidemiological model.
	:type basic_parameters: bool

	:param populations: whether to output parameters of the population model.
	:type populations: bool

	:param immunity_model: whether to output parameters of the immunity model.
	:type immunity_model: bool
	"""

def set_transmission_rate(self, rate, haplotype=None):
	"""
	.. _Wikipedia: https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

	Transmission rate is the the expected number of new infections from a single infected individual per time unit in the beginning of the epidemics when all but one of the hosts are susceptible.
See `Wikipedia`_ - that is parameter beta in SIR model.

	:param rate: transmission rate value.
	:type rate: float

	:param haplotype: haplotypes for which the new value is being set. See :ref:`Haplotypes` for details.
	:type haplotype: int or string or None
	"""


def set_recovery_rate(self, rate, haplotype=None):
	"""
	.. _Wikipedia: https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

	Recovery rate is the inverse of the expected recovery time - non-negative float value. See `Wikipedia`_ - that is parameter gamma in SIR model.

	:param rate: recovery rate value.
	:type rate: float

	:param haplotype: haplotypes for which the new value is being set. See :ref:`Haplotypes` for details.
	:type haplotype: int or string or None
	"""

def set_sampling_rate(self, rate, haplotype=None):
	"""
	Sampling rate: the rate at which infected individuals are sampled. The genealogy will be generated only for these samples. Alternatively, if *sampling_probability* is **True**, it will specify the fraction of recovered individuals which are sampled.

	:param rate: sampling rate value.
	:type rate: float

	:param haplotype: haplotypes for which the new value is being set. See :ref:`Haplotypes` for details.
	:type haplotype: int or string or None
	"""

def set_mutation_rate(self, rate=None, probabilities=None, haplotype=None, mutation=None):
	"""
	This method allows setting the mutation rate and the weights for each single nucleotide substitution (given the mutation happened).

	:param rate: mutation rate value (None would not change the old value).
	:type rate: float or None

	:param probabilities: weights of each single nucleotide substitution given mutation occured (None would not change the old values). See :ref:`Haplotypes` for example.
	:type probabilities: list of four non-negative integers or None

	:param haplotype: the ancestral haplotype for the mutation. See :ref:`Haplotypes` for details.
	:type haplotype: int or string or None

	:param mutation: the id of position on the haplotype where mutation arises.
	:type mutation: int or None
	"""

def set_susceptibility_type(self, susceptibility_type, haplotype=None):
	"""
	The type of immunity (or the susceptibility group) which an individual gets after being infected with a pathogen of a haplotype.

	:param susceptibility_type: immunity group id.
	:type susceptibility_type: int

	:param haplotype: haplotypes for which the new value is being set. See :ref:`Haplotypes` for details.
	:type haplotype: int or string or None
	"""

def set_susceptibility(self, rate, haplotype=None, susceptibility_type=None):
	"""
	Susceptibility is a multiplicative modifier of the transmission rate based on the susceptible host immunity with a certain *susceptibility_type*. It can decrease or increase the transmission rate of a particular haplotype to the individuals of said *susceptibility_type*.

	:param rate: susceptibility value.
	:type rate: float

	:param haplotype: haplotypes for which the new value is being set. See :ref:`Haplotypes` for details.
	:type haplotype: int or string or None

	:param susceptibility_type: immunity group id for which the new susceptibility value is being set.
	:type susceptibility_type: int or None
	"""

def set_immunity_transition(self, rate, source=None, target=None):
	"""
	The change of immunity without infection (e.g. due to vaccination or immunity loss with time).

	:param rate: transition rate value.
	:type rate: float

	:param source: source immunity group id (None means that the new value will be set to all immunity groups as source).
	:type source: int or None

	:param target: target immunity group id (None means that the new value will be set to all immunity groups as source).
	:type target: int or None
	"""


def set_population_size(self, size, population=None):
	"""
	Set the number of individuals in the population.

	:param size: total number of individuals in the population.
	:type size: int

	:param population: population for which the new population size is being set (in case of None the value will be updated for all populations).
	:type population: int or None
	"""

def set_contact_density(self, value, population=None):
	"""
	The relative number of contacts per time units.

	:param value: contact density value.
	:type value: float

	:param population: population for which the new contact density size is being set (in case of None the value will be updated for all populations).
	:type population: int or None
	"""

def set_lockdown(self, parameters, population=None):
	"""
	Setting conditions when lockdown in a population is imposed and lifted with the certain contact density during the lockdown.

	:param parameters: list with three elements: contact density value during lockdown, fraction of infectious population when the lockdown is set, fraction of infectious population when lockdown is lifted.
	:type parameters: list of length 3 with float values

	:param population: population for which the lockdown parameters are being set (in case of None the value will be updated for all populations).
	:type population: int or None
	"""

def set_sampling_multiplier(self, multiplier, population=None):
	"""
	The relative sampling in the population (multiplicative modifier). Sampling rate of each haplotype is modified by this factor.

	:param value: sampling multiplier value.
	:type value: float

	:param population: population for which the new sampling multiplier is being set (in case of None the value will be updated for all populations).
	:type population: int or None
	"""

def set_migration_probability(self, probability=None, total_probability=None, source=None, target=None):
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
	The method allows to move susceptible individuals between susceptibility groups.

    :param amount: the number of individuals to be moved.
    :param source_type: source susceptibility group from which individuals are moved.
    :param target_type: target susceptibility group into which individuals are moved.
    :param population: population for which the new value is being set (in case of None the value will be updated for all populations).
    :type amount: int
    :type source_type: int or None
    :type target_type: float
    :type population: int or None
	"""


def simulate(iterations=1000, sample_size=None, time=None, method='direct'):
	"""
	This methods starts the simulation. The simulation interrupts when either one of the conditions is satisfied: the number of iterations is iterations, the number of collected samples is sample_size, ot the total virtual time of the epidemic exceeds time
	It can be called multiple times, changes of most parameters are allowed between simulation runs.

    :param iterations: maximal number of iterations.
    :param sample_size: desired sample size.
    :param time: virtual (model) time to be simulated.
    :param method: 'direct' for the exact algorithm, 'tau' for tau-leaping approximated algorithm.
    :type iterations: int
    :type sample_size: int or None
    :type time: float or None
    :type method: string
	"""

def genealogy(seed=None):
	"""
	Generating a genealogy based on the chain of events generated during all the instances of simulate() method.

    :param seed: seed value (None for random seed).
    :type seed: int or None
	"""

def output_newick(file_template="newick_output", file_path = ''):
	"""
	Record format of binary trees

    :param file_template: template for file name
    :param file_path: path to output file
    :type file_template: str
    :type file_path: str
	"""

def output_mutations(file_template="mutation_output", file_path = ''):
	"""
	Information about all mutations

    :param file_template: template for file name
    :param file_path: path to output file
    :type file_template: str
    :type file_path: str
	"""

def output_migrations(file_template="migrations", file_path = ''):
	"""
	Information about all migrations

    :param file_template: template for file name
    :param file_path: path to output file
    :type file_template: str
    :type file_path: str
	"""

def output_parameters(name_file="parameters"):
	"""
	Make the directory with files for launching via console

	:param: name_file="parameters"
	:type: name_file = str
	"""

def sample_data():
	"""
	Show all information about sampling, time, place, etc.
	"""

def epidemiology_timelines(step=1000, output_file=False):
	"""
	Records simulation state changes over some period of time. step - a number of parts epidemiology_timelines is split on.
	"""


def add_plot_infectious(population, haplotype, step_num=100, label=None):
	"""
	Add to plot the trajectories of the change of the number of infectious individuals over time.  label allows to add the label to plot legend.

	:param: population, haplotype, step_num=100, label=None
	:type: population = int, haplotype = int or str, step_num = int, label = str or None
	"""

def add_plot_susceptible(population, susceptibility_type, step_num=100, label=None):
	"""
	Add to plot the trajectories of the change of the number of susceptible individuals over time.

	:param: population, susceptibility_type, step_num=100, label=None
	:type: population = int, susceptibility_type = int, step_num = int, label = str or None
	"""


def plot():
	"""
	Show plot with all the trajectories added by the previous two methods.
	"""

def add_title(name="Plot"):
	"""
	Plot title.

	:param: name="Plot"
	:type: name = str
	"""

def add_legend():
	"""
	Add legend to the plot.
	"""
