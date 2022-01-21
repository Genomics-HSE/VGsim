


class Simulator():
	"""
	This is the class which creates the simulation.
		sites_number: the number of mutable sites with strong phenotypic effect
		populations_number: the number of populations (demes)
		susceptibility_types: the number of susceptible groups (groups with different immunity response)

	:param: sites_number=0, populations_number=1, susceptibility_types=2, seed=None, sampling_probability=False, strong_migration=False
	:type: sites_number = int, populations_number = int, susceptibility_types = int, seed = float or None, sampling_probability = True or False, strong_migration = True or False
	"""
	def __init__(self):
		pass

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

	:param: basic_parameters=True, populations=True, immunity_model=True
	:type: basic_parameters = True or False, populations = True or False, immunity_model = True or False
	"""

def set_transmission_rate(rate, haplotype=None):
	"""
	.. _Wikipedia: https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

	Transmission rate: the expected number of new infections from a single infected individual per time unit in the beginning of the epidemics when (almost) all hosts are susceptible. See `Wikipedia`_ - parameter beta in SIR model.

    :param rate: transmission rate value.
    :param haplotype: haplotypes for which the new value is being set. See XXX for details.
    :type rate: float
    :type haplotype: int or string or None
	"""


def set_recovery_rate(rate, haplotype=None):
	"""
	.. _Wikipedia: https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

	Recovery rate: the inverse of the expected recovery time. `Wikipedia`_ - parameter gamma in SIR model.

	:param rate: recovery rate value.
    :param haplotype: haplotypes for which the new value is being set. See XXX for details.
    :type rate: float
    :type haplotype: int or string or None
	"""

def set_sampling_rate(rate, haplotype=None):
	"""
	Sampling rate: the rate at which infected individuals are sampled. The genealogy will be generated only for these samples. Alternatively, one can set probability=True in order to specify the fraction of recovered individuals which are sampled.

	:param rate: sampling rate value.
    :param haplotype: haplotypes for which the new value is being set. See XXX for details.
    :type rate: float
    :type haplotype: int or string or None
	"""

def set_mutation_rate(rate=None, substitution_weights=None, haplotype=None, site_id=None):
	"""
	This method allows setting the mutation rate and the weights for each single nucleotide substitution (given the mutation happened). haplotype is the ancestral haplotype for the mutation. site_id is the position on the haplotype where mutation arises. substitution_weights are given in the order of derived variant ATCG. The derived variant cannot be the same as the ancestral state, so the corresponding array entry will be ignored.

	:param rate: mutation rate value (None would not change the old value).
	:param substitution_weights: weights of each single nucleotide substitution given mutation occured (None would not change the old values). See XXX for example.
    :param haplotype: haplotypes for which the new value is being set. See XXX for details.
    :type rate: float or None
	:type substitution_weights: list of four non-negative integers or None
    :type haplotype: int or string or None
	"""

def set_susceptibility_type(susceptibility_type, haplotype=None):
	"""
	The type of immunity (or the susceptibility group) which an individual gets after being infected with a pathogen of a haplotype.

	:param susceptibility_type: immunity group id.
    :param haplotype: haplotypes for which the new value is being set. See XXX for details.
    :type susceptibility_type: int
    :type haplotype: int or string or None
	"""

def set_susceptibility(susceptibility, haplotype=None, susceptibility_type=None):
	"""
	Susceptibility is a multiplicative modifier of the transmission rate based on the susceptible host immunity with a susceptibility_type. It can decrease or increase the transmission rate of a particular haplotype to the individuals of susceptibility_type.

	:param susceptibility: susceptibility value.
    :param haplotype: haplotypes for which the new susceptibility value is being set. See XXX for details.
	:param susceptibility_type: immunity group id for which the new susceptibility value is being set.
    :type susceptibility: int
    :type haplotype: int or string or None
	:type susceptibility_type: int

	:param: susceptibility, haplotype=None, susceptibility_type=None
	:type: susceptibility = float, haplotype = int or str or None, susceptibility_type = int or None
	"""

def set_immunity_transition(rate, source=None, target=None):
	"""
	The change of immunity without infection (e.g. due to vaccination or immunity loss with time).

	:param rate: transition rate value.
    :param source: source immunity group id (None means that the new value will be set to all immunity groups as source).
	:param target: target immunity group id (None means that the new value will be set to all immunity groups as target).
    :type rate: float
    :type source: int or None
	:type target: int or None

	:param: rate, source=None, target=None
	:type: rate = float, source = int or None, target = int or None
	"""


def set_population_size(amount, population=None):
	"""
	The number of individuals in the population.

	:param amount: total number of individuals in the population.
    :param population: population for which the new population size is set (in case of None the value will be updated for all populations).
    :type amount: int
    :type population: int or None
	"""

def set_contact_density(value, population=None):
	"""
	The relative number of contacts per time units.

	:param value: contact density value.
    :param population: population for which the new contact density size is being set (in case of None the value will be updated for all populations).
    :type value: float
    :type population: int or None
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

def set_lockdown(parameters, population=None):
	"""
	Setting conditions when lockdown in a population is imposed and lifted with the contact density during the lockdown.

	:param parameters: list with three elements: contact density value during lockdown, fraction of infectious population when the lockdown is set, fraction of infectious population when lockdown is lifted.
    :param population: population for which the new lockdown parameters are being set (in case of None the value will be updated for all populations).
    :type value: list of length 3
    :type population: int or None
	"""

def set_sampling_multiplier(value, population=None):
	"""
	The relative sampling in the population (multiplicative modifier). Sampling rate of each haplotype is modified by this factor.

	:param: value, population=None
	:type: value = float, population = int or None
	"""

def set_migration_probability(probability=None, total_probability=None, source=None, target=None):
	"""
	The probability that an individual from the population source is travelling to the population target. cumulative=True means that probability is the total probability to find an individual from population source outside of its population. All the entries (except for the diagonal) of the corresponding row of the migration probability matrix will be set to probability/(K-1), where K is the number of populations in the simulation.

	:param: probability=None, total_probability=None, source=None, target=None
	:type:	probability = float or None, total_probability = float or None, source = int or None, target = int or None
	"""


def simulate(iterations, sample_size, time):
	"""
	This methods starts the simulation. The simulation interrupts when either one of the conditions is satisfied: the number of iterations is iterations, the number of collected samples is sample_size, ot the total virtual time of the epidemic exceeds time
	It can be called multiple times, changes of most parameters are allowed between simulation runs.

	:param: iterations, sample_size, time
	:type: iterations = int, sample_size = int, time = float
	"""

def genealogy(seed=None):
	"""
	Generating a genealogy based on the chain of events generated during all the instances of simulate() method.

	:param: seed=None
	:type: seed = float or None
	"""

def output_newick(name_file="newick_output"):
	"""
	Record format of binary trees

	:param: name_file="newick_output"
	:type: name_file = str
	"""

def output_mutations(name_file="mutation_output"):
	"""
	Information about all mutations

	:param: name_file="mutation_output"
	:type: name_file = str
	"""

def output_migrations(name_file="migrations"):
	"""
	Information about all migrations

	:param: name_file="migrations"
	:type: name_file = str
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
