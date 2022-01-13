class Simulator():
	"""
	#TODO

	:param: sites_number=0, populations_number=1, susceptibility_types=2, seed=None, sampling_probability=False, strong_migration=False
	:type: sites_number = int, populations_number = int, susceptibility_types = int, seed = float or None, sampling_probability = True or False, strong_migration = True or False
	"""
	def __init__(self):
		pass

def print_basic_parameters():
	"""
	#TODO
	"""

def print_populations():
	"""
	#TODO
	"""

def print_immunity_model():
	"""
	#TODO
	"""

def print_all(basic_parameters=True, populations=True, immunity_model=True):
	"""
	#TODO

	:param: basic_parameters=True, populations=True, immunity_model=True
	:type: basic_parameters = True or False, populations = True or False, immunity_model = True or False
	"""


def set_transmission_rate(rate, haplotype=None):
	"""
	.. _Wikipedia: https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

	Transmission rate: the expected number of new infections from a single infected individual per time unit in the beginning of the epidemics when (almost) all hosts are susceptible. See `Wikipedia`_ - parameter beta in SIR model.

    :param: rate, haplotype=None
    :type: rate = float, haplotype = int or str or None
	"""

def set_recovery_rate(rate, haplotype=None):
	"""
	.. _Wikipedia: https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

	Recovery rate: the inverse of the expected recovery time. `Wikipedia`_ - parameter gamma in SIR model.

    :param: rate, haplotype=None
    :type: rate = float, haplotype = int or str or None
	"""

def set_sampling_rate(rate, haplotype=None):
	"""
	Sampling rate: the rate at which infected individuals are sampled. The genealogy will be generated only for these samples. Alternatively, one can set probability=True in order to specify the fraction of recovered individuals which are sampled.

    :param: rate, haplotype=None
    :type: rate = float, haplotype = int or str or None
	"""

def set_mutation_rate(rate=None, substitution_weights=None, haplotype=None, site_id=None):
	"""
	This method allows setting the mutation rate and the weights for each single nucleotide substitution (given the mutation happened). haplotype is the ancestral haplotype for the mutation. site_id is the position on the haplotype where mutation arises. substitution_weights are given in the order of derived variant ATCG. The derived variant cannot be the same as the ancestral state, so the corresponding array entry will be ignored. 
	
    :param: rate=None, substitution_weights=None, haplotype=None, site_id=None
    :type: rate = float, substitution_weights = list of 4 elements, haplotype = int or str or None, site_id = int or None
	"""


def set_susceptibility_type(susceptibility_type, haplotype=None):
	"""
	The type of immunity (or the susceptibility group) which an individual gets after being infected with a pathogen of a haplotype.

	:param: susceptibility_type, haplotype=None
	:type: susceptibility_type = int, haplotype = int or str or None
	"""

def set_susceptibility(susceptibility, haplotype=None, susceptibility_type=None):
	"""
	Susceptibility is a multiplicative modifier of the transmission rate based on the susceptible host immunity with a susceptibility_type. It can decrease or increase the transmission rate of a particular haplotype to the individuals of susceptibility_type.

	:param: susceptibility, haplotype=None, susceptibility_type=None
	:type: susceptibility = float, haplotype = int or str or None, susceptibility_type = int or None
	"""

def set_immunity_transition(rate, source=None, target=None):
	"""
	The change of immunity without infection (e.g. due to vaccination or immunity loss with time).

	:param: rate, source=None, target=None
	:type: rate = float, source = int or None, target = int or None
	"""


def set_population_size(amount, population=None):
	"""
	The number of individuals in the population.

	:param: amount, population=None
	:type: amount = int, population = int or None
	"""

def set_contact_density(value, population=None):
	"""
	The relative number of contacts per time units.

	:param: value, population=None
	:type: value = float, population = int or None
	"""

def set_susceptible_individuals(amount, source_type, target_type, population=None):
	"""
	The method allows to move susceptible individuals between susceptibility groups. 

	:param: amount, source_type, target_type, population=None
	:type: amount = int, source_type = int, target_type = int, population = int or None
	"""

def set_lockdown(value, population=None):
	"""
	Setting conditions when lockdown in a population is imposed and lifted with the contact density during the lockdown.

	:param: value, population=None
	:type: valut = list, population = int or None
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

def add_legend():
	"""
	Add legend to the plot.
	"""

def add_title(name="Plot"):
	"""
	Plot title.

	:param: name="Plot"
	:type: name = str
	"""

def plot():
	"""
	Show plot with all the trajectories added by the previous two methods.
	"""


def simulate(iterations, sample_size, time):
	"""
	#TODO

	:param: iterations, sample_size, time
	:type: iterations = int, sample_size = int, time = float
	"""

def genealogy(seed=None):
	"""
	#TODO

	:param: seed=None
	:type: seed = float or None
	"""
