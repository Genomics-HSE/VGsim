from ._BirthDeath import BirthDeathModel
from .IO import writeGenomeNewick, writeMutations
from random import randrange
import sys
import matplotlib.pyplot as plt
import numpy as np
import time

class Simulator:
    """
    The class which creates and runs simulations.

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

    :param genome_length:
    :type genome_length: int

    :param recombination_probability:
    :type recombination_probability: float
    """
    def __init__(self, number_of_sites=0, populations_number=1, number_of_susceptible_groups=1, seed=None, \
        sampling_probability=False, memory_optimization=False, genome_length=int(1e6), recombination_probability=0.0, \
        number_of_states_allele=4):
        self.fig = None
        if seed == None:
            seed = int(randrange(sys.maxsize))
        print('User seed:', seed)

        self.simulation = BirthDeathModel(number_of_sites=number_of_sites, populations_number=populations_number, \
            number_of_susceptible_groups=number_of_susceptible_groups, seed=seed, sampling_probability=sampling_probability, \
            memory_optimization=memory_optimization, genome_length=genome_length, recombination_probability=recombination_probability, \
            number_of_states_allele=number_of_states_allele)

    def set_proportion_type_infections(self, proportion):
        self.simulation.set_proportion_type_infections(proportion)

    def print_basic_parameters(self):
        """
        This methods prints the basic parameters of the epidemiological model.
        """
        self.simulation.print_basic_parameters()

    def print_populations(self, population=True, susceptibles=True, infectious=True, migration=True):
        """
        This methods prints parameters of the population model.

        :param population: print informations about populations.
        :type population: bool

        :param susceptibles: print informations about susceptibles.
        :type susceptibles: bool

        :param infectious: print informations about infectious.
        :type infectious: bool

        :param migration: print informations about migration.
        :type migration: bool
        """
        self.simulation.print_populations(population=population, susceptibles=susceptibles, infectious=infectious, \
            migration=migration)

    def print_immunity_model(self, immunity=True, transition=True):
        """
        This methods prints the basic parameters of the immunity model.

        :param immunity: print informations about immunity.
        :type immunity: bool

        :param transition: print informations about immunity transition.
        :type transition: bool
        """
        self.simulation.print_immunity_model(immunity, transition)

    def print_mutations(self):
        self.simulation.print_mutations()

    def print_migrations(self):
        self.simulation.print_migrations()

    def print_all(self, basic_parameters=False, population=False, susceptible=False, infectious=False, migration=False, \
        immunity_model=False, immunity=False, transition=False):
        """
        This methods prints all the parameters of the simulation.

        :param basic_parameters: whether to output the basic parameters of the epidemiological model.
        :type basic_parameters: bool

        :param population: print informations about populations.
        :type population: bool

        :param susceptibles: print informations about susceptibles.
        :type susceptibles: bool

        :param infectious: print informations about infectious.
        :type infectious: bool

        :param migration: print informations about migration.
        :type migration: bool

        :param immunity: print informations about immunity.
        :type immunity: bool

        :param transition: print informations about immunity transition.
        :type transition: bool
        """
        if basic_parameters:
            self.simulation.print_basic_parameters()
        if population or susceptible or infectious or migration:
            self.simulation.print_populations(population=population, susceptibles=susceptibles, infectious=infectious, \
            migration=migration)
        if immunity or transition:
            self.simulation.print_immunity_model(immunity=immunity, transition=transition)


    def get_indexes_from_haplotype(self, haplotype):
        return np.array(self.simulation.create_list_for_cycles(haplotype, self.simulation.haplotype_number))

    @property
    def seed(self):
        return self.simulation.seed

    @property
    def sampling_probability(self):
        return self.simulation.sampling_probability

    @property
    def memory_optimization(self):
        return self.simulation.memory_optimization

    @property
    def number_of_sites(self):
        return self.simulation.number_of_sites

    @property
    def number_of_states_allele(self):
        return self.simulation.number_of_states_allele

    @property
    def haplotypes_number(self):
        return self.simulation.haplotypes_number

    @property
    def populations_number(self):
        return self.simulation.populations_number

    @property
    def number_of_susceptible_groups(self):
        return self.simulation.number_of_susceptible_groups

    @property
    def initial_haplotype(self):
        return self.simulation.initial_haplotype

    def set_initial_haplotype(self, amount):
        self.simulation.set_initial_haplotype(amount)

    @property
    def step_haplotype(self):
        return self.simulation.step_haplotype

    def set_step_haplotype(self, amount):
        self.simulation.set_step_haplotype(amount)

    @property
    def genome_length(self):
        return self.simulation.genome_length

    def set_genome_length(self, genome_length):
        self.simulation.set_genome_length(genome_length)

    @property
    def coinfection_parameters(self):
        return self.simulation.coinfection_parameters

    def set_coinfection_parameters(self, recombination):
        self.simulation.set_coinfection_parameters(recombination)

    @property
    def transmission_rate(self):
        return self.simulation.transmission_rate

    def set_transmission_rate(self, rate, haplotype=None):
        """
        .. _Wikipedia: https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

        Transmission rate is the the expected number of new infections from a single infected individual per time unit in the beginning of the epidemics when all but one of the hosts are susceptible. See `Wikipedia`_ - that is parameter beta in SIR model.

        :param rate: transmission rate value.
        :type rate: float

        :param haplotype: haplotypes for which the new value is being set. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
        :type haplotype: int(0 or 4) or string('T*' or 'AC') or int and string list([0, 4, 'T*']) or None
        """
        self.simulation.set_transmission_rate(rate, haplotype)

    @property
    def recovery_rate(self):
        return self.simulation.recovery_rate

    def set_recovery_rate(self, rate, haplotype=None):
        """
        .. _Wikipedia: https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

        Recovery rate is the inverse of the expected recovery time - non-negative float value. See `Wikipedia`_ - that is parameter gamma in SIR model.

        :param rate: recovery rate value.
        :type rate: float

        :param haplotype:haplotypes for which the new value is being set. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
        :type haplotype: int(0 or 4) or string('T*' or 'AC') or int and string list([0, 4, 'T*']) or None
        """
        self.simulation.set_recovery_rate(rate, haplotype)

    @property
    def sampling_rate(self):
        return self.simulation.sampling_rate

    def set_sampling_rate(self, rate, haplotype=None):
        """
        Sampling rate: the rate at which infected individuals are sampled. The genealogy will be generated only for these samples. Alternatively, if *sampling_probability* is **True**, it will specify the fraction of recovered individuals which are sampled.

        :param rate: sampling rate value.
        :type rate: float

        :param haplotype: haplotypes for which the new value is being set. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
        :type haplotype: int(0 or 4) or string('T*' or 'AC') or int and string list([0, 4, 'T*']) or None
        """
        self.simulation.set_sampling_rate(rate, haplotype)

    @property
    def mutation_rate(self):
        return self.simulation.mutation_rate
    
    def set_mutation_rate(self, rate, haplotype=None, mutation=None):
        """
        This method allows setting the mutation rate.

        :param rate: mutation rate value (None would not change the old value).
        :type rate: float or None

        :param haplotype: haplotypes for which the new value is being set. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
        :type haplotype: int(0 or 4) or string('T*' or 'AC') or int and string list([0, 4, 'T*']) or None

        :param mutation: the id of position on the haplotype where mutation arises.
        :type mutation: int or None
        """
        self.simulation.set_mutation_rate(rate, haplotype, mutation)

    @property
    def mutation_probabilities(self):
        return self.simulation.mutation_probabilities

    def set_mutation_probabilities(self, probabilities, haplotype=None, mutation=None):
        """
        This method allows setting the weights for each single nucleotide substitution (given the mutation happened).

        :param probabilities: weights of each single nucleotide substitution given mutation occured (None would not change the old values). `See for example <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
        :type probabilities: list of four non-negative integers or None

        :param haplotype: haplotypes for which the new value is being set. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
        :type haplotype: int(0 or 4) or string('T*' or 'AC') or int and string list([0, 4, 'T*']) or None

        :param mutation: the id of position on the haplotype where mutation arises.
        :type mutation: int or None
        """
        self.simulation.set_mutation_probabilities(probabilities, haplotype, mutation)

    @property
    def mutation_position(self):
        return self.simulation.mutation_position

    def set_mutation_position(self, mutation, position):
        self.simulation.set_mutation_position(mutation, position)


    @property
    def susceptibility_type(self):
        return self.simulation.susceptibility_type

    def set_susceptibility_type(self, susceptibility_type, haplotype=None):
        """
        The type of immunity (or the susceptibility group) which an individual gets after being infected with a pathogen of a haplotype.

        :param susceptibility_type: immunity group id.
        :type susceptibility_type: int

        :param haplotype: haplotypes for which the new value is being set. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
        :type haplotype: int(0 or 4) or string('T*' or 'AC') or int and string list([0, 4, 'T*']) or None
        """
        self.simulation.set_susceptibility_type(susceptibility_type, haplotype)

    @property
    def susceptibility(self):
        return self.simulation.susceptibility

    def set_susceptibility(self, rate, haplotype=None, susceptibility_type=None):
        """
        Susceptibility is a multiplicative modifier of the transmission rate based on the susceptible host immunity with a certain *susceptibility_type*. It can decrease or increase the transmission rate of a particular haplotype to the individuals of said *susceptibility_type*.

        :param rate: susceptibility value.
        :type rate: float

        :param haplotype: haplotypes for which the new value is being set. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
        :type haplotype: int(0 or 4) or string('T*' or 'AC') or int and string list([0, 4, 'T*']) or None

        :param susceptibility_type: immunity group id for which the new susceptibility value is being set.
        :type susceptibility_type: int or None
        """
        self.simulation.set_susceptibility(rate, haplotype, susceptibility_type)

    @property
    def immunity_transition(self):
        return self.simulation.immunity_transition

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
        self.simulation.set_immunity_transition(rate, source, target)


    @property
    def population_size(self):
        return self.simulation.population_size

    def set_population_size(self, size, population=None):
        """
        Set the number of individuals in the population.

        :param size: total number of individuals in the population.
        :type size: int

        :param population: population for which the new population size is being set (in case of None the value will be updated for all populations).
        :type population: int or None
        """
        self.simulation.set_population_size(size, population)

    @property
    def contact_density(self):
        return self.simulation.contact_density

    def set_contact_density(self, value, population=None):
        """
        The relative number of contacts per time units.

        :param value: contact density value.
        :type value: float

        :param population: population for which the new contact density size is being set (in case of None the value will be updated for all populations).
        :type population: int or None
        """
        self.simulation.set_contact_density(value, population)

    @property
    def npi(self):
        return self.simulation.npi

    def set_npi(self, parameters, population=None):
        """
        Setting conditions when lockdown in a population is imposed and lifted with the certain contact density during the lockdown.

        :param parameters: list with three elements: contact density value during lockdown, fraction of infectious population when the lockdown is set, fraction of infectious population when lockdown is lifted.
        :type parameters: list of length 3 with float values

        :param population: population for which the lockdown parameters are being set (in case of None the value will be updated for all populations).
        :type population: int or None
        """
        self.simulation.set_npi(parameters, population)

    @property
    def sampling_multiplier(self):
        return self.simulation.sampling_multiplier

    def set_sampling_multiplier(self, multiplier, population=None):
        """
        The relative sampling in the population (multiplicative modifier). Sampling rate of each haplotype is modified by this factor.

        :param value: sampling multiplier value.
        :type value: float

        :param population: population for which the new sampling multiplier is being set (in case of None the value will be updated for all populations).
        :type population: int or None
        """
        self.simulation.set_sampling_multiplier(multiplier, population)

    @property
    def migration_probability(self):
        return self.simulation.migration_probability

    def set_migration_probability(self, probability, source=None, target=None):
        """
        The probability that an individual from the population source is travelling to the population target.

        :param probability: probability of migration value.
        :type probability: float

        :param source: source population with infectious individual (None means that the new value will be set to all populations as source).
        :type source: int or None

        :param target: target population with infectious individual (None means that the new value will be set to all populations as source).
        :type target: int or None
        """
        self.simulation.set_migration_probability(probability, source, target)

    def set_total_migration_probability(self, total_probability):
        """
        The probability that an individual from the population source is travelling to the population target.

        :param total_probability: if True, all the entries (except for the diagonal element) of row corresponding to source population is filled with probabilitie/(K-1). That represents individual travelling out of its population to a random destination.
        :type total_probability: bool
        """
        self.simulation.set_total_migration_probability(total_probability)

    @property
    def susceptible(self):
        return self.simulation.susceptible

    def set_susceptible(self, amount, source_type, target_type, population=None):
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
        self.simulation.set_susceptible(amount, source_type, target_type, population)

    @property
    def infectious(self):
        return self.simulation.infectious

    def set_infectious(self, amount, source_type, target_haplotype, population=None):
        """
        The method allows moving individuals from susceptibility group to target haplotype.

        :param amount: the number of individuals to be moved.
        :type amount: int

        :param source_type: source susceptibility group from which individuals are moved.
        :type source_type: int

        :param target_haplotype: target haplotype into which individuals are moved.
        :type target_haplotype: int(0 or 4) or string('T*' or 'AC')

        :param population: population for which the new value is being set (in case of None the value will be updated for all populations).
        :type population: int or None
        """
        self.simulation.set_infectious(amount, source_type, target_haplotype, population)

    def set_chain_events(self, file_name):
        """
        Allows to import chain of events directly from the .npy file

        :param file_name: a name of the file (without the .npy extension)
        :type file_name: str
        """
        self.simulation.set_chain_events(file_name)

    def set_settings(self, file_template):
        """
        Currently is a placeholder and does nothing. Do not include in actual docs until it does something.
        """
        self.simulation.set_settings(file_template)

    def set_state(self, file_template):
        """
        Calls set_chain_events and set_setting. Now identical to set_events. Do not include in docs until it starts to differ.
        """
        self.set_chain_events(file_template)
        self.set_settings(file_template)

    # def export_newick всё что в файлы это export
    def export_newick(self, file_template=None, file_path = None):
        """
        Exports the result of simulation in the newick format.

        :param file_template: template for the file name
        :type file_template: str

        :param file_path: path to output file
        :type file_path: str
        """
        pruferSeq, times, mut, populations = self.simulation.output_tree_mutations()
        writeGenomeNewick(pruferSeq, times, populations, file_template, file_path)

    def export_mutations(self, file_template=None, file_path = None):
        """
        Exports the information about all mutations

        :param file_template: template for the file name
        :type file_template: str

        :param file_path: path to output file
        :type file_path: str
        """
        pruferSeq, times, mut, populations = self.simulation.output_tree_mutations()
        writeMutations(mut, len(pruferSeq), file_template, file_path)

    def export_migrations(self, file_template=None, file_path = None):
        """
        Exports the information about all migrations

        :param file_template: template for the file name
        :type file_template: str

        :param file_path: path to output file
        :type file_path: str
        """
        self.simulation.export_migrations(file_template, file_path)

    def output_sample_data(self, output_print=False):
        """
        Show all information about sampling, time, place, etc.

        :param output_print: if *Thue*, it outputs data about each sample to the console. Data includes time of sampling, population and haplotype of sampled individual.
        :type output_print: bool

        :return: 3 lists with times, population ids and haplotypes of the sampled individuals
        """
        time, pop, hap = self.simulation.output_sample_data()
        if output_print:
            return time, pop, hap
        else:
            print(time)
            print(pop)
            print(hap)

    def output_epidemiology_timelines(self, step=1000, output_file=False):
        """
        Records simulation state changes over some period of time.

        :param step: a number of steps which split the epidemiology timelines.
        :type step: str

        :param output_file: if **True**, the timeline is written in the output file
        :type output_file: bool
        """
        if output_file:
            self.simulation.output_epidemiology_timelines(step, output_file)
        else:
            return self.simulation.output_epidemiology_timelines(step, output_file)

    def export_chain_events(self, file_name="chain_events"):
        """
        Exports event chain.

        :param file_name: file name for the output
        :param file_name: str
        """
        self.simulation.export_chain_events(file_name)

    def export_settings(self, file_template="parameters"):
        """
        Exports settings.

        :param file_template: file name for the output
        :param file_template: str
        """
        self.simulation.export_settings(file_template)

    def export_state(self):
        """
        Exports event chain and settings.
        """
        self.output_chain_events()
        self.output_settings()

    def export_ts(self):
        return self.simulation.export_ts()

    def get_tree(self):
        return self.simulation.get_tree()

    def get_data_susceptible(self, population, susceptibility_type,
                             step_num):  # returns susceptible, time_points, lockdowns
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
        return self.simulation.get_data_susceptible(population, susceptibility_type, step_num)

    def get_data_infectious(self, population, haplotype,
                            step_num):  # returns infections, sample, time_points, lockdowns
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
        return self.simulation.get_data_infectious(population, haplotype, step_num)

    def add_plot_infectious(self, population, haplotype, step_num=100, label_infectious=None, label_samples=None):
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
        if self.fig == None:
            self.fig, self.ax = plt.subplots(figsize=(8, 6))
            self.ax.set_ylabel('Number of samples')
            self.ax.set_xlabel('Time')
            self.ax_2 = self.ax.twinx()
            self.ax_2.set_ylabel('Number of individuals')

        if isinstance(haplotype, int) == True:
            self.plot_infectious(population, haplotype, step_num, label_infectious, label_samples)
        elif isinstance(haplotype, str) == True:
            haplotypes = self.simulation.create_list_for_cycles(haplotype, self.simulation.haplotype_number)
            for hi in haplotypes:
                self.plot_infectious(population, hi, step_num, label_infectious, label_samples)
        else:
            print("#TODO")

    def plot_infectious(self, population, haplotype, step_num, label_infectious, label_samples):
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
        infections, sample, time_points, lockdowns = self.simulation.get_data_infectious(population, haplotype, step_num)
        if label_infectious == None:
            self.ax_2.plot(time_points, infections, label='Infectious pop:' + str(population) + ' hap:' + self.simulation.calculate_string(haplotype))
        elif isinstance(label_infectious, str) == True:
            self.ax_2.plot(time_points, infections, label=label_infectious)
        else:
            print("#TODO")

        if label_samples == None:
            self.ax.plot(time_points, sample, "--", label='Samples pop:' + str(population) + ' hap:' + self.simulation.calculate_string(haplotype))
        elif isinstance(label_samples, str) == True:
            self.ax.plot(time_points, sample, "--", label=label_samples)
        else:
            print("#TODO")

        if len(lockdowns) != 0:
            point = 0
            pointEnd = 0
            for ld in range(0, len(lockdowns), 2):
                if ld+1 == len(lockdowns):
                    while time_points[point] < lockdowns[ld][1]:
                        point += 1
                    plt.fill_between(time_points[point:], infections[point:], alpha=0.2)
                else:
                    while time_points[point] < lockdowns[ld][1]:
                        point += 1
                    while time_points[pointEnd] < lockdowns[ld+1][1]:
                        pointEnd += 1
                    if pointEnd == point:
                        continue
                    plt.fill_between(time_points[point:pointEnd+1], infections[point:pointEnd+1], alpha=0.2)

    def add_plot_susceptible(self, population, susceptibility_type, step_num=100, label_susceptible=None):
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
        if self.fig == None:
            self.fig, self.ax = plt.subplots(figsize=(8, 6))
            self.ax.set_ylabel('Number of samples')
            self.ax.set_xlabel('Time')
            self.ax_2 = self.ax.twinx()
            self.ax_2.set_ylabel('Number of individuals')

        susceptible, time_points, lockdowns = self.simulation.get_data_susceptible(population, susceptibility_type, step_num)
        if label_susceptible == None:
            self.ax_2.plot(time_points, susceptible, label='Susceptible pop:' + str(population) + ' sus:' + str(susceptibility_type))
        elif isinstance(label_susceptible, str) == True:
            self.ax_2.plot(time_points, susceptible, label=label_susceptible)
        else:
            print("#TODO")

        if len(lockdowns) != 0:
            point = 0
            pointEnd = 0
            for ld in range(0, len(lockdowns), 2):
                if ld+1 == len(lockdowns):
                    while time_points[point] < lockdowns[ld][1]:
                        point += 1
                    plt.fill_between(time_points[point:], susceptible[point:], alpha=0.2)
                else:
                    while time_points[point] < lockdowns[ld][1]:
                        point += 1
                    while time_points[pointEnd] < lockdowns[ld+1][1]:
                        pointEnd += 1
                    if pointEnd == point:
                        continue
                    plt.fill_between(time_points[point:pointEnd+1], susceptible[point:pointEnd+1], alpha=0.2)

    def add_legend(self):
        """
        Add legend to the plot.
        """
        lines_1, labels_1 = self.ax.get_legend_handles_labels()
        lines_2, labels_2 = self.ax_2.get_legend_handles_labels()
        lines = lines_1 + lines_2
        labels = labels_1 + labels_2
        self.ax.legend(lines, labels, loc=0)

    def add_title(self, name="Plot"):
        """
        Add plot title.

        :param name: the title for the plot
        :type name: str
        """
        self.ax.set_title(name)

    def plot(self, name_file=None):
        """
        Show plot with all the trajectories added by other methods.
        """
        if name_file:
            plt.savefig(name_file)
        else:
            plt.show()
        self.fig = None


    def simulate(self, iterations=1000, sample_size=None, epidemic_time=-1, method='direct', attempts=200):
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
        if sample_size is None:
            sample_size = iterations
        if epidemic_time is None:
            epidemic_time = -1

        start_time = time.time()
        if method == 'direct':
            self.simulation.SimulatePopulation(iterations, sample_size, epidemic_time, attempts)
            self.simulation.Stats(time.time() - start_time)
        elif method == 'tau':
            self.simulation.SimulatePopulation_tau(iterations, sample_size, epidemic_time, attempts)
            self.simulation.Stats(time.time() - start_time)
        else:
            print("Unknown method. Choose between 'direct' and 'tau'.")

    def genealogy(self, seed=None):
        """
        Generates a genealogy based on the chain of events generated during all the instances of simulate() method.

        :param seed: seed value (None for random seed).
        :type seed: int or None
        """
        start_time = time.time()
        self.simulation.GetGenealogy(seed)
        print(f"Getting genealogy time: {time.time() - start_time}")

    def print_recomb(self, left, right):
        self.simulation.print_recomb(left, right)

    def print_chain(self):
        self.simulation.print_chain()

    def print_tree(self):
        self.simulation.print_tree()


    def citation(self):
        """
        Prints a citation of the paper in the console
        """
        print("VGsim: scalable viral genealogy simulator for global pandemic")
        print("Vladimir Shchur, Vadim Spirin, Dmitry Sirotkin, EvgeniBurovski, Nicola De Maio, Russell Corbett-Detig")
        print("medRxiv 2021.04.21.21255891; doi: https://doi.org/10.1101/2021.04.21.21255891")

    def debug(self):
        """
        Returns all of the available information about simulation. Mainly used for the debug purposes. Not documented.
        """

        self.simulation.Debug()

    def get_proportion(self):
        """
        Used for paper preparation. Not documented.
        """
        return self.simulation.get_proportion()

    def print_counters(self):
        """
        Used for paper preparation. Not documented.
        """
        self.simulation.PrintCounters()

    def print_propensities(self):
        """
        Used for paper preparation. Not documented.
        """
        self.simulation.PrintPropensities()
