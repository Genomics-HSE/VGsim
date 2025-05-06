from random import randrange
import sys
import numpy as np
from tabulate import tabulate
import matplotlib.pyplot as plt

import source_VGsim

class Simulator:
    def __init__(self,
                 number_of_sites=0,
                 number_of_allele_states=4,
                 number_of_populations=1,
                 number_of_susceptible_groups=1,
                 seed=None,
                 sampling_probability=False):
        self._check_amount(number_of_sites, 'number of sites', zero=False)
        self._check_amount(number_of_allele_states, 'number of allele states')
        self._check_amount(number_of_populations, 'populations number')
        self._check_amount(number_of_susceptible_groups, 'number of susceptible groups')
        if not isinstance(sampling_probability, bool):
            raise ValueError('Incorrect value of sampling probability. Value of sampling probability should be True or False.')
        if seed is None:
            seed = int(randrange(sys.maxsize))
        self._check_amount(seed, 'seed', zero=False)

        self._simulator = source_VGsim.Simulator(number_of_sites, number_of_allele_states,
                                                 number_of_populations, number_of_susceptible_groups,
                                                 seed)
        self._number_of_sites = number_of_sites
        self._number_of_allele_states = number_of_allele_states
        self._number_of_haplotypes = number_of_allele_states**number_of_sites
        self._number_of_populations = number_of_populations
        self._number_of_susceptible_groups = number_of_susceptible_groups
        self._seed = seed
        self._sampling_probability = sampling_probability
        self._first_simulation = False

        self.fig = None

    def simulate(self, iterations=1000, sample_size=None, epidemic_time=0.0, method='direct', attempts=200):
        self._first_simulation = True
        if sample_size is None:
            sample_size = iterations
        self._simulator.simulate(iterations, sample_size, epidemic_time, method, attempts)
        self._simulator.Debug()

    def genealogy(self, seed):
        self._simulator.genealogy()

    def get_flat_chain(self):
        return self._simulator.get_flat_chain()

    def get_tree(self):
        return self._simulator.get_tree()

    def export_chain_events(self, file_name="chain_events"):
        chain = self._simulator.export_chain_events()
        np.save(file_name, chain)


    def get_number_of_sites(self):
        return self._number_of_sites

    def get_number_of_allele_states(self):
        return self._number_of_allele_states

    def get_number_of_haplotypes(self):
        return self._number_of_haplotypes

    def get_number_of_populations(self):
        return self._number_of_populations

    def get_number_of_susceptible_groups(self):
        return self._number_of_susceptible_groups

    def get_seed(self):
        return self._seed

    def get_sampling_probability(self):
        return self._sampling_probability

    def set_susceptibility_group(self, susceptibility_group, haplotype=None):
        """
        The type of immunity (or the susceptibility group) which an individual gets after being infected with a pathogen of a haplotype.

        :param susceptibility_group: immunity group id.
        :type susceptibility_group: int

        :param haplotype: haplotypes for which the new value is being set. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
        :type haplotype: int(0 or 4) or string('T*' or 'AC') or int and string list([0, 4, 'T*']) or None
        """
        if not isinstance(susceptibility_group, int):
            raise TypeError('Incorrect type of susceptibility group. Type should be int.')
        elif susceptibility_group < 0 or susceptibility_group >= self._number_of_susceptible_groups:
            raise IndexError('There are no such susceptibility group!')
        self._check_indexes(haplotype, self._number_of_haplotypes, 'haplotype', True)
        haplotypes = self._calculate_indexes(haplotype, self._number_of_haplotypes)

        for hn in haplotypes:
            self._simulator.set_susceptibility_group(susceptibility_group, hn)

    def get_susceptibility_group(self):
        return self._simulator.get_susceptibility_group()

    def set_transmission_rate(self, rate, haplotype=None):
        """
        .. _Wikipedia: https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

        Transmission rate is the the expected number of new infections from a single infected individual per time unit in the beginning of the epidemics when all but one of the hosts are susceptible. See `Wikipedia`_ - that is parameter beta in SIR model.

        :param rate: transmission rate value.
        :type rate: float

        :param haplotype: haplotypes for which the new value is being set. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
        :type haplotype: int(0 or 4) or string('T*' or 'AC') or int and string list([0, 4, 'T*']) or None
        """
        self._check_value(rate, 'transmission rate')
        self._check_indexes(haplotype, self._number_of_haplotypes, 'haplotype', True)
        haplotypes = self._calculate_indexes(haplotype, self._number_of_haplotypes)
        for hn in haplotypes:
            self._simulator.set_transmission_rate(rate, hn)

    def get_transmission_rate(self):
        return self._simulator.get_transmission_rate()

    def set_recovery_rate(self, rate, haplotype=None):
        """
        .. _Wikipedia: https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

        Transmission rate is the the expected number of new infections from a single infected individual per time unit in the beginning of the epidemics when all but one of the hosts are susceptible. See `Wikipedia`_ - that is parameter beta in SIR model.

        :param rate: transmission rate value.
        :type rate: float

        :param haplotype: haplotypes for which the new value is being set. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
        :type haplotype: int(0 or 4) or string('T*' or 'AC') or int and string list([0, 4, 'T*']) or None
        """
        self._check_value(rate, 'recovery rate')
        self._check_indexes(haplotype, self._number_of_haplotypes, 'haplotype', True)
        haplotypes = self._calculate_indexes(haplotype, self._number_of_haplotypes)
        for hn in haplotypes:
            self._simulator.set_recovery_rate(rate, hn)

    def get_recovery_rate(self):
        return self._simulator.get_recovery_rate()

    def set_sampling_rate(self, rate, haplotype=None):
        """
        .. _Wikipedia: https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

        Transmission rate is the the expected number of new infections from a single infected individual per time unit in the beginning of the epidemics when all but one of the hosts are susceptible. See `Wikipedia`_ - that is parameter beta in SIR model.

        :param rate: transmission rate value.
        :type rate: float

        :param haplotype: haplotypes for which the new value is being set. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
        :type haplotype: int(0 or 4) or string('T*' or 'AC') or int and string list([0, 4, 'T*']) or None
        """
        self._check_indexes(haplotype, self._number_of_haplotypes, 'haplotype', True)
        haplotypes = self._calculate_indexes(haplotype, self._number_of_haplotypes)

        if self._sampling_probability == True:
            self._check_value(rate, 'sampling probability', edge=1)
            recovery_rate = self._simulator.get_recovery_rate()
            sampling_rate = self._simulator.get_sampling_rate()

            for hn in haplotypes:
                deathRate = recovery_rate[hn] + sampling_rate[hn]
                self._simulator.set_recovery_rate((1 - rate) * deathRate, hn)
                self._simulator.set_sampling_rate(rate * deathRate, hn)
        else:
            self._check_value(rate, 'sampling rate')

            for hn in haplotypes:
                self._simulator.set_sampling_rate(rate, hn)

    def get_sampling_rate(self):
        return self._simulator.get_sampling_rate()

    def set_mutation_rate(self, rate, haplotype=None, mutation=None):
        """
        .. _Wikipedia: https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

        Transmission rate is the the expected number of new infections from a single infected individual per time unit in the beginning of the epidemics when all but one of the hosts are susceptible. See `Wikipedia`_ - that is parameter beta in SIR model.

        :param rate: transmission rate value.
        :type rate: float

        :param haplotype: haplotypes for which the new value is being set. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
        :type haplotype: int(0 or 4) or string('T*' or 'AC') or int and string list([0, 4, 'T*']) or None
        """
        self._check_value(rate, 'mutation rate')
        self._check_indexes(haplotype, self._number_of_haplotypes, 'haplotype', True)
        self._check_indexes(mutation, self._number_of_sites, 'mutation site')
        haplotypes = self._calculate_indexes(haplotype, self._number_of_haplotypes)
        sites = self._calculate_indexes(mutation, self._number_of_sites)

        for hn in haplotypes:
            for s in sites:
                self._simulator.set_mutation_rate(rate, hn, s)

    def get_mutation_rate(self):
        return self._simulator.get_mutation_rate()

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
        self._check_list(probabilities, 'probabilities list', self._number_of_allele_states)
        if isinstance(probabilities, list):
            for probability in probabilities:
                self._check_value(probability, 'mutation probabilities')
        self._check_indexes(haplotype, self._number_of_haplotypes, 'haplotype', True)
        self._check_indexes(mutation, self._number_of_sites, 'mutation site')
        haplotypes = self._calculate_indexes(haplotype, self._number_of_haplotypes)
        sites = self._calculate_indexes(mutation, self._number_of_sites)

        for hn in haplotypes:
            for s in sites:
                probabilities_base = list(probabilities)
                del probabilities_base[self._calculate_base(hn, s)]
                if sum(probabilities_base) == 0:
                    raise ValueError('Incorrect probabilities list. The sum of three elements without mutation allele should be more 0.')
                for i in range(3):
                    self._simulator.set_mutation_probabilities(probabilities_base[i], hn, s, i)


    def get_mutation_probabilities(self):
        return self._simulator.get_mutation_probabilities()

    def set_susceptibility(self, rate, haplotype=None, susceptibility_group=None):
        """
        Susceptibility is a multiplicative modifier of the transmission rate based on the susceptible host immunity with a certain *susceptibility_group*. It can decrease or increase the transmission rate of a particular haplotype to the individuals of said *susceptibility_group*.

        :param rate: susceptibility value.
        :type rate: float

        :param haplotype: haplotypes for which the new value is being set. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
        :type haplotype: int(0 or 4) or string('T*' or 'AC') or int and string list([0, 4, 'T*']) or None

        :param susceptibility_group: immunity group id for which the new susceptibility value is being set.
        :type susceptibility_group: int or None
        """
        self._check_value(rate, 'susceptibility rate')
        self._check_indexes(haplotype, self._number_of_haplotypes, 'haplotype', True)
        self._check_indexes(susceptibility_group, self._number_of_susceptible_groups, 'susceptibility group')
        haplotypes = self._calculate_indexes(haplotype, self._number_of_haplotypes)
        sus_types = self._calculate_indexes(susceptibility_group, self._number_of_susceptible_groups)

        for hn in haplotypes:
            for sn in sus_types:
                self._simulator.set_susceptibility(rate, hn, sn)

    def get_susceptibility(self):
        return self._simulator.get_susceptibility()

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
        self._check_value(rate, 'immunity transition rate')
        self._check_indexes(source, self._number_of_susceptible_groups, 'susceptibility group')
        self._check_indexes(target, self._number_of_susceptible_groups, 'susceptibility group')
        sus_types_1 = self._calculate_indexes(source, self._number_of_susceptible_groups)
        sus_types_2 = self._calculate_indexes(target, self._number_of_susceptible_groups)

        for sn1 in sus_types_1:
            for sn2 in sus_types_2:
                if sn1 == sn2:
                    continue
                self._simulator.set_immunity_transition(rate, sn1, sn2)

    def get_immunity_transition(self):
        return self._simulator.get_immunity_transition()

    def set_population_size(self, size, population=None):
        """
        Set the number of individuals in the population.

        :param size: total number of individuals in the population.
        :type size: int

        :param population: population for which the new population size is being set (in case of None the value will be updated for all populations).
        :type population: int or None
        """
        if self._first_simulation == True:
            raise ValueError('Changing population size is available only before first simulation!')
        self._check_amount(size, 'population size')
        self._check_index(population, self._number_of_populations, 'population')
        populations = self._calculate_index(population, self._number_of_populations)

        for pn in populations:
            self._simulator.set_population_size(size, pn)

    def get_population_size(self):
        return self._simulator.get_population_size()

    def set_contact_density(self, value, population=None):
        """
        The relative number of contacts per time units.

        :param value: contact density value.
        :type value: float

        :param population: population for which the new contact density size is being set (in case of None the value will be updated for all populations).
        :type population: int or None
        """
        self._check_value(value, 'contact density')
        self._check_indexes(population, self._number_of_populations, 'population')
        populations = self._calculate_indexes(population, self._number_of_populations)

        for pn in populations:
            self._simulator.set_contact_density(value, pn)

    def get_contact_density(self):
        return self._simulator.get_contact_density()

    def set_npi(self, parameters, population=None):
        """
        Setting conditions when lockdown in a population is imposed and lifted with the certain contact density during the lockdown.

        :param parameters: list with three elements: contact density value during lockdown, fraction of infectious population when the lockdown is set, fraction of infectious population when lockdown is lifted.
        :type parameters: list of length 3 with float values

        :param population: population for which the lockdown parameters are being set (in case of None the value will be updated for all populations).
        :type population: int or None
        """
        self._check_list(parameters, 'npi parameters', 3)
        self._check_value(parameters[0], 'first npi parameter')
        self._check_value(parameters[1], 'second npi parameter', edge=1)
        self._check_value(parameters[2], 'third npi parameter', edge=1)
        self._check_indexes(population, self._number_of_populations, 'population')
        populations = self._calculate_indexes(population, self._number_of_populations)

        for pn in populations:
            self._simulator.set_npi(parameters[0], parameters[1], parameters[2], pn)

    def get_npi(self):
        return self._simulator.get_npi()

    def set_sampling_multiplier(self, multiplier, population=None):
        """
        The relative sampling in the population (multiplicative modifier). Sampling rate of each haplotype is modified by this factor.

        :param value: sampling multiplier value.
        :type value: float

        :param population: population for which the new sampling multiplier is being set (in case of None the value will be updated for all populations).
        :type population: int or None
        """
        self._check_value(multiplier, 'sampling multiplier')
        self._check_indexes(population, self._number_of_populations, 'population')
        populations = self._calculate_indexes(population, self._number_of_populations)

        for pn in populations:
            self._simulator.set_sampling_multiplier(multiplier, pn)

    def get_sampling_multiplier(self):
        return self._simulator.get_sampling_multiplier()

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
        self._check_value(probability, 'migration probability', edge=1)
        self._check_indexes(source, self._number_of_populations, 'population')
        self._check_indexes(target, self._number_of_populations, 'population')
        populations_1 = self._calculate_indexes(source, self._number_of_populations)
        populations_2 = self._calculate_indexes(target, self._number_of_populations)

        for pn1 in populations_1:
            for pn2 in populations_2:
                if pn1 == pn2:
                    continue
                self._simulator.set_migration_probability(probability, pn1, pn2)

        self._check_migration_probability()

    def set_total_migration_probability(self, total_probability):
        """
        The probability that an individual from the population source is travelling to the population target.

        :param total_probability: if True, all the entries (except for the diagonal element) of row corresponding to source population is filled with probabilitie/(K-1). That represents individual travelling out of its population to a random destination.
        :type total_probability: bool
        """
        self._check_value(total_probability, 'total migration probability', edge=1)
        source_rate = 1.0 - total_probability
        target_rate = total_probability / (self._number_of_populations - 1)

        for pn1 in range(self._number_of_populations):
            for pn2 in range(self._number_of_populations):
                if pn1 == pn2:
                    self._simulator.set_migration_probability(source_rate, pn1, pn2)
                else:
                    self._simulator.set_migration_probability(target_rate, pn1, pn2)

        self._check_migration_probability()

    def get_migration_probability(self):
        return self._simulator.get_migration_probability()


    def print_basic_parameters(self):
        """
        This methods prints the basic parameters of the epidemiological model.
        """
        field = ["H", "TR", "RR", "SR", "ST"]
        for s in range(self._number_of_sites):
            field.append(f"M{s}")
            field.append(f"MW{s}")

        transmission_rate = self.get_transmission_rate()
        recovery_rate = self.get_recovery_rate()
        sampling_rate = self.get_sampling_rate()
        susceptibility_group = self.get_susceptibility_group()
        mutation_rate = self.get_mutation_rate()
        mutation_probabilities = self.get_mutation_probabilities()

        data = []
        for hn in range(self._number_of_haplotypes):
            row = [
                    self._calculate_string_from_number_haplotype(hn),
                    transmission_rate[hn],
                    recovery_rate[hn],
                    sampling_rate[hn],
                    susceptibility_group[hn]
                   ]
            for s in range(self._number_of_sites):
                row.append(mutation_rate[hn][s])
                row.append(self._calculate_colored_haplotype(mutation_probabilities, hn, s))
            data.append(row)

        print(tabulate(data, headers=field, tablefmt="fancy_grid", numalign="center", stralign="center"))
        print("Legend:")
        print("H - haplotype")
        print("TR - transmission rate")
        print("RR - recovery rate")
        print("SR - sampling rate")
        print("ST - susceptibility type")
        for s in range(self._number_of_sites):
            print(f"M{s} - {s} mutation rate")
            print(f"MW{s} - {s} mutation weights")
        print()

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
        if susceptibles or infectious:
            current_susceptible, current_infectious = self._simulator.get_current_individuals()

        if population:            
            field = ["ID", "Size", 'Actual size', "CD",'CDBLC', "CDALD", "SLD", "ELD", "SM"]

            sizes = self.get_population_size()
            actual_sizes = self._simulator.get_actual_size()
            contact_density = self.get_contact_density()
            contact_density_before_lockdown = self._simulator.get_contact_density_before_lockdown()
            contact_density_after_lockdown, start_lockdown, end_lockdown = self.get_npi()
            sampling_multiplier = self.get_sampling_multiplier()

            data = []
            for pn in range(self._number_of_populations):
                data.append([pn,
                             sizes[pn],
                             actual_sizes[pn],
                             contact_density[pn],
                             contact_density_before_lockdown[pn],
                             contact_density_after_lockdown[pn],
                             start_lockdown[pn],
                             end_lockdown[pn],
                             sampling_multiplier[pn]])

            print(tabulate(data, headers=field, tablefmt="fancy_grid", numalign="center", stralign="center"))
            print("Legend:")
            print("ID - number of population")
            print("Size - size of population")
            print("Actual size - actual size of population")
            print("CD - contact density")
            print("CDBLD - contact density without lockdown")
            print("CDALD - contact density at lockdown")
            print("SLD - start of lockdown")
            print("ELD - end of lockdown")
            print("SM - sampling multiplier")
            print()

        if susceptibles:
            field = ["ST\\ID"]
            for pn in range(self._number_of_populations):
                field.append(pn)

            data = []
            for sn in range(self._number_of_susceptible_groups):
                row = [sn]
                for pn in range(self._number_of_populations):
                    row.append(current_susceptible[pn][sn])
                data.append(row)

            print(tabulate(data, headers=field, tablefmt="fancy_grid", numalign="center", stralign="center"))
            print("Legend:")
            print("ID - ID population")
            print("ST - susceptibility type")
            print()

        if infectious:
            field = ["H\\ID"]
            for pn in range(self._number_of_populations):
                field.append(pn)

            data = []
            for hn in range(self._number_of_haplotypes):
                row = [self._calculate_string_from_number_haplotype(hn)]
                for pn in range(self._number_of_populations):
                    row.append(current_infectious[pn][hn])
                data.append(row)

            print(tabulate(data, headers=field, tablefmt="fancy_grid", numalign="center", stralign="center"))
            print("Legend:")
            print("ID - ID population")
            print("H - haplotype")
            print()

        if migration:
            field = ["S\\T"]
            for pn in range(self._number_of_populations):
                field.append(pn)

            migration_probability = self.get_migration_probability()

            data = []
            for pn1 in range(self._number_of_populations):
                row = [pn1]
                for pn2 in range(self._number_of_populations):
                    row.append(migration_probability[pn1][pn2])
                data.append(row)

            print(tabulate(data, headers=field, tablefmt="fancy_grid", numalign="center", stralign="center"))
            print("Legend:")
            print("S - ID source population")
            print("T - ID target population")
            print()

    def print_immunity_model(self, immunity=True, transition=True):
        """
        This methods prints the basic parameters of the immunity model.

        :param immunity: print informations about immunity.
        :type immunity: bool

        :param transition: print informations about immunity transition.
        :type transition: bool
        """
        if immunity:
            field = ["H\\ST"]
            for sn in range(self._number_of_susceptible_groups):
                field.append(f"S{sn}")

            susceptibility = self.get_susceptibility()

            data = []
            for hn in range(self._number_of_haplotypes):
                row = [self._calculate_string_from_number_haplotype(hn)]
                for sn in range(self._number_of_susceptible_groups):
                    row.append(susceptibility[hn][sn])
                data.append(row)

            print(tabulate(data, headers=field, tablefmt="fancy_grid", numalign="center", stralign="center"))
            print("Legend:")
            print("H - haplotype")
            print("ST - susceptibility type")
            print()

        if transition:
            field = ["ID"]
            for sn in range(self._number_of_susceptible_groups):
                field.append(sn)

            immunity_transition = self.get_immunity_transition()

            data = []
            for sn1 in range(self._number_of_susceptible_groups):
                row = [sn1]
                for sn2 in range(self._number_of_susceptible_groups):
                    row.append(immunity_transition[sn1][sn2])
                data.append(row)

            print(tabulate(data, headers=field, tablefmt="fancy_grid", numalign="center", stralign="center"))
            print("Legend:")
            print("ID - ID susceptibility type")
            print()

    # def print_mutations(self):
    #     self.simulation.print_mutations()

    # def print_migrations(self):
    #     self.simulation.print_migrations()


    def print_all(self, basic_parameters=False,
                        population=False,
                        susceptible=False,
                        infectious=False,
                        migration=False,
                        immunity_model=False,
                        immunity=False,
                        transition=False):
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
            self.print_basic_parameters()
        if population or susceptible or infectious or migration:
            self.print_populations(population=population,
                                   susceptibles=susceptibles,
                                   infectious=infectious,
                                   migration=migration)
        if immunity or transition:
            self.print_immunity_model(immunity=immunity,
                                      transition=transition)





    def get_data_susceptible(self, population=None, susceptibility_group=None, step_number=100):
        """
        Returns the list with the amount of susceptible individuals for the populations and the haplotypes for the entire period of time.

        :param population: Population to retrieve data from. By default - None.
        :type population: int or list of int or None

        :param susceptibility_group: Susceptibility group to retrieve data from. By default - None.
        :type susceptibility_group: int or list of int or None

        :param step_number: Number of steps which split the epidemiology timeline. By default - 100.
        :type step_number: int

        :return: A list of a 3-tuple: population, susceptibility group and list of amount of susceptible individuals.
        """
        self._check_indexes(population, self._number_of_populations, 'population')
        self._check_indexes(susceptibility_group, self._number_of_susceptible_groups, 'susceptibility group')
        self._check_amount(step_number, 'step number', False)

        populations = self._calculate_indexes(population, self._number_of_populations)
        susceptibility_groups = self._calculate_indexes(susceptibility_group, self._number_of_susceptible_groups)

        data = []
        for pn in populations:
            for sn in susceptibility_groups:
                data.append((pn, sn, self._simulator.get_data_susceptible(pn, sn, step_number)))

        return data

    def get_data_infected(self, population=None, haplotype=None, step_number=100):
        """
        Returns the list with the amount of infected individuals for the populations and the haplotypes for the entire period of time.

        :param population: Population to retrieve data from. By default - None.
        :type population: int or list of int or None

        :param haplotype: Haplotype to retrieve data from. By default - None.
        :type haplotype: int or str or list of int and str or None

        :param step_number: Number of steps which split the epidemiology timeline. By default - 100.
        :type step_number: int

        :return: A list of a 3-tuple: population, haplotype and list of amount of infected individuals.
        """
        self._check_indexes(population, self._number_of_populations, 'population')
        self._check_indexes(haplotype, self._number_of_haplotypes, 'haplotype')
        self._check_amount(step_number, 'step number', False)

        populations = self._calculate_indexes(population, self._number_of_populations)
        haplotypes = self._calculate_indexes(haplotype, self._number_of_haplotypes)

        data = []
        for pn in populations:
            for hn in haplotypes:
                data.append((pn, hn, self._simulator.get_data_infected(pn, hn, step_number)))

        return data

    def get_data_sample(self, population=None, haplotype=None, step_number=100):
        """
        Returns the list with the amount of sample individuals for the populations and the haplotypes for the entire period of time.

        :param population: Population to retrieve data from. By default - None.
        :type population: int or list of int or None

        :param haplotype: Haplotype to retrieve data from. By default - None.
        :type haplotype: int or str or list of int and str or None

        :param step_number: Number of steps which split the epidemiology timeline. By default - 100.
        :type step_number: int

        :return: A list of a 3-tuple: population, haplotype and list of amount of sample individuals.
        """
        self._check_indexes(population, self._number_of_populations, 'population')
        self._check_indexes(haplotype, self._number_of_haplotypes, 'haplotype')
        self._check_amount(step_number, 'step number', False)

        populations = self._calculate_indexes(population, self._number_of_populations)
        haplotypes = self._calculate_indexes(haplotype, self._number_of_haplotypes)

        data = []
        for pn in populations:
            for hn in haplotypes:
                data.append((pn, hn, self._simulator.get_data_sample(pn, hn, step_number)))

        return data

    def get_time_points(self, step_number=100):
        """
        Returns the list with the time points for the entire period of time.

        :param step_number: Number of steps which split the epidemiology timeline. By default - 100.
        :type step_number: int

        :return: List with the time points.
        """
        self._check_amount(step_number, 'step number', False)

        return self._simulator.get_time_points(step_number)

    # def get_lockdowns(self, population):
    #     """
    #     Returns the list with the information about lockdowns over some period of time.

    #     :param population: Population to retrieve data from. By default - None.
    #     :type population: int or list of int or None

    #     :return: 3 lists - amounts of sample individuals on the time points, time points themselves at which we retrieve amounts of individuals, and lockdowns data for the current epidemiology
    #     """
    #     return self._simulator.get_lockdowns(population)

    def add_plot_susceptible(self, population=None, susceptibility_group=None, step_number=100, label=None):
        """
        Add to plot the trajectories of the change of the amount of susceptible individuals over time.

        :param population: Population id for the plot. By default - None.
        :type population: int or list of int or None

        :param susceptibility_group: Susceptibility group for the plot. By default - None.
        :type susceptibility_group: int or list of int or None

        :param step_number: Number of steps which split the epidemiology timeline. By default - 100.
        :type step_number: int

        :param label: The label for the susceptible plot. If label is None then to use base label. By default - None.
        :type label: str or None
        """
        if label is not None and not isinstance(label, str):
            raise TypeError('Incorrect type of lable. Type should be string or None.')

        self._create_base_plot()

        time_points = self.get_time_points(step_number)

        for data in self.get_data_susceptible(population, susceptibility_group, step_number):
            if label is None:
                label_susceptible = f'Susceptible pop: {data[0]}, sus: {data[1]}'
            else:
                label_susceptible = label

            self.ax_2.plot(time_points, data[2], label=label_susceptible)

    def add_plot_infected(self, population=None, haplotype=None, step_number=100, label=None):
        """
        Add to plot the trajectories of the change of the amount of infectious individuals over time.

        :param population: Population id for the plot. By default - None.
        :type population: int or list of int or None

        :param haplotype: Haplotype for the plot. By default - None.
        :type haplotype: int or str or list of int and str or None

        :param step_number: Number of steps which split the epidemiology timeline. By default - 100.
        :type step_number: int

        :param label: The label for the infected plot. If label is None then to use base label. By default - None.
        :type label: str or None
        """
        if label is not None and not isinstance(label, str):
            raise TypeError('Incorrect type of lable. Type should be string or None.')

        self._create_base_plot()

        time_points = self.get_time_points(step_number)

        for data in self.get_data_infected(population, haplotype, step_number):
            if label is None:
                label_infected = f'Infected pop: {data[0]}, hap: {self._calculate_string_from_number_haplotype(data[1])}'
            else:
                label_infected = label

            self.ax_2.plot(time_points, data[2], label=label_infected)

    def add_plot_sample(self, population=None, haplotype=None, step_number=100, label=None):
        """
        Add to plot the trajectories of the change of the amount of sample individuals over time.

        :param population: Population id for the plot. By default - None.
        :type population: int or list of int or None

        :param haplotype: Haplotype for the plot. By default - None.
        :type haplotype: int or str or list of int and str or None

        :param step_number: Number of steps which split the epidemiology timeline. By default - 100.
        :type step_number: int

        :param label: The label for the sample plot. If label is None then to use base label. By default - None.
        :type label: str or None
        """
        if label is not None and not isinstance(label, str):
            raise TypeError('Incorrect type of lable. Type should be string or None.')

        self._create_base_plot()

        time_points = self.get_time_points(step_number)

        for data in self.get_data_sample(population, haplotype, step_number):
            if label is None:
                label_sample = f'Sample pop: {data[0]}, hap: {self._calculate_string_from_number_haplotype(data[1])}'
            else:
                label_sample = label

            self.ax.plot(time_points, data[2], label=label_sample)

    # def add_plot_lockdowns(self, population):
    #     """
    #     Add information about lockdowns for populations

    #     :param population: populaiton id for the plot
    #     :type population: int
    #     """
    #     self._create_base_plot()
    #     lockdowns = self.get_lockdowns(population)
    #     if len(lockdowns) != 0:
    #         point = 0
    #         pointEnd = 0
    #         for ld in range(0, len(lockdowns), 2):
    #             if ld + 1 == len(lockdowns):
    #                 while time_points[point] < lockdowns[ld][1]:
    #                     point += 1
    #                 plt.fill_between(time_points[point:], infections[point:], alpha=0.2)
    #             else:
    #                 while time_points[point] < lockdowns[ld][1]:
    #                     point += 1
    #                 while time_points[pointEnd] < lockdowns[ld + 1][1]:
    #                     pointEnd += 1
    #                 if pointEnd == point:
    #                     continue
    #                 plt.fill_between(time_points[point:pointEnd + 1], infections[point:pointEnd + 1], alpha=0.2)

    def add_legend(self):
        """
        Add legend to the plot.
        """
        lines_1, labels_1 = self.ax.get_legend_handles_labels()
        lines_2, labels_2 = self.ax_2.get_legend_handles_labels()
        lines = lines_1 + lines_2
        labels = labels_1 + labels_2
        self.ax.legend(lines, labels, loc=0)

    def add_title(self, title="Plot"):
        """
        Add plot title.

        :param name: The title for the plot. By default - Plot.
        :type title: str
        """
        if not isinstance(filename, str):
            raise TypeError('Incorrect type of title. Type should be string.')

        self.ax.set_title(title)

    def plot(self, filename=None):
        """
        Show plot with all the trajectories added by other methods.

        :param filename: The filename for file with plot. If filename is None then show plot. By default - None.
        :type filename: str or None
        """
        if filename is not None and not isinstance(filename, str):
            raise TypeError('Incorrect type of filename. Type should be string or None.')

        if filename is None:
            plt.show()
        else:
            plt.savefig(filename)
        self.fig = None

    def _check_amount(self, amount, type, zero=True):
        if not isinstance(amount, int):
            raise TypeError(f'Incorrect type of {type}. Type should be int.')
        elif amount <= 0 and zero:
            raise ValueError(f'Incorrect value of {type}. Value should be more 0.')
        elif amount < 0 and not zero:
            raise ValueError(f'Incorrect value of {type}. Value should be more or equal 0.')

    def _check_migration_probability(self):
        err = self._simulator.check_migration_probability()
        if err == 1:
            raise ValueError('Incorrect the sum of migration probabilities. The sum of migration probabilities from each population should be equal or less 1.')
        elif err == 2:
            raise ValueError('Incorrect value of migration probability. Value of migration probability from source population to target population should be more 0.')

    def _check_indexes(self, indexes, edge, text, hap=False, none=True):
        if isinstance(indexes, list):
            for index in indexes:
                self._check_index(index, edge, text, hap=hap, none=none)
        else:
            self._check_index(indexes, edge, text, hap=hap, none=none)

    def _check_index(self, index, edge, text, hap=False, none=True):
        if not none and index is None:
            raise TypeError(f'Incorrect type of {text}. Type should be int.')
        elif isinstance(index, int):
            if index < 0 or index >= edge:
                raise IndexError(f'There are no such {text}!')
        elif isinstance(index, str) and hap:
            if index.count("A") + index.count("T") + index.count("C") + index.count("G") + index.count("*") \
            != self._number_of_sites:
                raise ValueError('Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"*\" and length of haplotype should be equal number of mutations sites.')
        elif index is not None:
            if hap:
                raise TypeError('Incorrect type of haplotype. Type should be int or str or None.')
            else:
                raise TypeError(f'Incorrect type of {text}. Type should be int or None.')

    def _check_list(self, data, text, length):
        if isinstance(data, list):
            if len(data) != length:
                raise ValueError(f'Incorrect length of {text}. Length should be equal {length}.')
        else:
            raise TypeError('Incorrect type of ' + text + '. Type should be list.')

    def _check_value(self, value, text, edge=None, none=False):
        if none and not isinstance(value, (int, float)) and value is not None:
            raise TypeError(f'Incorrect type of {text}. Type should be int or float or None.')
        elif not none and not isinstance(value, (int, float)):
            raise TypeError(f'Incorrect type of {text}. Type should be int or float.')

        if edge is None and isinstance(value, (int, float)) and value < 0:
            raise ValueError(f'Incorrect value of {text}. Value should be more or equal 0.')
        elif edge is not None and isinstance(value, (int, float)) and (value < 0 or value > edge):
            raise ValueError(f'Incorrect value of {text}. Value should be more or equal 0 and equal or less {edge}.')

    def _calculate_base(self, haplotype, site):
        for _ in range(self._number_of_sites-site):
            base = haplotype % 4
            haplotype = haplotype // 4
        return base

    def _calculate_indexes(self, indexes_list, edge):
        if isinstance(indexes_list, list):
            indexes = set()
            for index in indexes_list:
                indexes.update(self._calculate_index(index, edge))
        else:
            indexes = set(self._calculate_index(indexes_list, edge))
        return indexes

    def _calculate_index(self, index, edge):
        if isinstance(index, str):
            haplotypes = [index]
            bases = ["A", "T", "C", "G"]
            for s in range(self._number_of_sites):
                for i in range(len(haplotypes)):
                    haplotype_old = haplotypes[i]
                    if haplotype_old[s] == "*":
                        for j in bases:
                            index = haplotype_old.replace("*", j, 1)
                            haplotypes.append(index)
            for i in range(len(haplotypes)-1, -1, -1):
                if haplotypes[i].count("*") != 0:
                    haplotypes.remove(haplotypes[i])
            for i in range(len(haplotypes)):
                haplotypes[i] = self._calculate_number_haplotype_from_string(haplotypes[i])

            return haplotypes
        elif isinstance(index, int):
            return [index]
        elif index is None:
            return range(edge)
        else:
            raise TypeError(f'Unknown index type for calculate index - {type(index)}!')

    def _calculate_string_from_number_haplotype(self, haplotype):
        if self._number_of_sites == 0:
            return 'no sites'

        bases = ["A", "T", "C", "G"]
        string = ""
        for _ in range(self._number_of_sites):
            string = bases[haplotype%4] + string
            haplotype = haplotype // 4
        return string

    def _calculate_number_haplotype_from_string(self, string):
        string = string[::-1]
        haplotype = 0
        for site in range(self._number_of_sites):
            if string[site] == "T":
                haplotype += (4**site)
            elif string[site] == "C":
                haplotype += 2*(4**site)
            elif string[site] == "G":
                haplotype += 3*(4**site)
        return haplotype

    def _calculate_colored_haplotype(self, mutation_probabilities, haplotype, site):
        hap = self._calculate_string_from_number_haplotype(haplotype)
        base = self._calculate_base(haplotype, site)
        bases = ["A", "T", "C", "G"]
        other_base = bases[:base] + bases[base + 1:]

        color_hap=[]
        for other in other_base:
            color_hap.append(f'{hap[:site]}\033[31m{other}\033[0m{hap[site+1:]}')
        color_hap.append(f'{hap[:site]}\033[31m{bases[base]}\033[0m{hap[site+1:]}')

        result = f'{color_hap[3]}->{color_hap[0]}: {mutation_probabilities[haplotype][site][0]}\n' + \
                 f'{color_hap[3]}->{color_hap[1]}: {mutation_probabilities[haplotype][site][1]}\n' + \
                 f'{color_hap[3]}->{color_hap[2]}: {mutation_probabilities[haplotype][site][2]}\n'
        return result

    def _create_base_plot(self):
        if self.fig is None:
            self.fig, self.ax = plt.subplots(figsize=(8, 6))
            self.ax.set_ylabel('Number of samples')
            self.ax.set_xlabel('Time')
            self.ax_2 = self.ax.twinx()
            self.ax_2.set_ylabel('Number of individuals')
