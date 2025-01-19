from random import randrange
import sys
import numpy as np
from tabulate import tabulate

import source_VGsim

class Simulator:
    def __init__(self, number_of_sites=0, number_of_populations=1, number_of_susceptible_groups=1, seed=None, sampling_probability=False):
        self._check_amount(number_of_sites, 'number of sites', zero=False)
        self._check_amount(number_of_populations, 'populations number')
        self._check_amount(number_of_susceptible_groups, 'number of susceptible groups')
        if not isinstance(sampling_probability, bool):
            raise ValueError('Incorrect value of sampling probability. Value of sampling probability should be True or False.')
        if seed == None:
            seed = int(randrange(sys.maxsize))
        self._check_amount(seed, 'seed', zero=False)

        self._simulator = source_VGsim.Simulator(number_of_sites, number_of_populations, number_of_susceptible_groups, seed)  # инициализация C++ класса
        self._number_of_sites = number_of_sites
        self._number_of_haplotypes = 4**number_of_sites
        self._number_of_populations = number_of_populations
        self._number_of_susceptible_groups = number_of_susceptible_groups
        self._seed = seed
        self._sampling_probability = sampling_probability
        self._first_simulation = False

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
        if isinstance(susceptibility_group, int) == False:
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
        self._check_list(probabilities, 'probabilities list', 4)
        if isinstance(probabilities, list):
            for i in range(4):
                self._check_value(probabilities[i], 'mutation probabilities')
        self._check_indexes(haplotype, self._number_of_haplotypes, 'haplotype', True)
        self._check_indexes(mutation, self._number_of_sites, 'mutation site')
        haplotypes = self._calculate_indexes(haplotype, self._number_of_haplotypes)
        sites = self._calculate_indexes(mutation, self._number_of_sites)

        for hn in haplotypes:
            for s in sites:
                probabilities_allele = list(probabilities)
                del probabilities_allele[self._calculate_allele(hn, s)]
                if sum(probabilities_allele) == 0:
                    raise ValueError('Incorrect probabilities list. The sum of three elements without mutation allele should be more 0.')
                for i in range(3):
                    self._simulator.set_mutation_probabilities(probabilities_allele[i], hn, s, i)


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
            # self.sizes[pn] = amount
            # self.susceptible[pn, 0] = amount
            # for sn in range(1, self.susNum):
            #     self.susceptible[pn, sn] = 0

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
            # self.contactDensity[pn] = value
            # self.contactDensityBeforeLockdown[pn] = value

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
            # self.contactDensityAfterLockdown[pn] = parameters[0]
            # self.startLD[pn] = parameters[1]
            # self.endLD[pn] = parameters[2]

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
        data = [field]

        transmission_rate = self.get_transmission_rate()
        recovery_rate = self.get_recovery_rate()
        sampling_rate = self.get_sampling_rate()
        susceptibility_group = self.get_susceptibility_group()
        mutation_rate = self.get_mutation_rate()
        mutation_probabilities = self.get_mutation_probabilities()

        for hn in range(self._number_of_haplotypes):
            row = [
                    self._calculate_string_from_haplotype(hn),
                    transmission_rate[hn],
                    recovery_rate[hn],
                    sampling_rate[hn],
                    susceptibility_group[hn]
                   ]
            for s in range(self._number_of_sites):
                row.append(mutation_rate[hn][s])
                row.append(self._calculate_colored_haplotype(mutation_probabilities, hn, s))
            data.append(row)

        print(tabulate(data, headers="firstrow", tablefmt="fancy_grid", numalign="center", stralign="center"))
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
            print(current_susceptible, current_infectious)

        # if population:
        #     table_populations = PrettyTable()
        #     table_populations.field_names = ["ID", "Size", 'Actual size', "CD",'CDBLC', "CDALD", "SLD", "ELD", "SM"]
        #     for pn in range(self.popNum):
        #         table_populations.add_row([pn, self.sizes[pn], self.actualSizes[pn], self.contactDensity[pn], \
        #             self.contactDensityBeforeLockdown[pn], self.contactDensityAfterLockdown[pn], self.startLD[pn], self.endLD[pn], \
        #             self.samplingMultiplier[pn]])

        #     print(table_populations)
        #     print("Legend:")
        #     print("ID - number of population")
        #     print("Size - size of population")
        #     print("Actual size - actual size of population")
        #     print("CD - contact density")
        #     print("CDBLD - contact density without lockdown")
        #     print("CDALD - contact density at lockdown")
        #     print("SLD - start of lockdown")
        #     print("ELD - end of lockdown")
        #     print("SM - sampling multiplier")
        #     print()

        # if susceptibles:
        #     table_susceptible = PrettyTable()
        #     field = ["ST\\ID"]
        #     for pn in range(self.popNum):
        #         field.append(pn)
        #     table_susceptible.field_names = field
        #     for sn in range(self.susNum):
        #         row = [sn]
        #         for pn in range(self.popNum):
        #             row.append(current_susceptible[pn][sn])
        #         table_susceptible.add_row(row)

        #     print(table_susceptible)
        #     print("Legend:")
        #     print("ID - ID population")
        #     print("ST - susceptibility type")
        #     print()

        # if infectious:
        #     table_infectious = PrettyTable()
        #     field = ["H\\ID"]
        #     for pn in range(self.popNum):
        #         field.append(pn)
        #     table_infectious.field_names = field
        #     for hn in range(self.hapNum):
        #         row = [self.calculate_string_from_haplotype(hn)]
        #         for pn in range(self.popNum):
        #             row.append(current_infectious[pn][hn])
        #         table_infectious.add_row(row)

        #     print(table_infectious)
        #     print("Legend:")
        #     print("ID - ID population")
        #     print("H - haplotype")
        #     print()

        # if migration:
        #     table_migration = PrettyTable()
        #     field = ["S\\T"]
        #     for pn in range(self.popNum):
        #         field.append(pn)
        #     table_migration.field_names = field
        #     for pn1 in range(self.popNum):
        #         row = [pn1]
        #         for pn2 in range(self.popNum):
        #             row.append(self.migrationRates[pn1, pn2])
        #         table_migration.add_row(row)

        #     print(table_migration)
        #     print("Legend:")
        #     print("S - ID source population")
        #     print("T - ID target population")
        #     print()

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
            data = [field]

            susceptibility = self.get_susceptibility()

            for hn in range(self._number_of_haplotypes):
                row = [self._calculate_string_from_haplotype(hn)]
                for sn in range(self._number_of_susceptible_groups):
                    row.append(susceptibility[hn][sn])
                data.append(row)

            print(tabulate(data, headers="firstrow", tablefmt="fancy_grid", numalign="center", stralign="center"))
            print("Legend:")
            print("H - haplotype")
            print("ST - susceptibility type")
            print()

        if transition:
            field = ["ID"]
            for sn in range(self._number_of_susceptible_groups):
                field.append(sn)
            data = [field]

            immunity_transition = self.get_immunity_transition()

            for sn1 in range(self._number_of_susceptible_groups):
                row = [sn1]
                for sn2 in range(self._number_of_susceptible_groups):
                    row.append(immunity_transition[sn1][sn2])
                data.append(row)

            print(tabulate(data, headers="firstrow", tablefmt="fancy_grid", numalign="center", stralign="center"))
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

    def _check_amount(self, amount, smth, zero=True):
        if isinstance(amount, int) == False:
            raise TypeError('Incorrect type of ' + smth + '. Type should be int.')
        elif amount <= 0 and zero:
            raise ValueError('Incorrect value of ' + smth + '. Value should be more 0.')
        elif amount < 0 and zero == False:
            raise ValueError('Incorrect value of ' + smth + '. Value should be more or equal 0.')

    def _check_migration_probability(self):
        err = self._simulator.check_migration_probability()
        if err == 1:
            raise ValueError('Incorrect the sum of migration probabilities. The sum of migration probabilities from each population should be equal or less 1.')
        elif err == 2:
            raise ValueError('Incorrect value of migration probability. Value of migration probability from source population to target population should be more 0.')

    def _check_indexes(self, index, edge, smth, hap=False, none=True):
        if isinstance(index, list):
            for i in index:
                self._check_index(i, edge, smth, hap=hap, none=none)
        else:
            self._check_index(index, edge, smth, hap=hap, none=none)

    def _check_index(self, index, edge, smth, hap=False, none=True):
        if none == False and index is None:
            raise TypeError('Incorrect type of ' + smth + '. Type should be int.')
        elif isinstance(index, int):
            if index < 0 or index >= edge:
                raise IndexError('There are no such ' + smth + '!')
        elif isinstance(index, str) and hap:
            if index.count("A") + index.count("T") + index.count("C") + index.count("G") + index.count("*") \
            != self._number_of_sites:
                raise ValueError('Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"*\" and length of haplotype should be equal number of mutations sites.')
        elif index is not None:
            if hap:
                raise TypeError('Incorrect type of haplotype. Type should be int or str or None.')
            else:
                raise TypeError('Incorrect type of ' + smth + '. Type should be int or None.')

    def _check_list(self, data, smth, length):
        if isinstance(data, list):
            if len(data) != length:
                raise ValueError('Incorrect length of ' + smth + '. Length should be equal ' + str(length) + '.')
        else:
            raise TypeError('Incorrect type of ' + smth + '. Type should be list.')

    def _check_value(self, value, smth, edge=None, none=False):
        if none and isinstance(value, (int, float)) == False and value is not None:
            raise TypeError('Incorrect type of ' + smth + '. Type should be int or float or None.')
        elif not none and isinstance(value, (int, float)) == False:
            raise TypeError('Incorrect type of ' + smth + '. Type should be int or float.')

        if edge is None and isinstance(value, (int, float)) and value < 0:
            raise ValueError('Incorrect value of ' + smth + '. Value should be more or equal 0.')
        elif edge is not None and isinstance(value, (int, float)) and (value < 0 or value > edge):
            raise ValueError('Incorrect value of ' + smth + '. Value should be more or equal 0 and equal or less ' + str(edge) + '.')

    def _calculate_allele(self, haplotype, site):
        for _ in range(self._number_of_sites-site):
            allele = haplotype % 4
            haplotype = haplotype // 4
        return allele

    def _calculate_haplotype_from_string(self, string):
        string = string[::-1]
        haplotype = 0
        for s in range(self._number_of_sites):
            if string[s] == "T":
                haplotype += (4**s)
            elif string[s] == "C":
                haplotype += 2*(4**s)
            elif string[s] == "G":
                haplotype += 3*(4**s)
        return haplotype

    def _calculate_indexes(self, indexes_list, edge):
        if isinstance(indexes_list, list):
            indexes = set()
            for i in indexes_list:
                indexes.update(self._calculate_index(i, edge))
        else:
            indexes = set(self._calculate_index(indexes_list, edge))
        return indexes

    def _calculate_index(self, index, edge):
        if isinstance(index, str):
            haplotypes = [index]
            letters = ["A", "T", "C", "G"]
            for s in range(self._number_of_sites):
                for i in range(len(haplotypes)):
                    haplotype_old = haplotypes[i]
                    if haplotype_old[s] == "*":
                        for j in letters:
                            index = haplotype_old.replace("*", j, 1)
                            haplotypes.append(index)
            for i in range(len(haplotypes)-1, -1, -1):
                if haplotypes[i].count("*") != 0:
                    haplotypes.remove(haplotypes[i])
            for i in range(len(haplotypes)):
                haplotypes[i] = self._calculate_haplotype_from_string(haplotypes[i])

            return haplotypes
        elif isinstance(index, int):
            return [index]
        else:
            return range(edge)

    def _calculate_colored_haplotype(self, mutation_probabilities, haplotype, site):
        hap = self._calculate_string_from_haplotype(haplotype)
        haplotypes = [hap[:site] + "A" + hap[site+1:], hap[:site] + "T" + hap[site+1:], hap[:site] + "C" + hap[site+1:], \
        hap[:site] + "G" + hap[site+1:]]
        for i in range(4):
            if haplotypes[i] == hap:
                haplotypes.remove(haplotypes[i])
                haplotypes.append(hap)
                break
        color_hap=[]
        for hapl in haplotypes:
            a = ""
            for s in range(self._number_of_sites):
                if s == site:
                    a = a + "\033[31m{}\033[0m" .format(hapl[s:s+1])
                else:
                    a = a + hapl[s:s+1]
            color_hap.append(a)
        string = color_hap[3] + "->" + color_hap[0] + ": " + str(mutation_probabilities[haplotype][site][0]) + "\n" + color_hap[3] + \
        "->" + color_hap[1] + ": " + str(mutation_probabilities[haplotype][site][1]) + "\n" + color_hap[3] + "->" + color_hap[2] + ": " + \
        str(mutation_probabilities[haplotype][site][2]) + "\n"
        return string

    def _calculate_string_from_haplotype(self, haplotype):
        letters = ["A", "T", "C", "G"]
        string = ""
        for s in range(self._number_of_sites):
            string = letters[haplotype%4] + string
            haplotype = haplotype // 4
        return string
