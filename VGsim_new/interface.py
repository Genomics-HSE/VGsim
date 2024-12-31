from random import randrange
import sys

import source_VGsim

class Simulator:
    def __init__(self, number_of_sites=0, number_of_populations=1, number_of_susceptible_groups=1, seed=None):
        if seed == None:
            seed = int(randrange(sys.maxsize))
        self.simulator = source_VGsim.Simulator(number_of_sites, number_of_populations, number_of_susceptible_groups, seed)  # инициализация C++ класса
        self.number_of_sites = number_of_sites
        self.number_of_haplotypes = 4**number_of_sites
        self.number_of_populations = number_of_populations
        self.number_of_susceptible_groups = number_of_susceptible_groups

    def simulate(self, iterations=1000, sample_size=None, epidemic_time=0.0, method='direct', attempts=200):
        if sample_size is None:
            sample_size = iterations
        self.simulator.simulate(iterations, sample_size, epidemic_time, method, attempts)
        self.simulator.Debug()

    def get_flat_chain(self):
        return self.simulator.get_flat_chain()

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
        elif susceptibility_group < 0 or susceptibility_group >= self.number_of_susceptible_groups:
            raise IndexError('There are no such susceptibility group!')
        self._check_indexes(haplotype, self.number_of_haplotypes, 'haplotype', True)
        haplotypes = self._calculate_indexes(haplotype, self.number_of_haplotypes)

        for hn in haplotypes:
            self.simulator.set_susceptibility_group(susceptibility_group, hn)

    def get_susceptibility_group(self):
        return self.simulator.get_susceptibility_group()

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
        self._check_indexes(haplotype, self.number_of_haplotypes, 'haplotype', True)
        haplotypes = self._calculate_indexes(haplotype, self.number_of_haplotypes)
        for hn in haplotypes:
            self.simulator.set_transmission_rate(rate, hn)

    def get_transmission_rate(self):
        return self.simulator.get_transmission_rate()

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
        self._check_indexes(haplotype, self.number_of_haplotypes, 'haplotype', True)
        haplotypes = self._calculate_indexes(haplotype, self.number_of_haplotypes)
        for hn in haplotypes:
            self.simulator.set_recovery_rate(rate, hn)

    def get_recovery_rate(self):
        return self.simulator.get_recovery_rate()

    def set_sampling_rate(self, rate, haplotype=None):
        """
        .. _Wikipedia: https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology

        Transmission rate is the the expected number of new infections from a single infected individual per time unit in the beginning of the epidemics when all but one of the hosts are susceptible. See `Wikipedia`_ - that is parameter beta in SIR model.

        :param rate: transmission rate value.
        :type rate: float

        :param haplotype: haplotypes for which the new value is being set. `See for details <https://vg-sim.readthedocs.io/en/latest/Haplotypes.html>`_.
        :type haplotype: int(0 or 4) or string('T*' or 'AC') or int and string list([0, 4, 'T*']) or None
        """
        self._check_value(rate, 'sampling rate')
        self._check_indexes(haplotype, self.number_of_haplotypes, 'haplotype', True)
        haplotypes = self._calculate_indexes(haplotype, self.number_of_haplotypes)
        for hn in haplotypes:
            self.simulator.set_sampling_rate(rate, hn)

    def get_sampling_rate(self):
        return self.simulator.get_sampling_rate()

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
        self._check_indexes(haplotype, self.number_of_haplotypes, 'haplotype', True)
        self._check_indexes(mutation, self.number_of_sites, 'mutation site')
        haplotypes = self._calculate_indexes(haplotype, self.number_of_haplotypes)
        sites = self._calculate_indexes(mutation, self.number_of_sites)

        for hn in haplotypes:
            for s in sites:
                self.simulator.set_mutation_rate(rate, hn, s)

    def get_mutation_rate(self):
        return self.simulator.get_mutation_rate()

    def set_mutation_probabilities(self, probabilities, haplotype, mutation):
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
        self._check_indexes(haplotype, self.number_of_haplotypes, 'haplotype', True)
        self._check_indexes(mutation, self.number_of_sites, 'mutation site')
        haplotypes = self._calculate_indexes(haplotype, self.number_of_haplotypes)
        sites = self._calculate_indexes(mutation, self.number_of_sites)

        for hn in haplotypes:
            for s in sites:
                probabilities_allele = list(probabilities)
                del probabilities_allele[self._calculate_allele(hn, s)]
                if sum(probabilities_allele) == 0:
                    raise ValueError('Incorrect probabilities list. The sum of three elements without mutation allele should be more 0.')
                for i in range(3):
                    self.simulator.set_mutation_probabilities(probabilities_allele[i], hn, s, i)


    def get_mutation_probabilities(self):
        return self.simulator.get_mutation_probabilities()

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
        self._check_indexes(haplotype, self.number_of_haplotypes, 'haplotype', True)
        self._check_indexes(susceptibility_group, self.number_of_susceptible_groups, 'susceptibility group')
        haplotypes = self._calculate_indexes(haplotype, self.number_of_haplotypes)
        sus_types = self._calculate_indexes(susceptibility_group, self.number_of_susceptible_groups)

        for hn in haplotypes:
            for sn in sus_types:
                self.simulator.set_susceptibility(rate, hn, sn)

    def get_susceptibility(self):
        return self.simulator.get_susceptibility()


    def _check_value(self, value, smth, edge=None, none=False):
        if none and isinstance(value, (int, float)) == False and value is not None:
            raise TypeError('Incorrect type of ' + smth + '. Type should be int or float or None.')
        elif not none and isinstance(value, (int, float)) == False:
            raise TypeError('Incorrect type of ' + smth + '. Type should be int or float.')

        if edge is None and isinstance(value, (int, float)) and value < 0:
            raise ValueError('Incorrect value of ' + smth + '. Value should be more or equal 0.')
        elif edge is not None and isinstance(value, (int, float)) and (value < 0 or value > edge):
            raise ValueError('Incorrect value of ' + smth + '. Value should be more or equal 0 and equal or less ' + str(edge) + '.')

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
            != self.number_of_sites:
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

    def _calculate_allele(self, haplotype, site):
        for _ in range(self.number_of_sites-site):
            allele = haplotype % 4
            haplotype = haplotype // 4
        return allele

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
            for s in range(self.number_of_sites):
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

    def _calculate_haplotype_from_string(self, string):
        string = string[::-1]
        haplotype = 0
        for s in range(self.number_of_sites):
            if string[s] == "T":
                haplotype += (4**s)
            elif string[s] == "C":
                haplotype += 2*(4**s)
            elif string[s] == "G":
                haplotype += 3*(4**s)
        return haplotype
