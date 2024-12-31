import source_VGsim

class Simulator:
    def __init__(self, number_of_sites=0, populations_number=1, number_of_susceptible_groups=1, seed=None):
        if seed == None:
            seed = int(randrange(sys.maxsize))
        self.simulator = source_VGsim.Simulator(number_of_sites, populations_number, number_of_susceptible_groups, seed)  # инициализация C++ класса
        self.hapNum = 4**number_of_sites
        self.sites = number_of_sites

    def simulate(self, iterations=1000, sample_size=None, epidemic_time=0.0, method='direct', attempts=200):
        if sample_size is None:
            sample_size = iterations
        self.simulator.simulate(iterations, sample_size, epidemic_time, method, attempts)
        self.simulator.Debug()

    def get_flat_chain(self):
        return self.simulator.get_flat_chain()

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
        self._check_indexes(haplotype, self.hapNum, 'haplotype', True)
        haplotypes = self._calculate_indexes(haplotype, self.hapNum)
        for hn in haplotypes:
            self.simulator.set_transmission_rate(rate, hn)

    def get_transmission_rate(self):
        return self.simulator.get_transmission_rate()

    def _check_value(self, value, smth, edge=None, none=False):
        if none and isinstance(value, (int, float)) == False and value is not None:
            raise TypeError('Incorrect type of ' + smth + '. Type should be int or float or None.')
        elif none and isinstance(value, (int, float)) == False:
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
            != self.sites:
                raise ValueError('Incorrect haplotype. Haplotype should contain only \"A\", \"T\", \"C\", \"G\", \"*\" and length of haplotype should be equal number of mutations sites.')
        elif index is not None:
            if hap:
                raise TypeError('Incorrect type of haplotype. Type should be int or str or None.')
            else:
                raise TypeError('Incorrect type of ' + smth + '. Type should be int or None.')

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
            for s in range(self.sites):
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
        for s in range(self.sites):
            if string[s] == "T":
                haplotype += (4**s)
            elif string[s] == "C":
                haplotype += 2*(4**s)
            elif string[s] == "G":
                haplotype += 3*(4**s)
        return haplotype
