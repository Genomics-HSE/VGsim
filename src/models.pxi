cdef class Mutations:
    cdef:
        vector[Py_ssize_t] nodeId, AS, DS, site
        vector[double] time

        Py_ssize_t number_of_sites, number_of_states_allele

    def __init__(self, number_of_sites, number_of_states_allele):#AS = ancestral state, DS = derived state
        self.number_of_sites = number_of_sites
        self.number_of_states_allele = number_of_states_allele

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void AddMutation(self, Py_ssize_t nodeId, Py_ssize_t haplotype, Py_ssize_t newHaplotype, double time):
        cdef:
            Py_ssize_t ASDSdigit4, site, digit4

        ASDSdigit4 = int(abs(newHaplotype - haplotype))
        site = self.number_of_sites - 1

        while ASDSdigit4 >= self.number_of_states_allele:
            ASDSdigit4 /= self.number_of_states_allele
            site -= 1

        digit4 = int(self.number_of_states_allele**(self.number_of_sites - site - 1))
        self.nodeId.push_back(nodeId)
        self.DS.push_back(newHaplotype // digit4 % self.number_of_states_allele)
        self.AS.push_back(haplotype // digit4 % self.number_of_states_allele)
        self.site.push_back(site)
        self.time.push_back(time)

    def get_mutation(self, id_mut):
        return self.nodeId[id_mut], self.DS[id_mut], self.AS[id_mut], self.site[id_mut], self.time[id_mut]


cdef class Migrations:
    cdef:
        vector[Py_ssize_t] nodeId, oldPop, newPop
        vector[double] time

    def __init__(self):
        pass

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void AddMigration(self, Py_ssize_t nodeId, double time, Py_ssize_t oldPop, Py_ssize_t newPop):
        self.nodeId.push_back(nodeId)
        self.time.push_back(time)
        self.oldPop.push_back(oldPop)
        self.newPop.push_back(newPop)

    def get_migration(self, id_mig):
        return self.nodeId[id_mig], self.time[id_mig], self.oldPop[id_mig], self.newPop[id_mig]


cdef class Lockdowns:
    cdef:
        vector[bint] states
        vector[Py_ssize_t] populationsId
        vector[double] times

    def __init__(self):
        pass

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void AddLockdown(self, bint state, Py_ssize_t populationId, double time):
        self.states.push_back(state)
        self.populationsId.push_back(populationId)
        self.times.push_back(time)


cdef class Recombination:
    cdef:
        vector[Py_ssize_t] positions
        # vector[Py_ssize_t] positions, hi1s, hi2s
        # vector[double] times

        vector[Py_ssize_t] idevents, his, hi2s, nhis, posRecombs

    def __init__(self):
        pass

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void AddRecombination_forward(self, Py_ssize_t Id, Py_ssize_t hi, Py_ssize_t hi2, Py_ssize_t posRecomb, Py_ssize_t nhi):
        self.idevents.push_back(Id)
        self.his.push_back(hi)
        self.hi2s.push_back(hi2)
        self.nhis.push_back(nhi)
        self.posRecombs.push_back(posRecomb)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void AddRecombination(self, Py_ssize_t position):
    # cdef void AddRecombination(self, Py_ssize_t position, Py_ssize_t hi1, Py_ssize_t hi2, double time):
        self.positions.push_back(position)
        # self.hi1s.push_back(hi1)
        # self.hi2s.push_back(hi2)
        # self.times.push_back(time)
