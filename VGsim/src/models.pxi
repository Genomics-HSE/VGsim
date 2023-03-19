cdef class Mutations:
    cdef:
        vector[Py_ssize_t] nodeId, AS, DS, site
        vector[double] time

    def __init__(self):#AS = ancestral state, DS = derived state
        pass

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void AddMutation(self, Py_ssize_t nodeId, Py_ssize_t haplotype, Py_ssize_t newHaplotype, double time):
        cdef:
            Py_ssize_t ASDSdigit4, site, digit4

        ASDSdigit4 = int(abs(newHaplotype - haplotype))
        site = 0
        while ASDSdigit4 >= 4:
            ASDSdigit4 = ASDSdigit4 / 4
            site += 1
        digit4 = int(4**site)
        self.nodeId.push_back(nodeId)
        self.DS.push_back(int(floor(newHaplotype/digit4) % 4))
        self.AS.push_back(int(floor(haplotype/digit4) % 4))
        self.site.push_back(int(site))
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
