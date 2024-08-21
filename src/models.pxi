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


cdef class RecombinationForward:
    cdef:
        double time
        Py_ssize_t position

    def __init__(self, double time, Py_ssize_t position):
        self.time = time
        self.position = position

# cdef class RecombinationBackward:
#     cdef:
#         Py_ssize_t child, parent, position, haplotype_left, haplotype_right

#     def __init__(self, Py_ssize_t child, Py_ssize_t parent, Py_ssize_t position, Py_ssize_t left, Py_ssize_t right):
#         self.child = child
#         self.parent = parent
#         self.position = position
#         self.haplotype_left = left
#         self.haplotype_right = right

cdef class Recombination:
    cdef:
        Py_ssize_t counter

        vector[double] time
        vector[Py_ssize_t] position

        vector[Py_ssize_t] child
        vector[Py_ssize_t] parent
        vector[Py_ssize_t] position_back
        vector[Py_ssize_t] left_haplotype
        vector[Py_ssize_t] right_haplotype

    def __init__(self):
        self.counter = 0

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void addForward(self, double time, Py_ssize_t position):
        self.time.push_back(time)
        self.position.push_back(position)
        self.counter += 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef Py_ssize_t getForward(self):
        self.counter -= 1
        # result = RecombinationForward(self.time[self.counter], self.position[self.counter])
        # return result
        return self.position[self.counter]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef float getPosition(self, double time):
        for i in range(self.time.size()):
            if self.time[i].time == time:
                return self.positions[i]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void addBackward(self, Py_ssize_t child, Py_ssize_t parent, Py_ssize_t position, Py_ssize_t left, Py_ssize_t right):
        self.child.push_back(child)
        self.parent.push_back(parent)
        self.position.push_back(position)
        self.left_haplotype.push_back(left)
        self.right_haplotype.push_back(right)

    def print_edges(self):
        for i in range(self.child.size()):
            print(f"Edge - {i}, child - {self.child[i]}, parent - {self.parent[i]}")



    # @cython.boundscheck(False)
    # @cython.wraparound(False)
    # cdef void AddRecombination_forward(self, Py_ssize_t Id, Py_ssize_t hi, Py_ssize_t hi2, Py_ssize_t posRecomb, Py_ssize_t nhi):
    #     self.idevents.push_back(Id)
    #     self.his.push_back(hi)
    #     self.hi2s.push_back(hi2)
    #     self.nhis.push_back(nhi)
    #     self.posRecombs.push_back(posRecomb)

    # @cython.boundscheck(False)
    # @cython.wraparound(False)
    # cdef void AddRecombination(self, Py_ssize_t position):
    # # cdef void AddRecombination(self, Py_ssize_t position, Py_ssize_t hi1, Py_ssize_t hi2, double time):
    #     self.positions.push_back(position)
    #     self.hi1s.push_back(hi1)
    #     self.hi2s.push_back(hi2)
    #     self.times.push_back(time)
