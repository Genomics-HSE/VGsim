# use named constants for event types
DEF BIRTH = 0
DEF DEATH = 1
DEF SAMPLING = 2
DEF MUTATION = 3
DEF SUSCCHANGE = 4
DEF MIGRATION = 5
DEF MULTITYPE = 6

cdef class Event:
    cdef:
        Py_ssize_t type_, haplotype, population, newHaplotype, newPopulation
        double time

    def __init__(self, double time, Py_ssize_t type_, Py_ssize_t haplotype, Py_ssize_t population, Py_ssize_t newHaplotype, Py_ssize_t newPopulation):
        self.time = time
        self.type_ = type_
        self.haplotype = haplotype
        self.population = population
        self.newHaplotype = newHaplotype
        self.newPopulation = newPopulation


cdef class Events:
    cdef:
        Py_ssize_t size, ptr

        Py_ssize_t[::1] types, haplotypes, populations, newHaplotypes, newPopulations
        double[::1] times

    def __init__(self):
        self.size = 0
        self.ptr = 0#pointer to the first empty cell

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void AddEvent(self, double time_, Py_ssize_t type_, Py_ssize_t haplotype, Py_ssize_t population, Py_ssize_t newHaplotype, Py_ssize_t newPopulation):
        self.times[ self.ptr ] = time_
        self.types[ self.ptr ] = type_
        self.haplotypes[ self.ptr ] = haplotype
        self.populations[ self.ptr ] = population
        self.newHaplotypes[ self.ptr ] = newHaplotype
        self.newPopulations[ self.ptr ] = newPopulation
        self.ptr += 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef Event GetEvent(self, Py_ssize_t e_id):
        ev = Event( self.times[ e_id ], self.types[ e_id ], self.haplotypes[ e_id ], self.populations[ e_id ], self.newHaplotypes[ e_id ], self.newPopulations[ e_id ])
        return( ev )

    cdef void CreateEvents(self, Py_ssize_t iterations):
        if self.ptr == 0:
            self.size += iterations
            self.times = np.zeros(self.size, dtype=float)
            self.types = np.zeros(self.size, dtype=np.int64)
            self.haplotypes = np.zeros(self.size, dtype=np.int64)
            self.populations = np.zeros(self.size, dtype=np.int64)
            self.newHaplotypes = np.zeros(self.size, dtype=np.int64)
            self.newPopulations = np.zeros(self.size, dtype=np.int64)
        else:
            self.times = np.concatenate((self.times, np.zeros(iterations + self.ptr - self.size, dtype=float)))
            self.types = np.concatenate((self.types, np.zeros(iterations + self.ptr - self.size, dtype=np.int64)))
            self.haplotypes = np.concatenate((self.haplotypes, np.zeros(iterations + self.ptr - self.size, dtype=np.int64)))
            self.populations = np.concatenate((self.populations, np.zeros(iterations + self.ptr - self.size, dtype=np.int64)))
            self.newHaplotypes = np.concatenate((self.newHaplotypes, np.zeros(iterations + self.ptr - self.size, dtype=np.int64)))
            self.newPopulations = np.concatenate((self.newPopulations, np.zeros(iterations + self.ptr - self.size, dtype=np.int64)))
            self.size = iterations + self.ptr

cdef class multiEvent:
    cdef:
        Py_ssize_t num, type_, haplotype, population, newHaplotype, newPopulation
        double time

    def __init__(self, Py_ssize_t num, double time, Py_ssize_t type_, Py_ssize_t haplotype, Py_ssize_t population, Py_ssize_t newHaplotype, Py_ssize_t newPopulation):
        self.num = num
        self.time = time
        self.type_ = type_
        self.haplotype = haplotype
        self.population = population
        self.newHaplotype = newHaplotype
        self.newPopulation = newPopulation

    def PrintEvent(self):
        if self.type_ == BIRTH:
            tn = "B"
        elif  self.type_ == DEATH:
            tn = "D"
        elif  self.type_ == SAMPLING:
            tn = "S"
        elif  self.type_ == MUTATION:
            tn = "MUT"
        elif  self.type_ == SUSCCHANGE:
            tn = "SUS"
        elif  self.type_ == MIGRATION:
            tn = "MIG"
        print("num=", self.num,
                      "  time=", self.time,
                      "  type=", tn,
                      "  hap=", self.haplotype,
                      "  pop=", self.population,
                      "  dest1=", self.newHaplotype,
                      "  dest2=", self.newPopulation)

cdef class multiEvents:
    cdef:
        Py_ssize_t size, ptr

        Py_ssize_t[::1] num, types, haplotypes, populations, newHaplotypes, newPopulations
        double[::1] times

    def __init__(self):
        self.size = 0
        self.ptr = 0#pointer to the first empty cell

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void AddEvents(self, Py_ssize_t num, double time_, Py_ssize_t type_, Py_ssize_t haplotype, Py_ssize_t population, Py_ssize_t newHaplotype, Py_ssize_t newPopulation):
        self.num[ self.ptr ] = num
        self.times[ self.ptr ] = time_
        self.types[ self.ptr ] = type_
        self.haplotypes[ self.ptr ] = haplotype
        self.populations[ self.ptr ] = population
        self.newHaplotypes[ self.ptr ] = newHaplotype
        self.newPopulations[ self.ptr ] = newPopulation
        self.ptr += 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef multiEvent GetEvent(self, Py_ssize_t e_id):
        ev = multiEvent( self.num[ e_id ], self.times[ e_id ], self.types[ e_id ], self.haplotypes[ e_id ], self.populations[ e_id ], self.newHaplotypes[ e_id ], self.newPopulations[ e_id ])
        return( ev )

    cdef void CreateEvents(self, Py_ssize_t iterations):
        if self.ptr == 0:
            self.size = iterations
            self.num = np.zeros(self.size, dtype=np.int64)
            self.times = np.zeros(self.size, dtype=float)
            self.types = np.zeros(self.size, dtype=np.int64)
            self.haplotypes = np.zeros(self.size, dtype=np.int64)
            self.populations = np.zeros(self.size, dtype=np.int64)
            self.newHaplotypes = np.zeros(self.size, dtype=np.int64)
            self.newPopulations = np.zeros(self.size, dtype=np.int64)
        else:
            self.num = np.concatenate((self.num, np.zeros(iterations + self.ptr - self.size, dtype=np.int64)))
            self.times = np.concatenate((self.times, np.zeros(iterations + self.ptr - self.size, dtype=float)))
            self.types = np.concatenate((self.types, np.zeros(iterations + self.ptr - self.size, dtype=np.int64)))
            self.haplotypes = np.concatenate((self.haplotypes, np.zeros(iterations + self.ptr - self.size, dtype=np.int64)))
            self.populations = np.concatenate((self.populations, np.zeros(iterations + self.ptr - self.size, dtype=np.int64)))
            self.newHaplotypes = np.concatenate((self.newHaplotypes, np.zeros(iterations + self.ptr - self.size, dtype=np.int64)))
            self.newPopulations = np.concatenate((self.newPopulations, np.zeros(iterations + self.ptr - self.size, dtype=np.int64)))
            self.size = iterations + self.ptr
