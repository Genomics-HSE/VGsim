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

        def GetData(self):
            return {'node': self.nodeId, 'AS': self.AS, 'DS': self.DS, 'site': self.site, 'time': self.time}
            # print("MutType, AS, DS: ", site, self.AS[self.AS.size()-1], self.DS[self.DS.size()-1])

cdef class Migrations:
    cdef:
        vector[Py_ssize_t] nodeId, oldPop, newPop
        vector[double] time

    def __init__(self):
        pass

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void AddMigration(self, Py_ssize_t nodeId, double time, Py_ssize_t oldPop, Py_ssize_t newPop):
        self.nodeId.push_back(nodeId)
        self.time.push_back(time)
        self.oldPop.push_back(oldPop)
        self.newPop.push_back(newPop)

class Population:
    def __init__(self, size=1000000, contactDensity=1.0):
        self.size = size
        self.contactDensity = contactDensity

class Lockdown:
    def __init__(self, conDenAfterLD=0.1, startLD=2, endLD=1):
        self.conDenAfterLD = conDenAfterLD
        self.startLD = startLD
        self.endLD = endLD


cdef class PopulationModel:
    cdef:
        Py_ssize_t globalInfectious, migPlus, migNonPlus, swapLockdown, migCounter, popNum

        Migrations mig

        long[::1] sizes, totalSusceptible, totalInfectious, lockdownON
        long[:,::1] susceptible

        double[::1] maxEffectiveMigration, contactDensity, contactDensityBeforeLockdown, contactDensityAfterLockdown, startLD, endLD, samplingMultiplier
        double[:,::1] migrationRates, effectiveMigration

    def __init__(self, population_sizes, susNum):
        self.popNum = len(population_sizes)
        self.sizes = np.asarray(population_sizes)
        self.susceptible = np.zeros((self.popNum, susNum), dtype=np.int64)
        self.totalSusceptible = np.zeros(self.popNum, dtype=np.int64)
        self.totalInfectious = np.zeros(self.popNum, dtype=np.int64)
        self.globalInfectious = 0

        self.contactDensity = np.ones(self.popNum, dtype=float)
        self.contactDensityBeforeLockdown = np.ones(self.popNum, dtype=float)
        self.contactDensityAfterLockdown = np.zeros(self.popNum, dtype=float)
        self.startLD = np.zeros(self.popNum, dtype=float)
        self.endLD = np.zeros(self.popNum, dtype=float)
        self.lockdownON = np.zeros(self.popNum, dtype=np.int64)
        self.swapLockdown = 0

        self.samplingMultiplier = np.ones(self.popNum)
        for pn in range(self.popNum):
            self.susceptible[pn, 0] = self.sizes[pn]
            self.totalSusceptible[pn] = self.sizes[pn]
            self.startLD[pn] = 1.01*self.sizes[pn]
            self.endLD[pn] = 1.0*self.sizes[pn]

        self.mig = Migrations()
        self.migrationRates = np.zeros((self.popNum, self.popNum), dtype=float)
        self.effectiveMigration = np.zeros((self.popNum, self.popNum), dtype=float)
        self.maxEffectiveMigration = np.zeros(self.popNum, dtype=float)
        self.migPlus = 0
        self.migNonPlus = 0
        self.migCounter = 0

        self.SetEffectiveMigration()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void SetEffectiveMigration(self):
        for pn in range(self.popNum):
            self.maxEffectiveMigration[pn] = 0.0
        for pn1 in range(self.popNum):
            for pn2 in range(self.popNum):
                self.effectiveMigration[pn1, pn2] = self.migrationRates[pn1, pn2]*self.contactDensity[pn2]/self.sizes[pn2]+self.migrationRates[pn2, pn1]*self.contactDensity[pn1]/self.sizes[pn1]
                if self.effectiveMigration[pn1, pn2] > self.maxEffectiveMigration[pn2]:
                    self.maxEffectiveMigration[pn2] = self.effectiveMigration[pn1, pn2]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef inline void NewInfection(self, Py_ssize_t pi, Py_ssize_t si):
        self.susceptible[pi, si] -= 1
        self.totalSusceptible[pi] -= 1
        self.totalInfectious[pi] += 1
        self.globalInfectious += 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef inline void NewRecovery(self, Py_ssize_t pi, Py_ssize_t si):
        self.susceptible[pi, si] += 1
        self.totalSusceptible[pi] += 1
        self.totalInfectious[pi] -= 1
        self.globalInfectious -= 1
