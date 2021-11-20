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


cdef class PopulationModel:
    cdef:
        bint strong_migration
        Py_ssize_t hapNum, popNum, susNum, swapLockdown, migPlus, migNonPlus, migCounter, globalInfectious

        Migrations mig

        long[::1] sizes, totalSusceptible, totalInfectious, lockdownON
        long[:,::1] susceptible, liveBranches

        double[::1] effectiveSizes, contactDensity, contactDensityBeforeLockdown, contactDensityAfterLockdown, startLD, endLD, samplingMultiplier, maxEffectiveMigration
        double[:,::1] migrationRates, effectiveMigration

    def __init__(self, popNum, susNum, hapNum, strong_migration):
        self.hapNum = hapNum
        self.popNum = popNum
        self.susNum = susNum
        self.strong_migration = strong_migration

        self.swapLockdown = 0
        self.migPlus = 0
        self.migNonPlus = 0
        self.migCounter = 0
        self.globalInfectious = 0

        self.mig = Migrations()

        self.sizes = np.zeros(self.popNum, dtype=np.int64)
        self.totalSusceptible = np.zeros(self.popNum, dtype=np.int64)
        self.totalInfectious = np.zeros(self.popNum, dtype=np.int64)
        self.lockdownON = np.zeros(self.popNum, dtype=np.int64)

        self.susceptible = np.zeros((self.popNum, self.susNum), dtype=np.int64)
        self.liveBranches = np.zeros((self.popNum, self.hapNum), dtype=np.int64)

        self.effectiveSizes = np.zeros(self.popNum, dtype=float)
        self.contactDensity = np.ones(self.popNum, dtype=float)
        self.contactDensityBeforeLockdown = np.ones(self.popNum, dtype=float)
        self.contactDensityAfterLockdown = np.zeros(self.popNum, dtype=float)
        self.startLD = np.ones(self.popNum, dtype=float)
        self.endLD = np.ones(self.popNum, dtype=float)
        self.samplingMultiplier = np.ones(self.popNum, dtype=float)
        self.maxEffectiveMigration = np.zeros(self.popNum, dtype=float)

        self.migrationRates = np.zeros((self.popNum, self.popNum), dtype=float)
        self.effectiveMigration = np.zeros((self.popNum, self.popNum), dtype=float)

        for pn in range(self.popNum):
            self.sizes[pn] = 1000000
            self.totalSusceptible[pn] = 1000000
            self.susceptible[pn, 0] = 1000000

        # self.NewInfection(0, 0, 0)
        self.UpdateAllRates()

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.cdivision(True)
    cdef void UpdateAllRates(self):
        for pn in range(self.popNum):
            self.maxEffectiveMigration[pn] = 0.0

        if self.strong_migration == True:
            for pn1 in range(self.popNum):
                for pn2 in range(self.popNum):
                    if pn1 == pn2:
                        summa = 1
                        for pn3 in range(self.popNum):
                            summa -= self.migrationRates[pn3, pn1]
                        self.effectiveSizes[pn1] = summa*self.sizes[pn1]
                    else:
                        self.effectiveSizes[pn1] = self.migrationRates[pn1, pn2]*self.sizes[pn1]
        else:
            for pn in range(self.popNum):
                self.effectiveSizes[pn] = self.sizes[pn]

        for pn1 in range(self.popNum):
            for pn2 in range(self.popNum):
                self.effectiveMigration[pn1, pn2] = self.migrationRates[pn1, pn2]*self.contactDensity[pn2]/self.effectiveSizes[pn2]+self.migrationRates[pn2, pn1]*self.contactDensity[pn1]/self.effectiveSizes[pn1]
                if self.effectiveMigration[pn1, pn2] > self.maxEffectiveMigration[pn2]:
                    self.maxEffectiveMigration[pn2] = self.effectiveMigration[pn1, pn2]

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef void FirstInfection(self):
        if self.globalInfectious == 0:
            for sn in range(self.susNum):
                if self.susceptible[0, sn] != 0:
                    self.NewInfection(0, sn, 0)
                    break

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef inline void NewInfection(self, Py_ssize_t pi, Py_ssize_t si, Py_ssize_t hi):
        self.susceptible[pi, si] -= 1
        self.totalSusceptible[pi] -= 1
        self.liveBranches[pi, hi] += 1
        self.totalInfectious[pi] += 1
        self.globalInfectious += 1

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef inline void NewRecovery(self, Py_ssize_t pi, Py_ssize_t si, Py_ssize_t hi):
        self.susceptible[pi, si] += 1
        self.totalSusceptible[pi] += 1
        self.liveBranches[pi, hi] -= 1
        self.totalInfectious[pi] -= 1
        self.globalInfectious -= 1

    def debug(self):
        print("Population model")
        print("Haplotypes number(const): ", self.hapNum)
        print("Populations number(const): ", self.popNum)
        print("Susceptible number(const): ", self.susNum)
        print("Strong migration(const): ", self.strong_migration)
        print("swapLockdown(mutable): ", self.swapLockdown)
        print("Migration plus(mutable): ", self.migPlus)
        print("Migration non plus(mutable): ", self.migNonPlus)
        print("Migration counter(mutable): ", self.migCounter)
        print("globalInfectious(mutable): ", self.globalInfectious)

        print("sizes(const): ", end="")
        for pn in range(self.popNum):
            print(self.sizes[pn], end=" ")
        print()
        print("totalSusceptible(mutable): ", end="")
        for pn in range(self.popNum):
            print(self.totalSusceptible[pn], end=" ")
        print()
        print("totalInfectious(mutable): ", end="")
        for pn in range(self.popNum):
            print(self.totalInfectious[pn], end=" ")
        print()
        print("lockdownON(mutable): ", end="")
        for pn in range(self.popNum):
            print(self.lockdownON[pn], end=" ")
        print()
        print("effectiveSizes(const): ", end=" ")
        for pn in range(self.popNum):
            print(self.effectiveSizes[pn], end=" ")
        print()
        print("contact density(const): ", end=" ")
        for pn in range(self.popNum):
            print(self.contactDensity[pn], end=" ")
        print()
        print("contactDensityBeforeLockdown(const): ", end=" ")
        for pn in range(self.popNum):
            print(self.contactDensityBeforeLockdown[pn], end=" ")
        print()
        print("contactDensityBeforeLockdown(const): ", end=" ")
        for pn in range(self.popNum):
            print(self.contactDensityAfterLockdown[pn], end=" ")
        print()
        print("startLD(const): ", end=" ")
        for pn in range(self.popNum):
            print(self.startLD[pn], end=" ")
        print()
        print("endLD(const): ", end=" ")
        for pn in range(self.popNum):
            print(self.endLD[pn], end=" ")
        print()
        print("samplingMultiplier(const): ", end=" ")
        for pn in range(self.popNum):
            print(self.samplingMultiplier[pn], end=" ")
        print()
        print("max effective migration(const): ", end=" ")
        for pn in range(self.popNum):
            print(self.maxEffectiveMigration[pn], end=" ")
        print()

        print("susceptible(mutable)----")
        for pn in range(self.popNum):
            for sn in range(self.susNum):
                print(self.susceptible[pn, sn], end=" ")
            print()
        print()
        print("liveBranches(mutable)----")
        for pn in range(self.popNum):
            for hn in range(self.hapNum):
                print(self.liveBranches[pn, hn], end=" ")
            print()
        print()
        print("Population model - migration rates(const)----")
        for pn1 in range(self.popNum):
            for pn2 in range(self.popNum):
                print(self.migrationRates[pn1, pn2], end=" ")
            print()
        print()
        print("Population model - effective migration(const)----")
        for pn1 in range(self.popNum):
            for pn2 in range(self.popNum):
                print(self.effectiveMigration[pn1, pn2], end=" ")
            print()
        print()
