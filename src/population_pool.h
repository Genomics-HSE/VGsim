#pragma once

class PopulationPool {
public:
    PopulationPool(uint64_t number_of_populations, uint64_t number_of_haplotypes, uint64_t number_of_susceptible_groups);
    ~PopulationPool();
    void Debug();
    void Update();

    void SaveInfections();
    void FirstInfections();
    void Restart();
    inline void NewInfections(uint64_t count, uint64_t haplotype, uint64_t group, uint64_t population);
    inline void NewRecoveries(uint64_t count, uint64_t haplotype, uint64_t group, uint64_t population);
    inline void NewMutation(uint64_t count, uint64_t old_haplotype, uint64_t new_haplotype, uint64_t population);
    inline void NewImmunity(uint64_t count, uint64_t old_group, uint64_t new_group, uint64_t population);
    void CheckLockdown(uint64_t population);
    std::vector<std::vector<uint64_t>> GetInfectious();

    inline uint64_t* GetSizeBegin() const;
    inline uint64_t GetInfected() const;
    inline uint64_t* GetInfectedPopBegin() const;
    inline uint64_t GetInfectedPop(uint64_t population) const;
    inline uint64_t* GetInfectedPopHapBegin(uint64_t population) const;
    inline uint64_t GetInfectedPopHap(uint64_t population, uint64_t haplotype) const;
    inline uint64_t GetSusceptibles() const;
    inline uint64_t* GetSusceptiblesPopBegin() const;
    inline uint64_t GetSusceptiblesPop(uint64_t population) const;
    inline uint64_t* GetSusceptiblesPopSusBegin(uint64_t population) const;
    inline uint64_t GetSusceptiblesPopSus(uint64_t population, uint64_t group) const;

    inline double GetContactDensity(uint64_t population) const;
    inline double* GetMaxEffectiveMigration() const;
    inline double GetMigrationProbability(uint64_t source, uint64_t target) const;
    inline double GetEffectiveMigration(uint64_t source, uint64_t target) const;
    inline double GetMultiplierMigrationProbability(uint64_t source, uint64_t target) const;
    inline double GetSamplingMultiplier(uint64_t population) const;
    inline double GetActualSize(uint64_t population) const;

private:
    void UpdateRates();

    inline uint64_t getNumberPopulations() const;
    inline uint64_t getNumberHaplotypes() const;
    inline uint64_t getNumberSusceptibleGroups() const;
    inline uint64_t getIndexPop(uint64_t first, uint64_t second) const;

    uint64_t number_of_haplotypes_;
    uint64_t number_of_populations_;
    uint64_t number_of_susceptible_groups_;

    uint64_t infected_;
    uint64_t susceptibles_;
    uint64_t* sizes_;
    uint64_t* infected_pop_;
    uint64_t* susceptibles_pop_;

    double* actual_sizes_;
    double* max_effective_migration_;
    double* migration_probability_;
    double* effective_migration_probability_;
    double* multiplier_migration_probability_;

    Population* populations_;
};
