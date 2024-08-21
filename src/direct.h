#pragma once

#include <iostream>
#include <random>

class Direct {
public:
    Direct(Counters* counters, PopulationPool* pool, Infectious* infectious_data, Susceptibles* susceptibles_data, Chain* chain, RandomGenerator* generator, uint64_t sites, uint64_t haplotypes, uint64_t populations, uint64_t susceptible_groups);
    ~Direct();
    void Debug();

    void Simulate(uint64_t iterations = 100'000, uint64_t number_attempts = 100);

private:
    void Update();
    void Restart(uint64_t index);
    void TimeStep();
    void GenerateEvent();
    void UpdateRates(uint64_t population, bool infection, bool immunity, bool migration);
    double BirthRate(uint64_t population, uint64_t haplotype);
    void ImmunityTransition(uint64_t population);
    void Transmission(uint64_t population, uint64_t haplotype);
    void Recovery(uint64_t population, uint64_t haplotype);
    void Sampling(uint64_t population, uint64_t haplotype);
    void Mutation(uint64_t population, uint64_t haplotype);
    void Migration();

    inline uint64_t getNumberSites() const;
    inline uint64_t getNumberHaplotypes() const;
    inline uint64_t getNumberPopulations() const;
    inline uint64_t getNumberSusceptibleGroups() const;
    inline uint64_t getIndexHap(uint64_t first, uint64_t second) const;
    inline uint64_t getIndexSus(uint64_t first, uint64_t second) const;
    inline uint64_t getIndexHap4(uint64_t first, uint64_t second, uint64_t third) const;
    inline uint64_t getIndexHapSus(uint64_t first, uint64_t second, uint64_t third) const;

    uint64_t number_of_sites_;
    uint64_t number_of_haplotypes_;
    uint64_t number_of_populations_;
    uint64_t number_of_susceptible_groups_;

    double rn_;
    double rate_;
    double rate_migration_;
    double* rate_pop_;
    double* rate_infection_;
    double* rate_immunity_;
    double* rate_pop_hap_;
    double* rate_pop_total_;
    double* rate_pop_sus_;
    double* rate_pop_hap_event_;
    double* suscept_hap_pop_rate_;
    double* rate_migration_pop_;
    double* max_effective_transmission_migration_;

    Counters* counters_;
    PopulationPool* pool_;
    Infectious* infectious_data_;
    Susceptibles* susceptibles_data_;
    Chain* chain_;
    RandomGenerator* generator_;
};
