#pragma once

#include <iostream>
#include <random>

class Direct {
public:
    Direct(Counters* counters, PopulationPool* pool, Infectious* infectious_data, Susceptibles* susceptibles_data, Chain* chain, RandomGenerator* generator, ConditionStop* stopper, Numbers numbers);
    void Debug();

    void Simulate();

private:
    void Update();
    void Restart();
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
    inline uint64_t getNumberAlleleStates() const;
    inline uint64_t getNumberHaplotypes() const;
    inline uint64_t getNumberPopulations() const;
    inline uint64_t getNumberSusceptibleGroups() const;
    inline uint64_t getIndexHap(uint64_t first, uint64_t second) const;
    inline uint64_t getIndexSus(uint64_t first, uint64_t second) const;
    inline uint64_t getIndexHap4(uint64_t first, uint64_t second, uint64_t third) const;
    inline uint64_t getIndexHapSus(uint64_t first, uint64_t second, uint64_t third) const;

    Numbers numbers_;

    double rn_;
    double rate_;
    double rate_migration_;
    ArrayBase<double> rate_pop_;
    ArrayBase<double> rate_infection_;
    ArrayBase<double> rate_immunity_;
    ArrayBase<double> rate_pop_hap_;
    ArrayBase<double> rate_pop_total_;
    ArrayBase<double> rate_pop_sus_;
    ArrayBase<double> rate_pop_hap_event_;
    ArrayBase<double> suscept_hap_pop_rate_;
    ArrayBase<double> rate_migration_pop_;
    ArrayBase<double> max_effective_transmission_migration_;

    Counters* counters_;
    PopulationPool* pool_;
    Infectious* infectious_data_;
    Susceptibles* susceptibles_data_;
    Chain* chain_;
    RandomGenerator* generator_;
    ConditionStop* stopper_;
};
