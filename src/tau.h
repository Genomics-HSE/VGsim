#pragma once

#include <iostream>
#include <random>

class Tau {
public:
    Tau(Counters* counters, PopulationPool* pool, Infectious* infectious_data, Susceptibles* susceptibles_data, Chain* chain, RandomGenerator* generator, Numbers numbers);
    ~Tau();
    void Debug();

    void Simulate(uint64_t iterations = 1'000, uint64_t number_attempts = 10);

private:
    void Update();
    void Propensities();
    void TimeStep();
    void ChooseTime();
    bool GenerateEvents();
    void UpdateCompartmentCounts();

    inline uint64_t getNumberSites() const;
    inline uint64_t getNumberHaplotypes() const;
    inline uint64_t getNumberPopulations() const;
    inline uint64_t getNumberSusceptibleGroups() const;
    inline uint64_t getIndexHap(uint64_t first, uint64_t second) const;
    inline uint64_t getIndexSus(uint64_t first, uint64_t second) const;
    inline uint64_t getIndexSit(uint64_t first, uint64_t second) const;

    Numbers numbers_;

    int64_t* infectious_delta_;
    int64_t* susceptibles_delta_;

    uint64_t* events_transmission_;
    uint64_t* events_recovery_;
    uint64_t* events_sampling_;
    uint64_t* events_mutation_;
    uint64_t* events_migration_;
    uint64_t* events_immunity_;

    double time_step_;

    double* infectious_tau_;
    double* susceptibles_tau_;
    double* propensities_transmission_;
    double* propensities_recovery_;
    double* propensities_sampling_;
    double* propensities_mutation_;
    double* propensities_migration_;
    double* propensities_immunity_;

    Counters* counters_;
    PopulationPool* pool_;
    Infectious* infectious_data_;
    Susceptibles* susceptibles_data_;
    Chain* chain_;
    RandomGenerator* generator_;
};
