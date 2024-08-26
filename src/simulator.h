#pragma once

#include "numbers.cpp"
#include "utils.cpp"
#include "chain.cpp"
#include "population_pool.cpp"
#include "infectious.cpp"
#include "susceptibles.cpp"
#include "counters.cpp"
#include "random_generator.cpp"
#include "condition_stop.cpp"
#include "direct.cpp"
#include "tau.cpp"
#include "arg.cpp"

#include <iostream>
#include <random>
#include <cmath>
#include <ctime>
#include <string>

class Simulator {
public:
    Simulator(uint64_t number_of_sites = 0, uint64_t number_of_populations = 1, uint64_t number_of_susceptible_groups = 1, uint64_t seed = 1234, uint64_t number_attempts = 100);
    void Debug();

    void Simulate(uint64_t iterations = 100'000, std::string type = "direct", uint64_t number_attempts = 100);
    void Genealogy();
    // void SetAttempts(uint64_t attempts);
    // void SetIterations(uint64_t iterations);


    // Setters
    // Susceptibles
    // void SetSusceptibilityTransition(double rate, int64_t source = -1, int64_t target = -1);

    // Getters

private:
    inline uint64_t getNumberSites() const;
    inline uint64_t getNumberHaplotypes() const;
    inline uint64_t getNumberPopulations() const;
    inline uint64_t getNumberSusceptibleGroups() const;

    Numbers numbers_;
    uint64_t seed_;

    Counters counters_;
    PopulationPool pool_;
    Infectious infectious_data_;
    Susceptibles susceptibles_data_;
    Chain chain_;
    RandomGenerator generator_;
    ConditionStop stopper_;
    Direct direct_;
    Tau tau_;
    ARG arg_;
};
