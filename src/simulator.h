#pragma once

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

    void SetAttempts(uint64_t attempts);
    void SetIterations(uint64_t iterations);


    // Setters
    // Susceptibles
    // void SetSusceptibilityTransition(double rate, int64_t source = -1, int64_t target = -1);

    // Getters

private:
    uint64_t number_of_sites_;
    uint64_t number_of_haplotypes_;
    uint64_t number_of_populations_;
    uint64_t number_of_susceptible_groups_;
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
};
