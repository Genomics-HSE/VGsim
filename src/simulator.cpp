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

#include "simulator.h"

constexpr double kTime = 1'000'000;

Simulator::Simulator(uint64_t number_of_sites, uint64_t number_of_populations, uint64_t number_of_susceptible_groups, uint64_t seed, uint64_t number_attempts)
    : numbers_({number_of_sites, static_cast<uint64_t>(std::pow(4, number_of_sites)), number_of_populations, number_of_susceptible_groups})
    , seed_(seed)
    
    , counters_(Counters())
    , pool_(PopulationPool(getNumberPopulations(), getNumberHaplotypes(), getNumberSusceptibleGroups()))
    , infectious_data_(Infectious(getNumberSites(), getNumberSusceptibleGroups()))
    , susceptibles_data_(Susceptibles(getNumberSusceptibleGroups())) 
    , chain_(Chain())
    , generator_(RandomGenerator(seed_))
    , stopper_(ConditionStop())
    , direct_(Direct(&counters_, &pool_, &infectious_data_, &susceptibles_data_, &chain_, &generator_, &stopper_, numbers_))
    , tau_(Tau(&counters_, &pool_, &infectious_data_, &susceptibles_data_, &chain_, &generator_, numbers_)) {
}

void Simulator::Debug() {
    std::cout << "Number of sites: " << getNumberSites() << std::endl;
    std::cout << "Number of haplotypes: " << getNumberHaplotypes() << std::endl;
    std::cout << "Number of populations: " << getNumberPopulations() << std::endl;
    std::cout << "Number of susceptible groups: " << getNumberSusceptibleGroups() << std::endl;
    std::cout << "Seed: " << seed_ << std::endl;
    counters_.Debug();
    pool_.Debug();
    infectious_data_.Debug();
    susceptibles_data_.Debug();
    direct_.Debug();
    tau_.Debug();
    chain_.Debug();
}

void Simulator::Simulate(uint64_t iterations, std::string type, uint64_t number_attempts) {
    uint64_t start_time = clock();
    if (type == "direct") {
        direct_.Simulate();
    } else if (type == "tau") {
        std::cout << "Tau!" << std::endl;
        tau_.Simulate(iterations, number_attempts);
    } else {
        std::cout << "Unknown type!" << std::endl;
    }
    uint64_t end_time = clock();
    std::cout << "Time: " << (end_time - start_time) / kTime << " s" << std::endl;
}

void Simulator::SetAttempts(uint64_t attempts) {
    stopper_.SetAttempts(attempts);
}

void Simulator::SetIterations(uint64_t iterations) {
    stopper_.SetIterations(iterations);
}

// void Simulator::SetSusceptibilityTransition(double rate, int64_t source, int64_t target) {
//     try {
//         CheckValue(rate, "immunity transition rate");
//         CheckIndex(source, getNumberSusceptibleGroups(), "source susceptibility type");
//         CheckIndex(target, getNumberSusceptibleGroups(), "target susceptibility type");
//     } catch (const std::exception& e) {
//         std::cout << e.what() << std::endl;
//         return;
//     }

//     for (uint64_t si : GetIndexes(source, getNumberSusceptibleGroups())) {
//         for (uint64_t ti : GetIndexes(target, getNumberSusceptibleGroups())) {
//             if (si != ti) {
//                 susceptibles_data_.SetSusceptibilityTransition(rate, si, ti);
//             }
//         }
//     }
// }

inline uint64_t Simulator::getNumberSites() const {
    return numbers_.sites;
}

inline uint64_t Simulator::getNumberHaplotypes() const {
    return numbers_.haplotypes;
}

inline uint64_t Simulator::getNumberPopulations() const {
    return numbers_.populations;
}

inline uint64_t Simulator::getNumberSusceptibleGroups() const {
    return numbers_.susceptible_groups;
}
