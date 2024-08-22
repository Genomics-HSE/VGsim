#pragma once

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
    : number_of_sites_(number_of_sites)
    , number_of_haplotypes_(std::pow(4, number_of_sites))
    , number_of_populations_(number_of_populations)
    , number_of_susceptible_groups_(number_of_susceptible_groups)
    , seed_(seed)
    
    , counters_(Counters())
    , pool_(PopulationPool(number_of_populations_, number_of_haplotypes_, number_of_susceptible_groups_))
    , infectious_data_(Infectious(number_of_sites_, number_of_susceptible_groups_))
    , susceptibles_data_(Susceptibles(number_of_susceptible_groups_)) 
    , chain_(Chain())
    , generator_(RandomGenerator(seed_))
    , stopper_(ConditionStop())
    , direct_(Direct(&counters_, &pool_, &infectious_data_, &susceptibles_data_, &chain_, &generator_, &stopper_, number_of_sites_, number_of_haplotypes_, number_of_populations_, number_of_susceptible_groups_))
    , tau_(Tau(&counters_, &pool_, &infectious_data_, &susceptibles_data_, &chain_, &generator_, number_of_sites_, number_of_haplotypes_, number_of_populations_, number_of_susceptible_groups_)) {
}

void Simulator::Debug() {
    std::cout << "Number of sites: " << number_of_sites_ << std::endl;
    std::cout << "Number of haplotypes: " << number_of_haplotypes_ << std::endl;
    std::cout << "Number of populations: " << number_of_populations_ << std::endl;
    std::cout << "Number of susceptible groups: " << number_of_susceptible_groups_ << std::endl;
    std::cout << "Seed: " << seed_ << std::endl;
    counters_.Debug();
    pool_.Debug();
    // infectious_data_.Debug();
    // susceptibles_data_.Debug();
    // direct_.Debug();
    tau_.Debug();
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
//         CheckIndex(source, number_of_susceptible_groups_, "source susceptibility type");
//         CheckIndex(target, number_of_susceptible_groups_, "target susceptibility type");
//     } catch (const std::exception& e) {
//         std::cout << e.what() << std::endl;
//         return;
//     }

//     for (uint64_t si : GetIndexes(source, number_of_susceptible_groups_)) {
//         for (uint64_t ti : GetIndexes(target, number_of_susceptible_groups_)) {
//             if (si != ti) {
//                 susceptibles_data_.SetSusceptibilityTransition(rate, si, ti);
//             }
//         }
//     }
// }
