#pragma once

#include "simulator.h"

constexpr double kTime = 1'000'000;

Simulator::Simulator(uint64_t number_of_sites, uint64_t number_of_populations, uint64_t number_of_susceptible_groups, uint64_t seed)
    : numbers_({number_of_sites, static_cast<uint64_t>(std::pow(4, number_of_sites)), number_of_populations, number_of_susceptible_groups})
    , seed_(seed)
    
    , counters_(Counters())
    , pool_(PopulationPool(getNumberPopulations(), getNumberHaplotypes(), getNumberSusceptibleGroups()))
    , infectious_data_(Infectious(getNumberSites(), getNumberSusceptibleGroups()))
    , susceptibles_data_(Susceptibles(getNumberSusceptibleGroups())) 
    , chain_(Chain(numbers_))
    , generator_(RandomGenerator(seed_))
    , stopper_(ConditionStop(&counters_, &pool_, &chain_))
    , direct_(Direct(&counters_, &pool_, &infectious_data_, &susceptibles_data_, &chain_, &generator_, &stopper_, numbers_))
    , tau_(Tau(&counters_, &pool_, &infectious_data_, &susceptibles_data_, &chain_, &generator_, &stopper_, numbers_))
    , arg_(ARG(numbers_, &chain_, &counters_, &pool_, &generator_)) {
}

void Simulator::Debug() {
    // std::cout << "Number of sites: " << getNumberSites() << std::endl;
    // std::cout << "Number of haplotypes: " << getNumberHaplotypes() << std::endl;
    // std::cout << "Number of populations: " << getNumberPopulations() << std::endl;
    // std::cout << "Number of susceptible groups: " << getNumberSusceptibleGroups() << std::endl;
    // std::cout << "Seed: " << seed_ << std::endl;
    counters_.Debug();
    // pool_.Debug();
    // infectious_data_.Debug();
    // susceptibles_data_.Debug();
    // chain_.Debug();
    chain_.LastTime();
    // direct_.Debug();
    // tau_.Debug();
    // arg_.Debug();
}

void Simulator::simulate(uint64_t iterations, uint64_t sampling, double epidemic_time, std::string type, uint64_t number_attempts) {
    stopper_.SetAttempts(number_attempts);
    stopper_.SetSampling(sampling);
    stopper_.SetIterations(iterations);
    stopper_.SetEpidemicTime(epidemic_time);
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

void Simulator::Genealogy() {
    arg_.CalculateGenealogy();
}

// void Simulator::SetAttempts(uint64_t attempts) {
//     stopper_.SetAttempts(attempts);
// }

// void Simulator::SetIterations(uint64_t iterations) {
//     stopper_.SetIterations(iterations);
// }

PyObject* Simulator::get_flat_chain() {
    boost::python::list ret;

    for (uint64_t index = 0; index < chain_.GetSize(); ++index) {
        Event event = chain_.GetEvent(index);
        ret.append(boost::python::object(static_cast<uint64_t>(event.type)));
        ret.append(boost::python::object(event.parameter1));
        ret.append(boost::python::object(event.parameter2));
        ret.append(boost::python::object(event.parameter3));
        ret.append(boost::python::object(event.parameter4));
    }

    return boost::python::incref(ret.ptr());
}

void Simulator::set_transmission_rate(double rate, uint64_t haplotype) {
    infectious_data_.set_transmission_rate(rate, haplotype);
}

PyObject* Simulator::get_transmission_rate() {
    return infectious_data_.get_transmission_rate();
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
