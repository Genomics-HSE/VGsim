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

void Simulator::set_susceptibility_group(uint64_t group, uint64_t haplotype) {
    infectious_data_.set_susceptibility_group(group, haplotype);
}

PyObject* Simulator::get_susceptibility_group() {
    return infectious_data_.get_susceptibility_group();
}

void Simulator::set_transmission_rate(double rate, uint64_t haplotype) {
    infectious_data_.set_transmission_rate(rate, haplotype);
}

PyObject* Simulator::get_transmission_rate() {
    return infectious_data_.get_transmission_rate();
}

void Simulator::set_recovery_rate(double rate, uint64_t haplotype) {
    infectious_data_.set_recovery_rate(rate, haplotype);
}

PyObject* Simulator::get_recovery_rate() {
    return infectious_data_.get_recovery_rate();
}

void Simulator::set_sampling_rate(double rate, uint64_t haplotype) {
    infectious_data_.set_sampling_rate(rate, haplotype);
}

PyObject* Simulator::get_sampling_rate() {
    return infectious_data_.get_sampling_rate();
}

void Simulator::set_mutation_rate(double rate, uint64_t haplotype, uint64_t mutation) {
    infectious_data_.set_mutation_rate(rate, haplotype, mutation);
}

PyObject* Simulator::get_mutation_rate() {
    return infectious_data_.get_mutation_rate();
}

void Simulator::set_mutation_probabilities(double rate, uint64_t haplotype, uint64_t mutation, uint64_t index) {
    infectious_data_.set_mutation_probabilities(rate, haplotype, mutation, index);
}

PyObject* Simulator::get_mutation_probabilities() {
    return infectious_data_.get_mutation_probabilities();
}

void Simulator::set_susceptibility(double rate, uint64_t haplotype, uint64_t group) {
    infectious_data_.set_susceptibility(rate, haplotype, group);
}

PyObject* Simulator::get_susceptibility() {
    return infectious_data_.get_susceptibility();
}

void Simulator::set_immunity_transition(double rate, uint64_t source_group, uint64_t target_group) {
    susceptibles_data_.set_immunity_transition(rate, source_group, target_group);
}

PyObject* Simulator::get_immunity_transition() {
    return susceptibles_data_.get_immunity_transition();
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
