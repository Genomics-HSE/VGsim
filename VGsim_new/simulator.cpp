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

void Simulator::genealogy() {
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

boost::python::tuple Simulator::get_tree() {
    std::vector<double> time_ = arg_.get_time();
    boost::python::list times;
    for (const double & time : time_) {
        times.append(boost::python::object(time));
    }

    std::vector<int64_t> tree_ = arg_.get_tree();
    boost::python::list tree;
    for (const int64_t & node : tree_) {
        tree.append(boost::python::object(node));
    }

    return boost::python::make_tuple(tree, times);
}

PyObject* Simulator::export_chain_events() {
    boost::python::list ret;

    boost::python::list types;
    boost::python::list haplotypes;
    boost::python::list populations;
    boost::python::list newHaplotypes;
    boost::python::list newPopulations;
    boost::python::list times;

    for (uint64_t index = 0; index < chain_.GetSize(); ++index) {
        Event event = chain_.GetEvent(index);
        types.append(boost::python::object(static_cast<uint64_t>(event.type)));
        haplotypes.append(boost::python::object(event.parameter1));
        populations.append(boost::python::object(event.parameter2));
        newHaplotypes.append(boost::python::object(event.parameter3));
        newPopulations.append(boost::python::object(event.parameter4));
        times.append(boost::python::object(event.time));
    }

    ret.append(times);
    ret.append(types);
    ret.append(haplotypes);
    ret.append(populations);
    ret.append(newHaplotypes);
    ret.append(newPopulations);

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

void Simulator::set_population_size(uint64_t size, uint64_t population) {
    pool_.set_population_size(size, population);
}

PyObject* Simulator::get_population_size() {
    return pool_.get_population_size();
}

void Simulator::set_contact_density(double value, uint64_t population) {
    pool_.set_contact_density(value, population);
}

PyObject* Simulator::get_contact_density() {
    return pool_.get_contact_density();
}

void Simulator::set_npi(double after, double start, double end, uint64_t population) {
    pool_.set_npi(after, start, end, population);
}

PyObject* Simulator::get_npi() {
    return pool_.get_npi();
}

void Simulator::set_sampling_multiplier(double multiplier, uint64_t population) {
    pool_.set_sampling_multiplier(multiplier, population);
}

PyObject* Simulator::get_sampling_multiplier() {
    return pool_.get_sampling_multiplier();
}

void Simulator::set_migration_probability(double probability, uint64_t source_population, uint64_t target_population) {
    pool_.set_migration_probability(probability, source_population, target_population);
}

PyObject* Simulator::get_migration_probability() {
    return pool_.get_migration_probability();
}

uint64_t Simulator::check_migration_probability() {
    return pool_.check_migration_probability();
}


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
