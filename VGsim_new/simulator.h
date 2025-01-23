#pragma once

#include "additional_information.cpp"
#include "utils.cpp"
#include "numbers.cpp"
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
    Simulator(uint64_t number_of_sites = 0, uint64_t number_of_populations = 1, uint64_t number_of_susceptible_groups = 1, uint64_t seed = 1234);
    void Debug();

    void simulate(uint64_t iterations, uint64_t sampling, double time, std::string type, uint64_t number_attempts);
    void genealogy();
    boost::python::tuple get_current_individuals();
    PyObject* get_actual_size();
    PyObject* get_contact_density_before_lockdown();

    PyObject* get_flat_chain();
    boost::python::tuple get_tree();
    PyObject* export_chain_events();

    // Infectious
    void set_susceptibility_group(uint64_t group, uint64_t haplotype);
    PyObject* get_susceptibility_group();
    void set_transmission_rate(double rate, uint64_t haplotype);
    PyObject* get_transmission_rate();
    void set_recovery_rate(double rate, uint64_t haplotype);
    PyObject* get_recovery_rate();
    void set_sampling_rate(double rate, uint64_t haplotype);
    PyObject* get_sampling_rate();
    void set_mutation_rate(double rate, uint64_t haplotype, uint64_t mutation);
    PyObject* get_mutation_rate();
    void set_mutation_probabilities(double rate, uint64_t haplotype, uint64_t mutation, uint64_t index);
    PyObject* get_mutation_probabilities();
    void set_susceptibility(double rate, uint64_t haplotype, uint64_t group);
    PyObject* get_susceptibility();

    // Susceptibles
    void set_immunity_transition(double rate, uint64_t source_group, uint64_t target_group);
    PyObject* get_immunity_transition();

    // Population pool
    void set_population_size(uint64_t size, uint64_t population);
    PyObject* get_population_size();
    void set_contact_density(double value, uint64_t population);
    PyObject* get_contact_density();
    void set_npi(double after, double start, double end, uint64_t population);
    PyObject* get_npi();
    void set_sampling_multiplier(double multiplier, uint64_t population);
    PyObject* get_sampling_multiplier();
    void set_migration_probability(double probability, uint64_t source_population, uint64_t target_population);
    PyObject* get_migration_probability();
    uint64_t check_migration_probability();

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
