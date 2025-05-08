#pragma once

#include "../utils/utils.cpp"

#include "susceptibles.h"

Susceptibles::Susceptibles(uint64_t number_of_susceptible_groups)
    : number_of_susceptible_groups_(number_of_susceptible_groups)
    , susceptibility_cumul_transition_(ArrayBase<double>(1, number_of_susceptible_groups_))
    , susceptibility_transition_(ArrayBase<double>(number_of_susceptible_groups_, number_of_susceptible_groups_)) {
    double transition = 0.0;
    for (uint64_t source = 0; source < getNumberSusceptibleGroups(); ++source) {
        susceptibility_cumul_transition_[source] = transition * (getNumberSusceptibleGroups() - 1);
        for (uint64_t target = 0; target < getNumberSusceptibleGroups(); ++target) {
            susceptibility_transition_[getIndexSus(source, target)] = source == target ? 0.0 : transition;
        }
    }
}

void Susceptibles::Debug() {
    std::cout << "SUSCEPTIBLES" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    PrintArray1nd("Susceptibility cumul transition:", susceptibility_cumul_transition_, getNumberSusceptibleGroups());
    PrintArray2nd("Susceptibility transition", susceptibility_transition_, getNumberSusceptibleGroups(), getNumberSusceptibleGroups());
    std::cout << std::endl;
}

void Susceptibles::Update() {
    for (uint64_t source = 0; source < getNumberSusceptibleGroups(); ++source) {
        susceptibility_cumul_transition_[source] = 0.0;
        for (uint64_t target = 0; target < getNumberSusceptibleGroups(); ++target) {
            if (source == target) {
                continue;
            }
            susceptibility_cumul_transition_[source] += susceptibility_transition_[getIndexSus(source, target)];
        }
    }
}


void Susceptibles::set_immunity_transition(double rate, uint64_t source_group, uint64_t target_group) {
    susceptibility_transition_[getIndexSus(source_group, target_group)] = rate;
}

boost::python::list Susceptibles::get_immunity_transition() {
    boost::python::list data;

    for (uint64_t source = 0; source < getNumberSusceptibleGroups(); ++source) {
        boost::python::list str_data;
        for (uint64_t target = 0; target < getNumberSusceptibleGroups(); ++target) {
            str_data.append(boost::python::object(susceptibility_transition_[getIndexSus(source, target)]));
        }
        data.append(str_data);
    }

    return data;
}


inline double Susceptibles::GetSusceptibilityTransition(uint64_t source, uint64_t target) const {
    return susceptibility_transition_[getIndexSus(source, target)];
}

inline double* Susceptibles::GetSusceptibilityTransitionBegin(uint64_t group) const {
    return susceptibility_transition_.getDataArray(group);
}

inline double Susceptibles::GetSusceptibilityCumulTransition(uint64_t group) const {
    return susceptibility_cumul_transition_[group];
}


inline uint64_t Susceptibles::getNumberSusceptibleGroups() const {
    return number_of_susceptible_groups_;
}

inline uint64_t Susceptibles::getIndexSus(uint64_t first, uint64_t second) const {
    return first * getNumberSusceptibleGroups() + second;
}
