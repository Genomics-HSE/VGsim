#pragma once

#include "utils.cpp"

#include "susceptibles.h"

Susceptibles::Susceptibles(uint64_t number_of_susceptible_groups)
    : number_of_susceptible_groups_(number_of_susceptible_groups)
    , susceptibility_cumul_transition_(new double[number_of_susceptible_groups_])
    , susceptibility_transition_(new double[number_of_susceptible_groups_ * number_of_susceptible_groups_]) {
    double transition = 0.00001;
    for (uint64_t source = 0; source < getNumberSusceptibleGroups(); ++source) {
        susceptibility_cumul_transition_[source] = transition * (getNumberSusceptibleGroups() - 1);
        for (uint64_t target = 0; target < getNumberSusceptibleGroups(); ++target) {
            susceptibility_transition_[getIndexSus(source, target)] = source == target ? 0.0 : transition;
        }
    }
}

Susceptibles::~Susceptibles() {
    delete[] susceptibility_cumul_transition_;
    delete[] susceptibility_transition_;
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


void Susceptibles::SetSusceptibilityTransition(double rate, uint64_t source, uint64_t target) {
    susceptibility_transition_[getIndexSus(source, target)] = rate;
}

inline double Susceptibles::GetSusceptibilityTransition(uint64_t source, uint64_t target) const {
    return susceptibility_transition_[getIndexSus(source, target)];
}

inline double* Susceptibles::GetSusceptibilityTransitionBegin(uint64_t group) const {
    return &susceptibility_transition_[getIndexSus(group, 0)];
}

inline double Susceptibles::GetSusceptibilityCumulTransition(uint64_t group) const {
    return susceptibility_cumul_transition_[group];
}

inline double* Susceptibles::GetSusceptibilityCumulTransitionBegin() const {
    return susceptibility_cumul_transition_;
}


inline uint64_t Susceptibles::getNumberSusceptibleGroups() const {
    return number_of_susceptible_groups_;
}

inline uint64_t Susceptibles::getIndexSus(uint64_t first, uint64_t second) const {
    return first * getNumberSusceptibleGroups() + second;
}
