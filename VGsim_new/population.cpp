#pragma once

#include "lockdown.cpp"

#include "population.h"

Population::Population() {
}

Population::~Population() {
    delete[] infected_;
    delete[] susceptibles_;
    delete[] infected_save_;
    delete[] susceptibles_save_;
}

void Population::SetParameters(uint64_t number_of_haplotypes, uint64_t number_of_susceptible_groups) {
    size_ = 0;
    number_of_haplotypes_ = number_of_haplotypes;
    number_of_susceptible_groups_ = number_of_susceptible_groups;
    infected_ = new uint64_t[number_of_haplotypes];
    susceptibles_ = new uint64_t[number_of_susceptible_groups];
    infected_save_ = new uint64_t[number_of_haplotypes];
    susceptibles_save_ = new uint64_t[number_of_susceptible_groups];
    lockdown_ = Lockdown(false, 1.0, 0.0, 1.0, 1.0);
}

inline void Population::NewInfections(uint64_t count, uint64_t haplotype, uint64_t group) {
    susceptibles_[group] -= count;
    infected_[haplotype] += count;
}

inline void Population::NewRecoveries(uint64_t count, uint64_t haplotype, uint64_t group) {
    susceptibles_[group] += count;
    infected_[haplotype] -= count;
}

inline void Population::NewMutation(uint64_t count, uint64_t old_haplotype, uint64_t new_haplotype) {
    infected_[old_haplotype] -= count;
    infected_[new_haplotype] += count;
}

inline void Population::NewImmunity(uint64_t count, uint64_t old_group, uint64_t new_group) {
    susceptibles_[old_group] -= count;
    susceptibles_[new_group] += count;
}

void Population::Save() {
    for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
        infected_save_[haplotype] = infected_[haplotype];
    }
    for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
        susceptibles_save_[group] = susceptibles_[group];
    }
}

void Population::Reset() {
    for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
        infected_[haplotype] = infected_save_[haplotype];
    }
    for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
        susceptibles_[group] = susceptibles_save_[group];
    }
}


// Population
void Population::SetSize(uint64_t size) {
    size_ = size;
    susceptibles_[0] = size_;
    for (uint64_t susceptible = 1; susceptible < getNumberSusceptibleGroups(); ++susceptible) {
        susceptibles_[susceptible] = 0;
    }
    for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
        infected_[haplotype] = 0;
    }
}

inline uint64_t Population::GetSize() const {
    return size_;
}

inline uint64_t Population::GetSusceptibles(uint64_t group) const {
    return susceptibles_[group];
}

inline uint64_t* Population::GetSusceptiblesBegin() const {
    return susceptibles_;
}

inline uint64_t Population::GetInfected(uint64_t haplotype) const {
    return infected_[haplotype];
}

inline uint64_t* Population::GetInfectedBegin() const {
    return infected_;
}


// Lockdown
inline void Population::Switch() {
    lockdown_.Switch();
}

inline double Population::GetContactDensity() const {
    return lockdown_.GetCurrentContactDensity();
}

inline void Population::SetOn(bool on) {
    lockdown_.SetOn(on);
}

inline bool Population::GetOn() const {
    return lockdown_.GetOn();
}

inline void Population::SetContactDensityBefore(double contact_density_before) {
    lockdown_.SetContactDensityBefore(contact_density_before);
}

inline double Population::GetContactDensityBefore() const {
    return lockdown_.GetContactDensityBefore();
}

inline void Population::SetContactDensityAfter(double contact_density_after) {
    lockdown_.SetContactDensityAfter(contact_density_after);
}

inline double Population::GetContactDensityAfter() const {
    return lockdown_.GetContactDensityAfter();
}

inline void Population::SetStart(double start) {
    lockdown_.SetStart(start);
}

inline double Population::GetStart() const {
    return lockdown_.GetStart();
}

inline void Population::SetEnd(double end) {
    lockdown_.SetEnd(end);
}

inline double Population::GetEnd() const {
    return lockdown_.GetEnd();
}


// Private
inline uint64_t Population::getNumberHaplotypes() const {
    return number_of_haplotypes_;
}

inline uint64_t Population::getNumberSusceptibleGroups() const {
    return number_of_susceptible_groups_;
}
