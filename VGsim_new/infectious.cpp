#pragma once

#include "utils.cpp"

#include "infectious.h"

Infectious::Infectious(uint64_t number_of_sites, uint64_t number_of_susceptible_groups)
    : number_of_sites_(number_of_sites)
    , number_of_haplotypes_(std::pow(4, number_of_sites))
    , number_of_susceptible_groups_(number_of_susceptible_groups)
    , susceptibility_groups_(new uint64_t[number_of_haplotypes_])
    , max_effective_transmission_(2.5)
    , transmission_rates_(new double[number_of_haplotypes_])
    , recovery_rates_(new double[number_of_haplotypes_])
    , sampling_rates_(new double[number_of_haplotypes_])
    , total_mutation_rates_(new double[number_of_haplotypes_])
    , mutation_rates_(new double[number_of_haplotypes_ * number_of_sites_])
    , total_sites_rates_(new double[number_of_haplotypes_ * number_of_sites_])
    , susceptibility_(new double[number_of_haplotypes_ * number_of_susceptible_groups_])
    , sites_rates_(new double[number_of_haplotypes_ * number_of_sites_ * 3]) {
    double mutation = 0.01;
    for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
        susceptibility_groups_[haplotype] = 0;
        transmission_rates_[haplotype] = 2.0;
        recovery_rates_[haplotype] = 1.0;
        sampling_rates_[haplotype] = 0.01;

        for (uint64_t site = 0; site < getNumberSites(); ++site) {
            mutation_rates_[getIndexSite(haplotype, site)] = mutation;
            for (uint64_t index = 0; index < 3; ++index) {
                sites_rates_[getIndexSite3(haplotype, site, index)] = 1.0;
            }
            total_sites_rates_[getIndexSite(haplotype, site)] = 3.0;
        }
        total_mutation_rates_[haplotype] = mutation * getNumberSites();
        susceptibility_[getIndexSus(haplotype, 0)] = 1.0;
        for (uint64_t group = 1; group < getNumberSusceptibleGroups(); ++group) {
            susceptibility_[getIndexSus(haplotype, group)] = 0.0;
        }
    }
}

Infectious::~Infectious() {
    delete[] susceptibility_groups_;
    delete[] transmission_rates_;
    delete[] recovery_rates_;
    delete[] sampling_rates_;
    delete[] total_mutation_rates_;
    delete[] mutation_rates_;
    delete[] total_sites_rates_;
    delete[] susceptibility_;
    delete[] sites_rates_;
}

void Infectious::Debug() {
    std::cout << "INFECTIOUS" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Max effective transmission: " << max_effective_transmission_ << std::endl;
    PrintArray1nd("Susceptibility types:", susceptibility_groups_, getNumberHaplotypes());
    PrintArray1nd("Transmission rates:", transmission_rates_, getNumberHaplotypes());
    PrintArray1nd("Recovery rates:", recovery_rates_, getNumberHaplotypes());
    PrintArray1nd("Sampling rates:", sampling_rates_, getNumberHaplotypes());
    PrintArray1nd("Total mutation rates:", total_mutation_rates_, getNumberHaplotypes());
    PrintArray2nd("Mutation rates", mutation_rates_, getNumberHaplotypes(), getNumberSites());
    PrintArray2nd("Total sites rates", total_sites_rates_, getNumberHaplotypes(), getNumberSites());
    PrintArray2nd("Susceptibility", susceptibility_, getNumberHaplotypes(), getNumberSusceptibleGroups());
    PrintArray3nd("Sites rates", sites_rates_, getNumberHaplotypes(), getNumberSites(), 3);
    std::cout << std::endl;
}


void Infectious::Update() {
    for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
        total_mutation_rates_[haplotype] = 0.0;
        for (uint64_t site = 0; site < getNumberSites(); ++site) {
            total_mutation_rates_[haplotype] += mutation_rates_[getIndexSite(haplotype, site)];
            total_sites_rates_[getIndexSite(haplotype, site)] = 0.0;
            for (uint64_t index = 0; index < 3; ++index) {
                total_sites_rates_[getIndexSite(haplotype, site)] += sites_rates_[getIndexSite3(haplotype, site, index)];
            }
        }
    }

    max_effective_transmission_ = 0.0;
    for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
        for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
            max_effective_transmission_ = std::max(max_effective_transmission_, transmission_rates_[haplotype] * susceptibility_[getIndexSus(haplotype, group)]);
        }
    }
}


void Infectious::set_susceptibility_group(uint64_t group, uint64_t haplotype) {
    susceptibility_groups_[haplotype] = group;
}

PyObject* Infectious::get_susceptibility_group() {
    boost::python::list data;

    for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
        data.append(boost::python::object(susceptibility_groups_[haplotype]));
    }

    return boost::python::incref(data.ptr());
}

void Infectious::set_transmission_rate(double rate, uint64_t haplotype) {
    transmission_rates_[haplotype] = rate;
}

PyObject* Infectious::get_transmission_rate() {
    boost::python::list data;

    for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
        data.append(boost::python::object(transmission_rates_[haplotype]));
    }

    return boost::python::incref(data.ptr());
}

void Infectious::set_recovery_rate(double rate, uint64_t haplotype) {
    recovery_rates_[haplotype] = rate;
}

PyObject* Infectious::get_recovery_rate() {
    boost::python::list data;

    for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
        data.append(boost::python::object(recovery_rates_[haplotype]));
    }

    return boost::python::incref(data.ptr());
}

void Infectious::set_sampling_rate(double rate, uint64_t haplotype) {
    sampling_rates_[haplotype] = rate;
}

PyObject* Infectious::get_sampling_rate() {
    boost::python::list data;

    for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
        data.append(boost::python::object(sampling_rates_[haplotype]));
    }

    return boost::python::incref(data.ptr());
}

void Infectious::set_mutation_rate(double rate, uint64_t haplotype, uint64_t mutation) {
    mutation_rates_[getIndexSite(haplotype, mutation)] = rate;
}

PyObject* Infectious::get_mutation_rate() {
    boost::python::list data;

    for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
        boost::python::list str_data;
        for (uint64_t site = 0; site < getNumberSites(); ++site) {
            str_data.append(boost::python::object(mutation_rates_[getIndexSite(haplotype, site)]));
        }
        data.append(str_data);
    }

    return boost::python::incref(data.ptr());
}

void Infectious::set_mutation_probabilities(double rate, uint64_t haplotype, uint64_t mutation, uint64_t index) {
    sites_rates_[getIndexSite3(haplotype, mutation, index)] = rate;
}

PyObject* Infectious::get_mutation_probabilities() {
    boost::python::list data;

    for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
        boost::python::list str_data;
        for (uint64_t site = 0; site < getNumberSites(); ++site) {
            boost::python::list str_data_2;
            for (uint64_t index = 0; index < 3; ++index) {
                str_data_2.append(boost::python::object(sites_rates_[getIndexSite3(haplotype, site, index)]));
            }
            str_data.append(str_data_2);
        }
        data.append(str_data);
    }

    return boost::python::incref(data.ptr());
}

void Infectious::set_susceptibility(double rate, uint64_t haplotype, uint64_t group) {
    susceptibility_[getIndexSus(haplotype, group)] = rate;
}

PyObject* Infectious::get_susceptibility() {
    boost::python::list data;

    for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
        boost::python::list str_data;
        for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
            str_data.append(boost::python::object(susceptibility_[getIndexSus(haplotype, group)]));
        }
        data.append(str_data);
    }

    return boost::python::incref(data.ptr());
}


inline double Infectious::GetMaxEffectiveTransmission() const {
    return max_effective_transmission_;
}

inline double Infectious::GetTransmissionSusceptibility(uint64_t haplotype, uint64_t group) const {
    return transmission_rates_[haplotype] * susceptibility_[getIndexSus(haplotype, group)];
}


inline uint64_t Infectious::GetSusceptibilityTypes(uint64_t haplotype) const {
    return susceptibility_groups_[haplotype];
}

inline double* Infectious::GetTransmissionRateBegin() const {
    return transmission_rates_;
}

inline double* Infectious::GetRecoveryRateBegin() const {
    return recovery_rates_;
}

inline double* Infectious::GetSamplingRateBegin() const {
    return sampling_rates_;
}

inline double* Infectious::GetMutationRateHapBegin(uint64_t haplotype) const {
    return &mutation_rates_[getIndexSite(haplotype, 0)];
}

inline double* Infectious::GetMutationRateBegin() const {
    return mutation_rates_;
}

inline double Infectious::GetTotalMutationRateHap(uint64_t haplotype) const {
    return total_mutation_rates_[haplotype];
}

inline double* Infectious::GetTotalMutationRateBegin() const {
    return total_mutation_rates_;
}

inline double Infectious::GetSitesRateHapSiteProbability(uint64_t haplotype, uint64_t site, uint64_t index) const {
    return sites_rates_[getIndexSite3(haplotype, site, index)] / GetTotalSitesRateHapSite(haplotype, site);
}

inline double* Infectious::GetSitesRateHapSiteBegin(uint64_t haplotype, uint64_t site) const {
    return &sites_rates_[getIndexSite3(haplotype, site, 0)];
}

inline double Infectious::GetTotalSitesRateHapSite(uint64_t haplotype, uint64_t site) const {
    return total_sites_rates_[getIndexSite(haplotype, site)];
}

inline double Infectious::GetSusceptibility(uint64_t haplotype, uint64_t group) const {
    return susceptibility_[getIndexSus(haplotype, group)];
}


inline uint64_t Infectious::getNumberSites() const {
    return number_of_sites_;
}

inline uint64_t Infectious::getNumberHaplotypes() const {
    return number_of_haplotypes_;
}

inline uint64_t Infectious::getNumberSusceptibleGroups() const {
    return number_of_susceptible_groups_;
}

inline uint64_t Infectious::getIndexSite(uint64_t first, uint64_t second) const {
    return first * getNumberSites() + second;
}

inline uint64_t Infectious::getIndexSus(uint64_t first, uint64_t second) const {
    return first * getNumberSusceptibleGroups() + second;
}

inline uint64_t Infectious::getIndexSite3(uint64_t first, uint64_t second, uint64_t third) const {
    return (first * getNumberSites() + second) * 3 + third;
}
