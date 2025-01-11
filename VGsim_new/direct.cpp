#pragma once

#include "utils.cpp"

#include "direct.h"

Direct::Direct(Counters* counters, PopulationPool* pool, Infectious* infectious_data, Susceptibles* susceptibles_data, Chain* chain, RandomGenerator* generator, ConditionStop* stopper, Numbers numbers)
    : numbers_(numbers)

    , rn_(0.0)
    , rate_(0.0)
    , rate_migration_(0.0)
    , rate_pop_(new double[getNumberPopulations()])
    , rate_infection_(new double[getNumberPopulations()])
    , rate_immunity_(new double[getNumberPopulations()])
    , rate_pop_hap_(new double[getNumberPopulations() * getNumberHaplotypes()])
    , rate_pop_total_(new double[getNumberPopulations() * getNumberHaplotypes()])
    , rate_pop_sus_(new double[getNumberPopulations() * getNumberSusceptibleGroups()])
    , rate_pop_hap_event_(new double[getNumberPopulations() * getNumberHaplotypes() * 4])
    , suscept_hap_pop_rate_(new double[getNumberPopulations() * getNumberHaplotypes() * getNumberSusceptibleGroups()])
    , rate_migration_pop_(new double[getNumberPopulations()])
    , max_effective_transmission_migration_(new double[getNumberPopulations()])
    
    , counters_(counters)
    , pool_(pool)
    , infectious_data_(infectious_data)
    , susceptibles_data_(susceptibles_data) 
    , chain_(chain) 
    , generator_(generator)
    , stopper_(stopper) {
}

Direct::~Direct() {
    delete[] rate_pop_;
    delete[] rate_infection_;
    delete[] rate_immunity_;
    delete[] rate_pop_hap_;
    delete[] rate_pop_total_;
    delete[] rate_pop_sus_;
    delete[] rate_pop_hap_event_;
    delete[] suscept_hap_pop_rate_;
    delete[] rate_migration_pop_;
    delete[] max_effective_transmission_migration_;
}

void Direct::Simulate() {
    pool_->SaveInfections();
    chain_->Reserve(stopper_->GetIterations());
    while (stopper_->CheckAttempt()) {
        pool_->FirstInfections();
        Update();
        stopper_->Restart();
        while (stopper_->CheckIteration()) {
            TimeStep();
            GenerateEvent();
        }
        if (stopper_->CheckRestart()) {
            Restart();
        } else {
            break;
        }
    }
}

void Direct::Debug() {
    std::cout << "DIRECT" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Rn: " << rn_ << std::endl;
    std::cout << "Rate: " << rate_ << std::endl;
    std::cout << "Rate migration: " << rate_migration_ << std::endl;
    PrintArray1nd("Rate pop:", rate_pop_, getNumberPopulations());
    PrintArray1nd("Rate infection:", rate_infection_, getNumberPopulations());
    PrintArray1nd("Rate immunity:", rate_immunity_, getNumberPopulations());
    PrintArray1nd("Rate migration pop:", rate_migration_pop_, getNumberPopulations());
    PrintArray1nd("Max effective transmission migration:", max_effective_transmission_migration_, getNumberPopulations());
    PrintArray2nd("Rate pop hap:", rate_pop_hap_, getNumberPopulations(), getNumberHaplotypes());
    PrintArray2nd("Rate pop total:", rate_pop_total_, getNumberPopulations(), getNumberHaplotypes());
    PrintArray2nd("Rate pop sus:", rate_pop_sus_, getNumberPopulations(), getNumberSusceptibleGroups());
    PrintArray3nd("Rate pop hap event:", rate_pop_hap_event_, getNumberPopulations(), getNumberHaplotypes(), 4);
    PrintArray3nd("Suscept hap pop rate:", suscept_hap_pop_rate_, getNumberPopulations(), getNumberHaplotypes(), getNumberSusceptibleGroups());
}


void Direct::Update() {
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        pool_->CheckLockdown(population);
    }

    susceptibles_data_->Update();
    infectious_data_->Update();
    pool_->Update();

    // #TODO структура для указателей
    double* transmission_begin = infectious_data_->GetTransmissionRateBegin();
    double* recovery_begin = infectious_data_->GetRecoveryRateBegin();
    double* sampling_begin = infectious_data_->GetSamplingRateBegin();
    double* total_mutation_begin = infectious_data_->GetTotalMutationRateBegin();
    
    rate_ = 0.0;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        rate_infection_[population] = 0.0;
        uint64_t* infected_begin = pool_->GetInfectedPopHapBegin(population);
        for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
            rate_pop_hap_event_[getIndexHap4(population, haplotype, 0)] = transmission_begin[haplotype] * BirthRate(population, haplotype);
            rate_pop_hap_event_[getIndexHap4(population, haplotype, 1)] = recovery_begin[haplotype];
            rate_pop_hap_event_[getIndexHap4(population, haplotype, 2)] = sampling_begin[haplotype] * pool_->GetSamplingMultiplier(population);
            rate_pop_hap_event_[getIndexHap4(population, haplotype, 3)] = total_mutation_begin[haplotype];

            rate_pop_total_[getIndexHap(population, haplotype)] = 0.0;
            for (int64_t i = 0; i < 4; ++i) {
                rate_pop_total_[getIndexHap(population, haplotype)] += rate_pop_hap_event_[getIndexHap4(population, haplotype, i)];
            }
            rate_pop_hap_[getIndexHap(population, haplotype)] = rate_pop_total_[getIndexHap(population, haplotype)] * infected_begin[haplotype];
            rate_infection_[population] += rate_pop_hap_[getIndexHap(population, haplotype)];
            for (uint64_t type = 0; type < getNumberSusceptibleGroups(); ++type) {
                suscept_hap_pop_rate_[getIndexHapSus(population, haplotype, type)] = pool_->GetSusceptiblesPopSus(population, type) * infectious_data_->GetSusceptibility(haplotype, type);
            }
        }
        rate_immunity_[population] = 0.0;
        for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
            rate_pop_sus_[getIndexSus(population, group)] = susceptibles_data_->GetSusceptibilityCumulTransition(group) * pool_->GetSusceptiblesPopSus(population, group);
            rate_immunity_[population] += rate_pop_sus_[getIndexSus(population, group)];
        }
        rate_pop_[population] = rate_infection_[population] + rate_immunity_[population];
        rate_ += rate_pop_[population];
    }

    uint64_t infected = pool_->GetInfected();
    uint64_t* infected_begin = pool_->GetInfectedPopBegin();
    uint64_t* susceptibles_begin = pool_->GetSusceptiblesPopBegin();
    double max_effective_transmission = infectious_data_->GetMaxEffectiveTransmission();
    double* max_effective_migration = pool_->GetMaxEffectiveMigration();

    rate_migration_ = 0.0;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        max_effective_transmission_migration_[population] = max_effective_migration[population] * max_effective_transmission;
        rate_migration_pop_[population] = max_effective_transmission_migration_[population] * susceptibles_begin[population] * (infected - infected_begin[population]);
        rate_migration_ += rate_migration_pop_[population];
    }
}

void Direct::Restart() {
    generator_->NextSeed();
    pool_->Restart();
    counters_->Restart();
    chain_->Restart();
}

double Direct::BirthRate(uint64_t source, uint64_t haplotype) {
    double rate = 0.0;

    for (uint64_t type = 0; type < getNumberSusceptibleGroups(); ++type) {
        suscept_hap_pop_rate_[getIndexHapSus(source, haplotype, type)] = pool_->GetSusceptiblesPopSus(source, type) * infectious_data_->GetSusceptibility(haplotype, type);
        for (uint64_t target = 0; target < getNumberPopulations(); ++target) {
            rate += suscept_hap_pop_rate_[getIndexHapSus(source, haplotype, type)] * pool_->GetMultiplierMigrationProbability(source, target);
        }
    }

    return rate;
}

void Direct::TimeStep() {
    chain_->AddTime(-(std::log(generator_->GetUniform()) / rate_));
}

void Direct::GenerateEvent() {
    rn_ = generator_->GetUniform();
    double choose = rn_ * (rate_ + rate_migration_);
    if (rate_ > choose) {
        rn_ = choose / rate_;
        uint64_t population = fastChoose(rate_pop_, rate_, &rn_);
        choose = rn_ * rate_pop_[population];
        if (rate_immunity_[population] > choose) {
            rn_ = choose / rate_immunity_[population];
            ImmunityTransition(population);
        } else {
            rn_ = (choose - rate_immunity_[population]) / rate_infection_[population];
            uint64_t haplotype = fastChoose(&rate_pop_hap_[getIndexHap(population, 0)], rate_infection_[population], &rn_);
            uint64_t event = fastChoose(&rate_pop_hap_event_[getIndexHap4(population, haplotype, 0)], rate_pop_total_[getIndexHap(population, haplotype)], &rn_);
            switch (static_cast<TypeEvents>(event)) {
                case TypeEvents::kTRANSMISSION: Transmission(population, haplotype); break;
                case TypeEvents::kRECOVERY: Recovery(population, haplotype); break;
                case TypeEvents::kSAMPLING: Sampling(population, haplotype); break;
                case TypeEvents::kMUTATION: Mutation(population, haplotype); break;
                case TypeEvents::kMIGRATION: break;
                case TypeEvents::kSUSCCHANGE: break;
                case TypeEvents::kMULTITYPE: break;
            }
        }
    } else {
        rn_ = (choose - rate_) / rate_migration_;
        Migration();
    }
}

void Direct::UpdateRates(uint64_t population, bool infection, bool immunity, bool migration) {
    if (infection || migration) {
        pool_->CheckLockdown(population);
    }

    if (infection) {
        uint64_t* infected_begin = pool_->GetInfectedPopHapBegin(population);
        double* transmission_begin = infectious_data_->GetTransmissionRateBegin();

        rate_infection_[population] = 0.0;
        for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
            rate_pop_hap_event_[getIndexHap4(population, haplotype, 0)] = BirthRate(population, haplotype) * transmission_begin[haplotype];
            rate_pop_total_[getIndexHap(population, haplotype)] = (rate_pop_hap_event_[getIndexHap4(population, haplotype, 0)] +
                                                                   rate_pop_hap_event_[getIndexHap4(population, haplotype, 1)] +
                                                                   rate_pop_hap_event_[getIndexHap4(population, haplotype, 2)] +
                                                                   rate_pop_hap_event_[getIndexHap4(population, haplotype, 3)]);
            rate_pop_hap_[getIndexHap(population, haplotype)] = rate_pop_total_[getIndexHap(population, haplotype)] * infected_begin[haplotype];
            rate_infection_[population] += rate_pop_hap_[getIndexHap(population, haplotype)];
        }
    }

    if (immunity) {
        uint64_t* susceptibles_begin = pool_->GetSusceptiblesPopSusBegin(population);

        rate_immunity_[population] = 0.0;
        for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
            rate_pop_sus_[getIndexSus(population, group)] = susceptibles_data_->GetSusceptibilityCumulTransition(group) * susceptibles_begin[group];
            rate_immunity_[population] += rate_pop_sus_[getIndexSus(population, group)];
        }
    }

    if (infection || immunity) {
        rate_pop_[population] = rate_infection_[population] + rate_immunity_[population];
        rate_ = 0.0;
        for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
            rate_ += rate_pop_[population];
        }
    }

    if (migration) {
        uint64_t infected = pool_->GetInfected();
        uint64_t* infected_begin = pool_->GetInfectedPopBegin();
        uint64_t* susceptibles_begin = pool_->GetSusceptiblesPopBegin();

        rate_migration_ = 0.0;
        for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
            rate_migration_pop_[population] = max_effective_transmission_migration_[population] * susceptibles_begin[population] * (infected - infected_begin[population]);
            rate_migration_ += rate_migration_pop_[population];
        }
    }
} 

void Direct::ImmunityTransition(uint64_t population) {
    uint64_t source_type = fastChoose(&rate_pop_sus_[getIndexSus(population, 0)], rate_immunity_[population], &rn_);
    uint64_t target_type = fastChoose(susceptibles_data_->GetSusceptibilityTransitionBegin(source_type), susceptibles_data_->GetSusceptibilityCumulTransition(source_type), &rn_);
    pool_->NewImmunity(1, source_type, target_type, population);

    UpdateRates(population, false, true, false);
    counters_->AddImmunity(1);
    chain_->AddEvent({TypeEvents::kSUSCCHANGE, source_type, population, target_type, 0});
}

void Direct::Transmission(uint64_t population, uint64_t haplotype) {
    double sum = 0.0;
    for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
        sum += suscept_hap_pop_rate_[getIndexHapSus(population, haplotype, group)];
    }
    uint64_t group = fastChoose(&suscept_hap_pop_rate_[getIndexHapSus(population, haplotype, 0)], sum, &rn_);
    pool_->NewInfections(1, haplotype, group, population);

    UpdateRates(population, true, true, true);
    counters_->AddTransmission(1);
    chain_->AddEvent({TypeEvents::kTRANSMISSION, haplotype, population, group, getNumberHaplotypes()});
}

void Direct::Recovery(uint64_t population, uint64_t haplotype) {
    uint64_t group = infectious_data_->GetSusceptibilityTypes(haplotype);
    pool_->NewRecoveries(1, haplotype, group, population);

    UpdateRates(population, true, true, true);
    counters_->AddRecovery(1);
    chain_->AddEvent({TypeEvents::kRECOVERY, haplotype, population, infectious_data_->GetSusceptibilityTypes(haplotype), 0});
}

void Direct::Sampling(uint64_t population, uint64_t haplotype) {
    uint64_t group = infectious_data_->GetSusceptibilityTypes(haplotype);
    pool_->NewRecoveries(1, haplotype, group, population);

    UpdateRates(population, true, true, true);
    counters_->AddSampling(1);
    chain_->AddEvent({TypeEvents::kSAMPLING, haplotype, population, infectious_data_->GetSusceptibilityTypes(haplotype), 0});
}

void Direct::Mutation(uint64_t population, uint64_t haplotype) {
    uint64_t site = fastChoose(infectious_data_->GetMutationRateHapBegin(haplotype), infectious_data_->GetTotalMutationRateHap(haplotype), &rn_);
    uint64_t allele = fastChoose(infectious_data_->GetSitesRateHapSiteBegin(haplotype, site), infectious_data_->GetTotalSitesRateHapSite(haplotype, site), &rn_);
    uint64_t new_haplotype = GetNewHaplotype(haplotype, site, allele, getNumberSites());
    pool_->NewMutation(1, haplotype, new_haplotype, population);

    UpdateRates(population, true, false, false);
    counters_->AddMutation(1);
    chain_->AddEvent({TypeEvents::kMUTATION, haplotype, population, new_haplotype, 0});
}

void Direct::Migration() {
    uint64_t target_population = fastChoose(rate_migration_pop_, rate_migration_, &rn_);
    uint64_t source_population = fastChooseSkip(pool_->GetInfectedPopBegin(), pool_->GetInfected() - pool_->GetInfectedPop(target_population), &rn_, target_population);
    uint64_t haplotype = fastChoose(pool_->GetInfectedPopHapBegin(source_population), pool_->GetInfectedPop(source_population), &rn_);
    uint64_t group = fastChoose(pool_->GetSusceptiblesPopSusBegin(target_population), pool_->GetSusceptiblesPop(target_population), &rn_);

    double accept = pool_->GetEffectiveMigration(source_population, target_population) * infectious_data_->GetTransmissionSusceptibility(haplotype, group) / max_effective_transmission_migration_[target_population];
    if (rn_ < accept) {
        pool_->NewInfections(1, haplotype, group, target_population);
        UpdateRates(target_population, true, true, true);
        counters_->AddMigrationAccept(1);
        chain_->AddEvent({TypeEvents::kMIGRATION, haplotype, source_population, group, target_population});
    } else {
        counters_->AddMigrationReject(1);
    }
}


inline uint64_t Direct::getNumberSites() const {
    return numbers_.sites;
}

inline uint64_t Direct::getNumberHaplotypes() const {
    return numbers_.haplotypes;
}

inline uint64_t Direct::getNumberPopulations() const {
    return numbers_.populations;
}

inline uint64_t Direct::getNumberSusceptibleGroups() const {
    return numbers_.susceptible_groups;
}

inline uint64_t Direct::getIndexHap(uint64_t first, uint64_t second) const {
    return first * getNumberHaplotypes() + second;
}

inline uint64_t Direct::getIndexSus(uint64_t first, uint64_t second) const {
    return first * getNumberSusceptibleGroups() + second;
}

inline uint64_t Direct::getIndexHap4(uint64_t first, uint64_t second, uint64_t third) const {
    return (first * getNumberHaplotypes() + second) * 4 + third;
}

inline uint64_t Direct::getIndexHapSus(uint64_t first, uint64_t second, uint64_t third) const {
    return (first * getNumberHaplotypes() + second) * getNumberSusceptibleGroups() + third;
}
