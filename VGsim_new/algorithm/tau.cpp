#pragma once

#include "../utils/utils.cpp"

#include "tau.h"

Tau::Tau(Counters* counters, PopulationPool* pool, Infectious* infectious_data, Susceptibles* susceptibles_data, Chain* chain, RandomGenerator* generator, ConditionStop* stopper, Numbers numbers)
    : numbers_(numbers)
    
    , infectious_delta_(new int64_t[getNumberPopulations() * getNumberHaplotypes()])
    , susceptibles_delta_(new int64_t[getNumberPopulations() * getNumberSusceptibleGroups()])
    , events_transmission_(new uint64_t[getNumberPopulations() * getNumberHaplotypes() * getNumberSusceptibleGroups()])
    , events_recovery_(new uint64_t[getNumberPopulations() * getNumberHaplotypes()])
    , events_sampling_(new uint64_t[getNumberPopulations() * getNumberHaplotypes()])
    , events_mutation_(new uint64_t[getNumberPopulations() * getNumberHaplotypes() * getNumberSites() * (getNumberAlleleStates() - 1)])
    , events_migration_(new uint64_t[getNumberPopulations() * getNumberPopulations() * getNumberHaplotypes() * getNumberSusceptibleGroups()])
    , events_immunity_(new uint64_t[getNumberPopulations() * getNumberSusceptibleGroups() * getNumberSusceptibleGroups()])

    , time_step_(0.0)

    , infectious_tau_(new double[getNumberPopulations() * getNumberHaplotypes()])
    , susceptibles_tau_(new double[getNumberPopulations() * getNumberSusceptibleGroups()])
    , propensities_transmission_(new double[getNumberPopulations() * getNumberHaplotypes() * getNumberSusceptibleGroups()])
    , propensities_recovery_(new double[getNumberPopulations() * getNumberHaplotypes()])
    , propensities_sampling_(new double[getNumberPopulations() * getNumberHaplotypes()])
    , propensities_mutation_(new double[getNumberPopulations() * getNumberHaplotypes() * getNumberSites() * (getNumberAlleleStates() - 1)])
    , propensities_migration_(new double[getNumberPopulations() * getNumberPopulations() * getNumberHaplotypes() * getNumberSusceptibleGroups()])
    , propensities_immunity_(new double[getNumberPopulations() * getNumberSusceptibleGroups() * getNumberSusceptibleGroups()])

    , counters_(counters)
    , pool_(pool)
    , infectious_data_(infectious_data)
    , susceptibles_data_(susceptibles_data) 
    , chain_(chain)
    , generator_(generator)
    , stopper_(stopper) {
}

Tau::~Tau() {
    delete[] infectious_delta_;
    delete[] susceptibles_delta_;
    delete[] events_transmission_;
    delete[] events_recovery_;
    delete[] events_sampling_;
    delete[] events_mutation_;
    delete[] events_migration_;
    delete[] events_immunity_;

    delete[] infectious_tau_;
    delete[] susceptibles_tau_;
    delete[] propensities_transmission_;
    delete[] propensities_recovery_;
    delete[] propensities_sampling_;
    delete[] propensities_mutation_;
    delete[] propensities_migration_;
    delete[] propensities_immunity_;
}

void Tau::Debug() {
    std::cout << "TAU" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    PrintArray2nd("Infectious delta", infectious_delta_, getNumberPopulations(), getNumberHaplotypes());
    PrintArray2nd("Susceptibles delta", susceptibles_delta_, getNumberPopulations(), getNumberSusceptibleGroups());
    PrintArray3nd("Events Transmission", events_transmission_, getNumberPopulations(), getNumberHaplotypes(), getNumberSusceptibleGroups());
    PrintArray2nd("Events Recovery", events_recovery_, getNumberPopulations(), getNumberHaplotypes());
    PrintArray2nd("Events Sampling", events_sampling_, getNumberPopulations(), getNumberHaplotypes());
    PrintArray4nd("Events Mutation", events_mutation_, getNumberPopulations(), getNumberHaplotypes(), getNumberSites(), getNumberAlleleStates() - 1);
    PrintArray4nd("Events Migration", events_migration_, getNumberPopulations(), getNumberPopulations(), getNumberHaplotypes(), getNumberSusceptibleGroups());
    PrintArray3nd("Events Immunity", events_immunity_, getNumberPopulations(), getNumberSusceptibleGroups(), getNumberSusceptibleGroups());

    PrintArray2nd("Infectious tau", infectious_tau_, getNumberPopulations(), getNumberHaplotypes());
    PrintArray2nd("Susceptibles tau", susceptibles_tau_, getNumberPopulations(), getNumberSusceptibleGroups());
    PrintArray3nd("Propensities Transmission", propensities_transmission_, getNumberPopulations(), getNumberHaplotypes(), getNumberSusceptibleGroups());
    PrintArray2nd("Propensities Recovery", propensities_recovery_, getNumberPopulations(), getNumberHaplotypes());
    PrintArray2nd("Propensities Sampling", propensities_sampling_, getNumberPopulations(), getNumberHaplotypes());
    PrintArray4nd("Propensities Mutation", propensities_mutation_, getNumberPopulations(), getNumberHaplotypes(), getNumberSites(), getNumberAlleleStates() - 1);
    PrintArray4nd("Propensities Migration", propensities_migration_, getNumberPopulations(), getNumberPopulations(), getNumberHaplotypes(), getNumberSusceptibleGroups());
    PrintArray3nd("Propensities Immunity", propensities_immunity_, getNumberPopulations(), getNumberSusceptibleGroups(), getNumberSusceptibleGroups());
}

void Tau::Simulate(uint64_t iterations, uint64_t number_attempts) {
    pool_->SaveInfections();
    chain_->Reserve(iterations);
    while (stopper_->CheckAttempt()) {
        pool_->FirstInfections();
        Update();
        stopper_->Restart();
        while (stopper_->CheckIteration()) {
            Propensities();
            TimeStep();
            UpdateCompartmentCounts();
        }
        // if (iterations >= 5 && iter < 5 && i != number_attempts) {
        //     // Restart(i);
        //     continue;
        // }
        break;
    }
}

void Tau::Update() {
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        pool_->CheckLockdown(population);
    }

    susceptibles_data_->Update();
    infectious_data_->Update();
    pool_->Update();
}

void Tau::Propensities() {
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
            infectious_tau_[getIndexHap(population, haplotype)] = 0.0;
        }
        for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
            susceptibles_tau_[getIndexSus(population, group)] = 0.0;
        }
    }

    uint64_t index_transmission = 0;
    double* transmission_rate = infectious_data_->GetTransmissionRateBegin();
    for (uint64_t target = 0; target < getNumberPopulations(); ++target) {
        uint64_t* infected = pool_->GetInfectedPopHapBegin(target);
        uint64_t* susceptibles = pool_->GetSusceptiblesPopSusBegin(target);
        for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
            for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
                propensities_transmission_[index_transmission] = 0.0;
                for (uint64_t source = 0; source < getNumberPopulations(); ++source) {
                    propensities_transmission_[index_transmission] += transmission_rate[haplotype] * infected[haplotype] *\
                    susceptibles[group] * infectious_data_->GetSusceptibility(haplotype, group) * pool_->GetMigrationProbability(source, target) *\
                    pool_->GetMigrationProbability(source, target) * pool_->GetContactDensity(source) / pool_->GetActualSize(source);
                }

                infectious_tau_[getIndexHap(target, haplotype)] += propensities_transmission_[index_transmission];
                susceptibles_tau_[getIndexSus(target, group)] -= propensities_transmission_[index_transmission];
                ++index_transmission;
            }
        }
    }

    uint64_t index_recovery_sampling = 0;
    double* recovery_rate = infectious_data_->GetRecoveryRateBegin();
    double* sampling_rate = infectious_data_->GetSamplingRateBegin();
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        uint64_t* infected = pool_->GetInfectedPopHapBegin(population);
        for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
            propensities_recovery_[index_recovery_sampling] = recovery_rate[haplotype] * infected[haplotype];
            propensities_sampling_[index_recovery_sampling] = sampling_rate[haplotype] * infected[haplotype];

            infectious_tau_[index_recovery_sampling] -= propensities_recovery_[index_recovery_sampling];
            susceptibles_tau_[getIndexSus(population, infectious_data_->GetSusceptibilityTypes(haplotype))] += propensities_recovery_[index_recovery_sampling];
            infectious_tau_[index_recovery_sampling] -= propensities_sampling_[index_recovery_sampling];
            susceptibles_tau_[getIndexSus(population, infectious_data_->GetSusceptibilityTypes(haplotype))] += propensities_sampling_[index_recovery_sampling];
            ++index_recovery_sampling;
        }
    }

    uint64_t index_mutation = 0;
    double* mutation_rate = infectious_data_->GetMutationRateBegin();
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        uint64_t* infected = pool_->GetInfectedPopHapBegin(population);
        for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
            for (uint64_t site = 0; site < getNumberSites(); ++site) {
                for (uint64_t i = 0; i < getNumberAlleleStates() - 1; ++i) {
                    propensities_mutation_[index_mutation] = mutation_rate[getIndexSit(haplotype, site)] * infected[haplotype] * infectious_data_->GetSitesRateHapSiteProbability(haplotype, site, i);

                    infectious_tau_[getIndexHap(population, GetNewHaplotype(haplotype, site, i, getNumberSites()))] += propensities_mutation_[index_mutation];
                    infectious_tau_[getIndexHap(population, haplotype)] -= propensities_mutation_[index_mutation];
                    ++index_mutation;
                }
            }
        }
    }

    uint64_t index_migration = 0;
    for (uint64_t source = 0; source < getNumberPopulations(); ++source) {
        uint64_t* infected = pool_->GetInfectedPopHapBegin(source);
        for (uint64_t target = 0; target < getNumberPopulations(); ++target) {
            if (source == target) {
                continue;
            }
            uint64_t* susceptibles = pool_->GetSusceptiblesPopSusBegin(target);
            for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
                for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
                    propensities_migration_[index_migration] = pool_->GetEffectiveMigration(source, target) * infected[haplotype] * susceptibles[group] * infectious_data_->GetSusceptibility(haplotype, group) * pool_->GetMigrationProbability(source, target);

                    infectious_tau_[getIndexHap(target, haplotype)] += propensities_migration_[index_migration];
                    susceptibles_tau_[getIndexSus(target, group)] -= propensities_migration_[index_migration];
                    ++index_migration;
                }
            }
        }
    }

    uint64_t index_immunity = 0;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        for (uint64_t source = 0; source < getNumberSusceptibleGroups(); ++source) {
            uint64_t* susceptibles = pool_->GetSusceptiblesPopSusBegin(population);
            for (uint64_t target = 0; target < getNumberSusceptibleGroups(); ++target) {
                if (source == target) {
                    continue;
                }
                propensities_immunity_[index_immunity] = susceptibles_data_->GetSusceptibilityTransition(source, target) * susceptibles[source];

                susceptibles_tau_[getIndexSus(population, target)] += propensities_immunity_[index_immunity];
                susceptibles_tau_[getIndexSus(population, source)] -= propensities_immunity_[index_immunity];
                ++index_immunity;
            }
        }
    }
}

void Tau::TimeStep() {
    ChooseTime();
    while (!GenerateEvents()) {
        time_step_ /= 2;
    }
    chain_->AddTime(time_step_);
}

void Tau::ChooseTime() {
    double epsilon = 0.03;
    double tmp;

    time_step_ = 1.0;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        uint64_t* infected = pool_->GetInfectedPopHapBegin(population);
        for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
            tmp = std::abs(infectious_tau_[getIndexHap(population, haplotype)]);
            if (tmp < 1e-8) {
                continue;
            }
            time_step_ = std::min(time_step_, std::max(epsilon * infected[haplotype], 1.0) / tmp);
        }

        uint64_t* susceptible = pool_->GetSusceptiblesPopSusBegin(population);
        for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
            tmp = std::abs(susceptibles_tau_[getIndexSus(population, group)]);
            if (tmp < 1e-8) {
                continue;
            }
            time_step_ = std::min(time_step_, std::max(epsilon * susceptible[group], 1.0) / tmp);
        }
    }
}

bool Tau::GenerateEvents() {
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
            infectious_delta_[getIndexHap(population, haplotype)] = 0;
        }
        for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
            susceptibles_delta_[getIndexSus(population, group)] = 0;
        }
    }

    uint64_t index_transmission = 0;
    for (uint64_t target = 0; target < getNumberPopulations(); ++target) {
        for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
            for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
                events_transmission_[index_transmission] = generator_->GetPoisson(time_step_ * propensities_transmission_[index_transmission]);

                infectious_delta_[getIndexHap(target, haplotype)] += events_transmission_[index_transmission];
                susceptibles_delta_[getIndexSus(target, group)] -= events_transmission_[index_transmission];
                ++index_transmission;
            }
        }
    }

    uint64_t index_recovery_sampling = 0;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
            events_recovery_[index_recovery_sampling] = generator_->GetPoisson(time_step_ * propensities_recovery_[index_recovery_sampling]);
            events_sampling_[index_recovery_sampling] = generator_->GetPoisson(time_step_ * propensities_sampling_[index_recovery_sampling]);
            
            infectious_delta_[index_recovery_sampling] -= events_recovery_[index_recovery_sampling];
            susceptibles_delta_[getIndexSus(population, infectious_data_->GetSusceptibilityTypes(haplotype))] += events_recovery_[index_recovery_sampling];
            infectious_delta_[index_recovery_sampling] -= events_sampling_[index_recovery_sampling];
            susceptibles_delta_[getIndexSus(population, infectious_data_->GetSusceptibilityTypes(haplotype))] += events_sampling_[index_recovery_sampling];
            ++index_recovery_sampling;
        }
    }

    uint64_t index_mutation = 0;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
            for (uint64_t site = 0; site < getNumberSites(); ++site) {
                for (uint64_t i = 0; i < getNumberAlleleStates() - 1; ++i) {
                    events_mutation_[index_mutation] = generator_->GetPoisson(time_step_ * propensities_mutation_[index_mutation]);

                    infectious_delta_[getIndexHap(population, GetNewHaplotype(haplotype, site, i, getNumberSites()))] += events_mutation_[index_mutation];
                    infectious_delta_[getIndexHap(population, haplotype)] -= events_mutation_[index_mutation];
                    ++index_mutation;
                }
            }
        }
    }

    uint64_t index_migration = 0;
    for (uint64_t source = 0; source < getNumberPopulations(); ++source) {
        for (uint64_t target = 0; target < getNumberPopulations(); ++target) {
            if (source == target) {
                continue;
            }
            for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
                for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
                    events_migration_[index_migration] = generator_->GetPoisson(time_step_ * propensities_migration_[index_migration]);

                    infectious_delta_[getIndexHap(target, haplotype)] += events_migration_[index_migration];
                    susceptibles_delta_[getIndexSus(target, group)] -= events_migration_[index_migration];
                    ++index_migration;
                }
            }
        }
    }

    uint64_t index_immunity = 0;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        for (uint64_t source = 0; source < getNumberSusceptibleGroups(); ++source) {
            for (uint64_t target = 0; target < getNumberSusceptibleGroups(); ++target) {
                if (source == target) {
                    continue;
                }
                events_immunity_[index_immunity] = generator_->GetPoisson(time_step_ * propensities_immunity_[index_immunity]);

                susceptibles_delta_[getIndexSus(population, target)] += events_immunity_[index_immunity];
                susceptibles_delta_[getIndexSus(population, source)] -= events_immunity_[index_immunity];
                ++index_immunity;
            }
        }
    }

    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        uint64_t* size = pool_->GetSizeBegin();
        uint64_t* infected = pool_->GetInfectedPopHapBegin(population);
        for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
            if (infectious_delta_[getIndexHap(population, haplotype)] + infected[haplotype] < 0 || infectious_delta_[getIndexHap(population, haplotype)] + infected[haplotype] > size[population]) {
                return false;
            }
        }

        uint64_t* susceptible = pool_->GetSusceptiblesPopSusBegin(population);
        for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
            if (susceptibles_delta_[getIndexSus(population, group)] + susceptible[group] < 0 || susceptibles_delta_[getIndexSus(population, group)] + susceptible[group] > size[population]) {
                return false;
            }
        }
    }

    return true;
}

void Tau::UpdateCompartmentCounts() {
    uint64_t count_event = 0;

    uint64_t index_transmission = 0;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
            for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
                if (events_transmission_[index_transmission] != 0) {
                    pool_->NewInfections(events_transmission_[index_transmission], haplotype, group, population);
                    chain_->AddMultievent({{TypeEvents::kTRANSMISSION, haplotype, population, group, 0}, events_transmission_[index_transmission]});
                    counters_->AddTransmission(events_transmission_[index_transmission]);
                    ++count_event;
                }
                ++index_transmission;
            }
        }
    }

    uint64_t index_recovery_sampling = 0;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
            if (events_recovery_[index_recovery_sampling] != 0) {
                pool_->NewRecoveries(events_recovery_[index_recovery_sampling], haplotype, infectious_data_->GetSusceptibilityTypes(haplotype), population);
                chain_->AddMultievent({{TypeEvents::kRECOVERY, haplotype, population, infectious_data_->GetSusceptibilityTypes(haplotype), 0}, events_recovery_[index_recovery_sampling]});
                counters_->AddRecovery(events_recovery_[index_recovery_sampling]);
                ++count_event;
            
            }
            if (events_sampling_[index_recovery_sampling] != 0) {
                pool_->NewRecoveries(events_sampling_[index_recovery_sampling], haplotype, infectious_data_->GetSusceptibilityTypes(haplotype), population);
                chain_->AddMultievent({{TypeEvents::kSAMPLING, haplotype, population, infectious_data_->GetSusceptibilityTypes(haplotype), 0}, events_sampling_[index_recovery_sampling]});
                counters_->AddSampling(events_sampling_[index_recovery_sampling]);
                ++count_event;
            }
            ++index_recovery_sampling;
        }
    }

    uint64_t index_mutation = 0;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
            for (uint64_t site = 0; site < getNumberSites(); ++site) {
                for (uint64_t i = 0; i < getNumberAlleleStates() - 1; ++i) {
                    if (events_mutation_[index_mutation] != 0) {
                        pool_->NewMutation(events_mutation_[index_mutation], haplotype, GetNewHaplotype(haplotype, site, i, getNumberSites()), population);
                        chain_->AddMultievent({{TypeEvents::kMUTATION, haplotype, population, GetNewHaplotype(haplotype, site, i, getNumberSites()), 0}, events_mutation_[index_mutation]});
                        counters_->AddMutation(events_mutation_[index_mutation]);
                        ++count_event;
                    }
                    ++index_mutation;
                }
            }
        }
    }

    uint64_t index_migration = 0;
    for (uint64_t source = 0; source < getNumberPopulations(); ++source) {
        for (uint64_t target = 0; target < getNumberPopulations(); ++target) {
            if (source == target) {
                continue;
            }
            for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
                for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
                    if (events_migration_[index_migration] != 0) {
                        pool_->NewInfections(events_migration_[index_migration], haplotype, group, target);
                        chain_->AddMultievent({{TypeEvents::kMIGRATION, haplotype, source, group, target}, events_migration_[index_migration]});
                        counters_->AddMigrationAccept(events_migration_[index_migration]);
                        ++count_event;
                    }
                    ++index_migration;
                }
            }
        }
    }

    uint64_t index_immunity = 0;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        for (uint64_t source = 0; source < getNumberSusceptibleGroups(); ++source) {
            for (uint64_t target = 0; target < getNumberSusceptibleGroups(); ++target) {
                if (source == target) {
                    continue;
                }
                if (events_immunity_[index_immunity] != 0) {
                    pool_->NewImmunity(events_immunity_[index_immunity], source, target, population);
                    chain_->AddMultievent({{TypeEvents::kSUSCCHANGE, source, population, target, 0}, events_immunity_[index_immunity]});
                    counters_->AddImmunity(events_immunity_[index_immunity]);
                    ++count_event;
                }
                ++index_immunity;
            }
        }
    }

    chain_->AddEvent({TypeEvents::kMULTITYPE, chain_->GetLastMultievent(), count_event, 0, 0});
}

inline uint64_t Tau::getNumberSites() const {
    return numbers_.sites;
}

inline uint64_t Tau::getNumberAlleleStates() const {
    return numbers_.allele_states;
}

inline uint64_t Tau::getNumberHaplotypes() const {
    return numbers_.haplotypes;
}

inline uint64_t Tau::getNumberPopulations() const {
    return numbers_.populations;
}

inline uint64_t Tau::getNumberSusceptibleGroups() const {
    return numbers_.susceptible_groups;
}

inline uint64_t Tau::getIndexHap(uint64_t first, uint64_t second) const {
    return first * getNumberHaplotypes() + second;
}

inline uint64_t Tau::getIndexSus(uint64_t first, uint64_t second) const {
    return first * getNumberSusceptibleGroups() + second;
}

inline uint64_t Tau::getIndexSit(uint64_t first, uint64_t second) const {
    return first * getNumberSites() + second;
}
