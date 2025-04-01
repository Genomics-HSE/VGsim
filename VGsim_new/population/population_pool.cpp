#pragma once

#include "../utils/utils.cpp"
#include "population.cpp"

#include "population_pool.h"

PopulationPool::PopulationPool(uint64_t number_of_populations, uint64_t number_of_haplotypes, uint64_t number_of_susceptible_groups)
    : number_of_haplotypes_(number_of_haplotypes)
    , number_of_populations_(number_of_populations)
    , number_of_susceptible_groups_(number_of_susceptible_groups)
    , infected_(0)
    , susceptibles_(1'000'000 * number_of_populations_)
    , sizes_(new uint64_t[number_of_populations_])
    , infected_pop_(new uint64_t[number_of_populations_])
    , susceptibles_pop_(new uint64_t[number_of_populations_])
    , actual_sizes_(new double[number_of_populations_])
    , max_effective_migration_(new double[number_of_populations_])
    , sampling_multiplier_(new double[number_of_populations_])
    , migration_probability_(new double[number_of_populations_ * number_of_populations])
    , effective_migration_probability_(new double[number_of_populations_ * number_of_populations])
    , multiplier_migration_probability_(new double[number_of_populations_ * number_of_populations_])
    , populations_(new Population[number_of_populations_]) {
    uint64_t base_size = 1'000'000;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        populations_[population].SetParameters(number_of_haplotypes_, number_of_susceptible_groups_);
        populations_[population].SetSize(base_size);
        actual_sizes_[population] = base_size;
        sizes_[population] = base_size;
        infected_pop_[population] = 0;
        susceptibles_pop_[population] = base_size;
        sampling_multiplier_[population] = 1.0;
    }
    double probability = 0.0;
    for (uint64_t source = 0; source < getNumberPopulations(); ++source) {
        for (uint64_t target = 0; target < getNumberPopulations(); ++target) {
            migration_probability_[getIndexPop(source, target)] = source == target ? 1.0 : probability;
        }
    }
}

PopulationPool::~PopulationPool() {
    delete[] sizes_;
    delete[] infected_pop_;
    delete[] susceptibles_pop_;

    delete[] actual_sizes_;
    delete[] max_effective_migration_;
    delete[] sampling_multiplier_;
    delete[] migration_probability_;
    delete[] effective_migration_probability_;
    delete[] multiplier_migration_probability_;
}

void PopulationPool::Debug() {
    std::cout << "POPULATION POOL" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Infected: " << infected_ << std::endl;
    std::cout << "Susceptibles: " << susceptibles_ << std::endl;
    PrintArray1nd("Infected population:", infected_pop_, getNumberPopulations());
    PrintArray1nd("Susceptibles population:", susceptibles_pop_, getNumberPopulations());
    PrintArray1nd("Actual sizes:", actual_sizes_, getNumberPopulations());
    PrintArray1nd("Max effective migration:", max_effective_migration_, getNumberPopulations());
    PrintArray2nd("Migration probability", migration_probability_, getNumberPopulations(), getNumberPopulations());
    PrintArray2nd("Effective migration probability", effective_migration_probability_, getNumberPopulations(), getNumberPopulations());
    PrintArray2nd("Multiplier migration rate", multiplier_migration_probability_, getNumberPopulations(), getNumberPopulations());

    std::cout << "Size:";
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        std::cout << " " << populations_[population].GetSize();
    }
    std::cout << std::endl;
    std::cout << "Sampling multiplier:";
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        std::cout << " " << GetSamplingMultiplier(population);
    }
    std::cout << std::endl;
    std::cout << "Susceptibles population" << std::endl;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        uint64_t* susceptibles = populations_[population].GetSusceptiblesBegin();
        for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
            std::cout << susceptibles[group] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Infected population" << std::endl;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        uint64_t* infected = populations_[population].GetInfectedBegin();
        for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
            std::cout << infected[haplotype] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Lockdown on/off:";
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        std::cout << " " << populations_[population].GetOn();
    }
    std::cout << std::endl;
    std::cout << "Contact density before:";
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        std::cout << " " << populations_[population].GetContactDensityBefore();
    }
    std::cout << std::endl;
    std::cout << "Contact density after:";
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        std::cout << " " << populations_[population].GetContactDensityAfter();
    }
    std::cout << std::endl;
    std::cout << "Start:";
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        std::cout << " " << populations_[population].GetStart();
    }
    std::cout << std::endl;
    std::cout << "End:";
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        std::cout << " " << populations_[population].GetEnd();
    }
    std::cout << std::endl;
    std::cout << std::endl;
}

boost::python::tuple PopulationPool::get_current_individuals() {
    boost::python::list susceptibles;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        uint64_t* susceptibles_data = populations_[population].GetSusceptiblesBegin();
        boost::python::list data;
        for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
            data.append(boost::python::object(susceptibles_data[group]));
        }
        susceptibles.append(data);
    }

    boost::python::list infected;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        uint64_t* infected_data = populations_[population].GetInfectedBegin();
        boost::python::list data;
        for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
            data.append(boost::python::object(infected_data[haplotype]));
        }
        infected.append(data);
    }

    return boost::python::make_tuple(susceptibles, infected);
}

boost::python::list PopulationPool::get_actual_size() {
    boost::python::list data;

    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        data.append(boost::python::object(actual_sizes_[population]));
    }

    return data;
}

boost::python::list PopulationPool::get_contact_density_before_lockdown() {
    boost::python::list data;

    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        data.append(boost::python::object(populations_[population].GetContactDensityBefore()));
    }

    return data;
}

void PopulationPool::set_population_size(uint64_t size, uint64_t population) {
    sizes_[population] = size;
    susceptibles_pop_[population] = size;
    populations_[population].SetSize(size);
}

boost::python::list PopulationPool::get_population_size() {
    boost::python::list data;

    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        data.append(boost::python::object(sizes_[population]));
    }

    return data;
}

void PopulationPool::set_contact_density(double value, uint64_t population) {
    populations_[population].SetContactDensityBefore(value);
}

boost::python::list PopulationPool::get_contact_density() {
    boost::python::list data;

    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        data.append(boost::python::object(populations_[population].GetContactDensity()));
    }

    return data;
}

void PopulationPool::set_npi(double after, double start, double end, uint64_t population) {
    populations_[population].SetContactDensityAfter(after);
    populations_[population].SetStart(start);
    populations_[population].SetEnd(end);
}

boost::python::list PopulationPool::get_npi() {
    boost::python::list data;

    boost::python::list str_data;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        str_data.append(boost::python::object(populations_[population].GetContactDensityAfter()));
    }
    data.append(str_data);

    boost::python::list str_data_2;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        str_data_2.append(boost::python::object(populations_[population].GetStart()));
    }
    data.append(str_data_2);

    boost::python::list str_data_3;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        str_data_3.append(boost::python::object(populations_[population].GetEnd()));
    }
    data.append(str_data_3);

    return data;
}

void PopulationPool::set_sampling_multiplier(double multiplier, uint64_t population) {
    sampling_multiplier_[population] = multiplier;
}

boost::python::list PopulationPool::get_sampling_multiplier() {
    boost::python::list data;

    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        data.append(boost::python::object(sampling_multiplier_[population]));
    }

    return data;
}

void PopulationPool::set_migration_probability(double probability, uint64_t source_population, uint64_t target_population) {
    migration_probability_[getIndexPop(source_population, target_population)] = probability;
}

boost::python::list PopulationPool::get_migration_probability() {
    boost::python::list data;

    for (uint64_t source = 0; source < getNumberPopulations(); ++source) {
        boost::python::list str_data;
        for (uint64_t target = 0; target < getNumberPopulations(); ++target) {
            str_data.append(boost::python::object(migration_probability_[getIndexPop(source, target)]));
        }
        data.append(str_data);
    }

    return data;
}

uint64_t PopulationPool::check_migration_probability() {
    for (uint64_t source = 0; source < getNumberPopulations(); ++source) {
        double summa = 0.0;
        migration_probability_[getIndexPop(source, source)] = 1.0;
        for (uint64_t target = 0; target < getNumberPopulations(); ++target) {
            if (source == target) {
                continue;
            }
            summa += migration_probability_[getIndexPop(source, target)];
            migration_probability_[getIndexPop(source, source)] -= migration_probability_[getIndexPop(source, target)];
        }
        if (summa > 1.0) {
            return 1;
        }
    }
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        if (migration_probability_[getIndexPop(population, population)] <= 1e-15) {
            return 2;
        }
    }
    return 0;
}


void PopulationPool::SaveInfections() {
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        populations_[population].Save();
    }
}

void PopulationPool::FirstInfections() {
    if (GetInfected() != 0) {
        return;
    }
    NewInfections(1, 0, 0, 0);
}

void PopulationPool::Restart() {
    infected_ = 0;
    susceptibles_ = 0;
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        populations_[population].Reset();
        uint64_t* infected_begin = populations_[population].GetInfectedBegin();
        infected_pop_[population] = 0;
        for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
            infected_pop_[population] += infected_begin[haplotype];
        }
        infected_ += infected_pop_[population];
        uint64_t* susceptibles_begin = populations_[population].GetSusceptiblesBegin();
        susceptibles_pop_[population] = 0;
        for (uint64_t group = 0; group < getNumberSusceptibleGroups(); ++group) {
            susceptibles_pop_[population] += susceptibles_begin[group];
        }
        susceptibles_ += susceptibles_pop_[population];
    }
}

void PopulationPool::Update() {
    for (uint64_t source = 0; source < getNumberPopulations(); ++source) {
        migration_probability_[getIndexPop(source, source)] = 1.0;
        actual_sizes_[source] = 0.0;
        for (uint64_t target = 0; target < getNumberPopulations(); ++target) {
            if (source == target) {
                continue;
            }
            migration_probability_[getIndexPop(source, source)] -= migration_probability_[getIndexPop(source, target)];
            actual_sizes_[source] +=
                migration_probability_[getIndexPop(target, source)] *
                populations_[target].GetSize();
        }
        actual_sizes_[source] +=
            migration_probability_[getIndexPop(source, source)] *
            populations_[source].GetSize();
    }

    UpdateRates();
}


inline void PopulationPool::NewInfections(uint64_t count, uint64_t haplotype, uint64_t group, uint64_t population) {
    infected_ += count;
    susceptibles_ -= count;
    infected_pop_[population] += count;
    susceptibles_pop_[population] -= count;
    populations_[population].NewInfections(count, haplotype, group);
}

inline void PopulationPool::NewRecoveries(uint64_t count, uint64_t haplotype, uint64_t group, uint64_t population) {
    infected_ -= count;
    susceptibles_ += count;
    infected_pop_[population] -= count;
    susceptibles_pop_[population] += count;
    populations_[population].NewRecoveries(count, haplotype, group);
}

inline void PopulationPool::NewMutation(uint64_t count, uint64_t old_haplotype, uint64_t new_haplotype, uint64_t population) {
    populations_[population].NewMutation(count, old_haplotype, new_haplotype);
}

inline void PopulationPool::NewImmunity(uint64_t count, uint64_t old_group, uint64_t new_group, uint64_t population) {
    populations_[population].NewImmunity(count, old_group, new_group);
}

void PopulationPool::CheckLockdown(uint64_t population) {
    if ((populations_[population].GetOn() && static_cast<double>(infected_pop_[population]) / populations_[population].GetSize() < populations_[population].GetEnd())
        || (!populations_[population].GetOn() && static_cast<double>(infected_pop_[population]) / populations_[population].GetSize() > populations_[population].GetStart())) {
        populations_[population].Switch();
        UpdateRates();
    }
}

std::vector<std::vector<uint64_t>> PopulationPool::GetInfectious() {
    std::vector<std::vector<uint64_t>> infectious(getNumberPopulations(), std::vector<uint64_t>(getNumberHaplotypes()));
    for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
        for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
            infectious[population][haplotype] = GetInfectedPopHap(population, haplotype);
        }
    }
    return infectious;
}

inline double PopulationPool::GetContactDensity(uint64_t population) const {
    return populations_[population].GetContactDensity();
}

inline double* PopulationPool::GetMaxEffectiveMigration() const {
    return max_effective_migration_;
}

inline double PopulationPool::GetMigrationProbability(uint64_t source, uint64_t target) const {
    return migration_probability_[getIndexPop(source, target)];
}

inline double PopulationPool::GetEffectiveMigration(uint64_t source, uint64_t target) const {
    return effective_migration_probability_[getIndexPop(source, target)];
}

inline double PopulationPool::GetMultiplierMigrationProbability(uint64_t source, uint64_t target) const {
    return multiplier_migration_probability_[getIndexPop(source, target)];
}

inline double PopulationPool::GetSamplingMultiplier(uint64_t population) const {
    return sampling_multiplier_[population];
}

inline double PopulationPool::GetActualSize(uint64_t population) const {
    return actual_sizes_[population];
}

inline uint64_t* PopulationPool::GetSizeBegin() const {
    return sizes_;
}

inline uint64_t PopulationPool::GetInfected() const {
    return infected_;
}

inline uint64_t* PopulationPool::GetInfectedPopBegin() const {
    return infected_pop_;
}

inline uint64_t PopulationPool::GetInfectedPop(uint64_t population) const {
    return infected_pop_[population];
}

inline uint64_t* PopulationPool::GetInfectedPopHapBegin(uint64_t population) const {
    return populations_[population].GetInfectedBegin();
}

inline uint64_t PopulationPool::GetInfectedPopHap(uint64_t population, uint64_t haplotype) const {
    return populations_[population].GetInfected(haplotype);
}

inline uint64_t PopulationPool::GetSusceptibles() const {
    return susceptibles_;
}

inline uint64_t* PopulationPool::GetSusceptiblesPopBegin() const {
    return susceptibles_pop_;
}

inline uint64_t PopulationPool::GetSusceptiblesPop(uint64_t population) const {
    return susceptibles_pop_[population];
}

inline uint64_t* PopulationPool::GetSusceptiblesPopSusBegin(uint64_t population) const {
    return populations_[population].GetSusceptiblesBegin();
}

inline uint64_t PopulationPool::GetSusceptiblesPopSus(uint64_t population, uint64_t group) const {
    return populations_[population].GetSusceptibles(group);
}

uint64_t PopulationPool::GetStartNumberSusceptible(uint64_t population, uint64_t group) const {
    return populations_[population].GetStartNumberSusceptible(group);
}

uint64_t PopulationPool::GetStartNumberInfected(uint64_t population, uint64_t haplotype) const {
    return populations_[population].GetStartNumberInfected(haplotype);
}


// Private
void PopulationPool::UpdateRates() {
    for (uint64_t source = 0; source < getNumberPopulations(); ++source) {
        for (uint64_t target = 0; target < getNumberPopulations(); ++target) {
            effective_migration_probability_[getIndexPop(source, target)] = 0.0;
            if (source == target) {
                continue;
            }
            for (uint64_t middle = 0; middle < getNumberPopulations(); ++middle) {
                effective_migration_probability_[getIndexPop(source, target)] += 
                    migration_probability_[getIndexPop(source, middle)] *
                    migration_probability_[getIndexPop(target, middle)] *
                    populations_[middle].GetContactDensity() /
                    actual_sizes_[middle];
            }
            max_effective_migration_[target] = std::max(effective_migration_probability_[getIndexPop(source, target)], max_effective_migration_[target]);
        }
    }

    for (uint64_t source = 0; source < getNumberPopulations(); ++source) {
        for (uint64_t target = 0; target < getNumberPopulations(); ++target) {
            multiplier_migration_probability_[getIndexPop(source, target)] =
                migration_probability_[getIndexPop(source, target)] *
                migration_probability_[getIndexPop(source, target)] *
                populations_[target].GetContactDensity() /
                actual_sizes_[target];
        }
    }
}

inline uint64_t PopulationPool::getNumberPopulations() const {
    return number_of_populations_;
}

inline uint64_t PopulationPool::getNumberHaplotypes() const {
    return number_of_haplotypes_;
}

inline uint64_t PopulationPool::getNumberSusceptibleGroups() const {
    return number_of_susceptible_groups_;
}

inline uint64_t PopulationPool::getIndexPop(uint64_t first, uint64_t second) const {
    return first * getNumberPopulations() + second;
}
