#pragma once

#include "arg.h"

ARG::ARG(Numbers numbers, Chain* chain, Counters* counters, PopulationPool* pool, RandomGenerator* generator)
    : numbers_(numbers)
    , chain_(chain)
    , counters_(counters)
    , pool_(pool)
    , generator_(generator) {
}

void ARG::CalculateGenealogy() {
    restart();
    std::vector<std::vector<std::vector<uint64_t>>> liveBranches(getNumberPopulations(), std::vector<std::vector<uint64_t>>(getNumberHaplotypes()));
    std::vector<std::vector<uint64_t>> infectious = pool_->GetInfectious();
    uint64_t currentNode = 0;

    std::vector<uint64_t> haplotypeToIndex(getNumberHaplotypes());
    for (uint64_t i = 0; i < getNumberHaplotypes(); ++i) {
        haplotypeToIndex[i] = i;
    }

    std::vector<std::vector<std::vector<uint64_t>>> newLineages(getNumberPopulations(), std::vector<std::vector<uint64_t>>(getNumberHaplotypes()));
    std::vector<uint64_t> infectiousDelta(getNumberPopulations() * getNumberHaplotypes(), 0);


    uint64_t haplotype;
    uint64_t population;
    uint64_t oldHaplotype;
    uint64_t newHaplotype;
    uint64_t source_population;
    uint64_t target_population;
    uint64_t start;
    uint64_t count;
    uint64_t countBranches;
    uint64_t countInfectious;
    uint64_t countBranchesTarget;
    uint64_t countBranchesSource;

    uint64_t leftChildIndex;
    uint64_t rightChildIndex;
    uint64_t leftChild;
    uint64_t rightChild;
    uint64_t parent;

    uint64_t node;
    uint64_t nodeIndex;
    uint64_t nodeTarget;
    uint64_t nodeSource;
    uint64_t targetChild;
    uint64_t sourceChild;

    uint64_t mt_ev_num;
    uint64_t mt_ev_num2;

    double probability;

    for (uint64_t index = chain_->GetSize(); index > 0; --index) {
        Event event = chain_->GetEvent(index - 1);
        double time = chain_->GetTime(index - 1);
        switch (static_cast<TypeEvents>(event.type)) {
            case TypeEvents::kTRANSMISSION: {
                haplotype = event.parameter1;
                population = event.parameter2;

                countBranches = liveBranches[population][haplotype].size();
                countInfectious = infectious[population][haplotypeToIndex[haplotype]];
                probability = static_cast<double>(countBranches) * (countBranches - 1) / countInfectious / (countInfectious - 1);
                if (generator_->GetUniform() < probability) {
                    leftChildIndex = static_cast<int64_t>(countBranches * generator_->GetUniform());
                    rightChildIndex = static_cast<int64_t>((countBranches - 1) * generator_->GetUniform());
                    if (rightChildIndex >= leftChildIndex) {
                        rightChildIndex += 1;
                    }
                    leftChild = liveBranches[population][haplotype][leftChildIndex];
                    rightChild = liveBranches[population][haplotype][rightChildIndex];
                    parent = currentNode++;
                    addNode(population, time, parent, leftChild, rightChild);
                    liveBranches[population][haplotype][leftChildIndex] = parent;
                    removeBranch(liveBranches, population, haplotype, rightChildIndex);
                }
                --infectious[population][haplotypeToIndex[haplotype]];
                break;
            }
            case TypeEvents::kRECOVERY: {
                haplotype = event.parameter1;
                population = event.parameter2;

                ++infectious[population][haplotypeToIndex[haplotype]];
                break;
            }
            case TypeEvents::kSAMPLING: {
                haplotype = event.parameter1;
                population = event.parameter2;

                addLeaf(population, time);
                liveBranches[population][haplotype].push_back(currentNode++);
                ++infectious[population][haplotypeToIndex[haplotype]];
                break;
            }
            case TypeEvents::kMUTATION: {
                oldHaplotype = event.parameter1;
                population = event.parameter2;
                newHaplotype = event.parameter3;

                countBranches = liveBranches[population][newHaplotype].size();
                countInfectious = infectious[population][haplotypeToIndex[newHaplotype]];
                probability = static_cast<double>(countBranches) / countInfectious;
                if (generator_->GetUniform() < probability) {
                    nodeIndex = static_cast<uint64_t>(countBranches * generator_->GetUniform());
                    node = liveBranches[population][newHaplotype][nodeIndex];
                    removeBranch(liveBranches, population, newHaplotype, nodeIndex);
                    liveBranches[population][oldHaplotype].push_back(node);
                    addMutation(node, oldHaplotype, newHaplotype, time);
                }
                --infectious[population][newHaplotype];
                ++infectious[population][oldHaplotype];
                break;
            }
            case TypeEvents::kMIGRATION: {
                haplotype = event.parameter1;
                source_population = event.parameter2;
                target_population = event.parameter4;

                countBranchesTarget = liveBranches[target_population][haplotype].size();
                countBranchesSource = liveBranches[source_population][haplotype].size();
                probability = static_cast<double>(countBranchesTarget) / infectious[target_population][haplotype];
                if (generator_->GetUniform() < probability) {
                    nodeTarget = static_cast<int64_t>(countBranchesTarget * generator_->GetUniform());
                    probability = static_cast<double>(countBranchesSource) / infectious[source_population][haplotype];
                    if (generator_->GetUniform() < probability) {
                        nodeSource = static_cast<int64_t>(countBranchesSource * generator_->GetUniform());
                        targetChild = liveBranches[target_population][haplotype][nodeTarget];
                        sourceChild = liveBranches[source_population][haplotype][nodeSource];
                        parent = currentNode++;
                        liveBranches[source_population][haplotype][nodeSource] = parent;
                        removeBranch(liveBranches, target_population, haplotype, nodeTarget);
                        addNode(source_population, time, parent, targetChild, sourceChild);
                        addMigration(targetChild, {source_population, target_population, time});
                    } else {
                        liveBranches[source_population][haplotype].push_back(liveBranches[target_population][haplotype][nodeTarget]);
                        removeBranch(liveBranches, target_population, haplotype, nodeTarget);
                    }
                }
                --infectious[target_population][haplotype];
                break;
            }
            case TypeEvents::kSUSCCHANGE:
                break;
            case TypeEvents::kMULTITYPE:
                start = event.parameter1;
                count = event.parameter2;

                for (uint64_t index_m = start; index_m < start + count; ++index_m) {
                    Multievent event_m = chain_->GetMultievent(index_m);
                    switch (event_m.type) {
                        case TypeEvents::kTRANSMISSION:
                            haplotype = event_m.parameter1;
                            population = event_m.parameter2;

                            countBranches = liveBranches[population][haplotype].size();
                            countInfectious = infectious[population][haplotype];
                            mt_ev_num = generator_->GetHypergeometric(event_m.count, countBranches * (countBranches - 1) / 2, countInfectious * (countInfectious - 1) / 2);
                            for (uint64_t _ = 0; _ < mt_ev_num; ++_) {
                                leftChildIndex = static_cast<int64_t>(countBranches * generator_->GetUniform());
                                rightChildIndex = static_cast<int64_t>((countBranches - 1) * generator_->GetUniform());
                                if (rightChildIndex >= leftChildIndex) {
                                    ++rightChildIndex;
                                }
                                leftChild = liveBranches[population][haplotype][leftChildIndex];
                                rightChild = liveBranches[population][haplotype][rightChildIndex];
                                parent = currentNode++;
                                addNode(population, time, parent, leftChild, rightChild);
                                removeBranch(liveBranches, population, haplotype, leftChildIndex);
                                removeBranch(liveBranches, population, haplotype, rightChildIndex);
                                newLineages[population][haplotype].push_back(parent);
                            }
                            infectiousDelta[getIndexHap(population, haplotype)] += event_m.count;
                            break;
                        case TypeEvents::kRECOVERY:
                            haplotype = event_m.parameter1;
                            population = event_m.parameter2;

                            infectiousDelta[getIndexHap(population, haplotype)] += event_m.count;
                            break;
                        case TypeEvents::kSAMPLING:
                            haplotype = event_m.parameter1;
                            population = event_m.parameter2;

                            for (uint64_t _ = 0; _ < event_m.count; ++_) {
                                addLeaf(population, time);
                                newLineages[population][haplotype].push_back(currentNode++);
                            }
                            infectiousDelta[getIndexHap(population, haplotype)] += event_m.count;
                            break;
                        case TypeEvents::kMUTATION:
                            oldHaplotype = event_m.parameter1;
                            population = event_m.parameter2;
                            newHaplotype = event_m.parameter3;

                            countBranches = liveBranches[population][newHaplotype].size();
                            mt_ev_num = generator_->GetHypergeometric(event_m.count, countBranches, infectious[population][newHaplotype]);
                            for (uint64_t _ = 0; _ < mt_ev_num; ++_) {
                                nodeIndex = static_cast<uint64_t>(countBranches * generator_->GetUniform());
                                node = liveBranches[population][newHaplotype][nodeIndex];
                                removeBranch(liveBranches, population, newHaplotype, nodeIndex);
                                newLineages[population][oldHaplotype].push_back(node);
                                addMutation(node, oldHaplotype, newHaplotype, time);
                                --countBranches;
                            }
                            infectiousDelta[getIndexHap(population, oldHaplotype)] += event_m.count;
                            infectiousDelta[getIndexHap(population, newHaplotype)] -= event_m.count;
                            break;
                        case TypeEvents::kMIGRATION:
                            haplotype = event_m.parameter1;
                            source_population = event_m.parameter2;
                            target_population = event_m.parameter4;

                            countBranchesTarget = liveBranches[target_population][haplotype].size();
                            mt_ev_num = generator_->GetHypergeometric(event_m.count, countBranchesTarget, infectious[target_population][haplotype]);
                            for (uint64_t _ = 0; _ < mt_ev_num; ++_) {
                                countBranchesSource = liveBranches[source_population][haplotype].size();
                                mt_ev_num2 = mt_ev_num = generator_->GetHypergeometric(event_m.count, countBranchesSource, infectious[source_population][haplotype]);
                                for (uint64_t _ = 0; _ < mt_ev_num2; ++_) {
                                    nodeTarget = static_cast<int64_t>(countBranchesTarget * generator_->GetUniform());
                                    nodeSource = static_cast<int64_t>(countBranchesSource * generator_->GetUniform());
                                    targetChild = liveBranches[target_population][haplotype][nodeTarget];
                                    sourceChild = liveBranches[source_population][haplotype][nodeSource];
                                    parent = currentNode++;
                                    addNode(parent, time, parent, leftChild, rightChild);
                                    removeBranch(liveBranches, source_population, haplotype, nodeSource);
                                    removeBranch(liveBranches, target_population, haplotype, nodeTarget);
                                    newLineages[population][haplotype].push_back(parent);
                                    addMigration(targetChild, {source_population, target_population, time});
                                    --countBranchesTarget;
                                    --countBranchesSource;
                                }
                                for (uint64_t _ = 0; _ < mt_ev_num - mt_ev_num2; ++_) {
                                    nodeTarget = static_cast<int64_t>(countBranchesTarget * generator_->GetUniform());
                                    newLineages[source_population][haplotype].push_back(nodeTarget);
                                    removeBranch(liveBranches, source_population, haplotype, nodeTarget);
                                    --countBranchesTarget;
                                }
                            }
                            break;
                        case TypeEvents::kSUSCCHANGE:
                            break;
                        case TypeEvents::kMULTITYPE:
                            break;
                    }
                }
                for (uint64_t population = 0; population < getNumberPopulations(); ++population) {
                    for (uint64_t haplotype = 0; haplotype < getNumberHaplotypes(); ++haplotype) {
                        infectious[population][haplotype] += infectiousDelta[getIndexHap(population, haplotype)];
                        infectiousDelta[getIndexHap(population, haplotype)] = 0;
                        liveBranches[population][haplotype].insert(liveBranches[population][haplotype].end(), newLineages[population][haplotype].begin(), newLineages[population][haplotype].end());
                        newLineages[population][haplotype].clear();
                    }
                }
                break;
        }
    }
    return;
}

void ARG::Debug() {
    std::cout << "Time:";
    for (double time : time_) {
        std::cout << " " << time;
    }
    std::cout << std::endl;
    std::cout << "Tree:";
    for (int64_t node : tree_) {
        std::cout << " " << node;
    }
    std::cout << std::endl;
    std::cout << "Population:";
    for (int64_t population : population_) {
        std::cout << " " << population;
    }
    std::cout << std::endl;
    std::cout << "Mutations\n";
    for (const Mutation& mutation : mutations_) {
        PrintMutation(mutation);
    }
    std::cout << std::endl;
    std::cout << "Mutations\n";
    for (auto migration : migrations_) {
        PrintMigration(migration.first, migration.second);
    }
    std::cout << std::endl;
}

std::vector<double> ARG::get_time() {
    return time_;
}

std::vector<int64_t> ARG::get_tree() {
    return tree_;
}

void ARG::restart() {
    uint64_t countSampling = counters_->GetSampling();
    time_.resize(0);
    time_.reserve(countSampling);
    tree_.resize(0);
    tree_.reserve(countSampling);
    population_.resize(0);
    population_.reserve(countSampling);
}

inline void ARG::removeBranch(std::vector<std::vector<std::vector<uint64_t>>>& branches, uint64_t population, uint64_t haplotype, uint64_t branch) {
    branches[population][haplotype][branch] = branches[population][haplotype].back();
    branches[population][haplotype].pop_back();
}

inline void ARG::addLeaf(int64_t population, double time) {
    tree_.push_back(-1);
    population_.push_back(population);
    time_.push_back(time);
}

inline void ARG::addNode(int64_t population, double time, int64_t parent, int64_t left, int64_t right) {
    addLeaf(population, time);

    tree_[left] = parent;
    tree_[right] = parent;
}

inline void ARG::addMutation(uint64_t node, uint64_t oldHaplotype, uint64_t newHaplotype, double time) {
    uint64_t site = 0;
    uint64_t derivedState = 0;
    uint64_t ancestralState = 0;

    uint64_t del = 1;
    for (uint64_t index = 0; index < getNumberSites(); ++index) {
        if (oldHaplotype % (4 * del) / del != newHaplotype % (4 * del) / del) {
            site = index;
            derivedState = oldHaplotype % (4 * del) / del;
            ancestralState = newHaplotype % (4 * del) / del;
            break;
        }
        del *= 4;
    }

    mutations_.push_back({node, site, derivedState, ancestralState, time});
}

inline void ARG::addMigration(uint64_t node, Migration migration) {
    migrations_.insert({node, migration});
}

inline uint64_t ARG::getNumberSites() const {
    return numbers_.sites;
}

inline uint64_t ARG::getNumberHaplotypes() const {
    return numbers_.haplotypes;
}

inline uint64_t ARG::getNumberPopulations() const {
    return numbers_.populations;
}

inline uint64_t ARG::getNumberSusceptibleGroups() const {
    return numbers_.susceptible_groups;
}

inline uint64_t ARG::getIndexHap(uint64_t first, uint64_t second) const {
    return first * getNumberHaplotypes() + second;
}
