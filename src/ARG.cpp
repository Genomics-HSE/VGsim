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
    uint64_t mt_ev_num;

    uint64_t size = chain_->GetSize();
    for (uint64_t index = 1; index <= size; ++index) {
        Event event = chain_->GetEvent(size - index);
        double time = chain_->GetTime(size - index);
        switch (static_cast<TypeEvents>(event.type)) {
            case TypeEvents::kTRANSMISSION: {
                haplotype = event.parameter1;
                population = event.parameter2;

                int64_t countBranches = liveBranches[population][haplotype].size();
                int64_t countInfectious = infectious[population][haplotypeToIndex[haplotype]];
                double probability = static_cast<double>(countBranches) * (countBranches - 1) / countInfectious / (countInfectious - 1);
                if (generator_->GetUniform() < probability) {
                    int64_t leftChildIndex = static_cast<int64_t>(countBranches * generator_->GetUniform());
                    int64_t rightChildIndex = static_cast<int64_t>((countBranches - 1) * generator_->GetUniform());
                    if (rightChildIndex >= leftChildIndex) {
                        rightChildIndex += 1;
                    }
                    int64_t leftChild = liveBranches[population][haplotype][leftChildIndex];
                    int64_t rightChild = liveBranches[population][haplotype][rightChildIndex];
                    int64_t parent = currentNode++;
                    addNode(population, time, parent, leftChild, rightChild);
                    liveBranches[population][haplotype][leftChildIndex] = parent;
                    liveBranches[population][haplotype][rightChildIndex] = liveBranches[population][haplotype][countBranches - 1];
                    liveBranches[population][haplotype].pop_back();
                }
                infectious[population][haplotypeToIndex[haplotype]] -= 1;
                break;
            }
            case TypeEvents::kRECOVERY: {
                haplotype = event.parameter1;
                population = event.parameter2;

                infectious[population][haplotypeToIndex[haplotype]] += 1;
                break;
            }
            case TypeEvents::kSAMPLING: {
                haplotype = event.parameter1;
                population = event.parameter2;

                addNode(population, time);
                liveBranches[population][haplotype].push_back(currentNode++);
                infectious[population][haplotypeToIndex[haplotype]] += 1;
                break;
            }
            case TypeEvents::kMUTATION: {
                oldHaplotype = event.parameter1;
                population = event.parameter2;
                newHaplotype = event.parameter3;

                uint64_t countBranches = liveBranches[population][newHaplotype].size();
                double probability = static_cast<double>(countBranches) / infectious[population][newHaplotype];
                if (generator_->GetUniform() < probability) {
                    uint64_t nodeIndex = static_cast<uint64_t>(countBranches * generator_->GetUniform());
                    uint64_t node = liveBranches[population][newHaplotype][nodeIndex];
                    liveBranches[population][newHaplotype][nodeIndex] = liveBranches[population][newHaplotype][countBranches - 1];
                    liveBranches[population][newHaplotype].pop_back();
                    liveBranches[population][oldHaplotype].push_back(node);
                    addMutation(node, oldHaplotype, newHaplotype, time);
                }
                infectious[population][newHaplotype] -= 1;
                infectious[population][oldHaplotype] += 1;
                break;
            }
            case TypeEvents::kMIGRATION: {
                haplotype = event.parameter1;
                source_population = event.parameter2;
                target_population = event.parameter4;

                int64_t countBranchesTarget = liveBranches[target_population][haplotype].size();
                int64_t countBranchesSource = liveBranches[source_population][haplotype].size();
                double probability = static_cast<double>(countBranchesTarget) / infectious[target_population][haplotype];
                if (generator_->GetUniform() < probability) {
                    int64_t nodeTarget = static_cast<int64_t>(countBranchesTarget * generator_->GetUniform());
                    probability = static_cast<double>(countBranchesSource) / infectious[source_population][haplotype];
                    if (generator_->GetUniform() < probability) {
                        int64_t nodeSource = static_cast<int64_t>(countBranchesSource * generator_->GetUniform());
                        int64_t targetChild = liveBranches[target_population][haplotype][nodeTarget];
                        int64_t sourceChild = liveBranches[source_population][haplotype][nodeSource];
                        int64_t parent = currentNode++;
                        liveBranches[source_population][haplotype][nodeSource] = parent;
                        liveBranches[target_population][haplotype][nodeTarget] = liveBranches[target_population][haplotype][countBranchesTarget - 1];
                        liveBranches[target_population][haplotype].pop_back();
                        addNode(source_population, time, parent, targetChild, sourceChild);
                        addMigration(targetChild, {source_population, target_population, time});
                    } else {
                        liveBranches[event.parameter2][haplotype].push_back(liveBranches[target_population][haplotype][nodeTarget]);
                        liveBranches[target_population][haplotype][nodeTarget] = liveBranches[target_population][haplotype][countBranchesTarget - 1];
                        liveBranches[target_population][haplotype].pop_back();
                    }
                }
                infectious[target_population][haplotype] -= 1;
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
                            // lbs = liveBranchesS[me_population][me_haplotype].size()
                            // lbs_e = self.infectious[me_population, me_haplotype]
                            // if me_num == 0 or lbs == 0:
                            //     mt_ev_num = 0
                            // else:
                            //     mt_ev_num = random_hypergeometric(self.seed.rng, int(lbs*(lbs-1.0)/2.0), int(lbs_e*(lbs_e-1)/2-lbs*(lbs-1)/2), me_num)
                            // for i in range(mt_ev_num):
                            //     n1 = int(floor( lbs*self.seed.uniform() ))
                            //     n2 = int(floor( (lbs-1)*self.seed.uniform() ))
                            //     if n2 >= n1:
                            //         n2 += 1
                            //     id1 = liveBranchesS[me_population][me_haplotype][n1]
                            //     id2 = liveBranchesS[me_population][me_haplotype][n2]
                            //     id3 = ptrTreeAndTime
                            //     newLineages[me_population][me_haplotype].push_back(id3)
                            //     if n1 == lbs-1:#TODO: need to check if we really need these if and elif. Most likely - no!
                            //         liveBranchesS[me_population][me_haplotype].pop_back()
                            //         liveBranchesS[me_population][me_haplotype][n2] = liveBranchesS[me_population][me_haplotype][lbs-2]
                            //         liveBranchesS[me_population][me_haplotype].pop_back()
                            //     elif n2 == lbs-1:
                            //         liveBranchesS[me_population][me_haplotype].pop_back()
                            //         liveBranchesS[me_population][me_haplotype][n1] = liveBranchesS[me_population][me_haplotype][lbs-2]
                            //         liveBranchesS[me_population][me_haplotype].pop_back()
                            //     else:
                            //         liveBranchesS[me_population][me_haplotype][n1] = liveBranchesS[me_population][me_haplotype][lbs-1]
                            //         liveBranchesS[me_population][me_haplotype].pop_back()
                            //         liveBranchesS[me_population][me_haplotype][n2] = liveBranchesS[me_population][me_haplotype][lbs-2]
                            //         liveBranchesS[me_population][me_haplotype].pop_back()
                            //     self.tree[id1] = id3
                            //     self.tree[id2] = id3
                            //     self.tree[ptrTreeAndTime] = -1
                            //     self.tree_pop[ptrTreeAndTime] = me_population
                            //     self.times[ptrTreeAndTime] = me_time
                            //     ptrTreeAndTime += 1
                            //     lbs -= 2
                            // self.infectiousDelta[me_population, me_haplotype] -= me_num


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
                                addNode(population, time);
                                newLineages[population][haplotype].push_back(currentNode++);
                            }
                            infectiousDelta[getIndexHap(population, haplotype)] += event_m.count;
                            break;
                        case TypeEvents::kMUTATION:
                            oldHaplotype = event_m.parameter1;
                            population = event_m.parameter2;
                            newHaplotype = event_m.parameter3;

                            countBranches = liveBranches[population][newHaplotype].size();
                            mt_ev_num = 0;
                            // if (countBranches != 0) {
                            //     mt_ev_num = generator_->GetHypergeometric(countBranches, infectious[population][newHaplotype] - countBranches, event_m.count);
                            // }

                            for (uint64_t _ = 0; _ < mt_ev_num; ++_) {
                                uint64_t nodeIndex = static_cast<uint64_t>(countBranches * generator_->GetUniform());
                                uint64_t node = liveBranches[population][newHaplotype][nodeIndex];
                                liveBranches[population][newHaplotype][nodeIndex] = liveBranches[population][newHaplotype][countBranches - 1];
                                liveBranches[population][newHaplotype].pop_back();
                                newLineages[population][oldHaplotype].push_back(node);
                                addMutation(node, oldHaplotype, newHaplotype, time);
                                --countBranches;
                            }
                            infectiousDelta[getIndexHap(population, oldHaplotype)] += event_m.count;
                            infectiousDelta[getIndexHap(population, newHaplotype)] -= event_m.count;
                            break;
                        case TypeEvents::kMIGRATION:
                            
                            break;
                        case TypeEvents::kSUSCCHANGE:
                            break;
                        case TypeEvents::kMULTITYPE:
                            break;
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

void ARG::restart() {
    uint64_t countSampling = counters_->GetSampling();
    time_.resize(0);
    time_.reserve(countSampling);
    tree_.resize(0);
    tree_.reserve(countSampling);
    population_.resize(0);
    population_.reserve(countSampling);
}

void ARG::addNode(int64_t population, double time, int64_t parent, int64_t left, int64_t right) {
    tree_.push_back(-1);
    population_.push_back(population);
    time_.push_back(time);

    if (left != -1 && right != -1) {
        tree_[left] = parent;
        tree_[right] = parent;
    }
}

void ARG::addMutation(uint64_t node, uint64_t oldHaplotype, uint64_t newHaplotype, double time) {
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

void ARG::addMigration(uint64_t node, Migration migration) {
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
