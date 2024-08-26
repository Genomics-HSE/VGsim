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
    std::vector<std::vector<std::vector<uint64_t>>> liveBranches(numbers_.populations, std::vector<std::vector<uint64_t>>(numbers_.haplotypes));
    std::vector<std::vector<uint64_t>> infectious = pool_->GetInfectious();
    uint64_t currentNode = 0;

    std::vector<uint64_t> haplotypeToIndex(numbers_.haplotypes);
    for (uint64_t i = 0; i < numbers_.haplotypes; ++i) {
        haplotypeToIndex[i] = i;
    }

    uint64_t size = chain_->GetSize();

    for (uint64_t index = 1; index <= size; ++index) {
        Event event = chain_->GetEvent(size - index);
        double time = chain_->GetTime(size - index);
        switch (event.type) {
            case kTRANSMISSION: {
                uint64_t haplotype = event.parameter1;
                uint64_t population = event.parameter2;

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
            case kRECOVERY: {
                uint64_t haplotype = event.parameter1;
                uint64_t population = event.parameter2;

                infectious[population][haplotypeToIndex[haplotype]] += 1;
                break;
            }
            case kSAMPLING: {
                uint64_t haplotype = event.parameter1;
                uint64_t population = event.parameter2;

                addNode(population, time);
                liveBranches[population][haplotype].push_back(currentNode++);
                infectious[population][haplotypeToIndex[haplotype]] += 1;
                break;
            }
            case kMUTATION: {
                uint64_t oldHaplotype = event.parameter1;
                uint64_t population = event.parameter2;
                uint64_t newHaplotype = event.parameter3;

                int64_t countBranches = liveBranches[population][newHaplotype].size();
                double probability = static_cast<double>(countBranches) / infectious[population][newHaplotype];
                if (generator_->GetUniform() < probability) {
                    int64_t nodeIndex = static_cast<int64_t>(countBranches * generator_->GetUniform());
                    int64_t node = liveBranches[population][newHaplotype][nodeIndex];
                    liveBranches[population][newHaplotype][nodeIndex] = liveBranches[population][newHaplotype][countBranches - 1];
                    liveBranches[population][newHaplotype].pop_back();
                    liveBranches[population][oldHaplotype].push_back(node);
                    // addMutation({node, oldHaplotype, newHaplotype, time});
                }
                infectious[population][newHaplotype] -= 1;
                infectious[population][oldHaplotype] += 1;
                break;
            }
            // case kMIGRATION: {
            //     int64_t lbs = liveBranches[event.parameter4][event.parameter1].size();
            //     double probability = static_cast<double>(lbs) / infectious[event.parameter4][event.parameter1];
            //     if (generator_->GetUniform() < probability) {
            //         int64_t nt = static_cast<int64_t>(lbs * generator_->GetUniform());
            //         int64_t lbss = liveBranches[event.parameter2][event.parameter1].size();
            //         double probability1 = static_cast<double>(lbss) / infectious[event.parameter2][event.parameter1];
            //         if (generator_->GetUniform() < probability1) {
            //             int64_t ns = static_cast<int64_t>(lbss * generator_->GetUniform());
            //             int64_t idt = liveBranches[event.parameter4][event.parameter1][nt];
            //             int64_t ids = liveBranches[event.parameter2][event.parameter1][ns];
            //             int64_t id3 = currentNode;
            //             liveBranches[event.parameter2][event.parameter1][ns] = id3;
            //             liveBranches[event.parameter4][event.parameter1][nt] =
            //                 liveBranches[event.parameter4][event.parameter1][lbs - 1];
            //             liveBranches[event.parameter4][event.parameter1].pop_back();
                        
            //             addNode(id3, event.parameter2, time, idt, ids);
            //             ++currentNode;
            //             // self.mig.AddMigration(idt, e_time, event.parameter2, event.parameter4)
            //         } else {
            //             liveBranches[event.parameter2][event.parameter1].push_back(liveBranches[event.parameter4][event.parameter1][nt]);
            //             liveBranches[event.parameter4][event.parameter1][nt] =
            //                 liveBranches[event.parameter4][event.parameter1][lbs - 1];
            //             liveBranches[event.parameter4][event.parameter1].pop_back();
            //         }
            //     }

            //     infectious[event.parameter4][event.parameter1] -= 1;
            //     break;
            // }
            // case kSUSCCHANGEN:
            //     break;
            // case kMULTITYPE:
            //     break;
        }
    }
    Debug();
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
    if (left != -1 && right != -1) {
        tree_[left] = parent;
        tree_[right] = parent;
    }
    tree_.push_back(-1);
    population_.push_back(population);
    time_.push_back(time);
}

    // void printTree() {
    //     for (uint64_t index = 0; index < tree_.size(); ++index) {
    //         std::cout << tree_[index] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // void addMutation() {
    //     mutation.push_back({node, });
    // }

    // Mutation getMutation() {
        
    // }

    // void addMigration(int64_t node, int64_t oldPopulation, int64_t newPopulation, double time) {
    //     migration.insert({node, {oldPopulation, newPopulation, time}});
    // }

    // Migration getMigration() {
    //     return migration[node];
    // }

    // bool isMigration(int64_t node) {
    //     return migration.contains(node);
    // }

    // void addRecombination(int64_t child, int64_t parent, int64_t position, int64_t left_haplotype, int64_t right_haplotype) {
    //     recombinations.insert({child, {parent, position, left_haplotype, right_haplotype}});
    // }

    // Recombination getRecombination(int64_t node) {
    //     return recombination[node];
    // }

    // bool isRecombination(int64_t node) {
    //     return recombination.contains(node);
    // }
