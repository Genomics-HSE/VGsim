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

    std::vector<uint64_t> hapToNum(numbers_.haplotypes);
    for (uint64_t i = 0; i < numbers_.haplotypes; ++i) {
        hapToNum[i] = i;
    }

    uint64_t size = chain_->GetSize();

    for (uint64_t index = 0; index < size; ++index) {
        for (double time : time_) {
            std::cout << time << " ";
        }
        std::cout << std::endl;

        for (uint64_t node : tree_) {
            std::cout << node << " ";
        }
        std::cout << std::endl;


        // std::cout << size - index - 1 << std::endl;
        Event event = chain_->GetEvent(size - index - 1);
        double time = chain_->GetTime(size - index - 1);
        std::cout << "Check!" << std::endl;
        switch (event.type) {
            case kTRANSMISSION: {
                int64_t countBranches = liveBranches[event.parameter2][event.parameter1].size();
                int64_t countInfectious = infectious[event.parameter2][hapToNum[event.parameter1]];
                double probability = static_cast<double>(countBranches) * (countBranches - 1) /
                    countInfectious / (countInfectious - 1);
                if (generator_->GetUniform() < probability) {
                    int64_t numberBranchKLeft = static_cast<int64_t>(countBranches * generator_->GetUniform());
                    int64_t numberBranchRight = static_cast<int64_t>((countBranches - 1) * generator_->GetUniform());
                    if (numberBranchRight >= numberBranchKLeft) {
                        numberBranchRight += 1;
                    }
                    int64_t idBranchLeft = liveBranches[event.parameter2][event.parameter1][numberBranchKLeft];
                    int64_t idBranchRight = liveBranches[event.parameter2][event.parameter1][numberBranchRight];
                    int64_t idBranchParent = currentNode;
                    addNode(event.parameter2, time, idBranchParent, idBranchLeft, idBranchRight);
                    ++currentNode;
                    liveBranches[event.parameter2][event.parameter1][numberBranchKLeft] = idBranchParent;
                    liveBranches[event.parameter2][event.parameter1][numberBranchRight] = liveBranches[event.parameter2][event.parameter1][countBranches - 1];
                    liveBranches[event.parameter2][event.parameter1].pop_back();
                }
                infectious[event.parameter2][hapToNum[event.parameter1]] -= 1;
                break;
            }
            case kRECOVERY: {
                infectious[event.parameter2][hapToNum[event.parameter1]] += 1;
                break;
            }
            case kSAMPLING: {
                addNode(event.parameter2, time);
                liveBranches[event.parameter2][event.parameter1].push_back(currentNode);
                ++currentNode;
                infectious[event.parameter2][hapToNum[event.parameter1]] += 1;
                break;
            }
        //     case kMUTATION: {
        //         int64_t countBranches = liveBranches[event.parameter2][event.parameter3].size();
        //         double probability = static_cast<double>(countBranches) /
        //             infectious[event.parameter2][event.parameter3];
        //         if (generator_->GetUniform() < probability) {
        //             int64_t numberBranch = static_cast<int64_t>(countBranches * generator_->GetUniform());
        //             int64_t idBranch = liveBranches[event.parameter2][event.parameter3][numberBranch];
        //             liveBranches[event.parameter2][event.parameter3][numberBranch] =
        //                 liveBranches[event.parameter2][event.parameter3][countBranches - 1];
        //             liveBranches[event.parameter2][event.parameter3].pop_back();
        //             liveBranches[event.parameter2][event.parameter1].push_back(idBranch);
        //             // tree.addMutation({idBranch, event.parameter1, event.parameter3, time});
        //         }
        //         infectious[event.parameter2][event.parameter3] -= 1;
        //         infectious[event.parameter2][event.parameter1] += 1;
        //         break;
        //     }
        //     case kMIGRATION: {
        //         int64_t lbs = liveBranches[event.parameter4][event.parameter1].size();
        //         double probability = static_cast<double>(lbs) / infectious[event.parameter4][event.parameter1];
        //         if (generator_->GetUniform() < probability) {
        //             int64_t nt = static_cast<int64_t>(lbs * generator_->GetUniform());
        //             int64_t lbss = liveBranches[event.parameter2][event.parameter1].size();
        //             double probability1 = static_cast<double>(lbss) / infectious[event.parameter2][event.parameter1];
        //             if (generator_->GetUniform() < probability1) {
        //                 int64_t ns = static_cast<int64_t>(lbss * generator_->GetUniform());
        //                 int64_t idt = liveBranches[event.parameter4][event.parameter1][nt];
        //                 int64_t ids = liveBranches[event.parameter2][event.parameter1][ns];
        //                 int64_t id3 = currentNode;
        //                 liveBranches[event.parameter2][event.parameter1][ns] = id3;
        //                 liveBranches[event.parameter4][event.parameter1][nt] =
        //                     liveBranches[event.parameter4][event.parameter1][lbs - 1];
        //                 liveBranches[event.parameter4][event.parameter1].pop_back();
                        
        //                 addNode(id3, event.parameter2, time, idt, ids);
        //                 ++currentNode;
        //                 // self.mig.AddMigration(idt, e_time, event.parameter2, event.parameter4)
        //             } else {
        //                 liveBranches[event.parameter2][event.parameter1].push_back(liveBranches[event.parameter4][event.parameter1][nt]);
        //                 liveBranches[event.parameter4][event.parameter1][nt] =
        //                     liveBranches[event.parameter4][event.parameter1][lbs - 1];
        //                 liveBranches[event.parameter4][event.parameter1].pop_back();
        //             }
        //         }

        //         infectious[event.parameter4][event.parameter1] -= 1;
        //         break;
        //     }
        //     case kSUSCCHANGEN:
        //         break;
        //     case kMULTITYPE:
        //         break;
        }
        std::cout << "Check - 2!" << std::endl;
    }

    return;
}

void ARG::Debug() {
    std::cout << "Some" << std::endl;
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
    tree_.push_back(0);
    population_.push_back(population);
    time_.push_back(time);

    // tree_[parent] = -1;
    // population_[parent] = population;
    // time_[parent] = time;
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
