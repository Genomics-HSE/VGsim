#pragma once

#include "chain.cpp"
#include "arg.cpp"

#include <vector>
#include <random>

ARG Genealogy(Chain& chain, std::vector<std::vector<int64_t>> infectious, const std::vector<int64_t>& hapToNum) {
    int64_t samplingCounter = 0;
    for (auto index = chain.begin(); index != chain.end(); ++index) {
        if (index->event_type == kSAMPLING) {
            ++samplingCounter;
        }
    }
    ARG tree{2 * samplingCounter - 1};
    std::vector<std::vector<std::vector<int64_t>>> liveBranches;
    liveBranches.resize(infectious.size());
    for (uint64_t i = 0; i < infectious.size(); ++i) {
        liveBranches[i].resize(infectious[i].size());
    }
    int64_t currentNode = 0;
    // for (uint64_t i = 0; i < infectious.size(); ++i) {
    //     for (uint64_t j = 0; j < infectious[i].size(); ++j) {
    //         std::cout << infectious[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> rn(0.0, 1.0);
    for (auto event = chain.end() - 1; event != chain.begin() - 1; --event) {
        switch (event->event_type) {
            case kTRANSMISSION: {
                int64_t countBranches = liveBranches[event->parameter2][event->parameter1].size();
                int64_t countInfectious = infectious[event->parameter2][hapToNum[event->parameter1]];
                double probability = static_cast<double>(countBranches) * (countBranches - 1.0) /
                    countInfectious / (countInfectious - 1.0);
                if (rn(gen) < probability) {
                    int64_t numberBranchKLeft = static_cast<int64_t>(countBranches * rn(gen));
                    int64_t numberBranchRight = static_cast<int64_t>((countBranches - 1) * rn(gen));
                    if (numberBranchRight >= numberBranchKLeft) {
                        numberBranchRight += 1;
                    }
                    int64_t idBranchLeft = liveBranches[event->parameter2][event->parameter1][numberBranchKLeft];
                    int64_t idBranchRight = liveBranches[event->parameter2][event->parameter1][numberBranchRight];
                    int64_t idBranchParent = currentNode;
                    tree.addNode(idBranchParent, event->parameter2, event->time, idBranchLeft, idBranchRight);
                    ++currentNode;
                    liveBranches[event->parameter2][event->parameter1][numberBranchKLeft] = idBranchParent;
                    liveBranches[event->parameter2][event->parameter1][numberBranchRight] = liveBranches[event->parameter2][event->parameter1][countBranches - 1];
                    liveBranches[event->parameter2][event->parameter1].pop_back();
                }
                infectious[event->parameter2][hapToNum[event->parameter1]] -= 1;
                break;
            }
            case kRECOVERY: {
                infectious[event->parameter2][hapToNum[event->parameter1]] += 1;
                break;
            }
            case kSAMPLING: {
                tree.addNode(currentNode, event->parameter2, event->time);
                liveBranches[event->parameter2][event->parameter1].push_back(currentNode);
                ++currentNode;
                infectious[event->parameter2][hapToNum[event->parameter1]] += 1;
                break;
            }
            case kMUTATION: {
                int64_t countBranches = liveBranches[event->parameter2][event->parameter3].size();
                double probability = static_cast<double>(countBranches) /
                    infectious[event->parameter2][event->parameter3];
                if (rn(gen) < probability) {
                    int64_t numberBranch = static_cast<int64_t>(countBranches * rn(gen));
                    int64_t idBranch = liveBranches[event->parameter2][event->parameter3][numberBranch];
                    liveBranches[event->parameter2][event->parameter3][numberBranch] =
                        liveBranches[event->parameter2][event->parameter3][countBranches - 1];
                    liveBranches[event->parameter2][event->parameter3].pop_back();
                    liveBranches[event->parameter2][event->parameter1].push_back(idBranch);
                    // tree.addMutation({idBranch, event->parameter1, event->parameter3, event->time});
                }
                infectious[event->parameter2][event->parameter3] -= 1;
                infectious[event->parameter2][event->parameter1] += 1;
                break;
            }
            case kMIGRATION: {
                int64_t lbs = liveBranches[event->parameter4][event->parameter1].size();
                double probability = static_cast<double>(lbs) / infectious[event->parameter4][event->parameter1];
                if (rn(gen) < probability) {
                    int64_t nt = static_cast<int64_t>(lbs * rn(gen));
                    int64_t lbss = liveBranches[event->parameter2][event->parameter1].size();
                    double probability1 = static_cast<double>(lbss) / infectious[event->parameter2][event->parameter1];
                    if (rn(gen) < probability1) {
                        int64_t ns = static_cast<int64_t>(lbss * rn(gen));
                        int64_t idt = liveBranches[event->parameter4][event->parameter1][nt];
                        int64_t ids = liveBranches[event->parameter2][event->parameter1][ns];
                        int64_t id3 = currentNode;
                        liveBranches[event->parameter2][event->parameter1][ns] = id3;
                        liveBranches[event->parameter4][event->parameter1][nt] =
                            liveBranches[event->parameter4][event->parameter1][lbs - 1];
                        liveBranches[event->parameter4][event->parameter1].pop_back();
                        
                        tree.addNode(id3, event->parameter2, event->time, idt, ids);
                        ++currentNode;
                        // self.mig.AddMigration(idt, e_time, event->parameter2, event->parameter4)
                    } else {
                        liveBranches[event->parameter2][event->parameter1].push_back(liveBranches[event->parameter4][event->parameter1][nt]);
                        liveBranches[event->parameter4][event->parameter1][nt] =
                            liveBranches[event->parameter4][event->parameter1][lbs - 1];
                        liveBranches[event->parameter4][event->parameter1].pop_back();
                    }
                }

                infectious[event->parameter4][event->parameter1] -= 1;
                break;
            }
            case kSUSCCHANGEN:
                break;
            case kMULTITYPE:
                break;
        }
        // std::cout << event->time << std::endl;
    }
    // std::cout << "Check!" << std::endl;
    // tree.printTree();
    // std::cout << "Check!" << std::endl;
    // for (uint64_t i = 0; i < infectious.size(); ++i) {
    //     for (uint64_t j = 0; j < infectious[i].size(); ++j) {
    //         std::cout << infectious[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // for (uint64_t i = 0; i < liveBranches.size(); ++i) {
    //     for (uint64_t j = 0; j < liveBranches[i].size(); ++j) {
    //         std::cout << "Размер: " << liveBranches[i][j].size() << std::endl;
    //         for (uint64_t k = 0; k < liveBranches[i][j].size(); ++k) {
    //             std::cout << liveBranches[i][j][k] << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << std::endl;
    // }
    return tree;
}
