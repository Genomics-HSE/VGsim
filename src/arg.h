#pragma once


#include <vector>
#include <unordered_map>
#include "additional_information.h"

class ARG {
public:
    ARG(int64_t size);

    void AddNode(int64_t parent, int64_t population, double time, int64_t left = -1, int64_t right = -1);

    void Debug();

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

private:
    std::vector<double> time_;
    std::vector<int64_t> tree_;
    std::vector<int64_t> population_;
    std::vector<Mutation> mutation_;
    std::unordered_map<int64_t, Migration> migration_;
    std::unordered_map<int64_t, Recombination> recombination_;
};
