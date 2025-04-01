#pragma once

#include <vector>
#include <unordered_map>
#include <iostream>

class ARG {
public:
    ARG(Numbers numbers, Chain* chain, Counters* counters, PopulationPool* pool, RandomGenerator* generator);

    void CalculateGenealogy();

    void Debug();

    std::vector<double> get_time();
    std::vector<int64_t> get_tree();

    // void printTree() {
    //     for (uint64_t index = 0; index < tree_.size(); ++index) {
    //         std::cout << tree_[index] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // Mutation getMutation() {
        
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
    void restart();
    inline void removeBranch(std::vector<std::vector<std::vector<uint64_t>>>& branches, uint64_t population, uint64_t haplotype, uint64_t branch);
    inline void addLeaf(int64_t population, double time);
    inline void addNode(int64_t population, double time, int64_t parent, int64_t left, int64_t right);
    inline void addMutation(uint64_t node, uint64_t oldHaplotype, uint64_t newHaplotype, double time);
    inline void addMigration(uint64_t node, Migration migration);

    inline uint64_t getNumberSites() const;
    inline uint64_t getNumberHaplotypes() const;
    inline uint64_t getNumberPopulations() const;
    inline uint64_t getNumberSusceptibleGroups() const;
    inline uint64_t getIndexHap(uint64_t first, uint64_t second) const;

    std::vector<double> time_;
    std::vector<int64_t> tree_;
    std::vector<int64_t> population_;
    std::vector<Mutation> mutations_;
    std::unordered_map<uint64_t, Migration> migrations_;
    // std::unordered_map<int64_t, Recombination> recombination_;

    Numbers numbers_;
    Chain* chain_;
    Counters* counters_;
    PopulationPool* pool_;
    RandomGenerator* generator_;
};
