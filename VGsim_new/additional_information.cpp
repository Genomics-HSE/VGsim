#pragma once

struct Mutation {
    uint64_t node;
    uint64_t site;
    uint64_t derivedState;
    uint64_t ancestralState;
    double time;
};

struct Migration {
    uint64_t oldPopulation;
    uint64_t newPopulation;
    double time;
};

struct Recombination {
    uint64_t parent;
    uint64_t position;
    uint64_t left_haplotype;
    uint64_t right_haplotype;
};
