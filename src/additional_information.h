#pragma once

struct Mutation {
    int64_t node;
    int64_t site;
    int64_t derivedState;
    int64_t ancestralState;
    double time;
};

struct Migration {
    int64_t oldPopulation;
    int64_t newPopulation;
    double time;
};

struct Recombination {
    int64_t parent;
    int64_t position;
    int64_t left_haplotype;
    int64_t right_haplotype;
};
