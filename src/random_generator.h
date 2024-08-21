#pragma once

#include <random>
#include <vector>

class RandomGenerator {
public:
    RandomGenerator(uint64_t seed);

    // void SetSeed(uint64_t new_seed);
    void NextSeed();
    double GetUniform();
    uint64_t GetPoisson(double lambda);

private:
    uint64_t seed_;
    std::mt19937 generator_;
    std::uniform_real_distribution<> uniform_;
};
