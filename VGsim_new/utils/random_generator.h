// #pragma once

// #include <random>
// #include <vector>

// class RandomGenerator {
// public:
//     RandomGenerator(uint64_t seed);

//     // void SetSeed(uint64_t new_seed);
//     void NextSeed();
//     double GetUniform();
//     uint64_t GetPoisson(double lambda);
//     uint64_t GetHypergeometric(uint64_t n, uint64_t K, uint64_t N);

// private:
//     uint64_t seed_;
//     std::mt19937 generator_;
//     std::uniform_real_distribution<> uniform_;
// };

#pragma once

class RandomGenerator {
public:
    RandomGenerator(uint32_t seed);

    void NextSeed();
    double GetUniform();
    uint64_t GetPoisson(double lambda);
    uint64_t GetHypergeometric(uint64_t n, uint64_t K, uint64_t N);

private:
    uint32_t seed_;
    uint32_t r[344];
    uint32_t n;
};
