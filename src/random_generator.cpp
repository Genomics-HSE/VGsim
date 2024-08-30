#pragma once

#include "random_generator.h"

RandomGenerator::RandomGenerator(uint64_t seed)
    : seed_(seed)
    , generator_(std::mt19937(seed_))
    , uniform_(std::uniform_real_distribution<>(0.0, 1.0)) {
}

// void RandomGenerator::SetSeed(uint64_t new_seed) {
//     seed_ = new_seed;
//     generator_ = std::mt19937(seed_);
// }

void RandomGenerator::NextSeed() {
    ++seed_;
    generator_ = std::mt19937(seed_);
}

double RandomGenerator::GetUniform() {
    return uniform_(generator_);
}

uint64_t RandomGenerator::GetPoisson(double lambda) {
    std::poisson_distribution<> poisson(lambda);
    return poisson(generator_);
}

uint64_t RandomGenerator::GetHypergeometric(uint64_t n, uint64_t K, uint64_t N) {
    uint64_t k = 0;
    double probability = GetUniform();
    double current_probability = 1.0;

    for (uint64_t i = 0; i < n; ++i) {
        current_probability *= (N - K - i) / (N - i);
    }

    while (current_probability < probability) {
        current_probability *= (K - k) * (n - k) / (k + 1) / (N - K - n + k + 1);
        ++k;
    }

    return k;
}
