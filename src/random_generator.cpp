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
