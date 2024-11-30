// #pragma once

// #include "random_generator.h"

// RandomGenerator::RandomGenerator(uint64_t seed)
//     : seed_(seed)
//     , generator_(std::mt19937(seed_))
//     , uniform_(std::uniform_real_distribution<>(0.0, 1.0)) {
// }

// // void RandomGenerator::SetSeed(uint64_t new_seed) {
// //     seed_ = new_seed;
// //     generator_ = std::mt19937(seed_);
// // }

// void RandomGenerator::NextSeed() {
//     ++seed_;
//     generator_ = std::mt19937(seed_);
// }

// double RandomGenerator::GetUniform() {
//     return uniform_(generator_);
// }

// uint64_t RandomGenerator::GetPoisson(double lambda) {
//     std::poisson_distribution<> poisson(lambda);
//     return poisson(generator_);
// }

// uint64_t RandomGenerator::GetHypergeometric(uint64_t n, uint64_t K, uint64_t N) {
//     uint64_t k = 0;
//     double probability = GetUniform();
//     double current_probability = 1.0;

//     for (uint64_t i = 0; i < n; ++i) {
//         current_probability *= (N - K - i) / (N - i);
//     }

//     while (current_probability < probability) {
//         current_probability *= (K - k) * (n - k) / (k + 1) / (N - K - n + k + 1);
//         ++k;
//     }

//     return k;
// }

#pragma once

#include "random_generator.h"
#include <limits>

RandomGenerator::RandomGenerator(uint32_t seed) {
    seed_ = seed;
    r[0] = seed_;
    for (int i = 1; i < 31; i++) {
        r[i] = (uint32_t)((16807 * (uint64_t)r[i - 1]) % 2147483647);
    }
    for (int i = 31; i < 34; i++) {
        r[i] = r[i - 31];
    }
    for (int i = 34; i < 344; i++) {
        r[i] = r[i - 31] + r[i - 3];
    }

    n = 0;
}

double RandomGenerator::GetUniform() {
    uint32_t x = r[n % 344] = r[(n + 313) % 344] + r[(n + 341) % 344];
    n = (n + 1) % 344;
    // std::cout << (double)(x >> 1) / std::numeric_limits<int32_t>::max() << std::endl;
    return (double)(x >> 1) / std::numeric_limits<int32_t>::max();
}

void RandomGenerator::NextSeed() {
    return;
}

uint64_t RandomGenerator::GetPoisson(double lambda) {
    return 0;
}

uint64_t RandomGenerator::GetHypergeometric(uint64_t n, uint64_t K, uint64_t N) {
    return 0;
}
