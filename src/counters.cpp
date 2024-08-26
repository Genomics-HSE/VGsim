#pragma once

#include <iostream>

#include "counters.h"

Counters::Counters() 
    : immunyty_(0)
    , transmission_(0)
    , recovery_(0)
    , sampling_(0)
    , mutation_(0)
    , migration_accept_(0)
    , migration_reject_(0) {}

void Counters::Debug() {
    std::cout << "COUNTERS\n";
    std::cout << "----------------------------------------\n";
    std::cout << "Immunity counter: " << immunyty_ << "\n";
    std::cout << "Transmission counter: " << transmission_ << "\n";
    std::cout << "Recovery counter: " << recovery_ << "\n";
    std::cout << "Sampling counter: " << sampling_ << "\n";
    std::cout << "Mutation counter: " << mutation_ << "\n";
    std::cout << "Migration accept counter: " << migration_accept_ << "\n";
    std::cout << "Migration reject counter: " << migration_reject_ << "\n";
    std::cout << std::endl;
}

void Counters::Restart() {
    immunyty_ = 0;
    transmission_ = 0;
    recovery_ = 0;
    sampling_ = 0;
    mutation_ = 0;
    migration_accept_ = 0;
    migration_reject_ = 0;
}

inline void Counters::AddImmunity(uint64_t count) {
    immunyty_ += count;
}

inline void Counters::AddTransmission(uint64_t count) {
    transmission_ += count;
}

inline void Counters::AddRecovery(uint64_t count) {
    recovery_ += count;
}

inline void Counters::AddSampling(uint64_t count) {
    sampling_ += count;
}

uint64_t Counters::GetSampling() {
    return sampling_;
}

inline void Counters::AddMutation(uint64_t count) {
    mutation_ += count;
}

inline void Counters::AddMigrationAccept(uint64_t count) {
    migration_accept_ += count;
}

inline void Counters::AddMigrationReject(uint64_t count) {
    migration_reject_ += count;
}
