#pragma once

#include "condition_stop.h"

ConditionStop::ConditionStop(Counters* counters, PopulationPool* pool, Chain* chain)
    : attempts_(0)
    , sampling_(0)
    , iterations_(0)
    , epidemic_time_(0.0)
    , current_attempts_(0)
    , current_iterations_(0)
    , counters_(counters)
    , pool_(pool)
    , chain_(chain) {
}

void ConditionStop::SetAttempts(uint64_t attempts) {
    attempts_ = attempts;
}

void ConditionStop::SetSampling(uint64_t sampling) {
    sampling_ = sampling;
}

void ConditionStop::SetIterations(uint64_t iterations) {
    iterations_ = iterations;
}

void ConditionStop::SetEpidemicTime(double epidemic_time) {
    epidemic_time_ = epidemic_time;
}

void ConditionStop::Restart() {
    current_attempts_ = 0;
    current_iterations_ = 0;
}

uint64_t ConditionStop::GetIterations() {
    return iterations_;
}

bool ConditionStop::CheckAttempt() {
    return current_attempts_++ != attempts_;
}

bool ConditionStop::CheckIteration() {
    return current_iterations_++ != iterations_
            && pool_->GetInfected() != 0
            && counters_->GetSampling() < sampling_
            && (epidemic_time_ == 0.0 || chain_->GetLastTime() < epidemic_time_);
}

bool ConditionStop::CheckRestart() {
    return iterations_ >= 100
            && current_iterations_ < 100
            && current_attempts_ != attempts_;
}

void ConditionStop::RevertIteration() {
    --current_iterations_;
}
