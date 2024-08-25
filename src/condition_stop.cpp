#pragma once

#include "condition_stop.h"

ConditionStop::ConditionStop()
    : attempts_(0)
    , iterations_(0)
    , current_attempts_(0)
    , current_iterations_(0) {
}

void ConditionStop::SetAttempts(uint64_t attempts) {
    attempts_ = attempts;
}

void ConditionStop::SetIterations(uint64_t iterations) {
    iterations_ = iterations;
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
    return current_iterations_++ != iterations_;
}

bool ConditionStop::CheckRestart() {
    return iterations_ >= 100 && current_iterations_ < 100 && current_attempts_ != attempts_;
}
