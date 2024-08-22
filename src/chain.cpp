#pragma once

#include "chain.h"

Chain::Chain() 
    : pointer_(0)
    , size_(0)
    , current_time_(0.0)
    , times_(new double[size_])
    , events_(new Event[size_]) {}

Chain::~Chain() {
    delete[] events_;
}

void Chain::Reserve(uint64_t add_size) {
    int64_t new_memory = add_size - size_ + pointer_;
    if (new_memory <= 0) {
        return;
    }
    double* new_times = new double[size_ + static_cast<uint64_t>(new_memory)];
    Event* new_events = new Event[size_ + static_cast<uint64_t>(new_memory)];
    for (uint64_t i = 0; i < size_; ++i) {
        new_times[i] = times_[i];
        new_events[i] = events_[i];
    }
    delete[] times_;
    delete[] events_;
    size_ += static_cast<uint64_t>(new_memory);
    times_ = new_times;
    events_ = new_events;
}

void Chain::Reset() {
    pointer_ = 0;
    size_ = 0;
    current_time_ = 0.0;
    delete[] times_;
    delete[] events_;
    times_ = new double[size_];
    events_ = new Event[size_];
}

void Chain::AddEvent(Event event) {
    times_[pointer_] = current_time_;
    events_[pointer_++] = event;
}

void Chain::AddTime(double time) {
    current_time_ += time;
}

uint64_t Chain::Size() {
    return size_;
}

void Chain::Debug() {
    std::cout << "Pointer: " << pointer_ << std::endl;
    std::cout << "Size: " << size_ << std::endl;
    for (uint64_t i = 0; i < size_; ++i) {
        switch (events_[i].type) {
            case kTRANSMISSION: DebugTransmission(times_[i], events_[i]); break;
            case kRECOVERY: DebugRecovery(times_[i], events_[i]); break;
            case kSAMPLING: DebugSampling(times_[i], events_[i]); break;
            case kMUTATION: DebugMutation(times_[i], events_[i]); break;
            case kMIGRATION: DebugMigration(times_[i], events_[i]); break;
            case kSUSCCHANGE: DebugSuscchange(times_[i], events_[i]); break;
        }
    }
}

void DebugTransmission(double time, Event& event) {
    std::cout << "Transmission, haplotype - " << event.parameter1
              << ", population - " << event.parameter2
              << ", group - " << event.parameter3
              << ", second haplotype - " << event.parameter4
              << ", time - " << time << "\n";
}

void DebugRecovery(double time, Event& event) {
    std::cout << "Recovery, haplotype - " << event.parameter1
              << ", population - " << event.parameter2
              << ", group - " << event.parameter3
              << ", time - " << time << "\n";
}

void DebugSampling(double time, Event& event) {
    std::cout << "Sampling, haplotype - " << event.parameter1
              << ", population - " << event.parameter2
              << ", group - " << event.parameter3
              << ", time - " << time << "\n";
}

void DebugMutation(double time, Event& event) {
    std::cout << "Mutation, source haplotype - " << event.parameter1
              << ", population - " << event.parameter2
              << ", target haplotype - " << event.parameter3
              << ", time - " << time << "\n";
}

void DebugMigration(double time, Event& event) {
    std::cout << "Migration, haplotype - " << event.parameter1
              << ", source population - " << event.parameter2
              << ", group - " << event.parameter3
              << ", target population - " << event.parameter4
              << ", time - " << time << "\n";
}

void DebugSuscchange(double time, Event& event) {
    std::cout << "Suscchange, source group - " << event.parameter1
              << ", population - " << event.parameter2
              << ", target group - " << event.parameter3
              << ", time - " << time << "\n";
}
