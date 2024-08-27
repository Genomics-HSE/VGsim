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

void Chain::Restart() {
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
    std::cout << "CHAIN" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Pointer: " << pointer_ << std::endl;
    std::cout << "Size: " << size_ << std::endl;
    for (uint64_t i = 0; i < pointer_; ++i) {
        DebugEvent(times_[i], events_[i]);
    }
    std::cout << std::endl;
}

uint64_t Chain::GetSize() {
    return pointer_;
}

Event Chain::GetEvent(uint64_t index) {
    return events_[index];
}

double Chain::GetTime(uint64_t index) {
    return times_[index];
}

void DebugEvent(double time, Event& event) {
    switch (event.type) {
        case kTRANSMISSION: DebugTransmission(time, event); break;
        case kRECOVERY: DebugRecovery(time, event); break;
        case kSAMPLING: DebugSampling(time, event); break;
        case kMUTATION: DebugMutation(time, event); break;
        case kMIGRATION: DebugMigration(time, event); break;
        case kSUSCCHANGE: DebugSuscchange(time, event); break;
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
