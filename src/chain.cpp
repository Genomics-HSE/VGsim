#pragma once

#include "chain.h"

Chain::Chain(Numbers numbers)
    : size_(0)
    , pointer_(0)
    , pointer_multievents_(0)
    , current_time_(0.0)
    , numbers_(numbers)
    , times_(new double[size_])
    , events_(new Event[size_])
    , multievents_(new Multievent[size_]) {
}

Chain::~Chain() {
    delete[] times_;
    delete[] events_;
    delete[] multievents_;
}

void Chain::Reserve(uint64_t add_size) {
    int64_t new_memory = add_size - size_ + pointer_;
    if (new_memory <= 0) {
        return;
    }
    double* new_times = new double[size_ + static_cast<uint64_t>(new_memory)];
    Event* new_events = new Event[size_ + static_cast<uint64_t>(new_memory)];
    Multievent* new_multievents = new Multievent[(size_ + static_cast<uint64_t>(new_memory)) * calculateNumberEvents()];
    for (uint64_t i = 0; i < size_; ++i) {
        new_times[i] = times_[i];
        new_events[i] = events_[i];
        new_multievents[i] = multievents_[i];
    }
    delete[] times_;
    delete[] events_;
    delete[] multievents_;
    size_ += static_cast<uint64_t>(new_memory);
    times_ = new_times;
    events_ = new_events;
    multievents_ = new_multievents;
}

void Chain::Restart() {
    pointer_ = 0;
    size_ = 0;
    current_time_ = 0.0;
    delete[] times_;
    delete[] events_;
    delete[] multievents_;
    times_ = new double[size_];
    events_ = new Event[size_];
    multievents_ = new Multievent[size_ * calculateNumberEvents()];
}

void Chain::AddTime(double time) {
    current_time_ += time;
}

void Chain::AddEvent(Event event) {
    times_[pointer_] = current_time_;
    events_[pointer_++] = event;
}

void Chain::AddMultievent(Multievent event) {
    multievents_[pointer_multievents_++] = event;
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
        std::cout << static_cast<int64_t>(events_[i].type) << ", "
                  << events_[i].parameter1 << ", "
                  << events_[i].parameter2 << ", "
                  << events_[i].parameter3 << ", "
                  << events_[i].parameter4 << ", "
                  << times_[i] << std::endl;
        // DebugEvent(times_[i], events_[i]);
    }
    std::cout << std::endl;
}

void Chain::LastTime() {
    std::cout << "Epidemic time: " << times_[pointer_ - 1] << std::endl;
}

uint64_t Chain::GetSize() {
    return pointer_;
}

uint64_t Chain::GetLastMultievent() {
    return pointer_multievents_;
}

double Chain::GetTime(uint64_t index) {
    return times_[index];
}

Event Chain::GetEvent(uint64_t index) {
    return events_[index];
}

Multievent Chain::GetMultievent(uint64_t index) {
    return multievents_[index];
}

uint64_t Chain::calculateNumberEvents() {
    uint64_t transmission = numbers_.populations * numbers_.haplotypes * numbers_.susceptible_groups;
    uint64_t recovery = numbers_.populations * numbers_.haplotypes;
    uint64_t sampling = numbers_.populations * numbers_.haplotypes;
    uint64_t mutation = numbers_.populations * numbers_.haplotypes * numbers_.sites * 3;
    uint64_t migration = numbers_.populations * (numbers_.populations - 1) * numbers_.haplotypes * numbers_.susceptible_groups;
    uint64_t suscchange = numbers_.populations * numbers_.susceptible_groups * (numbers_.susceptible_groups - 1);
    return transmission + recovery + sampling + mutation + migration + suscchange;
}

void DebugEvent(double time, Event& event) {
    switch (event.type) {
        case TypeEvents::kTRANSMISSION: DebugTransmission(time, event); break;
        case TypeEvents::kRECOVERY: DebugRecovery(time, event); break;
        case TypeEvents::kSAMPLING: DebugSampling(time, event); break;
        case TypeEvents::kMUTATION: DebugMutation(time, event); break;
        case TypeEvents::kMIGRATION: DebugMigration(time, event); break;
        case TypeEvents::kSUSCCHANGE: DebugSuscchange(time, event); break;
        case TypeEvents::kMULTITYPE: DebugMultitype(time, event); break;
    }
}

void DebugTransmission(double time, Event& event) {
    std::cout << "Transmission, haplotype - " << event.parameter1
              << ", population - " << event.parameter2
              << ", group - " << event.parameter3
              << ", second haplotype - " << event.parameter4
              << ", time - " << time
              << "\n";
}

void DebugRecovery(double time, Event& event) {
    std::cout << "Recovery, haplotype - " << event.parameter1
              << ", population - " << event.parameter2
              << ", group - " << event.parameter3
              << ", time - " << time
              << "\n";
}

void DebugSampling(double time, Event& event) {
    std::cout << "Sampling, haplotype - " << event.parameter1
              << ", population - " << event.parameter2
              << ", group - " << event.parameter3
              << ", time - " << time
              << "\n";
}

void DebugMutation(double time, Event& event) {
    std::cout << "Mutation, source haplotype - " << event.parameter1
              << ", population - " << event.parameter2
              << ", target haplotype - " << event.parameter3
              << ", time - " << time
              << "\n";
}

void DebugMigration(double time, Event& event) {
    std::cout << "Migration, haplotype - " << event.parameter1
              << ", source population - " << event.parameter2
              << ", group - " << event.parameter3
              << ", target population - " << event.parameter4
              << ", time - " << time
              << "\n";
}

void DebugSuscchange(double time, Event& event) {
    std::cout << "Suscchange, source group - " << event.parameter1
              << ", population - " << event.parameter2
              << ", target group - " << event.parameter3
              << ", time - " << time
              << "\n";
}

void DebugMultitype(double time, Event& event) {
    std::cout << "Multitype, start - " << event.parameter1
              << ", end - " << event.parameter2
              << ", time - " << time
              << "\n";
}
