#pragma once

#include "chain.h"

Chain::Chain(Numbers numbers)
    : size_(0)
    , pointer_(0)
    , pointer_multievents_(0)
    , current_time_(0.0)
    , numbers_(numbers)
    , events_(new Event[size_])
    , multievents_(new Multievent[size_]) {
}

Chain::~Chain() {
    delete[] events_;
    delete[] multievents_;
}

void Chain::Reserve(uint64_t add_size) {
    int64_t new_memory = add_size - size_ + pointer_;
    if (new_memory <= 0) {
        return;
    }
    Event* new_events = new Event[size_ + static_cast<uint64_t>(new_memory)];
    for (uint64_t i = 0; i < size_; ++i) {
        new_events[i] = events_[i];
    }
    delete[] events_;
    size_ += static_cast<uint64_t>(new_memory);
    events_ = new_events;
}

void Chain::Restart() {
    pointer_ = 0;
    size_ = 0;
    current_time_ = 0.0;
    delete[] events_;
    events_ = new Event[size_];
}

void Chain::AddTime(double time) {
    current_time_ += time;
}

void Chain::AddEvent(Event event) {
    event.time = current_time_;
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
                  << events_[i].time << std::endl;
        // DebugEvent(times_[i], events_[i]);
    }
    std::cout << std::endl;
}

double Chain::GetLastTime() {
    return current_time_;
}

void Chain::LastTime() {
    std::cout << "Epidemic time: " << current_time_ << std::endl;
}

uint64_t Chain::GetSize() {
    return pointer_;
}

uint64_t Chain::GetLastMultievent() {
    return pointer_multievents_;
}

double Chain::GetTime(uint64_t index) {
    return events_[index].time;
}

Event Chain::GetEvent(uint64_t index) {
    return events_[index];
}

Multievent Chain::GetMultievent(uint64_t index) {
    return multievents_[index];
}

// Utility
boost::python::list Chain::get_data_susceptible(uint64_t population, uint64_t group, uint64_t step_number, uint64_t start_amount) {
    boost::python::list data;

    double part_time = current_time_ / step_number;
    uint64_t counter = 0;
    for (uint64_t i = 0; i <= pointer_; ++i) {
        Event current_event = events_[i];
        while (part_time * counter < current_event.time) {
            ++counter;
            data.append(boost::python::object(start_amount));
        }

        switch (current_event.type) {
            case TypeEvents::kTRANSMISSION:
                if (current_event.parameter2 == population && current_event.parameter3 == group) {
                    --start_amount;
                }
                break;
            case TypeEvents::kRECOVERY:
            case TypeEvents::kSAMPLING:
                if (current_event.parameter2 == population && current_event.parameter3 == group) {
                    ++start_amount;
                }
                break;
            case TypeEvents::kMUTATION:
                break;
            case TypeEvents::kMIGRATION:
                if (current_event.parameter4 == population && current_event.parameter3 == group) {
                    --start_amount;
                }
                break;
            case TypeEvents::kSUSCCHANGE:
                if (current_event.parameter2 == population && current_event.parameter3 == group) {
                    ++start_amount;
                } else if (current_event.parameter2 == population && current_event.parameter1 == group) {
                    --start_amount;
                }
                break;
            case TypeEvents::kMULTITYPE:
                break;
        }
    }
    data.append(boost::python::object(start_amount));

    return data;
}

boost::python::list Chain::get_data_infected(uint64_t population, uint64_t haplotype, uint64_t step_number, uint64_t start_amount) {
    boost::python::list data;

    double part_time = current_time_ / step_number;
    uint64_t counter = 0;
    for (uint64_t i = 0; i <= pointer_; ++i) {
        Event current_event = events_[i];
        while (part_time * counter < current_event.time) {
            ++counter;
            data.append(boost::python::object(start_amount));
        }

        switch (current_event.type) {
            case TypeEvents::kTRANSMISSION:
                if (current_event.parameter2 == population && current_event.parameter1 == haplotype) {
                    ++start_amount;
                }
                break;
            case TypeEvents::kRECOVERY:
            case TypeEvents::kSAMPLING:
                if (current_event.parameter2 == population && current_event.parameter1 == haplotype) {
                    --start_amount;
                }
                break;
            case TypeEvents::kMUTATION:
                if (current_event.parameter2 == population && current_event.parameter3 == haplotype) {
                    ++start_amount;
                } else if (current_event.parameter2 == population && current_event.parameter1 == haplotype) {
                    --start_amount;
                }
                break;
            case TypeEvents::kMIGRATION:
                if (current_event.parameter4 == population && current_event.parameter1 == haplotype) {
                    ++start_amount;
                }
                break;
            case TypeEvents::kSUSCCHANGE:
                break;
            case TypeEvents::kMULTITYPE:
                break;
        }
    }
    data.append(boost::python::object(start_amount));

    return data;
}

boost::python::list Chain::get_data_sample(uint64_t population, uint64_t haplotype, uint64_t step_number) {
    boost::python::list data;
    uint64_t start_amount = 0;

    double part_time = current_time_ / step_number;
    uint64_t counter = 0;
    for (uint64_t i = 0; i <= pointer_; ++i) {
        Event current_event = events_[i];
        while (part_time * counter < current_event.time) {
            ++counter;
            data.append(boost::python::object(start_amount));
        }

        switch (current_event.type) {
            case TypeEvents::kTRANSMISSION:
                break;
            case TypeEvents::kRECOVERY:
                break;
            case TypeEvents::kSAMPLING:
                if (current_event.parameter2 == population && current_event.parameter1 == haplotype) {
                    ++start_amount;
                }
                break;
            case TypeEvents::kMUTATION:
                break;
            case TypeEvents::kMIGRATION:
                break;
            case TypeEvents::kSUSCCHANGE:
                break;
            case TypeEvents::kMULTITYPE:
                break;
        }
    }
    data.append(boost::python::object(start_amount));

    return data;
}

boost::python::list Chain::get_time_points(uint64_t step_number) {
    boost::python::list data;

    double part_time = current_time_ / step_number;
    for (uint64_t step = 0; step <= step_number; ++step) {
        data.append(boost::python::object(step * part_time));
    }

    return data;
}


// Private
uint64_t Chain::calculateNumberEvents() {
    uint64_t transmission = numbers_.populations * numbers_.haplotypes * numbers_.susceptible_groups;
    uint64_t recovery = numbers_.populations * numbers_.haplotypes;
    uint64_t sampling = numbers_.populations * numbers_.haplotypes;
    uint64_t mutation = numbers_.populations * numbers_.haplotypes * numbers_.sites * (numbers_.allele_states - 1);
    uint64_t migration = numbers_.populations * (numbers_.populations - 1) * numbers_.haplotypes * numbers_.susceptible_groups;
    uint64_t suscchange = numbers_.populations * numbers_.susceptible_groups * (numbers_.susceptible_groups - 1);
    return transmission + recovery + sampling + mutation + migration + suscchange;
}

void DebugEvent(Event& event) {
    switch (event.type) {
        case TypeEvents::kTRANSMISSION: DebugTransmission(event); break;
        case TypeEvents::kRECOVERY: DebugRecovery(event); break;
        case TypeEvents::kSAMPLING: DebugSampling(event); break;
        case TypeEvents::kMUTATION: DebugMutation(event); break;
        case TypeEvents::kMIGRATION: DebugMigration(event); break;
        case TypeEvents::kSUSCCHANGE: DebugSuscchange(event); break;
        case TypeEvents::kMULTITYPE: DebugMultitype(event); break;
    }
}

void DebugTransmission(Event& event) {
    std::cout << "Transmission, haplotype - " << event.parameter1
              << ", population - " << event.parameter2
              << ", group - " << event.parameter3
              << ", second haplotype - " << event.parameter4
              << ", time - " << event.time
              << "\n";
}

void DebugRecovery(Event& event) {
    std::cout << "Recovery, haplotype - " << event.parameter1
              << ", population - " << event.parameter2
              << ", group - " << event.parameter3
              << ", time - " << event.time
              << "\n";
}

void DebugSampling(Event& event) {
    std::cout << "Sampling, haplotype - " << event.parameter1
              << ", population - " << event.parameter2
              << ", group - " << event.parameter3
              << ", time - " << event.time
              << "\n";
}

void DebugMutation(Event& event) {
    std::cout << "Mutation, source haplotype - " << event.parameter1
              << ", population - " << event.parameter2
              << ", target haplotype - " << event.parameter3
              << ", time - " << event.time
              << "\n";
}

void DebugMigration(Event& event) {
    std::cout << "Migration, haplotype - " << event.parameter1
              << ", source population - " << event.parameter2
              << ", group - " << event.parameter3
              << ", target population - " << event.parameter4
              << ", time - " << event.time
              << "\n";
}

void DebugSuscchange(Event& event) {
    std::cout << "Suscchange, source group - " << event.parameter1
              << ", population - " << event.parameter2
              << ", target group - " << event.parameter3
              << ", time - " << event.time
              << "\n";
}

void DebugMultitype(Event& event) {
    std::cout << "Multitype, start - " << event.parameter1
              << ", end - " << event.parameter2
              << ", time - " << event.time
              << "\n";
}
