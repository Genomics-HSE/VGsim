#pragma once

enum class TypeEvents {kTRANSMISSION = 0,
                       kRECOVERY,
                       kSAMPLING,
                       kMUTATION,
                       kSUSCCHANGE,
                       kMIGRATION,
                       kMULTITYPE
                      };

struct Event {
    TypeEvents type;
    uint64_t parameter1;
    uint64_t parameter2;
    uint64_t parameter3;
    uint64_t parameter4;
    double time = 0.0;
};

struct Multievent : public Event {
    uint64_t count;
};


class Chain {
public:
    Chain(Numbers numbers);
    ~Chain();

    void Reserve(uint64_t add_size);
    void Restart();
    void AddTime(double time);
    void AddEvent(Event event);
    void AddMultievent(Multievent event);
    uint64_t Size();
    void Debug();
    double GetLastTime();
    void LastTime();

    uint64_t GetSize();
    uint64_t GetLastMultievent();
    double GetTime(uint64_t index);
    Event GetEvent(uint64_t index);
    Multievent GetMultievent(uint64_t index);

    // Utility
    boost::python::list get_data_susceptible(uint64_t population, uint64_t group, uint64_t step_number, uint64_t start_amount);
    boost::python::list get_data_infected(uint64_t population, uint64_t haplotype, uint64_t step_number, uint64_t start_amount);
    boost::python::list get_data_sample(uint64_t population, uint64_t haplotype, uint64_t step_number);
    boost::python::list get_time_points(uint64_t step_number);

private:
    uint64_t calculateNumberEvents();

    uint64_t size_;
    uint64_t pointer_;
    uint64_t pointer_multievents_;
    double current_time_;
    
    Numbers numbers_;

    Event* events_;
    Multievent* multievents_;
};

void DebugEvent(Event& event);
void DebugTransmission(Event& event);
void DebugRecovery(Event& event);
void DebugSampling(Event& event);
void DebugMutation(Event& event);
void DebugMigration(Event& event);
void DebugSuscchange(Event& event);
void DebugMultitype(Event& event);
