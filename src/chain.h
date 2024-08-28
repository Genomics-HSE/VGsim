#pragma once

enum class TypeEvents {kTRANSMISSION, kRECOVERY, kSAMPLING, kMUTATION, kMIGRATION, kSUSCCHANGE, kMULTITYPE};

struct Event {
    TypeEvents type;
    uint64_t parameter1;
    uint64_t parameter2;
    uint64_t parameter3;
    uint64_t parameter4;
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

    uint64_t GetSize();
    uint64_t GetLastMultievent();
    double GetTime(uint64_t index);
    Event GetEvent(uint64_t index);
    Multievent GetMultievent(uint64_t index);

private:
    uint64_t calculateNumberEvents();

    uint64_t size_;
    uint64_t pointer_;
    uint64_t pointer_multievents_;
    double current_time_;
    
    Numbers numbers_;

    double* times_;
    Event* events_;
    Multievent* multievents_;
};

void DebugEvent(double time, Event& event);
void DebugTransmission(double time, Event& event);
void DebugRecovery(double time, Event& event);
void DebugSampling(double time, Event& event);
void DebugMutation(double time, Event& event);
void DebugMigration(double time, Event& event);
void DebugSuscchange(double time, Event& event);
void DebugMultitype(double time, Event& event);
