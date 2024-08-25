#pragma once

constexpr uint64_t kTRANSMISSION = 0;
constexpr uint64_t kRECOVERY = 1;
constexpr uint64_t kSAMPLING = 2;
constexpr uint64_t kMUTATION = 3;
constexpr uint64_t kMIGRATION = 4;
constexpr uint64_t kSUSCCHANGE = 5;
constexpr uint64_t kMULTITYPE = 6;

struct Event {
    uint64_t type;
    uint64_t parameter1;
    uint64_t parameter2;
    uint64_t parameter3;
    uint64_t parameter4;
};

class Chain {
public:
    Chain();
    ~Chain();

    void Reserve(uint64_t add_size);
    void Restart();
    void AddTime(double time);
    void AddEvent(Event event);
    uint64_t Size();
    void Debug();

private:
    uint64_t pointer_;
    uint64_t size_;
    double current_time_;
    double* times_;
    Event* events_;
};

void DebugTransmission(double time, Event& event);
void DebugRecovery(double time, Event& event);
void DebugSampling(double time, Event& event);
void DebugMutation(double time, Event& event);
void DebugMigration(double time, Event& event);
void DebugSuscchange(double time, Event& event);
