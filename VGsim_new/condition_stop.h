#pragma once

class ConditionStop {
public:
    ConditionStop(Counters* counters, PopulationPool* pool, Chain* chain);

    void SetAttempts(uint64_t attempts);
    void SetSampling(uint64_t sampling);
    void SetIterations(uint64_t iterations);
    void SetEpidemicTime(double epidemic_time);
    void Restart();
    uint64_t GetIterations();

    bool CheckAttempt();
    bool CheckIteration();
    bool CheckRestart();

private:
    uint64_t attempts_;
    uint64_t sampling_;
    uint64_t iterations_;
    double epidemic_time_;
    uint64_t current_attempts_;
    uint64_t current_iterations_;

    Counters* counters_;
    PopulationPool* pool_;
    Chain* chain_;
};
