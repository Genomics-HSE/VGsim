#pragma once

class ConditionStop {
public:
    ConditionStop(PopulationPool* pool);

    void SetAttempts(uint64_t attempts);
    void SetIterations(uint64_t iterations);
    void Restart();
    uint64_t GetIterations();

    bool CheckAttempt();
    bool CheckIteration();
    bool CheckRestart();

private:
    uint64_t attempts_;
    uint64_t iterations_;
    uint64_t current_attempts_;
    uint64_t current_iterations_;

    PopulationPool* pool_;
};
