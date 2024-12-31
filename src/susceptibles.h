#pragma once

class Susceptibles {
public:
    Susceptibles(uint64_t number_of_susceptible_groups);
    ~Susceptibles();
    void Debug();
    void Update();

    void set_immunity_transition(double rate, uint64_t source_group, uint64_t target_group);
    PyObject* get_immunity_transition();

    void SetSusceptibilityTransition(double rate, uint64_t source, uint64_t target);
    inline double GetSusceptibilityTransition(uint64_t source, uint64_t target) const;
    inline double* GetSusceptibilityTransitionBegin(uint64_t group) const;
    inline double GetSusceptibilityCumulTransition(uint64_t group) const;
    inline double* GetSusceptibilityCumulTransitionBegin() const;
    
private:
    inline uint64_t getNumberSusceptibleGroups() const;
    inline uint64_t getIndexSus(uint64_t first, uint64_t second) const;

    uint64_t number_of_susceptible_groups_;

    double* susceptibility_cumul_transition_;
    double* susceptibility_transition_;
};
