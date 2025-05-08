#pragma once

class Susceptibles {
public:
    Susceptibles(uint64_t number_of_susceptible_groups);
    void Debug();
    void Update();

    void set_immunity_transition(double rate, uint64_t source_group, uint64_t target_group);
    boost::python::list get_immunity_transition();

    inline double GetSusceptibilityTransition(uint64_t source, uint64_t target) const;
    inline double* GetSusceptibilityTransitionBegin(uint64_t group) const;
    inline double GetSusceptibilityCumulTransition(uint64_t group) const;
    
private:
    inline uint64_t getNumberSusceptibleGroups() const;
    inline uint64_t getIndexSus(uint64_t first, uint64_t second) const;

    uint64_t number_of_susceptible_groups_;

    ArrayBase<double> susceptibility_cumul_transition_;
    ArrayBase<double> susceptibility_transition_;
};
