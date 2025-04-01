#pragma once

class Population {
public:
    Population();
    ~Population();


    void SetParameters(uint64_t number_of_haplotypes, uint64_t number_of_susceptible_groups);
    inline void NewInfections(uint64_t count, uint64_t haplotype, uint64_t group);
    inline void NewRecoveries(uint64_t count, uint64_t haplotype, uint64_t group);
    inline void NewMutation(uint64_t count, uint64_t old_haplotype, uint64_t new_haplotype);
    inline void NewImmunity(uint64_t count, uint64_t old_group, uint64_t new_group);
    void Save();
    void Reset();

    // Population
    void SetSize(uint64_t size);
    inline uint64_t GetSize() const;
    inline uint64_t GetSusceptibles(uint64_t group) const;
    inline uint64_t* GetSusceptiblesBegin() const;
    inline uint64_t GetInfected(uint64_t haplotype) const;
    inline uint64_t* GetInfectedBegin() const;

    // Lockdown
    inline void Switch();
    inline double GetContactDensity() const;
    inline void SetOn(bool on);
    inline bool GetOn() const;
    inline void SetContactDensityBefore(double contact_density_before);
    inline double GetContactDensityBefore() const;
    inline void SetContactDensityAfter(double contact_density_after);
    inline double GetContactDensityAfter() const;
    inline void SetStart(double start);
    inline double GetStart() const;
    inline void SetEnd(double end);
    inline double GetEnd() const;

    // Utility
    uint64_t GetStartNumberSusceptible(uint64_t group) const;
    uint64_t GetStartNumberInfected(uint64_t haplotype) const;

private:
    inline uint64_t getNumberHaplotypes() const;
    inline uint64_t getNumberSusceptibleGroups() const;

    uint64_t size_;
    uint64_t number_of_haplotypes_;
    uint64_t number_of_susceptible_groups_;
    uint64_t* susceptibles_;
    uint64_t* infected_;
    uint64_t* susceptibles_save_;
    uint64_t* infected_save_;
    Lockdown lockdown_;
};
