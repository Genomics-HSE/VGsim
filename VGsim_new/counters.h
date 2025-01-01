#pragma once

class Counters {
public:
    Counters();
    void Debug();
    void Restart();

    inline void AddImmunity(uint64_t count);
    inline void AddTransmission(uint64_t count);
    inline void AddRecovery(uint64_t count);
    inline void AddSampling(uint64_t count);
    uint64_t GetSampling();
    inline void AddMutation(uint64_t count);
    inline void AddMigrationAccept(uint64_t count);
    inline void AddMigrationReject(uint64_t count);

private:
    uint64_t immunyty_;
    uint64_t transmission_;
    uint64_t recovery_;
    uint64_t sampling_;
    uint64_t mutation_;
    uint64_t migration_accept_;
    uint64_t migration_reject_;
};
