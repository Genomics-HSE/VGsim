#pragma once

class Lockdown {
public:
    Lockdown();
    Lockdown(bool on, double contact_density_before, double contact_density_after, double start, double end);

    void Switch();
    inline double GetCurrentContactDensity() const;

    void SetOn(bool on);
    inline bool GetOn() const;
    void SetContactDensityBefore(double contact_density_before);
    inline double GetContactDensityBefore() const;
    void SetContactDensityAfter(double contact_density_after);
    inline double GetContactDensityAfter() const;
    void SetStart(double start);
    inline double GetStart() const;
    void SetEnd(double end);
    inline double GetEnd() const;

private:
    bool on_;
    double contact_density_before_;
    double contact_density_after_;
    double start_;
    double end_;
};
