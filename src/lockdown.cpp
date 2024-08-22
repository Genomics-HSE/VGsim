#pragma once

#include "lockdown.h"

Lockdown::Lockdown() {
}

Lockdown::Lockdown(bool on, double contact_density_before, double contact_density_after, double start, double end)
    : on_(on)
    , contact_density_before_(contact_density_before)
    , contact_density_after_(contact_density_after)
    , start_(start)
    , end_(end) {
}


void Lockdown::Switch() {
    on_ = !on_;
}

inline double Lockdown::GetCurrentContactDensity() const {
    return on_ ? contact_density_after_ : contact_density_before_;
}


void Lockdown::SetOn(bool on) {
    on_ = on;
}

inline bool Lockdown::GetOn() const {
    return on_;
}

void Lockdown::SetContactDensityBefore(double contact_density_before) {
    contact_density_before_ = contact_density_before;
}

inline double Lockdown::GetContactDensityBefore() const {
    return contact_density_before_;
}

void Lockdown::SetContactDensityAfter(double contact_density_after) {
    contact_density_after_ = contact_density_after;
}

inline double Lockdown::GetContactDensityAfter() const {
    return contact_density_after_;
}

void Lockdown::SetStart(double start) {
    start_ = start;
}

inline double Lockdown::GetStart() const {
    return start_;
}

void Lockdown::SetEnd(double end) {
    end_ = end;
}

inline double Lockdown::GetEnd() const {
    return end_;
}
