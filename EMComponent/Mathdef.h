//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "Custom.h"

namespace component {

const value_t MAX_VALUE = std::numeric_limits<value_t>::max();
const value_t MIN_VALUE = std::numeric_limits<value_t>::lowest();

const value_t PI = arma::Datum<value_t>::pi;
const value_t PI2 = 2.0f * PI;
const value_t PI4 = 4.0f * PI;
const value_t PI9 = 9.0f * PI;
const value_t PI16 = 16.0f * PI;
const value_t PISq = PI * PI;
const value_t PIhalves = 0.5f * PI;
const value_t PIthirds = PI / 3.0f;
const value_t PItwothirds = PI2 / 3.0f;
const value_t PIfourths = 0.25f * PI;

const value_t RadiansToDegrees = 180.f / PI;
const value_t DegreesToRadians = PI / 180.0f;

const value_t cc = arma::Datum<value_t>::c_0;
const value_t mu0 = arma::Datum<value_t>::mu_0;
const value_t eps0 = arma::Datum<value_t>::eps_0;

const Complex J0(0, 1);
const value_t Z0 = 120.f * PI;

template<class T>
inline T Min(const T& x, const T& y)
{
    return (x < y ? x : y);
}

template<class T>
inline T Max(const T& x, const T& y)
{
    return (y < x ? x : y);
}

} // namespace component