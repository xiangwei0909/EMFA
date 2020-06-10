//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "Custom.h"

namespace component {

class VectorR4 {
public:
    VectorR4() :x(0.0), y(0.0), z(0.0), w(0.0) {}
    VectorR4(const VectorR4&) = default;
    VectorR4(VectorR4&&) = default;
    ~VectorR4() = default;
public:
    VectorR4& operator=(const VectorR4&) = default;
    VectorR4& operator=(VectorR4&&) = default;
public:
    value_t x, y, z, w;
};

} // namespace component