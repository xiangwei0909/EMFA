//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "Custom.h"
#include "VectorR3.h"

namespace component {

class Mesh;

class FMM_Boundary {
public:
    FMM_Boundary();
    ~FMM_Boundary();

    void        init(std::shared_ptr<Mesh> pmesh, value_t lambda, value_t delta);

    int         inWhichBox(const VectorR3& mid) const;
    int         boxNum() const;
    void        getXYZBox(int& x, int& y, int& z) const;
    VectorR3    getCenter(int idx) const;
    value_t     getDelta() const;
    void        reportInfo(Qostream& strm) const;
private:
    int         x_box_, y_box_, z_box_;
    value_t     delta_;
    VectorR3    min_pos_, max_pos_;
};

inline int FMM_Boundary::boxNum() const
{
    return x_box_ * y_box_ * z_box_;
}

inline void FMM_Boundary::getXYZBox(int & x, int & y, int & z) const
{
    x = x_box_;
    y = y_box_;
    z = z_box_;
}

inline value_t FMM_Boundary::getDelta() const
{
    return delta_;
}

}   // namespace component