//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include"VectorR3.h"

namespace component {

class Boundary {
public:
    void        setBoundary(const VectorR3& _min, const VectorR3& _max);
    VectorR3    getBoundaryMin() const;
    VectorR3    getBoundaryMax() const;
    VectorR3    getBoundarySize() const;
    bool        inBoundary(const VectorR3& v) const;

private:
    VectorR3 minBox, maxBox;
};

inline void Boundary::setBoundary(const VectorR3& _min, const VectorR3& _max)
{
    minBox = _min;
    maxBox = _max;
}

inline VectorR3 Boundary::getBoundaryMin() const
{
    return minBox;
}

inline VectorR3 Boundary::getBoundaryMax() const
{
    return maxBox;
}

inline VectorR3 Boundary::getBoundarySize() const
{
    return (maxBox - minBox);
}

inline bool Boundary::inBoundary(const VectorR3& v) const
{
    if (v.x<minBox.x || v.y<minBox.y || v.z<minBox.z || v.x>maxBox.x || v.y>maxBox.y || v.z>maxBox.z)
        return false;
    return true;
}

} // namespace component