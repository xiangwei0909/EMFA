//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "VectorR4.h"

namespace component {

class VectorR3 {

public:
    VectorR3() :x(0.0f), y(0.0f), z(0.0f) {}
    VectorR3(value_t _x, value_t _y, value_t _z)
        :x(_x), y(_y), z(_z)
    {
    }
    VectorR3(const VectorR3&) = default;
    VectorR3(VectorR3&&) = default;
    ~VectorR3() {}

public:
    VectorR3& operator=(const VectorR3&);
    VectorR3& operator+=(const VectorR3& v)
    {
        x += v.x; y += v.y; z += v.z;
        return *this;
    }
    VectorR3& operator-=(const VectorR3& v)
    {
        x -= v.x; y -= v.y; z -= v.z;
        return *this;
    }

    VectorR3& operator*=(const VectorR3& v)
    {
        value_t tx = x, ty = y;
        x = y*v.z - z*v.y;
        y = z*v.x - tx*v.z;
        z = tx*v.y - ty*v.x;
        return *this;
    }
    VectorR3& operator*=(const value_t d)
    {
        x *= d; y *= d; z *= d;
        return *this;
    }
    VectorR3& operator/=(const value_t d)
    {
        assert(d != 0);
        register value_t mInv = value_t(1.0) / d; x *= mInv; y *= mInv; z *= mInv;
        return *this;
    }
    VectorR3& VectorR3::SetFromHg(const VectorR4& v)
    {
        value_t wInv = value_t(1.0) / v.w;
        assert(wInv != 0);
        x = v.x*wInv;
        y = v.y*wInv;
        z = v.z*wInv;
        return *this;
    }
    value_t Norm() const
    {
        return (std::sqrt(x*x + y*y + z*z));
    }
    VectorR3& Normalize()
    {
        *this /= Norm();
        return *this;
    }
public:
    value_t x, y, z;
};
inline VectorR3& VectorR3::operator=(const VectorR3& v1)
{
    if (&v1 != this)
    {
        this->x = v1.x;
        this->y = v1.y;
        this->z = v1.z;
    }
    return *this;
}

inline VectorR3 operator+(const VectorR3& v1, const VectorR3& v2);
inline VectorR3 operator-(const VectorR3& v1, const VectorR3& v2);
inline value_t  operator^(const VectorR3& v1, const VectorR3& v2);
inline VectorR3 operator*(const VectorR3& v1, const VectorR3& v2);
inline VectorR3 operator*(const VectorR3& v1, const value_t d);
inline VectorR3 operator*(const value_t d, const VectorR3& v2);
inline VectorR3 operator/(const VectorR3& v1, const value_t d);
inline value_t Area(const VectorR3& v1, const VectorR3& v2, const VectorR3& v3);
inline value_t Volume(const VectorR3& v1, const VectorR3& v2, const VectorR3& v3, const VectorR3& v4);
inline VectorR3 Normal(const VectorR3& v1, const VectorR3& v2, const VectorR3& v3);
inline VectorR3 operator+(const VectorR3& v1, const VectorR3& v2)
{
    auto tmp = v1;
    return tmp += v2;
}

inline VectorR3 operator-(const VectorR3& v1, const VectorR3& v2)
{
    auto tmp = v1;
    return tmp -= v2;
}

inline value_t operator^(const VectorR3& v1, const VectorR3& v2)
{
    return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

inline VectorR3 operator*(const VectorR3& v1, const VectorR3& v2)
{
    auto tmp = v1;
    return tmp *= v2;
}

inline VectorR3 operator*(const VectorR3& v1, const value_t d)
{
    auto tmp = v1;
    return tmp *= d;
}

inline VectorR3 operator*(const value_t d, const VectorR3& v1)
{
    return v1 * d;
}

inline VectorR3 operator/(const VectorR3& v1, const value_t d)
{
    auto tmp = v1;
    return tmp /= d;
}

inline bool operator==(const VectorR3& v1, const VectorR3& v2)
{
    return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
}

inline bool operator!=(const VectorR3& v1, const VectorR3& v2)
{
    return !(v1 == v2);
}

inline value_t Area(const VectorR3& v1, const VectorR3& v2, const VectorR3& v3)
{
    value_t la = (v1 - v2).Norm();
    value_t lb = (v2 - v3).Norm();
    value_t lc = (v3 - v1).Norm();
    value_t ls = (la + lb + lc) / 2.0f;
    return sqrt(ls*(ls - la)*(ls - lb)*(ls - lc));
}

inline value_t Volume(const VectorR3& v1, const VectorR3& v2, const VectorR3& v3, const VectorR3& v4)
{
	VectorR3 Edge12 = v2 - v1;
	VectorR3 Edge13 = v3 - v1;
	VectorR3 Edge14 = v4 - v1;

	VectorR3 Buff = Edge12*Edge13;

	return Buff^Edge14 / 6.0f;
}

inline VectorR3 Normal(const VectorR3& v1, const VectorR3& v2, const VectorR3& v3)
{
	VectorR3 Edge12 = v2 - v1;
	VectorR3 Edge23 = v3 - v2;

	VectorR3 Buff = (Edge12*Edge23).Normalize();

	return Buff;
}
} // namespace component