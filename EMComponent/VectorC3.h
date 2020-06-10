//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "VectorR3.h"

namespace component {

class VectorC3 {
public:
    VectorC3() : x(0, 0), y(0, 0), z(0, 0) { }
    VectorC3(value_t xVal, value_t yVal, value_t zVal)
        :x(Complex(xVal, 0)), y(Complex(yVal, 0)), z(Complex(zVal, 0))
    {
    }
    VectorC3(const Complex& xVal, const Complex& yVal, const Complex& zVal)
        : x(xVal), y(yVal), z(zVal)
    {
    }
    VectorC3(const VectorC3&) = default;
    VectorC3(VectorC3&&) = default;
    ~VectorC3() = default;

    VectorC3& operator=(const VectorC3& v) = default;
    VectorC3& operator=(VectorC3&&) = default;
    VectorC3& operator+=(const VectorC3& v);
	VectorC3& operator+=(const VectorR3& v);
    VectorC3& operator-=(const VectorC3& v);
    VectorC3& operator*=(const VectorC3& v);
    VectorC3& operator*=(const value_t d);
    VectorC3& operator/=(const value_t d);
    value_t norm() const
    {
        return std::norm(x) + std::norm(y) + std::norm(z);
    }
public:
    Complex x, y, z;
};
inline VectorC3& VectorC3::operator*=(const VectorC3& v)
{
    Complex tempx = x, tempy = y;
    x = y*v.z - z*v.y;
    y = z*v.x - tempx*v.z;
    z = tempx*v.y - tempy*v.x;
    return *this;
}

inline VectorC3& VectorC3::operator+=(const VectorC3& v)
{
    x += v.x, y += v.y, z += v.z;
    return *this;
}

inline VectorC3& VectorC3::operator+=(const VectorR3& v)
{
	x += v.x, y += v.y, z += v.z;
	return *this;
}

inline VectorC3& VectorC3::operator-=(const VectorC3& v)
{
    x -= v.x, y -= v.y, z -= v.z;
    return *this;
}

inline VectorC3 & VectorC3::operator*=(const value_t d)
{
    x *= d; y *= d; z *= d;
    return *this;
}

inline VectorC3 & VectorC3::operator/=(const value_t d)
{
    assert(d != 0);
    x /= d; y /= d; z /= d;
    return *this;
}

//非成员函数
inline VectorC3 operator+(const VectorC3 &v1, const VectorC3 &v2);
inline VectorC3 operator+(const VectorR3 &v1, const VectorC3 &v2);
inline VectorC3 operator+(const VectorC3 &v1, const VectorR3 &v2);
inline VectorC3 operator-(const VectorC3 &v1, const VectorC3 &v2);
inline VectorC3 operator*(const VectorC3 &v1, const VectorC3 &v2);
inline VectorC3 operator*(const VectorC3 &vc, const VectorR3 &vr);
inline VectorC3 operator*(const VectorR3 &vr, const VectorC3 &vc);
inline VectorC3 operator*(const VectorC3 &v, const value_t d);
inline VectorC3 operator*(const value_t d, const VectorC3 &v);
inline VectorC3 operator/(const VectorC3 &v, const value_t d);
inline Complex  operator^(const VectorC3 &v1, const VectorC3 &v2);
inline Complex  operator^(const VectorC3 &vc, const VectorR3 &vr);
inline Complex  operator^(const VectorR3 &vr, const VectorC3 &vc);
inline VectorC3 operator*(const VectorR3& u, const Complex m);
inline VectorC3 operator*(const Complex m, const VectorR3& u);
inline VectorC3 operator*(const VectorC3& u, const Complex& c);
inline VectorC3 operator*(const Complex& c, const VectorC3& u);
inline value_t  norm(const VectorC3& vc);

inline VectorC3 operator+(const VectorC3 &v1, const VectorC3 &v2)
{
    auto tmp = v1;
    return tmp += v2;
}

inline VectorC3 operator+(const VectorR3 &v1, const VectorC3 &v2)
{
	return(VectorC3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z));
}

inline VectorC3 operator+(const VectorC3 &v1, const VectorR3 &v2)
{
	return(VectorC3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z));
}

inline VectorC3 operator-(const VectorC3 &v1, const VectorC3 &v2)
{
    auto tmp = v1;
    return tmp -= v2;
}

inline VectorC3 operator*(const VectorC3 &v1, const VectorC3 &v2)
{
    auto tmp = v1;
    return tmp *= v2;
}

inline VectorC3 operator*(const VectorC3 &vc, const VectorR3 &vr)
{
    return (VectorC3(
        vc.y*vr.z - vc.z*vr.y,
        vc.z*vr.x - vc.x*vr.z,
        vc.x*vr.y - vc.y*vr.x
    ));
}

inline VectorC3 operator*(const VectorR3 &vr, const VectorC3 &vc)
{
    return (VectorC3(
        vr.y*vc.z - vr.z*vc.y,
        vr.z*vc.x - vr.x*vc.z,
        vr.x*vc.y - vr.y*vc.x
    ));
}

inline VectorC3 operator*(const VectorC3 &v, const value_t d)
{
    auto tmp = v;
    return tmp *= d;
}

inline VectorC3 operator*(const value_t d, const VectorC3 &v)
{
    auto tmp = v;
    return tmp *= d;
}

inline VectorC3 operator/(const VectorC3 &v, const value_t d)
{
    auto tmp = v;
    return tmp /= d;
}

inline Complex  operator^(const VectorC3 &v1, const VectorC3 &v2)
{
    return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
}

inline Complex  operator^(const VectorC3 &vc, const VectorR3 &vr)
{
    return (vc.x*vr.x + vc.y*vr.y + vc.z*vr.z);
}

inline Complex  operator^(const VectorR3 &vr, const VectorC3 &vc)
{
    return (vc.x*vr.x + vc.y*vr.y + vc.z*vr.z);
}

inline VectorC3 operator*(const VectorR3& u, const Complex m)
{
    return VectorC3(u.x*m, u.y*m, u.z*m);
}

inline VectorC3 operator*(const Complex m, const VectorR3& u)
{
    return VectorC3(u.x*m, u.y*m, u.z*m);
}

inline VectorC3 operator*(const VectorC3& u, const Complex& c)
{
    return VectorC3(u.x*c, u.y*c, u.z*c);
}

inline VectorC3 operator*(const Complex& c, const VectorC3& u)
{
    return u * c;
}

inline VectorC3 conj(const VectorC3& vc)
{
    return VectorC3(std::conj(vc.x), std::conj(vc.y), std::conj(vc.z));
}

inline value_t norm(const VectorC3& vc)
{
    return vc.norm();
}

} // namespace component