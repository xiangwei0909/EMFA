#include "Triangle.h"

using namespace component;

void Triangle::setVertex(const int _A, const int _B, const int _C)
{
    nA = _A;
    nB = _B;
    nC = _C;
}

void Triangle::setVertex(const VectorR3& _VA, const VectorR3& _VB, const VectorR3& _VC)
{
    VertexA = _VA;
    VertexB = _VB;
    VertexC = _VC;
    getInfo();
}

void Triangle::setVertex(const VectorR3& _VA, const VectorR3& _VB, const VectorR3& _VC, const VectorR3 & _Nml)
{
    VertexA = _VA;
    VertexB = _VB;
    VertexC = _VC;
    Normal = _Nml;
}

void Triangle::getVertex(int& _A, int& _B, int& _C) const
{
    _A = nA;
    _B = nB;
    _C = nC;
}

void Triangle::getVertex(VectorR3& _VA, VectorR3& _VB, VectorR3& _VC)const
{
    _VA = VertexA;
    _VB = VertexB;
    _VC = VertexC;
}

VectorR3 Triangle::getCenter() const
{
    return (VertexA + VertexB + VertexC) / 3;
}

void Triangle::getInfo()
{
    auto EdgeAB = VertexB - VertexA;
    auto EdgeBC = VertexC - VertexB;
    auto EdgeCA = VertexA - VertexC;

    if ((EdgeAB^EdgeBC) < (EdgeBC^EdgeCA))
        Normal = EdgeAB*EdgeBC;
    else
        Normal = EdgeBC*EdgeCA;

    auto mag = Normal.Norm();
    assert(mag >= 0);
    Normal /= mag;
}
