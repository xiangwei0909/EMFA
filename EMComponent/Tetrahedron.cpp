#include "Tetrahedron.h"

using namespace component;

Tetrahedron::Tetrahedron(void)
{
}


Tetrahedron::~Tetrahedron(void)
{
}

void Tetrahedron::SetVertices(const int a, const int b, const int c, const int d)
{
	nA = a;
	nB = b;
	nC = c;
	nD = d;
}

void Tetrahedron::SetVertices(const VectorR3& vA, const VectorR3& vB, const VectorR3& vC, const VectorR3& vD)
{
	VertexA = vA;
	VertexB = vB;
	VertexC = vC;
	VertexD = vD;
	ComputeTetraHedronVolume();
}

void Tetrahedron::GetVertices(int& a, int& b, int& c, int& d) const
{
	a = nA;
	b = nB;
	c = nC;
	d = nD;
}

void Tetrahedron::GetVertices(VectorR3 & vertA, VectorR3 & vertB, VectorR3 & vertC, VectorR3 & vertD)  const
{
	vertA = VertexA;
	vertB = VertexB;
	vertC = VertexC;
	vertD = VertexD;
}

void Tetrahedron::ComputeTetraHedronVolume()
{
	VectorR3 EdgeAB = VertexB - VertexA;
	VectorR3 EdgeAC = VertexC - VertexA;
	VectorR3 EdgeAD = VertexD - VertexA;

	VectorR3 Buff = EdgeAB*EdgeAC;
	Volume = Buff^EdgeAD/6.0f;
}

// PreCalcInfo takes the vertex values and computes information to
//		help with intersections with rays.
// void Triangle::PreCalcInfo()
// {
// 	VectorR3 EdgeAB = VertexB - VertexA;
// 	VectorR3 EdgeBC = VertexC - VertexB;
// 	VectorR3 EdgeCA = VertexA - VertexC;
// 
// 	if ( (EdgeAB^EdgeBC) < (EdgeBC^EdgeCA) ) {
// 		Normal = EdgeAB*EdgeBC;
// 	}
// 	else {
// 		Normal = EdgeBC*EdgeCA;
// 	}
// 	float mag = Normal.Norm();
// 	if ( mag>0.0 ) {
// 		Normal /= mag;		// Unit vector to triangle's plane
// 	}
// }