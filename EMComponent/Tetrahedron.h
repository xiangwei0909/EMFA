//**********************************************************
// Author: Xiang Wei
// License: MIT
//**********************************************************
#pragma once
#include"VectorR3.h"


namespace component {
class Tetrahedron
{
public:
	Tetrahedron(void);
	~Tetrahedron(void);

public:
	void SetVertices(const int a, const int b, const int c, const int d);
	void SetVertices(const VectorR3& vA, const VectorR3& vB, const VectorR3& vC, const VectorR3& vD);
	void SetEpsilon(const float a) { m_ep = a; };
	void GetVertices(int& a, int& b, int& c, int& d) const;
	void GetVertices(VectorR3 & vertA, VectorR3 & vertB, VectorR3 & vertC, VectorR3 & vertD) const;

	void ComputeTetraHedronVolume();
	float GetVolume() { return Volume; };
	void GetEpsilon(float& a) { a = m_ep; };

public:
	void SetIsVisible(bool is) { isVisible = is; };
	bool GetIsVisible() { return isVisible; };

private:
	int nA;
	int nB;
	int nC;
	int nD;
	VectorR3 VertexA;
	VectorR3 VertexB;
	VectorR3 VertexC;
	VectorR3 VertexD;
	bool isVisible;
	float Volume;
	float m_ep;
	//VectorR3 Normal;	// Unit normal to the plane of triangle
};
}