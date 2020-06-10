//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include"VectorR3.h"

namespace component {

class Triangle {
public:
    void        setIsCommonFace(int _No);
    bool        getIsCommonFace() const;
    void        setID(int _id);
    int         getID() const;
    void        setVertex(const int _A, const int _B, const int _C);
    void        setVertex(const VectorR3& _A, const VectorR3& _B, const VectorR3& _C);
    void        setVertex(const VectorR3& _A, const VectorR3& _B, const VectorR3& _C, const VectorR3& _Nml);
	void		setEps(const value_t& _eps) { eps = _eps;}
	void		setFla(const Qstring& _fla) { fla = _fla; }
	void		setgroupID(const int _groupID) { groupID = _groupID; }
    void        getVertex(int& _A, int& _B, int& _C)const;
    void        getVertex(VectorR3& _VA, VectorR3& _VB, VectorR3& _VC)const;

    void        setVisible(bool _is) { isVisible = _is; }
	void		getEps(value_t& _eps) { _eps = eps; }
	void		getFla(Qstring& _fla) { _fla = fla; }
    bool        getVisible()const { return isVisible; }
    VectorR3            getCenter() const;
    const VectorR3&     getNormal()const { return Normal; }
	int			getgroupID()const { return groupID; }
private:
    void        getInfo(); //求解面法向量
private:
    int         nA;
    int         nB;
    int         nC;
    bool        isVisible;
    bool        isCommonFace;
    int         id;
    VectorR3    VertexA;
    VectorR3    VertexB;
    VectorR3    VertexC;
    VectorR3    Normal;
	value_t		eps;
	Qstring		fla;
	int			groupID;

};

inline void Triangle::setIsCommonFace(int _No)
{
    isCommonFace = _No == 1 ? true : false;
}

inline bool Triangle::getIsCommonFace() const
{
    return isCommonFace;
}

inline void Triangle::setID(int _id)
{
    id = _id;
}

inline int Triangle::getID() const
{
    return id;
}

} // namespace component