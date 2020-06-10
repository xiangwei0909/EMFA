//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "Custom.h"
#include "VectorR3.h"

namespace component {

class Mesh;

class AIMRect {
public:
    AIMRect();
    ~AIMRect();

public:
    void        init(std::shared_ptr<Mesh> _pMesh, value_t _xspacing, value_t _yspacing, value_t _zspacing, int _points);
    void        getPointInfo(std::vector<int>& _pointNo, std::vector<value_t>& _xLength,
                             std::vector<value_t>& _yLength, std::vector<value_t>& _zLength, const VectorR3& _mid);
    int         getPointsNum() const;
    void        getXYZPointNum(int& _xNum, int& _yNum, int& _zNum) const;
    value_t     getDistance(int fx, int fy, int fz, int tx, int ty, int tz) const;
    void        getGroupDist(const VectorR3& _f, const VectorR3& _s, std::vector<value_t>& _dist,
                             std::vector<int>& _fPoint, std::vector<int>& _sPoint) const;
    void        getGroupDist(const VectorR3& _f, const VectorR3& _s, std::vector<int>& _dist,
                             std::vector<int>& _fPoint, std::vector<int>& _sPoint) const;

    void        reportInfo(Qostream& strm) const;
private:
    VectorR3    minBox, maxBox;
    value_t     x_spacing, y_spacing, z_spacing;
    int         xPointNum, yPointNum, zPointNum;
    int         aim_order_;
};

inline int AIMRect::getPointsNum() const
{
    return xPointNum * yPointNum * zPointNum;
}

inline void AIMRect::getXYZPointNum(int & _xNum, int & _yNum, int & _zNum) const
{
    _xNum = xPointNum;
    _yNum = yPointNum;
    _zNum = zPointNum;
}

}