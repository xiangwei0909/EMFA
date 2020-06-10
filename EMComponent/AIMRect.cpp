#include "AIMRect.h"
#include "Mesh.h"

using namespace component;

AIMRect::AIMRect()
{
}


AIMRect::~AIMRect()
{
}

void AIMRect::init(std::shared_ptr<Mesh> _pMesh, value_t _xspacing, value_t _yspacing, value_t _zspacing,int _points)
{
    x_spacing = _xspacing;
    y_spacing = _yspacing;
    z_spacing = _zspacing;
    aim_order_   = _points;

    _pMesh->getBoundary(minBox, maxBox);
    minBox = minBox - VectorR3(_xspacing, _yspacing, _zspacing);    //  all direction extend spacing
    maxBox = maxBox + VectorR3(_xspacing, _yspacing, _zspacing);

    xPointNum = static_cast<int>((maxBox.x - minBox.x) / _xspacing + 0.5) + 1;
    yPointNum = static_cast<int>((maxBox.y - minBox.y) / _yspacing + 0.5) + 1;
    zPointNum = static_cast<int>((maxBox.z - minBox.z) / _zspacing + 0.5) + 1;
}

void AIMRect::getPointInfo(std::vector<int>& _pointNo, std::vector<value_t>& _xLength,
                            std::vector<value_t>& _yLength, std::vector<value_t>& _zLength,
                                const VectorR3 & _mid)
{
    _pointNo.reserve(aim_order_ * aim_order_ * aim_order_);
    _xLength.reserve(aim_order_);
    _yLength.reserve(aim_order_);
    _zLength.reserve(aim_order_);
    auto tmpMid = _mid - minBox;
    int xPos = static_cast<int>(tmpMid.x / x_spacing + 0.5);
    int yPos = static_cast<int>(tmpMid.y / y_spacing + 0.5);
    int zPos = static_cast<int>(tmpMid.z / z_spacing + 0.5);

    int range = aim_order_ / 2;
    int indexX = yPointNum * zPointNum;
    int indexY = zPointNum;
    for (int x = -range; x <= range; ++x)
    {
        auto idx = (xPos + x) * indexX;
        for (int y = -range; y <= range; ++y)
        {
            auto idy = (yPos + y) * indexY;
            for (int z = -range; z <= range; ++z)
            {
                auto idz = (zPos + z);
                _pointNo.push_back(idx + idy + idz);
            }
        }
    }

    for (int i = -range; i <= range; ++i)
    {
        auto idx = xPos + i;
        auto idy = yPos + i;
        auto idz = zPos + i;

        _xLength.push_back(idx * x_spacing - tmpMid.x);
        _yLength.push_back(idy * y_spacing - tmpMid.y);
        _zLength.push_back(idz * z_spacing - tmpMid.z);
    }
}

value_t AIMRect::getDistance(int fx, int fy, int fz, int tx, int ty, int tz) const
{
    auto x_dist = (fx - tx) * x_spacing;
    auto y_dist = (fy - ty) * y_spacing;
    auto z_dist = (fz - tz) * z_spacing;

    return std::sqrt(x_dist * x_dist + y_dist * y_dist + z_dist * z_dist);
}

void AIMRect::getGroupDist(const VectorR3 & _f, const VectorR3 & _s, std::vector<value_t>& _dist,
                                std::vector<int>& _fPoint, std::vector<int>& _sPoint) const
{
    auto cap = aim_order_ * aim_order_ * aim_order_;
    _dist.reserve(cap * cap);
    _fPoint.reserve(cap);
    _sPoint.reserve(cap);

    auto fMid = _f - minBox;
    int fxPos = static_cast<int>(fMid.x / x_spacing + 0.5);
    int fyPos = static_cast<int>(fMid.y / y_spacing + 0.5);
    int fzPos = static_cast<int>(fMid.z / z_spacing + 0.5);

    auto sMid = _s - minBox;
    int sxPos = static_cast<int>(sMid.x / x_spacing + 0.5);
    int syPos = static_cast<int>(sMid.y / y_spacing + 0.5);
    int szPos = static_cast<int>(sMid.z / z_spacing + 0.5);

    int range = aim_order_ / 2;

    for (int fx = -range; fx <= range; ++fx)
    {
        for (int fy = -range; fy <= range; ++fy)
        {
            for (int fz = -range; fz <= range; ++fz)
            {
                for (int sx = -range; sx <= range; ++sx)
                {
                    auto x_sp = fxPos + fx - (sxPos + sx);
                    auto x_dist = x_sp * x_spacing;
                    x_dist = x_dist * x_dist;
                    for (int sy = -range; sy <= range; ++sy)
                    {
                        auto y_sp = (fyPos + fy - (syPos + sy));
                        auto y_dist = y_sp * y_spacing;
                        y_dist = y_dist * y_dist;
                        for (int sz = -range; sz <= range; ++sz)
                        {
                            auto z_sp = (fzPos + fz - (szPos + sz));
                            auto z_dist = z_sp * z_spacing;
                            z_dist = z_dist * z_dist;
                            _dist.push_back(std::sqrt(x_dist + y_dist + z_dist));
                        }
                    }
                }
            }
        }
    }

    int indexX = yPointNum * zPointNum;
    int indexY = zPointNum;
    for (int x = -range; x <= range; ++x)
    {
        auto idfx = (fxPos + x) * indexX;
        auto idsx = (sxPos + x) * indexX;
        for (int y = -range; y <= range; ++y)
        {
            auto idfy = (fyPos + y) * indexY;
            auto idsy = (syPos + y) * indexY;
            for (int z = -range; z <= range; ++z)
            {
                auto idfz = (fzPos + z);
                auto idsz = (szPos + z);
                _fPoint.push_back(idfx + idfy + idfz);
                _sPoint.push_back(idsx + idsy + idsz);
            }
        }
    }
}

void AIMRect::getGroupDist(const VectorR3 & _f, const VectorR3 & _s, std::vector<int>& _dist,
                            std::vector<int>& _fPoint, std::vector<int>& _sPoint) const
{
    auto cap = aim_order_ * aim_order_ * aim_order_;
    _dist.reserve(cap * cap);
    _fPoint.reserve(cap);
    _sPoint.reserve(cap);

    auto fMid = _f - minBox;
    int fxPos = static_cast<int>(fMid.x / x_spacing + 0.5);
    int fyPos = static_cast<int>(fMid.y / y_spacing + 0.5);
    int fzPos = static_cast<int>(fMid.z / z_spacing + 0.5);

    auto sMid = _s - minBox;
    int sxPos = static_cast<int>(sMid.x / x_spacing + 0.5);
    int syPos = static_cast<int>(sMid.y / y_spacing + 0.5);
    int szPos = static_cast<int>(sMid.z / z_spacing + 0.5);

    int range = aim_order_ / 2;

    for (int fx = -range; fx <= range; ++fx)
    {
        for (int fy = -range; fy <= range; ++fy)
        {
            for (int fz = -range; fz <= range; ++fz)
            {
                for (int sx = -range; sx <= range; ++sx)
                {
                    auto x_sp = std::abs(fxPos + fx - (sxPos + sx));
                    for (int sy = -range; sy <= range; ++sy)
                    {
                        auto y_sp = std::abs(fyPos + fy - (syPos + sy));
                        for (int sz = -range; sz <= range; ++sz)
                        {
                            auto z_sp = std::abs(fzPos + fz - (szPos + sz));
                            _dist.push_back((x_sp * yPointNum + y_sp) * zPointNum + z_sp);
                        }
                    }
                }
            }
        }
    }

    int indexX = yPointNum * zPointNum;
    int indexY = zPointNum;
    for (int x = -range; x <= range; ++x)
    {
        auto idfx = (fxPos + x) * indexX;
        auto idsx = (sxPos + x) * indexX;
        for (int y = -range; y <= range; ++y)
        {
            auto idfy = (fyPos + y) * indexY;
            auto idsy = (syPos + y) * indexY;
            for (int z = -range; z <= range; ++z)
            {
                auto idfz = (fzPos + z);
                auto idsz = (szPos + z);
                _fPoint.push_back(idfx + idfy + idfz);
                _sPoint.push_back(idsx + idsy + idsz);
            }
        }
    }
}


void AIMRect::reportInfo(Qostream & strm) const
{
    auto box = maxBox - minBox;
    auto flagState = strm.flags();
    strm << std::fixed << std::setprecision(2);
    strm << HEADING "AIM Information:\n" TRAILING;
    strm << LEVEL1 "Rect min:" << '(' << minBox.x << ',' << minBox.y << ',' << minBox.z << ")\n"
        << LEVEL1 "Rect max:" << '(' << maxBox.x << ',' << maxBox.y << ',' << maxBox.z << ")\n"
        << LEVEL1 "Rect size:" << '(' << box.x << ',' << box.y << ',' << box.z << ")\n"
        << LEVEL1 "X|Y|Z spacing:" << x_spacing << "|" << y_spacing << "|" << z_spacing << "\n"
        << LEVEL1 "X|Y|Z points:" << xPointNum << "|" << yPointNum << "|" << zPointNum << "\n";
    strm.flags(flagState);
}
