//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "Custom.h"
#include "VectorR3.h"

namespace component {
class Mesh;

class ContourEdge {
public:
    enum EdgeType { CONFORMAL, NONCONFORMAL };
    struct RWG {
        RWG(int _v1, int _v2, int _vp, int _fp)
            :v1ID(_v1), v2ID(_v2), vpID(_vp),
            vmID(-1), fpID(_fp), fmID(-1) { } 
        int     v1ID;
        int     v2ID;
        int     vpID;
        int     vmID;
        int     fpID;
        int     fmID;
        value_t length;
    };
    struct HRWG {
        HRWG(int nv1, int nv2, int nvx, int fid, int id, EdgeType tp, value_t len)
        {
            v1 = nv1;
            v2 = nv2;
            vx = nvx;
            tid = fid;
            idx = id;
            type = tp;
            length = len;
        }
        int     v1;     //  start point
        int     v2;     //  end point
        int     vx;
        int     tid;    // triangle id
        int     idx;     // common edge id
        EdgeType    type;
        value_t     length;
    };
    struct Line {
        VectorR3 start, end;
        value_t length() const { return (end - start).Norm(); }
    };
    using ContourLine = std::vector<std::pair<int, Line>>;
    //using ContourLine = std::vector<std::pair<int, value_t>>;

public:
    ContourEdge();
    ~ContourEdge();

    ContourEdge(const ContourEdge&) = delete;
    ContourEdge& operator=(const ContourEdge&) = delete;

public:
    void buildEdge(std::shared_ptr<Mesh> pmodel, policy::BasisFunction type)
    {
        basis_type_ = type;
        if (basis_type_ == policy::HALF_RWG)
            buildContourEdge(pmodel);
        else
            buildMixedEdge(pmodel);
    }
    int getUnknownsNum() const
    {
        return basis_type_  == policy::HALF_RWG ? getUnknownsNum1() : getUnknownsNum2();
    }
    const HRWG& getHRWGRef(size_t index) const
    {
        if (basis_type_ == policy::HALF_RWG)
            return getHRWGRef1(index);
        else
            return getHRWGRef2(index);
    }
    // IP_IE_DDM
    bool getContourLine(size_t index, int trangle_id, Line* ln) const
    {
        for (const auto& line : line_list_[index])
        {
            if (line.first == trangle_id)
            {
                *ln = line.second;
                return true;
            }
        }
        return false;
    }
    // IEDG
    bool getContourLineLength(size_t index, int trangle_id, value_t* len) const
    {
        for (const auto& line : line_list_[index])
        {
            if (line.first == trangle_id)
            {
                *len = (line.second.end - line.second.start).Norm();
                return true;
            }
        }
        return false;
    }

    void clear();
    void reportInfo(Qostream& strm) const;
    int getRWGNum() const
    {
        return static_cast<int>(fce_list_.size());
    }
    int getHRWGNum() const
    {
        return static_cast<int>(hce_list_.size());
    }
    const RWG& getRWGRef(size_t index) const
    {
        return fce_list_[index];
    }
    void Debug(Qostream& strm, std::shared_ptr<Mesh> pmode) const
    {
        if (basis_type_ == policy::HALF_RWG)
            Debug1(strm, pmode);
        else
            Debug2(strm, pmode);
    }
private:
    void buildContourEdge(std::shared_ptr<Mesh> pmodel);
    void buildMixedEdge(std::shared_ptr<Mesh> pmodel);
    void buildCommonEdge(std::vector<RWG>& celist, std::shared_ptr<Mesh> pmode);

    int getUnknownsNum1() const
    {
        return static_cast<int>(ce_list_.size());
    }
    int getUnknownsNum2() const
    {
        return static_cast<int>(fce_list_.size() + hce_list_.size());
    }
    const HRWG& getHRWGRef1(size_t index) const
    {
        return ce_list_[index];
    }
    const HRWG& getHRWGRef2(size_t index) const
    {
        return hce_list_[index];
    }
    void Debug1(Qostream& strm, std::shared_ptr<Mesh> pmode) const;

    void Debug2(Qostream& strm, std::shared_ptr<Mesh> pmode) const;

private:
    policy::BasisFunction basis_type_;
    size_t total_edges_;
    std::vector<HRWG> ce_list_;
    std::vector<ContourLine> line_list_;
    std::vector<RWG> fce_list_;
    std::vector<HRWG> hce_list_;
};

} // namespace component