#include "ContourEdge.h"
#include "Mesh.h"
#include "tools.h"

using namespace component;

ContourEdge::ContourEdge()
{
}


ContourEdge::~ContourEdge()
{
}

void ContourEdge::buildContourEdge(std::shared_ptr<Mesh> pmodel)
{
    std::vector<RWG> edge_list, left_comedge, left_list;
    buildCommonEdge(edge_list, pmodel);
    total_edges_ = edge_list.size();
    ce_list_.reserve(3 * pmodel->getTriangleNum());
    int ce_edge_no = 0;
    for (const auto& edge : edge_list)
    {
        if (edge.fmID != -1)
        {
            ce_list_.push_back(HRWG(edge.v1ID, edge.v2ID, edge.vpID, edge.fpID,
                                    ce_edge_no, CONFORMAL, edge.length));
            ce_list_.push_back(HRWG(edge.v2ID, edge.v1ID, edge.vmID, edge.fmID,
                                    ce_edge_no, CONFORMAL, edge.length));
            ++ce_edge_no;
        }
        else
            left_comedge.push_back(edge);
    }

    for (auto cur = left_comedge.begin(); cur != left_comedge.end(); ++cur)
    {
        auto v1 = pmodel->getVertex(cur->v1ID);
        auto v2 = pmodel->getVertex(cur->v2ID);
        auto iter = cur + 1;
        for (; iter != left_comedge.end(); ++iter)
        {
            auto v3 = pmodel->getVertex(iter->v1ID);
            auto v4 = pmodel->getVertex(iter->v2ID);
            if (v1 == v4 && v2 == v3)
            {
                ce_list_.push_back(HRWG(cur->v1ID, cur->v2ID, cur->vpID, cur->fpID,
                                        ce_edge_no, CONFORMAL, cur->length));
                ce_list_.push_back(HRWG(iter->v1ID, iter->v2ID, iter->vpID, iter->fpID,
                                        ce_edge_no, CONFORMAL, iter->length));
                ++ce_edge_no;
                break;
            }
        }
        if (iter == left_comedge.end())
            left_list.push_back(*cur);
    }
  
    const double threshold = 0.95f; // threshold depend on difference of the scale between adjacent subdomains
    ce_edge_no = 0;
    for (auto cur = left_list.begin(); cur != left_list.end(); ++cur)
    {
        std::vector<std::pair<int, Line>> lines;
        auto v1 = pmodel->getVertex(cur->v1ID);
        auto v2 = pmodel->getVertex(cur->v2ID);
        auto v1v2 = v2 - v1;
        double nv1v2 = v1v2.Norm();
        for (auto iter = left_list.begin(); iter != left_list.end(); ++iter)
        {
            if (cur == iter) continue;
            auto v3 = pmodel->getVertex(iter->v1ID);
            auto v4 = pmodel->getVertex(iter->v2ID);
            VectorR3 t1 = v4 - v1, t2 = v3 - v1, t3 = v4 - v2, t4 = v3 - v2;
            if ((v1v2 ^ t2) <= 0 || (v1v2 ^ t3) >= 0)    continue;

            //  avoid to divide zero if v1-v4 and v2-v3 overlap
            double t1sign = (v1 == v4) ? 1.0 : (v1v2 ^ t1) / (nv1v2 * t1.Norm());
            double t4sign = (v2 == v3) ? 1.0 : (v1v2 ^ t4) / (nv1v2 * t4.Norm());
            if (abs(t1sign) > threshold && abs(t4sign) > threshold)
            {
                Line segment;
                int start = t1sign < 0 ? cur->v1ID : iter->v2ID;
                int end = t4sign < 0 ? iter->v1ID : cur->v2ID;
                segment.start = pmodel->getVertex(start);
                segment.end = pmodel->getVertex(end);
                lines.push_back(std::make_pair(iter->fpID, segment));
            }
        }
        if (!lines.empty())
        {
            ce_list_.push_back(HRWG(cur->v1ID, cur->v2ID, cur->vpID, cur->fpID,
                                    ce_edge_no, NONCONFORMAL, cur->length));
            ++ce_edge_no;
            line_list_.push_back(std::move(lines));
        }
    }
}

void ContourEdge::clear()
{
    ce_list_.clear();
    line_list_.clear();
}

void ContourEdge::reportInfo(Qostream & strm) const
{
    auto old_state = strm.flags();
    if (basis_type_ == policy::HALF_RWG)
    {
        size_t mem = ce_list_.size() * ce_list_.size() * 2 * sizeof(value_t) / (1024 * 1024);
        size_t rwgmem = sizeof(HRWG) * ce_list_.size() / (1024 * 1024);
        strm << HEADING "ContourEdge Information\n" TRAILING;
        strm << LEVEL1 "Total edges:" << total_edges_ << '\n'
            << LEVEL1 "Total half-RWGs:" << ce_list_.size() << '\n'
            << LEVEL2 "conformal: " << ce_list_.size() - line_list_.size() << '\n'
            << LEVEL2 "nonconformal: " << line_list_.size() << '\n'
            << LEVEL1 "Memory prediction: " << mem << " MB" << FORMAT_MEMORY(mem) << '\n'
            << LEVEL1 "Memory detail: " << '\n'
            << LEVEL2 "RWG:" << rwgmem << " MB" << FORMAT_MEMORY(rwgmem) << '\n';
    }
    else
    {
        size_t mem = (fce_list_.size() + hce_list_.size()) * (fce_list_.size() + hce_list_.size()) * 2 * sizeof(value_t) / (1024 * 1024);
        size_t rwgmem = (sizeof(HRWG) * hce_list_.size() + sizeof(RWG) * fce_list_.size()) / (1024 * 1024);
        strm << HEADING "CommonEdge Information\n" TRAILING;
        strm << LEVEL1 "Total edges:" << total_edges_ << '\n'
            << LEVEL1 "Total valid RWGs:" << fce_list_.size() + hce_list_.size() << '\n'
            << LEVEL2 "full-RWGs:" << fce_list_.size() << '\n'
            << LEVEL2 "half-RWGs: " << hce_list_.size() << '\n'
            << LEVEL1 "Memory prediction: " << mem << " MB" << FORMAT_MEMORY(mem) << '\n'
            << LEVEL1 "Memory detail: " << '\n'
            << LEVEL2 "RWG:" << rwgmem << " MB" << FORMAT_MEMORY(rwgmem) << '\n';
    }
    strm.flags(old_state);
}

void ContourEdge::Debug1(Qostream & strm, std::shared_ptr<Mesh> pmode) const
{
    auto old_state = strm.flags();
    strm << std::left << "ContourEdge Information(Debug)\n";
    strm << " conformal edges: " << ce_list_.size() - line_list_.size() << '\n'
        << " nonconformal edges: " << line_list_.size() << std::endl;
    for (size_t i = 0; i < line_list_.size(); ++i)
    {
        const auto& cur = line_list_[i];
        strm << i << ": " << cur.size() << "    [";
        for (const auto& e : cur)
            strm << e.first << " (" << (e.second.end - e.second.start).Norm() << ") ";
        strm << "]" << std::endl;
    }
    strm.flags(old_state);
}

void ContourEdge::buildCommonEdge(std::vector<RWG>& edge_list, std::shared_ptr<Mesh> pmode)
{
    int tri_num = pmode->getTriangleNum();
    edge_list.reserve((tri_num >> 1) + tri_num);
    int nA, nB, nC;
    for (int i = 0; i < tri_num; ++i)
    {
        const auto& tri = pmode->getTriangleRef(i);
        tri.getVertex(nA, nB, nC);
        bool new_edge1 = true, new_edge2 = true, new_edge3 = true;
        int nv1, nv2;
        for (auto& edge : edge_list)
        {
            if (edge.fmID != -1)
                continue;
            nv1 = edge.v1ID;
            nv2 = edge.v2ID;
            if (nA == nv2 && nB == nv1)
            {
                new_edge1 = false;
                edge.vmID = nC;
                edge.fmID = i;
            }
            else if (nB == nv2 && nC == nv1)
            {
                new_edge2 = false;
                edge.vmID = nA;
                edge.fmID = i;
            }
            else if (nC == nv2 && nA == nv1)
            {
                new_edge3 = false;
                edge.vmID = nB;
                edge.fmID = i;
            }
        }
        if (new_edge1)   edge_list.push_back(RWG(nA, nB, nC, i));
        if (new_edge2)   edge_list.push_back(RWG(nB, nC, nA, i));
        if (new_edge3)   edge_list.push_back(RWG(nC, nA, nB, i));
    }
    VectorR3 v1, v2;
    for (auto& edge : edge_list)
    {
        v1 = pmode->getVertex(edge.v1ID);
        v2 = pmode->getVertex(edge.v2ID);
        edge.length = (v2 - v1).Norm();
    }
}



void component::ContourEdge::buildMixedEdge(std::shared_ptr<Mesh> pmodel)
{
    std::vector<RWG> edge_list, left_list;
    buildCommonEdge(edge_list, pmodel);
    total_edges_ = edge_list.size();
    for (const auto& edge : edge_list)
    {
        if (edge.fmID != -1)
            fce_list_.push_back(edge);
        else
            left_list.push_back(edge);
    }

    const double threshold = 0.95f; // threshold depend on difference of the scale between adjacent subdomains
    int ce_edge_no = 0;
    for (auto cur = left_list.begin(); cur != left_list.end(); ++cur)
    {
        std::vector<std::pair<int, Line>> lines;
        auto v1 = pmodel->getVertex(cur->v1ID);
        auto v2 = pmodel->getVertex(cur->v2ID);
        auto v1v2 = v2 - v1;
        double nv1v2 = v1v2.Norm();
        for (auto iter = left_list.begin(); iter != left_list.end(); ++iter)
        {
            if (cur == iter) continue;
            auto v3 = pmodel->getVertex(iter->v1ID);
            auto v4 = pmodel->getVertex(iter->v2ID);
            VectorR3 t1 = v4 - v1, t2 = v3 - v1, t3 = v4 - v2, t4 = v3 - v2;
            if ((v1v2 ^ t2) <= 0 || (v1v2 ^ t3) >= 0)    continue;

            //  avoid to divide zero if v1-v4 and v2-v3 overlap
            double t1sign = (v1 == v4) ? 1.0 : (v1v2 ^ t1) / (nv1v2 * t1.Norm());
            double t4sign = (v2 == v3) ? 1.0 : (v1v2 ^ t4) / (nv1v2 * t4.Norm());
            if (abs(t1sign) > threshold && abs(t4sign) > threshold)
            {
                Line segment;
                int start = t1sign < 0 ? cur->v1ID : iter->v2ID;
                int end = t4sign < 0 ? iter->v1ID : cur->v2ID;
                segment.start = pmodel->getVertex(start);
                segment.end = pmodel->getVertex(end);
                lines.push_back(std::make_pair(iter->fpID, segment));
            }
        }
        if (!lines.empty())
        {
            hce_list_.push_back(HRWG(cur->v1ID, cur->v2ID, cur->vpID, cur->fpID,
                ce_edge_no, NONCONFORMAL, cur->length));
            ++ce_edge_no;
            line_list_.push_back(std::move(lines));
        }
    }
}

void component::ContourEdge::Debug2(Qostream & strm, std::shared_ptr<Mesh> pmode) const
{
    auto old_state = strm.flags();
    strm << std::left << "ContourEdge Information(Debug2):\n";
    strm << " RWGs: " << fce_list_.size() << '\n'
        << " HRWGs: " << hce_list_.size() << std::endl;
    for (size_t i = 0; i < line_list_.size(); ++i)
    {
        const auto& cur = line_list_[i];
        strm << i + 1 << ": " << cur.size() << "    [";
        for (const auto& e : cur)
            strm << e.first << " (" << (e.second.end - e.second.start).Norm() << ") ";
        strm << "]" << std::endl;
    }
    strm.flags(old_state);
}