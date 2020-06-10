//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#ifndef _COMMONEDGE_H_
#define _COMMONEDGE_H_
#include "Mesh.h"

namespace component {

class CommonEdge {
public:
    struct SCE {
        SCE(int _v1, int _v2, int _vp, int _fp)
            :v1(_v1), v2(_v2), vxp(_vp),
            vxm(-1), tpid(_fp), tmid(-1),bedge(-1),sedid(0)
        {
        }
		SCE(int _v1, int _v2, int _vp, int _fp, value_t _ep)
			:v1(_v1), v2(_v2), vxp(_vp), ep(_ep),
			vxm(-1), tpid(_fp), tmid(-1), bedge(-1), sedid(0), em(1.0f)
		{
		}
		SCE(int _v1, int _v2, int _vp, int _fp,Qstring _fla)
			:v1(_v1), v2(_v2), vxp(_vp), fla_p(_fla),
			vxm(-1), tpid(_fp), tmid(-1), bedge(-1), sedid(0), em(1.0f)
		{
		}
        int v1;
        int v2;
        int vxp;
        int vxm;
        int tpid;
        int tmid;
		int bedge;
		int sedid;
		value_t ep;
		value_t em;
        value_t length;
		Qstring fla_m;
		Qstring fla_p;
    };
	struct SPT {
		SPT(int _v1,int _v2,int _v3,int _tid, value_t _ep)
			:v1(_v1),v2(_v2),v3(_v3),tid(_tid),ep(_ep){}
		int v1;
		int v2;
		int v3;
		int tid;
		value_t ep;
	};
    CommonEdge();
    ~CommonEdge();
    CommonEdge(const CommonEdge&) = delete;
    CommonEdge& operator=(const CommonEdge&) = delete;
public:
    int         getPositiveFace(int index) const;
    value_t     getCommonEdgeLength(int index) const;
	void		getBEDGEandSEDID(int index,int &_bedge,int &_sedid) const;
	void		getSEDCommonEdgeSize(int index,int &_start,int &_end) const;
    void        getCommonEdgeTri(int index, int &_fp, int &_fm) const;
	void        getCommonEdge(int index, int &_v1, int &_v2) const;
	void		getCommonEdge(int index, int &_v1, int &_v2, int &_sedid) const;
	void        getCommonEdge(int index, int &_vp, int &_vm, int &_fp, int &_fm, value_t &_len) const;
    void        getCommonEdge(int index, int &_vp, int &_vm, int &_fp, int &_fm, value_t &_len,int& _bedge) const;
	void		getCommonEdge(int index, int &_vp, int &_vm, int &_fp, int &_fm, value_t &_len, int &_sedid, int &_bedge) const;
	void		getCommonEdge(int index, int &_vp, int &_vm, int &_fp, int &_fm, value_t &_len, value_t &_ep, value_t &_em,int &_bedge) const;
	void		getCommonEdge(int index, int &_v1, int &_v2, int &_fp, int &_fm, value_t &_len, Qstring& _fla_p, Qstring& _fla_m) const;
    const SCE&  getRWGRef(size_t i) const;
    bool        isCommonEdge(int index)const;
	int			getCommon_ce() const;
	void		getBoundary_ce(int &_b_ce1, int &_b_ce2) const;

    void        buildCommonEdge(std::shared_ptr<Mesh> mesh_ptr_);
	void		buildContinuousEdges(std::shared_ptr<Mesh> mesh_ptr_);

	void		buildSEDCommonEdge_com(std::shared_ptr<Mesh> mesh_ptr_);
	void		buildSEDCommonEdge_boundary(std::shared_ptr<Mesh> mesh_ptr_,value_t _Dx,value_t _Dy);//old
	void		buildSEDCommonEdge_boundary(std::shared_ptr<Mesh> mesh_ptr_, value_t _Dx, value_t _Dy,value_t _lamda);//new
	//void		buildSPINCommonEdge_boundary(std::shared_ptr<Mesh> mesh_ptr_,value_t _angle);
	void		combineSEDCommonEdge(std::shared_ptr<Mesh> mesh_ptr_);
	
	void		buildIBCCommonEdge(std::shared_ptr<Mesh> mesh_ptr_);

    int         getCommonEdgeNum()const;
	int			getPluseTriNum()const;
	void		getPluseTri(int index, int &fp, value_t &ep) const;
    void        clear();
    bool        writeLengthData(const std::string& filename, value_t lamda) const;
    void        reportInfo(Qostream& strm) const;
private:
    void        updateCommonEdge(int nA, int nB, int nC, int faceID);
	void		updateIBCCommonEdge(int nA, int nB, int nC, int faceID, Qstring fla);
    void        setEdgeTempAt(int index, int _vm, int _fp);
	void		setEdgeTempAt(int index, int _vm, int _fp, Qstring _eps);
    void        getEdgeTempAt(int index, int &_v1, int &_v2) const;
	void		getEdgeTempAt(int index, int &_v1, int &_v2, Qstring& _eps) const;
    bool        isCommonEdgeTemp(int index) const;
    int         getCommonEdgeTempNum()const;

private:
    std::vector<SCE> m_ceListTemp;
    std::vector<SCE> m_ceList;
	std::vector<SCE> m_beList;

	std::vector<SCE> d_ceListTemp;
	std::vector<SCE> d_ceList;
	std::vector<SPT> d_ptList;

private:
	std::vector<SCE> m_ceList_SED0;
	std::vector<SCE> m_ceList_SED1;
	std::vector<SCE> m_ceList_SED2;
	std::vector<SCE> m_ceList_SED3;    
	std::vector<SCE> m_ceList_SED4;
	std::vector<SCE> m_ceList_SED5;
	std::vector<SCE> m_ceList_SED6;
	std::vector<SCE> m_ceList_SED7;
	std::vector<SCE> m_ceList_SED8;
	int				 m_ce_size[10];
private:
    int              totalEdges;
	int				 boundaryEdges;
	int				 common_ce;
	int              b_ce1;
	int              b_ce2;
	
};

inline int CommonEdge::getPositiveFace(int index) const
{
    return m_ceList[index].tpid;
}

inline value_t CommonEdge::getCommonEdgeLength(int index) const
{
    return m_ceList[index].length;
}

inline void CommonEdge::getBEDGEandSEDID(int index,int &_bedge,int &_sedid) const
{
	_bedge = m_ceList[index].bedge;
	_sedid = m_ceList[index].sedid;
}

inline void CommonEdge::getCommonEdgeTri(int index, int & _fp, int & _fm) const
{
    auto& edge = m_ceList[index];
    _fp = edge.tpid;
    _fm = edge.tmid;
}

inline void CommonEdge::getSEDCommonEdgeSize(int index, int &_start,int &_end) const
{
	_start = m_ce_size[index];
	_end = m_ce_size[index + 1];
}

inline void CommonEdge::getCommonEdge(int index, int & _v1, int & _v2) const
{
    _v1 = m_ceList[index].v1;
    _v2 = m_ceList[index].v2;
}

inline void CommonEdge::getCommonEdge(int index, int &_v1, int &_v2, int &_sedid) const
{
	_v1 = m_ceList[index].v1;
	_v2 = m_ceList[index].v2;
	_sedid = m_ceList[index].sedid;
}

inline const CommonEdge::SCE & component::CommonEdge::getRWGRef(size_t i) const
{
    return m_ceList[i];
}

inline bool CommonEdge::isCommonEdge(int index) const
{
    return m_ceList[index].tmid != -1;
}

inline int CommonEdge::getCommonEdgeNum() const
{
    return static_cast<int>(m_ceList.size());
}

inline int CommonEdge::getPluseTriNum() const
{
	return static_cast<int>(d_ptList.size());
}

inline void CommonEdge::getPluseTri(int index, int &fp, value_t &ep) const
{
	fp = d_ptList[index].tid;
	ep = d_ptList[index].ep;
}

inline void CommonEdge::clear()
{
    m_ceList.clear();
}

inline bool CommonEdge::isCommonEdgeTemp(int index) const
{
    return m_ceListTemp[index].tmid != -1;
}

inline int CommonEdge::getCommonEdgeTempNum() const
{
    return static_cast<int>(m_ceListTemp.size());
}

inline int CommonEdge::getCommon_ce() const
{
	return common_ce;
}

inline void CommonEdge::getBoundary_ce(int &_b_ce1, int &_b_ce2) const
{
	_b_ce1 = b_ce1;
	_b_ce2 = b_ce2;
}
} // namespace component
#endif