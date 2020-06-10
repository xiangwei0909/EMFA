//**********************************************************
// Author: Xiang Wei
// License: MIT
//**********************************************************

#ifndef _COMMONTRIANGLE_H_
#define _COMMONTRIANGLE_H_
#include "Mesh.h"

namespace component {

class CommonTriangle {
public:
	struct SCT
	{
		SCT(int _v1, int _v2, int _v3, int _vp, int _fp, value_t _eps)
			:v1(_v1), v2(_v2), v3(_v3), vxp(_vp), vxm(-1), tpid(_fp), tmid(-1), bedge(-1), eps_p(_eps),eps_m(1.0f),sedid(0)
		{}

		int v1;
		int v2;
		int v3;
		int vxp;
		int vxm;
		int tpid;
		int tmid;
		int bedge;
		int sedid;
		value_t area;
		value_t eps_p;
		value_t eps_m;
	};

	struct h_SWG
	{
		h_SWG(int _v1,int _v2,int _v3,int _vx,int _t,value_t _eps)
			:v1(_v1),v2(_v2),v3(_v3),vx(_vx),tid(_t),eps(_eps)
		{}

		int v1;
		int v2;
		int v3;
		int vx;
		int tid;
		value_t eps;
		value_t area;
	};

	CommonTriangle();
	~CommonTriangle();
	CommonTriangle(const CommonTriangle&) = delete;
	CommonTriangle& operator=(const CommonTriangle&) = delete;

public:
	int         getPositiveTriangle(int index) const;
	value_t     getCommonTriangleArea(int index) const;
	int			getbedge(int index) const;
	void	    getBEDGEandSEDID(int index,int &_bedge,int &_sedid) const;
	void		getSEDCommonTriangleSize(int index, int &_start,int &_end) const;
	void        getCommonTriangleTet(int index, int &_fp, int &_fm) const;
	void		getTet_eps(int index, value_t &_eps) const;
	void        getCommonTriangle(int index, int &_v1, int &_v2,int &_v3)const;
	void        getCommonTriangle(int index, int &_vp, int &_vm, int &_fp, int &_fm, value_t &area) const;
	void		getCommonTriangle(int index, int &_vp, int &_vm, int &_fp, int &_fm, value_t &area, int &_sedid, int &_bedge) const;
	void        getCommonTriangle(int index, int &_vp, int &_vm, int &_fp, int &_fm, value_t &area,value_t &_epsp,value_t &_epsm) const;
	void		getCommonTriangle(int index, int &_vp, int &_vm, int &_fp, int &_fm, value_t &area, value_t &_epsp, value_t &_epsm, int &_sedid, int &_bedge) const;
	void		gethSWG(int index, int &_v1, int &_v2, int &_v3, int &_vx, int &_tid, value_t &_eps, value_t &_area) const;
	const SCT&  getSWGRef(size_t i) const;
	bool        isCommonTriangle(int index)const;
	int			getCommon_ct() const;

	void        buildCommonTriangle(std::shared_ptr<Mesh> mesh_ptr_);
	void		buildCommonTriangleContinuous(std::shared_ptr<Mesh> mesh_ptr_);
	void		buildSEDCommonTriangle(std::shared_ptr<Mesh> mesh_ptr_);
	void		combineSEDCommonTriangle(std::shared_ptr<Mesh> mesh_ptr_);
	void		buildSEDCommonTriangle_com(std::shared_ptr<Mesh> mesh_ptr_,int &isCon);
	void		buildSEDCommonTriangle_boundary(std::shared_ptr<Mesh> mesh_ptr_,value_t _Dx,value_t _Dy);
	void		buildSEDCommonTriangle_boundary(std::shared_ptr<Mesh> mesh_ptr_, value_t _Dx, value_t _Dy,value_t _lamda);
	void		buildhSWG(std::shared_ptr<Mesh> mesh_ptr_);
	int         getCommonTriangleNum()const;
	int			gethSWGNum()const;
	void        clear();
	bool        writeAreaData(const std::string& filename, value_t lamda) const;
	void        reportInfo(Qostream& strm) const;
private:
	void        updateCommonTriangle(int nA, int nB, int nC, int nD, int tetID, value_t _eps);
	void        setTriangleTempAt(int index, int _vm, int _fp, value_t _eps);
	void        getTriangleTempAt(int index, int &_v1, int &_v2,int &_v3) const;
	bool        isCommonTriangleTemp(int index) const;
	int         getCommonTriangleTempNum()const;

private:
	std::vector<SCT> m_ctListTemp;
	std::vector<SCT> m_ctList;
	std::vector<SCT> m_btList;
private:
	std::vector<SCT> m_ctList_SED0;
	std::vector<SCT> m_ctList_SED1;
	std::vector<SCT> m_ctList_SED2;
	std::vector<SCT> m_ctList_SED3;
	std::vector<SCT> m_ctList_SED4;
	std::vector<SCT> m_ctList_SED5;
	std::vector<SCT> m_ctList_SED6;
	std::vector<SCT> m_ctList_SED7;
	std::vector<SCT> m_ctList_SED8;
	std::vector<std::vector<SCT>> m_ctList_SED;

private:
	std::vector<h_SWG> m_hSWGList;
	int              totalTriangles;
	int				 bounaryEdges;
	int				 m_ct_size[10];
	int				 common_ct;
};

inline int CommonTriangle::getPositiveTriangle(int index) const
{
	return m_ctList[index].tpid;
}

inline value_t CommonTriangle::getCommonTriangleArea(int index) const
{
	return m_ctList[index].area;
}

inline int CommonTriangle::getbedge(int index) const
{
	return m_ctList[index].bedge;
}

inline void CommonTriangle::getSEDCommonTriangleSize(int index, int &_start,int &_end)const
{
	_start = m_ct_size[index];
	_end = m_ct_size[index + 1];

}

inline void CommonTriangle::getBEDGEandSEDID(int index,int &_bedge,int &_sedid) const
{
	_bedge = m_ctList[index].bedge;
	_sedid = m_ctList[index].sedid;
}

inline void CommonTriangle::getCommonTriangleTet(int index, int & _fp, int & _fm) const
{
	auto& triangle = m_ctList[index];
	_fp = triangle.tpid;
	_fm = triangle.tmid;
}

inline void CommonTriangle::getCommonTriangle(int index, int & _v1, int & _v2,int & _v3) const
{
	_v1 = m_ctList[index].v1;
	_v2 = m_ctList[index].v2;
	_v3 = m_ctList[index].v3;
}

inline void CommonTriangle::getTet_eps(int index, value_t &_eps) const
{
	_eps = m_ctList[index].eps_p;
}

inline const CommonTriangle::SCT & component::CommonTriangle::getSWGRef(size_t i) const
{
	return m_ctList[i];
}

inline bool CommonTriangle::isCommonTriangle(int index) const
{
	return m_ctList[index].tmid != -1;
}

inline int CommonTriangle::getCommonTriangleNum() const
{
	return static_cast<int>(m_ctList.size());
}

inline int CommonTriangle::gethSWGNum() const
{
	return static_cast<int>(m_hSWGList.size());
}

inline void CommonTriangle::clear()
{
	m_ctList.clear();
}

inline bool CommonTriangle::isCommonTriangleTemp(int index) const
{
	return m_ctListTemp[index].tmid != -1;
}

inline int CommonTriangle::getCommonTriangleTempNum() const
{
	return static_cast<int>(m_ctListTemp.size());
}

inline int CommonTriangle::getCommon_ct() const
{
	return common_ct;
}

}





#endif