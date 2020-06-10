#include "CommonTriangle.h"
#include "ErrorHandle.h"
#include "tools.h"

using namespace component;

CommonTriangle::CommonTriangle()
{
}

CommonTriangle::~CommonTriangle()
{
}

void CommonTriangle::getCommonTriangle(int index, int & _vp, int & _vm, int & _fp, int & _fm, value_t & area,value_t &_epsp,value_t &_epsm) const
{
	auto& triangle = m_ctList[index];
	_vp = triangle.vxp;
	_vm = triangle.vxm;
	_fp = triangle.tpid;
	_fm = triangle.tmid;
	area = triangle.area;
	_epsp = triangle.eps_p;
	_epsm = triangle.eps_m;
}

void CommonTriangle::getCommonTriangle(int index, int & _vp, int & _vm, int & _fp, int & _fm, value_t & area, value_t &_epsp, value_t &_epsm,int &_sedid,int &_bedge) const
{
	auto& triangle = m_ctList[index];
	_vp = triangle.vxp;
	_vm = triangle.vxm;
	_fp = triangle.tpid;
	_fm = triangle.tmid;
	area = triangle.area;
	_epsp = triangle.eps_p;
	_epsm = triangle.eps_m;
	_sedid = triangle.sedid;
	_bedge = triangle.bedge;
}

void CommonTriangle::gethSWG(int index, int &_v1, int &_v2, int &_v3, int &_vx, int &_tid, value_t &_eps, value_t &_area) const
{
	auto &temp = m_hSWGList[index];
	_v1 = temp.v1;
	_v2 = temp.v2;
	_v3 = temp.v3;
	_vx = temp.vx;
	_tid = temp.tid;
	_eps = temp.eps;
	_area = temp.area;
}

void CommonTriangle::getCommonTriangle(int index, int &_vp, int &_vm, int &_fp, int &_fm, value_t &area) const
{
	auto& triangle = m_ctList[index];
	_vp = triangle.vxp;
	_vm = triangle.vxm;
	_fp = triangle.tpid;
	_fm = triangle.tmid;
	area = triangle.area;
}

void CommonTriangle::getCommonTriangle(int index, int &_vp, int &_vm, int &_fp, int &_fm, value_t &area,int &_sedid,int &_bedge) const
{
	auto& triangle = m_ctList[index];
	_vp = triangle.vxp;
	_vm = triangle.vxm;
	_fp = triangle.tpid;
	_fm = triangle.tmid;
	area = triangle.area;
	_sedid = triangle.sedid;
	_bedge = triangle.bedge;
}

void CommonTriangle::buildCommonTriangle(std::shared_ptr<Mesh> mesh_ptr_)
{
	auto tetNum = mesh_ptr_->getTetrahedronNum();
	int  nAID, nBID, nCID, nDID;
	value_t eps;
	for (int i = 0; i < tetNum; ++i)
	{
		auto& tet = mesh_ptr_->getTetrahedronRef(i);
		tet.GetVertices(nAID, nBID, nCID, nDID);
		tet.GetEpsilon(eps);
		updateCommonTriangle(nAID, nBID, nCID, nDID, i,eps);
	}

	totalTriangles = static_cast<int>(m_ctListTemp.size());
	for (int i = 0; i < totalTriangles; ++i)
	{
		m_ctList.push_back(std::move(m_ctListTemp[i]));
	}
	m_ctListTemp.clear();
	m_ctListTemp.shrink_to_fit();

	m_ctList.shrink_to_fit();

	VectorR3 vertex1, vertex2, vertex3;
	std::stringstream tmpcelist;
	
	for (auto &ce : m_ctList)
	{
		vertex1 = mesh_ptr_->getVertex(ce.v1);
		vertex2 = mesh_ptr_->getVertex(ce.v2);
		vertex3 = mesh_ptr_->getVertex(ce.v3);
		ce.area = Area(vertex1, vertex2, vertex3);
		tmpcelist << ce.v1+1 << '\t' << ce.v2+1 << '\t' << ce.v3+1 << '\t' << ce.vxp+1 << '\t' << ce.vxm+1 <<'\t'<<ce.tpid+1<<'\t'<<ce.tmid+1<< '\t' << ce.eps_p << '\t' << ce.eps_m << '\n';
	}
	Qofstream tmpfile("temp_celist.dat", std::ios::out);
	tmpfile << tmpcelist.rdbuf();

	Qcout << "It is already output the commonTriangle!" << std::endl;
}

void CommonTriangle::buildCommonTriangleContinuous(std::shared_ptr<Mesh> mesh_ptr_)
{
	auto tetNum = mesh_ptr_->getTetrahedronNum();
	int  nAID, nBID, nCID, nDID;
	value_t eps;
	for (int i = 0; i < tetNum; ++i)
	{
		auto& tet = mesh_ptr_->getTetrahedronRef(i);
		tet.GetVertices(nAID, nBID, nCID, nDID);
		tet.GetEpsilon(eps);
		updateCommonTriangle(nAID, nBID, nCID, nDID, i, eps);
	}

	totalTriangles = static_cast<int>(m_ctListTemp.size());
	for (int i = 0; i < totalTriangles; ++i)
	{
		if (!isCommonTriangleTemp(i))
		{
			m_btList.push_back(std::move(m_ctListTemp[i]));
			continue;
		}
		m_ctList.push_back(std::move(m_ctListTemp[i]));
	}

	m_ctListTemp.clear();
	m_ctListTemp.shrink_to_fit();

	/////////////////////////////////////////
	bounaryEdges = m_btList.size();
	VectorR3 minBox, maxBox;
	mesh_ptr_->getBoundary(minBox, maxBox);
	std::vector<SCT> maxEdges_x, maxEdges_y;
	std::vector<SCT> minEdges_x, minEdges_y;

	for (int i = 0; i < bounaryEdges; i++)
	{
		VectorR3 v1 = mesh_ptr_->getVertex(m_btList[i].v1);
		VectorR3 v2 = mesh_ptr_->getVertex(m_btList[i].v2);
		VectorR3 v3 = mesh_ptr_->getVertex(m_btList[i].v3);
		VectorR3 v_cen = (v1 + v2 + v3) / 3.0f;

		if (abs(v_cen.x - minBox.x) < 1e-6)
			minEdges_x.push_back(std::move(m_btList[i]));
		else if (abs(v_cen.x - maxBox.x) < 1e-6)
			maxEdges_x.push_back(std::move(m_btList[i]));
		else if (abs(v_cen.y - minBox.y) < 1e-6)
			minEdges_y.push_back(std::move(m_btList[i]));
		else if (abs(v_cen.y - maxBox.y) < 1e-6)
			maxEdges_y.push_back(std::move(m_btList[i]));
		else
			m_ctList.push_back(std::move(m_btList[i]));
	}

	if ((minEdges_x.size() == maxEdges_x.size()) && (minEdges_y.size() == maxEdges_y.size()))
		Qcout << "The boundary number is match..." << std::endl;
	else
		Qcout << "The boundary number isn't match...." << std::endl;

	int bnum1 = maxEdges_x.size();
	int bnum2 = maxEdges_y.size();
	for (int i = 0; i < bnum1; ++i)
	{
		VectorR3 vmax1 = mesh_ptr_->getVertex(maxEdges_x[i].v1);
		VectorR3 vmax2 = mesh_ptr_->getVertex(maxEdges_x[i].v2);
		VectorR3 vmax3 = mesh_ptr_->getVertex(maxEdges_x[i].v3);
		for (int j = 0; j < bnum1; ++j)
		{
			VectorR3 vmin1 = mesh_ptr_->getVertex(minEdges_x[j].v1);
			VectorR3 vmin2 = mesh_ptr_->getVertex(minEdges_x[j].v2);
			VectorR3 vmin3 = mesh_ptr_->getVertex(minEdges_x[j].v3);
			if ((abs(vmax1.y + vmax2.y + vmax3.y - vmin1.y - vmin2.y - vmin3.y) < 1e-6)
				&& (abs(vmax1.z + vmax2.z + vmax3.z - vmin1.z - vmin2.z - vmin3.z) < 1e-6))
			{
				maxEdges_x[i].vxm = minEdges_x[j].vxp;
				maxEdges_x[i].tmid = minEdges_x[j].tpid;
				maxEdges_x[i].eps_m = minEdges_x[j].eps_p;
				maxEdges_x[i].bedge = 1;
				//maxEdges_x[i].bedge = 2;
				//maxEdges_x[i].length = (vmax1 - vmax2).Norm();
				m_ctList.push_back(std::move(maxEdges_x[i]));
			}
		}
	}

	for (int i = 0; i < bnum2; ++i)
	{
		VectorR3 vmax1 = mesh_ptr_->getVertex(maxEdges_y[i].v1);
		VectorR3 vmax2 = mesh_ptr_->getVertex(maxEdges_y[i].v2);
		VectorR3 vmax3 = mesh_ptr_->getVertex(maxEdges_y[i].v3);
		for (int j = 0; j < bnum2; ++j)
		{
			VectorR3 vmin1 = mesh_ptr_->getVertex(minEdges_y[j].v1);
			VectorR3 vmin2 = mesh_ptr_->getVertex(minEdges_y[j].v2);
			VectorR3 vmin3 = mesh_ptr_->getVertex(minEdges_y[j].v3);
			if ((abs(vmax1.x + vmax2.x + vmax3.x - vmin1.x - vmin2.x - vmin3.x) < 1e-6)
				&& (abs(vmax1.z + vmax2.z + vmax3.z - vmin1.z - vmin2.z - vmin3.z) < 1e-6))
			{
				maxEdges_y[i].vxm = minEdges_y[j].vxp;
				maxEdges_y[i].tmid = minEdges_y[j].tpid;
				maxEdges_y[i].eps_m = minEdges_y[j].eps_p;
				maxEdges_y[i].bedge = 2;
				//maxEdges_x[i].bedge = 2;
				//maxEdges_x[i].length = (vmax1 - vmax2).Norm();
				m_ctList.push_back(std::move(maxEdges_y[i]));
			}
		}
	}

	m_ctList.shrink_to_fit();

	VectorR3 vertex1, vertex2, vertex3;
	std::stringstream tmpcelist;

	for (auto &ce : m_ctList)
	{
		vertex1 = mesh_ptr_->getVertex(ce.v1);
		vertex2 = mesh_ptr_->getVertex(ce.v2);
		vertex3 = mesh_ptr_->getVertex(ce.v3);
		ce.area = Area(vertex1, vertex2, vertex3);
		tmpcelist << ce.v1 + 1 << '\t' << ce.v2 + 1 << '\t' << ce.v3 + 1 << '\t' << ce.vxp + 1 << '\t' << ce.vxm + 1 << '\t' << ce.tpid + 1 << '\t' << ce.tmid + 1 << '\t' << ce.eps_p << '\t' << ce.eps_m << '\n';
	}
	Qofstream tmpfile("temp_celist.dat", std::ios::out);
	tmpfile << tmpcelist.rdbuf();

	Qcout << "It is already output the commonTriangle!" << std::endl;

}

void CommonTriangle::buildSEDCommonTriangle(std::shared_ptr<Mesh> mesh_ptr_)
{
}

void CommonTriangle::combineSEDCommonTriangle(std::shared_ptr<Mesh> mesh_ptr_)
{
	int unk_SED;

	unk_SED = m_ctList_SED0.size();
	m_ct_size[0] = 0;
	for (int i = 0; i < unk_SED; i++)
	{
		m_ctList_SED0[i].sedid = 0;
		m_ctList.push_back(std::move(m_ctList_SED0[i]));
	}
	m_ct_size[1] = m_ctList.size();

	unk_SED = m_ctList_SED1.size();
	for (int i = 0; i < unk_SED; i++)
	{
		m_ctList_SED1[i].sedid = 1;
		m_ctList.push_back(std::move(m_ctList_SED1[i]));
	}
	m_ct_size[2] = m_ctList.size();

	unk_SED = m_ctList_SED2.size();
	for (int i = 0; i < unk_SED; i++)
	{
		m_ctList_SED2[i].sedid = 2;
		m_ctList.push_back(std::move(m_ctList_SED2[i]));
	}
	m_ct_size[3] = m_ctList.size();

	unk_SED = m_ctList_SED3.size();
	for (int i = 0; i < unk_SED; i++)
	{
		m_ctList_SED3[i].sedid = 3;
		m_ctList.push_back(std::move(m_ctList_SED3[i]));
	}
	m_ct_size[4] = m_ctList.size();

	unk_SED = m_ctList_SED4.size();
	for (int i = 0; i < unk_SED; i++)
	{
		m_ctList_SED4[i].sedid = 4;
		m_ctList.push_back(std::move(m_ctList_SED4[i]));
	}
	m_ct_size[5] = m_ctList.size();

	unk_SED = m_ctList_SED5.size();
	for (int i = 0; i < unk_SED; i++)
	{
		m_ctList_SED5[i].sedid = 5;
		m_ctList.push_back(std::move(m_ctList_SED5[i]));
	}
	m_ct_size[6] = m_ctList.size();

	unk_SED = m_ctList_SED6.size();
	for (int i = 0; i < unk_SED; i++)
	{
		m_ctList_SED6[i].sedid = 6;
		m_ctList.push_back(std::move(m_ctList_SED6[i]));
	}
	m_ct_size[7] = m_ctList.size();

	unk_SED = m_ctList_SED7.size();
	for (int i = 0; i < unk_SED; i++)
	{
		m_ctList_SED7[i].sedid = 7;
		m_ctList.push_back(std::move(m_ctList_SED7[i]));
	}
	m_ct_size[8] = m_ctList.size();

	unk_SED = m_ctList_SED8.size();
	for (int i = 0; i < unk_SED; i++)
	{
		m_ctList_SED8[i].sedid = 8;
		m_ctList.push_back(std::move(m_ctList_SED8[i]));
	}
	m_ct_size[9] = m_ctList.size();

	VectorR3 vertex1, vertex2, vertex3;
	std::stringstream tmpctlist;

	for (auto &ct : m_ctList)
	{
		vertex1 = mesh_ptr_->getVertex(ct.v1);
		vertex2 = mesh_ptr_->getVertex(ct.v2);
		vertex3 = mesh_ptr_->getVertex(ct.v3);
		ct.area = Area(vertex1, vertex2, vertex3);
		tmpctlist << ct.v1 + 1 << '\t' << ct.v2 + 1 << '\t' << ct.v3 + 1 << '\t' << ct.vxp + 1 << '\t' << ct.vxm + 1 << '\t' << ct.tpid + 1 << '\t' << ct.tmid + 1 << '\t' << ct.eps_p << '\t' << ct.eps_m << '\t' << ct.bedge << '\t' << ct.sedid << '\n';
	}
	Qofstream tmpfile("temp_ctlist.dat", std::ios::out);
	tmpfile << tmpctlist.rdbuf();

	Qcout << "It is already output the commonTriangle!" << std::endl;
}

void component::CommonTriangle::buildSEDCommonTriangle_com(std::shared_ptr<Mesh> mesh_ptr_,int& isCon)
{
	auto tetNum = mesh_ptr_->getTetrahedronNum();
	int  nAID, nBID, nCID, nDID;
	value_t eps;
	for (int i = 0; i < tetNum; ++i)
	{
		auto& tet = mesh_ptr_->getTetrahedronRef(i);
		tet.GetVertices(nAID, nBID, nCID, nDID);
		tet.GetEpsilon(eps);
		updateCommonTriangle(nAID, nBID, nCID, nDID, i, eps);
	}

	totalTriangles = static_cast<int>(m_ctListTemp.size());
	if (isCon == 0)
	{
		for (int i = 0; i < totalTriangles; ++i)
		{
			m_ctList_SED0.push_back(std::move(m_ctListTemp[i]));
			m_ctList_SED1.push_back(std::move(m_ctListTemp[i]));
			m_ctList_SED2.push_back(std::move(m_ctListTemp[i]));
			m_ctList_SED3.push_back(std::move(m_ctListTemp[i]));
			m_ctList_SED4.push_back(std::move(m_ctListTemp[i]));
			m_ctList_SED5.push_back(std::move(m_ctListTemp[i]));
			m_ctList_SED6.push_back(std::move(m_ctListTemp[i]));
			m_ctList_SED7.push_back(std::move(m_ctListTemp[i]));
			m_ctList_SED8.push_back(std::move(m_ctListTemp[i]));
		}
	}
	else
	{
		for (int i = 0; i < totalTriangles; ++i)
		{
			if (!isCommonTriangleTemp(i))
			{
				m_btList.push_back(std::move(m_ctListTemp[i]));
				continue;
			}
			m_ctList_SED0.push_back(std::move(m_ctListTemp[i]));
			m_ctList_SED1.push_back(std::move(m_ctListTemp[i]));
			m_ctList_SED2.push_back(std::move(m_ctListTemp[i]));
			m_ctList_SED3.push_back(std::move(m_ctListTemp[i]));
			m_ctList_SED4.push_back(std::move(m_ctListTemp[i]));
			m_ctList_SED5.push_back(std::move(m_ctListTemp[i]));
			m_ctList_SED6.push_back(std::move(m_ctListTemp[i]));
			m_ctList_SED7.push_back(std::move(m_ctListTemp[i]));
			m_ctList_SED8.push_back(std::move(m_ctListTemp[i]));
		}
	}

	//m_ctListTemp.clear();
	m_ctListTemp.shrink_to_fit();
	common_ct = m_ctList_SED0.size();
}

void component::CommonTriangle::buildSEDCommonTriangle_boundary(std::shared_ptr<Mesh> mesh_ptr_, value_t _Dx,value_t _Dy)
{
	bounaryEdges = m_btList.size();
	VectorR3 minBox, maxBox;
	mesh_ptr_->getBoundary(minBox, maxBox);
	value_t dis_x, dis_y, dis_z;
	mesh_ptr_->getSize(dis_x, dis_y, dis_z);
	//maxBox.x = minBox.x + _Dx;
	//maxBox.y = minBox.y + _Dy;
	std::vector<SCT> maxEdges_x, maxEdges_y;
	std::vector<SCT> minEdges_x, minEdges_y;

	for (int i = 0; i < bounaryEdges; i++)
	{
		VectorR3 v1 = mesh_ptr_->getVertex(m_btList[i].v1);
		VectorR3 v2 = mesh_ptr_->getVertex(m_btList[i].v2);
		VectorR3 v3 = mesh_ptr_->getVertex(m_btList[i].v3);
		VectorR3 v_cen = (v1 + v2 + v3) / 3.0f;

		if ((abs(v_cen.x - minBox.x) < 1e-6) && (abs(_Dx - dis_x)<1e-6))
		{
			
		}
		else if ((abs(v_cen.x - maxBox.x) < 1e-6) && (abs(_Dx - dis_x)<1e-6))
		{
			
		}
		else if ((abs(v_cen.y - minBox.y) < 1e-6) && (abs(_Dy - dis_y)<1e-6))
		{
			
		}
		else if ((abs(v_cen.y - maxBox.y) < 1e-6) && (abs(_Dy - dis_y)<1e-6))
		{
			
		}
		else
		{
			m_ctList_SED0.push_back(std::move(m_btList[i]));
			m_ctList_SED1.push_back(std::move(m_btList[i]));
			m_ctList_SED2.push_back(std::move(m_btList[i]));
			m_ctList_SED3.push_back(std::move(m_btList[i]));
			m_ctList_SED4.push_back(std::move(m_btList[i]));
			m_ctList_SED5.push_back(std::move(m_btList[i]));
			m_ctList_SED6.push_back(std::move(m_btList[i]));
			m_ctList_SED7.push_back(std::move(m_btList[i]));
			m_ctList_SED8.push_back(std::move(m_btList[i]));
		}
	}
	common_ct = m_ctList_SED0.size();

	for (int i = 0; i < bounaryEdges; i++)
	{
		VectorR3 v1 = mesh_ptr_->getVertex(m_btList[i].v1);
		VectorR3 v2 = mesh_ptr_->getVertex(m_btList[i].v2);
		VectorR3 v3 = mesh_ptr_->getVertex(m_btList[i].v3);
		VectorR3 v_cen = (v1 + v2 + v3) / 3.0f;

		if ((abs(v_cen.x - minBox.x) < 1e-6)&&(abs(_Dx-dis_x)<1e-6))
		{
			minEdges_x.push_back(std::move(m_btList[i]));
			m_ctList_SED0.push_back(std::move(m_btList[i]));
			m_ctList_SED1.push_back(std::move(m_btList[i]));
			m_ctList_SED2.push_back(std::move(m_btList[i]));
		}
		else if ((abs(v_cen.x - maxBox.x) < 1e-6)&& (abs(_Dx - dis_x)<1e-6))
		{
			maxEdges_x.push_back(std::move(m_btList[i]));
			m_ctList_SED6.push_back(std::move(m_btList[i]));
			m_ctList_SED7.push_back(std::move(m_btList[i]));
			m_ctList_SED8.push_back(std::move(m_btList[i]));
		}
		else if ((abs(v_cen.y - minBox.y) < 1e-6)&&(abs(_Dy-dis_y)<1e-6))
		{
			minEdges_y.push_back(std::move(m_btList[i]));
			m_ctList_SED0.push_back(std::move(m_btList[i]));
			m_ctList_SED3.push_back(std::move(m_btList[i]));
			m_ctList_SED6.push_back(std::move(m_btList[i]));
		}
		else if ((abs(v_cen.y - maxBox.y) < 1e-6)&& (abs(_Dy - dis_y)<1e-6))
		{
			maxEdges_y.push_back(std::move(m_btList[i]));
			m_ctList_SED2.push_back(std::move(m_btList[i]));
			m_ctList_SED5.push_back(std::move(m_btList[i]));
			m_ctList_SED8.push_back(std::move(m_btList[i]));
		}
		else
		{
			/*m_ctList_SED0.push_back(std::move(m_btList[i]));
			m_ctList_SED1.push_back(std::move(m_btList[i]));
			m_ctList_SED2.push_back(std::move(m_btList[i]));
			m_ctList_SED3.push_back(std::move(m_btList[i]));
			m_ctList_SED4.push_back(std::move(m_btList[i]));
			m_ctList_SED5.push_back(std::move(m_btList[i]));
			m_ctList_SED6.push_back(std::move(m_btList[i]));
			m_ctList_SED7.push_back(std::move(m_btList[i]));
			m_ctList_SED8.push_back(std::move(m_btList[i]));*/
		}
	}

	if ((minEdges_x.size() == maxEdges_x.size()) && (minEdges_y.size() == maxEdges_y.size()))
		Qcout << "The boundary number is match..." << std::endl;
	else
		Qcout << "The boundary number isn't match...." << std::endl;

	int bnum1 = maxEdges_x.size();
	int bnum2 = maxEdges_y.size();
	for (int i = 0; i < bnum1; ++i)
	{
		VectorR3 vmax1 = mesh_ptr_->getVertex(maxEdges_x[i].v1);
		VectorR3 vmax2 = mesh_ptr_->getVertex(maxEdges_x[i].v2);
		VectorR3 vmax3 = mesh_ptr_->getVertex(maxEdges_x[i].v3);
		for (int j = 0; j < bnum1; ++j)
		{
			VectorR3 vmin1 = mesh_ptr_->getVertex(minEdges_x[j].v1);
			VectorR3 vmin2 = mesh_ptr_->getVertex(minEdges_x[j].v2);
			VectorR3 vmin3 = mesh_ptr_->getVertex(minEdges_x[j].v3);
			if ((abs(vmax1.y + vmax2.y + vmax3.y - vmin1.y - vmin2.y - vmin3.y) < 1e-6)
				&& (abs(vmax1.z + vmax2.z + vmax3.z - vmin1.z - vmin2.z - vmin3.z) < 1e-6))
			{
				maxEdges_x[i].vxm = minEdges_x[j].vxp;
				maxEdges_x[i].tmid = minEdges_x[j].tpid;
				maxEdges_x[i].eps_m = minEdges_x[j].eps_p;
				maxEdges_x[i].bedge = 1;

				m_ctList_SED0.push_back(std::move(maxEdges_x[i]));
				m_ctList_SED1.push_back(std::move(maxEdges_x[i]));
				m_ctList_SED2.push_back(std::move(maxEdges_x[i]));
				m_ctList_SED3.push_back(std::move(maxEdges_x[i]));
				m_ctList_SED4.push_back(std::move(maxEdges_x[i]));
				m_ctList_SED5.push_back(std::move(maxEdges_x[i]));
			}
		}
	}

	for (int i = 0; i < bnum2; ++i)
	{
		VectorR3 vmax1 = mesh_ptr_->getVertex(maxEdges_y[i].v1);
		VectorR3 vmax2 = mesh_ptr_->getVertex(maxEdges_y[i].v2);
		VectorR3 vmax3 = mesh_ptr_->getVertex(maxEdges_y[i].v3);
		for (int j = 0; j < bnum2; ++j)
		{
			VectorR3 vmin1 = mesh_ptr_->getVertex(minEdges_y[j].v1);
			VectorR3 vmin2 = mesh_ptr_->getVertex(minEdges_y[j].v2);
			VectorR3 vmin3 = mesh_ptr_->getVertex(minEdges_y[j].v3);
			if ((abs(vmax1.x + vmax2.x + vmax3.x - vmin1.x - vmin2.x - vmin3.x) < 1e-6)
				&& (abs(vmax1.z + vmax2.z + vmax3.z - vmin1.z - vmin2.z - vmin3.z) < 1e-6))
			{
				maxEdges_y[i].vxm = minEdges_y[j].vxp;
				maxEdges_y[i].tmid = minEdges_y[j].tpid;
				maxEdges_y[i].eps_m = minEdges_y[j].eps_p;
				maxEdges_y[i].bedge = 2;

				m_ctList_SED0.push_back(std::move(maxEdges_y[i]));
				m_ctList_SED1.push_back(std::move(maxEdges_y[i]));
				m_ctList_SED3.push_back(std::move(maxEdges_y[i]));
				m_ctList_SED4.push_back(std::move(maxEdges_y[i]));
				m_ctList_SED6.push_back(std::move(maxEdges_y[i]));
				m_ctList_SED7.push_back(std::move(maxEdges_y[i]));
			}
		}
	}
}

void component::CommonTriangle::buildSEDCommonTriangle_boundary(std::shared_ptr<Mesh> mesh_ptr_, value_t _Dx, value_t _Dy,value_t _lamda)
{
	bounaryEdges = m_btList.size();
	VectorR3 minBox, maxBox;
	mesh_ptr_->getBoundary(minBox, maxBox);
	value_t dis_x, dis_y, dis_z;
	mesh_ptr_->getSize(dis_x, dis_y, dis_z);
	//maxBox.x = minBox.x + _Dx;
	//maxBox.y = minBox.y + _Dy;
	std::vector<SCT> maxEdges_x, maxEdges_y;
	std::vector<SCT> minEdges_x, minEdges_y;

	for (int i = 0; i < bounaryEdges; i++)
	{
		VectorR3 v1 = mesh_ptr_->getVertex(m_btList[i].v1);
		VectorR3 v2 = mesh_ptr_->getVertex(m_btList[i].v2);
		VectorR3 v3 = mesh_ptr_->getVertex(m_btList[i].v3);
		VectorR3 v_cen = (v1 + v2 + v3) / 3.0f;

		if ((abs(v_cen.x - minBox.x) < 1e-3*_lamda) && (abs(_Dx - dis_x)<1e-3*_lamda))
		{

		}
		else if ((abs(v_cen.x - maxBox.x) < 1e-3*_lamda) && (abs(_Dx - dis_x)<1e-3*_lamda))
		{

		}
		else if ((abs(v_cen.y - minBox.y) < 1e-3*_lamda) && (abs(_Dy - dis_y)<1e-3*_lamda))
		{

		}
		else if ((abs(v_cen.y - maxBox.y) < 1e-3*_lamda) && (abs(_Dy - dis_y)<1e-3*_lamda))
		{

		}
		else
		{
			m_ctList_SED0.push_back(std::move(m_btList[i]));
			m_ctList_SED1.push_back(std::move(m_btList[i]));
			m_ctList_SED2.push_back(std::move(m_btList[i]));
			m_ctList_SED3.push_back(std::move(m_btList[i]));
			m_ctList_SED4.push_back(std::move(m_btList[i]));
			m_ctList_SED5.push_back(std::move(m_btList[i]));
			m_ctList_SED6.push_back(std::move(m_btList[i]));
			m_ctList_SED7.push_back(std::move(m_btList[i]));
			m_ctList_SED8.push_back(std::move(m_btList[i]));
		}
	}
	common_ct = m_ctList_SED0.size();

	for (int i = 0; i < bounaryEdges; i++)
	{
		VectorR3 v1 = mesh_ptr_->getVertex(m_btList[i].v1);
		VectorR3 v2 = mesh_ptr_->getVertex(m_btList[i].v2);
		VectorR3 v3 = mesh_ptr_->getVertex(m_btList[i].v3);
		VectorR3 v_cen = (v1 + v2 + v3) / 3.0f;

		if ((abs(v_cen.x - minBox.x) < 1e-3*_lamda) && (abs(_Dx - dis_x)<1e-3*_lamda))
		{
			minEdges_x.push_back(std::move(m_btList[i]));
			m_ctList_SED0.push_back(std::move(m_btList[i]));
			m_ctList_SED1.push_back(std::move(m_btList[i]));
			m_ctList_SED2.push_back(std::move(m_btList[i]));
		}
		else if ((abs(v_cen.x - maxBox.x) < 1e-3*_lamda) && (abs(_Dx - dis_x)<1e-3*_lamda))
		{
			maxEdges_x.push_back(std::move(m_btList[i]));
			m_ctList_SED6.push_back(std::move(m_btList[i]));
			m_ctList_SED7.push_back(std::move(m_btList[i]));
			m_ctList_SED8.push_back(std::move(m_btList[i]));
		}
		else if ((abs(v_cen.y - minBox.y) < 1e-3*_lamda) && (abs(_Dy - dis_y)<1e-3*_lamda))
		{
			minEdges_y.push_back(std::move(m_btList[i]));
			m_ctList_SED0.push_back(std::move(m_btList[i]));
			m_ctList_SED3.push_back(std::move(m_btList[i]));
			m_ctList_SED6.push_back(std::move(m_btList[i]));
		}
		else if ((abs(v_cen.y - maxBox.y) < 1e-3*_lamda) && (abs(_Dy - dis_y)<1e-3*_lamda))
		{
			maxEdges_y.push_back(std::move(m_btList[i]));
			m_ctList_SED2.push_back(std::move(m_btList[i]));
			m_ctList_SED5.push_back(std::move(m_btList[i]));
			m_ctList_SED8.push_back(std::move(m_btList[i]));
		}
		else
		{
			/*m_ctList_SED0.push_back(std::move(m_btList[i]));
			m_ctList_SED1.push_back(std::move(m_btList[i]));
			m_ctList_SED2.push_back(std::move(m_btList[i]));
			m_ctList_SED3.push_back(std::move(m_btList[i]));
			m_ctList_SED4.push_back(std::move(m_btList[i]));
			m_ctList_SED5.push_back(std::move(m_btList[i]));
			m_ctList_SED6.push_back(std::move(m_btList[i]));
			m_ctList_SED7.push_back(std::move(m_btList[i]));
			m_ctList_SED8.push_back(std::move(m_btList[i]));*/
		}
	}

	if ((minEdges_x.size() == maxEdges_x.size()) && (minEdges_y.size() == maxEdges_y.size()))
		Qcout << "The boundary number is match..." << std::endl;
	else
		Qcout << "The boundary number isn't match...." << std::endl;

	int bnum1 = maxEdges_x.size();
	int bnum2 = maxEdges_y.size();
	for (int i = 0; i < bnum1; ++i)
	{
		VectorR3 vmax1 = mesh_ptr_->getVertex(maxEdges_x[i].v1);
		VectorR3 vmax2 = mesh_ptr_->getVertex(maxEdges_x[i].v2);
		VectorR3 vmax3 = mesh_ptr_->getVertex(maxEdges_x[i].v3);
		for (int j = 0; j < bnum1; ++j)
		{
			VectorR3 vmin1 = mesh_ptr_->getVertex(minEdges_x[j].v1);
			VectorR3 vmin2 = mesh_ptr_->getVertex(minEdges_x[j].v2);
			VectorR3 vmin3 = mesh_ptr_->getVertex(minEdges_x[j].v3);
			if ((abs(vmax1.y + vmax2.y + vmax3.y - vmin1.y - vmin2.y - vmin3.y) < 1e-3*_lamda)
				&& (abs(vmax1.z + vmax2.z + vmax3.z - vmin1.z - vmin2.z - vmin3.z) < 1e-3*_lamda))
			{
				maxEdges_x[i].vxm = minEdges_x[j].vxp;
				maxEdges_x[i].tmid = minEdges_x[j].tpid;
				maxEdges_x[i].eps_m = minEdges_x[j].eps_p;
				maxEdges_x[i].bedge = 1;

				m_ctList_SED0.push_back(std::move(maxEdges_x[i]));
				m_ctList_SED1.push_back(std::move(maxEdges_x[i]));
				m_ctList_SED2.push_back(std::move(maxEdges_x[i]));
				m_ctList_SED3.push_back(std::move(maxEdges_x[i]));
				m_ctList_SED4.push_back(std::move(maxEdges_x[i]));
				m_ctList_SED5.push_back(std::move(maxEdges_x[i]));
			}
		}
	}

	for (int i = 0; i < bnum2; ++i)
	{
		VectorR3 vmax1 = mesh_ptr_->getVertex(maxEdges_y[i].v1);
		VectorR3 vmax2 = mesh_ptr_->getVertex(maxEdges_y[i].v2);
		VectorR3 vmax3 = mesh_ptr_->getVertex(maxEdges_y[i].v3);
		for (int j = 0; j < bnum2; ++j)
		{
			VectorR3 vmin1 = mesh_ptr_->getVertex(minEdges_y[j].v1);
			VectorR3 vmin2 = mesh_ptr_->getVertex(minEdges_y[j].v2);
			VectorR3 vmin3 = mesh_ptr_->getVertex(minEdges_y[j].v3);
			if ((abs(vmax1.x + vmax2.x + vmax3.x - vmin1.x - vmin2.x - vmin3.x) < 1e-3*_lamda)
				&& (abs(vmax1.z + vmax2.z + vmax3.z - vmin1.z - vmin2.z - vmin3.z) < 1e-3*_lamda))
			{
				maxEdges_y[i].vxm = minEdges_y[j].vxp;
				maxEdges_y[i].tmid = minEdges_y[j].tpid;
				maxEdges_y[i].eps_m = minEdges_y[j].eps_p;
				maxEdges_y[i].bedge = 2;

				m_ctList_SED0.push_back(std::move(maxEdges_y[i]));
				m_ctList_SED1.push_back(std::move(maxEdges_y[i]));
				m_ctList_SED3.push_back(std::move(maxEdges_y[i]));
				m_ctList_SED4.push_back(std::move(maxEdges_y[i]));
				m_ctList_SED6.push_back(std::move(maxEdges_y[i]));
				m_ctList_SED7.push_back(std::move(maxEdges_y[i]));
			}
		}
	}
}

void CommonTriangle::buildhSWG(std::shared_ptr<Mesh> model)
{
	auto tetNum = model->getTetrahedronNum();
	int nA, nB, nC, nD;
	value_t eps;
	for (int i = 0; i < tetNum; i++)
	{
		auto& tet = model->getTetrahedronRef(i);
		tet.GetVertices(nA, nB, nC, nD);
		tet.GetEpsilon(eps);
		m_hSWGList.push_back(h_SWG(nA, nC, nB, nD, i, eps));
		m_hSWGList.push_back(h_SWG(nA, nD, nC, nB, i, eps));
		m_hSWGList.push_back(h_SWG(nA, nB, nD, nC, i, eps));
		m_hSWGList.push_back(h_SWG(nB, nC, nD, nA, i, eps));
	}

	VectorR3 vertex1, vertex2, vertex3;
	std::stringstream tmpcelist;

	for (auto &ce : m_hSWGList)
	{
		vertex1 = model->getVertex(ce.v1);
		vertex2 = model->getVertex(ce.v2);
		vertex3 = model->getVertex(ce.v3);
		ce.area = Area(vertex1, vertex2, vertex3);
		tmpcelist << ce.v1 + 1 << '\t' << ce.v2 + 1 << '\t' << ce.v3 + 1 << '\t' << ce.vx + 1 << '\t' << ce.tid + 1 << '\t' << ce.eps << '\n';
	}
	Qofstream tmpfile("temp_celist.dat", std::ios::out);
	tmpfile << tmpcelist.rdbuf();

	Qcout << "It is already output the commonTriangle!" << std::endl;
}

bool CommonTriangle::writeAreaData(const std::string & filename, value_t lamda) const
{
	Qofstream output(filename, std::ios::out);
	if (output.fail())
		throw FileError("fail opening file: " + filename);

	auto trianglesNum = m_ctList.size();

	output << std::fixed << std::setprecision(4);
	output << "area data: " << output.widen('\n');
	for (size_t i = 0; i < trianglesNum; ++i)
	{
		output << std::setw(10) << i + 1 << std::setw(15) << m_ctList[i].area / lamda << output.widen('\n');
	}

	return true;
}

void CommonTriangle::reportInfo(Qostream & strm) const
{
	size_t mem = m_ctList.size() * m_ctList.size() * 2 * sizeof(value_t) / (1024 * 1024);
	size_t cemem = sizeof(SCT) * m_ctList.size() / (1024 * 1024);
	auto oldState = strm.flags();
	strm << HEADING "CommonTriangle Information\n" TRAILING;
	strm << LEVEL1 "Total triangles: " << totalTriangles << strm.widen('\n')
		<< LEVEL1 "Total SWGs: " << m_ctList.size() << strm.widen('\n')
		<< LEVEL1 "Memory prediction: " << mem << " MB" << FORMAT_MEMORY(mem) << '\n'
		<< LEVEL1 "Memory detail:" << '\n'
		<< LEVEL2 "SWG:" << cemem << " MB" << FORMAT_MEMORY(cemem) << '\n';
	strm.flush();
	strm.flags(oldState);
}

void CommonTriangle::updateCommonTriangle(int nA, int nB, int nC, int nD, int tetID, value_t _eps)
{
	bool newTriangle1 = true, newTriangle2 = true, newTriangle3 = true, newTriangle4 = true;
	int  v1ID, v2ID, v3ID;

	int triNum = static_cast<int>(m_ctListTemp.size());
	for (int i = 0; i < triNum; ++i)
	{
		if (isCommonTriangleTemp(i))
			continue;
		getTriangleTempAt(i, v1ID, v2ID, v3ID);
		if ((nA == v1ID&&nB == v2ID&&nC == v3ID) || (nA == v1ID&&nC == v2ID&&nB == v3ID)|| (nB == v1ID&&nA == v2ID&&nC == v3ID) || (nB == v1ID&&nC == v2ID&&nA == v3ID)|| (nC == v1ID&&nA == v2ID&&nB == v3ID) || (nC == v1ID&&nB == v2ID&&nA == v3ID))
		{
			newTriangle1 = false;
			setTriangleTempAt(i, nD, tetID, _eps);
		}
		else if ((nB == v1ID&&nC == v2ID&&nD == v3ID) || (nB == v1ID&&nD == v2ID&&nC == v3ID)|| (nC == v1ID&&nB == v2ID&&nD == v3ID) || (nC == v1ID&&nD == v2ID&&nB == v3ID)|| (nD == v1ID&&nB == v2ID&&nC == v3ID) || (nD == v1ID&&nC == v2ID&&nB == v3ID))
		{
			newTriangle2 = false;
			setTriangleTempAt(i, nA, tetID, _eps);
		}
		else if ((nA == v1ID&&nC == v2ID&&nD == v3ID) || (nA == v1ID&&nD == v2ID&&nC == v3ID)|| (nC == v1ID&&nA == v2ID&&nD == v3ID) || (nC == v1ID&&nD == v2ID&&nA == v3ID)|| (nD == v1ID&&nA == v2ID&&nC == v3ID) || (nD == v1ID&&nC == v2ID&&nA == v3ID))
		{
			newTriangle3 = false;
			setTriangleTempAt(i, nB, tetID, _eps);
		}
		else if ((nA == v1ID&&nB == v2ID&&nD == v3ID) || (nA == v1ID&&nD == v2ID&&nB == v3ID)|| (nB == v1ID&&nA == v2ID&&nD == v3ID) || (nB == v1ID&&nD == v2ID&&nA == v3ID)|| (nD == v1ID&&nA == v2ID&&nB == v3ID) || (nD == v1ID&&nB == v2ID&&nA == v3ID))
		{
			newTriangle4 = false;
			setTriangleTempAt(i, nC, tetID, _eps);
		}
	}
	if (newTriangle1)
	{
		m_ctListTemp.push_back(SCT(nA, nC, nB, nD, tetID, _eps));
	}
	if (newTriangle2)
	{
		m_ctListTemp.push_back(SCT(nB, nC, nD, nA, tetID, _eps));
	}
	if (newTriangle3)
	{
		m_ctListTemp.push_back(SCT(nA, nD, nC, nB, tetID, _eps));
	}
	if (newTriangle4)
	{
		m_ctListTemp.push_back(SCT(nA, nB, nD, nC, tetID, _eps));
	}
}

void CommonTriangle::setTriangleTempAt(int index, int _vm, int _fp, value_t _eps)
{
	m_ctListTemp[index].vxm = _vm;
	m_ctListTemp[index].tmid = _fp;
	m_ctListTemp[index].eps_m = _eps;
}

void CommonTriangle::getTriangleTempAt(int index, int & _v1, int & _v2,int & _v3) const
{
	_v1 = m_ctListTemp[index].v1;
	_v2 = m_ctListTemp[index].v2;
	_v3 = m_ctListTemp[index].v3;
}


