#include "CommonEdge.h"
#include "ErrorHandle.h"
#include "tools.h"
#include "Mathdef.h"

using namespace component;

CommonEdge::CommonEdge()
{
}

CommonEdge::~CommonEdge()
{
}

void CommonEdge::getCommonEdge(int index, int & _vp, int & _vm, int & _fp, int & _fm, value_t & _len) const
{
	auto& edge = m_ceList[index];
	_vp = edge.vxp;
	_vm = edge.vxm;
	_fp = edge.tpid;
	_fm = edge.tmid;
	_len = edge.length;
}

void CommonEdge::getCommonEdge(int index, int & _vp, int & _vm, int & _fp, int & _fm, value_t & _len,int &_bedge) const
{
    auto& edge = m_ceList[index];
    _vp = edge.vxp;
    _vm = edge.vxm;
    _fp = edge.tpid;
    _fm = edge.tmid;
    _len = edge.length;
	_bedge = edge.bedge;
}

void CommonEdge::getCommonEdge(int index, int & _vp, int & _vm, int & _fp, int & _fm, value_t & _len, int &_sedid, int &_bedge) const
{
	auto& edge = m_ceList[index];
	_vp = edge.vxp;
	_vm = edge.vxm;
	_fp = edge.tpid;
	_fm = edge.tmid;
	_len = edge.length;
	_sedid = edge.sedid;
	_bedge = edge.bedge;
}

void CommonEdge::getCommonEdge(int index, int & _vp, int & _vm, int & _fp, int & _fm, value_t & _len, value_t &_ep,value_t &_em,int &_bedge) const
{
	auto& edge = m_ceList[index];
	_vp = edge.vxp;
	_vm = edge.vxm;
	_fp = edge.tpid;
	_fm = edge.tmid;
	_len = edge.length;
	_ep = edge.ep;
	_em = edge.em;
	_bedge = edge.bedge;
}

void CommonEdge::getCommonEdge(int index, int & _vp, int & _vm, int & _fp, int & _fm, value_t & _len, Qstring &_fla_p,Qstring &_fla_m) const
{
	auto& edge = m_ceList[index];
	_vp = edge.vxp;
	_vm = edge.vxm;
	_fp = edge.tpid;
	_fm = edge.tmid;
	_len = edge.length;
	_fla_p = edge.fla_p;
	_fla_m = edge.fla_m;
}

void CommonEdge::buildCommonEdge(std::shared_ptr<Mesh> mesh_ptr_)
{
    auto triNum = mesh_ptr_->getTriangleNum();
    int  nAID, nBID, nCID;
    for (int i = 0; i < triNum; ++i)
    {
        auto& tri = mesh_ptr_->getTriangleRef(i);
        tri.getVertex(nAID, nBID, nCID);
        updateCommonEdge(nAID, nBID, nCID, i);
    }

    totalEdges = static_cast<int>(m_ceListTemp.size());
    for (int i = 0; i < totalEdges; ++i)
    {
		if (!isCommonEdgeTemp(i))
		{
			m_beList.push_back(std::move(m_ceListTemp[i]));
			continue;
		}
        m_ceList.push_back(std::move(m_ceListTemp[i]));
    }
    m_ceListTemp.clear();
    m_ceListTemp.shrink_to_fit();
    m_ceList.shrink_to_fit();

    VectorR3 vertex1, vertex2;
    for (auto &ce : m_ceList)
    {
        vertex1     = mesh_ptr_->getVertex(ce.v1);
        vertex2     = mesh_ptr_->getVertex(ce.v2);
        ce.length = (vertex1 - vertex2).Norm();
    }
}

void CommonEdge::buildContinuousEdges(std::shared_ptr<Mesh> mesh_ptr_)
{
	buildCommonEdge(mesh_ptr_);
	boundaryEdges = m_beList.size();
	VectorR3 minBox, maxBox;
	mesh_ptr_->getBoundary(minBox, maxBox);
	std::vector<SCE> maxEdges_x, maxEdges_y;
	std::vector<SCE> minEdges_x, minEdges_y;
	for (int i = 0; i < boundaryEdges; ++i)
	{
		VectorR3 v1 = mesh_ptr_->getVertex(m_beList[i].v1);
		VectorR3 v2 = mesh_ptr_->getVertex(m_beList[i].v2);
		if (abs(v1.x - v2.x) < 1e-4)
		{
			if (abs(v1.x - minBox.x) < 1e-4)
				minEdges_x.push_back(std::move(m_beList[i]));
			else if(abs(v1.x - maxBox.x) < 1e-4)
				maxEdges_x.push_back(std::move(m_beList[i]));
		}
		else if (abs(v1.y - v2.y) < 1e-4)
		{
			if (abs(v1.y - minBox.y) < 1e-4)
				minEdges_y.push_back(std::move(m_beList[i]));
			else if (abs(v1.y - maxBox.y) < 1e-4)
				maxEdges_y.push_back(std::move(m_beList[i]));
		}
	}
	if ((minEdges_x.size() == maxEdges_x.size())&&(minEdges_y.size() == maxEdges_y.size()))
		Qcout << "The boundary number is match..." << std::endl;
	else
		Qcout << "The boundary number isn't match...." << std::endl;

	int bnum1 = maxEdges_x.size();
	int bnum2 = maxEdges_y.size();
	auto triNum = mesh_ptr_->getTriangleNum();
	for (int i = 0; i < bnum1; ++i)
	{
		VectorR3 vmax1 = mesh_ptr_->getVertex(maxEdges_x[i].v1);
		VectorR3 vmax2 = mesh_ptr_->getVertex(maxEdges_x[i].v2);
		for (int j = 0; j < bnum1; ++j)
		{
			VectorR3 vmin1 = mesh_ptr_->getVertex(minEdges_x[j].v1);
			VectorR3 vmin2 = mesh_ptr_->getVertex(minEdges_x[j].v2);
			if (abs(vmax1.y + vmax2.y - vmin1.y - vmin2.y) < 1e-4)
			{
				maxEdges_x[i].vxm = minEdges_x[j].vxp;
				maxEdges_x[i].tmid = minEdges_x[j].tpid;
				maxEdges_x[i].bedge = 1;
				maxEdges_x[i].length = (vmax1 - vmax2).Norm();
				m_ceList.push_back(std::move(maxEdges_x[i]));
			}
		}
	}

	for (int i = 0; i < bnum2; ++i)
	{
		VectorR3 vmax1 = mesh_ptr_->getVertex(maxEdges_y[i].v1);
		VectorR3 vmax2 = mesh_ptr_->getVertex(maxEdges_y[i].v2);
		for (int j = 0; j < bnum2; ++j)
		{
			VectorR3 vmin1 = mesh_ptr_->getVertex(minEdges_y[j].v1);
			VectorR3 vmin2 = mesh_ptr_->getVertex(minEdges_y[j].v2);
			if (abs(vmax1.x + vmax2.x - vmin1.x - vmin2.x) < 1e-4)
			{
				maxEdges_y[i].vxm = minEdges_y[j].vxp;
				maxEdges_y[i].tmid = minEdges_y[j].tpid;
				maxEdges_y[i].bedge = 2;
				maxEdges_y[i].length = (vmax1 - vmax2).Norm();
				m_ceList.push_back(std::move(maxEdges_y[i]));
			}
		}
	}
	minEdges_x.clear();
	minEdges_y.clear();
	maxEdges_x.clear();
	maxEdges_y.clear();
	m_ceList.shrink_to_fit();

	std::stringstream tmpcelist;
	for (auto &ce : m_ceList)
	{
		tmpcelist << ce.v1 + 1 << '\t' << ce.v2 + 1 << '\t' << ce.vxp + 1 << '\t' << ce.vxm + 1 << '\t' << ce.tpid + 1 << '\t' << ce.tmid + 1 <<'\t'<<ce.bedge <<'\n';
	}
	Qofstream tmpfile("temp_celist.dat", std::ios::out);
	tmpfile << tmpcelist.rdbuf();
	tmpfile.flush();
	Qcout << "It is already output the commonEdge!" << std::endl;
}

void CommonEdge::buildSEDCommonEdge_com(std::shared_ptr<Mesh> mesh_ptr_)
{
	auto triNum = mesh_ptr_->getTriangleNum();
	int  nAID, nBID, nCID;
	for (int i = 0; i < triNum; ++i)
	{
		auto& tri = mesh_ptr_->getTriangleRef(i);
		tri.getVertex(nAID, nBID, nCID);
		updateCommonEdge(nAID, nBID, nCID, i);
	}

	totalEdges = static_cast<int>(m_ceListTemp.size());
	for (int i = 0; i < totalEdges; ++i)
	{
		if (!isCommonEdgeTemp(i))
		{
			m_beList.push_back(std::move(m_ceListTemp[i]));
			continue;
		}
		m_ceList_SED0.push_back(std::move(m_ceListTemp[i]));
		m_ceList_SED1.push_back(std::move(m_ceListTemp[i]));
		m_ceList_SED2.push_back(std::move(m_ceListTemp[i]));
		m_ceList_SED3.push_back(std::move(m_ceListTemp[i]));
		m_ceList_SED4.push_back(std::move(m_ceListTemp[i]));
		m_ceList_SED5.push_back(std::move(m_ceListTemp[i]));
		m_ceList_SED6.push_back(std::move(m_ceListTemp[i]));
		m_ceList_SED7.push_back(std::move(m_ceListTemp[i]));
		m_ceList_SED8.push_back(std::move(m_ceListTemp[i]));
	}
	common_ce = m_ceList_SED0.size();
	m_ceListTemp.clear();
	m_ceListTemp.shrink_to_fit();
	//m_ceList.shrink_to_fit();

	VectorR3 vertex1, vertex2;
	for (auto &ce : m_ceList)
	{
		vertex1 = mesh_ptr_->getVertex(ce.v1);
		vertex2 = mesh_ptr_->getVertex(ce.v2);
		ce.length = (vertex1 - vertex2).Norm();
	}
}

void CommonEdge::buildSEDCommonEdge_boundary(std::shared_ptr<Mesh> mesh_ptr_,value_t _Dx,value_t _Dy)
{
	//buildSEDCommonEdge_com(mesh_ptr_);
	boundaryEdges = m_beList.size();
	VectorR3 minBox, maxBox;
	mesh_ptr_->getBoundary(minBox, maxBox);
	value_t dis_x, dis_y, dis_z;
	mesh_ptr_->getSize(dis_x, dis_y, dis_z);
	//maxBox.x = minBox.x + _Dx;
	//maxBox.y = minBox.y + _Dy;
	std::vector<SCE> maxEdges_x, maxEdges_y;
	std::vector<SCE> minEdges_x, minEdges_y;
	for (int i = 0; i < boundaryEdges; ++i)
	{
		VectorR3 v1 = mesh_ptr_->getVertex(m_beList[i].v1);
		VectorR3 v2 = mesh_ptr_->getVertex(m_beList[i].v2);
		if ((abs(v1.x - v2.x) < 1e-4) && (abs(_Dx - dis_x)<1e-4))
		{
			if ((abs(v1.x - minBox.x) < 1e-4)&&(abs(_Dx-dis_x)<1e-4))
				minEdges_x.push_back(std::move(m_beList[i]));
			else if ((abs(v1.x - maxBox.x) < 1e-4)&& (abs(_Dx - dis_x)<1e-4))
				maxEdges_x.push_back(std::move(m_beList[i]));
		}
		else if ((abs(v1.y - v2.y) < 1e-4) && (abs(_Dy - dis_y)<1e-4))
		{
			if ((abs(v1.y - minBox.y) < 1e-4)&& (abs(_Dy - dis_y)<1e-4))
				minEdges_y.push_back(std::move(m_beList[i]));
			else if ((abs(v1.y - maxBox.y) < 1e-4)&& (abs(_Dy - dis_y)<1e-4))
				maxEdges_y.push_back(std::move(m_beList[i]));
		}
	}
	if ((minEdges_x.size() == maxEdges_x.size()) && (minEdges_y.size() == maxEdges_y.size()))
		Qcout << "The boundary number is match..." << std::endl;
	else
		Qcout << "The boundary number isn't match...." << std::endl;

	int bnum1 = maxEdges_x.size();
	int bnum2 = maxEdges_y.size();
	b_ce1 = bnum1;
	b_ce2 = bnum2;
	auto triNum = mesh_ptr_->getTriangleNum();
	for (int i = 0; i < bnum1; ++i)
	{
		VectorR3 vmax1 = mesh_ptr_->getVertex(maxEdges_x[i].v1);
		VectorR3 vmax2 = mesh_ptr_->getVertex(maxEdges_x[i].v2);
		for (int j = 0; j < bnum1; ++j)
		{
			VectorR3 vmin1 = mesh_ptr_->getVertex(minEdges_x[j].v1);
			VectorR3 vmin2 = mesh_ptr_->getVertex(minEdges_x[j].v2);
			if ((abs(vmax1.y + vmax2.y - vmin1.y - vmin2.y) < 1e-3)
				&&(abs(vmax1.z + vmax2.z - vmin1.z - vmin2.z)<1e-3))
			{
				maxEdges_x[i].vxm = minEdges_x[j].vxp;
				maxEdges_x[i].tmid = minEdges_x[j].tpid;
				maxEdges_x[i].bedge = 1;
				maxEdges_x[i].length = (vmax1 - vmax2).Norm();
				m_ceList_SED0.push_back(std::move(maxEdges_x[i]));
				m_ceList_SED1.push_back(std::move(maxEdges_x[i]));
				m_ceList_SED2.push_back(std::move(maxEdges_x[i]));
				m_ceList_SED3.push_back(std::move(maxEdges_x[i]));
				m_ceList_SED4.push_back(std::move(maxEdges_x[i]));
				m_ceList_SED5.push_back(std::move(maxEdges_x[i]));
			}
		}
	}

	for (int i = 0; i < bnum2; ++i)
	{
		VectorR3 vmax1 = mesh_ptr_->getVertex(maxEdges_y[i].v1);
		VectorR3 vmax2 = mesh_ptr_->getVertex(maxEdges_y[i].v2);
		for (int j = 0; j < bnum2; ++j)
		{
			VectorR3 vmin1 = mesh_ptr_->getVertex(minEdges_y[j].v1);
			VectorR3 vmin2 = mesh_ptr_->getVertex(minEdges_y[j].v2);
			if ((abs(vmax1.x + vmax2.x - vmin1.x - vmin2.x) < 1e-3) 
				&& (abs(vmax1.z + vmax2.z - vmin1.z - vmin2.z)<1e-3))
			{
				maxEdges_y[i].vxm = minEdges_y[j].vxp;
				maxEdges_y[i].tmid = minEdges_y[j].tpid;
				maxEdges_y[i].bedge = 2;
				maxEdges_y[i].length = (vmax1 - vmax2).Norm();
				m_ceList_SED0.push_back(std::move(maxEdges_y[i]));
				m_ceList_SED1.push_back(std::move(maxEdges_y[i]));
				m_ceList_SED3.push_back(std::move(maxEdges_y[i]));
				m_ceList_SED4.push_back(std::move(maxEdges_y[i]));
				m_ceList_SED6.push_back(std::move(maxEdges_y[i]));
				m_ceList_SED7.push_back(std::move(maxEdges_y[i]));
			}
		}
	}
	minEdges_x.clear();
	minEdges_y.clear();
	maxEdges_x.clear();
	maxEdges_y.clear();
	//m_ceList.shrink_to_fit();


	
}

void CommonEdge::buildSEDCommonEdge_boundary(std::shared_ptr<Mesh> mesh_ptr_, value_t _Dx, value_t _Dy,value_t _lamda)
{
	//buildSEDCommonEdge_com(mesh_ptr_);
	boundaryEdges = m_beList.size();
	VectorR3 minBox, maxBox;
	mesh_ptr_->getBoundary(minBox, maxBox);
	value_t dis_x, dis_y, dis_z;
	mesh_ptr_->getSize(dis_x, dis_y, dis_z);
	//maxBox.x = minBox.x + _Dx;
	//maxBox.y = minBox.y + _Dy;
	std::vector<SCE> maxEdges_x, maxEdges_y;
	std::vector<SCE> minEdges_x, minEdges_y;
	for (int i = 0; i < boundaryEdges; ++i)
	{
		VectorR3 v1 = mesh_ptr_->getVertex(m_beList[i].v1);
		VectorR3 v2 = mesh_ptr_->getVertex(m_beList[i].v2);
		if ((abs(v1.x - v2.x) < 1e-4*_lamda) && (abs(_Dx - dis_x)<1e-4*_lamda))
		{
			if ((abs(v1.x - minBox.x) < 1e-4*_lamda) && (abs(_Dx - dis_x)<1e-4*_lamda))
				minEdges_x.push_back(std::move(m_beList[i]));
			else if ((abs(v1.x - maxBox.x) < 1e-4*_lamda) && (abs(_Dx - dis_x)<1e-4*_lamda))
				maxEdges_x.push_back(std::move(m_beList[i]));
		}
		else if ((abs(v1.y - v2.y) < 1e-4*_lamda) && (abs(_Dy - dis_y)<1e-4*_lamda))
		{
			if ((abs(v1.y - minBox.y) < 1e-4*_lamda) && (abs(_Dy - dis_y)<1e-4*_lamda))
				minEdges_y.push_back(std::move(m_beList[i]));
			else if ((abs(v1.y - maxBox.y) < 1e-4*_lamda) && (abs(_Dy - dis_y)<1e-4*_lamda))
				maxEdges_y.push_back(std::move(m_beList[i]));
		}
	}
	if ((minEdges_x.size() == maxEdges_x.size()) && (minEdges_y.size() == maxEdges_y.size()))
		Qcout << "The boundary number is match..." << std::endl;
	else
		Qcout << "The boundary number isn't match...." << std::endl;

	int bnum1 = maxEdges_x.size();
	int bnum2 = maxEdges_y.size();
	b_ce1 = bnum1;
	b_ce2 = bnum2;
	auto triNum = mesh_ptr_->getTriangleNum();
	for (int i = 0; i < bnum1; ++i)
	{
		VectorR3 vmax1 = mesh_ptr_->getVertex(maxEdges_x[i].v1);
		VectorR3 vmax2 = mesh_ptr_->getVertex(maxEdges_x[i].v2);
		for (int j = 0; j < bnum1; ++j)
		{
			VectorR3 vmin1 = mesh_ptr_->getVertex(minEdges_x[j].v1);
			VectorR3 vmin2 = mesh_ptr_->getVertex(minEdges_x[j].v2);
			if ((abs(vmax1.y + vmax2.y - vmin1.y - vmin2.y) < 1e-4*_lamda)
				&& (abs(vmax1.z + vmax2.z - vmin1.z - vmin2.z)<1e-4*_lamda))
			{
				maxEdges_x[i].vxm = minEdges_x[j].vxp;
				maxEdges_x[i].tmid = minEdges_x[j].tpid;
				maxEdges_x[i].bedge = 1;
				maxEdges_x[i].length = (vmax1 - vmax2).Norm();
				m_ceList_SED0.push_back(std::move(maxEdges_x[i]));
				m_ceList_SED1.push_back(std::move(maxEdges_x[i]));
				m_ceList_SED2.push_back(std::move(maxEdges_x[i]));
				m_ceList_SED3.push_back(std::move(maxEdges_x[i]));
				m_ceList_SED4.push_back(std::move(maxEdges_x[i]));
				m_ceList_SED5.push_back(std::move(maxEdges_x[i]));
			}
		}
	}

	for (int i = 0; i < bnum2; ++i)
	{
		VectorR3 vmax1 = mesh_ptr_->getVertex(maxEdges_y[i].v1);
		VectorR3 vmax2 = mesh_ptr_->getVertex(maxEdges_y[i].v2);
		for (int j = 0; j < bnum2; ++j)
		{
			VectorR3 vmin1 = mesh_ptr_->getVertex(minEdges_y[j].v1);
			VectorR3 vmin2 = mesh_ptr_->getVertex(minEdges_y[j].v2);
			if ((abs(vmax1.x + vmax2.x - vmin1.x - vmin2.x) < 1e-4*_lamda)
				&& (abs(vmax1.z + vmax2.z - vmin1.z - vmin2.z)<1e-4*_lamda))
			{
				maxEdges_y[i].vxm = minEdges_y[j].vxp;
				maxEdges_y[i].tmid = minEdges_y[j].tpid;
				maxEdges_y[i].bedge = 2;
				maxEdges_y[i].length = (vmax1 - vmax2).Norm();
				m_ceList_SED0.push_back(std::move(maxEdges_y[i]));
				m_ceList_SED1.push_back(std::move(maxEdges_y[i]));
				m_ceList_SED3.push_back(std::move(maxEdges_y[i]));
				m_ceList_SED4.push_back(std::move(maxEdges_y[i]));
				m_ceList_SED6.push_back(std::move(maxEdges_y[i]));
				m_ceList_SED7.push_back(std::move(maxEdges_y[i]));
			}
		}
	}
	minEdges_x.clear();
	minEdges_y.clear();
	maxEdges_x.clear();
	maxEdges_y.clear();
	//m_ceList.shrink_to_fit();



}

/*void CommonEdge::buildSPINCommonEdge_boundary(std::shared_ptr<Mesh> mesh_ptr_,value_t _angle)
{
	//buildSEDCommonEdge_com(mesh_ptr_);
	boundaryEdges = m_beList.size();
	//VectorR3 minBox, maxBox;
	//mesh_ptr_->getBoundary(minBox, maxBox);
	//value_t dis_x, dis_y, dis_z;
	//mesh_ptr_->getSize(dis_x, dis_y, dis_z);
	//maxBox.x = minBox.x + _Dx;
	//maxBox.y = minBox.y + _Dy;
	std::vector<SCE> maxEdges;
	std::vector<SCE> minEdges;
	value_t theta = _angle * DegreesToRadians;
	for (int i = 0; i < boundaryEdges; ++i)
	{
		VectorR3 v1 = mesh_ptr_->getVertex(m_beList[i].v1);
		VectorR3 v2 = mesh_ptr_->getVertex(m_beList[i].v2);
		if (abs((v1.y*v2.x) - (v2.y*v1.x)) < 1e-4)
		{
			if (abs(v1.y -(v1.x*tan(theta))) < 1e-4)
				maxEdges.push_back(std::move(m_beList[i]));
			else if ((abs(v1.x - maxBox.x) < 1e-4) && (abs(_Dx - dis_x)<1e-4))
				maxEdges_x.push_back(std::move(m_beList[i]));
		}
		else if (abs(v1.y - v2.y) < 1e-4)
		{
			if ((abs(v1.y - minBox.y) < 1e-4) && (abs(_Dy - dis_y)<1e-4))
				minEdges_y.push_back(std::move(m_beList[i]));
			else if ((abs(v1.y - maxBox.y) < 1e-4) && (abs(_Dy - dis_y)<1e-4))
				maxEdges_y.push_back(std::move(m_beList[i]));
		}
	}
	if ((minEdges_x.size() == maxEdges_x.size()) && (minEdges_y.size() == maxEdges_y.size()))
		Qcout << "The boundary number is match..." << std::endl;
	else
		Qcout << "The boundary number isn't match...." << std::endl;

	int bnum1 = maxEdges_x.size();
	int bnum2 = maxEdges_y.size();
	auto triNum = mesh_ptr_->getTriangleNum();
	for (int i = 0; i < bnum1; ++i)
	{
		VectorR3 vmax1 = mesh_ptr_->getVertex(maxEdges_x[i].v1);
		VectorR3 vmax2 = mesh_ptr_->getVertex(maxEdges_x[i].v2);
		for (int j = 0; j < bnum1; ++j)
		{
			VectorR3 vmin1 = mesh_ptr_->getVertex(minEdges_x[j].v1);
			VectorR3 vmin2 = mesh_ptr_->getVertex(minEdges_x[j].v2);
			if ((abs(vmax1.y + vmax2.y - vmin1.y - vmin2.y) < 1e-4)
				&& (abs(vmax1.z + vmax2.z - vmin1.z - vmin2.z)<1e-4))
			{
				maxEdges_x[i].vxm = minEdges_x[j].vxp;
				maxEdges_x[i].tmid = minEdges_x[j].tpid;
				maxEdges_x[i].bedge = 1;
				maxEdges_x[i].length = (vmax1 - vmax2).Norm();
				m_ceList_SED0.push_back(std::move(maxEdges_x[i]));
				m_ceList_SED1.push_back(std::move(maxEdges_x[i]));
				m_ceList_SED2.push_back(std::move(maxEdges_x[i]));
				m_ceList_SED3.push_back(std::move(maxEdges_x[i]));
				m_ceList_SED4.push_back(std::move(maxEdges_x[i]));
				m_ceList_SED5.push_back(std::move(maxEdges_x[i]));
			}
		}
	}

	for (int i = 0; i < bnum2; ++i)
	{
		VectorR3 vmax1 = mesh_ptr_->getVertex(maxEdges_y[i].v1);
		VectorR3 vmax2 = mesh_ptr_->getVertex(maxEdges_y[i].v2);
		for (int j = 0; j < bnum2; ++j)
		{
			VectorR3 vmin1 = mesh_ptr_->getVertex(minEdges_y[j].v1);
			VectorR3 vmin2 = mesh_ptr_->getVertex(minEdges_y[j].v2);
			if ((abs(vmax1.x + vmax2.x - vmin1.x - vmin2.x) < 1e-4)
				&& (abs(vmax1.z + vmax2.z - vmin1.z - vmin2.z)<1e-4))
			{
				maxEdges_y[i].vxm = minEdges_y[j].vxp;
				maxEdges_y[i].tmid = minEdges_y[j].tpid;
				maxEdges_y[i].bedge = 2;
				maxEdges_y[i].length = (vmax1 - vmax2).Norm();
				m_ceList_SED0.push_back(std::move(maxEdges_y[i]));
				m_ceList_SED1.push_back(std::move(maxEdges_y[i]));
				m_ceList_SED3.push_back(std::move(maxEdges_y[i]));
				m_ceList_SED4.push_back(std::move(maxEdges_y[i]));
				m_ceList_SED6.push_back(std::move(maxEdges_y[i]));
				m_ceList_SED7.push_back(std::move(maxEdges_y[i]));
			}
		}
	}
	minEdges_x.clear();
	minEdges_y.clear();
	maxEdges_x.clear();
	maxEdges_y.clear();
	//m_ceList.shrink_to_fit();

}*/

void CommonEdge::buildIBCCommonEdge(std::shared_ptr<Mesh> mesh_ptr_)
{
	auto triNum = mesh_ptr_->getTriangleNum();
	int  nAID, nBID, nCID;
	//value_t eps;
	Qstring fla;
	for (int i = 0; i < triNum; ++i)
	{
		auto& tri = mesh_ptr_->getTriangleRef(i);
		tri.getVertex(nAID, nBID, nCID);
		tri.getFla(fla);
		updateIBCCommonEdge(nAID, nBID, nCID, i, fla);
		//d_ptList.push_back(SCE(nAID, nBID, nCID, i, fla));       //将冲击基函数的三角形保存
	}

	totalEdges = static_cast<int>(m_ceListTemp.size());
	for (int i = 0; i < totalEdges; ++i)
	{
		/*if (!isCommonEdgeTemp(i))
		{
			m_beList.push_back(std::move(m_ceListTemp[i]));
			continue;
		}*/
		m_ceList.push_back(std::move(m_ceListTemp[i]));
	}
	m_ceListTemp.clear();
	m_ceListTemp.shrink_to_fit();
	m_ceList.shrink_to_fit();

	VectorR3 vertex1, vertex2;
	for (auto &ce : m_ceList)
	{
		vertex1 = mesh_ptr_->getVertex(ce.v1);
		vertex2 = mesh_ptr_->getVertex(ce.v2);
		ce.length = (vertex1 - vertex2).Norm();
	}
}

void CommonEdge::combineSEDCommonEdge(std::shared_ptr<Mesh> mesh_ptr_)
{
	int unk_SED;

	m_ce_size[0] = 0;
	unk_SED = m_ceList_SED0.size();
	for (int i = 0; i < unk_SED; i++)
	{
		m_ceList_SED0[i].sedid = 0;
		m_ceList.push_back(std::move(m_ceList_SED0[i]));
	}
	m_ce_size[1] = m_ceList.size();

	unk_SED = m_ceList_SED1.size();
	for (int i = 0; i < unk_SED; i++)
	{
		m_ceList_SED1[i].sedid = 1;
		m_ceList.push_back(std::move(m_ceList_SED1[i]));
	}
	m_ce_size[2] = m_ceList.size();

	unk_SED = m_ceList_SED2.size();
	for (int i = 0; i < unk_SED; i++)
	{
		m_ceList_SED2[i].sedid = 2;
		m_ceList.push_back(std::move(m_ceList_SED2[i]));
	}
	m_ce_size[3] = m_ceList.size();

	unk_SED = m_ceList_SED3.size();
	for (int i = 0; i < unk_SED; i++)
	{
		m_ceList_SED3[i].sedid = 3;
		m_ceList.push_back(std::move(m_ceList_SED3[i]));
	}
	m_ce_size[4] = m_ceList.size();

	unk_SED = m_ceList_SED4.size();
	for (int i = 0; i < unk_SED; i++)
	{
		m_ceList_SED4[i].sedid = 4;
		m_ceList.push_back(std::move(m_ceList_SED4[i]));
	}
	m_ce_size[5] = m_ceList.size();

	unk_SED = m_ceList_SED5.size();
	for (int i = 0; i < unk_SED; i++)
	{
		m_ceList_SED5[i].sedid = 5;
		m_ceList.push_back(std::move(m_ceList_SED5[i]));
	}
	m_ce_size[6] = m_ceList.size();

	unk_SED = m_ceList_SED6.size();
	for (int i = 0; i < unk_SED; i++)
	{
		m_ceList_SED6[i].sedid = 6;
		m_ceList.push_back(std::move(m_ceList_SED6[i]));
	}
	m_ce_size[7] = m_ceList.size();

	unk_SED = m_ceList_SED7.size();
	for (int i = 0; i < unk_SED; i++)
	{
		m_ceList_SED7[i].sedid = 7;
		m_ceList.push_back(std::move(m_ceList_SED7[i]));
	}
	m_ce_size[8] = m_ceList.size();

	unk_SED = m_ceList_SED8.size();
	for (int i = 0; i < unk_SED; i++)
	{
		m_ceList_SED8[i].sedid = 8;
		m_ceList.push_back(std::move(m_ceList_SED8[i]));
	}
	m_ce_size[9] = m_ceList.size();

	VectorR3 vertex1, vertex2;
	for (auto &ce : m_ceList)
	{
		vertex1 = mesh_ptr_->getVertex(ce.v1);
		vertex2 = mesh_ptr_->getVertex(ce.v2);
		ce.length = (vertex1 - vertex2).Norm();
	}

	std::stringstream tmpcelist;
	for (auto &ce : m_ceList)
	{
		tmpcelist << ce.v1 + 1 << '\t' << ce.v2 + 1 << '\t' << ce.vxp + 1 << '\t' << ce.vxm + 1 << '\t' << ce.tpid + 1 << '\t' << ce.tmid + 1 << '\t' << ce.bedge << '\t' << ce.sedid << '\t' << ce.length << '\n';
	}
	Qofstream tmpfile("temp_celist.dat", std::ios::out);
	tmpfile << tmpcelist.rdbuf();
	tmpfile.flush();
	Qcout << "It is already output the commonEdge!" << std::endl;
}

bool CommonEdge::writeLengthData(const std::string & filename, value_t lamda) const
{
    Qofstream output(filename, std::ios::out);
    if (output.fail())
        throw FileError("fail opening file: " + filename);

    auto edgesNum = m_ceList.size();

    output << std::fixed << std::setprecision(4);
    output << "length data: " << output.widen('\n');
    for (size_t i = 0; i < edgesNum; ++i)
    {
        output << std::setw(10) << i + 1 << std::setw(15) << m_ceList[i].length/lamda << output.widen('\n');
    }
    
    return true;
}

void CommonEdge::reportInfo(Qostream & strm) const
{
    size_t mem = m_ceList.size() * m_ceList.size() * 2 * sizeof(value_t) / (1024 * 1024);
    size_t cemem = sizeof(SCE) * m_ceList.size() / (1024 * 1024);
    auto oldState = strm.flags();
    strm << HEADING "CommonEdge Information\n" TRAILING;
    strm << LEVEL1 "Total edges: " << totalEdges << strm.widen('\n')
        << LEVEL1 "Total RWGs: " << m_ceList.size() << strm.widen('\n')
        << LEVEL1 "Memory prediction: " << mem << " MB" << FORMAT_MEMORY(mem) << '\n'
        << LEVEL1 "Memory detail:" << '\n'
        << LEVEL2 "RWG:" << cemem << " MB" << FORMAT_MEMORY(cemem) << '\n';
    strm.flush();
    strm.flags(oldState);
}

void CommonEdge::updateCommonEdge(int nA, int nB, int nC, int faceID)
{
    bool newEdge1 = true, newEdge2 = true, newEdge3 = true;
    int  v1ID, v2ID;

    int edgeNum = static_cast<int>(m_ceListTemp.size());
    for (int i = 0; i < edgeNum; ++i)
    {
        if (isCommonEdgeTemp(i))
            continue;
        getEdgeTempAt(i, v1ID, v2ID);
        if (nA == v2ID&&nB == v1ID)
        {
            newEdge1 = false;
            setEdgeTempAt(i, nC, faceID);
        }
        else if (nB == v2ID&&nC == v1ID)
        {
            newEdge2 = false;
            setEdgeTempAt(i, nA, faceID);
        }
        else if (nC == v2ID&&nA == v1ID)
        {
            newEdge3 = false;
            setEdgeTempAt(i, nB, faceID);
        }
    }
    if (newEdge1)
    {
        m_ceListTemp.push_back(SCE(nA, nB, nC, faceID));
    }
    if (newEdge2)
    {
        m_ceListTemp.push_back(SCE(nB, nC, nA, faceID));
    }
    if (newEdge3)
    {
        m_ceListTemp.push_back(SCE(nC, nA, nB, faceID));
    }
}

void CommonEdge::updateIBCCommonEdge(int nA, int nB, int nC, int faceID, Qstring fla)
{
	bool newEdge1 = true, newEdge2 = true, newEdge3 = true;
	int  v1ID, v2ID;
	//value_t eps_temp;
	Qstring fla_temp;
	int edgeNum = static_cast<int>(m_ceListTemp.size());
	for (int i = 0; i < edgeNum; ++i)
	{
		if (isCommonEdgeTemp(i))
			continue;
		getEdgeTempAt(i, v1ID, v2ID);
		if (nA == v2ID && nB == v1ID)
		{
			newEdge1 = false;
			setEdgeTempAt(i, nC, faceID, fla);
		}
		else if (nB == v2ID && nC == v1ID)
		{
			newEdge2 = false;
			setEdgeTempAt(i, nA, faceID,fla);
		}
		else if (nC == v2ID && nA == v1ID)
		{
			newEdge3 = false;
			setEdgeTempAt(i, nB, faceID, fla);
		}
	}
	if (newEdge1)
	{
		m_ceListTemp.push_back(SCE(nA, nB, nC, faceID, fla_temp));
	}
	if (newEdge2)
	{
		m_ceListTemp.push_back(SCE(nB, nC, nA, faceID, fla_temp));
	}
	if (newEdge3)
	{
		m_ceListTemp.push_back(SCE(nC, nA, nB, faceID, fla_temp));
	}
}

void CommonEdge::setEdgeTempAt(int index, int _vm, int _fp)
{
    m_ceListTemp[index].vxm = _vm;
    m_ceListTemp[index].tmid = _fp;
}

void CommonEdge::setEdgeTempAt(int index, int _vm, int _fp, Qstring _fla)
{
	m_ceListTemp[index].vxm = _vm;
	m_ceListTemp[index].tmid = _fp;
	m_ceListTemp[index].fla_m = _fla;
	m_ceListTemp[index].bedge = 1;
}

void CommonEdge::getEdgeTempAt(int index, int & _v1, int & _v2) const
{
    _v1 = m_ceListTemp[index].v1;
    _v2 = m_ceListTemp[index].v2;
}

void CommonEdge::getEdgeTempAt(int index, int & _v1, int & _v2,Qstring& _fla) const
{
	//_v1 = m_ceListTemp[index].v1;
	//_v2 = m_ceListTemp[index].v2;
	//_fla = m_ceListTemp[index].fla;
}


