#include "DGVSIE.h"
#include "ConfigLoader.h"
#include "Mathdef.h"
#include "Quadrature.h"
#include "VectorC3.h"
#include "tools.h"
#include "CommonEdge.h"
#include "CommonTriangle.h"
#include "Tetrahedron.h"
#include "iml.h"

using namespace component;
using namespace mom;
using namespace math;
using std::setw;


DGVSIE::DGVSIE()
{

}

DGVSIE::~DGVSIE()
{
}

void DGVSIE::init(component::ConfigLoader * ploader)
{
	EM::init(ploader);
	if (rtype_ == policy::ResultType::Rad)
		readExcEdges(dir_ + '/' + folder_name_ + ".rad");

	dir_ = tool::creatFolder(dir_, folder_name_ + "_DGVSIE");

	auto log_name = dir_ + (rtype_ == policy::ResultType::Rad ? "/rad_" : "/sca_") + "runtime.log";
	logger_.open(log_name);

	mesh_ptr_ = std::make_shared<Mesh>();
	ce_ptr_ = std::make_shared<CommonEdge>();
	ct_ptr_ = std::make_shared<CommonTriangle>();

	mesh_ptr_->loadVSIEMesh(ploader->getMeshPath());
	mesh_ptr_->getBoundary(min_boundary, max_boundary);

	//HALF SWG
	ct_ptr_->buildhSWG(mesh_ptr_);
	unknowns_v = ct_ptr_->gethSWGNum();
	//RWG
	ce_ptr_->buildCommonEdge(mesh_ptr_);
	unknowns_s = ce_ptr_->getCommonEdgeNum();

	unknowns_ = unknowns_s + unknowns_v;

	isContinuous = ploader->getisContinuous();
	tet_num = mesh_ptr_->getTetrahedronNum();

	iter_num = ploader->getMaxIterationNum();
	iter_threshold = ploader->getIterationThreshold();

	Isfast = ploader->getIsfast();
	k_ = PI2 * incidence_.freq / cc;
	omiga = PI2*incidence_.freq;
	d_x = (ploader->getd_x())*cc / incidence_.freq;
	d_y = (ploader->getd_y())*cc / incidence_.freq;
	d_z = (ploader->getd_z())*cc / incidence_.freq;
	prepareFSPGF(ploader->getDx(), ploader->getDy(), ploader->getArray_x(), ploader->getArray_y(), ploader->gett_sum());
	//prepareVIE(ploader->getEps1(), ploader->getMu());

	logger_ << ploader->getReport();
	reportInfo(logger_);
	TIME_LOG("init");
}

void DGVSIE::prepareFSPGF(value_t _Dx, value_t _Dy, int _Array_x, int _Array_y, int _t_sum)
{
	auto lamda = cc / incidence_.freq;
	Dx = _Dx*lamda;
	Dy = _Dy*lamda;
	Array_x = _Array_x;
	Array_y = _Array_y;
	t_sum = _t_sum;
}

void DGVSIE::solve()
{
	SEGMENT("Solve");
	Z.zeros(unknowns_, unknowns_);

	Z2K1.zeros(tet_num, tet_num);
	Z2K2.set_size(tet_num, tet_num);
	Z2K3.set_size(tet_num, tet_num);
	Z2K4.zeros(tet_num, tet_num);
	Z5K.zeros(tet_num, tet_num);

	RIVP.zeros(tet_num);
	RRIV.zeros(tet_num);
	IV.zeros(tet_num);

	RIV.resize(tet_num);
	IVP.resize(tet_num);

	LOG(fillZ(), "Filling impedance matrix");

	/*if (Isfast != 1)
	{
	LOG(FillHZ(), "Filling impedance matrix");
	}
	else
	{
	LOG(CalculatePGFgrid(), "CalculatePGFgrid");
	LOG(KernelIntegral(), "Processing the Kernel Integral");
	LOG(FastfillZ(), "Fast filling impedance matrix");
	}*/


	TIME_LOG("fillZ");

	I.zeros(unknowns_);
	V.zeros(unknowns_);

	if (rtype_ == policy::ResultType::Rad)
	{
		LOG(radiateV(), "Filling rad voltage vector");
	}
	else
	{
		LOG(fillV(), "Filling sca voltage vector");
	}
	TIME_LOG("fillV");

	//Qcout << "Iterative solve:" << std::endl;
	//int iterNum = 0;
	//if(!iml::BICGSTAB(Z, V, I, 0.01, 500, iterNum))
	//{
	//  Qcerr << "fatal error: fail solving Z*I = V" << std::endl;
	//  return;
	//}
	//Qcout << "  -> Iterative number: " << iterNum << std::endl;
	int useless = 0;
	Qcout << setw(30) << "Solving matrix equation:";
	if (!arma::solve(I, Z, V, arma::solve_opts::fast))
	//if (!iml::BICGSTAB(Z, V, I, iter_threshold, iter_num, useless))
		throw std::runtime_error("fail solving matrix equation");
	if (writeZIVData())
		Qcout << "The data is already saved!" << std::endl;
	Qcout << "success" << std::endl;
	TIME_LOG("solve");

#ifdef _DEBUG
	//if (!writeZIVData())
	//  Qcerr << "fail writing ZIV data" << std::endl;
#endif
}

void DGVSIE::output()
{
	SEGMENT("Result");
	Result result;
	if (rtype_ == policy::ResultType::Rad)
	{
		LOG(result.getRadiationField(this, &DGVSIE::getEFiled), "Calculating Radiation");
		TIME_LOG("getRadiation");
	}
	else
	{
		LOG(result.getBistaticRCS(this, &DGVSIE::getBiRCS), "Calculating BistaticRCS");
		TIME_LOG("getBistaticRCS");
	}
}

void DGVSIE::clear()
{
	Z.reset();
	I.reset();
	V.reset();
	exc_edge_.clear();

	mesh_ptr_->clear();
	ct_ptr_->clear();
}


void DGVSIE::reportInfo(Qostream& strm) const
{
	mesh_ptr_->reportInfo(strm);
	ct_ptr_->reportInfo(strm);
}

Complex DGVSIE::vZppKernel(VectorR3 * vf4, VectorR3 * vs4, VectorR3 & vf, VectorR3 & vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, value_t &epsp, value_t &epsm)
{
	value_t r, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	//四面体不重合，所以Z1=0

	//z2,z5,z6
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//z2,z5
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg4f[p] - vg4s[q]).Norm();
			G = exp(-J0 * k_ * r) / r;
			Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s  *  G);//a2
			Z5 += wt4[p] * wt4[q] * G;
		}

		//z6
		if (abs(epsp - epsm) > 1e-4)
		{
			for (int q = 0; q < 3; q++)
			{
				r = (vg4f[p] - vg3s[q]).Norm();
				G = exp(-J0*k_*r) / r;
				Z6 += wt4[p] * w3[q] * G;
			}
		}
	}

	//z3,z4
	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 4; q++)
			{
				r = (vg3f[p] - vg4s[q]).Norm();
				G = exp(-J0*k_*r) / r;
				Z3 += w3[p] * wt4[q] * G;
			}

			if (abs((epsm - epsp)) > 1e-4)
			{
				for (int q = 0; q < 3; q++)
				{
					r = (vg3f[p] - vg3s[q]).Norm();
					G = exp(-J0*k_*r) / r;
					Z4 += w3[p] * w3[q] * G;
				}
			}
		}
	}


	Z = (-1.0f*k_* kp / 9.0f) * Z2 + (1.0f / k_) * (-kp*Z3 - (km - kp)*Z4 + kp*Z5 + (km - kp)*Z6);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex DGVSIE::vZpmKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, value_t &epsp, value_t &epsm)
{
	value_t r, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	//四面体不重合，所以Z1=0

	//z2,z5
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg4f[p] - vg4s[q]).Norm();
			G = exp(-J0 * k_ * r) / r;
			Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s  *  G);//a2
			Z5 += wt4[p] * wt4[q] * G;//a3
		}
	}

	//z3
	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 4; q++)
			{
				r = (vg3f[p] - vg4s[q]).Norm();
				G = exp(-J0*k_*r) / r;
				Z3 += w3[p] * wt4[q] * G;
			}
		}
	}


	Z = (k_*km / 9.0f)*Z2 + (1.0f / k_) * (km* Z3 - km*Z5);
	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex DGVSIE::vZmpKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, value_t &epsp, value_t &epsm)
{
	value_t r, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	//四面体不重合，所以Z1=0

	//z2,z5,z6
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg4f[p] - vg4s[q]).Norm();
			G = exp(-J0 * k_ * r) / r;
			Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s  *  G);//a2
			Z5 += wt4[p] * wt4[q] * G;
		}

		if (abs((epsm - epsp)) > 1e-4)
		{
			for (int q = 0; q < 3; q++)
			{
				r = (vg4f[p] - vg3s[q]).Norm();
				G = exp(-J0*k_*r) / r;
				Z6 += wt4[p] * w3[q] * G;
			}
		}
	}

	Z = (k_*kp / 9.0f)*Z2 + (1.0f / k_) *(-kp*Z5 - (km - kp)*Z6);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex DGVSIE::vZmmKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, value_t &epsp, value_t &epsm)
{
	value_t r, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	//四面体不重合，所以Z1=0

	//z2,z5
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg4f[p] - vg4s[q]).Norm();
			G = exp(-J0 * k_ * r) / r;
			Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s  *  G);//a2
			Z5 += wt4[p] * wt4[q] * G;
		}
	}

	Z = (-k_*km / 9.0f)*Z2 + (1.0f / k_) *(km*Z5);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex DGVSIE::eZppSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, value_t &epsp, value_t &epsm)
{
	value_t r, volume, area, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);
	area = Area(vs3[0], vs3[1], vs3[2]);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//四面体不重合，所以Z1=0
	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4s[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}
	//z2,z5,z6
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//z2,z5非奇异部分
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg4f[p] - vg4s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s  *  G);
			Z5 += wt4[p] * wt4[q] * G;
		}
		//z2,z5奇异部分
		VectorR3 Ivp(0.0f, 0.0f, 0.0f);
		value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
			Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
			Iv += 0.5f*IvSingular(v4_3[j], vg4f[p]);
		}
		Z2 += wt4[p] * (rou_f ^ (Ivp + Iv*(vg4f[p] - vs))) / volume;
		Z5 += wt4[p] * (Iv / volume);

		//z6
		if (abs(epsp - epsm) > 1e-4)
		{
			//z6非奇异部分
			for (int q = 0; q < 3; q++)
			{
				r = (vg4f[p] - vg3s[q]).Norm();
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
				Z6 += wt4[p] * w3[q] * G;
			}
			//z6奇异部分
			Z6 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
		}
	}
	//z3,z4非奇异部分
	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 4; q++)
			{
				r = (vg3f[p] - vg4s[q]).Norm();
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
				Z3 += w3[p] * wt4[q] * G;
			}
			//z3奇异部分
			value_t Iv = 0.0f;
			for (int j = 0; j < 4; j++)
			{
				Iv += 0.5f*IvSingular(v4_3[j], vg3f[p]);
			}
			Z3 += w3[p] * Iv / volume;

			//z4非奇异
			if (abs(epsp - epsm) > 1e-4)
			{
				for (int q = 0; q < 3; q++)
				{
					r = (vg3f[p] - vg3s[q]).Norm();
					G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
					Z4 += w3[p] * w3[q] * G;
				}
				//z4奇异
				Z4 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
			}
		}
	}
	Z = (1.0f / k_) *(PI4 / (9.0f*volume*epsp))*Z1 - (k_* kp / 9.0f)*Z2 + (1.0f / k_) *(-kp*Z3 - (km - kp)*Z4 + kp*Z5 + (km - kp)*Z6);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex DGVSIE::eZpmSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, value_t &epsp, value_t &epsm)
{
	value_t r, volume, area, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);
	area = Area(vs3[0], vs3[1], vs3[2]);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//四面体不重合，所以Z1=0
	//z1
	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4s[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}
	//z2,z5
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//z2,z5非奇异部分
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg4f[p] - vg4s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s  *  G);
			Z5 += wt4[p] * wt4[q] * G;
		}
		//z2,z5奇异部分
		VectorR3 Ivp(0, 0, 0); value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
			Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
			Iv += 0.5f*IvSingular(v4_3[j], vg4f[p]);
		}
		Z2 += wt4[p] * rou_f ^ (Ivp + Iv*(vg4f[p] - vs)) / volume;
		Z5 += wt4[p] * Iv / volume;

		//z6
	}
	//z3,z4非奇异部分
	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 4; q++)
			{
				r = (vg3f[p] - vg4s[q]).Norm();
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
				Z3 += w3[p] * wt4[q] * G;
			}
			//z3奇异部分
			value_t Iv = 0.0f;
			for (int j = 0; j < 4; j++)
			{
				Iv += 0.5f*IvSingular(v4_3[j], vg3f[p]);
			}
			Z3 += w3[p] * Iv / volume;

			//z4
		}
	}

	Z = (-1.0f / k_) *(PI4 / (9.0f * volume * epsm))*Z1 + (k_*km / 9.0f)*Z2 + (1.0f / k_) *(km*Z3 - km*Z5);
	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex DGVSIE::eZmpSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, value_t &epsp, value_t &epsm)
{
	value_t r, volume, area, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);
	area = Area(vs3[0], vs3[1], vs3[2]);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//四面体不重合，所以Z1=0
	//z1
	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4s[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}
	//z2,z5,z6
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//z2,z5非奇异部分
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg4f[p] - vg4s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s  *  G);
			Z5 += wt4[p] * wt4[q] * G;
		}
		//z2,z5奇异部分
		VectorR3 Ivp(0, 0, 0); value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
			Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
			Iv += 0.5f*IvSingular(v4_3[j], vg4f[p]);
		}
		Z2 += wt4[p] * rou_f ^ (Ivp + Iv * (vg4f[p] - vs)) / volume;
		Z5 += wt4[p] * Iv / volume;


		if (abs(epsp - epsm) > 1e-4)
		{
			//z6非奇异部分
			for (int q = 0; q < 3; q++)
			{
				r = (vg4f[p] - vg3s[q]).Norm();
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
				Z6 += wt4[p] * w3[q] * G;
			}
			//z6奇异部分
			value_t Is = 0.0f;
			Z6 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
		}
	}
	//z3,z4非奇异部分

	Z = (-1.0f / k_) *(PI4 / (9.0f*volume*epsp))*Z1 + (k_*kp / 9.0f)*Z2 + (1.0f / k_) *(-kp*Z5 - (km - kp)*Z6);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex DGVSIE::eZmmSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, value_t &epsp, value_t &epsm)
{
	value_t r, volume, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//四面体不重合，所以Z1=0
	//z1
	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4s[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}
	//z2,z5,z6
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//z2,z5非奇异部分
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg4f[p] - vg4s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s  *  G);
			Z5 += wt4[p] * wt4[q] * G;
		}
		//z2,z5奇异部分
		VectorR3 Ivp(0, 0, 0); value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
			Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
			Iv += 0.5f*IvSingular(v4_3[j], vg4f[p]);
		}
		Z2 += wt4[p] * rou_f ^ (Ivp + Iv*(vg4f[p] - vs)) / volume;
		Z5 += wt4[p] * Iv / volume;
	}
	//z3,z4非奇异部分

	Z = (1.0f / k_) *(PI4 / (9.0f*volume*epsm))*Z1 - (k_* km / 9.0f)*Z2 + (1.0f / k_) *(km*Z5);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

value_t DGVSIE::IsSingular(VectorR3 * vs3, VectorR3 & vf)
{
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[1]);
	VectorR3 ni = n.Normalize();
	VectorR3 rou = vf - ni*(ni^vf);

	value_t Is = 0.0f;
	value_t d = (vf - vs3[0]) ^ ni;
	for (int i = 0; i < 3; i++)
	{
		VectorR3 li = (vs3[(i + 1) % 3] - vs3[i]).Normalize();//not normalize
		VectorR3 r_f = vf - vs3[i];
		VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
		VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
		VectorR3 ui = li*ni;
		value_t p0 = abs((vs3[i] - vf) ^ ui);
		value_t lp = (vs3[(i + 1) % 3] - vf) ^ li;
		value_t lm = (vs3[i] - vf) ^ li;
		value_t lp1 = (vf - vs3[(i + 1) % 3]) ^ li;
		value_t lm1 = (vf - vs3[i]) ^ li;
		value_t rp = sqrt(p0*p0 + d*d + lp*lp);
		value_t rm = sqrt(p0*p0 + d*d + lm*lm);

		if (p0 < 1e-8)
			continue;

		VectorR3 pi = ((rou_p - rou) - lp * li).Normalize();
		value_t r02 = p0*p0 + d*d;
		Is += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
	}
	return Is;
}

value_t DGVSIE::IvSingular(VectorR3 * vs3, VectorR3 & vf)
{
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[1]);
	VectorR3 ni = n.Normalize();
	VectorR3 rou = vf - ni*(ni^vf);

	value_t Iv = 0.0f;
	value_t d = (vf - vs3[0]) ^ ni;
	for (int i = 0; i < 3; i++)
	{
		VectorR3 li = (vs3[(i + 1) % 3] - vs3[i]).Normalize();//not normalize
		VectorR3 r_f = vf - vs3[i];
		VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
		VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
		VectorR3 ui = li*ni;
		value_t p0 = abs((vs3[i] - vf) ^ ui);
		value_t lp = (vs3[(i + 1) % 3] - vf) ^ li;
		value_t lm = (vs3[i] - vf) ^ li;
		value_t rp = sqrt(p0*p0 + d*d + lp*lp);
		value_t rm = sqrt(p0*p0 + d*d + lm*lm);

		if (p0 < 1e-8)
			continue;

		VectorR3 pi = ((rou_p - rou) - lp * li).Normalize();
		value_t r02 = p0*p0 + d*d;
		Iv += (pi^ui)*(abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))) - p0*log((rp + lp) / (rm + lm)));
	}
	return d*Iv;
}

VectorR3 DGVSIE::IvpSingular(VectorR3 * vs3, VectorR3 & vf)
{
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[1]);
	VectorR3 ni = n.Normalize();
	VectorR3 rou = vf - ni*(ni^vf);
	value_t Ivp(0.0f);
	value_t d = (vf - vs3[0]) ^ ni;
	for (int i = 0; i < 3; i++)
	{
		VectorR3 li = (vs3[(i + 1) % 3] - vs3[i]).Normalize();//not normalize

		VectorR3 r_f = vf - vs3[i];
		VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
		VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
		VectorR3 ui = li*ni;
		value_t p0 = abs((vs3[i] - vf) ^ ui);
		value_t lp = (vs3[(i + 1) % 3] - vf) ^ li;
		value_t lm = (vs3[i] - vf) ^ li;
		value_t rp = sqrt(p0*p0 + d*d + lp*lp);
		value_t rm = sqrt(p0*p0 + d*d + lm*lm);

		if (p0 < 1e-8)
			continue;

		VectorR3 pi = ((rou_p - rou) - lp * li).Normalize();
		value_t r02 = p0*p0 + d*d;
		value_t d3 = abs(d*d*d);
		value_t d2 = d*d;
		Ivp += (pi^ui)*(((p0*(r02 + 2.0f*d2)) / 2.0f)*log((rp + lp) / (rm + lm)) + (p0 / 2.0f)*(lp*rp - lm*rm) - d3*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
	}
	return Ivp*ni;
}

void DGVSIE::TetToTri(VectorR3 * v4, VectorR3 ** v4_3)
{
	v4_3[0][0] = v4[0];
	v4_3[0][1] = v4[2];
	v4_3[0][2] = v4[1];

	v4_3[1][0] = v4[0];
	v4_3[1][1] = v4[3];
	v4_3[1][2] = v4[2];

	v4_3[2][0] = v4[0];
	v4_3[2][1] = v4[1];
	v4_3[2][2] = v4[3];

	v4_3[3][0] = v4[1];
	v4_3[3][1] = v4[2];
	v4_3[3][2] = v4[3];
}

Complex DGVSIE::z1(VectorR3 * vf4, VectorR3 * vs4, VectorR3 & vf, VectorR3 & vs, value_t & epsp, value_t & epsm)
{

	return Complex();
}

Complex DGVSIE::eVKernel(VectorR3 *vp4, VectorR3 *vm4, VectorR3 &vp, VectorR3 &vm, int &_nvm)
{
	VectorR3 vgp[5], vgm[5];
	VectorR3 rou_p, rou_m; //ρ+-
	Complex Ve;
	value_t kp = 0.0f, km = 0.0f;
	//kp = (epsp - 1.0f) / epsp;
	//km = (epsm - 1.0f) / epsm;

	Gauss5Point(vp4[0], vp4[1], vp4[2], vp4[3], vgp);
	Gauss5Point(vm4[0], vm4[1], vm4[2], vm4[3], vgm);

	VectorC3 Ei(0, 0, 0);
	Complex Vep(0, 0), Vem(0, 0), G(0, 0);

	for (int a = 0; a < 5; a++)
	{
		//////////////////////////////////////////////////////////////////////////
		rou_p = vgp[a] - vp;
		G = exp(-J0 * k_ * (vgp[a] ^ inc_k_));

		Ei = G * inc_e_;
		Vep += wt5[a] * (rou_p ^ Ei);

		//////////////////////////////////////////////////////////////////////////
		if (_nvm != -1)
		{
			rou_m = vm - vgm[a];
			G = exp(-J0 * k_ * (vgm[a] ^ inc_k_));

			Ei = G * inc_e_;
			Vem += wt5[a] * (rou_m ^ Ei);
		}

	}

	Ve = Vep + Vem;

	return Ve;

}

Complex DGVSIE::eHVKernel(VectorR3 *v4, VectorR3 &vx)
{
	VectorR3 vg4[4];
	VectorR3 rou;
	VectorC3 Ei(0, 0, 0);
	Complex Ve(0, 0), G(0, 0);
	Gauss4Point(v4[0], v4[1], v4[2], v4[3], vg4);

	for (int p = 0; p < 4; p++)
	{
		rou = vg4[p] - vx;
		G = exp(-J0*k_*(vg4[p] ^ inc_k_));
		Ei = G*inc_e_;

		Ve += wt4[p] * (Ei^rou);
	}

	return Ve;
}

Complex DGVSIE::eMVKernel(VectorR3 *vp3, VectorR3 *vm3, VectorR3 &vp, VectorR3 &vm)
{
	VectorR3 vgp[7], vgm[7];
	VectorR3 rou_p, rou_m; //ρ+-
	Complex Ve;

	Gauss7Point(vp3[0], vp3[1], vp3[2], vgp);
	Gauss7Point(vm3[0], vm3[1], vm3[2], vgm);

	VectorC3 Ei(0, 0, 0)/*, Hi(0, 0, 0)*/;
	Complex Vep(0, 0), Vem(0, 0)/*, Vmp(0, 0), Vmm(0, 0)*/, G(0, 0);

	for (int a = 0; a < 7; a++)
	{
		//////////////////////////////////////////////////////////////////////////
		rou_p = vgp[a] - vp;
		G = exp(-J0 * k_ * (vgp[a] ^ inc_k_));

		Ei = G * inc_e_;
		Vep += w7[a] * (rou_p ^ Ei);

		//////////////////////////////////////////////////////////////////////////
		rou_m = vm - vgm[a];
		G = exp(-J0 * k_ * (vgm[a] ^ inc_k_));

		Ei = G * inc_e_;
		Vem += w7[a] * (rou_m ^ Ei);
	}

	Ve = (Vep + Vem);

	return Ve;
	//后面要乘以 （0.5*L）
}

void DGVSIE::fillZ()
{

	FillZVV();
	FillZVS();
	FillZSV();
	FillZSS();

	/*int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3];
	Tetrahedron tet_plu, tet_min;
	value_t a_fld, a_src;//场与源的公共面面积
	value_t epsp, epsm;//源处介电常数

					   //Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm_f, vm_s, v_ct1, v_ct2, v_ct3;
	Complex coef = omiga*Z0 / PI4;

	tool::BarAndPercent bar_perc;
	for (int f = 0; f < unknowns_; f++)
	{
		bar_perc(f + 1, unknowns_);

		ct_ptr_->getCommonTriangle(f, vp, vm_f, t_fld_plu, t_fld_min, a_fld);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);
		if (vm_f != -1)
		{
			v_fld_plu = mesh_ptr_->getVertex(vp);//场顶点+
			v_fld_min = mesh_ptr_->getVertex(vm_f);//场顶点-

			tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//场四面体-

			tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);
		}
		else
		{
			v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
			tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
			tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);
		}

		for (int s = 0; s < unknowns_; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ct_ptr_->getCommonTriangle(s, vp, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm);
			ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
			v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);
			if (vm_s != -1)
			{
				v_src_plu = mesh_ptr_->getVertex(vp);
				v_src_min = mesh_ptr_->getVertex(vm_s);

				tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
				tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);

				tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);
				tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);
			}
			else
			{
				v_src_plu = mesh_ptr_->getVertex(vp);
				tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
				tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);
			}
			//field tetrahedron+ <-->source tetrahedron+ 
			if (t_fld_plu == t_src_plu)
			{
				Zpp = eZppSingular(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, epsp, epsm);
			}
			else
			{
				Zpp = vZppKernel(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, epsp, epsm);
			}
			//field tetrahedron+ <-->source tetrahedron-
			if (t_fld_plu == t_src_min)
			{
				Zpm = eZpmSingular(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, epsp, epsm);
			}
			else if (t_src_min != -1)
			{
				Zpm = vZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, epsp, epsm);
			}
			//field tetrahedron- <-->source tetrahedron+
			if (t_fld_min == t_src_plu)
			{
				Zmp = eZmpSingular(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, epsp, epsm);
			}
			else if (t_fld_min != -1)
			{
				Zmp = vZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, epsp, epsm);
			}
			//field tetrahedron- <-->source tetrahedron-
			if ((t_fld_min != -1) && (t_src_min != -1))
			{
				if (t_fld_min == t_src_min)
				{
					Zmm = eZmmSingular(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, epsp, epsm);
				}
				else
					Zmm = vZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, epsp, epsm);
			}

			Z(f, s) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}*/
}

Complex DGVSIE::FZppKernel(VectorR3 * vf4, VectorR3 * vs4, VectorR3 & vf, VectorR3 & vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_p, value_t &epsp, value_t &epsm)
{
	value_t kp, km;
	VectorR3 r;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	//四面体不重合，所以Z1=0

	//z2,z5
	Z2 = Z2K1(t_f_p, t_s_p) - (vf^Z2K2(t_f_p, t_s_p)) - (vs^Z2K3(t_f_p, t_s_p)) + ((vs^vf)*Z2K4(t_f_p, t_s_p));
	//z5
	Z5 = Z5K(t_f_p, t_s_p);
	//z6
	for (int p = 0; p < 4; p++)
	{
		if (abs(epsp - epsm) > 1e-4)
		{
			for (int q = 0; q < 3; q++)
			{
				//r = (vg4f[p] - vg3s[q]).Norm();
				//G = exp(-J0*k_*r) / r;
				r = vg4f[p] - vg3s[q];

				if ((abs(r.x) < d_x) && (abs(r.y) < d_y) && (abs(r.z) < d_z))
				{
					G = FSPGF(r, t_sum, Dx, Dy);
				}
				else
					G = InterpolarPGF(r);

				Z6 += wt4[p] * w3[q] * G;
			}
		}
	}

	//z3,z4
	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 4; q++)
			{
				//r = (vg3f[p] - vg4s[q]).Norm();
				//G = exp(-J0*k_*r) / r;
				r = vg3f[p] - vg4s[q];

				if ((abs(r.x) < d_x) && (abs(r.y) < d_y) && (abs(r.z) < d_z))
				{
					G = FSPGF(r, t_sum, Dx, Dy);
				}
				else
					G = InterpolarPGF(r);

				Z3 += w3[p] * wt4[q] * G;
			}

			if (abs((epsm - epsp)) > 1e-4)
			{
				for (int q = 0; q < 3; q++)
				{
					//r = (vg3f[p] - vg3s[q]).Norm();
					//G = exp(-J0*k_*r) / r;
					r = vg3f[p] - vg3s[q];

					if ((abs(r.x) < d_x) && (abs(r.y) < d_y) && (abs(r.z) < d_z))
					{
						G = FSPGF(r, t_sum, Dx, Dy);
					}
					else
						G = InterpolarPGF(r);

					Z4 += w3[p] * w3[q] * G;
				}
			}
		}
	}


	Z = (-1.0f*k_* kp / 9.0f) * Z2 + (1.0f / k_) * (-kp*Z3 - (km - kp)*Z4 + kp*Z5 + (km - kp)*Z6);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex DGVSIE::FZpmKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_m, value_t &epsp, value_t &epsm)
{
	value_t kp, km;
	VectorR3 r;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	//四面体不重合，所以Z1=0

	//z2,z5
	Z2 = Z2K1(t_f_p, t_s_m) - (vf^Z2K2(t_f_p, t_s_m)) - (vs^Z2K3(t_f_p, t_s_m)) + ((vs^vf)*Z2K4(t_f_p, t_s_m));
	Z5 = Z5K(t_f_p, t_s_m);
	//z3
	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 4; q++)
			{
				//r = (vg3f[p] - vg4s[q]).Norm();
				//G = exp(-J0*k_*r) / r;
				r = vg3f[p] - vg4s[q];

				if ((abs(r.x) < d_x) && (abs(r.y) < d_y) && (abs(r.z) < d_z))
				{
					G = FSPGF(r, t_sum, Dx, Dy);
				}
				else
					G = InterpolarPGF(r);

				Z3 += w3[p] * wt4[q] * G;
			}
		}
	}


	Z = (k_*km / 9.0f)*Z2 + (1.0f / k_) * (km* Z3 - km*Z5);
	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex DGVSIE::FZmpKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_p, value_t &epsp, value_t &epsm)
{
	value_t kp, km;
	VectorR3 r;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	//四面体不重合，所以Z1=0

	//z2,z5
	Z2 = Z2K1(t_f_m, t_s_p) - (vf^Z2K2(t_f_m, t_s_p)) - (vs^Z2K3(t_f_m, t_s_p)) + ((vs^vf)*Z2K4(t_f_m, t_s_p));
	Z5 = Z5K(t_f_m, t_s_p);

	//z6
	for (int p = 0; p < 4; p++)
	{
		if (abs((epsm - epsp)) > 1e-4)
		{
			for (int q = 0; q < 3; q++)
			{
				//r = (vg4f[p] - vg3s[q]).Norm();
				//G = exp(-J0*k_*r) / r;
				r = vg4f[p] - vg3s[q];

				if ((abs(r.x) < d_x) && (abs(r.y) < d_y) && (abs(r.z) < d_z))
				{
					G = FSPGF(r, t_sum, Dx, Dy);
				}
				else
					G = InterpolarPGF(r);

				Z6 += wt4[p] * w3[q] * G;
			}
		}
	}

	Z = (k_*kp / 9.0f)*Z2 + (1.0f / k_) *(-kp*Z5 - (km - kp)*Z6);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex DGVSIE::FZmmKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_m, value_t &epsp, value_t &epsm)
{
	value_t kp, km;
	VectorR3 r;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	//四面体不重合，所以Z1=0

	//z2,z5
	Z2 = Z2K1(t_f_m, t_s_m) - (vf^Z2K2(t_f_m, t_s_m)) - (vs^Z2K3(t_f_m, t_s_m)) + ((vs^vf)*Z2K4(t_f_m, t_s_m));
	Z5 = Z5K(t_f_m, t_s_m);

	Z = (-k_*km / 9.0f)*Z2 + (1.0f / k_) *(km*Z5);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex DGVSIE::FZppSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_p, value_t &epsp, value_t &epsm)
{
	value_t r, volume, area, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);
	area = Area(vs3[0], vs3[1], vs3[2]);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//四面体不重合，所以Z1=0
	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4s[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}
	//z2,z5,z6
	Z2 = Z2K1(t_f_p, t_s_p) - (vf^Z2K2(t_f_p, t_s_p)) - (vs^Z2K3(t_f_p, t_s_p)) + ((vs^vf)*Z2K4(t_f_p, t_s_p)) +
		((RIVP(t_f_p) + RRIV(t_f_p) - (vs^RIV[t_f_p]) - (vf^IVP[t_f_p]) - (vf^RIV[t_f_p]) + ((vf^vs)*IV(t_f_p))) / volume);
	Z5 = Z5K(t_f_p, t_s_p) + (IV(t_f_p) / volume);

	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//z6
		if (abs(epsp - epsm) > 1e-4)
		{
			//z6非奇异部分
			for (int q = 0; q < 3; q++)
			{
				r = (vg4f[p] - vg3s[q]).Norm();
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
				Z6 += wt4[p] * w3[q] * G;
			}
			//z6奇异部分
			Z6 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
		}
	}

	//z3,z4非奇异部分
	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 4; q++)
			{
				r = (vg3f[p] - vg4s[q]).Norm();
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
				Z3 += w3[p] * wt4[q] * G;
			}
			//z3奇异部分
			value_t Iv = 0.0f;
			for (int j = 0; j < 4; j++)
			{
				Iv += 0.5f*IvSingular(v4_3[j], vg3f[p]);
			}
			Z3 += w3[p] * Iv / volume;

			//z4非奇异
			if (abs(epsp - epsm) > 1e-4)
			{
				for (int q = 0; q < 3; q++)
				{
					r = (vg3f[p] - vg3s[q]).Norm();
					G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
					Z4 += w3[p] * w3[q] * G;
				}
				//z4奇异
				Z4 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
			}
		}
	}


	Z = (1.0f / k_) *(PI4 / (9.0f*volume*epsp))*Z1 - (k_* kp / 9.0f)*Z2 + (1.0f / k_) *(-kp*Z3 - (km - kp)*Z4 + kp*Z5 + (km - kp)*Z6);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex DGVSIE::FZpmSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_m, value_t &epsp, value_t &epsm)
{
	value_t r, volume, area, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);
	area = Area(vs3[0], vs3[1], vs3[2]);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4s[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}
	//z2,z5
	Z2 = Z2K1(t_f_p, t_s_m) - (vf^Z2K2(t_f_p, t_s_m)) - (vs^Z2K3(t_f_p, t_s_m)) + ((vs^vf)*Z2K4(t_f_p, t_s_m)) +
		((RIVP(t_f_p) + RRIV(t_f_p) - (vs^RIV[t_f_p]) - (vf^IVP[t_f_p]) - (vf^RIV[t_f_p]) + ((vf^vs)*IV(t_f_p))) / volume);
	Z5 = Z5K(t_f_p, t_s_m) + (IV(t_f_p) / volume);
	//z3,z4非奇异部分
	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 4; q++)
			{
				r = (vg3f[p] - vg4s[q]).Norm();
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
				Z3 += w3[p] * wt4[q] * G;
			}
			//z3奇异部分
			value_t Iv = 0.0f;
			for (int j = 0; j < 4; j++)
			{
				Iv += 0.5f*IvSingular(v4_3[j], vg3f[p]);
			}
			Z3 += w3[p] * Iv / volume;

			//z4
		}
	}

	Z = (-1.0f / k_) *(PI4 / (9.0f * volume * epsm))*Z1 + (k_*km / 9.0f)*Z2 + (1.0f / k_) *(km*Z3 - km*Z5);
	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex DGVSIE::FZmpSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_p, value_t &epsp, value_t &epsm)
{
	value_t r, volume, area, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);
	area = Area(vs3[0], vs3[1], vs3[2]);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//四面体不重合，所以Z1=0
	//z1
	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4s[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}
	//z2,z5
	Z2 = Z2K1(t_f_m, t_s_p) - (vf^Z2K2(t_f_m, t_s_p)) - (vs^Z2K3(t_f_m, t_s_p)) + ((vs^vf)*Z2K4(t_f_m, t_s_p)) +
		((RIVP(t_f_m) + RRIV(t_f_m) - (vs^RIV[t_f_m]) - (vf^IVP[t_f_m]) - (vf^RIV[t_f_m]) + ((vf^vs)*IV(t_f_m))) / volume);
	Z5 = Z5K(t_f_m, t_s_p) + (IV(t_f_m) / volume);

	//z6
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		if (abs(epsp - epsm) > 1e-4)
		{
			//z6非奇异部分
			for (int q = 0; q < 3; q++)
			{
				r = (vg4f[p] - vg3s[q]).Norm();
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
				Z6 += wt4[p] * w3[q] * G;
			}
			//z6奇异部分
			value_t Is = 0.0f;
			Z6 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
		}
	}

	Z = (-1.0f / k_) *(PI4 / (9.0f*volume*epsp))*Z1 + (k_*kp / 9.0f)*Z2 + (1.0f / k_) *(-kp*Z5 - (km - kp)*Z6);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex DGVSIE::FZmmSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_m, value_t &epsp, value_t &epsm)
{
	value_t r, volume, kp, km;
	kp = (epsp - 1.0f) / epsp;
	km = (epsm - 1.0f) / epsm;
	Complex Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), Z(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//四面体不重合，所以Z1=0
	//z1
	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4s[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}
	//z2,z5
	Z2 = Z2K1(t_f_m, t_s_m) - (vf^Z2K2(t_f_m, t_s_m)) - (vs^Z2K3(t_f_m, t_s_m)) + ((vs^vf)*Z2K4(t_f_m, t_s_m)) +
		((RIVP(t_f_m) + RRIV(t_f_m) - (vs^RIV[t_f_m]) - (vf^IVP[t_f_m]) - (vf^RIV[t_f_m]) + ((vf^vs)*IV(t_f_m))) / volume);
	Z5 = Z5K(t_f_m, t_s_m) + (IV(t_f_m) / volume);

	Z = (1.0f / k_) *(PI4 / (9.0f*volume*epsm))*Z1 - (k_* km / 9.0f)*Z2 + (1.0f / k_) *(km*Z5);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex DGVSIE::HZppKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &t_f, int &t_s, value_t eps_f, value_t eps_s)
{
	value_t kn, r;
	//VectorR3 r;
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	kn = 1.0f - (1.0f / eps_s);

	//z2,z3
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;

			r = (vg4f[p] - vg4s[q]).Norm();
			//r = vg4f[p] - vg4s[q];
			G = exp(-J0*k_*r) / r;
			//G = FSPGF(r, t_sum, Dx, Dy);

			Z2 += wt4[p] * wt4[q] * (rou_f^rou_s)*G;
			Z3 += wt4[p] * wt4[q] * G;
		}
		//z4
		for (int q = 0; q < 3; q++)
		{
			r = (vg4f[p] - vg3s[q]).Norm();
			//r = vg4f[p] - vg3s[q];
			G = exp(-J0*k_*r) / r;
			//G = FSPGF(r, t_sum, Dx, Dy);

			Z4 += wt4[p] * w3[q] * G;
		}
	}

	//z5,z6
	for (int p = 0; p < 3; p++)
	{
		for (int q = 0; q < 4; q++)
		{
			r = (vg3f[p] - vg4s[q]).Norm();
			//r = vg3f[p] - vg4s[q];
			G = exp(-J0*k_*r) / r;
			//G = FSPGF(r, t_sum, Dx, Dy); 

			Z5 += w3[p] * wt4[q] * G;
		}

		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			//r = vg3f[p] - vg3s[q];
			G = exp(-J0*k_*r) / r;
			//G = FSPGF(r, t_sum, Dx, Dy);

			Z6 += w3[p] * w3[q] * G;
		}
	}

	Z = -(kn / 9.0f)*Z2 + (kn / (k_*k_))*(Z3 - Z4 - Z5 + Z6);

	return Z;
}

Complex DGVSIE::HZ1Singular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &t_f, int &t_s, value_t eps_f, value_t eps_s)
{
	value_t r, volume, area, kn;
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);

	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);
	area = Area(vs3[0], vs3[1], vs3[2]);
	kn = 1.0f - (1.0f / eps_s);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4f[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}

	//z2,z3
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//z2,z3非奇异部分
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg4f[p] - vg4s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);

			Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s)  *  G;
			Z3 += wt4[p] * wt4[q] * G;
		}
		//z2,z3奇异部分
		VectorR3 Ivp(0.0f, 0.0f, 0.0f);
		value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
			Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
			Iv += 0.5f*IvSingular(v4_3[j], vg4f[p]);
		}
		Z2 += wt4[p] * (rou_f ^ (Ivp + Iv*(vg4f[p] - vs))) / volume;
		Z3 += wt4[p] * (Iv / volume);

		//z4
		//z4非奇异部分
		for (int q = 0; q < 3; q++)
		{
			r = (vg4f[p] - vg3s[q]).Norm();
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G = exp(-J0*k_*r) / r;

			Z4 += wt4[p] * w3[q] * G;
		}
		//z4奇异部分
		//Z4 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
	}

	//z5,z6
	for (int p = 0; p < 3; p++)
	{
		//z5,z6非奇异
		for (int q = 0; q < 4; q++)
		{
			r = (vg3f[p] - vg4s[q]).Norm();
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G = exp(-J0*k_*r) / r;
			Z5 += w3[p] * wt4[q] * G;
		}
		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z6 += w3[p] * w3[q] * G;
		}

		//z6奇异部分
		/*value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
		Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
		Iv += 0.5f*IvSingular(v4_3[j], vg3f[p]);
		}
		Z5 += w3[p] * (Iv / volume);*/
		Z6 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}

	Z = (PI4 / (eps_s*9.0f*volume*k_*k_))*Z1 - (kn / 9.0f)*Z2 + (kn / (k_*k_))*(Z3 - Z4 - Z5 + Z6);

	return Z;
}

Complex DGVSIE::HZ2Singular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &t_f, int &t_s, value_t eps_f, value_t eps_s)
{
	value_t r, volume, area, kn;
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);

	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);
	area = Area(vs3[0], vs3[1], vs3[2]);
	kn = 1.0f - (1.0f / eps_s);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4f[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}

	//z2,z3
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//z2,z3非奇异部分
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg4f[p] - vg4s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s)  *  G;
			Z3 += wt4[p] * wt4[q] * G;
		}
		//z2,z3奇异部分
		VectorR3 Ivp(0.0f, 0.0f, 0.0f);
		value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
			Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
			Iv += 0.5f*IvSingular(v4_3[j], vg4f[p]);
		}
		Z2 += wt4[p] * (rou_f ^ (Ivp + Iv*(vg4f[p] - vs))) / volume;
		Z3 += wt4[p] * (Iv / volume);

		//z4
		//z4非奇异部分
		for (int q = 0; q < 3; q++)
		{
			r = (vg4f[p] - vg3s[q]).Norm();
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G = exp(-J0*k_*r) / r;
			Z4 += wt4[p] * w3[q] * G;
		}
		//z4奇异部分
		//Z4 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
	}

	//z5,z6
	for (int p = 0; p < 3; p++)
	{
		//z5,z6非奇异
		for (int q = 0; q < 4; q++)
		{
			r = (vg3f[p] - vg4s[q]).Norm();
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G = exp(-J0*k_*r) / r;
			Z5 += w3[p] * wt4[q] * G;
		}
		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G = exp(-J0*k_*r) / r;
			Z6 += w3[p] * w3[q] * G;
		}

		//z5,z6奇异部分
		//value_t Iv = 0.0f;
		//for (int j = 0; j < 4; j++)
		//{
		//Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
		//Iv += 0.5f*IvSingular(v4_3[j], vg3f[p]);
		//}
		//Z5 += w3[p] * (Iv / volume);
		//Z6 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}

	Z = (PI4 / (eps_s*9.0f*volume*k_*k_))*Z1 - (kn / 9.0f)*Z2 + (kn / (k_*k_))*(Z3 - Z4 - Z5 + Z6);

	return Z;
}

Complex DGVSIE::HZ3Singular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &t_f, int &t_s, value_t eps_f, value_t eps_s)
{
	value_t  volume, area, kn, r;
	//VectorR3 r;
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], vg1f[1], vg1s[1], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	Gauss1Point(vf3[0], vf3[1], vf3[2], vg1f);
	Gauss1Point(vs3[0], vs3[1], vs3[2], vg1s);

	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);
	area = Area(vs3[0], vs3[1], vs3[2]);
	kn = 1.0f - (1.0f / eps_s);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//z1
	//for (int p = 0; p < 4; p++)
	//{
	//	rou_f = vg4f[p] - vf;
	//	rou_s = vg4s[p] - vs;
	//	Z1 += wt4[p] * (rou_f^rou_s);
	//}

	//z2,z3
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//z2,z3非奇异部分
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg4f[p] - vg4s[q]).Norm();
			//r = vg4f[p] - vg4s[q];
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G = exp(-J0*k_*r) / r;
			//G = FSPGF(r, t_sum, Dx, Dy);
			Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s)  *  G;
			Z3 += wt4[p] * wt4[q] * G;
		}
		//z2,z3奇异部分
		/*VectorR3 Ivp(0.0f, 0.0f, 0.0f);
		value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
		Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
		Iv += 0.5f*IvSingular(v4_3[j], vg4f[p]);
		}
		Z2 += wt4[p] * (rou_f ^ (Ivp + Iv*(vg4f[p] - vs))) / volume;
		Z3 += wt4[p] * (Iv / volume);*/

		//z4
		//z4非奇异部分
		for (int q = 0; q < 3; q++)
		{
			r = (vg4f[p] - vg3s[q]).Norm();
			//r = vg4f[p] - vg3s[q];
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G = exp(-J0*k_*r) / r;
			//G = FSPGF(r, t_sum, Dx, Dy);
			Z4 += wt4[p] * w3[q] * G;
		}
		//z4奇异部分
		//Z4 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
	}

	//z5,z6
	for (int p = 0; p < 3; p++)
	{
		//z5非奇异
		for (int q = 0; q < 4; q++)
		{
			r = (vg3f[p] - vg4s[q]).Norm();
			//r = vg3f[p] - vg4s[q];
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G = exp(-J0*k_*r) / r;
			//G = FSPGF(r, t_sum, Dx, Dy);
			Z5 += w3[p] * wt4[q] * G;
		}
		//for (int q = 0; q < 3; q++)
		//{
		//	r = (vg3f[p] - vg3s[q]).Norm();
		//	G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
		//	//G = exp(-J0*k_*r) / r;
		//	Z6 += w3[p] * w3[q] * G;
		//}

		//z5,z6奇异部分
		//value_t Iv = 0.0f;
		//for (int j = 0; j < 4; j++)
		//{
		//Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
		//Iv += 0.5f*IvSingular(v4_3[j], vg3f[p]);
		//}
		//Z5 += w3[p] * (Iv / volume);
		//Z6 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}
	value_t r1;
	for (int p = 0; p < 3; p++)
	{
		//z6非奇异部分
		for (int q = 0; q < 3; q++)
		{
			r1 = (vg3f[p] - vg3s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r1 + (J0 * k_ / 6.0f) * (k_ * k_ * r1 * r1 - 6);
			//G = exp(-J0*k_*r) / r;
			Z6 += w3[p] * w3[q] * G;
		}
		//z6奇异部分
		Z6 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}

	Z = (PI4 / (eps_s*9.0f*volume*k_*k_))*Z1 - (kn / 9.0f)*Z2 + (kn / (k_*k_))*(Z3 - Z4 - Z5 + Z6);

	return Z;
}


Complex DGVSIE::FHZppKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &t_f, int &t_s, value_t eps_f, value_t eps_s)
{
	value_t kn, r;
	//VectorR3 r;
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	kn = 1.0f - (1.0f / eps_s);

	//z2,z3
	Z2 = Z2K1(t_f, t_s) - (vf^Z2K2(t_f, t_s)) - (vs^Z2K3(t_f, t_s)) + ((vs^vf)*Z2K4(t_f, t_s));
	Z3 = Z5K(t_f, t_s);

	//z4
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		/*for (int q = 0; q < 4; q++)
		{
		rou_s = vg4s[q] - vs;

		r = (vg4f[p] - vg4s[q]).Norm();
		//r = vg4f[p] - vg4s[q];
		G = exp(-J0*k_*r) / r;
		//G = FSPGF(r, t_sum, Dx, Dy);

		Z2 += wt4[p] * wt4[q] * (rou_f^rou_s)*G;
		Z3 += wt4[p] * wt4[q] * G;
		}*/
		//z4
		for (int q = 0; q < 3; q++)
		{
			r = (vg4f[p] - vg3s[q]).Norm();
			//r = vg4f[p] - vg3s[q];
			G = exp(-J0*k_*r) / r;
			//G = FSPGF(r, t_sum, Dx, Dy);

			Z4 += wt4[p] * w3[q] * G;
		}
	}

	//z5,z6
	for (int p = 0; p < 3; p++)
	{
		for (int q = 0; q < 4; q++)
		{
			r = (vg3f[p] - vg4s[q]).Norm();
			//r = vg3f[p] - vg4s[q];
			G = exp(-J0*k_*r) / r;
			//G = FSPGF(r, t_sum, Dx, Dy); 

			Z5 += w3[p] * wt4[q] * G;
		}

		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			//r = vg3f[p] - vg3s[q];
			G = exp(-J0*k_*r) / r;
			//G = FSPGF(r, t_sum, Dx, Dy);

			Z6 += w3[p] * w3[q] * G;
		}
	}

	Z = -(kn / 9.0f)*Z2 + (kn / (k_*k_))*(Z3 - Z4 - Z5 + Z6);

	return Z;
}

Complex DGVSIE::FHZ1Singular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &t_f, int &t_s, value_t eps_f, value_t eps_s)
{
	value_t r, volume, area, kn;
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);

	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);
	area = Area(vs3[0], vs3[1], vs3[2]);
	kn = 1.0f - (1.0f / eps_s);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4f[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}

	//z2,z3
	Z2 = Z2K1(t_f, t_s) - (vf^Z2K2(t_f, t_s)) - (vs^Z2K3(t_f, t_s)) + ((vs^vf)*Z2K4(t_f, t_s)) +
		((RIVP(t_f) + RRIV(t_f) - (vs^RIV[t_f]) - (vf^IVP[t_f]) - (vf^RIV[t_f]) + ((vf^vs)*IV(t_f))) / volume);
	Z3 = Z5K(t_f, t_s) + (IV(t_f) / volume);

	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//z2,z3非奇异部分
		/*for (int q = 0; q < 4; q++)
		{
		rou_s = vg4s[q] - vs;
		r = (vg4f[p] - vg4s[q]).Norm();
		G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);

		Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s)  *  G;
		Z3 += wt4[p] * wt4[q] * G;
		}
		//z2,z3奇异部分
		VectorR3 Ivp(0.0f, 0.0f, 0.0f);
		value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
		Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
		Iv += 0.5f*IvSingular(v4_3[j], vg4f[p]);
		}
		Z2 += wt4[p] * (rou_f ^ (Ivp + Iv*(vg4f[p] - vs))) / volume;
		Z3 += wt4[p] * (Iv / volume);*/

		//z4
		//z4非奇异部分
		for (int q = 0; q < 3; q++)
		{
			r = (vg4f[p] - vg3s[q]).Norm();
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G = exp(-J0*k_*r) / r;

			Z4 += wt4[p] * w3[q] * G;
		}
		//z4奇异部分
		//Z4 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
	}

	//z5,z6
	for (int p = 0; p < 3; p++)
	{
		//z5,z6非奇异
		for (int q = 0; q < 4; q++)
		{
			r = (vg3f[p] - vg4s[q]).Norm();
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G = exp(-J0*k_*r) / r;
			Z5 += w3[p] * wt4[q] * G;
		}
		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z6 += w3[p] * w3[q] * G;
		}

		//z6奇异部分
		/*value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
		Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
		Iv += 0.5f*IvSingular(v4_3[j], vg3f[p]);
		}
		Z5 += w3[p] * (Iv / volume);*/
		Z6 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}

	Z = (PI4 / (eps_s*9.0f*volume*k_*k_))*Z1 - (kn / 9.0f)*Z2 + (kn / (k_*k_))*(Z3 - Z4 - Z5 + Z6);

	return Z;
}

Complex DGVSIE::FHZ2Singular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &t_f, int &t_s, value_t eps_f, value_t eps_s)
{
	value_t r, volume, area, kn;
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);

	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);
	area = Area(vs3[0], vs3[1], vs3[2]);
	kn = 1.0f - (1.0f / eps_s);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//z1
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		rou_s = vg4f[p] - vs;
		Z1 += wt4[p] * (rou_f^rou_s);
	}

	//z2,z3
	Z2 = Z2K1(t_f, t_s) - (vf^Z2K2(t_f, t_s)) - (vs^Z2K3(t_f, t_s)) + ((vs^vf)*Z2K4(t_f, t_s)) +
		((RIVP(t_f) + RRIV(t_f) - (vs^RIV[t_f]) - (vf^IVP[t_f]) - (vf^RIV[t_f]) + ((vf^vs)*IV(t_f))) / volume);
	Z3 = Z5K(t_f, t_s) + (IV(t_f) / volume);

	//z4
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//z2,z3非奇异部分
		/*for (int q = 0; q < 4; q++)
		{
		rou_s = vg4s[q] - vs;
		r = (vg4f[p] - vg4s[q]).Norm();
		G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
		Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s)  *  G;
		Z3 += wt4[p] * wt4[q] * G;
		}
		//z2,z3奇异部分
		VectorR3 Ivp(0.0f, 0.0f, 0.0f);
		value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
		Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
		Iv += 0.5f*IvSingular(v4_3[j], vg4f[p]);
		}
		Z2 += wt4[p] * (rou_f ^ (Ivp + Iv*(vg4f[p] - vs))) / volume;
		Z3 += wt4[p] * (Iv / volume);*/

		//z4
		//z4非奇异部分
		for (int q = 0; q < 3; q++)
		{
			r = (vg4f[p] - vg3s[q]).Norm();
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G = exp(-J0*k_*r) / r;
			Z4 += wt4[p] * w3[q] * G;
		}
		//z4奇异部分
		//Z4 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
	}

	//z5,z6
	for (int p = 0; p < 3; p++)
	{
		//z5,z6非奇异
		for (int q = 0; q < 4; q++)
		{
			r = (vg3f[p] - vg4s[q]).Norm();
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G = exp(-J0*k_*r) / r;
			Z5 += w3[p] * wt4[q] * G;
		}
		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G = exp(-J0*k_*r) / r;
			Z6 += w3[p] * w3[q] * G;
		}

		//z5,z6奇异部分
		//value_t Iv = 0.0f;
		//for (int j = 0; j < 4; j++)
		//{
		//Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
		//Iv += 0.5f*IvSingular(v4_3[j], vg3f[p]);
		//}
		//Z5 += w3[p] * (Iv / volume);
		//Z6 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}

	Z = (PI4 / (eps_s*9.0f*volume*k_*k_))*Z1 - (kn / 9.0f)*Z2 + (kn / (k_*k_))*(Z3 - Z4 - Z5 + Z6);

	return Z;
}

Complex DGVSIE::FHZ3Singular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &t_f, int &t_s, value_t eps_f, value_t eps_s)
{
	value_t  volume, area, kn, r;
	//VectorR3 r;
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), Z4(0.0f, 0.0f), Z5(0.0f, 0.0f), Z6(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg4s[4], vg3f[3], vg3s[3], vg1f[1], vg1s[1], rou_f, rou_s;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	Gauss1Point(vf3[0], vf3[1], vf3[2], vg1f);
	Gauss1Point(vs3[0], vs3[1], vs3[2], vg1s);

	volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);
	area = Area(vs3[0], vs3[1], vs3[2]);
	kn = 1.0f - (1.0f / eps_s);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);

	//z1
	//for (int p = 0; p < 4; p++)
	//{
	//	rou_f = vg4f[p] - vf;
	//	rou_s = vg4s[p] - vs;
	//	Z1 += wt4[p] * (rou_f^rou_s);
	//}

	//z2,z3
	Z2 = Z2K1(t_f, t_s) - (vf^Z2K2(t_f, t_s)) - (vs^Z2K3(t_f, t_s)) + ((vs^vf)*Z2K4(t_f, t_s));
	Z3 = Z5K(t_f, t_s);

	//z4
	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//z2,z3非奇异部分
		/*for (int q = 0; q < 4; q++)
		{
		rou_s = vg4s[q] - vs;
		r = (vg4f[p] - vg4s[q]).Norm();
		//r = vg4f[p] - vg4s[q];
		//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
		G = exp(-J0*k_*r) / r;
		//G = FSPGF(r, t_sum, Dx, Dy);
		Z2 += wt4[p] * wt4[q] * (rou_f ^ rou_s)  *  G;
		Z3 += wt4[p] * wt4[q] * G;
		}
		//z2,z3奇异部分
		VectorR3 Ivp(0.0f, 0.0f, 0.0f);
		value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
		Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
		Iv += 0.5f*IvSingular(v4_3[j], vg4f[p]);
		}
		Z2 += wt4[p] * (rou_f ^ (Ivp + Iv*(vg4f[p] - vs))) / volume;
		Z3 += wt4[p] * (Iv / volume);*/

		//z4
		//z4非奇异部分
		for (int q = 0; q < 3; q++)
		{
			r = (vg4f[p] - vg3s[q]).Norm();
			//r = vg4f[p] - vg3s[q];
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G = exp(-J0*k_*r) / r;
			//G = FSPGF(r, t_sum, Dx, Dy);
			Z4 += wt4[p] * w3[q] * G;
		}
		//z4奇异部分
		//Z4 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
	}

	//z5,z6
	for (int p = 0; p < 3; p++)
	{
		//z5非奇异
		for (int q = 0; q < 4; q++)
		{
			r = (vg3f[p] - vg4s[q]).Norm();
			//r = vg3f[p] - vg4s[q];
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G = exp(-J0*k_*r) / r;
			//G = FSPGF(r, t_sum, Dx, Dy);
			Z5 += w3[p] * wt4[q] * G;
		}
		//for (int q = 0; q < 3; q++)
		//{
		//	r = (vg3f[p] - vg3s[q]).Norm();
		//	G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
		//	//G = exp(-J0*k_*r) / r;
		//	Z6 += w3[p] * w3[q] * G;
		//}

		//z5,z6奇异部分
		//value_t Iv = 0.0f;
		//for (int j = 0; j < 4; j++)
		//{
		//Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg4f[p]);
		//Iv += 0.5f*IvSingular(v4_3[j], vg3f[p]);
		//}
		//Z5 += w3[p] * (Iv / volume);
		//Z6 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}
	value_t r1;
	for (int p = 0; p < 3; p++)
	{
		//z6非奇异部分
		for (int q = 0; q < 3; q++)
		{
			r1 = (vg3f[p] - vg3s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r1 + (J0 * k_ / 6.0f) * (k_ * k_ * r1 * r1 - 6);
			//G = exp(-J0*k_*r) / r;
			Z6 += w3[p] * w3[q] * G;
		}
		//z6奇异部分
		Z6 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}

	Z = (PI4 / (eps_s*9.0f*volume*k_*k_))*Z1 - (kn / 9.0f)*Z2 + (kn / (k_*k_))*(Z3 - Z4 - Z5 + Z6);

	return Z;
}


Complex DGVSIE::ZSVppKernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, value_t eps_s)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg3f[3], vg3s[3], vg4s[4], rou_f, rou_s;
	value_t r, kn;
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	kn = 1.0f - (1.0f / eps_s);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vg3f[p] - vf;
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg3f[p] - vg4s[q]).Norm();
			G = exp(-J0*k_*r) / r;
			Z1 += w3[p] * wt4[q] * (rou_f^rou_s)*G;
			Z2 += w3[p] * wt4[q] * G;
		}

		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			G=exp(-J0*k_*r) / r;
			Z3 += w3[p] * w3[q] * G;
		}
	}

	Z = (-kn / 6.0f)*Z1 + (kn / (k_*k_))*(Z2 - Z3);

	return Z;
}

Complex DGVSIE::ZSVmpKernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, value_t eps_s)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg3f[3], vg3s[3], vg4s[4], rou_f, rou_s;
	value_t r, kn;
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	kn = 1.0f - (1.0f / eps_s);
	for (int p = 0; p < 3; p++)
	{
		rou_f = vf - vg3f[p];
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg3f[p] - vg4s[q]).Norm();
			G = exp(-J0*k_*r) / r;
			Z1 += w3[p] * wt4[q] * (rou_f^rou_s)*G;
			Z2 += w3[p] * wt4[q] * G;
		}

		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			G = exp(-J0*k_*r) / r;
			Z3 += w3[p] * w3[q] * G;
		}
	}

	Z = (-kn / 6.0f)*Z1 - (kn / (k_*k_))*(Z2 - Z3);

	return Z;
}

Complex DGVSIE::ZSVppSingular(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, value_t eps_s)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg3f[3], vg3s[3], vg4s[4], rou_f, rou_s;
	value_t r, kn;
	value_t area = Area(vs3[0], vs3[1], vs3[2]);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	kn = 1.0f - (1.0f / eps_s);
	for (int p = 0; p < 3; p++)
	{
		rou_f = vg3f[p] - vf;
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg3f[p] - vg4s[q]).Norm();
			G = exp(-J0*k_*r) / r;
			Z1 += w3[p] * wt4[q] * (rou_f^rou_s)*G;
			Z2 += w3[p] * wt4[q] * G;
		}

		//z3非奇异部分
		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z3 += w3[p] * w3[q] * G;
		}
		//z3奇异部分
		Z3+= w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}

	Z = (-kn / 6.0f)*Z1 + (kn / (k_*k_))*(Z2 - Z3);

	return Z;
}

Complex DGVSIE::ZSVmpSingular(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, value_t eps_s)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg3f[3], vg3s[3], vg4s[4], rou_f, rou_s;
	value_t r, kn;
	value_t area = Area(vs3[0], vs3[1], vs3[2]);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	kn = 1.0f - (1.0f / eps_s);
	for (int p = 0; p < 3; p++)
	{
		rou_f = vg3f[p] - vf;
		for (int q = 0; q < 4; q++)
		{
			rou_s = vs - vg4s[q];
			r = (vg3f[p] - vg4s[q]).Norm();
			G = exp(-J0*k_*r) / r;
			Z1 += w3[p] * wt4[q] * (rou_f^rou_s)*G;
			Z2 += w3[p] * wt4[q] * G;
		}

		//z3非奇异部分
		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z3 += w3[p] * w3[q] * G;
		}
		//z3奇异部分
		Z3 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}

	Z = (-kn / 6.0f)*Z1 - (kn / (k_*k_))*(Z2 - Z3);

	return Z;
}


Complex DGVSIE::ZVSppKernel(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, value_t eps_s)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg3f[3], vg3s[3], rou_f, rou_s;
	value_t r;
	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);

	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vg3s[q] - vs;
			r = (vg4f[p] - vg3s[q]).Norm();
			G = exp(-J0*k_*r) / r;
			Z1 += wt4[p] * w3[q] * (rou_f^rou_s)*G;
			Z2 += wt4[p] * w3[q] * G;
		}
	}

	for (int p = 0; p < 3; p++)
	{
		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			G = exp(-J0*k_*r) / r;
			Z3 += w3[p] * w3[q] * G;
		}
	}

	Z = (1.0f / 6.0f)*Z1 - (1.0f / (k_*k_))*(Z2 - Z3);//外面乘以J0*k_*Z0/PI4

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}

	return Z;
}

Complex DGVSIE::ZVSpmKernel(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, value_t eps_s)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg3f[3], vg3s[3], rou_f, rou_s;
	value_t r;
	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);

	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vs - vg3s[q];
			r = (vg4f[p] - vg3s[q]).Norm();
			G = exp(-J0*k_*r) / r;
			Z1 += wt4[p] * w3[q] * (rou_f^rou_s)*G;
			Z2 += wt4[p] * w3[q] * G;
		}
	}

	for (int p = 0; p < 3; p++)
	{
		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			G = exp(-J0*k_*r) / r;
			Z3 += w3[p] * w3[q] * G;
		}
	}

	Z = (1.0f / 6.0f)*Z1 + (1.0f / (k_*k_))*(Z2 - Z3);//外面乘以J0*k_*Z0/PI4

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}

	return Z;
}

Complex DGVSIE::ZVSppSingular(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, value_t eps_s)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg3f[3], vg3s[3], rou_f, rou_s;
	value_t r;
	value_t area = Area(vs3[0], vs3[1], vs3[2]);
	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);

	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vg3s[q] - vs;
			r = (vg4f[p] - vg3s[q]).Norm();
			G = exp(-J0*k_*r) / r;
			Z1 += wt4[p] * w3[q] * (rou_f^rou_s)*G;
			Z2 += wt4[p] * w3[q] * G;
		}
	}

	for (int p = 0; p < 3; p++)
	{
		//z3非奇异部分
		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z3 += w3[p] * w3[q] * G;
		}
		//z3奇异部分
		Z3 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}

	Z = (1.0f / 6.0f)*Z1 - (1.0f / (k_*k_))*(Z2 - Z3);//外面乘以J0*k_*Z0/PI4

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}

	return Z;
}

Complex DGVSIE::ZVSpmSingular(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, value_t eps_s)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg3f[3], vg3s[3], rou_f, rou_s;
	value_t r;
	value_t area = Area(vs3[0], vs3[1], vs3[2]);
	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);

	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vs - vg3s[q];
			r = (vg4f[p] - vg3s[q]).Norm();
			G = exp(-J0*k_*r) / r;
			Z1 += wt4[p] * w3[q] * (rou_f^rou_s)*G;
			Z2 += wt4[p] * w3[q] * G;
		}
	}

	for (int p = 0; p < 3; p++)
	{
		//z3非奇异部分
		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z3 += w3[p] * w3[q] * G;
		}
		//z3奇异部分
		Z3 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}

	Z = (1.0f / 6.0f)*Z1 + (1.0f / (k_*k_))*(Z2 - Z3);//外面乘以J0*k_*Z0/PI4

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}

	return Z;
}


Complex DGVSIE::ZSSppKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	value_t r;
	Complex Ze(0, 0), Zm(0, 0), G(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs;
			r = (vgf[p] - vgs[q]).Norm();
			G = exp(-J0 * k_ * r) / r;
			Ze += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * G;
		}
	}
	return Ze;
}

Complex DGVSIE::ZSSpmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	value_t r;
	Complex Ze(0, 0), Zm(0, 0), G(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vs - vgs[q];
			r = (vgf[p] - vgs[q]).Norm();
			G = exp(-J0 * k_ * r) / r;
			Ze += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * G;
		}
	}
	return Ze;
}

Complex DGVSIE::ZSSmpKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	value_t r;
	Complex Ze(0, 0), Zm(0, 0), G(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vf - vgf[p];
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs;
			r = (vgf[p] - vgs[q]).Norm();
			G = exp(-J0 * k_ * r) / r;
			Ze += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * G;
		}
	}
	return Ze;
}

Complex DGVSIE::ZSSmmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	value_t r;
	Complex Ze(0, 0), Zm(0, 0), G(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vf - vgf[p];
		for (int q = 0; q < 3; q++)
		{
			rou_s = vs - vgs[q];
			r = (vgf[p] - vgs[q]).Norm();
			G = exp(-J0 * k_ * r) / r;
			Ze += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * G;
		}
	}
	return Ze;
}

Complex DGVSIE::ZSSppSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	/*Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0);// Ze3(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_s = Area(vs3[0], vs3[1], vs3[2]);

	VectorR3 vsc = (vs3[0] + vs3[1] + vs3[2]) / 3.0;
	VectorR3 rou_sc = vsc - vs; //.1.

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//////////////////////////////////////////////////////////////////////////
	for (int p = 0; p < 3; p++)
	{
	//电场  非奇异部分
	rou_f = vgf[p] - vf;
	for (int q = 0; q < 3; q++)
	{
	rou_s = vgs[q] - vs; //.2.
	value_t r = (vgf[p] - vgs[q]).Norm();
	Complex coff = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
	Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * coff; //.3.
	}
	//电场  奇异部分
	Complex Ztemp(0, 0);
	value_t Is2(0.0f);
	VectorR3 Is1(0.0f, 0.0f, 0.0f);
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[1]);
	VectorR3 ni = n.Normalize();
	VectorR3 rou = vgf[p] - ni*(ni^vgf[p]);
	VectorR3 rou_n = vs - ni*(ni^vs);
	for (int i = 0; i < 3; i++)
	{
	VectorR3 li = vs3[(i + 1) % 3] - vs3[i];//not normalize
	VectorR3 lj = vs3[(i + 2) % 3] - vs3[(i + 1) % 3];//not normalize
	VectorR3 r_f = vgf[p] - vs3[i];
	VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
	VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
	VectorR3 u = li*n;
	VectorR3 ui = u.Normalize();//normalize
	VectorR3 li_1 = li.Normalize();//normalize
	value_t p0 = abs((vs3[i] - vgf[p]) ^ ui);
	value_t d = (vgf[p] - vs3[i]) ^ ni;
	//value_t lp = (rou_p - rou) ^ li_1;
	//value_t lm = (rou_m - rou) ^ li_1;
	value_t lp = (vs3[(i + 1) % 3] - vgf[p]) ^ li_1;
	value_t lm = (vs3[i] - vgf[p]) ^ li_1;
	//value_t rp = sqrt(p0*p0 + d*d + lp*lp);
	//value_t rm = sqrt(p0*p0 + d*d + lm*lm);
	value_t rp = (vgf[p] - vs3[(i + 1) % 3]).Norm();
	value_t rm = (vgf[p] - vs3[i]).Norm();

	if (p0 < 1e-10)
	continue;
	VectorR3 pi = ((rou_p - rou) - lp * li_1).Normalize();
	value_t r02 = p0*p0 + d*d;
	Is1 += (0.5f*ui)*(r02*log((rp + lp) / (rm + lm)) + (lp*rp) - (lm*rm));
	Is2 += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
	}
	Ze2 += w3[p] * (0.25f * (rou_f ^ (Is1 + Is2*(rou - rou_n))) - (Is2 / (k_ * k_)));
	}*/

	/*if (abs(Sigma) > 1e-8)
	{
	for (int i = 0; i < 3; i++)
	{
	rou_f = vgf[i] - vf;
	rou_s = vgf[i] - vs;
	Ze3 += w3[i] * rou_f^rou_s;
	}
	Ze3 = (PI / (S_s*J0*k_*Z0*Sigma))*Ze3;
	}*/
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_s = Area(vs3[0], vs3[1], vs3[2]);

	VectorR3 vsc = (vs3[0] + vs3[1] + vs3[2]) / 3.0;
	VectorR3 rou_sc = vsc - vs; //.1.

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//////////////////////////////////////////////////////////////////////////
	for (int p = 0; p < 3; p++)
	{
		//电场  非奇异部分
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs; //.2.
			value_t r = (vgf[p] - vgs[q]).Norm();
			Complex coff = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * coff; //.3.
		}
		//电场  奇异部分
		Complex Ztemp(0, 0);
		for (int i = 0; i < 3; i++)
		{
			VectorR3 Li = vs3[(i + 1) % 3] - vs3[i];
			Li.Normalize();//边向量

			value_t Lpi = (vs3[(i + 1) % 3] - vgf[p]) ^ Li; //(边正顶点 - 场点)·(边向量)
			value_t Lmi = (vs3[i] - vgf[p]) ^ Li; //(边负顶点 - 场点)·(边向量)
			value_t Rpi = (vs3[(i + 1) % 3] - vgf[p]).Norm(); //边正顶点 到 场点 的距离
			value_t Rmi = (vs3[i % 3] - vgf[p]).Norm(); //边负顶点 到 场点 的距离
			VectorR3 sj0 = vs3[i] + Li * ((vgf[p] - vs3[i]) ^ Li); //垂足在u方向位置矢量
			value_t P0i = (sj0 - vgf[p]).Norm(); //垂足 到 场点 的距离
			if (P0i > 1e-10)
				Ztemp += P0i * log((Rpi + Lpi) / (Rmi + Lmi)); //.4.
		}
		Ze2 += w3[p] * (0.25f * (rou_f ^ rou_sc) - 1.0f / (k_ * k_))*Ztemp;
	}

	Ze = (Ze1 + 1.0f / S_s * Ze2);// +Ze3;

	return Ze;
}

Complex DGVSIE::ZSSpmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	/*Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0);// Ze3(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_s = Area(vs3[0], vs3[1], vs3[2]);

	VectorR3 vsc = (vs3[0] + vs3[1] + vs3[2]) / 3.0;
	VectorR3 rou_sc = vs - vsc; //.1.

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//////////////////////////////////////////////////////////////////////////
	for (int p = 0; p < 3; p++)
	{
	//电场  非奇异部分
	rou_f = vgf[p] - vf;
	for (int q = 0; q < 3; q++)
	{
	rou_s = vs - vgs[q]; //.2.
	value_t r = (vgf[p] - vgs[q]).Norm();
	Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
	Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * coff; //.3.
	}
	//电场  奇异部分

	Complex Ztemp(0.0, 0.0);
	value_t Is2(0.0f);
	VectorR3 Is1(0.0f, 0.0f, 0.0f);
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[1]);
	VectorR3 ni = n.Normalize();
	VectorR3 rou = vgf[p] - ni*(ni^vgf[p]);
	VectorR3 rou_n = vs - ni*(ni^vs);
	for (int i = 0; i < 3; i++)
	{
	VectorR3 li = vs3[(i + 1) % 3] - vs3[i];//not normalize
	VectorR3 lj = vs3[(i + 2) % 3] - vs3[(i + 1) % 3];//not normalize
	VectorR3 r_f = vgf[p] - vs3[i];
	VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
	VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
	VectorR3 u = li*n;
	VectorR3 ui = u / u.Norm();//normalize
	VectorR3 li_1 = li / li.Norm();//normalize
	value_t p0 = abs((vs3[i] - vgf[p]) ^ ui);
	value_t d = (vgf[p] - vs3[i]) ^ ni;
	value_t lp = (rou_p - rou) ^ li_1;
	value_t lm = (rou_m - rou) ^ li_1;
	value_t rp = sqrt(p0*p0 + d*d + lp*lp);
	value_t rm = sqrt(p0*p0 + d*d + lm*lm);

	if (p0 < 1e-10)
	continue;
	VectorR3 pi = ((rou_p - rou) - lp * li_1).Normalize();
	value_t r02 = p0*p0 + d*d;
	Is1 += (0.5f*ui)*(r02*log((rp + lp) / (rm + lm)) + (lp*rp) - (lm*rm));
	Is2 += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
	}
	Ze2 += w3[p] * (0.25f * (-1.0f*rou_f ^ (Is1 + Is2*(rou - rou_n))) + (Is2 / (k_ * k_)));
	}*/

	/*if (abs(Sigma) > 1e-8)
	{
	for (int i = 0; i < 3; i++)
	{
	rou_f = vgf[i] - vf;
	rou_s = vgf[i] - vs;
	Ze3 += w3[i] * rou_f^rou_s;
	}
	Ze3 = (PI / (S_s*J0*k_*Z0*Sigma))*Ze3;
	}*/
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_s = Area(vs3[0], vs3[1], vs3[2]);

	VectorR3 vsc = (vs3[0] + vs3[1] + vs3[2]) / 3.0;
	VectorR3 rou_sc = vs - vsc; //.1.

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//////////////////////////////////////////////////////////////////////////
	for (int p = 0; p < 3; p++)
	{
		//电场  非奇异部分
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vs - vgs[q]; //.2.
			value_t r = (vgf[p] - vgs[q]).Norm();
			Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
			Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * coff; //.3.
		}
		//电场  奇异部分

		Complex Ztemp(0.0, 0.0);
		for (int i = 0; i < 3; i++)
		{
			VectorR3 Li = vs3[(i + 1) % 3] - vs3[i];
			Li.Normalize();//边向量

			value_t Lpi = (vs3[(i + 1) % 3] - vgf[p]) ^ Li; //(边正顶点 - 场点)·(边向量)
			value_t Lmi = (vs3[i] - vgf[p]) ^ Li; //(边负顶点 - 场点)·(边向量)
			value_t Rpi = (vs3[(i + 1) % 3] - vgf[p]).Norm(); //边正顶点 到 场点 的距离
			value_t Rmi = (vs3[i % 3] - vgf[p]).Norm(); //边负顶点 到 场点 的距离
			VectorR3 sj0 = vs3[i] + Li * ((vgf[p] - vs3[i]) ^ Li); //垂足在u方向位置矢量
			value_t P0i = (sj0 - vgf[p]).Norm(); //垂足 到 场点 的距离
			if (P0i > 1e-10)
				Ztemp += P0i * log((Rpi + Lpi) / (Rmi + Lmi)); //.4.
		}
		Ze2 += w3[p] * (0.25f * (rou_f ^ rou_sc) + 1.0f / (k_ * k_))*Ztemp;
	}

	Ze = (Ze1 + 1.0f / S_s * Ze2);// +Ze3;

	return Ze;
}

Complex DGVSIE::ZSSmpSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	/*Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0);// Ze3(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_s = Area(vs3[0], vs3[1], vs3[2]);

	VectorR3 vsc = (vs3[0] + vs3[1] + vs3[2]) / 3.0;
	VectorR3 rou_sc = vsc - vs; //.1.

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//////////////////////////////////////////////////////////////////////////
	for (int p = 0; p < 3; p++)
	{
	//电场  非奇异部分
	rou_f = vf - vgf[p];
	for (int q = 0; q < 3; q++)
	{
	rou_s = vgs[q] - vs; //.2.
	value_t r = (vgf[p] - vgs[q]).Norm();
	Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
	Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * coff; //.3.
	}
	//电场  奇异部分

	Complex Ztemp(0.0, 0.0);
	value_t Is2(0.0f);
	VectorR3 Is1(0.0f, 0.0f, 0.0f);
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[1]);
	VectorR3 ni = n.Normalize();
	VectorR3 rou = vgf[p] - ni*(ni^vgf[p]);
	VectorR3 rou_n = vs - ni*(ni^vs);
	for (int i = 0; i < 3; i++)
	{
	VectorR3 li = vs3[(i + 1) % 3] - vs3[i];//not normalize
	VectorR3 lj = vs3[(i + 2) % 3] - vs3[(i + 1) % 3];//not normalize
	VectorR3 r_f = vgf[p] - vs3[i];
	VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
	VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
	VectorR3 u = li*n;
	VectorR3 ui = u / u.Norm();//normalize
	VectorR3 li_1 = li / li.Norm();//normalize
	value_t p0 = abs((vs3[i] - vgf[p]) ^ ui);
	value_t d = (vgf[p] - vs3[i]) ^ ni;
	value_t lp = (rou_p - rou) ^ li_1;
	value_t lm = (rou_m - rou) ^ li_1;
	value_t rp = sqrt(p0*p0 + d*d + lp*lp);
	value_t rm = sqrt(p0*p0 + d*d + lm*lm);

	if (p0 < 1e-10)
	continue;
	VectorR3 pi = ((rou_p - rou) - lp * li_1).Normalize();
	value_t r02 = p0*p0 + d*d;
	Is1 += (0.5f*ui)*(r02*log((rp + lp) / (rm + lm)) + (lp*rp) - (lm*rm));
	Is2 += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
	}
	Ze2 += w3[p] * (0.25f * (rou_f ^ (Is1 + Is2*(rou - rou_n))) + (Is2 / (k_ * k_)));
	}*/


	/*if (abs(Sigma) > 1e-8)
	{
	for (int i = 0; i < 3; i++)
	{
	rou_f = vgf[i] - vf;
	rou_s = vgf[i] - vs;
	Ze3 += w3[i] * rou_f^rou_s;
	}
	Ze3 = (PI / (S_s*J0*k_*Z0*Sigma))*Ze3;
	}*/
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_s = Area(vs3[0], vs3[1], vs3[2]);

	VectorR3 vsc = (vs3[0] + vs3[1] + vs3[2]) / 3.0;
	VectorR3 rou_sc = vsc - vs; //.1.

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//////////////////////////////////////////////////////////////////////////
	for (int p = 0; p < 3; p++)
	{
		//电场  非奇异部分
		rou_f = vf - vgf[p];
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs; //.2.
			value_t r = (vgf[p] - vgs[q]).Norm();
			Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
			Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * coff; //.3.
		}
		//电场  奇异部分

		Complex Ztemp(0.0, 0.0);
		for (int i = 0; i < 3; i++)
		{
			VectorR3 Li = vs3[(i + 1) % 3] - vs3[i];
			Li.Normalize();//边向量

			value_t Lpi = (vs3[(i + 1) % 3] - vgf[p]) ^ Li; //(边正顶点 - 场点)·(边向量)
			value_t Lmi = (vs3[i] - vgf[p]) ^ Li; //(边负顶点 - 场点)·(边向量)
			value_t Rpi = (vs3[(i + 1) % 3] - vgf[p]).Norm(); //边正顶点 到 场点 的距离
			value_t Rmi = (vs3[i % 3] - vgf[p]).Norm(); //边负顶点 到 场点 的距离
			VectorR3 sj0 = vs3[i] + Li * ((vgf[p] - vs3[i]) ^ Li); //垂足在u方向位置矢量
			value_t P0i = (sj0 - vgf[p]).Norm(); //垂足 到 场点 的距离
			if (P0i > 1e-10)
				Ztemp += P0i * log((Rpi + Lpi) / (Rmi + Lmi)); //.4.
		}
		Ze2 += w3[p] * (0.25f * (rou_f ^ rou_sc) + 1.0f / (k_ * k_))*Ztemp;
	}

	Ze = (Ze1 + 1.0f / S_s * Ze2);// +Ze3;

	return Ze;
}

Complex DGVSIE::ZSSmmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	/*Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0);// Ze3(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_s = Area(vs3[0], vs3[1], vs3[2]);

	VectorR3 vsc = (vs3[0] + vs3[1] + vs3[2]) / 3.0;
	VectorR3 rou_sc = vs - vsc; //.1.

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//////////////////////////////////////////////////////////////////////////
	for (int p = 0; p < 3; p++)
	{
	//电场  非奇异部分
	rou_f = vf - vgf[p];
	for (int q = 0; q < 3; q++)
	{
	rou_s = vs - vgs[q]; //.2.
	value_t r = (vgf[p] - vgs[q]).Norm();
	Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
	Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * coff; //.3.
	}
	//电场  奇异部分

	Complex Ztemp(0.0, 0.0);
	value_t Is2(0.0f);
	VectorR3 Is1(0.0f, 0.0f, 0.0f);
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[1]);
	VectorR3 ni = n.Normalize();
	VectorR3 rou = vgf[p] - (ni^vgf[p])*ni;
	VectorR3 rou_n = vs - (ni^vs)*ni;
	for (int i = 0; i < 3; i++)
	{
	VectorR3 li = vs3[(i + 1) % 3] - vs3[i];//not normalize
	VectorR3 lj = vs3[(i + 2) % 3] - vs3[(i + 1) % 3];//not normalize
	VectorR3 r_f = vgf[p] - vs3[i];
	VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
	VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
	VectorR3 u = li*n;
	VectorR3 ui = u / u.Norm();//normalize
	VectorR3 li_1 = li / li.Norm();//normalize
	value_t p0 = abs((vs3[i] - vgf[p]) ^ ui);
	value_t d = (vgf[p] - vs3[i]) ^ ni;
	value_t lp = (rou_p - rou) ^ li_1;
	value_t lm = (rou_m - rou) ^ li_1;
	value_t rp = sqrt(p0*p0 + d*d + lp*lp);
	value_t rm = sqrt(p0*p0 + d*d + lm*lm);

	if (p0 < 1e-10)
	continue;
	VectorR3 pi = ((rou_p - rou) - lp * li_1).Normalize();
	value_t r02 = p0*p0 + d*d;
	Is1 += (0.5f*ui)*(r02*log((rp + lp) / (rm + lm)) + (lp*rp) - (lm*rm));
	Is2 += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
	}
	Ze2 += w3[p] * (0.25f * (-0.1f*rou_f ^ (Is1 + Is2*(rou - rou_n))) - (Is2 / (k_ * k_)));
	}*/

	/*if (abs(Sigma) > 1e-8)
	{
	for (int i = 0; i < 3; i++)
	{
	rou_f = vgf[i] - vf;
	rou_s = vgf[i] - vs;
	Ze3 += w3[i] * rou_f^rou_s;
	}
	Ze3 = (PI / (S_s*J0*k_*Z0*Sigma))*Ze3;
	}*/
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_s = Area(vs3[0], vs3[1], vs3[2]);

	VectorR3 vsc = (vs3[0] + vs3[1] + vs3[2]) / 3.0;
	VectorR3 rou_sc = vs - vsc; //.1.

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//////////////////////////////////////////////////////////////////////////
	for (int p = 0; p < 3; p++)
	{
		//电场  非奇异部分
		rou_f = vf - vgf[p];
		for (int q = 0; q < 3; q++)
		{
			rou_s = vs - vgs[q]; //.2.
			value_t r = (vgf[p] - vgs[q]).Norm();
			Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
			Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * coff; //.3.
		}
		//电场  奇异部分

		Complex Ztemp(0.0, 0.0);
		for (int i = 0; i < 3; i++)
		{
			VectorR3 Li = vs3[(i + 1) % 3] - vs3[i];
			Li.Normalize();//边向量

			value_t Lpi = (vs3[(i + 1) % 3] - vgf[p]) ^ Li; //(边正顶点 - 场点)·(边向量)
			value_t Lmi = (vs3[i] - vgf[p]) ^ Li; //(边负顶点 - 场点)·(边向量)
			value_t Rpi = (vs3[(i + 1) % 3] - vgf[p]).Norm(); //边正顶点 到 场点 的距离
			value_t Rmi = (vs3[i % 3] - vgf[p]).Norm(); //边负顶点 到 场点 的距离
			VectorR3 sj0 = vs3[i] + Li * ((vgf[p] - vs3[i]) ^ Li); //垂足在u方向位置矢量
			value_t P0i = (sj0 - vgf[p]).Norm(); //垂足 到 场点 的距离
			if (P0i > 1e-10)
				Ztemp += P0i * log((Rpi + Lpi) / (Rmi + Lmi)); //.4.
		}
		Ze2 += w3[p] * (0.25f * (rou_f ^ rou_sc) - 1.0f / (k_ * k_))*Ztemp;
	}

	Ze = (Ze1 + 1.0f / S_s * Ze2);// +Ze3;

	return Ze;
}


Complex DGVSIE::FSPGF(VectorR3 &r, int _t_sum, value_t _Dx, value_t _Dy)
{
	Complex PG(0.0, 0.0), PG1(0.0, 0.0), PG2(0.0, 0.0);
	VectorR3 R = r;
	value_t H;
	double delta = 0.0;
	Complex a(0.0, 0.0);
	H = sqrt(PI / (Dx*Dy));
	for (int i = -_t_sum; i <= _t_sum; i++)
	{
		for (int j = -_t_sum; j <= _t_sum; j++)
		{

			delta = ((float)i*PI / Dx)*((float)i*PI / Dx) + ((float)j*PI / Dy) *((float)j*PI / Dy);
			//R = ((d_x - m*Dx) ^ 2 + (d_y - n*Dy) ^ 2 + (d_z) ^ 2);
			R.x = r.x - i*Dx;
			R.y = r.y - j*Dy;
			value_t R_N = R.Norm();

			if (delta >= ((k_ * k_) / 4.0f))
				a = sqrt(delta - ((k_ * k_) / 4.0f));
			else
				a = J0*(value_t)sqrt(((k_ * k_) / 4.0f - delta));

			PG1 += ((exp(-2.0f*J0 * PI*(i*r.x / Dx + j*r.y / Dy))) / a)*((exp(2.0f * a*r.z))*(erfc(r.z*H + a / H)) + (exp(-2.0f * a*r.z))*(erfc((a / H) - r.z*H)));
			PG2 += (exp(-J0*k_*(inc_k_.x*i*Dx + inc_k_.y*j*Dy)) / R_N)*(erfc(R_N*H + (-J0*k_) / (2.0f * H))*exp(-J0*k_*R_N)).real();
			/*if (isinf(PG1.real()) || isinf(PG1.imag()) || isinf(PG2.imag()) || isinf(PG2.real()))
			{
			Qcout << "This number isn't a number!" << std::endl;
			Qcout << "please cheak this part!" << std::endl;
			}*/
			if (_isnan(PG1.real()) || _isnan(PG2.real()))
			{
				Qcout << "This element is not a number!" << std::endl;
				Qcout << "Please check this part!" << std::endl;
			}
		}
	}

	PG = (PI*PG1) / (2.0f * Dx*Dy) + PG2;
	return PG;
}

Complex DGVSIE::erfz(Complex z)
{
	/*std::complex<double> j0(0.0, 1.0);
	double pi = 3.1415926;
	value_t sqrtpi = sqrt(pi);
	value_t abs_z = abs(z);



	Complex f;

	if (abs_z <= 8)
	{
	int n = 32;
	double x = z.real(), y = z.imag();
	std::complex<double> k1 = 2 / pi * exp(-x*x), k2 = exp(-2.0 * j0  * x*y);
	std::complex<double> s1 = erf(x);
	std::complex<double> s2(0.0, 0.0), s3(0.0, 0.0), s4(0.0, 0.0), s5(0.0, 0.0), s6(0.0, 0.0);
	if (abs(x) > 1e-5)
	{
	s2 = k1 / (4.0 * x)*(1.0 - k2);
	}
	else
	{
	s2 = (j0 / pi) *y;
	}

	f = s1 + s2;

	double xk=0.0, yk=0.0;

	if (abs(y) > 1e-5)
	{
	xk = x;
	yk = y;
	}

	s5 = 0;
	for (int i = 1; i <= n; i++)
	{
	s3 = exp(-(1.0*(i*i)) / 4.0) / (1.0*(i*i + 4 * xk*xk));
	s4 = 2.0 * xk - k2*(2.0 * xk*cosh(i*yk) - 1.0*i*j0*sinh(i*yk));
	s5 = s5 + s3*s4;
	}
	s6 = k1*s5;
	f = f + (arma::cx_float)s6;

	return f;
	}
	else
	{
	if (z.real() < 0)
	z = -z;
	int nmax = 193;
	Complex s = (1.0, 0.0);
	Complex y = 2.0f * z*z;
	for (int i = nmax; i >= 1; i -= 2)
	{
	s = 1.0f - 1.0f*i*(s / y);
	}

	f = 1.0f - s*exp(-z*z) / (sqrtpi*z);

	if (z.real() < 0)
	f = -f;

	if (abs(z.real()) < 1e-8)
	f = f - 1.0f;

	return f;
	}*/

	double x = z.real(), y = z.imag();
	std::complex<double> k1 = 2.0 / PI * exp(-x*x), k2 = exp(-2.0 * (std::complex<double>)J0  * x*y);
	std::complex<double> s1 = erf(x);
	std::complex<double> s2(0.0, 0.0), s3(0.0, 0.0), s4(0.0, 0.0), s5(0.0, 0.0), s6(0.0, 0.0);
	int n = 25;
	Complex f;
	if (abs(x) > 1e-5)
	{
		s2 = k1 / (4.0 * x)*(1.0 - k2);
	}
	else
	{
		s2 = ((std::complex<double>)J0 / (double)PI) *y;
	}


	s5 = 0;
	for (int i = 1; i <= n; i++)
	{
		s3 = exp(-(1.0*(i*i)) / 4.0) / (1.0*(i*i + 4.0 * x*x));
		s4 = 2.0 * x - k2*(2.0 * x*cosh(i*y) - 1.0*i*(std::complex<double>)J0*sinh(i*y));
		s5 = s5 + s3*s4;
	}
	s6 = k1*s5;
	f = (arma::cx_float)(s1 + s2 + s6);

	return f;
}

void DGVSIE::CalculatePGFgrid()
{
	grid_x = std::ceil((max_boundary.x - min_boundary.x) / d_x) + 1;
	grid_y = std::ceil((max_boundary.y - min_boundary.y) / d_y) + 1;
	grid_z = std::ceil((max_boundary.z - min_boundary.z) / d_z) + 1;

	//int grid_t = grid_x*grid_y*grid_z;
	tool::BarAndPercent bar;
	for (int z = 0; z < grid_z; z++)
	{
		bar(z + 1, grid_z);
		for (int y = 0; y < grid_y; y++)
		{
			for (int x = 0; x < grid_x; x++)
			{
				if ((x == 0) && (y == 0) && (z == 0))
				{
					PGFgrid.push_back((0.0f, 0.0f));
					continue;
				}
				VectorR3 r;
				r.x = x*d_x;
				r.y = y*d_y;
				r.z = z*d_z;
				PGFgrid.push_back(FSPGF(r, t_sum, Dx, Dy));
			}
		}
	}
	Qcout << "The PGFgrid's filling is accomplished!" << std::endl;
}

Complex DGVSIE::InterpolarPGF(VectorR3 &r)
{
	int x_seq, y_seq, z_seq;
	x_seq = floor(abs(r.x) / d_x);
	y_seq = floor(abs(r.y) / d_y);
	z_seq = floor(abs(r.z) / d_z);

	int n000 = z_seq*grid_x*grid_y + y_seq*grid_x + x_seq;
	int n100 = n000 + 1;
	int n010 = n000 + grid_x;
	int n110 = n010 + 1;
	int n001 = n000 + grid_x*grid_y;
	int n101 = n001 + 1;
	int n011 = n001 + grid_x;
	int n111 = n011 + 1;

	Complex G000 = PGFgrid[n000];
	Complex G100 = PGFgrid[n100];
	Complex G010 = PGFgrid[n010];
	Complex G110 = PGFgrid[n110];
	Complex G001 = PGFgrid[n001];
	Complex G101 = PGFgrid[n101];
	Complex G011 = PGFgrid[n011];
	Complex G111 = PGFgrid[n111];

	value_t effx = (abs(r.x) - x_seq*d_x) / d_x;
	value_t effy = (abs(r.y) - y_seq*d_y) / d_y;
	value_t effz = (abs(r.z) - z_seq*d_z) / d_z;

	Complex G = (1.0f - effx)*(1.0f - effy)*(1.0f - effz)*G000
		+ (1.0f - effx)*(1.0f - effy)*effz*G001
		+ (1.0f - effx)*effy*(1.0f - effz)*G010
		+ (1.0f - effx)*effy*effz*G011
		+ effx*(1.0f - effy)*(1.0f - effz)*G100
		+ effx*(1.0f - effy)*effz*G101
		+ effx*effy*(1.0f - effz)*G110
		+ effx*effy*effz*G111;

	return G;
}

void DGVSIE::FastfillZ()
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3];
	Tetrahedron tet_plu, tet_min;
	value_t a_fld, a_src;//场与源的公共面面积
	value_t epsp, epsm;//源处介电常数

					   //Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp_f, vp_s, vm_f, vm_s, v_ct1, v_ct2, v_ct3;
	Complex coef = omiga*Z0 / PI4;

	tool::BarAndPercent bar_perc;
	for (int f = 0; f < unknowns_v; f++)
	{
		bar_perc(f + 1, unknowns_v);

		ct_ptr_->getCommonTriangle(f, vp_f, vm_f, t_fld_plu, t_fld_min, a_fld);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);
		//VectorR3 nor_temp = Normal(v_fld_ct3[0], v_fld_ct3[1], v_fld_ct3[2]);

		if (vm_f != -1)
		{
			v_fld_plu = mesh_ptr_->getVertex(vp_f);//场顶点+
			v_fld_min = mesh_ptr_->getVertex(vm_f);//场顶点-

			tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//场四面体-

			tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);
		}
		else
		{
			v_fld_plu = mesh_ptr_->getVertex(vp_f);//场点+
			tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//场四面体+
			tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);

			VectorR3 nor_temp = Normal(v_fld_ct3[0], v_fld_ct3[1], v_fld_ct3[2]);
			if (isContinuous == 1)
			{
				if (abs(abs(nor_temp.z) - 1.0f) > 1e-5)
					vm_f = -2;
			}

		}
		for (int s = 0; s < unknowns_; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ct_ptr_->getCommonTriangle(s, vp_s, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm);
			ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
			v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);
			//VectorR3 nor_temp = Normal(v_fld_ct3[0], v_fld_ct3[1], v_fld_ct3[2]);

			if (vm_s != -1)
			{
				v_src_plu = mesh_ptr_->getVertex(vp_s);
				v_src_min = mesh_ptr_->getVertex(vm_s);

				tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
				tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);

				tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);
				tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);
			}
			else
			{
				v_src_plu = mesh_ptr_->getVertex(vp_s);
				tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
				tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);
				VectorR3 nor_temp = Normal(v_fld_ct3[0], v_fld_ct3[1], v_fld_ct3[2]);
				if (isContinuous == 1)
				{
					if (abs(abs(nor_temp.z) - 1.0f) > 1e-5)
					{
						epsm = epsp;
						vm_s = -2;
					}

				}
			}

			//field tetrahedron+ <-->source tetrahedron+ 
			if (t_fld_plu == t_src_plu)
			{
				Zpp = FZppSingular(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
			}
			else
			{
				Zpp = FZppKernel(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
			}
			//field tetrahedron+ <-->source tetrahedron-
			if (t_fld_plu == t_src_min)
			{
				Zpm = FZpmSingular(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			}
			else if (t_src_min != -1)
			{
				Zpm = FZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			}
			//field tetrahedron- <-->source tetrahedron+
			if (t_fld_min == t_src_plu)
			{
				Zmp = FZmpSingular(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			}
			else if (t_fld_min != -1)
			{
				Zmp = FZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			}
			//field tetrahedron- <-->source tetrahedron-
			if ((t_fld_min != -1) && (t_src_min != -1))
			{
				if (t_fld_min == t_src_min)
				{
					Zmm = FZmmSingular(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
				}
				else
					Zmm = FZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
			}

			Z(f, s) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void DGVSIE::KernelIntegral()
{
	Tetrahedron tet_i, tet_j;
	VectorR3 v_tet_i[4], v_tet_j[4];
	VectorR3 vgi[4], vgj[4];
	tool::BarAndPercent bar;

	for (int i = 0; i < tet_num; i++)
	{
		bar(i + 1, tet_num);
		tet_i = mesh_ptr_->getTetrahedronRef(i);
		tet_i.GetVertices(v_tet_i[0], v_tet_i[1], v_tet_i[2], v_tet_i[3]);
		Gauss4Point(v_tet_i[0], v_tet_i[1], v_tet_i[2], v_tet_i[3], vgi);
		for (int j = i; j < tet_num; j++)
		{
			tet_j = mesh_ptr_->getTetrahedronRef(j);
			tet_j.GetVertices(v_tet_j[0], v_tet_j[1], v_tet_j[2], v_tet_j[3]);
			Gauss4Point(v_tet_j[0], v_tet_j[1], v_tet_j[2], v_tet_j[3], vgj);
			FillZKernel(vgi, vgj, i, j, v_tet_j);
		}
	}
	Qcout << "test" << std::endl;
}

void DGVSIE::FillZKernel(VectorR3 *vm4, VectorR3 *vn4, int &m, int &n, VectorR3 *v4)
{
	Complex G(0, 0);
	value_t r;
	//VectorR3 r;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			r = (vm4[i] - vn4[j]).Norm();
			//r = vm4[i] - vn4[j];

			if (m != n)
			{
				G = exp(-J0*k_*r) / r;
				/*if ((abs(r.x) < d_x) && (abs(r.y) < d_y) && (abs(r.z) < d_z))
				{
				G = FSPGF(r, t_sum, Dx, Dy);
				}
				else
				G = InterpolarPGF(r);*/

				Z2K1(m, n) += wt4[i] * wt4[j] * (vm4[i] ^ vn4[j])*G;
				Z2K1(n, m) = Z2K1(m, n);
				//Qcout << Z2K1(m, n) << std::endl;
				Z2K2(m, n) += wt4[i] * wt4[j] * vn4[j] * G;
				Z2K3(m, n) += wt4[i] * wt4[j] * vm4[i] * G;

				Z2K2(n, m) = Z2K3(m, n);
				Z2K3(n, m) = Z2K2(m, n);

				Z2K4(m, n) += wt4[i] * wt4[j] * G;
				Z2K4(n, m) = Z2K4(m, n);

				Z5K(m, n) = Z2K4(m, n);
				Z5K(n, m) = Z2K4(m, n);
			}
			else
			{
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r* r - 6);

				Z2K1(m, m) += wt4[i] * wt4[j] * (vm4[i] ^ vn4[j])*G;
				Z2K2(m, m) += wt4[i] * wt4[j] * vn4[j] * G;
				Z2K3(m, m) += wt4[i] * wt4[j] * vm4[i] * G;
				Z2K4(m, m) += wt4[i] * wt4[j] * G;
				Z5K(m, m) = Z2K4(m, m);
			}

		}
		if (m == n)
		{
			VectorR3 Ivp(0.0f, 0.0f, 0.0f);
			value_t Iv = 0.0f;

			Pv4_3[0] = v4_3[0];
			Pv4_3[1] = v4_3[1];
			Pv4_3[2] = v4_3[2];
			Pv4_3[3] = v4_3[3];
			TetToTri(v4, Pv4_3);

			for (int k = 0; k < 4; k++)
			{
				Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[k], vm4[i]);
				Iv += 0.5f*IvSingular(v4_3[k], vm4[i]);
			}
			RIVP(m) += wt4[i] * (vm4[i] ^ Ivp);
			RRIV(m) += wt4[i] * (vm4[i] ^ vm4[i])*Iv;
			RIV[m] += wt4[i] * Iv*vm4[i];
			IVP[m] += wt4[i] * Ivp;
			IV(m) += wt4[i] * Iv;
		}
	}
}

void DGVSIE::FillZVV()
{
	int t_fld, t_src;
	VectorR3 vx_fld, vx_src, v_fld_tet4[4], v_src_tet4[4], v_fld_tri3[3], v_src_tri3[3], v_fld_cen, v_src_cen;
	Tetrahedron tet_fld, tet_src;
	value_t a_fld, a_src;
	value_t eps_f, eps_s;

	KernelIntegral();

	Complex coef = Z0*Z0*k_ / (J0*PI4);
	int vx, v_ftri1, v_ftri2, v_ftri3, v_stri1, v_stri2, v_stri3;
	tool::BarAndPercent bar;
	for (int f = 0; f < unknowns_v; f++)
	{
		bar(f + 1, unknowns_v);
		ct_ptr_->gethSWG(f, v_ftri1, v_ftri2, v_ftri3, vx, t_fld, eps_f, a_fld);
		tet_fld = mesh_ptr_->getTetrahedronRef(t_fld);
		tet_fld.GetVertices(v_fld_tet4[0], v_fld_tet4[1], v_fld_tet4[2], v_fld_tet4[3]);
		v_fld_tri3[0] = mesh_ptr_->getVertex(v_ftri1);
		v_fld_tri3[1] = mesh_ptr_->getVertex(v_ftri2);
		v_fld_tri3[2] = mesh_ptr_->getVertex(v_ftri3);
		vx_fld = mesh_ptr_->getVertex(vx);
		v_fld_cen = (v_fld_tri3[0] + v_fld_tri3[1] + v_fld_tri3[2]) / 3.0f;

		for (int s = 0; s < unknowns_v; s++)
		{
			ct_ptr_->gethSWG(s, v_stri1, v_stri2, v_stri3, vx, t_src, eps_s, a_src);
			tet_src = mesh_ptr_->getTetrahedronRef(t_src);
			tet_src.GetVertices(v_src_tet4[0], v_src_tet4[1], v_src_tet4[2], v_src_tet4[3]);
			v_src_tri3[0] = mesh_ptr_->getVertex(v_stri1);
			v_src_tri3[1] = mesh_ptr_->getVertex(v_stri2);
			v_src_tri3[2] = mesh_ptr_->getVertex(v_stri3);
			vx_src = mesh_ptr_->getVertex(vx);
			v_src_cen = (v_src_tri3[0] + v_src_tri3[1] + v_src_tri3[2]) / 3.0f;
			Complex z_temp(0.0f, 0.0f);
			if ((t_fld == t_src) && ((v_fld_cen - v_src_cen).Norm() < 1e-6))
			{
				z_temp = FHZ1Singular(v_fld_tet4, v_src_tet4, vx_fld, vx_src, v_fld_tri3, v_src_tri3, t_fld, t_src, eps_f, eps_s);
			}
			else if (t_fld == t_src)
			{
				z_temp = FHZ2Singular(v_fld_tet4, v_src_tet4, vx_fld, vx_src, v_fld_tri3, v_src_tri3, t_fld, t_src, eps_f, eps_s);
			}
			else if ((v_fld_cen - v_src_cen).Norm() < 1e-6)
			{
				z_temp = FHZ3Singular(v_fld_tet4, v_src_tet4, vx_fld, vx_src, v_fld_tri3, v_src_tri3, t_fld, t_src, eps_f, eps_s);
			}
			else
			{
				z_temp = FHZppKernel(v_fld_tet4, v_src_tet4, vx_fld, vx_src, v_fld_tri3, v_src_tri3, t_fld, t_src, eps_f, eps_s);
			}

			Z(f, s) = a_fld * a_src * coef * z_temp;
		}
	}
}

void DGVSIE::FillZVS()
{
	int t_fld;
	int f_src_plu, f_src_min;
	VectorR3 v_fld, v_fld4[4], v_fld_ct3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
	Tetrahedron tet;
	Triangle tri_plu, tri_min;
	value_t a_fld, l_src,eps_f;

	int n_fld, n_ct1, n_ct2, n_ct3;
	int n_src_plu, n_src_min;
	int tri_plu1, tri_plu2, tri_plu3, tri_min1, tri_min2, tri_min3;
	Complex coef(0.0f, 0.0f);

	//填充ZVS矩阵
	tool::BarAndPercent bar_perc2;
	coef = J0*k_*Z0 / PI4;
	for (int f = 0; f < unknowns_v; f++)
	{
		bar_perc2(f + 1, unknowns_v);
		ct_ptr_->gethSWG(f, n_ct1, n_ct2, n_ct3, n_fld, t_fld, eps_f, a_fld);
		v_fld_ct3[0] = mesh_ptr_->getVertex(n_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(n_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(n_ct3);

		v_fld = mesh_ptr_->getVertex(n_fld);
		tet = mesh_ptr_->getTetrahedronRef(t_fld);
		tet.GetVertices(v_fld4[0], v_fld4[1], v_fld4[2], v_fld4[3]);

		for (int s = 0; s < unknowns_s; s++)
		{
			ce_ptr_->getCommonEdge(s, n_src_plu, n_src_min, f_src_plu, f_src_min, l_src);

			v_src_plu = mesh_ptr_->getVertex(n_src_plu);
			v_src_min = mesh_ptr_->getVertex(n_src_min);

			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
			tri_min = mesh_ptr_->getTriangleRef(f_src_min);

			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);
			tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
			tri_min.getVertex(tri_min1, tri_min2, tri_min3);
			//fill
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);

			if ((tri_plu1 == n_ct1&&tri_plu2 == n_ct2&&tri_plu3 == n_ct3)
				|| (tri_plu1 == n_ct1&&tri_plu2 == n_ct3&&tri_plu3 == n_ct2)
				|| (tri_plu1 == n_ct2&&tri_plu2 == n_ct1&&tri_plu3 == n_ct3)                  
				|| (tri_plu1 == n_ct2&&tri_plu2 == n_ct3&&tri_plu3 == n_ct1)
				|| (tri_plu1 == n_ct3&&tri_plu2 == n_ct1&&tri_plu3 == n_ct2)
				|| (tri_plu1 == n_ct3&&tri_plu2 == n_ct2&&tri_plu3 == n_ct1))
			{
				Zpp = ZVSppSingular(v_fld4, v_src_plu3, v_fld, v_src_plu, v_fld_ct3, eps_f);
			}
			else if ((tri_min1 == n_ct1&&tri_min2 == n_ct2&&tri_min3 == n_ct3)
				|| (tri_min1 == n_ct1&&tri_min2 == n_ct3&&tri_min3 == n_ct2)
				|| (tri_min1 == n_ct2&&tri_min2 == n_ct1&&tri_min3 == n_ct3)
				|| (tri_min1 == n_ct2&&tri_min2 == n_ct3&&tri_min3 == n_ct1)
				|| (tri_min1 == n_ct3&&tri_min2 == n_ct1&&tri_min3 == n_ct2)
				|| (tri_min1 == n_ct3&&tri_min2 == n_ct2&&tri_min3 == n_ct1))
			{
				Zpm = ZVSpmSingular(v_fld4, v_src_min3, v_fld, v_src_min, v_fld_ct3, eps_f);
			}
			else
			{
				Zpp = ZVSppKernel(v_fld4, v_src_plu3, v_fld, v_src_plu, v_fld_ct3, eps_f);
				Zpm = ZVSpmKernel(v_fld4, v_src_min3, v_fld, v_src_min, v_fld_ct3, eps_f);
			}
			Z(f, s + unknowns_v) = coef*a_fld*l_src*(Zpp + Zpm);
		}

	}

}

void DGVSIE::FillZSV()
{
	int t_src;
	int f_fld_plu, f_fld_min;
	VectorR3 v_fld_plu, v_fld_plu3[3], v_fld_min, v_fld_min3[3];
	VectorR3 v_src, v_src4[4], v_src_ct3[3];
	Tetrahedron tet;
	Triangle tri_plu, tri_min;
	value_t a_src, l_fld, eps_s;

	int n_src, n_ct1, n_ct2, n_ct3;
	int n_fld_plu, n_fld_min;
	int tri_plu1, tri_plu2, tri_plu3, tri_min1, tri_min2, tri_min3;
	Complex coef(0.0f, 0.0f);

	//填充ZSV
	tool::BarAndPercent bar_perc3;
	coef = Z0*Z0*k_ / (J0*PI4);

	for (int f = 0; f < unknowns_s; f++)
	{
		bar_perc3(f + 1, unknowns_s);
		ce_ptr_->getCommonEdge(f, n_fld_plu, n_fld_min, f_fld_plu, f_fld_min, l_fld);

		v_fld_plu = mesh_ptr_->getVertex(n_fld_plu);//场点+
		v_fld_min = mesh_ptr_->getVertex(n_fld_min);//场点-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);
		tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
		tri_min.getVertex(tri_min1, tri_min2, tri_min3);
		for (int s = 0; s < unknowns_v; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ct_ptr_->gethSWG(s, n_ct1, n_ct2, n_ct3, n_src, t_src, eps_s, a_src);
			tet = mesh_ptr_->getTetrahedronRef(t_src);
			tet.GetVertices(v_src4[0], v_src4[1], v_src4[2],v_src4[3]);
			v_src_ct3[0] = mesh_ptr_->getVertex(n_ct1);
			v_src_ct3[1] = mesh_ptr_->getVertex(n_ct2);
			v_src_ct3[2] = mesh_ptr_->getVertex(n_ct3);

			if ((tri_plu1 == n_ct1&&tri_plu2 == n_ct2&&tri_plu3 == n_ct3)
				|| (tri_plu1 == n_ct1&&tri_plu2 == n_ct3&&tri_plu3 == n_ct2)
				|| (tri_plu1 == n_ct2&&tri_plu2 == n_ct1&&tri_plu3 == n_ct3)
				|| (tri_plu1 == n_ct2&&tri_plu2 == n_ct3&&tri_plu3 == n_ct1)
				|| (tri_plu1 == n_ct3&&tri_plu2 == n_ct1&&tri_plu3 == n_ct2)
				|| (tri_plu1 == n_ct3&&tri_plu2 == n_ct2&&tri_plu3 == n_ct1))
			{
				Zpp = ZSVppSingular(v_fld_plu3, v_src4, v_fld_plu, v_src, v_src_ct3, eps_s);
			}
			else if ((tri_min1 == n_ct1&&tri_min2 == n_ct2&&tri_min3 == n_ct3)
				|| (tri_min1 == n_ct1&&tri_min2 == n_ct3&&tri_min3 == n_ct2)
				|| (tri_min1 == n_ct2&&tri_min2 == n_ct1&&tri_min3 == n_ct3)
				|| (tri_min1 == n_ct2&&tri_min2 == n_ct3&&tri_min3 == n_ct1)
				|| (tri_min1 == n_ct3&&tri_min2 == n_ct1&&tri_min3 == n_ct2)
				|| (tri_min1 == n_ct3&&tri_min2 == n_ct2&&tri_min3 == n_ct1))
			{
				Zmp = ZSVmpSingular(v_fld_plu3, v_src4, v_fld_plu, v_src, v_src_ct3, eps_s);
			}
			else
			{
				Zpp = ZSVppKernel(v_fld_plu3, v_src4, v_fld_plu, v_src, v_src_ct3, eps_s);
				Zmp = ZSVmpKernel(v_fld_plu3, v_src4, v_fld_plu, v_src, v_src_ct3, eps_s);
			}

			Z(f + unknowns_v, s) = coef*l_fld*a_src*(Zpp + Zmp);
		}
	}
}

void DGVSIE::FillZSS()
{
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu3[3], v_fld_min3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
	Triangle tri_plu, tri_min;
	value_t l_fld, l_src;//场与源的公共边长度

	Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm;
	Complex coef = (J0 * k_ * Z0) / PI4;
	//填充ZMM矩阵
	tool::BarAndPercent bar_perc4;
	//coef = (J0*k_*Z0) / PI4;
	for (int f = 0; f < unknowns_s; f++)
	{
		bar_perc4(f + 1, unknowns_s);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		v_fld_min = mesh_ptr_->getVertex(vm);//场点-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);
		for (int s = 0; s < unknowns_s; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src);

			v_src_plu = mesh_ptr_->getVertex(vp);
			v_src_min = mesh_ptr_->getVertex(vm);

			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
			tri_min = mesh_ptr_->getTriangleRef(f_src_min);

			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);

			//field face+ <-->  source face+
			if (f_fld_plu == f_src_plu)
			{
				Zpp = ZSSppSingular(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			else
			{
				Zpp = ZSSppKernel(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			//field face+ <--> source face-
			if (f_fld_plu == f_src_min)
			{
				Zpm = ZSSpmSingular(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
			}
			else
			{
				Zpm = ZSSpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
			}
			//field face- <--> source face+
			if (f_fld_min == f_src_plu)
			{
				Zmp = ZSSmpSingular(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
			}
			else
			{
				Zmp = ZSSmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
			}
			//field face- <--> source face-
			if (f_fld_min == f_src_min)
			{
				Zmm = ZSSmmSingular(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
			}
			else
			{
				Zmm = ZSSmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
			}
			Z(f + unknowns_v, s + unknowns_v) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void DGVSIE::fillV()
{
	//填充SWG的右边项
	VectorR3 v4[4], vx;
	int nv1, nv2, nv3, nvx, tid;
	value_t eps, area;
	Tetrahedron tet;

	tool::BarAndPercent bar_perc1;
	for (int u = 0; u < unknowns_v; u++)
	{
		bar_perc1(u + 1, unknowns_v);

		ct_ptr_->gethSWG(u, nv1, nv2, nv3, nvx, tid, eps, area);
		vx = mesh_ptr_->getVertex(nvx);
		tet = mesh_ptr_->getTetrahedronRef(tid);
		tet.GetVertices(v4[0], v4[1], v4[2], v4[3]);

		Complex vk = eHVKernel(v4, vx);
		V(u) = area*vk / 3.0f;
	}

	//填充RWG的右边项
	int nvp, nvm, fp, fm;
	VectorR3 vp, vm, vp3[3], vm3[3];
	value_t ln;
	Triangle tri_plu, tri_min;

	tool::BarAndPercent bar_perc2;   //
	for (int u = 0; u < unknowns_s; u++)
	{
		bar_perc2(u + 1, unknowns_s);    //

		ce_ptr_->getCommonEdge(u, nvp, nvm, fp, fm, ln);
		vp = mesh_ptr_->getVertex(nvp);
		vm = mesh_ptr_->getVertex(nvm);
		tri_plu = mesh_ptr_->getTriangleRef(fp);
		tri_min = mesh_ptr_->getTriangleRef(fm);
		tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);
		tri_min.getVertex(vm3[0], vm3[1], vm3[2]);

		Complex vk = eMVKernel(vp3, vm3, vp, vm);
		V(u) = 0.5f * ln * vk;
	}

	/*int nvp, nvm, fp, fm;
	VectorR3 vp, vm, vp4[4], vm4[4];
	value_t area, epsp = 1.0f, epsm = 1.0f;
	Tetrahedron tet_plu, tet_min;

	tool::BarAndPercent bar_perc;   //
	for (int u = 0; u < unknowns_; u++)
	{
		bar_perc(u + 1, unknowns_);    //

		ct_ptr_->getCommonTriangle(u, nvp, nvm, fp, fm, area, epsp, epsm);
		vp = mesh_ptr_->getVertex(nvp);
		tet_plu = mesh_ptr_->getTetrahedronRef(fp);
		tet_plu.GetVertices(vp4[0], vp4[1], vp4[2], vp4[3]);
		if (nvm != -1)
		{
			vm = mesh_ptr_->getVertex(nvm);
			tet_min = mesh_ptr_->getTetrahedronRef(fm);
			tet_min.GetVertices(vm4[0], vm4[1], vm4[2], vm4[3]);
		}

		Complex vk = eVKernel(vp4, vm4, vp, vm, nvm);
		V(u) = area * vk / 3.0f;
	}*/
}

void DGVSIE::FillHV()
{
	VectorR3 v4[4], vx;
	int nv1, nv2, nv3, nvx, tid;
	value_t eps, area;
	Tetrahedron tet;

	tool::BarAndPercent bar_perc;
	for (int u = 0; u < unknowns_; u++)
	{
		bar_perc(u + 1, unknowns_);

		ct_ptr_->gethSWG(u, nv1, nv2, nv3, nvx, tid, eps, area);
		vx = mesh_ptr_->getVertex(nvx);
		tet = mesh_ptr_->getTetrahedronRef(tid);
		tet.GetVertices(v4[0], v4[1], v4[2], v4[3]);

		Complex vk = eHVKernel(v4, vx);
		V(u) = area*vk / 3.0f;
	}
}

void DGVSIE::prepareVIE(value_t _eps1, value_t _mu1)
{
	eps1_ = _eps1;
	mu1_ = _mu1;
}

bool DGVSIE::readExcEdges(const Qstring & rad_file, Qstring& stateInfo)
{
	Qifstream radStream(rad_file, std::ios::in);
	if (radStream.fail())
	{
		stateInfo.assign("read " + rad_file + " failed");
		return false;
	}
	std::pair<int, int> edge;
	while (radStream >> edge.first >> edge.second)
	{
		--edge.first;
		--edge.second;
		exc_edge_.push_back(edge);
	}

	if (radStream.eof())
		return true;
	else
		return false;
}

void DGVSIE::readExcEdges(const Qstring & rad_file)
{
	Qifstream rad_stream(rad_file, std::ios::in);
	if (rad_stream.fail())
		throw component::FileError("fail loading rad file: " + rad_file);

	std::pair<int, int> edge;
	while (rad_stream >> edge.first >> edge.second)
	{
		--edge.first;
		--edge.second;
		exc_edge_.push_back(edge);
	}
	rad_stream.close();
}

void DGVSIE::radiateV()
{
	int nv1, nv2;

	tool::BarAndPercent bar_perc;
	/*for (int u = 0; u < unknowns_; ++u)
	{
	bar_perc(u + 1, unknowns_);

	ce_ptr_->getCommonEdge(u, nv1, nv2);
	auto length = ce_ptr_->getCommonEdgeLength(u);

	for (const auto& elem : exc_edge_)
	{
	if ((elem.first == nv1 && elem.second == nv2) || (elem.first == nv2 && elem.second == nv1))
	{
	if (elem.first == nv1)
	V(u) = length;
	else
	V(u) = -length;
	break;
	}
	}
	}*/
}

value_t DGVSIE::getBiRCS(const VectorR3& sca_k) const
{
	//介质体
	int nv1, nv2, nv3, nvx, tid;
	VectorR3 vx, v4[4];
	value_t S, eps;

	VectorR3 vg4[4];

	VectorR3 rou;
	Complex G0(0, 0);
	VectorC3 Es(0, 0, 0), Es_v(0, 0, 0), Es_s(0, 0, 0);

	for (int u = 0; u < unknowns_v; u++)
	{
		ct_ptr_->gethSWG(u, nv1, nv2, nv3, nvx, tid, eps, S);

		vx = mesh_ptr_->getVertex(nvx);
		const auto& tet = mesh_ptr_->getTetrahedronRef(tid);
		tet.GetVertices(v4[0], v4[1], v4[2], v4[3]);
		Gauss4Point(v4[0], v4[1], v4[2], v4[3], vg4);

		VectorC3 Esp(0, 0, 0);
		for (int g = 0; g < 4; g++)
		{
			rou = vg4[g] - vx;
			G0 = exp(J0 * k_ * (vg4[g] ^ sca_k));
			Esp += rou * wt4[g] * G0*(1.0f - (1.0f / eps));
		}
		Es_v += Esp * (I[u] * S);
	}

	//金属表面
	int nvp, nvm, fp, fm;
	value_t L;
	VectorR3 vp, vm, vp3[3], vm3[3], rou_p, rou_m;
	VectorR3 vgp_m[7], vgm_m[7];
	for (int u = 0; u < unknowns_s; u++)
	{
		ce_ptr_->getCommonEdge(u, nvp, nvm, fp, fm, L);
		vp = mesh_ptr_->getVertex(nvp);
		vm = mesh_ptr_->getVertex(nvm);
		const auto& tri_plu = mesh_ptr_->getTriangleRef(fp);
		const auto& tri_min = mesh_ptr_->getTriangleRef(fm);
		tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);
		tri_min.getVertex(vm3[0], vm3[1], vm3[2]);

		Gauss7Point(vp3[0], vp3[1], vp3[2], vgp_m);
		Gauss7Point(vm3[0], vm3[1], vm3[2], vgm_m);

		VectorC3 Esp(0, 0, 0), Esm(0, 0, 0);
		for (int g = 0; g < 7; g++)
		{
			rou_p = vgp_m[g] - vp;
			G0 = exp(J0 * k_ * (vgp_m[g] ^ sca_k));
			Esp = Esp + rou_p * w7[g] * G0;

			rou_m = vm - vgm_m[g];
			G0 = exp(J0 * k_ * (vgm_m[g] ^ sca_k));
			Esm = Esm + rou_m * w7[g] * G0;
		}
		Es_s += (Esp + Esm) * (I[u + unknowns_v] * L);
	}


	auto es = sca_k*(sca_k*((Z0*Es_v / 3.0f) + (Es_s / 2.0f)));
	value_t rcs = k_ * k_ * Z0*Z0* es.norm() / PI4;

	return 10 * log10(rcs);
}

value_t DGVSIE::getBiRCS_h(const VectorR3& sca_k) const
{
	/*int nvp, nvm, fp, fm;
	VectorR3 vp, vm, vp4[4], vm4[4];
	value_t S, epsp = 0.0f, epsm = 0.0f;

	VectorR3 vgp[5], vgm[5];

	VectorR3 rou_p, rou_m;
	Complex G0(0, 0);
	VectorC3 Es(0, 0, 0);

	for (int u = 0; u < unknowns_; u++)
	{
	ct_ptr_->getCommonTriangle(u, nvp, nvm, fp, fm, S, epsp, epsm);

	vp = mesh_ptr_->getVertex(nvp);
	const auto& tet_plu = mesh_ptr_->getTetrahedronRef(fp);
	tet_plu.GetVertices(vp4[0], vp4[1], vp4[2], vp4[3]);
	Gauss5Point(vp4[0], vp4[1], vp4[2], vp4[3], vgp);

	if (nvm != -1)
	{
	vm = mesh_ptr_->getVertex(nvm);
	const auto& tet_min = mesh_ptr_->getTetrahedronRef(fm);
	tet_min.GetVertices(vm4[0], vm4[1], vm4[2], vm4[3]);
	Gauss5Point(vm4[0], vm4[1], vm4[2], vm4[3], vgm);
	}

	VectorC3 Esp(0, 0, 0), Esm(0, 0, 0);
	for (int g = 0; g < 5; g++)
	{
	rou_p = vgp[g] - vp;
	G0 = exp(J0 * k_ * (vgp[g] ^ sca_k));
	Esp += rou_p * wt5[g] * G0*(1.0f - (1.0f / epsp));
	if (nvm != -1)
	{
	rou_m = vm - vgm[g];
	G0 = exp(J0 * k_ * (vgm[g] ^ sca_k));
	Esm += rou_m * wt5[g] * G0*(1.0f - (1.0f / epsm));
	}

	}
	Es += (Esp + Esm) * (I[u] * S);
	}*/
	int nv1, nv2, nv3, nvx, tid;
	VectorR3 vx, v4[4];
	value_t S, eps;

	VectorR3 vg4[4];

	VectorR3 rou;
	Complex G0(0, 0);
	VectorC3 Es(0, 0, 0);

	for (int u = 0; u < unknowns_; u++)
	{
		ct_ptr_->gethSWG(u, nv1, nv2, nv3, nvx, tid, eps, S);

		vx = mesh_ptr_->getVertex(nvx);
		const auto& tet = mesh_ptr_->getTetrahedronRef(tid);
		tet.GetVertices(v4[0], v4[1], v4[2], v4[3]);
		Gauss4Point(v4[0], v4[1], v4[2], v4[3], vg4);

		VectorC3 Esp(0, 0, 0);
		for (int g = 0; g < 4; g++)
		{
			rou = vg4[g] - vx;
			G0 = exp(J0 * k_ * (vg4[g] ^ sca_k));
			Esp += rou * wt4[g] * G0*(1.0f - (1.0f / eps));
		}
		Es += Esp * (I[u] * S);
	}

	auto es = sca_k * (sca_k * Es) / 6.0f;
	value_t rcs = k_ * k_ * omiga * omiga *Z0*Z0* es.norm() / PI;

	return 10 * log10(rcs);
}

void DGVSIE::getEFiled(component::FieldData * pdata, const VectorR3 & rad_k, const VectorR3 & rad_ev, const VectorR3 & rad_eh) const
{
	int nvp, nvm, fp, fm;
	VectorR3 vp, vm, vp3[3], vm3[3];
	value_t L;
	Triangle tri_plu, tri_min;

	VectorR3 vgp[7], vgm[7];

	VectorR3 rou_p, rou_m;
	Complex  G0(0, 0);
	VectorC3 Es(0, 0, 0);

	for (int u = 0; u < unknowns_; u++)
	{
		ct_ptr_->getCommonTriangle(u, nvp, nvm, fp, fm, L);
		vp = mesh_ptr_->getVertex(nvp);
		vm = mesh_ptr_->getVertex(nvm);
		tri_plu = mesh_ptr_->getTriangleRef(fp);
		tri_min = mesh_ptr_->getTriangleRef(fm);
		tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);
		tri_min.getVertex(vm3[0], vm3[1], vm3[2]);

		Gauss7Point(vp3[0], vp3[1], vp3[2], vgp);
		Gauss7Point(vm3[0], vm3[1], vm3[2], vgm);

		VectorC3 Esp(0, 0, 0), Esm(0, 0, 0);
		for (int g = 0; g < 7; g++)
		{
			rou_p = vgp[g] - vp;
			G0 = exp(J0 * k_ * (vgp[g] ^ rad_k));
			Esp = Esp + rou_p * w7[g] * G0;

			rou_m = vm - vgm[g];
			G0 = exp(J0 * k_ * (vgm[g] ^ rad_k));
			Esm = Esm + rou_m * w7[g] * G0;
		}
		Es = Es + (Esp + Esm) * (I[u] * L);
	}
	//auto dir_theta = sca_k_ * Es;
	const auto coeff = -J0 * Z0 * k_ / (2.0f * PI4);
	pdata->ftheta = coeff * (Es ^ rad_ev);
	pdata->fphi = coeff * (Es ^ rad_eh);
}

bool DGVSIE::writeZIVData()
{
	Qofstream outputZ(dir_ + "/matrix_Z.txt");
	if (outputZ.fail())
		return false;
	Z.save(outputZ, arma::arma_ascii);

	Qofstream outputV(dir_ + "/matrix_V.txt");
	if (outputV.fail())
		return false;
	V.save(outputV, arma::arma_ascii);

	Qofstream outputI(dir_ + "/martix_I.txt");
	if (outputI.fail())
		return false;
	I.save(outputI, arma::arma_ascii);

	return true;
}

