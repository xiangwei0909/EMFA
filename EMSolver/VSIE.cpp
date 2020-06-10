#include "VSIE.h"
#include "ConfigLoader.h"
#include "Mathdef.h"
#include "Quadrature.h"
#include "tools.h"
#include "CommonEdge.h"
#include "CommonTriangle.h"
#include "Triangle.h"
#include "Tetrahedron.h"
#include "iml.h"


using namespace component;
using namespace mom;
using namespace math;
using namespace Eigen;
using std::setw;


VSIE::VSIE()
{

}

VSIE::~VSIE()
{
}

void VSIE::init(component::ConfigLoader * ploader)
{
	EM::init(ploader);
	
	if (rtype_ == policy::ResultType::Rad)
		readExcEdges(dir_ + '/' + folder_name_ + ".rad");

	dir_ = tool::creatFolder(dir_, folder_name_ + "_VSIE");

	auto log_name = dir_ + (rtype_ == policy::ResultType::Rad ? "/rad_" : "/sca_") + "runtime.log";
	logger_.open(log_name);

	mesh_ptr_ = std::make_shared<Mesh>();
	ce_ptr_ = std::make_shared<CommonEdge>();
	ct_ptr_ = std::make_shared<CommonTriangle>();

	mesh_ptr_->loadVSIEMesh(ploader->getMeshPath());
	ce_ptr_->buildCommonEdge(mesh_ptr_);
	ct_ptr_->buildCommonTriangle(mesh_ptr_);

	unknowns_t = ct_ptr_->getCommonTriangleNum();
	unknowns_e = ce_ptr_->getCommonEdgeNum();
	unknowns_ = unknowns_t + unknowns_e;
	tet_num = mesh_ptr_->getTetrahedronNum();

	k_ = PI2 * incidence_.freq / cc;
	omiga = PI2*incidence_.freq;
	Sigma = ploader->getSigma();

	iter_num = ploader->getMaxIterationNum();
	iter_threshold = ploader->getIterationThreshold();

	//prepareVIE(ploader->getEps1(), ploader->getMu());

	logger_ << ploader->getReport();
	reportInfo(logger_);
	TIME_LOG("init");
}

void VSIE::solve()
{
	SEGMENT("Solve");
	Z.zeros(unknowns_, unknowns_);
	Z_eigen.resize(unknowns_, unknowns_);

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
	TIME_LOG("fillZ");

	I.zeros(unknowns_);
	I_eigen.resize(unknowns_);
	V.zeros(unknowns_);
	V_eigen.resize(unknowns_);

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
	Qcout << setw(30) << "Solving matrix equation:" << std::endl;

	/*for (int m = 0; m < unknowns_; m++)
	{
		for (int n = 0; n < unknowns_; n++)
		{
			Z_eigen(m, n) = Z(m, n);
		}
		V_eigen[m] = V[m];
	}*/

	I_eigen = (Z_eigen.lu().solve(V_eigen));
	
	//if (!arma::solve(I, Z, V, arma::solve_opts::fast))
		Qcout << setw(30) << "Finish the solving" << std::endl;
	//if(!iml::CGNE(Z,V,I,iter_threshold,iter_num,useless))
	
	Qcout << "success" << std::endl;
	TIME_LOG("solve");
	if (!writeZIVData())
		Qcout << "The data is already saved!" << std::endl;

#ifdef _DEBUG
	//if (!writeZIVData())
	//  Qcerr << "fail writing ZIV data" << std::endl;
#endif
}

void VSIE::output()
{
	SEGMENT("Result");
	Result result;
	if (rtype_ == policy::ResultType::Rad)
	{
		LOG(result.getRadiationField(this, &VSIE::getEFiled), "Calculating Radiation");
		TIME_LOG("getRadiation");
	}
	else
	{
		LOG(result.getBistaticRCS(this, &VSIE::getBiRCS), "Calculating BistaticRCS");
		TIME_LOG("getBistaticRCS");
	}
	LOG(result.getNearField(this, &VSIE::getNearEField), "Calculating Near field");
	LOG(result.getCurrentDistribution(this, &VSIE::calculateSurfaceCurrent), "Calculating surface current");
}

void VSIE::clear()
{
	Z.reset();
	I.reset();
	V.reset();
	exc_edge_.clear();

	mesh_ptr_->clear();
	ct_ptr_->clear();
}


void VSIE::reportInfo(Qostream& strm) const
{
	mesh_ptr_->reportInfo(strm);
	ct_ptr_->reportInfo(strm);
}

Complex VSIE::DDZppKernel(VectorR3 * vf4, VectorR3 * vs4, VectorR3 & vf, VectorR3 & vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_p, value_t &epsp, value_t &epsm)
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
	//�����岻�غϣ�����Z1=0

	//z2,z5
	Z2 = Z2K1(t_f_p, t_s_p) - (vf^Z2K2(t_f_p, t_s_p)) - (vs^Z2K3(t_f_p, t_s_p)) + ((vs^vf)*Z2K4(t_f_p, t_s_p));
	Z5 = Z5K(t_f_p, t_s_p);

	//z6
	for (int p = 0; p < 4; p++)
	{
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

Complex VSIE::DDZpmKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_m, value_t &epsp, value_t &epsm)
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
	//�����岻�غϣ�����Z1=0

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

Complex VSIE::DDZmpKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_p, value_t &epsp, value_t &epsm)
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
	//�����岻�غϣ�����Z1=0

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

Complex VSIE::DDZmmKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_m, value_t &epsp, value_t &epsm)
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
	//�����岻�غϣ�����Z1=0

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

Complex VSIE::DDZppSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_p, value_t &epsp, value_t &epsm)
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

	//�����岻�غϣ�����Z1=0
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
			//z6�����첿��
			for (int q = 0; q < 3; q++)
			{
				r = (vg4f[p] - vg3s[q]).Norm();
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
				Z6 += wt4[p] * w3[q] * G;
			}
			//z6���첿��
			Z6 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
		}
	}
	//z3,z4�����첿��
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
			//z3���첿��
			value_t Iv = 0.0f;
			for (int j = 0; j < 4; j++)
			{
				Iv += 0.5f*IvSingular(v4_3[j], vg3f[p]);
			}
			Z3 += w3[p] * Iv / volume;

			//z4������
			if (abs(epsp - epsm) > 1e-4)
			{
				for (int q = 0; q < 3; q++)
				{
					r = (vg3f[p] - vg3s[q]).Norm();
					G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
					Z4 += w3[p] * w3[q] * G;
				}
				//z4����
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

Complex VSIE::DDZpmSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_m, value_t &epsp, value_t &epsm)
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

	//z3,z4�����첿��
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
			//z3���첿��
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

Complex VSIE::DDZmpSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_p, value_t &epsp, value_t &epsm)
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

	//�����岻�غϣ�����Z1=0
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
			//z6�����첿��
			for (int q = 0; q < 3; q++)
			{
				r = (vg4f[p] - vg3s[q]).Norm();
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
				Z6 += wt4[p] * w3[q] * G;
			}
			//z6���첿��
			value_t Is = 0.0f;
			Z6 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
		}
	}
	//z3,z4�����첿��

	Z = (-1.0f / k_) *(PI4 / (9.0f*volume*epsp))*Z1 + (k_*kp / 9.0f)*Z2 + (1.0f / k_) *(-kp*Z5 - (km - kp)*Z6);

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex VSIE::DDZmmSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_m, value_t &epsp, value_t &epsm)
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

	//�����岻�غϣ�����Z1=0
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


Complex VSIE::DMZppKernel(VectorR3 * vf4, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 * vf3, int & _vm_f)
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
			Z1 += wt4[p]*w3[q] * (rou_f^rou_s)*G;
			Z2 += wt4[p]*w3[q] * G;
		}
	}

	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 3; q++)
			{
				r = (vg3f[p] - vg3s[q]).Norm();
				G = exp(-J0*k_*r) / r;
				Z3 += w3[p] * w3[q] * G;
			}
		}
	}
	
	Z = (1.0f / 6.0f)*Z1 - (1.0f / (k_*k_))*(Z2 - Z3);//�������J0*k_*Z0/PI4

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}

	return Z;
}

Complex VSIE::DMZpmKernel(VectorR3 * vf4, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 * vf3, int & _vm_f)
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
			Z1 += wt4[p]*w3[q] * (rou_f^rou_s)*G;
			Z2 += wt4[p]*w3[q] * G;
		}
	}

	if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 3; q++)
			{
				r = (vg3f[p] - vg3s[q]).Norm();
				G = exp(-J0*k_*r) / r;
				Z3 += w3[p] * w3[q] * G;
			}
		}
	}


	Z = -(1.0f / 6.0f)*Z1 + (1.0f / (k_*k_))*(Z2 - Z3);//�������J0*k_*Z0

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex VSIE::DMZmpKernel(VectorR3 * vf4, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 * vf3, int & _vm_f)
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
			Z1 += wt4[p]*w3[q] * (rou_f^rou_s)*G;
			Z2 += wt4[p]*w3[q] * G;
		}
	}

	/*if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 3; q++)
			{
				r = (vg1f[p] - vg3s[q]).Norm();
				G = exp(-J0*k_*r) / r;
				Z3 += w3[p] * w3[q] * G;
			}
		}
	}*/

	Z = -(1.0f / 6.0f)*Z1 + (1.0f / (k_*k_))*(Z2);//�������J0*k_*Z0

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex VSIE::DMZmmKernel(VectorR3 * vf4, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 * vf3, int & _vm_f)
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
			Z1 += wt4[p]*w3[q] * (rou_f^rou_s)*G;
			Z2 += wt4[p]*w3[q] * G;
		}
	}

	/*if (_vm_f == -1)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 3; q++)
			{
				r = (vg1f[p] - vg3s[q]).Norm();
				G = exp(-J0*k_*r) / r;
				Z3 += w3[p] * w3[q] * G;
			}
		}
	}*/

	Z = (1.0f / 6.0f)*Z1 - (1.0f / (k_*k_))*(Z2);//�������J0*k_*Z0

	return Z;
}


Complex	VSIE::DMZppSingular(VectorR3 * vf4, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 * vf3, int & _vm_f)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg3f[3], vg3s[3], rou_f, rou_s;
	value_t r,area;
	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	area = Area(vs3[0], vs3[1], vs3[2]);

	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//�����첿��
		for (int q = 0; q < 3; q++)
		{
			rou_s = vg3s[q] - vs;
			r = (vg4f[p] - vg3s[q]).Norm();
			G = exp(-J0*k_*r) / r;
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z1 += wt4[p]*w3[q] * (rou_f^rou_s)*G;
			Z2 += wt4[p]*w3[q] * G;
		}
		//���첿��
		//Z1 += wt4[p] * (rou_f ^ (IspSingular(vs3, vg4f[p]) + (vg4f[p] - vs)*IsSingular(vs3, vg4f[p]))) / area;
		//Z2 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
	}

	for (int p = 0; p < 3; p++)
	{
		//Z3�����첿��
		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z3 += w3[p] * w3[q] * G;
		}
		//Z3���첿��
		Z3 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}

	Z = (1.0f / 6.0f)*Z1 - (1.0f / (k_*k_))*(Z2 - Z3);//�������J0*k_*Z0/PI4

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}

	return Z;
}

Complex	VSIE::DMZpmSingular(VectorR3 * vf4, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 * vf3, int & _vm_f)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg4f[4], vg3f[3], vg3s[3], rou_f, rou_s;
	value_t r, area;

	Gauss4Point(vf4[0], vf4[1], vf4[2], vf4[3], vg4f);
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	area = Area(vs3[0], vs3[1], vs3[2]);

	for (int p = 0; p < 4; p++)
	{
		rou_f = vg4f[p] - vf;
		//�����첿��
		for (int q = 0; q < 3; q++)
		{
			rou_s = vg3s[q] - vs;
			r = (vg4f[p] - vg3s[q]).Norm();
			G = exp(-J0*k_*r) / r;
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z1 += wt4[p] * w3[q] * (rou_f^rou_s)*G;
			Z2 += wt4[p] * w3[q] * G;
		}
		//���첿��
		//Z1 += wt4[p] * (rou_f ^ (IspSingular(vs3, vg4f[p]) + (vg4f[p] - vs)*IsSingular(vs3, vg4f[p]))) / area;
		//Z2 += wt4[p] * IsSingular(vs3, vg4f[p]) / area;
	}

	for (int p = 0; p < 3; p++)
	{
		//Z3�����첿��
		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z3 += w3[p] * w3[q] * G;
		}
		//Z3���첿��
		Z3 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}

	Z = -(1.0f / 6.0f)*Z1 + (1.0f / (k_*k_))*(Z2 - Z3);//�������J0*k_*Z0/PI4

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}

	return Z;
}


Complex VSIE::MDZppKernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg3f[3], vg3s[3], vg4s[4], rou_f, rou_s;
	value_t r, kp, km;
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	kp = 1.0f - (1.0f / epsp);
	km = 1.0f - (1.0f / epsm);
	//kp = epsp - 1;
	//km = epsm - 1;
	for (int p = 0; p < 3; p++)
	{
		rou_f = vg3f[p] - vf;
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg3f[p] - vg4s[q]).Norm();
			G = exp(-J0*k_*r) / r;
			Z1 += w3[p] * wt4[q] * (rou_f^rou_s*G);
			Z2 += w3[p] * wt4[q] * G;//in the end need *kp
		}
	}

	if (abs(epsm-epsp)>1e-6)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 3; q++)
			{
				r = (vg3f[p] - vg3s[q]).Norm();
				G = exp(-J0*k_*r) / r;
				Z3 += w3[p] * w3[q] * G;//in the end need *(km-kp)
			}
		}
	}
	
	Z = (1.0f / 6.0f)*(-kp*Z1 + (6.0f / (k_*k_))*(kp*Z2 + (km-kp)*Z3));

	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex VSIE::MDZpmKernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg3f[3], vg3s[3], vg4s[4], rou_f, rou_s;
	value_t r, kp, km;
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	kp = 1.0f - (1.0f / epsp);
	km = 1.0f - (1.0f / epsm);
	//kp = epsp - 1;
	//km = epsm - 1;
	for (int p = 0; p < 3; p++)
	{
		rou_f = vg3f[p] - vf;
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg3f[p] - vg4s[q]).Norm();
			G = exp(-J0*k_*r) / r;
			Z1 += w3[p] * wt4[q] * (rou_f^rou_s*G);
			Z2 += w3[p] * wt4[q] * G;//in the end need *kp
		}
	}

	/*if (abs(epsm - epsp)<1e-6)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 3; q++)
			{
				r = (vg3f[p] - vg5s[q]).Norm();
				G = exp(-J0*k_*r) / r;
				Z3 += w3[p] * w3[q] * G;//in the end need *(km-kp)
			}
		}
	}*/

	Z = (1.0f / 6.0f)*(km*Z1 - (6.0f / (k_*k_))*(km*Z2));
	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex VSIE::MDZmpKernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg3f[3], vg3s[3], vg4s[4], rou_f, rou_s;
	value_t r, kp, km;
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	kp = 1.0f - (1.0f / epsp);
	km = 1.0f - (1.0f / epsm);
	//kp = epsp - 1;
	//km = epsm - 1;
	for (int p = 0; p < 3; p++)
	{
		rou_f = vg3f[p] - vf;
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg3f[p] - vg4s[q]).Norm();
			G = exp(-J0*k_*r) / r;
			Z1 += w3[p] * wt4[q] * (rou_f^rou_s*G);
			Z2 += w3[p] * wt4[q] * G;//in the end need *kp
		}
	}

	if (abs(epsm - epsp)>1e-6)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 3; q++)
			{
				r = (vg3f[p] - vg3s[q]).Norm();
				G = exp(-J0*k_*r) / r;
				Z3 += w3[p] * w3[q] * G;//in the end need *(km-kp)
			}
		}
	}

	Z = (1.0f / 6.0f)*(kp*Z1 - (6.0f / (k_*k_))*(kp*Z2 + (km-kp)*Z3));
	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex VSIE::MDZmmKernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg3f[3], vg3s[3], vg4s[4], rou_f, rou_s;
	value_t r, kp, km;
	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	kp = 1.0f - (1.0f / epsp);
	km = 1.0f - (1.0f / epsm);
	//kp = epsp - 1;
	//km = epsm - 1;
	for (int p = 0; p < 3; p++)
	{
		rou_f = vg3f[p] - vf;
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg3f[p] - vg4s[q]).Norm();
			G = exp(-J0*k_*r) / r;
			Z1 += w3[p] * wt4[q] * (rou_f^rou_s*G);
			Z2 += w3[p] * wt4[q] * G;//in the end need *kp
		}
	}

	/*if (abs(epsm - epsp)<1e-6)
	{
		for (int p = 0; p < 3; p++)
		{
			for (int q = 0; q < 3; q++)
			{
				r = (vg3f[p] - vg5s[q]).Norm();
				G = exp(-J0*k_*r) / r;
				Z3 += w3[p] * w3[q] * G;//in the end need *(km-kp)
			}
		}
	}*/


	Z = (1.0f / 6.0f)*(-km*Z1 + (6.0f / (k_*k_))*(km*Z2));
	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}


Complex VSIE::MDZppSingular(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg3f[3], vg3s[3], vg4s[4], rou_f, rou_s;
	value_t r, kp, km;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	kp = 1.0f - (1.0f / epsp);
	km = 1.0f - (1.0f / epsm);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);
	value_t volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vg3f[p] - vf;
		//�����첿��
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg3f[p] - vg4s[q]).Norm();
			G = exp(-J0 * k_*r) / r;
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z1 += w3[p] * wt4[q] * (rou_f^rou_s)*G;
			Z2 += w3[p] * wt4[q] * G;//in the end need *kp
		}
		//���첿��
		/*VectorR3 Ivp(0.0f, 0.0f, 0.0f);
		value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
			Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg3f[p]);
			Iv += 0.5f*IvSingular(v4_3[j], vg3f[p]);
		}
		Z1 += w3[p] * (rou_f ^ (Ivp + Iv*(vg3f[p] - vs))) / volume;
		Z2 += w3[p] * (Iv / volume);*/
	}

	value_t area = Area(vs3[0], vs3[1], vs3[2]);
	for (int p = 0; p < 3; p++)
	{
		//Z3�����첿��
		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z3 += w3[p] * w3[q] * G;//in the end need *(km-kp)
		}
		//Z3���첿��
		Z3 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}

	/*Complex Zi(0, 0);
	if (abs(Sigma) > 1e-8)
	{
		for (int i = 0; i < 3; i++)
		{
			rou_f = vg3s[i] - vf;
			rou_s = vg3s[i] - vs;
			Zi += w3[i] * rou_f^rou_s;
		}
		Zi = (Sigma*PI4 / (volume*Z0*k_*omiga))*Zi;
	}*/

	Z = (1.0f / 6.0f)*(-kp * Z1 + (6.0f / (k_*k_))*(kp*Z2 + (km - kp)*Z3));
	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}

Complex VSIE::MDZmpSingular(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm)
{
	Complex Z(0.0f, 0.0f), Z1(0.0f, 0.0f), Z2(0.0f, 0.0f), Z3(0.0f, 0.0f), G(0.0f, 0.0f);
	VectorR3 vg3f[3], vg3s[3], vg4s[4], rou_f, rou_s;
	value_t r, kp, km;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	Gauss3Point(vf3[0], vf3[1], vf3[2], vg3f);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vg3s);
	Gauss4Point(vs4[0], vs4[1], vs4[2], vs4[3], vg4s);
	kp = 1.0f - (1.0f / epsp);
	km = 1.0f - (1.0f / epsm);

	Pv4_3[0] = v4_3[0];
	Pv4_3[1] = v4_3[1];
	Pv4_3[2] = v4_3[2];
	Pv4_3[3] = v4_3[3];
	TetToTri(vs4, Pv4_3);
	value_t volume = Volume(vs4[0], vs4[1], vs4[2], vs4[3]);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vg3f[p] - vf;
		//�����첿��
		for (int q = 0; q < 4; q++)
		{
			rou_s = vg4s[q] - vs;
			r = (vg3f[p] - vg4s[q]).Norm();
			G = exp(-J0*k_*r) / r;
			//G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z1 += w3[p] * wt4[q] * (rou_f^rou_s)*G;
			Z2 += w3[p] * wt4[q] * G;//in the end need *kp
		}
		//���첿��
		/*VectorR3 Ivp(0.0f, 0.0f, 0.0f);
		value_t Iv = 0.0f;
		for (int j = 0; j < 4; j++)
		{
			Ivp += (1.0f / 3.0f)*IvpSingular(v4_3[j], vg3f[p]);
			Iv += 0.5f*IvSingular(v4_3[j], vg3f[p]);
		}
		Z1 += w3[p] * (rou_f ^ (Ivp + Iv*(vg3f[p] - vs))) / volume;
		Z2 += w3[p] * (Iv / volume);*/
	}

	value_t area = Area(vs3[0], vs3[1], vs3[2]);
	for (int p = 0; p < 3; p++)
	{
		//Z3�����첿��
		for (int q = 0; q < 3; q++)
		{
			r = (vg3f[p] - vg3s[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Z3 += w3[p] * w3[q] * G;//in the end need *(km-kp)
		}
		//Z3���첿��
		Z3 += w3[p] * IsSingular(vs3, vg3f[p]) / area;
	}

	/*Complex Zi(0, 0);
	if (abs(Sigma) > 1e-8)
	{
		for (int i = 0; i < 3; i++)
		{
			rou_f = vg3s[i] - vf;
			rou_s = vg3s[i] - vs;
			Zi += w3[i] * rou_f^rou_s;
		}
		Zi = (Sigma*PI4 / (volume*Z0*k_*omiga))*Zi;
	}*/

	Z = (1.0f / 6.0f)*(kp*Z1 - (6.0f / (k_*k_))*(kp*Z2 + (km - kp)*Z3));
	if (_isnan(Z.real()))
	{
		Qcout << "This number is not a number!" << std::endl;
	}
	return Z;
}


Complex VSIE::MMZppKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex VSIE::MMZpmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex VSIE::MMZmpKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex VSIE::MMZmmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex VSIE::MMZppSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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
		//�糡  �����첿��
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs; //.2.
			value_t r = (vgf[p] - vgs[q]).Norm();
			Complex coff = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * coff; //.3.
		}
		//�糡  ���첿��
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
		//�糡  �����첿��
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs; //.2.
			value_t r = (vgf[p] - vgs[q]).Norm();
			Complex coff = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * coff; //.3.
		}
		//�糡  ���첿��
		Complex Ztemp(0, 0);
		value_t Is2(0.0f);
		VectorR3 Is1(0.0f, 0.0f, 0.0f);
		VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[0]);
		VectorR3 ni = n.Normalize();
		VectorR3 rou = vgf[p] - ni * (ni^vgf[p]);
		VectorR3 rou_n = vs - ni * (ni^vs);
		for (int i = 0; i < 3; i++)
		{
			VectorR3 li = vs3[(i + 1) % 3] - vs3[i];//not normalize
			VectorR3 lj = vs3[(i + 2) % 3] - vs3[(i + 1) % 3];//not normalize
			VectorR3 r_f = vgf[p] - vs3[i];
			VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
			VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
			VectorR3 u = li * n;
			VectorR3 ui = u.Normalize();//normalize
			VectorR3 li_1 = li.Normalize();//normalize
			value_t p0 = abs((vs3[i] - vgf[p]) ^ ui);
			value_t d = (vgf[p] - vs3[i]) ^ ni;
			value_t lp = (vs3[(i + 1) % 3] - vgf[p]) ^ li_1;
			value_t lm = (vs3[i] - vgf[p]) ^ li_1;
			value_t rp = (vgf[p] - vs3[(i + 1) % 3]).Norm();
			value_t rm = (vgf[p] - vs3[i]).Norm();

			if (p0 < 1e-10)
				continue;
			VectorR3 pi = ((rou_p - rou) - lp * li_1).Normalize();
			value_t r02 = p0 * p0 + d * d;
			Is1 += (0.5f*ui)*(r02*log((rp + lp) / (rm + lm)) + (lp*rp) - (lm*rm));
			Is2 += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
		}
		Ze2 += w3[p] * (0.25f * (rou_f ^ (Is1 + Is2 * (vgf[p] - vs))) - (Is2 / (k_ * k_)));
	}

	Complex Ze3(0, 0);
	if (abs(Sigma) > 1e-8)
	{
		for (int i = 0; i < 3; i++)
		{
			rou_f = vgf[i] - vf;
			rou_s = vgf[i] - vs;
			Ze3 += w3[i] * rou_f^rou_s;
		}
		Ze3 = (PI *Sigma / (S_s*J0*k_*Z0))*Ze3;
	}
	Ze = (Ze1 + 1.0f / S_s * Ze2) +Ze3;

	return Ze;
}

Complex VSIE::MMZpmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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
		//�糡  �����첿��
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vs - vgs[q]; //.2.
			value_t r = (vgf[p] - vgs[q]).Norm();
			Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
			Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * coff; //.3.
		}
		//�糡  ���첿��

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
		//�糡  �����첿��
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vs - vgs[q]; //.2.
			value_t r = (vgf[p] - vgs[q]).Norm();
			Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
			Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * coff; //.3.
		}
		//�糡  ���첿��

		Complex Ztemp(0.0, 0.0);
		value_t Is2(0.0f);
		VectorR3 Is1(0.0f, 0.0f, 0.0f);
		VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[0]);
		VectorR3 ni = n.Normalize();
		VectorR3 rou = vgf[p] - ni * (ni^vgf[p]);
		VectorR3 rou_n = vs - ni * (ni^vs);
		for (int i = 0; i < 3; i++)
		{
			VectorR3 li = vs3[(i + 1) % 3] - vs3[i];//not normalize
			VectorR3 lj = vs3[(i + 2) % 3] - vs3[(i + 1) % 3];//not normalize
			VectorR3 r_f = vgf[p] - vs3[i];
			VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
			VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
			//VectorR3 n = (li*lj);
			VectorR3 u = li * n;
			//VectorR3 ni = n / n.Norm();//normalize
			VectorR3 ui = u / u.Norm();//normalize
			VectorR3 li_1 = li / li.Norm();//normalize
			value_t p0 = abs((vs3[i] - vgf[p]) ^ ui);
			value_t d = (vgf[p] - vs3[i]) ^ ni;
			value_t lp = (rou_p - rou) ^ li_1;
			value_t lm = (rou_m - rou) ^ li_1;
			value_t rp = sqrt(p0*p0 + d * d + lp * lp);
			value_t rm = sqrt(p0*p0 + d * d + lm * lm);

			if (p0 < 1e-10)
				continue;
			VectorR3 pi = ((rou_p - rou) - lp * li_1).Normalize();
			value_t r02 = p0 * p0 + d * d;
			Is1 += (0.5f*ui)*(r02*log((rp + lp) / (rm + lm)) + (lp*rp) - (lm*rm));
			Is2 += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
		}
		Ze2 += w3[p] * (0.25f * (-1.0f*rou_f ^ (Is1 + Is2 * (vgf[p] - vs))) + (Is2 / (k_ * k_)));
	}

	Complex Ze3(0, 0);
	if (abs(Sigma) > 1e-8)
	{
		for (int i = 0; i < 3; i++)
		{
			rou_f = vgf[i] - vf;
			rou_s = vgf[i] - vs;
			Ze3 += w3[i] * rou_f^rou_s;
		}
		Ze3 = (PI *Sigma / (S_s*J0*k_*Z0))*Ze3;
	}

	Ze = (Ze1 + 1.0f / S_s * Ze2) - Ze3;

	return Ze;
}

Complex VSIE::MMZmpSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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
		//�糡  �����첿��
		rou_f = vf - vgf[p];
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs; //.2.
			value_t r = (vgf[p] - vgs[q]).Norm();
			Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
			Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * coff; //.3.
		}
		//�糡  ���첿��

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
		//�糡  �����첿��
		rou_f = vf - vgf[p];
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs; //.2.
			value_t r = (vgf[p] - vgs[q]).Norm();
			Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
			Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * coff; //.3.
		}
		//�糡  ���첿��

		Complex Ztemp(0.0, 0.0);
		value_t Is2(0.0f);
		VectorR3 Is1(0.0f, 0.0f, 0.0f);
		VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[0]);
		VectorR3 ni = n.Normalize();
		VectorR3 rou = vgf[p] - ni * (ni^vgf[p]);
		VectorR3 rou_n = vs - ni * (ni^vs);
		for (int i = 0; i < 3; i++)
		{
			VectorR3 li = vs3[(i + 1) % 3] - vs3[i];//not normalize
			VectorR3 lj = vs3[(i + 2) % 3] - vs3[(i + 1) % 3];//not normalize
			VectorR3 r_f = vgf[p] - vs3[i];
			VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
			VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
			VectorR3 u = li * n;
			VectorR3 ui = u / u.Norm();//normalize
			VectorR3 li_1 = li / li.Norm();//normalize
			value_t p0 = abs((vs3[i] - vgf[p]) ^ ui);
			value_t d = (vgf[p] - vs3[i]) ^ ni;
			value_t lp = (rou_p - rou) ^ li_1;
			value_t lm = (rou_m - rou) ^ li_1;
			value_t rp = sqrt(p0*p0 + d * d + lp * lp);
			value_t rm = sqrt(p0*p0 + d * d + lm * lm);

			if (p0 < 1e-10)
				continue;
			VectorR3 pi = ((rou_p - rou) - lp * li_1).Normalize();
			value_t r02 = p0 * p0 + d * d;
			Is1 += (0.5f*ui)*(r02*log((rp + lp) / (rm + lm)) + (lp*rp) - (lm*rm));
			Is2 += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
		}
		Ze2 += w3[p] * (0.25f * (rou_f ^ (Is1 + Is2 * (vgf[p] - vs))) + (Is2 / (k_ * k_)));
	}

	Complex Ze3(0, 0);
	if (abs(Sigma) > 1e-8)
	{
		for (int i = 0; i < 3; i++)
		{
			rou_f = vgf[i] - vf;
			rou_s = vgf[i] - vs;
			Ze3 += w3[i] * rou_f^rou_s;
		}
		Ze3 = (PI *Sigma / (S_s*J0*k_*Z0))*Ze3;
	}

	Ze = (Ze1 + 1.0f / S_s * Ze2) - Ze3;

	return Ze;
}

Complex VSIE::MMZmmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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
		//�糡  �����첿��
		rou_f = vf - vgf[p];
		for (int q = 0; q < 3; q++)
		{
			rou_s = vs - vgs[q]; //.2.
			value_t r = (vgf[p] - vgs[q]).Norm();
			Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
			Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * coff; //.3.
		}
		//�糡  ���첿��

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
		//�糡  �����첿��
		rou_f = vf - vgf[p];
		for (int q = 0; q < 3; q++)
		{
			rou_s = vs - vgs[q]; //.2.
			value_t r = (vgf[p] - vgs[q]).Norm();
			Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
			Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * coff; //.3.
		}
		//�糡  ���첿��

		Complex Ztemp(0.0, 0.0);
		value_t Is2(0.0f);
		VectorR3 Is1(0.0f, 0.0f, 0.0f);
		VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[0]);
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
			VectorR3 u = li * n;
			VectorR3 ui = u / u.Norm();//normalize
			VectorR3 li_1 = li / li.Norm();//normalize
			value_t p0 = abs((vs3[i] - vgf[p]) ^ ui);
			value_t d = (vgf[p] - vs3[i]) ^ ni;
			value_t lp = (rou_p - rou) ^ li_1;
			value_t lm = (rou_m - rou) ^ li_1;
			value_t rp = sqrt(p0*p0 + d * d + lp * lp);
			value_t rm = sqrt(p0*p0 + d * d + lm * lm);

			if (p0 < 1e-10)
				continue;
			VectorR3 pi = ((rou_p - rou) - lp * li_1).Normalize();
			value_t r02 = p0 * p0 + d * d;
			Is1 += (0.5f*ui)*(r02*log((rp + lp) / (rm + lm)) + (lp*rp) - (lm*rm));
			Is2 += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
		}
		Ze2 += w3[p] * (0.25f * (-1.0f*rou_f ^ (Is1 + Is2 * (vgf[p] - vs))) - (Is2 / (k_ * k_)));
	}

	Complex Ze3(0, 0);
	if (abs(Sigma) > 1e-8)
	{
		for (int i = 0; i < 3; i++)
		{
			rou_f = vgf[i] - vf;
			rou_s = vgf[i] - vs;
			Ze3 += w3[i] * rou_f^rou_s;
		}
		Ze3 = (PI *Sigma / (S_s*J0*k_*Z0))*Ze3;
	}

	Ze = (Ze1 + 1.0f / S_s * Ze2) +Ze3;

	return Ze;
}

value_t VSIE::IsSingular(VectorR3 * vs3, VectorR3 & vf)
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

VectorR3 VSIE::IspSingular(VectorR3 * vs3, VectorR3 & vf)
{
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[1]);
	VectorR3 ni = n.Normalize();
	VectorR3 rou = vf - ni*(ni^vf);

	VectorR3 Isp(0, 0, 0);
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
		//Isp += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
		Isp += (r02*(log((rp + lp) / (rm + lm))) + rp*lp - rm*lm)*ui;
	}
	return 0.5f*Isp;
}

value_t VSIE::IvSingular(VectorR3 * vs3, VectorR3 & vf)
{
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[0]);
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

VectorR3 VSIE::IvpSingular(VectorR3 * vs3, VectorR3 & vf)
{
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[0]);
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

void VSIE::TetToTri(VectorR3 * v4, VectorR3 ** v4_3)
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

Complex VSIE::z1(VectorR3 * vf4, VectorR3 * vs4, VectorR3 & vf, VectorR3 & vs, value_t & epsp, value_t & epsm)
{

	return Complex();
}

Complex VSIE::DVKernel(VectorR3 *vp4, VectorR3 *vm4, VectorR3 &vp, VectorR3 &vm, int &_nvm)
{
	VectorR3 vgp[4], vgm[4];
	VectorR3 rou_p, rou_m; //��+-
	Complex Ve;
	value_t kp = 0.0f, km = 0.0f;
	//kp = (epsp - 1.0f) / epsp;
	//km = (epsm - 1.0f) / epsm;

	Gauss4Point(vp4[0], vp4[1], vp4[2], vp4[3], vgp);
	Gauss4Point(vm4[0], vm4[1], vm4[2], vm4[3], vgm);

	VectorC3 Ei(0, 0, 0);
	Complex Vep(0, 0), Vem(0, 0), G(0, 0);

	for (int a = 0; a < 4; a++)
	{
		//////////////////////////////////////////////////////////////////////////
		rou_p = vgp[a] - vp;
		G = exp(-J0 * k_ * (vgp[a] ^ inc_k_));

		Ei = G * inc_e_;
		Vep += wt4[a] * (rou_p ^ Ei);

		//////////////////////////////////////////////////////////////////////////
		if (_nvm != -1)
		{
			rou_m = vm - vgm[a];
			G = exp(-J0 * k_ * (vgm[a] ^ inc_k_));

			Ei = G * inc_e_;
			Vem += wt4[a] * (rou_m ^ Ei);
		}

	}

	Ve = Vep + Vem;

	return Ve;

}

Complex VSIE::MVKernel(VectorR3 *vp3, VectorR3 *vm3, VectorR3 &vp, VectorR3 &vm)
{
	VectorR3 vgp[7], vgm[7];
	VectorR3 rou_p, rou_m; //��+-
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
	//����Ҫ���� ��0.5*L��
}

void VSIE::KernelIntegral()
{
	Tetrahedron tet_i, tet_j;
	VectorR3 v_tet_i[4], v_tet_j[4];
	VectorR3 vgi[4], vgj[4];

	Qcout << '\n' << "Fill the Kernelintegral:";
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

void VSIE::FillZKernel(VectorR3 *vm4, VectorR3 *vn4, int &m, int &n, VectorR3 *v4)
{
	Complex G(0, 0);
	value_t r;
	VectorR3 v4_3[4][3], *Pv4_3[4];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			r = (vm4[i] - vn4[j]).Norm();

			if (m != n)
			{
				G = exp(-J0*k_*r) / r;

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
				G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);

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

void VSIE::fillZ()
{
	KernelIntegral();
	fillZDD();
	fillZDM();
	fillZMD();
	fillZMM();
	/*int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle 
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3], v_fld_plu3[3], v_fld_min3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3], v_src_plu3[3], v_src_min3[3];
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;
	value_t a_fld, a_src;//SWG����Դ�Ĺ��������
	value_t l_fld, l_src;//RWG����Դ�Ĺ����߳���
	value_t epsp, epsm;//Դ����糣��

	//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm_f, vm_s,vm, v_ct1, v_ct2, v_ct3;
	int tri_plu1, tri_plu2, tri_plu3, tri_min1, tri_min2, tri_min3;
	Complex coef(0.0f,0.0f);

	//���ZDD����
	tool::BarAndPercent bar_perc1;
	coef = Z0*Z0 / (J0* PI4);
	for (int f = 0; f < unknowns_t; f++)
	{
		bar_perc1(f + 1, unknowns_t);

		ct_ptr_->getCommonTriangle(f, vp, vm_f, t_fld_plu, t_fld_min, a_fld);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);
		if (vm_f != -1)
		{
			v_fld_plu = mesh_ptr_->getVertex(vp);//������+
			v_fld_min = mesh_ptr_->getVertex(vm_f);//������-

			tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//��������+
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//��������-

			tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);
		}
		else
		{
			v_fld_plu = mesh_ptr_->getVertex(vp);//����+
			tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//��������+
			tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);
		}

		for (int s = 0; s < unknowns_t; s++)
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
				Zpp = DDZppSingular(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, epsp, epsm);
			}
			else
			{
				Zpp = DDZppKernel(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, epsp, epsm);
			}
			//field tetrahedron+ <-->source tetrahedron-
			if (t_fld_plu == t_src_min)
			{
				Zpm = DDZpmSingular(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, epsp, epsm);
			}
			else if (t_src_min != -1)
			{
				Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, epsp, epsm);
			}
			//field tetrahedron- <-->source tetrahedron+
			if (t_fld_min == t_src_plu)
			{
				Zmp = DDZmpSingular(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, epsp, epsm);
			}
			else if (t_fld_min != -1)
			{
				Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, epsp, epsm);
			}
			//field tetrahedron- <-->source tetrahedron-
			if ((t_fld_min != -1) && (t_src_min != -1))
			{
				if (t_fld_min == t_src_min)
				{
					Zmm = DDZmmSingular(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, epsp, epsm);
				}
				else
					Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, epsp, epsm);
			}

			Z(f, s) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}

	//���ZDM����
	tool::BarAndPercent bar_perc2;
	coef = J0*k_*Z0 / PI4;
	for (int f = 0; f < unknowns_t; f++)
	{
		bar_perc2(f + 1, unknowns_t);
		ct_ptr_->getCommonTriangle(f, vp, vm_f, t_fld_plu, t_fld_min, a_fld);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);
		if (vm_f != -1)
		{
			v_fld_plu = mesh_ptr_->getVertex(vp);//������+
			v_fld_min = mesh_ptr_->getVertex(vm_f);//������-

			tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//��������+
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//��������-

			tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);
		}
		else
		{
			v_fld_plu = mesh_ptr_->getVertex(vp);//����+
			tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//��������+
			tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);
		}
		for (int s = 0; s < unknowns_e; s++)
		{
			ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src);

			v_src_plu = mesh_ptr_->getVertex(vp);
			v_src_min = mesh_ptr_->getVertex(vm);

			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
			tri_min = mesh_ptr_->getTriangleRef(f_src_min);

			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);
			tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
			tri_min.getVertex(tri_min1, tri_min2, tri_min3);
			//fill
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			
			if (vm_f != -1)
			{
				Zmp = DMZmpKernel(v_fld_min4, v_src_plu3, v_fld_min, v_src_plu, v_fld_ct3, vm_f);
				Zmm = DMZmmKernel(v_fld_min4, v_src_min3, v_fld_min, v_src_min, v_fld_ct3, vm_f);
			}
			else if ((tri_plu1 == v_ct1&&tri_plu2 == v_ct2&&tri_plu3 == v_ct3) 
			      || (tri_plu1 == v_ct1&&tri_plu2 == v_ct3&&tri_plu3 == v_ct2)
		          || (tri_plu1 == v_ct2&&tri_plu2 == v_ct1&&tri_plu3 == v_ct3)
		          || (tri_plu1 == v_ct2&&tri_plu2 == v_ct3&&tri_plu3 == v_ct1)
			      || (tri_plu1 == v_ct3&&tri_plu2 == v_ct1&&tri_plu3 == v_ct2)
			      || (tri_plu1 == v_ct3&&tri_plu2 == v_ct2&&tri_plu3 == v_ct1))
			{
				Zpp = DMZppSingular(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
			}
			else if ((tri_min1 == v_ct1&&tri_min2 == v_ct2&&tri_min3 == v_ct3)
				  || (tri_min1 == v_ct1&&tri_min2 == v_ct3&&tri_min3 == v_ct2)
				  || (tri_min1 == v_ct2&&tri_min2 == v_ct1&&tri_min3 == v_ct3)
				  || (tri_min1 == v_ct2&&tri_min2 == v_ct3&&tri_min3 == v_ct1)
				  || (tri_min1 == v_ct3&&tri_min2 == v_ct1&&tri_min3 == v_ct2)
				  || (tri_min1 == v_ct3&&tri_min2 == v_ct2&&tri_min3 == v_ct1))
			{
				Zpm = DMZpmSingular(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else
			{
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			
			Z(f, s + unknowns_t) = coef*a_fld*l_src*(Zpp + Zpm + Zmp + Zmm);
		}
		
	}

	//���ZDM����
	tool::BarAndPercent bar_perc3;
	coef = Z0*Z0*k_ / (J0*PI4);
	for (int f = 0; f < unknowns_e; f++)
	{
		bar_perc3(f + 1, unknowns_e);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld);

		v_fld_plu = mesh_ptr_->getVertex(vp);//����+
		v_fld_min = mesh_ptr_->getVertex(vm);//����-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//��������+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//��������-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);
		tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
		tri_min.getVertex(tri_min1, tri_min2, tri_min3);
		for (int s = 0; s < unknowns_t; s++)
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
			//fill
			if ((tri_plu1 == v_ct1&&tri_plu2 == v_ct2&&tri_plu3 == v_ct3)
			 || (tri_plu1 == v_ct1&&tri_plu2 == v_ct3&&tri_plu3 == v_ct2)
			 || (tri_plu1 == v_ct2&&tri_plu2 == v_ct1&&tri_plu3 == v_ct3)
			 || (tri_plu1 == v_ct2&&tri_plu2 == v_ct3&&tri_plu3 == v_ct1)
			 || (tri_plu1 == v_ct3&&tri_plu2 == v_ct1&&tri_plu3 == v_ct2)
			 || (tri_plu1 == v_ct3&&tri_plu2 == v_ct2&&tri_plu3 == v_ct1))
			{
				Zpp = DMZppSingular(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
			}
			else
				Zpp = DMZppKernel(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);

			
			if ((tri_min1 == v_ct1&&tri_min2 == v_ct2&&tri_min3 == v_ct3)
			 || (tri_min1 == v_ct1&&tri_min2 == v_ct3&&tri_min3 == v_ct2)
			 || (tri_min1 == v_ct2&&tri_min2 == v_ct1&&tri_min3 == v_ct3)
			 || (tri_min1 == v_ct2&&tri_min2 == v_ct3&&tri_min3 == v_ct1)
			 || (tri_min1 == v_ct3&&tri_min2 == v_ct1&&tri_min3 == v_ct2)
			 || (tri_min1 == v_ct3&&tri_min2 == v_ct2&&tri_min3 == v_ct1))
			{
				Zmp = DMZmpSingular(v_fld_min3, v_src_plu4, v_fld_min, v_fld_plu, v_src_ct3, vm_s, epsp, epsm);
			}
			else
				Zmp = DMZmpKernel(v_fld_min3, v_src_plu4, v_fld_min, v_fld_plu, v_src_ct3, vm_s, epsp, epsm);

			if (vm_s != -1)
			{
				Zpm = DMZpmKernel(v_fld_plu3, v_src_min4, v_fld_plu, v_src_min, v_src_ct3, vm_s, epsp, epsm);
				Zmm = DMZmmKernel(v_fld_min3, v_src_min4, v_fld_min, v_fld_min, v_src_ct3, vm_s, epsp, epsm);
			}
			
			Z(f+unknowns_t,s) = coef*l_fld*a_src*(Zpp + Zpm + Zmp + Zmm);
		}
	}

	//���ZMM����
	tool::BarAndPercent bar_perc4;
	coef = J0*k_*Z0 / PI4;
	for (int f = 0; f < unknowns_e; f++)
	{
		bar_perc4(f + 1, unknowns_e);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld);

		v_fld_plu = mesh_ptr_->getVertex(vp);//����+
		v_fld_min = mesh_ptr_->getVertex(vm);//����-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//��������+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//��������-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);
		for (int s = 0; s < unknowns_e; s++)
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
				Zpp = MMZppSingular(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			else
			{
				Zpp = MMZppKernel(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			//field face+ <--> source face-
			if (f_fld_plu == f_src_min)
			{
				Zpm = MMZpmSingular(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
			}
			else
			{
				Zpm = MMZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
			}
			//field face- <--> source face+
			if (f_fld_min == f_src_plu)
			{
				Zmp = MMZmpSingular(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
			}
			else
			{
				Zmp = MMZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
			}
			//field face- <--> source face-
			if (f_fld_min == f_src_min)
			{
				Zmm = MMZmmSingular(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
			}
			else
			{
				Zmm = MMZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
			}
			Z(f+unknowns_t, s+unknowns_t) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}*/

}

void VSIE::fillZDD()
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3];
	Tetrahedron tet_plu, tet_min;
	value_t a_fld, a_src;//����Դ�Ĺ��������
	value_t epsp, epsm;//Դ����糣��

	//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp_f, vp_s, vm_f, vm_s, v_ct1, v_ct2, v_ct3;

	Complex coef(0.0f, 0.0f);
	Qcout << '\n' << "Fill the matrix ZDD:";
	//���ZDD����
	tool::BarAndPercent bar_perc1;
	coef = Z0 / (J0* PI4);
	//coef = omiga * Z0 / PI4;
	for (int f = 0; f < unknowns_t; f++)
	{
		bar_perc1(f + 1, unknowns_t);

		ct_ptr_->getCommonTriangle(f, vp_f, vm_f, t_fld_plu, t_fld_min, a_fld);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);
		if (vm_f != -1)//���ڿ����߽��������
		{
			v_fld_plu = mesh_ptr_->getVertex(vp_f);//������+
			v_fld_min = mesh_ptr_->getVertex(vm_f);//������-

			tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//��������+
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//��������-

			tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);
		}
		else//�ڿ����߽��������
		{
			v_fld_plu = mesh_ptr_->getVertex(vp_f);//����+
			tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//��������+
			tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);
		}

		for (int s = 0; s < unknowns_t; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ct_ptr_->getCommonTriangle(s, vp_s, vm_s, t_src_plu, t_src_min, a_src, epsp, epsm);
			ct_ptr_->getCommonTriangle(s, v_ct1, v_ct2, v_ct3);
			v_src_ct3[0] = mesh_ptr_->getVertex(v_ct1);
			v_src_ct3[1] = mesh_ptr_->getVertex(v_ct2);
			v_src_ct3[2] = mesh_ptr_->getVertex(v_ct3);
			if (vm_s != -1)//���ڿ����߽��������
			{
				v_src_plu = mesh_ptr_->getVertex(vp_s);
				v_src_min = mesh_ptr_->getVertex(vm_s);

				tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
				tet_min = mesh_ptr_->getTetrahedronRef(t_src_min);

				tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);
				tet_min.GetVertices(v_src_min4[0], v_src_min4[1], v_src_min4[2], v_src_min4[3]);
			}
			else//�ڿ����߽��������
			{
				v_src_plu = mesh_ptr_->getVertex(vp_s);
				tet_plu = mesh_ptr_->getTetrahedronRef(t_src_plu);
				tet_plu.GetVertices(v_src_plu4[0], v_src_plu4[1], v_src_plu4[2], v_src_plu4[3]);
			}

			//field tetrahedron+ <-->source tetrahedron+ 
			if (t_fld_plu == t_src_plu)
			{
				Zpp = DDZppSingular(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
			}
			else
			{
				Zpp = DDZppKernel(v_fld_plu4, v_src_plu4, v_fld_plu, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_plu, epsp, epsm);
			}
			//field tetrahedron+ <-->source tetrahedron-
			if (t_fld_plu == t_src_min)
			{
				Zpm = DDZpmSingular(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			}
			else if (t_src_min != -1)
			{
				Zpm = DDZpmKernel(v_fld_plu4, v_src_min4, v_fld_plu, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_plu, t_src_min, epsp, epsm);
			}
			//field tetrahedron- <-->source tetrahedron+
			if (t_fld_min == t_src_plu)
			{
				Zmp = DDZmpSingular(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			}
			else if (t_fld_min != -1)
			{
				Zmp = DDZmpKernel(v_fld_min4, v_src_plu4, v_fld_min, v_src_plu, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_plu, epsp, epsm);
			}
			//field tetrahedron- <-->source tetrahedron-
			if ((t_fld_min != -1) && (t_src_min != -1))
			{
				if (t_fld_min == t_src_min)
				{
					Zmm = DDZmmSingular(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
				}
				else
					Zmm = DDZmmKernel(v_fld_min4, v_src_min4, v_fld_min, v_src_min, v_fld_ct3, v_src_ct3, vm_f, vm_s, t_fld_min, t_src_min, epsp, epsm);
			}

			Z_eigen(f, s) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
			//Z(f, s) = coef * a_fld * a_src * (Zpp + Zpm + Zmp + Zmm);
			
		}
	}
}

void VSIE::fillZDM()
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle 
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;
	value_t a_fld, a_src;//SWG����Դ�Ĺ��������
	value_t l_fld, l_src;//RWG����Դ�Ĺ����߳���
	value_t epsp, epsm;//Դ����糣��

	//Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm_f, vm_s, vm, v_ct1, v_ct2, v_ct3;
	int tri_plu1, tri_plu2, tri_plu3, tri_min1, tri_min2, tri_min3;
	Complex coef(0.0f, 0.0f);
	Qcout << '\n' << "Fill the matrix ZDM:";
	//���ZDM����
	tool::BarAndPercent bar_perc2;
	coef = J0*k_*Z0 / PI4;
	for (int f = 0; f < unknowns_t; f++)
	{
		bar_perc2(f + 1, unknowns_t);
		ct_ptr_->getCommonTriangle(f, vp, vm_f, t_fld_plu, t_fld_min, a_fld);
		ct_ptr_->getCommonTriangle(f, v_ct1, v_ct2, v_ct3);
		v_fld_ct3[0] = mesh_ptr_->getVertex(v_ct1);
		v_fld_ct3[1] = mesh_ptr_->getVertex(v_ct2);
		v_fld_ct3[2] = mesh_ptr_->getVertex(v_ct3);
		if (vm_f != -1)//���ڿ����߽��������
		{
			v_fld_plu = mesh_ptr_->getVertex(vp);//������+
			v_fld_min = mesh_ptr_->getVertex(vm_f);//������-

			tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//��������+
			tet_min = mesh_ptr_->getTetrahedronRef(t_fld_min);//��������-

			tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);
			tet_min.GetVertices(v_fld_min4[0], v_fld_min4[1], v_fld_min4[2], v_fld_min4[3]);
		}
		else//�ڿ����߽��������
		{
			v_fld_plu = mesh_ptr_->getVertex(vp);//����+
			tet_plu = mesh_ptr_->getTetrahedronRef(t_fld_plu);//��������+
			tet_plu.GetVertices(v_fld_plu4[0], v_fld_plu4[1], v_fld_plu4[2], v_fld_plu4[3]);
		}
		for (int s = 0; s < unknowns_e; s++)
		{
			ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src);

			v_src_plu = mesh_ptr_->getVertex(vp);
			v_src_min = mesh_ptr_->getVertex(vm);

			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
			tri_min = mesh_ptr_->getTriangleRef(f_src_min);

			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);
			tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
			tri_min.getVertex(tri_min1, tri_min2, tri_min3);
			//fill
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);

			if (vm_f != -1)//�������������Ӵ�
			{
				Zmp = DMZmpKernel(v_fld_min4, v_src_plu3, v_fld_min, v_src_plu, v_fld_ct3, vm_f);
				Zmm = DMZmmKernel(v_fld_min4, v_src_min3, v_fld_min, v_src_min, v_fld_ct3, vm_f);
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else if ((tri_plu1 == v_ct1&&tri_plu2 == v_ct2&&tri_plu3 == v_ct3)
				  || (tri_plu1 == v_ct1&&tri_plu2 == v_ct3&&tri_plu3 == v_ct2)
				  || (tri_plu1 == v_ct2&&tri_plu2 == v_ct1&&tri_plu3 == v_ct3)
				  || (tri_plu1 == v_ct2&&tri_plu2 == v_ct3&&tri_plu3 == v_ct1)
				  || (tri_plu1 == v_ct3&&tri_plu2 == v_ct1&&tri_plu3 == v_ct2)
				  || (tri_plu1 == v_ct3&&tri_plu2 == v_ct2&&tri_plu3 == v_ct1))//������������������Ԫ�غ�
			{
				Zpp = DMZppSingular(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else if ((tri_min1 == v_ct1&&tri_min2 == v_ct2&&tri_min3 == v_ct3)
				  || (tri_min1 == v_ct1&&tri_min2 == v_ct3&&tri_min3 == v_ct2)
				  || (tri_min1 == v_ct2&&tri_min2 == v_ct1&&tri_min3 == v_ct3)
				  || (tri_min1 == v_ct2&&tri_min2 == v_ct3&&tri_min3 == v_ct1)
				  || (tri_min1 == v_ct3&&tri_min2 == v_ct1&&tri_min3 == v_ct2)
				  || (tri_min1 == v_ct3&&tri_min2 == v_ct2&&tri_min3 == v_ct1))//�������������θ���Ԫ�غ�
			{
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmSingular(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}
			else
			{
				Zpp = DMZppKernel(v_fld_plu4, v_src_plu3, v_fld_plu, v_src_plu, v_fld_ct3, vm_f);
				Zpm = DMZpmKernel(v_fld_plu4, v_src_min3, v_fld_plu, v_src_min, v_fld_ct3, vm_f);
			}

			
			Z_eigen(f, s + unknowns_t) = coef*a_fld*l_src*(Zpp + Zpm + Zmp + Zmm);
			//Z(f, s + unknowns_t) = coef * a_fld*l_src*(Zpp + Zpm + Zmp + Zmm);
		}

	}
}

void VSIE::fillZMD()
{
	int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle 
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3], v_fld_plu3[3], v_fld_min3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3], v_src_plu3[3], v_src_min3[3];
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;
	value_t a_fld, a_src;//SWG����Դ�Ĺ��������
	value_t l_fld, l_src;//RWG����Դ�Ĺ����߳���
	value_t epsp, epsm;//Դ����糣��

					   //Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm_f, vm_s, vm, v_ct1, v_ct2, v_ct3;
	int tri_plu1, tri_plu2, tri_plu3, tri_min1, tri_min2, tri_min3;
	Complex coef(0.0f, 0.0f);
	Qcout << '\n' << "Fill the matrix ZDM:";
	//���ZMD����
	tool::BarAndPercent bar_perc3;
	coef = Z0*k_ / (J0*PI4);
	//coef = k_ * Z0*omiga / PI4;
	for (int f = 0; f < unknowns_e; f++)
	{
		bar_perc3(f + 1, unknowns_e);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld);

		v_fld_plu = mesh_ptr_->getVertex(vp);//����+
		v_fld_min = mesh_ptr_->getVertex(vm);//����-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//��������+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//��������-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);
		tri_plu.getVertex(tri_plu1, tri_plu2, tri_plu3);
		tri_min.getVertex(tri_min1, tri_min2, tri_min3);
		for (int s = 0; s < unknowns_t; s++)
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
			//fill
			if (vm_s != -1)
			{
				Zpp = MDZppKernel(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
				Zpm = MDZpmKernel(v_fld_plu3, v_src_min4, v_fld_plu, v_src_min, v_src_ct3, vm_s, epsp, epsm);
				Zmp = MDZmpKernel(v_fld_min3, v_src_plu4, v_fld_min, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
				Zmm = MDZmmKernel(v_fld_min3, v_src_min4, v_fld_min, v_src_min, v_src_ct3, vm_s, epsp, epsm);
			}
			else if ((tri_plu1 == v_ct1&&tri_plu2 == v_ct2&&tri_plu3 == v_ct3)
			    || (tri_plu1 == v_ct1&&tri_plu2 == v_ct3&&tri_plu3 == v_ct2)
				|| (tri_plu1 == v_ct2&&tri_plu2 == v_ct1&&tri_plu3 == v_ct3)
				|| (tri_plu1 == v_ct2&&tri_plu2 == v_ct3&&tri_plu3 == v_ct1)
				|| (tri_plu1 == v_ct3&&tri_plu2 == v_ct1&&tri_plu3 == v_ct2)
				|| (tri_plu1 == v_ct3&&tri_plu2 == v_ct2&&tri_plu3 == v_ct1))
			{
				Zpp = MDZppSingular(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
				Zmp = MDZmpKernel(v_fld_min3, v_src_plu4, v_fld_min, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
			}
			else if ((tri_min1 == v_ct1&&tri_min2 == v_ct2&&tri_min3 == v_ct3)
				|| (tri_min1 == v_ct1&&tri_min2 == v_ct3&&tri_min3 == v_ct2)
				|| (tri_min1 == v_ct2&&tri_min2 == v_ct1&&tri_min3 == v_ct3)
				|| (tri_min1 == v_ct2&&tri_min2 == v_ct3&&tri_min3 == v_ct1)
				|| (tri_min1 == v_ct3&&tri_min2 == v_ct1&&tri_min3 == v_ct2)
				|| (tri_min1 == v_ct3&&tri_min2 == v_ct2&&tri_min3 == v_ct1))
			{
				Zpp = MDZppKernel(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
				Zmp = MDZmpSingular(v_fld_min3, v_src_plu4, v_fld_min, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
			}
			else
			{
				Zpp = MDZppKernel(v_fld_plu3, v_src_plu4, v_fld_plu, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
				Zmp = MDZmpKernel(v_fld_min3, v_src_plu4, v_fld_min, v_src_plu, v_src_ct3, vm_s, epsp, epsm);
			}

			Z_eigen(f + unknowns_t, s) = coef*l_fld*a_src*(Zpp + Zpm + Zmp + Zmm);
			//Z(f + unknowns_t, s) = coef * l_fld*a_src*(Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void VSIE::fillZMM()
{
	/*int t_fld_plu, t_fld_min, t_src_plu, t_src_min; //tetrahedron +- of source and field triangle
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle 
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu4[4], v_fld_min4[4], v_fld_ct3[3], v_fld_plu3[3], v_fld_min3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu4[4], v_src_min4[4], v_src_ct3[3], v_src_plu3[3], v_src_min3[3];
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;
	value_t a_fld, a_src;//SWG����Դ�Ĺ��������
	value_t l_fld, l_src;//RWG����Դ�Ĺ����߳���
	value_t epsp, epsm;//Դ����糣��

					   //Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm_f, vm_s, vm, v_ct1, v_ct2, v_ct3;
	int tri_plu1, tri_plu2, tri_plu3, tri_min1, tri_min2, tri_min3;
	Complex coef(0.0f, 0.0f);*/

	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu3[3], v_fld_min3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
	Triangle tri_plu, tri_min;
	value_t l_fld, l_src;//����Դ�Ĺ����߳���

	Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm;
	Complex coef = (J0 * k_ * Z0) / PI4;
	Qcout << '\n' << "Fill the matrix ZMM:";
	//���ZMM����
	tool::BarAndPercent bar_perc4;
	//coef = (J0*k_*Z0) / PI4;
	for (int f = 0; f < unknowns_e; f++)
	{
		bar_perc4(f + 1, unknowns_e);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld);

		v_fld_plu = mesh_ptr_->getVertex(vp);//����+
		v_fld_min = mesh_ptr_->getVertex(vm);//����-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//��������+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//��������-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);
		for (int s = 0; s < unknowns_e; s++)
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
				Zpp = MMZppSingular(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			else
			{
				Zpp = MMZppKernel(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			//field face+ <--> source face-
			if (f_fld_plu == f_src_min)
			{
				Zpm = MMZpmSingular(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
			}
			else
			{
				Zpm = MMZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
			}
			//field face- <--> source face+
			if (f_fld_min == f_src_plu)
			{
				Zmp = MMZmpSingular(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
			}
			else
			{
				Zmp = MMZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
			}
			//field face- <--> source face-
			if (f_fld_min == f_src_min)
			{
				Zmm = MMZmmSingular(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
			}
			else
			{
				Zmm = MMZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
			}

			Z_eigen(f + unknowns_t, s + unknowns_t) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
			//Z(f + unknowns_t, s + unknowns_t) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void VSIE::fillV()
{
	int nvp, nvm, fp, fm;
	VectorR3 vp, vm, vp4[4], vm4[4];
	VectorR3 vp3[3], vm3[3];
	value_t area, ln;
	Tetrahedron tet_plu, tet_min;
	Triangle tri_plu, tri_min;

	tool::BarAndPercent bar_perc1;   //
	for (int u = 0; u < unknowns_t; u++)
	{
		bar_perc1(u + 1, unknowns_t);    //

		ct_ptr_->getCommonTriangle(u, nvp, nvm, fp, fm, area);
		vp = mesh_ptr_->getVertex(nvp);
		tet_plu = mesh_ptr_->getTetrahedronRef(fp);
		tet_plu.GetVertices(vp4[0], vp4[1], vp4[2], vp4[3]);
		if (nvm != -1)
		{
			vm = mesh_ptr_->getVertex(nvm);
			tet_min = mesh_ptr_->getTetrahedronRef(fm);
			tet_min.GetVertices(vm4[0], vm4[1], vm4[2], vm4[3]);
		}

		Complex vk = DVKernel(vp4, vm4, vp, vm, nvm);

		V_eigen(u) = area * vk / 3.0f;
		//V(u) = area * vk / 3.0f;
	}

	tool::BarAndPercent bar_perc2;   //
	for (int u = 0; u < unknowns_e; u++)
	{
		bar_perc2(u + 1, unknowns_e);    //

		ce_ptr_->getCommonEdge(u, nvp, nvm, fp, fm, ln);
		vp = mesh_ptr_->getVertex(nvp);
		vm = mesh_ptr_->getVertex(nvm);
		tri_plu = mesh_ptr_->getTriangleRef(fp);
		tri_min = mesh_ptr_->getTriangleRef(fm);
		tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);
		tri_min.getVertex(vm3[0], vm3[1], vm3[2]);

		Complex vk = MVKernel(vp3, vm3, vp, vm);

		V_eigen(u + unknowns_t) = 0.5f * ln * vk;
		//V(u + unknowns_t) = 0.5f * ln * vk;
	}
}

void VSIE::prepareVIE(value_t _eps1, value_t _mu1)
{
	eps1_ = _eps1;
	mu1_ = _mu1;
}

bool VSIE::readExcEdges(const Qstring & rad_file, Qstring& stateInfo)
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

void VSIE::readExcEdges(const Qstring & rad_file)
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

void VSIE::radiateV()
{
	int nv1, nv2;

	tool::BarAndPercent bar_perc;
	for (int u = 0; u < unknowns_e; ++u)
	{
		bar_perc(u + 1, unknowns_e);

		ce_ptr_->getCommonEdge(u, nv1, nv2);
		auto length = ce_ptr_->getCommonEdgeLength(u);

		for (const auto& elem : exc_edge_)
		{
			if ((elem.first == nv1 && elem.second == nv2) || (elem.first == nv2 && elem.second == nv1))
			{
				if (elem.first == nv1)
					V_eigen(u+unknowns_t) = length;
				else
					V_eigen(u+unknowns_t) = -length;
			break;
			}
		}
	}
}

value_t VSIE::getBiRCS(const VectorR3& sca_k) const
{
	int nvp, nvm, fp, fm;
	VectorR3 vp, vm, vp4[4], vm4[4], vp3[3], vm3[3];
	value_t S, L, epsp = 0.0f, epsm = 0.0f;

	VectorR3 vgp_d[5], vgm_d[5];
	VectorR3 vgp_m[7], vgm_m[7];

	VectorR3 rou_p, rou_m;
	Complex G0(0, 0);
	VectorC3 Es_d(0, 0, 0), Es_m(0, 0, 0);

	for (int u = 0; u < unknowns_t; u++)
	{
		ct_ptr_->getCommonTriangle(u, nvp, nvm, fp, fm, S, epsp, epsm);

		vp = mesh_ptr_->getVertex(nvp);
		const auto& tet_plu = mesh_ptr_->getTetrahedronRef(fp);
		tet_plu.GetVertices(vp4[0], vp4[1], vp4[2], vp4[3]);
		Gauss5Point(vp4[0], vp4[1], vp4[2], vp4[3], vgp_d);

		if (nvm != -1)
		{
			vm = mesh_ptr_->getVertex(nvm);
			const auto& tet_min = mesh_ptr_->getTetrahedronRef(fm);
			tet_min.GetVertices(vm4[0], vm4[1], vm4[2], vm4[3]);
			Gauss5Point(vm4[0], vm4[1], vm4[2], vm4[3], vgm_d);
		}

		VectorC3 Esp(0, 0, 0), Esm(0, 0, 0);
		for (int g = 0; g < 5; g++)
		{
			rou_p = vgp_d[g] - vp;
			G0 = exp(J0 * k_ * (vgp_d[g] ^ sca_k));
			Esp += rou_p * wt5[g] * G0*(1.0f - (1.0f / epsp));
			if (nvm != -1)
			{
				rou_m = vm - vgm_d[g];
				G0 = exp(J0 * k_ * (vgm_d[g] ^ sca_k));
				Esm += rou_m * wt5[g] * G0*(1.0f - (1.0f / epsm));
			}

		}
		Es_d += (Esp + Esm) * (I_eigen[u] * S);
	}

	for (int u = 0; u < unknowns_e; u++)
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
		Es_m += (Esp + Esm) * (I_eigen[u+unknowns_t] * L);
	}

	auto es = sca_k * (sca_k*((Es_d / 3.0f) + (Es_m / 2.0f)));
	//auto es = sca_k*(sca_k*((J0*omiga*Es_d / 3.0f) + (Es_m / 2.0f)));
	value_t rcs = k_ * k_ * Z0*Z0* es.norm() / PI4;

	return 10 * log10(rcs);
}

void VSIE::getEFiled(component::FieldData * pdata, const VectorR3 & rad_k, const VectorR3 & rad_ev, const VectorR3 & rad_eh) const
{
	int nvp, nvm, fp, fm;
	VectorR3 vp, vm, vp4[4], vm4[4], vp3[3], vm3[3];
	value_t S, L, epsp = 0.0f, epsm = 0.0f;

	VectorR3 vgp_d[5], vgm_d[5];
	VectorR3 vgp_m[7], vgm_m[7];

	VectorR3 rou_p, rou_m;
	Complex G0(0, 0);
	VectorC3 Es_d(0, 0, 0), Es_m(0, 0, 0);

	for (int u = 0; u < unknowns_t; u++)
	{
		ct_ptr_->getCommonTriangle(u, nvp, nvm, fp, fm, S, epsp, epsm);

		vp = mesh_ptr_->getVertex(nvp);
		const auto& tet_plu = mesh_ptr_->getTetrahedronRef(fp);
		tet_plu.GetVertices(vp4[0], vp4[1], vp4[2], vp4[3]);
		Gauss5Point(vp4[0], vp4[1], vp4[2], vp4[3], vgp_d);

		if (nvm != -1)
		{
			vm = mesh_ptr_->getVertex(nvm);
			const auto& tet_min = mesh_ptr_->getTetrahedronRef(fm);
			tet_min.GetVertices(vm4[0], vm4[1], vm4[2], vm4[3]);
			Gauss5Point(vm4[0], vm4[1], vm4[2], vm4[3], vgm_d);
		}

		VectorC3 Esp(0, 0, 0), Esm(0, 0, 0);
		for (int g = 0; g < 5; g++)
		{
			rou_p = vgp_d[g] - vp;
			G0 = exp(J0 * k_ * (vgp_d[g] ^ rad_k));
			Esp += rou_p * wt5[g] * G0*(1.0f - (1.0f / epsp));
			if (nvm != -1)
			{
				rou_m = vm - vgm_d[g];
				G0 = exp(J0 * k_ * (vgm_d[g] ^ rad_k));
				Esm += rou_m * wt5[g] * G0*(1.0f - (1.0f / epsm));
			}

		}
		Es_d += (Esp + Esm) * (I_eigen[u] * S);
	}

	for (int u = 0; u < unknowns_e; u++)
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
			G0 = exp(J0 * k_ * (vgp_m[g] ^ rad_k));
			Esp = Esp + rou_p * w7[g] * G0;

			rou_m = vm - vgm_m[g];
			G0 = exp(J0 * k_ * (vgm_m[g] ^ rad_k));
			Esm = Esm + rou_m * w7[g] * G0;
		}
		Es_m += (Esp + Esm) * (I_eigen[u + unknowns_t] * L);
	}

	auto es = (Es_d / 3.0f) + (Es_m / 2.0f);
	const auto coeff = -J0 * Z0 * k_ / (PI4);
	pdata->ftheta = coeff * (es ^ rad_ev);
	pdata->fphi = coeff * (es ^ rad_eh);
}

void VSIE::getNearEField(std::vector<component::NearFieldData>* data) const
{
	VectorR3 ori = nearfield_.origin_point;
	VectorR3 end = nearfield_.end_point;
	int samp_x_num = nearfield_.sampling_x;
	int samp_y_num = nearfield_.sampling_y;
	int samp_z_num = nearfield_.sampling_z;

	value_t delta_x, delta_y, delta_z;
	std::vector<VectorR3> pointArr;

	if ((end - ori).x < 1e-6)
		delta_x = 0;
	else if (samp_x_num > 1)
		delta_x = (end - ori).x / (samp_x_num - 1);
	else
		delta_x = 0;

	if ((end - ori).y < 1e-6)
		delta_y = 0;
	else if (samp_y_num > 1)
		delta_y = (end - ori).y / (samp_y_num - 1);
	else
		delta_y = 0;

	if ((end - ori).z < 1e-6)
		delta_z = 0;
	else if (samp_z_num > 1)
		delta_z = (end - ori).z / (samp_z_num - 1);
	else
		delta_z = 0;

	VectorR3 point = nearfield_.origin_point;
	VectorR3 point_temp;
	for (int i = 0; i < samp_z_num; i++)
	{
		point_temp.z = point.z + i * delta_z;
		for (int j = 0; j < samp_y_num; j++)
		{
			point_temp.y = point.y + j * delta_y;
			for (int k = 0; k < samp_x_num; k++)
			{
				point_temp.x = point.x + k * delta_x;
				pointArr.push_back(point_temp);
			}
		}
	}
	int point_num = pointArr.size();
	
	VectorC3 E_near;
	Complex coff = -J0 * k_*Z0 / PI4;

	for (int p = 0; p < point_num; p++)
	{
		VectorC3 E_vol, E_surf,E_inc, e_vol_temp, e_surf_temp;
		VectorR3 ob_point = pointArr[p];
		for (int t = 0; t < unknowns_t; t++)
		{
			NEFkernel_vol(ob_point, t, e_vol_temp);
			E_vol += I(t)*e_vol_temp;
		}
		for (int e = 0; e < unknowns_e; e++)
		{
			NEFkernel_surf(ob_point, e, e_surf_temp);
			E_surf += I(e + unknowns_t)*e_surf_temp;
		}
		E_near = coff * (E_vol + E_surf);
		E_inc = inc_e_*exp(-J0 * k_*(ob_point^inc_k_));
		NearFieldData EFdata(ob_point, E_near);
		EFdata.sca_inc = E_near +E_inc;
		data->push_back(EFdata);
	}
}

void VSIE::NEFkernel_vol(VectorR3& ob, int& unk,VectorC3& e_vol) const
{
	int nvp, nvm, tp, tm, v1, v2, v3;
	value_t area, epsp, epsm;
	VectorR3 vp, vm, vp4[4], vm4[4], vt3[3];
	VectorR3 vgp4[4], vgm4[4],vg3[3];

	ct_ptr_->getCommonTriangle(unk, nvp, nvm, tp, tm, area, epsp, epsm);
	ct_ptr_->getCommonTriangle(unk, v1, v2, v3);

	auto& tet_plu = mesh_ptr_->getTetrahedronRef(tp);
	tet_plu.GetVertices(vp4[0], vp4[1], vp4[2], vp4[3]);
	Gauss4Point(vp4[0], vp4[1], vp4[2], vp4[3], vgp4);
	vp = mesh_ptr_->getVertex(nvp);

	if (tm != -1)
	{
		auto& tet_min = mesh_ptr_->getTetrahedronRef(tm);
		tet_min.GetVertices(vm4[0], vm4[1], vm4[2], vm4[3]);
		Gauss4Point(vm4[0], vm4[1], vm4[2], vm4[3], vgm4);
		vm = mesh_ptr_->getVertex(nvm);
	}

	vt3[0] = mesh_ptr_->getVertex(v1);
	vt3[1] = mesh_ptr_->getVertex(v2);
	vt3[2] = mesh_ptr_->getVertex(v3);
	Gauss3Point(vt3[0], vt3[1], vt3[2], vg3);

	VectorC3 E1, E2, E3;
	VectorR3 rou_p,rou_m,r;
	Complex G0,G_grad;
	value_t kp, km;
	kp = 1.0f - (1.0f / epsp);
	km = 1.0f - (1.0f / epsm);

	for (int i = 0; i < 4; i++)
	{
		rou_p = vgp4[i] - vp;
		r = ob - vgp4[i];

		G0 = (exp(-J0 * k_*r.Norm())) / r.Norm();
		G_grad = G0 * (-J0 * k_ - (1.0f / r.Norm()));

		E1 += wt4[i] * rou_p*kp*G0;
		E2 += wt4[i] * kp*r.Normalize()*G_grad;
	}

	if (abs(epsp - epsm) > 1e-5)
	{
		for (int i = 0; i < 3; i++)
		{
			r = ob - vg3[i];
			G0 = (exp(-J0 * k_*r.Norm())) / r.Norm();
			G_grad = G0 * (-J0 * k_ - (1.0f / r.Norm()));
			E3 += w3[i] * (km - kp) *r.Normalize()*G_grad;
		}
	}
	if (tm != -1)
	{
		for (int i = 0; i < 4; i++)
		{
			rou_m = vm - vgm4[i];
			r = vgm4[i] - ob;

			G0 = (exp(-J0 * k_*r.Norm())) / r.Norm();
			G_grad = G0 * (-J0 * k_ - (1.0f / r.Norm()));

			E1 += wt4[i] * rou_m * km * G0;
			E2 += wt4[i] * km * r.Normalize() * G_grad;
		}
	}

	e_vol = (area*((E1/3.0f) + (E2 + E3) / (k_*k_)));
}

void VSIE::NEFkernel_surf(VectorR3& ob, int& unk,VectorC3& e_surf) const
{
	int nvp, nvm, fp, fm, v1, v2;
	value_t len;
	VectorR3 vp, vm, vp3[3], vm3[3];
	VectorR3 vgp3[3], vgm3[3];

	ce_ptr_->getCommonEdge(unk, nvp, nvm, fp, fm, len);

	auto &tri_plu = mesh_ptr_->getTriangleRef(fp);
	auto &tri_min = mesh_ptr_->getTriangleRef(fm);

	tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);
	tri_min.getVertex(vm3[0], vm3[1], vm3[2]);

	vp = mesh_ptr_->getVertex(nvp);
	vm = mesh_ptr_->getVertex(nvm);

	Gauss3Point(vp3[0], vp3[1], vp3[2], vgp3);
	Gauss3Point(vm3[0], vm3[1], vm3[2], vgm3);

	VectorR3 r, rou_p, rou_m;
	VectorC3 E1, E2;
	Complex G0, G_grad;

	for (int i = 0; i < 3; i++)
	{
		r = ob - vgp3[i];
		rou_p = vgp3[i] - vp;

		G0 = (exp(-J0 * k_*r.Norm())) / r.Norm();
		G_grad = G0 * (-J0 * k_ - (1.0f / r.Norm()));

		E1 += w3[i] * rou_p*G0;
		E2 += w3[i] * r.Normalize()*G_grad;
	}

	for (int i = 0; i < 3; i++)
	{
		r = vgm3[i] - ob;
		rou_m = vm - vgm3[i];

		G0 = (exp(-J0 * k_*r.Norm())) / r.Norm();
		G_grad = G0 * (-J0 * k_ - (1.0f / r.Norm()));

		E1 += w3[i] * rou_m*G0;
		E2 += w3[i] * r.Normalize() *G_grad;
	}

	e_surf = (len*((E1/2.0f) + E2 / (k_*k_)));
}

bool VSIE::writeZIVData()
{
	Qofstream outputZ(dir_ + "/matrix_Z.txt");
	if (outputZ.fail())
		return false;
	outputZ <<Z_eigen;

	Qofstream outputV(dir_ + "/matrix_V.txt");
	if (outputV.fail())
		return false;
	outputV <<V_eigen;

	Qofstream outputI(dir_ + "/martix_I.txt");
	if (outputI.fail())
		return false;
	outputI << I_eigen;

	return true;
}

void VSIE::calculateSurfaceCurrent(std::vector<component::CurrentData>* currents) const
{
	Qcx_vec cur_vec = I;
	const size_t tri_num = mesh_ptr_->getTriangleNum();
	const size_t node_num = mesh_ptr_->getNodeNum();
	std::vector<std::vector<value_t>> node(node_num);
	std::vector<std::vector<int>> triangle(tri_num);
	for (int u = 0; u < unknowns_e; ++u)
	{
		auto& rwg = ce_ptr_->getRWGRef(u);
		//node[rwg.v1].push_back(u);
		//node[rwg.v2].push_back(u);
		triangle[rwg.tpid].push_back(u);
		triangle[rwg.tmid].push_back(u);
	}
	currents->resize(tri_num);
	VectorR3 vx[3];
	int nx[3];

	value_t len[3];
	for (size_t t = 0; t < tri_num; ++t)
	{
		VectorC3 mag[3];
		auto& unks = triangle[t];
		for (size_t i = 0; i < unks.size(); ++i)
		{
			auto& rwg = ce_ptr_->getRWGRef(unks[i]);
			vx[i] = (t == rwg.tpid) ? mesh_ptr_->getVertex(rwg.vxp) : mesh_ptr_->getVertex(rwg.vxm);
			nx[i] = (t == rwg.tpid) ? rwg.vxp : rwg.vxm;
			len[i] = rwg.length;
		}
		//������ԵRWG
		if (unks.size() == 1)
		{
			auto& rwg = ce_ptr_->getRWGRef(unks[0]);
			vx[1] = mesh_ptr_->getVertex(rwg.v1);
			vx[2] = mesh_ptr_->getVertex(rwg.v2);
			nx[1] = rwg.v1;
			nx[2] = rwg.v2;

		}
		if (unks.size() == 2)
		{
			auto& rwg1 = ce_ptr_->getRWGRef(unks[0]);
			auto& rwg2 = ce_ptr_->getRWGRef(unks[1]);
			if ((rwg1.v1 != rwg2.v1) && (rwg1.v1 != rwg2.v2))
			{
				vx[2] = mesh_ptr_->getVertex(rwg1.v2);
				nx[2] = rwg1.v2;
			}
			else
			{
				vx[2] = mesh_ptr_->getVertex(rwg1.v1);
				nx[2] = rwg1.v1;
			}
		}
		auto center = (vx[0] + vx[1] + vx[2]) / 3;
		auto dominator = 2 * Area(vx[0], vx[1], vx[2]);
		VectorC3 cen_cur(0, 0, 0);
		for (size_t i = 0; i < unks.size(); ++i)
		{
			auto& rwg = ce_ptr_->getRWGRef(unks[i]);
			cen_cur += (t == rwg.tpid ? 1.0f : -1.0f) * cur_vec(unks[i]) * (len[i] / dominator) * (center - vx[i]);
			//VectorC3 ncur(0, 0, 0);
			for (size_t j = 0; j < 3; ++j)
			{
				if (i == j) continue;
				//auto& srwg = ce_ptr_->getRWGRef(unks[j]);
				mag[j] += (t == rwg.tpid ? 1.0f : -1.0f) * cur_vec(unks[i]) * (len[i] / dominator) * (vx[j] - vx[i]);
			}
			//mag[i] = ncur;
			//node[nx[i]].push_back(std::sqrt(ncur.norm()));
		}

		/*value_t Ap(0.0f), Am(0.0f);
		for (int i = 0; i < 3; ++i)
		{
		//auto& _nx = nx[i];
		auto& _node = node[nx[i]];
		VectorC3 ncur(0, 0, 0);
		for (int j = 0; j < _node.size(); ++j)
		{
		auto& rwg = ce_ptr_->getRWGRef(_node[j]);
		Ap = Area(mesh_ptr_->getVertex(rwg.v1), mesh_ptr_->getVertex(rwg.v2), mesh_ptr_->getVertex(rwg.vxp));
		Am = Area(mesh_ptr_->getVertex(rwg.v1), mesh_ptr_->getVertex(rwg.v2), mesh_ptr_->getVertex(rwg.vxm));
		ncur += cur_vec(_node[j])*rwg.length*(((vx[i] - mesh_ptr_->getVertex(rwg.vxp)) / Ap) +((mesh_ptr_->getVertex(rwg.vxm) - vx[i]) / Am));
		}
		if (_node.size() != 0)
		{
		mag[i] = 0.5f*ncur/((value_t)_node.size());
		}
		else
		mag[i] = ncur;

		}*/
		auto& data = (*currents)[t];
		data.v1 = vx[0];
		data.v2 = vx[1];
		data.v3 = vx[2];
		data.n[0] = nx[0];
		data.n[1] = nx[1];
		data.n[2] = nx[2];
		data.magnc = std::sqrt(cen_cur.norm());
		/*data.magn1 = std::sqrt(mag[0].norm());
		data.magn2 = std::sqrt(mag[1].norm());
		data.magn3 = std::sqrt(mag[2].norm());*/
		node[nx[0]].push_back(std::sqrt(mag[0].norm()));
		node[nx[1]].push_back(std::sqrt(mag[1].norm()));
		node[nx[2]].push_back(std::sqrt(mag[2].norm()));
	}
	for (int t = 0; t < tri_num; ++t)
	{
		//auto &tri = mesh_ptr_->getTriangleRef(t);
		auto &data = (*currents)[t];
		value_t magni[3] = { 0.0f,0.0f,0.0f };
		for (int i = 0; i < 3; ++i)
		{
			auto &_node = node[data.n[i]];
			int k = _node.size();
			if (k != 0)
			{
				for (int j = 0; j < k; ++j)
				{
					magni[i] += _node[j];
				}
				magni[i] = magni[i] / ((value_t)k);
			}
		}
		data.magn1 = magni[0];
		data.magn2 = magni[1];
		data.magn3 = magni[2];
	}
	Qcout << "test" << std::endl;
}