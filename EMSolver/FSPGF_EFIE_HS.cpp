#include "FSPGF_EFIE_HS.h"
#include "ConfigLoader.h"
#include "Mathdef.h"
#include "Quadrature.h"
#include "VectorC3.h"
#include "tools.h"
#include "CommonEdge.h"
#include "iml.h"

using namespace component;
using namespace mom;
using namespace math;
using std::setw;


FSPGF_EFIE_HS::FSPGF_EFIE_HS()
{

}

FSPGF_EFIE_HS::~FSPGF_EFIE_HS()
{
}

void FSPGF_EFIE_HS::init(component::ConfigLoader * ploader)
{
	EM::init(ploader);
	if (rtype_ == policy::ResultType::Rad)
		readExcEdges(dir_ + '/' + folder_name_ + ".rad");

	dir_ = tool::creatFolder(dir_, folder_name_ + "_FSPGF_EFIE_HS");

	auto log_name = dir_ + (rtype_ == policy::ResultType::Rad ? "/rad_" : "/sca_") + "runtime.log";
	logger_.open(log_name);

	mesh_ptr_ = std::make_shared<Mesh>();
	ce_ptr_ = std::make_shared<CommonEdge>();

	mesh_ptr_->loadMeshFile(ploader->getMeshPath());
	isContinuous = ploader->getisContinuous();
	if (isContinuous == 1)
	{
		ce_ptr_->buildContinuousEdges(mesh_ptr_);
	}
	else
	{
		ce_ptr_->buildCommonEdge(mesh_ptr_);
	}

	unknowns_ = ce_ptr_->getCommonEdgeNum();
	k_ = PI2 * incidence_.freq / cc;


	prepareFSPGF(ploader->getDx(), ploader->getDy(), ploader->getArray_x(), ploader->getArray_y(), ploader->gett_sum());
	Sigma = ploader->getSigma();

	logger_ << ploader->getReport();
	reportInfo(logger_);
	TIME_LOG("init");
}

void FSPGF_EFIE_HS::prepareFSPGF(value_t _Dx, value_t _Dy, int _Array_x, int _Array_y, int _t_sum)
{
	auto lamda = cc / incidence_.freq;
	Dx = _Dx*lamda;
	Dy = _Dy*lamda;
	Array_x = _Array_x;
	Array_y = _Array_y;
	t_sum = _t_sum;
}

void FSPGF_EFIE_HS::solve()
{
	SEGMENT("Solve");
	Z.zeros(unknowns_, unknowns_);
	LOG(fillZ(), "Filling impedance matrix");
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

	Qcout << setw(30) << "Solving matrix equation:";
	if (!arma::solve(I, Z, V, arma::solve_opts::fast))
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

void FSPGF_EFIE_HS::output()
{
	SEGMENT("Result");
	Result result;
	if (rtype_ == policy::ResultType::Rad)
	{
		LOG(result.getRadiationField(this, &FSPGF_EFIE_HS::getEFiled), "Calculating Radiation");
		TIME_LOG("getRadiation");
	}
	else
	{
		LOG(result.getBistaticRCS(this, &FSPGF_EFIE_HS::getBiRCS), "Calculating BistaticRCS");
		TIME_LOG("getBistaticRCS");
	}
	LOG(result.getCurrentDistribution(this, &FSPGF_EFIE_HS::calculateSurfaceCurrent), "Calculating surface current");
	TIME_LOG("getSurfaceCurrent");
}

void FSPGF_EFIE_HS::clear()
{
	Z.reset();
	I.reset();
	V.reset();
	exc_edge_.clear();

	mesh_ptr_->clear();
	ce_ptr_->clear();
}


void FSPGF_EFIE_HS::reportInfo(Qostream& strm) const
{
	mesh_ptr_->reportInfo(strm);
	ce_ptr_->reportInfo(strm);
}

Complex FSPGF_EFIE_HS::eZppKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	VectorR3 r,r1;
	//value_t r;
	Complex Ze(0, 0), Zm(0, 0), G(0.0, 0.0),G1(0.0,0.0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs;
			/*r = (vgf[p] - vgs[q]).Norm();
			G = exp(-J0 * k_ * r) / r;*/
			r = vgf[p] - vgs[q];
			r1 = r;
			r1.z = r1.z + 2.0f*vgs[q].z;
			G = FSPGF(r, t_sum, Dx, Dy);
			G1 = FSPGF(r1, t_sum, Dx, Dy);
			Ze += w3[p] * w3[q] * ((0.25f * (rou_f ^ rou_s) - (1.0f / (k_ * k_))) * (G - G1) + 0.5f*rou_f.z*rou_s.z*G1);
		}
	}
	return Ze;
}

Complex FSPGF_EFIE_HS::eZpmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	VectorR3 r, r1;
	//value_t r;
	Complex Ze(0, 0), Zm(0, 0), G(0.0, 0.0), G1(0.0, 0.0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vs - vgs[q];
			/*r = (vgf[p] - vgs[q]).Norm();
			G = exp(-J0 * k_ * r) / r;*/
			r = vgf[p] - vgs[q];
			r1 = r;
			r1.z = r1.z + 2.0f*vgs[q].z;
			G = FSPGF(r, t_sum, Dx, Dy);
			G1 = FSPGF(r1, t_sum, Dx, Dy);
			Ze += w3[p] * w3[q] * ((0.25f * (rou_f ^ rou_s) + (1.0f / (k_ * k_)))* (G - G1) + 0.5f*rou_f.z*rou_s.z*G1);
		}
	}
	return Ze;
}

Complex FSPGF_EFIE_HS::eZmpKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	VectorR3 r, r1;
	//value_t r;
	Complex Ze(0, 0), Zm(0, 0), G(0.0, 0.0), G1(0.0, 0.0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vf - vgf[p];
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs;
			/*r = (vgf[p] - vgs[q]).Norm();
			G = exp(-J0 * k_ * r) / r;*/
			r = vgf[p] - vgs[q];
			r1 = r;
			r1.z = r1.z + 2.0f*vgs[q].z;
			G = FSPGF(r, t_sum, Dx, Dy);
			G1 = FSPGF(r1, t_sum, Dx, Dy);
			Ze += w3[p] * w3[q] * ((0.25f * (rou_f ^ rou_s) + (1.0f / (k_ * k_)))* (G - G1) + 0.5f*rou_f.z*rou_s.z*G1);
		}
	}
	return Ze;
}

Complex FSPGF_EFIE_HS::eZmmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	VectorR3 r, r1;
	//value_t r;
	Complex Ze(0, 0), Zm(0, 0), G(0.0, 0.0), G1(0.0, 0.0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vf - vgf[p];
		for (int q = 0; q < 3; q++)
		{
			rou_s = vs - vgs[q];
			/*r = (vgf[p] - vgs[q]).Norm();
			G = exp(-J0 * k_ * r) / r;*/
			r = vgf[p] - vgs[q];
			r1 = r;
			r1.z = r1.z + 2.0f*vgs[q].z;
			G = FSPGF(r, t_sum, Dx, Dy);
			G1 = FSPGF(r1, t_sum, Dx, Dy);
			Ze += w3[p] * w3[q] * ((0.25f * (rou_f ^ rou_s) - (1.0f / (k_ * k_)))* (G - G1) + 0.5f*rou_f.z*rou_s.z*G1);
		}
	}
	return Ze;
}

Complex FSPGF_EFIE_HS::eZppSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	Complex Ze(0, 0), Zen(0, 0), Zes(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0);
	Complex G1(0.0, 0.0);
	VectorR3 r1;

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

	Zes = (Ze1 + 1.0f / S_s * Ze2);
	
	/*VectorR3 r;
	Complex G(0, 0);
	for (int p = 0; p < 3; p++)
	{
	rou_f = vgf[p] - vf;
	for (int q = 0; q < 3; q++)
	{
	rou_s = vgs[q] - vs;
	r = vgf[p] - vgs[q];
	G = FSPGFSingular(r, t_sum, Dx, Dy);
	Zen += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * G;
	}
	}*/
	//镜像部分
	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs;
			/*r = (vgf[p] - vgs[q]).Norm();
			G = exp(-J0 * k_ * r) / r;*/
			r1 = vgf[p] - vgs[q];
			r1.z = r1.z + 2.0f*vgs[q].z;
			G1 = FSPGF(r1, t_sum, Dx, Dy);
			Zen += w3[p] * w3[q] * ((0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * (-1.0f* G1) + 0.5f*rou_f.z*rou_s.z*G1);
		}
	}

	if (abs(Sigma) > 1e-6)
	{
		for (int i = 0; i < 3; i++)
		{
			rou_f = vgf[i] - vf;
			rou_s = vgf[i] - vs;
			Ze3 += w3[i] * rou_f^rou_s;
		}
		Ze3 = (PI / (Sigma*S_s*J0*k_*Z0))*Ze3;
	}

	Ze = Zen + Zes + Ze3;

	return Ze;
}

Complex FSPGF_EFIE_HS::eZpmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	Complex Ze(0, 0), Zen(0, 0), Zes(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;
	Complex G1(0.0, 0.0);
	VectorR3 r1;

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

	Zes = (Ze1 + 1.0f / S_s * Ze2);

	//其他周期单元部分
	/*VectorR3 r;
	Complex G(0, 0);
	for (int p = 0; p < 3; p++)
	{
	rou_f = vgf[p] - vf;
	for (int q = 0; q < 3; q++)
	{
	rou_s = vgs[q] - vs;
	r = vgf[p] - vgs[q];
	G = FSPGFSingular(r, t_sum, Dx, Dy);
	Zen += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * G;
	}
	}*/
	//镜像部分
	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vs - vgs[q];
			/*r = (vgf[p] - vgs[q]).Norm();
			G = exp(-J0 * k_ * r) / r;*/
			r1 = vgf[p] - vgs[q];
			r1.z = r1.z + 2.0f*vgs[q].z;
			G1 = FSPGF(r1, t_sum, Dx, Dy);
			Zen += w3[p] * w3[q] * ((0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * (-1.0f* G1) + 0.5f*rou_f.z*rou_s.z*G1);
		}
	}

	if (abs(Sigma) > 1e-6)
	{
		for (int i = 0; i < 3; i++)
		{
			rou_f = vgf[i] - vf;
			rou_s = vgf[i] - vs;
			Ze3 += w3[i] * rou_f^rou_s;
		}
		Ze3 = (PI / (Sigma*S_s*J0*k_*Z0))*Ze3;
	}

	Ze = Zen + Zes + Ze3;

	return Ze;
}

Complex FSPGF_EFIE_HS::eZmpSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	Complex Ze(0, 0), Zen(0, 0), Zes(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	Complex G1(0.0, 0.0);
	VectorR3 r1;

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

	Zes = (Ze1 + 1.0f / S_s * Ze2);
	//其他周期单元部分
	/*VectorR3 r;
	Complex G(0, 0);
	for (int p = 0; p < 3; p++)
	{
	rou_f = vgf[p] - vf;
	for (int q = 0; q < 3; q++)
	{
	rou_s = vgs[q] - vs;
	r = vgf[p] - vgs[q];
	G = FSPGFSingular(r, t_sum, Dx, Dy);
	Zen += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * G;
	}
	}*/
	//镜像部分
	for (int p = 0; p < 3; p++)
	{
		rou_f = vf - vgf[p];
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs;
			/*r = (vgf[p] - vgs[q]).Norm();
			G = exp(-J0 * k_ * r) / r;*/
			r1 = vgf[p] - vgs[q];
			r1.z = r1.z + 2.0f*vgs[q].z;
			G1 = FSPGF(r1, t_sum, Dx, Dy);
			Zen += w3[p] * w3[q] * ((0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * (-1.0f* G1) + 0.5f*rou_f.z*rou_s.z*G1);
		}
	}

	if (abs(Sigma) > 1e-6)
	{
		for (int i = 0; i < 3; i++)
		{
			rou_f = vgf[i] - vf;
			rou_s = vgf[i] - vs;
			Ze3 += w3[i] * rou_f^rou_s;
		}
		Ze3 = (PI / (Sigma*S_s*J0*k_*Z0))*Ze3;
	}

	Ze = Zen + Zes + Ze3;

	return Ze;
}

Complex FSPGF_EFIE_HS::eZmmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	Complex Ze(0, 0), Zen(0, 0), Zes(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	Complex G1(0.0, 0.0);
	VectorR3 r1;

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

	Zes = (Ze1 + 1.0f / S_s * Ze2);
	//其他周期单元部分
	/*VectorR3 r;
	Complex G(0, 0);
	for (int p = 0; p < 3; p++)
	{
	rou_f = vgf[p] - vf;
	for (int q = 0; q < 3; q++)
	{
	rou_s = vgs[q] - vs;
	r = vgf[p] - vgs[q];
	G = FSPGFSingular(r, t_sum, Dx, Dy);
	Zen += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * G;
	}
	}*/
	//镜像部分
	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vs - vgs[q];
			/*r = (vgf[p] - vgs[q]).Norm();
			G = exp(-J0 * k_ * r) / r;*/
			r1 = vgf[p] - vgs[q];
			r1.z = r1.z + 2.0f*vgs[q].z;
			G1 = FSPGF(r1, t_sum, Dx, Dy);
			Zen += w3[p] * w3[q] * ((0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * (-1.0f* G1) + 0.5f*rou_f.z*rou_s.z*G1);
		}
	}

	if (abs(Sigma) > 1e-6)
	{
		for (int i = 0; i < 3; i++)
		{
			rou_f = vgf[i] - vf;
			rou_s = vgf[i] - vs;
			Ze3 += w3[i] * rou_f^rou_s;
		}
		Ze3 = (PI / (Sigma*S_s*J0*k_*Z0))*Ze3;
	}

	Ze = Zen + Zes + Ze3;

	return Ze;
}

Complex FSPGF_EFIE_HS::eVKernel(VectorR3 *vp3, VectorR3 *vm3, VectorR3 &vp, VectorR3 &vm)
{
	VectorR3 vgp[7], vgm[7];
	VectorR3 rou_p, rou_m; //ρ+-
	Complex Ve;

	Gauss7Point(vp3[0], vp3[1], vp3[2], vgp);
	Gauss7Point(vm3[0], vm3[1], vm3[2], vgm);

	VectorC3 Ei(0, 0, 0), Er(0, 0, 0)/*, Hi(0, 0, 0)*/;
	Complex Vep(0, 0), Vem(0, 0)/*, Vmp(0, 0), Vmm(0, 0)*/, G(0, 0),Gr(0,0);

	for (int a = 0; a < 7; a++)
	{
		//////////////////////////////////////////////////////////////////////////
		rou_p = vgp[a] - vp;
		G = exp(-J0 * k_ * (vgp[a] ^ inc_k_));
		Gr = exp(-J0 * k_ * (vgp[a] ^ r_inc_k_));

		Ei = G * inc_e_;
		Er = Gr * r_inc_e_;
		Vep += w7[a] * (rou_p ^ (Ei + Er));

		//////////////////////////////////////////////////////////////////////////
		rou_m = vm - vgm[a];
		G = exp(-J0 * k_ * (vgm[a] ^ inc_k_));
		Gr = exp(-J0 * k_ * (vgm[a] ^ r_inc_k_));

		Ei = G * inc_e_;
		Er = Gr * r_inc_e_;
		Vem += w7[a] * (rou_m ^ (Ei + Er));
	}

	Ve = (Vep + Vem);

	return Ve;
	//后面要乘以 （0.5*L）
}

Complex FSPGF_EFIE_HS::FSPGF(VectorR3 &r, int _t_sum, value_t _Dx, value_t _Dy)
{
	Complex PG(0.0, 0.0), PG1(0.0, 0.0), PG2(0.0, 0.0);
	VectorR3 R = r;
	value_t H;
	value_t delta = 0.0f;
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
			if (delta >= ((k_ * k_) / 4))
				a = sqrt(delta - (k_ * k_) / 4);
			else
				a = J0*sqrt(((k_ * k_) / 4 - delta));

			PG1 += ((exp(-2.0f*J0 * PI*(i*r.x / Dx + j*r.y / Dy))) / a)*((exp(2.0f * a*r.z))*(erfc(r.z*H + a / H)) + (exp(-2.0f * a*r.z))*(erfc((a / H) - r.z*H)));
			PG2 += (exp(-J0*k_*(inc_k_.x*i*Dx + inc_k_.y*j*Dy)) / R_N)*(erfc(R_N*H + (-J0*k_) / (2.0f * H))*exp(-J0*k_*R_N)).real();
		}
	}
	PG = (PI*PG1) / (2.0f * Dx*Dy) + PG2;
	return PG;
}

Complex FSPGF_EFIE_HS::FSPGFSingular(VectorR3 &r, int _t_sum, value_t _Dx, value_t _Dy)
{
	Complex PG(0.0, 0.0), PG1(0.0, 0.0), PG2(0.0, 0.0);
	VectorR3 R = r;
	value_t H;
	H = sqrt(PI / (Dx*Dy));
	//PG = FSPGF(r, _t_sum, _Dx, _Dy) - exp(-J0*k_*r.Norm()) / (PI4*r.Norm());

	for (int i = -_t_sum; i <= _t_sum; ++i)
	{
		for (int j = -_t_sum; j <= _t_sum; ++j)
		{
			if ((i == 0) && (j == 0))
				continue;
			value_t delta;
			Complex a;
			delta = ((float)i*PI / Dx)*((float)i*PI / Dx) + ((float)j*PI / Dy) *((float)j*PI / Dy);
			//R = ((d_x - m*Dx) ^ 2 + (d_y - n*Dy) ^ 2 + (d_z) ^ 2);
			R.x = r.x - i*Dx;
			R.y = r.y - j*Dy;
			if (delta >= ((k_ * k_) / 4))
				a = sqrt(delta - (k_ * k_) / 4);
			else
				a = J0*sqrt(((k_ * k_) / 4 - delta));

			PG1 += ((exp(-2.0f*J0 * PI*(i*r.x / Dx + j*r.y / Dy))) / a)*((exp(2.0f * a*r.z))*(erfc(r.z*H + a / H)) + (exp(-2.0f * a*r.z))*(erfc(-r.z*H + a / H)));
			PG2 += (1.0 / R.Norm())*(erfc(abs(R.Norm()*H + (-J0*k_) / (2 * H)))*exp(-J0*k_*R.Norm())).real();
		}
	}
	PG = PG1 / (8 * Dx*Dy) + PG2 / (PI4);
	return PG;
}

VectorR3 FSPGF_EFIE_HS::Translateoff(VectorR3 r, int _x, int _y, int _z) const
{
	VectorR3 r_t;
	r_t.x = r.x + Dx*_x;
	r_t.y = r.y + Dy*_y;
	r_t.z = r.z;

	return r_t;
}

Complex FSPGF_EFIE_HS::erfz(Complex z)
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

		double xk = 0.0, yk = 0.0;

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

		if (abs(z.real()) < 1e-5)
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

void FSPGF_EFIE_HS::fillZ()
{
	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu3[3], v_fld_min3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
	Triangle tri_plu, tri_min;
	value_t l_fld, l_src;//场与源的公共边长度

	Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm;
	Complex coef = (J0 * k_ * Z0) / PI4;

	tool::BarAndPercent bar_perc;
	for (int f = 0; f < unknowns_; f++)
	{
		bar_perc(f + 1, unknowns_);

		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		v_fld_min = mesh_ptr_->getVertex(vm);//场点-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);

		for (int s = 0; s < unknowns_; s++)
		{
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
				Zpp = eZppSingular(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			else
			{
				Zpp = eZppKernel(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
			}
			//field face+ <--> source face-
			if (f_fld_plu == f_src_min)
			{
				Zpm = eZpmSingular(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
			}
			else
			{
				Zpm = eZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
			}
			//field face- <--> source face+
			if (f_fld_min == f_src_plu)
			{
				Zmp = eZmpSingular(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
			}
			else
			{
				Zmp = eZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
			}
			//field face- <--> source face-
			if (f_fld_min == f_src_min)
			{
				Zmm = eZmmSingular(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
			}
			else
			{
				Zmm = eZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
			}
			Z(f, s) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void FSPGF_EFIE_HS::fillV()
{
	int nvp, nvm, fp, fm, bedge;
	VectorR3 vp, vm, vp3[3], vm3[3];
	value_t ln;
	Triangle tri_plu, tri_min;

	const auto &trinum = mesh_ptr_->getTriangleNum();
	VectorR3 trans_vector = mesh_ptr_->getSize();

	tool::BarAndPercent bar_perc;   //
	for (int u = 0; u < unknowns_; u++)
	{
		bar_perc(u + 1, unknowns_);    //

		ce_ptr_->getCommonEdge(u, nvp, nvm, fp, fm, ln, bedge);
		vp = mesh_ptr_->getVertex(nvp);
		vm = mesh_ptr_->getVertex(nvm);
		tri_plu = mesh_ptr_->getTriangleRef(fp);
		//tri_min = mesh_ptr_->getTriangleRef(fm);
		tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);
		//tri_min.getVertex(vm3[0], vm3[1], vm3[2]);

		if (bedge == 1)
		{
			tri_min = mesh_ptr_->getTriangleRef(fm);
			tri_min.getVertex(vm3[0], vm3[1], vm3[2]);
			vm.y += trans_vector.y;
			vm3[0].y += trans_vector.y;
			vm3[1].y += trans_vector.y;
			vm3[2].y += trans_vector.y;
		}
		else if (bedge == 2)
		{
			tri_min = mesh_ptr_->getTriangleRef(fm);
			tri_min.getVertex(vm3[0], vm3[1], vm3[2]);
			vm.x += trans_vector.x;
			vm3[0].x += trans_vector.x;
			vm3[1].x += trans_vector.x;
			vm3[2].x += trans_vector.x;
		}
		else if (bedge == 0)
		{
			tri_min = mesh_ptr_->getTriangleRef(fm);
			tri_min.getVertex(vm3[0], vm3[1], vm3[2]);
		}

		Complex vk = eVKernel(vp3, vm3, vp, vm);
		V(u) = 0.5f * ln * vk;
	}
}

bool FSPGF_EFIE_HS::readExcEdges(const Qstring & rad_file, Qstring& stateInfo)
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

void FSPGF_EFIE_HS::readExcEdges(const Qstring & rad_file)
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

void FSPGF_EFIE_HS::radiateV()
{
	int nv1, nv2;

	tool::BarAndPercent bar_perc;
	for (int u = 0; u < unknowns_; ++u)
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
	}
}

value_t FSPGF_EFIE_HS::getBiRCS(const VectorR3 & sca_k) const
{
	int nvp, nvm, fp, fm;
	VectorR3 vp, vm, vp3[3], vm3[3], r_tp, r_tm, r_tp1, r_tm1;
	value_t L;
	value_t phaseoff_x, phaseoff_y;
	value_t arrfactor_x, arrfactor_y;

	VectorR3 vgp[7], vgm[7];

	VectorR3 rou_p, rou_m, rou_p1, rou_m1;
	Complex G0(0, 0), G01(0, 0);
	VectorC3 Es(0, 0, 0);

	for (int u = 0; u < unknowns_; u++)
	{
		ce_ptr_->getCommonEdge(u, nvp, nvm, fp, fm, L);
		vp = mesh_ptr_->getVertex(nvp);
		vm = mesh_ptr_->getVertex(nvm);
		const auto& tri_plu = mesh_ptr_->getTriangleRef(fp);
		const auto& tri_min = mesh_ptr_->getTriangleRef(fm);
		tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);
		tri_min.getVertex(vm3[0], vm3[1], vm3[2]);

		Gauss7Point(vp3[0], vp3[1], vp3[2], vgp);
		Gauss7Point(vm3[0], vm3[1], vm3[2], vgm);

		VectorC3 Esp(0, 0, 0), Esm(0, 0, 0), Esp1(0, 0, 0), Esm1(0, 0, 0);
		phaseoff_x = sca_k.x*Dx*k_;
		phaseoff_y = sca_k.y*Dy*k_;
		if (abs((sin(phaseoff_x / 2.0f))) < 1e-5)
		{
			arrfactor_x = (value_t)Array_x;
		}
		else
		{
			arrfactor_x = (sin((value_t)Array_x*phaseoff_x / 2.0f)) / (sin(phaseoff_x / 2.0f));
		}

		if (abs(sin(phaseoff_y / 2.0f)) < 1e-5)
		{
			arrfactor_y = (value_t)Array_y;
		}
		else
		{
			arrfactor_y = (sin((value_t)Array_y*phaseoff_y / 2.0f)) / (sin(phaseoff_y / 2.0f));
		}
		for (int x = 0; x < 1; x++)
		{
			for (int y = 0; y < 1; y++)
			{
				for (int g = 0; g < 7; g++)
				{
					//+face
					rou_p = vgp[g] - vp;
					r_tp = Translateoff(vgp[g], x, y, 0);
					G0 = exp(J0 * k_ * (r_tp ^ sca_k));
					Esp += rou_p * w7[g] * G0;
					//+face mirror
					r_tp1 = r_tp;
					r_tp1.z =  - 1.0f*r_tp.z;
					G01 = exp(J0 * k_ * (r_tp1 ^ sca_k));
					rou_p1 = -1.0f*rou_p;
					rou_p1.z = rou_p.z;
					Esp1 += rou_p1 * w7[g] * G01;

					//-face
					rou_m = vm - vgm[g];
					r_tm = Translateoff(vgm[g], x, y, 0);
					G0 = exp(J0 * k_ * (r_tm ^ sca_k));
					Esm += rou_m * w7[g] * G0;
					//-face mirror
					r_tm1 = r_tm;
					r_tm1.z =-1.0f*r_tm.z;
					G01 = exp(J0 * k_ * (r_tm1 ^ sca_k));
					rou_m1 = -1.0f*rou_m;
					rou_m1.z = rou_m.z;
					Esm1 += rou_m1 * w7[g] * G01;
					//rou_p = vgp[g] - vp;
					////r_tp = Translateoff(vgp[g], x, y, 0);
					//G0 = exp(J0 * k_ * (vgp[g] ^ sca_k))*4.0f*cos((value_t)x*Dx*sca_k.x)*cos((value_t)y*Dy*sca_k.y);
					//Esp = Esp + rou_p * w7[g] * G0;

					//rou_m = vm - vgm[g];
					////r_tm = Translateoff(vgm[g], x, y, 0);
					//G0 = exp(J0 * k_ * (vgm[g] ^ sca_k))*4.0f*cos((value_t)x*Dx*sca_k.x)*cos((value_t)y*Dy*sca_k.y);
					//Esm = Esm + rou_m * w7[g] * G0;
				}
			}
		}
		Es += (Esp + Esm + Esp1 + Esm1) * (I[u] * L);
	}

	auto es = 0.5f * sca_k * (sca_k * (Es*arrfactor_x*arrfactor_y));
	value_t rcs = Z0 * Z0 * k_ * k_ * es.norm() / PI4;

	return 10 * log10(rcs);
}

void FSPGF_EFIE_HS::getEFiled(component::FieldData * pdata, const VectorR3 & rad_k, const VectorR3 & rad_ev, const VectorR3 & rad_eh) const
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
		ce_ptr_->getCommonEdge(u, nvp, nvm, fp, fm, L);
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

void FSPGF_EFIE_HS::calculateSurfaceCurrent(std::vector<component::CurrentData>* currents) const
{
	Qcx_vec cur_vec = I;
	const size_t tri_num = mesh_ptr_->getTriangleNum();
	std::vector<std::vector<int>> triangles(tri_num);
	for (int u = 0; u < unknowns_; ++u)
	{
		auto& rwg = ce_ptr_->getRWGRef(u);
		triangles[rwg.tpid].push_back(u);
		triangles[rwg.tmid].push_back(u);
	}
	currents->resize(tri_num);
	VectorR3 vx[3];
	VectorC3 mag[3];
	value_t len[3];
	for (size_t t = 0; t < tri_num; ++t)
	{
		auto& unks = triangles[t];
		for (size_t i = 0; i < unks.size(); ++i)
		{
			auto& rwg = ce_ptr_->getRWGRef(unks[i]);
			vx[i] = (t == rwg.tpid) ? mesh_ptr_->getVertex(rwg.vxp) : mesh_ptr_->getVertex(rwg.vxm);
			len[i] = rwg.length;
		}
		//处理边缘RWG
		if (unks.size() == 1)
		{
			auto& rwg = ce_ptr_->getRWGRef(unks[0]);
			vx[1] = mesh_ptr_->getVertex(rwg.v1);
			vx[2] = mesh_ptr_->getVertex(rwg.v2);
		}
		if (unks.size() == 2)
		{
			auto& rwg1 = ce_ptr_->getRWGRef(unks[0]);
			auto& rwg2 = ce_ptr_->getRWGRef(unks[1]);
			if ((rwg1.v1 != rwg2.v1) && (rwg1.v1 != rwg2.v2))
			{
				vx[2] = mesh_ptr_->getVertex(rwg1.v2);
			}
			else
				vx[2] = mesh_ptr_->getVertex(rwg1.v1);
		}
		auto center = (vx[0] + vx[1] + vx[2]) / 3;
		auto dominator = 2 * Area(vx[0], vx[1], vx[2]);
		VectorC3 cen_cur(0, 0, 0);
		for (size_t i = 0; i < unks.size(); ++i)
		{
			auto& rwg = ce_ptr_->getRWGRef(unks[i]);
			cen_cur += (t == rwg.tpid ? 1.0f : -1.0f) * cur_vec(unks[i]) * (len[i] / dominator) * (center - vx[i]);
			VectorC3 ncur(0, 0, 0);
			for (size_t j = 0; j < unks.size(); ++j)
			{
				if (i == j) continue;
				auto& srwg = ce_ptr_->getRWGRef(unks[j]);
				ncur += (t == srwg.tpid ? 1.0f : -1.0f) * cur_vec(unks[j]) * (len[j] / dominator) * (vx[i] - vx[j]);
			}
			mag[i] = ncur;
		}
		auto& data = (*currents)[t];
		data.v1 = vx[0];
		data.v2 = vx[1];
		data.v3 = vx[2];
		data.magnc = std::sqrt(cen_cur.norm());
		data.magn1 = std::sqrt(mag[0].norm());
		data.magn2 = std::sqrt(mag[1].norm());
		data.magn3 = std::sqrt(mag[2].norm());
	}
}

bool FSPGF_EFIE_HS::writeZIVData()
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