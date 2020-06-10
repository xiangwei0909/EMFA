#include "IBC.h"
#include "Miscellaneous.h"
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


IBC::IBC()
{

}

IBC::~IBC()
{
}

void IBC::init(component::ConfigLoader * ploader)
{
	EM::init(ploader);
	if (rtype_ == policy::ResultType::Rad)
		readExcEdges(dir_ + '/' + folder_name_ + ".rad");

	dir_ = tool::creatFolder(dir_, folder_name_ + "_IBC");

	auto log_name = dir_ + (rtype_ == policy::ResultType::Rad ? "/rad_" : "/sca_") + "runtime.log";
	logger_.open(log_name);

	mesh_ptr_ = std::make_shared<Mesh>();
	ce_ptr_ = std::make_shared<CommonEdge>();

	mesh_ptr_->loadIBCMesh(ploader->getMeshPath());
	ce_ptr_->buildIBCCommonEdge(mesh_ptr_);
	Sigma = ploader->getSigma();
	//////////////////////////////////////////////////////////
	eps = ploader->getEpsilon2();
	thickness = ploader->getThickness_sheet();
	Zs = -J0 * sqrt(1.0f / eps)*tan(sqrt(eps)*k_*thickness);
	//////////////////////////////////////////////////////////
	unknowns_ = ce_ptr_->getCommonEdgeNum();
	k_ = PI2 * incidence_.freq / cc;
	logger_ << ploader->getReport();
	reportInfo(logger_);
	TIME_LOG("init");
}

void IBC::solve()
{
	SEGMENT("Solve");
	Z.setZero(unknowns_, unknowns_);
	LOG(fillZ(), "Filling impedance matrix");
	TIME_LOG("fillZ");

	I.setZero(unknowns_);
	V.setZero(unknowns_);

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
	I = Z.partialPivLu().solve(V);
	//if (!arma::solve(I, Z, V, arma::solve_opts::fast))
		//throw std::runtime_error("fail solving matrix equation");
	if (writeZIVData())
		Qcout << "The data is already saved!" << std::endl;
	Qcout << "success" << std::endl;
	TIME_LOG("solve");

#ifdef _DEBUG
	//if (!writeZIVData())
	//  Qcerr << "fail writing ZIV data" << std::endl;
#endif
}

void IBC::output()
{
	SEGMENT("Result");
	Result result;
	if (rtype_ == policy::ResultType::Rad)
	{
		LOG(result.getRadiationField(this, &IBC::getEFiled), "Calculating Radiation");
		TIME_LOG("getRadiation");
		LOG(getCurrentOnFeed(), "Calculating the Current on the feed");
	}
	else
	{
		LOG(result.getBistaticRCS(this, &IBC::getBiRCS), "Calculating BistaticRCS");
		TIME_LOG("getBistaticRCS");
	}
	LOG(result.getCurrentDistribution(this, &IBC::calculateSurfaceCurrent), "Calculating surface current");
}

void IBC::clear()
{
	Z.resize(0, 0);
	I.resize(0, 0);
	V.resize(0, 0);
	
	exc_edge_.clear();

	mesh_ptr_->clear();
	ce_ptr_->clear();
}


void IBC::reportInfo(Qostream& strm) const
{
	mesh_ptr_->reportInfo(strm);
	ce_ptr_->reportInfo(strm);
}

Complex IBC::eZppKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex IBC::eZpmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex IBC::eZmpKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex IBC::eZmmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex IBC::eZppSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0);
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
			value_t r02 = p0 * p0 + d * d;
			Is1 += (0.5f*ui)*(r02*log((rp + lp) / (rm + lm)) + (lp*rp) - (lm*rm));
			Is2 += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
		}
		Ze2 += w3[p] * (0.25f * (rou_f ^ (Is1 + Is2 * (vgf[p] - vs))) - (Is2 / (k_ * k_)));
	}

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
	Ze = (Ze1 + (1.0f / S_s) * Ze2) + Ze3;

	return Ze;
}

Complex IBC::eZpmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0);
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
		Ze2 += w3[p] * (0.25f * (-1.0f*rou_f ^ (Is1 + Is2 * (vgf[p] - vs))) + (Is2 / (k_ * k_)));
	}

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

	Ze = (Ze1 + (1.0f / S_s) * Ze2) + Ze3;

	return Ze;
}

Complex IBC::eZmpSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0);
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


	if (abs(Sigma) > 1e-8)
	{
		for (int i = 0; i < 3; i++)
		{
			rou_f = vgf[i] - vf;
			rou_s = vgf[i] - vs;
			Ze3 += w3[i] * rou_f^rou_s;
		}
		Ze3 = (PI*Sigma / (S_s*J0*k_*Z0))*Ze3;
	}

	Ze = (Ze1 + (1.0f / S_s) * Ze2) + Ze3;

	return Ze;
}

Complex IBC::eZmmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0);
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

	if (abs(Sigma) > 1e-8)
	{
		for (int i = 0; i < 3; i++)
		{
			rou_f = vgf[i] - vf;
			rou_s = vgf[i] - vs;
			Ze3 += w3[i] * rou_f^rou_s;
		}
		Ze3 = (PI*Sigma / (S_s*J0*k_*Z0))*Ze3;
	}

	Ze = (Ze1 + (1.0f / S_s) * Ze2) + Ze3;

	return Ze;
}

Complex IBC::mZppKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf)
{
	value_t r;
	Complex  Zm(0, 0), G(0, 0);
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
			Zm += w3[p] * w3[q] * (rou_f ^ (nf * ((vgf[p] - vgs[q]) * rou_s))) * (value_t(1.0) + J0 * k_ * r) * G / (r * r);
		}
	}
	Zm = Zm / (16 * PI);
	return Zm;
}

Complex IBC::mZpmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf)
{
	value_t r;
	Complex Zm(0, 0), G(0, 0);
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
			Zm += w3[p] * w3[q] * (rou_f ^ (nf * ((vgf[p] - vgs[q]) * rou_s))) * (value_t(1.0) + J0 * k_ * r) * G / (r * r);
		}
	}
	Zm = Zm / (16 * PI);
	return Zm;
}

Complex IBC::mZmpKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf)
{
	value_t r;
	Complex Zm(0, 0), G(0, 0);
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
			Zm += w3[p] * w3[q] * (rou_f ^ (nf * ((vgf[p] - vgs[q]) * rou_s))) * (value_t(1.0) + J0 * k_ * r) * G / (r * r);
		}
	}
	Zm = Zm / (16 * PI);
	return  Zm;
}

Complex IBC::mZmmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf)
{
	value_t r;
	Complex Zm(0, 0), G(0, 0);
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

			Zm += w3[p] * w3[q] * (rou_f ^ (nf * ((vgf[p] - vgs[q]) * rou_s))) * (value_t(1.0) + J0 * k_ * r) * G / (r * r);
		}
	}
	Zm = Zm / (16 * PI);
	return Zm;
}

Complex IBC::mZppSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	Complex Zm(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_f = Area(vf3[0], vf3[1], vf3[2]);

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	for (int p = 0; p < 3; p++)
	{
		//磁场 奇异部分
		rou_f = vgf[p] - vf;
		rou_s = vgs[p] - vs;
		Zm += w3[p] * (rou_f ^ rou_s);
	}

	Zm = Zm / (value_t(8.0) * S_f);

	return Zm;
}

Complex IBC::mZpmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	Complex Zm(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_f = Area(vf3[0], vf3[1], vf3[2]);

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	for (int p = 0; p < 3; p++)
	{
		//磁场 奇异部分
		rou_f = vgf[p] - vf;
		rou_s = vs - vgs[p];
		Zm += w3[p] * (rou_f ^ rou_s);
	}

	Zm = Zm / (value_t(8.0) * S_f);

	return Zm;
}

Complex IBC::mZmpSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	Complex Zm(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_f = Area(vf3[0], vf3[1], vf3[2]);

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//////////////////////////////////////////////////////////////////////////
	for (int p = 0; p < 3; p++)
	{
		//磁场 奇异部分
		rou_f = vf - vgf[p];
		rou_s = vgs[p] - vs;
		Zm += w3[p] * (rou_f ^ rou_s);
	}

	Zm = Zm / (value_t(8.0) * S_f);

	return Zm;
}

Complex IBC::mZmmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
	Complex Zm(0, 0);
	VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	value_t S_f = Area(vf3[0], vf3[1], vf3[2]);

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	for (int p = 0; p < 3; p++)
	{
		//磁场 奇异部分
		rou_f = vf - vgf[p];
		rou_s = vs - vgs[p];
		Zm += w3[p] * (rou_f ^ rou_s);
	}

	Zm = Zm / (value_t(8.0) * S_f);

	return Zm;
}

Complex IBC::cZppKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf, value_t &msign, value_t& dsign)
{
	return integral::cfieZKernel(k_, vf3, vs3, vf, vs, nf, 1, 1, msign, dsign);
	//value_t r;
	//Complex Ze(0, 0), Zm(0, 0), G(0, 0);
	//VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	//Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	//Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//for (int p = 0; p < 3; p++)
	//{
	//    rou_f = vgf[p] - vf;
	//    for (int q = 0; q < 3; q++)
	//    {
	//        rou_s = vgs[q] - vs;
	//        r = (vgf[p] - vgs[q]).Norm();
	//        G = exp(-J0 * k_ * r) / r;
	//        Ze += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * G;
	//        Zm += w3[p] * w3[q] * (rou_f ^ (nf * ((vgf[p] - vgs[q]) * rou_s))) * (1.0f + J0 * k_ * r) * G / (r * r);
	//    }
	//}
	//Ze = (J0 * k_ * Z0 / PI4) * Ze;
	//Zm = Zm / (16 * PI);
	//return alpha_ * Ze + (1 - alpha_) * Z0 * Zm;
}

Complex IBC::cZpmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf, value_t &msign, value_t& dsign)
{
	return integral::cfieZKernel(k_, vf3, vs3, vf, vs, nf, 1, -1, msign, dsign);
	//value_t r;
	//Complex Ze(0, 0), Zm(0, 0), G(0, 0);
	//VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	//Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	//Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//for (int p = 0; p < 3; p++)
	//{
	//    rou_f = vgf[p] - vf;
	//    for (int q = 0; q < 3; q++)
	//    {
	//        rou_s = vs - vgs[q];
	//        r = (vgf[p] - vgs[q]).Norm();
	//        G = exp(-J0 * k_ * r) / r;
	//        Ze += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * G;
	//        Zm += w3[p] * w3[q] * (rou_f ^ (nf * ((vgf[p] - vgs[q]) * rou_s))) * (1.0f + J0 * k_ * r) * G / (r * r);
	//    }
	//}
	//Ze = (J0 * k_ * Z0 / PI4) * Ze;
	//Zm = Zm / (16 * PI);
	//return alpha_ * Ze + (1 - alpha_) * Z0 * Zm;
}

Complex IBC::cZmpKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf, value_t &msign, value_t& dsign)
{
	return integral::cfieZKernel(k_, vf3, vs3, vf, vs, nf, -1, 1, msign, dsign);
	//value_t r;
	//Complex Ze(0, 0), Zm(0, 0), G(0, 0);
	//VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	//Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	//Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//for (int p = 0; p < 3; p++)
	//{
	//    rou_f = vf - vgf[p];
	//    for (int q = 0; q < 3; q++)
	//    {
	//        rou_s = vgs[q] - vs;
	//        r = (vgf[p] - vgs[q]).Norm();
	//        G = exp(-J0 * k_ * r) / r;
	//        Ze += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * G;
	//        Zm += w3[p] * w3[q] * (rou_f ^ (nf * ((vgf[p] - vgs[q]) * rou_s))) * (1.0f + J0 * k_ * r) * G / (r * r);
	//    }
	//}
	//Ze = (J0 * k_ * Z0 / PI4) * Ze;
	//Zm = Zm / (16 * PI);
	//return alpha_ * Ze + (1 - alpha_) * Z0 * Zm;
}

Complex IBC::cZmmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf, value_t &msign, value_t& dsign)
{
	return integral::cfieZKernel(k_, vf3, vs3, vf, vs, nf, -1, -1, msign, dsign);
	//value_t r;
	//Complex Ze(0, 0), Zm(0, 0), G(0, 0);
	//VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	//Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	//Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	//for (int p = 0; p < 3; p++)
	//{
	//    rou_f = vf - vgf[p];
	//    for (int q = 0; q < 3; q++)
	//    {
	//        rou_s = vs - vgs[q];
	//        r = (vgf[p] - vgs[q]).Norm();
	//        G = exp(-J0 * k_ * r) / r;

	//        Ze += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * G;
	//        Zm += w3[p] * w3[q] * (rou_f ^ (nf * ((vgf[p] - vgs[q]) * rou_s))) * (1.0f + J0 * k_ * r) * G / (r * r);
	//    }
	//}
	//Ze = (J0 * k_ * Z0 / PI4) * Ze;
	//Zm = Zm / (16 * PI);
	//return alpha_ * Ze + (1 - alpha_) * Z0 * Zm;
}

Complex IBC::cZppSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, value_t &msign, value_t& dsign)
{
	return integral::cfieZSingular(k_, vf3, vs3, vf, vs, 1, 1, msign, dsign);
	//Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Zm(0, 0);
	//VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	//value_t S_f = Area(vf3[0], vf3[1], vf3[2]);
	//value_t S_s = Area(vs3[0], vs3[1], vs3[2]);

	//Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	//Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	////////////////////////////////////////////////////////////////////////////
	//for (int p = 0; p < 3; p++)
	//{
	//    //电场  非奇异部分
	//    rou_f = vgf[p] - vf;
	//    for (int q = 0; q < 3; q++)
	//    {
	//        rou_s = vgs[q] - vs; //.2.
	//        value_t r = (vgf[p] - vgs[q]).Norm();
	//        Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
	//        Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * coff; //.3.
	//    }
	//    //电场  奇异部分
	//    for (int i = 0; i < 3; i++)
	//    {
	//        VectorR3 Li = vs3[(i + 1) % 3] - vs3[i];
	//        Li.Normalize();//边向量

	//        value_t Lpi = (vs3[(i + 1) % 3] - vgf[p]) ^ Li; //(边正顶点 - 场点)·(边向量)
	//        value_t Lmi = (vs3[i] - vgf[p]) ^ Li; //(边负顶点 - 场点)·(边向量)
	//        value_t Rpi = (vs3[(i + 1) % 3] - vgf[p]).Norm(); //边正顶点 到 场点 的距离
	//        value_t Rmi = (vs3[i % 3] - vgf[p]).Norm(); //边负顶点 到 场点 的距离
	//        VectorR3 sj0 = vs3[i] + Li * ((vgf[p] - vs3[i]) ^ Li); //垂足在u方向位置矢量
	//        value_t P0i = (sj0 - vgf[p]).Norm(); //垂足 到 场点 的距离
	//        if (P0i > 1e-10f)
	//            Ze2 += w3[p] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * P0i * log((Rpi + Lpi) / (Rmi + Lmi)); //.4.
	//    }

	//    //磁场 奇异部分
	//    rou_f = vgf[p] - vf;
	//    rou_s = vgs[p] - vs;
	//    Zm += w3[p] * (rou_f ^ rou_s);
	//}

	//Ze = (J0 * k_ * Z0 / PI4) * (Ze1 + 1.0f / S_s * Ze2);
	//Zm = Zm / (8.0f * S_f);

	//return alpha_ * Ze + (1 - alpha_) * Z0 * Zm;
}

Complex IBC::cZpmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, value_t &msign, value_t& dsign)
{
	return integral::cfieZSingular(k_, vf3, vs3, vf, vs, 1, -1, msign, dsign);
	//Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Zm(0, 0);
	//VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	//value_t S_f = Area(vf3[0], vf3[1], vf3[2]);
	//value_t S_s = Area(vs3[0], vs3[1], vs3[2]);

	//Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	//Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	////////////////////////////////////////////////////////////////////////////
	//for (int p = 0; p < 3; p++)
	//{
	//    //电场  非奇异部分
	//    rou_f = vgf[p] - vf; //.2.
	//    for (int q = 0; q < 3; q++)
	//    {
	//        rou_s = vs - vgs[q]; //.2.
	//        value_t r = (vgf[p] - vgs[q]).Norm();
	//        Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
	//        Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * coff; //.3.
	//    }
	//    //电场  奇异部分
	//    for (int i = 0; i < 3; i++)
	//    {
	//        VectorR3 Li = vs3[(i + 1) % 3] - vs3[i];
	//        Li.Normalize();//边向量

	//        value_t Lpi = (vs3[(i + 1) % 3] - vgf[p]) ^ Li; //(边正顶点 - 场点)·(边向量)
	//        value_t Lmi = (vs3[i] - vgf[p]) ^ Li; //(边负顶点 - 场点)·(边向量)
	//        value_t Rpi = (vs3[(i + 1) % 3] - vgf[p]).Norm(); //边正顶点 到 场点 的距离
	//        value_t Rmi = (vs3[i % 3] - vgf[p]).Norm(); //边负顶点 到 场点 的距离
	//        VectorR3 sj0 = vs3[i] + Li * ((vgf[p] - vs3[i]) ^ Li); //垂足在u方向位置矢量
	//        value_t P0i = (sj0 - vgf[p]).Norm(); //垂足 到 场点 的距离
	//        if (P0i > 1e-10)
	//            Ze2 += w3[p] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * P0i * log((Rpi + Lpi) / (Rmi + Lmi)); //.4.
	//    }

	//    //磁场 奇异部分
	//    rou_f = vgf[p] - vf;
	//    rou_s = vs - vgs[p];
	//    Zm += w3[p] * (rou_f ^ rou_s);
	//}

	//Ze = (J0 * k_ * Z0 / PI4) * (Ze1 + 1.0f / S_s * Ze2);//fix by wangjian 2015
	//Zm = Zm / (8.0f * S_f);

	//return alpha_ * Ze + (1 - alpha_) * Z0 * Zm;
}

Complex IBC::cZmpSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, value_t &msign, value_t& dsign)
{
	return integral::cfieZSingular(k_, vf3, vs3, vf, vs, -1, 1, msign, dsign);
	//Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Zm(0, 0);
	//VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	//value_t S_f = Area(vf3[0], vf3[1], vf3[2]);
	//value_t S_s = Area(vs3[0], vs3[1], vs3[2]);

	//Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	//Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	////////////////////////////////////////////////////////////////////////////
	//for (int p = 0; p < 3; p++)
	//{
	//    //电场  非奇异部分
	//    rou_f = vf - vgf[p];
	//    for (int q = 0; q < 3; q++)
	//    {
	//        rou_s = vgs[q] - vs; //.2.
	//        value_t r = (vgf[p] - vgs[q]).Norm();
	//        Complex coff = -0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
	//        Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * coff; //.3.
	//    }
	//    //电场  奇异部分
	//    for (int i = 0; i < 3; i++)
	//    {
	//        VectorR3 Li = vs3[(i + 1) % 3] - vs3[i];
	//        Li.Normalize();//边向量

	//        value_t Lpi = (vs3[(i + 1) % 3] - vgf[p]) ^ Li; //(边正顶点 - 场点)·(边向量)
	//        value_t Lmi = (vs3[i] - vgf[p]) ^ Li; //(边负顶点 - 场点)·(边向量)
	//        value_t Rpi = (vs3[(i + 1) % 3] - vgf[p]).Norm(); //边正顶点 到 场点 的距离
	//        value_t Rmi = (vs3[i % 3] - vgf[p]).Norm(); //边负顶点 到 场点 的距离
	//        VectorR3 sj0 = vs3[i] + Li * ((vgf[p] - vs3[i]) ^ Li); //垂足在u方向位置矢量
	//        value_t P0i = (sj0 - vgf[p]).Norm(); //垂足 到 场点 的距离
	//        if (P0i > 1e-10)
	//            Ze2 += w3[p] * (0.25f * (rou_f ^ rou_s) + 1.0f / (k_ * k_)) * P0i * log((Rpi + Lpi) / (Rmi + Lmi)); //.4.
	//    }

	//    //磁场 奇异部分
	//    rou_f = vf - vgf[p];
	//    rou_s = vgs[p] - vs;
	//    Zm += w3[p] * (rou_f ^ rou_s);
	//}

	////Ze = Ze1 + 0.5 / S_s * Ze2;
	//Ze = (J0 * k_ * Z0 / PI4) * (Ze1 + 1.0f / S_s * Ze2);
	//Zm = Zm / (8.0f * S_f);

	//return alpha_ * Ze + (1 - alpha_) * Z0 * Zm;
}

Complex IBC::cZmmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, value_t &msign, value_t& dsign)
{
	return integral::cfieZSingular(k_, vf3, vs3, vf, vs, -1, -1, msign, dsign);
	//Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Zm(0, 0);
	//VectorR3 vgf[3], vgs[3], rou_f, rou_s;

	//value_t S_f = Area(vf3[0], vf3[1], vf3[2]);
	//value_t S_s = Area(vs3[0], vs3[1], vs3[2]);


	//Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	//Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	////////////////////////////////////////////////////////////////////////////
	//for (int p = 0; p < 3; p++)
	//{
	//    //电场  非奇异部分
	//    rou_f = vf - vgf[p];
	//    for (int q = 0; q < 3; q++)
	//    {
	//        rou_s = vs - vgs[q]; //.2.
	//        value_t r = (vgf[p] - vgs[q]).Norm();
	//        Complex coff = - 0.5f * k_ * k_ * r + J0 * k_ / 6.0f * (k_ * k_ * r * r - 6);
	//        Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * coff; //.3.
	//    }
	//    //电场  奇异部分
	//    for (int i = 0; i < 3; i++)
	//    {
	//        VectorR3 Li = vs3[(i + 1) % 3] - vs3[i];
	//        Li.Normalize();//边向量

	//        value_t Lpi = (vs3[(i + 1) % 3] - vgf[p]) ^ Li; //(边正顶点 - 场点)·(边向量)
	//        value_t Lmi = (vs3[i] - vgf[p]) ^ Li; //(边负顶点 - 场点)·(边向量)
	//        value_t Rpi = (vs3[(i + 1) % 3] - vgf[p]).Norm(); //边正顶点 到 场点 的距离
	//        value_t Rmi = (vs3[i % 3] - vgf[p]).Norm(); //边负顶点 到 场点 的距离
	//        VectorR3 sj0 = vs3[i] + Li * ((vgf[p] - vs3[i]) ^ Li); //垂足在u方向位置矢量
	//        value_t P0i = (sj0 - vgf[p]).Norm(); //垂足 到 场点 的距离
	//        if (P0i > 1e-10)
	//            Ze2 += w3[p] * (0.25f * (rou_f ^ rou_s) - 1.0f / (k_ * k_)) * P0i * log((Rpi + Lpi) / (Rmi + Lmi)); //.4.
	//    }

	//    //磁场 奇异部分
	//    rou_f = vf - vgf[p];
	//    rou_s = vs - vgs[p];
	//    Zm += w3[p] * (rou_f ^ rou_s);
	//}

	////Ze = Ze1 + 0.5 / S_s * Ze2;
	//Ze = (J0 * k_ * Z0 / PI4) * (Ze1 + 1.0f / S_s * Ze2);
	//Zm = Zm / (8.0f * S_f);

	//return alpha_ * Ze + (1 - alpha_) * Z0 * Zm;
}

Complex IBC::eVKernel(VectorR3 *vp3, VectorR3 *vm3, VectorR3 &vp, VectorR3 &vm)
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

void IBC::fillZ()
{
	Z.setZero(unknowns_, unknowns_);

	int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle
	VectorR3 v_fld_plu, v_fld_min, v_fld_plu3[3], v_fld_min3[3];
	VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
	Triangle tri_plu, tri_min;
	value_t l_fld, l_src;//场与源的公共边长度
	value_t msign, dsign;
	Qstring fla_fld_plu, fla_fld_min, fla_src_plu, fla_src_min;
	Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
	int vp, vm;

	tool::BarAndPercent bar_perc;
	for (int f = 0; f < unknowns_; f++)
	{
		bar_perc(f + 1, unknowns_);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld, fla_fld_plu, fla_fld_min);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		v_fld_min = mesh_ptr_->getVertex(vm);//场点-

		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_min = mesh_ptr_->getTriangleRef(f_fld_min);//场三角面-

		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);

		VectorR3 nml_plu = tri_plu.getNormal();
		VectorR3 nml_min = tri_min.getNormal();

		for (int s = 0; s < unknowns_; s++)
		{
			ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src, fla_src_plu, fla_src_min);

			v_src_plu = mesh_ptr_->getVertex(vp);
			v_src_min = mesh_ptr_->getVertex(vm);

			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
			tri_min = mesh_ptr_->getTriangleRef(f_src_min);

			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);

			msign = 1;
			dsign = 1;
			///////////////////////////////////////////////////////////////////////
			//if(fla_src_plu)
			//field face+ <-->  source face+
			if (f_fld_plu == f_src_plu)
			{
				Zpp = cZppSingular(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu, msign, dsign);
			}
			else
			{
				Zpp = cZppKernel(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu, nml_plu, msign, dsign);
			}
			//field face+ <--> source face-
			if (f_fld_plu == f_src_min)
			{
				Zpm = cZpmSingular(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min, msign, dsign);
			}
			else
			{
				Zpm = cZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min, nml_plu, msign, dsign);
			}
			//field face- <--> source face+
			if (f_fld_min == f_src_plu)
			{
				Zmp = cZmpSingular(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu, msign, dsign);
			}
			else
			{
				Zmp = cZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu, nml_min, msign, dsign);
			}
			//field face- <--> source face-
			if (f_fld_min == f_src_min)
			{
				Zmm = cZmmSingular(v_fld_min3, v_src_min3, v_fld_min, v_src_min, msign, dsign);
			}
			else
			{
				Zmm = cZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min, nml_min, msign, dsign);
			}
			Z(f, s) = l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
		}
	}
}

void IBC::fillV()
{
	int nvp, nvm, fp, fm;
	VectorR3 vp, vm, vp3[3], vm3[3];
	value_t ln;
	Triangle tri_plu, tri_min;

	tool::BarAndPercent bar_perc;   //
	for (int u = 0; u < unknowns_; u++)
	{
		bar_perc(u + 1, unknowns_);    //

		ce_ptr_->getCommonEdge(u, nvp, nvm, fp, fm, ln);
		vp = mesh_ptr_->getVertex(nvp);
		vm = mesh_ptr_->getVertex(nvm);
		tri_plu = mesh_ptr_->getTriangleRef(fp);
		tri_min = mesh_ptr_->getTriangleRef(fm);
		tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);
		tri_min.getVertex(vm3[0], vm3[1], vm3[2]);

		Complex vk = eVKernel(vp3, vm3, vp, vm);
		V(u) = 0.5f * ln * vk;
	}
}

bool IBC::readExcEdges(const Qstring & rad_file, Qstring& stateInfo)
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

void IBC::readExcEdges(const Qstring & rad_file)
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

void IBC::radiateV()
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

value_t IBC::getBiRCS(const VectorR3 & sca_k) const
{
	int nvp, nvm, fp, fm;
	VectorR3 vp, vm, vp3[3], vm3[3];
	value_t L;

	VectorR3 vgp[7], vgm[7];

	VectorR3 rou_p, rou_m;
	Complex G0(0, 0);
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

		VectorC3 Esp(0, 0, 0), Esm(0, 0, 0);
		for (int g = 0; g < 7; g++)
		{
			rou_p = vgp[g] - vp;
			G0 = exp(J0 * k_ * (vgp[g] ^ sca_k));
			Esp = Esp + rou_p * w7[g] * G0;

			rou_m = vm - vgm[g];
			G0 = exp(J0 * k_ * (vgm[g] ^ sca_k));
			Esm = Esm + rou_m * w7[g] * G0;
		}
		Es = Es + (Esp + Esm) * (I[u] * L);
	}

	auto es = 0.5f * sca_k * (sca_k * Es);
	value_t rcs = Z0 * Z0 * k_ * k_ * es.norm() / PI4;

	//return rcs;
	return 10 * log10(rcs);
}

void IBC::getEFiled(component::FieldData * pdata, const VectorR3 & rad_k, const VectorR3 & rad_ev, const VectorR3 & rad_eh) const
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

void IBC::getCurrentOnFeed()
{
	int nv1, nv2;
	Qofstream outputI(dir_ + "/I_feed.txt");
	for (int u = 0; u < unknowns_; u++)
	{
		ce_ptr_->getCommonEdge(u, nv1, nv2);
		for (const auto& elem : exc_edge_)
		{
			if ((elem.first == nv1 && elem.second == nv2) || (elem.first == nv2 && elem.second == nv1))
			{
				Complex i_temp;
				auto length = ce_ptr_->getCommonEdgeLength(u);
				if (elem.first == nv1)
					i_temp = I(u)*length;
				else
					i_temp = I(u)*length;
				outputI << std::setw(18) << i_temp.real() << std::setw(18) << i_temp.imag() << '\n';
			}
		}
	}
	outputI.flush();
	outputI.close();
}

bool IBC::writeZIVData()
{
	Qofstream outputZ(dir_ + "/matrix_Z.txt");
	if (outputZ.fail())
		return false;
	outputZ << Z;
	outputZ.flush();
	outputZ.close();

	Qofstream outputV(dir_ + "/matrix_V.txt");
	if (outputV.fail())
		return false;
	outputV << V;
	outputV.flush();
	outputV.close();

	Qofstream outputI(dir_ + "/martix_I.txt");
	if (outputI.fail())
		return false;
	outputI << I;
	outputI.flush();
	outputI.close();

	return true;
}

void IBC::calculateSurfaceCurrent(std::vector<component::CurrentData>* currents) const
{
	CVector cur_vec = I;
	const size_t tri_num = mesh_ptr_->getTriangleNum();
	const size_t node_num = mesh_ptr_->getNodeNum();
	std::vector<std::vector<value_t>> node(node_num);
	std::vector<std::vector<int>> triangle(tri_num);
	for (int u = 0; u < unknowns_; ++u)
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
		//处理边缘RWG
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
