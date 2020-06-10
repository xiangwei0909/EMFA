#include "TDS.h"
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


TDS::TDS()
{

}

TDS::~TDS()
{
}

void TDS::init(component::ConfigLoader * ploader)
{
	EM::init(ploader);
	if (rtype_ == policy::ResultType::Rad)
		readExcEdges(dir_ + '/' + folder_name_ + ".rad");

	dir_ = tool::creatFolder(dir_, folder_name_ + "_TDS");

	auto log_name = dir_ + (rtype_ == policy::ResultType::Rad ? "/rad_" : "/sca_") + "runtime.log";
	logger_.open(log_name);

	mesh_ptr_ = std::make_shared<Mesh>();
	ce_ptr_ = std::make_shared<CommonEdge>();

	mesh_ptr_->loadTDSMesh(ploader->getMeshPath());
	ce_ptr_->buildTDSCommonEdge(mesh_ptr_);
	Sigma = ploader->getSigma();
	omiga = PI2 * incidence_.freq;
	Thickness = ploader->getThickness_sheet();
	unknowns_t = ce_ptr_->getCommonEdgeNum();
	unknowns_n = ce_ptr_->getPluseTriNum();
	k_ = PI2 * incidence_.freq / cc;
	logger_ << ploader->getReport();
	reportInfo(logger_);
	TIME_LOG("init");
}

void TDS::solve()
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

void TDS::output()
{
	SEGMENT("Result");
	Result result;
	if (rtype_ == policy::ResultType::Rad)
	{
		LOG(result.getRadiationField(this, &TDS::getEFiled), "Calculating Radiation");
		TIME_LOG("getRadiation");
	}
	else
	{
		LOG(result.getBistaticRCS(this, &TDS::getBiRCS), "Calculating BistaticRCS");
		TIME_LOG("getBistaticRCS");
	}
	LOG(result.getCurrentDistribution(this, &TDS::calculateSurfaceCurrent), "Calculating surface current");
}

void TDS::clear()
{
	Z.reset();
	I.reset();
	V.reset();
	exc_edge_.clear();

	mesh_ptr_->clear();
	ce_ptr_->clear();
}


void TDS::reportInfo(Qostream& strm) const
{
	mesh_ptr_->reportInfo(strm);
	ce_ptr_->reportInfo(strm);
}

Complex TDS::eZppKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex TDS::eZpmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex TDS::eZmpKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex TDS::eZmmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex TDS::eZppSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex TDS::eZpmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex TDS::eZmpSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex TDS::eZmmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex TDS::TTZppKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &f_nor, VectorR3 &s_nor, VectorR3 *vf2, VectorR3 *vs2, value_t &eps_p, value_t &eps_m, int &bedge_f, int &bedge_s)
{
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0), Ze4(0, 0), Ze5(0, 0), Ze6(0, 0);
	Complex  G(0, 0), G_up(0, 0), G_c(0, 0), G_c_up(0, 0), G_c_c(0, 0);
	VectorR3 vf3_tds[3], vs3_tds[3],vs3_up[3], vf2_c[2], vs2_c[2],vf_tds,v;
	VectorR3 vgf[3], vgs[3], vgf_tds[3], vgs_tds[3], vgs_up[3], vgf_c[3], vgs_c[3];
	VectorR3 v2_m[2];
	VectorR3 rou_f, rou_s;
	value_t r,r_up,r_c,r_c_up,r_c_c;
	value_t k_p, k_m;

	TransTri(vf3, vs3_up, s_nor, 1.0f);
	TransTri(vf3, vf3_tds, s_nor, 0.5f);
	TransTri(vs3, vs3_tds, s_nor, 0.5f);
	TransLine(vf2, vf2_c, s_nor, 0.5f);
	TransLine(vs2, vs2_c, s_nor, 0.5f);

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	Gauss3Point(vf3_tds[0], vf3_tds[1], vf3_tds[2], vgf_tds);
	GaussLine3Point(vf2_c[0], vf2_c[1], vgf_c);

	Gauss3Point(vs3_tds[0], vs3_tds[1], vs3_tds[2], vgs_tds);
	Gauss3Point(vs3_up[0], vs3_up[1], vs3_up[2], vgs_up);
	GaussLine3Point(vs2_c[0], vs2_c[1], vgs_c);

	

	k_p = (1.0f / eps_p) - 1.0f;
	k_m = (1.0f / eps_m) - 1.0f;

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs;

			r = (vgf_tds[p] - vgs_tds[q]).Norm();
			r_up = (vgf_tds[p] - vgs_up[q]).Norm();
			r_c = (vgf_tds[p] - vgs_c[q]).Norm();
			r_c_up = (vgf_c[p] - vgs_up[q]).Norm();
			r_c_c = (vgf_c[p] - vgs_c[q]).Norm();

			G = exp(-J0 * k_*r) / r;
			G_up = exp(-J0 * k_*r_up) / r_up;
			G_c = exp(-J0 * k_*r_c) / r_c;
			G_c_up = exp(-J0 * k_*r_c_up) / r_c_up;
			G_c_c = exp(-J0 * k_*r_c_c) / r_c_c;

			Ze2 += w3[p] * w3[q] * (rou_f^(rou_s-(Thickness*s_nor)))*G;
			Ze3 += w3[p] * w3[q] * G_up;
			Ze4 += w3[p] * glw3[q] * G_c;

			if (bedge_f == 1)
			{
				Ze5 += glw3[p] * w3[q] * G_c_up;
				Ze6 += glw3[p] * glw3[q] * G_c_c;
			}
		}
	}

	Ze = (k_p * Thickness*Ze2/4.0f) + ((k_p*Thickness) / (k_*k_))*(Ze5 - Ze3) + ((k_m - k_p)*Thickness / (k_*k_))*(Ze6 - Ze4);

	return Ze;
}

Complex TDS::TTZpmKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &f_nor, VectorR3 &s_nor, VectorR3 *vf2, VectorR3 *vs2, value_t &eps_p, value_t &eps_m, int &bedge_f, int &bedge_s)
{
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0), Ze4(0, 0), Ze5(0, 0), Ze6(0, 0);
	Complex  G(0, 0), G_up(0, 0), G_c(0, 0), G_c_up(0, 0), G_c_c(0, 0);
	VectorR3 vf3_tds[3], vs3_tds[3], vs3_up[3], vf2_c[2], vs2_c[2], vf_tds, v;
	VectorR3 vgf[3], vgs[3], vgf_tds[3], vgs_tds[3], vgs_up[3], vgf_c[3], vgs_c[3];
	VectorR3 v2_m[2];
	VectorR3 rou_f, rou_s;
	value_t r, r_up, r_c, r_c_up, r_c_c;
	value_t k_p, k_m;

	TransTri(vf3, vs3_up, s_nor, 1.0f);
	TransTri(vf3, vf3_tds, s_nor, 0.5f);
	TransTri(vs3, vs3_tds, s_nor, 0.5f);
	TransLine(vf2, vf2_c, s_nor, 0.5f);
	TransLine(vs2, vs2_c, s_nor, 0.5f);

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	Gauss3Point(vf3_tds[0], vf3_tds[1], vf3_tds[2], vgf_tds);
	GaussLine3Point(vf2_c[0], vf2_c[1], vgf_c);

	Gauss3Point(vs3_tds[0], vs3_tds[1], vs3_tds[2], vgs_tds);
	Gauss3Point(vs3_up[0], vs3_up[1], vs3_up[2], vgs_up);
	GaussLine3Point(vs2_c[0], vs2_c[1], vgs_c);



	k_p = (1.0f / eps_p) - 1.0f;
	k_m = (1.0f / eps_m) - 1.0f;

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs;

			r = (vgf_tds[p] - vgs_tds[q]).Norm();
			r_up = (vgf_tds[p] - vgs_up[q]).Norm();
			r_c = (vgf_tds[p] - vgs_c[q]).Norm();
			r_c_up = (vgf_c[p] - vgs_up[q]).Norm();
			r_c_c = (vgf_c[p] - vgs_c[q]).Norm();

			G = exp(-J0 * k_*r) / r;
			G_up = exp(-J0 * k_*r_up) / r_up;
			G_c = exp(-J0 * k_*r_c) / r_c;
			G_c_up = exp(-J0 * k_*r_c_up) / r_c_up;
			G_c_c = exp(-J0 * k_*r_c_c) / r_c_c;

			Ze2 += w3[p] * w3[q] * (rou_f ^ (rou_s - (Thickness*s_nor)))*G;
			Ze3 += w3[p] * w3[q] * G_up;
			Ze4 += w3[p] * glw3[q] * G_c;

			if (bedge_f == 1)
			{
				Ze5 += glw3[p] * w3[q] * G_c_up;
				Ze6 += glw3[p] * glw3[q] * G_c_c;
			}
		}
	}

	Ze = (-k_m * Thickness*Ze2 / 4.0f) + ((k_m*Thickness) / (k_*k_))*(Ze5 - Ze3);

	return Ze;
}

Complex TDS::TTZmpKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &f_nor, VectorR3 &s_nor, VectorR3 *vf2, VectorR3 *vs2, value_t &eps_p, value_t &eps_m, int &bedge_f, int &bedge_s)
{
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0), Ze4(0, 0), Ze5(0, 0), Ze6(0, 0);
	Complex  G(0, 0), G_up(0, 0), G_c(0, 0), G_c_up(0, 0), G_c_c(0, 0);
	VectorR3 vf3_tds[3], vs3_tds[3], vs3_up[3], vf2_c[2], vs2_c[2], vf_tds, v;
	VectorR3 vgf[3], vgs[3], vgf_tds[3], vgs_tds[3], vgs_up[3], vgf_c[3], vgs_c[3];
	VectorR3 v2_m[2];
	VectorR3 rou_f, rou_s;
	value_t r, r_up, r_c, r_c_up, r_c_c;
	value_t k_p, k_m;

	TransTri(vf3, vs3_up, s_nor, 1.0f);
	TransTri(vf3, vf3_tds, s_nor, 0.5f);
	TransTri(vs3, vs3_tds, s_nor, 0.5f);
	TransLine(vf2, vf2_c, s_nor, 0.5f);
	TransLine(vs2, vs2_c, s_nor, 0.5f);

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	Gauss3Point(vf3_tds[0], vf3_tds[1], vf3_tds[2], vgf_tds);
	GaussLine3Point(vf2_c[0], vf2_c[1], vgf_c);

	Gauss3Point(vs3_tds[0], vs3_tds[1], vs3_tds[2], vgs_tds);
	Gauss3Point(vs3_up[0], vs3_up[1], vs3_up[2], vgs_up);
	GaussLine3Point(vs2_c[0], vs2_c[1], vgs_c);



	k_p = (1.0f / eps_p) - 1.0f;
	k_m = (1.0f / eps_m) - 1.0f;

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs;

			r = (vgf_tds[p] - vgs_tds[q]).Norm();
			r_up = (vgf_tds[p] - vgs_up[q]).Norm();
			r_c = (vgf_tds[p] - vgs_c[q]).Norm();
			r_c_up = (vgf_c[p] - vgs_up[q]).Norm();
			r_c_c = (vgf_c[p] - vgs_c[q]).Norm();

			G = exp(-J0 * k_*r) / r;
			G_up = exp(-J0 * k_*r_up) / r_up;
			G_c = exp(-J0 * k_*r_c) / r_c;
			G_c_up = exp(-J0 * k_*r_c_up) / r_c_up;
			G_c_c = exp(-J0 * k_*r_c_c) / r_c_c;

			Ze2 += w3[p] * w3[q] * (rou_f ^ (rou_s - (Thickness*s_nor)))*G;
			Ze3 += w3[p] * w3[q] * G_up;
			Ze4 += w3[p] * glw3[q] * G_c;

			if (bedge_f == 1)
			{
				Ze5 += glw3[p] * w3[q] * G_c_up;
				Ze6 += glw3[p] * glw3[q] * G_c_c;
			}
		}
	}

	Ze = (-k_p * Thickness*Ze2 / 4.0f) + ((k_p*Thickness) / (k_*k_))*(Ze5 - Ze3);

	return Ze;
}

Complex TDS::TTZmmKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &f_nor, VectorR3 &s_nor, VectorR3 *vf2, VectorR3 *vs2, value_t &eps_p, value_t &eps_m, int &bedge_f, int &bedge_s)
{
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0), Ze4(0, 0), Ze5(0, 0), Ze6(0, 0);
	Complex  G(0, 0), G_up(0, 0), G_c(0, 0), G_c_up(0, 0), G_c_c(0, 0);
	VectorR3 vf3_tds[3], vs3_tds[3], vs3_up[3], vf2_c[2], vs2_c[2], vf_tds, v;
	VectorR3 vgf[3], vgs[3], vgf_tds[3], vgs_tds[3], vgs_up[3], vgf_c[3], vgs_c[3];
	VectorR3 v2_m[2];
	VectorR3 rou_f, rou_s;
	value_t r, r_up, r_c, r_c_up, r_c_c;
	value_t k_p, k_m;

	TransTri(vf3, vs3_up, s_nor, 1.0f);
	TransTri(vf3, vf3_tds, s_nor, 0.5f);
	TransTri(vs3, vs3_tds, s_nor, 0.5f);
	TransLine(vf2, vf2_c, s_nor, 0.5f);
	TransLine(vs2, vs2_c, s_nor, 0.5f);

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	Gauss3Point(vf3_tds[0], vf3_tds[1], vf3_tds[2], vgf_tds);
	GaussLine3Point(vf2_c[0], vf2_c[1], vgf_c);

	Gauss3Point(vs3_tds[0], vs3_tds[1], vs3_tds[2], vgs_tds);
	Gauss3Point(vs3_up[0], vs3_up[1], vs3_up[2], vgs_up);
	GaussLine3Point(vs2_c[0], vs2_c[1], vgs_c);



	k_p = (1.0f / eps_p) - 1.0f;
	k_m = (1.0f / eps_m) - 1.0f;

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs;

			r = (vgf_tds[p] - vgs_tds[q]).Norm();
			r_up = (vgf_tds[p] - vgs_up[q]).Norm();
			r_c = (vgf_tds[p] - vgs_c[q]).Norm();
			r_c_up = (vgf_c[p] - vgs_up[q]).Norm();
			r_c_c = (vgf_c[p] - vgs_c[q]).Norm();

			G = exp(-J0 * k_*r) / r;
			G_up = exp(-J0 * k_*r_up) / r_up;
			G_c = exp(-J0 * k_*r_c) / r_c;
			G_c_up = exp(-J0 * k_*r_c_up) / r_c_up;
			G_c_c = exp(-J0 * k_*r_c_c) / r_c_c;

			Ze2 += w3[p] * w3[q] * (rou_f ^ (rou_s - (Thickness*s_nor)))*G;
			Ze3 += w3[p] * w3[q] * G_up;
			Ze4 += w3[p] * glw3[q] * G_c;

			if (bedge_f == 1)
			{
				Ze5 += glw3[p] * w3[q] * G_c_up;
				Ze6 += glw3[p] * glw3[q] * G_c_c;
			}
		}
	}

	Ze = (k_m * Thickness*Ze2 / 4.0f) + ((k_m*Thickness) / (k_*k_))*(Ze5 - Ze3) + ((k_m - k_p)*Thickness / (k_*k_))*(Ze6 - Ze4);

	return Ze;
}

Complex TDS::TTZppSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &f_nor, VectorR3 &s_nor, VectorR3 *vf2, VectorR3 *vs2, value_t &eps_p, value_t &eps_m, int &bedge_f, int &bedge_s)
{
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0), Ze4(0, 0), Ze5(0, 0), Ze6(0, 0);
	Complex  G(0, 0), G_up(0, 0), G_c(0, 0), G_c_up(0, 0), G_c_c(0, 0);
	VectorR3 vf3_tds[3], vs3_tds[3], vs3_up[3], vf2_c[2], vs2_c[2], vf_tds, v;
	VectorR3 vgf[3], vgs[3], vgf_tds[3], vgs_tds[3], vgs_up[3], vgf_c[3], vgs_c[3];
	VectorR3 v2_m[2];
	VectorR3 rou_f, rou_s;
	value_t r, r_up, r_c, r_c_up, r_c_c;
	value_t k_p, k_m;

	TransTri(vf3, vs3_up, s_nor, 1.0f);
	TransTri(vf3, vf3_tds, s_nor, 0.5f);
	TransTri(vs3, vs3_tds, s_nor, 0.5f);
	TransLine(vf2, vf2_c, s_nor, 0.5f);
	TransLine(vs2, vs2_c, s_nor, 0.5f);

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	Gauss3Point(vf3_tds[0], vf3_tds[1], vf3_tds[2], vgf_tds);
	GaussLine3Point(vf2_c[0], vf2_c[1], vgf_c);

	Gauss3Point(vs3_tds[0], vs3_tds[1], vs3_tds[2], vgs_tds);
	Gauss3Point(vs3_up[0], vs3_up[1], vs3_up[2], vgs_up);
	GaussLine3Point(vs2_c[0], vs2_c[1], vgs_c);



	k_p = (1.0f / eps_p) - 1.0f;
	k_m = (1.0f / eps_m) - 1.0f;

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			rou_s = vgs[q] - vs;

			r = (vgf_tds[p] - vgs_tds[q]).Norm();
			r_up = (vgf_tds[p] - vgs_up[q]).Norm();
			r_c = (vgf_tds[p] - vgs_c[q]).Norm();
			r_c_up = (vgf_c[p] - vgs_up[q]).Norm();
			r_c_c = (vgf_c[p] - vgs_c[q]).Norm();

			G = exp(-J0 * k_*r) / r;
			G_up = exp(-J0 * k_*r_up) / r_up;
			G_c = exp(-J0 * k_*r_c) / r_c;
			G_c_up = exp(-J0 * k_*r_c_up) / r_c_up;
			G_c_c = exp(-J0 * k_*r_c_c) / r_c_c;

			Ze2 += w3[p] * w3[q] * (rou_f ^ (rou_s - (Thickness*s_nor)))*G;
			Ze3 += w3[p] * w3[q] * G_up;
			Ze4 += w3[p] * glw3[q] * G_c;

			if (bedge_f == 1)
			{
				Ze5 += glw3[p] * w3[q] * G_c_up;
				Ze6 += glw3[p] * glw3[q] * G_c_c;
			}
		}
	}

	Ze = (k_p * Thickness*Ze2 / 4.0f) + ((k_p*Thickness) / (k_*k_))*(Ze5 - Ze3) + ((k_m - k_p)*Thickness / (k_*k_))*(Ze6 - Ze4);

	return Ze;
}

Complex TDS::TTZpmSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &f_nor, VectorR3 &s_nor, VectorR3 *vf2, VectorR3 *vs2, value_t &eps_p, value_t &eps_m, int &bedge_f, int &bedge_s)
{
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0), Ze4(0, 0), G(0, 0), G_m(0, 0), G_grad(0, 0);
	VectorR3 vgf[3], vgs[3], vgs_m[3], vgs_l[3], vgs_lm[3];
	VectorR3 v2_m[2];
	VectorR3 rou_f, rou_s;
	value_t r, r_m;
	value_t k_p, k_m;
	value_t area;

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);
	GaussLine3Point(v2[0], v2[1], vgs_l);

	TransTri(vgs, vgs_m, s_nor, 1.0f);
	TransTri(vgs_l, vgs_lm, s_nor, 0.5f);
	//TransLine(v2, v2_m, s_nor, 0.5);
	//GaussLine3Point(v2_m[0], v2_m[1], vgs_l);

	k_p = (1.0f / eps_p) - 1.0f;
	k_m = (1.0f / eps_m) - 1.0f;
	area = Area(vs3[0], vs3[1], vs3[2]);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)//非奇异部分
		{
			rou_s = vgs[q] - vs;
			r = (vgf[p] - vgs[q]).Norm();
			r_m = (vgf[p] - vgs_m[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G_m = value_t(-0.5) * k_ * k_ * r_m + (J0 * k_ / 6.0f) * (k_ * k_ * r_m * r_m - 6);
			//G = exp(-J0 * k_*r) / r;
			//G_m = exp(-J0 * k_*r_m) / r_m;
			//G_grad = (-J0 * k_ - (1.0f / r))*G;

			Ze2 += w3[p] * w3[q] * (rou_f^rou_s)*G;
			Ze3 += w3[p] * w3[q] * G_m;
		}
		//奇异部分
		Ze2 += w3[p] * (rou_f ^ (IspSingular(vs3, vgf[p]) + (vgf[p] - vs)*IsSingular(vs3, vgf[p]))) / area;
		Ze3 += w3[p] * IsSingular(vgs_m, vgf[p]) / area;

		for (int q = 0; q < 3; q++)
		{
			r = (vgf[p] - vgs_lm[q]).Norm();
			G_grad = (-J0 * k_ - (1.0f / r))*(exp(-J0 * k_*r) / r);
			Ze4 += w3[p] * glw3[q] * G_grad;
		}
		if (bedge == 1)
		{

		}
	}

	Ze = (-k_p * Ze2 / 4.0f) + (Thickness*(k_m - k_p) / (k_*k_))*(Ze3 + Ze4 / 2.0f);

	return Ze;
}

Complex TDS::TTZmpSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &f_nor, VectorR3 &s_nor, VectorR3 *vf2, VectorR3 *vs2, value_t &eps_p, value_t &eps_m, int &bedge_f, int &bedge_s)
{
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0), Ze4(0, 0), G(0, 0), G_m(0, 0), G_grad(0, 0);
	VectorR3 vgf[3], vgs[3], vgs_m[3], vgs_l[3], vgs_lm[3];
	VectorR3 v2_m[2];
	VectorR3 rou_f, rou_s;
	value_t r, r_m;
	value_t k_p, k_m;
	value_t area;

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);
	GaussLine3Point(v2[0], v2[1], vgs_l);

	TransTri(vgs, vgs_m, s_nor, 1.0f);
	TransTri(vgs_l, vgs_lm, s_nor, 0.5f);
	//TransLine(v2, v2_m, s_nor, 0.5);
	//GaussLine3Point(v2_m[0], v2_m[1], vgs_l);

	k_p = (1.0f / eps_p) - 1.0f;
	k_m = (1.0f / eps_m) - 1.0f;
	area = Area(vs3[0], vs3[1], vs3[2]);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)//非奇异部分
		{
			rou_s = vgs[q] - vs;
			r = (vgf[p] - vgs[q]).Norm();
			r_m = (vgf[p] - vgs_m[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G_m = value_t(-0.5) * k_ * k_ * r_m + (J0 * k_ / 6.0f) * (k_ * k_ * r_m * r_m - 6);
			//G = exp(-J0 * k_*r) / r;
			//G_m = exp(-J0 * k_*r_m) / r_m;
			//G_grad = (-J0 * k_ - (1.0f / r))*G;

			Ze2 += w3[p] * w3[q] * (rou_f^rou_s)*G;
			Ze3 += w3[p] * w3[q] * G_m;
		}
		//奇异部分
		Ze2 += w3[p] * (rou_f ^ (IspSingular(vs3, vgf[p]) + (vgf[p] - vs)*IsSingular(vs3, vgf[p]))) / area;
		Ze3 += w3[p] * IsSingular(vgs_m, vgf[p]) / area;

		for (int q = 0; q < 3; q++)
		{
			r = (vgf[p] - vgs_lm[q]).Norm();
			G_grad = (-J0 * k_ - (1.0f / r))*(exp(-J0 * k_*r) / r);
			Ze4 += w3[p] * glw3[q] * G_grad;
		}
		if (bedge == 1)
		{

		}
	}

	Ze = (-k_p * Ze2 / 4.0f) + (Thickness*(k_m - k_p) / (k_*k_))*(Ze3 + Ze4 / 2.0f);

	return Ze;
}

Complex TDS::TTZmmSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &f_nor, VectorR3 &s_nor, VectorR3 *vf2, VectorR3 *vs2, value_t &eps_p, value_t &eps_m, int &bedge_f, int &bedge_s)
{
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0), Ze4(0, 0), G(0, 0), G_m(0, 0), G_grad(0, 0);
	VectorR3 vgf[3], vgs[3], vgs_m[3], vgs_l[3], vgs_lm[3];
	VectorR3 v2_m[2];
	VectorR3 rou_f, rou_s;
	value_t r, r_m;
	value_t k_p, k_m;
	value_t area;

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);
	GaussLine3Point(v2[0], v2[1], vgs_l);

	TransTri(vgs, vgs_m, s_nor, 1.0f);
	TransTri(vgs_l, vgs_lm, s_nor, 0.5f);
	//TransLine(v2, v2_m, s_nor, 0.5);
	//GaussLine3Point(v2_m[0], v2_m[1], vgs_l);

	k_p = (1.0f / eps_p) - 1.0f;
	k_m = (1.0f / eps_m) - 1.0f;
	area = Area(vs3[0], vs3[1], vs3[2]);

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)//非奇异部分
		{
			rou_s = vgs[q] - vs;
			r = (vgf[p] - vgs[q]).Norm();
			r_m = (vgf[p] - vgs_m[q]).Norm();
			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G_m = value_t(-0.5) * k_ * k_ * r_m + (J0 * k_ / 6.0f) * (k_ * k_ * r_m * r_m - 6);
			//G = exp(-J0 * k_*r) / r;
			//G_m = exp(-J0 * k_*r_m) / r_m;
			//G_grad = (-J0 * k_ - (1.0f / r))*G;

			Ze2 += w3[p] * w3[q] * (rou_f^rou_s)*G;
			Ze3 += w3[p] * w3[q] * G_m;
		}
		//奇异部分
		Ze2 += w3[p] * (rou_f ^ (IspSingular(vs3, vgf[p]) + (vgf[p] - vs)*IsSingular(vs3, vgf[p]))) / area;
		Ze3 += w3[p] * IsSingular(vgs_m, vgf[p]) / area;

		for (int q = 0; q < 3; q++)
		{
			r = (vgf[p] - vgs_lm[q]).Norm();
			G_grad = (-J0 * k_ - (1.0f / r))*(exp(-J0 * k_*r) / r);
			Ze4 += w3[p] * glw3[q] * G_grad;
		}
		if (bedge == 1)
		{

		}
	}

	Ze = (k_p * Ze2 / 4.0f) + (Thickness*(k_m - k_p) / (k_*k_))*(Ze3 + Ze4 / 2.0f);

	return Ze;
}

Complex TDS::TNZppKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &s_nor, value_t &eps)
{
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0);
	VectorR3 vgf[3], vgs[3], vgs_m[3];
	VectorR3 rou_f;
	value_t r, r_m;
	value_t k_p;
	Complex G(0, 0), G_m(0, 0);

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	TransTri(vgs, vgs_m, s_nor, 1.0f);
	k_p = (1.0f / eps) - 1.0f;

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			r = (vgf[p] - vgs[q]).Norm();
			r_m = (vgf[p] - vgs_m[q]).Norm();

			G = exp(-J0 * k_*r) / r;
			G_m = exp(-J0 * k_*r_m) / r_m;

			Ze1 += (rou_f^s_nor)*G;
			Ze2 += G;
			Ze3 += G_m;
		}
	}

	Ze = (Thickness * k_p*Ze1/2.0f) + k_p * (Ze2 - Ze3);
}

Complex TDS::TNZpmKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &s_nor, value_t &eps)
{
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0);
	VectorR3 vgf[3], vgs[3], vgs_m[3];
	VectorR3 rou_f;
	value_t r, r_m;
	value_t k_p;
	Complex G(0, 0), G_m(0, 0);

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	TransTri(vgs, vgs_m, s_nor, 1.0f);
	k_p = (1.0f / eps) - 1.0f;

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			r = (vgf[p] - vgs[q]).Norm();
			r_m = (vgf[p] - vgs_m[q]).Norm();

			G = exp(-J0 * k_*r) / r;
			G_m = exp(-J0 * k_*r_m) / r_m;

			Ze1 += (rou_f^s_nor)*G;
			Ze2 += G;
			Ze3 += G_m;
		}
	}

	Ze = (Thickness * k_p*Ze1 / 2.0f) + k_p * (Ze2 - Ze3);
}

Complex TDS::TNZmpKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &s_nor, value_t &eps)
{
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0);
	VectorR3 vgf[3], vgs[3], vgs_m[3];
	VectorR3 rou_f;
	value_t r, r_m;
	value_t k_p;
	Complex G(0, 0), G_m(0, 0);

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	TransTri(vgs, vgs_m, s_nor, 1.0f);
	k_p = (1.0f / eps) - 1.0f;

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			r = (vgf[p] - vgs[q]).Norm();
			r_m = (vgf[p] - vgs_m[q]).Norm();

			G = exp(-J0 * k_*r) / r;
			G_m = exp(-J0 * k_*r_m) / r_m;

			Ze1 += (rou_f^s_nor)*G;
			Ze2 += G;
			Ze3 += G_m;
		}
	}

	Ze = (-Thickness * k_p*Ze1 / 2.0f) + k_p * (Ze2 - Ze3);
}

Complex TDS::TNZmmKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &s_nor, value_t &eps)
{
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0);
	VectorR3 vgf[3], vgs[3], vgs_m[3];
	VectorR3 rou_f;
	value_t r, r_m;
	value_t k_p;
	Complex G(0, 0), G_m(0, 0);

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	TransTri(vgs, vgs_m, s_nor, 1.0f);
	k_p = (1.0f / eps) - 1.0f;

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			r = (vgf[p] - vgs[q]).Norm();
			r_m = (vgf[p] - vgs_m[q]).Norm();

			G = exp(-J0 * k_*r) / r;
			G_m = exp(-J0 * k_*r_m) / r_m;

			Ze1 += (rou_f^s_nor)*G;
			Ze2 += G;
			Ze3 += G_m;
		}
	}

	Ze = (Thickness * k_p*Ze1 / 2.0f) + k_p * (Ze2 - Ze3);
}

Complex TDS::TNZppSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &s_nor, value_t &eps)
{
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0);
	VectorR3 vgf[3], vgs[3], vgs_m[3];
	VectorR3 rou_f;
	value_t r, r_m;
	value_t k_p;
	Complex G(0, 0), G_m(0, 0);

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	TransTri(vgs, vgs_m, s_nor, 1.0f);
	k_p = (1.0f / eps) - 1.0f;

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			r = (vgf[p] - vgs[q]).Norm();
			r_m = (vgf[p] - vgs_m[q]).Norm();

			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G_m = exp(-J0 * k_*r_m) / r_m;

			Ze1 += (rou_f^s_nor)*G;
			Ze2 += G;
			Ze3 += G_m;
		}
		Ze1 += w3[p] * (rou_f^s_nor)*IsSingular(vgs, vgf[p]);
		Ze2 += w3[p] * IsSingular(vgs, vgf[p]);
	}

	Ze = (Thickness * k_p*Ze1 / 2.0f) + k_p * (Ze2 - Ze3);
}

Complex TDS::TNZpmSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &s_nor, value_t &eps)
{

}

Complex TDS::TNZmpSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &s_nor, value_t &eps)
{
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Ze3(0, 0);
	VectorR3 vgf[3], vgs[3], vgs_m[3];
	VectorR3 rou_f;
	value_t r, r_m;
	value_t k_p;
	Complex G(0, 0), G_m(0, 0);

	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	TransTri(vgs, vgs_m, s_nor, 1.0f);
	k_p = (1.0f / eps) - 1.0f;

	for (int p = 0; p < 3; p++)
	{
		rou_f = vgf[p] - vf;
		for (int q = 0; q < 3; q++)
		{
			r = (vgf[p] - vgs[q]).Norm();
			r_m = (vgf[p] - vgs_m[q]).Norm();

			G = value_t(-0.5) * k_ * k_ * r + (J0 * k_ / 6.0f) * (k_ * k_ * r * r - 6);
			G_m = exp(-J0 * k_*r_m) / r_m;

			Ze1 += (rou_f^s_nor)*G;
			Ze2 += G;
			Ze3 += G_m;
		}
		Ze1 += w3[p] * (rou_f^s_nor)*IsSingular(vgs, vgf[p]);
		Ze2 += w3[p] * IsSingular(vgs, vgf[p]);
	}

	Ze = (-Thickness * k_p*Ze1 / 2.0f) + k_p * (Ze2 - Ze3);
}

Complex TDS::TNZmmSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &s_nor, value_t &eps)
{

}

Complex TDS::eVKernel(VectorR3 *vp3, VectorR3 *vm3, VectorR3 &vp, VectorR3 &vm)
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

void TDS::fillZ()
{
	fillZTT();
	fillZTN();
	fillZNT();
	fillZNN();
}

void TDS::fillZTT()
{
	int vp, vm, f_fld_plu, f_fld_min, f_src_plu, f_src_min, bedge;
	int v1, v2;
	value_t l_fld, l_src, eps_p, eps_m;
	VectorR3 v_fld_plu, v_fld_plu3[3], v_fld_min, v_fld_min3[3];
	VectorR3 v_src_plu, v_src_plu3[3], v_src_min, v_src_min3[3];
	VectorR3 v_src_ce2[2];
	VectorR3 v_fld_pluNor, v_fld_minNor, v_src_pluNor, v_src_minNor;
	Triangle tri_plu, tri_min;

	Complex coff;//待求解
	coff = k_ * Z0*omiga / PI4;
	tool::BarAndPercent bar;
	for (int f = 0; f < unknowns_t; f++)
	{
		bar(f + 1, unknowns_t);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		v_fld_pluNor = tri_plu.getNormal();

		if (vm != -1)
		{
			v_fld_min = mesh_ptr_->getVertex(vp);//场点+
			tri_min = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
			tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);
			v_fld_minNor = tri_min.getNormal();
		}

		for (int s = 0; s < unknowns_t; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src, eps_p, eps_m, bedge);

			v_src_plu = mesh_ptr_->getVertex(vp);//场点+
			tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);//场三角面+
			tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
			v_src_pluNor = tri_plu.getNormal();

			if (vm != -1)
			{
				v_src_min = mesh_ptr_->getVertex(vp);//场点+
				tri_min = mesh_ptr_->getTriangleRef(f_src_plu);//场三角面+
				tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);
				v_src_minNor = tri_min.getNormal();
			}
			else
			{
				ce_ptr_->getCommonEdge(s, v1, v2);
				v_src_ce2[0] = mesh_ptr_->getVertex(v1);
				v_src_ce2[1] = mesh_ptr_->getVertex(v2);
			}

			if (f_fld_plu == f_src_plu)
			{
				Zpp = TTZppSingular(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu, v_fld_pluNor, v_src_pluNor,v_src_ce2, eps_p, eps_m, bedge);
			}
			else
			{
				Zpp = TTZppKernel(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu, v_fld_pluNor, v_src_pluNor, v_src_ce2, eps_p, eps_m, bedge);
			}

			if(f_fld_min!=-1)
				if (f_fld_min == f_src_plu)
				{
					Zmp = TTZmpSingular(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu, v_fld_minNor, v_src_pluNor, v_src_ce2, eps_p, eps_m, bedge);
				}
				else
				{
					Zmp = TTZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu, v_fld_minNor, v_src_pluNor, v_src_ce2, eps_p, eps_m, bedge);
				}
			
			if(f_src_min!=-1)
				if (f_fld_plu == f_src_min)
				{
					Zpm = TTZpmSingular(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min, v_fld_pluNor, v_src_minNor, v_src_ce2, eps_p, eps_m, bedge);
				}
				else
				{
					Zpm = TTZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min, v_fld_pluNor, v_src_minNor, v_src_ce2, eps_p, eps_m, bedge);
				}

			if (f_fld_min != -1 && f_src_min != -1)
				if (f_fld_min == f_src_min)
				{
					Zmm = TTZmmSingular(v_fld_min3, v_src_min3, v_fld_min, v_src_min, v_fld_minNor, v_src_minNor, v_src_ce2, eps_p, eps_m, bedge);
				}
				else
				{
					Zmm = TTZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min, v_fld_minNor, v_src_minNor, v_src_ce2, eps_p, eps_m, bedge);
				}

			Z(f, s) = coff * l_fld*l_src* (Zpp + Zpm + Zmp + Zmm); //需要修改
		}
	}
}

void TDS::fillZTN()
{
	int vp, vm, f_fld_plu, f_fld_min, f_src_n, bedge;
	int v1, v2;
	value_t l_fld, l_src, eps;
	VectorR3 v_fld_plu, v_fld_plu3[3], v_fld_min, v_fld_min3[3];
	VectorR3 v_src_plu, v_src_n3[3], v_src_min, v_src_min3[3];
	VectorR3 v_src_ce2[2];
	VectorR3 v_fld_pluNor, v_fld_minNor, v_src_Nor;
	Triangle tri_plu, tri_min, tri_n;

	Complex coff;//待求解
	coff = k_ * Z0*omiga / PI4;
	tool::BarAndPercent bar;
	for (int f = 0; f < unknowns_t; f++)
	{
		bar(f + 1, unknowns_t);
		ce_ptr_->getCommonEdge(f, vp, vm, f_fld_plu, f_fld_min, l_fld);

		v_fld_plu = mesh_ptr_->getVertex(vp);//场点+
		tri_plu = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
		tri_plu.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
		v_fld_pluNor = tri_plu.getNormal();

		if (vm != -1)
		{
			v_fld_min = mesh_ptr_->getVertex(vp);//场点+
			tri_min = mesh_ptr_->getTriangleRef(f_fld_plu);//场三角面+
			tri_min.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);
			v_fld_minNor = tri_min.getNormal();
		}
		for (int s = 0; s < unknowns_n; s++)
		{
			Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
			ce_ptr_->getPluseTri(s, f_src_n, eps);
			tri_n = mesh_ptr_->getTriangleRef(f_src_n);
			tri_n.getVertex(v_src_n3[0], v_src_n3[1], v_src_n3[2]);
			v_src_Nor = tri_n.getNormal();
			if (f_fld_plu == f_src_n)
			{
				Zpp = TNZppSingular(v_fld_plu3, v_src_n3, v_fld_plu, v_src_Nor, eps);
			}
			else
			{
				Zpp = TNZppKernel(v_fld_plu3, v_src_n3, v_fld_plu, v_src_Nor, eps);
			}

			if (f_fld_min != -1)
			{
				if (f_fld_min == f_src_n)
				{
					Zmp = TNZmpSingular(v_fld_min3, v_src_n3, v_fld_min, v_src_Nor, eps);
				}
				else
				{
					Zmp = TNZmpKernel(v_fld_min3, v_src_n3, v_fld_min, v_src_Nor, eps);
				}
			}
			Z(f, s + unknowns_t) = coff * l_fld*(Zpp + Zmp);
		}
	}
}

void TDS::fillZNT()
{

}

void TDS::fillZNN()
{

}

void TDS::fillV()
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

bool TDS::readExcEdges(const Qstring & rad_file, Qstring& stateInfo)
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

void TDS::readExcEdges(const Qstring & rad_file)
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

void TDS::radiateV()
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

value_t TDS::getBiRCS(const VectorR3 & sca_k) const
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

void TDS::getEFiled(component::FieldData * pdata, const VectorR3 & rad_k, const VectorR3 & rad_ev, const VectorR3 & rad_eh) const
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

void TDS::TransTri(VectorR3 *v3_ori, VectorR3 *v3_end,VectorR3& Nor,value_t wgt)
{
	v3_end[0] = v3_ori[0] + wgt * Thickness * Nor;
	v3_end[1] = v3_ori[1] + wgt * Thickness * Nor;
	v3_end[2] = v3_ori[2] + wgt * Thickness * Nor;
}

void TDS::TransLine(VectorR3 *v2_ori, VectorR3 *v2_end, VectorR3& Nor, value_t wgt)
{
	v2_end[0] = v2_ori[0] + wgt * Thickness*Nor;
	v2_end[1] = v2_ori[1] + wgt * Thickness*Nor;
}

void TDS::TransPoi(VectorR3 &v_ori, VectorR3 &v_end, VectorR3 &Nor, value_t wgt)
{
	v_end = v_ori + wgt * Thickness*Nor;
}

value_t TDS::IsSingular(VectorR3 * vs3, VectorR3 & vf)
{
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[1]);
	VectorR3 ni = n.Normalize();
	VectorR3 rou = vf - ni * (ni^vf);

	value_t Is = 0.0f;
	value_t d = (vf - vs3[0]) ^ ni;
	for (int i = 0; i < 3; i++)
	{
		VectorR3 li = (vs3[(i + 1) % 3] - vs3[i]).Normalize();//not normalize
		VectorR3 r_f = vf - vs3[i];
		VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
		VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
		VectorR3 ui = li * ni;
		value_t p0 = abs((vs3[i] - vf) ^ ui);
		value_t lp = (vs3[(i + 1) % 3] - vf) ^ li;
		value_t lm = (vs3[i] - vf) ^ li;
		value_t lp1 = (vf - vs3[(i + 1) % 3]) ^ li;
		value_t lm1 = (vf - vs3[i]) ^ li;
		value_t rp = sqrt(p0*p0 + d * d + lp * lp);
		value_t rm = sqrt(p0*p0 + d * d + lm * lm);

		if (p0 < 1e-8)
			continue;

		VectorR3 pi = ((rou_p - rou) - lp * li).Normalize();
		value_t r02 = p0 * p0 + d * d;
		Is += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
	}
	return Is;
}

VectorR3 TDS::IspSingular(VectorR3 * vs3, VectorR3 & vf)
{
	VectorR3 n = (vs3[1] - vs3[0])*(vs3[2] - vs3[1]);
	VectorR3 ni = n.Normalize();
	VectorR3 rou = vf - ni * (ni^vf);

	VectorR3 Isp(0, 0, 0);
	value_t d = (vf - vs3[0]) ^ ni;
	for (int i = 0; i < 3; i++)
	{
		VectorR3 li = (vs3[(i + 1) % 3] - vs3[i]).Normalize();//not normalize
		VectorR3 r_f = vf - vs3[i];
		VectorR3 rou_p = vs3[(i + 1) % 3] - (ni ^ vs3[(i + 1) % 3]) * ni;
		VectorR3 rou_m = vs3[i] - (ni ^ vs3[i]) * ni;
		VectorR3 ui = li * ni;
		value_t p0 = abs((vs3[i] - vf) ^ ui);
		value_t lp = (vs3[(i + 1) % 3] - vf) ^ li;
		value_t lm = (vs3[i] - vf) ^ li;
		value_t lp1 = (vf - vs3[(i + 1) % 3]) ^ li;
		value_t lm1 = (vf - vs3[i]) ^ li;
		value_t rp = sqrt(p0*p0 + d * d + lp * lp);
		value_t rm = sqrt(p0*p0 + d * d + lm * lm);

		if (p0 < 1e-8)
			continue;

		VectorR3 pi = ((rou_p - rou) - lp * li).Normalize();
		value_t r02 = p0 * p0 + d * d;
		//Isp += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
		Isp += (r02*(log((rp + lp) / (rm + lm))) + rp * lp - rm * lm)*ui;
	}
	return 0.5f*Isp;
}

bool TDS::writeZIVData()
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

void TDS::calculateSurfaceCurrent(std::vector<component::CurrentData>* currents) const
{
	Qcx_vec cur_vec = I;
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
