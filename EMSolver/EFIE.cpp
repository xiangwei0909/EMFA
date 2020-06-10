#include "EFIE.h"
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


EFIE::EFIE()
{

}

EFIE::~EFIE()
{
}

void EFIE::init(component::ConfigLoader * ploader)
{
    EM::init(ploader);
    if(rtype_ == policy::ResultType::Rad)
        readExcEdges(dir_ + '/' + folder_name_ + ".rad");

    dir_ = tool::creatFolder(dir_, folder_name_ + "_EFIE");

    auto log_name = dir_ + (rtype_ == policy::ResultType::Rad ? "/rad_" : "/sca_") + "runtime.log";
    logger_.open(log_name);

    mesh_ptr_ = std::make_shared<Mesh>();
    ce_ptr_ = std::make_shared<CommonEdge>();

    mesh_ptr_->loadMeshFile(ploader->getMeshPath());
    ce_ptr_->buildCommonEdge(mesh_ptr_);
	Sigma = ploader->getSigma();
	Dx = ploader->getDx();
	Dy = ploader->getDy();
	Array_x = ploader->getArray_x();
	Array_y = ploader->getArray_y();
    unknowns_ = ce_ptr_->getCommonEdgeNum();
    k_ = PI2 * incidence_.freq / cc;
	threshold_edm = 3.0f*cc / incidence_.freq;
	logger_ << ploader->getReport();
	reportInfo(logger_);
    TIME_LOG("init");
}

void EFIE::solve()
{
    SEGMENT("Solve");
    Z.zeros(unknowns_, unknowns_);
	Z_eigen.setZero(unknowns_, unknowns_);
    LOG(fillZ(), "Filling impedance matrix");
    TIME_LOG("fillZ");

    I.zeros(unknowns_);
    V.zeros(unknowns_);
	I_eigen.setZero(unknowns_);
	V_eigen.setZero(unknowns_);

    if(rtype_ == policy::ResultType::Rad)
    {
        LOG(radiateV(), "Filling rad voltage vector");
    }
    else
    {
        LOG(fillV(), "Filling sca voltage vector");
    }
    TIME_LOG("fillV");

    Qcout << "Iterative solve:" << std::endl;
    int iterNum = 0;
    if(!iml::GMRES(Z_eigen, V_eigen, I_eigen, 0.001, 30, 5))
    {
      Qcerr << "fatal error: fail solving Z*I = V" << std::endl;
      return;
    }
    Qcout << "  -> Iterative number: " << iterNum << std::endl;

    //Qcout << setw(30) << "Solving matrix equation:";
    //if (!arma::solve(I, Z, V, arma::solve_opts::fast))
       //throw std::runtime_error("fail solving matrix equation");
    //I_eigen = Z_eigen.partialPivLu().solve(V_eigen);
	if (writeZIVData())
		Qcout << "The data is already saved!" << std::endl;
    Qcout << "success" << std::endl;
    TIME_LOG("solve");

#ifdef _DEBUG
    //if (!writeZIVData())
    //  Qcerr << "fail writing ZIV data" << std::endl;
#endif
}

void EFIE::output()
{
    SEGMENT("Result");
    Result result;
    if(rtype_ == policy::ResultType::Rad)
    {
        LOG(result.getRadiationField(this, &EFIE::getEFiled), "Calculating Radiation");
        TIME_LOG("getRadiation");
		LOG(getCurrentOnFeed(), "Calculating the Current on the feed");
    }
    else
    {
        LOG(result.getBistaticRCS(this, &EFIE::getBiRCS), "Calculating BistaticRCS");
        TIME_LOG("getBistaticRCS");
    }

	//LOG(result.getNearField(this, &EFIE::getNearEField), "Calculating Near field");
	LOG(result.getCurrentDistribution(this, &EFIE::calculateSurfaceCurrent), "Calculating surface current");
}

void EFIE::clear()
{
    Z.reset();
    I.reset();
    V.reset();
    exc_edge_.clear();

    mesh_ptr_->clear();
    ce_ptr_->clear();
}


void EFIE::reportInfo(Qostream& strm) const
{
    mesh_ptr_->reportInfo(strm);
    ce_ptr_->reportInfo(strm);
}

Complex EFIE::eZppKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex EFIE::eZpmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex EFIE::eZmpKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex EFIE::eZmmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex EFIE::eZppSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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
		VectorR3 n = (vs3[1]-vs3[0])*(vs3[2]-vs3[0]);
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
			//VectorR3 n = (li*lj);
			VectorR3 u = li*n;
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
			value_t r02 = p0*p0 + d*d;
			Is1 += (0.5f*ui)*(r02*log((rp + lp) / (rm + lm)) + (lp*rp) - (lm*rm));
			Is2 += (pi^ui)*((p0*log((rp + lp) / (rm + lm))) - abs(d)*(atan((p0*lp) / (r02 + abs(d)*rp)) - atan((p0*lm) / (r02 + abs(d)*rm))));
		}
        Ze2 += w3[p] * (0.25f * (rou_f ^ (Is1+Is2*(vgf[p]-vs))) - (Is2 / (k_ * k_)));
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

Complex EFIE::eZpmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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
		Ze2 += w3[p] * (0.25f * (-1.0f*rou_f ^ (Is1 + Is2*(vgf[p]-vs))) + (Is2 / (k_ * k_)));
    }

	if (abs(Sigma) > 1e-8)
	{
		for (int i = 0; i < 3; i++)
		{
			rou_f = vgf[i] - vf;
			rou_s = vgf[i] - vs;
			Ze3 += w3[i] * rou_f^rou_s;
		}
		Ze3 = (PI *Sigma/ (S_s*J0*k_*Z0))*Ze3;
	}

    Ze = (Ze1 + (1.0f / S_s) * Ze2)+Ze3;

    return Ze;
}

Complex EFIE::eZmpSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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
		Ze2 += w3[p] * (0.25f * (rou_f ^ (Is1 + Is2*(vgf[p]-vs))) + (Is2 / (k_ * k_)));
    }


	if (abs(Sigma) > 1e-8)
	{
		for (int i = 0; i < 3; i++)
		{
			rou_f = vgf[i] - vf;
			rou_s = vgf[i] - vs;
			Ze3 += w3[i] * rou_f^rou_s;
		}
		Ze3 = (PI*Sigma/ (S_s*J0*k_*Z0))*Ze3;
	}

    Ze =(Ze1 + (1.0f / S_s) * Ze2) + Ze3;

    return Ze;
}

Complex EFIE::eZmmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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
			Is2 += (pi^ui)*((p0*log((rp + lp) / (rm + lm)))-abs(d)*(atan((p0*lp)/(r02+abs(d)*rp))-atan((p0*lm)/(r02+abs(d)*rm))));
		}
		Ze2 += w3[p] * (0.25f * (-1.0f*rou_f ^ (Is1 + Is2*(vgf[p]-vs))) - (Is2 / (k_ * k_)));
    }

	if (abs(Sigma) > 1e-8)
	{
		for (int i = 0; i < 3; i++)
		{
			rou_f = vgf[i] - vf;
			rou_s = vgf[i] - vs;
			Ze3 += w3[i] * rou_f^rou_s;
		}
		Ze3 = (PI*Sigma/ (S_s*J0*k_*Z0))*Ze3;
	}

    Ze = (Ze1 + (1.0f / S_s) * Ze2) + Ze3;

    return Ze;
}

Complex EFIE::EDMKernel(VectorR3 *vf_p3, VectorR3 *vf_m3, VectorR3 *vs_p3, VectorR3 *vs_m3, value_t &l_fld, value_t &l_src)
{
	VectorR3 vf_cen_plu = center(vf_p3);
	VectorR3 vf_cen_min = center(vf_m3);

	VectorR3 vs_cen_plu = center(vs_p3);
	VectorR3 vs_cen_min = center(vs_m3);

	VectorR3 r_fld = (vf_cen_plu + vf_cen_min) / 2.0f;
	VectorR3 r_src = (vs_cen_plu + vs_cen_min) / 2.0f;

	VectorR3 R = r_fld - r_src;
	value_t R_norm = R.Norm();

	VectorR3 m_s = l_src * (vs_cen_min - vs_cen_plu);
	Complex C = (1.0f / (R_norm*R_norm))*(1.0f + 1.0f / (J0*k_*R_norm));
	VectorR3 M_s = (R^m_s)*R / (R_norm*R_norm);

	VectorC3 E;
	E = ((M_s - m_s)*((J0*k_ / R_norm) + C) + 2.0f*M_s*C)*exp(-J0 * k_*R_norm)*Z0 / PI4;

	return (-l_fld *(E ^ (vf_cen_min - vf_cen_plu)));
}

Complex EFIE::eVKernel(VectorR3 *vp3, VectorR3 *vm3, VectorR3 &vp, VectorR3 &vm)
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

void EFIE::fillZ()
{
    int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle
    VectorR3 v_fld_plu, v_fld_min, v_fld_plu3[3], v_fld_min3[3];
    VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
	VectorR3 v_fld_cen, v_src_cen;
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

		v_fld_cen = (center(v_fld_plu3)+center(v_fld_min3))/2.0f;

        for (int s = 0; s < unknowns_; s++)
        {
            ce_ptr_->getCommonEdge(s, vp, vm, f_src_plu, f_src_min, l_src);

            v_src_plu = mesh_ptr_->getVertex(vp);
            v_src_min = mesh_ptr_->getVertex(vm);

            tri_plu = mesh_ptr_->getTriangleRef(f_src_plu);
            tri_min = mesh_ptr_->getTriangleRef(f_src_min);

            tri_plu.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
            tri_min.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);
			
			v_src_cen = (center(v_src_plu3) + center(v_src_min3)) / 2.0f;
			if ((v_fld_cen - v_src_cen).Norm() > threshold_edm)
			{
				Z(f, s) = EDMKernel(v_fld_plu3, v_fld_min3, v_src_plu3, v_src_min3, l_fld, l_src);
			}
			else
			{
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
					Zpm = -eZppKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
				}
				//field face- <--> source face+
				if (f_fld_min == f_src_plu)
				{
					Zmp = eZmpSingular(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
				}
				else
				{
					Zmp = -eZppKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
				}
				//field face- <--> source face-
				if (f_fld_min == f_src_min)
				{
					Zmm = eZmmSingular(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
				}
				else
				{
					Zmm = eZppKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
				}
				Z(f, s) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
				Z_eigen(f, s) = coef * l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
			}
            
        }
    }
}

void EFIE::fillV()
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
		V_eigen(u) = 0.5f * ln * vk;
    }
}

bool EFIE::readExcEdges(const Qstring & rad_file, Qstring& stateInfo)
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

void EFIE::readExcEdges(const Qstring & rad_file)
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

void EFIE::radiateV()
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

value_t EFIE::getBiRCS(const VectorR3 & sca_k) const
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
        Es = Es + (Esp + Esm) * (I_eigen[u] * L);
    }

    auto es = 0.5f * sca_k * (sca_k * Es);
    value_t rcs = Z0 * Z0 * k_ * k_ * es.norm() / PI4;

	//return rcs;
    return 10 * log10(rcs);
}

void EFIE::getEFiled(component::FieldData * pdata, const VectorR3 & rad_k, const VectorR3 & rad_ev, const VectorR3 & rad_eh) const
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
		//////////////////////////
		/*int nv1, nv2;
		VectorR3 v1, v2;
		VectorR3 ce_cen;
		ce_ptr_->getCommonEdge(u, nv1, nv2);
		v1 = mesh_ptr_->getVertex(nv1);
		v2 = mesh_ptr_->getVertex(nv2);
		ce_cen = (v1 + v2) / 2.0f;
		if (ce_cen.x > 0.225f || ce_cen.y > 0.252f)
			continue;*/
		//////////////////////////
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

void EFIE::getCurrentOnFeed()
{
	int nv1, nv2;
	CVector ifed;
	ifed.setZero(Array_x*Array_y);
	Qofstream outputI(dir_ + "/I_feed.txt");
	for (int u = 0; u < unknowns_; u++)
	{
		ce_ptr_->getCommonEdge(u, nv1, nv2);
		for (const auto& elem : exc_edge_)
		{
			if ((elem.first == nv1 && elem.second == nv2) || (elem.first == nv2 && elem.second == nv1))
			{
				Complex i_temp;
				VectorR3 va, vb,vmid;
				int xid, yid;
				va = mesh_ptr_->getVertex(nv1);
				vb = mesh_ptr_->getVertex(nv2);
				vmid = (va + vb) / 2.0f;

				auto length = ce_ptr_->getCommonEdgeLength(u);
				if (elem.first == nv1)
					i_temp = I(u)*length;
				else
					i_temp = I(u)*length;

				xid = round(vmid.x / Dx);
				yid = round(vmid.y / Dy);

				ifed(xid*Array_y + yid) = i_temp;
				//outputI << std::setw(18) << i_temp.real() << std::setw(18) << i_temp.imag() << '\n';
			}
		}
	}
	for (int i = 0; i < Array_x*Array_y; i++)
	{
		outputI << std::setw(18) << ifed(i).real() << std::setw(18) << ifed(i).imag() << '\n';
	}
	outputI.flush();
	outputI.close();
}

bool EFIE::writeZIVData()
{
    Qofstream outputZ(dir_ + "/matrix_Z.txt");
    if (outputZ.fail())
        return false;
    Z.save(outputZ,arma::arma_ascii);

    Qofstream outputV(dir_ + "/matrix_V.txt");
    if (outputV.fail())
        return false;
    V.save(outputV,arma::arma_ascii);

    Qofstream outputI(dir_ + "/martix_I.txt");
    if (outputI.fail())
        return false;
    I.save(outputI,arma::arma_ascii);

    return true;
}

void EFIE::calculateSurfaceCurrent(std::vector<component::CurrentData>* currents) const
{
	/*Qcx_vec cur_vec = I;
	const size_t tri_num = mesh_ptr_->getTriangleNum();
	const size_t node_num = mesh_ptr_->getNodeNum();
	std::vector<std::vector<value_t>> node(node_num);
	std::vector<std::vector<int>> triangle(tri_num);
	currents->resize(tri_num);
	for (int u = 0; u < unknowns_; ++u)
	{
		auto& rwg = ce_ptr_->getRWGRef(u);
		triangle[rwg.tpid].push_back(u);
		triangle[rwg.tmid].push_back(u);
	}

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
			for (size_t j = 0; j < 3; ++j)
			{
				if (i == j) continue;
				mag[j] += (t == rwg.tpid ? 1.0f : -1.0f) * cur_vec(unks[i]) * (len[i] / dominator) * (vx[j] - vx[i]);
			}

		}

		auto& data = (*currents)[t];
		data.v1 = vx[0];
		data.v2 = vx[1];
		data.v3 = vx[2];
		data.n[0] = nx[0];
		data.n[1] = nx[1];
		data.n[2] = nx[2];
		data.magnc = std::sqrt(cen_cur.norm());
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
	}*/

	int tri_num,start,end;
	int nv1, nv2, nvp, nvm, fp, fm;
	value_t len;
	VectorR3 vf[3], vs[3];
	tri_num = mesh_ptr_->getTriangleNum();
	currents->resize(tri_num);

	for (int x = 0; x < 1; x++)
	{
		for (int y = 0; y <1; y++)
		{
			int bedge;
			//int sed_id = ArrayToSED(x, y);
			int ver_num = mesh_ptr_->getNodeNum();

			std::vector<std::vector<value_t>> cur_ver(ver_num);
			std::vector<std::vector<int>> tri_rwg(tri_num);
			std::vector<value_t> cur_ver_mag(ver_num);
			std::vector<value_t> cur_cen(tri_num);

			//ce_ptr_->getSEDCommonEdgeSize(sed_id, start, end);
			//CVector I_vec = I_tot_eigen[y + x * Array_y];

			//遍历公共边，将公共边压入至三角形中，之后可知每个三角形的公共边
			for (int u = 0; u < unknowns_; u++)
			{
				ce_ptr_->getCommonEdge(u, nvp, nvm, fp, fm, len, bedge);

				tri_rwg[fp].push_back(u);
				tri_rwg[fm].push_back(u);
			}

			//遍历三角形，将每个三角形上的所有公共边对三个顶点以及中心点的贡献求和压入顶点电流栈
			for (int t = 0; t < tri_num; t++)
			{
				int n[3];
				VectorR3 v[3];

				Triangle tri = mesh_ptr_->getTriangleRef(t);
				tri.getVertex(n[0], n[1], n[2]);
				v[0] = mesh_ptr_->getVertex(n[0]);
				v[1] = mesh_ptr_->getVertex(n[1]);
				v[2] = mesh_ptr_->getVertex(n[2]);

				value_t area = Area(v[0], v[1], v[2]);
				VectorR3 v_cen = (v[0] + v[1] + v[2]) / 3.0f;

				auto& ce = tri_rwg[t];
				
				for (int k = 0; k < 3; k++)//遍历该三角形的三个顶点
				{
					VectorC3 cur_node_vec(0.0f, 0.0f, 0.0f);
					//VectorC3 cur_cen_vec(0.0f, 0.0f, 0.0f);
					for (int e = 0; e < ce.size(); e++)//遍历该三角形的所有公共边，计算公共边对该顶点的贡献
					{
						ce_ptr_->getCommonEdge(ce[e], nv1, nv2);
						ce_ptr_->getCommonEdge(ce[e], nvp, nvm, fp, fm, len);
						value_t sign = 1.0f;
						VectorR3 vx;

						if (t == fp)
							vx = mesh_ptr_->getVertex(nvp);
						else
						{
							vx = mesh_ptr_->getVertex(nvm);
							sign = -1.0f;
						}

						if (n[k] == nv1)
						{
							VectorR3 v1 = mesh_ptr_->getVertex(n[k]);
							cur_node_vec += sign * I(ce[e])*(v1 - vx)*len / (2.0f * area);
						}
						else if (n[k] == nv2)
						{
							VectorR3 v2 = mesh_ptr_->getVertex(n[k]);
							cur_node_vec += sign * I(ce[e])*(v2 - vx)*len / (2.0f * area);
						}
					}
					cur_ver[n[k]].push_back(std::sqrt(cur_node_vec.norm()));
				}

				VectorC3 cur_cen_vec(0.0f, 0.0f, 0.0f);
				for (int e = 0; e < ce.size(); e++)//计算所有公共边对中心点的贡献
				{
					ce_ptr_->getCommonEdge(ce[e], nv1, nv2);
					ce_ptr_->getCommonEdge(ce[e], nvp, nvp, fp, fm, len);
					
					if (t == fp)
					{
						VectorR3 vp = mesh_ptr_->getVertex(nvp);
						cur_cen_vec += I(ce[e])*(v_cen - vp)*len / (2.0f * area);
					}
					else
					{
						VectorR3 vm = mesh_ptr_->getVertex(nvm);
						cur_cen_vec += I(ce[e])*(vm - v_cen)*len / (2.0f * area);
					}
				}
				cur_cen[t]=(std::sqrt(cur_cen_vec.norm())); 
			}

			//遍历所有顶点，将顶点的电流值求和取平均
			for (int n = 0; n < ver_num; n++)
			{
				cur_ver_mag[n] = cur_ver[n].empty() ? 0.0f : std::accumulate(cur_ver[n].begin(), cur_ver[n].end(), 0.0f) / cur_ver[n].size();
			}

			//遍历三角形，将三角形的相关信息以及电流信息写入文件中
			for (int t = 0; t < tri_num; t++)
			{
				auto& data = (*currents)[t];
				auto& tri = mesh_ptr_->getTriangleRef(t);
				int n[3];
				VectorR3 v[3];
				tri.getVertex(n[0], n[1], n[2]);
				v[0] = mesh_ptr_->getVertex(n[0]);
				v[1] = mesh_ptr_->getVertex(n[1]);
				v[2] = mesh_ptr_->getVertex(n[2]);

				//shiftcoordv3(x, y, v);
				data.magn1 = cur_ver_mag[n[0]];
				data.magn2 = cur_ver_mag[n[1]];
				data.magn3 = cur_ver_mag[n[2]];
				data.magnc = cur_cen[t];
				data.v1 = v[0];
				data.v2 = v[1];
				data.v3 = v[2];
			}

			cur_ver.clear();
			tri_rwg.clear();
			cur_ver_mag.clear();
			cur_cen.clear();
		}
	}


 	Qcout << "test" << std::endl;
}
