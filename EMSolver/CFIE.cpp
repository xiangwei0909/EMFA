#include "CFIE.h"
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

CFIE::CFIE()
{
}

CFIE::~CFIE()
{
}

void CFIE::init(component::ConfigLoader * ploader)
{
    EM::init(ploader);
    alpha_ = ploader->getAlpha();
    if(rtype_ == policy::ResultType::Rad)
    {
        readExcEdges(dir_ + '/' + folder_name_ + ".rad");
        alpha_ = 1;
    }
    dir_ = tool::creatFolder(dir_, folder_name_ + "_CFIE_" + std::to_string(alpha_).substr(0, 3));

    auto log_name = dir_ + (rtype_ == policy::ResultType::Rad ? "/rad_" : "/sca_") + "runtime.log";
    logger_.open(log_name);

    mesh_ptr_ = std::make_shared<Mesh>();
    ce_ptr_ = std::make_shared<CommonEdge>();

    mesh_ptr_->loadMeshFile(ploader->getMeshPath());
    ce_ptr_->buildCommonEdge(mesh_ptr_);

    unknowns_ = ce_ptr_->getCommonEdgeNum();
    k_ = PI2 * incidence_.freq / cc;

	logger_ << ploader->getReport();
    reportInfo(logger_);
    TIME_LOG("init");
}

void CFIE::solve()
{
    SEGMENT("Solve");
    LOG(fillZ(), "Filling impedance matrix");
    TIME_LOG("fillZ");

    if(rtype_ == policy::ResultType::Rad)
    {
        LOG(radiateV(), "Filling rad voltage vector");
    }
    else
    {
        LOG(fillV(), "Filling sca voltage vector");
    }
    TIME_LOG("fillV");

    Qcout << setw(30) << "Solving matrix equation:";
    //if (!arma::solve(I, Z, V, arma::solve_opts::fast))
    //    throw std::runtime_error("fail solving matrix equation");
    //Qcout << "success" << std::endl;

    Qcout << '\n';
    int iter_num = 0;
    if (!iml::BICGSTAB(Z, V, I, 0.01f, 500, iter_num))
        logger_ << LEVEL1 "WARNING: Iteration cannot converge" << std::endl;
    logger_ << LEVEL1 "Iterative number: " << iter_num << '\n';
    Qcout << LEVEL1 "Iterative number: " << iter_num << std::endl;

    TIME_LOG("solve");

#ifdef _DEBUG
    if (!writeZIVData())
        Qcerr << "fail writing ZIV data" << std::endl;
    TIME_LOG("writeZIV");
#endif
}

void CFIE::output()
{
    SEGMENT("Result");
    Result result;
    if(rtype_ == policy::ResultType::Rad)
    {
        LOG(result.getRadiationField(this, &CFIE::getEFiled), "Calculating Radiation");
        TIME_LOG("getRadiation");
    }
    else
    {
        LOG(result.getBistaticRCS(this, &CFIE::getBiRCS), "Calculating BistaticRCS");
        TIME_LOG("getBistaticRCS");
    }
    LOG(result.getCurrentDistribution(this, &CFIE::calculateSurfaceCurrent), "Calculating surface current");
    TIME_LOG("getSurfaceCurrent");
}

void CFIE::clear()
{
    Z.reset();
    I.reset();
    V.reset();
    exc_edge_.clear();
    mesh_ptr_->clear();
    ce_ptr_->clear();
}

void CFIE::reportInfo(Qostream & strm) const
{
    mesh_ptr_->reportInfo(strm);
    ce_ptr_->reportInfo(strm);
}

Complex CFIE::cZppKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf)
{
    return integral::cfieZKernel(alpha_, k_, vf3, vs3, vf, vs, nf, 1, 1);
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

Complex CFIE::cZpmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf)
{
    return integral::cfieZKernel(alpha_, k_, vf3, vs3, vf, vs, nf, 1, -1);
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

Complex CFIE::cZmpKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf)
{
    return integral::cfieZKernel(alpha_, k_, vf3, vs3, vf, vs, nf, -1, 1);
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

Complex CFIE::cZmmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf)
{
    return integral::cfieZKernel(alpha_, k_, vf3, vs3, vf, vs, nf, -1, -1);
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

Complex CFIE::cZppSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
    return integral::cfieZSingular(alpha_, k_, vf3, vs3, vf, vs, 1, 1);
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

Complex CFIE::cZpmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
    return integral::cfieZSingular(alpha_, k_, vf3, vs3, vf, vs, 1, -1);
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

Complex CFIE::cZmpSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
    return integral::cfieZSingular(alpha_, k_, vf3, vs3, vf, vs, -1, 1);
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

Complex CFIE::cZmmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
    return integral::cfieZSingular(alpha_, k_, vf3, vs3, vf, vs, -1, -1);
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

Complex CFIE::cVKernel(VectorR3 * vp3, VectorR3 * vm3, VectorR3 & np, VectorR3 & nm, VectorR3 & vp, VectorR3 & vm)
{
    return integral::cfieVKernel(alpha_, k_, inc_k_, inc_e_, inc_h_, vp3, vm3, np, nm, vp, vm);
    //VectorR3 vgp[3], vgm[3];
    //VectorR3 rou_p, rou_m; //ρ+-
    //Complex Ve, Vm;

    //Gauss3Point(vp3[0], vp3[1], vp3[2], vgp);
    //Gauss3Point(vm3[0], vm3[1], vm3[2], vgm);

    //VectorC3 Ei(0, 0, 0), Hi(0, 0, 0);
    //Complex Vep(0, 0), Vem(0, 0), Vmp(0, 0), Vmm(0, 0), G(0, 0);

    //for (int a = 0; a < 3; a++)
    //{
    //    //////////////////////////////////////////////////////////////////////////
    //    rou_p = vgp[a] - vp;
    //    G = exp(-J0 * k_ * (vgp[a] ^ inc_k_));

    //    Ei = G * inc_e_;
    //    Vep += w3[a] * (rou_p ^ Ei);
    //    Hi = G * inc_h_;//这里多了个Z0
    //    Vmp += w3[a] * (rou_p ^ (np * Hi));

    //    //////////////////////////////////////////////////////////////////////////
    //    rou_m = vm - vgm[a];
    //    G = exp(-J0 * k_ * (vgm[a] ^ inc_k_));

    //    Ei = G * inc_e_;
    //    Vem += w3[a] * (rou_m ^ Ei);
    //    Hi = G * inc_h_;//这里多了个Z0
    //    Vmm += w3[a] * (rou_m ^ (nm * Hi));
    //}

    //Ve = (Vep + Vem);
    //Vm = (Vmp + Vmm);

    //return alpha_ * Ve + (1 - alpha_) * Vm;   //这里就不需要Z0
    //                                            //后面要乘以 （0.5*L）
}

void CFIE::fillZ()
{
    Z.zeros(unknowns_, unknowns_);

    int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle
    VectorR3 v_fld_plu, v_fld_min, v_fld_plu3[3], v_fld_min3[3];
    VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
    Triangle tri_plu, tri_min;
    value_t l_fld, l_src;//场与源的公共边长度

    Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
    int vp, vm;

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

        VectorR3 nml_plu = tri_plu.getNormal();
        VectorR3 nml_min = tri_min.getNormal();

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
                Zpp = cZppSingular(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
            }
            else
            {
                Zpp = cZppKernel(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu, nml_plu);
            }
            //field face+ <--> source face-
            if (f_fld_plu == f_src_min)
            {
                Zpm = cZpmSingular(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
            }
            else
            {
                Zpm = cZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min, nml_plu);
            }
            //field face- <--> source face+
            if (f_fld_min == f_src_plu)
            {
                Zmp = cZmpSingular(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
            }
            else
            {
                Zmp = cZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu, nml_min);
            }
            //field face- <--> source face-
            if (f_fld_min == f_src_min)
            {
                Zmm = cZmmSingular(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
            }
            else
            {
                Zmm = cZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min, nml_min);
            }
            Z(f, s) = l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
        }
    }
}

void CFIE::fillV()
{
    V.zeros(unknowns_);

    int nvp, nvm, fp, fm;
    VectorR3 vp, vm, vp3[3], vm3[3];
    value_t ln;

    tool::BarAndPercent bar_perc;
    for (int u = 0; u < unknowns_; u++)
    {
        bar_perc(u + 1, unknowns_);

        ce_ptr_->getCommonEdge(u, nvp, nvm, fp, fm, ln);
        vp = mesh_ptr_->getVertex(nvp);
        vm = mesh_ptr_->getVertex(nvm);
        auto& tri_plu = mesh_ptr_->getTriangleRef(fp);
        auto& tri_min = mesh_ptr_->getTriangleRef(fm);
        tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);
        tri_min.getVertex(vm3[0], vm3[1], vm3[2]);

        VectorR3 np = mesh_ptr_->getTriangleRef(fp).getNormal();
        VectorR3 nm = mesh_ptr_->getTriangleRef(fm).getNormal();

        Complex vk = cVKernel(vp3, vm3, np, nm, vp, vm);
        V(u) = 0.5f * ln * vk;
    }
}

void CFIE::readExcEdges(const Qstring & rad_file)
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

void CFIE::radiateV()
{
    V.zeros(unknowns_);

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

value_t CFIE::getBiRCS(const VectorR3& sca_k) const
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
    return 10 * log10(rcs);
}

void CFIE::getEFiled(component::FieldData * pdata, const VectorR3 & rad_k, const VectorR3 & rad_ev, const VectorR3 & rad_eh) const
{
    int nvp, nvm, fp, fm;
    VectorR3 vp, vm, vp3[3], vm3[3];
    value_t L;
    VectorR3 vgp[7], vgm[7], rou_p, rou_m;
    Complex  G0(0, 0);
    VectorC3 Es(0, 0, 0);

    for (int u = 0; u < unknowns_; u++)
    {
        ce_ptr_->getCommonEdge(u, nvp, nvm, fp, fm, L);
        vp = mesh_ptr_->getVertex(nvp);
        vm = mesh_ptr_->getVertex(nvm);
        auto& tri_plu = mesh_ptr_->getTriangleRef(fp);
        auto& tri_min = mesh_ptr_->getTriangleRef(fm);
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
    const auto coeff = -J0 * Z0 * k_ / (2.0f * PI4);
    pdata->ftheta = coeff * (Es ^ rad_ev);
    pdata->fphi = coeff * (Es ^ rad_eh);
}

void CFIE::calculateSurfaceCurrent(std::vector<component::CurrentData>* currents) const
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

bool CFIE::writeZIVData()
{
    Qofstream outputZ(dir_ + '/' + "matrix_Z.txt");
    if (outputZ.fail())
        return false;
    Z.save(outputZ, arma::arma_ascii);

    Qofstream outputV(dir_ + '/' + "matrix_V.txt");
    if (outputV.fail())
        return false;
    V.save(outputV, arma::arma_ascii);

    Qofstream outputI(dir_ + '/' + "martix_I.txt");
    if (outputI.fail())
        return false;
    I.save(outputI, arma::arma_ascii);

    return true;
}
