#include "MFIE.h"
#include "ConfigLoader.h"
#include "tools.h"
#include "Mathdef.h"
#include "Quadrature.h"
#include "VectorC3.h"
#include "Mesh.h"
#include "CommonEdge.h"

using namespace component;
using namespace mom;
using namespace math;
using std::setw;

MFIE::MFIE()
{
}


MFIE::~MFIE()
{
}

void MFIE::init(component::ConfigLoader * ploader)
{
    EM::init(ploader);
    if(rtype_ == policy::ResultType::Rad)
        throw std::runtime_error("MFIE cannot solve radiation problem");

    dir_ = tool::creatFolder(dir_, folder_name_ + "_MFIE");

    auto log_name = dir_ + "/sca_runtime.log";
    logger_.open(log_name);

    mesh_ptr_ = std::make_shared<Mesh>();
    ce_ptr_ = std::make_shared<CommonEdge>();

    mesh_ptr_->loadMeshFile(ploader->getMeshPath());
    ce_ptr_->buildCommonEdge(mesh_ptr_);

    unknowns_ = ce_ptr_->getCommonEdgeNum();

    auto m_w = PI2 * incidence_.freq;
    k_ = m_w / cc;
	logger_ << ploader->getReport();
	reportInfo(logger_);
    TIME_LOG("init");
}

void MFIE::solve()
{
    SEGMENT("Solve");
    Z.zeros(unknowns_, unknowns_);
    I.zeros(unknowns_);
    V.zeros(unknowns_);

    LOG(fillZ(), "Filling impedance matrix");
    TIME_LOG("fillZ");

    LOG(fillV(), "Filling voltage vector");
    TIME_LOG("fillV");

    Qcout << setw(30) << "Solving matrix equation:";
    if (!arma::solve(I, Z, V, arma::solve_opts::fast))
        throw std::runtime_error("fail solving matrix equation");
    Qcout << "success" << std::endl;
    TIME_LOG("solve");

#ifdef _DEBUG
    if (!writeZIVData())
        Qcerr << "fail writing ZIV data" << std::endl;
#endif
}

void MFIE::output()
{
    SEGMENT("Result");
    Result result;
    LOG(result.getBistaticRCS(this, &MFIE::getBiRCS), "Calculating BistaticRCS");
    TIME_LOG("getBistaticRCS");
}

void MFIE::clear()
{
    Z.reset();
    I.reset();
    V.reset();

    mesh_ptr_->clear();
    ce_ptr_->clear();
}

void MFIE::reportInfo(Qostream & strm) const
{
    mesh_ptr_->reportInfo(strm);
    ce_ptr_->reportInfo(strm);
}

Complex MFIE::mZppKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf)
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

Complex MFIE::mZpmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf)
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

Complex MFIE::mZmpKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf)
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

Complex MFIE::mZmmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf)
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

Complex MFIE::mZppSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex MFIE::mZpmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex MFIE::mZmpSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex MFIE::mZmmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
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

Complex MFIE::mVKernel(VectorR3 * vp3, VectorR3 * vm3, VectorR3 & np, VectorR3 & nm, VectorR3 & vp, VectorR3 & vm)
{
    VectorR3 vgp[7], vgm[7];
    VectorR3 rou_p, rou_m; //ρ+-
    Complex Vm;

    Gauss7Point(vp3[0], vp3[1], vp3[2], vgp);
    Gauss7Point(vm3[0], vm3[1], vm3[2], vgm);

    VectorC3 /*Ei(0, 0, 0), */Hi(0, 0, 0);
    Complex /*Vep(0, 0), Vem(0, 0), */Vmp(0, 0), Vmm(0, 0), G(0, 0);

    for (int a = 0; a < 7; a++)
    {
        //////////////////////////////////////////////////////////////////////////
        rou_p = vgp[a] - vp;
        G = exp(-J0 * k_ * (vgp[a] ^ inc_k_));

        Hi = G * inc_h_;
        Vmp += w7[a] * (rou_p ^ (np * Hi));

        //////////////////////////////////////////////////////////////////////////
        rou_m = vm - vgm[a];
        G = exp(-J0 * k_ * (vgm[a] ^ inc_k_));

        Hi = G * inc_h_;
        Vmm += w7[a] * (rou_m ^ (nm * Hi));
    }

    Vm = (Vmp + Vmm);

    return Vm;
    //后面要乘以 （0.5*L）}
}

void MFIE::fillZ()
{
    int f_fld_plu, f_fld_min, f_src_plu, f_src_min; //face +- of source and field triangle
    VectorR3 v_fld_plu, v_fld_min, v_fld_plu3[3], v_fld_min3[3];
    VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
    Triangle tri_plu, tri_min;
    value_t l_fld, l_src;//场与源的公共边长度

    Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
    int vp, vm;

    tool::BarAndPercent bar_perc;   //
    for (int f = 0; f < unknowns_; f++)
    {
        bar_perc(f + 1, unknowns_);    //

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
                Zpp = mZppSingular(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu);
            }
            else
            {
                Zpp = mZppKernel(v_fld_plu3, v_src_plu3, v_fld_plu, v_src_plu,nml_plu);
            }
            //field face+ <--> source face-
            if (f_fld_plu == f_src_min)
            {
                Zpm = mZpmSingular(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min);
            }
            else
            {
                Zpm = mZpmKernel(v_fld_plu3, v_src_min3, v_fld_plu, v_src_min,nml_plu);
            }
            //field face- <--> source face+
            if (f_fld_min == f_src_plu)
            {
                Zmp = mZmpSingular(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu);
            }
            else
            {
                Zmp = mZmpKernel(v_fld_min3, v_src_plu3, v_fld_min, v_src_plu,nml_min);
            }
            //field face- <--> source face-
            if (f_fld_min == f_src_min)
            {
                Zmm = mZmmSingular(v_fld_min3, v_src_min3, v_fld_min, v_src_min);
            }
            else
            {
                Zmm = mZmmKernel(v_fld_min3, v_src_min3, v_fld_min, v_src_min,nml_min);
            }
            Z(f, s) = l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
        }
    }
}

void MFIE::fillV()
{
    int nvp, nvm, fp, fm;
    VectorR3 vp, vm, vp3[3], vm3[3];
    value_t ln;
    Triangle tri_plu, tri_min;

    VectorR3 vgp[7], vgm[7];
    VectorR3 rou_p, rou_m; //ρ+-
    Complex Ve(0, 0);

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

        VectorR3 np = mesh_ptr_->getTriangleRef(fp).getNormal();
        VectorR3 nm = mesh_ptr_->getTriangleRef(fm).getNormal();

        Complex vk = mVKernel(vp3, vm3, np, nm, vp, vm);
        V(u) = value_t(0.5) * ln * vk / Z0;
    }
}

value_t MFIE::getBiRCS(const VectorR3 & sca_k) const
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

bool MFIE::writeZIVData()
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
