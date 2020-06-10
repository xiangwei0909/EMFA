#include "FMM.h"
#include "Miscellaneous.h"
#include "ConfigLoader.h"
#include "Mathdef.h"
#include "Quadrature.h"
#include "tools.h"
#include "CommonEdge.h"

using namespace component;
using namespace mom;
using namespace math;
using std::setw;

FMM::FMM()
{
}


FMM::~FMM()
{
}

void FMM::init(component::ConfigLoader * ploader)
{
    EM::init(ploader);
    alpha_ = ploader->getAlpha();
    if(rtype_ == policy::ResultType::Rad)
    {
        readExcEdges(dir_ + '/' + folder_name_ + ".rad");
        alpha_ = 1;
    }
    dir_ = tool::creatFolder(dir_, folder_name_ + "_FMM_CFIE_" + std::to_string(alpha_).substr(0, 3));

    auto log_name = dir_ + (rtype_ == policy::ResultType::Rad ? "/rad_" : "/sca_") + "runtime.log";
    logger_.open(log_name);

    mesh_ptr_ = std::make_shared<Mesh>();
    ce_ptr_ = std::make_shared<CommonEdge>();

    mesh_ptr_->loadMeshFile(ploader->getMeshPath());
    ce_ptr_->buildCommonEdge(mesh_ptr_);
    unknowns_ = ce_ptr_->getCommonEdgeNum();

    fmm_box_.init(mesh_ptr_, cc / incidence_.freq, ploader->getFMMBoxLength());
    iteration_threshold_ = ploader->getIterationThreshold();
    max_iteration_num_ = ploader->getMaxIterationNum();

    k_ = PI2 * incidence_.freq / cc;
    auto kd = 1.732051f * k_ * fmm_box_.getDelta();
    sam_theta_ = static_cast<int>(kd + 6 * cbrtf(kd));
    sam_phi_ = 2 * sam_theta_;

    setBoxParameter();
	logger_ << ploader->getReport();
	reportInfo(logger_);
    reportBoxInfo();    //  write box infomation to file

    TIME_LOG("init");
}

void FMM::solve()
{
    SEGMENT("Solve");
    preCalculateKArray();

    LOG(preCalculateTrans(), "Calculating transfer");
    TIME_LOG("CalculateTrans");

    LOG(preCalculateRadAndRecv(), "Calculating rad and recv");
    TIME_LOG("CalculateRadAndRecv");

    LOG(calculateNearZ(), "Filling near impedance");
    TIME_LOG("calculateNearZ");

    LOG(fillVbyBox(), "Filling voltage vector");
    TIME_LOG("fillVbyBox");

    Qcout << setw(30) << "Iterative solving:" << std::endl;
    int iter_num = iterativeSolve();
    if (iter_num == max_iteration_num_)
        logger_ << LEVEL1 "WARNING: Iteration cannot converge" << std::endl;
    TIME_LOG("iterativeSolve");

    logger_ << LEVEL1 "Iterative number:" << iter_num << '\n';
    Qcout << LEVEL1 "Iterative number: " << iter_num << std::endl;
}

void FMM::output()
{
    SEGMENT("Result");
    Result result;
    if(rtype_ == policy::ResultType::Rad)
    {
        LOG(result.getRadiationField(this, &FMM::getEFiled), "Calculating Radiation");
        TIME_LOG("getRadiation");
    }
    else
    {
        LOG(result.getBistaticRCS(this, &FMM::getBiRCS), "Calculating BistaticRCS");
        TIME_LOG("getBistaticRCS");
    }
}

void FMM::clear()
{
    rad_.clear();
    rec_.clear();
    coeff_.reset();
    sk_.clear();
    gk_.clear();
    rcs.clear();
    I_.reset();
    V_.reset();
    nearZ_.reset();
    transfer_.clear();
    arr_k_.clear();
    box_array_.clear();
    ce_ptr_->clear();
    mesh_ptr_->clear();
}

void FMM::reportInfo(Qostream & strm) const
{
    mesh_ptr_->reportInfo(strm);
    ce_ptr_->reportInfo(strm);
    fmm_box_.reportInfo(strm);
    reportFFMInfo(strm);
}

void FMM::Test()
{
    //Legendre legendre;
    //std::cout << "Legendre:\n";
    //std::cout << "n = 3, x = 0.25: " << legendre(3, 0.25) << std::endl;
    //std::cout << "n = 4, x = 0.25: " << legendre(4, 0.25) << std::endl;
    //SphHankel sph_hankel;
    //value_t x = 1.2345f;
    //std::cout << "sphere_hankel:\n";
    //std::cout << "n = 0, x = " << x << ": " << sph_hankel(0, x) << std::endl;
    //std::cout << "n = 1, x = " << x << ": " << sph_hankel(1, x) << std::endl;
    //std::cout << "n = 2, x = " << x << ": " << sph_hankel(2, x) << std::endl;
    //std::cout << "n = 3, x = " << x << ": " << sph_hankel(3, x) << std::endl;
    //Complex h0, h1, h2, h3;
    //h0 = J0 * exp(-J0 * x) / x;
    //h1 = Complex(sin(x) / x - cos(x), cos(x) / x + sin(x)) / x;
    //h2 = ((3.0f / (x * x) - 1.0f) * sin(x) / x - (3.0f * cos(x) / (x * x))) - J0 * ((-3.0f / (x * x) + 1.0f) * cos(x) / x - 3 * sin(x) / (x * x));
    //h3 = ((15.0f / (x * x * x) - 6.0f / x) * sin(x) / x - (15.0f / (x * x) - 1.0f) * cos(x) / x) - J0 * ((-15.0f / (x * x * x) + 6.0f / x) * cos(x) / x - (15 / (x * x) - 1) * sin(x) / x);
    //std::cout << "sky book: \n";
    //std::cout << "n = 0, x = " << x << ": " << h0 << std::endl;
    //std::cout << "n = 1, x = " << x << ": " << h1 << std::endl;
    //std::cout << "n = 2, x = " << x << ": " << h2 << std::endl;
    //std::cout << "n = 3, x = " << x << ": " << h3 << std::endl;
}

///////////////////////////////////////////////////////////////////////////////

void FMM::setBoxParameter()
{
    int num = fmm_box_.boxNum();
    std::vector<Box> arr_box(num);
    for (int u = 0; u < unknowns_; ++u)
    {
        int nv1, nv2;
        ce_ptr_->getCommonEdge(u, nv1, nv2);
        auto v1 = mesh_ptr_->getVertex(nv1);
        auto v2 = mesh_ptr_->getVertex(nv2);
        int idx = fmm_box_.inWhichBox((v1 + v2) / 2.0f);
        arr_box[idx].unknown.push_back(u);
    }

    for (size_t i = 0; i < arr_box.size(); ++i)
    {
        if (!arr_box[i].unknown.empty())
        {
            arr_box[i].center = fmm_box_.getCenter(static_cast<int>(i));
            box_array_.push_back(std::move(arr_box[i]));
        }
    }

    for (size_t idx = 0; idx < box_array_.size(); ++idx)
        box_array_[idx].index = idx;

    const value_t near_limit = 1.75f * fmm_box_.getDelta();
    for (auto& cur_box : box_array_)
    {
        for (auto& nbox : box_array_)
        {
            if ((cur_box.center - nbox.center).Norm() < near_limit)
                cur_box.near.push_back(&nbox);
            else
                cur_box.far.push_back(&nbox);
        }
    }

    unk_near_ = 0;
    size_t aux = 0;
    for (auto& cur_box : box_array_)
    {
        cur_box.offset = aux;
        const size_t rwg_num = cur_box.unknown.size();
        aux += rwg_num;
        for (const auto& nbox : cur_box.near)
            unk_near_ += rwg_num * nbox->unknown.size();
    }
    assert(aux == unknowns_);
}

void FMM::prepareSolve()
{
    const size_t knum = sam_phi_ * sam_theta_;
    sk_.set_size(knum, box_array_.size());
    gk_.set_size(knum, box_array_.size());
    coeff_.set_size(knum);
    size_t cur_k = 0;
    for (int p = 0; p < sam_phi_; ++p)
        for (int t = 0; t < sam_theta_; ++t)
            coeff_(cur_k++) = gl_phi_.weight(p) * gl_theta_.weight(t) * sin(PIhalves * (gl_theta_.root(t) + 1));
    coeff_ = (k_ * k_ * Z0 / 32) * coeff_;
    assert(cur_k == knum);
}

int FMM::iterativeSolve()
{
    prepareSolve();

    //  BICGSTAB iterative algorithm
    int iter_num = 0;
    Qcx_vec x(unknowns_, arma::fill::zeros);
    Qcx_vec r0 = V_ - matrixVectorMultiply(x);
    Qcx_vec r1 = r0;
    Complex rou0(1, 0), alph(1, 0), omiga(1, 0);
    Qcx_vec v, p0, p1, s, t;
    v.zeros(arma::size(x));
    p0.zeros(arma::size(x));
    Complex rou1, beta;
    const auto b_normal = arma::norm(V_);

    auto oldState = Qcout.flags();
    Qcout << std::fixed << std::setprecision(10) << std::right;
    while(iter_num < max_iteration_num_)
    {
        ++iter_num;
        rou1 = arma::cdot(r1, r0);
        beta = (rou1 / rou0) * (alph / omiga);
        p1 = r0 + beta * (p0 - omiga * v);
        v = matrixVectorMultiply(p1);
        alph = rou1 / (arma::cdot(r1, v));
        s = r0 - alph * v;
        x = x + alph * p1;
        auto resdual = arma::norm(s) / b_normal;
        Qcout << std::setw(7)  << iter_num << " of " << max_iteration_num_
              << std::setw(15) << resdual  << std::endl;
        if (resdual < iteration_threshold_)
            break;
        t = matrixVectorMultiply(s);
        omiga = arma::cdot(t, s) / arma::cdot(t, t);
        x = x + omiga * s;
        r0 = s - omiga * t;
        resdual = arma::norm(s) / b_normal;
        if (resdual < iteration_threshold_)
            break;
        rou0 = rou1;
        p0 = p1;
    }
    I_ = x;
    Qcout.flags(oldState);
    return iter_num;
}

Qcx_vec FMM::matrixVectorMultiply(const Qcx_vec& b)
{
    const size_t knum = sk_.row();
    for (const auto& cur_box : box_array_)
    {
        const size_t idx = cur_box.index;
        const size_t rwg_num = cur_box.unknown.size();
        const size_t base = cur_box.offset;
        const auto& rad_mat = rad_[idx];
        for (size_t k = 0; k < knum; ++k)
        {
            VectorC3 val(0, 0, 0);
            for (size_t u = 0; u < rwg_num; ++u)
                val += rad_mat(u, k) * b(base + u);
            sk_(k, idx) = val;
        }
    }

    gk_.reset(VectorC3(0, 0, 0));
    for (const auto& cur_box : box_array_)
    {
        const size_t idx = cur_box.index;
        const auto& trans = transfer_[idx];
        size_t cur_k = 0;
        for (const auto fptr : cur_box.far)
        {
            const size_t fidx = fptr->index;
            for (size_t k = 0; k < knum; ++k)
                gk_(k, idx) += trans(cur_k++) * sk_(k, fidx);
        }
        assert(cur_k == trans.n_elem);
    }

    Qcx_vec far(unknowns_);
    size_t cur_u = 0;
    for (const auto& cur_box : box_array_)
    {
        const size_t idx = cur_box.index;
        const size_t rwg_num = cur_box.unknown.size();
        const auto& rec_mat = rec_[idx];
        for (size_t u = 0; u < rwg_num; ++u)
        {
            Complex val(0, 0);
            for (size_t k = 0; k < knum; ++k)
                val += coeff_(k) * (rec_mat(k, u) ^ gk_(k, idx));
            far(cur_u++) = val;
        }
    }
    assert(cur_u == unknowns_);
    return nearZ_ * b + far;
}

void FMM::preCalculateKArray()
{
    gl_theta_.set(sam_theta_);
    gl_phi_.set(sam_phi_);

    std::vector<value_t> sint(sam_theta_), cost(sam_theta_);
    for (int i = 0; i < sam_theta_; ++i)
    {
        value_t theta = PIhalves * (gl_theta_.root(i) + 1);
        sint[i] = sin(theta);
        cost[i] = cos(theta);
    }

    std::vector<value_t> sinp(sam_phi_), cosp(sam_phi_);
    for (int i = 0; i < sam_phi_; ++i)
    {
        value_t phi = PI * (gl_phi_.root(i) + 1);
        sinp[i] = sin(phi);
        cosp[i] = cos(phi);
    }

    arr_k_.set_size(sam_theta_, sam_phi_);
    for (int p = 0; p < sam_phi_; ++p)
        for (int t = 0; t < sam_theta_; ++t)
            arr_k_(t, p) = VectorR3(sint[t] * cosp[p], sint[t] * sinp[p], cost[t]);
}

void FMM::preCalculateTrans()
{
    transfer_.resize(box_array_.size());
    const size_t knum = sam_theta_ * sam_phi_;
    const Complex j_l[4] = { 1, -J0, -1, J0 };
    SphHankel sph_hankel;
    Legendre legendre;
    tool::BarAndPercent bar_perc;                       //  
    for (size_t idx = 0; idx < box_array_.size(); ++idx)
    {
        bar_perc(idx + 1, box_array_.size());                //  

        const auto& box = box_array_[idx];
        const auto& far = box.far;
        const auto& center = box.center;
        transfer_[idx].set_size(knum, far.size());  //  K x Nfar
        for (size_t fid = 0; fid < far.size(); ++fid)
        {
            for (size_t k = 0; k < knum; ++k)
            {
                const auto rvec = center - far[fid]->center;
                const value_t vec_norm = rvec.Norm();
                const value_t kd_norm = k_ * vec_norm;
                const value_t kd_unit = (arr_k_(k) ^ rvec) / vec_norm;
                int minL = static_cast<int>(sam_theta_ < kd_norm ? sam_theta_ : kd_norm);
                Complex val(0, 0);
                for (int cl = 0; cl <= minL; ++cl)
                    val += (2 * cl + 1.0f) * j_l[cl % 4] * sph_hankel(kd_norm, cl) * legendre(kd_unit, cl);
                transfer_[idx](k, fid) = val;
            }
        } // end for far boxes
    }
}

void FMM::preCalculateRadAndRecv()
{
    const size_t box_num = box_array_.size();
    const size_t knum = sam_phi_ * sam_theta_;
    rad_.resize(box_num);
    rec_.resize(box_num);

    tool::BarAndPercent bar_perc;   //
    for (size_t idx = 0; idx < box_num; ++idx)
    {
        bar_perc(idx + 1, box_num);     //

        const auto& cur_box = box_array_[idx];
        const auto& edges = cur_box.unknown;
        const auto& center = cur_box.center;
        const size_t rwg_num = edges.size();
        auto& rad_mat = rad_[idx];
        auto& rec_mat = rec_[idx];
        rad_mat.set_size(rwg_num, knum);    //  RWG x K
        rec_mat.set_size(knum, rwg_num);    //  K x RWG

        VectorR3 v_plu, v_min, nml_plu, nml_min;
        VectorR3 v_plu3[3], v_min3[3], vpg3[3], vmg3[3];

        for (size_t i = 0; i < rwg_num; ++i)
        {
            int u = edges[i];
            int nvp, nvm, nfp, nfm;
            value_t length;
            ce_ptr_->getCommonEdge(u, nvp, nvm, nfp, nfm, length);
            v_plu = mesh_ptr_->getVertex(nvp);
            v_min = mesh_ptr_->getVertex(nvm);

            auto& tri_plu = mesh_ptr_->getTriangleRef(nfp);
            auto& tri_min = mesh_ptr_->getTriangleRef(nfm);
            tri_plu.getVertex(v_plu3[0], v_plu3[1], v_plu3[2]);
            tri_min.getVertex(v_min3[0], v_min3[1], v_min3[2]);
            nml_plu = tri_plu.getNormal();
            nml_min = tri_min.getNormal();

            Gauss3Point(v_plu3[0], v_plu3[1], v_plu3[2], vpg3);
            Gauss3Point(v_min3[0], v_min3[1], v_min3[2], vmg3);

            for (size_t k = 0; k < knum; ++k)
            {
                const auto tmp_rad = radiationFunction(vpg3, vmg3, v_plu, v_min, center, arr_k_(k));
                const auto tmp_rec = halfRecvFunction(vpg3, vmg3, v_plu, v_min, nml_plu, nml_min, center, arr_k_(k));
                rad_mat(i, k) = (length / 2) * tmp_rad;
                rec_mat(k, i) = (length / 2) * (alpha_ * conj(tmp_rad) + (1 - alpha_) * tmp_rec);
            }
        } 
    } // end for box
}

VectorC3 FMM::radiationFunction(const VectorR3* vpgs3, const VectorR3* vmgs3, const VectorR3& vp,
                                const VectorR3& vm, const VectorR3& center, const VectorR3& sam_k)
{
    VectorC3 val(0, 0, 0);
    VectorR3 rou_plu, rou_min;
    for (int p = 0; p < 3; ++p)
    {
        rou_plu = vpgs3[p] - vp;
        rou_min = vm - vmgs3[p];
        val += w3[p] * rou_plu * std::exp(J0 * k_ * (sam_k ^ (vpgs3[p] - center)));
        val += w3[p] * rou_min * std::exp(J0 * k_ * (sam_k ^ (vmgs3[p] - center)));
    }
    return val - ((sam_k ^ val) * sam_k);
}

VectorC3 FMM::halfRecvFunction(const VectorR3 * vpgs3, const VectorR3 * vmgs3, const VectorR3 & vp,
                               const VectorR3 & vm, const VectorR3 & nml_plu, const VectorR3 & nml_min,
                               const VectorR3 & center, const VectorR3 & sam_k)
{
    VectorC3 val(0, 0, 0);
    VectorR3 rou_plu, rou_min;
    for (int p = 0; p < 3; ++p)
    {
        rou_plu = vpgs3[p] - vp;
        rou_min = vm - vmgs3[p];
        val += w3[p] * (nml_plu * rou_plu) * std::exp(-J0 * k_ * (sam_k ^ (vpgs3[p] - center)));
        val += w3[p] * (nml_min * rou_min) * std::exp(-J0 * k_ * (sam_k ^ (vmgs3[p] - center)));
    }
    return sam_k * val;
}

void FMM::calculateNearZ()
{
    Qumat   location(2, unk_near_);
    Qcx_vec elements(unk_near_);
    size_t  aux = 0, cur_rwg = 0;

    tool::BarAndPercent bar_perc;       //  
    for (size_t idx = 0; idx < box_array_.size(); ++idx)
    {
        bar_perc(idx + 1, box_array_.size());   //

        const auto& sbox = box_array_[idx];
        for (int s : sbox.unknown)
        {
            for (const auto nptr : sbox.near)
                fillNearZbyColumn(nptr->unknown, s, nptr->offset, cur_rwg, location, elements, aux);
            ++cur_rwg;
        }
    }
    assert(cur_rwg == unknowns_);
    assert(aux == unk_near_);
    nearZ_ = Qsp_cx_mat(location, elements, false);
}

void FMM::fillNearZbyColumn(const std::vector<int>& flds, int sid, size_t fpos, size_t spos,
                            Qumat& loc, Qcx_vec& elems, size_t& aux)
{
    int f_fld_plu, f_fld_min, f_src_plu, f_src_min; 
    VectorR3 v_fld_plu, v_fld_min, v_fld_plu3[3], v_fld_min3[3];
    VectorR3 v_src_plu, v_src_min, v_src_plu3[3], v_src_min3[3];
    value_t l_fld, l_src;

    Complex Zpp(0, 0), Zpm(0, 0), Zmp(0, 0), Zmm(0, 0);
    int nvp, nvm;

    ce_ptr_->getCommonEdge(sid, nvp, nvm, f_src_plu, f_src_min, l_src);

    v_src_plu = mesh_ptr_->getVertex(nvp);
    v_src_min = mesh_ptr_->getVertex(nvm);

    auto& tri_plu_s = mesh_ptr_->getTriangleRef(f_src_plu);
    auto& tri_min_s = mesh_ptr_->getTriangleRef(f_src_min);

    tri_plu_s.getVertex(v_src_plu3[0], v_src_plu3[1], v_src_plu3[2]);
    tri_min_s.getVertex(v_src_min3[0], v_src_min3[1], v_src_min3[2]);

    auto fRWG = flds.size();
    for (size_t f = 0; f < fRWG; ++f)
    {
        ce_ptr_->getCommonEdge(flds[f], nvp, nvm, f_fld_plu, f_fld_min, l_fld);

        v_fld_plu = mesh_ptr_->getVertex(nvp);
        v_fld_min = mesh_ptr_->getVertex(nvm);

        auto& tri_plu_f = mesh_ptr_->getTriangleRef(f_fld_plu);
        auto& tri_min_f = mesh_ptr_->getTriangleRef(f_fld_min);

        tri_plu_f.getVertex(v_fld_plu3[0], v_fld_plu3[1], v_fld_plu3[2]);
        tri_min_f.getVertex(v_fld_min3[0], v_fld_min3[1], v_fld_min3[2]);

        VectorR3 nml_plu = tri_plu_f.getNormal();
        VectorR3 nml_min = tri_min_f.getNormal();

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
        //  inserted by column-major order
        loc(0, aux) = fpos + f;
        loc(1, aux) = spos;
        elems(aux) = l_fld * l_src * (Zpp + Zpm + Zmp + Zmm);
        ++aux;
    }
}

void mom::FMM::fillVbyBox()
{
    V_.set_size(unknowns_);
    size_t curpos = 0;

    tool::BarAndPercent bar_perc;
    for (size_t idx = 0; idx < box_array_.size(); ++idx)
    {
        bar_perc(idx + 1, box_array_.size());

        const auto& box = box_array_[idx];
        assert(curpos == box.offset);

        for (int eid : box.unknown)
        {
            int nvp, nvm, fp, fm;
            VectorR3 vp, vm, vp3[3], vm3[3];
            value_t length;

            ce_ptr_->getCommonEdge(eid, nvp, nvm, fp, fm, length);
            vp = mesh_ptr_->getVertex(nvp);
            vm = mesh_ptr_->getVertex(nvm);
            auto& tri_plu = mesh_ptr_->getTriangleRef(fp);
            auto& tri_min = mesh_ptr_->getTriangleRef(fm);
            tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);
            tri_min.getVertex(vm3[0], vm3[1], vm3[2]);

            VectorR3 nml_plu = tri_plu.getNormal();
            VectorR3 nml_min = tri_min.getNormal();

            V_(curpos++) = 0.5f * length * cVKernel(vp3, vm3, nml_plu, nml_min, vp, vm);
        }
    }
    assert(curpos == unknowns_);
}

bool FMM::readExcEdges(const Qstring & rad_file, Qstring & stateInfo)
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

void FMM::readExcEdges(const Qstring & rad_file)
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

void FMM::radiateV()
{
    int nv1, nv2;

    tool::BarAndPercent bar_perc;   //
    for (int u = 0; u < unknowns_; ++u)
    {
        bar_perc(u + 1, unknowns_); //

        ce_ptr_->getCommonEdge(u, nv1, nv2);
        const value_t length = ce_ptr_->getCommonEdgeLength(u);

        for (const auto& elem : exc_edge_)
        {
            if ((elem.first == nv1 && elem.second == nv2) || (elem.first == nv2 && elem.second == nv1))
            {
                if (elem.first == nv1)
                    V_(u) = length;
                else
                    V_(u) = -length;
                break;
            }
        }
    }
}

value_t FMM::getBiRCS(const VectorR3 & sca_k) const
{
    size_t curpos = 0;

    int nvp, nvm, fp, fm;
    VectorR3 vp, vm, vp3[3], vm3[3];
    value_t L;
    VectorR3 vgp[7], vgm[7];
    VectorR3 rou_p, rou_m;
    Complex G0(0, 0);
    VectorC3 Es(0, 0, 0);

    for (const auto& cur_box : box_array_)
    {
        auto& edges = cur_box.unknown;
        for (int u : edges)
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
                G0 = exp(J0 * k_ * (vgp[g] ^ sca_k));
                Esp = Esp + rou_p * w7[g] * G0;

                rou_m = vm - vgm[g];
                G0 = exp(J0 * k_ * (vgm[g] ^ sca_k));
                Esm = Esm + rou_m * w7[g] * G0;
            }
            Es = Es + (Esp + Esm) * (I_(curpos++) * L);
        }
    }
    assert(curpos == unknowns_);

    auto es = 0.5f * sca_k * (sca_k * Es);
    value_t rcs = Z0 * Z0 * k_ * k_ * es.norm() / PI4;
    return 10 * log10(rcs);
}

void FMM::getEFiled(component::FieldData * pdata, const VectorR3 & rad_k, const VectorR3 & rad_ev, const VectorR3 & rad_eh) const
{

}

void FMM::reportFFMInfo(Qostream & strm) const
{
    auto old_flag = strm.flags();
    size_t elems = 0;
    for (const auto& cur_box : box_array_)
        elems += cur_box.far.size();
    elems = elems * sam_phi_ * sam_theta_;
    auto near_mem = unk_near_ * 2 * sizeof(value_t) / 1024 / 1024;
    auto rdcv_mem = 12 * unknowns_ * sam_phi_ * sam_theta_ * sizeof(value_t) / 1024 / 1024;
    auto tran_mem = 2 * elems * sizeof(value_t) / 1024 / 1024;
    auto orig_mem = 2 * sizeof(value_t) * unknowns_ * unknowns_ / 1024 / 1024;
    auto total_mem = near_mem + rdcv_mem + tran_mem;
    bool overflow = total_mem > orig_mem;
    auto saved_mem = overflow ? (total_mem - orig_mem) : (orig_mem - total_mem);
    strm << HEADING "FMM Information:\n" TRAILING
        << LEVEL1 "Total boxes: " << fmm_box_.boxNum() << '\n'
        << LEVEL1 "Valid boxes: " << box_array_.size() << '\n'
        << LEVEL1 "Average per box: " << unknowns_ / box_array_.size() << '\n'
        << LEVEL1 "Number of near: " << unk_near_ << '\n'
        << LEVEL1 "Sparse rate: " << 100.0f * unk_near_ / unknowns_ / unknowns_ << "%\n"
        << LEVEL1 "Number of multipoles: " << sam_theta_ << '\n'
        << LEVEL1 "Near impedance: " << near_mem << " MB" << FORMAT_MEMORY(near_mem) << '\n'
        << LEVEL1 "Transfer: " << tran_mem << " MB" << FORMAT_MEMORY(tran_mem) << '\n'
        << LEVEL1 "Rad and Recv: " << rdcv_mem << " MB" << FORMAT_MEMORY(rdcv_mem) << '\n'
        << LEVEL1 "Total memory" << total_mem << " MB" << FORMAT_MEMORY(total_mem) << '\n'
        << LEVEL1 "Saved memory: " << (overflow ? '-' : '\0') << saved_mem << " MB" << FORMAT_MEMORY(saved_mem) << '\n';
    strm.flush();
    strm.flags(old_flag);
}

void mom::FMM::reportBoxInfo() const
{
    Qofstream ffm_log(dir_ + "/ffm_box_log.txt", std::ios::out);
    ffm_log << std::left;
    reportFFMInfo(ffm_log);
    ffm_log << std::right;
    for (const auto& cur_box : box_array_)
    {
        ffm_log << "    -> Box " << setw(4) << cur_box.index + 1 << ":    near: " << setw(3) << cur_box.near.size()
                << "    far: " << setw(5) << cur_box.far.size()
                << "    unknowns: " << setw(5) << cur_box.unknown.size() << '\n';
    }
    ffm_log.close();
}

Complex FMM::cZppKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf)
{
    return integral::cfieZKernel(alpha_, k_, vf3, vs3, vf, vs, nf, 1, 1);
}

Complex FMM::cZpmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf)
{
    return integral::cfieZKernel(alpha_, k_, vf3, vs3, vf, vs, nf, 1, -1);
}

Complex FMM::cZmpKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf)
{
    return integral::cfieZKernel(alpha_, k_, vf3, vs3, vf, vs, nf, -1, 1);
}

Complex FMM::cZmmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 & nf)
{
    return integral::cfieZKernel(alpha_, k_, vf3, vs3, vf, vs, nf, -1, -1);
}

Complex FMM::cZppSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
    return integral::cfieZSingular(alpha_, k_, vf3, vs3, vf, vs, 1, 1);
}

Complex FMM::cZpmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
    return integral::cfieZSingular(alpha_, k_, vf3, vs3, vf, vs, 1, -1);
}

Complex FMM::cZmpSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
    return integral::cfieZSingular(alpha_, k_, vf3, vs3, vf, vs, -1, 1);
}

Complex FMM::cZmmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs)
{
    return integral::cfieZSingular(alpha_, k_, vf3, vs3, vf, vs, -1, -1);
}

Complex FMM::cVKernel(VectorR3 * vp3, VectorR3 * vm3, VectorR3 & np, VectorR3 & nm, VectorR3 & vp, VectorR3 & vm)
{
    return integral::cfieVKernel(alpha_, k_, inc_k_, inc_e_, inc_h_, vp3, vm3, np, nm, vp, vm);
}
