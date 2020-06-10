#include "MLFMA.h"
#include "Miscellaneous.h"
#include "ConfigLoader.h"
#include "Mathdef.h"
#include "Quadrature.h"
#include "tools.h"
#include "Mesh.h"
#include "CommonEdge.h"

using namespace component;
using namespace mom;
using namespace math;
using std::setw;

MLFMA::MLFMA()
{
}

MLFMA::~MLFMA()
{
}

void MLFMA::init(component::ConfigLoader * ploader)
{
    EM::init(ploader);
    alpha_ = ploader->getAlpha();
    if(rtype_ == policy::ResultType::Rad)
    {
        //readExcEdges(dir_ + '/' + folder_name_ + ".rad");
        alpha_ = 1;
    }
    dir_ = tool::creatFolder(dir_, folder_name_ + "_MLFMA_CFIE_" + std::to_string(alpha_).substr(0, 3));

    mesh_ptr_ = std::make_shared<Mesh>();
    ce_ptr_ = std::make_shared<CommonEdge>();

    mesh_ptr_->loadMeshFile(ploader->getMeshPath());
    ce_ptr_->buildCommonEdge(mesh_ptr_);
    unknowns_ = ce_ptr_->getCommonEdgeNum();

    iteration_threshold_ = ploader->getIterationThreshold();
    max_iteration_num_ = ploader->getMaxIterationNum();
    preconditioning_ = ploader->getPreconditioningSwitch();
    row_threshold_ = ploader->getRowThreshold();
    col_threshold_ = ploader->getColThreshold();

    k_ = PI2 * incidence_.freq / cc;
    boundary_.init(mesh_ptr_, cc / incidence_.freq, ploader->getMLFMABoxThreshold());
    level_ = boundary_.getLevel();
    if (level_ < 2)
        throw std::runtime_error("level number of MLFMA is too small (MLFMA box threshold is too large)");
    
    pool_.init(ploader->getThreadNumber());
    auto task_num = pool_.getThreadNum() * ploader->getTaskFactor();

    buildBox();
    prepareAssignment(task_num);

    reportLevelDetail();
    auto log_name = dir_ + (rtype_ == policy::ResultType::Rad ? "/rad_" : "/sca_") + "parallel_runtime.log";
    logger_.open(log_name);
	logger_ << ploader->getReport();
    reportInfo(logger_);
    TIME_LOG("init");
}

void MLFMA::solve()
{
    SEGMENT("Solve");
    //if (!Validate())
    //{
    //    Qcout << "Validiation failed" << std::endl;
    //    return;
    //}

    levelRootsAndWeights();
    calculateInterpolationMtr();
    //Test();
    calculateLevelK();

    LOG(calculateLevelTrans(), "Calculating transfer");
    TIME_LOG("calculateLevelTrans");

    LOG(calculateRadAndRecv(), "Calculating rad and recv");
    TIME_LOG("calculateRadAndRecv");

    LOG(calculateNearZ(), "Filling near impedance");
    TIME_LOG("calculateNearZ");

    LOG(fillVbyBox(), "Filling voltage vector");
    TIME_LOG("fillVbyBox");

    if (preconditioning_)
    {
        LOG(calculatePreconditioner(), "Calculating preconditioner");
        TIME_LOG("Preconditioner");
    }

    Qcout << setw(30) << "Iterative solving:" << std::endl;
    int iter_num = iterativeSolve();
    if (iter_num == max_iteration_num_)
        logger_ << LEVEL1 "WARNING: Iteration cannot converge" << std::endl;
    logger_ << LEVEL1 "Iterative number: " << iter_num << '\n';
    Qcout << LEVEL1 "Iterative number: " << iter_num << std::endl;
    TIME_LOG("iterativeSolve");
}

void MLFMA::output()
{
    SEGMENT("Result");
    Result result;
    if(rtype_ == policy::ResultType::Rad)
    {
        //LOG(result.getRadiationField(this, &MLFMA::getEFiled), "Calculating Radiation");
        TIME_LOG("getRadiation");
    }
    else
    {
        LOG(result.getBistaticRCS(this, &MLFMA::getBiRCS), "Calculating BistaticRCS");
        TIME_LOG("getBistaticRCS");
    }
    LOG(result.getCurrentDistribution(this, &MLFMA::calculateSurfaceCurrent), "Calculating surface current");
    TIME_LOG("getSurfaceCurrent");
}

void MLFMA::clear()
{
    pool_.clear();
    assignment_.reset();
    precond_.reset();
    nearZ_.reset();
    V_.reset();
    I_.reset();
    coeff_.reset();
    level_rad_.clear();
    level_rec_.clear();
    rad_.clear();
    rec_.clear();
    level_trans_.clear();
    level_k_.clear();
    arr_ip_.clear();
    glxp_.clear();
    glxt_.clear();
    glw_.clear();
    box_len_.clear();
    level_box_.clear();
    ce_ptr_->clear();
    mesh_ptr_->clear();
}

void MLFMA::reportInfo(Qostream & strm) const
{
    mesh_ptr_->reportInfo(strm);
    ce_ptr_->reportInfo(strm);
    boundary_.reportInfo(strm);
    pool_.reportInfo(strm);
    this->reportLevelInfo(strm);
    this->reportTaskAssigmentInfo(strm);
    this->reportMLFMAInfo(strm);
}

void MLFMA::Test()
{
    auto func = [](value_t x, value_t y) { return 2 * x * x + 3 * y * y; };

    auto& samx1 = glxt_.back();
    auto& samy1 = glxp_.back();
    auto& gw1 = glw_.back();
    auto& samx2 = *(glxt_.end() - 2);
    auto& samy2 = *(glxp_.end() - 2);
    auto& gw2 = *(glw_.end() - 2);

    size_t x1_num = samx1.n_rows;
    size_t y1_num = samy1.n_rows;
    size_t x2_num = samx2.n_rows;
    size_t y2_num = samy2.n_rows;

    Qcout << std::setprecision(10);

    const auto& ip_mat = arr_ip_.back();
    Qvec f1(x1_num * y1_num, arma::fill::zeros);
    for (int y = 0; y < y1_num; ++y)
        for (int x = 0; x < x1_num; ++x)
            f1(y * x1_num + x) = func(samx1(x), samy1(y));

    Qvec f2(x2_num * y2_num, arma::fill::zeros);

    size_t cur_num = 0;
    for (size_t cur_row = 0; cur_row < x2_num * y2_num; ++cur_row)
    {
        for (size_t c = 0; c < 9; ++c)
        {
            assert(cur_row == ip_mat.row(cur_num));
            f2(cur_row) += ip_mat(cur_num) * f1(ip_mat.col(cur_num));
            ++cur_num;
        }
    }
    assert(cur_num == 9 * x2_num * y2_num);

    for (size_t y = 0; y < y2_num; ++y)
    {
        for (size_t x = 0; x < x2_num; ++x)
        {
            size_t cur_id = y * x2_num + x;
            Qcout << cur_id << ": " << f2(cur_id) << "    " << func(samx2(x), samy2(y)) << '\n';
        }
    }

    Qvec f3(x1_num * y1_num, arma::fill::zeros);
    for (cur_num = 0; cur_num < ip_mat.size(); ++cur_num)
    {
        const size_t ck = ip_mat.col(cur_num);
        const size_t pk = ip_mat.row(cur_num);
        f3(ck) += ip_mat(cur_num) * f2(pk) * gw2(pk) / gw1(ck);
    }

    for (size_t y = 0; y < y1_num; ++y)
    {
        for (size_t x = 0; x < x1_num; ++x)
        {
            size_t cur_id = y * x1_num + x;
            Qcout << cur_id << " : " << f3(cur_id) << "    " << f1(cur_id) << '\n';
        }
    }
}

bool MLFMA::Validate()
{
    Qcout << setw(30) << "Validation Information:";
    //  parent validation
    for (int cur_level = 0; cur_level < level_; ++cur_level)
    {
        for (auto& cur_box : level_box_[cur_level])
        {
            bool well = true;
            for (auto cptr : cur_box.child)
            {
                if (&cur_box != cptr->parent)
                {
                    well = false;
                    break;
                }
            }
            if (!well)
            {
                Qcout << "parent -> child -> parent error" << std::endl;
                return false;
            }
        }
    }
    for (int cur_level = 1; cur_level <= level_; ++cur_level)
    {
        size_t total = 0;
        for (auto& cur_box : level_box_[cur_level])
        {
            total += cur_box.unknown.size();

            //  child validation
            bool well = false;
            for (auto sptr : cur_box.parent->child)
            {
                if (&cur_box == sptr)
                {
                    well = true;
                    break;
                }
            }
            if (!well)
            {
                Qcout << "child -> parent -> child error" << std::endl;
                return false;
            }

            //  near validation
            for (auto nptr : cur_box.near)
            {
                well = false;
                for (auto sptr : nptr->near)
                {
                    if (&cur_box == sptr)
                    {
                        well = true;
                        break;
                    }
                }
                if (!well)
                {
                    Qcout << "near <-> near error" << std::endl;
                    return false;
                }
            }

            // far validation
            for (auto fptr : cur_box.far)
            {
                well = false;
                for (auto sptr : fptr->far)
                {
                    if (&cur_box == sptr)
                    {
                        well = true;
                        break;
                    }
                }
                if (!well)
                {
                    Qcout << "far <-> far error" << std::endl;
                    return false;
                }
            }
            
        } // box
        if (total != unknowns_)
        {
            Qcout << "expect != acutal\nacutal: " << total << ", expect: " << unknowns_ << std::endl;
            return false;
        }
    } // level

    //  offset validation
    auto iter = level_box_.back().begin();
    for (; iter != level_box_.back().end() - 1; ++iter)
    {
        if (iter->offset + iter->unknown.size() != (iter + 1)->offset)
        {
            Qcout << "current offset + current size != next offest" << std::endl;
            return false;
        }
    }
    if (iter->offset + iter->unknown.size() != unknowns_)
    {
        Qcout << "last offset + last size != unknowns" << std::endl;
        return false;
    }

    Qcout << "pass (^_^)" << std::endl;
    return true;
}

//////////////////////////////////////////////////////////////

void MLFMA::prepareAssignment(size_t thread_num)
{
    assignment_.zeros(thread_num + 1, level_ - 1);
    for (int cur_level = 2; cur_level <= level_; ++cur_level)
    {
        size_t idx = cur_level - 2;
        size_t box_num = level_box_[cur_level].size();
        size_t average = box_num / thread_num;
        size_t remaining = box_num % thread_num;
        size_t cur_num = 0;
        for (size_t i = 0; i <= thread_num; ++i)
        {
            assignment_(i, idx) = cur_num;
            cur_num += average;
            if (remaining)
            {
                ++cur_num;
                --remaining;
            }
        }
        assert(cur_num = box_num + average);
    }
}

void MLFMA::allocateForSolve()
{
    level_rad_.resize(level_ - 1);
    level_rec_.resize(level_ - 1);
    for (int cur_level = 2; cur_level <= level_; ++cur_level)
    {
        int idx = cur_level - 2;
        auto& up_mat = level_rad_[idx];
        auto& down_mat = level_rec_[idx];
        auto knum = level_k_[idx].size();
        auto bnum = level_box_[cur_level].size();
        up_mat.set_size(knum, bnum);
        down_mat.set_size(knum, bnum);
    }
}

void MLFMA::calculateCoeff()
{
    const auto& gw = glw_.back();
    const auto& gt = glxt_.back();
    const size_t knum = gw.n_elem;
    const size_t pnum = glxp_.back().n_elem;
    const size_t tnum = gt.n_elem;
    coeff_.set_size(knum);

    size_t cur_k = 0;
    for (size_t p = 0; p < pnum; ++p)
        for (size_t t = 0; t < tnum; ++t, ++cur_k)
            coeff_(cur_k) = gw(cur_k) * sin(PIhalves * (gt(t) + 1));

    assert(cur_k == level_k_.back().size());

    const auto const_num = k_ * k_ * Z0 / 32;
    coeff_ = const_num * coeff_;
}

int MLFMA::iterativeSolve()
{
    allocateForSolve();
    calculateCoeff();

    if (preconditioning_)
        V_ = precond_ * V_;
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

    const int old_flag = Qcout.flags();
    Qcout << std::fixed << std::setprecision(10) << std::right;
    while (iter_num < max_iteration_num_)
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
        Qcout << std::setw(7) << iter_num << " of " << max_iteration_num_
              << std::setw(15) << resdual << std::endl;
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
    Qcout.flags(old_flag);
    return iter_num;
}

Qcx_vec MLFMA::matrixVectorMultiply(const Qcx_vec & b)
{
    for (int cur_level = 2; cur_level <= level_; ++cur_level)
        level_rad_[cur_level - 2].reset(VectorC3(0, 0, 0));

    const size_t task_num = assignment_.n_rows - 1;
    size_t idx = level_ - 2;
    for (size_t t = 0; t < task_num; ++t)
    {
        size_t begin = assignment_(t, idx);
        size_t end = assignment_(t + 1, idx);
        if (begin < end)
        {
            auto task = std::bind(&MLFMA::parallelFinestLayerUpward, this, begin, end, std::cref(b));
            pool_.submit(task);
        }
    }
    pool_.run();

    upwardPass();

    for (int cur_level = 2; cur_level <= level_; ++cur_level)
        level_rec_[cur_level - 2].reset(VectorC3(0, 0, 0));
    idx = 0;
    for (size_t t = 0; t < task_num; ++t)
    {
        size_t begin = assignment_(t, idx);
        size_t end = assignment_(t + 1, idx);
        if (begin < end)
        {
            auto task = std::bind(&MLFMA::parallelFinestLayerDownward, this, begin, end);
            pool_.submit(task);
        }
    }
    pool_.run();

    downwardPass();

    Qcx_vec far(unknowns_, arma::fill::zeros);
    idx = level_ - 2;
    for (size_t t = 0; t < task_num; ++t)
    {
        size_t begin = assignment_(t, idx);
        size_t end = assignment_(t + 1, idx);
        if (begin < end)
        {
            auto task = std::bind(&MLFMA::parallelFar, this, begin, end, &far);
            pool_.submit(task);
        }
    }
    pool_.run();

    if (preconditioning_)
        return precond_ * (nearZ_ * b + far);
    return nearZ_ * b + far;
}

void MLFMA::parallelFinestLayerUpward(size_t begin, size_t end, const Qcx_vec & b)
{
    auto& up_mat = level_rad_.back();
    size_t knum = up_mat.row();
    for (size_t ib = begin; ib < end; ++ib)
    {
        auto& cur_box = level_box_.back()[ib];
        size_t rwg_num = cur_box.unknown.size();
        size_t base = cur_box.offset;
        size_t bidx = cur_box.index;
        const auto& rad_mat = rad_[bidx];
        for (size_t k = 0; k < knum; ++k)
        {
            VectorC3 val(0, 0, 0);
            for (size_t u = 0; u < rwg_num; ++u)
                val += b(base + u) * rad_mat(u, k);
            up_mat(k, bidx) = val;
        }
    }
}

void MLFMA::parallelFinestLayerDownward(size_t begin, size_t end)
{
    auto& down_mat = level_rec_.front();
    const auto& u_mat = level_rad_.front();
    const auto& arr_trans = level_trans_.front();
    const auto& box_array = level_box_[2];
    const size_t knum = u_mat.row();
    for (size_t ib = begin; ib < end; ++ib)
    {
        size_t bidx = box_array[ib].index;
        const auto& trans = arr_trans[bidx];
        size_t aux = 0;
        for (const auto fptr : box_array[ib].far)
        {
            size_t fidx = fptr->index;
            for (size_t k = 0; k < knum; ++k)
                down_mat(k, bidx) += trans(aux++) * u_mat(k, fidx);
        }
        assert(aux == trans.n_elem);
    }
}

void MLFMA::parallelFar(size_t begin, size_t end, Qcx_vec * pfar) const
{

    const auto& d_mat = level_rec_.back();
    const auto& box_array = level_box_.back();
    const size_t knum = d_mat.row();
    for (size_t ib = begin; ib < end; ++ib)
    {
        const size_t rwg_num = box_array[ib].unknown.size();
        const size_t bidx = box_array[ib].index;
        const size_t offset = box_array[ib].offset;
        const auto& rec_mat = rec_[bidx];
        for (size_t u = 0; u < rwg_num; ++u)
        {
            Complex val(0, 0);
            for (size_t k = 0; k < knum; ++k)
                val += coeff_(k) * (rec_mat(k, u) ^ d_mat(k, bidx));
            (*pfar)(offset + u) = val;
        }
    }
}

void MLFMA::upwardPass()
{
    for (int cur_level = level_ - 1; cur_level >= 2; --cur_level)
    {
        const size_t task_num = assignment_.n_rows - 1;
        size_t idx = cur_level - 1;
        for (size_t t = 0; t < task_num; ++t)
        {
            size_t begin = assignment_(t, idx);
            size_t end = assignment_(t + 1, idx);
            if (begin < end)
            {
                auto task = std::bind(&MLFMA::parallelUpwardPass, this, begin, end, cur_level);
                pool_.submit(task);
            }
        }
        pool_.run();
    } // previous level
}

void MLFMA::parallelUpwardPass(size_t begin, size_t end, int cur_level)
{
    int idx = cur_level - 2;
    const auto& arr_cbox = level_box_[cur_level + 1];
    const auto& ip_mat = arr_ip_[idx];
    const auto& in_mat = level_rad_[idx + 1];
    const auto& arr_k = level_k_[idx];
    const auto knum = arr_k.size();
    wt::QMat<VectorC3> out_mat;
    out_mat.set_size(level_rad_[idx].row(), level_rad_[idx].col(), VectorC3());

    for (size_t cidx = begin; cidx < end; ++cidx)
    {
        VectorR3 rvec = arr_cbox[cidx].parent->center - arr_cbox[cidx].center;
        size_t bid = arr_cbox[cidx].parent->index;
        size_t cur_num = 0;
        for (size_t k = 0; k < knum; ++k)
        {
            const Complex phase = exp(-J0 * k_ * (arr_k(k) ^ rvec));
            VectorC3 val(0, 0, 0);
            for (size_t p = 0; p < 9; ++p)
            {
                val += ip_mat(cur_num) * in_mat(ip_mat.col(cur_num), cidx);
                ++cur_num;
            }
            out_mat(k, bid) += phase * val;
        }
        assert(cur_num == 9 * out_mat.row());
    }
    MutexGuard guard(mutex_);
    level_rad_[idx] += out_mat;
}

void MLFMA::downwardPass()
{
    for (int cur_level = 3; cur_level <= level_; ++cur_level)
    {
        const size_t task_num = assignment_.n_rows - 1;
        size_t idx = cur_level - 2;
        for (size_t t = 0; t < task_num; ++t)
        {
            size_t begin = assignment_(t, idx);
            size_t end = assignment_(t + 1, idx);
            if (begin < end)
            {
                auto task = std::bind(&MLFMA::parallelDownwardPass, this, begin, end, idx);
                pool_.submit(task);
            }
        }
        pool_.run();
    } // next level
}

void MLFMA::parallelDownwardPass(size_t begin, size_t end, size_t idx)
{
    const auto& arr_trans = level_trans_[idx];
    const auto& up_mat = level_rad_[idx];
    auto& down_mat = level_rec_[idx];
    const auto& pdown = level_rec_[idx - 1];
    const auto& ip_mat = arr_ip_[idx - 1];
    const auto& arr_k = level_k_[idx - 1];
    const auto& gw = glw_[idx];
    const auto& pgw = glw_[idx - 1];
    const auto& box_array = level_box_[idx + 2];
    const size_t knum = up_mat.row();
    for (size_t ib = begin; ib < end; ++ib)
    {
        auto& cur_box = box_array[ib];
        size_t bidx = cur_box.index;
        const auto& trans = arr_trans[bidx];
        size_t aux = 0;
        for (const auto fptr : cur_box.far)
        {
            size_t fidx = fptr->index;
            for (size_t k = 0; k < knum; ++k)
                down_mat(k, bidx) += trans(aux++) * up_mat(k, fidx);
        }
        assert(aux == trans.n_elem);

        size_t pidx = cur_box.parent->index;
        VectorR3 rvec = cur_box.parent->center - cur_box.center;
        for (size_t cur_num = 0; cur_num < ip_mat.size(); ++cur_num)
        {
            const size_t ck = ip_mat.col(cur_num);
            const size_t pk = ip_mat.row(cur_num);
            const value_t coeff = ip_mat(cur_num) * pgw(pk) / gw(ck);
            const Complex phase = exp(J0 * k_ * (arr_k(pk) ^ rvec));
            down_mat(ck, bidx) += coeff * phase * pdown(pk, pidx);
        }
    }
}

void MLFMA::buildBox()
{
    std::vector<VectorR3> center_points(unknowns_);
    fillCenterPoints(center_points);//公共边中心点

    level_box_.resize(level_ + 1);
    box_len_.resize(level_ + 1);
    value_t cur_box_len = boundary_.getMaxLength() / 2;
    buildFirstLevel(center_points, cur_box_len);
    box_len_[1] = cur_box_len;

    for (int cur_level = 2; cur_level <= level_; ++cur_level)
    {
        cur_box_len /= 2;
        box_len_[cur_level] = cur_box_len;
        buildCurLevel(center_points, cur_level, cur_box_len);
    }

    setBoxParameter();

    size_t total_near = 0;
    for (const auto& cur_box : level_box_.back())
    {
        size_t srwg = cur_box.unknown.size();
        size_t crwg = 0;
        for (auto nptr : cur_box.near)
            crwg += nptr->unknown.size();
        total_near += srwg * crwg;
    }
    unk_near_ = total_near;

    col_threshold_ *= box_len_.back();
    row_threshold_ *= box_len_.back();
    std::vector<size_t> cols;
    precond_num_ = 0;
    if (preconditioning_)
    {
        for (const auto& cur_box : level_box_.back())
        {
            for (const auto nptr : cur_box.near)
            {
                auto& unks = nptr->unknown;
                for (size_t i = 0; i < unks.size(); ++i)
                {
                    auto dist = (center_points[unks[i]] - cur_box.center).Norm();
                    if (dist <= col_threshold_) cols.push_back(nptr->offset + i);
                }
            }
            precond_num_ += cols.size() * cur_box.unknown.size();
            cols.clear();
        }
    }
}

void MLFMA::fillCenterPoints(std::vector<VectorR3>& centers)
{
    int nv1, nv2;
    for (int u = 0; u < unknowns_; ++u)
    {
        ce_ptr_->getCommonEdge(u, nv1, nv2);
        auto v1 = mesh_ptr_->getVertex(nv1);
        auto v2 = mesh_ptr_->getVertex(nv2);
        centers[u] = (v1 + v2) / 2;
    }
}

void MLFMA::buildFirstLevel(const std::vector<VectorR3>& cps, value_t box_len)
{
    Box root;
    root.center = boundary_.getCenter();
    level_box_[0].push_back(root);

    auto center = level_box_[0][0].center;
    Box childbox[8];
    for (int x = 0; x < 2; ++x)
    {
        for (int y = 0; y < 2; ++y)
        {
            for (int z = 0; z < 2; ++z)
            {
                int index = (x * 2 + y) * 2 + z;
                childbox[index].center = center + box_len * VectorR3(x - 0.5f, y - 0.5f, z - 0.5f);
                childbox[index].parent = &level_box_[0][0];
            }
        }
    }
    for (int u = 0; u < unknowns_; ++u)
    {
        int x_flag = cps[u].x < center.x ? 0 : 1;
        int y_flag = cps[u].y < center.y ? 0 : 1;
        int z_flag = cps[u].z < center.z ? 0 : 1;
        int index = (x_flag * 2 + y_flag) * 2 + z_flag;
        childbox[index].unknown.push_back(u);
    }
    size_t cur_pos = 0;
    for (int c = 0; c < 8; ++c)
    {
        if (!childbox[c].unknown.empty())
        {
            childbox[c].offset = cur_pos;
            cur_pos += childbox[c].unknown.size();
            level_box_[1].push_back(childbox[c]);
        }
    }
    assert(cur_pos == unknowns_);
}

void MLFMA::buildCurLevel(const std::vector<VectorR3>& cps, int cur_level, value_t box_len)
{
    auto& parents = level_box_[cur_level - 1];
    for (auto& pbox : parents)
    {
        Box childbox[8];
        auto center = pbox.center;
        for (int x = 0; x < 2; ++x)
        {
            for (int y = 0; y < 2; ++y)
            {
                for (int z = 0; z < 2; ++z)
                {
                    int index = (x * 2 + y) * 2 + z;
                    childbox[index].center = center + box_len * VectorR3(x - 0.5f, y - 0.5f, z - 0.5f);
                    childbox[index].parent = &pbox;
                }
            }
        }
        for (int u : pbox.unknown)
        {
            int x_flag = cps[u].x < center.x ? 0 : 1;
            int y_flag = cps[u].y < center.y ? 0 : 1;
            int z_flag = cps[u].z < center.z ? 0 : 1;
            int index = (x_flag * 2 + y_flag) * 2 + z_flag;
            childbox[index].unknown.push_back(u);
        }

        size_t cur_pos = pbox.offset;
        for (int c = 0; c < 8; ++c)
        {
            if (!childbox[c].unknown.empty())
            {
                childbox[c].offset = cur_pos;
                cur_pos += childbox[c].unknown.size();
                level_box_[cur_level].push_back(childbox[c]);
            }
        }
    }
}

void MLFMA::setBoxParameter()
{
    int cur_level = 1;
    for (; cur_level <= level_; ++cur_level)
        for (auto& cur_box : level_box_[cur_level])
            cur_box.parent->child.push_back(&cur_box);

    // level = 1
    for (size_t f = 0; f < level_box_[1].size(); ++f)
        for (size_t s = 0; s < level_box_[1].size(); ++s)
            level_box_[1][f].near.push_back(&level_box_[1][s]);

    value_t cur_box_len = boundary_.getMaxLength() / 4;
    // level > 1
    for (cur_level = 2; cur_level <= level_; ++cur_level, cur_box_len /= 2)
    {
        value_t len_limit = 1.75f * cur_box_len;
        for (auto& cur_box : level_box_[cur_level])
        {
            auto& parent = cur_box.parent;
            for (auto pnear : parent->near)
            {
                for (auto near_cptr : pnear->child)
                {
                    value_t dist = (near_cptr->center - cur_box.center).Norm();
                    if (dist < len_limit)
                        cur_box.near.push_back(near_cptr);
                    else
                        cur_box.far.push_back(near_cptr);
                }
            }
        } // current box
    } // current level
    for (cur_level = 2; cur_level <= level_; ++cur_level)
    {
        auto& arr_box = level_box_[cur_level];
        for (size_t id = 0; id < arr_box.size(); ++id)
            arr_box[id].index = id;
    }
}

int MLFMA::getMultipoleNumber(value_t box_len) const
{
    value_t kd = 1.732051f * k_ * box_len;
    return static_cast<int>(kd + 6 * cbrtf(kd));
}

void MLFMA::levelRootsAndWeights()
{
    glxt_.resize(level_ - 1);
    glxp_.resize(level_ - 1);
    glw_.resize(level_ - 1);
    for (int cur_level = 2; cur_level <= level_; ++cur_level)
    {
        int cur_poles = getMultipoleNumber(box_len_[cur_level]);
        auto& theta_vec = glxt_[cur_level - 2];
        auto& phi_vec = glxp_[cur_level - 2];
        auto& weight_vec = glw_[cur_level - 2];
        theta_vec.set_size(cur_poles);
        phi_vec.set_size(2 * cur_poles);
        weight_vec.set_size(2 * cur_poles * cur_poles);

        GaussLegendre gl_theta, gl_phi;
        gl_theta.set(cur_poles);
        for (int t = 0; t < cur_poles; ++t)
            theta_vec(t) = gl_theta.root(t);
        gl_phi.set(2 * cur_poles);
        for (int p = 0; p < 2 * cur_poles; ++p)
            phi_vec(p) = gl_phi.root(p);

        size_t cur_gw_idx = 0;
        for (int p = 0; p < 2 * cur_poles; ++p)
        {
            value_t phiw = gl_phi.weight(p);
            for (int t = 0; t < cur_poles; ++t, ++cur_gw_idx)
                weight_vec(cur_gw_idx) = phiw * gl_theta.weight(t);
        }
        assert(cur_gw_idx == weight_vec.n_elem);
    }
}

void MLFMA::calculateInterpolationMtr()
{
    arr_ip_.resize(level_ - 2);

    Qumat   tloc, ploc;
    Qmat    telems, pelems;
    for (int cur_level = 2; cur_level < level_; ++cur_level)
    {
        const auto& ct_vec = glxt_[cur_level - 2];
        const auto& cp_vec = glxp_[cur_level - 2];
        const auto& nt_vec = glxt_[cur_level - 1];
        const auto& np_vec = glxp_[cur_level - 1];
        const size_t ct_num = ct_vec.n_rows;
        const size_t nt_num = nt_vec.n_rows;

        auto& ip_mat = arr_ip_[cur_level - 2];
        ip_mat.set_size(2 * ct_num * ct_num, 9);

        lagrangeIPWeights(nt_vec, ct_vec, tloc, telems);
        lagrangeIPWeights(np_vec, cp_vec, ploc, pelems);

        size_t theta_num = telems.n_cols, phi_num = pelems.n_cols;
        size_t cur_num = 0, cur_row = 0;
        for (size_t p = 0; p < phi_num; ++p)
        {
            for (size_t t = 0; t < theta_num; ++t, ++cur_row)
            {
                for (size_t ip = 0; ip < 3; ++ip)
                {
                    size_t base = ploc(ip, p) * nt_num;
                    for (size_t it = 0; it < 3; ++it)
                        ip_mat.set_value(cur_row, base + tloc(it, t), cur_num++, pelems(ip, p) * telems(it, t));
                } // fill 9 nonzero elements
            }
        }
        assert(cur_num == 18 * ct_num * ct_num);
        assert(cur_row == glw_[cur_level - 2].n_elem);
    } // end for level
}

void MLFMA::lagrangeIPWeights(const Qvec & nxt_sam, const Qvec & cur_sam, Qumat & loc, Qmat & elems) const
{
    size_t cols = cur_sam.n_rows;
    loc.set_size(3, cols);
    elems.set_size(3, cols);
    size_t nlast = nxt_sam.n_rows - 1;
    size_t nid = 2, cid = 0;
    size_t pos1, pos2, pos3;
    value_t low = nxt_sam(1), high = nxt_sam(nlast - 1);
    while (cid < cols)
    {
        const value_t xx = cur_sam(cid);
        if (xx < low)
        {
            pos1 = 0;
            pos2 = 1;
            pos3 = 2;
        }
        else if (xx > high)
        {
            pos1 = nlast - 2;
            pos2 = nlast - 1;
            pos3 = nlast;
        }
        else
        {
            nid += nxt_sam(nid) < xx ? 1 : 0;
            size_t offset = (nxt_sam(nid) - xx < xx - nxt_sam(nid - 1)) ? 1 : 0;
            pos1 = nid + offset - 2;
            pos2 = nid + offset - 1;
            pos3 = nid + offset;
        }
        loc(0, cid) = pos1;
        loc(1, cid) = pos2;
        loc(2, cid) = pos3;

        const auto x0 = nxt_sam(pos1);
        const auto x1 = nxt_sam(pos2);
        const auto x2 = nxt_sam(pos3);
        //  three-points lagrange interpolation weights
        elems(0, cid) = (xx - x1) * (xx - x2) / ((x0 - x1) * (x0 - x2));
        elems(1, cid) = (xx - x0) * (xx - x2) / ((x1 - x0) * (x1 - x2));
        elems(2, cid) = (xx - x0) * (xx - x1) / ((x2 - x0) * (x2 - x1));
        ++cid;
    }
}

void MLFMA::calculateNearZ()
{
    Qumat   location(2, unk_near_);
    Qcx_vec elements(unk_near_);
    const auto& box_array = level_box_.back();
    size_t cur_base = 0;
    size_t task_num = assignment_.n_rows - 1;
    size_t idx = level_ - 2;
    for (size_t t = 0; t < task_num; ++t)
    {
        size_t begin = assignment_(t, idx);
        size_t end = assignment_(t + 1, idx);
        if (begin < end)
        {
            size_t task_unk = 0;
            for (size_t ib = begin; ib < end; ++ib)
            {
                size_t srwg = box_array[ib].unknown.size();
                size_t frwg = 0;
                for (const auto nptr : box_array[ib].near)
                    frwg += nptr->unknown.size();
                task_unk += srwg * frwg;
            }
            auto task = std::bind(&MLFMA::parallelNearZ, this, begin, end, &location, &elements, cur_base);
            pool_.submit(task, task_unk);
            cur_base += task_unk;
        }
    }
    pool_.interactiveRun();
    nearZ_ = Qsp_cx_mat(location, elements, false);
}

void MLFMA::parallelNearZ(size_t begin, size_t end, Qumat * ploc, Qcx_vec * pelems, size_t ebase) const
{
    for (size_t sb = begin; sb < end; ++sb)
    {
        auto& cur_box = level_box_.back()[sb];
        auto& srcs = cur_box.unknown;
        for (size_t s = 0; s < srcs.size(); ++s)
        {
            size_t spos = cur_box.offset + s;
            for (const auto nptr : cur_box.near)
                fillNearZbyColumn(nptr->unknown, srcs[s], nptr->offset, spos, *ploc, *pelems, ebase);
        }
    }
}

void MLFMA::fillNearZbyColumn(const std::vector<int>& flds, int sid, size_t fpos, size_t spos, Qumat & loc, Qcx_vec & elems, size_t & aux) const
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

void MLFMA::calculateLevelK()
{
    level_k_.resize(level_ - 1);
    Qvec    tsin, tcos, psin, pcos;
    for (int cur_level = 2; cur_level <= level_; ++cur_level)
    {
        int idx = cur_level - 2;
        const auto& tvec = glxt_[idx];
        const auto& pvec = glxp_[idx];
        const size_t tnum = tvec.n_rows;
        const size_t pnum = pvec.n_rows;
        tsin.set_size(tnum);
        tcos.set_size(tnum);
        psin.set_size(pnum);
        pcos.set_size(pnum);
        for (size_t t = 0; t < tnum; ++t)
        {
            value_t theta = PIhalves * (tvec(t) + 1);
            tsin(t) = sin(theta);
            tcos(t) = cos(theta);
        }
        for (size_t p = 0; p < pnum; ++p)
        {
            value_t phi = PI * (pvec(p) + 1);
            psin(p) = sin(phi);
            pcos(p) = cos(phi);
        }
        auto& arr_k = level_k_[idx];
        arr_k.set_size(tnum, pnum);
        for (size_t p = 0; p < pnum; ++p)
            for (size_t t = 0; t < tnum; ++t)
                arr_k(t, p) = VectorR3(tsin(t) * pcos(p), tsin(t) * psin(p), tcos(t));
    }
}

void MLFMA::calculateLevelTrans()
{
    level_trans_.resize(level_ - 1);
    size_t task_num = assignment_.n_rows - 1;
    for (int cur_level = 2; cur_level <= level_; ++cur_level)
    {
        size_t idx = cur_level - 2;
        size_t box_num = assignment_(task_num, idx);
        level_trans_[idx].resize(box_num);
        for (size_t t = 0; t < task_num; ++t)
        {
            size_t begin = assignment_(t, idx);
            size_t end = assignment_(t + 1, idx);
            if (begin < end)
            {
                auto task = std::bind(&MLFMA::parallelLevelTrans, this, begin, end, idx);
                pool_.submit(task);
            }
        }
        pool_.run();
    }
}

void MLFMA::parallelLevelTrans(size_t begin, size_t end, size_t idx)
{
    const Complex j_l[4] = { 1, -J0, -1, J0 };
    Legendre legendre;
    SphHankel sph_hankel;

    const auto& arr_k = level_k_[idx];
    const auto& box_vec = level_box_[idx + 2];
    auto& arr_trans = level_trans_[idx];
    for (size_t ib = begin; ib < end; ++ib)
    {
        const auto& cur_box = box_vec[ib];
        arr_trans[ib].set_size(arr_k.size(), cur_box.far.size()); // K x Nfar
        size_t aux = 0;
        for (const auto fptr : cur_box.far)
        {
            auto vec = cur_box.center - fptr->center;
            value_t vec_norm = vec.Norm();
            value_t kd_norm = k_ * vec_norm;
            for (size_t kth = 0; kth < arr_k.size(); ++kth)
            {
                value_t kd_unit = (arr_k(kth) ^ vec) / vec_norm;
                int minL = static_cast<int>(arr_k.row() < kd_norm ? arr_k.row() : kd_norm);
                Complex val(0, 0);
                for (int cl = 0; cl <= minL; ++cl)
                    val += (2 * cl + 1.0f) * j_l[cl % 4] * sph_hankel(kd_norm, cl) * legendre(kd_unit, cl);
                arr_trans[ib](aux++) = val;
            }
        }
        assert(aux == arr_trans[ib].n_elem);
    }
}

void MLFMA::calculateRadAndRecv()
{
    size_t task_num = assignment_.n_rows - 1;
    size_t idx = level_ - 2;
    size_t box_num = assignment_(task_num, idx);
    rad_.resize(box_num);
    rec_.resize(box_num);
    for (size_t t = 0; t < task_num; ++t)
    {
        size_t begin = assignment_(t, idx);
        size_t end = assignment_(t + 1, idx);
        if (begin < end)
        {
            auto task = std::bind(&MLFMA::parallelRadAndRecv, this, begin, end);
            pool_.submit(task);
        }
    }
    pool_.interactiveRun();
}

void MLFMA::parallelRadAndRecv(size_t begin, size_t end)
{
    const auto& box_array = level_box_.back();
    const auto& arr_k = level_k_.back();
    const size_t knum = arr_k.size();
    for (size_t ib = begin; ib < end; ++ib)
    {
        auto& center = box_array[ib].center;
        auto& edges = box_array[ib].unknown;
        const size_t rwg_num = edges.size();
        auto& rad_mat = rad_[ib];
        auto& rec_mat = rec_[ib];
        rad_mat.set_size(rwg_num, knum);    // RWG x K
        rec_mat.set_size(knum, rwg_num);    // K x RWG

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
                const auto& sam_k = arr_k(k);
                auto tmp_rad = radiationFunction(vpg3, vmg3, v_plu, v_min, center, sam_k);
                auto tmp_rec = halfRecvFunction(vpg3, vmg3, v_plu, v_min, nml_plu, nml_min, center, sam_k);
                rad_mat(i, k) = (length / 2) * tmp_rad;
                rec_mat(k, i) = (length / 2) * ((alpha_ * conj(tmp_rad)) + ((1 - alpha_) * tmp_rec));
            }
        } // end for rwg_num
    } // end for box_num

}

void MLFMA::calculatePreconditioner()
{
    std::vector<VectorR3> center_points(unknowns_);
    fillCenterPoints(center_points);

    std::vector<size_t> base(assignment_.n_rows);
    preparepreconditioner(center_points, base);

    Qumat location(2, precond_num_);
    Qcx_vec elements(precond_num_);
    size_t task_num = assignment_.n_rows - 1;
    size_t idx = level_ - 2;
    for (size_t t = 0; t < task_num; ++t)
    {
        size_t begin = assignment_(t, idx);
        size_t end = assignment_(t + 1, idx);
        if (begin < end)
        {
            auto task = std::bind(&MLFMA::parallelPreconditioner, this, begin, end, std::cref(center_points), &location, &elements, base[t]);
            pool_.submit(task, base[t + 1] - base[t]);
        }
    }
    pool_.interactiveRun();
    precond_ = Qsp_cx_mat(location, elements, false);
}

void MLFMA::preparepreconditioner(const std::vector<VectorR3>& centers, std::vector<size_t>& base) const
{
    size_t cur_base = 0;
    size_t level_id = level_ - 2, idx = 0;
    for (const auto& cur_box : level_box_.back())
    {
        if (cur_box.index == assignment_(idx, level_id))
        {
            base[idx] = cur_base;
            ++idx;
        }
        size_t cols = 0;
        for (const auto nptr : cur_box.near)
        {
            auto& unks = nptr->unknown;
            for (size_t i = 0; i < unks.size(); ++i)
            {
                auto dist = (centers[unks[i]] - cur_box.center).Norm();
                if (dist <= col_threshold_) ++cols;
            }
        }
        cur_base += cols * cur_box.unknown.size();
    }
    base.back() = cur_base;
    assert(cur_base == precond_num_);
}

void MLFMA::parallelPreconditioner(size_t begin, size_t end, const std::vector<VectorR3>& centers, Qumat * ploc, Qcx_vec * pelems, size_t ebase) const
{
    const auto& box_array = level_box_.back();
    std::vector<size_t> rows, cols;
    Qcx_mat near_imp, matQ, matR;
    size_t cur_col = 0, cur_num = ebase;
    for (size_t ib = begin; ib < end; ++ib)
    {
        auto& cur_box = box_array[ib];
        size_t offset = 0;
        for (const auto nptr : cur_box.near)
        {
            if (nptr == &cur_box) offset = rows.size();
            auto& unks = nptr->unknown;
            for (size_t i = 0; i < unks.size(); ++i)
            {
                auto dist = (centers[unks[i]] - cur_box.center).Norm();
                if (dist <= row_threshold_) rows.push_back(nptr->offset + i);
                if (dist <= col_threshold_) cols.push_back(nptr->offset + i);
            }
        }
        near_imp.set_size(rows.size(), cols.size());
        for (size_t s = 0; s < cols.size(); ++s)
            for (size_t f = 0; f < rows.size(); ++f)
                near_imp(f, s) = nearZ_(rows[f], cols[s]);
        {
            MutexGuard guard(mutex_);
            arma::qr_econ(matQ, matR, near_imp);
        }
        Qcx_mat res = arma::inv(matR) * matQ.submat(offset, 0, arma::size(cur_box.unknown.size(), cols.size())).t();
        cur_col = cur_box.offset;
        for (size_t u = 0; u < cur_box.unknown.size(); ++u)
        {
            for (size_t c = 0; c < cols.size(); ++c, ++cur_num)
            {
                (*ploc)(0, cur_num) = cols[c];
                (*ploc)(1, cur_num) = cur_col;
                (*pelems)(cur_num) = res(c, u);
            }
            ++cur_col;
        }
        rows.clear();
        cols.clear();
    }
}

VectorC3 MLFMA::radiationFunction(const VectorR3* vpgs3, const VectorR3* vmgs3, const VectorR3& vp,
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

VectorC3 MLFMA::halfRecvFunction(const VectorR3 * vpgs3, const VectorR3 * vmgs3, const VectorR3 & vp,
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

void MLFMA::reportLevelInfo(Qostream & strm) const
{
    int old_flag = strm.flags();
    strm << HEADING "MLFMA Level Information\n" TRAILING
        << LEVEL1 "Number of level: " << boundary_.getLevel() + 1 << '\n'
        << LEVEL1 "Number of valid boxes:" << '\n';
    for (size_t cur_level = 1; cur_level < level_box_.size(); ++cur_level)
        strm << "    => level " << setw(2) << cur_level << setw(15) << ':' << level_box_[cur_level].size() << '\n';
    strm << std::flush;
    strm.flags(old_flag);
}

void MLFMA::reportLevelDetail() const
{
    Qofstream level_log(dir_ + "/level_detail_info.txt", std::ios::out);
    level_log << std::left;
    reportLevelInfo(level_log);
    level_log << HEADING "Detail of level box: \n" TRAILING;
    for (size_t cur_level = 1; cur_level < level_box_.size(); ++cur_level)
    {
        level_log << "  -> level " << cur_level << " :\n";
        int cur_idx = 0;
        for (const auto& cur_box : level_box_[cur_level])
        {
            level_log << "    => box " << setw(8) << ++cur_idx << ":    child: " << cur_box.child.size()
                << "    near: " << setw(3) << cur_box.near.size()
                << "    far: " << setw(3) << cur_box.far.size()
                << "    unknown: " << cur_box.unknown.size() << '\n';
        }
    }
    level_log.close();
}

void MLFMA::reportTaskAssigmentInfo(Qostream & strm) const
{
    int old_state = strm.flags();
    size_t task_num = assignment_.n_rows - 1;
    strm << HEADING "Task Assignment Information\n" TRAILING
        << LEVEL1 "Task packet number: " << task_num << '\n'
        << LEVEL1 "Task distribution:" << "Boxes in per task packet\n";
    for (size_t cur_level = 2; cur_level <= level_; ++cur_level)
        strm << "    => level " << setw(2) << cur_level << setw(15) << ':'
        << "max: " << setw(3) << assignment_(1, cur_level - 2) - assignment_(0, cur_level - 2)
        << " min: " << setw(3) << assignment_(task_num, cur_level - 2) - assignment_(task_num - 1, cur_level - 2)
        << " total: " << assignment_(task_num, cur_level - 2) << '\n';
    strm.flush();
    strm.flags(old_state);
}

void MLFMA::reportMLFMAInfo(Qostream & strm) const
{
    std::vector<int> level_pole(level_ + 1, 0);
    for (int cur_level = 2; cur_level <= level_; ++cur_level)
        level_pole[cur_level] = getMultipoleNumber(box_len_[cur_level]);

    size_t ip_ap_elems = 0;
    for (int cur_level = 2; cur_level < level_; ++cur_level)
    {
        int ck_num = 2 * level_pole[cur_level] * level_pole[cur_level];
        int nk_num = 2 * level_pole[cur_level + 1] * level_pole[cur_level + 1];
        ip_ap_elems += ck_num + nk_num;
    }
    ip_ap_elems *= 9;

    size_t trans_elems = 0;
    for (int cur_level = 2; cur_level <= level_; ++cur_level)
    {
        int ck_num = 2 * level_pole[cur_level] * level_pole[cur_level];
        size_t far_box_num = 0;
        for (const auto& cur_box : level_box_[cur_level])
            far_box_num += cur_box.far.size();
        trans_elems += ck_num * far_box_num;
    }
    
    auto ip_ap_mem = ip_ap_elems * sizeof(value_t) / 1024 / 1024;
    auto ip_ap_aux_mem = ip_ap_elems * 2 * sizeof(size_t) / 1024 / 1024;
    auto near_mem = unk_near_ * 2 * sizeof(value_t) / 1024 / 1024;
    auto near_aux_mem = unk_near_ * 2 * sizeof(size_t) / 1024 / 1024;
    auto trans_mem = trans_elems * 2 * sizeof(value_t) / 1024 / 1024;
    auto rdcv_mem = 24 * sizeof(value_t) * unknowns_ * level_pole[level_] * level_pole[level_] / 1024 / 1024;
    auto precond_mem = precond_num_ * 2 * sizeof(value_t) / 1024 / 1024;
    auto precond_aux_mem = precond_num_ * 2 * sizeof(size_t) / 1024 / 1024;
    auto total_mem = ip_ap_aux_mem + ip_ap_aux_mem + near_mem + near_aux_mem + trans_mem + rdcv_mem + precond_mem + precond_aux_mem;
    auto mom_mem = 2 * sizeof(value_t) * unknowns_ * unknowns_ / 1024 / 1024;
    bool overflow = mom_mem < total_mem;
    auto save_mem = overflow ? total_mem - mom_mem : mom_mem - total_mem;

    int old_flag = strm.flags();
    strm << HEADING "MLFMA Information\n" TRAILING;
	strm << LEVEL1 "Number of multipoles:" << '\n';
    for (int cur_level = 2; cur_level <= level_; ++cur_level)
        strm << "    => level " << setw(2) << cur_level << setw(15) << ':' << level_pole[cur_level] << '\n';
    strm << LEVEL1 "Average per box:" << unknowns_ / level_box_.back().size() << '\n'
        << LEVEL1 "Number of near:" << unk_near_ << '\n';
    if (preconditioning_)
        strm << LEVEL1 "Number of precondition:" << precond_num_ << '\n';
    strm << LEVEL1 "Sparse rate:" << unk_near_ * 100.0f / unknowns_ / unknowns_ << "%\n"
        << LEVEL1 "Interpolation:" << ip_ap_mem << " MB" << FORMAT_MEMORY(ip_ap_mem) << '\n'
        << LEVEL2 "Auxiliary data:" << ip_ap_aux_mem << " MB" << FORMAT_MEMORY(ip_ap_aux_mem) << '\n'
        << LEVEL1 "Near impedance:" << near_mem << " MB" << FORMAT_MEMORY(near_mem) << '\n'
        << LEVEL2 "Auxiliary data:" << near_aux_mem << " MB" << FORMAT_MEMORY(near_aux_mem) << '\n'
        << LEVEL1 "Transfer:" << trans_mem << " MB" << FORMAT_MEMORY(trans_mem) << '\n'
        << LEVEL1 "Rad and Recv:" << rdcv_mem << " MB" << FORMAT_MEMORY(rdcv_mem) << '\n';
    if (preconditioning_)
        strm << LEVEL1 "Preconditioner:" << precond_mem << " MB" << FORMAT_MEMORY(precond_mem) << '\n'
            << LEVEL2 "Auxiliary data:" << precond_aux_mem << " MB" << FORMAT_MEMORY(precond_aux_mem) << '\n';
    strm << LEVEL1 "Total memory:" << total_mem << " MB" << FORMAT_MEMORY(total_mem) << '\n'
        << LEVEL1 "Saved memory: " << (overflow ? "-" : "") << save_mem << " MB" << FORMAT_MEMORY(save_mem) << '\n';
    strm.flush();
    strm.flags(old_flag);
}

void MLFMA::fillVbyBox()
{
    V_.set_size(unknowns_);
    size_t curpos = 0;

    size_t cur_bp = 0, total_bp = level_box_.back().size(); //
    tool::BarAndPercent bar_perc;   //
    for (const auto& cur_box : level_box_.back())
    {
        bar_perc(++cur_bp, total_bp);  //

        auto& edges = cur_box.unknown;
        assert(curpos == cur_box.offset);

        for (int eid : edges)
        {
            int nvp, nvm, fp, fm;
            VectorR3 vp, vm, vp3[3], vm3[3];
            value_t length;

            ce_ptr_->getCommonEdge(eid, nvp, nvm, fp, fm, length);
            vp = mesh_ptr_->getVertex(nvp);
            vm = mesh_ptr_->getVertex(nvm);
            const auto& tri_plu = mesh_ptr_->getTriangleRef(fp);
            const auto& tri_min = mesh_ptr_->getTriangleRef(fm);
            tri_plu.getVertex(vp3[0], vp3[1], vp3[2]);
            tri_min.getVertex(vm3[0], vm3[1], vm3[2]);

            VectorR3 nml_plu = tri_plu.getNormal();
            VectorR3 nml_min = tri_min.getNormal();

            V_(curpos++) = 0.5f * length * cVKernel(vp3, vm3, nml_plu, nml_min, vp, vm);
        }
    }
    assert(curpos == unknowns_);
}

value_t MLFMA::getBiRCS(const VectorR3 & sca_k) const
{
    size_t curpos = 0;

    int nvp, nvm, fp, fm;
    VectorR3 vp, vm, vp3[3], vm3[3];
    value_t length;
    VectorR3 vgp[7], vgm[7];
    VectorR3 rou_p, rou_m;
    Complex G0(0, 0);
    VectorC3 Es(0, 0, 0);

    for (const auto& cur_box : level_box_.back())
    {
        auto& edges = cur_box.unknown;
        for (int u : edges)
        {
            ce_ptr_->getCommonEdge(u, nvp, nvm, fp, fm, length);
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
            Es = Es + (Esp + Esm) * (I_(curpos++) * length);
        }
    }
    assert(curpos == unknowns_);

    auto es = 0.5f * sca_k * (sca_k * Es);
    value_t rcs = Z0 * Z0 * k_ * k_ * es.norm() / PI4;
    return 10 * log10(rcs);
}

void MLFMA::getEFiled(component::FieldData * pdata, const VectorR3 & rad_k, const VectorR3 & rad_ev, const VectorR3 & rad_eh) const
{

}

void MLFMA::getNearEField(std::vector<component::NearFieldData>* data) const
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
	size_t curpos = 0;
	for (int p = 0; p < point_num; p++)
	{
		VectorC3 E_vol, E_surf, e_vol_temp, e_surf_temp;
		VectorR3 ob_point = pointArr[p];
		
		/*for (int e = 0; e < unknowns_e; e++)
		{
			NEFkernel_surf(ob_point, e, e_surf_temp);
			E_surf += I(e + unknowns_t)*e_surf_temp;
		}*/
		for (const auto& cur_box : level_box_.back())
		{
			auto& edges = cur_box.unknown;
			for (int u : edges)
			{
				NEFkernel_surf(ob_point, u, e_surf_temp);
				E_surf += I_(curpos++)*e_surf_temp;
			}
		}

		E_near = coff * (E_vol + E_surf);
		NearFieldData EFdata(ob_point, E_near);
		data->push_back(EFdata);
	}
}

void MLFMA::NEFkernel_surf(VectorR3& ob, int& unk, VectorC3& e_surf) const
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

	e_surf = (len*((E1 / 2.0f) + E2 / (k_*k_)));
}

void MLFMA::calculateSurfaceCurrent(std::vector<component::CurrentData>* currents) const
{
    Qcx_vec cur_vec(unknowns_, arma::fill::zeros);
    size_t cur_idx = 0;
    for (const auto& box : level_box_.back())
        for (int u : box.unknown)
            cur_vec(u) = I_(cur_idx++);
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

///////////////////////////////////////////////////////////////////////////////////////////////////////

Complex MLFMA::cZppKernel(const VectorR3 * vf3, const VectorR3 * vs3, const VectorR3 & vf, const VectorR3 & vs, const VectorR3 & nf) const
{
    return integral::cfieZKernel(alpha_, k_, vf3, vs3, vf, vs, nf, 1, 1);
}

Complex MLFMA::cZpmKernel(const VectorR3 * vf3, const VectorR3 * vs3, const VectorR3 & vf, const VectorR3 & vs, const VectorR3 & nf) const
{
    return integral::cfieZKernel(alpha_, k_, vf3, vs3, vf, vs, nf, 1, -1);
}

Complex MLFMA::cZmpKernel(const VectorR3 * vf3, const VectorR3 * vs3, const VectorR3 & vf, const VectorR3 & vs, const VectorR3 & nf) const
{
    return integral::cfieZKernel(alpha_, k_, vf3, vs3, vf, vs, nf, -1, 1);
}

Complex MLFMA::cZmmKernel(const VectorR3 * vf3, const VectorR3 * vs3, const VectorR3 & vf, const VectorR3 & vs, const VectorR3 & nf) const
{
    return integral::cfieZKernel(alpha_, k_, vf3, vs3, vf, vs, nf, -1, -1);
}

Complex MLFMA::cZppSingular(const VectorR3 * vf3, const VectorR3 * vs3, const VectorR3 & vf, const VectorR3 & vs) const
{
    return integral::cfieZSingular(alpha_, k_, vf3, vs3, vf, vs, 1, 1);
}

Complex MLFMA::cZpmSingular(const VectorR3 * vf3, const VectorR3 * vs3, const VectorR3 & vf, const VectorR3 & vs) const
{
    return integral::cfieZSingular(alpha_, k_, vf3, vs3, vf, vs, 1, -1);
}

Complex MLFMA::cZmpSingular(const VectorR3 * vf3, const VectorR3 * vs3, const VectorR3 & vf, const VectorR3 & vs) const
{
    return integral::cfieZSingular(alpha_, k_, vf3, vs3, vf, vs, -1, 1);
}

Complex MLFMA::cZmmSingular(const VectorR3 * vf3, const VectorR3 * vs3, const VectorR3 & vf, const VectorR3 & vs) const
{
    return integral::cfieZSingular(alpha_, k_, vf3, vs3, vf, vs, -1, -1);
}

Complex MLFMA::cVKernel(const VectorR3 * vp3, const VectorR3 * vm3, const VectorR3 & np, const VectorR3 & nm, const VectorR3 & vp, const VectorR3 & vm) const
{
    return integral::cfieVKernel(alpha_, k_, inc_k_, inc_e_, inc_h_, vp3, vm3, np, nm, vp, vm);
}
