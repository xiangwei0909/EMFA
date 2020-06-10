#include "FastWT.h"
#include "Miscellaneous.h"
#include "ConfigLoader.h"
#include "Mathdef.h"
#include "Quadrature.h"
#include "tools.h"
#include "Mesh.h"
#include "ContourEdge.h"

using namespace component;
using namespace mom;
using namespace math;

using std::setw;


FastWT::FastWT()
{
}


FastWT::~FastWT()
{
}

void FastWT::init(component::ConfigLoader * ploader)
{
    EM::init(ploader);
    alpha_ = ploader->getAlpha();
    if(rtype_ == policy::ResultType::Rad)
    {
        //readExcEdges(dir_ + '/' + folder_name_ + ".rad");
        alpha_ = 1;
    }
    dir_ = tool::creatFolder(dir_, folder_name_ + "_FastWT_" + std::to_string(alpha_).substr(0, 3));

    mesh_ptr_ = std::make_shared<Mesh>();
    ce_ptr_ = std::make_shared<ContourEdge>();

    mesh_ptr_->loadMeshFile(ploader->getMeshPath());
    ce_ptr_->buildEdge(mesh_ptr_, policy::HALF_RWG);
    unknowns_ = ce_ptr_->getUnknownsNum();

    iteration_threshold_ = ploader->getIterationThreshold();
    max_iteration_num_ = ploader->getMaxIterationNum();
    preconditioning_ = ploader->getPreconditioningSwitch();
    row_threshold_ = ploader->getRowThreshold();
    col_threshold_ = ploader->getColThreshold();
    beta_ = ploader->getStabilizationFactor() * abs(log10(ploader->getAverageSize()));

    auto wave_length = cc / incidence_.freq;
    k_ = PI2 / wave_length;

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

void FastWT::solve()
{
    SEGMENT("Solve");

    levelRootsAndWeights();
    calculateInterpolationMtr();
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

void FastWT::output()
{
    SEGMENT("Result");
    Result result;
    if(rtype_ == policy::ResultType::Rad)
    {
        //LOG(result.getRadiationField(this, &FastWT::getEFiled), "Calculating Radiation");
        TIME_LOG("getRadiation");
    }
    else
    {
        LOG(result.getBistaticRCS(this, &FastWT::getBiRCS), "Calculating BistaticRCS");
        TIME_LOG("getBistaticRCS");
    }
    LOG(result.getCurrentDistribution(this, &FastWT::calculateSurfaceCurrent), "Calculating surface current");
    TIME_LOG("getSurfaceCurrent");
}

void FastWT::clear()
{
    pool_.clear();
    assignment_.reset();
    nearZ_.reset();
    precond_.reset();
    V_.reset();
    I_.reset();
    coeff_.reset();
    pcoeff_.reset();
    rad_.clear();
    rec_.clear();
    prad_.clear();
    prec_.clear();
    level_rad_.clear();
    level_rec_.clear();
    plevel_rad_.clear();
    plevel_rec_.clear();
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

void FastWT::reportInfo(Qostream & strm) const
{
    mesh_ptr_->reportInfo(strm);
    ce_ptr_->reportInfo(strm);
    boundary_.reportInfo(strm);
    pool_.reportInfo(strm);
    this->reportLevelInfo(strm);
    this->reportTaskAssigmentInfo(strm);
    this->reportFastWTInfo(strm);
}

void FastWT::buildBox()
{
    std::vector<VectorR3> center_points(unknowns_);
    fillCenterPoints(center_points);

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
    precond_num_ = 0;
    if (preconditioning_)
    {
        for (const auto& cur_box : level_box_.back())
        {
            size_t cols = 0;
            for (const auto nptr : cur_box.near)
            {
                auto& unks = nptr->unknown;
                for (size_t i = 0; i < unks.size(); ++i)
                {
                    auto dist = (center_points[unks[i]] - cur_box.center).Norm();
                    if (dist <= col_threshold_) ++cols;
                }
            }
            precond_num_ += cols * cur_box.unknown.size();
        }
    }
}

void FastWT::fillCenterPoints(std::vector<VectorR3>& centers)
{
    for (int u = 0; u < unknowns_; ++u)
    {
        auto& hrwg = ce_ptr_->getHRWGRef(u);
        auto v1 = mesh_ptr_->getVertex(hrwg.v1);
        auto v2 = mesh_ptr_->getVertex(hrwg.v2);
        centers[u] = (v1 + v2) / 2;
    }
}

void FastWT::buildFirstLevel(const std::vector<VectorR3>& cps, value_t box_len)
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

void FastWT::buildCurLevel(const std::vector<VectorR3>& cps, int cur_level, value_t box_len)
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

void FastWT::setBoxParameter()
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

int FastWT::getMultipoleNumber(value_t box_len) const
{
    value_t kd = 1.732051f * k_ * box_len;
    return static_cast<int>(kd + 6 * cbrtf(kd));
}

void FastWT::prepareAssignment(size_t thread_num)
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

void FastWT::levelRootsAndWeights()
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

void FastWT::calculateInterpolationMtr()
{
    arr_ip_.resize(level_ - 2);

    Qumat tloc, ploc;
    Qmat telems, pelems;
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

void FastWT::lagrangeIPWeights(const Qvec & nxt_sam, const Qvec & cur_sam, Qumat & loc, Qmat & elems) const
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

void FastWT::calculateLevelK()
{
    level_k_.resize(level_ - 1);
    Qvec tsin, tcos, psin, pcos;
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

void FastWT::calculateLevelTrans()
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
                auto task = std::bind(&FastWT::parallelLevelTrans, this, begin, end, idx);
                pool_.submit(task);
            }
        }
        pool_.run();
    } 
}

void FastWT::parallelLevelTrans(size_t begin, size_t end, size_t idx)
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

void FastWT::calculateRadAndRecv()
{
    size_t task_num = assignment_.n_rows - 1;
    size_t idx = level_ - 2;
    size_t box_num = assignment_(task_num, idx);
    rad_.resize(box_num);
    rec_.resize(box_num);
    prad_.resize(box_num);
    prec_.resize(box_num);
    for (size_t t = 0; t < task_num; ++t)
    {
        size_t begin = assignment_(t, idx);
        size_t end = assignment_(t + 1, idx);
        if (begin < end)
        {
            auto task = std::bind(&FastWT::parallelRadAndRecv, this, begin, end);
            pool_.submit(task);
        }
    }
    pool_.interactiveRun();
}

void FastWT::parallelRadAndRecv(size_t begin, size_t end)
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
        auto& prad_mat = prad_[ib];
        auto& prec_mat = prec_[ib];
        rad_mat.set_size(rwg_num, knum);    // RWG x K
        rec_mat.set_size(knum, rwg_num);    // K x RWG
        prad_mat.set_size(knum, rwg_num);   // K x RWG
        prec_mat.set_size(knum, rwg_num);   // K x RWG
        VectorR3 vx3[3], vg3[3], vp, nml, v1, v2, vgl3[3];
        for (size_t i = 0; i < rwg_num; ++i)
        {
            auto& hrwg = ce_ptr_->getHRWGRef(edges[i]);
            v1 = mesh_ptr_->getVertex(hrwg.v1);
            v2 = mesh_ptr_->getVertex(hrwg.v2);
            vp = mesh_ptr_->getVertex(hrwg.vx);
            auto& triangle = mesh_ptr_->getTriangleRef(hrwg.tid);
            triangle.getVertex(vx3[0], vx3[1], vx3[2]);
            nml = triangle.getNormal();
            Gauss3Point(vx3[0], vx3[1], vx3[2], vg3);
            GaussLine3Point(v1, v2, vgl3);
            for (size_t k = 0; k < knum; ++k)
            {
                const auto& sam_k = arr_k(k);
                auto tmp_rad = radiationFunction(vg3, vp, center, sam_k);
                auto tmp_rec = halfRecvFunction(vg3, vp, nml, center, sam_k);
                rad_mat(i, k) = (hrwg.length / 2) * tmp_rad;
                rec_mat(k, i) = (hrwg.length / 2) * ((alpha_ * conj(tmp_rad)) + ((1 - alpha_) * tmp_rec));
                auto line_rad = lineRadiationFunction(vgl3, vp, center, sam_k);
                prad_mat(k, i) = hrwg.length * line_rad;
                prec_mat(k, i) = hrwg.length * conj(line_rad);
            }
        } // end for rwg_num
    } // end for box_num
}

VectorC3 FastWT::radiationFunction(const VectorR3 * vpgs3, const VectorR3 & vp, const VectorR3 & center, const VectorR3 & sam_k) const
{
    VectorC3 val(0, 0, 0);
    for (int p = 0; p < 3; ++p)
        val += w3[p] * (vpgs3[p] - vp) * std::exp(J0 * k_ * (sam_k ^ (vpgs3[p] - center)));
    return val - ((sam_k ^ val) * sam_k);
}

VectorC3 FastWT::halfRecvFunction(const VectorR3 * vpgs3, const VectorR3 & vp, const VectorR3 & nml_plu, const VectorR3 & center, const VectorR3 & sam_k) const
{
    VectorC3 val(0, 0, 0);
    for (int p = 0; p < 3; ++p)
        val += w3[p] * (nml_plu * (vpgs3[p] - vp)) * std::exp(-J0 * k_ * (sam_k ^ (vpgs3[p] - center)));
    return sam_k * val;
}

Complex FastWT::lineRadiationFunction(const VectorR3 * vpgs3, const VectorR3 & vp, const VectorR3 & center, const VectorR3 & sam_k) const
{
    Complex val(0, 0);
    for (int p = 0; p < 3; ++p)
        val += glw3[p] * std::exp(J0 * k_ * (sam_k ^ (vpgs3[p] - center)));
    return val;
}

void FastWT::calculateNearZ()
{
    Qumat   location(2, unk_near_);
    Qcx_vec elements(unk_near_);
    const auto& box_array = level_box_.back();
    size_t  cur_base = 0;
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
            auto task = std::bind(&FastWT::parallelNearZ, this, begin, end, &location, &elements, cur_base);
            pool_.submit(task, task_unk);
            cur_base += task_unk;
        }
    }
    assert(cur_base == unk_near_);
    pool_.interactiveRun();
    nearZ_ = Qsp_cx_mat(location, elements, false);
}

void FastWT::parallelNearZ(size_t begin, size_t end, Qumat * ploc, Qcx_vec * pelems, size_t ebase) const
{
    for (size_t sb = begin; sb < end; ++sb)
    {
        auto& cur_box = level_box_.back()[sb];
        auto& srcs = cur_box.unknown;
        for (size_t s = 0; s < srcs.size(); ++s)
        {
            size_t spos = cur_box.offset + s;
            for (const auto& nptr : cur_box.near)
                fillNearZbyColumn(nptr->unknown, srcs[s], nptr->offset, spos, *ploc, *pelems, ebase);
        }
    }
}

void FastWT::fillNearZbyColumn(const std::vector<int>& flds, int sid, size_t fpos, size_t spos, Qumat & loc, Qcx_vec & elems, size_t & aux) const
{
    VectorR3 vs3[3], vf3[3], vfx, vsx, vf1, vf2, vs1, vs2;
    auto& shrwg = ce_ptr_->getHRWGRef(sid);
    auto& striangle = mesh_ptr_->getTriangleRef(shrwg.tid);
    vs1 = mesh_ptr_->getVertex(shrwg.v1);
    vs2 = mesh_ptr_->getVertex(shrwg.v2);
    vsx = mesh_ptr_->getVertex(shrwg.vx);
    striangle.getVertex(vs3[0], vs3[1], vs3[2]);

    for (size_t f = 0; f < flds.size(); ++f)
    {
        auto& fhrwg = ce_ptr_->getHRWGRef(flds[f]);
        auto& ftriangle = mesh_ptr_->getTriangleRef(fhrwg.tid);
        vf1 = mesh_ptr_->getVertex(fhrwg.v1);
        vf2 = mesh_ptr_->getVertex(fhrwg.v2);
        vfx = mesh_ptr_->getVertex(fhrwg.vx);
        ftriangle.getVertex(vf3[0], vf3[1], vf3[2]);
        auto nml = ftriangle.getNormal();
        Complex z_fs(0, 0), ip(0, 0), line_term(0, 0);
        if (shrwg.tid == fhrwg.tid)
        {
            z_fs = shrwg.length * fhrwg.length * cZSingular(vf3, vs3, vfx, vsx);
            if (sid == flds[f])
                ip = stabilizationTerm(fhrwg.length);
        }
        else
        {
            z_fs = shrwg.length * fhrwg.length * cZKernel(vf3, vs3, vfx, vsx, nml);
            value_t line_len;
            if (shrwg.type == ContourEdge::EdgeType::CONFORMAL &&
                fhrwg.type == ContourEdge::EdgeType::CONFORMAL &&
                shrwg.idx == fhrwg.idx)
            {
                ip = stabilizationTerm(fhrwg.length);
            }
            else if (shrwg.type == ContourEdge::EdgeType::NONCONFORMAL &&
                fhrwg.type == ContourEdge::EdgeType::NONCONFORMAL &&
                ce_ptr_->getContourLineLength(fhrwg.idx, shrwg.tid, &line_len))
            {
                ip = stabilizationTerm(line_len);
            }
        }
        line_term = -(fhrwg.length * consistencyTerm(vf3, vs1, vs2) + shrwg.length * consistencyTerm(vs3, vf1, vf2));
        //  inserted by column-major order
        loc(0, aux) = fpos + f;
        loc(1, aux) = spos;
        elems(aux) = z_fs + alpha_ * (line_term + ip);
        ++aux;
    }
}

void FastWT::fillVbyBox()
{
    V_.set_size(unknowns_);
    size_t curpos = 0;
    VectorR3 vx3[3], vp, nml;

    TOOL size_t cur_bp = 0, total_bp = level_box_.back().size(); //
    TOOL tool::BarAndPercent bar_perc;   //
    for (const auto& cur_box : level_box_.back())
    {
        TOOL bar_perc(++cur_bp, total_bp);  //

        auto& edges = cur_box.unknown;
        assert(curpos == cur_box.offset);

        for (int fld : edges)
        {
            auto& hrwg = ce_ptr_->getHRWGRef(fld);
            auto& triangle = mesh_ptr_->getTriangleRef(hrwg.tid);
            triangle.getVertex(vx3[0], vx3[1], vx3[2]);
            vp = mesh_ptr_->getVertex(hrwg.vx);
            nml = triangle.getNormal();
            V_(curpos++) = 0.5f * hrwg.length * cVKernel(vx3, vp, nml);
        }
    }
    assert(curpos == unknowns_);
}

void FastWT::calculatePreconditioner()
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
            auto task = std::bind(&FastWT::parallelPreconditioner, this, begin, end, std::cref(center_points), &location, &elements, base[t]);
            pool_.submit(task, base[t + 1] - base[t]);
        }
    }
    pool_.run();
    precond_ = Qsp_cx_mat(location, elements, false);
}

void FastWT::preparepreconditioner(const std::vector<VectorR3>& centers, std::vector<size_t>& base) const
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

void FastWT::parallelPreconditioner(size_t begin, size_t end, const std::vector<VectorR3> &centers, Qumat * ploc, Qcx_vec * pelems, size_t ebase) const
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
            // thread-unsafe function
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

int FastWT::iterativeSolve()
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

void FastWT::allocateForSolve()
{
    level_rad_.resize(level_ - 1);
    level_rec_.resize(level_ - 1);
    plevel_rad_.resize(level_ - 1);
    plevel_rec_.resize(level_ - 1);
    for (int cur_level = 2; cur_level <= level_; ++cur_level)
    {
        int idx = cur_level - 2;
        auto& up_mat = level_rad_[idx];
        auto& down_mat = level_rec_[idx];
        auto knum = level_k_[idx].size();
        auto bnum = level_box_[cur_level].size();
        up_mat.set_size(knum, bnum);
        down_mat.set_size(knum, bnum);
        plevel_rad_[idx].set_size(knum, bnum);
        plevel_rec_[idx].set_size(knum, bnum);
    }
}

void FastWT::calculateCoeff()
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
    pcoeff_ = (alpha_ * Z0 / 32) * coeff_;
    coeff_ = (k_ * k_ * Z0 / 32) * coeff_;
}

Qcx_vec FastWT::matrixVectorMultiply(const Qcx_vec & b)
{
    for (int cur_level = 2; cur_level <= level_; ++cur_level)
    {
        level_rad_[cur_level - 2].reset(VectorC3(0, 0, 0));
        plevel_rad_[cur_level - 2].fill(Complex(0, 0));
    }
    const size_t task_num = assignment_.n_rows - 1;
    size_t idx = level_ - 2;
    for (size_t t = 0; t < task_num; ++t)
    {
        size_t begin = assignment_(t, idx);
        size_t end = assignment_(t + 1, idx);
        if (begin < end)
        {
            auto task = std::bind(&FastWT::parallelFinestLayerUpward, this, begin, end, std::cref(b));
            pool_.submit(task);
        }
    }
    pool_.run();

    upwardPass();

    for (int cur_level = 2; cur_level <= level_; ++cur_level)
    {
        level_rec_[cur_level - 2].reset(VectorC3(0, 0, 0));
        plevel_rec_[cur_level - 2].fill(Complex(0, 0));
    }
    idx = 0;
    for (size_t t = 0; t < task_num; ++t)
    {
        size_t begin = assignment_(t, idx);
        size_t end = assignment_(t + 1, idx);
        if (begin < end)
        {
            auto task = std::bind(&FastWT::parallelFinestLayerDownward, this, begin, end);
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
            auto task = std::bind(&FastWT::parallelFar, this, begin, end, &far);
            pool_.submit(task);
        }
    }
    pool_.run();

    if (preconditioning_)
        return precond_ * (nearZ_ * b + far);
    return nearZ_ * b + far;
}

void FastWT::parallelFinestLayerUpward(size_t begin, size_t end, const Qcx_vec & b)
{
    auto& up_mat = level_rad_.back();
    auto& pup_mat = plevel_rad_.back();

    size_t knum = up_mat.row();
    for (size_t ib = begin; ib < end; ++ib)
    {
        auto& cur_box = level_box_.back()[ib];
        size_t rwg_num = cur_box.unknown.size();
        size_t base = cur_box.offset;
        size_t bidx = cur_box.index;
        const auto& rad_mat = rad_[bidx];
        const auto& prad_mat = prad_[bidx];
        for (size_t k = 0; k < knum; ++k)
        {
            VectorC3 val(0, 0, 0);
            for (size_t u = 0; u < rwg_num; ++u)
                val += b(base + u) * rad_mat(u, k);
            up_mat(k, bidx) = val;
        }
        pup_mat.col(bidx) = prad_mat * b.subvec(base, arma::size(rwg_num, 1));
    }
}

void FastWT::parallelFinestLayerDownward(size_t begin, size_t end)
{
    auto& down_mat = level_rec_.front();
    auto& pdown_mat = plevel_rec_.front();
    const auto& u_mat = level_rad_.front();
    const auto& pu_mat = plevel_rad_.front();
    const auto& arr_trans = level_trans_.front();
    const auto& box_array = level_box_[2];
    const size_t knum = u_mat.row();
    for (size_t ib = begin; ib < end; ++ib)
    {
        size_t bidx = box_array[ib].index;
        const auto& trans = arr_trans[bidx];
        size_t aux = 0, cur_col = 0;
        for (const auto fptr : box_array[ib].far)
        {
            size_t fidx = fptr->index;
            for (size_t k = 0; k < knum; ++k)
                down_mat(k, bidx) += trans(aux++) * u_mat(k, fidx);
            pdown_mat.col(bidx) += trans.col(cur_col++) % pu_mat.col(fidx);
        }
        assert(aux == trans.n_elem);
        assert(cur_col == trans.n_cols);
    }
}

void FastWT::parallelFar(size_t begin, size_t end, Qcx_vec * pfar) const
{
    const auto& d_mat = level_rec_.back();
    const auto& pd_mat = plevel_rec_.back();
    const auto& box_array = level_box_.back();
    const size_t knum = d_mat.row();
    for(size_t ib = begin; ib < end; ++ib)
    {
        const size_t rwg_num = box_array[ib].unknown.size();
        const size_t bidx = box_array[ib].index;
        const size_t offset = box_array[ib].offset;
        const auto& rec_mat = rec_[bidx];
        const auto& prec_mat = prec_[bidx];
        for (size_t u = 0; u < rwg_num; ++u)
        {
            Complex val(0, 0);
            for (size_t k = 0; k < knum; ++k)
                val += coeff_(k) * (rec_mat(k, u) ^ d_mat(k, bidx));
            val += arma::dot(pcoeff_ % prec_mat.col(u), pd_mat.col(bidx));
            (*pfar)(offset + u) = val;
        }
    }
}

void FastWT::upwardPass()
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
                auto task = std::bind(&FastWT::parallelUpwardPass, this, begin, end, cur_level);
                pool_.submit(task);
            }
        }
        pool_.run();
    } // previous level
}

void FastWT::parallelUpwardPass(size_t begin, size_t end, int cur_level)
{
    size_t idx = cur_level - 2;
    const auto& arr_cbox = level_box_[cur_level + 1];
    const auto& ip_mat = arr_ip_[idx];
    const auto& in_mat = level_rad_[idx + 1];
    const auto& pin_mat = plevel_rad_[idx + 1];
    const auto& arr_k = level_k_[idx];
    const auto knum = arr_k.size();
    wt::QMat<VectorC3> out_mat;
    Qcx_mat pout_mat;
    out_mat.set_size(level_rad_[idx].row(), level_rad_[idx].col(), VectorC3());
    pout_mat.zeros(arma::size(plevel_rad_[idx]));
    for(size_t cidx = begin; cidx < end; ++cidx)
    {
        VectorR3 rvec = arr_cbox[cidx].parent->center - arr_cbox[cidx].center;
        size_t bid = arr_cbox[cidx].parent->index;
        size_t cur_num = 0;
        for (size_t k = 0; k < knum; ++k)
        {
            const Complex phase = exp(-J0 * k_ * (arr_k(k) ^ rvec));
            VectorC3 val(0, 0, 0);
            Complex line_val(0, 0);
            for (size_t p = 0; p < 9; ++p)
            {
                val += ip_mat(cur_num) * in_mat(ip_mat.col(cur_num), cidx);
                line_val += ip_mat(cur_num) * pin_mat(ip_mat.col(cur_num), cidx);
                ++cur_num;
            }
            out_mat(k, bid) += phase * val;
            pout_mat(k, bid) += phase * line_val;
        }
        assert(cur_num == 9 * out_mat.row());
    }
    MutexGuard guard(mutex_);
    level_rad_[idx] += out_mat;
    plevel_rad_[idx] += pout_mat;
}

void FastWT::downwardPass()
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
                auto task = std::bind(&FastWT::parallelDownwardPass, this, begin, end, idx);
                pool_.submit(task);
            }
        }
        pool_.run();
    } // next level
}

void FastWT::parallelDownwardPass(size_t begin, size_t end, size_t idx)
{
    const auto& arr_trans = level_trans_[idx];
    const auto& up_mat = level_rad_[idx];
    const auto& pup_mat = plevel_rad_[idx];
    auto& down_mat = level_rec_[idx];
    auto& pdown_mat = plevel_rec_[idx];
    const auto& pdown = level_rec_[idx - 1];
    const auto& ppdown = plevel_rec_[idx - 1];
    const auto& ip_mat = arr_ip_[idx - 1];
    const auto& arr_k = level_k_[idx - 1];
    const auto& gw = glw_[idx];
    const auto& pgw = glw_[idx - 1];
    const auto& box_array = level_box_[idx + 2];
    const size_t knum = up_mat.row();
    for(size_t ib = begin; ib < end; ++ib)
    {
        auto& cur_box = box_array[ib];
        size_t bidx = cur_box.index;
        const auto& trans = arr_trans[bidx];
        size_t aux = 0, cur_col = 0;
        for (const auto fptr : cur_box.far)
        {
            size_t fidx = fptr->index;
            for (size_t k = 0; k < knum; ++k)
                down_mat(k, bidx) += trans(aux++) * up_mat(k, fidx);
            pdown_mat.col(bidx) += trans.col(cur_col++) % pup_mat.col(fidx);
        }
        assert(aux == trans.n_elem);
        assert(cur_col == trans.n_cols);

        size_t pidx = cur_box.parent->index;
        VectorR3 rvec = cur_box.parent->center - cur_box.center;
        for (size_t cur_num = 0; cur_num < ip_mat.size(); ++cur_num)
        {
            const size_t ck = ip_mat.col(cur_num);
            const size_t pk = ip_mat.row(cur_num);
            const value_t coeff = ip_mat(cur_num) * pgw(pk) / gw(ck);
            const Complex phase = exp(J0 * k_ * (arr_k(pk) ^ rvec));
            down_mat(ck, bidx) += coeff * phase * pdown(pk, pidx);
            pdown_mat(ck, bidx) += coeff * phase * ppdown(pk, pidx);
        }
    }
}

void FastWT::reportLevelInfo(Qostream & strm) const
{
    int old_flag = strm.flags();
    strm << HEADING "MLFMA Level Information\n" TRAILING
        << LEVEL1 "Number of level: " << boundary_.getLevel() + 1 << '\n'
        << LEVEL1 "Number of valid boxes:" << '\n';
    for (size_t cur_level = 1; cur_level < level_box_.size(); ++cur_level)
        strm << "    => level " << setw(2) << cur_level << setw(15) << ':' << level_box_[cur_level].size() << '\n';
    strm.flush();
    strm.flags(old_flag);
}

void FastWT::reportTaskAssigmentInfo(Qostream & strm) const
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

void FastWT::reportLevelDetail() const
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

void FastWT::reportFastWTInfo(Qostream & strm) const
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
    auto prdcv_mem = rdcv_mem / 3;
    auto precond_mem = precond_num_ * 2 * sizeof(value_t) / 1024 / 1024;
    auto precond_aux_mem = precond_num_ * 2 * sizeof(size_t) / 1024 / 1024;
    auto total_mem = ip_ap_aux_mem + ip_ap_aux_mem + near_mem + near_aux_mem + trans_mem + rdcv_mem + prdcv_mem + precond_mem + precond_aux_mem;
    auto mom_mem = 2 * sizeof(value_t) * unknowns_ * unknowns_ / 1024 / 1024;
    bool overflow = mom_mem < total_mem;
    auto save_mem = overflow ? total_mem - mom_mem : mom_mem - total_mem;

    int old_flag = strm.flags();
    strm << HEADING "FastWT Information\n" TRAILING;
    strm << LEVEL1 "Stabilization parameter:" << beta_ << '\n'
        << LEVEL1 "Number of multipoles:" << '\n';
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
        << LEVEL1 "Rad and Recv:" << rdcv_mem << " MB" << FORMAT_MEMORY(rdcv_mem) << '\n'
        << LEVEL1 "Interior penalty:" << prdcv_mem << " MB" << FORMAT_MEMORY(prdcv_mem) << '\n';
    if (preconditioning_)
    {
        strm << LEVEL1 "Preconditioner:" << precond_mem << " MB" << FORMAT_MEMORY(precond_mem) << '\n'
            << LEVEL2 "Auxiliary data:" << precond_aux_mem << " MB" << FORMAT_MEMORY(precond_aux_mem) << '\n';
    }
    strm << LEVEL1 "Total memory:" << total_mem << " MB" << FORMAT_MEMORY(total_mem) << '\n'
        << LEVEL1 "Saved memory:" << (overflow ? "-" : "") << save_mem << " MB" << FORMAT_MEMORY(save_mem) << '\n';
    strm.flush();
    strm.flags(old_flag);
}

value_t FastWT::getBiRCS(const VectorR3 & sca_k) const
{
    VectorR3 vx3[3], vp, vg3[3], rou;
    VectorC3 Es(0, 0, 0);
    int cur_row = 0;

    for (const auto& cur_box : level_box_.back())
    {
        auto& edges = cur_box.unknown;
        for (int u : edges)
        {
            auto& hrwg = ce_ptr_->getHRWGRef(u);
            mesh_ptr_->getTriangleRef(hrwg.tid).getVertex(vx3[0], vx3[1], vx3[2]);
            vp = mesh_ptr_->getVertex(hrwg.vx);
            Gauss3Point(vx3[0], vx3[1], vx3[2], vg3);

            VectorC3 u_es(0, 0, 0);
            for (int g = 0; g < 3; ++g)
            {
                rou = vg3[g] - vp;
                u_es = u_es + rou * w3[g] * exp(J0 * k_ * (vg3[g] ^ sca_k));
            }
            Es = Es + u_es * (I_(cur_row++) * hrwg.length);
        }
    }
    assert(cur_row == unknowns_);

    auto es = 0.5f * sca_k * (sca_k * Es);
    value_t rcs = Z0 * Z0 * k_ * k_ * es.norm() / PI4;
    return 10 * log10(rcs);
}

void FastWT::getEFiled(component::FieldData * pdata, const VectorR3 & rad_k, const VectorR3 & rad_ev, const VectorR3 & rad_eh) const
{
}

void FastWT::calculateSurfaceCurrent(std::vector<component::CurrentData>* currents) const
{
    Qcx_vec cur_vec(unknowns_, arma::fill::zeros);
    size_t cur_idx = 0;
    for (const auto& box : level_box_.back())
        for (int u : box.unknown)
            cur_vec(u) = I_(cur_idx++);
    size_t tri_num = mesh_ptr_->getTriangleNum();
    std::vector<std::vector<int>> triangles(tri_num);

    for (int u = 0; u < unknowns_; ++u)
    {
        auto& hrwg = ce_ptr_->getHRWGRef(u);
        triangles[hrwg.tid].push_back(u);
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
            auto& hrwg = ce_ptr_->getHRWGRef(unks[i]);
            vx[i] = mesh_ptr_->getVertex(hrwg.vx);
            len[i] = hrwg.length;
        }

        auto center = (vx[0] + vx[1] + vx[2]) / 3;
        auto dominator = 2 * Area(vx[0], vx[1], vx[2]);
        VectorC3 cen_cur(0, 0, 0);
        for (size_t i = 0; i < unks.size(); ++i)
        {
            cen_cur += cur_vec(unks[i]) * (len[i] / dominator) * (center - vx[i]);
            VectorC3 ncur(0, 0, 0);
            for (size_t j = 0; j < unks.size(); ++j)
            {
                if (i == j) continue;
                ncur += cur_vec(unks[j]) * (len[j] / dominator) * (vx[i] - vx[j]);
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

Complex FastWT::cZKernel(const VectorR3 * vf3, const VectorR3 * vs3, const VectorR3 & vf, const VectorR3 & vs, const VectorR3 & nml) const
{
    return integral::cfieZKernel(alpha_, k_, vf3, vs3, vf, vs, nml, 1, 1);
}

Complex FastWT::cZSingular(const VectorR3 * vf3, const VectorR3 * vs3, const VectorR3 & vf, const VectorR3 & vs) const
{
    return integral::cfieZSingular(alpha_, k_, vf3, vs3, vf, vs, 1, 1);
}

Complex FastWT::cVKernel(const VectorR3 * vx3, const VectorR3 & vp, const VectorR3 & nml) const
{
    VectorR3 vg[3], rou;
    Complex Ve, Vm, G;
    Gauss3Point(vx3[0], vx3[1], vx3[2], vg);
    VectorC3 Ei(0, 0, 0), Hi(0, 0, 0);
    for (int a = 0; a < 3; ++a)
    {
        rou = vg[a] - vp;
        G = exp(-J0 * k_ * (vg[a] ^ inc_k_));
        Ei = G * inc_e_;
        Ve += w3[a] * (rou ^ Ei);
        Hi = G * inc_h_;    //  include Z0
        Vm += w3[a] * (rou ^ (nml * Hi));
    }
    return alpha_ * Ve + (1 - alpha_) * Vm;   // Note: multiple (0.5 * length)
}

Complex FastWT::lineIntegral(const VectorR3 * vx3, const VectorR3 & v1, const VectorR3 & v2) const
{
    VectorR3 vg3[3];
    Gauss3Point(vx3[0], vx3[1], vx3[2], vg3);
    auto v1v2 = v2 - v1;
    Complex line_term(0, 0);
    for (int p = 0; p < 3; ++p)
    {
        auto vi = v1 + gl3[p] * v1v2;
        for (int q = 0; q < 3; ++q)
        {
            value_t r = (vi - vg3[q]).Norm();
            auto green = exp(-J0 * k_ * r) / r;
            line_term += glw3[p] * w3[q] * green;
        }
    }
    return line_term * (v1 - v2).Norm();
}

Complex FastWT::consistencyTerm(const VectorR3 * vx3, const VectorR3 & v1, const VectorR3 & v2) const
{
    return (Z0 / (k_ * PI4 * J0)) * lineIntegral(vx3, v1, v2);
}

value_t FastWT::stabilizationTerm(value_t len) const
{
    return 2.0f * len * beta_ * Z0 / k_;
}
