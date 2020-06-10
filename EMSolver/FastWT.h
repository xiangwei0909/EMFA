//  Description:    Fast Electromagnetic Solver Based on Preconditioning IP-SIE (FastWT)
//  Author:         ZhangQiang
//  Date:           2017-05-15
//  License:        MIT
//  Note:           Multithreading

#pragma once
#include "EM.h"
#include "ContourEdge.h"
#include "MLFMA_Boundary.h"
#include "VectorC3.h"
#include "ThreadPool.h"

namespace mom {

using component::VectorC3;
using component::MLFMA_Boundary;

class FastWT : public EM {
    struct Box {
        size_t              offset;
        size_t              index;
        VectorR3            center;
        Box*                parent = nullptr;
        std::vector<int>    unknown;
        std::vector<Box*>   child;
        std::vector<Box*>   near;
        std::vector<Box*>   far;
    };
    using MeshPtr = std::shared_ptr<component::Mesh>;
    using ComEdgePtr = std::shared_ptr<component::ContourEdge>;
    using BoxArray = std::vector<std::vector<Box>>;
    using KArray = wt::QMat<VectorR3>;
    using VC3Array = wt::QMat<VectorC3>;
    using IPMat = wt::QSpMat<value_t>;
    using TransArray = std::vector<Qcx_mat>;
    using MutexGuard = std::lock_guard<std::mutex>;
public:
    FastWT();
    ~FastWT();

public:
    void init(component::ConfigLoader* ploader) override;
    void solve() override;
    void output() override;
    void clear() override;
    void reportInfo(Qostream& strm) const override;

private:
    void buildBox();
    void fillCenterPoints(std::vector<VectorR3>& centers);
    void buildFirstLevel(const std::vector<VectorR3>& cps, value_t box_len);
    void buildCurLevel(const std::vector<VectorR3>& cps, int cur_level, value_t box_len);
    void setBoxParameter();
    int getMultipoleNumber(value_t box_len) const;
    void prepareAssignment(size_t thread_num);

    void levelRootsAndWeights();
    void calculateInterpolationMtr();
    void lagrangeIPWeights(const Qvec& nxt_sam, const Qvec& cur_sam, Qumat& loc, Qmat& elems) const;
    void calculateLevelK();

    void calculateLevelTrans();
    void parallelLevelTrans(size_t begin, size_t end, size_t idx);

    void calculateRadAndRecv();
    void parallelRadAndRecv(size_t begin, size_t end);
    VectorC3 radiationFunction(const VectorR3* vpgs3, const VectorR3& vp, const VectorR3& center, const VectorR3& sam_k) const;
    VectorC3 halfRecvFunction(const VectorR3* vpgs3, const VectorR3& vp, const VectorR3& nml_plu, const VectorR3& center, const VectorR3& sam_k) const;
    Complex lineRadiationFunction(const VectorR3* vpgs3, const VectorR3& vp, const VectorR3& center, const VectorR3& sam_k) const;

    void calculateNearZ();
    void parallelNearZ(size_t begin, size_t end, Qumat* ploc, Qcx_vec* pelems, size_t ebase) const;
    void fillNearZbyColumn(const std::vector<int>& flds, int sid, size_t fpos, size_t spos,
                           Qumat& loc, Qcx_vec& elems, size_t& aux) const;
    void fillVbyBox();

    void calculatePreconditioner();
    void preparepreconditioner(const std::vector<VectorR3>& centers, std::vector<size_t>& base) const;
    void parallelPreconditioner(size_t begin, size_t end, const std::vector<VectorR3>& centers, Qumat *ploc, Qcx_vec *pelems, size_t ebase) const;

    int iterativeSolve();
    void allocateForSolve();
    void calculateCoeff();
    Qcx_vec matrixVectorMultiply(const Qcx_vec& b);
    void parallelFinestLayerUpward(size_t begin, size_t end, const Qcx_vec& b);
    void parallelFinestLayerDownward(size_t begin, size_t end);
    void parallelFar(size_t begin, size_t end, Qcx_vec *pfar) const;
    void upwardPass();
    void parallelUpwardPass(size_t begin, size_t end, int cur_level);
    void downwardPass();
    void parallelDownwardPass(size_t begin, size_t end, size_t idx);

    void reportLevelInfo(Qostream& strm) const;
    void reportTaskAssigmentInfo(Qostream& strm) const;
    void reportLevelDetail() const;
    void reportFastWTInfo(Qostream& strm) const;

    value_t getBiRCS(const VectorR3& sca_k) const;
    void getEFiled(component::FieldData* pdata, const VectorR3& rad_k, const VectorR3& rad_ev, const VectorR3& rad_eh) const;
    void calculateSurfaceCurrent(std::vector<component::CurrentData> *currents) const;

    Complex cZKernel(const VectorR3 *vf3, const VectorR3 *vs3, const VectorR3 &vf, const VectorR3 &vs, const VectorR3 &nml) const;
    Complex cZSingular(const VectorR3 *vf3, const VectorR3 *vs3, const VectorR3 &vf, const VectorR3 &vs) const;
    Complex cVKernel(const VectorR3 *vx3, const VectorR3 &vp, const VectorR3& nml) const;

    Complex lineIntegral(const VectorR3* vx3, const VectorR3& v1, const VectorR3& v2) const;
    Complex consistencyTerm(const VectorR3* vx3, const VectorR3& v1, const VectorR3& v2) const;
    value_t stabilizationTerm(value_t len) const;

private:
    bool preconditioning_;
    int unknowns_;
    int level_;
    int max_iteration_num_;

    size_t unk_near_, precond_num_;
    value_t alpha_;
    value_t k_;
    value_t beta_;
    value_t iteration_threshold_;
    value_t col_threshold_, row_threshold_;

    MeshPtr mesh_ptr_;
    ComEdgePtr ce_ptr_;

    MLFMA_Boundary boundary_;
    BoxArray level_box_;

    Qsp_cx_mat nearZ_, precond_;
    Qvec coeff_, pcoeff_;
    Qcx_vec V_;
    Qcx_vec I_;

    std::vector<value_t> box_len_;
    std::vector<Qvec> glxt_, glxp_, glw_;
    std::vector<IPMat> arr_ip_;
    std::vector<KArray> level_k_;
    std::vector<TransArray> level_trans_;
    std::vector<VC3Array> rad_, rec_;
    std::vector<VC3Array> level_rad_, level_rec_;
    std::vector<Qcx_mat> prad_, prec_;
    std::vector<Qcx_mat> plevel_rad_, plevel_rec_;

    mutable std::mutex mutex_;
    arma::Mat<size_t> assignment_;
    component::ThreadPool pool_;
};

} // namespace mom