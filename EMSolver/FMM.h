//  Description:    Fast Multipole Method (FMM)
//  Author:         ZhangQiang
//  Date:           2016-12-17
//  License:        MIT
//  Note:           Optimize at 2016-12-30 

#pragma once
#include "EM.h"
#include "VectorC3.h"
#include "FMM_Boundary.h"
#include "Quadrature.h"

namespace mom {

using component::VectorC3;
using component::FMM_Boundary;

class FMM : public EM {
    struct Box {
        size_t              index;
        size_t              offset;
        VectorR3            center;
        std::vector<Box*>   near;
        std::vector<Box*>   far;
        std::vector<int>    unknown;
    };
    using MeshPtr       = std::shared_ptr<component::Mesh>;
    using ComEdgePtr    = std::shared_ptr<component::CommonEdge>;
    using BoxArray      = std::vector<Box>;
    using VC3Array      = wt::QMat<VectorC3>;
    using KArray        = wt::QMat<VectorR3>;
    using ResultArray   = std::vector<std::vector<value_t>>;
    using TransArray    = std::vector<Qcx_mat>;
    using RadVector     = std::vector<std::pair<int, int>>;

public:
    FMM();
    ~FMM();

    void        init(component::ConfigLoader* ploader) override;
    void        solve() override;
    void        output() override;
    void        clear() override;
    void        reportInfo(Qostream& strm) const override;
    void        Test();

private:
    void        setBoxParameter();
    void        prepareSolve();
    int         iterativeSolve();
    Qcx_vec     matrixVectorMultiply(const Qcx_vec& b);
    void        preCalculateKArray();
    void        preCalculateTrans();
    void        preCalculateRadAndRecv();
    VectorC3    radiationFunction(const VectorR3* vpgs3, const VectorR3* vmgs3, const VectorR3& vp,
                                  const VectorR3& vm, const VectorR3& center, const VectorR3& sam_k);
    VectorC3    halfRecvFunction(const VectorR3* vpgs3, const VectorR3* vmgs3, const VectorR3& vp,
                                 const VectorR3& vm, const VectorR3& nml_plu, const VectorR3& nml_min,
                                 const VectorR3& center, const VectorR3& sam_k);
    void        calculateNearZ();
    void        fillNearZbyColumn(const std::vector<int>& flds, int sid, size_t fpos, size_t spos,
                                  Qumat& loc, Qcx_vec& elems, size_t& aux);
    void        fillVbyBox();
    bool        readExcEdges(const Qstring& rad_file, Qstring& stateInfo);
    void        readExcEdges(const Qstring& rad_file);
    void        radiateV();

    value_t     getBiRCS(const VectorR3& sca_k) const;
    void        getEFiled(component::FieldData* pdata, const VectorR3& rad_k, const VectorR3& rad_ev, const VectorR3& rad_eh) const;

    void        reportFFMInfo(Qostream& strm) const;
    void        reportBoxInfo() const;

    Complex     cZppKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &nf);
    Complex     cZpmKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &nf);
    Complex     cZmpKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &nf);
    Complex     cZmmKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &nf);

    Complex     cZppSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
    Complex     cZpmSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
    Complex     cZmpSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
    Complex     cZmmSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
    
    Complex     cVKernel(VectorR3 *vp3, VectorR3 *vm3, VectorR3 &np, VectorR3 &nm, VectorR3 &vp, VectorR3 &vm);
private:
    int         unknowns_;
    int         sam_theta_;
    int         sam_phi_;
    int         max_iteration_num_;

    size_t      unk_near_;
    value_t     alpha_;
    value_t     k_;
    value_t     iteration_threshold_;

    MeshPtr         mesh_ptr_;
    ComEdgePtr      ce_ptr_;
    FMM_Boundary    fmm_box_;
    BoxArray        box_array_;
    KArray          arr_k_;

    TransArray      transfer_;
    Qsp_cx_mat      nearZ_;
    Qcx_vec         V_;
    Qcx_vec         I_;
    RadVector       exc_edge_;
    ResultArray     rcs;
    Qvec            coeff_;
    VC3Array        sk_, gk_;

    std::vector<VC3Array>   rad_, rec_;
    math::GaussLegendre     gl_theta_, gl_phi_;
};

} // namespace mom