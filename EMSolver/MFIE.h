//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "EM.h"

namespace mom {

class MFIE :public EM {
    using MeshPointer   = std::shared_ptr<component::Mesh>;
    using CEPointer     = std::shared_ptr<component::CommonEdge>;

public:
    MFIE();
    ~MFIE();
public:
    void        init(component::ConfigLoader* ploader) override;
    void        solve() override;
    void        output() override;
    void        clear() override;
    void        reportInfo(Qostream& strm) const override;

protected:
    Complex     mZppKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &nf);
    Complex     mZpmKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &nf);
    Complex     mZmpKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &nf);
    Complex     mZmmKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &nf);

    Complex     mZppSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
    Complex     mZpmSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
    Complex     mZmpSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
    Complex     mZmmSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);

    Complex     mVKernel(VectorR3 *vp3, VectorR3 *vm3, VectorR3 &np, VectorR3 &nm, VectorR3 &vp, VectorR3 &vm);
private:

    void        fillZ();
    void        fillV();
    value_t     getBiRCS(const VectorR3& sca_k) const;
    bool        writeZIVData();
private:
    int         unknowns_;
    value_t     k_;
    Qcx_mat     Z;
    Qcx_vec     I;
    Qcx_vec     V;

    MeshPointer mesh_ptr_;
    CEPointer   ce_ptr_;
};

} // namespace mom
