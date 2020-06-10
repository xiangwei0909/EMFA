//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "EM.h"
#include "Mathdef.h"

namespace mom {

class CFIE :public EM {
    using MeshPointer   = std::shared_ptr<component::Mesh>;
    using CEPointer     = std::shared_ptr<component::CommonEdge>;
    using RadVector     = std::vector<std::pair<int, int>>;

public:
    CFIE();
    ~CFIE();

public:
    void        init(component::ConfigLoader * ploader) override;
    void        solve() override;
    void        output() override;
    void        clear() override;
    void        reportInfo(Qostream& strm) const override;

protected:
    Complex     cZppKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &nf);
    Complex     cZpmKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &nf);
    Complex     cZmpKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &nf);
    Complex     cZmmKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &nf);

    Complex     cZppSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
    Complex     cZpmSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
    Complex     cZmpSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
    Complex     cZmmSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);

    Complex     cVKernel(VectorR3 *vp3, VectorR3 *vm3, VectorR3 &np, VectorR3 &nm, VectorR3 &vp, VectorR3 &vm);

    void        readExcEdges(const Qstring& rad_file);
    void        radiateV();
    value_t     getBiRCS(const VectorR3& sca_k) const;
    void        getEFiled(component::FieldData* pdata, const VectorR3& rad_k, const VectorR3& rad_ev, const VectorR3& rad_eh) const;
    void        calculateSurfaceCurrent(std::vector<component::CurrentData> *currents) const;
    bool        writeZIVData();

    virtual void fillZ();
    virtual void fillV();
protected:
    int         unknowns_;
    value_t     k_;
    value_t     alpha_;

    Qcx_mat     Z;
    Qcx_vec     I;
    Qcx_vec     V;

    RadVector   exc_edge_;
    MeshPointer mesh_ptr_;
    CEPointer   ce_ptr_;
};

} // namespace mom

