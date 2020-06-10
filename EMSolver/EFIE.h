//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "EM.h"

namespace mom {

class EFIE :public EM {
    using MeshPointer   = std::shared_ptr<component::Mesh>;
    using CEPointer     = std::shared_ptr<component::CommonEdge>;
    using RadVector     = std::vector<std::pair<int, int>>;
public:
    EFIE();
    ~EFIE();
public:
    void        init(component::ConfigLoader* ploader) override;
    void        solve() override;
    void        output() override;
    void        clear() override;
    void        reportInfo(Qostream& strm) const override;

protected:
    Complex     eZppKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
    Complex     eZpmKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
    Complex     eZmpKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
    Complex     eZmmKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);

    Complex     eZppSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
    Complex     eZpmSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
    Complex     eZmpSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
    Complex     eZmmSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);

    Complex     eVKernel(VectorR3 *vp3, VectorR3 *vm3, VectorR3 &vp, VectorR3 &vm);
	Complex		EDMKernel(VectorR3 *vf_p3, VectorR3 *vf_m3, VectorR3 *vs_p3, VectorR3 *vs_m3, value_t &l_fld, value_t &l_src);
	VectorR3	center(VectorR3 *v3);
protected:
    bool        readExcEdges(const Qstring& rad_file, Qstring& stateInfo);
    void        readExcEdges(const Qstring& rad_file);
    void        fillZ();
    void        fillV();
    void        radiateV();
    value_t     getBiRCS(const VectorR3& sca_k) const;
    void        getEFiled(component::FieldData* pdata, const VectorR3& rad_k, const VectorR3& rad_ev, const VectorR3& rad_eh) const;
	void        getCurrentOnFeed();

    bool        writeZIVData();
	void        calculateSurfaceCurrent(std::vector<component::CurrentData> *currents) const;
protected:
    int         unknowns_;
	value_t		threshold_edm;
    value_t     k_;
    Qcx_mat     Z;
	CMatrix		Z_eigen;
    Qcx_vec     I;
	CVector		I_eigen;
    Qcx_vec     V;
	CVector		V_eigen;
	Complex		Sigma;
    RadVector   exc_edge_;

	int			Array_x;
	int			Array_y;
	value_t		Dx;
	value_t		Dy;

    MeshPointer mesh_ptr_;
    CEPointer   ce_ptr_;
};
inline VectorR3 EFIE::center(VectorR3 *v3)
{
	return ((v3[0] + v3[1] + v3[2]) / 3.0f);
}
} // namespace mom