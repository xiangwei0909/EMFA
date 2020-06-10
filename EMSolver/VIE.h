//**********************************************************
// Author: Xiang Wei
// License: MIT
//**********************************************************

#pragma once
#include "EM.h"
#include "VectorC3.h"
namespace mom {

typedef wt::QMat<component::VectorC3> Qc3_mat;
typedef wt::QMat<component::VectorR3> Qr3_mat;
typedef std::vector<VectorR3> Vec_R3;
class VIE :public EM {
	using MeshPointer = std::shared_ptr<component::Mesh>;
	using CEPointer = std::shared_ptr<component::CommonTriangle>;
	using RadVector = std::vector<std::pair<int, int>>;
public:
	VIE();
	~VIE();
public:
	void        init(component::ConfigLoader* ploader) override;
	void        solve() override;
	void        output() override;
	void        clear() override;
	void        reportInfo(Qostream& strm) const override;

protected:
	Complex     vZppKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, value_t &epsp, value_t &epsm);
	Complex     vZpmKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, value_t &epsp, value_t &epsm);
	Complex     vZmpKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, value_t &epsp, value_t &epsm);
	Complex     vZmmKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, value_t &epsp, value_t &epsm);

	Complex     eZppSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, value_t &epsp, value_t &epsm);
	Complex     eZpmSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, value_t &epsp, value_t &epsm);
	Complex     eZmpSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, value_t &epsp, value_t &epsm);
	Complex     eZmmSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, value_t &epsp, value_t &epsm);
	
	Complex     FZppKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_p, value_t &epsp, value_t &epsm);
	Complex     FZpmKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_m, value_t &epsp, value_t &epsm);
	Complex     FZmpKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_p, value_t &epsp, value_t &epsm);
	Complex     FZmmKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_m, value_t &epsp, value_t &epsm);

	Complex     FZppSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_p, value_t &epsp, value_t &epsm);
	Complex     FZpmSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_m, value_t &epsp, value_t &epsm);
	Complex     FZmpSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_p, value_t &epsp, value_t &epsm);
	Complex     FZmmSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_m, value_t &epsp, value_t &epsm);

	
	value_t		IsSingular(VectorR3 *vs3, VectorR3 &vf);
	value_t		IvSingular(VectorR3 *vs3, VectorR3 &vf);
	VectorR3	IvpSingular(VectorR3 *vs3, VectorR3 &vf);
	void		TetToTri(VectorR3 *v4, VectorR3 **v4_3);

	Complex		z1(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, value_t &epsp, value_t &epsm);

	Complex     eVKernel(VectorR3 *vp4, VectorR3 *vm4, VectorR3 &vp, VectorR3 &vm, int &_nvm);
protected:
	void		prepareVIE(value_t _eps1, value_t _mu1);
	bool        readExcEdges(const Qstring& rad_file, Qstring& stateInfo);
	void        readExcEdges(const Qstring& rad_file);
	void        fillZ();
	void		FastfillZ();
	void		KernelIntegral();
	void		FillZKernel(VectorR3 *vm4, VectorR3 *vn4, int &m, int &n, VectorR3 *v4);

	void        fillV();
	void        radiateV();
	value_t     getBiRCS(const VectorR3& sca_k) const;
	void        getEFiled(component::FieldData* pdata, const VectorR3& rad_k, const VectorR3& rad_ev, const VectorR3& rad_eh) const;

	bool        writeZIVData();
protected:
	int         unknowns_;
	int			tet_num;
	int			Isfast;
	value_t		mu1_, eps1_;
	value_t     k_;
	value_t		omiga;
	Qcx_mat     Z;
	
	Qcx_mat		Z2K1, Z2K4;
	Qc3_mat		Z2K2, Z2K3;

	Qcx_mat		Z2EK1, Z2EK4;
	Qc3_mat		Z2EK2, Z2EK3;

	Qcx_mat		Z5K;

	Qcx_vec		RIVP, RRIV, IV;
	Vec_R3		RIV, IVP;

	Qcx_vec     I;
	Qcx_vec     V;

	RadVector   exc_edge_;

	MeshPointer mesh_ptr_;
	CEPointer   ct_ptr_;
};

} // namespace mom