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
class DGVSIE :public EM {
	using MeshPointer = std::shared_ptr<component::Mesh>;
	using CEPointer = std::shared_ptr<component::CommonEdge>;
	using CTPointer = std::shared_ptr<component::CommonTriangle>;
	using RadVector = std::vector<std::pair<int, int>>;
public:
	DGVSIE();
	~DGVSIE();
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

	Complex		HZppKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &t_f, int &t_s, value_t eps_f, value_t eps_s);
	Complex		HZ1Singular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &t_f, int &t_s, value_t eps_f, value_t eps_s);
	Complex		HZ2Singular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &t_f, int &t_s, value_t eps_f, value_t eps_s);
	Complex		HZ3Singular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &t_f, int &t_s, value_t eps_f, value_t eps_s);

	Complex		FHZppKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &t_f, int &t_s, value_t eps_f, value_t eps_s);
	Complex		FHZ1Singular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &t_f, int &t_s, value_t eps_f, value_t eps_s);
	Complex		FHZ2Singular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &t_f, int &t_s, value_t eps_f, value_t eps_s);
	Complex		FHZ3Singular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &t_f, int &t_s, value_t eps_f, value_t eps_s);

	Complex		ZSVppKernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, value_t eps_s);
	Complex		ZSVmpKernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, value_t eps_s);
	Complex		ZSVppSingular(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, value_t eps_s);
	Complex		ZSVmpSingular(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, value_t eps_s);

	Complex		ZVSppKernel(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, value_t eps_s);
	Complex		ZVSpmKernel(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, value_t eps_s);
	Complex		ZVSppSingular(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, value_t eps_s);
	Complex		ZVSpmSingular(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, value_t eps_s);

	Complex		ZSSppKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs);
	Complex		ZSSpmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs);
	Complex		ZSSmpKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs);
	Complex		ZSSmmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs);
	Complex		ZSSppSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs);
	Complex		ZSSpmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs);
	Complex		ZSSmpSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs);
	Complex		ZSSmmSingular(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs);

	Complex		FSPGF(VectorR3 &r, int _t_sum, value_t _Dx, value_t _Dy);
	Complex		erfz(Complex z);
	Complex		erfc(Complex z) { return 1.0f - erfz(z); }

	value_t		IsSingular(VectorR3 *vs3, VectorR3 &vf);
	value_t		IvSingular(VectorR3 *vs3, VectorR3 &vf);
	VectorR3	IvpSingular(VectorR3 *vs3, VectorR3 &vf);
	void		TetToTri(VectorR3 *v4, VectorR3 **v4_3);

	Complex		z1(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, value_t &epsp, value_t &epsm);

	Complex     eVKernel(VectorR3 *vp4, VectorR3 *vm4, VectorR3 &vp, VectorR3 &vm, int &_nvm);
	Complex		eHVKernel(VectorR3 *v4, VectorR3 &vx);
	Complex		eMVKernel(VectorR3 *vp3, VectorR3 *vm3, VectorR3 &vp, VectorR3 &vm);
protected:
	void		prepareVIE(value_t _eps1, value_t _mu1);
	void		prepareFSPGF(value_t _Dx, value_t _Dy, int _Array_x, int _Array_y, int _t_sum);
	bool        readExcEdges(const Qstring& rad_file, Qstring& stateInfo);
	void        readExcEdges(const Qstring& rad_file);
	void        fillZ();
	void		FastfillZ();
	void		KernelIntegral();
	void		CalculatePGFgrid();
	Complex		InterpolarPGF(VectorR3 &r);
	void		FillZKernel(VectorR3 *vm4, VectorR3 *vn4, int &m, int &n, VectorR3 *v4);
	void		FillZVV();
	void		FillZVS();
	void		FillZSV();
	void		FillZSS();

	void        fillV();
	void		FillHV();
	void        radiateV();
	value_t     getBiRCS(const VectorR3& sca_k) const;
	value_t		getBiRCS_h(const VectorR3& sca_k) const;
	void        getEFiled(component::FieldData* pdata, const VectorR3& rad_k, const VectorR3& rad_ev, const VectorR3& rad_eh) const;

	bool        writeZIVData();
protected:
	int         unknowns_;
	int			unknowns_v;
	int			unknowns_s;
	int			tet_num;
	int			Isfast;
	int			Array_x;
	int			Array_y;
	int			t_sum;
	int			grid_x;
	int			grid_y;
	int			grid_z;
	int			isContinuous;
	int			iter_num;
	value_t		iter_threshold;
	VectorR3	min_boundary, max_boundary;

	value_t		mu1_, eps1_;
	value_t     k_;
	value_t		Dx;
	value_t		Dy;
	value_t		d_x;
	value_t		d_y;
	value_t		d_z;
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
	std::vector<Complex> PGFgrid;
	RadVector   exc_edge_;

	MeshPointer mesh_ptr_;
	CEPointer	ce_ptr_;
	CTPointer   ct_ptr_;
};

} // namespace mom


