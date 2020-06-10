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

	class PGF_VSIE :public EM {
		using MeshPointer = std::shared_ptr<component::Mesh>;
		using CEPointer = std::shared_ptr<component::CommonEdge>;
		using CTPointer = std::shared_ptr<component::CommonTriangle>;
		using RadVector = std::vector<std::pair<int, int>>;
	public:
		PGF_VSIE();
		~PGF_VSIE();
	public:
		void        init(component::ConfigLoader* ploader) override;
		void        solve() override;
		void        output() override;
		void        clear() override;
		void        reportInfo(Qostream& strm) const override;

	protected:
		Complex     DDZppKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_p, value_t &epsp, value_t &epsm);
		Complex     DDZpmKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_m, value_t &epsp, value_t &epsm);
		Complex     DDZmpKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_p, value_t &epsp, value_t &epsm);
		Complex     DDZmmKernel(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_m, value_t &epsp, value_t &epsm);

		Complex     DDZppSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_p, value_t &epsp, value_t &epsm);
		Complex     DDZpmSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_m, value_t &epsp, value_t &epsm);
		Complex     DDZmpSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_p, value_t &epsp, value_t &epsm);
		Complex     DDZmmSingular(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_m, value_t &epsp, value_t &epsm);

		Complex		DMZppKernel(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, int &_vm_f);
		Complex		DMZpmKernel(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, int &_vm_f);
		Complex		DMZmpKernel(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, int &_vm_f);
		Complex		DMZmmKernel(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, int &_vm_f);

		Complex		DMZppSingular(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, int &_vm_f);
		Complex		DMZpmSingular(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, int &_vm_f);


		Complex		MDZppKernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm);
		Complex		MDZmpKernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm);
		Complex		MDZpmKernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm);
		Complex		MDZmmKernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm);

		Complex		MDZppSingular(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm);
		Complex		MDZmpSingular(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm);

		Complex     MMZppKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
		Complex     MMZpmKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
		Complex     MMZmpKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
		Complex     MMZmmKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);

		Complex     MMZppSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
		Complex     MMZpmSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
		Complex     MMZmpSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
		Complex     MMZmmSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);

		value_t		IsSingular(VectorR3 *vs3, VectorR3 &vf);
		VectorR3    IspSingular(VectorR3 * vs3, VectorR3 & vf);
		value_t		IvSingular(VectorR3 *vs3, VectorR3 &vf);
		VectorR3	IvpSingular(VectorR3 *vs3, VectorR3 &vf);
		void		TetToTri(VectorR3 *v4, VectorR3 **v4_3);

		Complex		z1(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, value_t &epsp, value_t &epsm);

		Complex     DVKernel(VectorR3 *vp4, VectorR3 *vm4, VectorR3 &vp, VectorR3 &vm, int &_nvm);
		Complex     MVKernel(VectorR3 *vp3, VectorR3 *vm3, VectorR3 &vp, VectorR3 &vm);

		Complex		FSPGF(VectorR3 &r, int _t_sum, value_t _Dx, value_t _Dy) const;
		Complex		FSPGF_singular(VectorR3 &r, int _t_sum, value_t _Dx, value_t _Dy);
		Complex		erfz(Complex z) const;
		Complex		erfc(Complex z) const{ return 1.0f - erfz(z); }

		void		KernelIntegral();
		void		FillZKernel(VectorR3 *vm4, VectorR3 *vn4, int &m, int &n, VectorR3 *v4);
		void		CalculatePGFgrid();
		Complex		InterpolarPGF(VectorR3 &r);
	protected:
		void		prepareFSPGF(value_t _Dx, value_t _Dy, int _Array_x, int _Array_y, int _t_sum);
		void		prepareVIE(value_t _eps1, value_t _mu1);
		bool        readExcEdges(const Qstring& rad_file, Qstring& stateInfo);
		void        readExcEdges(const Qstring& rad_file);
		void        fillZ();
		void		fillZDD();
		void		fillZDM();
		void		fillZMD();
		void		fillZMM();
		void        fillV();
		void        radiateV();
		value_t     getBiRCS(const VectorR3& sca_k) const;
		void        getEFiled(component::FieldData* pdata, const VectorR3& rad_k, const VectorR3& rad_ev, const VectorR3& rad_eh) const;
		value_t		getReflectField(const VectorR3& sca_k) const;
		void		getNearEField(std::vector<component::NearFieldData>* data) const;
		void		NEFkernel_vol(VectorR3& ob, int& unk, component::VectorC3& e_vol) const;
		void		NEFkernel_surf(VectorR3& ob, int& unk, component::VectorC3& e_surf) const;
		void		calculateSurfaceCurrent(std::vector<component::CurrentData>* currents) const;

		bool        writeZIVData();
	protected:
		int			tet_num;
		int			iter_num;
		value_t		iter_threshold;

		Complex		Sigma;

		int			grid_x;
		int			grid_y;
		int			grid_z;
		value_t		d_x;
		value_t		d_y;
		value_t		d_z;
		VectorR3	min_boundary, max_boundary;
		std::vector<Complex> PGFgrid;

		int			isContinuous;
		value_t		Dx;
		value_t		Dy;
		int			Array_x;
		int			Array_y;
		int			t_sum;

		int         unknowns_;
		int			unknowns_t;
		int			unknowns_e;
		
		value_t		mu1_, eps1_;
		value_t     k_;
		value_t		omiga;
		Qcx_mat     Z;
		Qcx_vec     I;
		Qcx_vec     V;

		Qcx_mat		Z2K1, Z2K4;
		Qc3_mat		Z2K2, Z2K3;

		Qcx_mat		Z2EK1, Z2EK4;
		Qc3_mat		Z2EK2, Z2EK3;

		Qcx_mat		Z5K;

		Qcx_vec		RIVP, RRIV, IV;
		Vec_R3		RIV, IVP;

		RadVector   exc_edge_;

		MeshPointer mesh_ptr_;
		CEPointer   ce_ptr_;
		CTPointer   ct_ptr_;
	};

} // namespace mom
