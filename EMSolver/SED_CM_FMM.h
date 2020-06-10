//**********************************************************
// Author: Xiang Wei
// License: MIT
//**********************************************************

#pragma once
#include "EM.h"
#include "VectorC3.h"
#include "FMM_Boundary.h"
#include "Quadrature.h"
namespace mom {

	typedef wt::QMat<component::VectorC3> Qc3_mat;
	typedef wt::QMat<component::VectorR3> Qr3_mat;
	typedef std::vector<VectorR3> Vec_R3;
	using component::VectorC3;

	class SED_CM_FMM :public EM {
		using MeshPointer = std::shared_ptr<component::Mesh>;
		using CEPointer = std::shared_ptr<component::CommonEdge>;
		using CTPointer = std::shared_ptr<component::CommonTriangle>;
		using RadVector = std::vector<std::pair<int, int>>;
		using MatArray = std::vector<Qcx_mat>;
		using MatArrayEi = std::vector<CMatrix>;
		using MatVec = std::vector<Qcx_vec>;
		using MatVecEi = std::vector<CVector>;
		using VC3Array = wt::QMat<VectorC3>;
		using KArray = wt::QMat<VectorR3>;
		using TransArray = std::vector<CMatrix>;
		using MutexGuard = std::lock_guard<std::mutex>;
	public:
		SED_CM_FMM();
		~SED_CM_FMM();
	public:
		struct KIM
		{
			KIM(int n)
			{
				Z2K1.zeros(n, n);
				Z2K2.set_size(n, n);
				Z2K3.set_size(n, n);
				Z2K4.zeros(n, n);
				Z5K.zeros(n, n);

				RIVP.zeros(n);
				RRIV.zeros(n);
				IV.zeros(n);

				RIV.resize(n);
				IVP.resize(n);
			}
			Qcx_mat		Z2K1, Z2K4;
			Qc3_mat		Z2K2, Z2K3;

			Qcx_mat		Z5K;

			Qcx_vec		RIVP, RRIV, IV;
			Vec_R3		RIV, IVP;
		};
		struct Blk
		{
			Blk(int idx, int idy)
			{
				_idx = idx;
				_idy = idy;
			}
			int _idx;
			int _idy;
			int _id_sed;
			int _id_near;
			int unknown;
			int offset=0;
			int preoffset = 0;
			std::vector<Blk*> near;
		};
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

		Complex     DDZppKernel_fast(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_p, value_t &epsp, value_t &epsm, KIM& ik);
		Complex     DDZpmKernel_fast(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_m, value_t &epsp, value_t &epsm, KIM& ik);
		Complex     DDZmpKernel_fast(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_p, value_t &epsp, value_t &epsm, KIM& ik);
		Complex     DDZmmKernel_fast(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_m, value_t &epsp, value_t &epsm, KIM& ik);

		Complex     DDZppSingular_fast(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_p, value_t &epsp, value_t &epsm, KIM& ik);
		Complex     DDZpmSingular_fast(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_p, int &t_s_m, value_t &epsp, value_t &epsm, KIM& ik);
		Complex     DDZmpSingular_fast(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_p, value_t &epsp, value_t &epsm, KIM& ik);
		Complex     DDZmmSingular_fast(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, VectorR3 *vs3, int &_vm_f, int &_vm_s, int &t_f_m, int &t_s_m, value_t &epsp, value_t &epsm, KIM& ik);

		Complex		DMZppKernel(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, int &_vm_f);
		Complex		DMZpmKernel(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, int &_vm_f);
		Complex		DMZmpKernel(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, int &_vm_f);
		Complex		DMZmmKernel(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, int &_vm_f);

		Complex		DMZppSingular(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, int &_vm_f);
		Complex		DMZpmSingular(VectorR3 *vf4, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 *vf3, int &_vm_f);


		Complex		MDZppkernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm);
		Complex		MDZmpkernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm);
		Complex		MDZpmkernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm);
		Complex		MDZmmkernel(VectorR3 *vf3, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, VectorR3 *vs3, int &_vm_s, value_t &epsp, value_t &epsm);

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

		Complex		MMEDMKernel(VectorR3 *vf_p3, VectorR3 *vf_m3, VectorR3 *vs_p3, VectorR3 *vs_m3, value_t &l_fld, value_t &l_src);

		value_t		IsSingular(VectorR3 *vs3, VectorR3 &vf);
		VectorR3    IspSingular(VectorR3 * vs3, VectorR3 & vf);
		value_t		IvSingular(VectorR3 *vs3, VectorR3 &vf);
		VectorR3	IvpSingular(VectorR3 *vs3, VectorR3 &vf);
		void		TetToTri(VectorR3 *v4, VectorR3 **v4_3);

		Complex		z1(VectorR3 *vf4, VectorR3 *vs4, VectorR3 &vf, VectorR3 &vs, value_t &epsp, value_t &epsm);

		Complex     DVKernel(VectorR3 *vp4, VectorR3 *vm4, VectorR3 &vp, VectorR3 &vm, int &_nvm);
		Complex     MVKernel(VectorR3 *vp3, VectorR3 *vm3, VectorR3 &vp, VectorR3 &vm);

		Complex		FSPGF(VectorR3 &r, int _t_sum, value_t _Dx, value_t _Dy);
		Complex		erfz(Complex z);
		Complex		erfc(Complex z) { return 1.0f - erfz(z); }

		void		KernelIntegral();
		KIM 		SED_KernelIntegral(int &sedid_f, int &sedid_s);
		KIM 		RED_KernelIntegral(int &xid_f, int &yid_f, int &xid_s, int &yid_s);
		void		FillZKernel(VectorR3 *vm4, VectorR3 *vn4, int &m, int &n, VectorR3 *v4);
		void		FillZKernel_SED_COM(VectorR3 *vm4, VectorR3 *vn4, int &m, int &n, VectorR3 *v4, KIM& ik);
		void		FillZKernel_SED_DIFF(VectorR3 *vm4, VectorR3 *vn4, int &m, int &n, VectorR3 *v4, KIM& ik);
		void		CalculatePGFgrid();
		Complex		InterpolarPGF(VectorR3 &r);

	protected:
		void		prepareFSPGF(value_t _Dx, value_t _Dy, int _Array_x, int _Array_y, int _t_sum);
		void		prepareVIE(value_t _eps1, value_t _mu1);
		bool        readExcEdges(const Qstring& rad_file, Qstring& stateInfo);
		void        readExcEdges(const Qstring& rad_file);
		void        fillZ();
		void		fillZDD(int &start_f, int &end_f, int &start_s, int &end_s);
		void		fillZMD(int &start_f, int &end_f, int &start_s, int &end_s);
		void		fillZDM(int &start_f, int &end_f, int &start_s, int &end_s);
		void		fillZMM(int &start_f, int &end_f, int &start_s, int &end_s);
		void        fillV();
		void        radiateV();
		value_t     getBiRCS(const VectorR3& sca_k);
		void        getEFiled(component::FieldData* pdata, const VectorR3& rad_k, const VectorR3& rad_ev, const VectorR3& rad_eh) const;
		void		getCurrentOnFeed();
		void		getNearEField(std::vector<component::NearFieldData>* data) const;
		void		NEFkernel_vol(VectorR3& ob, int& unk, component::VectorC3& e_vol) const;
		void		NEFkernel_surf(VectorR3& ob, int& unk, component::VectorC3& e_surf) const;
		void		NEFkernel_surf(VectorR3& ob, int& _x, int& _y, int& unk, component::VectorC3& e_surf) const;
		value_t		getReflectField(const VectorR3& sca_k) const;
		void		calculateSurfaceCurrent(std::vector<component::CurrentData>* currents) const;
		void		readFEKOcurrent();
		bool        writeZIVData();
	protected:
		void		fillSEDZ();
		void		fillSEDZDD(int &_f_sed, int &_s_sed, int &_start_tf, int &_end_tf, int &_start_ts, int &_end_ts, int &_start_ef, int &_end_ef, int &_start_es, int &_end_es);
		void		fillSEDZDM(int &_f_sed, int &_s_sed, int &_start_tf, int &_end_tf, int &_start_ts, int &_end_ts, int &_start_ef, int &_end_ef, int &_start_es, int &_end_es);
		void		fillSEDZMD(int &_f_sed, int &_s_sed, int &_start_tf, int &_end_tf, int &_start_ts, int &_end_ts, int &_start_ef, int &_end_ef, int &_start_es, int &_end_es);
		void		fillSEDZMM(int &_f_sed, int &_s_sed, int &_start_tf, int &_end_tf, int &_start_ts, int &_end_ts, int &_start_ef, int &_end_ef, int &_start_es, int &_end_es);
		void		fillCOUPZDD(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n);
		void		fillCOUPZDM(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n);
		void		fillCOUPZMD(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n);
		void		fillCOUPZMM(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n);
		void		fillSEDV();
		void		fillSEDV_rad();
		void		solveIsed();
		void		fillREDZ();
		void		fillREDV();

		void		fillREDV_rad();
		void		REDZKernel(int &_x_m, int &_y_m, int &_x_n, int &_y_n);
		void		REDVKernel(int &_x, int &_y, int &_id);
		void		REDVKernel_rad(int &_X, int &_y, int &_id);
		void		fillSVDI();
		void		constructI();
		value_t		testRCS(const VectorR3& sca_k);
		Complex		SEDDVKernel(VectorR3 *vp4, VectorR3 *vm4, VectorR3 &vp, VectorR3 &vm, int &_nvm, VectorR3 &inc_k, VectorR3 &inc_e);
		Complex		SEDMVKernel(VectorR3 *vp3, VectorR3 *vm3, VectorR3 &vp, VectorR3 &vm, VectorR3 &inc_k, VectorR3 &inc_e);

	protected:
		void		Task_construct();
		void		Task_assign(size_t task_num);
		void		fillSEDZ_MP();
		void		SEDZ_MP_Kernel(int begin, int end);
		void		fillREDZ_MP();
		void		REDZ_MP_Kernel(int begin, int end);
		void		fillREDV_MP();
		void		REDV_MP_Kernel(int begin, int end);
		//value_t     getBiRCS_MP(const VectorR3& sca_k);
		//void        getEFiled_MP(component::FieldData* pdata, const VectorR3& rad_k, const VectorR3& rad_ev, const VectorR3& rad_eh) const;
		void		fillSEDZDD_MP(int &_f_sed, int &_s_sed, int &_start_tf, int &_end_tf, int &_start_ts, int &_end_ts, int &_start_ef, int &_end_ef, int &_start_es, int &_end_es);

		void		fillCOUPZDD_MP(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup);
		void		fillCOUPZDM_MP(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup);
		void		fillCOUPZMD_MP(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup);
		void		fillCOUPZMM_MP(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup);

		void		fillCOUPZDD_MP_C(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup);
		void		fillCOUPZDM_MP_C(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup);
		void		fillCOUPZMD_MP_C(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup);
		void		fillCOUPZMM_MP_C(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup);

		void		fillCOUPZDD_MP_B(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup);
		void		fillCOUPZDM_MP_B(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup);
		void		fillCOUPZMD_MP_B(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup);
		void		fillCOUPZMM_MP_B(int &_x_m, int &_y_m, int &_x_n, int &_y_n, int &_k_m, int &_k_n, CMatrix& Z_coup);

		CMatrix		coupZ_MP(int &_x_m, int &_y_m, int &_x_n, int &_y_n);
		CMatrix		coupZ_MP_COM(int &_x_m, int &_y_m, int &_x_n, int &_y_n);

	protected:
		void		getJCM();
	protected:
		void		buildblk();
		void		findnearblk(Blk& _blk);
		void        preCalculateKArray();
		void        preCalculateTrans();
		void		preCalculateTrans_v2();
		void		CalculateTransKernel(int begin, int end);
		void        preCalculateTrans(int &_x, int &_y,CVector& _transfer) const;
		void        preCalculateRadAndRecv();
		void		preCalculateRadAndRecv_v2();
		void		CalculateNearZ();
		void		CalculateNearZ_v2();
		void		preCalculateNearZ();
		void		fillprenearZ_kernel(int &_idx, int &_idy);
		void		fillnearZ_kernel(int &_idx, int &_idy);
		void		fillnearZ_kernelMP(int begin, int end);
		void		fastfillnearZ_kernelMP(int begin, int end);
		void		fillnearZ_byblkarr(int begin, int end);
		void        prepareSolve();
		void		calculatePreconditioner();
		CVector     matrixVectorMultiply(const int &_x,const CVector& b,const CVector&_transfer) const;
		CVector     matrixVectorMultiply(const CVector& b);
		CVector		matrixVectorMultiply_v2(const CVector& b);
		void		fillsk_kernel(int begin, int end,CVector &b);
		void		fillgk_kernel(int begin, int end);
		void		fillfarZb_kernel(int begin, int end);
		void		fillUfarZb_kernel(int begin, int end);
		int			iterativeSolve_BICG();
		int			iterativeSolve_CGN();
		int			iterativeSolve_GMRES();
		void		fillZbyFMM();
		
	private:
		VectorC3    radiationFunction(const VectorR3* vpgs3, const VectorR3* vmgs3, const VectorR3& vp, 
									  const VectorR3& vm, const VectorR3& center, const VectorR3& sam_k);
		VectorC3    halfRecvFunction(const VectorR3* vpgs3, const VectorR3* vmgs3, const VectorR3& vp,
			const VectorR3& vm, const VectorR3& nml_plu, const VectorR3& nml_min,
			const VectorR3& center, const VectorR3& sam_k);
		VectorR3	dipoleFunction(VectorR3 *v_p3, VectorR3 *v_m3, value_t &len);
		Complex		dpcouplFunction(VectorC3& dp_f, VectorC3& dp_s, VectorR3& _R);

	protected:
		void		shiftcoordv4(int &sedid, VectorR3 *v4) const;
		void		shiftcoordv1(int &sedid, VectorR3 &v) const;
		void		shiftcoordv3(int &sedid, VectorR3 *v3) const;

		void		shiftcoordv4(int &xid, int &yid, VectorR3 *v4) const;
		void		shiftcoordv3(int &xid, int &yid, VectorR3 *v3) const;
		void		shiftcoordv1(int &xid, int &yid, VectorR3 &v) const;

		void		shiftbedgecoordv4(int &bedge, VectorR3 *v4) const;
		void		shiftbedgecoordv3(int &bedge, VectorR3 *v3) const;
		void		shiftbedgecoordv1(int &bedge, VectorR3 &v) const;

		void		shiftcoordv4(int &sedid, VectorR3 *v4, value_t &scale) const;
		void		shiftcoordv1(int &sedid, VectorR3 &v, value_t &scale) const;
		void		shiftcoordv3(int &sedid, VectorR3 *v3, value_t &scale) const;

		void		shiftcoordv4(int &xid, int &yid, VectorR3 *v4, value_t &scale) const;
		void		shiftcoordv3(int &xid, int &yid, VectorR3 *v3, value_t &scale) const;
		void		shiftcoordv1(int &xid, int &yid, VectorR3 &v, value_t &scale) const;

		void		scalecoordv1(int &xid, int &yid, VectorR3 &v, value_t &scale) const;
		void		scalecoordv3(int &xid, int &yid, VectorR3 *v3, value_t &scale) const;

		void		shiftbedgecoordv4(int &bedge, VectorR3 *v4, value_t &scale) const;
		void		shiftbedgecoordv3(int &bedge, VectorR3 *v3, value_t &scale) const;
		void		shiftbedgecoordv1(int &bedge, VectorR3 &v, value_t &scale) const;

		VectorR3	centerv4(VectorR3 *v4);
		VectorR3	centerv3(VectorR3 *v3);
		int			ArrayToSED(int _x, int _y) const;
		int			ArrayToNear(int _x, int _y) const;
		int			NearArrayToSED(int _x, int _y) const;
		void        getBoundaryNum(int &_sedid, int &_b1, int &_b2) const;
	protected:
		int			tet_num;
		int			max_iter_num_;
		int			unk_t[9];
		int			unk_e[9];
		int			svd_k;
		value_t		iter_threshold;
		value_t		svd_threshold;
		value_t		threshold_edm;
		Complex		Sigma;

		int			grid_x;
		int			grid_y;
		int			grid_z;
		value_t		d_x;
		value_t		d_y;
		value_t		d_z;
		VectorR3	min_boundary, max_boundary;
		VectorR3	cen_box;
		std::vector<Complex> PGFgrid;

		int			isContinuous;
		value_t		Dx;
		value_t		Dy;
		value_t		Phase_0;
		value_t		Phase_x;
		value_t		Phase_y;
		int			Array_x;
		int			Array_y;
		int			Array_num;
		value_t		scale_x;
		value_t		scale_y;
		int			t_sum;

		int         unknowns_;
		int			unknowns_t;
		int			unknowns_e;
		int			com_ce;
		int			com_ct;
		int         bou_ce1;
		int         bou_ce2;

		value_t		mu1_, eps1_;
		value_t     k_;
		value_t		omiga;
		Qcx_mat     Z;
		Qcx_mat		Z_sed;
		CMatrix		Z_sed_eigen;
		Qcx_mat		Z_coup;
		CMatrix		Z_coup_eigen;
		Qcx_mat		Z_red;
		CMatrix		Z_red_eigen;
		Qcx_vec		I;
		MatVec      I_tot;
		MatArrayEi	I_tot_eigen;
		Qcx_mat		I_sed;
		CMatrix		I_sed_eigen;
		CMatrix		I_sed_eigen_rad;
		MatArray	I_svd;
		MatArrayEi	I_svd_eigen;
		Qcx_vec		I_red;
		CVector		I_red_eigen;
		Qcx_vec     V;
		Qcx_mat		V_sed;
		CMatrix		V_sed_eigen;
		CMatrix		V_sed_eigen_rad;
		Qcx_vec		V_red;
		CVector		V_red_eigen;
		Qcx_vec		V_temp;
		CVector		V_temp_eigen;
		CVector		I_feed;

		Qcx_mat		Z2K1, Z2K4;
		Qc3_mat		Z2K2, Z2K3;

		Qcx_mat		Z2EK1, Z2EK4;
		Qc3_mat		Z2EK2, Z2EK3;

		Qcx_mat		Z5K;

		Qcx_vec		RIVP, RRIV, IV;
		Vec_R3		RIV, IVP;

		RadVector   exc_edge_;
		component::MultipleIncidence multiInc;
		VectorR3	PW_inc_k;
		VectorR3	PW_inc_e;

		MeshPointer mesh_ptr_;
		MeshPointer feko_mesh_ptr_;
		CEPointer   ce_ptr_;
		CTPointer   ct_ptr_;

		mutable std::mutex mutex_;
		std::vector<int> Task_assign1, Task_assign2, Task_assign3, rad_rec_assign,trans_assign,dp_assign;
		component::ThreadPool pool_;
		std::vector<std::pair<int, int>> SED_task1, SED_task2,rad_rec_task,trans_task,dp_task;
		std::vector<int> REDV_task;

		int		ex_num;
		std::vector<VectorR3> Centercoord;
		std::vector<Complex> J_center;

		std::vector<Blk> BlockArray;
		int         sam_theta_;
		int         sam_phi_;
		Qr3_mat      arr_k_;
		bool            preconditioning_;
		//RadVector       exc_edge_;
		CVector         coeff_;
		CMatrix			coeffMatrix_;
		//C3Matrix        sk_, gk_;
		VC3Array		sk_, gk_, dp_,dpk_;
		SpCMat		nearZ_, nearZ_0inv, nearZ_1,nearZ_01;
		CMatrix		farZ_, preNearZ_,preLU_;
		CVector		farZb_,UfarZb_;
		std::vector<LocSpE>	Loc_near_ele_,Loc_nearZ1_ele_;
		//std::vector<C3Matrix>   rad_, rec_;
		//std::vector<C3Matrix>	sed_rad_, sed_rec_;
		std::vector<VC3Array> rad_, rec_;
		std::vector<VC3Array> sed_rad_, sed_rec_;
		math::GaussLegendre     gl_theta_, gl_phi_;
		TransArray		transfer_;
	};

	inline VectorR3 SED_CM_FMM::centerv4(VectorR3 *v4)
	{
		return (v4[0] + v4[1] + v4[2] + v4[3]) / 4.0f;
	}

	inline VectorR3 SED_CM_FMM::centerv3(VectorR3 *v3)
	{
		return (v3[0] + v3[1] + v3[2]) / 3.0f;
	}
} // namespace mom





