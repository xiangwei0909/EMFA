//**********************************************************
// Author: Xiang Wei
// License: MIT
//**********************************************************

#pragma once
#include "EM.h"

namespace mom {

	class TDS :public EM {
		using MeshPointer = std::shared_ptr<component::Mesh>;
		using CEPointer = std::shared_ptr<component::CommonEdge>;
		using RadVector = std::vector<std::pair<int, int>>;
	public:
		TDS();
		~TDS();
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

		Complex		TTZppKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &f_nor, VectorR3 &s_nor, VectorR3 *vf2, VectorR3 *vs2, value_t &eps_p, value_t &eps_m, int &bedge_f, int &bedge_s);
		Complex		TTZpmKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &f_nor, VectorR3 &s_nor, VectorR3 *vf2, VectorR3 *vs2, value_t &eps_p, value_t &eps_m, int &bedge_f, int &bedge_s);
		Complex		TTZmpKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &f_nor, VectorR3 &s_nor, VectorR3 *vf2, VectorR3 *vs2, value_t &eps_p, value_t &eps_m, int &bedge_f, int &bedge_s);
		Complex		TTZmmKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &f_nor, VectorR3 &s_nor, VectorR3 *vf2, VectorR3 *vs2, value_t &eps_p, value_t &eps_m, int &bedge_f, int &bedge_s);

		Complex		TTZppSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &f_nor, VectorR3 &s_nor, VectorR3 *vf2, VectorR3 *vs2, value_t &eps_p, value_t &eps_m, int &bedge_f, int &bedge_s);
		Complex		TTZpmSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &f_nor, VectorR3 &s_nor, VectorR3 *vf2, VectorR3 *vs2, value_t &eps_p, value_t &eps_m, int &bedge_f, int &bedge_s);
		Complex		TTZmpSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &f_nor, VectorR3 &s_nor, VectorR3 *vf2, VectorR3 *vs2, value_t &eps_p, value_t &eps_m, int &bedge_f, int &bedge_s);
		Complex		TTZmmSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs, VectorR3 &f_nor, VectorR3 &s_nor, VectorR3 *vf2, VectorR3 *vs2, value_t &eps_p, value_t &eps_m, int &bedge_f, int &bedge_s);

		Complex		TNZppKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &s_nor, value_t &eps);
		Complex		TNZpmKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &s_nor, value_t &eps);
		Complex		TNZmpKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &s_nor, value_t &eps);
		Complex		TNZmmKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &s_nor, value_t &eps);

		Complex		TNZppSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &s_nor, value_t &eps);
		Complex		TNZpmSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &s_nor, value_t &eps);
		Complex		TNZmpSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &s_nor, value_t &eps);
		Complex		TNZmmSingular(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &s_nor, value_t &eps);

		Complex     eVKernel(VectorR3 *vp3, VectorR3 *vm3, VectorR3 &vp, VectorR3 &vm);
	protected:
		bool        readExcEdges(const Qstring& rad_file, Qstring& stateInfo);
		void        readExcEdges(const Qstring& rad_file);
		void        fillZ();
		void		fillZTT();
		void		fillZTN();
		void		fillZNT();
		void		fillZNN();
		void        fillV();
		void        radiateV();
		value_t     getBiRCS(const VectorR3& sca_k) const;
		void        getEFiled(component::FieldData* pdata, const VectorR3& rad_k, const VectorR3& rad_ev, const VectorR3& rad_eh) const;

		void		TransTri(VectorR3 *v3_ori, VectorR3 *v3_end, VectorR3& Nor, value_t wgt);
		void		TransLine(VectorR3 *v2_ori, VectorR3 *v2_end, VectorR3& Nor, value_t wgt);
		void		TransPoi(VectorR3 &v_ori, VectorR3 &v_end, VectorR3 &Nor, value_t wgt);
		value_t		IsSingular(VectorR3 *vs3, VectorR3 &vf);
		VectorR3    IspSingular(VectorR3 * vs3, VectorR3 & vf);
		bool        writeZIVData();
		void        calculateSurfaceCurrent(std::vector<component::CurrentData> *currents) const;
	protected:
		int         unknowns_;
		int			unknowns_n, unknowns_t;
		value_t     k_;
		Qcx_mat     Z;
		Qcx_vec     I;
		Qcx_vec     V;
		Complex		Sigma;
		value_t		omiga;
		value_t		Thickness;
		RadVector   exc_edge_;

		MeshPointer mesh_ptr_;
		CEPointer   ce_ptr_;
	};

} // namespace mom