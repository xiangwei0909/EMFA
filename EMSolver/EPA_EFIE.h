//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "EM.h"

namespace mom {

	class EPA_EFIE :public EM {
		using MeshPointer = std::shared_ptr<component::Mesh>;
		using CEPointer = std::shared_ptr<component::CommonEdge>;
		using RadVector = std::vector<std::pair<int, int>>;
	public:
		EPA_EFIE();
		~EPA_EFIE();
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

		Complex		ZuKernel(VectorR3 *vf3, VectorR3 *vs3, VectorR3 &vf, VectorR3 &vs);
		Complex		ZioEJmKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs);
		Complex		ZoiHJeKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 &nf);
		Complex		ZoiEJeKernel(VectorR3 * vf3, VectorR3 * vs3, VectorR3 & vf, VectorR3 & vs, VectorR3 &nf);

		Complex     eVKernel(VectorR3 *vp3, VectorR3 *vm3, VectorR3 &vp, VectorR3 &vm);
		void		heVKernel(Complex &vh,Complex &ve,VectorR3 *vp3, VectorR3 *vm3, VectorR3 &vp, VectorR3 &vm, VectorR3 &nplu, VectorR3&nmin);

		Complex		EDMKernel(VectorR3 *vf_p3, VectorR3 *vf_m3, VectorR3 *vs_p3, VectorR3 *vs_m3, value_t &l_fld, value_t &l_src);
		VectorR3	center(VectorR3 *v3);
	protected:
		bool        readExcEdges(const Qstring& rad_file, Qstring& stateInfo);
		void        readExcEdges(const Qstring& rad_file);
		void        fillZ();
		void        fillV();
		void        radiateV();
		value_t     getBiRCS(const VectorR3& sca_k) const;
		value_t     getEFIEBiRCS(const VectorR3& sca_k) const;
		void        getEFiled(component::FieldData* pdata, const VectorR3& rad_k, const VectorR3& rad_ev, const VectorR3& rad_eh) const;
		void        getCurrentOnFeed();

		bool        writeZIVData();
		void        calculateSurfaceCurrent(std::vector<component::CurrentData> *currents) const;

		void		fillZ_EPA();
		void		fillZ_u();
		void		fillZ_efie();
		void		fillZ_io();
		void		fillZ_oi();
		void		fillV_inc();
	protected:
		int         unknowns_;
		int			unknowns_epf;
		value_t		threshold_edm;
		value_t     k_;
		Qcx_mat     Z;
		Qcx_vec     I;
		Qcx_vec     V;
		CMatrix		Z_i;
		CMatrix		Z_u;
		CMatrix		Z_oi;
		CMatrix		Z_io;
		CMatrix		Z_efie;
		CMatrix		Z_final;
		CVector		V_inc;
		CVector		I_sca;
		CVector		I_test;
		CVector		V_test;
		Complex		Sigma;
		RadVector   exc_edge_;


		MeshPointer mesh_ptr_;
		MeshPointer epf_ptr_;
		CEPointer   ce_ptr_;
		CEPointer	epf_ce_;
		Qstring		dir_sur_;
	};
	inline VectorR3 EPA_EFIE::center(VectorR3 *v3)
	{
		return ((v3[0] + v3[1] + v3[2]) / 3.0f);
	}
} // namespace mom
