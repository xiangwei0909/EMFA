#pragma once
#include "EM.h"

namespace mom {

class FSPGF_EFIE :public EM {
	using MeshPointer = std::shared_ptr<component::Mesh>;
	using CEPointer = std::shared_ptr<component::CommonEdge>;
	using RadVector = std::vector<std::pair<int, int>>;
public:
	FSPGF_EFIE();
	~FSPGF_EFIE();
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

	Complex		FSPGF(VectorR3 &r, int _t_sum, value_t _Dx, value_t _Dy);
	Complex		FSPGFSingular(VectorR3 &r, int _t_sum, value_t _Dx, value_t _Dy);

	VectorR3	Translateoff(VectorR3 &r, int _x, int _y, int _z) const;
	Complex		erfz(Complex z);
	Complex		erfc(Complex z) {return 1.0f - erfz(z); }
protected:
	void		prepareFSPGF(value_t _Dx, value_t _Dy, int _Array_x, int _Array_y, int _t_sum);
	bool        readExcEdges(const Qstring& rad_file, Qstring& stateInfo);
	void        readExcEdges(const Qstring& rad_file);
	void        fillZ();
	void        fillV();
	void        radiateV();
	value_t     getBiRCS(const VectorR3& sca_k) const;
	void        getEFiled(component::FieldData* pdata, const VectorR3& rad_k, const VectorR3& rad_ev, const VectorR3& rad_eh) const;

	void        calculateSurfaceCurrent(std::vector<component::CurrentData> *currents) const;
	bool        writeZIVData();
protected:
	int         unknowns_;
	int			Array_x;
	int			Array_y;
	int			t_sum;
	int			isContinuous;
	value_t     k_;
	value_t		Dx;
	value_t		Dy;
	Complex     Sigma;
	Qcx_mat     Z;
	Qcx_vec     I;
	Qcx_vec     V;

	RadVector   exc_edge_;

	MeshPointer mesh_ptr_;
	CEPointer   ce_ptr_;
};

} // namespace mom