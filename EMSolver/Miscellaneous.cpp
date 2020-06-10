#include "Miscellaneous.h"

using namespace component;
using namespace math;

Complex mom::integral::cfieZKernel(value_t alpha, value_t k, const VectorR3 * vf3, const VectorR3 * vs3, const VectorR3 & vf,
                                   const VectorR3 & vs, const VectorR3 & nf, value_t fsign, value_t ssign)
{
    VectorR3 vgf[3], vgs[3];
    math::Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
    math::Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

    Complex Ze(0, 0), Zm(0, 0), G(0, 0);
    VectorR3 rou_f, rou_s;
    for (int p = 0; p < 3; p++)
    {
        rou_f = fsign * (vgf[p] - vf);
        for (int q = 0; q < 3; q++)
        {
            rou_s = ssign * (vgs[q] - vs);
            value_t r = (vgf[p] - vgs[q]).Norm();
            G = exp(-J0 * k * r) / r;
            Ze += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - (fsign * ssign) / (k * k)) * G;
            Zm += w3[p] * w3[q] * (rou_f ^ (nf * ((vgf[p] - vgs[q]) * rou_s))) * (1.0f + J0 * k * r) * G / (r * r);
        }
    }
    Ze = (J0 * k * Z0 / PI4) * Ze;
    Zm = Zm / (16 * PI);
    return alpha * Ze + (1 - alpha) * Z0 * Zm;
}

Complex mom::integral::cfieZSingular(value_t alpha, value_t k, const VectorR3 * vf3, const VectorR3 *vs3, const VectorR3 & vf, 
                                     const VectorR3 & vs, value_t fsign, value_t ssign)
{
    VectorR3 vgf[3], vgs[3];
    Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
    Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);
    VectorR3 rou_f, rou_s;
    Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Zm(0, 0);
    value_t area = Area(vf3[0], vf3[1], vf3[2]);

    for (int p = 0; p < 3; p++)
    {
        rou_f = fsign * (vgf[p] - vf);
        for (int q = 0; q < 3; q++)
        {
            rou_s = ssign * (vgs[q] - vs); 
            value_t r = (vgf[p] - vgs[q]).Norm();
            Complex coff = -0.5f * k * k * r + J0 * k / 6.0f * (k * k * r * r - 6);
            Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - (fsign * ssign) / (k * k)) * coff;
        }

        for (int i = 0; i < 3; i++)
        {
            VectorR3 Li = vs3[(i + 1) % 3] - vs3[i];
            Li.Normalize();

            value_t Lpi = (vs3[(i + 1) % 3] - vgf[p]) ^ Li; //(边正顶点 - 场点)·(边向量)
            value_t Lmi = (vs3[i] - vgf[p]) ^ Li; //(边负顶点 - 场点)·(边向量)
            value_t Rpi = (vs3[(i + 1) % 3] - vgf[p]).Norm(); //边正顶点 到 场点 的距离
            value_t Rmi = (vs3[i % 3] - vgf[p]).Norm(); //边负顶点 到 场点 的距离
            VectorR3 sj0 = vs3[i] + Li * ((vgf[p] - vs3[i]) ^ Li); //垂足在u方向位置矢量
            value_t P0i = (sj0 - vgf[p]).Norm(); //垂足 到 场点 的距离
            if (P0i > 1e-10f)
                Ze2 += w3[p] * (0.25f * (rou_f ^ rou_s) - (fsign * ssign) / (k * k)) * P0i * log((Rpi + Lpi) / (Rmi + Lmi)); 
        }

        rou_f = fsign * (vgf[p] - vf);
        rou_s = ssign * (vgs[p] - vs);
        Zm += w3[p] * (rou_f ^ rou_s);
    }

    Ze = (J0 * k * Z0 / PI4) * (Ze1 + 1.0f / area * Ze2);
    Zm = Zm / (8.0f * area);

    return alpha * Ze + (1 - alpha) * Z0 * Zm;
}

Complex mom::integral::cfieZKernel(value_t k, const VectorR3 * vf3, const VectorR3 * vs3, const VectorR3 & vf,
	const VectorR3 & vs, const VectorR3 & nf, value_t fsign, value_t ssign,value_t msign,value_t dsign)
{
	VectorR3 vgf[3], vgs[3];
	math::Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	math::Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);

	Complex Ze(0, 0), Zm(0, 0), G(0, 0);
	VectorR3 rou_f, rou_s;
	for (int p = 0; p < 3; p++)
	{
		rou_f = fsign * (vgf[p] - vf);
		for (int q = 0; q < 3; q++)
		{
			rou_s = ssign * (vgs[q] - vs);
			value_t r = (vgf[p] - vgs[q]).Norm();
			G = exp(-J0 * k * r) / r;
			Ze += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - (fsign * ssign) / (k * k)) * G;
			Zm += w3[p] * w3[q] * (rou_f ^ (nf * ((vgf[p] - vgs[q]) * rou_s))) * (1.0f + J0 * k * r) * G / (r * r);
		}
	}
	Ze = (J0 * k * Z0 / PI4) * Ze;
	Zm = Zm / (16 * PI);
	return msign*Ze + dsign*Z0 * Zm;
}

Complex mom::integral::cfieZSingular(value_t k, const VectorR3 * vf3, const VectorR3 *vs3, const VectorR3 & vf,
	const VectorR3 & vs, value_t fsign, value_t ssign, value_t msign, value_t dsign)
{
	VectorR3 vgf[3], vgs[3];
	Gauss3Point(vf3[0], vf3[1], vf3[2], vgf);
	Gauss3Point(vs3[0], vs3[1], vs3[2], vgs);
	VectorR3 rou_f, rou_s;
	Complex Ze(0, 0), Ze1(0, 0), Ze2(0, 0), Zm(0, 0);
	value_t area = Area(vf3[0], vf3[1], vf3[2]);

	for (int p = 0; p < 3; p++)
	{
		rou_f = fsign * (vgf[p] - vf);
		for (int q = 0; q < 3; q++)
		{
			rou_s = ssign * (vgs[q] - vs);
			value_t r = (vgf[p] - vgs[q]).Norm();
			Complex coff = -0.5f * k * k * r + J0 * k / 6.0f * (k * k * r * r - 6);
			Ze1 += w3[p] * w3[q] * (0.25f * (rou_f ^ rou_s) - (fsign * ssign) / (k * k)) * coff;
		}

		for (int i = 0; i < 3; i++)
		{
			VectorR3 Li = vs3[(i + 1) % 3] - vs3[i];
			Li.Normalize();

			value_t Lpi = (vs3[(i + 1) % 3] - vgf[p]) ^ Li; //(边正顶点 - 场点)·(边向量)
			value_t Lmi = (vs3[i] - vgf[p]) ^ Li; //(边负顶点 - 场点)·(边向量)
			value_t Rpi = (vs3[(i + 1) % 3] - vgf[p]).Norm(); //边正顶点 到 场点 的距离
			value_t Rmi = (vs3[i % 3] - vgf[p]).Norm(); //边负顶点 到 场点 的距离
			VectorR3 sj0 = vs3[i] + Li * ((vgf[p] - vs3[i]) ^ Li); //垂足在u方向位置矢量
			value_t P0i = (sj0 - vgf[p]).Norm(); //垂足 到 场点 的距离
			if (P0i > 1e-10f)
				Ze2 += w3[p] * (0.25f * (rou_f ^ rou_s) - (fsign * ssign) / (k * k)) * P0i * log((Rpi + Lpi) / (Rmi + Lmi));
		}

		rou_f = fsign * (vgf[p] - vf);
		rou_s = ssign * (vgs[p] - vs);
		Zm += w3[p] * (rou_f ^ rou_s);
	}

	Ze = (J0 * k * Z0 / PI4) * (Ze1 + 1.0f / area * Ze2);
	Zm = Zm / (8.0f * area);

	return msign*Ze + dsign * Z0 * Zm;
}

Complex mom::integral::cfieVKernel(value_t alpha, value_t k, const VectorR3 & vec_k, const VectorR3 & vec_ei, const VectorR3 & vec_hi, 
                                    const VectorR3 * vp3, const VectorR3 * vm3, const VectorR3 & np, const VectorR3 & nm, 
                                    const VectorR3 & vp, const VectorR3 & vm)
{
    VectorR3 vgp[3], vgm[3];
    VectorR3 rou_p, rou_m; 
    Complex Ve, Vm;

    Gauss3Point(vp3[0], vp3[1], vp3[2], vgp);
    Gauss3Point(vm3[0], vm3[1], vm3[2], vgm);

    VectorC3 Ei(0, 0, 0), Hi(0, 0, 0);
    Complex Vep(0, 0), Vem(0, 0), Vmp(0, 0), Vmm(0, 0), G(0, 0);

    for (int a = 0; a < 3; a++)
    {
        rou_p = vgp[a] - vp;
        G = exp(-J0 * k * (vgp[a] ^ vec_k));

        Ei = G * vec_ei;
        Vep += w3[a] * (rou_p ^ Ei);
        Hi = G * vec_hi;    //  include Z0
        Vmp += w3[a] * (rou_p ^ (np * Hi));

        rou_m = vm - vgm[a];
        G = exp(-J0 * k * (vgm[a] ^ vec_k));

        Ei = G * vec_ei;
        Vem += w3[a] * (rou_m ^ Ei);
        Hi = G * vec_hi;    // include Z0
        Vmm += w3[a] * (rou_m ^ (nm * Hi));
    }

    Ve = (Vep + Vem);
    Vm = (Vmp + Vmm);

    return alpha * Ve + (1 - alpha) * Vm;   // Note: multiple (0.5 * length)
}
