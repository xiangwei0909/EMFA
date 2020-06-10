#include "Quadrature.h"

void math::Gauss1Point(const VectorR3& v1, const VectorR3& v2, const VectorR3& v3, VectorR3* p)
{
    p[0] = g1*(v1 + v2 + v3);
}

void math::Gauss1Point(const VectorR3& v1, const VectorR3& v2, const VectorR3& v3, const VectorR3& v4, VectorR3* p)
{
	p[0] = (v1 + v2 + v3 + v4)/4.0f;
}

void math::Gauss3Point(const VectorR3& v1, const VectorR3& v2, const VectorR3& v3, VectorR3* p)
{
    p[0] = g3_1[0] * v1 + g3_1[1] * v2 + g3_1[2] * v3;
    p[1] = g3_2[0] * v1 + g3_2[1] * v2 + g3_2[2] * v3;
    p[2] = g3_3[0] * v1 + g3_3[1] * v2 + g3_3[2] * v3;
}

void math::Gauss4Point(const VectorR3& v1, const VectorR3& v2, const VectorR3& v3,const VectorR3& v4, VectorR3* p)
{
	p[0] = gt4_1[0] * v1 + gt4_1[1] * v2 + gt4_1[2] * v3 + gt4_1[3] * v4;
	p[1] = gt4_2[0] * v1 + gt4_2[1] * v2 + gt4_2[2] * v3 + gt4_2[3] * v4;
	p[2] = gt4_3[0] * v1 + gt4_3[1] * v2 + gt4_3[2] * v3 + gt4_3[3] * v4;
	p[3] = gt4_4[0] * v1 + gt4_4[1] * v2 + gt4_4[2] * v3 + gt4_4[3] * v4;
}

void math::Gauss4Point(const VectorR3& v1, const VectorR3& v2, const VectorR3& v3, VectorR3* p)
{
    p[0] = g4_1[0] * v1 + g4_1[1] * v2 + g4_1[2] * v3;
    p[1] = g4_2[0] * v1 + g4_2[1] * v2 + g4_2[2] * v3;
    p[2] = g4_3[0] * v1 + g4_3[1] * v2 + g4_3[2] * v3;
}

void math::Gauss5Point(VectorR3& v1, VectorR3& v2, VectorR3& v3, VectorR3& v4, VectorR3* p)
{
	p[0] = gt5_1[0] * v1 + gt5_1[1] * v2 + gt5_1[2] * v3 + gt5_1[3] * v4;
	p[1] = gt5_2[0] * v1 + gt5_2[1] * v2 + gt5_2[2] * v3 + gt5_2[3] * v4;
	p[2] = gt5_3[0] * v1 + gt5_3[1] * v2 + gt5_3[2] * v3 + gt5_3[3] * v4;
	p[3] = gt5_4[0] * v1 + gt5_4[1] * v2 + gt5_4[2] * v3 + gt5_4[3] * v4;
	p[4] = gt5_5[0] * v1 + gt5_5[1] * v2 + gt5_5[2] * v3 + gt5_5[3] * v4;
}

void math::Gauss7Point(const VectorR3& v1, const VectorR3& v2, const VectorR3& v3, VectorR3* p)
{
    p[0] = g7_1[0] * v1 + g7_1[1] * v2 + (1 - g7_1[0] - g7_1[1])*v3;
    p[1] = g7_2[0] * v1 + g7_2[1] * v2 + (1 - g7_2[0] - g7_2[1])*v3;
    p[2] = g7_3[0] * v1 + g7_3[1] * v2 + (1 - g7_3[0] - g7_3[1])*v3;
    p[3] = g7_4[0] * v1 + g7_4[1] * v2 + (1 - g7_4[0] - g7_4[1])*v3;
    p[4] = g7_5[0] * v1 + g7_5[1] * v2 + (1 - g7_5[0] - g7_5[1])*v3;
    p[5] = g7_6[0] * v1 + g7_6[1] * v2 + (1 - g7_6[0] - g7_6[1])*v3;
    p[6] = g7_7[0] * v1 + g7_7[1] * v2 + (1 - g7_7[0] - g7_7[1])*v3;
}

void math::GaussLine3Point(const VectorR3 & v1, const VectorR3 & v2, VectorR3 * p)
{
    auto v1v2 = v2 - v1;
    for (int g = 0; g < 3; ++g)
        p[g] = v1 + gl3[g] * v1v2;
}

/////////////////////////////////////////////////////////////////////////////////////////

math::GaussLegendre::Evaluation::Evaluation(value_t x, int n)
: _x(x), _v(1), _d(0) 
{
    this->evaluate(x, n);
}

void math::GaussLegendre::Evaluation::evaluate(value_t x, int n)
{
    this->_x = x;

    value_t vsub1 = x;
    value_t vsub2 = 1;
    value_t f = 1 / (x * x - 1);

    for (int i = 2; i <= n; ++i)
    {
        this->_v = ((2 * i - 1) * x * vsub1 - (i - 1) * vsub2) / i;
        this->_d = i * f * (x * this->_v - vsub1);

        vsub2 = vsub1;
        vsub1 = this->_v;
    }
}

math::GaussLegendre::GaussLegendre()
{

}

void math::GaussLegendre::set(int n)
{
    _r.resize(n);
    _w.resize(n);

    // Solve roots and weights
    for (int i = 0; i < n; ++i)
    {
        value_t dr = 1;

        Evaluation eval(cos(component::PI * (i + 0.75f) / (n + 0.5f)), n);
        do {
            dr = eval.v() / eval.d();
            eval.evaluate(eval.x() - dr, n);
        } while (fabs(dr) > 2e-7f);

        this->_r[n - 1 - i] = eval.x();
        this->_w[n - 1 - i] = 2 / ((1 - eval.x() * eval.x()) * eval.d() * eval.d());
    }
}

void math::GaussLegendre::printRootsAndWeights() const
{
    std::cout << std::setprecision(10) << "    weights          root" << '\n';
    for (size_t i = 0; i < _r.size(); ++i)
        std::cout << i + 1 << ": " << _w[i] << "    " << _r[i] << '\n';
}

