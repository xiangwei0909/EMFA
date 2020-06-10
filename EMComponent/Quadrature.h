//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "VectorR3.h"
#include "Mathdef.h"

namespace math {

using component::VectorR3;

const value_t w1 = 1.0f;    //一点高斯采样 权值
const value_t g1 = 1.0f / 3.0f; //一点高斯采样 系数

const value_t gl3[3] = { 0.1127016654f, 0.5f, 0.8872983346f }; // 一维高斯积分系数 [0, 1]
const value_t glw3[3] = { 0.2777777777775f, 0.444444444444445f, 0.2777777777775f }; // 一维高斯积分权值[0, 1]

const value_t w3[3] = {
    1.0f / 3.0f,
    1.0f / 3.0f,
    1.0f / 3.0f
};      //三点高斯采样 权值
const value_t g3_1[3] = { 2.0f / 3.0f, 1.0f / 6.0f, 1.0f / 6.0f };  //三点高斯采样 系数
const value_t g3_2[3] = { 1.0f / 6.0f, 2.0f / 3.0f, 1.0f / 6.0f };  //三点高斯采样 系数
const value_t g3_3[3] = { 1.0f / 6.0f, 1.0f / 6.0f, 2.0f / 3.0f };  //三点高斯采样 系数


const value_t w4[4] = {
    -0.562500f,
    0.5208333333f,
    0.5208333333f,
    0.5208333333f
};      //三点高斯采样 权值
const value_t g4_1[3] = { 1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f };  //三点高斯采样 系数
const value_t g4_2[3] = { 0.6f, 0.2f, 0.2f };   //三点高斯采样 系数
const value_t g4_3[3] = { 0.2f, 0.6f, 0.2f };   //三点高斯采样 系数
const value_t g4_4[3] = { 0.2f, 0.2f, 0.6f };   //三点高斯采样 系数


const value_t w7[7] = {
    0.112500000000000f * 2,
    0.062969590272413f * 2,
    0.062969590272413f * 2,
    0.062969590272413f * 2,
    0.066197076394253f * 2,
    0.066197076394253f * 2,
    0.066197076394253f * 2
};

static const float gt5_1[4] = { 1.0f / 4.0f,1.0f / 4.0f,1.0f / 4.0f,1.0f / 4.0f };	//四点高斯采样 系数
static const float gt5_2[4] = { 1.0f / 2.0f,1.0f / 6.0f,1.0f / 6.0f,1.0f / 6.0f };	//四点高斯采样 系数
static const float gt5_3[4] = { 1.0f / 6.0f,1.0f / 2.0f,1.0f / 6.0f,1.0f / 6.0f };	//四点高斯采样 系数
static const float gt5_4[4] = { 1.0f / 6.0f,1.0f / 6.0f,1.0f / 2.0f,1.0f / 6.0f };	//四点高斯采样 系数
static const float gt5_5[4] = { 1.0f / 6.0f,1.0f / 6.0f,1.0f / 6.0f,1.0f / 2.0f };	//四点高斯采样 系数

static const float wt5[5] = {
	-4.0f / 5.0f,
	9.0f / 20.0f,
	9.0f / 20.0f,
	9.0f / 20.0f,
	9.0f / 20.0f
};

const value_t g7_1[2] = { 0.333333333333333f, 0.333333333333333f };
const value_t g7_2[2] = { 0.797426985353087f, 0.101286507323456f };
const value_t g7_3[2] = { 0.101286507323456f, 0.797426985353087f };
const value_t g7_4[2] = { 0.101286507323456f, 0.101286507323456f };
const value_t g7_5[2] = { 0.470142064105115f, 0.470142064105115f };
const value_t g7_6[2] = { 0.470142064105115f, 0.059715871789770f };
const value_t g7_7[2] = { 0.059715871789770f, 0.470142064105115f };//第三个分量为1-

const value_t gt4_1[4] = { 0.554370564668483f,0.148543145110506f,0.148543145110506f ,0.148543145110506f };
const value_t gt4_2[4] = { 0.148543145110506f,0.554370564668483f,0.148543145110506f ,0.148543145110506f };
const value_t gt4_3[4] = { 0.148543145110506f,0.148543145110506f,0.554370564668483f ,0.148543145110506f };
const value_t gt4_4[4] = { 0.148543145110506f,0.148543145110506f,0.148543145110506f ,0.554370564668483f };

const value_t wt4[4] = { 0.25f,0.25f,0.25f,0.25f };

void Gauss1Point(const VectorR3& v1, const VectorR3& v2, const VectorR3& v3, VectorR3* p);

void Gauss1Point(const VectorR3& v1, const VectorR3& v2, const VectorR3& v3, const VectorR3& v4, VectorR3* p);

void Gauss3Point(const VectorR3& v1, const VectorR3& v2, const VectorR3& v3, VectorR3* p);

void Gauss4Point(const VectorR3& v1, const VectorR3& v2, const VectorR3& v3, const VectorR3& v4, VectorR3* p);

void Gauss4Point(const VectorR3& v1, const VectorR3& v2, const VectorR3& v3, VectorR3* p);

void Gauss5Point(VectorR3& v1, VectorR3& v2, VectorR3& v3, VectorR3& v4, VectorR3* p);

void Gauss7Point(const VectorR3& v1, const VectorR3& v2, const VectorR3& v3, VectorR3* p);

void GaussLine3Point(const VectorR3& v1, const VectorR3& v2, VectorR3*p);

class Legendre {
public:
    value_t operator()(value_t x, int n)
    {
        value_t val;
        if (n == 0)
        {
            prev2 = 1.0f;
            val = 1.0f;
        }
        else if (n == 1)
        {
            prev1 = x;
            val = x;
        }
        else
        {
            val = ((2 * n - 1) * x * prev1 - (n - 1) * prev2) / n;
            prev2 = prev1;
            prev1 = val;
        }
        return val;
    }
private:
    value_t prev1, prev2;
};

class SphHankel {
public:
    Complex operator()(value_t x, int n)
    {
        Complex val;
        if (n == 0)
        {
            prev2 = Complex(sin(x), cos(x)) / x;
            val = prev2;
        }
        else if (n == 1)
        {
            prev1 = Complex(sin(x) / x - cos(x), cos(x) / x + sin(x)) / x;
            val = prev1;
        }
        else
        {
            val = ((2 * n - 1) / x) * prev1 - prev2;
            prev2 = prev1;
            prev1 = val;
        }
        return val;
    }
private:
    Complex prev1, prev2;
};

//  Gauss-Lengendre roots and weights 
//  Ref: rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature
class GaussLegendre {
    class Evaluation {
    public:
        explicit Evaluation(value_t x, int n);
        void evaluate(value_t x, int n);

        value_t v() const { return this->_v; }
        value_t d() const { return this->_d; }
        value_t x() const { return this->_x; }

    private:
        value_t _x;
        value_t _v;
        value_t _d;
    };
public:
    GaussLegendre();

    void set(int n);
    value_t root(size_t i) const { return _r[i]; }
    value_t weight(size_t i) const { return _w[i]; }

    void printRootsAndWeights() const;

private:
    std::vector<value_t> _r;
    std::vector<value_t> _w;
};

} // namespace math