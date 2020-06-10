//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "Custom.h"

namespace iml {

bool    BICGSTAB(const Qcx_mat& _A, const Qcx_vec& _b, Qcx_vec& _x, value_t _eps, int _max_iter, int& _iterNum);
bool    BICGSTAB(const CMatrix& _A, const CVector& _b, CVector& _x, value_t _eps, int _max_iter, int& _iterNum);
bool	BICGSTAB(CVector& x, MulVct mul, const CVector & b, int iterMax, value_t iterEps, int &iterNum);
bool	GMRES(const CMatrix& _A, const CVector& _b, CVector& _x, value_t _eps, int _max_iter,int _out_iter);

bool    CG(const Qcx_mat& _A, const Qcx_vec& _b, Qcx_vec& _x, value_t _eps, int _max_iter, int& _iterNum);
bool    CG(const CMatrix& _A, const CVector& _b, CVector& _x, value_t _eps, int _max_iter, int& _iterNum);
bool    BICG(const Qcx_mat& _A, const Qcx_vec& _b, Qcx_vec& _x, value_t _eps, int _max_iter, int& _iterNum);
bool    CGNE(const Qcx_mat& _A, const Qcx_vec& _b, Qcx_vec& _x, value_t _eps, int _max_iter, int& _iterNum);
bool    GCR(const Qcx_mat& A, const Qcx_vec& b, Qcx_vec& x, value_t eps, int max_iter, int& iter_num);

} // namespace iml

namespace test {

void GCR();

} // namespace test