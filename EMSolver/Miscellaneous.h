//**********************************************************
// Author: Zhang Qiang
// License: MIT
//**********************************************************

#pragma once
#include "Custom.h"
#include "Mathdef.h"
#include "Quadrature.h"
#include "VectorR3.h"
#include "VectorC3.h"

using component::VectorR3;

namespace mom {

namespace integral {

Complex cfieZKernel(value_t alpha, value_t k, const VectorR3 * vf3, const VectorR3 * vs3, const VectorR3 & vf,
                    const VectorR3 & vs, const VectorR3 & nf, value_t fsign, value_t ssign);

Complex cfieZSingular(value_t alpha, value_t k, const VectorR3 * vf3, const VectorR3 * vs3, const VectorR3 & vf,
                      const VectorR3 & vs, value_t fsign, value_t ssign);

Complex cfieZKernel(value_t k, const VectorR3 * vf3, const VectorR3 * vs3, const VectorR3 & vf,
	const VectorR3 & vs, const VectorR3 & nf, value_t fsign, value_t ssign,value_t msign,value_t dsign);

Complex cfieZSingular(value_t k, const VectorR3 * vf3, const VectorR3 * vs3, const VectorR3 & vf,
	const VectorR3 & vs, value_t fsign, value_t ssign, value_t msign, value_t dsign);

Complex cfieVKernel(value_t alpha, value_t k, const VectorR3& vec_k, const VectorR3& vec_ei, const VectorR3& vec_hi,
                    const VectorR3 * vp3, const VectorR3 * vm3, const VectorR3 & np, const VectorR3 & nm,
                    const VectorR3 & vp, const VectorR3 & vm);
} // namespace integral

} // namespace mom