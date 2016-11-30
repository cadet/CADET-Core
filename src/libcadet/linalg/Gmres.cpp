// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2016: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "cadet/cadetCompilerInfo.hpp"
#include "linalg/Gmres.hpp"

#include <sundials/sundials_spgmr.h>
#include "SundialsVector.hpp"

#include <type_traits>

namespace cadet
{

namespace linalg
{

// Wrapper function that calls the user function with the supplied user data
int gmresCallback(void* userData, N_Vector v, N_Vector z)
{
	Gmres* const g = static_cast<Gmres*>(userData);
	Gmres::MatrixVectorMultFun callback = g->matrixVectorMultiplier();
	return callback(g->userData(), NVEC_DATA(v), NVEC_DATA(z));
}

Gmres::Gmres() CADET_NOEXCEPT : _mem(nullptr), _ortho(Orthogonalization::ModifiedGramSchmidt), _maxRestarts(0), _matrixSize(0), _matVecMul(nullptr), _userData(nullptr)
{
}

Gmres::~Gmres() CADET_NOEXCEPT
{
	if (_mem)
		SpgmrFree(_mem);
}

void Gmres::initialize(unsigned int matrixSize, unsigned int maxKrylov)
{
	initialize(matrixSize, maxKrylov, Orthogonalization::ModifiedGramSchmidt, 10);
}

void Gmres::initialize(unsigned int matrixSize, unsigned int maxKrylov, Orthogonalization om, unsigned int maxRestarts)
{
	_matrixSize = matrixSize;
	if (maxKrylov == 0)
		maxKrylov = _matrixSize;

	// Create a template vector for the malloc routine of SPGMR
	N_Vector NV_tmpl = NVec_New(matrixSize);
	NVec_Const(0.0, NV_tmpl);

	// Size of allocated memory is either _maxKrylov or _cc.neq_bnd()
	_mem = SpgmrMalloc(maxKrylov, NV_tmpl);

	NVec_Destroy(NV_tmpl);
}

int Gmres::solve(double tolerance, double const* weight, double const* rhs, double* sol)
{
	// Create init-guess/solution vector by bending pointer
	N_Vector NV_sol = NVec_NewEmpty(_matrixSize);
	NVEC_DATA(NV_sol) = sol;

	// Create weight vector by bending pointer
	N_Vector NV_weight = NVec_NewEmpty(_matrixSize);
	NVEC_DATA(NV_weight) = const_cast<double*>(weight);

	// Create right hand side vector by pointer bending
	N_Vector NV_rhs = NVec_NewEmpty(_matrixSize);
	NVEC_DATA(NV_rhs) = const_cast<double*>(rhs);

//	double tolerance = _cc.sqrt_neq() * IDA_mem->ida_epsNewt * _schurSafety;

	const int gsType = static_cast<typename std::underlying_type<Orthogonalization>::type>(_ortho);
	int nIter = 0;
	int nPrecondSolve = 0;
	double res_norm = -1.0;

	const int flag = SpgmrSolve(_mem, this, NV_sol, NV_rhs,
			PREC_NONE, gsType, tolerance, _maxRestarts, NULL,
			NV_weight, NV_weight, &gmresCallback, NULL, 
			&res_norm, &nIter, &nPrecondSolve);

	// Free NVector memory space
	NVec_Destroy(NV_rhs);
	NVec_Destroy(NV_weight);
	NVec_Destroy(NV_sol);

	return flag;
}

const char* Gmres::getReturnFlagName(int flag) const CADET_NOEXCEPT
{
	switch (flag)
	{
	case  0: return "SPGMR_SUCCESS";            // Converged
	case  1: return "SPGMR_RES_REDUCED";        // Did not converge, but reduced
												// norm of residual

	case  2: return "SPGMR_CONV_FAIL";          // Failed to converge
	case  3: return "SPGMR_QRFACT_FAIL";        // QRfact found singular matrix
	case  4: return "SPGMR_PSOLVE_FAIL_REC";    // psolve failed recoverably
	case  5: return "SPGMR_ATIMES_FAIL_REC";    // atimes failed recoverably
	case  6: return "SPGMR_PSET_FAIL_REC";      // pset faild recoverably

	case -1: return "SPGMR_MEM_NULL";           // mem argument is NULL
	case -2: return "SPGMR_ATIMES_FAIL_UNREC";  // atimes returned failure flag
	case -3: return "SPGMR_PSOLVE_FAIL_UNREC";  // psolve failed unrecoverably
	case -4: return "SPGMR_GS_FAIL";            // Gram-Schmidt routine faiuled
	case -5: return "SPGMR_QRSOL_FAIL";         // QRsol found singular R
	case -6: return "SPGMR_PSET_FAIL_UNREC";    // pset failed unrecoverably
	default: return "NO_VALID_FLAG";
	}
}

}  // namespace linalg

}  // namespace cadet
