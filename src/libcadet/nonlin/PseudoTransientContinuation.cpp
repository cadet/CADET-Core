// =============================================================================
//  CADET
//  
//  Copyright © The CADET Authors
//            Please see the CONTRIBUTORS.md file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "nonlin/PseudoTransientContinuation.hpp"
#include "cadet/ParameterProvider.hpp"
#include "linalg/DenseMatrix.hpp"

namespace cadet
{

namespace nonlin
{

PseudoTransientContinuationSolver::PseudoTransientContinuationSolver() : _tau(20.0), _scale(nullptr), _maxIter(100), _numNonMonotone(5), _variant(false) { }
PseudoTransientContinuationSolver::~PseudoTransientContinuationSolver() { }

bool PseudoTransientContinuationSolver::configure(IParameterProvider& paramProvider)
{
	if (paramProvider.exists("INIT_PSEUDOTIMESTEP"))
		_tau = paramProvider.getDouble("INIT_PSEUDOTIMESTEP");
	if (paramProvider.exists("MAX_NONMONOTONE_ITERATIONS"))
		_numNonMonotone = paramProvider.getInt("MAX_NONMONOTONE_ITERATIONS");
	if (paramProvider.exists("MAX_ITERATIONS"))
		_maxIter = paramProvider.getInt("MAX_ITERATIONS");
	if (paramProvider.exists("NUMERIC_VARIANT"))
		_variant = paramProvider.getBool("NUMERIC_VARIANT");
	return true;
}

bool PseudoTransientContinuationSolver::solve(std::function<bool(double const* const, double* const)> residual, std::function<bool(double const* const, linalg::detail::DenseMatrixBase& jac)> jacobian,
		double tol, double* const point, double* const workingMemory, linalg::detail::DenseMatrixBase& jacMatrix, unsigned int size) const
{
	return pseudoTransientContinuation(residual, jacobian, _maxIter, tol, _numNonMonotone,
		_tau, _scale, _variant, point, jacMatrix, workingMemory, size);
}

} // namespace nonlin

} // namespace cadet
