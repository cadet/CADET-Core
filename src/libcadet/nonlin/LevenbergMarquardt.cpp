// =============================================================================
//  CADET
//  
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "nonlin/LevenbergMarquardt.hpp"
#include "cadet/ParameterProvider.hpp"
#include "cadet/Exceptions.hpp"
#include "linalg/DenseMatrix.hpp"

#include <string>

namespace cadet
{

namespace nonlin
{

LevenbergMarquardtSolver::LevenbergMarquardtSolver() : _initDamping(1e-2), _maxIter(50) { }
LevenbergMarquardtSolver::~LevenbergMarquardtSolver() { }

bool LevenbergMarquardtSolver::configure(IParameterProvider& paramProvider)
{
	if (paramProvider.exists("INIT_DAMPING"))
	{
		_initDamping = paramProvider.getDouble("INIT_DAMPING");
		if (_initDamping < 0.0)
			throw InvalidParameterException("Field INIT_DAMPING has to be non-negative, but is " + std::to_string(_initDamping));
	}
	if (paramProvider.exists("MAX_ITERATIONS"))
	{
		const int maxIter = paramProvider.getInt("MAX_ITERATIONS");
		if (maxIter < 1)
			throw InvalidParameterException("Field MAX_ITERATIONS has to be at least 1, but is " + std::to_string(maxIter));

		_maxIter = static_cast<unsigned int>(maxIter);
	}
	return true;
}

bool LevenbergMarquardtSolver::solve(std::function<bool(double const* const, double* const)> residual, std::function<bool(double const* const, linalg::detail::DenseMatrixBase& jac)> jacobian,
		double tol, double* const point, double* const workingMemory, linalg::detail::DenseMatrixBase& jacMatrix, unsigned int size) const
{
	return levenbergMarquardt(residual, jacobian, _maxIter, tol, _initDamping, point, workingMemory, jacMatrix, size);
}


} // namespace nonlin

} // namespace cadet
