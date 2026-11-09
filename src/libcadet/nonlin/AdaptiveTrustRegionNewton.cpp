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

#include "nonlin/AdaptiveTrustRegionNewton.hpp"
#include "cadet/ParameterProvider.hpp"
#include "cadet/Exceptions.hpp"
#include "linalg/DenseMatrix.hpp"

#include <string>

namespace cadet
{

namespace nonlin
{

AdaptiveTrustRegionNewtonSolver::AdaptiveTrustRegionNewtonSolver() : _initDamping(1e-2), _minDamping(1e-4), _maxIter(50) { }
AdaptiveTrustRegionNewtonSolver::~AdaptiveTrustRegionNewtonSolver() { }

bool AdaptiveTrustRegionNewtonSolver::configure(IParameterProvider& paramProvider)
{
	if (paramProvider.exists("INIT_DAMPING"))
	{
		_initDamping = paramProvider.getDouble("INIT_DAMPING");
		if (_initDamping < 0.0)
			throw InvalidParameterException("Field INIT_DAMPING has to be non-negative, but is " + std::to_string(_initDamping));
	}
	if (paramProvider.exists("MIN_DAMPING"))
	{
		_minDamping = paramProvider.getDouble("MIN_DAMPING");
		if (_minDamping < 0.0)
			throw InvalidParameterException("Field MIN_DAMPING has to be non-negative, but is " + std::to_string(_minDamping));
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

bool AdaptiveTrustRegionNewtonSolver::solve(std::function<bool(double const* const, double* const)> residual, std::function<bool(double const* const, linalg::detail::DenseMatrixBase& jac)> jacobian,
		double tol, double* const point, double* const workingMemory, linalg::detail::DenseMatrixBase& jacMatrix, unsigned int size) const
{
	double* const scaleFactors = workingMemory + 4 * size;
	return adaptiveTrustRegionNewtonMethod(residual, [&](double const* const x, double* const y) -> bool {
			if (!jacobian(x, jacMatrix))
				return false;

			jacMatrix.rowScaleFactors(scaleFactors);
			jacMatrix.scaleRows(scaleFactors);

			return jacMatrix.factorize() && jacMatrix.solve(scaleFactors, y);
		},
		_maxIter, tol, _initDamping, _minDamping, point, workingMemory, size);
}


RobustAdaptiveTrustRegionNewtonSolver::RobustAdaptiveTrustRegionNewtonSolver() : _initDamping(1e-2), _minDamping(1e-4), _maxIter(50) { }
RobustAdaptiveTrustRegionNewtonSolver::~RobustAdaptiveTrustRegionNewtonSolver() { }

bool RobustAdaptiveTrustRegionNewtonSolver::configure(IParameterProvider& paramProvider)
{
	if (paramProvider.exists("INIT_DAMPING"))
	{
		_initDamping = paramProvider.getDouble("INIT_DAMPING");
		if (_initDamping < 0.0)
			throw InvalidParameterException("Field INIT_DAMPING has to be non-negative, but is " + std::to_string(_initDamping));
	}
	if (paramProvider.exists("MIN_DAMPING"))
	{
		_minDamping = paramProvider.getDouble("MIN_DAMPING");
		if (_minDamping < 0.0)
			throw InvalidParameterException("Field MIN_DAMPING has to be non-negative, but is " + std::to_string(_minDamping));
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

bool RobustAdaptiveTrustRegionNewtonSolver::solve(std::function<bool(double const* const, double* const)> residual, std::function<bool(double const* const, linalg::detail::DenseMatrixBase& jac)> jacobian,
		double tol, double* const point, double* const workingMemory, linalg::detail::DenseMatrixBase& jacMatrix, unsigned int size) const
{
	double* const scaleFactors = workingMemory + 4 * size;
	return robustAdaptiveTrustRegionNewtonMethod(residual, [&](double const* const x, double* const y) -> bool {
			if (!jacobian(x, jacMatrix))
				return false;

			jacMatrix.rowScaleFactors(scaleFactors);
			jacMatrix.scaleRows(scaleFactors);

			return jacMatrix.factorize() && jacMatrix.solve(scaleFactors, y);
		},
		[&](double* const y) -> bool {
			return jacMatrix.solve(y);
		},
		_maxIter, tol, _initDamping, _minDamping, point, workingMemory, size);
}


} // namespace nonlin

} // namespace cadet
