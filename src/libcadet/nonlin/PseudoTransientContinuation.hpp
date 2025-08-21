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

/**
 * @file 
 * Provides pseudo-transient continuation method for solving nonlinear equation systems
 */

#ifndef LIBCADET_PSEUDOTRANSIENTCONTINUATION_HPP_
#define LIBCADET_PSEUDOTRANSIENTCONTINUATION_HPP_

#include "common/CompilerSpecific.hpp"
#include "nonlin/Solver.hpp"

#include <cmath>
#include <functional>
#include <algorithm>
#include <vector>

#include "linalg/Norms.hpp"
#include "linalg/DenseMatrix.hpp"

namespace cadet
{

namespace nonlin
{

	/**
	 * @brief Iterate output policy that does nothing
	 */
	struct VoidPTCIterateOutputPolicy
	{
		inline static void iteration(unsigned int idxIter, double const* const x, double const* const fx, double residualNorm, double stepSize, unsigned int size) { }
	};

	/**
	 * @brief Solves nonlinear equations using a pseudo-transient continuation method
	 * @details The nonlinear equation system @f$ f(x) = 0 @f$ is solved using pseudo-transient
     *          continuation (PTC), which introduce pseudo time and computes the steady
     *          state of the ODE @f$ dx / dt = f(x) @f$.
     *
     *          After resolving the transient regime at the beginning of the time
     *          integration, PTC quickly increases the time step size to arrive at the
     *          equilibrium. It can be interpreted as a linear-implicit Euler scheme for
     *          the ODE, which tends to a Newton method as the time step size increases.
     *
     *          The step size is adaptively determined by switched evolution relaxation
     *          (SER).
     *
     *          See \cite Deuflhard2011.
	 *          
	 *          This algorithm is designed to solve the nonlinear equation system
	 *          @f[\begin{align} f(x) = 0, \qquad f: \mathbb{R}^n \to \mathbb{R}^n. \end{align}@f]
	 *          A solution is indicated by the error test
	 *          @f[\begin{align} \left\lVert f(x) \right\rVert_{\ell^2} \leq \text{tol} \end{align}@f]
	 *          with a given tolerance (@p resTol).
	 *          
	 *          The residuals @f$ f(x) @f$ are evaluated using the given function @p residual.
	 *          The Jacobian @f$ \frac{\mathrm{d}f}{\mathrm{d}x}(x) @f$ are evaluated using the given
	 *          function @p jacobian.
	 *          
	 *          The initial step size of the linear-implicit Euler scheme is provided in @p tau.
	 * @param [in] residual Function providing the residual @f$ f(x) @f$ at position @f$ x @f$ of the nonlinear equation 
	 *             system @f$ f(x) = 0 @f$ to be solved.
	 *             The signature of the function is
	 *             `bool residual(double const* const x, double* const r)`
	 *             where the return value communicates whether the evaluation has been successful.
	 *             On exit, the residual is returned in-place in @f$ r @f$.
	 *             If an error is indicated by returning @c false, the algorithm aborts with failure.
	 * @param [in] jacobian Function calculating the Jacobian matrix @f$ J_f(x) @f$ at position @f$ x @f$. 
	 *             The signature of the function is
	 *             `bool jacobian(double const* const x, linalg::detail::DenseMatrixBase& J)`
	 *             where the return value communicates whether the evaluation has been successful.
	 *             On exit, the Jacobian is returned in-place in @f$ J @f$.
	 *             If an error is indicated by returning @c false, the algorithm aborts with failure.
	 * @param [in] maxIter Maximum number of iterations
	 * @param [in] resTol Termination criterion on the residual @f$\ell^2@f$-norm
	 * @param [in] maxNonMonotone Maximum number of non-monotone residual updates before the algorithm aborts
	 * @param [in] tau Initial time step for the Euler scheme
	 * @param [in] scale Scaling factors, optional (use @c nullptr to disable)
	 * @param [in] variant Use a different variant for updating the position, which is possibly numerically more stable
	 * @param [in,out] point On entry initial guess, on exit solution or last iterate
	 * @param [in] jacMat Memory for Jacobian matrix
	 * @param [in] workingMemory Additional memory of size @f$ 4n @f$ required for performing the iterations, where @f$ n @f$ is the problem @p size
	 * @param [in] size Size of the problem (i.e., number of equations, length of residual, columns of Jacobian etc.)
	 * @tparam IterateOutputPolicy Policy that handles output of intermediate values (useful for debugging), see VoidPTCIterateOutputPolicy
	 * @return @c true if a solution meeting the residual tolerance was found, @c false otherwise
	 */
	template <typename IterateOutputPolicy = VoidPTCIterateOutputPolicy>
	bool pseudoTransientContinuation(std::function<bool(double const* const, double* const)> residual, std::function<bool(double const* const, linalg::detail::DenseMatrixBase&)> jacobian,
		int maxIter, double resTol, int maxNonMonotone, double tau, double const* const scale, bool variant, double* const point, linalg::detail::DenseMatrixBase& jacMat, double* const workingMemory, int size)
	{
		// Split working memory into parts
		double* const residualMem = workingMemory;
		double* const dx = workingMemory + size;
		double* const qrMem = workingMemory + 2 * size;

		residual(point, residualMem);
		double normFk = 0.0;
		if (scale)
			normFk = cadet::linalg::l2normWeighted(residualMem, scale, size);
		else
			normFk = cadet::linalg::l2Norm(residualMem, size);

		int k = 0;
		int nNonMonotoneResidual = 0;
		double inv_tau = 1.0 / tau;

		IterateOutputPolicy::iteration(0, point, residualMem, normFk, tau, size);

		while (k <= maxIter)
		{
			if (normFk <= resTol)
				return true;
			
			jacobian(point, jacMat);
			if (scale)
				jacMat.scaleRows(scale);

			for (int i = 0; i < size; ++i)
			{
				if (scale)
					dx[i] = residualMem[i] / scale[i];
				else
					dx[i] = residualMem[i];
			}

			if (variant)
			{
				// (1/tau * I - F'(x)) * dx = F(x), update x_new = x + dx
				// -(1/tau * I - F'(x)) * dx = -F(x), update x_new = x + dx
				// (F'(x) - 1/tau * I) * dx = -F(x), update x_new = x + dx
				// (F'(x) - 1/tau * I) * (-dx) = F(x), update x_new = x - dx

				if (scale)
				{
					for (int i = 0; i < size; ++i)
						jacMat.native(i, i) -= inv_tau / scale[i];
				}
				else
				{
					for (int i = 0; i < size; ++i)
						jacMat.native(i, i) -= inv_tau;
				}

				jacMat.robustFactorize(qrMem);
				jacMat.robustSolve(dx, qrMem);

				for (int i = 0; i < size; ++i)
					point[i] = point[i] - dx[i];
			}
			else
			{
				// (I - tau * F'(x)) * dx = F(x), update x_new = x + tau * dx
				for (int i = 0; i < size; ++i)
				{
					for (int j = 0; j < size; ++j)
						jacMat.native(i, j) *= -tau;

					if (scale)
						jacMat.native(i, i) += 1.0 / scale[i];
					else
						jacMat.native(i, i) += 1.0;
				}

				jacMat.robustFactorize(qrMem);
				jacMat.robustSolve(dx, qrMem);

				for (int i = 0; i < size; ++i)
					point[i] = point[i] + tau * dx[i];
			}

			residual(point, residualMem);

			double normFkp1 = 0.0;
			if (scale)
				normFkp1 = cadet::linalg::l2normWeighted(residualMem, scale, size);
			else
				normFkp1 = cadet::linalg::l2Norm(residualMem, size);

			// Guard against division by zero
			if (normFkp1 == 0.0)
				return true;

			// Update step size
			tau = tau * normFk / normFkp1;
			inv_tau = 1.0 / tau;

			// Check monotonicity of residual norm
			if (normFkp1 < normFk)
				nNonMonotoneResidual = 0;
			else
			{
				++nNonMonotoneResidual;
				if (nNonMonotoneResidual >= maxNonMonotone)
					return false;
			}

			normFk = normFkp1;
			++k;

			IterateOutputPolicy::iteration(k, point, residualMem, normFk, tau, size);
		}

		return true;
	}

	/**
	 * @brief Uses a pseudo-transient continuation method for solving nonlinear equations
	 * @details Wraps pseudoTransientContinuation() function.
	 */
	class PseudoTransientContinuationSolver : public Solver
	{
	public:
		PseudoTransientContinuationSolver();
		virtual ~PseudoTransientContinuationSolver();

		static const char* identifier() { return "PTC"; }
		virtual const char* name() const { return PseudoTransientContinuationSolver::identifier(); }
		virtual bool configure(IParameterProvider& paramProvider);

		virtual unsigned int workspaceSize(unsigned int problemSize) const
		{
			return 4 * problemSize;
		}
		
		virtual unsigned int numTuningParameters() const { return 5; }

		virtual bool solve(std::function<bool(double const* const, double* const)> residual, std::function<bool(double const* const, linalg::detail::DenseMatrixBase& jac)> jacobian,
			double tol, double* const point, double* const workingMemory, linalg::detail::DenseMatrixBase& jacMatrix, unsigned int size) const;
	
	protected:
		double _tau; //!< Initial step size
		double* _scale; //!< Scaling factors
		int _maxIter; //!< Maximum number of iterations
		int _numNonMonotone; //!< Maximum number of non-monotone error trajectory
		bool _variant; //!< Different numeric calculation, but mathematically equivalent
	};

} // namespace nonlin

} // namespace cadet

#endif  // LIBCADET_PSEUDOTRANSIENTCONTINUATION_HPP_
