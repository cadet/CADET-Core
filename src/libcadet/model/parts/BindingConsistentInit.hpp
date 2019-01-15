// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2019: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file 
 * Defines the operator responsible for binding.
 */

#ifndef LIBCADET_BINDINGCONSISTENTINIT_HPP_
#define LIBCADET_BINDINGCONSISTENTINIT_HPP_

//#include "common/CompilerSpecific.hpp"
#include "linalg/BandMatrix.hpp"
#include "linalg/DenseMatrix.hpp"

#include <algorithm>

namespace cadet
{

namespace model
{

namespace parts
{

/**
 * @brief Handles consistent initialization of the solid phase of cells with binding
 */
class BindingConsistentInitializer
{
public:

	template <typename StateType, typename ResidualType, typename ParamType, typename BindingType>
	static inline void residual(const ParamType& t, double z, double r, unsigned int secIdx, const ParamType& timeFactor, 
		StateType const* y, double const* yDot, ResidualType* res, BindingType* const binding, double* const buffer)
	{
		binding->residual(t, z, r, secIdx, timeFactor, y, yDot, res, buffer);
	}

	template <typename StateType, typename ResidualType, typename ParamType, typename BindingType, typename MatrixRowIterator>
	static inline void residualWithJacobian(const ParamType& t, double z, double r, unsigned int secIdx, const ParamType& timeFactor, 
		StateType const* y, double const* yDot, ResidualType* res, BindingType* const binding, double* const buffer, MatrixRowIterator jac)
	{
		binding->residual(t, z, r, secIdx, timeFactor, y, yDot, res, buffer);
		binding->analyticJacobian(static_cast<double>(t), z, r, secIdx, reinterpret_cast<double const*>(y), jac, buffer);
	}

	/**
	 * @brief Sets up the linear equations required for consistent initialization of state vector time derivatives in cells with binding
	 * @details State vector time derivatives are consistently initialized by solving a linear system.
	 *          In case of dynamic equations, the matrix row corresponds to the Jacobian of the binding model
	 *          equations with respect to @f$ \dot{y} @f$. For algebraic equations, we need to take the derivatve
	 *          of the equation with respect to time, which is the Jacobian with respect to the state vector @f$ y @f$.
	 *          The right hand side is in both cases given by @f$ 0 @f$. An exception is presented by dynamic equations
	 *          that depend on time (e.g., through use of external functions), where the right hand side is 
	 *          @f[ -\frac{\partial \mathrm{res}(t, y, \dot{y})}{\partial t} @f]
	 *          
	 *          Note that the Jacobian row iterators are not advanced by this function.
	 * 
	 * @param [in] binding Binding model
	 * @param [in] timeFactor Used to compute parameter derivatives with respect to section length,
	 *             originates from time transformation and is premultiplied to time derivatives
	 * @param [in] jac Row iterator to first row of bound phase in the factorizable Jacobian matrix
	 * @param [in] origJacobian Row iterator to first row of bound phase in original Jacobian matrix
	 * @param [in,out] qShellDot Pointer to bound phase in state vector time derivatives
	 * @param [in] t Current time point
	 * @param [in] z Axial position in normalized coordinates (column inlet = 0, column outlet = 1)
	 * @param [in] r Radial position in normalized coordinates (outer shell = 1, inner center = 0)
	 * @param [in] secIdx Index of the current section
	 * @param [in,out] workSpace Memory work space
	 * @tparam BindingType Type of binding model
	 * @tparam MatrixIteratorType Type of original matrix row iterator 
	 * @tparam FactorizableMatrixIteratorType Type of factorizable matrix row iterator 
	 */
	// CADET_ALWAYS_INLINE
	template <typename BindingType, typename MatrixIteratorType, typename FactorizableMatrixIteratorType>
	static inline void consistentInitialTimeDerivative(BindingType* const binding, double timeFactor, FactorizableMatrixIteratorType jac, MatrixIteratorType origJacobian, double* qShellDot,
		double t, double z, double r, unsigned int secIdx, void* workSpace)
	{
		// Overwrite rows corresponding to algebraic equations with the Jacobian and set right hand side to 0

		// Get start and length of algebraic block
		unsigned int algStart = 0;
		unsigned int algLen = 0;
		binding->getAlgebraicBlock(algStart, algLen);

		// Get row iterators to algebraic block
		jac += algStart;
		origJacobian += algStart;

		// Pointer to right hand side of algebraic block
		double* const qShellDotAlg = qShellDot + algStart;

		// Copy matrix rows
		for (unsigned int algRow = 0; algRow < algLen; ++algRow, ++jac, ++origJacobian)
			jac.copyRowFrom(origJacobian);

		// Right hand side is -\frac{\partial res(t, y, \dot{y})}{\partial t}
		// If the residual is not explicitly depending on time, this expression is 0
		std::fill(qShellDotAlg, qShellDotAlg + algLen, 0.0);
		if (binding->dependsOnTime())
		{
			// TODO: Memory does alias since qShellDotAlg array is subarray of qShellDot array
			binding->timeDerivativeAlgebraicResidual(t, z, r, secIdx, qShellDot, qShellDotAlg, workSpace);
			for (unsigned int algRow = 0; algRow < algLen; ++algRow)
				qShellDotAlg[algRow] *= -1.0;
		}
	}

	/**
	 * @brief Sets up the linear equations required for consistent initialization of state vector time derivatives in cells with binding
	 * @details State vector time derivatives are consistently initialized by solving a linear system.
	 *          In case of dynamic equations, the matrix row corresponds to the Jacobian of the binding model
	 *          equations with respect to @f$ \dot{y} @f$. For algebraic equations, we need to take the derivatve
	 *          of the equation with respect to time, which is the Jacobian with respect to the state vector @f$ y @f$.
	 *          The right hand side is in both cases given by @f$ 0 @f$. An exception is presented by dynamic equations
	 *          that depend on time (e.g., through use of external functions), where the right hand side is 
	 *          @f[ -\frac{\partial \mathrm{res}(t, y, \dot{y})}{\partial t} @f]
	 *          
	 *          Note that the Jacobian row iterators are not advanced by this function.
	 * 
	 * @param [in] algStart Start index of the algebraic block in the bound phase
	 * @param [in] algLen Length of the algebraic block
	 * @param [in] binding Binding model
	 * @param [in] timeFactor Used to compute parameter derivatives with respect to section length,
	 *             originates from time transformation and is premultiplied to time derivatives
	 * @param [in] jacobianMatrix Factorizable matrix of size @p algLen x @p algLen
	 * @param [in] origJacobian Original Jacobian matrix with respect to state vector
	 * @param [in] localCpOffset Offset to the current particle shell's liquid phase from the beginning of the state vector
	 * @param [in] jacRowOffset Offset to the current particle shell's liquid phase from the beginning of the particle block
	 * @param [in] strideParLiquid Number of liquid phase DOFs in a particle shell
	 * @param [in] strideBound Number of bound phase DOFs in a particle shell
	 * @param [in,out] sensY Sensitivity state vector, will be modified with consistent initial values for the algebraic variables in the current particle shell
	 * @param [in] sensYdot Negative residual (right hand side for algebraic equation system)
	 * @param [in,out] workSpace Memory work space for preconditioning (@c nullptr to disable preconditioning)
	 * @tparam MatrixType Type of original matrix
	 * @tparam FactorizableMatrixType Type of factorizable matrix
	 */
	// CADET_ALWAYS_INLINE
	template <typename MatrixType, typename FactorizableMatrixType>
	static inline void consistentInitialSensitivityState(const unsigned int algStart, const unsigned int algLen, FactorizableMatrixType& jacobianMatrix, 
		MatrixType& origJacobian, const unsigned int localCpOffset, const unsigned int jacRowOffset, const int strideParLiquid, const unsigned int strideBound,
		double* const sensY, double* const sensYdot, double* const workSpace)
	{
		const int localOffset = localCpOffset + strideParLiquid;

		// Get pointer to q variables in a shell of particle pblk
		double* const qShell = sensY + localOffset;
		// Pointer to -dF / dp
		double* const dFdP = sensYdot + localOffset;
		// Pointer to c_p variables in this shell
		double* const cpShell = sensY + localCpOffset;

		// In general, the linear system looks like this
		// [c_p | q_diff | q_alg | q_diff ] * state + dF /dp = 0
		// We want to solve the q_alg block, which means we have to solve
		// [q_alg] * state = -[c_p | q_diff | 0 | q_diff ] * state - dF / dp
		// Note that we do not have to worry about fluxes since we are dealing
		// with bound states here.

		// Overwrite state with right hand side

		// Copy -dF / dp to state
		std::copy(dFdP + algStart, dFdP + algStart + algLen, qShell + algStart);

		// Subtract [c_p | q_diff] * state
		consistentInitialSensitivityStateLeftSubmatrixMultiply(origJacobian, cpShell, qShell, jacRowOffset, algStart, algLen, strideParLiquid);

		// Subtract [q_diff] * state (potential differential block behind q_alg block)
		if (algStart + algLen < strideBound)
			consistentInitialSensitivityStateRightSubmatrixMultiply(origJacobian, qShell, jacRowOffset, algStart, algLen, strideBound);

		// Copy main block to dense matrix
		consistentInitialSensitivityStateCopySubmatrix(jacobianMatrix, origJacobian, jacRowOffset, algStart, algLen);

		// Solve for algebraic variables, possibly using preconditioning
		if (workSpace)
			solve(jacobianMatrix, qShell + algStart, workSpace);
		else
			solve(jacobianMatrix, qShell + algStart);
	}

	/**
	 * @brief Sets up the linear equations required for consistent initialization of state vector time derivatives in cells with binding
	 * @details State vector time derivatives are consistently initialized by solving a linear system.
	 *          In case of dynamic equations, the matrix row corresponds to the Jacobian of the binding model
	 *          equations with respect to @f$ \dot{y} @f$. For algebraic equations, we need to take the derivatve
	 *          of the equation with respect to time, which is the Jacobian with respect to the state vector @f$ y @f$.
	 *          The right hand side is in both cases given by @f$ 0 @f$. An exception is presented by dynamic equations
	 *          that depend on time (e.g., through use of external functions), where the right hand side is 
	 *          @f[ -\frac{\partial \mathrm{res}(t, y, \dot{y})}{\partial t} @f]
	 *          
	 *          Note that the Jacobian row iterators are not advanced by this function.
	 * 
	 * @param [in] binding Binding model
	 * @param [in] jac Row iterator to first row of bound phase in the factorizable Jacobian matrix
	 * @param [in] origJacobian Row iterator to first row of bound phase in original Jacobian matrix
	 * @param [in,out] qShellDot Pointer to bound phase in state vector time derivatives
	 */
	// CADET_ALWAYS_INLINE
	template <typename BindingType, typename MatrixIteratorType, typename FactorizableMatrixIteratorType>
	static inline void consistentInitialSensitivityTimeDerivative(BindingType* const binding, FactorizableMatrixIteratorType jac, MatrixIteratorType origJacobian, double* qShellDot)
	{
		// Overwrite rows corresponding to algebraic equations with the Jacobian and set right hand side to 0

		// Get start and length of algebraic block
		unsigned int algStart = 0;
		unsigned int algLen = 0;
		binding->getAlgebraicBlock(algStart, algLen);

		// Get row iterators to algebraic block
		jac += algStart;
		origJacobian += algStart;

		// Pointer to right hand side of algebraic block
		double* const qShellDotAlg = qShellDot + algStart;

		// Copy rows and reset right hand side
		for (unsigned int algRow = 0; algRow < algLen; ++algRow, ++jac, ++origJacobian)
		{
			jac.copyRowFrom(origJacobian);

			// Right hand side is -\frac{\partial^2 res(t, y, \dot{y})}{\partial p \partial t}
			// If the residual is not explicitly depending on time, this expression is 0
			// @todo This is wrong if external functions are used. Take that into account!
			qShellDotAlg[algRow] = 0.0;
		}
	}

private:

	/**
	 * @brief Multiplies the submatrix in front of the algebraic block with a given vector
	 * @details Due to different interfaces, banded matrices and dense matrices have to be handled differently.
	 *          Dense matrices use absolute column indexing, whereas banded matrices use diagonal based indexing.
	 *          This is achieved by specializing this function template for each matrix type.
	 * 
	 * @param [in] origJacobian Original Jacobian
	 * @param [in] cpShell Pointer to first component in liquid phase of state vector
	 * @param [in,out] qShell Pointer to first bound state in state vector
	 * @param [in] jacRowOffset Row offset of the algebraic block in the original Jacobian
	 * @param [in] algStart Offset of the algebraic block among the variables
	 * @param [in] algLen Length of the algebraic block
	 * @param [in] strideBound Total number of bound states per cell
	 * @tparam MatrixType Type of original Jacobian matrix
	 */
	static inline void consistentInitialSensitivityStateLeftSubmatrixMultiply(const linalg::BandMatrix& origJacobian, double const* const cpShell, double* const qShell,
		const unsigned int jacRowOffset, const unsigned int algStart, const unsigned int algLen, const int strideParLiquid)
	{
		// Subtract [c_p | q_diff] * state
		// Note that band matrices use start diagonal index
		origJacobian.submatrixMultiplyVector(cpShell, jacRowOffset + algStart, -strideParLiquid - static_cast<int>(algStart), 
			algLen, static_cast<unsigned int>(strideParLiquid) + algStart, -1.0, 1.0, qShell + algStart);
	}

	static inline void consistentInitialSensitivityStateLeftSubmatrixMultiply(const linalg::detail::DenseMatrixBase& origJacobian, double const* const cpShell, double* const qShell,
		const unsigned int jacRowOffset, const unsigned int algStart, const unsigned int algLen, const int strideParLiquid)
	{
		// Subtract [c_p | q_diff] * state
		// Note that dense matrices use absolute column start index
		origJacobian.submatrixMultiplyVector(cpShell, jacRowOffset + algStart, static_cast<int>(jacRowOffset) - strideParLiquid, 
			algLen, static_cast<unsigned int>(strideParLiquid) + algStart, -1.0, 1.0, qShell + algStart);
	}

	/**
	 * @brief Multiplies the submatrix behind the algebraic block with a given vector
	 * @details Due to different interfaces, banded matrices and dense matrices have to be handled differently.
	 *          Dense matrices use absolute column indexing, whereas banded matrices use diagonal based indexing.
	 *          This is achieved by specializing this function template for each matrix type.
	 * 
	 * @param [in] origJacobian Original Jacobian
	 * @param [in,out] qShell Pointer to first bound state in state vector
	 * @param [in] jacRowOffset Row offset of the algebraic block in the original Jacobian
	 * @param [in] algStart Offset of the algebraic block among the variables
	 * @param [in] algLen Length of the algebraic block
	 * @param [in] strideBound Total number of bound states per cell
	 * @tparam MatrixType Type of original Jacobian matrix
	 */
	static inline void consistentInitialSensitivityStateRightSubmatrixMultiply(const linalg::BandMatrix& origJacobian, double* const qShell,
		const unsigned int jacRowOffset, const unsigned int algStart, const unsigned int algLen, const unsigned int strideBound)
	{
		// Subtract [q_diff] * state (potential differential block behind q_alg block)
		origJacobian.submatrixMultiplyVector(qShell + algStart + algLen, jacRowOffset + algStart, algStart + algLen, 
			algLen, strideBound - algStart - algLen, -1.0, 1.0, qShell + algStart);
	}

	static inline void consistentInitialSensitivityStateRightSubmatrixMultiply(const linalg::detail::DenseMatrixBase& origJacobian, double* const qShell,
		const unsigned int jacRowOffset, const unsigned int algStart, const unsigned int algLen, const unsigned int strideBound)
	{
		// Subtract [q_diff] * state (potential differential block behind q_alg block)
		origJacobian.submatrixMultiplyVector(qShell + algStart + algLen, jacRowOffset + algStart, static_cast<unsigned int>(jacRowOffset) + algStart + algLen, 
			algLen, strideBound - algStart - algLen, -1.0, 1.0, qShell + algStart);
	}

	/**
	 * @brief Copies the submatrix from the original Jacobian that refers to the algebraic block into another matrix
	 * @details Due to different interfaces, banded matrices and dense matrices have to be handled differently.
	 *          Dense matrices use absolute column indexing, whereas banded matrices use diagonal based indexing.
	 *          This is achieved by specializing this function template for each matrix type.
	 * 
	 * @param [out] jacobianMatrix Target matrix the algebraic block is copied into
	 * @param [in] origJacobian Original Jacobian block referring to algebraic variables
	 * @param [in] jacRowOffset Row offset of the algebraic block in the original Jacobian
	 * @param [in] algStart Offset of the algebraic block among the variables
	 * @param [in] algLen Length of the algebraic block
	 * @tparam MatrixType Type of the source matrix
	 */
	static inline void consistentInitialSensitivityStateCopySubmatrix(linalg::detail::DenseMatrixBase& jacobianMatrix, const linalg::BandMatrix& origJacobian, 
		const unsigned int jacRowOffset, const unsigned int algStart, const unsigned int algLen)
	{
		jacobianMatrix.copySubmatrixFromBanded(origJacobian, jacRowOffset + algStart, 0, algLen, algLen);
	}

	static inline void consistentInitialSensitivityStateCopySubmatrix(linalg::detail::DenseMatrixBase& jacobianMatrix, const linalg::detail::DenseMatrixBase& origJacobian, 
		const unsigned int jacRowOffset, const unsigned int algStart, const unsigned int algLen)
	{
		jacobianMatrix.copySubmatrix(origJacobian, jacRowOffset + algStart, jacRowOffset + algStart, algLen, algLen);
	}

	/**
	 * @brief Solves a linear system with diagonal preconditioning
	 * @param [in,out] mat On entry matrix to be solved, on exit factorized matrix
	 * @param [in,out] rhs On entry right hand side of equation system, on exit solution
	 * @param [in,out] workSpace Scratch space of the same size as @p rhs for storing diagonal preconditioner
	 * @tparam MatrixType Type of factorizable matrix
	 */
	template <typename MatrixType>
	static inline void solve(MatrixType& mat, double* const rhs, double* const workSpace)
	{
		// Precondition
		mat.rowScaleFactors(workSpace);
		mat.scaleRows(workSpace);

		// Solve
		mat.factorize();
		mat.solve(workSpace, rhs);
	}

	/**
	 * @brief Solves a linear system (without diagonal preconditioning)
	 * @param [in,out] mat On entry matrix to be solved, on exit factorized matrix
	 * @param [in,out] rhs On entry right hand side of equation system, on exit solution
	 * @tparam MatrixType Type of factorizable matrix
	 */
	template <typename MatrixType>
	static inline void solve(MatrixType& mat, double* const rhs)
	{
		mat.factorize();
		mat.solve(rhs);
	}

};

} // namespace operators
} // namespace model
} // namespace cadet

#endif  // LIBCADET_BINDINGCONSISTENTINIT_HPP_
