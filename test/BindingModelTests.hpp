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
 * Defines tests for binding models.
 */

#ifndef CADETTEST_BINDINGMODELTEST_HPP_
#define CADETTEST_BINDINGMODELTEST_HPP_

#include <limits>

namespace cadet
{

namespace test
{

namespace binding
{

	/**
	 * @brief Checks the analytic Jacobian of the binding model against AD
	 * @param [in] modelName Name of the binding model
	 * @param [in] nComp Number of components
	 * @param [in] nBound Array with number of bound states for each component
	 * @param [in] isKinetic Determines whether kinetic or quasi-stationary binding mode is applied
	 * @param [in] config JSON string with binding model parameters
	 * @param [in] point Liquid phase and solid phase values to check Jacobian at
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testJacobianAD(const char* modelName, unsigned int nComp, unsigned int const* nBound, bool isKinetic, const char* config, double const* point, double absTol = 0.0, double relTol = std::numeric_limits<float>::epsilon() * 100.0);

	/**
	 * @brief Checks the time derivative Jacobian of the binding model against finite differences
	 * @details Uses centered finite differences. The analytic Jacobian is provided by jacobianAddDiscretized().
	 * @param [in] modelName Name of the binding model
	 * @param [in] nComp Number of components
	 * @param [in] nBound Array with number of bound states for each component
	 * @param [in] isKinetic Determines whether kinetic or quasi-stationary binding mode is applied
	 * @param [in] config JSON string with binding model parameters
	 * @param [in] h Step size of centered finite differences
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testTimeDerivativeJacobianFD(const char* modelName, unsigned int nComp, unsigned int const* nBound, bool isKinetic, const char* config, double h = 1e-6, double absTol = 0.0, double relTol = std::numeric_limits<float>::epsilon() * 100.0);

	/**
	 * @brief Checks the multiplyWithDerivativeJacobian() function against jacobianAddDiscretized()
	 * @param [in] modelName Name of the binding model
	 * @param [in] nComp Number of components
	 * @param [in] nBound Array with number of bound states for each component
	 * @param [in] isKinetic Determines whether kinetic or quasi-stationary binding mode is applied
	 * @param [in] config JSON string with binding model parameters
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testTimeDerivativeJacobianMultiplyFunction(const char* modelName, unsigned int nComp, unsigned int const* nBound, bool isKinetic, const char* config, double absTol = 0.0, double relTol = std::numeric_limits<float>::epsilon() * 100.0);

	/**
	 * @brief Checks consistency of getAlgebraicBlock() and hasAlgebraicEquations() against structure of time derivative Jacobian
	 * @details Algebraic equations are detected in the time derivative Jacobian by checking for all-zero rows.
	 *          The Jacobian is obtained from jacobianAddDiscretized().
	 * @param [in] modelName Name of the binding model
	 * @param [in] nComp Number of components
	 * @param [in] nBound Array with number of bound states for each component
	 * @param [in] isKinetic Determines whether kinetic or quasi-stationary binding mode is applied
	 * @param [in] config JSON string with binding model parameters
	 */
	void testConsistencyOfAlgebraicEquationFunctions(const char* modelName, unsigned int nComp, unsigned int const* nBound, bool isKinetic, const char* config);

	/**
	 * @brief Checks consistent initialization of algebraic equations
	 * @details Calls consistentInitialState() and checks if the residual of the algebraic equations is (almost) 0.
	 * @param [in] modelName Name of the binding model
	 * @param [in] nComp Number of components
	 * @param [in] nBound Array with number of bound states for each component
	 * @param [in] isKinetic Determines whether kinetic or quasi-stationary binding mode is applied
	 * @param [in] config JSON string with binding model parameters
	 * @param [in] useAD Determines whether the Jacobian is computed via AD or analytically
	 * @param [in] point Liquid phase and initial solid phase values to start the process from
	 * @param [in] consTol Error tolerance for consistent initialization solver
	 * @param [in] absTol Error tolerance for checking whether algebraic residual is 0
	 */
	 void testConsistentInitialization(const char* modelName, unsigned int nComp, unsigned int const* nBound, bool isKinetic, const char* config, bool useAD, double const* point, double consTol, double absTol);

	/**
	 * @brief Checks residual and analytic Jacobian of normal model variant against externally dependent ones
	 * @param [in] modelName Name of the binding model
	 * @param [in] modelNameExt Name of the externally dependent binding model variant
	 * @param [in] nComp Number of components
	 * @param [in] nBound Array with number of bound states for each component
	 * @param [in] isKinetic Determines whether kinetic or quasi-stationary binding mode is applied
	 * @param [in] config JSON string with binding model parameters for both variants
	 * @param [in] point Liquid phase and solid phase values to evaluate residual at
	 */
	void testNormalExternalConsistency(const char* modelName, const char* modelNameExt, unsigned int nComp, unsigned int const* nBound, bool isKinetic, const char* config, double const* point);

	/**
	 * @brief Checks whether Jacobian columns of non-binding liquid phase components are all zero
	 * @param [in] modelName Name of the binding model
	 * @param [in] nComp Number of components
	 * @param [in] nBound Array with number of bound states for each component
	 * @param [in] isKinetic Determines whether kinetic or quasi-stationary binding mode is applied
	 * @param [in] config JSON string with binding model parameters
	 * @param [in] useAD Determines whether the Jacobian is computed via AD or analytically
	 * @param [in] point Liquid phase and solid phase values to evaluate Jacobian at
	 */
	void testNonBindingConsistency(const char* modelName, unsigned int nComp, unsigned int const* nBound, bool isKinetic, const char* config, bool useAD, double const* point);

	/**
	 * @brief Checks whether Jacobian and residual of variants with all-binding and some non-binding components match
	 * @param [in] modelName Name of the binding model
	 * @param [in] nCompBnd Number of components in all binding variant
	 * @param [in] nCompNonBnd Number of components in non-binding variant
	 * @param [in] nBound Array with number of bound states for each component in all binding variant
	 * @param [in] nBoundNonBnd Array with number of bound states for each component in non-binding variant
	 * @param [in] isKinetic Determines whether kinetic or quasi-stationary binding mode is applied
	 * @param [in] configBnd JSON string with binding model parameters for all binding variant
	 * @param [in] configNonBnd JSON string with binding model parameters for non-binding variant
	 * @param [in] useAD Determines whether the Jacobian is computed via AD or analytically
	 * @param [in] pointBnd Liquid phase and solid phase values to evaluate Jacobian at in all binding variant
	 * @param [in] pointNonBnd Liquid phase and solid phase values to evaluate Jacobian at in non-binding variant
	 */
	void testNonbindingBindingConsistency(const char* modelName, unsigned int nCompBnd, unsigned int nCompNonBnd, unsigned int const* nBound, unsigned int const* nBoundNonBnd, bool isKinetic, const char* configBnd, const char* configNonBnd, bool useAD, double const* pointBnd, double const* pointNonBnd);

} // namespace binding
} // namespace test
} // namespace cadet

#endif  // CADETTEST_BINDINGMODELTEST_HPP_
