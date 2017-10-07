// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2017: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file 
 * Defines simulation tests for column unit operations.
 */

#ifndef CADETTEST_COLUMNSIMTEST_HPP_
#define CADETTEST_COLUMNSIMTEST_HPP_

#include <limits>

namespace cadet
{

class JsonParameterProvider;

namespace test
{

namespace column
{

	/**
	 * @brief Sets the number of axial cells in a configuration of a column-like unit operation
	 * @details Overwrites the NCOL field in the discretization group of the given ParameterProvider.
	 * @param [in,out] jpp ParameterProvider to change the number of axial cells in
	 * @param [in] nCol Number of axial cells
	 */
	void setNumAxialCells(cadet::JsonParameterProvider& jpp, unsigned int nCol);

	/**
	 * @brief Sets the WENO order in a configuration of a column-like unit operation
	 * @details Overwrites the WENO_ORDER field in the weno group of the given ParameterProvider.
	 * @param [in,out] jpp ParameterProvider to change the WENO order in
	 * @param [in] order Target order
	 */
	void setWenoOrder(cadet::JsonParameterProvider& jpp, int order);

	/**
	 * @brief Reverses the flow of a column-like unit operation
	 * @param [in,out] jpp ParameterProvider to change the flow direction in
	 */
	void reverseFlow(cadet::JsonParameterProvider& jpp);

	/**
	 * @brief Sets the binding mode
	 * @details Overwrites the IS_KINETIC field in the adsorption group of the given ParameterProvider.
	 * @param [in,out] jpp ParameterProvider to change the binding mode in
	 * @param [in] isKinetic Determines whether kinetic or quasi-stationary binding is used
	 */
	void setBindingMode(cadet::JsonParameterProvider& jpp, bool isKinetic);

	/**
	 * @brief Runs a simulation test comparing against (semi-)analytic single component pulse injection reference data
	 * @details Linear binding model is used in the column-like unit operation.
	 * @param [in] uoType Unit operation type
	 * @param [in] refFileRelPath Path to the reference data file from the directory of this file
	 * @param [in] forwardFlow Determines whether the unit operates in forward flow (@c true) or backwards flow (@c false)
	 * @param [in] dynamicBinding Determines whether dynamic binding (@c true) or rapid equilibrium (@c false) is used
	 * @param [in] nCol Number of axial cells
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testAnalyticBenchmark(const char* uoType, const char* refFileRelPath, bool forwardFlow, bool dynamicBinding, unsigned int nCol, double absTol, double relTol);

	/**
	 * @brief Runs a simulation test comparing against (semi-)analytic single component pulse injection reference data
	 * @details The component is assumed to be non-binding.
	 * @param [in] uoType Unit operation type
	 * @param [in] refFileRelPath Path to the reference data file from the directory of this file
	 * @param [in] forwardFlow Determines whether the unit operates in forward flow (@c true) or backwards flow (@c false)
	 * @param [in] nCol Number of axial cells
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testAnalyticNonBindingBenchmark(const char* uoType, const char* refFileRelPath, bool forwardFlow, unsigned int nCol, double absTol, double relTol);

	/**
	 * @brief Runs a simulation test comparing forward and backwards flow in the load-wash-elution example
	 * @param [in] uoType Unit operation type
	 * @param [in] wenoOrder Order of the WENO method
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testWenoForwardBackward(const char* uoType, int wenoOrder, double absTol, double relTol);

	/**
	 * @brief Checks the full Jacobian against AD and FD pattern switching from forward to backward flow and back
	 * @details Checks the analytic Jacobian against the AD Jacobian and checks both against the FD pattern.
	 *          Checks both forward and backward flow mode as well as switching between them.
	 * 
	 * @param [in] uoType Unit operation type
	 * @param [in] wenoOrder Order of the WENO method
	 */
	void testJacobianWenoForwardBackward(const std::string& uoType, int wenoOrder);

	/**
	 * @brief Checks the (analytic) time derivative Jacobian against FD
	 * @details Uses centered finite differences.
	 * @param [in] uoType Unit operation type
	 * @param [in] h Step size of centered finite differences
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testTimeDerivativeJacobianFD(const std::string& uoType, double h = 1e-6, double absTol = 0.0, double relTol = std::numeric_limits<float>::epsilon() * 100.0);

	/**
	 * @brief Checks the forward sensitivity residual using analytic Jacobians
	 * @details Uses centered finite differences.
	 * @param [in] uoType Unit operation type
	 * @param [in] h Step size of centered finite differences
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testFwdSensJacobians(const std::string& uoType, double h = 1e-6, double absTol = 0.0, double relTol = std::numeric_limits<float>::epsilon() * 100.0);

} // namespace column
} // namespace test
} // namespace cadet

#endif
