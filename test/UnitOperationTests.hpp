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
 * Defines tests for general unit operations.
 */

#ifndef CADETTEST_UNITOPTESTS_HPP_
#define CADETTEST_UNITOPTESTS_HPP_

namespace cadet
{

class IUnitOperation;

namespace test
{

namespace unitoperation
{

	/**
	 * @brief Checks the consistent initialization of a generic unit operation
	 * @param [in] unit Configured unit operation
	 * @param [in] adEnabled Determines whether AD is used for calculating the Jacobian
	 * @param [in] y Initial state to start consistent initialization from
	 * @param [in] consTol Tolerance for the consistent initialization algorithm
	 * @param [in] absTol Allowed maximum residual error after consistent initialization
	 */
	void testConsistentInitialization(IUnitOperation* const unit, bool adEnabled, double* const y, double consTol, double absTol);

	/**
	 * @brief Checks the consistent initialization of the sensitivies of a generic unit operation
	 * @param [in] unit Configured unit operation
	 * @param [in] adEnabled Determines whether AD is used for calculating the Jacobian
	 * @param [in] y State of original system
	 * @param [in] yDot Time derivative of state of original system
	 * @param [in] absTol Allowed maximum residual error after consistent initialization
	 */
	void testConsistentInitializationSensitivity(cadet::IUnitOperation* const unit, bool adEnabled, double const* const y, double const* const yDot, double absTol);

} // namespace unitoperation
} // namespace test
} // namespace cadet

#endif  // CADETTEST_UNITOPTESTS_HPP_
