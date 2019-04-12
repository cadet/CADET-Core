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
 * Defines tests for reaction models.
 */

#ifndef CADETTEST_REACTIONODELTEST_HPP_
#define CADETTEST_REACTIONODELTEST_HPP_

#include <limits>

namespace cadet
{

namespace test
{

namespace reaction
{

	/**
	 * @brief Checks the analytic Jacobians of the dynamic reaction model against AD
	 * @param [in] modelName Name of the reaction model
	 * @param [in] nComp Number of components
	 * @param [in] nBound Array with number of bound states for each component
	 * @param [in] config JSON string with reaction model parameters
	 * @param [in] point Liquid phase and solid phase values to check Jacobian at
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testDynamicJacobianAD(const char* modelName, unsigned int nComp, unsigned int const* nBound, const char* config, double const* point, double absTol = 0.0, double relTol = std::numeric_limits<float>::epsilon() * 100.0);

} // namespace reaction
} // namespace test
} // namespace cadet

#endif  // CADETTEST_REACTIONODELTEST_HPP_
