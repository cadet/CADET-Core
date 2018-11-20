// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2018: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file 
 * Defines utility functions for tests.
 */

#ifndef CADETTEST_TESTUTILS_HPP_
#define CADETTEST_TESTUTILS_HPP_

#include <functional>
#include <algorithm>

namespace cadet
{

namespace test
{

namespace util
{

	/**
	 * @brief Fills the vector with the values of a given function
	 * @details The function @p f uses the current index to assign a value.
	 * @param [out] y Populated vector
	 * @param [in] f Function for computing the content of the vector
	 * @param [in] numDofs Size of the vector
	 */
	inline void populate(double* y, std::function<double(unsigned int)> f, unsigned int numDofs)
	{
		for (unsigned int i = 0; i < numDofs; ++i)
			y[i] = f(i);
	}

	/**
	 * @brief Sequentially places a sequence multiple times into an array
	 * @details Writes a given sequence multiple times into an array by concatenation.
	 *          Writes a total of @c size * times elements.
	 * @param [out] dest Target array
	 * @param [in] src Source sequence
	 * @param [in] size Length of the source sequence
	 * @param [in] times Number of repetitions of the source sequence
	 */
	inline void repeat(double* dest, double const* src, unsigned int size, unsigned int times)
	{
		for (unsigned int i = 0; i < times; ++i, dest += size)
		{
			std::copy(src, src + size, dest);
		}
	}

} // namespace util
} // namespace test
} // namespace cadet

#endif  // CADETTEST_TESTUTILS_HPP_
