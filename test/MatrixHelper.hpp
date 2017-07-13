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
 * Defines helper functions for checking matrices.
 */

#ifndef CADETTEST_MATRIXHELPER_HPP_
#define CADETTEST_MATRIXHELPER_HPP_

#include <vector>

namespace cadet
{

namespace test
{

/**
 * @brief Checks a given matrix against a linear array
 * @details The matrix is expected as dense matrix in row-major format.
 * @param [in] mat Dense matrix to be checked
 * @param [in] matRef Reference values
 */
inline void checkMatrixAgainstLinearArray(double const* mat, const std::vector<double>& matRef)
{
	for (unsigned int i = 0; i < matRef.size(); ++i)
		CHECK(matRef[i] == Approx(mat[i]));
}

} // namespace test
} // namespace cadet

#endif  // CADETTEST_MATRIXHELPER_HPP_
