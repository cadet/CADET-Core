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

#ifndef LIBCADET_MATHUTIL_HPP_
#define LIBCADET_MATHUTIL_HPP_

#include <cmath>

#include "cadet/cadetCompilerInfo.hpp"

namespace cadet 
{

#if defined(ACTIVE_ADOLC)
	/**
	 * @brief Squares the given value
	 * @details Calculates @f$ x * x @f$
	 * @param [in] x Value that is squared
	 * @return Squared value
	 */
    template <typename real_t> inline real_t sqr(const real_t& x) CADET_NOEXCEPT { return x * x; }
#else
	/**
	 * @brief Squares the given value
	 * @details Calculates @f$ x * x @f$
	 * @param [in] x Value that is squared
	 * @return Squared value
	 */
	inline double sqr(const double x) CADET_NOEXCEPT { return x * x; }
#endif

} // namespace cadet

#endif  // LIBCADET_MATHUTIL_HPP_
