// =============================================================================
//  CADET
//  
//  Copyright © 2008-2020: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file Automatic differentiation (AD) library integration
 */

#ifndef LIBCADET_AUTODIFF_HPP_
#define LIBCADET_AUTODIFF_HPP_

#include "cadet/cadetCompilerInfo.hpp"
#include "common/CompilerSpecific.hpp"

#if defined(ACTIVE_SFAD) || defined(ACTIVE_SETFAD)

	#define SFAD_DEFAULT_DIR 80

	#if defined(ACTIVE_SFAD)
		#include "sfad.hpp"
	#else
		#include "setfad.hpp"
	#endif

	#define ACTIVE_INIT SFAD_GLOBAL_GRAD_SIZE

	namespace cadet
	{
		
		#if defined(ACTIVE_SFAD)
			typedef sfad::Fwd<double> active;
		#else
			typedef sfad::FwdET<double> active;
		#endif

		namespace ad
		{
			/**
			 * @brief Returns the maximum number of allowed AD directions (seed vectors)
			 * @return Maximum number of allowed AD directions
			 */
			inline size_t getMaxDirections() CADET_NOEXCEPT { return SFAD_DEFAULT_DIR; }

			/**
			 * @brief Returns the current number of AD directions (seed vectors)
			 * @return Current number of AD directions
			 */
			inline size_t getDirections() CADET_NOEXCEPT { return sfad::getGradientSize(); }

			/**
			 * @brief Sets the current number of AD directions (seed vectors)
			 * @details The number of AD directions must not exceed the value returned by getMaxDirections().
			 * 
			 * @param [in] numDir Number of required AD directions
			 */
			inline void setDirections(size_t n)
			{
				cadet_assert(n <= SFAD_DEFAULT_DIR);
				sfad::setGradientSize(n);
			}
		}
	}

#else

	#error No active data type defined!

#endif  // #if defined ACTIVE_


namespace cadet
{
	/**
	 * @brief Selects the @c active type between @c double and @c active
	 * @tparam A Type A
	 * @tparam B Type B
	 */
	template <typename A, typename B>
	struct DoubleActivePromoter { };

	template <>
	struct DoubleActivePromoter<cadet::active, cadet::active> { typedef cadet::active type; };

	template <>
	struct DoubleActivePromoter<cadet::active, double> { typedef cadet::active type; };

	template <>
	struct DoubleActivePromoter<double, cadet::active> { typedef cadet::active type; };

	template <>
	struct DoubleActivePromoter<double, double> { typedef double type; };

	template <typename A>
	using ActivePromoter = DoubleActivePromoter<A, double>;

	/**
	 * @brief Selects the @c double type between @c double and @c active
	 * @tparam A Type A
	 * @tparam B Type B
	 */
	template <typename A, typename B>
	struct DoubleActiveDemoter { };

	template <>
	struct DoubleActiveDemoter<cadet::active, cadet::active> { typedef cadet::active type; };

	template <>
	struct DoubleActiveDemoter<cadet::active, double> { typedef double type; };

	template <>
	struct DoubleActiveDemoter<double, cadet::active> { typedef double type; };

	template <>
	struct DoubleActiveDemoter<double, double> { typedef double type; };

	template <typename A>
	using DoubleDemoter = DoubleActiveDemoter<A, cadet::active>;
}

#endif  // LIBCADET_AUTODIFF_HPP_
