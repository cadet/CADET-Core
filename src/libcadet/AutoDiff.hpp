// =============================================================================
//  CADET
//  
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
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

#include "CompileTimeConfig.hpp"
#include "cadet/cadetCompilerInfo.hpp"
#include "common/CompilerSpecific.hpp"

#ifdef ENABLE_DG
	#include <Eigen/Core>
#endif

#include <type_traits>

#if defined(ACTIVE_SFAD) || defined(ACTIVE_SETFAD)

	#if !defined(SFAD_DEFAULT_DIR) && defined(NUM_MAX_AD_DIRS)
		#define SFAD_DEFAULT_DIR NUM_MAX_AD_DIRS
	#endif

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
			inline std::size_t getMaxDirections() CADET_NOEXCEPT { return SFAD_DEFAULT_DIR; }

			/**
			 * @brief Returns the current number of AD directions (seed vectors)
			 * @return Current number of AD directions
			 */
			inline std::size_t getDirections() CADET_NOEXCEPT { return sfad::getGradientSize(); }

			/**
			 * @brief Sets the current number of AD directions (seed vectors)
			 * @details The number of AD directions must not exceed the value returned by getMaxDirections().
			 * 
			 * @param [in] numDir Number of required AD directions
			 */
			inline void setDirections(std::size_t n)
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
	struct DoubleActivePromoterImpl { };

	template <>
	struct DoubleActivePromoterImpl<cadet::active, cadet::active> { typedef cadet::active type; };

	template <>
	struct DoubleActivePromoterImpl<cadet::active, double> { typedef cadet::active type; };

	template <>
	struct DoubleActivePromoterImpl<double, cadet::active> { typedef cadet::active type; };

	template <>
	struct DoubleActivePromoterImpl<double, double> { typedef double type; };

	/**
	 * @brief Selects the @c active type between @c double and @c active
	 * @tparam A Type A
	 * @tparam B Type B
	 */
	template <typename A, typename B>
	struct DoubleActivePromoter { typedef typename DoubleActivePromoterImpl<std::decay_t<A>, std::decay_t<B>>::type type; };

	template <typename A>
	using ActivePromoter = DoubleActivePromoter<A, double>;

	/**
	 * @brief Selects the @c double type between @c double and @c active
	 * @tparam A Type A
	 * @tparam B Type B
	 */
	template <typename A, typename B>
	struct DoubleActiveDemoterImpl { };

	template <>
	struct DoubleActiveDemoterImpl<cadet::active, cadet::active> { typedef cadet::active type; };

	template <>
	struct DoubleActiveDemoterImpl<cadet::active, double> { typedef double type; };

	template <>
	struct DoubleActiveDemoterImpl<double, cadet::active> { typedef double type; };

	template <>
	struct DoubleActiveDemoterImpl<double, double> { typedef double type; };

	/**
	 * @brief Selects the @c double type between @c double and @c active
	 * @tparam A Type A
	 * @tparam B Type B
	 */
	template <typename A, typename B>
	struct DoubleActiveDemoter { typedef typename DoubleActiveDemoterImpl<std::decay_t<A>, std::decay_t<B>>::type type; };

	template <typename A>
	using DoubleDemoter = DoubleActiveDemoter<A, cadet::active>;

	/**
	 * @brief Selects value type @c double or reference type @c active&
	 * @tparam A Base type (i.e., @c double or @c active)
	 */
	template <typename A>
	struct ActiveRefOrDouble { };

	template <>
	struct ActiveRefOrDouble<cadet::active> { typedef cadet::active& type; };

	template <>
	struct ActiveRefOrDouble<const cadet::active> { typedef const cadet::active& type; };

	template <>
	struct ActiveRefOrDouble<double> { typedef double type; };

	template <>
	struct ActiveRefOrDouble<const double> { typedef const double type; };
}

#ifdef ENABLE_DG
	namespace Eigen {

		template<> struct NumTraits<cadet::active>
			: NumTraits<double> // permits to get the epsilon, dummy_precision, lowest, highest functions
		{
			typedef cadet::active Real;
			typedef cadet::active NonInteger;
			typedef cadet::active Nested;

			enum {
				IsComplex = 0,
				IsInteger = 0,
				IsSigned = 1,
				RequireInitialization = 1,
				ReadCost = 1,
				AddCost = 3,
				MulCost = 3
			};
		};

			// specify return types concerning active double scalar operations for eigen
			template<>
			struct ScalarBinaryOpTraits<cadet::active, cadet::active, internal::scalar_product_op<cadet::active, cadet::active>> {
				typedef cadet::active ReturnType;
			};
			template<>
			struct ScalarBinaryOpTraits<cadet::active, double, internal::scalar_product_op<cadet::active, double>> {
				typedef cadet::active ReturnType;
			};
			template<>
			struct ScalarBinaryOpTraits<double, cadet::active, internal::scalar_product_op<double, cadet::active>> {
				typedef cadet::active ReturnType;
			};
	}
#endif

#endif  // LIBCADET_AUTODIFF_HPP_
