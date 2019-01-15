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
 * @file Automatic differentiation (AD) library integration
 */

#ifndef LIBCADET_AUTODIFF_HPP_
#define LIBCADET_AUTODIFF_HPP_

#include "cadet/cadetCompilerInfo.hpp"
#include "common/CompilerSpecific.hpp"

#if defined(ACTIVE_ADOLC)

	#define ADOLC_TAPELESS
	#define NUMBER_DIRECTIONS 80
	#include <adolc/adouble.h>

	#define ACTIVE_INIT ADOLC_TAPELESS_UNIQUE_INTERNALS

	namespace cadet
	{

		/**
		 * @brief Automatic differentiation (AD) type used as replacement for double
		 */
		class active: public adtl::adouble
		{
		public:

			typedef size_t idx_t;

			// constructors
			active() :
				adtl::adouble()
			{
			}

			active(const double& d) :
				adtl::adouble(d)
			{
			}

			active(const adtl::adouble& ad) :
				adtl::adouble(ad)
			{
			}

			// operators
			const active& operator=(const adtl::adouble& rhs)
			{
				*(static_cast<adtl::adouble*> (this)) = rhs;
				return *this;
			}

			// operators
			const active& operator=(double rhs)
			{
				*(static_cast<adtl::adouble*> (this)) = rhs;
				return *this;
			}

			inline size_t gradientSize() const { return adtl::ADOLC_numDir; }

			///todo remove this cast operator - dirty hack to cast actives to doubles in some residual functions
			// that are actually not called.
			explicit operator double() const
			{
				return getValue();
			}

			inline void fillADValue(const real_t v)
			{
				fillADValue(0, adtl::ADOLC_numDir, v);
			}
			inline void fillADValue(const size_t start, const real_t v)
			{
				fillADValue(start, adtl::ADOLC_numDir, v);
			}
			inline void fillADValue(const size_t start, const size_t end, const real_t v)
			{
				for (size_t i = start; i < end; ++i)
					setADValue(i, v);
			}

		};

		namespace ad
		{
			/**
			 * @brief Returns the maximum number of allowed AD directions (seed vectors)
			 * @return Maximum number of allowed AD directions
			 */
			inline size_t getMaxDirections() CADET_NOEXCEPT { return NUMBER_DIRECTIONS; }

			/**
			 * @brief Returns the current number of AD directions (seed vectors)
			 * @return Current number of AD directions
			 */
			inline size_t getDirections() CADET_NOEXCEPT { return adtl::ADOLC_numDir; }

			/**
			 * @brief Sets the current number of AD directions (seed vectors)
			 * @details The number of AD directions must not exceed the value returned by getMaxDirections().
			 * 
			 * @param [in] numDir Number of required AD directions
			 */
			inline void setDirections(size_t numDir)
			{
				cadet_assert(numDir <= NUMBER_DIRECTIONS);
				adtl::ADOLC_numDir = numDir; 
			}
		}
	}

#elif defined(ACTIVE_SFAD) || defined(ACTIVE_SETFAD)

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
			typedef sfad::Fwd<double, sfad::StackStorage> active;
		#else
			typedef sfad::FwdET<double, sfad::StackStorage> active;
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


#endif  // LIBCADET_AUTODIFF_HPP_
