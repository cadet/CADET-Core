// =============================================================================
//  SFAD - Simple Forward Automatic Differentiation
//  
//  Copyright © 2015-2018: Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#ifndef _SFAD_COMMON_HPP_
#define _SFAD_COMMON_HPP_

#ifndef SFAD_DEFAULT_DIR
	#define SFAD_DEFAULT_DIR 80
#endif

#ifndef SFAD_GLOBAL_GRAD_SIZE
	#define SFAD_GLOBAL_GRAD_SIZE std::size_t sfad::detail::globalGradSize = SFAD_DEFAULT_DIR;
#endif

// Improve branch prediction by marking likely and unlikely execution paths
// Example:
//     if (some_unlikely_condition) { ... }
// transforms to
//     if (sfad_unlikely(some_unlikely_condition)) { ... }
#ifdef __GNUC__
	#define sfad_likely(x) __builtin_expect(!!(x), 1)
	#define sfad_unlikely(x) __builtin_expect(!!(x), 0)
#endif

#ifndef sfad_likely
	#define sfad_likely(x) (x)
	#define sfad_unlikely(x) (x)
#endif

#if defined(__clang__) && defined(__apple_build_version__)
	// Apple Clang
	#if ((__clang_major__ * 100) + __clang_minor__) >= 400 && __has_feature(cxx_noexcept)
		#define SFAD_COMPILER_CXX_NOEXCEPT 1
	#else
		#define SFAD_COMPILER_CXX_NOEXCEPT 0
	#endif

#elif defined(__clang__)
	// Clang
	#if ((__clang_major__ * 100) + __clang_minor__) >= 304 && __has_feature(cxx_noexcept)
		#define SFAD_COMPILER_CXX_NOEXCEPT 1
	#else
		#define SFAD_COMPILER_CXX_NOEXCEPT 0
	#endif

#elif defined(__INTEL_COMPILER) || defined(__ICC)
	// Intel icpc
	#if (__INTEL_COMPILER >= 1400)
		#define SFAD_COMPILER_CXX_NOEXCEPT 1
	#else
		#define SFAD_COMPILER_CXX_NOEXCEPT 0
	#endif

#elif defined(__GNUC__)
	// GNU GCC
	#if (__GNUC__ * 100 + __GNUC_MINOR__) >= 406 && (__cplusplus >= 201103L || (defined(__GXX_EXPERIMENTAL_CXX0X__) && __GXX_EXPERIMENTAL_CXX0X__))
		#define SFAD_COMPILER_CXX_NOEXCEPT 1
	#else
		#define SFAD_COMPILER_CXX_NOEXCEPT 0
	#endif

#elif defined(_MSC_VER)
	// MS Visual C++
	#if _MSC_VER >= 1900
		#define SFAD_COMPILER_CXX_NOEXCEPT 1
	#else
		#define SFAD_COMPILER_CXX_NOEXCEPT 0
	#endif

#endif

#if SFAD_COMPILER_CXX_NOEXCEPT
	#define SFAD_NOEXCEPT noexcept
	#define SFAD_NOEXCEPT_EXPR(X) noexcept(X)
#else
	#define SFAD_NOEXCEPT
	#define SFAD_NOEXCEPT_EXPR(X)
#endif


#include <algorithm>

namespace sfad
{
	namespace detail
	{
		extern std::size_t globalGradSize;
	}

	inline void setGradientSize(const std::size_t n) SFAD_NOEXCEPT
	{
		detail::globalGradSize = n;
	}

	inline std::size_t getGradientSize() SFAD_NOEXCEPT
	{
		return detail::globalGradSize;
	}	


	template <typename real_t>
	class HeapStorage
	{
	public:
		HeapStorage() : _grad(new real_t[detail::globalGradSize]) { }
		HeapStorage(HeapStorage<real_t>&& other) SFAD_NOEXCEPT : _grad(other._grad)
		{
			other._grad = 0;
		}
		HeapStorage(const HeapStorage<real_t>& cpy) : _grad(new real_t[detail::globalGradSize])
		{
			copyGradient(cpy._grad);
		}

		~HeapStorage() SFAD_NOEXCEPT
		{
			delete[] _grad;
		}

		void resizeGradient()
		{
			delete[] _grad;

			_grad = new real_t[detail::globalGradSize];
			std::fill(_grad, _grad + detail::globalGradSize, real_t(0));
		}

	protected:
		real_t* _grad;

		void moveAssign(HeapStorage&& other) SFAD_NOEXCEPT
		{
			delete[] _grad;

			_grad = other._grad;
			other._grad = 0;
		}

		void copyGradient(real_t const* const cpy)
		{
			std::copy(cpy, cpy + detail::globalGradSize, _grad);
		}
	};


	template <typename real_t>
	class StackStorage
	{
	public:
		StackStorage() SFAD_NOEXCEPT { }
		StackStorage(const StackStorage<real_t>& cpy) SFAD_NOEXCEPT { copyGradient(cpy._grad); }
//		StackStorage(StackStorage<real_t>&& other) : _grad(std::move(other._grad)) { }
		StackStorage(StackStorage<real_t>&& other) SFAD_NOEXCEPT { copyGradient(other._grad); }

		~StackStorage() SFAD_NOEXCEPT { }

		void resizeGradient() SFAD_NOEXCEPT { }

	protected:
		real_t _grad[SFAD_DEFAULT_DIR];

		void moveAssign(StackStorage&& other) SFAD_NOEXCEPT
		{
//			_grad = std::move(other._grad);
			copyGradient(other._grad);
		}

		void copyGradient(real_t const* const cpy)
		{
			std::copy(cpy, cpy + detail::globalGradSize, _grad);
		}
	};

}

#endif
