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
 * Provides a very simple memory pool for arrays of different objects.
 */

#ifndef LIBCADET_MEMORYPOOL_HPP_
#define LIBCADET_MEMORYPOOL_HPP_

#include "common/CompilerSpecific.hpp"
#include <type_traits>

namespace cadet
{
	/**
	 * @brief Manages a memory pool for arrays of various types
	 * @details An ArrayPool manages a block of memory which can be used to create arrays of
	 *          various types. It is initialized to allocate a block of memory big enough to
	 *          hold the maximum number of elements of the biggest type. The user can then
	 *          create arrays of different types and sizes without actually allocating new
	 *          memory (saves perfomance). Note that the elements have to be default constructible.
	 */
	class ArrayPool
	{
	public:

		/**
		 * @brief Creates an empty ArrayPool
		 */
		ArrayPool() : _mem(nullptr), _numElements(0)
		{
#ifdef DEBUG
			_capacity = 0;
#endif
		}

		/**
		 * @brief Creates an ArrayPool with the given size in bytes
		 * @details The size of the pool should be computed by taking @c sizeof(BiggestType) * maxNumElements
		 * 
		 * @param [in] maxBytes Size of the pool in bytes
		 */
		ArrayPool(unsigned int maxBytes) : _mem(new char[maxBytes]), _numElements(0)
		{
#ifdef DEBUG
			_capacity = maxBytes;
#endif
		}

		~ArrayPool() CADET_NOEXCEPT { delete[] _mem; }

		/**
		 * @brief Resizes the memory pool to the given size
		 * @details Already created arrays are invalid as the underlying memory pool is reallocated.
		 *          The user has to make sure that the created array is destroyed prior to calling resize().
		 * 
		 * @param [in] maxBytes Size of the pool in bytes
		 */
		void resize(unsigned int maxBytes)
		{
			delete[] _mem;
			_mem = new char[maxBytes];
			_numElements = 0;

#ifdef DEBUG
			_capacity = maxBytes;
#endif
		}

		/**
		 * @brief Creates an array in the pool
		 * @details Since the array is always created in the same memory, previous arrays will be overwritten.
		 *          Array elements are default constructed.
		 * 
		 * @param [in] n Number of elements
		 * @tparam T Type of elements
		 * @return Pointer to array of the requested size and type
		 */
		template <typename T>
		inline typename std::enable_if<!std::is_arithmetic<T>::value, T* const>::type create(unsigned int n)
		{
			cadet_assert(sizeof(T) * n <= _capacity);

			// Call constructor on every element
			T* cur = reinterpret_cast<T*>(_mem);
			for (unsigned int i = 0; i < n; ++i)
			{
				// Default constructor
				new(cur) T;
				++cur;
			}

			_numElements = n;
			return reinterpret_cast<T*>(_mem);
		}

		/**
		 * @brief Destroys a currently active array
		 * @details The destructor of the array elements is called
		 * 
		 * @tparam T Type of elements
		 */
		template <typename T>
		inline typename std::enable_if<!std::is_arithmetic<T>::value, void>::type destroy()
		{
			// Call destructor on every element
			T* cur = reinterpret_cast<T*>(_mem);
			for (unsigned int i = 0; i < _numElements; ++i)
			{
				cur->~T();
				++cur;
			}
			_numElements = 0;
		}

		/**
		 * @brief Creates an array of integral / arithmetic types (e.g., double, int) in the pool
		 * @details Since the array is always created in the same memory, previous arrays will be overwritten.
		 *          Array elements are initialized to zero.
		 * 
		 * @param [in] n Number of elements
		 * @tparam T Type of elements
		 * @return Pointer to array of the requested size and type
		 */
		template <typename T>
		inline typename std::enable_if<std::is_arithmetic<T>::value, T* const>::type create(unsigned int n)
		{
			cadet_assert(sizeof(T) * n <= _capacity);

			// Arithmetic types do not need constructor, but initialization to zero
			T* cur = reinterpret_cast<T*>(_mem);
			for (unsigned int i = 0; i < n; ++i)
			{
				*cur = T(0);
				++cur;
			}

			_numElements = n;
			return reinterpret_cast<T*>(_mem);
		}

		/**
		 * @brief Destroys a currently active array
		 * @tparam T Type of elements
		 */
		template <typename T>
		inline typename std::enable_if<std::is_arithmetic<T>::value, void>::type destroy()
		{
			// Call to destructor is not required
			_numElements = 0;
		}

		inline const unsigned int numElements() const { return _numElements; }

	protected:
		char* _mem; //<! Memory block
		unsigned int _numElements; //<! Current number of created elements
#ifdef DEBUG
		unsigned int _capacity; //<! Capacity of the pool, only available in debug builds
#endif
	};

} // namespace cadet

#endif  // LIBCADET_MEMORYPOOL_HPP_
