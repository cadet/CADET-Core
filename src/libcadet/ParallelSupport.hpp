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
 * Helper functions and macros for parallelization.
 */

#ifndef LIBCADET_PARALLEL_SUPPORT_HPP_
#define LIBCADET_PARALLEL_SUPPORT_HPP_

#include "cadet/cadetCompilerInfo.hpp"
#include "common/CompilerSpecific.hpp"

#ifdef CADET_PARALLELIZE
	#define CADET_PARFOR_END )
	#define CADET_PARNODE_END )

	#include "tbb/task_arena.h"
	#include "tbb/cache_aligned_allocator.h"
	#include <vector>

	namespace cadet
	{
	namespace util
	{

		/**
		 * @brief Provides a distinct array of @c item_t for each thread (thread local storage)
		 * @details For each thread, a dynamic array of certain size is allocated.
		 *          The arrays are aligned to cache lines in order to avoid false sharing.
		 *          Note that it is not guaranteed that the exact same thread always obtains
		 *          the same array as the thread index of a thread may change over time.
		 * @tparam item_t Array item type
		 */
		template <typename item_t>
		class ThreadLocalStorage
		{
		public:
			ThreadLocalStorage() : _data(0, nullptr) { }
			~ThreadLocalStorage()
			{
				tbb::cache_aligned_allocator<item_t> allocator;
				for (item_t* it : _data)
					// Size of allocated array is not used by TBB
					allocator.deallocate(it, 0);
			}

			ThreadLocalStorage(const ThreadLocalStorage&) = delete;
			ThreadLocalStorage(ThreadLocalStorage&& mv) CADET_NOEXCEPT = default;

			ThreadLocalStorage& operator=(const ThreadLocalStorage&) = delete;

			#ifdef COMPILER_SUPPORT_NOEXCEPT_DEFAULTED_MOVE
				ThreadLocalStorage& operator=(ThreadLocalStorage&&) CADET_NOEXCEPT = default;
			#else
				ThreadLocalStorage& operator=(ThreadLocalStorage&&) = default;
			#endif

			/**
			 * @brief Allocates arrays for each thread with the given number of elements
			 * @details The number of threads can expand or shrink and the arrays are
			 *          created or released accordingly.
			 * @param [in] numThreads Number of threads
			 * @param [in] arraySize Number of elements in each array
			 */
			inline void resize(unsigned int numThreads, unsigned int arraySize)
			{
				tbb::cache_aligned_allocator<item_t> allocator;

				if (numThreads < _data.size())
				{
					// Remove superfluous arrays at the end
					for (unsigned int i = 0; i < numThreads - _data.size(); ++i)
					{
						// Size of allocated array is not used by TBB
						allocator.deallocate(_data.back(), 0);
						_data.pop_back();
					}
				}
				else if (numThreads > _data.size())
				{
					// Add some more arrays
					const unsigned int oldSize = _data.size();
					_data.resize(numThreads, nullptr);
					for (unsigned int i = oldSize; i < _data.size(); ++i)
						_data[i] = allocator.allocate(arraySize);
				}
			}

			/**
			 * @brief Allocates arrays for each thread with the given number of elements
			 * @details The number of threads can expand or shrink and the arrays are
			 *          created or released accordingly.
			 *
			 *          The number of threads is inferred from TBB.
			 * @param [in] arraySize Number of elements in each array
			 */
			inline void resize(unsigned int arraySize)
			{
				resize(tbb::this_task_arena::max_concurrency(), arraySize);
			}

			/**
			 * @brief Return pointer to array of the current thread
			 * @return Pointer to array of the current thread
			 */
			inline item_t* get() const
			{
				#ifdef CADET_DEBUG
			        const int threadIdx = tbb::this_task_arena::current_thread_index();
					cadet_assert(threadIdx != tbb::task_arena::not_initialized);
					cadet_assert(threadIdx >= 0);
					cadet_assert(threadIdx < _data.size());
				#endif
				return _data[tbb::this_task_arena::current_thread_index()];
			}

			/**
			 * @brief Return pointer to array of the given thread index
			 * @param [in] idx Thread index
			 * @return Pointer to array of the given thread index
			 */
			inline item_t* get(unsigned int idx) const
			{
				#ifdef CADET_DEBUG
					cadet_assert(idx < _data.size());
				#endif
				return _data[idx];
			}

		private:
			std::vector<item_t*> _data;
		};

		typedef ThreadLocalStorage<double> ThreadLocalArray;

		/**
		 * @brief Returns the maximum number of threads at this point in the code
		 * @details The current maximum number of threads may change from point in code
		 *          to point in code. It is governed by several TBB settings, most importantly
		 *          the first instantiation of tbb::task_scheduler_init.
		 * @return Maximum number of threads
		 */
		inline unsigned int getMaxThreads() { return tbb::this_task_arena::max_concurrency(); }

	} // namespace util
	} // namespace cadet

#else
	#define CADET_PARFOR_END
	#define CADET_PARNODE_END

	namespace cadet
	{
	namespace util
	{

		/**
		 * @brief Provides a distinct array of @c item_t for each thread (thread local storage)
		 * @details As there is only one thread, this class simply wraps a pointer to an array.
		 * @tparam item_t Array item type
		 */
		template <typename item_t>
		class ThreadLocalStorage
		{
		public:
			ThreadLocalStorage() : _data(nullptr) { }
			~ThreadLocalStorage()
			{
				delete[] _data;
			}

			ThreadLocalStorage(const ThreadLocalStorage&) = delete;
			ThreadLocalStorage(ThreadLocalStorage&& mv) CADET_NOEXCEPT = default;

			ThreadLocalStorage& operator=(const ThreadLocalStorage&) = delete;

			#ifdef COMPILER_SUPPORT_NOEXCEPT_DEFAULTED_MOVE
				ThreadLocalStorage& operator=(ThreadLocalStorage&&) CADET_NOEXCEPT = default;
			#else
				ThreadLocalStorage& operator=(ThreadLocalStorage&&) = default;
			#endif

			/**
			 * @brief Allocates arrays for each thread with the given number of elements
			 * @param [in] numThreads Number of threads
			 * @param [in] arraySize Number of elements in each array
			 */
			inline void resize(unsigned int numThreads, unsigned int arraySize)
			{
				if (!_data)
					_data = new item_t[arraySize];
			}

			/**
			 * @brief Allocates arrays for each thread with the given number of elements
			 * @param [in] numThreads Number of threads
			 * @param [in] arraySize Number of elements in each array
			 */
			inline void resize(unsigned int arraySize)
			{
				resize(1u, arraySize);
			}

			/**
			 * @brief Return pointer to array of the current thread
			 * @return Pointer to array of the current thread
			 */
			inline item_t* get() const
			{
				return _data;
			}

			/**
			 * @brief Return pointer to array of the given thread index
			 * @param [in] idx Thread index
			 * @return Pointer to array of the given thread index
			 */
			inline item_t* get(unsigned int idx) const
			{
				return _data;
			}

		private:
			item_t* _data;
		};

		typedef ThreadLocalStorage<double> ThreadLocalArray;

		/**
		 * @brief Returns the maximum number of threads at this point in the code
		 * @details This is just one in the single threaded case.
		 * @return Maximum number of threads
		 */
		inline unsigned int getMaxThreads() { return 1; }

	} // namespace util
	} // namespace cadet

#endif

#endif  // LIBCADET_PARALLEL_SUPPORT_HPP_
