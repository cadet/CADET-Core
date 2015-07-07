// =============================================================================
//  SFAD - Simple Forward Automatic Differentiation
//  
//  Copyright © 2015: Samuel Leweke¹
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
	#define SFAD_GLOBAL_GRAD_SIZE size_t sfad::detail::globalGradSize = SFAD_DEFAULT_DIR;
#endif

namespace sfad
{
	namespace detail
	{
		extern size_t globalGradSize;
	}

	static void setGradientSize(const size_t n)
	{
		detail::globalGradSize = n;
	}

	static size_t getGradientSize()
	{
		return detail::globalGradSize;
	}	


	template <typename real_t>
	class HeapStorage
	{
	public:
		HeapStorage() : _grad(new real_t[detail::globalGradSize]) { }
		HeapStorage(HeapStorage<real_t>&& other) : _grad(other._grad)
		{
			other._grad = 0;
		}
		HeapStorage(const HeapStorage<real_t>& cpy) : _grad(new real_t[detail::globalGradSize])
		{
			copyGradient(cpy._grad);
		}

		~HeapStorage()
		{
			if (_grad)
				delete[] _grad;
		}

		void resizeGradient()
		{
			if (_grad)
				delete[] _grad;

			_grad = new real_t[detail::globalGradSize];
			for (size_t i = 0; i < detail::globalGradSize; ++i)
				_grad[i] = real_t(0);
		}

	protected:
		real_t* _grad;

		void moveAssign(HeapStorage&& other)
		{
			if (_grad)
				delete[] _grad;

			_grad = other._grad;
			other._grad = 0;
		}

		void copyGradient(real_t const* const cpy)
		{
			for (size_t i = 0; i < detail::globalGradSize; ++i)
			{
				_grad[i] = cpy[i];
			}
		}
	};


	template <typename real_t>
	class StackStorage
	{
	public:
		StackStorage() { }
		StackStorage(const StackStorage<real_t>& cpy) { copyGradient(cpy._grad); }
//		StackStorage(StackStorage<real_t>&& other) : _grad(std::move(other._grad)) { }
		StackStorage(StackStorage<real_t>&& other) { copyGradient(other._grad); }

		~StackStorage() { }

		void resizeGradient() { }

	protected:
		real_t _grad[SFAD_DEFAULT_DIR];

		void moveAssign(StackStorage&& other)
		{
//			_grad = std::move(other._grad);
			copyGradient(other._grad);
		}

		void copyGradient(real_t const* const cpy)
		{
			for (size_t i = 0; i < detail::globalGradSize; ++i)
			{
				_grad[i] = cpy[i];
			}
		}
	};

}

#endif
