// =============================================================================
//  SFAD - Simple Forward Automatic Differentiation
//  
//  Copyright © 2015-2017: Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#ifndef _SFAD_MAIN_HPP_
#define _SFAD_MAIN_HPP_

#include <cmath>
#include <algorithm>
#include <limits>
#include <utility>

#include "sfad-common.hpp"

namespace sfad
{
	template <typename real_t, template <class T> class storage_t>
	class Fwd : public storage_t<real_t>
	{
	public:
		typedef std::size_t idx_t;

		Fwd() SFAD_NOEXCEPT_EXPR(noexcept(storage_t<real_t>())) : storage_t<real_t>(), _val(0)
		{
			setADValue(real_t(0));
		}
		Fwd(const real_t val) SFAD_NOEXCEPT_EXPR(noexcept(storage_t<real_t>())) : storage_t<real_t>(), _val(val)
		{
			setADValue(real_t(0));
		}
		Fwd(const real_t val, real_t const* const grad) SFAD_NOEXCEPT_EXPR(noexcept(storage_t<real_t>())) : storage_t<real_t>(), _val(val)
		{
			storage_t<real_t>::copyGradient(grad);
		}
		Fwd(const Fwd<real_t, storage_t>& cpy) SFAD_NOEXCEPT_EXPR(noexcept(storage_t<real_t>(cpy))) : storage_t<real_t>(cpy), _val(cpy._val) { }
		Fwd(Fwd<real_t, storage_t>&& other) SFAD_NOEXCEPT_EXPR(noexcept(storage_t<real_t>(std::move(other)))) : storage_t<real_t>(std::move(other)), _val(std::move(other._val)) { }

		~Fwd() { }

		Fwd<real_t, storage_t>& operator=(Fwd<real_t, storage_t>&& other) SFAD_NOEXCEPT
		{
			_val = std::move(other._val);
			storage_t<real_t>::moveAssign(std::move(other));

			return *this;
		}

		Fwd<real_t, storage_t>& operator=(const Fwd<real_t, storage_t>& other)
		{
			if (sfad_likely(this != &other))
			{
				_val = other._val;
				storage_t<real_t>::copyGradient(other._grad);
			}

			return *this;
		}

		const idx_t gradientSize() const SFAD_NOEXCEPT { return detail::globalGradSize; }

		template<typename T, template <class T2> class s_t> friend void swap (Fwd<T, s_t>& x, Fwd<T, s_t>& y) SFAD_NOEXCEPT;

		// ADOL-C compatibility

		inline real_t getValue() SFAD_NOEXCEPT { return _val; }
		inline const real_t getValue() const SFAD_NOEXCEPT { return _val; }
		inline void setValue(const real_t v) SFAD_NOEXCEPT { _val = v; }
		
		inline real_t getADValue(const idx_t idx) { return storage_t<real_t>::_grad[idx]; }
		inline const real_t getADValue(const idx_t idx) const { return storage_t<real_t>::_grad[idx]; }
		inline void setADValue(const idx_t idx, const real_t v) { storage_t<real_t>::_grad[idx] = v; }
		inline void setADValue(const real_t v)
		{
			fillADValue(v);
		}

		inline void fillADValue(const real_t v)
		{
			fillADValue(0, detail::globalGradSize, v);
		}
		inline void fillADValue(const idx_t start, const real_t v)
		{
			fillADValue(start, detail::globalGradSize, v);
		}
		inline void fillADValue(const idx_t start, const idx_t end, const real_t v)
		{
			std::fill(storage_t<real_t>::_grad + start, storage_t<real_t>::_grad + end, v);
		}

		// Modern C++ accessor

		inline real_t& operator[](const idx_t idx) { return storage_t<real_t>::_grad[idx]; }
		inline const real_t operator[](const idx_t idx) const { return storage_t<real_t>::_grad[idx]; }

		explicit operator real_t() const SFAD_NOEXCEPT { return _val; }

		// Operators with non-temporary results
		
		// Assignment
		inline Fwd<real_t, storage_t>& operator=(const real_t v)
		{
			_val = v;
			setADValue(real_t(0));

			return *this;
		}

		// Addition
		inline Fwd<real_t, storage_t>& operator+=(const real_t v)
		{
			_val += v;
			return *this;
		}

		inline Fwd<real_t, storage_t>& operator+=(const Fwd<real_t, storage_t>& a)
		{
			_val += a._val;
			for (idx_t i = 0; i < detail::globalGradSize; ++i)
				storage_t<real_t>::_grad[i] += a._grad[i];

			return *this;
		}

		// Substraction
		inline Fwd<real_t, storage_t>& operator-=(const real_t v)
		{
			_val -= v;
			return *this;
		}

		inline Fwd<real_t, storage_t>& operator-=(const Fwd<real_t, storage_t>& a)
		{
			_val -= a._val;
			for (idx_t i = 0; i < detail::globalGradSize; ++i)
				storage_t<real_t>::_grad[i] -= a._grad[i];

			return *this;
		}

		// Multiplication
		inline Fwd<real_t, storage_t>& operator*=(const real_t v)
		{
			_val *= v;
			for (idx_t i = 0; i < detail::globalGradSize; ++i)
				storage_t<real_t>::_grad[i] *= v;
			return *this;
		}

		inline Fwd<real_t, storage_t>& operator*=(const Fwd<real_t, storage_t>& a)
		{
			for (idx_t i = 0; i < detail::globalGradSize; ++i)
				storage_t<real_t>::_grad[i] = a._val * storage_t<real_t>::_grad[i] + _val * a._grad[i];

			_val *= a._val;
			return *this;
		}
		
		// Division
		inline Fwd<real_t, storage_t>& operator/=(const real_t v)
		{
			_val /= v;
			for (idx_t i = 0; i < detail::globalGradSize; ++i)
				storage_t<real_t>::_grad[i] /= v;
			return *this;
		}

		inline Fwd<real_t, storage_t>& operator/=(const Fwd<real_t, storage_t>& a)
		{
			for (idx_t i = 0; i < detail::globalGradSize; ++i)
//				_grad[i] = (_grad[i] - _val / a._val * a._grad[i]) / a._val;
				storage_t<real_t>::_grad[i] = (storage_t<real_t>::_grad[i] * a._val - _val * a._grad[i]) / (a._val * a._val);

			_val /= a._val;
			return *this;
		}

		// Comparisons
		inline bool operator!=(const Fwd<real_t, storage_t>& v) const SFAD_NOEXCEPT { return v != _val; }
		inline bool operator!=(const real_t v) const SFAD_NOEXCEPT { return v != _val; }
		inline friend bool operator!=(const real_t v, const Fwd<real_t, storage_t>& a) SFAD_NOEXCEPT { return v != a._val; }

		inline bool operator==(const Fwd<real_t, storage_t>& v) const SFAD_NOEXCEPT { return v == _val; }
		inline bool operator==(const real_t v) const SFAD_NOEXCEPT { return v == _val; }
		inline friend bool operator==(const real_t v, const Fwd<real_t, storage_t>& a) SFAD_NOEXCEPT { return v == a._val; }

		inline bool operator<=(const Fwd<real_t, storage_t>& v) const SFAD_NOEXCEPT { return _val <= v._val; }
		inline bool operator<=(const real_t v) const SFAD_NOEXCEPT { return _val <= v; }
		inline friend bool operator<=(const real_t v, const Fwd<real_t, storage_t>& a) SFAD_NOEXCEPT { return v <= a._val; }

		inline bool operator>=(const Fwd<real_t, storage_t>& v) const SFAD_NOEXCEPT { return _val >= v._val; }
		inline bool operator>=(const real_t v) const SFAD_NOEXCEPT { return _val >= v; }
		inline friend bool operator>= (const real_t v, const Fwd<real_t, storage_t>& a) SFAD_NOEXCEPT { return v >= a._val; }

		inline bool operator>(const Fwd<real_t, storage_t>& v) const SFAD_NOEXCEPT { return _val > v._val; }
		inline bool operator>(const real_t v) const SFAD_NOEXCEPT { return _val > v; }
		inline friend bool operator>(const real_t v, const Fwd<real_t, storage_t>& a) SFAD_NOEXCEPT { return v > a._val; }

		inline bool operator<(const Fwd<real_t, storage_t>& v) const SFAD_NOEXCEPT { return _val < v._val; }
		inline bool operator<(const real_t v) const SFAD_NOEXCEPT { return _val < v; }
		inline friend bool operator<(const real_t v, const Fwd<real_t, storage_t>& a) SFAD_NOEXCEPT { return v < a._val; }

		// Operators with temporary results
		
		// Unary sign
		inline Fwd<real_t, storage_t> operator-() const
		{
			Fwd<real_t, storage_t> cpy(-_val, false);
			for (idx_t i = 0; i < detail::globalGradSize; ++i)
				cpy._grad[i] = -storage_t<real_t>::_grad[i];

			return cpy;
		}

		inline Fwd<real_t, storage_t> operator+() const { return *this; }

		// Addition
		inline Fwd<real_t, storage_t> operator+(const real_t v) const
		{
			return Fwd<real_t, storage_t>(_val + v, storage_t<real_t>::_grad);
		}

		inline Fwd<real_t, storage_t> operator+(const Fwd<real_t, storage_t>& a) const
		{
			Fwd<real_t, storage_t> cpy(_val + a._val, false);
			for (idx_t i = 0; i < detail::globalGradSize; ++i)
				cpy._grad[i] = storage_t<real_t>::_grad[i] + a._grad[i];
			return cpy;
		}

		inline friend Fwd<real_t, storage_t> operator+(const real_t v, const Fwd<real_t, storage_t>& a)
		{
			return Fwd<real_t, storage_t>(v + a._val, a._grad);
		}
		
		// Substraction
		inline Fwd<real_t, storage_t> operator-(const real_t v) const
		{
			return Fwd<real_t, storage_t>(_val - v, storage_t<real_t>::_grad);
		}

		inline Fwd<real_t, storage_t> operator-(const Fwd<real_t, storage_t>& a) const
		{
			Fwd<real_t, storage_t> cpy(_val - a._val, false);
			for (idx_t i = 0; i < detail::globalGradSize; ++i)
				cpy._grad[i] = storage_t<real_t>::_grad[i] - a._grad[i];
			return cpy;
		}

		inline friend Fwd<real_t, storage_t> operator-(const real_t v, const Fwd<real_t, storage_t>& a)
		{
			Fwd<real_t, storage_t> res(v - a._val, false);
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = -a._grad[i];
			return res;
		}
		
		// Multiplication
		inline Fwd<real_t, storage_t> operator*(const real_t v) const
		{
			Fwd<real_t, storage_t> res(_val * v);
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = v * storage_t<real_t>::_grad[i];
			return res;
		}

		inline Fwd<real_t, storage_t> operator*(const Fwd<real_t, storage_t>& a) const
		{
			Fwd<real_t, storage_t> cpy(_val * a._val, false);
			for (idx_t i = 0; i < detail::globalGradSize; ++i)
				cpy._grad[i] = a._val * storage_t<real_t>::_grad[i] + _val * a._grad[i];
			return cpy;
		}

		inline friend Fwd<real_t, storage_t> operator*(const real_t v, const Fwd<real_t, storage_t>& a)
		{
			Fwd<real_t, storage_t> res(v * a._val, false);
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = v * a._grad[i];
			return res;
		}
	
		// Division
		inline Fwd<real_t, storage_t> operator/(const real_t v) const
		{
			Fwd<real_t, storage_t> res(_val / v, false);
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = storage_t<real_t>::_grad[i] / v;
			return res;
		}

		inline Fwd<real_t, storage_t> operator/(const Fwd<real_t, storage_t>& a) const
		{
			Fwd<real_t, storage_t> res(_val / a._val, false);
			for (idx_t i = 0; i < detail::globalGradSize; ++i)
//				res._grad[i] = (storage_t<real_t>::_grad[i] - _val / a._val * a._grad[i]) / a._val;
				res._grad[i] = (storage_t<real_t>::_grad[i] * a._val - _val * a._grad[i]) / (a._val * a._val);
			return res;
		}

		inline friend Fwd<real_t, storage_t> operator/(const real_t v, const Fwd<real_t, storage_t>& a)
		{
			Fwd<real_t, storage_t> res(v / a._val, false);
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
//				res._grad[i] = -(v / (a._val * a._val) * a._grad[i]);
				res._grad[i] = -v * a._grad[i] / (a._val * a._val);
			return res;
		}

		// Math functions
		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> exp(const Fwd<T, s_t> &a);
		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> log(const Fwd<T, s_t> &a);
		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> log10(const Fwd<T, s_t> &a);
		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> sqrt(const Fwd<T, s_t> &a);
		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> sqr(const Fwd<T, s_t> &a);

		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> sin(const Fwd<T, s_t> &a);
		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> cos(const Fwd<T, s_t> &a);
		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> tan(const Fwd<T, s_t> &a);
		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> asin(const Fwd<T, s_t> &a);
		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> acos(const Fwd<T, s_t> &a);
		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> atan(const Fwd<T, s_t> &a);

		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> pow(const Fwd<T, s_t> &a, T v);
		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> pow(T v, const Fwd<T, s_t> &a);
		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> pow(const Fwd<T, s_t> &a, const Fwd<T, s_t> &b);

		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> sinh(const Fwd<T, s_t> &a);
		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> cosh(const Fwd<T, s_t> &a);
		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> tanh(const Fwd<T, s_t> &a);

		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> fabs(const Fwd<T, s_t> &a);

		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> ceil(const Fwd<T, s_t> &a);
		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> floor(const Fwd<T, s_t> &a);

		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> fmax(const Fwd<T, s_t> &a, const Fwd<T, s_t> &b);
		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> fmax(T v, const Fwd<T, s_t> &a);
		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> fmax(const Fwd<T, s_t> &a, T v);

		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> fmin(const Fwd<T, s_t> &a, const Fwd<T, s_t> &b);
		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> fmin(T v, const Fwd<T, s_t> &a);
		template<typename T, template <class T2> class s_t> inline friend Fwd<T, s_t> fmin(const Fwd<T, s_t> &a, T v);

	protected:
		Fwd(const real_t val, bool dummy) : storage_t<real_t>(), _val(val) { }

		real_t _val;
	};

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> exp(const Fwd<real_t, storage_t> &a)
	{
		Fwd<real_t, storage_t> res(std::exp(a._val), false);
		for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
			res._grad[i] = a._grad[i] * res._val;
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> log(const Fwd<real_t, storage_t> &a)
	{
//		using std::copysign;

		Fwd<real_t, storage_t> res(std::log(a._val), false);
		if (sfad_likely(a._val > real_t(0)))
		{
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = a._grad[i] / a._val;
		}
		else if (a._val == real_t(0))
		{
			const real_t inf = std::numeric_limits<real_t>::infinity();
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = copysign(inf, -a._grad[i]);
		}
		else
		{
			const real_t nAn = std::numeric_limits<real_t>::quiet_NaN();
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = nAn;
		}

		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> log10(const Fwd<real_t, storage_t> &a)
	{
//		using std::copysign;

		Fwd<real_t, storage_t> res(std::log10(a._val), false);
		if (sfad_likely(a._val > real_t(0)))
		{
			const real_t tmp = std::log(real_t(10)) * a._val;
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = a._grad[i] / tmp;
		}
		else if (a._val == real_t(0))
		{
			const real_t inf = std::numeric_limits<real_t>::infinity();
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = copysign(inf, -a._grad[i]);
		}
		else
		{
			const real_t nAn = std::numeric_limits<real_t>::quiet_NaN();
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = nAn;
		}

		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> sqrt(const Fwd<real_t, storage_t> &a)
	{
//		using std::copysign;

		Fwd<real_t, storage_t> res(std::sqrt(a._val), false);
		if (sfad_likely(a._val > real_t(0)))
		{
			const real_t tmp = real_t(2) * res._val;
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = a._grad[i] / tmp;
		}
		else if (a._val == real_t(0))
		{
			const real_t inf = std::numeric_limits<real_t>::infinity();
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = copysign(inf, a._grad[i]);
		}
		else
		{
			const real_t nAn = std::numeric_limits<real_t>::quiet_NaN();
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = nAn;
		}

		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> sqr(const Fwd<real_t, storage_t> &a)
	{
		Fwd<real_t, storage_t> res(a._val * a._val, false);
		const real_t tmp = real_t(2) * a._val;
		for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
			res._grad[i] = tmp * a._grad[i];
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> sin(const Fwd<real_t, storage_t> &a)
	{
		Fwd<real_t, storage_t> res(std::sin(a._val), false);
		const real_t tmp = std::cos(a._val);
		for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
			res._grad[i] = a._grad[i] * tmp;
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> cos(const Fwd<real_t, storage_t> &a)
	{
		Fwd<real_t, storage_t> res(std::cos(a._val), false);
		const real_t tmp = -std::sin(a._val);
		for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
			res._grad[i] = a._grad[i] * tmp;
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> tan(const Fwd<real_t, storage_t> &a)
	{
		Fwd<real_t, storage_t> res(std::tan(a._val), false);

		const real_t tmpCos = std::cos(a._val);
		const real_t tmp = tmpCos * tmpCos;
		for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
			res._grad[i] = a._grad[i] / tmp;
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> asin(const Fwd<real_t, storage_t> &a)
	{
		Fwd<real_t, storage_t> res(std::asin(a._val), false);
		const real_t tmp = std::sqrt(real_t(1) - a._val * a._val);
		for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
			res._grad[i] = a._grad[i] / tmp;
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> acos(const Fwd<real_t, storage_t> &a)
	{
		Fwd<real_t, storage_t> res(std::acos(a._val), false);
		const real_t tmp = std::sqrt(real_t(1) - a._val * a._val);
		for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
			res._grad[i] = -a._grad[i] / tmp;
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> atan(const Fwd<real_t, storage_t> &a)
	{
		Fwd<real_t, storage_t> res(std::atan(a._val), false);
		const real_t tmp = real_t(1) + a._val * a._val;
		for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
			res._grad[i] = a._grad[i] / tmp;
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> pow(const Fwd<real_t, storage_t> &a, real_t v)
	{
		Fwd<real_t, storage_t> res(std::pow(a._val, v), false);
		const real_t tmp = v * std::pow(a._val, v - real_t(1));
		for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
			res._grad[i] = a._grad[i] * tmp;
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> pow(real_t v, const Fwd<real_t, storage_t> &a)
	{
		Fwd<real_t, storage_t> res(std::pow(v, a._val), false);
		const real_t tmp = res._val * std::log(v);
		for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
			res._grad[i] = a._grad[i] * tmp;
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> pow(const Fwd<real_t, storage_t> &a, const Fwd<real_t, storage_t> &b)
	{
		Fwd<real_t, storage_t> res(std::pow(a._val, b._val), false);
		const real_t tmp1 = b._val * std::pow(a._val, b._val - real_t(1));
		const real_t tmp2 = res._val * std::log(a._val);
		for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
			res._grad[i] = a._grad[i] * tmp1 + b._grad[i] * tmp2;
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> sinh (const Fwd<real_t, storage_t> &a)
	{
		Fwd<real_t, storage_t> res(std::sinh(a._val), false);
		const real_t tmp = std::cosh(a._val);
		for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
			res._grad[i] = a._grad[i] * tmp;
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> cosh (const Fwd<real_t, storage_t> &a)
	{
		Fwd<real_t, storage_t> res(std::cosh(a._val), false);
		const real_t tmp = std::sinh(a._val);
		for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
			res._grad[i] = a._grad[i] * tmp;
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> tanh (const Fwd<real_t, storage_t> &a)
	{
		Fwd<real_t, storage_t> res(std::tanh(a._val), false);
/*
		const real_t tmp = real_t(1) - res._val * res._val;
		for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
			res._grad[i] = a._grad[i] * tmp;
*/
		const real_t tmp = std::cosh(a._val);
		const real_t tmp2 = tmp * tmp;
		for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
			res._grad[i] = a._grad[i] / tmp2;
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> fabs (const Fwd<real_t, storage_t> &a)
	{
		Fwd<real_t, storage_t> res(std::abs(a._val), false);
		
		if (a._val > real_t(0))
		{
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = a._grad[i];
		}
		else if (a._val < real_t(0))
		{
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = -a._grad[i];
		}
		else
		{
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
			{
				if (a._grad[i] > real_t(0))
					res._grad[i] = a._grad[i];
				else if (a._grad[i] < real_t(0))
					res._grad[i] = -a._grad[i];
				else
					res._grad[i] = a._grad[i];
			}
				
		}

		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> ceil (const Fwd<real_t, storage_t> &a)
	{
		Fwd<real_t, storage_t> res(std::ceil(a._val), false);
		const real_t tmp(0);
		for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
			res._grad[i] = tmp;
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> floor (const Fwd<real_t, storage_t> &a)
	{
		Fwd<real_t, storage_t> res(std::floor(a._val), false);
		const real_t tmp(0);
		for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
			res._grad[i] = tmp;
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> fmax (const Fwd<real_t, storage_t> &a, const Fwd<real_t, storage_t> &b)
	{
		Fwd<real_t, storage_t> res(real_t(0), false);
		const real_t diff = a._val - b._val;
		if (diff > real_t(0))
		{
			res._val = a._val;
			res.copyGradient(a._grad);
		}
		else if (diff < real_t(0))
		{
			res._val = b._val;
			res.copyGradient(b._grad);
		}
		else
		{
			res._val = b._val;
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = std::max(a._grad[i], b._grad[i]);
		}
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> fmax (real_t v, const Fwd<real_t, storage_t> &a)
	{
		Fwd<real_t, storage_t> res(real_t(0), false);
		const real_t diff = v - a._val;
		if (diff > real_t(0))
		{
			res._val = v;
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = real_t(0);
		}
		else if (diff < real_t(0))
		{
			res._val = a._val;
			res.copyGradient(a._grad);
		}
		else
		{
			res._val = a._val;
			const real_t tmp(0);
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = std::max(tmp, a._grad[i]);
		}
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> fmax (const Fwd<real_t, storage_t> &a, real_t v)
	{
		Fwd<real_t, storage_t> res(real_t(0), false);
		const real_t diff = a._val - v;
		if (diff > real_t(0))
		{
			res._val = a._val;
			res.copyGradient(a._grad);
		}
		else if (diff < real_t(0))
		{
			res._val = v;
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = real_t(0);
		}
		else
		{
			res._val = a._val;
			const real_t tmp(0);
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = std::max(tmp, a._grad[i]);
		}
		return res;
	}
	
	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> fmin (const Fwd<real_t, storage_t> &a, const Fwd<real_t, storage_t> &b)
	{
		Fwd<real_t, storage_t> res(real_t(0), false);
		const real_t diff = a._val - b._val;
		if (diff < real_t(0))
		{
			res._val = a._val;
			res.copyGradient(a._grad);
		}
		else if (diff > real_t(0))
		{
			res._val = b._val;
			res.copyGradient(b._grad);
		}
		else
		{
			res._val = b._val;
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = std::min(a._grad[i], b._grad[i]);
		}
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> fmin (real_t v, const Fwd<real_t, storage_t> &a)
	{
		Fwd<real_t, storage_t> res(real_t(0), false);
		const real_t diff = v - a._val;
		if (diff < real_t(0))
		{
			res._val = v;
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = real_t(0);
		}
		else if (diff > real_t(0))
		{
			res._val = a._val;
			res.copyGradient(a._grad);
		}
		else
		{
			res._val = a._val;
			const real_t tmp(0);
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = std::min(tmp, a._grad[i]);
		}
		return res;
	}

	template <typename real_t, template <class T> class storage_t>
	inline Fwd<real_t, storage_t> fmin (const Fwd<real_t, storage_t> &a, real_t v)
	{
		Fwd<real_t, storage_t> res(real_t(0), false);
		const real_t diff = a._val - v;
		if (diff < real_t(0))
		{
			res._val = a._val;
			res.copyGradient(a._grad);
		}
		else if (diff > real_t(0))
		{
			res._val = v;
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = real_t(0);
		}
		else
		{
			res._val = a._val;
			const real_t tmp(0);
			for (typename Fwd<real_t, storage_t>::idx_t i = 0; i < detail::globalGradSize; ++i)
				res._grad[i] = std::min(tmp, a._grad[i]);
		}
		return res;
	}

	template <typename real_t, template <class T> class storage_t> inline Fwd<real_t, storage_t> max (const Fwd<real_t, storage_t> &a, const Fwd<real_t, storage_t> &b) { return fmax(a, b); }
	template <typename real_t, template <class T> class storage_t> inline Fwd<real_t, storage_t> max (real_t v, const Fwd<real_t, storage_t> &a) { return fmax(v, a); }
	template <typename real_t, template <class T> class storage_t> inline Fwd<real_t, storage_t> max (const Fwd<real_t, storage_t> &a, real_t v) { return fmax(a, v); }
	template <typename real_t, template <class T> class storage_t> inline Fwd<real_t, storage_t> min (const Fwd<real_t, storage_t> &a, const Fwd<real_t, storage_t> &b) { return fmin(a, b); }
	template <typename real_t, template <class T> class storage_t> inline Fwd<real_t, storage_t> min (real_t v, const Fwd<real_t, storage_t> &a) { return fmin(v, a); }
	template <typename real_t, template <class T> class storage_t> inline Fwd<real_t, storage_t> min (const Fwd<real_t, storage_t> &a, real_t v) { return fmin(a, v); }

	template <typename real_t, template <class T> class storage_t> inline Fwd<real_t, storage_t> abs (const Fwd<real_t, storage_t> &a) { return fabs(a); }

	template <typename real_t, template <class T> class storage_t>
	void swap(Fwd<real_t, storage_t>& x, Fwd<real_t, storage_t>& y) SFAD_NOEXCEPT
	{
		using std::swap;
		swap(x._val, y._val);
		swap(x._grad, y._grad);
	}

}

#endif
