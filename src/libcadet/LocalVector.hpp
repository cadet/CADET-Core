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
 * Defines local memory versions of some data structures
 */

#ifndef LIBCADET_LOCALVECTOR_HPP_
#define LIBCADET_LOCALVECTOR_HPP_

#include "common/CompilerSpecific.hpp"
#include "SlicedVector.hpp"

#include <vector>
#include <algorithm>
#include <cstring>
#include <type_traits>

namespace cadet
{

namespace util
{

/**
 * @brief Moves a pointer by a given number of bytes
 * @param [in] p Pointer to move
 * @param [in] delta Number of bytes to move the pointer
 * @return Pointer at new location
 */
template <typename T>
inline T* advancePointer(T* p, std::ptrdiff_t delta) CADET_NOEXCEPT
{
	return reinterpret_cast<T*>(reinterpret_cast<char*>(p) + delta);
}

template <typename T>
inline T* advancePointer(void* p, std::ptrdiff_t delta) CADET_NOEXCEPT
{
	return reinterpret_cast<T*>(reinterpret_cast<char*>(p) + delta);
}

/**
 * @brief Represents a fixed size std::vector<T> in a contiguous buffer
 * @details The contents of the std::vector<T> are appended to the memory of the LocalVector<T>
 *          in a contiguous buffer. Thus, the full data is represented in a linear slab of memory.
 *          
 *          The LocalVector cannot be expanded or shrunk, but its contents can be changed.
 *          Memory is not owned by this class and must be created or destroyed elsewhere.
 * @tparam T Type of data stored in the vector
 */
template <class T>
class LocalVector
{
public:
	LocalVector() CADET_NOEXCEPT : _size(0), _data(nullptr) { }
	~LocalVector() CADET_NOEXCEPT { }

	// Default copy and move mechanisms
	LocalVector(const LocalVector<T>& cpy) = default;
	LocalVector(LocalVector<T>&& cpy) CADET_NOEXCEPT = default;

	inline LocalVector<T>& operator=(const LocalVector<T>& cpy) = default;

#ifdef COMPILER_SUPPORT_NOEXCEPT_DEFAULTED_MOVE
	inline LocalVector<T>& operator=(LocalVector<T>&& cpy) CADET_NOEXCEPT = default;
#else
	inline LocalVector<T>& operator=(LocalVector<T>&& cpy) = default;
#endif

	inline T& operator[](std::size_t idx) { return _data[idx]; }
	inline const T& operator[](std::size_t idx) const { return _data[idx]; }

	inline T* data() CADET_NOEXCEPT { return _data; }
	inline T const* data() const CADET_NOEXCEPT { return _data; }

	inline std::size_t size() const CADET_NOEXCEPT { return _size; }

	inline T& front() { return _data[0]; }
	inline const T& front() const { return _data[0]; }

	inline T& back() { return _data[_size - 1]; }
	inline const T& back() const { return _data[_size - 1]; }

	/**
	 * @brief Assigns the data of the given vector to this instance
	 * @details Copies the data into the local buffer.
	 * @param [in] v Data source
	 */
	inline void assign(const std::vector<T>& v)
	{
		std::memcpy(_data, v.data(), v.size() * sizeof(T));
	}

	/**
	 * @brief Prepares internal data structures for representing the given vector in some buffer
	 * @details Places pointers inside the buffer given by @p ptrData such that there
	 *          is enough room for copying the given vector @p v.
	 * @param [in,out] ptrData Pointer to buffer where the data is stored
	 * @param [in] v Template
	 */
	inline void fromTemplate(void* ptrData, const std::vector<T>& v)
	{
		_data = reinterpret_cast<T*>(ptrData);
		_size = v.size();
	}

	/**
	 * @brief Stores the given data in the given buffer
	 * @details Combines fromTemplate() and assign().
	 * @param [in,out] ptrData Pointer to buffer where the data is stored
	 * @param [in] v Data source
	 */
	inline void from(void* ptrData, const std::vector<T>& v)
	{
		fromTemplate(ptrData, v);
		assign(v);
	}

	/**
	 * @brief Returns the size in bytes for representing a vector<T>
	 * @param [in] v Data to represent in a linearized local buffer
	 * @return Size in bytes required for representing the given data
	 */
	static inline std::size_t requiredMemoryFor(const std::vector<T>& v) CADET_NOEXCEPT { return requiredMemoryFor(v.size()); }
	static inline std::size_t requiredMemoryFor(std::size_t n) CADET_NOEXCEPT { return sizeof(LocalVector<T>) + sizeof(T) * n; }

	/**
	 * @brief Constructs a copy of the given data in the given memory buffer
	 * @details Constructs the LocalVector<T> at the beginning of the given buffer and
	 *          appends the actual data of the vector<T>. The required size of the buffer
	 *          is given by requiredMemoryFor().
	 * 
	 * @param [in,out] ptr Buffer
	 * @param [in] v Data to copy
	 * @return LocalVector resembling a copy of the given vector
	 */
	static inline LocalVector* constructIn(void* ptr, const std::vector<T>& v)
	{
		return constructIn(ptr, v.size(), v.data());
	}

	/**
	 * @brief Constructs a copy of the given data in the given memory buffer
	 * @details Constructs the LocalVector<T> at the beginning of the given buffer and
	 *          appends the actual data of the array. The required size of the buffer
	 *          is given by requiredMemoryFor().
	 * 
	 * @param [in,out] ptr Buffer
	 * @param [in] n Number of elements in the data array @p v
	 * @param [in] v Data to copy
	 * @return LocalVector resembling a copy of the given data array
	 */
	static inline LocalVector* constructIn(void* ptr, std::size_t n, T const* v)
	{
		return constructIn(ptr, advancePointer<T>(ptr, sizeof(LocalVector<T>)), n, v);
	}

	/**
	 * @brief Constructs a copy of the given data in the given memory buffer
	 * @details Constructs the LocalVector<T> at one place in the given buffer and
	 *          puts the actual data of the vector<T> at some other place. The required
	 *          size of the buffer is given by requiredMemoryFor().
	 * 
	 * @param [in,out] ptrBase Pointer to buffer where class is to be instantiated
	 * @param [in,out] ptrData Pointer to buffer where data is placed
	 * @param [in] v Data to copy
	 * @return LocalVector resembling a copy of the given vector
	 */
	static inline LocalVector* constructIn(void* ptrBase, void* ptrData, const std::vector<T>& v)
	{
		return constructIn(ptrBase, ptrData, v.size(), v.data());
	}

	/**
	 * @brief Constructs a copy of the given data in the given memory buffer
	 * @details Constructs the LocalVector<T> at one place in the given buffer and
	 *          puts the actual data of the vector<T> at some other place. The required
	 *          size of the buffer is given by requiredMemoryFor().
	 * 
	 * @param [in,out] ptrBase Pointer to buffer where class is to be instantiated
	 * @param [in,out] ptrData Pointer to buffer where data is placed
	 * @param [in] n Number of elements in the data array @p v
	 * @param [in] v Data to copy
	 * @return LocalVector resembling a copy of the given data array
	 */
	static inline LocalVector* constructIn(void* ptrBase, void* ptrData, std::size_t n, T const* v)
	{
		static_assert(std::is_trivially_copyable<T>::value || std::is_same<T, active>::value, "memcpy requires is_trivially_copyable<T>");

		// Construct LocalVector<T>
		LocalVector<T>* const lv = reinterpret_cast<LocalVector<T>*>(ptrBase);
		new(lv) LocalVector<T>(n);

		// Copy data and assign
		lv->_data = reinterpret_cast<T*>(ptrData);
		std::memcpy(lv->_data, v, n);

		return lv;
	}

protected:
	LocalVector(std::size_t size) CADET_NOEXCEPT : _size(size), _data(nullptr) { }

	std::size_t _size;
	T* _data;
};


/**
 * @brief Represents a fixed size SlicedVector<T> in a contiguous buffer
 * @details The contents of the SlicedVector<T> are appended to the memory of the LocalSlicedVector<T>
 *          in a contiguous buffer. Thus, the full data is represented in a linear slab of memory.
 *          
 *          The LocalSlicedVector cannot be expanded or shrunk, but its contents can be changed.
 *          Memory is not owned by this class and must be created or destroyed elsewhere.
 * @tparam T Type of data stored in the sliced vector
 */
template <class T>
class LocalSlicedVector
{
public:
	typedef typename SlicedVector<T>::size_type size_type;

	LocalSlicedVector() CADET_NOEXCEPT : _indexSize(0), _values(nullptr), _index(nullptr) { }
	~LocalSlicedVector() CADET_NOEXCEPT { }

	// Default copy and move mechanisms
	LocalSlicedVector(const LocalSlicedVector<T>& cpy) = default;
	LocalSlicedVector(LocalSlicedVector<T>&& cpy) CADET_NOEXCEPT = default;

	inline LocalSlicedVector<T>& operator=(const LocalSlicedVector<T>& cpy) = default;

#ifdef COMPILER_SUPPORT_NOEXCEPT_DEFAULTED_MOVE
	inline LocalSlicedVector<T>& operator=(LocalSlicedVector<T>&& cpy) CADET_NOEXCEPT = default;
#else
	inline LocalSlicedVector<T>& operator=(LocalSlicedVector<T>&& cpy) = default;
#endif

	/**
	 * @brief Checks whether this LocalSlicedVector is empty
	 * @return @c true if it is empty, otherwise @c false
	 */
	inline bool empty() const { return _indexSize == 1; }

	/**
	 * @brief Returns the number of slices
	 * @return Number of slices
	 */
	inline size_type slices() const CADET_NOEXCEPT { return _indexSize - 1; }

	/**
	 * @brief Returns a pointer to the first element of a given slice
	 * @param [in] idxSlice Index of the slice
	 * @return Pointer to the first element of the given slice
	 */
	inline T const* operator[](size_type idxSlice) const { return _values + _index[idxSlice]; }

	/**
	 * @brief Returns a pointer to the first element of a given slice
	 * @param [in] idxSlice Index of the slice
	 * @return Pointer to the first element of the given slice
	 */
	inline T* operator[](size_type idxSlice) { return _values + _index[idxSlice]; }

	/**
	 * @brief Returns the element at the given position in the given slice
	 * @param [in] idxSlice Index of the slice
	 * @param [in] idxElem Index of the element within the slice
	 * @return Element at the given index of the given slice
	 */
	inline const T& operator()(size_type idxSlice, size_type idxElem) const { return _values[_index[idxSlice] + idxElem]; }

	/**
	 * @brief Returns the element at the given position in the given slice
	 * @param [in] idxSlice Index of the slice
	 * @param [in] idxElem Index of the element within the slice
	 * @return Element at the given index of the given slice
	 */
	inline T& operator()(size_type idxSlice, size_type idxElem) { return _values[_index[idxSlice] + idxElem]; }

	/**
	 * @brief Returns the element at the given linear index
	 * @param [in] idx Linear index of the requested element
	 * @return Element at the given linear index
	 */
	inline T& native(size_type idx) { return _values[idx]; }

	/**
	 * @brief Returns the element at the given linear index
	 * @param [in] idx Linear index of the requested element
	 * @return Element at the given linear index
	 */
	inline const T& native(size_type idx) const { return _values[idx]; }

	/**
	 * @brief Returns the element at the given position in the given slice
	 * @param [in] idxSlice Index of the slice
	 * @param [in] idxElem Index of the element within the slice
	 * @return Element at the given index of the given slice
	 */
	inline const T& at(size_type idxSlice, size_type idxElem) const { return _values[_index[idxSlice] + idxElem]; }

	/**
	 * @brief Returns the element at the given position in the given slice
	 * @param [in] idxSlice Index of the slice
	 * @param [in] idxElem Index of the element within the slice
	 * @return Element at the given index of the given slice
	 */
	inline T& at(size_type idxSlice, size_type idxElem) { return _values[_index[idxSlice] + idxElem]; }

	/**
	 * @brief Returns a pointer to the first element of a given slice
	 * @param [in] idxSlice Index of the slice
	 * @return Pointer to the first element of the given slice
	 */
	inline T const* at(size_type idxSlice) const { return _values + _index[idxSlice]; }

	/**
	 * @brief Returns a pointer to the first element of a given slice
	 * @param [in] idxSlice Index of the slice
	 * @return Pointer to the first element of the given slice
	 */
	inline T* at(size_type idxSlice) { return _values + _index[idxSlice]; }

	/**
	 * @brief Returns a pointer to the first element of the last slice
	 * @return Pointer to the first element of the last slice
	 */
	inline T const* back() const 
	{
		cadet_assert(!empty());
		return &_values[_index[_indexSize - 2]];
	}
	
	inline T* back()
	{
		cadet_assert(!empty());
		return &_values[_index[_indexSize - 2]];
	}

	/**
	 * @brief Returns the number of elements in a slice
	 * @param [in] idxSlice Index of the slice
	 * @return Number of elmenets in the slice
	 */
	inline size_type sliceSize(size_type idxSlice) const { return _index[idxSlice + 1] - _index[idxSlice]; }

	/**
	 * @brief Returns the total number of stored items
	 * @details All items in every slice are counted.
	 * @return Total number of stored items
	 */
	inline size_type size() const CADET_NOEXCEPT { return _index[_indexSize - 1]; }

	/**
	 * @brief Returns the offset of the slice with the given index in the linear array
	 * @param [in] idx Index of the slice
	 * @return Offset of the given slice in the linearized storage
	 */
	inline size_type sliceOffset(size_type idx) const { return _index[idx]; }

	/**
	 * @brief Returns a pointer to the first element of the underlying linear array
	 * @return Pointer to first element of underlying linear array
	 */
	inline T* data() CADET_NOEXCEPT
	{
		return _values;
	}

	inline T const* data() const CADET_NOEXCEPT
	{
		return _values;
	}

	/**
	 * @brief Returns a pointer to the first element of the slice index array
	 * @return Pointer to first element of slice index array
	 */
	inline size_type* indices() CADET_NOEXCEPT
	{
		return _index;
	}

	inline size_type const* indices() const CADET_NOEXCEPT
	{
		return _index;
	}

	/**
	 * @brief Assigns the data of the given SlicedVector to this instance
	 * @details Copies the data into the local buffer.
	 * @param [in] v Data source
	 */
	inline void assign(const SlicedVector<T>& v)
	{
		std::memcpy(_values, v.data(), v.size() * sizeof(T));
	}

	/**
	 * @brief Assigns the slice indices of the given SlicedVector to this instance
	 * @details Copies the indices into the local buffer.
	 * @param [in] v Data source
	 */
	inline void assignIndices(const SlicedVector<T>& v)
	{
		std::memcpy(_index, v.indices(), (v.slices() + 1) * sizeof(typename SlicedVector<T>::size_type));
	}

	/**
	 * @brief Prepares internal data structures for representing the given SlicedVector in some buffer
	 * @details Places pointers inside the buffer given by @p ptrData such that there
	 *          is enough room for copying the given SlicedVector @p v.
	 * @param [in,out] ptrData Pointer to buffer where the data is stored
	 * @param [in] v Template
	 */
	inline void fromTemplate(void* ptrData, const SlicedVector<T>& v)
	{
		_values = reinterpret_cast<T*>(ptrData);
		_index = advancePointer<size_type>(_values, sizeof(T) * v.size());
		_indexSize = v.slices() + 1;
		assignIndices(v);
	}

	/**
	 * @brief Stores the given data in the given buffer
	 * @details Combines fromTemplate() and assign().
	 * @param [in,out] ptrData Pointer to buffer where the data is stored
	 * @param [in] v Data source
	 */
	inline void from(void* ptrData, const SlicedVector<T>& v)
	{
		fromTemplate(ptrData, v);
		assign(v);
	}

	/**
	 * @brief Returns the size in bytes for representing a SlicedVector<T>
	 * @param [in] v Data to represent in a linearized local buffer
	 * @return Size in bytes required for representing the given data
	 */
	static inline std::size_t requiredMemoryFor(const SlicedVector<T>& v) CADET_NOEXCEPT { return sizeof(LocalSlicedVector<T>) + (v.slices() + 1) * sizeof(size_type) + v.size() * sizeof(T); }

	/**
	 * @brief Constructs a copy of the given data in the given memory buffer
	 * @details Constructs the LocalSlicedVector<T> at the beginning of the given buffer and
	 *          appends the actual data of the SlicedVector<T>. The required size of the buffer
	 *          is given by requiredMemoryFor().
	 * 
	 * @param [in,out] ptr Buffer
	 * @param [in] v Data to copy
	 * @return LocalSlicedVector resembling a copy of the given SlicedVector
	 */
	static inline LocalSlicedVector* constructIn(void* ptr, const SlicedVector<T>& v)
	{
		return constructIn(ptr, advancePointer<T>(ptr, sizeof(LocalSlicedVector<T>)), v);
	}

	/**
	 * @brief Constructs a copy of the given data in the given memory buffer
	 * @details Constructs the LocalSlicedVector<T> at one place in the given buffer and
	 *          puts the actual data of the SlicedVector<T> at another. The required size
	 *          of the buffer is given by requiredMemoryFor().
	 * 
	 * @param [in,out] ptrBase Pointer to buffer where class is to be instantiated
	 * @param [in,out] ptrData Pointer to buffer where data is placed
	 * @param [in] v Data to copy
	 * @return LocalSlicedVector resembling a copy of the given SlicedVector
	 */
	static inline LocalSlicedVector* constructIn(void* ptrBase, void* ptrData, const SlicedVector<T>& v)
	{
		static_assert(std::is_trivially_copyable<T>::value || std::is_same<T, active>::value, "memcpy requires is_trivially_copyable<T>");

		// Construct LocalSlicedVector<T>
		LocalSlicedVector<T>* const lv = reinterpret_cast<LocalSlicedVector<T>*>(ptrBase);
		new(lv) LocalSlicedVector<T>(v.slices() + 1);

		// Copy data and assign
		lv->_values = reinterpret_cast<T*>(ptrData);
		std::memcpy(lv->_values, v.data(), v.size());
		lv->_index = advancePointer<size_type>(lv->_values, sizeof(T) * v.size());
		std::memcpy(lv->_index, v.indices(), v.slices() + 1);

		return lv;
	}

protected:
	LocalSlicedVector(std::size_t nSlices) CADET_NOEXCEPT : _indexSize(nSlices), _values(nullptr), _index(nullptr) { }

	std::size_t _indexSize; //!< Size of the _index array
	T* _values; //!< Data
	size_type* _index; //!< Index array with start and end indices of the slices in _values
};


/**
 * @brief Returns the active data of the container
 * @param [in] v Container
 * @return Active data of the container
 */
inline active* dataOfLocalVersion(active& v) CADET_NOEXCEPT
{
	return &v;
}

inline active* dataOfLocalVersion(LocalVector<active>& v) CADET_NOEXCEPT
{
	return v.data();
}

inline active* dataOfLocalVersion(LocalSlicedVector<active>& v) CADET_NOEXCEPT
{
	return v.data();
}

/**
 * @brief Returns a pointer behind the data of the given container
 * @param [in] v Container
 * @return Pointer behind the data of the given container
 */
inline void* ptrToEndOfData(active& v) CADET_NOEXCEPT
{
	return &v;
}

inline void* ptrToEndOfData(LocalVector<active>& v) CADET_NOEXCEPT
{
	return v.data() + v.size();
}

inline void* ptrToEndOfData(LocalSlicedVector<active>& v) CADET_NOEXCEPT
{
	return v.indices() + v.slices() + 1;
}

/**
 * @brief Returns the number of bytes required for representing the given data with a local equivalent
 * @param [in] v Data to represent
 * @return Number of bytes required for representing the data with a local equivalent
 */
inline std::size_t memoryForLocalVersionOf(double v) CADET_NOEXCEPT
{
	return sizeof(double);
}

inline std::size_t memoryForLocalVersionOf(const active& v) CADET_NOEXCEPT
{
	return sizeof(active);
}

template <class T>
inline std::size_t memoryForLocalVersionOf(const std::vector<T>& v) CADET_NOEXCEPT
{
	return LocalVector<T>::requiredMemoryFor(v);
}

template <class T>
inline std::size_t memoryForLocalVersionOf(const SlicedVector<T>& v) CADET_NOEXCEPT
{
	return LocalSlicedVector<T>::requiredMemoryFor(v);
}

/**
 * @brief Returns the size of the given data in bytes
 * @param [in] v Data
 * @return Size of the data in bytes
 */
inline std::size_t memoryForDataOf(double v) CADET_NOEXCEPT
{
	return 0;
}

inline std::size_t memoryForDataOf(const active& v) CADET_NOEXCEPT
{
	return 0;
}

template <class T>
inline std::size_t memoryForDataOf(const std::vector<T>& v) CADET_NOEXCEPT
{
	return v.size() * sizeof(T);
}

template <class T>
inline std::size_t memoryForDataOf(const SlicedVector<T>& v) CADET_NOEXCEPT
{
	return (v.slices() + 1) * sizeof(typename SlicedVector<T>::size_type) + v.size() * sizeof(T);
}

template <class T>
inline std::size_t memoryForDataOf(const LocalVector<T>& v) CADET_NOEXCEPT
{
	return v.size() * sizeof(T);
}

template <class T>
inline std::size_t memoryForDataOf(const LocalSlicedVector<T>& v) CADET_NOEXCEPT
{
	return (v.slices() + 1) * sizeof(typename LocalSlicedVector<T>::size_type) + v.size() * sizeof(T);
}

/**
 * @brief Provides the type of the local version of a given data type
 */
template <typename T>
struct localVersionOf { };

template <>
struct localVersionOf<double> { typedef double type; };

template <>
struct localVersionOf<active> { typedef active type; };

template <typename T>
struct localVersionOf<std::vector<T>> { typedef LocalVector<T> type; };

template <typename T>
struct localVersionOf<SlicedVector<T>> { typedef LocalSlicedVector<T> type; };

template <typename T>
struct localVersionOf<LocalVector<T>> { typedef LocalVector<T> type; };

template <typename T>
struct localVersionOf<LocalSlicedVector<T>> { typedef LocalSlicedVector<T> type; };

} // namespace util
} // namespace cadet

#endif  // LIBCADET_LOCALVECTOR_HPP_
