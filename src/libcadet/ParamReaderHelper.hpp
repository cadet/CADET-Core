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
 * Provides helper functions for handling parameters
 */

#ifndef LIBCADET_PARAMREADERHELPER_HPP_
#define LIBCADET_PARAMREADERHELPER_HPP_

#include "cadet/ParameterId.hpp"
#include "cadet/ParameterProvider.hpp"
#include "cadet/Exceptions.hpp"

#include "common/CompilerSpecific.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace cadet
{
	/**
	 * @brief Reads a scalar, possibly vectorial, parameter into a vector
	 * @details This function automatically detects if the parameter is vectorial.
	 *          If the user always wants a vectorial result, the single scalar value is
	 *          replicated to reach the requested array size.
	 * 
	 * @param [in] dest Destination vector in which the data is saved
	 * @param [in] paramProvider Parameter provider from which is read
	 * @param [in] dataSet Name of the dataset
	 * @param [in] nExpand How often a single scalar parameter is replicated (@c 1 if scalars should remain scalar)
	 * @tparam ValType Type of the parameter, such as @c active or @c double
	 */
	template <typename ValType>
	inline void readScalarParameterOrArray(std::vector<ValType>& dest, IParameterProvider& paramProvider, const std::string& dataSet, unsigned int nExpand)
	{
		dest.clear();
		if (paramProvider.isArray(dataSet))
		{
			// Copy all values from provider to destination
			const std::vector<double> vals = paramProvider.getDoubleArray(dataSet);
			dest.resize(vals.size());
			for (unsigned int i = 0; i < vals.size(); ++i)
				dest[i] = vals[i];
		}
		else
		{
			// Only scalar value in provider, expand to given size
			const double val = paramProvider.getDouble(dataSet);
			dest.resize(nExpand);
			for (unsigned int i = 0; i < nExpand; ++i)
				dest[i] = val;
		}
	}

	
	/**
	 * @brief Reads a vector valued parameter, that may also be a matrix
	 * @details It is automatically detected if the parameter is matrixoid or vectorial.
	 *          If the user always wants a matrixoid result, the vectorial parameter
	 *          is replicated to reach the requested size. The ordering is
	 *          vec0rep0, vec1rep0, vec2rep0, vec0rep1, vec1rep1, vec2rep1, etc.
	 * 
	 * @param [in] dest Destination vector in which the data is saved
	 * @param [in] paramProvider Parameter provider from which is read
	 * @param [in] dataSet Name of the dataset
	 * @param [in] nExpect Expected number of elements in the vectorial parameter
	 * @param [in] nExpand How often the vectorial parameter is replicated (@c 1 if vectorials should remain vectorial)
	 * @tparam ValType Type of the parameter, such as @c active or @c double
	 */
	template <typename ValType>
	inline void readParameterMatrix(std::vector<ValType>& dest, IParameterProvider& paramProvider, const std::string& dataSet, unsigned int nExpect, unsigned int nExpand)
	{
		dest.clear();
		const std::vector<double> vals = paramProvider.getDoubleArray(dataSet);

		if (vals.size() == nExpect)
		{
			// Read expected number of values, so copy them
			dest.resize(nExpect);
			for (unsigned int i = 0; i < nExpect; ++i)
				dest[i] = vals[i];
		}
		else
		{
			// We copy the read values nExpand times
			dest.resize(vals.size() * nExpand);
			for (unsigned int j = 0; j < nExpand; ++j)
			{
				for (unsigned int i = 0; i < vals.size(); ++i)
					dest[i + j * nExpand] = vals[i];
			}
		}
	}


	/**
	 * @brief Reads a parameter that depends on components and bound states available in a state-major matrix
	 * @details The ordering is expected to be state-major, that is, for a 3 component and 2 bound state parameter
	 *          the parameters @f$p_{\text{State},\text{Comp}} @f$ are stored like this:
	 *          @f[ \begin{pmatrix} p_{11} & p_{12} & p_{13} \\ p_{21} & p_{22} & p_{23} \\ p_{31} & p_{32} & p_{33}
	 *              \Leftrightarrow \left( p_{11}, p_{12}, p_{13}, p_{21}, \dots, p_{23} \right). @f]
	 *          The rows of the matrix which contain the parameters for all components of one specific bound state are
	 *          appended to one another to form a long array.
	 * 
	 * @param [in] dest Destination vector in which the data is saved
	 * @param [in] paramProvider Parameter provider from which is read
	 * @param [in] dataSet Name of the dataset
	 * @param [nComp] Number of required components
	 * @param [nStates] Number of required bound states
	 * @tparam SliceContainer_t Type of the slice container, such as @c cadet::util::SlicedVector
	 * @tparam ValType Type of the parameter, such as @c active or @c double
	 */
	template <typename SliceContainer_t, typename ValType>
	inline void readBoundStateDependentParameter(SliceContainer_t& dest, IParameterProvider& paramProvider, const std::string& dataSet, unsigned int nComp, unsigned int nStates)
	{
		dest.clear();
		const std::vector<double> vals = paramProvider.getDoubleArray(dataSet);

		if (vals.size() < nComp * nStates)
		{
			throw InvalidParameterException("Not enough elements in dataset " + dataSet + " (expected " + std::to_string(nComp * nStates)
					+ " but got only " + std::to_string(vals.size()) + ")");
		}

		for (unsigned int i = 0; i < nStates; ++i)
		{
			dest.pushBackSlice(nComp);
			ValType* const destCur = dest[i];
			for (unsigned int j = 0; j < nComp; ++j)
				destCur[j] = vals[i * nComp + j];
		}
	}


	/**
	 * @brief Reads a parameter that depends on components and bound states available in a component-major matrix
	 * @details The ordering is expected to be component-major, that is, for a 3 component and 2 bound state parameter
	 *          the parameters @f$p_{\text{State},\text{Comp}} @f$ are stored like this:
	 *          @f[ \begin{pmatrix} p_{11} & p_{12} & p_{13} \\ p_{21} & p_{22} & p_{23} \\ p_{31} & p_{32} & p_{33}
	 *              \Leftrightarrow \left( p_{11}, p_{21}, p_{12}, p_{22}, \dots, p_{23} \right). @f]
	 *          The rows of the matrix which contain the parameters for all components of one specific bound state are
	 *          appended to one another to form a long array.
	 * 
	 * @param [in] dest Destination vector in which the data is saved
	 * @param [in] paramProvider Parameter provider from which is read
	 * @param [in] dataSet Name of the dataset
	 * @param [nComp] Number of required components
	 * @param [nStates] Array with number of bound states for each component
	 * @tparam SliceContainer_t Type of the slice container, such as @c cadet::util::SlicedVector
	 * @tparam ValType Type of the parameter, such as @c active or @c double
	 */
	template <typename SliceContainer_t, typename ValType>
	inline void readBoundStateDependentParameter(SliceContainer_t& dest, IParameterProvider& paramProvider, const std::string& dataSet, unsigned int nComp, unsigned int const* const nStates)
	{
		dest.clear();
		const std::vector<double> vals = paramProvider.getDoubleArray(dataSet);

		unsigned int curIdx = 0;
		for (unsigned int i = 0; i < nComp; ++i)
		{
			if (vals.size() < curIdx + nStates[i])
				throw InvalidParameterException("Not enough elements in dataset " + dataSet + " (expected at least " + std::to_string(curIdx + nStates[i])
					+ " but got only " + std::to_string(vals.size()) + ")");

			dest.pushBackSlice(nStates[i]);
			ValType* const destCur = dest[i];
			for (unsigned int j = 0; j < nStates[i]; ++j, ++curIdx)
				destCur[j] = vals[curIdx];
		}
	}


	/**
	 * @brief Reads a matrix valued parameter that depends on components and bound states available in a component-row-major 3D array
	 * @details The ordering is expected to be component-row-major, that is, for a 3 component and 2 bound state matrix valued parameter
	 *          the parameters @f$p_{\text{Comp},\text{State},\text{State}} @f$ are stored like this:
	 *          @f[ p_{i,\cdot,\cdot} = \begin{pmatrix} p_{i,1,1} & p_{i,1,2} \\ p_{i,2,1} & p_{i,2,2} \end{pmatrix} 
	 *              \Leftrightarrow \left( p_{111}, p_{112}, p_{121}, p_{122}, p_{211}, p_{212}, p_{221}, p_{222}, \dots, p_{322} \right). @f]
	 *          The parameters are stored in the same ordering as they are read.
	 * 
	 * @param [in] dest Destination vector in which the data is saved
	 * @param [in] paramProvider Parameter provider from which is read
	 * @param [in] dataSet Name of the dataset
	 * @param [nComp] Number of required components
	 * @param [nStates] Array with number of bound states for each component
	 * @tparam SliceContainer_t Type of the slice container, such as @c cadet::util::SlicedVector
	 * @tparam ValType Type of the parameter, such as @c active or @c double
	 */
	template <typename SliceContainer_t, typename ValType>
	inline void readMatrixValuedBoundStateDependentParameter(SliceContainer_t& dest, IParameterProvider& paramProvider, const std::string& dataSet, unsigned int nComp, unsigned int const* const nStates)
	{
		dest.clear();
		const std::vector<double> vals = paramProvider.getDoubleArray(dataSet);

		unsigned int curIdx = 0;
		for (unsigned int i = 0; i < nComp; ++i)
		{
			const unsigned int curSize = nStates[i] * nStates[i];
			if (vals.size() < curIdx + curSize)
				throw InvalidParameterException("Not enough elements in dataset " + dataSet + " (expected at least " + std::to_string(curIdx + curSize)
					+ " but got only " + std::to_string(vals.size()) + ")");

			dest.pushBackSlice(curSize);
			ValType* const destCur = dest[i];
			for (unsigned int j = 0; j < curSize; ++j, ++curIdx)
				destCur[j] = vals[curIdx];
		}
	}


	/**
	 * @brief Selects an array of possible section dependent vectorial parameters from a list
	 * @details The parameters may be section dependent, but do not have to. The function is given
	 *          a vector with all parameter values for all indices (vectorial parameter) and sections 
	 *          (if section dependent). The vector contains only @p nElements elements if the
	 *          parameter is not section dependent. Otherwise, it should contain as much elements
	 *          as @c nSec * nElements. The correct element array is selected according to the index of
	 *          the current section and the section dependency of the parameter.
	 * 
	 * @param [in] data Vector with all parameter values
	 * @param [in] nElements Number of elements of the vectorial parameter (without taking section dependency into account)
	 * @param [in] secIdx Index of the current section
	 * @return Pointer to the first element of the correct vectorial parameter array
	 */
	inline active const* getSectionDependentSlice(const std::vector<active>& data, unsigned int nElements, unsigned int secIdx)
	{
		// Check if the vector contains enough elements
		if (cadet_unlikely(data.size() >= nElements * (secIdx + 1)))
		{
			// Advance data pointer to correct position
			active const* out = data.data();
			return out + nElements * secIdx;
		}
		else
			return data.data();
	}


	/**
	 * @brief Selects a possibly section dependent scalar parameter from a list of parameters
	 * @details The parameter may be section dependent, but does not have to. The function is given
	 *          a vector with all parameter values. The vector contains only one element if the
	 *          parameter is not section dependent. Otherwise, it should contain as much elements
	 *          as there are sections. The correct element is selected according to the index of
	 *          the current section and the section dependency of the parameter.
	 * 
	 * @param [in] data Vector with all parameter values (for every section, or just one)
	 * @param [in] secIdx Index of the current section
	 * @return The correct scalar parameter
	 */
	inline const active& getSectionDependentScalar(const std::vector<active>& data, unsigned int secIdx)
	{
		if (cadet_unlikely(data.size() > 1))
		{
			cadet_assert(data.size() > secIdx);
			return data[secIdx];
		}
		else
		{
			cadet_assert(data.size() == 1);
			return data[0];
		}
	}


	/**
	 * @brief Registers scalar parameters that may be section dependent
	 * 
	 * @param [in] nameHash Hashed name of the parameter
	 * @param [in,out] map Map to which the parameters are added
	 * @param [in] params Vector with parameters to be registered (only 1 if not section dependent)
	 * @param [in] unitOpIdx Index of the unit operation
	 */
	inline void registerScalarSectionDependentParam(const StringHash nameHash, std::unordered_map<ParameterId, active*>& map, std::vector<active>& params, UnitOpIdx unitOpIdx)
	{
		if (params.size() > 1)
		{
			for (unsigned int i = 0; i < params.size(); ++i)
			{
				map[makeParamId(nameHash, unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, i)] = &params[i];
			}
		}
		else
		{
			map[makeParamId(nameHash, unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &params[0];
		}
	}


	/**
	 * @brief Registers scalar parameters that may be bound phase dependent
	 * 
	 * @param [in] nameHash Hashed name of the parameter
	 * @param [in,out] map Map to which the parameters are added
	 * @param [in] params Vector with parameters to be registered (only 1 if not bound state dependent)
	 * @param [in] unitOpIdx Index of the unit operation
	 */
	inline void registerScalarBoundStateDependentParam(const StringHash nameHash, std::unordered_map<ParameterId, active*>& map, std::vector<active>& params, UnitOpIdx unitOpIdx)
	{
		for (unsigned int i = 0; i < params.size(); ++i)
		{
			map[makeParamId(nameHash, unitOpIdx, CompIndep, i, ReactionIndep, SectionIndep)] = &params[i];
		}
	}


	/**
	 * @brief Registers vector parameters (component dependent)
	 * @param [in] nameHash Hashed name of the parameter
	 * @param [in,out] map Map to which the parameters are added
	 * @param [in] params Vector with parameters to be registered
	 * @param [in] unitOpIdx Index of the unit operation
	 */
	inline void registerComponentDependentParam(const StringHash nameHash, std::unordered_map<ParameterId, active*>& map, std::vector<active>& params, UnitOpIdx unitOpIdx)
	{
		for (unsigned int j = 0; j < params.size(); ++j)
		{
			map[makeParamId(nameHash, unitOpIdx, j, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &params[j];
		}
	}


	/**
	 * @brief Registers vector parameters (component dependent) that may be section dependent
	 * @details The ordering is section-major, that is comp0sec0, comp1sec0, comp2sec0, comp0sec1, comp1sec1, comp2sec1, ...
	 *          The number of components is given in @p compStride.
	 *
	 * @param [in] nameHash Hashed name of the parameter
	 * @param [in,out] map Map to which the parameters are added
	 * @param [in] params Vector with parameters to be registered
	 * @param [in] unitOpIdx Index of the unit operation
	 * @param [in] compStride Distance between two parameters with consecutive sections
	 */
	inline void registerComponentSectionDependentParam(const StringHash nameHash, std::unordered_map<ParameterId, active*>& map, std::vector<active>& params, UnitOpIdx unitOpIdx, unsigned int compStride)
	{
		if (params.size() > compStride)
		{
			for (unsigned int i = 0; i < params.size() / compStride; ++i)
			{
				for (unsigned int j = 0; j < compStride; ++j)
				{
					map[makeParamId(nameHash, unitOpIdx, j, BoundPhaseIndep, ReactionIndep, i)] = &params[i * compStride + j];
				}
			}
		}
		else
		{
			for (unsigned int j = 0; j < compStride; ++j)
			{
				map[makeParamId(nameHash, unitOpIdx, j, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &params[j];
			}
		}
	}


	/**
	 * @brief Registers vector parameters (component dependent) that are also dependent on the bound phase
	 * @details The ordering is state-major, that is comp0state0, comp1state0, comp2state0, comp0state1, comp1state1, comp2state1, ...
	 *          The number of components is given in @p compStride.
	 *
	 * @param [in] nameHash Hashed name of the parameter
	 * @param [in,out] map Map to which the parameters are added
	 * @param [in] params Vector with parameters to be registered
	 * @param [in] unitOpIdx Index of the unit operation
	 * @param [in] compStride Distance between two parameters with consecutive bound phases
	 */
	inline void registerComponentBoundStateDependentParam(const StringHash nameHash, std::unordered_map<ParameterId, active*>& map, std::vector<active>& params, UnitOpIdx unitOpIdx, unsigned int compStride)
	{
		for (unsigned int i = 0; i < params.size() / compStride; ++i)
		{
			for (unsigned int j = 0; j < compStride; ++j)
			{
				map[makeParamId(nameHash, unitOpIdx, j, i, ReactionIndep, SectionIndep)] = &params[i * compStride + j];
			}
		}
	}


	/**
	 * @brief Registers vector parameters (component dependent) that are also dependent on the bound phase
	 * @details The ordering is state-major, that is comp0state0, comp1state0, comp2state0, comp0state1, comp1state1, comp2state1, ...
	 *          It is assumed that only one bound state is present.
	 *
	 * @param [in] nameHash Hashed name of the parameter
	 * @param [in,out] map Map to which the parameters are added
	 * @param [in] params Vector with parameters to be registered
	 * @param [in] unitOpIdx Index of the unit operation
	 */
	inline void registerComponentBoundStateDependentParam(const StringHash nameHash, std::unordered_map<ParameterId, active*>& map, std::vector<active>& params, UnitOpIdx unitOpIdx)
	{
		registerComponentBoundStateDependentParam(nameHash, map, params, unitOpIdx, params.size());
	}


	/**
	 * @brief Registers parameters that depend on component and bound state stored in state-major ordering
	 * @details The ordering is state-major, that is comp0state0, comp1state0, comp2state0, comp0state1, comp1state1, comp2state1, ...
	 *
	 * @param [in] nameHash Hashed name of the parameter
	 * @param [in,out] map Map to which the parameters are added
	 * @param [in] params Sliced vector with parameters to be registered (each slice corresponds to one bound state and contains all components)
	 * @param [in] unitOpIdx Index of the unit operation
	 * @tparam SliceContainer_t Type of the slice container, such as @c cadet::util::SlicedVector
	 */
	template <typename SliceContainer_t>
	inline void registerComponentBoundStateDependentParam(const StringHash nameHash, std::unordered_map<ParameterId, active*>& map, SliceContainer_t& params, UnitOpIdx unitOpIdx)
	{
		const unsigned int numComp = params.sliceSize(0);
		const unsigned int numStates = params.slices();
		for (unsigned int state = 0; state < numStates; ++state)
		{
			active* const slice = params[state];
			for (unsigned int comp = 0; comp < numComp; ++comp)
			{
				map[makeParamId(nameHash, unitOpIdx, comp, state, ReactionIndep, SectionIndep)] = slice + comp;
			}
		}
	}


	/**
	 * @brief Registers parameters that depend on component and bound state stored in component-major ordering
	 * @details The ordering is component-major, that is comp0state0, comp0state1, comp1state0, comp1state1, comp2state0, comp2state1, ...
	 *
	 * @param [in] nameHash Hashed name of the parameter
	 * @param [in,out] map Map to which the parameters are added
	 * @param [in] params Sliced vector with parameters to be registered (each slice corresponds to one component and contains all bound states)
	 * @param [in] unitOpIdx Index of the unit operation
	 * @tparam SliceContainer_t Type of the slice container, such as @c cadet::util::SlicedVector
	 */
	template <typename SliceContainer_t>
	inline void registerComponentBoundStateDependentParamCompMajor(const StringHash nameHash, std::unordered_map<ParameterId, active*>& map, SliceContainer_t& params, UnitOpIdx unitOpIdx)
	{
		const unsigned int numComp = params.slices();
		for (unsigned int comp = 0; comp < numComp; ++comp)
		{
			const unsigned int numStates = params.sliceSize(comp);
			active* const slice = params[comp];
			for (unsigned int state = 0; state < numStates; ++state)
			{
				map[makeParamId(nameHash, unitOpIdx, comp, state, ReactionIndep, SectionIndep)] = slice + state;
			}
		}
	}

} // namespace cadet

#endif  // LIBCADET_PARAMREADERHELPER_HPP_
