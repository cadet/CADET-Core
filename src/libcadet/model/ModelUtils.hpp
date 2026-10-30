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
 * @file 
 * Defines helper functions for models.
 */

#ifndef LIBCADET_MODELUTILS_HPP_
#define LIBCADET_MODELUTILS_HPP_

#include "cadet/ParameterId.hpp"
#include "ParamIdUtil.hpp"
#include "AutoDiff.hpp"
#include "SensParamUtil.hpp"

#include <cstddef>
#include <unordered_set>
#include <unordered_map>

namespace cadet
{

namespace model
{

/**
 * @brief Bulk side operator of the film transfer term
 * @details The film transfer flux of the particle at discrete bulk point @c j enters the bulk residual
 *          at the discrete points @c j + offset(j)[k] with the weights @c weight(j, comp)[k], for
 *          @c k = 0 ... stencilSize - 1.
 *
 *          A stencil of one applies whenever the mass matrix of the bulk balance and the weight of the
 *          film transfer source are proportional. That is the case for every discretization and geometry
 *          except a DGSEM discretized column that represents its flow through area by the nodal
 *          interpolant, i.e. one with a prescribed cross section or an axially varying porosity, and
 *          whose porosity or film transfer coefficient varies within an element. The single weight is
 *          then the familiar prefactor
 *          \f$ \left(1 - \varepsilon^{\mathrm{b}}\right) / \varepsilon^{\mathrm{b}} \f$ and the single
 *          offset is zero. The operator is assembled by the bulk transport operator, see
 *          parts::IConvectionDispersionOperatorBase1D::assembleFilmTransferOperator().
 */
struct bulkFilmTransfer
{
	active const* weights = nullptr; //!< Weights, laid out as described by the strides below
	int const* offsets = nullptr; //!< Signed offsets in discrete points, stencilSize per discrete point; nullptr for the node local operator, whose only offset is zero
	int stencilSize = 1; //!< Number of discrete points that the flux of one discrete point is spread over
	int strideWeightPoint = 1; //!< Stride between the weights of two discrete points
	int strideWeightComp = 0; //!< Stride between the weights of two components, zero if the weights do not depend on the component
	int strideBulkPoint = 0; //!< Stride between two discrete bulk points in the global state vector

	//! @brief Returns whether the flux of a discrete point enters the bulk residual at that point only
	inline bool nodeLocal() const { return stencilSize == 1; }

	//! @brief Returns the stencilSize weights of the flux at discrete point @p point and component @p comp
	inline active const* weight(int point, int comp) const
	{
		return weights + static_cast<std::size_t>(point) * strideWeightPoint + static_cast<std::size_t>(comp) * strideWeightComp;
	}

	//! @brief Returns the stencilSize offsets of the flux at discrete point @p point, which are all zero if offsets is nullptr
	inline int const* offset(int point) const { return offsets + static_cast<std::size_t>(point) * stencilSize; }
};

/**
 * @brief Returns whether the associated model has multiple bound states
 * @details Components with multiple bound states have several entries in the solid phase.
 * @param [in] nBound Array with number of bound states for each component
 * @param [in] nComp Number of components
 * @return @c true if multiple bound states are present, otherwise @c false
 */
inline bool hasMultipleBoundStates(unsigned int const* const nBound, unsigned int nComp)
{
	for (unsigned int i = 0; i < nComp; ++i)
	{
		if (nBound[i] > 1)
			return true;
	}
	return false;
}

/**
 * @brief Returns whether the associated model has non-binding components
 * @details Non-binding components do not have a solid phase representation.
 * @param [in] nBound Array with number of bound states for each component
 * @param [in] nComp Number of components
 * @return @c true if non-binding components are present, otherwise @c false
 */
inline bool hasNonBindingComponents(unsigned int const* const nBound, unsigned int nComp)
{
	for (unsigned int i = 0; i < nComp; ++i)
	{
		if (nBound[i] == 0)
			return true;
	}
	return false;
}

/**
 * @brief Returns the number of bound states
 * @param [in] nBound Array with number of bound states for each component
 * @param [in] nComp Number of components
 * @return Number of bound states
 */
inline unsigned int numBoundStates(unsigned int const* const nBound, unsigned int nComp)
{
	unsigned int totalBound = 0;
	for (unsigned int i = 0; i < nComp; ++i)
	{
		totalBound += nBound[i];
	}
	return totalBound;
}

/**
 * @brief Returns the number of channel states for each component
 * @param [in] nChannel Number of channels
 * @param [in] nComp Number of components
 * @return Number of channel states of each component
 */
inline unsigned int numExchangeStates(unsigned int nChannel, unsigned int nComp)
{
	return nChannel * nComp;
}

/**
 * @brief Returns the number of binding components (i.e., components that have at least @c 1 bound state)
 * @param [in] nBound Array with number of bound states for each component
 * @param [in] nComp Number of components
 * @return Number of binding components
 */
inline unsigned int numBindingComponents(unsigned int const* const nBound, unsigned int nComp)
{
	unsigned int totalBound = 0;
	for (unsigned int i = 0; i < nComp; ++i)
	{
		if (nBound[i] > 0)
			++totalBound;
	}
	return totalBound;
}

/**
 * @brief Returns the number of bound states of the first component with non-zero bound states
 * @param [in] nBound Array with number of bound states for each component
 * @param [in] nComp Number of components
 * @return Number of first non-zero bound states
 */
inline unsigned int firstNonEmptyBoundStates(unsigned int const* const nBound, unsigned int nComp)
{
	unsigned int nStates = 0;
	for (unsigned int i = 0; i < nComp; ++i)
	{
		if (nBound[i] == 0)
			continue;

		nStates = nBound[i];
		break;
	}
	return nStates;
}




template <typename T>
void getAllParameterValues(std::unordered_map<ParameterId, double>& data, const std::vector<T*>& items, bool singleItem)
{
	if (items.empty())
		return;

	if (singleItem)
	{
		if (!items[0])
			return;

		const std::unordered_map<ParameterId, double> localData = items[0]->getAllParameterValues();
		for (const std::pair<const ParameterId, double>& val : localData)
			data[val.first] = val.second;
	}
	else
	{
		for (T const* bm : items)
		{
            if (!bm)
                continue;

            const std::unordered_map<ParameterId, double> localData = bm->getAllParameterValues();
			for (const std::pair<const ParameterId, double>& val : localData)
				data[val.first] = val.second;
		}
	}
}

template <typename T>
bool getParameterDouble(const ParameterId& pId, const std::vector<T*>& items, bool singleItem, double& out)
{
	// Check binding model parameters
	if (items.empty())
		return false;

	if (singleItem)
	{
		if (!items[0])
			return false;

		active const* const val = items[0]->getParameter(pId);
		if (val)
		{
			out = static_cast<double>(*val);
			return true;
		}
	}
	else
	{
		for (T* bm : items)
		{
            if (!bm)
                continue;

            active const* const val = bm->getParameter(pId);
			if (val)
			{
				out = static_cast<double>(*val);
				return true;
			}
		}
	}

	// Not found
	return false;
}

template <typename T>
bool hasParameter(const ParameterId& pId, const std::vector<T*>& items, bool singleItem)
{
	if (items.empty())
		return false;

	if (singleItem)
	{
		if (items[0] && items[0]->hasParameter(pId))
			return true;
	}
	else
	{
		for (T const* bm : items)
		{
			if (bm && bm->hasParameter(pId))
				return true;
		}
	}

	return false;
}

template <typename T, typename param_t>
bool setParameter(const ParameterId& pId, param_t value, const std::vector<T*>& items, bool singleItem)
{
	if (items.empty())
		return false;

	if (singleItem)
	{
		if (items[0] && items[0]->setParameter(pId, value))
			return true;
	}
	else
	{
		for (T* bm : items)
		{
			if (bm && bm->setParameter(pId, value))
				return true;
		}
	}

	return false;
}

template <typename T>
bool setSensitiveParameterValue(const ParameterId& pId, double value, const std::unordered_set<active*>& sensParams, const std::vector<T*>& items, bool singleItem)
{
	if (items.empty())
		return false;

	if (singleItem)
	{
		if (!items[0])
			return false;

		active* const val = items[0]->getParameter(pId);
		if (val && contains(sensParams, val))
		{
			val->setValue(value);
			return true;
		}
	}
	else
	{
		for (T* bm : items)
		{
            if (!bm)
                continue;

			active* const val = bm->getParameter(pId);
			if (val && contains(sensParams, val))
			{
				val->setValue(value);
				return true;
			}
		}
	}

	return false;
}

template <typename T>
bool setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue, std::unordered_set<active*>& sensParams, const std::vector<T*>& items, bool singleItem)
{
	if (items.empty())
		return false;

	if (singleItem)
	{
		if (!items[0])
			return false;

		active* const paramBinding = items[0]->getParameter(pId);
		if (paramBinding)
		{
			// Register parameter and set AD seed / direction
			sensParams.insert(paramBinding);
			paramBinding->setADValue(adDirection, adValue);
			return true;
		}
	}
	else
	{
		for (T* bm : items)
		{
			if (!bm)
				continue;

			active* const paramBinding = bm->getParameter(pId);
			if (paramBinding)
			{
				// Register parameter and set AD seed / direction
				sensParams.insert(paramBinding);
				paramBinding->setADValue(adDirection, adValue);
				return true;
			}
		}
	}

	return false;
}


} // namespace model
} // namespace cadet

#endif  // LIBCADET_MODELUTILS_HPP_
