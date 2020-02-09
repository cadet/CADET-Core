// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2020: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
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

namespace cadet
{

namespace model
{

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

} // namespace model
} // namespace cadet

#endif  // LIBCADET_MODELUTILS_HPP_
