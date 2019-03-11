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
 * Provides helper functions for parameter multiplexing
 */

#ifndef LIBCADET_PARAMMULTIPLEXING_HPP_
#define LIBCADET_PARAMMULTIPLEXING_HPP_

#include "cadet/ParameterId.hpp"
#include "AutoDiff.hpp"

#include <vector>
#include <algorithm>

namespace cadet
{

namespace model
{

	enum class MultiplexMode : int
	{
		Independent,
		Radial,
		RadialSection,
		Component,
		ComponentRadial,
		ComponentRadialSection,
		ComponentSection,
		Section,
		Axial,
		AxialRadial,
		Type,
		ComponentType,
		ComponentSectionType
	};

} // namespace model

} // namespace cadet

#endif  // LIBCADET_PARAMMULTIPLEXING_HPP_
