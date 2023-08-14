// =============================================================================
//  CADET
//  
//  Copyright © 2008-2022: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file
 * Defines a discontinuous Galerkin library that implements required basic functionality and discrete operators.
 */

#ifndef LIBCADET_DGLIBRARY_HPP_
#define LIBCADET_DGLIBRARY_HPP_

 //#include "ParamIdUtil.hpp"
#include "AutoDiff.hpp"
//#include "Memory.hpp"
//#include "SimulationTypes.hpp"
//#include "linalg/BandedEigenSparseRowIterator.hpp"

//#include <unordered_map>
//#include <unordered_set>
#include <Eigen/Dense>
//#include <Eigen/Sparse>
//#include <vector>

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <math.h>
#endif


namespace cadet
{

	//class IParameterProvider;
	//class IConfigHelper;
	//struct AdJacobianParams;
	//struct SimulationTime;
	//class IModel;

	namespace model
	{

		//class IParameterParameterDependence;

		namespace parts
		{




		}
	}
}

#endif  // LIBCADET_DGLIBRARY_HPP_