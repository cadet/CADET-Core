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
 * Helper functions and macros for parallelization.
 */

#ifndef LIBCADET_PARALLEL_SUPPORT_HPP_
#define LIBCADET_PARALLEL_SUPPORT_HPP_

#ifdef CADET_PARALLELIZE
	#define CADET_PARFOR_END )
	#define CADET_PARNODE_END )
#else
	#define CADET_PARFOR_END
	#define CADET_PARNODE_END
#endif

#endif  // LIBCADET_PARALLEL_SUPPORT_HPP_
