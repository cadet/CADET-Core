// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2016: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file 
 * Declares some OpenMP helper functions and types
 */

#ifndef LIBCADET_OMPSUPPORT_HPP_
#define LIBCADET_OMPSUPPORT_HPP_

#ifdef _OPENMP
	#include <omp.h>

	#ifdef _MSC_VER
		// MS Visual Studio only implements OpenMP 2.0
 		// Thus, it does not support loop variables that are not of type int
		typedef int ompuint_t;
	#else
		// All other (sufficiently recent) compilers support OpenMP 3.0
		typedef unsigned int ompuint_t;
	#endif
#else

	inline int omp_get_max_threads() { return 1; }
	inline int omp_get_thread_num() { return 0; }

	typedef unsigned int ompuint_t;
#endif

#endif  // LIBCADET_OMPSUPPORT_HPP_
