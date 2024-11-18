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
 * Defines preprocessor macros for importing and exporting symbols from and to dynamic libraries.
 */

#ifndef LIBCADET_LIBEXPORT_HPP_
#define LIBCADET_LIBEXPORT_HPP_


// Export and import classes when using MS Visual Studio compiler
#ifndef CADET_API
	#ifdef _MSC_VER
		#if defined(libcadet_shared_EXPORTS) || defined(libcadet_static_EXPORTS) || defined(libcadet_EXPORTS)
			#define CADET_API _declspec(dllexport)
		#else
			#define CADET_API _declspec(dllimport)
		#endif
	#else
		#define CADET_API __attribute__((visibility("default")))
	#endif
#endif

#endif  // LIBCADET_LIBEXPORT_HPP_
