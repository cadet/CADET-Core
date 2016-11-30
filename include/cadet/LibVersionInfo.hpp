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
 * Provides version info.
 */

#ifndef LIBCADET_LIBVERSIONINFO_HPP_
#define LIBCADET_LIBVERSIONINFO_HPP_

#include "cadet/LibExportImport.hpp"
#include "cadet/cadetCompilerInfo.hpp"

extern "C"
{
	/**
	 * @brief Returns the version string of the libcadet library
	 * @return Version string
	 */
	CADET_API const char* cadetGetLibraryVersion();

	/**
	 * @brief Returns the git commit hash of the source which was used to build the binaries
	 * @return Git commit hash as string
	 */
	CADET_API const char* cadetGetLibraryCommitHash();

	/**
	 * @brief Returns the git refspec of the source which was used to build the binaries
	 * @return Git refspec
	 */
	CADET_API const char* cadetGetLibraryBranchRefspec();
}

namespace cadet
{

	/**
	 * @brief Returns the version string of the libcadet library
	 * @sa cadetGetLibraryVersion()
	 * @return Version string
	 */
	CADET_API const char* getLibraryVersion() CADET_NOEXCEPT;

	/**
	 * @brief Returns the git commit hash of the source which was used to build the binaries
	 * @sa cadetGetLibraryCommitHash()
	 * @return Git commit hash as string
	 */
	CADET_API const char* getLibraryCommitHash() CADET_NOEXCEPT;

	/**
	 * @brief Returns the git refspec of the source which was used to build the binaries
	 * @sa cadetGetLibraryBranchRefspec()
	 * @return Git refspec
	 */
	CADET_API const char* getLibraryBranchRefspec() CADET_NOEXCEPT;

} // namespace cadet

#endif  // LIBCADET_LIBVERSIONINFO_HPP_
