# =============================================================================
#  CADET - The Chromatography Analysis and Design Toolkit
#  
#  Copyright © 2008-2018: The CADET Authors
#            Please see the AUTHORS and CONTRIBUTORS file.
#  
#  All rights reserved. This program and the accompanying materials
#  are made available under the terms of the GNU Public License v3.0 (or, at
#  your option, any later version) which accompanies this distribution, and
#  is available at http://www.gnu.org/licenses/gpl.html
# =============================================================================
 
# Require a fairly new cmake version
cmake_minimum_required (VERSION 3.6)

# Prohibit in-source build
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
	message (FATAL_ERROR "In-source build prohibited.")
endif ()

# Set module path in order to use custom CMake modules
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

find_package(Git)

# Write the current version number to variable
if (GIT_FOUND)
	if (EXISTS "${CMAKE_SOURCE_DIR}/.git")
		execute_process(COMMAND ${GIT_EXECUTABLE} describe --abbrev=0 HEAD
		                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		                OUTPUT_VARIABLE CADET_VERSION
		                OUTPUT_STRIP_TRAILING_WHITESPACE)

		if (NOT "${CADET_VERSION}" STREQUAL "")
			message(STATUS "Get version from git")

			# Remove first character ("v")
			string(LENGTH "${CADET_VERSION}" CADET_VERSION_STRLEN)
			math(EXPR CADET_VERSION_STRLEN "${CADET_VERSION_STRLEN}-1")
			string(SUBSTRING "${CADET_VERSION}" 1 ${CADET_VERSION_STRLEN}  CADET_VERSION)
		endif ()
	endif ()
endif ()

# In case of missing tags, default to versions.txt file
if ("${CADET_VERSION}" STREQUAL "")
	message(STATUS "Get version from file")
	file(STRINGS "${CMAKE_SOURCE_DIR}/version.txt" CADET_VERSION)
endif ()

message(STATUS "CADET version: ${CADET_VERSION}")

# Get current commit hash from git
if (GIT_FOUND)
	include(GetGitRevisionDescription)
	get_git_head_revision(GIT_REFSPEC GIT_SHA1)
endif()
if (NOT DEFINED GIT_SHA1)
	set(GIT_SHA1 "NO-COMMIT-HASH")
	set(GIT_REFSPEC "NO-REFSPEC")
endif()
message(STATUS "Current git HEAD: ${GIT_REFSPEC} SHA1 ${GIT_SHA1}")

# Name of the current project
project (CadetFramework CXX C)

# Set some variables with info on the CADET project
set (PACKAGE_VERSION   "${CADET_VERSION}")
set (PACKAGE_BUGREPORT "cadet@fz-juelich.de")
set (PACKAGE_NAME      "CADET")
set (PACKAGE_STRING    "CADET ${PACKAGE_VERSION}")
set (PACKAGE_TARNAME   "cadet")


# ---------------------------------------------------
#   Other configuration options
# ---------------------------------------------------

option (STANDALONE "Build binaries with as much static linking as possible" OFF)
if (STANDALONE)
	set(BLA_STATIC ON)
	set(SUNDIALS_USE_STATIC_LIBRARIES ON)
	set(HDF5_USE_STATIC_LIBRARIES ON)
endif()

option (LOGGING "Enables or disables the logging framework" ON)
if (NOT LOGGING)
	add_definitions(-DCADET_LOGGING_DISABLE)
endif ()

option (BENCHMARK_MODE "In benchmark mode fine-grained timing is enabled" OFF)
if (BENCHMARK_MODE)
	add_definitions(-DCADET_BENCHMARK_MODE)
endif ()

option (PLATFORM_TIMER "Use a platform-dependent timer implementation" OFF)
if (PLATFORM_TIMER)
	add_definitions(-DCADET_USE_PLATFORM_TIMER)
endif()

# Enable multi-threading in debug builds
option (DEBUG_THREADING "Enable multi-threading in debug builds" OFF)

# Option that allows users to build release or debug version
if (NOT CMAKE_BUILD_TYPE)
	set (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel" FORCE)
	message(STATUS "Build type: ${CMAKE_BUILD_TYPE} (default)")
	add_definitions(-DNDEBUG)
endif ()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	add_definitions(-DCADET_LOGLEVEL_MIN=Trace -DDEBUG -D_DEBUG)
else()
	add_definitions(-DCADET_LOGLEVEL_MIN=Warning)
endif()

# ---------------------------------------------------
#   Check build environment
# ---------------------------------------------------
include(WriteCompilerDetectionHeader)

find_package(TBB COMPONENTS tbb OPTIONAL_COMPONENTS tbb_preview)
set(CADET_PARALLEL_FLAG "")
if (TBB_FOUND)
	# Use tbb_preview instead of tbb if it has been found
	if (TBB_tbb_preview_FOUND)
		if (TBB_LIBRARIES_DEBUG)
			set(TBB_LIBRARIES_DEBUG ${TBB_tbb_preview_LIBRARY_DEBUG})
		endif()
		if (TBB_LIBRARIES_RELEASE)
			set(TBB_LIBRARIES_RELEASE ${TBB_tbb_preview_LIBRARY_RELEASE})
		endif()
	endif()

	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
		if (TBB_LIBRARIES_DEBUG)
			set(TBB_LIBRARIES ${TBB_LIBRARIES_DEBUG})
		else()
			# Fallback to release if debug libs not available
			set(TBB_LIBRARIES ${TBB_LIBRARIES_RELEASE})
		endif()

		if (DEBUG_THREADING)
			set(CADET_PARALLEL_FLAG "-DCADET_PARALLELIZE")
		endif()
	else()
		if (TBB_LIBRARIES_RELEASE)
			set(TBB_LIBRARIES ${TBB_LIBRARIES_RELEASE})
		else()
			# Fallback to debug if release libs not available
			set(TBB_LIBRARIES ${TBB_LIBRARIES_DEBUG})
		endif()

		set(CADET_PARALLEL_FLAG "-DCADET_PARALLELIZE")
	endif()
else()
	set(TBB_INCLUDE_DIRS "")
	set(TBB_LIBRARIES "")
	set(TBB_LIBRARIES_DEBUG "")
	set(TBB_DEFINITIONS_DEBUG "")
	set(TBB_LIBRARIES_RELEASE "")
	set(TBB_DEFINITIONS_RELEASE "")
endif()

find_package(LAPACK)

# SUNDIALS_ROOT environment variable can be used to find SUNDIALS package
find_package(SUNDIALS REQUIRED COMPONENTS sundials_idas sundials_nvecserial OPTIONAL_COMPONENTS sundials_nvecopenmp)
option (SUNDIALS_PREFER_OPENMP "Prefer OpenMP vector implementation of SUNDIALS if available (for large problems)" OFF)

# Check whether OpenMP is available in SUNDIAL'S NVECTOR module
if (SUNDIALS_sundials_nvecopenmp_LIBRARY)
	# If we are not in Debug mode, prefer OpenMP over serial version
	if ((NOT CMAKE_BUILD_TYPE STREQUAL "Debug") AND SUNDIALS_PREFER_OPENMP)
		list(REMOVE_ITEM SUNDIALS_LIBRARIES ${SUNDIALS_sundials_nvecserial_LIBRARY})
		add_definitions(-DCADET_SUNDIALS_OPENMP)
	else()
		list(REMOVE_ITEM SUNDIALS_LIBRARIES ${SUNDIALS_sundials_nvecopenmp_LIBRARY})
	endif()
endif()

# Try to find static libs
set(SUNDIALS_STATIC_LIBRARIES )
if (SUNDIALS_FOUND)
	if (NOT WIN32)
		# On Linux and Mac OS we can look for .a archive files
		foreach(_LIB IN LISTS SUNDIALS_LIBRARIES)
			get_filename_component(_LIBROOT ${_LIB} DIRECTORY)
			get_filename_component(_LIBNAME ${_LIB} NAME_WE)
			find_file(_LIBSTATIC "${_LIBNAME}.a" PATHS ${_LIBROOT} NO_DEFAULT_PATH)
			if (_LIBSTATIC)
				list(APPEND SUNDIALS_STATIC_LIBRARIES ${_LIBSTATIC})
			endif()
			unset(_LIBROOT CACHE)
			unset(_LIBNAME CACHE)
			unset(_LIBSTATIC CACHE)
			unset(_LIB CACHE)
			unset(_LIBROOT)
			unset(_LIBNAME)
			unset(_LIBSTATIC)
			unset(_LIB)
		endforeach()

		# Fall back to dynamic libs if we haven't found static libs
		if (NOT SUNDIALS_STATIC_LIBRARIES)
			set(SUNDIALS_STATIC_LIBRARIES ${SUNDIALS_LIBRARIES})  
		endif()
	else()
		# On Windows we have no way to distinguish .lib files that import dlls or represent static libs
		set(SUNDIALS_STATIC_LIBRARIES ${SUNDIALS_LIBRARIES})  
	endif()
endif()

# HDF5_ROOT environment variable can be used to find HDF5 package
find_package(HDF5)
if (HDF5_FOUND)

	# Make sure HDF5_LIBRARY_DIRS is defined
	if ((NOT DEFINED HDF5_LIBRARY_DIRS) OR (NOT HDF5_LIBRARY_DIRS) OR ("${HDF5_LIBRARY_DIRS}" STREQUAL ""))
		list(GET HDF5_LIBRARIES 0 HDF5_LIB_TEMP)
		get_filename_component(HDF5_LIBRARY_DIRS ${HDF5_LIB_TEMP} DIRECTORY)
		unset(HDF5_LIB_TEMP)
	endif()

	# Check if we need additional libraries for linking (i.e., zlib, szip)
	include(${CMAKE_ROOT}/Modules/CheckCXXSourceCompiles.cmake)
	include(${CMAKE_ROOT}/Modules/CMakePushCheckState.cmake)

	cmake_push_check_state(RESET)

	# Set libs and includes
	set(CMAKE_REQUIRED_LIBRARIES ${HDF5_LIBRARIES})
	set(CMAKE_REQUIRED_INCLUDES ${HDF5_INCLUDE_DIRS})

	CHECK_CXX_SOURCE_COMPILES("#include <hdf5.h>\nint main(int argc, char** argv){\n H5Zfilter_avail(H5Z_FILTER_SZIP);\nH5Zfilter_avail(H5Z_FILTER_DEFLATE);\nreturn 0;\n}\n" HDF5_DONT_NEED_ZLIBS)

	# Reset libs and includes
	cmake_pop_check_state()

	# Find szip and zlib libs if we need them
	if (NOT HDF5_DONT_NEED_ZLIBS)
		find_library(HDF5_SZLIB NAMES szip libszip PATHS ${HDF5_LIBRARY_DIRS})
		find_library(HDF5_ZLIB NAMES zlib libzlib PATHS ${HDF5_LIBRARY_DIRS})

		if (HDF5_SZLIB)
			list(APPEND HDF5_LIBRARIES ${HDF5_SZLIB})
		endif()
		if (HDF5_ZLIB)
			list(APPEND HDF5_LIBRARIES ${HDF5_ZLIB})
		endif()
		unset(HDF5_SZLIB)
		unset(HDF5_ZLIB)
	endif()
endif()

# Add old and new Matlab versions
set(MATLAB_ADDITIONAL_VERSIONS)

# Set Matlab_ROOT_DIR from MATLAB_ROOT environment variable
if (DEFINED ENV{MATLAB_ROOT})
	set(Matlab_ROOT_DIR $ENV{MATLAB_ROOT})
endif()
find_package(Matlab COMPONENTS MX_LIBRARY LAPACK_LIBRARY OPTIONAL_COMPONENTS TBB_LIBRARY)

# Check TBB <-> Matlab compatibility
option(MATLAB_FORCE_ENABLE_TBB "Force Matlab TBB support" OFF)
if (Matlab_FOUND)
	if (TBB_tbb_preview_FOUND)
		# Linking to tbb_preview usually works
		set(Matlab_TBB_COMPATIBLE ON)
		set(Matlab_TBB_LIBRARY ${TBB_LIBRARIES})
	else()
		set(Matlab_TBB_COMPATIBLE OFF)
		if (Matlab_TBB_LIBRARY_FOUND AND TBB_FOUND)
			if ((TBB_INTERFACE_VERSION LESS Matlab_TBB_INTERFACE_VERSION) OR (TBB_INTERFACE_VERSION EQUAL Matlab_TBB_INTERFACE_VERSION))
				# System TBB is older than or same age as Matlab's TBB
				set(Matlab_TBB_COMPATIBLE ON)
			else()
				if (Matlab_FORCE_ENABLE_TBB)
					set(Matlab_TBB_COMPATIBLE ON)
				else()
					message(WARNING "Matlab's TBB version is older (${Matlab_TBB_INTERFACE_VERSION}) than the system headers (${TBB_INTERFACE_VERSION}). TBB parallelization is disabled in MEX interface. Set -DMatlab_FORCE_ENABLE_TBB=ON to enforce TBB support.")
				endif()
			endif()
		elseif(TBB_FOUND AND Matlab_FORCE_ENABLE_TBB)
			set(Matlab_TBB_COMPATIBLE ON)
			set(Matlab_TBB_LIBRARY ${TBB_LIBRARIES})
		elseif(NOT TBB_FOUND AND NOT Matlab_TBB_LIBRARY_FOUND)
			set(Matlab_TBB_LIBRARY "")
		endif()
	endif()
endif()

#find_package(LATEX COMPONENTS PDFLATEX BIBER XINDY)

# ---------------------------------------------------
#   Which modules to build?
# ---------------------------------------------------
# CADET-CLI module
if ((IS_DIRECTORY "${PROJECT_SOURCE_DIR}/src/cadet-cli") AND HDF5_FOUND)
	option (BUILD_CADET_CLI  "Build the CADET-CLI binary" ON)
else ()
	set (BUILD_CADET_CLI OFF)
endif ()

# CADET-MI module
if ((IS_DIRECTORY "${PROJECT_SOURCE_DIR}/src/cadet-mex") AND Matlab_FOUND)
	option (BUILD_CADET_MEX "Build the Matlab interface" ON)
else ()
	set (BUILD_CADET_MEX OFF)
endif ()

# TOOLS module
if ((IS_DIRECTORY "${PROJECT_SOURCE_DIR}/src/tools") AND HDF5_FOUND)
	option (BUILD_TOOLS  "Build tools" ON)
else ()
	set (BUILD_TOOLS OFF)
endif ()

# TEST module
if (IS_DIRECTORY "${PROJECT_SOURCE_DIR}/test")
	option (BUILD_TESTS "Build tests" ON)
else ()
	set (BUILD_TESTS OFF)
endif ()
# ---------------------------------------------------



# ---------------------------------------------------
#   Add selected modules to the build system and add the targets to the list of all targets
# ---------------------------------------------------

include("${CMAKE_SOURCE_DIR}/cmake/TargetFunctions.cmake")

# Shared components
add_subdirectory (src/libcadet)
# Add LIBCADET targets to the list of all targets
set (ALL_TARGETS ${LIB_TARGETS})

# CADET-CLI module
if (BUILD_CADET_CLI)
	add_subdirectory (src/cadet-cli)
	# Add CADET-CLI target
	list (APPEND ALL_TARGETS ${CADET_CLI_TARGETS})
endif ()

# CADET-MEX module
if (BUILD_CADET_MEX)
	add_subdirectory (src/cadet-mex)
	# Add CADET-MEX target
endif ()

# CADET-MEX module
if (BUILD_TOOLS)
	add_subdirectory (src/tools)
	# Add CADET-TOOLS target
	list (APPEND ALL_TARGETS ${TOOLS_TARGETS})
endif ()

# TESTS module
if (BUILD_TESTS)
	add_subdirectory (test)
endif ()
# ---------------------------------------------------



# ---------------------------------------------------
#   Set properties, definitions, install target etc.
# ---------------------------------------------------

foreach(_TARGET IN LISTS ALL_TARGETS)
	cadet_apply_compile_options(${_TARGET})
endforeach()
# ---------------------------------------------------

message("")
message("------------------------------- Summary -------------------------------")
message("C++ compiler name: ${CMAKE_CXX_COMPILER_ID} at ${CMAKE_CXX_COMPILER}")
message("Build type: ${CMAKE_BUILD_TYPE}")
message("Build shared libs: ${BUILD_SHARED_LIBS}")
message("Build static libs: ${BUILD_STATIC_LIBS}")
message("Source dir: ${CMAKE_SOURCE_DIR}")
message("Binary dir: ${CMAKE_BINARY_DIR}")
message("Install dir: ${CMAKE_INSTALL_PREFIX}")
message("C Flags: ${CMAKE_C_FLAGS}")
message("C++ Flags: ${CMAKE_CXX_FLAGS}")
message("------------------------------- Modules -------------------------------")
message("CADET-CLI: ${BUILD_CADET_CLI}")
message("CADET-MEX: ${BUILD_CADET_MEX}")
message("TOOLS: ${BUILD_TOOLS}")
message("Tests: ${BUILD_TESTS}")
message("------------------------------- Options -------------------------------")
message("Logging: ${LOGGING}")
message("Benchmark mode: ${BENCHMARK_MODE}")
message("Platform-dependent timer: ${PLATFORM_TIMER}")
message("Standalone mode: ${STANDALONE}")
message("AD library: ${ADLIB}")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	message("Check analytic Jacobian: ${CHECK_ANALYTIC_JACOBIAN}")
endif()
message("----------------------------- Dependencies ----------------------------")

message("Found LAPACK: ${LAPACK_FOUND}")
if (LAPACK_FOUND)
	message("  Linker flags ${LAPACK_LINKER_FLAGS}")
	message("  Libs ${LAPACK_LIBRARIES}")
endif()

message("Found TBB: ${TBB_FOUND}")
if (TBB_FOUND)
	message("  Version ${TBB_VERSION} (Interface ${TBB_INTERFACE_VERSION})")
	message("  Include ${TBB_INCLUDE_DIRS}")
	message("  Libs ${TBB_LIBRARIES}")
	message("  Release Libs ${TBB_LIBRARIES_RELEASE}")
	message("  Def release ${TBB_DEFINITIONS_RELEASE}")
	if ((CMAKE_BUILD_TYPE STREQUAL "Debug") AND DEBUG_THREADING)
		message("  Debug Libs ${TBB_LIBRARIES_DEBUG}")
		message("  Def debug ${TBB_DEFINITIONS_DEBUG}")
	endif()
endif()

message("Found SUNDIALS: ${SUNDIALS_FOUND}")
if (SUNDIALS_FOUND)
	message("  Version ${SUNDIALS_VERSION}")
	message("  Includes ${SUNDIALS_INCLUDE_DIRS}")
	message("  Libs ${SUNDIALS_LIBRARIES}")
	message("  Static libs ${SUNDIALS_STATIC_LIBRARIES}")
endif()

message("Found HDF5: ${HDF5_FOUND}")
if (HDF5_FOUND)
	message("  Version ${HDF5_VERSION}")
	message("  Includes ${HDF5_INCLUDE_DIRS}")
	message("  Libs ${HDF5_LIBRARIES}")
endif()

message("Found Matlab: ${Matlab_FOUND}")
if (Matlab_FOUND)
	message("  Version ${Matlab_VERSION_STRING}")
	message("  In ${Matlab_ROOT_DIR}")
	message("  Includes ${Matlab_INCLUDE_DIRS}")
	message("  Libs ${Matlab_LIBRARIES}")
	message("  MKL ${Matlab_LAPACK_LIBRARY} ${Matlab_BLAS_LIBRARY}")
	message("  LAPACK functions have trailing underscore: ${Matlab_LAPACK_TRAILING_UNDERSCORE}")
	message("  Declares mxCreateUninitNumericMatrix: ${Matlab_DECLARES_CREATEUNINITNUMERICMATRIX}")
	message("  Declares mxCreateUninitNumericArray: ${Matlab_DECLARES_CREATEUNINITNUMERICARRAY}")
	message("  TBB: ${Matlab_TBB_LIBRARY} (interface version ${Matlab_TBB_INTERFACE_VERSION})")
	message("  Parallelization (TBB compatibility): ${Matlab_TBB_COMPATIBLE}")
endif()
message("-----------------------------------------------------------------------")
message("")
