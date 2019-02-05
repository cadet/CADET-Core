# =============================================================================
#  CADET - The Chromatography Analysis and Design Toolkit
#  
#  Copyright © 2008-2019: The CADET Authors
#            Please see the AUTHORS and CONTRIBUTORS file.
#  
#  All rights reserved. This program and the accompanying materials
#  are made available under the terms of the GNU Public License v3.0 (or, at
#  your option, any later version) which accompanies this distribution, and
#  is available at http://www.gnu.org/licenses/gpl.html
# =============================================================================

# Find UMFPACK, the direct solver for sparse linear systems.
#
# To provide the module with a hint about where to find your UMFPACK installation,
# you can set the environment variable UMFPACK_ROOT. The FindUMFPACK module will
# then look in this path when searching for UMFPACK paths and libraries.
#
# This module will define the following variables:
#  UMFPACK_FOUND - true if UMFPACK was UMFPACK on the system
#  UMFPACK_INCLUDE_DIRS - Location of the UMFPACK includes
#  UMFPACK_LIBRARIES - Required libraries for all requested components
#  UMFPACK_VERSION_MAJOR - Major version
#  UMFPACK_VERSION_MINOR - Minro version
#  UMFPACK_VERSION_PATCH - Patch level
#  UMFPACK_VERSION - Full version string


include(FindPackageHandleStandardArgs)

# find the UMFPACK include directories
find_path( UMFPACK_INCLUDE_DIRS umfpack.h
    PATHS
        ${UMFPACK_ROOT}
        ${UMFPACK_ROOT_DIR}
    ENV
        UMFPACK_ROOT
        UMFPACK_ROOT_DIR
    PATH_SUFFIXES
        include
        include/UMFPACK
        include/suitesparse
        suitesparse
        suitesparse/include
)

if(UMFPACK_INCLUDE_DIRS)
    # extract version
    file(READ "${UMFPACK_INCLUDE_DIRS}/umfpack.h" _UMFPACK_VERSION_FILE)

    string(REGEX REPLACE ".*#define UMFPACK_MAIN_VERSION[ \t]*([0-9]+)[ \t\r\n]*.*" "\\1" UMFPACK_VERSION_MAJOR "${_UMFPACK_VERSION_FILE}")
    string(REGEX REPLACE ".*#define UMFPACK_SUB_VERSION[ \t]*([0-9]+)[ \t\r\n]*.*" "\\1" UMFPACK_VERSION_MINOR "${_UMFPACK_VERSION_FILE}")
    string(REGEX REPLACE ".*#define UMFPACK_SUBSUB_VERSION[ \t]*([0-9]+)[ \t\r\n]*.*" "\\1" UMFPACK_VERSION_PATCH "${_UMFPACK_VERSION_FILE}")
    set(UMFPACK_VERSION "${UMFPACK_VERSION_MAJOR}.${UMFPACK_VERSION_MINOR}.${UMFPACK_VERSION_PATCH}")

endif()


# find the UMFPACK libraries
find_library(UMFPACK_LIBRARIES
    NAMES umfpack
        umfpack.4
        umfpack.5
    PATHS
        ${UMFPACK_INCLUDE_DIR}
        ${UMFPACK_INCLUDE_DIR}/..
        ${UMFPACK_ROOT}
        ${UMFPACK_ROOT_DIR}
    ENV
        UMFPACK_ROOT
        UMFPACK_ROOT_DIR
    PATH_SUFFIXES 
        lib
        lib64
        Lib
        Lib64
)

# find UMFPACK dependency libraries
find_library(CHOLMOD_LIBRARIES
    NAMES cholmod
    PATHS
        ${UMFPACK_INCLUDE_DIR}
        ${UMFPACK_INCLUDE_DIR}/..
        ${UMFPACK_ROOT}
        ${UMFPACK_ROOT_DIR}
    ENV
        UMFPACK_ROOT
        UMFPACK_ROOT_DIR
    PATH_SUFFIXES 
        lib
        lib64
        Lib
        Lib64
)

find_library(COLAMD_LIBRARIES
    NAMES colamd
    PATHS
        ${UMFPACK_INCLUDE_DIR}
        ${UMFPACK_INCLUDE_DIR}/..
        ${UMFPACK_ROOT}
        ${UMFPACK_ROOT_DIR}
    ENV
        UMFPACK_ROOT
        UMFPACK_ROOT_DIR
    PATH_SUFFIXES 
        lib
        lib64
        Lib
        Lib64
)

find_library(AMD_LIBRARIES
    NAMES amd
    PATHS
        ${UMFPACK_INCLUDE_DIR}
        ${UMFPACK_INCLUDE_DIR}/..
        ${UMFPACK_ROOT}
        ${UMFPACK_ROOT_DIR}
    ENV
        UMFPACK_ROOT
        UMFPACK_ROOT_DIR
    PATH_SUFFIXES 
        lib
        lib64
        Lib
        Lib64
)

find_library(CONF_LIBRARIES
    NAMES suitesparseconfig
    PATHS
        ${UMFPACK_INCLUDE_DIR}
        ${UMFPACK_INCLUDE_DIR}/..
        ${UMFPACK_ROOT}
        ${UMFPACK_ROOT_DIR}
    ENV
        UMFPACK_ROOT
        UMFPACK_ROOT_DIR
    PATH_SUFFIXES 
        lib
        lib64
        Lib
        Lib64
)

list(APPEND UMFPACK_LIBRARIES ${COLAMD_LIBRARIES} ${CHOLMOD_LIBRARIES} ${AMD_LIBRARIES} ${CONF_LIBRARIES})

find_package_handle_standard_args( UMFPACK
    REQUIRED_VARS BLAS_FOUND UMFPACK_LIBRARIES UMFPACK_INCLUDE_DIRS
    VERSION_VAR UMFPACK_VERSION
)

mark_as_advanced(
    UMFPACK_LIBRARIES
    UMFPACK_INCLUDE_DIRS
)
