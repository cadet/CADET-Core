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

# Find SuperLU, the direct solver for sparse linear systems.
#
# To provide the module with a hint about where to find your SuperLU installation,
# you can set the environment variable SUPERLU_ROOT. The FindSuperLU module will
# then look in this path when searching for SuperLU paths and libraries.
#
# This module will define the following variables:
#  SUPERLU_FOUND - true if SuperLU was SuperLU on the system
#  SUPERLU_INCLUDE_DIRS - Location of the SuperLU includes
#  SUPERLU_LIBRARIES - Required libraries for all requested components
#  SUPERLU_VERSION_MAJOR - Major version
#  SUPERLU_VERSION_MINOR - Minro version
#  SUPERLU_VERSION_PATCH - Patch level
#  SUPERLU_VERSION - Full version string
#  SUPERLU_INT_TYPE - Integer type for indices used by SuperLU


include(FindPackageHandleStandardArgs)

# find the SuperLU include directories
find_path( SUPERLU_INCLUDE_DIRS supermatrix.h
    PATHS
        ${SUPERLU_ROOT}
        ${SUPERLU_ROOT_DIR}
    ENV
        SUPERLU_ROOT
        SUPERLU_ROOT_DIR
    PATH_SUFFIXES
        include
        include/superlu
        superlu
        SuperLU
)

if(SUPERLU_INCLUDE_DIRS)
    # extract version
    file(READ "${SUPERLU_INCLUDE_DIRS}/slu_util.h" _SUPERLU_VERSION_FILE)

    string(REGEX REPLACE ".*#define SUPERLU_MAJOR_VERSION [ \t\r\n]* ([0-9]+).*" "\\1" SUPERLU_VERSION_MAJOR "${_SUPERLU_VERSION_FILE}")
    string(REGEX REPLACE ".*#define SUPERLU_MINOR_VERSION [ \t\r\n]* ([0-9]+).*" "\\1" SUPERLU_VERSION_MINOR "${_SUPERLU_VERSION_FILE}")
    string(REGEX REPLACE ".*#define SUPERLU_PATCH_VERSION [ \t\r\n]* ([0-9]+).*" "\\1" SUPERLU_VERSION_PATCH "${_SUPERLU_VERSION_FILE}")
    set(SUPERLU_VERSION "${SUPERLU_VERSION_MAJOR}.${SUPERLU_VERSION_MINOR}.${SUPERLU_VERSION_PATCH}")


    # extract int type
    file(READ "${SUPERLU_INCLUDE_DIRS}/slu_sdefs.h" _SUPERLU_INTTYPE_FILE)
    string(REGEX REPLACE ".*typedef (.*) int_t;.*.*" "\\1" SUPERLU_INT_TYPE "${_SUPERLU_INTTYPE_FILE}")

endif()


# find the SuperLU libraries
find_library(SUPERLU_LIBRARIES
    NAMES superlu
        superlu_5.2.1
        superlu_5.2
        superlu_5.1.1
        superlu_5.1
        superlu_5.0
        superlu_4.3
        superlu_4.2
        superlu_4.1
        superlu_4.0
    PATHS
        ${SUPERLU_INCLUDE_DIR}
        ${SUPERLU_ROOT}
        ${SUPERLU_ROOT_DIR}
    ENV
        SUPERLU_ROOT
        SUPERLU_ROOT_DIR
    PATH_SUFFIXES 
        lib
        lib64
        Lib
        Lib64
)

find_package_handle_standard_args( SUPERLU
    REQUIRED_VARS BLAS_FOUND SUPERLU_LIBRARIES SUPERLU_INCLUDE_DIRS
    VERSION_VAR SUPERLU_VERSION
)

mark_as_advanced(
    SUPERLU_LIBRARIES
    SUPERLU_INCLUDE_DIRS
)
