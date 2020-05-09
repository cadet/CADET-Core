# =============================================================================
#  CADET - The Chromatography Analysis and Design Toolkit
#  
#  Copyright © 2008-2020: The CADET Authors
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
# Static libraries can be preferred by setting UMFPACK_PREFER_STATIC_LIBS to TRUE.
#
# This module will define the following variables:
#  UMFPACK_FOUND - true if UMFPACK was UMFPACK on the system
#  UMFPACK_INCLUDE_DIRS - Location of the UMFPACK includes
#  UMFPACK_LIBRARIES - Required libraries for all requested components
#  UMFPACK_VERSION_MAJOR - Major version
#  UMFPACK_VERSION_MINOR - Minro version
#  UMFPACK_VERSION_PATCH - Patch level
#  UMFPACK_VERSION - Full version string
#
# This module will also create the UMFPACK::UMFPACK target.


if (NOT DEFINED UMFPACK_PREFER_STATIC_LIBS)
    set(UMFPACK_PREFER_STATIC_LIBS OFF)
endif()

find_package(PkgConfig QUIET)
if (PKG_CONFIG_FOUND)
    pkg_check_modules(PKGCONFIG_UMFPACK QUIET SuiteSparse)
else()
    set(PKGCONFIG_UMFPACK_INCLUDE_DIRS "")
endif()

# find the UMFPACK include directories
find_path(UMFPACK_INCLUDE_DIRS umfpack.h
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

if (UMFPACK_INCLUDE_DIRS)
    # extract version
    file(READ "${UMFPACK_INCLUDE_DIRS}/umfpack.h" _UMFPACK_VERSION_FILE)

    string(REGEX REPLACE ".*#define UMFPACK_MAIN_VERSION[ \t]*([0-9]+)[ \t\r\n]*.*" "\\1" UMFPACK_VERSION_MAJOR "${_UMFPACK_VERSION_FILE}")
    string(REGEX REPLACE ".*#define UMFPACK_SUB_VERSION[ \t]*([0-9]+)[ \t\r\n]*.*" "\\1" UMFPACK_VERSION_MINOR "${_UMFPACK_VERSION_FILE}")
    string(REGEX REPLACE ".*#define UMFPACK_SUBSUB_VERSION[ \t]*([0-9]+)[ \t\r\n]*.*" "\\1" UMFPACK_VERSION_PATCH "${_UMFPACK_VERSION_FILE}")
    set(UMFPACK_VERSION "${UMFPACK_VERSION_MAJOR}.${UMFPACK_VERSION_MINOR}.${UMFPACK_VERSION_PATCH}")

endif()

# prefer static libs by prioritizing .lib and .a suffixes in CMAKE_FIND_LIBRARY_SUFFIXES
if (UMFPACK_PREFER_STATIC_LIBS)
    set(_UMFPACK_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
    if (WIN32)
        list(INSERT CMAKE_FIND_LIBRARY_SUFFIXES 0 .lib .a)
    else()
        list(INSERT CMAKE_FIND_LIBRARY_SUFFIXES 0 .a)
    endif()
endif()

# find the UMFPACK libraries
find_library(UMFPACK_LIBRARY
    NAMES
        umfpack
        umfpack.4
        umfpack.5
        libumfpack
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
find_library(CHOLMOD_LIBRARY
    NAMES
        cholmod
        libcholmod
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

find_library(COLAMD_LIBRARY
    NAMES
        colamd
        libcolamd
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

find_library(CCOLAMD_LIBRARY
    NAMES
        ccolamd
        libccolamd
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

find_library(CAMD_LIBRARY
    NAMES
        camd
        libcamd
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

find_library(AMD_LIBRARY
    NAMES
        amd
        libamd
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

find_library(CONF_LIBRARY
    NAMES
        suitesparseconfig
        libsuitesparseconfig
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

if (NOT DEFINED BLAS_FOUND)
    find_package(BLAS QUIET)
endif()

list(APPEND UMFPACK_LIBRARIES ${UMFPACK_LIBRARY} ${COLAMD_LIBRARY} ${CCOLAMD_LIBRARY} ${AMD_LIBRARY} ${CAMD_LIBRARY} ${CHOLMOD_LIBRARY} ${CONF_LIBRARY} ${BLAS_LIBRARIES})

find_library(METIS_LIBRARY
    NAMES
        metis
        libmetis
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

if (METIS_LIBRARY)
    list(APPEND UMFPACK_LIBRARIES ${METIS_LIBRARY})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(UMFPACK
    REQUIRED_VARS BLAS_FOUND UMFPACK_LIBRARIES UMFPACK_INCLUDE_DIRS
    VERSION_VAR UMFPACK_VERSION
)

mark_as_advanced(
    UMFPACK_LIBRARIES
    UMFPACK_INCLUDE_DIRS
)

# restore original find_library suffixes
if (UMFPACK_PREFER_STATIC_LIBS)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${_UMFPACK_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
endif()

include(FeatureSummary)
set_package_properties(UMFPACK PROPERTIES
    URL "http://faculty.cse.tamu.edu/davis/suitesparse.html"
    DESCRIPTION "UMFPACK sparse direct solver from SuiteSparse"
)

if (UMFPACK_FOUND AND NOT TARGET UMFPACK::UMFPACK)
    # it is unknown whether we have caught shared or static library

    add_library(UMFPACK::COLAMD UNKNOWN IMPORTED)
    set_target_properties(UMFPACK::COLAMD PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${UMFPACK_INCLUDE_DIRS}"
        IMPORTED_LOCATION ${COLAMD_LIBRARY}
    )

    add_library(UMFPACK::CCOLAMD UNKNOWN IMPORTED)
    set_target_properties(UMFPACK::CCOLAMD PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${UMFPACK_INCLUDE_DIRS}"
        IMPORTED_LOCATION ${CCOLAMD_LIBRARY}
    )

    add_library(UMFPACK::AMD UNKNOWN IMPORTED)
    set_target_properties(UMFPACK::AMD PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${UMFPACK_INCLUDE_DIRS}"
        IMPORTED_LOCATION ${AMD_LIBRARY}
    )

    add_library(UMFPACK::CAMD UNKNOWN IMPORTED)
    set_target_properties(UMFPACK::CAMD PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${UMFPACK_INCLUDE_DIRS}"
        IMPORTED_LOCATION ${CAMD_LIBRARY}
    )

    add_library(UMFPACK::CHOLMOD UNKNOWN IMPORTED)
    set_target_properties(UMFPACK::CHOLMOD PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${UMFPACK_INCLUDE_DIRS}"
        IMPORTED_LOCATION ${CHOLMOD_LIBRARY}
    )

    add_library(UMFPACK::config UNKNOWN IMPORTED)
    set_target_properties(UMFPACK::config PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${UMFPACK_INCLUDE_DIRS}"
        IMPORTED_LOCATION ${CONF_LIBRARY}
    )

    add_library(UMFPACK::UMFPACK UNKNOWN IMPORTED)
    set_target_properties(UMFPACK::UMFPACK PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${UMFPACK_INCLUDE_DIRS}"
        IMPORTED_LOCATION ${UMFPACK_LIBRARY}
    )

    target_link_libraries(UMFPACK::UMFPACK INTERFACE 
        UMFPACK::COLAMD
        UMFPACK::CCOLAMD
        UMFPACK::AMD
        UMFPACK::CAMD
        UMFPACK::CHOLMOD
        UMFPACK::config
        ${BLAS_LIBRARIES}
    )

    if (METIS_LIBRARY)
        add_library(UMFPACK::METIS UNKNOWN IMPORTED)
        set_target_properties(UMFPACK::METIS PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${UMFPACK_INCLUDE_DIRS}"
            IMPORTED_LOCATION ${METIS_LIBRARY}
        )
        target_link_libraries(UMFPACK::UMFPACK INTERFACE UMFPACK::METIS)
    endif()

endif()
