# =============================================================================
#  CADET - The Chromatography Analysis and Design Toolkit
#  
#  Copyright © 2008-2015: Eric von Lieres¹, Joel Andersson¹,
#                         Andreas Puettmann¹, Sebastian Schnittert¹,
#                         Samuel Leweke¹
#                                      
#    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
#  
#  All rights reserved. This program and the accompanying materials
#  are made available under the terms of the GNU Public License v3.0 (or, at
#  your option, any later version) which accompanies this distribution, and
#  is available at http://www.gnu.org/licenses/gpl.html
# =============================================================================

# Find CADET, the Chromatography Analysis and DEsign Toolkit.
#
# On UNIX systems, this module will read the variable CADET_USE_STATIC_LIBRARIES
# to determine whether or not to prefer a static link over a dynamic link for CADET
# and all of it's dependencies. To use this feature, make sure that the
# CADET_USE_STATIC_LIBRARIES variable is set before the call to find_package.
#
# To provide the module with a hint about where to find your CADET installation,
# you can set the environment variable CADET_ROOT. The FindCADET module will then
# look in this path when searching for CADET paths and libraries.
#
# This module will define the following variables:
#  CADET_INCLUDE_DIRS - Location of the CADET includes
#  CADET_FOUND - true if CADET was found on the system
#  CADET_LIBRARIES - Required libraries for all requested components

include(FindPackageHandleStandardArgs)

set (CADET_LIB_NAME cadet)

# Option that allows users to specify custom CADET path
if (NOT "$ENV{CADET_ROOT}" STREQUAL "")
    list (APPEND CMAKE_INCLUDE_PATH "$ENV{CADET_ROOT}")
    list (APPEND CMAKE_LIBRARY_PATH "$ENV{CADET_ROOT}")
#else ()
#    message (STATUS "No valid environment variable 'CADET_ROOT'     found, you can specify custom path in 'PATH_CADET_ROOT' cache variable")
endif ()


#if (NOT PATH_CADET_ROOT)
#    set (PATH_CADET_ROOT "" CACHE STRING "Optional path to the CADET library and include direrctory" FORCE)
#else ()
#    list (APPEND CMAKE_INCLUDE_PATH "${PATH_CADET_ROOT}")
#    list (APPEND CMAKE_LIBRARY_PATH "${PATH_CADET_ROOT}")
#endif ()



# List of user definable search paths
set (CADET_USER_PATHS
    $ENV{HOME}/.local
    $ENV{HOME}
    $ENV{HOME}/libs/cadet
    $ENV{HOME}/opt/cadet
)

#if( CADET_INCLUDE_DIRS AND CADET_LIBRARIES )
#    # Do nothing: we already have CADET_INCLUDE_PATH and CADET_LIBRARIES in the
#    # cache, it would be a shame to override them
#else()

    # find the CADET include directories
    find_path( CADET_INCLUDE_DIRS Simulator.hpp
        ENV
            CADET_ROOT
        PATHS
            ${CADET_USER_PATHS}
        PATH_SUFFIXES
            include
            Include
    )

    # find the CADET libraries
    if( UNIX AND CADET_USE_STATIC_LIBRARIES )
        # According to bug 1643 on the CMake bug tracker, this is the
        # preferred method for searching for a static library.
        # See http://www.cmake.org/Bug/view.php?id=1643.  We search
        # first for the full static library name, but fall back to a
        # generic search on the name if the static search fails.
        set( THIS_LIBRARY_SEARCH lib${CADET_LIB_NAME}.a ${CADET_LIB_NAME} )
    else()
        set( CMAKE_FIND_LIBRARY_SUFFIXES ".lib;.dll" )
        set( THIS_LIBRARY_SEARCH ${CADET_LIB_NAME} )
    endif()

    find_library( CADET_LIBRARIES
        NAMES ${THIS_LIBRARY_SEARCH}
        ENV
            CADET_ROOT
        PATHS
            ${CADET_USER_PATHS}
        PATH_SUFFIXES
            lib
            Lib
    )

#endif()

find_package_handle_standard_args( CADET DEFAULT_MSG
    CADET_LIBRARIES
    CADET_INCLUDE_DIRS
)

mark_as_advanced(
    CADET_LIBRARIES
    CADET_INCLUDE_DIRS
)
