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

# Name of the current project
project (${CADET_MI_PRJ} CXX)

# Set variables defining source and include directories
set (CADET_MI_SRC_DIR "${CMAKE_SOURCE_DIR}/src/cadet-mi")

set (COMMON_INC_DIR "${CMAKE_SOURCE_DIR}/include/common")
set (COMMON_SRC_DIR "${CMAKE_SOURCE_DIR}/src/common")

# Add include directories for access to exported LIBCADET header files.
include_directories (${LIBCADET_INC_DIR} ${COMMON_INC_DIR})

# Info message
message (STATUS "Added CADET-MI module")

# Reduce Warning level in VC++
if (CMAKE_BUILD_TOOL MATCHES "(msdev|devenv|nmake)")
    add_definitions (/W2)
endif()


# ------------------------------------------------------------------------------------------------------------------------
#   Create mex file 
# ------------------------------------------------------------------------------------------------------------------------
if (MATLAB_FOUND)

	# Add include directories for Matlab header files
	include_directories (${MATLAB_INCLUDE_DIR})

    # Add the build target for the static chromsim library
    add_library (cadetmex SHARED CadetMex.cpp)

    # Set the mex file name
    set_target_properties (cadetmex
        PROPERTIES
            OUTPUT_NAME "CadetMex" PREFIX "" SUFFIX ".${MATLAB_MEX_EXT}")

	# Add dependency on the LIBCADET library
	add_dependencies(cadetmex ${LIBCADET_NAME_MEX})

	# Link
	target_link_libraries(cadetmex ${LIBCADET_NAME_MEX} ${MATLAB_MEX_LIB} ${MATLAB_MX_LIB} ${MATLAB_MAT_LIB})

    # use, i.e. don't skip the full RPATH for the build tree
    set_target_properties (cadetmex PROPERTIES SKIP_BUILD_RPATH  FALSE)

    # when building, don't use the install RPATH already
    # (but later on when installing)
    set_target_properties (cadetmex PROPERTIES BUILD_WITH_INSTALL_RPATH FALSE) 

    set_target_properties (cadetmex PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    set_target_properties (cadetmex PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)


    # the RPATH to be used when installing, but only if it's not a system directory
    list (FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" ISSYSTEMDIR)
    if ("${ISSYSTEMDIR}" STREQUAL "-1")
       set_target_properties (cadetmex PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    endif ("${ISSYSTEMDIR}" STREQUAL "-1")

    if (CXX_COMPILER_NAME STREQUAL "icpc")
        set_target_properties (cadetmex PROPERTIES LINK_FLAGS "-pthread -openmp")
        set_target_properties (cadetmex PROPERTIES COMPILE_FLAGS "-std=c++0x -fPIC -fno-omit-frame-pointer -pthread -openmp -mkl")
        add_definitions(-DMATLAB_MEX_FILE -D_GNU_SOURCE)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set_target_properties (cadetmex PROPERTIES LINK_FLAGS "-fno-common -fexceptions")
        set_target_properties (cadetmex PROPERTIES COMPILE_FLAGS "-fPIC -fno-omit-frame-pointer -pthread -fopenmp -std=c++0x")
        add_definitions(-DMATLAB_MEX_FILE) 
    elseif (CMAKE_COMPILER_IS_GNUCXX OR CXX_COMPILER_NAME STREQUAL "g++" OR CXX_COMPILER_NAME STREQUAL "c++" OR CXX_COMPILER_NAME STREQUAL "g++.exe")
        set_target_properties (cadetmex PROPERTIES LINK_FLAGS "-pthread -fopenmp")
        set_target_properties (cadetmex PROPERTIES COMPILE_FLAGS "-fPIC -fno-omit-frame-pointer -pthread -fopenmp -std=c++0x")
        add_definitions(-DMATLAB_MEX_FILE -D_GNU_SOURCE) 
    elseif (CXX_COMPILER_NAME STREQUAL "cl")
        set_target_properties (cadetmex PROPERTIES COMPILE_FLAGS "/openmp")
        add_definitions(-DMATLAB_MEX_FILE) 
    endif ()

    # Install mex file to source directory
    install (TARGETS cadetmex DESTINATION ${CADET_MI_SRC_DIR})

endif ()
# ------------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------------
#   Set installation related stuff
# ------------------------------------------------------------------------------------------------------------------------

# Install the cadet-cs executable
install (CODE "MESSAGE(\"\nInstall CADET-MI\n\")")
#install (TARGETS ${EXPORT_ENUM_NAME} DESTINATION bin)

# ------------------------------------------------------------------------------------------------------------------------
