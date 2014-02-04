# =============================================================================
#  CADET - The Chromatography Analysis and Design Toolkit
#  
#  Copyright © 2008-2014: Eric von Lieres¹, Joel Andersson,
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
project (${LIBCADET_PRJ})

# Set variables defining source and include directories
set (LIBCADET_SRC_DIR "${CMAKE_SOURCE_DIR}/src/libcadet")
set (LIBCADET_INC_DIR "${CMAKE_SOURCE_DIR}/include/libcadet" CACHE INTERNAL "Path to exported libcadet headers")

set (COMMON_SRC_DIR "${CMAKE_SOURCE_DIR}/src/common")
set (COMMON_INC_DIR "${CMAKE_SOURCE_DIR}/include/common")

# Git integration
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Cadet.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/Cadet.cpp" @ONLY)

# Add variable LIBCADET_SOURCES with the sources for LIBCADET
set (LIBCADET_SOURCES
     ${LIBCADET_SRC_DIR}/CadetConvenience.cpp
     ${LIBCADET_SRC_DIR}/GeneralRateModel.cpp
     ${LIBCADET_SRC_DIR}/JacobianData.cpp
     ${LIBCADET_SRC_DIR}/ParticleDiscretization.cpp
     ${LIBCADET_SRC_DIR}/SchurSolver.cpp
     ${LIBCADET_SRC_DIR}/SimulatorPImpl.cpp
     ${LIBCADET_SRC_DIR}/TimeIntegrator.cpp
     ${LIBCADET_SRC_DIR}/WenoScheme.cpp
     ${LIBCADET_SRC_DIR}/mylapack.cpp
     ${CMAKE_CURRENT_BINARY_DIR}/Cadet.cpp
    )

# Add variable LIBCADET_HEADERS with the exported LIBCADET header files
set (LIBCADET_HEADERS
     ${LIBCADET_INC_DIR}/CadetEnumeration.hpp
     ${LIBCADET_INC_DIR}/CadetException.hpp
     ${LIBCADET_INC_DIR}/Simulator.hpp
    )

# Add include directories for access to exported header files.
include_directories (${LIBCADET_INC_DIR} ${COMMON_INC_DIR})



# Specify the location of additional CMAKE modules
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")



# ------------------------------------------------------------------------------------------------------------------------
#   Other configuration options
# ------------------------------------------------------------------------------------------------------------------------
# Option that allows users to build static and/or shared libraries
option (BUILD_SHARED_LIBS "Build shared libraries" ON)
option (BUILD_STATIC_LIBS "Build static libraries" OFF)
option (BUILD_MEX_LIBS    "Build MATLAB mex interface" ON)


# Make sure we build at least one type of libraries
if (NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)
    print_warning ("Both static and shared library generation were disabled"
        "Building shared libraries was re-enabled")
    set (BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries" FORCE)
endif ()


# Specify the VERSION and SOVERSION for shared libraries
set (LIBCADET_VERSION   "0.0.0")
set (LIBCADET_SOVERSION "0")


# Option that allows users to switch on verification of analytical Jacobian by AD
option (VERIFY_ANALYTICAL_JAC "Toggle verification of analytical Jacobian by AD (only DEBUG builds)" OFF)
if (VERIFY_ANALYTICAL_JAC AND CMAKE_BUILD_TYPE STREQUAL "debug")
    add_definitions (-DVERIFY_ANALYTICAL_JAC)
endif ()


# Option that allows users to specify the datatype used as active type for AD code parts
set (USE_ACTIVE_TYPE "adouble" CACHE STRING "Defines the type used as active type for AD code parts, options are 'adouble'" FORCE)
if (USE_ACTIVE_TYPE STREQUAL "adouble")
    add_definitions (-DACTIVE_ADOLC)             # Use ADOLC as active type

    # Add ADOLC headers
    include_directories ("${CMAKE_SOURCE_DIR}/ThirdParty/ADOL-C/include")
endif ()
# ------------------------------------------------------------------------------------------------------------------------





# ------------------------------------------------------------------------------------------------------------------------
#   Build the static library
# ------------------------------------------------------------------------------------------------------------------------
if (BUILD_STATIC_LIBS)

    set (LIBCADET_NAME_STATIC "${LIBCADET_NAME}_static" CACHE INTERNAL "Name of static cadet library target")

    # Add the build target for the static chromsim library
    add_library (${LIBCADET_NAME_STATIC} STATIC ${LIBCADET_SOURCES})

    # Set the library name and make sure it is not deleted
    set_target_properties (${LIBCADET_NAME_STATIC}
        PROPERTIES OUTPUT_NAME ${CADET_NAME} CLEAN_DIRECT_OUTPUT 1)

    list (APPEND LIBCADET_TARGETS ${LIBCADET_NAME_STATIC})

endif ()
# ------------------------------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------------------------------------------
#   Build the shared library
# ------------------------------------------------------------------------------------------------------------------------
if (BUILD_SHARED_LIBS)

    set (LIBCADET_NAME_SHARED "${LIBCADET_NAME}_shared" CACHE INTERNAL "Name of shared cadet library target")

    # Add the build target for the static chromsim library
    add_library (${LIBCADET_NAME_SHARED} SHARED ${LIBCADET_SOURCES})

    # Set the library name and make sure it is not deleted
    set_target_properties (${LIBCADET_NAME_SHARED}
        PROPERTIES
            OUTPUT_NAME ${CADET_NAME} CLEAN_DIRECT_OUTPUT 1)

    # Set VERSION and SOVERSION for shared libraries
    set_target_properties(${LIBCADET_NAME_SHARED}
        PROPERTIES
            VERSION ${LIBCADET_VERSION} SOVERSION ${LIBCADET_SOVERSION})

    list (APPEND LIBCADET_TARGETS ${LIBCADET_NAME_SHARED})
endif ()
# ------------------------------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------------------------------------------
#   Build the MATLAB mex interface
# ------------------------------------------------------------------------------------------------------------------------
find_package (MATLAB)
if (MATLAB_FOUND AND BUILD_MEX_LIBS)

    set (LIBCADET_NAME_MEX "${LIBCADET_NAME}_mex" CACHE INTERNAL "Name of shared cadet library target")

    # Add include directories for Matlab header files
    include_directories (${MATLAB_INCLUDE_DIR})

    # Add the build target for the static chromsim library
    add_library (${LIBCADET_NAME_MEX} SHARED ${LIBCADET_SOURCES})

    # Set the library name and make sure it is not deleted
    set_target_properties (${LIBCADET_NAME_MEX}
        PROPERTIES
            OUTPUT_NAME "${CADET_NAME}-mex" CLEAN_DIRECT_OUTPUT 1)

    # Set VERSION and SOVERSION for shared libraries
    set_target_properties (${LIBCADET_NAME_MEX}
        PROPERTIES
            VERSION ${LIBCADET_VERSION} SOVERSION ${LIBCADET_SOVERSION})

    list (APPEND LIBCADET_TARGETS ${LIBCADET_NAME_MEX})

    # Make variables visible in parent scope
    set (LIBCADET_NAME_MEX "${LIBCADET_NAME_MEX}" PARENT_SCOPE)
    set (MATLAB_FOUND TRUE PARENT_SCOPE)
    set (MATLAB_LIBRARIES ${MATLAB_LIBRARIES} PARENT_SCOPE)
    set (MATLAB_ROOT ${MATLAB_ROOT} PARENT_SCOPE)
    set (MATLAB_MEX_EXT ${MATLAB_MEX_EXT} PARENT_SCOPE)
    set (MATLAB_ROOT ${MATLAB_ROOT} PARENT_SCOPE)
    set (MATLAB_LIB_DIR ${MATLAB_LIB_DIR} PARENT_SCOPE)
    set (MATLAB_MEX_LIB ${MATLAB_MEX_LIB} PARENT_SCOPE)
    set (MATLAB_MX_LIB ${MATLAB_MX_LIB} PARENT_SCOPE)
    set (MATLAB_MAT_LIB ${MATLAB_MAT_LIB} PARENT_SCOPE)
    set (MATLAB_INCLUDE_DIR ${MATLAB_INCLUDE_DIR} PARENT_SCOPE)

else ()

    # Make variable visible in parent scope
    set (MATLAB_FOUND FALSE PARENT_SCOPE)

endif ()
# ------------------------------------------------------------------------------------------------------------------------


# Make variable visible in parent scope
set (LIBCADET_TARGETS ${LIBCADET_TARGETS} PARENT_SCOPE)


# Link to SUNDIALS
find_package (SUNDIALS REQUIRED COMPONENTS sundials_idas sundials_nvecserial)
foreach (_TARGET ${LIBCADET_TARGETS})
    target_link_libraries (${_TARGET} ${SUNDIALS_LIBRARIES})
endforeach ()
include_directories (${SUNDIALS_INCLUDE_DIRS})


# Link to native BLAS and LAPACK on all other than WINDOWS machines
if (NOT WIN32)
    find_package (LAPACK)
    foreach (_TARGET ${LIBCADET_TARGETS})
        target_link_libraries (${_TARGET} ${LAPACK_LIBRARIES})
    endforeach ()
else () # On WINDOWS environment var LAPACK_ROOT must be set!
    find_library (F2C_LIBRARIES
                      NAMES         libf2c
                      PATHS         ENV LAPACK_ROOT
                      PATH_SUFFIXES lib
                 )
    find_library (BLAS_LIBRARIES
                      NAMES         blas
                      PATHS         ENV LAPACK_ROOT
                      PATH_SUFFIXES lib
                 )
    find_library (CLAPACK_LIBRARIES
                      NAMES         lapack
                      PATHS         ENV LAPACK_ROOT
                      PATH_SUFFIXES lib
                 )
    foreach (_TARGET ${LIBCADET_TARGETS})
        target_link_libraries (${_TARGET} ${F2C_LIBRARIES} ${BLAS_LIBRARIES} ${CLAPACK_LIBRARIES})
    endforeach ()
    mark_as_advanced(FORCE F2C_LIBRARIES BLAS_LIBRARIES CLAPACK_LIBRARIES)
endif ()

# Set CADET_VERSION 
set_target_properties (${LIBCADET_TARGETS}
                       PROPERTIES COMPILE_DEFINITIONS CADET_VERSION="${CADET_VERSION}"
                       COMPILE_DEFINITIONS_DEBUG CADET_VERSION="${CADET_VERSION}_DEBUG")

# Reduce Warning level in VC++
if(CMAKE_BUILD_TOOL MATCHES "(msdev|devenv|nmake)")
    add_definitions(/W2)
endif()


# ------------------------------------------------------------------------------------------------------------------------
#   Use full RPATH on the mex library
# ------------------------------------------------------------------------------------------------------------------------
if (MATLAB_FOUND AND BUILD_MEX_LIBS)

    # link to MEX libraries
    target_link_libraries(${LIBCADET_NAME_MEX} ${MATLAB_MEX_LIB} ${MATLAB_MX_LIB} ${MATLAB_MAT_LIB})

    # use, i.e. don't skip the full RPATH for the build tree
    set_target_properties (${LIBCADET_NAME_MEX} PROPERTIES SKIP_BUILD_RPATH  FALSE)

    # when building, don't use the install RPATH already
    # (but later on when installing)
    set_target_properties (${LIBCADET_NAME_MEX} PROPERTIES BUILD_WITH_INSTALL_RPATH FALSE) 

    set_target_properties (${LIBCADET_NAME_MEX} PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    set_target_properties (${LIBCADET_NAME_MEX} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)


    # the RPATH to be used when installing, but only if it's not a system directory
    list (FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" ISSYSTEMDIR)
    if ("${ISSYSTEMDIR}" STREQUAL "-1")
       set_target_properties (${LIBCADET_NAME_MEX} PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    endif ("${ISSYSTEMDIR}" STREQUAL "-1")

    # Add MATLABMEX and MATLAB_MEX_FILE define to compile flags
    get_target_property (MEX_DEFINES ${LIBCADET_NAME_MEX} COMPILE_DEFINITIONS)
    
    if (MEX_DEFINES)
        list (APPEND MEX_DEFINES "MATLABMEX=1")
        list (APPEND MEX_DEFINES "MATLAB_MEX_FILE=1")
    else ()
        set (MEX_DEFINES "MATLABMEX=1" "MATLAB_MEX_FILE=1")
    endif ()
    set_target_properties (${LIBCADET_NAME_MEX} PROPERTIES COMPILE_DEFINITIONS "${MEX_DEFINES}")

endif ()
# ------------------------------------------------------------------------------------------------------------------------



# Info message
message (STATUS "Added LIBCADET module")


# ------------------------------------------------------------------------------------------------------------------------
#   Set installation related stuff
# ------------------------------------------------------------------------------------------------------------------------

# Install the LIBCADET library and header files
install (CODE "MESSAGE(\"\nInstall LIBCADET\n\")")
install (TARGETS ${LIBCADET_TARGETS} DESTINATION lib)
install (FILES ${LIBCADET_HEADERS} DESTINATION include)

# ------------------------------------------------------------------------------------------------------------------------
