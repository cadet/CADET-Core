# ---------------------------------------------------------------
# Programmer: Radu Serban, David J. Gardner, Cody J. Balos,
#             and Slaven Peles @ LLNL
# ---------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department
# of Energy by Lawrence Livermore National Laboratory in part under
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# ---------------------------------------------------------------
# Top level CMakeLists.txt for SUNDIALS (for cmake build system)
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Initial commands
# ---------------------------------------------------------------

# Libraries linked via full path no longer produce linker search paths
# Allows examples to build
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

# MACOSX_RPATH is enabled by default
# Fixes dynamic loading on OSX
if(POLICY CMP0042)
  cmake_policy(SET CMP0042 NEW) # Added in CMake 3.0
else()
  if(APPLE)
    set(CMAKE_MACOSX_RPATH 1)
  endif()
endif()

# Project SUNDIALS (initially only C supported)
# sets PROJECT_SOURCE_DIR and PROJECT_BINARY_DIR variables
PROJECT(sundials C)

# Set some variables with info on the SUNDIALS project
SET(PACKAGE_BUGREPORT "woodward6@llnl.gov")
SET(PACKAGE_NAME "SUNDIALS")
SET(PACKAGE_STRING "SUNDIALS 3.2.1")
SET(PACKAGE_TARNAME "sundials")

# set SUNDIALS version numbers
# (use "" for the version label if none is needed)
SET(PACKAGE_VERSION_MAJOR "3")
SET(PACKAGE_VERSION_MINOR "2")
SET(PACKAGE_VERSION_PATCH "1")
SET(PACKAGE_VERSION_LABEL "")

IF(PACKAGE_VERSION_LABEL)
  SET(PACKAGE_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}-${PACKAGE_VERSION_LABEL}")
ELSE()
  SET(PACKAGE_VERSION "${PACKAGE_VERSION_MAJOR}.${PACKAGE_VERSION_MINOR}.${PACKAGE_VERSION_PATCH}")
ENDIF()

# Hide some cache variables
MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)

# Always show the C compiler and flags
MARK_AS_ADVANCED(CLEAR
  CMAKE_C_COMPILER
  CMAKE_C_FLAGS)

# Specify the VERSION and SOVERSION for shared libraries

SET(idaslib_VERSION "2.2.1")
SET(idaslib_SOVERSION "2")

SET(nveclib_VERSION "3.2.1")
SET(nveclib_SOVERSION "3")

# Specify the location of additional CMAKE modules
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

# Get correct build paths automatically, but expose CMAKE_INSTALL_LIBDIR
# as a regular cache variable so that a user can more easily see what
# the library dir was set to be by GNUInstallDirs.
INCLUDE(GNUInstallDirs)
MARK_AS_ADVANCED(CLEAR CMAKE_INSTALL_LIBDIR)

# ---------------------------------------------------------------
# MACRO definitions
# ---------------------------------------------------------------
INCLUDE(SundialsCMakeMacros)

# ---------------------------------------------------------------
# Check for deprecated SUNDIALS CMake options/variables
# ---------------------------------------------------------------
INCLUDE(SundialsDeprecated)

# ---------------------------------------------------------------
# Which modules to build?
# ---------------------------------------------------------------

# For each SUNDIALS solver available (i.e. for which we have the
# sources), give the user the option of enabling/disabling it.

OPTION(BUILD_IDAS "Build the IDAS library" ON)

# ---------------------------------------------------------------
# Option to specify precision (realtype)
# ---------------------------------------------------------------

SET(SUNDIALS_PRECISION "double")
SET(FPRECISION_LEVEL "8")

# ---------------------------------------------------------------
# Option to specify index type
# ---------------------------------------------------------------

SET(DOCSTR "Signed 64-bit (64) or signed 32-bit (32) integer")
SHOW_VARIABLE(SUNDIALS_INDEX_SIZE STRING "${DOCSTR}" "32")
SET(DOCSTR "Integer type to use for indices in SUNDIALS")
SHOW_VARIABLE(SUNDIALS_INDEX_TYPE STRING "${DOCSTR}" "")
MARK_AS_ADVANCED(SUNDIALS_INDEX_TYPE)
include(SundialsIndexSize)

# ---------------------------------------------------------------
# Options to build static and/or shared libraries
# ---------------------------------------------------------------

SET(SUNDIALS_EXPORT "#define SUNDIALS_EXPORT")

# ---------------------------------------------------------------
# Option to use the generic math libraries (UNIX only)
# ---------------------------------------------------------------

IF(UNIX)
  OPTION(USE_GENERIC_MATH "Use generic (std-c) math libraries" ON)
  IF(USE_GENERIC_MATH)
    # executables will be linked against -lm
    SET(EXTRA_LINK_LIBS -lm)
    # prepare substitution variable for sundials_config.h
    SET(SUNDIALS_USE_GENERIC_MATH TRUE)
  ENDIF(USE_GENERIC_MATH)
ENDIF(UNIX)

## clock-monotonic, see if we need to link with rt
include(CheckSymbolExists)
set(CMAKE_REQUIRED_LIBRARIES_SAVE ${CMAKE_REQUIRED_LIBRARIES})
set(CMAKE_REQUIRED_LIBRARIES rt)
CHECK_SYMBOL_EXISTS(_POSIX_TIMERS "unistd.h;time.h" SUNDIALS_POSIX_TIMERS)
set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES_SAVE})
if(SUNDIALS_POSIX_TIMERS)
  find_library(SUNDIALS_RT_LIBRARY NAMES rt)
  mark_as_advanced(SUNDIALS_RT_LIBRARY)
  if(SUNDIALS_RT_LIBRARY)
    # sundials_config.h symbol
    SET(SUNDIALS_HAVE_POSIX_TIMERS TRUE)
    set(EXTRA_LINK_LIBS ${EXTRA_LINK_LIBS} ${SUNDIALS_RT_LIBRARY})
  endif()
endif()


# ===============================================================
# Options for Parallelism
# ===============================================================

OPTION(MPI_ENABLE "Enable MPI support" OFF)
OPTION(OPENMP_ENABLE "Enable OpenMP support" OFF)
OPTION(PTHREAD_ENABLE "Enable Pthreads support" OFF)
OPTION(CUDA_ENABLE "Enable CUDA support" OFF)
OPTION(RAJA_ENABLE "Enable RAJA support" OFF)

# ===============================================================
# Options for external packages
# ===============================================================

# ---------------------------------------------------------------
# Enable BLAS support?
# ---------------------------------------------------------------
OPTION(BLAS_ENABLE "Enable BLAS support" OFF)

# ===============================================================
# Add any platform specifc settings
# ===============================================================

# Under Windows, add compiler directive to inhibit warnings
# about use of unsecure functions

IF(WIN32)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)
ENDIF(WIN32)

IF(APPLE)
  SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -undefined dynamic_lookup")
ENDIF(APPLE)

# ===============================================================
# Add source and configuration files
# ===============================================================

# ---------------------------------------------------------------
# Configure the header file sundials_config.h
# ---------------------------------------------------------------

# All required substitution variables should be available at this point.
# Generate the header file and place it in the binary dir.
CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/include/sundials/sundials_config.in
  ${PROJECT_BINARY_DIR}/include/sundials/sundials_config.h
  )

# Add the include directory in the source tree and the one in
# the binary tree (for the header file sundials_config.h)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include)

# ---------------------------------------------------------------
# Add selected modules to the build system
# ---------------------------------------------------------------

ADD_SUBDIRECTORY(src)

