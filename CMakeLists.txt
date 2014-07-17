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
 
# Require a fairly new cmake version
cmake_minimum_required (VERSION 2.6)

# Prohibit in-source build
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message (FATAL_ERROR "In-source build prohibited.")
endif ()

# Set module path in order to use custom CMake modules
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules")

find_package(Git)

# Write the current version number to variable and file
if (GIT_FOUND)
  if (EXISTS "${CMAKE_SOURCE_DIR}/.git")
    execute_process(COMMAND ${GIT_EXECUTABLE} describe --abbrev=0 HEAD
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                    OUTPUT_VARIABLE CADET_VERSION
                    OUTPUT_STRIP_TRAILING_WHITESPACE)

    # In case of missing tags, default to versions.txt file
    if (NOT "${CADET_VERSION}" STREQUAL "")
      message(STATUS "Get version from git")
      
      # Remove first character ("v")
      string(LENGTH "${CADET_VERSION}" CADET_VERSION_STRLEN)
      math(EXPR CADET_VERSION_STRLEN "${CADET_VERSION_STRLEN}-1")
      string(SUBSTRING "${CADET_VERSION}" 1 ${CADET_VERSION_STRLEN}  CADET_VERSION)
      
      # Save to file; Use shell since writing to source directory is forbidden by CMAKE_DISABLE_SOURCE_CHANGES
      if (UNIX)
        execute_process (COMMAND sh -c "echo \"${CADET_VERSION}\" > ${CMAKE_SOURCE_DIR}/version.txt")
      else ()
        execute_process (COMMAND cmd /C "echo ${CADET_VERSION} > ${CMAKE_SOURCE_DIR}/version.txt")
      endif ()
    endif ()
  endif ()
endif ()

if ("${CADET_VERSION}" STREQUAL "")
  message(STATUS "Get version from file")
  file(STRINGS "${CMAKE_SOURCE_DIR}/version.txt" CADET_VERSION)
endif ()

message(STATUS "CADET version is ${CADET_VERSION}")

# Get current commit hash from git
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
if (NOT DEFINED GIT_SHA1)
  set(GIT_SHA1 "NO-COMMIT-HASH")
  set(GIT_REFSPEC "NO-REFSPEC")
endif()
message(STATUS "Current git HEAD: ${GIT_REFSPEC} SHA1 ${GIT_SHA1}")


# Names of the executables
set (CADET_NAME       "cadet")
set (LIBCADET_NAME    "libcadet")
set (CADET_CS_NAME    "cadet-cs")
set (CADET_CS_MICRO_NAME  "cadet-cs_micro")

set (CADET_MI_NAME    "cadet-mi")

set (CONVERTER_NAME   "converter")
set (H5TOXML_NAME     "h5toxml")
set (XMLTOH5_NAME     "xmltoh5")

set (BENCHMARK_NAME   "benchmarks")
set (SCL_NAME         "createSCL")
set (VAR_SCL_NAME     "var_createSCL")
set (LWE_NAME         "createLWE")
set (VAR_LWE_NAME     "var_createLWE")
set (C14_NAME         "createC14")
set (EXTL_NAME        "createEXTL")

set (TESTS_NAME       "tests")
set (RESULTCOMP_NAME  "rescomp")

# Names of the projects
set (CADET_PRJ        "prj_${CADET_NAME}")
set (LIBCADET_PRJ     "prj_${LIBCADET_NAME}")
set (CADET_CS_PRJ     "prj_${CADET_CS_NAME}")
set (CADET_MI_PRJ     "prj_${CADET_MI_NAME}")

set (CONVERTER_PRJ    "prj_${CONVERTER_NAME}")
set (BENCHMARK_PRJ    "prj_${BENCHMARK_NAME}")
set (TESTS_PRJ        "prj_${TESTS_NAME}")

# Name of the current project
project (${CADET_PRJ} CXX)

string (TOUPPER CADET_NAME CADET_NAME_UP)
# Set some variables with info on the CADET project
set (PACKAGE_VERSION   "${CADET_VERSION}")
set (PACKAGE_BUGREPORT "cadet@fz-juelich.de")
set (PACKAGE_NAME      "${CADET_NAME_UP}")
set (PACKAGE_STRING    "${CADET_NAME_UP} ${PACKAGE_VERSION}")
set (PACKAGE_TARNAME   "${CADET_NAME}")


# ------------------------------------------------------------------------------------------------------------------------
#   Other configuration options
# ------------------------------------------------------------------------------------------------------------------------
# Option that allows users to build release or debug version
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE "release" CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif ()

# Set default compiler compiler flags
if (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel") OR (CXX_COMPILER_NAME STREQUAL "icpc") OR (CXX_COMPILER_NAME STREQUAL "icc"))
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -xHOST -DNDEBUG -ipo -ip" CACHE STRING "Default flags for optimized release builds" FORCE)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native" CACHE STRING "Default flags for optimized release builds" FORCE)
elseif (CMAKE_COMPILER_IS_GNUCXX OR CXX_COMPILER_NAME STREQUAL "g++" OR CXX_COMPILER_NAME STREQUAL "c++" OR CXX_COMPILER_NAME STREQUAL "g++.exe")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native" CACHE STRING "Default flags for optimized release builds" FORCE)
elseif (CXX_COMPILER_NAME STREQUAL "cl")
endif ()

option (BENCHMARK_MODE "In benchmark mode logging framework is disabled, only timing is printed" OFF)
if (BENCHMARK_MODE)
  add_definitions(-DBENCHMARK_MODE)
endif ()
# ------------------------------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------------------------------------------
#   Which modules to build?
# ------------------------------------------------------------------------------------------------------------------------
# CADET-CS module
if (IS_DIRECTORY "${PROJECT_SOURCE_DIR}/src/cadet-cs")
  option (BUILD_CADET_CS  "Build the CADET-CS binary"  ON)
else ()
  set (BUILD_CADET_CS OFF)
endif ()

# CADET-MI module
if (IS_DIRECTORY "${PROJECT_SOURCE_DIR}/src/cadet-mi")
    option (BUILD_CADET_MI "Build the matlab interface tools" ON)
else ()
    set (BUILD_CADET_MI OFF)
endif ()

# CONVERTER module
if (IS_DIRECTORY "${PROJECT_SOURCE_DIR}/src/common/convert")
    option (BUILD_CONVERTER "Build the HDF5 and XML converter" ON)
else ()
    set (BUILD_CONVERTER OFF)
endif ()

# BENCHMARK module
if (IS_DIRECTORY "${PROJECT_SOURCE_DIR}/src/common/benchmark")
    option (BUILD_BENCHMARK "Build the benchmark case creation tools" ON)
else ()
    set (BUILD_BENCHMARK OFF)
endif ()

# TEST module
if (IS_DIRECTORY "${PROJECT_SOURCE_DIR}/test")
    option (BUILD_TESTS "Build tests" ON)
else ()
    set (BUILD_TESTS OFF)
endif ()
# ------------------------------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------------------------------------------
#   Add selected modules to the build system and add the targets to the list of all targets
# ------------------------------------------------------------------------------------------------------------------------

# Shared components
add_subdirectory (src/libcadet)
# Add LIBCADET targets to the list of all targets
set (CADET_TARGETS ${LIBCADET_TARGETS})

# CADET-CS module
if (BUILD_CADET_CS)
    add_subdirectory (src/cadet-cs)
    # Add CADET-CS target
    list (APPEND CADET_TARGETS ${CADET_CS_NAME})
    list (APPEND CADET_TARGETS ${CADET_CS_MICRO_NAME})
endif ()

# CADET-MI module
if (BUILD_CADET_MI)
    add_subdirectory (src/cadet-mi)
    # Add CADET-MI target
endif ()

# CONVERTER module
if (BUILD_CONVERTER)
    add_subdirectory (src/common/convert)
    # Add CONVERTER target
    list (APPEND CADET_TARGETS ${H5TOXML_NAME} ${XMLTOH5_NAME})
endif ()

# BENCHMARK module
if (BUILD_BENCHMARK)
    add_subdirectory (src/common/benchmark)
    # Add BENCHMARK target
    list (APPEND CADET_TARGETS ${SCL_NAME} ${LWE_NAME})
endif ()

# TESTS module
if (BUILD_TESTS)
    add_subdirectory (test)
    # Add TESTS target
    list (APPEND CADET_TARGETS ${RESULTCOMP_NAME})
endif ()
# ------------------------------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------------------------------------------
#   Set properties, definitions, install target etc.
# ------------------------------------------------------------------------------------------------------------------------

# Depending on the used compiler, we enable warnings and OpenMP support: -fopenmp -Wall -Wvla
# Disable line-wrapping in gcc/gxx to support eclipse CDT error parsing: -fmessage-length=0
# Enable C++11 - Standard compiler support: -std=c++0x 
get_filename_component (CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)
message (STATUS "C++ compiler name: ${CXX_COMPILER_NAME}")
foreach (_TARGET ${CADET_TARGETS})
    if (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel") OR (CXX_COMPILER_NAME STREQUAL "icpc") OR (CXX_COMPILER_NAME STREQUAL "icc"))
        set_target_properties (${_TARGET} PROPERTIES COMPILE_FLAGS "-openmp -w2 -std=c++0x -mkl -fno-rtti" LINK_FLAGS "-mkl -openmp")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        set_target_properties (${_TARGET} PROPERTIES COMPILE_FLAGS "-Wall -Wvla -std=c++0x -fno-rtti" LINK_FLAGS "")
    elseif (CMAKE_COMPILER_IS_GNUCXX OR CXX_COMPILER_NAME STREQUAL "g++" OR CXX_COMPILER_NAME STREQUAL "c++" OR CXX_COMPILER_NAME STREQUAL "g++.exe")
        set_target_properties (${_TARGET} PROPERTIES COMPILE_FLAGS "-fopenmp -Wall -Wvla -std=c++0x -fmessage-length=0 -fno-rtti" LINK_FLAGS "-fopenmp")
    elseif (CXX_COMPILER_NAME STREQUAL "cl")
        set_target_properties (${_TARGET} PROPERTIES COMPILE_FLAGS "/openmp")
    endif ()
endforeach ()

# ------------------------------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------------------------------------------
#   Set installation related stuff
# ------------------------------------------------------------------------------------------------------------------------

# Set install directory to current directory
#set (CMAKE_INSTALL_PREFIX "$ENV{HOME}/opt/cadet" CACHE
#  STRING "Install path prefix, prepended onto install directories." FORCE)

# ------------------------------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------------------------------------------
#   Misc stuff
# ------------------------------------------------------------------------------------------------------------------------
# Make HDF5_DIR variable invisible
mark_as_advanced(FORCE HDF5_DIR)
# ------------------------------------------------------------------------------------------------------------------------

