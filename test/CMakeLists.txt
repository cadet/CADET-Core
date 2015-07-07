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
project (${TESTS_PRJ})

# Set variables defining include directories
set (LIBCADET_INC_DIR "${CMAKE_SOURCE_DIR}/include/libcadet")
set (COMMON_INC_DIR "${CMAKE_SOURCE_DIR}/include/common")

# Add include directories for access to exported header files.
include_directories (${LIBCADET_INC_DIR} ${COMMON_INC_DIR} "${CMAKE_SOURCE_DIR}/ThirdParty/tclap/include" "${CMAKE_SOURCE_DIR}/ThirdParty/ADOL-C/include" "${CMAKE_SOURCE_DIR}/src/libcadet")

# Find HDF5
if (NOT WIN32)
    find_package (HDF5 REQUIRED CXX)
    include_directories (${HDF5_INCLUDE_DIRS} "${HDF5_INCLUDE_DIRS}/cpp")
    add_definitions (${HDF5_DEFINITIONS})
else () # On WINDOWS environment var HDF5_ROOT must be set!
    set (HDF5_DIR "$ENV{HDF5_ROOT}/cmake/hdf5")
    find_package(HDF5 NO_MODULE)
    include_directories (${HDF5_INCLUDE_DIR} ${HDF5_INCLUDE_DIR_CPP})
endif ()

# ------------------------------------------------------------------------------------------------------------------------
#   Build result comparer
# ------------------------------------------------------------------------------------------------------------------------

add_executable (${RESULTCOMP_NAME} rescomp.cpp)
target_link_libraries (${RESULTCOMP_NAME} ${HDF5_LIBRARIES})

# ------------------------------------------------------------------------------------------------------------------------
#   Build regression test creators
# ------------------------------------------------------------------------------------------------------------------------

add_executable (createRegTestC14 createRegTestC14.cpp)
target_link_libraries (createRegTestC14 ${HDF5_LIBRARIES})

add_executable (createRegTestLWE createRegTestLWE.cpp)
target_link_libraries (createRegTestLWE ${HDF5_LIBRARIES})

add_executable (createRegTestSCL createRegTestSCL.cpp)
target_link_libraries (createRegTestSCL ${HDF5_LIBRARIES})

# ------------------------------------------------------------------------------------------------------------------------
#   Build SFAD tests
# ------------------------------------------------------------------------------------------------------------------------

add_executable (compareSFADadolC compareSFADadolC.cpp)
target_link_libraries (compareSFADadolC)

add_executable (compareSETFADadolC compareSETFADadolC.cpp)
target_link_libraries (compareSETFADadolC)

if (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel") OR (CXX_COMPILER_NAME STREQUAL "icpc") OR (CXX_COMPILER_NAME STREQUAL "icc"))
  set_target_properties (compareSFADadolC PROPERTIES COMPILE_FLAGS "-w2 -std=c++0x")
  set_target_properties (compareSETFADadolC PROPERTIES COMPILE_FLAGS "-w2 -std=c++0x")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set_target_properties (compareSFADadolC PROPERTIES COMPILE_FLAGS "-Wall -Wvla -std=c++0x -fno-rtti")
  set_target_properties (compareSETFADadolC PROPERTIES COMPILE_FLAGS "-Wall -Wvla -std=c++0x -fno-rtti")
elseif (CMAKE_COMPILER_IS_GNUCXX OR CXX_COMPILER_NAME STREQUAL "g++" OR CXX_COMPILER_NAME STREQUAL "c++" OR CXX_COMPILER_NAME STREQUAL "g++.exe")
  set_target_properties (compareSFADadolC PROPERTIES COMPILE_FLAGS "-Wall -Wvla -std=c++0x")
  set_target_properties (compareSETFADadolC PROPERTIES COMPILE_FLAGS "-Wall -Wvla -std=c++0x")
endif ()
