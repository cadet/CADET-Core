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
project (${TESTS_PRJ})

# Set variables defining include directories
set (LIBCADET_INC_DIR "${CMAKE_SOURCE_DIR}/include/libcadet")
set (COMMON_INC_DIR "${CMAKE_SOURCE_DIR}/include/common")

# Add include directories for access to exported header files.
include_directories (${LIBCADET_INC_DIR} ${COMMON_INC_DIR} "${CMAKE_SOURCE_DIR}/ThirdParty/tclap/include")

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
