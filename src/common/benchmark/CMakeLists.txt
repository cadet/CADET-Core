# =============================================================================
#  CADET - The Chromatography Analysis and Design Toolkit
#  
#  Copyright © 2008-2015: Eric von Lieres¹, Joel Andersson,
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
project (${BENCHMARK_PRJ} CXX)

# Set variables defining source and include directories
set (BENCHMARK_SRC_DIR "${CMAKE_SOURCE_DIR}/src/common/benchmark")

set (COMMON_INC_DIR "${CMAKE_SOURCE_DIR}/include/common")
set (COMMON_SRC_DIR "${CMAKE_SOURCE_DIR}/src/common")

# Add variable SCL_SOURCES with the sources for SCL
set (SCL_SOURCES
     ${BENCHMARK_SRC_DIR}/createSCL.cpp
    )

# Add variable LWE_SOURCES with the sources for LWE
set (LWE_SOURCES
     ${BENCHMARK_SRC_DIR}/createLWE.cpp
    )

# Add variable C14_SOURCES with the sources for C14
set (C14_SOURCES
     ${BENCHMARK_SRC_DIR}/createC14.cpp
    )

# Add variable EXTL_SOURCES with the sources for EXTL
set (EXTL_SOURCES
     ${BENCHMARK_SRC_DIR}/createEXTL.cpp
    )


# Add include directories for access to exported LIBCADET header files.
include_directories (${LIBCADET_INC_DIR} ${COMMON_INC_DIR})

# Add the executable SCL
add_executable (${SCL_NAME} ${SCL_SOURCES})

# Add the executable LWE
add_executable (${LWE_NAME} ${LWE_SOURCES})

# Add the executable C14
add_executable (${C14_NAME} ${C14_SOURCES})

# Add the executable EXTL
add_executable (${EXTL_NAME} ${EXTL_SOURCES})


# ------------------------------------------------------------------------------------------------------------------------
#   Linking to external libraries
# ------------------------------------------------------------------------------------------------------------------------

# Link to HDF5
if (NOT WIN32)
    find_package (HDF5 REQUIRED CXX)
    include_directories (${HDF5_INCLUDE_DIRS} "${HDF5_INCLUDE_DIRS}/cpp")
    add_definitions (${HDF5_DEFINITIONS})
else () # On WINDOWS environment var HDF5_ROOT must be set!
    set (HDF5_DIR "$ENV{HDF5_ROOT}/cmake/hdf5")
    find_package(HDF5 NO_MODULE)
    include_directories (${HDF5_INCLUDE_DIR} ${HDF5_INCLUDE_DIR_CPP})
endif ()
target_link_libraries(${SCL_NAME}  ${HDF5_LIBRARIES})
target_link_libraries(${LWE_NAME}  ${HDF5_LIBRARIES})
target_link_libraries(${C14_NAME}  ${HDF5_LIBRARIES})
target_link_libraries(${EXTL_NAME} ${HDF5_LIBRARIES})
# ------------------------------------------------------------------------------------------------------------------------

# Info message
message (STATUS "Added BENCHMARK module")


# Reduce Warning level in VC++
if(CMAKE_BUILD_TOOL MATCHES "(msdev|devenv|nmake)")
    add_definitions(/W2)
endif()


# ------------------------------------------------------------------------------------------------------------------------
#   Set installation related stuff
# ------------------------------------------------------------------------------------------------------------------------

# Install the cadet-cs executable
install (CODE "MESSAGE(\"\nInstall BENCHMARK\n\")")
install (TARGETS ${SCL_NAME} ${LWE_NAME} ${C14_NAME} ${THM_NAME} DESTINATION bin)

# ------------------------------------------------------------------------------------------------------------------------
