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
project (${CADET_CS_PRJ} CXX)

# Set variables defining source and include directories
set (CADET_CS_SRC_DIR "${CMAKE_SOURCE_DIR}/src/cadet-cs")

set (COMMON_SRC_DIR "${CMAKE_SOURCE_DIR}/src/common")
set (COMMON_INC_DIR "${CMAKE_SOURCE_DIR}/include/common")
set (PUGIXML_DIR "${CMAKE_SOURCE_DIR}/ThirdParty/pugixml")

# Add variable CADET_CS_SOURCES with the sources for CADET-CS
set (CADET_CS_SOURCES
     ${CADET_CS_SRC_DIR}/cadet-cs.cpp
    )

# Add variable CADET_CS_MICRO_SOURCES with the sources for CADET-CS-MICRO
set (CADET_CS_MICRO_SOURCES
     ${CADET_CS_SRC_DIR}/cadet-cs_micro.cpp
    )

# Add variable COMMON_SOURCES with commonly used sources
set (COMMON_SOURCES
     ${PUGIXML_DIR}/pugixml.cpp
    )
    
if (NOT WIN32)
  list (APPEND COMMON_SOURCES  ${COMMON_SRC_DIR}/gnuplot/Gnuplot.cpp)
endif ()
  
# Add include directories for access to exported LIBCADET header files.
include_directories (${LIBCADET_INC_DIR} ${COMMON_INC_DIR} ${PUGIXML_DIR})



# ------------------------------------------------------------------------------------------------------------------------
#   Linking to LIBCADET and add dependencies
# ------------------------------------------------------------------------------------------------------------------------

# Add the executable CADET-CS
add_executable (${CADET_CS_NAME} ${CADET_CS_SOURCES} ${COMMON_SOURCES})

# Add the executable CADET-CS-MICRO
add_executable (${CADET_CS_MICRO_NAME} ${CADET_CS_MICRO_SOURCES} ${COMMON_SOURCES})

# Add dependency on the LIBCADET library
foreach (_TARGET ${LIBCADET_TARGETS})
    add_dependencies(${CADET_CS_NAME} ${_TARGET})
    add_dependencies(${CADET_CS_MICRO_NAME} ${_TARGET})
endforeach ()

if (BUILD_SHARED_LIBS)
    target_link_libraries(${CADET_CS_NAME} ${LIBCADET_NAME_SHARED})
    target_link_libraries(${CADET_CS_MICRO_NAME} ${LIBCADET_NAME_SHARED})
else ()
    target_link_libraries(${CADET_CS_NAME} ${LIBCADET_NAME_STATIC})
    target_link_libraries(${CADET_CS_MICRO_NAME} ${LIBCADET_NAME_STATIC})
endif ()
# ------------------------------------------------------------------------------------------------------------------------



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
target_link_libraries(${CADET_CS_NAME} ${HDF5_LIBRARIES})
# ------------------------------------------------------------------------------------------------------------------------

# Info message
message (STATUS "Added CADET-CS module")


# Reduce Warning level in VC++ and prevent WinDef.h from defining min and max macros
if(CMAKE_BUILD_TOOL MATCHES "(msdev|devenv|nmake)")
    add_definitions(/W2 /DNOMINMAX)
endif()


# ------------------------------------------------------------------------------------------------------------------------
#   Set installation related stuff
# ------------------------------------------------------------------------------------------------------------------------

# Install the cadet-cs executable
install (CODE "MESSAGE(\"\nInstall CADET-CS\n\")")
install (TARGETS ${CADET_CS_NAME} DESTINATION bin)
install (TARGETS ${CADET_CS_MICRO_NAME} DESTINATION bin)

# ------------------------------------------------------------------------------------------------------------------------
