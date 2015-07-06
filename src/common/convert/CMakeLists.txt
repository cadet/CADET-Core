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
project (${CONVERTER_PRJ} CXX)

# Set variables defining source and include directories
set (CONVERTER_SRC_DIR "${CMAKE_SOURCE_DIR}/src/common/convert")
set (PUGIXML_DIR "${CMAKE_SOURCE_DIR}/ThirdParty/pugixml")

# Add variable H5TOXML_SOURCES with the sources for H5TOXML
set (H5TOXML_SOURCES
     ${CONVERTER_SRC_DIR}/h5toxml.cpp
     ${PUGIXML_DIR}/pugixml.cpp
    )

# Add variable XMLTOH5_SOURCES with the sources for XMLTOH5
set (XMLTOH5_SOURCES
     ${CONVERTER_SRC_DIR}/xmltoh5.cpp
     ${PUGIXML_DIR}/pugixml.cpp
    )

# Add include directories for access to exported LIBCADET header files.
include_directories (${LIBCADET_INC_DIR} ${PUGIXML_DIR} ${CMAKE_SOURCE_DIR}/include/common)

# Add the executable H5TOXML
add_executable (${H5TOXML_NAME} ${H5TOXML_SOURCES})

# Add the executable XMLTOH5
add_executable (${XMLTOH5_NAME} ${XMLTOH5_SOURCES})



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
target_link_libraries(${H5TOXML_NAME} ${HDF5_LIBRARIES})
target_link_libraries(${XMLTOH5_NAME} ${HDF5_LIBRARIES})
# ------------------------------------------------------------------------------------------------------------------------

# Info message
message (STATUS "Added CONVERTER module")


# Reduce Warning level in VC++
if(CMAKE_BUILD_TOOL MATCHES "(msdev|devenv|nmake)")
    add_definitions(/W2)
endif()



# ------------------------------------------------------------------------------------------------------------------------
#   Set installation related stuff
# ------------------------------------------------------------------------------------------------------------------------

# Install the cadet-cs executable
install (CODE "MESSAGE(\"\nInstall CONVERTER\n\")")
install (TARGETS ${H5TOXML_NAME} ${XMLTOH5_NAME} DESTINATION bin)

# ------------------------------------------------------------------------------------------------------------------------
