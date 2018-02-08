# =============================================================================
#  CADET - The Chromatography Analysis and Design Toolkit
#  
#  Copyright © 2008-2018: The CADET Authors
#            Please see the AUTHORS and CONTRIBUTORS file.
#  
#  All rights reserved. This program and the accompanying materials
#  are made available under the terms of the GNU Public License v3.0 (or, at
#  your option, any later version) which accompanies this distribution, and
#  is available at http://www.gnu.org/licenses/gpl.html
# =============================================================================

# Name of the current project
project (CadetCLI CXX C)

# Add the executable CADET-CLI
add_executable (cadet-cli ${CMAKE_SOURCE_DIR}/ThirdParty/pugixml/pugixml.cpp
    ${CMAKE_SOURCE_DIR}/src/cadet-cli/cadet-cli.cpp
    ${CMAKE_SOURCE_DIR}/src/io/JsonParameterProvider.cpp)

# ---------------------------------------------------
#   Linking to LIBCADET and add dependencies
# ---------------------------------------------------

# Prefer static link
if (BUILD_STATIC_LIBS OR STANDALONE)
	target_link_libraries(cadet-cli PRIVATE libcadet_static)
else ()
	target_link_libraries(cadet-cli PRIVATE libcadet_shared)
endif ()
# ---------------------------------------------------
  
# Add include directories for access to exported LIBCADET header files.
target_include_directories(cadet-cli PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/ThirdParty/pugixml ${CMAKE_SOURCE_DIR}/ThirdParty/json ${CMAKE_SOURCE_DIR}/ThirdParty/tclap/include ${CMAKE_BINARY_DIR})

cadet_target_compile_features(cadet-cli)


# ---------------------------------------------------
#   Linking to external libraries
# ---------------------------------------------------

# Link to HDF5
target_include_directories(cadet-cli PRIVATE ${HDF5_INCLUDE_DIRS})
target_compile_definitions(cadet-cli PRIVATE ${HDF5_DEFINITIONS})
target_link_libraries(cadet-cli PRIVATE ${HDF5_LIBRARIES})

# Link to TBB for timer
if (BENCHMARK_MODE AND CADET_PARALLELIZE)
	target_include_directories(cadet-cli PRIVATE ${TBB_INCLUDE_DIRS})
	target_link_libraries(cadet-cli PRIVATE ${TBB_LIBRARIES})
endif()

# ---------------------------------------------------

# ---------------------------------------------------
#   Setup installation
# ---------------------------------------------------

# Install the cadet-cli executable
install (CODE "MESSAGE(\"\nInstall CADET-CLI\n\")")
install (TARGETS cadet-cli RUNTIME DESTINATION bin)

if (BUILD_CADET_MEX)
	# Also install into matlab/bin
	install (TARGETS cadet-cli RUNTIME DESTINATION matlab/bin)
endif()

# ---------------------------------------------------

set (CADET_CLI_TARGETS cadet-cli PARENT_SCOPE)

# Info message
message (STATUS "Added CADET-CLI module")
