# =============================================================================
#  CADET
#  
#  Copyright © 2008-2021: The CADET Authors
#            Please see the AUTHORS and CONTRIBUTORS file.
#  
#  All rights reserved. This program and the accompanying materials
#  are made available under the terms of the GNU Public License v3.0 (or, at
#  your option, any later version) which accompanies this distribution, and
#  is available at http://www.gnu.org/licenses/gpl.html
# =============================================================================

# Name of the current project
project(CadetCLI CXX C)

# Add the executable CADET-CLI
add_executable(cadet-cli
	${CMAKE_SOURCE_DIR}/ThirdParty/pugixml/pugixml.cpp
	${CMAKE_SOURCE_DIR}/src/cadet-cli/cadet-cli.cpp
	${CMAKE_SOURCE_DIR}/src/io/JsonParameterProvider.cpp
	${CMAKE_SOURCE_DIR}/src/cadet-cli/ProgressBar.cpp
	${CMAKE_SOURCE_DIR}/src/cadet-cli/SignalHandler.cpp
)

# ---------------------------------------------------
#   Linking to LIBCADET and add dependencies
# ---------------------------------------------------

if (ENABLE_STATIC_LINK_CLI)
	target_link_libraries(cadet-cli PRIVATE libcadet_static)
else()
	target_link_libraries(cadet-cli PRIVATE libcadet_shared)

	# On Windows, copy libcadet.dll file to directory of cadet-cli.exe in order to ease debugging
	if (WIN32)
		add_custom_command(
			TARGET cadet-cli
			POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy
				$<TARGET_FILE:libcadet_shared>
				$<TARGET_FILE_DIR:cadet-cli>/$<TARGET_FILE_NAME:libcadet_shared>
		)
	endif()
endif()
  
# Add include directories for access to exported LIBCADET header files.
target_include_directories(cadet-cli PRIVATE ${CMAKE_SOURCE_DIR}/ThirdParty/pugixml ${CMAKE_SOURCE_DIR}/ThirdParty/json ${CMAKE_SOURCE_DIR}/ThirdParty/tclap/include ${CMAKE_BINARY_DIR})

# Link to HDF5
target_link_libraries(cadet-cli PRIVATE HDF5::HDF5)

# Link to TBB for timer
if (ENABLE_BENCHMARK OR CADET_PARALLEL_FLAG)
	target_link_libraries(cadet-cli PRIVATE ${TBB_TARGET})
endif()

# ---------------------------------------------------
#   Setup installation
# ---------------------------------------------------

# Install the cadet-cli executable
install(CODE "MESSAGE(\"\nInstall CADET-CLI\n\")")
install(TARGETS cadet-cli RUNTIME)

if (BUILD_CADET_MEX)
	# Also install into matlab/bin
	install(TARGETS cadet-cli RUNTIME DESTINATION matlab/bin)
endif()

# ---------------------------------------------------

# Info message
message(STATUS "Added CADET-CLI module")
