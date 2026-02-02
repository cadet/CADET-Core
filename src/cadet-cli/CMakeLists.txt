# =============================================================================
#  CADET
#  
#  Copyright © 2008-present: The CADET-Core Authors
#            Please see the AUTHORS.md file.
#  
#  All rights reserved. This program and the accompanying materials
#  are made available under the terms of the GNU Public License v3.0 (or, at
#  your option, any later version) which accompanies this distribution, and
#  is available at http://www.gnu.org/licenses/gpl.html
# =============================================================================

# Name of the current project
project(CadetCLI CXX C)

# Wheel build
option(CADET_WHEEL "Option for building wheels" OFF)

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

if (CADET_WHEEL)
  install(CODE "MESSAGE(\"\nInstall CADET-CLI for wheels\n\")")

  if (APPLE)
    set_target_properties(cadet-cli PROPERTIES
      BUILD_WITH_INSTALL_RPATH ON
      INSTALL_RPATH "@loader_path/../lib;@loader_path/../.dylibs"
    )
  elseif (UNIX)  # UNIX AND NOT APPLE
    set_target_properties(cadet-cli PROPERTIES
      BUILD_WITH_INSTALL_RPATH ON
      INSTALL_RPATH "$ORIGIN/../lib"
    )
  endif()

  install(TARGETS cadet-cli RUNTIME DESTINATION cadet_core/bin)
else()
  install(CODE "MESSAGE(\"\nInstall CADET-CLI\n\")")
  install(TARGETS cadet-cli RUNTIME)
endif()



# ---------------------------------------------------

# Info message
message(STATUS "Added CADET-CLI module")
