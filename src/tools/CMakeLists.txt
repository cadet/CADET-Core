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
project(CadetTools CXX C)

set(TOOLS_TARGETS)

# Wheel build
option(CADET_WHEEL "Option for building wheels" OFF)

# Add the executables
if (HDF5_FOUND)

	add_executable(createLWE createLWE.cpp)
	list(APPEND TOOLS_TARGETS createLWE)

	add_executable(convertFile convertFile.cpp ${CMAKE_SOURCE_DIR}/src/io/FileIO.cpp FormatConverter.cpp ${CMAKE_SOURCE_DIR}/ThirdParty/pugixml/pugixml.cpp)
	target_include_directories(convertFile PRIVATE ${CMAKE_SOURCE_DIR}/ThirdParty/pugixml ${CMAKE_SOURCE_DIR}/ThirdParty/json)
	list(APPEND TOOLS_TARGETS convertFile)
endif()

foreach(_TARGET IN LISTS TOOLS_TARGETS)
	# Add include directories for access to exported LIBCADET header files.
	target_include_directories(${_TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR})
	# Add include directories for third party components
	target_include_directories(${_TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/ThirdParty/tclap/include)

	target_link_libraries(${_TARGET} PRIVATE CADET::CompileOptions HDF5::HDF5)
endforeach()

# ---------------------------------------------------

if (CADET_WHEEL)
  install(CODE "MESSAGE(\"\nInstall CADET Tools for wheels\n\")")

  if (APPLE)
    set_target_properties(${TOOLS_TARGETS} PROPERTIES
      BUILD_WITH_INSTALL_RPATH ON
      INSTALL_RPATH "@loader_path/../lib;@loader_path/../.dylibs"
    )
  elseif (UNIX)  # UNIX AND NOT APPLE
    set_target_properties(${TOOLS_TARGETS} PROPERTIES
      BUILD_WITH_INSTALL_RPATH ON
      INSTALL_RPATH "$ORIGIN/../lib"
    )
  endif()

  install(TARGETS ${TOOLS_TARGETS} RUNTIME DESTINATION cadet_core/bin)
else()
	install(CODE "MESSAGE(\"\nInstall CADET Tools\n\")")
	install(TARGETS ${TOOLS_TARGETS} RUNTIME)
endif()

# ---------------------------------------------------

unset(TOOLS_TARGETS)

# Info message
message(STATUS "Added tools")

