# =============================================================================
#  CADET - The Chromatography Analysis and Design Toolkit
#  
#  Copyright © 2008-2016: The CADET Authors
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
    ${CMAKE_SOURCE_DIR}/src/cadet-cli/cadet-cli.cpp)

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
target_include_directories (cadet-cli PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/ThirdParty/pugixml ${CMAKE_SOURCE_DIR}/ThirdParty/tclap/include ${CMAKE_BINARY_DIR})

target_compile_features(cadet-cli PRIVATE cxx_alias_templates cxx_defaulted_functions cxx_delegating_constructors
        cxx_explicit_conversions cxx_generalized_initializers cxx_inheriting_constructors cxx_rvalue_references
        cxx_lambdas cxx_nullptr cxx_auto_type cxx_range_for cxx_right_angle_brackets cxx_deleted_functions cxx_nullptr 
        cxx_strong_enums cxx_uniform_initialization cxx_template_template_parameters cxx_defaulted_move_initializers)



# ---------------------------------------------------
#   Linking to external libraries
# ---------------------------------------------------

# Link to HDF5
target_include_directories (cadet-cli PRIVATE ${HDF5_INCLUDE_DIRS})
target_compile_definitions (cadet-cli PRIVATE ${HDF5_DEFINITIONS})
target_link_libraries(cadet-cli PRIVATE ${HDF5_LIBRARIES})
# ---------------------------------------------------

# ---------------------------------------------------
#   Setup installation
# ---------------------------------------------------

# Install the cadet-cli executable
install (CODE "MESSAGE(\"\nInstall CADET-CLI\n\")")
install (TARGETS cadet-cli RUNTIME DESTINATION bin)

# ---------------------------------------------------

set (CADETCLI_TARGETS cadet-cli PARENT_SCOPE)

# Info message
message (STATUS "Added CADET-CLI module")
