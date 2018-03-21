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
project (CadetBuildTools CXX C)

# Add the template code generator
add_executable(templateCodeGen ${CMAKE_SOURCE_DIR}/src/build-tools/templateCodeGen.cpp)
target_include_directories(templateCodeGen PRIVATE ${CMAKE_SOURCE_DIR}/ThirdParty/json ${CMAKE_SOURCE_DIR}/ThirdParty/inja)
cadet_target_compile_features(templateCodeGen)

# Info message
message (STATUS "Added build tools")

