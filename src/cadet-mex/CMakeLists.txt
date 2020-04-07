# =============================================================================
#  CADET - The Chromatography Analysis and Design Toolkit
#  
#  Copyright © 2008-2020: The CADET Authors
#            Please see the AUTHORS and CONTRIBUTORS file.
#  
#  All rights reserved. This program and the accompanying materials
#  are made available under the terms of the GNU Public License v3.0 (or, at
#  your option, any later version) which accompanies this distribution, and
#  is available at http://www.gnu.org/licenses/gpl.html
# =============================================================================

# Name of the current project
project(CadetMex CXX C)


# ---------------------------------------------------
#   Create mex file 
# ---------------------------------------------------
matlab_add_mex(NAME cadet-mex SRC cadet-mex.cpp MatlabCommands.cpp MatlabExternalFunction.cpp OUTPUT_NAME CadetMex LINK_TO libcadet_mex)

# Add include directories for access to exported LIBCADET header files.
target_include_directories(cadet-mex PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src/cadet-mex ${CMAKE_BINARY_DIR})
target_compile_definitions(cadet-mex PRIVATE CADET_MATLABMEX MATLAB_MEX_FILE)

# Check if we are running on a 32 or 64 bit system
if (CMAKE_SIZEOF_VOID_P EQUAL 8)
	# 64 bit
	target_compile_definitions(cadet-mex PRIVATE CADET_MEX_64BIT)
else()
	# 32 bit
	target_compile_definitions(cadet-mex PRIVATE CADET_MEX_32BIT)
endif()

if (Matlab_DECLARES_CREATEUNINITNUMERICMATRIX)
	target_compile_definitions(cadet-mex PRIVATE MATLAB_HAVE_CREATEUNINITNUMERICMATRIX)
endif()

if (Matlab_DECLARES_CREATEUNINITNUMERICARRAY)
	target_compile_definitions(cadet-mex PRIVATE MATLAB_HAVE_CREATEUNINITNUMERICARRAY)
endif()

# ---------------------------------------------------
#   Setup installation
# ---------------------------------------------------

# Install the cadet-mex executable
install(CODE "MESSAGE(\"\nInstall CADET-MEX\n\")")
install(DIRECTORY ${CMAKE_SOURCE_DIR}/src/cadet-mi/ DESTINATION matlab FILES_MATCHING REGEX ".*\.(m|md)$")
install(FILES ${CMAKE_SOURCE_DIR}/AUTHORS ${CMAKE_SOURCE_DIR}/CONTRIBUTORS ${CMAKE_SOURCE_DIR}/LICENSE.txt DESTINATION matlab)
install(TARGETS cadet-mex RUNTIME DESTINATION matlab/bin LIBRARY DESTINATION matlab/bin)

# ---------------------------------------------------

# Info message
message(STATUS "Added CADET-MEX module")
