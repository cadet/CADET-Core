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
project (CadetMex CXX C)


# ---------------------------------------------------
#   Create mex file 
# ---------------------------------------------------
matlab_add_mex(NAME cadet-mex SRC cadet-mex.cpp MatlabCommands.cpp MatlabExternalFunction.cpp OUTPUT_NAME CadetMex LINK_TO libcadet_mex ${Matlab_LAPACK_LIBRARY} ${Matlab_BLAS_LIBRARY})

# Add include directories for access to exported LIBCADET header files.
target_include_directories (cadet-mex PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src/cadet-mex ${CMAKE_BINARY_DIR})

target_compile_features(cadet-mex PRIVATE cxx_alias_templates cxx_defaulted_functions cxx_delegating_constructors
        cxx_explicit_conversions cxx_generalized_initializers cxx_inheriting_constructors cxx_rvalue_references
        cxx_lambdas cxx_nullptr cxx_auto_type cxx_range_for cxx_right_angle_brackets cxx_deleted_functions
        cxx_strong_enums cxx_uniform_initialization cxx_template_template_parameters cxx_nullptr)

target_compile_definitions(cadet-mex PRIVATE -DCADET_MATLABMEX -DMATLAB_MEX_FILE)

# Check if we are running on a 32 or 64 bit system
if (CMAKE_SIZEOF_VOID_P EQUAL 8)
	# 64 bit
	target_compile_definitions(cadet-mex PRIVATE -DCADET_MEX_64BIT)
else()
	# 32 bit
	target_compile_definitions(cadet-mex PRIVATE -DCADET_MEX_32BIT)
endif()

if (Matlab_DECLARES_CREATEUNINITNUMERICMATRIX)
	target_compile_definitions(cadet-mex PRIVATE -DMATLAB_HAVE_CREATEUNINITNUMERICMATRIX)
endif()

if (Matlab_DECLARES_CREATEUNINITNUMERICARRAY)
	target_compile_definitions(cadet-mex PRIVATE -DMATLAB_HAVE_CREATEUNINITNUMERICARRAY)
endif()

# ---------------------------------------------------

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  if (DEBUG_OPENMP AND OPENMP_FOUND)
	  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	      target_compile_options(cadet-mex PRIVATE ${OpenMP_CXX_FLAGS} -w2 -std=c++0x -fno-rtti)
	  elseif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang"))
	      set_target_properties (cadet-mex PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
	      target_compile_options(cadet-mex PRIVATE ${OpenMP_CXX_FLAGS} -Wall -Wvla -pedantic -fno-rtti)
	  elseif (CMAKE_COMPILER_IS_GNUCXX OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
	      set_target_properties (cadet-mex PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
	      target_compile_options(cadet-mex PRIVATE ${OpenMP_CXX_FLAGS} -Wall -Wvla -pedantic -fno-rtti)
	  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	      # GR-= Disable RTTI
	      target_compile_options(cadet-mex PRIVATE "/GR-" ${OpenMP_CXX_FLAGS})
	  endif()
  else()
	  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	      target_compile_options(cadet-mex PRIVATE -w2 -std=c++0x -fno-rtti)
	  elseif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang"))
	      target_compile_options(cadet-mex PRIVATE -Wall -Wvla -pedantic -fno-rtti)
	  elseif (CMAKE_COMPILER_IS_GNUCXX OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
	      target_compile_options(cadet-mex PRIVATE -Wall -Wvla -pedantic -fno-rtti)
	  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	      # GR-= Disable RTTI
	      target_compile_options(cadet-mex PRIVATE "/GR-")
	  endif()
  endif()
else()
  if (OPENMP_FOUND)
	  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	      target_compile_options(cadet-mex PRIVATE ${OpenMP_CXX_FLAGS} -w2 -std=c++0x -fno-rtti)
	  elseif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang"))
	      set_target_properties (cadet-mex PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
	      target_compile_options(cadet-mex PRIVATE ${OpenMP_CXX_FLAGS} -Wall -Wvla -pedantic -fno-rtti)
	  elseif (CMAKE_COMPILER_IS_GNUCXX OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
	      set_target_properties (cadet-mex PROPERTIES LINK_FLAGS "${OpenMP_CXX_FLAGS}")
	      target_compile_options(cadet-mex PRIVATE ${OpenMP_CXX_FLAGS} -Wall -Wvla -pedantic -std=c++0x -fno-rtti)
	  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	      # GR-= Disable RTTI, LTCG = Link time code generation, /GL = Whole program optimization
	      set_target_properties (cadet-mex PROPERTIES LINK_FLAGS "/LTCG")
	      target_compile_options(cadet-mex PRIVATE "/GL" "/GR-" ${OpenMP_CXX_FLAGS})
	  endif ()
  else()
	  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	      target_compile_options(cadet-mex PRIVATE -w2 -std=c++0x -fno-rtti)
	  elseif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang"))
	      target_compile_options(cadet-mex PRIVATE -Wall -Wvla -pedantic -fno-rtti)
	  elseif (CMAKE_COMPILER_IS_GNUCXX OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
	      target_compile_options(cadet-mex PRIVATE -Wall -Wvla -pedantic -std=c++0x -fno-rtti)
	  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
	      # GR-= Disable RTTI, LTCG = Link time code generation, /GL = Whole program optimization
	      set_target_properties (cadet-mex PROPERTIES LINK_FLAGS "/LTCG")
	      target_compile_options(cadet-mex PRIVATE "/GL" "/GR-")
	  endif ()
  endif()
endif()

# ---------------------------------------------------
#   Setup installation
# ---------------------------------------------------

# Install the cadet-cs executable
install (CODE "MESSAGE(\"\nInstall CADET-MEX\n\")")

# ---------------------------------------------------

# Info message
message (STATUS "Added CADET-MEX module")
