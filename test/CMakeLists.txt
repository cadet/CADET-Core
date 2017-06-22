# =============================================================================
#  CADET - The Chromatography Analysis and Design Toolkit
#  
#  Copyright © 2008-2017: The CADET Authors
#            Please see the AUTHORS and CONTRIBUTORS file.
#  
#  All rights reserved. This program and the accompanying materials
#  are made available under the terms of the GNU Public License v3.0 (or, at
#  your option, any later version) which accompanies this distribution, and
#  is available at http://www.gnu.org/licenses/gpl.html
# =============================================================================

# Name of the current project
project (CadetTest CXX C)

set(TEST_HDF5_TARGETS )
set(TEST_NONLINALG_TARGETS )
set(TEST_LIBCADET_TARGETS )

# Add the executables
if (LAPACK_FOUND)
	add_executable (testAdaptiveTRNewton testAdaptiveTRNewton.cpp)
	list(APPEND TEST_NONLINALG_TARGETS testAdaptiveTRNewton)

	add_executable (testSMANonlinearSolve testSMANonlinearSolve.cpp)
	list(APPEND TEST_NONLINALG_TARGETS testSMANonlinearSolve)
endif()

add_executable (testLogging testLogging.cpp)


# CATCH unit tests
add_executable (testRunner testRunner.cpp JsonParameterProvider.cpp GRM-Residual.cpp GRM-Simulation.cpp BandMatrix.cpp DenseMatrix.cpp StringHashing.cpp AD.cpp)

list(APPEND TEST_LIBCADET_TARGETS testRunner)
list(APPEND TEST_NONLINALG_TARGETS testRunner)

cadet_choose_ad_lib(testRunner)

list(APPEND TEST_TARGETS ${TEST_NONLINALG_TARGETS} ${TEST_LIBCADET_TARGETS} ${TEST_HDF5_TARGETS} testLogging)

foreach(_TARGET IN LISTS TEST_TARGETS)
	# Add include directories for access to exported LIBCADET header files.
	target_include_directories (${_TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src/libcadet ${CMAKE_BINARY_DIR})
	# Add include directories for third party components
	target_include_directories (${_TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/ThirdParty/Catch ${CMAKE_SOURCE_DIR}/ThirdParty/json ${CMAKE_SOURCE_DIR}/ThirdParty/pugixml ${CMAKE_SOURCE_DIR}/ThirdParty/tclap/include)
	# Add SFAD and SETFAD AD libraries
	target_include_directories (${_TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/include/ad)

	cadet_target_compile_features(${_TARGET})
endforeach()

# ---------------------------------------------------
#   Linking to external libraries
# ---------------------------------------------------

# Link to HDF5
foreach(_TARGET IN LISTS TEST_HDF5_TARGETS)
	target_include_directories (${_TARGET} PRIVATE ${HDF5_INCLUDE_DIRS})
	target_compile_definitions (${_TARGET} PRIVATE ${HDF5_DEFINITIONS})
	target_link_libraries(${_TARGET} PRIVATE ${HDF5_LIBRARIES})
endforeach()

# Link to nonlinalg lib
foreach(_TARGET IN LISTS TEST_NONLINALG_TARGETS)
	target_link_libraries(${_TARGET} PRIVATE libcadet_nonlinalg_static)
endforeach()

# Link to libcadet
foreach(_TARGET IN LISTS TEST_LIBCADET_TARGETS)
	if (BUILD_SHARED_LIBS)
		target_link_libraries(${_TARGET} PRIVATE libcadet_shared)
	else ()
		target_link_libraries(${_TARGET} PRIVATE libcadet_static)
	endif ()
endforeach()
# ---------------------------------------------------

set (TEST_TARGETS ${TEST_TARGETS} PARENT_SCOPE)

# Info message
message (STATUS "Added tests")

