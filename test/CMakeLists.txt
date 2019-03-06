# =============================================================================
#  CADET - The Chromatography Analysis and Design Toolkit
#  
#  Copyright © 2008-2019: The CADET Authors
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
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Paths.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/Paths.cpp" @ONLY)

set(TEST_ADDITIONAL_SOURCES "")
if (ENABLE_GRM_2D)
	list(APPEND TEST_ADDITIONAL_SOURCES SparseFactorizableMatrix.cpp TwoDimConvectionDispersionOperator.cpp)
endif()

add_executable (testRunner testRunner.cpp JsonTestModels.cpp ColumnTests.cpp UnitOperationTests.cpp SimHelper.cpp ParticleHelper.cpp
	GeneralRateModel.cpp GeneralRateModel2D.cpp LumpedRateModelWithPores.cpp LumpedRateModelWithoutPores.cpp
	CSTR-Residual.cpp CSTR-Simulation.cpp 
	ConvectionDispersionOperator.cpp
	BindingModelTests.cpp BindingModels.cpp
	ModelSystem.cpp
	BandMatrix.cpp DenseMatrix.cpp SparseMatrix.cpp StringHashing.cpp LogUtils.cpp AD.cpp 
	"${CMAKE_CURRENT_BINARY_DIR}/Paths.cpp" "${CMAKE_SOURCE_DIR}/src/io/JsonParameterProvider.cpp"
	${TEST_ADDITIONAL_SOURCES}
	$<TARGET_OBJECTS:libcadet_object>)

target_link_libraries(testRunner PRIVATE ${SUNDIALS_LIBRARIES})
if (CADET_PARALLEL_FLAG)
	target_compile_definitions(testRunner PRIVATE ${CADET_PARALLEL_FLAG})
	target_include_directories(testRunner PRIVATE ${TBB_INCLUDE_DIRS})
	target_link_libraries(testRunner PRIVATE ${TBB_LIBRARIES})
endif()
if (ENABLE_GRM_2D)
	target_include_directories(testRunner PRIVATE ${CMAKE_BINARY_DIR}/src/libcadet)
	if (SUPERLU_FOUND)
		target_include_directories (testRunner PRIVATE ${SUPERLU_INCLUDE_DIRS})
		target_link_libraries(testRunner PRIVATE ${SUPERLU_LIBRARIES})
	endif()
	if (UMFPACK_FOUND)
		target_include_directories (testRunner PRIVATE ${UMFPACK_INCLUDE_DIRS})
		target_link_libraries(testRunner PRIVATE ${UMFPACK_LIBRARIES})
	endif()
endif()

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

# ---------------------------------------------------

set (TEST_TARGETS ${TEST_TARGETS} PARENT_SCOPE)

# Info message
message (STATUS "Added tests")

