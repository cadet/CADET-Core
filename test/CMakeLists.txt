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
project(CadetTest CXX C)

set(TEST_HDF5_TARGETS )
set(TEST_NONLINALG_TARGETS )
set(TEST_LIBCADET_TARGETS )

# Add the executables
if (LAPACK_FOUND)
	add_executable(testAdaptiveTRNewton testAdaptiveTRNewton.cpp)
	target_include_directories(testAdaptiveTRNewton PRIVATE ${CMAKE_SOURCE_DIR}/src/libcadet ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/ThirdParty/tclap/include)
	target_link_libraries(testAdaptiveTRNewton PRIVATE CADET::CompileOptions libcadet_nonlinalg_static ${LAPACK_LIBRARIES})

	add_executable(testSMANonlinearSolve testSMANonlinearSolve.cpp)
	target_include_directories(testSMANonlinearSolve PRIVATE ${CMAKE_SOURCE_DIR}/src/libcadet ${CMAKE_SOURCE_DIR}/include ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/ThirdParty/tclap/include)
	target_link_libraries(testSMANonlinearSolve PRIVATE CADET::CompileOptions libcadet_nonlinalg_static ${LAPACK_LIBRARIES})
endif()

add_executable(testLogging testLogging.cpp)
target_link_libraries(testLogging PRIVATE CADET::CompileOptions)

if (NOT WIN32)
	add_executable(testRadialKernel testRadialKernel.cpp TimeIntegrator.cpp ${CMAKE_SOURCE_DIR}/src/libcadet/Logging.cpp ${CMAKE_SOURCE_DIR}/src/libcadet/AutoDiff.cpp ${CMAKE_SOURCE_DIR}/src/libcadet/model/paramdep/ParameterDependenceBase.cpp)
	target_link_libraries(testRadialKernel PRIVATE CADET::CompileOptions CADET::AD SUNDIALS::sundials_idas ${SUNDIALS_NVEC_TARGET} ${EIGEN_TARGET} )
	target_include_directories(testRadialKernel PRIVATE ${CMAKE_SOURCE_DIR}/test ${CMAKE_BINARY_DIR}/src/libcadet)

	if (ENABLE_GRM_2D)
		target_include_directories(testRadialKernel PRIVATE ${CMAKE_BINARY_DIR}/src/libcadet)
		if (SUPERLU_FOUND)
			target_link_libraries(testRadialKernel PRIVATE SuperLU::SuperLU)
		endif()
		if (UMFPACK_FOUND)
			target_link_libraries(testRadialKernel PRIVATE UMFPACK::UMFPACK)
		endif()
	endif()

	list(APPEND TEST_NONLINALG_TARGETS testRadialKernel)
	list(APPEND TEST_HDF5_TARGETS testRadialKernel)
endif()


# CATCH unit tests
set(LIBCADET_SHARED_LIBFILE "$<TARGET_FILE:libcadet_shared>")
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/Paths.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/Paths.cpp.in" @ONLY)
file(GENERATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/Paths_$<CONFIG>.cpp" INPUT "${CMAKE_CURRENT_BINARY_DIR}/Paths.cpp.in")

set(TEST_ADDITIONAL_SOURCES "")
if (ENABLE_2D_MODELS)
	list(APPEND TEST_ADDITIONAL_SOURCES SparseFactorizableMatrix.cpp TwoDimConvectionDispersionOperator.cpp)
endif()

if (ENABLE_DG)
	list(APPEND TEST_ADDITIONAL_SOURCES GeneralRateModelDG.cpp LumpedRateModelWithPoresDG.cpp LumpedRateModelWithoutPoresDG.cpp)
	if (ENABLE_2D_MODELS)
		list(APPEND TEST_ADDITIONAL_SOURCES LumpedRateModelWithPoresDG2D.cpp)
	endif()
endif()

add_executable(testRunner testRunner.cpp JsonTestModels.cpp ColumnTests.cpp UnitOperationTests.cpp SimHelper.cpp ParticleHelper.cpp
	GeneralRateModel.cpp GeneralRateModel2D.cpp LumpedRateModelWithPores.cpp LumpedRateModelWithoutPores.cpp
	RadialGeneralRateModel.cpp RadialLumpedRateModelWithPores.cpp RadialLumpedRateModelWithoutPores.cpp
	MultiChannelTransportModel.cpp
	Crystallization.cpp
	CSTR-Residual.cpp CSTR-Simulation.cpp
	ConvectionDispersionOperator.cpp
	CellKernelTests.cpp
	BindingModelTests.cpp BindingModels.cpp BindingModelAutoJacobian.cpp
	ReactionModelTests.cpp ReactionModels.cpp
	ParamDepTests.cpp ParameterDependencies.cpp
	ModelSystem.cpp
	BandMatrix.cpp DenseMatrix.cpp SparseMatrix.cpp StringHashing.cpp LogUtils.cpp AD.cpp Subset.cpp Graph.cpp
	"${CMAKE_CURRENT_BINARY_DIR}/Paths_$<CONFIG>.cpp" "${CMAKE_SOURCE_DIR}/src/io/JsonParameterProvider.cpp"
	${TEST_ADDITIONAL_SOURCES}
	$<TARGET_OBJECTS:libcadet_object>)

target_include_directories(testRunner PRIVATE ${CMAKE_BINARY_DIR}/src/libcadet)
target_link_libraries(testRunner PRIVATE CADET::CompileOptions CADET::AD SUNDIALS::sundials_idas ${SUNDIALS_NVEC_TARGET} ${TBB_TARGET} ${EIGEN_TARGET})

if (ENABLE_2D_MODELS)
	if (SUPERLU_FOUND)
		target_link_libraries(testRunner PRIVATE SuperLU::SuperLU)
	endif()
	if (UMFPACK_FOUND)
		target_link_libraries(testRunner PRIVATE UMFPACK::UMFPACK)
	endif()
endif()

list(APPEND TEST_LIBCADET_TARGETS testRunner)
list(APPEND TEST_NONLINALG_TARGETS testRunner)
list(APPEND TEST_HDF5_TARGETS testRunner)

list(APPEND TEST_TARGETS ${TEST_NONLINALG_TARGETS} ${TEST_LIBCADET_TARGETS} ${TEST_HDF5_TARGETS} testLogging)

foreach(_TARGET IN LISTS TEST_TARGETS)
	# Add include directories for access to exported LIBCADET header files.
	target_include_directories(${_TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src/libcadet ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}/src/libcadet)
	# Add include directories for third party components
	target_include_directories(${_TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/ThirdParty/Catch ${CMAKE_SOURCE_DIR}/ThirdParty/json ${CMAKE_SOURCE_DIR}/ThirdParty/pugixml ${CMAKE_SOURCE_DIR}/ThirdParty/tclap/include)
endforeach()

# ---------------------------------------------------
#   Linking to external libraries
# ---------------------------------------------------

# Link to HDF5
foreach(_TARGET IN LISTS TEST_HDF5_TARGETS)
	target_link_libraries(${_TARGET} PRIVATE HDF5::HDF5)
endforeach()

# Link to nonlinalg lib
foreach(_TARGET IN LISTS TEST_NONLINALG_TARGETS)
	target_link_libraries(${_TARGET} PRIVATE libcadet_nonlinalg_static)
endforeach()

# ---------------------------------------------------

add_executable(testCAPIv1 "${CMAKE_CURRENT_BINARY_DIR}/Paths_$<CONFIG>.cpp" testCAPIv1.cpp)
target_include_directories(testCAPIv1 PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/ThirdParty/json)

# Set C++ standard for the target
target_compile_features(testCAPIv1 PRIVATE cxx_std_23)

# ---------------------------------------------------

unset(TEST_HDF5_TARGETS)
unset(TEST_NONLINALG_TARGETS)
unset(TEST_TARGETS)
unset(TEST_ADDITIONAL_SOURCES)

# Info message
message(STATUS "Added tests")

