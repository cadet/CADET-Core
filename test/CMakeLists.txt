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
project (CadetTest CXX C)

set(TEST_HDF5_TARGETS )
set(TEST_NONLINALG_TARGETS )
set(TEST_LIBCADET_TARGETS )

# Add the executables
if (HDF5_FOUND)
    add_executable (createLWE createLWE.cpp)
    list(APPEND TEST_HDF5_TARGETS createLWE)

    add_executable (createSCLin createSCLin.cpp)
    list(APPEND TEST_HDF5_TARGETS createSCLin)

    add_executable (createSCLinStep createSCLinStep.cpp)
    list(APPEND TEST_HDF5_TARGETS createSCLinStep)

    add_executable (createMCLin createMCLin.cpp)
    list(APPEND TEST_HDF5_TARGETS createMCLin)

    add_executable (createSCLang createSCLang.cpp)
    list(APPEND TEST_HDF5_TARGETS createSCLang)

    add_executable (createConvBenchmark createConvBenchmark.cpp)
    list(APPEND TEST_HDF5_TARGETS createConvBenchmark)
endif()

if (LAPACK_FOUND)
    add_executable (bandMatrixOutput bandMatrixOutput.cpp)
    list(APPEND TEST_NONLINALG_TARGETS bandMatrixOutput)

    add_executable (testAdaptiveTRNewton testAdaptiveTRNewton.cpp)
    list(APPEND TEST_NONLINALG_TARGETS testAdaptiveTRNewton)

    add_executable (testSMANonlinearSolve testSMANonlinearSolve.cpp)
    list(APPEND TEST_NONLINALG_TARGETS testSMANonlinearSolve)

    add_executable (testDenseFactorizations testDenseFactorizations.cpp)
    list(APPEND TEST_NONLINALG_TARGETS testDenseFactorizations)

    add_executable (testDenseBandedCopy testDenseBandedCopy.cpp)
    list(APPEND TEST_NONLINALG_TARGETS testDenseBandedCopy)

    add_executable (testBandSubmatrixMultiply testBandSubmatrixMultiply.cpp)
    list(APPEND TEST_NONLINALG_TARGETS testBandSubmatrixMultiply)

    add_executable (testDenseSubmatrixFromAD testDenseSubmatrixFromAD.cpp)
    list(APPEND TEST_NONLINALG_TARGETS testDenseSubmatrixFromAD)
    list(APPEND TEST_LIBCADET_TARGETS testDenseSubmatrixFromAD)
endif()

add_executable (testRowColIndexConverter testRowColIndexConverter.cpp)
add_executable (testLogging testLogging.cpp)

list(APPEND TEST_TARGETS ${TEST_NONLINALG_TARGETS} ${TEST_LIBCADET_TARGETS} ${TEST_HDF5_TARGETS} testRowColIndexConverter testLogging)

foreach(_TARGET IN LISTS TEST_TARGETS)
    # Add include directories for access to exported LIBCADET header files.
    target_include_directories (${_TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src/libcadet ${CMAKE_BINARY_DIR})
    #${CMAKE_BINARY_DIR}/src/libcadet 
    # Add include directories for third party components
    target_include_directories (${_TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/ThirdParty/pugixml ${CMAKE_SOURCE_DIR}/ThirdParty/tclap/include)
    # Add SFAD and SETFAD AD libraries
    target_include_directories (${_TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/include/ad)

    target_compile_features(${_TARGET} PUBLIC cxx_alias_templates cxx_defaulted_functions cxx_delegating_constructors
        cxx_explicit_conversions cxx_generalized_initializers cxx_inheriting_constructors cxx_rvalue_references
        cxx_lambdas cxx_nullptr cxx_auto_type cxx_range_for cxx_right_angle_brackets cxx_deleted_functions cxx_nullptr 
        cxx_strong_enums cxx_uniform_initialization cxx_template_template_parameters cxx_defaulted_move_initializers)
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

