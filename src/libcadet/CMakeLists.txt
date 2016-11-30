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
project (CadetLibrary CXX C)

# Option that allows users to build static and/or shared libraries
option (BUILD_SHARED_LIBS "Build shared libraries" ON)
option (BUILD_STATIC_LIBS "Build static libraries" OFF)

if (STANDALONE)
    set(BUILD_SHARED_LIBS OFF)
    set(BUILD_STATIC_LIBS ON)
endif()

# Make sure we build at least one type of libraries
if (NOT BUILD_STATIC_LIBS AND NOT BUILD_SHARED_LIBS)
    print_warning ("Both static and shared library generation were disabled"
        "Building shared libraries was re-enabled")
    set (BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries")
endif ()


# Git integration
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/VersionInfo.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/VersionInfo.cpp" @ONLY)

# LAPACK interface
#configure_file("${CMAKE_CURRENT_SOURCE_DIR}/LapackInterface.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/LapackInterface.hpp" @ONLY)

# Compiler features
write_compiler_detection_header(
    FILE "${CMAKE_BINARY_DIR}/cadet/cadetCompilerInfo.hpp"
    PREFIX CADET
    COMPILERS GNU Clang AppleClang MSVC #Intel (unsupported at the moment)
    FEATURES cxx_noexcept cxx_user_literals cxx_constexpr cxx_variadic_templates
    PROLOG "\
// =============================================================================\n\
//  CADET - The Chromatography Analysis and Design Toolkit\n\
//  \n\
//  Copyright © 2008-2016: The CADET Authors\n\
//            Please see the AUTHORS and CONTRIBUTORS file.\n\
//  \n\
//  All rights reserved. This program and the accompanying materials\n\
//  are made available under the terms of the GNU Public License v3.0 (or, at\n\
//  your option, any later version) which accompanies this distribution, and\n\
//  is available at http://www.gnu.org/licenses/gpl.html\n\
// =============================================================================\n"
    EPILOG "\
#if CADET_COMPILER_CXX_USER_LITERALS && CADET_COMPILER_CXX_CONSTEXPR\n\
    #define CADET_COMPILETIME_HASH 1\n\
#else\n\
    #define CADET_COMPILETIME_HASH 0\n\
#endif\n"
)

# Check for some compiler bugs
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)
#  include(${CMAKE_ROOT}/Modules/CheckCXXSourceCompiles.cmake)
#  include(${CMAKE_ROOT}/Modules/CMakePushCheckState.cmake)

cmake_push_check_state(RESET)

# Enable C++11
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_REQUIRED_FLAGS "-std=c++0x")
elseif (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang"))
    set(CMAKE_REQUIRED_FLAGS "-std=c++11")
elseif (CMAKE_COMPILER_IS_GNUCXX OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
    set(CMAKE_REQUIRED_FLAGS "-std=c++11")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(CMAKE_REQUIRED_FLAGS "")
endif ()

set(CMAKE_REQUIRED_INCLUDES "${CMAKE_BINARY_DIR}/cadet/")
CHECK_CXX_SOURCE_COMPILES("#include \"cadetCompilerInfo.hpp\" \n class Test { public: Test() { }\n Test& operator=(Test&& cpy) CADET_NOEXCEPT = default; \n private: \n int a; }; int main(int argc, char** argv) { return 0; }" COMPILER_SUPPORT_NOEXCEPT_DEFAULTED_MOVE)

cmake_pop_check_state()

# Add variable LIBCADET_SOURCES with the sources for LIBCADET
set (LIBCADET_SOURCES
     ${CMAKE_CURRENT_BINARY_DIR}/VersionInfo.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/Logging.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/FactoryFuncs.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/ModelBuilderImpl.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/SimulatorImpl.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/AutoDiff.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/AdUtils.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/Weno.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/BindingModelFactory.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/ModelSystemImpl.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/InletModel.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/OutletModel.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/GeneralRateModel.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/GeneralRateModel-LinearSolver.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/GeneralRateModel-InitialConditions.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/binding/BindingModelBase.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/binding/LinearBinding.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/binding/StericMassActionBinding.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/binding/LangmuirBinding.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/binding/AntiLangmuirBinding.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/binding/BiLangmuirBinding.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/binding/KumarLangmuirBinding.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/binding/MobilePhaseModulatorLangmuirBinding.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/binding/SelfAssociationBinding.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/binding/BiStericMassActionBinding.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/binding/SaskaBinding.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/inlet/PiecewiseCubicPoly.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/extfun/LinearInterpolationExternalFunction.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/model/extfun/PiecewiseCubicPolyExternalFunction.cpp
    )

# Add variable LIBCADET_NONLINALG_SOURCES with the sources for LIBCADET_NONLINALG
set (LIBCADET_NONLINALG_SOURCES
     ${CMAKE_SOURCE_DIR}/src/libcadet/linalg/BandMatrix.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/linalg/DenseMatrix.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/linalg/SparseMatrix.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/linalg/Gmres.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/nonlin/AdaptiveTrustRegionNewton.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/nonlin/LevenbergMarquardt.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/nonlin/CompositeSolver.cpp
     ${CMAKE_SOURCE_DIR}/src/libcadet/nonlin/Solver.cpp
    )


set (LIBCADET_TARGETS )
set (LIBNONLINALG_TARGETS )

if (LAPACK_FOUND)
    set(LIB_LAPACK_DEFINE "-DCADET_LAPACK_TRAILING_UNDERSCORE")

    # Add the build target for the static nonlinalg library
    add_library(libcadet_nonlinalg_static STATIC ${LIBCADET_NONLINALG_SOURCES})
    set_target_properties(libcadet_nonlinalg_static PROPERTIES OUTPUT_NAME cadet_nonlinalg_static)
    target_compile_definitions(libcadet_nonlinalg_static PRIVATE -Dlibcadet_nonlinalg_static_EXPORTS ${LIB_LAPACK_DEFINE})

    list (APPEND LIBNONLINALG_TARGETS libcadet_nonlinalg_static)

    # ---------------------------------------------------
    #   Build the static library
    # ---------------------------------------------------
    if (BUILD_STATIC_LIBS AND LAPACK_FOUND)

        # Add the build target for the static CADET library
        add_library(libcadet_static STATIC ${LIBCADET_SOURCES})
        set_target_properties(libcadet_static PROPERTIES OUTPUT_NAME cadet_static)
        target_compile_definitions(libcadet_static PRIVATE -Dlibcadet_static_EXPORTS ${LIB_LAPACK_DEFINE})
        target_link_libraries(libcadet_static PRIVATE libcadet_nonlinalg_static)

        list (APPEND LIBCADET_TARGETS libcadet_static)

    endif ()
    # ---------------------------------------------------

    # ---------------------------------------------------
    #   Build the shared library
    # ---------------------------------------------------
    if (BUILD_SHARED_LIBS AND LAPACK_FOUND)

        # Add the build target for the shared CADET library
        add_library(libcadet_shared SHARED ${LIBCADET_SOURCES})
        set_target_properties(libcadet_shared PROPERTIES OUTPUT_NAME cadet)
        target_compile_definitions(libcadet_shared PRIVATE ${LIB_LAPACK_DEFINE})
        target_link_libraries(libcadet_shared PRIVATE libcadet_nonlinalg_static)

        list (APPEND LIBCADET_TARGETS libcadet_shared)

    endif ()
    # ---------------------------------------------------
endif ()

set(ALL_LOCAL_TARGETS ${LIBNONLINALG_TARGETS} ${LIBCADET_TARGETS})

foreach(_TARGET IN LISTS ALL_LOCAL_TARGETS)
    set_target_properties(${_TARGET} PROPERTIES DEBUG_POSTFIX _d VERSION ${CADET_VERSION} SOVERSION "0")
    set_target_properties(${_TARGET} PROPERTIES POSITION_INDEPENDENT_CODE ON)

    target_compile_features(${_TARGET} PUBLIC cxx_alias_templates cxx_defaulted_functions cxx_delegating_constructors
        cxx_explicit_conversions cxx_generalized_initializers cxx_inheriting_constructors cxx_rvalue_references
        cxx_lambdas cxx_nullptr cxx_auto_type cxx_range_for cxx_right_angle_brackets cxx_deleted_functions cxx_nullptr 
        cxx_strong_enums cxx_uniform_initialization cxx_template_template_parameters cxx_defaulted_move_initializers)

    target_link_libraries(${_TARGET} PRIVATE ${SUNDIALS_LIBRARIES})
endforeach()

# Link to BLAS and LAPACK (except for Intel compiler which uses Intel MKL)
if (NOT ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel"))
    foreach(_TARGET IN LISTS LIBNONLINALG_TARGETS)
         target_link_libraries (${_TARGET} PRIVATE ${LAPACK_LIBRARIES})
         set_target_properties (${_TARGET} PROPERTIES LINK_FLAGS "${LAPACK_LINKER_FLAGS}")
    endforeach()
endif()

# Add include directories
foreach(_TARGET IN LISTS LIBNONLINALG_TARGETS)
    target_include_directories (${_TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src/libcadet ${CMAKE_BINARY_DIR} ${SUNDIALS_INCLUDE_DIRS})
endforeach()
foreach(_TARGET IN LISTS LIBCADET_TARGETS)
    target_include_directories (${_TARGET} PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include> 
        PRIVATE ${CMAKE_SOURCE_DIR}/src/libcadet ${CMAKE_BINARY_DIR} ${SUNDIALS_INCLUDE_DIRS})
endforeach()

# ---------------------------------------------------
#   Build the Matlab library
# ---------------------------------------------------
if (BUILD_CADET_MEX)

    set(MATLAB_LAPACK_DEFINE)
    if (Matlab_LAPACK_TRAILING_UNDERSCORE)
        set(MATLAB_LAPACK_DEFINE "-DCADET_LAPACK_TRAILING_UNDERSCORE")
    endif()

    # Add the build target for the mex nonlinalg library
    add_library(libcadet_nonlinalg_mex STATIC ${LIBCADET_NONLINALG_SOURCES})
    set_target_properties(libcadet_nonlinalg_mex PROPERTIES OUTPUT_NAME cadet_nonlinalg_mex)
    set_target_properties(libcadet_nonlinalg_mex PROPERTIES LINK_SEARCH_START_STATIC ON)
    set_target_properties(libcadet_nonlinalg_mex PROPERTIES LINK_SEARCH_END_STATIC ON)
    target_compile_definitions(libcadet_nonlinalg_mex PRIVATE -DCADET_MATLABMEX -DMATLAB_MEX_FILE ${MATLAB_LAPACK_DEFINE} -Dlibcadet_mex_EXPORTS)
    target_include_directories(libcadet_nonlinalg_mex PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src/libcadet ${CMAKE_BINARY_DIR} ${SUNDIALS_INCLUDE_DIRS})

    set_target_properties(libcadet_nonlinalg_mex PROPERTIES DEBUG_POSTFIX _d VERSION ${CADET_VERSION} SOVERSION "0")
    set_target_properties(libcadet_nonlinalg_mex PROPERTIES POSITION_INDEPENDENT_CODE ON)

    target_compile_features(libcadet_nonlinalg_mex PUBLIC cxx_alias_templates cxx_defaulted_functions cxx_delegating_constructors
        cxx_explicit_conversions cxx_generalized_initializers cxx_inheriting_constructors cxx_rvalue_references
        cxx_lambdas cxx_nullptr cxx_auto_type cxx_range_for cxx_right_angle_brackets cxx_deleted_functions cxx_nullptr 
        cxx_strong_enums cxx_uniform_initialization cxx_template_template_parameters cxx_defaulted_move_initializers)

    target_link_libraries (libcadet_nonlinalg_mex PRIVATE ${SUNDIALS_STATIC_LIBRARIES})

    list(APPEND ALL_LOCAL_TARGETS libcadet_nonlinalg_mex)

    # Add the build target for the mex CADET library
    add_library(libcadet_mex STATIC ${LIBCADET_SOURCES})
    set_target_properties(libcadet_mex PROPERTIES OUTPUT_NAME cadet_mex)
    set_target_properties(libcadet_mex PROPERTIES LINK_SEARCH_START_STATIC ON)
    set_target_properties(libcadet_mex PROPERTIES LINK_SEARCH_END_STATIC ON)
    target_compile_definitions(libcadet_mex PRIVATE -DCADET_MATLABMEX -DMATLAB_MEX_FILE ${MATLAB_LAPACK_DEFINE} -Dlibcadet_mex_EXPORTS)
    target_link_libraries(libcadet_mex PRIVATE libcadet_nonlinalg_mex ${SUNDIALS_STATIC_LIBRARIES})
    target_include_directories (libcadet_mex PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include> 
        PRIVATE ${CMAKE_SOURCE_DIR}/src/libcadet ${CMAKE_BINARY_DIR} ${SUNDIALS_INCLUDE_DIRS})

    set_target_properties(libcadet_mex PROPERTIES DEBUG_POSTFIX _d VERSION ${CADET_VERSION} SOVERSION "0")
    set_target_properties(libcadet_mex PROPERTIES POSITION_INDEPENDENT_CODE ON)

    target_compile_features(libcadet_mex PUBLIC cxx_alias_templates cxx_defaulted_functions cxx_delegating_constructors
        cxx_explicit_conversions cxx_generalized_initializers cxx_inheriting_constructors cxx_rvalue_references
        cxx_lambdas cxx_nullptr cxx_auto_type cxx_range_for cxx_right_angle_brackets cxx_deleted_functions cxx_nullptr 
        cxx_strong_enums cxx_uniform_initialization cxx_template_template_parameters cxx_defaulted_move_initializers)

    list(APPEND LIBCADET_TARGETS libcadet_mex)
    list(APPEND ALL_LOCAL_TARGETS libcadet_mex)

    unset(MATLAB_LAPACK_DEFINE)
endif ()
# ---------------------------------------------------

unset(LIBCADET_NONLINALG_SOURCES)
unset(LIBCADET_SOURCES)

# ---------------------------------------------------
#   Other configuration options
# ---------------------------------------------------

# Option that allows users to switch on verification of analytical Jacobian by AD
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    option (CHECK_ANALYTIC_JACOBIAN "Toggle verification of analytical Jacobian by AD (only Debug builds)" OFF)
endif ()

# Option that allows users to specify the AD library
set (ADLIB "sfad" CACHE STRING "Selects the AD library, options are 'adolc', 'sfad', 'setfad'")
string(TOLOWER ${ADLIB} ADLIB)
if (ADLIB STREQUAL "adolc")
    message (STATUS "AD library: ADOL-C")
elseif (ADLIB STREQUAL "sfad")
    message (STATUS "AD library: SFAD")
elseif (ADLIB STREQUAL "setfad")
    message (STATUS "AD library: SETFAD")
else ()
    message (FATAL_ERROR "Unkown AD library ${ADLIB} (options are 'adolc', 'sfad', 'setfad')")
endif ()

foreach(_TARGET IN LISTS LIBCADET_TARGETS)
    if (CHECK_ANALYTIC_JACOBIAN AND CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_definitions(${_TARGET} PRIVATE -DCADET_CHECK_ANALYTIC_JACOBIAN)
    endif ()

    if (ADLIB STREQUAL "adolc")
        target_compile_definitions(${_TARGET} PRIVATE -DACTIVE_ADOLC)
        target_include_directories(${_TARGET} PRIVATE "${CMAKE_SOURCE_DIR}/ThirdParty/ADOL-C/include")
    elseif (ADLIB STREQUAL "sfad")
        target_compile_definitions(${_TARGET} PRIVATE -DACTIVE_SFAD)
        target_include_directories(${_TARGET} PRIVATE "${CMAKE_SOURCE_DIR}/include/ad")
    elseif (ADLIB STREQUAL "setfad")
        target_compile_definitions(${_TARGET} PRIVATE -DACTIVE_SETFAD)
        target_include_directories(${_TARGET} PRIVATE "${CMAKE_SOURCE_DIR}/include/ad")
    endif ()
endforeach()
# ---------------------------------------------------

# ---------------------------------------------------
#   Set installation related stuff
# ---------------------------------------------------

install (CODE "MESSAGE(\"\nInstall LIBCADET\n\")")
install(TARGETS ${ALL_LOCAL_TARGETS} EXPORT libcadet-targets
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
install(FILES ${CMAKE_BINARY_DIR}/cadet/cadetCompilerInfo.hpp DESTINATION include/cadet)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/cadet DESTINATION include)
install(EXPORT libcadet-targets DESTINATION lib)

# ---------------------------------------------------

# Make variable visible in parent scope
set (LIB_TARGETS ${ALL_LOCAL_TARGETS} PARENT_SCOPE)
unset(LIBCADET_TARGETS)
unset(LIBNONLINALG_TARGETS)
unset(ALL_LOCAL_TARGETS)

# Info message
message (STATUS "Added LIBCADET module")
