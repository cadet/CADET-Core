# =============================================================================
#  CADET - The Chromatography Analysis and Design Toolkit
#  
#  Copyright © 2008-2015: Eric von Lieres¹, Joel Andersson¹,
#                         Andreas Puettmann¹, Sebastian Schnittert¹,
#                         Samuel Leweke¹
#                                      
#    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
#  
#  All rights reserved. This program and the accompanying materials
#  are made available under the terms of the GNU Public License v3.0 (or, at
#  your option, any later version) which accompanies this distribution, and
#  is available at http://www.gnu.org/licenses/gpl.html
# =============================================================================

# This module looks for MATLAB
# Defines:
#  MATLAB_ROOT:        root directory of most recent Matlab installation
#  MATLAB_INCLUDE_DIR: include path for mex.h, engine.h
#  MATLAB_LIB_DIR:     directory containing the Matlab libraries
#  MATLAB_MEX_EXT:     file extension of mex files on this platform
#  MATLAB_LIBRARIES:   required libraries: libmex, etc
#  MATLAB_MEX_LIBRARY: path to libmex.lib
#  MATLAB_MX_LIBRARY:  path to libmx.lib
#  MATLAB_ENG_LIBRARY: path to libeng.lib


include(FindPackageHandleStandardArgs)

set (MATLAB_FOUND FALSE)

# Detect architecture
set (ARCHITECTURE i386)
set (BITNESS 32)
if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    set (BITNESS 64)
    set (ARCHITECTURE amd64)
endif()

# Step 1: Find Matlab root
if (WIN32)

    # Set mex file extension
    if (BITNESS EQUAL 32)
        set (MATLAB_MEX_EXT "mexw32")
    else ()
        set (MATLAB_MEX_EXT "mexw64")
    endif ()

    # Find the most recent Matlab version in the registry
    if (NOT DEFINED MATLAB_ROOT)

        set (ROOT_CANDIDATES    "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\8.3;MATLABROOT]"     # Future
                                "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\8.2;MATLABROOT]"     # R2013b
                                "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\8.1;MATLABROOT]"     # R2013a
                                "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\8.0;MATLABROOT]"     # R2012b
                                "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.14;MATLABROOT]"    # R2012a
                                "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.13;MATLABROOT]"    # R2011b
                                "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.12;MATLABROOT]"    # R2011a
                                "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.11.1;MATLABROOT]"  # R2010bSP1
                                "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.11;MATLABROOT]"    # R2010b
                                "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.10;MATLABROOT]"    # R2010a
                                "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.9.1;MATLABROOT]"   # R2009bSP1
                                "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.9;MATLABROOT]"     # R2009b
                                "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MathWorks\\MATLAB\\7.8;MATLABROOT]"     # R2009a
            )

        # Look in candidates
        find_path (MATLAB_INCLUDE_DIR "mex.h" PATHS ENV MATLAB_ROOT ${ROOT_CANDIDATES} PATH_SUFFIXES extern/include)

        if (NOT IS_DIRECTORY ${MATLAB_INCLUDE_DIR})
            return ()
        else ()
            get_filename_component(MATLAB_ROOT ${MATLAB_INCLUDE_DIR}/../.. ABSOLUTE)
        endif ()

    endif ()
    set (MATLAB_LIB_DIR_CANDIDATES ${MATLAB_ROOT}/extern/lib/win${BITNESS}/microsoft ${MATLAB_ROOT}/extern/lib ${MATLAB_ROOT})

else (WIN32)

    if (APPLE)
        set (MATLAB_MEX_EXT "mexmaci64")
        set (MATLAB_ARCH "maci64")
    else (APPLE)
        set (MATLAB_MEX_EXT "mexa64")
        set (MATLAB_ARCH "glnxa64")
    endif (APPLE)

    # Find the most recent Matlab version in the filesystem
    if (NOT DEFINED MATLAB_ROOT)

        if (APPLE)
            set (ROOT_CANDIDATES    "/Applications/MATLAB\ R2014a"
                                    "/Applications/MATLAB\ R2013b"
                                    "/Applications/MATLAB\ R2013a"
                                    "/Applications/MATLAB\ R2012b"
                                    "/Applications/MATLAB\ R2012a"
                                    "/Applications/MATLAB\ R2011b"
                                    "/Applications/MATLAB\ R2011a"
                                    "/Applications/MATLAB\ R2010bSP1"
                                    "/Applications/MATLAB\ R2010b"
                                    "/Applications/MATLAB\ R2010a"
                                    "/Applications/MATLAB\ R2009bSP1"
                                    "/Applications/MATLAB\ R2009b"
                                    "/Applications/MATLAB\ R2009a"
                )
        else (APPLE)
            set (ROOT_CANDIDATES    "/usr/local/MATLAB/R2014a"
                                    "/usr/local/MATLAB/R2013b"
                                    "/usr/local/MATLAB/R2013a"
                                    "/usr/local/MATLAB/R2012b"
                                    "/usr/local/MATLAB/R2012a"
                                    "/usr/local/MATLAB/R2011b"
                                    "/usr/local/MATLAB/R2011a"
                                    "/usr/local/MATLAB/R2010bSP1"
                                    "/usr/local/MATLAB/R2010b"
                                    "/usr/local/MATLAB/R2010a"
                                    "/usr/local/MATLAB/R2009bSP1"
                                    "/usr/local/MATLAB/R2009b"
                                    "/usr/local/MATLAB/R2009a"
                                    "/usr/local/matlab/R2014a"
                                    "/usr/local/matlab/R2013b"
                                    "/usr/local/matlab/R2013a"
                                    "/usr/local/matlab/R2012b"
                                    "/usr/local/matlab/R2012a"
                                    "/usr/local/matlab/R2011b"
                                    "/usr/local/matlab/R2011a"
                                    "/usr/local/matlab/R2010bSP1"
                                    "/usr/local/matlab/R2010b"
                                    "/usr/local/matlab/R2010a"
                                    "/usr/local/matlab/R2009bSP1"
                                    "/usr/local/matlab/R2009b"
                                    "/usr/local/matlab/R2009a"
                                    "/opt/MATLAB/R2014a"
                                    "/opt/MATLAB/R2013b"
                                    "/opt/MATLAB/R2013a"
                                    "/opt/MATLAB/R2012b"
                                    "/opt/MATLAB/R2012a"
                                    "/opt/MATLAB/R2011b"
                                    "/opt/MATLAB/R2011a"
                                    "/opt/MATLAB/R2010bSP1"
                                    "/opt/MATLAB/R2010b"
                                    "/opt/MATLAB/R2010a"
                                    "/opt/MATLAB/R2009bSP1"
                                    "/opt/MATLAB/R2009b"
                                    "/opt/MATLAB/R2009a"
                                    "/opt/matlab/R2014a"
                                    "/opt/matlab/R2013b"
                                    "/opt/matlab/R2013a"
                                    "/opt/matlab/R2012b"
                                    "/opt/matlab/R2012a"
                                    "/opt/matlab/R2011b"
                                    "/opt/matlab/R2011a"
                                    "/opt/matlab/R2010bSP1"
                                    "/opt/matlab/R2010b"
                                    "/opt/matlab/R2010a"
                                    "/opt/matlab/R2009bSP1"
                                    "/opt/matlab/R2009b"
                                    "/opt/matlab/R2009a"
                )
        endif (APPLE)

        # Look in candidates
        find_path (MATLAB_INCLUDE_DIR "mex.h" PATHS ENV MATLAB_ROOT ${ROOT_CANDIDATES} PATH_SUFFIXES extern/include)

        if (NOT IS_DIRECTORY ${MATLAB_INCLUDE_DIR})
            MESSAGE(STATUS "Matlab root: Not found")
            return ()
        else ()
            get_filename_component(MATLAB_ROOT ${MATLAB_INCLUDE_DIR}/../.. ABSOLUTE)
            MESSAGE(STATUS "Matlab root: ${MATLAB_ROOT}")
        endif ()

    endif ()

    set (MATLAB_LIB_DIR_CANDIDATES ${MATLAB_ROOT}/bin/${MATLAB_ARCH} ${MATLAB_ROOT})

endif (WIN32)

MESSAGE(STATUS "Looking for Matlab libs in ${MATLAB_LIB_DIR_CANDIDATES}")

# Step 2: Find library directory and libs
#set (MATLAB_INCLUDE_DIR "${MATLAB_ROOT}/extern/include")
find_library (MATLAB_MEX_LIB    NAMES libmex      mex       PATHS ${MATLAB_LIB_DIR_CANDIDATES} NO_DEFAULT_PATH)
find_library (MATLAB_MX_LIB     NAMES libmx       mx        PATHS ${MATLAB_LIB_DIR_CANDIDATES} NO_DEFAULT_PATH)
find_library (MATLAB_MAT_LIB    NAMES libmat      mat       PATHS ${MATLAB_LIB_DIR_CANDIDATES} NO_DEFAULT_PATH)
find_library (MATLAB_ENG_LIB    NAMES libeng      eng       PATHS ${MATLAB_LIB_DIR_CANDIDATES} NO_DEFAULT_PATH)
find_library (MATLAB_LAPACK_LIB NAMES libmwlapack mwlapack  PATHS ${MATLAB_LIB_DIR_CANDIDATES} NO_DEFAULT_PATH)
find_library (MATLAB_BLAS_LIB   NAMES libmwblas   mwblas    PATHS ${MATLAB_LIB_DIR_CANDIDATES} NO_DEFAULT_PATH)

set (MATLAB_LIBRARIES ${MATLAB_MEX_LIB} ${MATLAB_MX_LIB} ${MATLAB_MAT_LIB} ${MATLAB_ENG_LIB} ${MATLAB_LAPACK_LIB} ${MATLAB_BLAS_LIB})
get_filename_component(MATLAB_LIB_DIR ${MATLAB_MEX_LIB} PATH)

if (MATLAB_INCLUDE_DIR AND MATLAB_MEX_LIB AND MATLAB_LIBRARIES)
    set (MATLAB_FOUND TRUE)
endif (MATLAB_INCLUDE_DIR AND MATLAB_MEX_LIB AND MATLAB_LIBRARIES)

find_package_handle_standard_args( MATLAB DEFAULT_MSG
    MATLAB_LIBRARIES
    MATLAB_INCLUDE_DIR
)

mark_as_advanced (
    MATLAB_FOUND
    MATLAB_ROOT
    MATLAB_LIB_DIR
    MATLAB_INCLUDE_DIR
    MATLAB_MEX_EXT
    MATLAB_LIBRARIES
    MATLAB_MEX_LIB
    MATLAB_MX_LIB
    MATLAB_MAT_LIB
    MATLAB_ENG_LIB
    MATLAB_LAPACK_LIB
    MATLAB_BLAS_LIB
)
