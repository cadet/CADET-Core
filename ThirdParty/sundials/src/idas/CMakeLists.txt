# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds @ SMU
#              Radu Serban @ LLNL
# ---------------------------------------------------------------
# LLNS/SMU Copyright Start
# Copyright (c) 2017, Southern Methodist University and 
# Lawrence Livermore National Security
#
# This work was performed under the auspices of the U.S. Department 
# of Energy by Southern Methodist University and Lawrence Livermore 
# National Laboratory under Contract DE-AC52-07NA27344.
# Produced at Southern Methodist University and the Lawrence 
# Livermore National Laboratory.
#
# All rights reserved.
# For details, see the LICENSE file.
# LLNS/SMU Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the IDAS library

# Add variable idas_SOURCES with the sources for the IDAS library
SET(idas_SOURCES
  idas.c
  idaa.c
  idas_io.c
  idas_ic.c
  idaa_io.c
  idas_direct.c
  idas_spils.c
  idas_bbdpre.c
  )

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the IDAS library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_matrix.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_band.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_dense.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_direct.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_iterative.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_version.c
  ${sundials_SOURCE_DIR}/src/nvec_ser/nvector_serial.c
  )

# Add variable sunlinsol_SOURCES with the common SUNLinearSolver sources which will
# also be included in the IDAS library
SET(sunlinsol_SOURCES
  ${sundials_SOURCE_DIR}/src/sunlinsol_spbcgs/sunlinsol_spbcgs.c
  ${sundials_SOURCE_DIR}/src/sunlinsol_spfgmr/sunlinsol_spfgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol_spgmr/sunlinsol_spgmr.c
  ${sundials_SOURCE_DIR}/src/sunlinsol_sptfqmr/sunlinsol_sptfqmr.c
  )

# Add variable idas_HEADERS with the exported IDAS header files
SET(idas_HEADERS
  idas.h
  idas_bbdpre.h
  idas_direct.h
  idas_spils.h
  )

# Add prefix with complete path to the IDAS header files
ADD_PREFIX(${sundials_SOURCE_DIR}/include/idas/ idas_HEADERS)

# Add source directories to include directories for access to
# implementation only header files.
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(../sundials)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Add the build target for the static IDAS library
ADD_LIBRARY(sundials_idas_static STATIC 
  ${idas_SOURCES} ${shared_SOURCES} ${sunlinsol_SOURCES})

# Set the library name and make sure it is not deleted
SET_TARGET_PROPERTIES(sundials_idas_static
  PROPERTIES OUTPUT_NAME sundials_idas CLEAN_DIRECT_OUTPUT 1)
