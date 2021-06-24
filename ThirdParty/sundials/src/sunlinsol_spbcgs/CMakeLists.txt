# ---------------------------------------------------------------
# Programmer:  Daniel R. Reynolds @ SMU
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
# CMakeLists.txt file for the SPBCGS SUNLinearSolver library

# Add variable sunlinsolspbcgs_SOURCES with the sources for the SUNLINSOLSPBCGS lib
SET(sunlinsolspbcgs_SOURCES sunlinsol_spbcgs.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLSPBCGS library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_iterative.c
  )

# Add variable sunlinsolspbcgs_HEADERS with the exported SUNLINSOLSPBCGS header files
SET(sunlinsolspbcgs_HEADERS 
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_spbcgs.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLSPBCGS library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLSPBCGS library
ADD_LIBRARY(sundials_sunlinsolspbcgs_static STATIC ${sunlinsolspbcgs_SOURCES} ${shared_SOURCES})
SET_TARGET_PROPERTIES(sundials_sunlinsolspbcgs_static
  PROPERTIES OUTPUT_NAME sundials_sunlinsolspbcgs CLEAN_DIRECT_OUTPUT 1)
