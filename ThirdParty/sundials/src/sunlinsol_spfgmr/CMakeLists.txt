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
# CMakeLists.txt file for the SPFGMR SUNLinearSolver library

# Add variable sunlinsolspfgmr_SOURCES with the sources for the SUNLINSOLSPFGMR lib
SET(sunlinsolspfgmr_SOURCES sunlinsol_spfgmr.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the SUNLINSOLSPFGMR library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_nvector.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_linearsolver.c
  ${sundials_SOURCE_DIR}/src/sundials/sundials_iterative.c
  )

# Add variable sunlinsolspfgmr_HEADERS with the exported SUNLINSOLSPFGMR header files
SET(sunlinsolspfgmr_HEADERS 
  ${sundials_SOURCE_DIR}/include/sunlinsol/sunlinsol_spfgmr.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the SUNLINSOLSPFGMR library
#  - Set the library name and make sure it is not deleted
#  - Install the SUNLINSOLSPFGMR library
ADD_LIBRARY(sundials_sunlinsolspfgmr_static STATIC ${sunlinsolspfgmr_SOURCES} ${shared_SOURCES})
SET_TARGET_PROPERTIES(sundials_sunlinsolspfgmr_static
  PROPERTIES OUTPUT_NAME sundials_sunlinsolspfgmr CLEAN_DIRECT_OUTPUT 1)
