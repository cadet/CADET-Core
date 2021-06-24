# ---------------------------------------------------------------
# Programmer:  Radu Serban @ LLNL
# ---------------------------------------------------------------
# LLNS Copyright Start
# Copyright (c) 2014, Lawrence Livermore National Security
# This work was performed under the auspices of the U.S. Department 
# of Energy by Lawrence Livermore National Laboratory in part under 
# Contract W-7405-Eng-48 and in part under Contract DE-AC52-07NA27344.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# LLNS Copyright End
# ---------------------------------------------------------------
# CMakeLists.txt file for the serial NVECTOR library

# Add variable nvecserial_SOURCES with the sources for the NVECSERIAL lib
SET(nvecserial_SOURCES nvector_serial.c)

# Add variable shared_SOURCES with the common SUNDIALS sources which will
# also be included in the NVECSERIAL library
SET(shared_SOURCES
  ${sundials_SOURCE_DIR}/src/sundials/sundials_math.c
  )

# Add variable nvecserial_HEADERS with the exported NVECSERIAL header files
SET(nvecserial_HEADERS
  ${sundials_SOURCE_DIR}/include/nvector/nvector_serial.h
  )

# Add source directory to include directories
INCLUDE_DIRECTORIES(.)

# Define C preprocessor flag -DBUILD_SUNDIALS_LIBRARY 
ADD_DEFINITIONS(-DBUILD_SUNDIALS_LIBRARY)

# Rules for building and installing the static library:
#  - Add the build target for the NVECSERIAL library
#  - Set the library name and make sure it is not deleted
#  - Install the NVECSERIAL library
ADD_LIBRARY(sundials_nvecserial_static STATIC ${nvecserial_SOURCES} ${shared_SOURCES})
SET_TARGET_PROPERTIES(sundials_nvecserial_static
  PROPERTIES OUTPUT_NAME sundials_nvecserial CLEAN_DIRECT_OUTPUT 1)
