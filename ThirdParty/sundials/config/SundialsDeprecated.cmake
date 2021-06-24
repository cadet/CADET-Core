# ---------------------------------------------------------------
# Programmer:  David J. Gardner @ LLNL
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
# Print warning is the user sets a deprecated CMake variable and
# copy the value into the correct CMake variable
# ---------------------------------------------------------------

# SUNDIALS_INDEX_TYPE got new behavior
if(SUNDIALS_INDEX_TYPE)
  string(TOUPPER ${SUNDIALS_INDEX_TYPE} tmp)

  if(tmp STREQUAL "INT32_T")
    PRINT_WARNING("SUNDIALS_INDEX_TYPE overrides the standard types SUNDIALS looks for."
    "Setting SUNDIALS_INDEX_SIZE to 32 and clearing SUNDIALS_INDEX_TYPE.")
    FORCE_VARIABLE(SUNDIALS_INDEX_SIZE STRING "SUNDIALS index size" 32)
    FORCE_VARIABLE(SUNDIALS_INDEX_TYPE STRING "SUNDIALS index type" "")
  elseif(tmp STREQUAL "INT64_T")
    PRINT_WARNING("SUNDIALS_INDEX_TYPE overrides the standard types SUNDIALS looks for."
    "Setting SUNDIALS_INDEX_SIZE to 64 and clearing SUNDIALS_INDEX_TYPE.")
    FORCE_VARIABLE(SUNDIALS_INDEX_SIZE STRING "SUNDIALS index size" 64)
    FORCE_VARIABLE(SUNDIALS_INDEX_TYPE STRING "SUNDIALS index type" "")
  else()
    PRINT_WARNING("SUNDIALS_INDEX_TYPE overrides the standard types SUNDIALS looks for." "")
  endif()
endif()
