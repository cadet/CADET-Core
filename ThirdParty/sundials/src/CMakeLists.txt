# ---------------------------------------------------------------
# Programmer:  David Gardner and Slaven Peles @ LLNL
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
# src level CMakeLists.txt for SUNDIALS (for cmake build system)
# ---------------------------------------------------------------

# Always add SUNDIALS provided serial modules
ADD_SUBDIRECTORY(nvec_ser)

# Always add SUNDIALS provided iterative linear solver modules
ADD_SUBDIRECTORY(sunlinsol_spgmr)
ADD_SUBDIRECTORY(sunlinsol_spfgmr)
ADD_SUBDIRECTORY(sunlinsol_spbcgs)
ADD_SUBDIRECTORY(sunlinsol_sptfqmr)

# IDAS library
ADD_SUBDIRECTORY(idas)
