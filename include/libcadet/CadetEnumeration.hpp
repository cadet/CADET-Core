// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2014: Eric von Lieres¹, Joel Andersson,
//                         Andreas Puettmann¹, Sebastian Schnittert¹,
//                         Samuel Leweke¹
//                                      
//    ¹ Forschungszentrum Juelich GmbH, IBG-1, Juelich, Germany.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#ifndef LIBCADET_CADETENUMERATION_HPP_
#define LIBCADET_CADETENUMERATION_HPP_

#include "CadetException.hpp"

namespace cadet
{


enum ChromatographyType {
    ChromatographyType_begin,   // Must be first

    GENERAL_RATE_MODEL,
    LUMPED_RATE_MODEL, // might be implemented later

    ChromatographyType_end      // Must be last
};

enum AdsorptionType {
    AdsorptionType_begin,       // Must be first

    MULTI_COMPONENT_LANGMUIR,
    MOBILE_PHASE_MODULATORS,
    STERIC_MASS_ACTION,
    SELF_ASSOCIATION,
    EXTERNAL_LANGMUIR,
    EXTERNAL_STERIC_MASS_ACTION,
    LINEAR,
    MULTI_COMPONENT_BILANGMUIR,

    AdsorptionType_end          // Must be last
};

enum ReconstructionType {
    ReconstructionType_begin,    // Must be first

    WENO_REC,

    ReconstructionType_end       // Must be last
};

enum ParticleDiscType {
    ParticleDiscType_begin,     // Must be first

    EQUIDISTANT_PAR,
    EQUIVOLUME_PAR,
    USER_DEFINED_PAR,

    ParticleDiscType_end        // Must be last
};

// This enumeration must contain all parameters that are potentially set sensitive
enum ParameterName {
    ParameterName_begin,        // Must be first

    COL_LENGTH,
    COL_POROSITY,
    COL_DISPERSION,
    VELOCITY,

    FILM_DIFFUSION,

    PAR_RADIUS,
    PAR_POROSITY,
    PAR_DIFFUSION,
    PAR_SURFDIFFUSION,

    MCL_KA,
    MCL_KD,
    MCL_QMAX,

    MPM_KA,
    MPM_KD,
    MPM_QMAX,
    MPM_BETA,
    MPM_GAMMA,

    SMA_KA,
    SMA_KD,
    SMA_NU,
    SMA_SIGMA,
    SMA_LAMBDA,

    SAI_KA1,
    SAI_KA2,
    SAI_KD,
    SAI_NU,
    SAI_SIGMA,
    SAI_LAMBDA,

    MCBL_KA1,
    MCBL_KD1,
    MCBL_KA2,
    MCBL_KD2,
    MCBL_QMAX1,
    MCBL_QMAX2,

    EXTL_KA,
    EXTL_KA_T,
    EXTL_KA_TT,
    EXTL_KA_TTT,
    EXTL_KD,
    EXTL_KD_T,
    EXTL_KD_TT,
    EXTL_KD_TTT,
    EXTL_QMAX,
    EXTL_QMAX_T,
    EXTL_QMAX_TT,
    EXTL_QMAX_TTT,

    EXTSMA_KA,
    EXTSMA_KA_T,
    EXTSMA_KA_TT,
    EXTSMA_KA_TTT,
    EXTSMA_KD,
    EXTSMA_KD_T,
    EXTSMA_KD_TT,
    EXTSMA_KD_TTT,
    EXTSMA_NU,
    EXTSMA_NU_T,
    EXTSMA_NU_TT,
    EXTSMA_NU_TTT,
    EXTSMA_SIGMA,
    EXTSMA_SIGMA_T,
    EXTSMA_SIGMA_TT,
    EXTSMA_SIGMA_TTT,
    EXTSMA_LAMBDA,
    EXTSMA_LAMBDA_T,
    EXTSMA_LAMBDA_TT,
    EXTSMA_LAMBDA_TTT,

    LIN_KA,
    LIN_KD,

    INLET_PARAMETER,

    ParameterName_end           // Must be last
};

// This enumeration contains other parameters that should be able to be represented as string
enum OtherName {
    OtherName_begin,        // Must be first

    CHROMATOGRAPHY_TYPE,

    // Model group
    IS_KINETIC,
    ADSORPTION_TYPE,
    NCOMP,
    INIT_C,
    INIT_Q,

    // Column inlet group
    NSEC,
    SECTION_TIMES,
    SECTION_CONTINUITY,

    // External group
    EXT_PROFILE,
    EXT_PROF_DELTA,
    EXT_VELOCITY,

    // Discretization group
    NCOL,
    NPAR,
    PAR_DISC_TYPE,
    PAR_DISC_VECTOR,
    RECONSTRUCTION,

    // Weno scheme group
    BOUNDARY_MODEL,
    WENO_EPS,
    WENO_ORDER,

    // Schur solver group
    GS_TYPE,
    MAX_KRYLOV,
    MAX_RESTARTS,
    SCHUR_SAFETY,

    // Time integrator group
    ABSTOL,
    RELTOL,
    INIT_STEP_SIZE,
    MAX_STEPS,

    // Solver group
    PRINT_PROGRESS,
    PRINT_STATISTICS,
    PRINT_TIMING,
    PRINT_PARAMLIST,
    PRINT_CONFIG,
    USE_ANALYTIC_JACOBIAN,
    WRITE_AT_USER_TIMES,
    WRITE_SOLUTION_TIMES,
    WRITE_SOLUTION_COLUMN_OUTLET,
    WRITE_SOLUTION_COLUMN_INLET,
    WRITE_SENS_COLUMN_OUTLET,
    WRITE_SOLUTION_ALL,
    WRITE_SENS_ALL,
    LOG_LEVEL,
    USER_SOLUTION_TIMES,
    NTHREADS,

    // Sensitivity group
    NSENS,
    SENS_METHOD,
    SENS_NAME,
    SENS_COMP,
    SENS_SECTION,
    SENS_ABSTOL,
    SENS_FD_DELTA,

    // Output group
    SOLUTION_TIMES,
    SOLUTION_COLUMN,
    SOLUTION_PARTICLE,
    SOLUTION_BOUNDARY,
    SOLUTION_COLUMN_OUTLET,
    SOLUTION_COLUMN_INLET,
    SENS_COLUMN,
    SENS_PARTICLE,
    SENS_BOUNDARY,
    SENS_COLUMN_OUTLET,

    OtherName_end           // Must be last
};


enum GroupName {
    GroupName_begin,        // Must be first

    GRP_IN,

    GRP_IN_MODEL,
    GRP_IN_ADSORPTION,
    GRP_IN_INLET,
    GRP_IN_EXTERNAL,

    GRP_IN_DISCRETIZATION,
    GRP_IN_WENO,

    GRP_IN_SOLVER,
    GRP_IN_SCHUR,
    GRP_IN_TIME,

    GRP_IN_SENSITIVITY,

    GRP_OUT,

    GRP_OUT_SOLUTION,
    GRP_OUT_SENS,

    GroupName_end           // Must be last
};


inline const char* e2s(ChromatographyType chromType)
{
    switch (chromType)
    {
    case GENERAL_RATE_MODEL:          return "GENERAL_RATE_MODEL";
    case LUMPED_RATE_MODEL:           return "LUMPED_RATE_MODEL";
    default:                          return "UNKNOWN_CHROMATOGRAPHY_MODEL";
    }
}


inline const char* e2s(AdsorptionType adsType)
{
    switch (adsType)
    {
    case MULTI_COMPONENT_LANGMUIR:     return "MULTI_COMPONENT_LANGMUIR";
    case MOBILE_PHASE_MODULATORS:      return "MOBILE_PHASE_MODULATORS";
    case STERIC_MASS_ACTION:           return "STERIC_MASS_ACTION";
    case SELF_ASSOCIATION:             return "SELF_ASSOCIATION";
    case EXTERNAL_LANGMUIR:            return "EXTERNAL_LANGMUIR";
    case EXTERNAL_STERIC_MASS_ACTION:  return "EXTERNAL_STERIC_MASS_ACTION";
    case LINEAR:                       return "LINEAR";
    case MULTI_COMPONENT_BILANGMUIR:   return "MULTI_COMPONENT_BILANGMUIR";
    default:                           return "UNKNOWN_ADSORPTION_TYPE";
    }
}

inline const char* e2s(ReconstructionType recType)
{
    switch (recType)
    {
    case WENO_REC:                    return "WENO";
    default:                          return "UNKNOWN_RECONSTRUCTION";
    }
}

inline const char* e2s(ParticleDiscType parDiscType)
{
    switch (parDiscType)
    {
    case EQUIDISTANT_PAR:             return "EQUIDISTANT_PAR";
    case EQUIVOLUME_PAR:              return "EQUIVOLUME_PAR";
    case USER_DEFINED_PAR:            return "USER_DEFINED_PAR";
    default:                          return "UNKNOWN_PARTICLE_DISC_TYPE";
    }
}


inline const char* e2s(ParameterName param)
{
    switch (param)
    {
    case COL_LENGTH:                  return "COL_LENGTH";
    case COL_POROSITY:                return "COL_POROSITY";
    case COL_DISPERSION:              return "COL_DISPERSION";
    case VELOCITY:                    return "VELOCITY";

    case FILM_DIFFUSION:              return "FILM_DIFFUSION";

    case PAR_RADIUS:                  return "PAR_RADIUS";
    case PAR_POROSITY:                return "PAR_POROSITY";
    case PAR_DIFFUSION:               return "PAR_DIFFUSION";
    case PAR_SURFDIFFUSION:           return "PAR_SURFDIFFUSION";

    case MCL_KA:                      return "MCL_KA";
    case MCL_KD:                      return "MCL_KD";
    case MCL_QMAX:                    return "MCL_QMAX";

    case MPM_KA:                      return "MPM_KA";
    case MPM_KD:                      return "MPM_KD";
    case MPM_QMAX:                    return "MPM_QMAX";
    case MPM_BETA:                    return "MPM_BETA";
    case MPM_GAMMA:                   return "MPM_GAMMA";

    case SMA_KA:                      return "SMA_KA";
    case SMA_KD:                      return "SMA_KD";
    case SMA_NU:                      return "SMA_NU";
    case SMA_SIGMA:                   return "SMA_SIGMA";
    case SMA_LAMBDA:                  return "SMA_LAMBDA";

    case SAI_KA1:                     return "SAI_KA1";
    case SAI_KA2:                     return "SAI_KA2";
    case SAI_KD:                      return "SAI_KD";
    case SAI_NU:                      return "SAI_NU";
    case SAI_SIGMA:                   return "SAI_SIGMA";
    case SAI_LAMBDA:                  return "SAI_LAMBDA";

    case MCBL_KA1:                    return "MCBL_KA1";
    case MCBL_KD1:                    return "MCBL_KD1";
    case MCBL_KA2:                    return "MCBL_KA2";
    case MCBL_KD2:                    return "MCBL_KD2";
    case MCBL_QMAX1:                  return "MCBL_QMAX1";
    case MCBL_QMAX2:                  return "MCBL_QMAX2";

    case EXTL_KA:                     return "EXTL_KA";
    case EXTL_KA_T:                   return "EXTL_KA_T";
    case EXTL_KA_TT:                  return "EXTL_KA_TT";
    case EXTL_KA_TTT:                 return "EXTL_KA_TTT";
    case EXTL_KD:                     return "EXTL_KD";
    case EXTL_KD_T:                   return "EXTL_KD_T";
    case EXTL_KD_TT:                  return "EXTL_KD_TT";
    case EXTL_KD_TTT:                 return "EXTL_KD_TTT";
    case EXTL_QMAX:                   return "EXTL_QMAX";
    case EXTL_QMAX_T:                 return "EXTL_QMAX_T";
    case EXTL_QMAX_TT:                return "EXTL_QMAX_TT";
    case EXTL_QMAX_TTT:               return "EXTL_QMAX_TTT";

    case EXTSMA_KA:                   return "EXTSMA_KA";
    case EXTSMA_KA_T:                 return "EXTSMA_KA_T";
    case EXTSMA_KA_TT:                return "EXTSMA_KA_TT";
    case EXTSMA_KA_TTT:               return "EXTSMA_KA_TTT";
    case EXTSMA_KD:                   return "EXTSMA_KD";
    case EXTSMA_KD_T:                 return "EXTSMA_KD_T";
    case EXTSMA_KD_TT:                return "EXTSMA_KD_TT";
    case EXTSMA_KD_TTT:               return "EXTSMA_KD_TTT";
    case EXTSMA_NU:                   return "EXTSMA_NU";
    case EXTSMA_NU_T:                 return "EXTSMA_NU_T";
    case EXTSMA_NU_TT:                return "EXTSMA_NU_TT";
    case EXTSMA_NU_TTT:               return "EXTSMA_NU_TTT";
    case EXTSMA_SIGMA:                return "EXTSMA_SIGMA";
    case EXTSMA_SIGMA_T:              return "EXTSMA_SIGMA_T";
    case EXTSMA_SIGMA_TT:             return "EXTSMA_SIGMA_TT";
    case EXTSMA_SIGMA_TTT:            return "EXTSMA_SIGMA_TTT";
    case EXTSMA_LAMBDA:               return "EXTSMA_LAMBDA";
    case EXTSMA_LAMBDA_T:             return "EXTSMA_LAMBDA_T";
    case EXTSMA_LAMBDA_TT:            return "EXTSMA_LAMBDA_TT";
    case EXTSMA_LAMBDA_TTT:           return "EXTSMA_LAMBDA_TTT";

    case LIN_KA:                      return "LIN_KA";
    case LIN_KD:                      return "LIN_KD";

    case INLET_PARAMETER:             return "INLET_PARAMETER";

    default:                          return "UNKNOWN_PARAMETER";
    }
}


inline const char* e2s(OtherName other)
{
    switch (other)
    {
    case CHROMATOGRAPHY_TYPE:         return "CHROMATOGRAPHY_TYPE";

    case ADSORPTION_TYPE:             return "ADSORPTION_TYPE";
    case NCOMP:                       return "NCOMP";
    case INIT_C:                      return "INIT_C";
    case INIT_Q:                      return "INIT_Q";
    case IS_KINETIC:                  return "IS_KINETIC";
    case EXT_PROFILE:                 return "EXT_PROFILE";
    case EXT_PROF_DELTA:              return "EXT_PROF_DELTA";
    case EXT_VELOCITY:                return "EXT_VELOCITY";

    case NSEC:                        return "NSEC";
    case SECTION_TIMES:               return "SECTION_TIMES";
    case SECTION_CONTINUITY:          return "SECTION_CONTINUITY";

    case NCOL:                        return "NCOL";
    case NPAR:                        return "NPAR";
    case PAR_DISC_TYPE:               return "PAR_DISC_TYPE";
    case PAR_DISC_VECTOR:             return "PAR_DISC_VECTOR";
    case RECONSTRUCTION:              return "RECONSTRUCTION";

    case BOUNDARY_MODEL:              return "BOUNDARY_MODEL";
    case WENO_EPS:                    return "WENO_EPS";
    case WENO_ORDER:                  return "WENO_ORDER";

    case ABSTOL:                      return "ABSTOL";
    case RELTOL:                      return "RELTOL";
    case INIT_STEP_SIZE:              return "INIT_STEP_SIZE";
    case MAX_STEPS:                   return "MAX_STEPS";

    case GS_TYPE:                     return "GS_TYPE";
    case MAX_KRYLOV:                  return "MAX_KRYLOV";
    case MAX_RESTARTS:                return "MAX_RESTARTS";
    case SCHUR_SAFETY:                return "SCHUR_SAFETY";

    case PRINT_PROGRESS:              return "PRINT_PROGRESS";
    case PRINT_STATISTICS:            return "PRINT_STATISTICS";
    case PRINT_TIMING:                return "PRINT_TIMING";
    case PRINT_PARAMLIST:             return "PRINT_PARAMLIST";
    case PRINT_CONFIG:                return "PRINT_CONFIG";
    case USE_ANALYTIC_JACOBIAN:       return "USE_ANALYTIC_JACOBIAN";
    case WRITE_AT_USER_TIMES:         return "WRITE_AT_USER_TIMES";
    case WRITE_SOLUTION_TIMES:        return "WRITE_SOLUTION_TIMES";
    case WRITE_SOLUTION_COLUMN_OUTLET:return "WRITE_SOLUTION_COLUMN_OUTLET";
    case WRITE_SOLUTION_COLUMN_INLET: return "WRITE_SOLUTION_COLUMN_INLET";
    case WRITE_SENS_COLUMN_OUTLET:    return "WRITE_SENS_COLUMN_OUTLET";
    case WRITE_SOLUTION_ALL:          return "WRITE_SOLUTION_ALL";
    case WRITE_SENS_ALL:              return "WRITE_SENS_ALL";
    case LOG_LEVEL:                   return "LOG_LEVEL";
    case USER_SOLUTION_TIMES:         return "USER_SOLUTION_TIMES";
    case NTHREADS:                    return "NTHREADS";

    case NSENS:                       return "NSENS";
    case SENS_METHOD:                 return "SENS_METHOD";
    case SENS_NAME:                   return "SENS_NAME";
    case SENS_COMP:                   return "SENS_COMP";
    case SENS_SECTION:                return "SENS_SECTION";
    case SENS_ABSTOL:                 return "SENS_ABSTOL";
    case SENS_FD_DELTA:               return "SENS_FD_DELTA";

    case SOLUTION_TIMES:              return "SOLUTION_TIMES";
    case SOLUTION_COLUMN:             return "SOLUTION_COLUMN";
    case SOLUTION_PARTICLE:           return "SOLUTION_PARTICLE";
    case SOLUTION_BOUNDARY:           return "SOLUTION_BOUNDARY";
    case SOLUTION_COLUMN_OUTLET:      return "SOLUTION_COLUMN_OUTLET";
    case SOLUTION_COLUMN_INLET:       return "SOLUTION_COLUMN_INLET";
    case SENS_COLUMN:                 return "SENS_COLUMN";
    case SENS_PARTICLE:               return "SENS_PARTICLE";
    case SENS_BOUNDARY:               return "SENS_BOUNDARY";
    case SENS_COLUMN_OUTLET:          return "SENS_COLUMN_OUTLET";

    default:                          return "UNKNOWN_PARAMETER";
    }
}

inline const char* e2s(GroupName group)
{
    switch (group)
    {
    case GRP_IN:                      return "/input";

    case GRP_IN_MODEL:                return "/input/model";
    case GRP_IN_ADSORPTION:           return "/input/model/adsorption";
    case GRP_IN_INLET:                return "/input/model/inlet";
    case GRP_IN_EXTERNAL:             return "/input/model/external";

    case GRP_IN_DISCRETIZATION:       return "/input/discretization";
    case GRP_IN_WENO:                 return "/input/discretization/weno";

    case GRP_IN_SOLVER:               return "/input/solver";
    case GRP_IN_SCHUR:                return "/input/solver/schur_solver";
    case GRP_IN_TIME:                 return "/input/solver/time_integrator";


    case GRP_IN_SENSITIVITY:          return "/input/sensitivity";

    case GRP_OUT:                     return "/output";

    case GRP_OUT_SOLUTION:            return "/output/solution";
    case GRP_OUT_SENS:                return "/output/sensitivity";

    default:                          return "UNKNOWN_GROUP";
    }
}


// the conversion is case sensitive !
inline void s2e(const std::string& name, ChromatographyType& enum_ret)
{
    for (int enum_it = ChromatographyType_begin; enum_it != ChromatographyType_end; ++enum_it )
    {
        enum_ret = static_cast<ChromatographyType>(enum_it);
        if (name == e2s(enum_ret)) return;
    }
    throw CadetException("No enum entry for given chromatography type \"" + name + "\"!");
}

inline void s2e(const std::string& name, AdsorptionType& enum_ret)
{
    for (int enum_it = AdsorptionType_begin; enum_it != AdsorptionType_end; ++enum_it )
    {
        enum_ret = static_cast<AdsorptionType>(enum_it);
        if (name == e2s(enum_ret)) return;
    }
    throw CadetException("No enum entry for given adsorption type \"" + name + "\"!");
}

inline void s2e(const std::string& name, ReconstructionType& enum_ret)
{
    for (int enum_it = ReconstructionType_begin; enum_it != ReconstructionType_end; ++enum_it )
    {
        enum_ret = static_cast<ReconstructionType>(enum_it);
        if (name == e2s(enum_ret)) return;
    }
    throw CadetException("No enum entry for given reconstruction type \"" + name + "\"!");
}

inline void s2e(const std::string& name, ParticleDiscType& enum_ret)
{
    for (int enum_it = ParticleDiscType_begin; enum_it != ParticleDiscType_end; ++enum_it )
    {
        enum_ret = static_cast<ParticleDiscType>(enum_it);
        if (name == e2s(enum_ret)) return;
    }
    throw CadetException("No enum entry for given particle discretization type \"" + name + "\"!");
}

inline void s2e(const std::string& name, ParameterName& enum_ret)
{
    for (int enum_it = ParameterName_begin; enum_it != ParameterName_end; ++enum_it )
    {
        enum_ret = static_cast<ParameterName>(enum_it);
        if (name == e2s(enum_ret)) return;
    }
    throw CadetException("No enum entry for given parameter name \"" + name + "\"!");
}

inline void s2e(const std::string& name, OtherName& enum_ret)
{
    for (int enum_it = OtherName_begin; enum_it != OtherName_end; ++enum_it )
    {
        enum_ret = static_cast<OtherName>(enum_it);
        if (name == e2s(enum_ret)) return;
    }
    throw CadetException("No enum entry for given parameter name \"" + name + "\"!");
}

inline void s2e(const std::string& name, GroupName& enum_ret)
{
    for (int enum_it = GroupName_begin; enum_it != GroupName_end; ++enum_it )
    {
        enum_ret = static_cast<GroupName>(enum_it);
        if (name == e2s(enum_ret)) return;
    }
    throw CadetException("No enum entry for given group name \"" + name + "\"!");
}



} // namespace cadet

#endif /* LIBCADET_CADETENUMERATION_HPP_ */
