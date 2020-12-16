// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2020: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include <catch.hpp>

#include "BindingModelTests.hpp"
#include "BindingModels.hpp"

CADET_BINDINGTEST("LINEAR", "EXT_LINEAR", (1,1), (1,0,1), (1.0, 2.0, 0.0, 0.0), (1.0, 3.0, 2.0, 0.0, 0.0), \
	R"json( "LIN_KA": [1.0, 2.0],
	        "LIN_KD": [0.1, 0.2]
	)json", \
	R"json( "LIN_KA": [1.0, 3.0, 2.0],
	        "LIN_KD": [0.1, 0.3, 0.2]
	)json", \
	R"json( "EXT_LIN_KA": [0.0, 0.0],
	        "EXT_LIN_KA_T": [1.0, 2.0],
	        "EXT_LIN_KA_TT": [0.0, 0.0],
	        "EXT_LIN_KA_TTT": [0.0, 0.0],
	        "EXT_LIN_KD": [0.0, 0.0],
	        "EXT_LIN_KD_T": [0.1, 0.2],
	        "EXT_LIN_KD_TT": [0.0, 0.0],
	        "EXT_LIN_KD_TTT": [0.0, 0.0]
	)json", \
	R"json( "EXT_LIN_KA": [0.0, 0.0, 0.0],
	        "EXT_LIN_KA_T": [1.0, 3.0, 2.0],
	        "EXT_LIN_KA_TT": [0.0, 0.0, 0.0],
	        "EXT_LIN_KA_TTT": [0.0, 0.0, 0.0],
	        "EXT_LIN_KD": [0.0, 0.0, 0.0],
	        "EXT_LIN_KD_T": [0.1, 0.3, 0.2],
	        "EXT_LIN_KD_TT": [0.0, 0.0, 0.0],
	        "EXT_LIN_KD_TTT": [0.0, 0.0, 0.0]
	)json", \
	1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING);


CADET_BINDINGTEST("MULTI_COMPONENT_LANGMUIR", "EXT_MULTI_COMPONENT_LANGMUIR", (1,1), (1,0,1), (1.0, 2.0, 0.0, 0.0), (1.0, 3.0, 2.0, 0.0, 0.0), \
	R"json( "MCL_KA": [1.14, 2.0],
	        "MCL_KD": [0.004, 0.008],
	        "MCL_QMAX": [4.88, 3.5]
	)json", \
	R"json( "MCL_KA": [1.14, 1.0, 2.0],
	        "MCL_KD": [0.004, 2.0, 0.008],
	        "MCL_QMAX": [4.88, 3.0, 3.5]
	)json", \
	R"json( "EXT_MCL_KA": [0.0, 0.0],
	        "EXT_MCL_KA_T": [1.14, 2.0],
	        "EXT_MCL_KA_TT": [0.0, 0.0],
	        "EXT_MCL_KA_TTT": [0.0, 0.0],
	        "EXT_MCL_KD": [0.0, 0.0],
	        "EXT_MCL_KD_T": [0.004, 0.008],
	        "EXT_MCL_KD_TT": [0.0, 0.0],
	        "EXT_MCL_KD_TTT": [0.0, 0.0],
	        "EXT_MCL_QMAX": [0.0, 0.0],
	        "EXT_MCL_QMAX_T": [4.88, 3.5],
	        "EXT_MCL_QMAX_TT": [0.0, 0.0],
	        "EXT_MCL_QMAX_TTT": [0.0, 0.0]
	)json", \
	R"json( "EXT_MCL_KA": [0.0, 0.0, 0.0],
	        "EXT_MCL_KA_T": [1.14, 1.0, 2.0],
	        "EXT_MCL_KA_TT": [0.0, 0.0, 0.0],
	        "EXT_MCL_KA_TTT": [0.0, 0.0, 0.0],
	        "EXT_MCL_KD": [0.0, 0.0, 0.0],
	        "EXT_MCL_KD_T": [0.004, 2.0, 0.008],
	        "EXT_MCL_KD_TT": [0.0, 0.0, 0.0],
	        "EXT_MCL_KD_TTT": [0.0, 0.0, 0.0],
	        "EXT_MCL_QMAX": [0.0, 0.0, 0.0],
	        "EXT_MCL_QMAX_T": [4.88, 3.0, 3.5],
	        "EXT_MCL_QMAX_TT": [0.0, 0.0, 0.0],
	        "EXT_MCL_QMAX_TTT": [0.0, 0.0, 0.0]
	)json", \
	1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING);


CADET_BINDINGTEST("MULTI_COMPONENT_ANTILANGMUIR", "EXT_MULTI_COMPONENT_ANTILANGMUIR", (1,1), (1,0,1), (1.0, 2.0, 0.0, 0.0), (1.0, 3.0, 2.0, 0.0, 0.0), \
	R"json( "MCAL_KA": [1.14, 2.0],
	        "MCAL_KD": [0.004, 0.008],
	        "MCAL_QMAX": [4.88, 3.5],
	        "MCAL_ANTILANGMUIR": [1.0, -1.0]
	)json", \
	R"json( "MCAL_KA": [1.14, 1.0, 2.0],
	        "MCAL_KD": [0.004, 2.0, 0.008],
	        "MCAL_QMAX": [4.88, 3.0, 3.5],
	        "MCAL_ANTILANGMUIR": [1.0, 0.0, -1.0]
	)json", \
	R"json( "EXT_MCAL_KA": [0.0, 0.0],
	        "EXT_MCAL_KA_T": [1.14, 2.0],
	        "EXT_MCAL_KA_TT": [0.0, 0.0],
	        "EXT_MCAL_KA_TTT": [0.0, 0.0],
	        "EXT_MCAL_KD": [0.0, 0.0],
	        "EXT_MCAL_KD_T": [0.004, 0.008],
	        "EXT_MCAL_KD_TT": [0.0, 0.0],
	        "EXT_MCAL_KD_TTT": [0.0, 0.0],
	        "EXT_MCAL_QMAX": [0.0, 0.0],
	        "EXT_MCAL_QMAX_T": [4.88, 3.5],
	        "EXT_MCAL_QMAX_TT": [0.0, 0.0],
	        "EXT_MCAL_QMAX_TTT": [0.0, 0.0],
	        "EXT_MCAL_ANTILANGMUIR": [1.0, -1.0],
	        "EXT_MCAL_ANTILANGMUIR_T": [0.0, 0.0],
	        "EXT_MCAL_ANTILANGMUIR_TT": [0.0, 0.0],
	        "EXT_MCAL_ANTILANGMUIR_TTT": [0.0, 0.0]
	)json", \
	R"json( "EXT_MCAL_KA": [0.0, 0.0, 0.0],
	        "EXT_MCAL_KA_T": [1.14, 1.0, 2.0],
	        "EXT_MCAL_KA_TT": [0.0, 0.0, 0.0],
	        "EXT_MCAL_KA_TTT": [0.0, 0.0, 0.0],
	        "EXT_MCAL_KD": [0.0, 0.0, 0.0],
	        "EXT_MCAL_KD_T": [0.004, 2.0, 0.008],
	        "EXT_MCAL_KD_TT": [0.0, 0.0, 0.0],
	        "EXT_MCAL_KD_TTT": [0.0, 0.0, 0.0],
	        "EXT_MCAL_QMAX": [0.0, 0.0, 0.0],
	        "EXT_MCAL_QMAX_T": [4.88, 3.0, 3.5],
	        "EXT_MCAL_QMAX_TT": [0.0, 0.0, 0.0],
	        "EXT_MCAL_QMAX_TTT": [0.0, 0.0, 0.0],
	        "EXT_MCAL_ANTILANGMUIR": [1.0, 0.0, -1.0],
	        "EXT_MCAL_ANTILANGMUIR_T": [0.0, 0.0, 0.0],
	        "EXT_MCAL_ANTILANGMUIR_TT": [0.0, 0.0, 0.0],
	        "EXT_MCAL_ANTILANGMUIR_TTT": [0.0, 0.0, 0.0]
	)json", \
	1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING);


CADET_BINDINGTEST("MOBILE_PHASE_MODULATOR", "EXT_MOBILE_PHASE_MODULATOR", (1,1,1), (1,1,0,1), (1.2, 1.5, 2.0, 0.5, 1.5, 1.8), (1.2, 1.5, 0.0, 2.0, 0.5, 1.5, 1.8), \
	R"json( "MPM_KA": [0.0, 1.14, 2.0],
	        "MPM_KD": [0.0, 0.004, 0.008],
	        "MPM_QMAX": [0.0, 4.88, 3.5],
	        "MPM_GAMMA": [0.0, 0.5, -1.0],
	        "MPM_BETA": [0.0, 1.5, 2.0]
	)json", \
	R"json( "MPM_KA": [0.0, 1.14, 1.0, 2.0],
	        "MPM_KD": [0.0, 0.004, 2.0, 0.008],
	        "MPM_QMAX": [0.0, 4.88, 3.0, 3.5],
	        "MPM_GAMMA": [0.0, 0.5, 0.0, -1.0],
	        "MPM_BETA": [0.0, 1.5, 0.0, 2.0]
	)json", \
	R"json( "EXT_MPM_KA": [0.0, 0.0, 0.0],
	        "EXT_MPM_KA_T": [0.0, 1.14, 2.0],
	        "EXT_MPM_KA_TT": [0.0, 0.0, 0.0],
	        "EXT_MPM_KA_TTT": [0.0, 0.0, 0.0],
	        "EXT_MPM_KD": [0.0, 0.0, 0.0],
	        "EXT_MPM_KD_T": [0.0, 0.004, 0.008],
	        "EXT_MPM_KD_TT": [0.0, 0.0, 0.0],
	        "EXT_MPM_KD_TTT": [0.0, 0.0, 0.0],
	        "EXT_MPM_QMAX": [0.0, 0.0, 0.0],
	        "EXT_MPM_QMAX_T": [0.0, 4.88, 3.5],
	        "EXT_MPM_QMAX_TT": [0.0, 0.0, 0.0],
	        "EXT_MPM_QMAX_TTT": [0.0, 0.0, 0.0],
	        "EXT_MPM_GAMMA": [0.0, 0.0, 0.0],
	        "EXT_MPM_GAMMA_T": [0.0, 0.5, -1.0],
	        "EXT_MPM_GAMMA_TT": [0.0, 0.0, 0.0],
	        "EXT_MPM_GAMMA_TTT": [0.0, 0.0, 0.0],
	        "EXT_MPM_BETA": [0.0, 0.0, 0.0],
	        "EXT_MPM_BETA_T": [0.0, 1.5, 2.0],
	        "EXT_MPM_BETA_TT": [0.0, 0.0, 0.0],
	        "EXT_MPM_BETA_TTT": [0.0, 0.0, 0.0]
	)json", \
	R"json( "EXT_MPM_KA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MPM_KA_T": [0.0, 1.14, 1.0, 2.0],
	        "EXT_MPM_KA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MPM_KA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MPM_KD": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MPM_KD_T": [0.0, 0.004, 2.0, 0.008],
	        "EXT_MPM_KD_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MPM_KD_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MPM_QMAX": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MPM_QMAX_T": [0.0, 4.88, 3.0, 3.5],
	        "EXT_MPM_QMAX_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MPM_QMAX_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MPM_GAMMA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MPM_GAMMA_T": [0.0, 0.5, 0.0, -1.0],
	        "EXT_MPM_GAMMA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MPM_GAMMA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MPM_BETA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MPM_BETA_T": [0.0, 1.5, 0.0, 2.0],
	        "EXT_MPM_BETA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MPM_BETA_TTT": [0.0, 0.0, 0.0, 0.0]
	)json", \
	1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING);


CADET_BINDINGTEST_MULTI("EXTENDED_MOBILE_PHASE_MODULATOR", "EXT_EXTENDED_MOBILE_PHASE_MODULATOR", " all linear", (1,1,1), (1,1,0,1), (1.2, 1.5, 2.0, 0.5, 1.5, 1.8), (1.2, 1.5, 0.0, 2.0, 0.5, 1.5, 1.8), \
	R"json( "EMPM_KA": [0.0, 1.14, 2.0],
	        "EMPM_KD": [0.0, 0.004, 0.008],
	        "EMPM_QMAX": [0.0, 4.88, 3.5],
	        "EMPM_GAMMA": [0.0, 0.5, -1.0],
	        "EMPM_BETA": [0.0, 1.5, 2.0],
	        "EMPM_COMP_MODE": [1,1,1]
	)json", \
	R"json( "EMPM_KA": [0.0, 1.14, 1.0, 2.0],
	        "EMPM_KD": [0.0, 0.004, 2.0, 0.008],
	        "EMPM_QMAX": [0.0, 4.88, 3.0, 3.5],
	        "EMPM_GAMMA": [0.0, 0.5, 0.0, -1.0],
	        "EMPM_BETA": [0.0, 1.5, 0.0, 2.0],
	        "EMPM_COMP_MODE": [1,1,1,1]
	)json", \
	R"json( "EXT_EMPM_KA": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_T": [0.0, 1.14, 2.0],
	        "EXT_EMPM_KA_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KD": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_T": [0.0, 0.004, 0.008],
	        "EXT_EMPM_KD_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX": [0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_T": [0.0, 4.88, 3.5],
	        "EXT_EMPM_QMAX_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA": [0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_T": [0.0, 0.5, -1.0],
	        "EXT_EMPM_GAMMA_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA": [0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_T": [0.0, 1.5, 2.0],
	        "EXT_EMPM_BETA_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_TTT": [0.0, 0.0, 0.0],
	        "EMPM_COMP_MODE": [1,1,1]
	)json", \
	R"json( "EXT_EMPM_KA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_T": [0.0, 1.14, 1.0, 2.0],
	        "EXT_EMPM_KA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KD": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_T": [0.0, 0.004, 2.0, 0.008],
	        "EXT_EMPM_KD_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_T": [0.0, 4.88, 3.0, 3.5],
	        "EXT_EMPM_QMAX_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_T": [0.0, 0.5, 0.0, -1.0],
	        "EXT_EMPM_GAMMA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_T": [0.0, 1.5, 0.0, 2.0],
	        "EXT_EMPM_BETA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EMPM_COMP_MODE": [1,1,1,1]
	)json", \
	1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING);


CADET_BINDINGTEST_MULTI("EXTENDED_MOBILE_PHASE_MODULATOR", "EXT_EXTENDED_MOBILE_PHASE_MODULATOR", " all Langmuir", (1,1,1), (1,1,0,1), (1.2, 1.5, 2.0, 0.5, 1.5, 1.8), (1.2, 1.5, 0.0, 2.0, 0.5, 1.5, 1.8), \
	R"json( "EMPM_KA": [1.6, 1.14, 2.0],
	        "EMPM_KD": [0.006, 0.004, 0.008],
	        "EMPM_QMAX": [1.0, 4.88, 3.5],
	        "EMPM_GAMMA": [1.5, 0.5, -1.0],
	        "EMPM_BETA": [0.8, 1.5, 2.0],
	        "EMPM_COMP_MODE": [2,2,2]
	)json", \
	R"json( "EMPM_KA": [1.6, 1.14, 1.0, 2.0],
	        "EMPM_KD": [0.006, 0.004, 2.0, 0.008],
	        "EMPM_QMAX": [1.0, 4.88, 3.0, 3.5],
	        "EMPM_GAMMA": [1.5, 0.5, 0.0, -1.0],
	        "EMPM_BETA": [0.8, 1.5, 0.0, 2.0],
	        "EMPM_COMP_MODE": [2,2,2,2]
	)json", \
	R"json( "EXT_EMPM_KA": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_T": [1.6, 1.14, 2.0],
	        "EXT_EMPM_KA_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KD": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_T": [0.006, 0.004, 0.008],
	        "EXT_EMPM_KD_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX": [0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_T": [1.0, 4.88, 3.5],
	        "EXT_EMPM_QMAX_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA": [0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_T": [1.5, 0.5, -1.0],
	        "EXT_EMPM_GAMMA_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA": [0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_T": [0.8, 1.5, 2.0],
	        "EXT_EMPM_BETA_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_TTT": [0.0, 0.0, 0.0],
	        "EMPM_COMP_MODE": [2,2,2]
	)json", \
	R"json( "EXT_EMPM_KA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_T": [1.6, 1.14, 1.0, 2.0],
	        "EXT_EMPM_KA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KD": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_T": [0.006, 0.004, 2.0, 0.008],
	        "EXT_EMPM_KD_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_T": [1.0, 4.88, 3.0, 3.5],
	        "EXT_EMPM_QMAX_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_T": [1.5, 0.5, 0.0, -1.0],
	        "EXT_EMPM_GAMMA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_T": [0.8, 1.5, 0.0, 2.0],
	        "EXT_EMPM_BETA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EMPM_COMP_MODE": [2,2,2,2]
	)json", \
	1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING);


CADET_BINDINGTEST_MULTI("EXTENDED_MOBILE_PHASE_MODULATOR", "EXT_EXTENDED_MOBILE_PHASE_MODULATOR", " mixed", (1,1,1), (1,1,0,1), (1.2, 1.5, 2.0, 0.5, 1.5, 1.8), (1.2, 1.5, 0.0, 2.0, 0.5, 1.5, 1.8), \
	R"json( "EMPM_KA": [1.6, 1.14, 2.0],
	        "EMPM_KD": [0.006, 0.004, 0.008],
	        "EMPM_QMAX": [1.0, 4.88, 3.5],
	        "EMPM_GAMMA": [1.5, 0.5, -1.0],
	        "EMPM_BETA": [0.8, 1.5, 2.0],
	        "EMPM_COMP_MODE": [2,1,2]
	)json", \
	R"json( "EMPM_KA": [1.6, 1.14, 1.0, 2.0],
	        "EMPM_KD": [0.006, 0.004, 2.0, 0.008],
	        "EMPM_QMAX": [1.0, 4.88, 3.0, 3.5],
	        "EMPM_GAMMA": [1.5, 0.5, 0.0, -1.0],
	        "EMPM_BETA": [0.8, 1.5, 0.0, 2.0],
	        "EMPM_COMP_MODE": [2,1,2,2]
	)json", \
	R"json( "EXT_EMPM_KA": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_T": [1.6, 1.14, 2.0],
	        "EXT_EMPM_KA_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KD": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_T": [0.006, 0.004, 0.008],
	        "EXT_EMPM_KD_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX": [0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_T": [1.0, 4.88, 3.5],
	        "EXT_EMPM_QMAX_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA": [0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_T": [1.5, 0.5, -1.0],
	        "EXT_EMPM_GAMMA_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA": [0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_T": [0.8, 1.5, 2.0],
	        "EXT_EMPM_BETA_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_TTT": [0.0, 0.0, 0.0],
	        "EMPM_COMP_MODE": [2,1,2]
	)json", \
	R"json( "EXT_EMPM_KA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_T": [1.6, 1.14, 1.0, 2.0],
	        "EXT_EMPM_KA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KD": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_T": [0.006, 0.004, 2.0, 0.008],
	        "EXT_EMPM_KD_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_T": [1.0, 4.88, 3.0, 3.5],
	        "EXT_EMPM_QMAX_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_T": [1.5, 0.5, 0.0, -1.0],
	        "EXT_EMPM_GAMMA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_T": [0.8, 1.5, 0.0, 2.0],
	        "EXT_EMPM_BETA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EMPM_COMP_MODE": [2,1,2,2]
	)json", \
	1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING);


CADET_BINDINGTEST_MULTI("EXTENDED_MOBILE_PHASE_MODULATOR", "EXT_EXTENDED_MOBILE_PHASE_MODULATOR", " Langmuir modified", (1,1,1), (1,1,0,1), (1.2, 1.5, 2.0, 0.5, 1.5, 1.8), (1.2, 1.5, 0.0, 2.0, 0.5, 1.5, 1.8), \
	R"json( "EMPM_KA": [1.6, 1.14, 2.0],
	        "EMPM_KD": [0.006, 0.004, 0.008],
	        "EMPM_QMAX": [1.0, 4.88, 3.5],
	        "EMPM_GAMMA": [1.5, 0.5, -1.0],
	        "EMPM_BETA": [0.8, 1.5, 2.0],
	        "EMPM_COMP_MODE": [2,2,0]
	)json", \
	R"json( "EMPM_KA": [1.6, 1.14, 1.0, 2.0],
	        "EMPM_KD": [0.006, 0.004, 2.0, 0.008],
	        "EMPM_QMAX": [1.0, 4.88, 3.0, 3.5],
	        "EMPM_GAMMA": [1.5, 0.5, 0.0, -1.0],
	        "EMPM_BETA": [0.8, 1.5, 0.0, 2.0],
	        "EMPM_COMP_MODE": [2,2,2,0]
	)json", \
	R"json( "EXT_EMPM_KA": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_T": [1.6, 1.14, 2.0],
	        "EXT_EMPM_KA_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KD": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_T": [0.006, 0.004, 0.008],
	        "EXT_EMPM_KD_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX": [0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_T": [1.0, 4.88, 3.5],
	        "EXT_EMPM_QMAX_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA": [0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_T": [1.5, 0.5, -1.0],
	        "EXT_EMPM_GAMMA_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA": [0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_T": [0.8, 1.5, 2.0],
	        "EXT_EMPM_BETA_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_TTT": [0.0, 0.0, 0.0],
	        "EMPM_COMP_MODE": [2,2,0]
	)json", \
	R"json( "EXT_EMPM_KA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_T": [1.6, 1.14, 1.0, 2.0],
	        "EXT_EMPM_KA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KD": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_T": [0.006, 0.004, 2.0, 0.008],
	        "EXT_EMPM_KD_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_T": [1.0, 4.88, 3.0, 3.5],
	        "EXT_EMPM_QMAX_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_T": [1.5, 0.5, 0.0, -1.0],
	        "EXT_EMPM_GAMMA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_T": [0.8, 1.5, 0.0, 2.0],
	        "EXT_EMPM_BETA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EMPM_COMP_MODE": [2,2,2,0]
	)json", \
	1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING);


CADET_BINDINGTEST_MULTI("EXTENDED_MOBILE_PHASE_MODULATOR", "EXT_EXTENDED_MOBILE_PHASE_MODULATOR", " mixed modified", (1,1,1), (1,1,0,1), (1.2, 1.5, 2.0, 0.5, 1.5, 1.8), (1.2, 1.5, 0.0, 2.0, 0.5, 1.5, 1.8), \
	R"json( "EMPM_KA": [1.6, 1.14, 2.0],
	        "EMPM_KD": [0.006, 0.004, 0.008],
	        "EMPM_QMAX": [1.0, 4.88, 3.5],
	        "EMPM_GAMMA": [1.5, 0.5, -1.0],
	        "EMPM_BETA": [0.8, 1.5, 2.0],
	        "EMPM_COMP_MODE": [2,1,0]
	)json", \
	R"json( "EMPM_KA": [1.6, 1.14, 1.0, 2.0],
	        "EMPM_KD": [0.006, 0.004, 2.0, 0.008],
	        "EMPM_QMAX": [1.0, 4.88, 3.0, 3.5],
	        "EMPM_GAMMA": [1.5, 0.5, 0.0, -1.0],
	        "EMPM_BETA": [0.8, 1.5, 0.0, 2.0],
	        "EMPM_COMP_MODE": [2,1,2,0]
	)json", \
	R"json( "EXT_EMPM_KA": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_T": [1.6, 1.14, 2.0],
	        "EXT_EMPM_KA_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KD": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_T": [0.006, 0.004, 0.008],
	        "EXT_EMPM_KD_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX": [0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_T": [1.0, 4.88, 3.5],
	        "EXT_EMPM_QMAX_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA": [0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_T": [1.5, 0.5, -1.0],
	        "EXT_EMPM_GAMMA_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_TTT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA": [0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_T": [0.8, 1.5, 2.0],
	        "EXT_EMPM_BETA_TT": [0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_TTT": [0.0, 0.0, 0.0],
	        "EMPM_COMP_MODE": [2,1,0]
	)json", \
	R"json( "EXT_EMPM_KA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_T": [1.6, 1.14, 1.0, 2.0],
	        "EXT_EMPM_KA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KD": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_T": [0.006, 0.004, 2.0, 0.008],
	        "EXT_EMPM_KD_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_KD_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_T": [1.0, 4.88, 3.0, 3.5],
	        "EXT_EMPM_QMAX_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_QMAX_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_T": [1.5, 0.5, 0.0, -1.0],
	        "EXT_EMPM_GAMMA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_GAMMA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_T": [0.8, 1.5, 0.0, 2.0],
	        "EXT_EMPM_BETA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_EMPM_BETA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EMPM_COMP_MODE": [2,1,2,0]
	)json", \
	1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING);


CADET_BINDINGTEST("KUMAR_MULTI_COMPONENT_LANGMUIR", "EXT_KUMAR_MULTI_COMPONENT_LANGMUIR", (0,1,1), (0,1,0,1), (1.2, 1.5, 2.0, 0.75, 1.5), (1.2, 1.5, 0.0, 2.0, 0.75, 1.5), \
	R"json( "KMCL_KA": [0.0, 1.14, 2.0],
	        "KMCL_KD": [0.0, 0.004, 0.008],
	        "KMCL_QMAX": [0.0, 4.88, 3.5],
	        "KMCL_TEMP": 0.5,
	        "KMCL_NU": [0.0, 1.2, 2.2],
	        "KMCL_KACT": [0.0, 1.5, 2.0]
	)json", \
	R"json( "KMCL_KA": [0.0, 1.14, 1.0, 2.0],
	        "KMCL_KD": [0.0, 0.004, 2.0, 0.008],
	        "KMCL_QMAX": [0.0, 4.88, 3.0, 3.5],
	        "KMCL_TEMP": 0.5,
	        "KMCL_NU": [0.0, 1.2, 0.0, 2.2],
	        "KMCL_KACT": [0.0, 1.5, 0.0, 2.0]
	)json", \
	R"json( "EXT_KMCL_KA": [0.0, 0.0, 0.0],
	        "EXT_KMCL_KA_T": [0.0, 1.14, 2.0],
	        "EXT_KMCL_KA_TT": [0.0, 0.0, 0.0],
	        "EXT_KMCL_KA_TTT": [0.0, 0.0, 0.0],
	        "EXT_KMCL_KD": [0.0, 0.0, 0.0],
	        "EXT_KMCL_KD_T": [0.0, 0.004, 0.008],
	        "EXT_KMCL_KD_TT": [0.0, 0.0, 0.0],
	        "EXT_KMCL_KD_TTT": [0.0, 0.0, 0.0],
	        "EXT_KMCL_TEMP": 0.0,
	        "EXT_KMCL_TEMP_T": 0.5,
	        "EXT_KMCL_TEMP_TT": 0.0,
	        "EXT_KMCL_TEMP_TTT": 0.0,
	        "EXT_KMCL_QMAX": [0.0, 0.0, 0.0],
	        "EXT_KMCL_QMAX_T": [0.0, 4.88, 3.5],
	        "EXT_KMCL_QMAX_TT": [0.0, 0.0, 0.0],
	        "EXT_KMCL_QMAX_TTT": [0.0, 0.0, 0.0],
	        "EXT_KMCL_NU": [0.0, 0.0, 0.0],
	        "EXT_KMCL_NU_T": [0.0, 1.2, 2.2],
	        "EXT_KMCL_NU_TT": [0.0, 0.0, 0.0],
	        "EXT_KMCL_NU_TTT": [0.0, 0.0, 0.0],
	        "EXT_KMCL_KACT": [0.0, 0.0, 0.0],
	        "EXT_KMCL_KACT_T": [0.0, 1.5, 2.0],
	        "EXT_KMCL_KACT_TT": [0.0, 0.0, 0.0],
	        "EXT_KMCL_KACT_TTT": [0.0, 0.0, 0.0]
	)json", \
	R"json( "EXT_KMCL_KA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_KMCL_KA_T": [0.0, 1.14, 1.0, 2.0],
	        "EXT_KMCL_KA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_KMCL_KA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_KMCL_KD": [0.0, 0.0, 0.0, 0.0],
	        "EXT_KMCL_KD_T": [0.0, 0.004, 2.0, 0.008],
	        "EXT_KMCL_KD_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_KMCL_KD_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_KMCL_QMAX": [0.0, 0.0, 0.0, 0.0],
	        "EXT_KMCL_QMAX_T": [0.0, 4.88, 3.0, 3.5],
	        "EXT_KMCL_QMAX_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_KMCL_QMAX_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_KMCL_TEMP": 0.0,
	        "EXT_KMCL_TEMP_T": 0.5,
	        "EXT_KMCL_TEMP_TT": 0.0,
	        "EXT_KMCL_TEMP_TTT": 0.0,
	        "EXT_KMCL_NU": [0.0, 0.0, 0.0, 0.0],
	        "EXT_KMCL_NU_T": [0.0, 1.2, 0.0, 2.2],
	        "EXT_KMCL_NU_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_KMCL_NU_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_KMCL_KACT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_KMCL_KACT_T": [0.0, 1.5, 0.0, 2.0],
	        "EXT_KMCL_KACT_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_KMCL_KACT_TTT": [0.0, 0.0, 0.0, 0.0]
	)json", \
	1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_USED, CADET_DONT_COMPARE_BINDING_VS_NONBINDING);
// Note that we cannot enable binding vs non-binding test since the first component has to be non-binding,
// the liquid phase component matters in the Jacobian is not a column to ignore.


CADET_BINDINGTEST_ALLBINDING("MULTI_COMPONENT_BILANGMUIR", "EXT_MULTI_COMPONENT_BILANGMUIR", (2,2), (1.0, 2.0, 0.0, 0.0, 0.0, 0.0), \
	R"json( "MCBL_KA": [1.14, 2.0, 2.28, 4.0],
	        "MCBL_KD": [0.004, 0.008, 0.002, 0.003],
	        "MCBL_QMAX": [4.88, 3.5, 3.88, 2.5]
	)json", \
	R"json( "EXT_MCBL_KA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCBL_KA_T": [1.14, 2.0, 2.28, 4.0],
	        "EXT_MCBL_KA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCBL_KA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCBL_KD": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCBL_KD_T": [0.004, 0.008, 0.002, 0.003],
	        "EXT_MCBL_KD_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCBL_KD_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCBL_QMAX": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCBL_QMAX_T": [4.88, 3.5, 3.88, 2.5],
	        "EXT_MCBL_QMAX_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCBL_QMAX_TTT": [0.0, 0.0, 0.0, 0.0]
	)json", \
	1e-10, 1e-10);


CADET_BINDINGTEST("STERIC_MASS_ACTION", "EXT_STERIC_MASS_ACTION", (1,1,1), (1,1,0,1), (1.2, 2.0, 1.5, 80.0, 3.5, 2.7), (1.2, 2.0, 3.0, 1.5, 80.0, 3.5, 2.7), \
	R"json( "SMA_KA": [0.0, 3.55, 1.59],
	        "SMA_KD": [0.0, 10.0, 10.0],
	        "SMA_NU": [1.5, 2.0, 1.5],
	        "SMA_SIGMA": [0.0, 11.83, 10.6],
	        "SMA_LAMBDA": 100.0
	)json", \
	R"json( "SMA_KA": [0.0, 3.55, 7.7, 1.59],
	        "SMA_KD": [0.0, 10.0, 10.0, 10.0],
	        "SMA_NU": [1.5, 2.0, 3.7, 1.5],
	        "SMA_SIGMA": [0.0, 11.83, 10.0, 10.6],
	        "SMA_LAMBDA": 100.0
	)json", \
	R"json( "EXT_SMA_KA": [0.0, 0.0, 0.0],
	        "EXT_SMA_KA_T": [0.0, 3.55, 1.59],
	        "EXT_SMA_KA_TT": [0.0, 0.0, 0.0],
	        "EXT_SMA_KA_TTT": [0.0, 0.0, 0.0],
	        "EXT_SMA_KD": [0.0, 0.0, 0.0],
	        "EXT_SMA_KD_T": [0.0, 10.0, 10.0],
	        "EXT_SMA_KD_TT": [0.0, 0.0, 0.0],
	        "EXT_SMA_KD_TTT": [0.0, 0.0, 0.0],
	        "EXT_SMA_NU": [1.5, 0.0, 0.0],
	        "EXT_SMA_NU_T": [0.0, 2.0, 1.5],
	        "EXT_SMA_NU_TT": [0.0, 0.0, 0.0],
	        "EXT_SMA_NU_TTT": [0.0, 0.0, 0.0],
	        "EXT_SMA_SIGMA": [0.0, 0.0, 0.0],
	        "EXT_SMA_SIGMA_T": [0.0, 11.83, 10.6],
	        "EXT_SMA_SIGMA_TT": [0.0, 0.0, 0.0],
	        "EXT_SMA_SIGMA_TTT": [0.0, 0.0, 0.0],
	        "EXT_SMA_LAMBDA": 0.0,
	        "EXT_SMA_LAMBDA_T": 100.0,
	        "EXT_SMA_LAMBDA_TT": 0.0,
	        "EXT_SMA_LAMBDA_TTT": 0.0
	)json", \
	R"json( "EXT_SMA_KA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SMA_KA_T": [0.0, 3.55, 7.7, 1.59],
	        "EXT_SMA_KA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SMA_KA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SMA_KD": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SMA_KD_T": [0.0, 10.0, 10.0, 10.0],
	        "EXT_SMA_KD_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SMA_KD_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SMA_NU": [1.5, 0.0, 0.0, 0.0],
	        "EXT_SMA_NU_T": [0.0, 2.0, 3.7, 1.5],
	        "EXT_SMA_NU_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SMA_NU_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SMA_SIGMA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SMA_SIGMA_T": [0.0, 11.83, 10.0, 10.6],
	        "EXT_SMA_SIGMA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SMA_SIGMA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SMA_LAMBDA": 0.0,
	        "EXT_SMA_LAMBDA_T": 100.0,
	        "EXT_SMA_LAMBDA_TT": 0.0,
	        "EXT_SMA_LAMBDA_TTT": 0.0
	)json", \
	1e-10, 1e-8, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING);


CADET_BINDINGTEST_ALLBINDING("BI_STERIC_MASS_ACTION", "EXT_BI_STERIC_MASS_ACTION", (2,2,2), (1.2, 2.0, 1.5, 80.0, 80.0, 3.5, 3.5, 2.7, 2.7), \
	R"json( "BISMA_KA": [0.0, 3.55, 1.59, 0.0, 3.55, 1.59],
	        "BISMA_KD": [0.0, 10.0, 10.0, 0.0, 10.0, 10.0],
	        "BISMA_NU": [1.2, 2.0, 1.5, 1.5, 2.0, 1.5],
	        "BISMA_SIGMA": [0.0, 11.83, 10.6, 0.0, 11.83, 10.6],
	        "BISMA_LAMBDA": [100.0, 100.0]
	)json", \
	R"json( "EXT_BISMA_KA": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_BISMA_KA_T": [0.0, 3.55, 1.59, 0.0, 3.55, 1.59],
	        "EXT_BISMA_KA_TT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_BISMA_KA_TTT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_BISMA_KD": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_BISMA_KD_T": [0.0, 10.0, 10.0, 0.0, 10.0, 10.0],
	        "EXT_BISMA_KD_TT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_BISMA_KD_TTT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_BISMA_NU": [1.2, 0.0, 0.0, 1.5, 0.0, 0.0],
	        "EXT_BISMA_NU_T": [0.0, 2.0, 1.5, 0.0, 2.0, 1.5],
	        "EXT_BISMA_NU_TT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_BISMA_NU_TTT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_BISMA_SIGMA": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_BISMA_SIGMA_T": [0.0, 11.83, 10.6, 0.0, 11.83, 10.6],
	        "EXT_BISMA_SIGMA_TT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_BISMA_SIGMA_TTT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_BISMA_LAMBDA": [0.0, 0.0],
	        "EXT_BISMA_LAMBDA_T": [100.0, 100.0],
	        "EXT_BISMA_LAMBDA_TT": [0.0, 0.0],
	        "EXT_BISMA_LAMBDA_TTT": [0.0, 0.0]
	)json", \
	1e-10, 1e-8);


CADET_BINDINGTEST("SELF_ASSOCIATION", "EXT_SELF_ASSOCIATION", (1,1,1), (1,1,0,1), (1.2, 2.0, 1.5, 80.0, 3.5, 2.7), (1.2, 2.0, 3.0, 1.5, 80.0, 3.5, 2.7), \
	R"json( "SAI_KA1": [0.0, 3.55, 1.59],
	        "SAI_KA2": [0.0, 1.5, 2.5],
	        "SAI_KD": [0.0, 10.0, 10.0],
	        "SAI_NU": [1.5, 2.0, 1.5],
	        "SAI_SIGMA": [0.0, 11.83, 10.6],
	        "SAI_LAMBDA": 100.0
	)json", \
	R"json( "SAI_KA1": [0.0, 3.55, 7.7, 1.59],
	        "SAI_KA2": [0.0, 1.5, 2.0, 2.5],
	        "SAI_KD": [0.0, 10.0, 10.0, 10.0],
	        "SAI_NU": [1.5, 2.0, 3.7, 1.5],
	        "SAI_SIGMA": [0.0, 11.83, 10.0, 10.6],
	        "SAI_LAMBDA": 100.0
	)json", \
	R"json( "EXT_SAI_KA1": [0.0, 0.0, 0.0],
	        "EXT_SAI_KA1_T": [0.0, 3.55, 1.59],
	        "EXT_SAI_KA1_TT": [0.0, 0.0, 0.0],
	        "EXT_SAI_KA1_TTT": [0.0, 0.0, 0.0],
	        "EXT_SAI_KA2": [0.0, 0.0, 0.0],
	        "EXT_SAI_KA2_T": [0.0, 1.5, 2.5],
	        "EXT_SAI_KA2_TT": [0.0, 0.0, 0.0],
	        "EXT_SAI_KA2_TTT": [0.0, 0.0, 0.0],
	        "EXT_SAI_KD": [0.0, 0.0, 0.0],
	        "EXT_SAI_KD_T": [0.0, 10.0, 10.0],
	        "EXT_SAI_KD_TT": [0.0, 0.0, 0.0],
	        "EXT_SAI_KD_TTT": [0.0, 0.0, 0.0],
	        "EXT_SAI_NU": [1.5, 0.0, 0.0],
	        "EXT_SAI_NU_T": [0.0, 2.0, 1.5],
	        "EXT_SAI_NU_TT": [0.0, 0.0, 0.0],
	        "EXT_SAI_NU_TTT": [0.0, 0.0, 0.0],
	        "EXT_SAI_SIGMA": [0.0, 0.0, 0.0],
	        "EXT_SAI_SIGMA_T": [0.0, 11.83, 10.6],
	        "EXT_SAI_SIGMA_TT": [0.0, 0.0, 0.0],
	        "EXT_SAI_SIGMA_TTT": [0.0, 0.0, 0.0],
	        "EXT_SAI_LAMBDA": 0.0,
	        "EXT_SAI_LAMBDA_T": 100.0,
	        "EXT_SAI_LAMBDA_TT": 0.0,
	        "EXT_SAI_LAMBDA_TTT": 0.0
	)json", \
	R"json( "EXT_SAI_KA1": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SAI_KA1_T": [0.0, 3.55, 7.7, 1.59],
	        "EXT_SAI_KA1_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SAI_KA1_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SAI_KA2": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SAI_KA2_T": [0.0, 1.5, 2.0, 2.5],
	        "EXT_SAI_KA2_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SAI_KA2_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SAI_KD": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SAI_KD_T": [0.0, 10.0, 10.0, 10.0],
	        "EXT_SAI_KD_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SAI_KD_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SAI_NU": [1.5, 0.0, 0.0, 0.0],
	        "EXT_SAI_NU_T": [0.0, 2.0, 3.7, 1.5],
	        "EXT_SAI_NU_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SAI_NU_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SAI_SIGMA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SAI_SIGMA_T": [0.0, 11.83, 10.0, 10.6],
	        "EXT_SAI_SIGMA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SAI_SIGMA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SAI_LAMBDA": 0.0,
	        "EXT_SAI_LAMBDA_T": 100.0,
	        "EXT_SAI_LAMBDA_TT": 0.0,
	        "EXT_SAI_LAMBDA_TTT": 0.0
	)json", \
	1e-10, 1e-4, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING);


CADET_BINDINGTEST("SASKA", "EXT_SASKA", (1,1), (1,0,1), (1.0, 2.0, 0.0, 0.0), (1.0, 3.0, 2.0, 0.0, 0.0), \
	R"json( "SASKA_H": [1.5, 2.0],
	        "SASKA_K": [0.1, 0.3, 0.4, 0.5]
	)json", \
	R"json( "SASKA_H": [1.5, 3.0, 2.0],
	        "SASKA_K": [0.1, 0.0, 0.3, 0.0, 0.0, 0.0, 0.4, 0.0, 0.5]
	)json", \
	R"json( "EXT_SASKA_H": [0.0, 0.0],
	        "EXT_SASKA_H_T": [1.5, 2.0],
	        "EXT_SASKA_H_TT": [0.0, 0.0],
	        "EXT_SASKA_H_TTT": [0.0, 0.0],
	        "EXT_SASKA_K": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SASKA_K_T": [0.1, 0.3, 0.4, 0.5],
	        "EXT_SASKA_K_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_SASKA_K_TTT": [0.0, 0.0, 0.0, 0.0]
	)json", \
	R"json( "EXT_SASKA_H": [0.0, 0.0, 0.0],
	        "EXT_SASKA_H_T": [1.5, 3.0, 2.0],
	        "EXT_SASKA_H_TT": [0.0, 0.0, 0.0],
	        "EXT_SASKA_H_TTT": [0.0, 0.0, 0.0],
	        "EXT_SASKA_K": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_SASKA_K_T": [0.1, 0.0, 0.3, 0.0, 0.0, 0.0, 0.4, 0.0, 0.5],
	        "EXT_SASKA_K_TT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_SASKA_K_TTT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	)json", \
	1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_USED, CADET_COMPARE_BINDING_VS_NONBINDING);


CADET_BINDINGTEST("MULTI_COMPONENT_SPREADING", "EXT_MULTI_COMPONENT_SPREADING", (2,2), (2,0,2), (1.2, 1.5, 0.1, 0.2, 0.3, 0.4), (1.2, 0.5, 1.5, 0.1, 0.2, 0.3, 0.4), \
	R"json( "MCSPR_KA": [1.14, 2.0, 1.5, 1.9],
	        "MCSPR_KD": [0.004, 0.008, 0.006, 0.002],
	        "MCSPR_QMAX": [4.88, 3.5, 4.5, 3.7],
	        "MCSPR_K12": [0.5, 1.1, 0.7, 1.2],
	        "MCSPR_K21": [0.6, 1.5, 2.0, 0.8]
	)json", \
	R"json( "MCSPR_KA": [1.14, 0.0, 2.0, 1.5, 0.0, 1.9],
	        "MCSPR_KD": [0.004, 0.0, 0.008, 0.006, 0.0, 0.002],
	        "MCSPR_QMAX": [4.88, 0.0, 3.5, 4.5, 0.0, 3.7],
	        "MCSPR_K12": [0.5, 0.0, 1.1, 0.7, 0.0, 1.2],
	        "MCSPR_K21": [0.6, 0.0, 1.5, 2.0, 0.0, 0.8]
	)json", \
	R"json( "EXT_MCSPR_KA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_KA_T": [1.14, 2.0, 1.5, 1.9],
	        "EXT_MCSPR_KA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_KA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_KD": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_KD_T": [0.004, 0.008, 0.006, 0.002],
	        "EXT_MCSPR_KD_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_KD_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_QMAX": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_QMAX_T": [4.88, 3.5, 4.5, 3.7],
	        "EXT_MCSPR_QMAX_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_QMAX_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_K12": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_K12_T": [0.5, 1.1, 0.7, 1.2],
	        "EXT_MCSPR_K12_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_K12_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_K21": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_K21_T": [0.6, 1.5, 2.0, 0.8],
	        "EXT_MCSPR_K21_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_K21_TTT": [0.0, 0.0, 0.0, 0.0]
	)json", \
	R"json( "EXT_MCSPR_KA": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_KA_T": [1.14, 0.0, 2.0, 1.5, 0.0, 1.9],
	        "EXT_MCSPR_KA_TT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_KA_TTT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_KD": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_KD_T": [0.004, 0.0, 0.008, 0.006, 0.0, 0.002],
	        "EXT_MCSPR_KD_TT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_KD_TTT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_QMAX": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_QMAX_T": [4.88, 0.0, 3.5, 4.5, 0.0, 3.7],
	        "EXT_MCSPR_QMAX_TT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_QMAX_TTT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_K12": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_K12_T": [0.5, 0.0, 1.1, 0.7, 0.0, 1.2],
	        "EXT_MCSPR_K12_TT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_K12_TTT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_K21": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_K21_T": [0.6, 0.0, 1.5, 2.0, 0.0, 0.8],
	        "EXT_MCSPR_K21_TT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MCSPR_K21_TTT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	)json", \
	1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING);


// @todo Find a better starting point for consistent initialization
CADET_BINDINGTEST("MULTISTATE_STERIC_MASS_ACTION", "EXT_MULTISTATE_STERIC_MASS_ACTION", (1,2,1), (1,2,0,1), (1.2, 2.0, 1.5, 80.0, 2.5, 2.7, 1.8), (1.2, 2.0, 3.0, 1.5, 80.0, 2.5, 2.7, 1.8), \
	R"json( "MSSMA_KA": [0.0, 3.55, 1.59, 4.0],
	        "MSSMA_KD": [0.0, 10.0, 10.0, 10.0],
	        "MSSMA_NU": [1.2, 1.5, 2.0, 1.9],
	        "MSSMA_SIGMA": [0.0, 11.83, 10.6, 9.8],
	        "MSSMA_RATES": [0.0, 0.9, 0.8, 1.2, 1.1, 1.4],
	        "MSSMA_LAMBDA": 100.0
	)json", \
	R"json( "MSSMA_KA": [0.0, 3.55, 1.59, 4.0],
	        "MSSMA_KD": [0.0, 10.0, 10.0, 10.0],
	        "MSSMA_NU": [1.2, 1.5, 2.0, 1.9],
	        "MSSMA_SIGMA": [0.0, 11.83, 10.6, 9.8],
	        "MSSMA_RATES": [0.0, 0.9, 0.8, 1.2, 1.1, 1.4],
	        "MSSMA_LAMBDA": 100.0
	)json", \
	R"json( "EXT_MSSMA_KA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_KA_T": [0.0, 3.55, 1.59, 4.0],
	        "EXT_MSSMA_KA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_KA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_KD": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_KD_T": [0.0, 10.0, 10.0, 10.0, 0.0],
	        "EXT_MSSMA_KD_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_KD_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_NU": [1.2, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_NU_T": [0.0, 1.5, 2.0, 1.9],
	        "EXT_MSSMA_NU_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_NU_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_SIGMA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_SIGMA_T": [0.0, 11.83, 10.6, 9.8],
	        "EXT_MSSMA_SIGMA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_SIGMA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_RATES": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_RATES_T": [0.0, 0.9, 0.8, 1.2, 1.1, 1.4],
	        "EXT_MSSMA_RATES_TT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_RATES_TTT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_LAMBDA": 0.0,
	        "EXT_MSSMA_LAMBDA_T": 100.0,
	        "EXT_MSSMA_LAMBDA_TT": 0.0,
	        "EXT_MSSMA_LAMBDA_TTT": 0.0
	)json", \
	R"json( "EXT_MSSMA_KA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_KA_T": [0.0, 3.55, 1.59, 4.0],
	        "EXT_MSSMA_KA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_KA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_KD": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_KD_T": [0.0, 10.0, 10.0, 10.0, 0.0],
	        "EXT_MSSMA_KD_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_KD_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_NU": [1.2, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_NU_T": [0.0, 1.5, 2.0, 1.9],
	        "EXT_MSSMA_NU_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_NU_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_SIGMA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_SIGMA_T": [0.0, 11.83, 10.6, 9.8],
	        "EXT_MSSMA_SIGMA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_SIGMA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_RATES": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_RATES_T": [0.0, 0.9, 0.8, 1.2, 1.1, 1.4],
	        "EXT_MSSMA_RATES_TT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_RATES_TTT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_MSSMA_LAMBDA": 0.0,
	        "EXT_MSSMA_LAMBDA_T": 100.0,
	        "EXT_MSSMA_LAMBDA_TT": 0.0,
	        "EXT_MSSMA_LAMBDA_TTT": 0.0
	)json", \
	1e-10, 5e-2, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING);


CADET_BINDINGTEST_SINGLE("SIMPLIFIED_MULTISTATE_STERIC_MASS_ACTION", (1,2,1), (1,2,0,1), (1.2, 2.0, 1.5, 80.0, 2.5, 2.7, 1.8), (1.2, 2.0, 3.0, 1.5, 80.0, 2.5, 2.7, 1.8), \
	R"json( "SMSSMA_KA": [0.0, 3.55, 1.59, 4.0],
	        "SMSSMA_KD": [0.0, 10.0, 10.0, 10.0],
	        "SMSSMA_NU_MIN": [0.0, 1.5, 1.9],
	        "SMSSMA_NU_MAX": [0.0, 2.0, 2.1],
	        "SMSSMA_NU_QUAD": [0.0, 0.1, 0.1],
	        "SMSSMA_SIGMA_MIN": [0.0, 10.6, 9.8],
	        "SMSSMA_SIGMA_MAX": [0.0, 11.83, 10.1],
	        "SMSSMA_SIGMA_QUAD": [0.0, -0.1, -0.8],
	        "SMSSMA_KSW": [0.0, 0.9, 1.4],
	        "SMSSMA_KSW_LIN": [0.0, 0.1, 0.1],
	        "SMSSMA_KSW_QUAD": [0.0, -0.2, -0.1],
	        "SMSSMA_KWS": [0.0, 0.8, 1.1],
	        "SMSSMA_KWS_LIN": [0.0, 0.2, 0.3],
	        "SMSSMA_KWS_QUAD": [0.0, -0.1, -0.1],
	        "SMSSMA_LAMBDA": 100.0,
	        "SMSSMA_NU_SALT": 1.2
	)json", \
	R"json( "SMSSMA_KA": [0.0, 3.55, 1.59, 4.0],
	        "SMSSMA_KD": [0.0, 10.0, 10.0, 10.0],
	        "SMSSMA_NU_MIN": [0.0, 1.5, 0.0, 1.9],
	        "SMSSMA_NU_MAX": [0.0, 2.0, 0.0, 2.1],
	        "SMSSMA_NU_QUAD": [0.0, 0.1, 0.0, 0.1],
	        "SMSSMA_SIGMA_MIN": [0.0, 10.6, 0.0, 9.8],
	        "SMSSMA_SIGMA_MAX": [0.0, 11.83, 0.0, 10.1],
	        "SMSSMA_SIGMA_QUAD": [0.0, -0.1, 0.0, -0.8],
	        "SMSSMA_KSW": [0.0, 0.9, 0.0, 1.4],
	        "SMSSMA_KSW_LIN": [0.0, 0.1, 0.0, 0.1],
	        "SMSSMA_KSW_QUAD": [0.0, -0.2, 0.0, -0.1],
	        "SMSSMA_KWS": [0.0, 0.8, 0.0, 1.1],
	        "SMSSMA_KWS_LIN": [0.0, 0.2, 0.0, 0.3],
	        "SMSSMA_KWS_QUAD": [0.0, -0.1, 0.0, -0.1],
	        "SMSSMA_LAMBDA": 100.0,
	        "SMSSMA_NU_SALT": 1.2
	)json", \
	1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING);

CADET_BINDINGTEST("GENERALIZED_ION_EXCHANGE", "EXT_GENERALIZED_ION_EXCHANGE", (1,0,1,1), (1,0,1,0,1), (1.2, 0.5, 2.0, 1.5, 80.0, 3.5, 2.7), (1.2, 0.5, 2.0, 3.0, 1.5, 80.0, 3.5, 2.7), \
	R"json( "GIEX_KA": [0.0, 0.0, 3.55, 1.59],
	        "GIEX_KA_LIN": [0.0, 0.0, 0.8, 0.9],
	        "GIEX_KA_QUAD": [0.0, 0.0, 0.5, 0.6],
	        "GIEX_KA_SALT": [0.0, 0.0, 1.2, 0.8],
	        "GIEX_KA_PROT": [0.0, 0.0, 1.1, 0.9],
	        "GIEX_KD": [0.0, 0.0, 10.0, 10.0],
	        "GIEX_KD_LIN": [0.0, 0.0, 0.6, 0.5],
	        "GIEX_KD_QUAD": [0.0, 0.0, 0.9, 0.7],
	        "GIEX_KD_SALT": [0.0, 0.0, 1.0, 1.2],
	        "GIEX_KD_PROT": [0.0, 0.0, 0.9, 1.1],
	        "GIEX_NU": [1.2, 0.0, 2.0, 1.5],
	        "GIEX_NU_LIN": [0.0, 0.0, 0.5, 1.1],
	        "GIEX_NU_QUAD": [0.0, 0.0, 0.8, 0.9],
	        "GIEX_SIGMA": [0.0, 0.0, 11.83, 10.6],
	        "GIEX_LAMBDA": 100.0
	)json", \
	R"json( "GIEX_KA": [0.0, 0.0, 3.55, 7.7, 1.59],
	        "GIEX_KA_LIN": [0.0, 0.0, 0.8, 1.2, 0.9],
	        "GIEX_KA_QUAD": [0.0, 0.0, 0.5, 3.2, 0.6],
	        "GIEX_KA_SALT": [0.0, 0.0, 1.2, 4.1, 0.8],
	        "GIEX_KA_PROT": [0.0, 0.0, 1.1, 6.1, 0.9],
	        "GIEX_KD": [0.0, 0.0, 10.0, 10.0, 10.0],
	        "GIEX_KD_LIN": [0.0, 0.0, 0.6, 0.9, 0.5],
	        "GIEX_KD_QUAD": [0.0, 0.0, 0.9, 2.1, 0.7],
	        "GIEX_KD_SALT": [0.0, 0.0, 1.0, 3.7, 1.2],
	        "GIEX_KD_PROT": [0.0, 0.0, 0.9, 4.8, 1.1],
	        "GIEX_NU": [1.2, 0.0, 2.0, 3.7, 1.5],
	        "GIEX_NU_LIN": [0.0, 0.0, 0.5, 1.3, 1.1],
	        "GIEX_NU_QUAD": [0.0, 0.0, 0.8, 4.7, 0.9],
	        "GIEX_SIGMA": [0.0, 0.0, 11.83, 10.0, 10.6],
	        "GIEX_LAMBDA": 100.0
	)json", \
	R"json( "EXT_GIEX_KA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_T": [0.0, 0.0, 3.55, 1.59],
	        "EXT_GIEX_KA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_LIN": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_LIN_T": [0.0, 0.0, 0.8, 0.9],
	        "EXT_GIEX_KA_LIN_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_LIN_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_QUAD": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_QUAD_T": [0.0, 0.0, 0.5, 0.6],
	        "EXT_GIEX_KA_QUAD_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_QUAD_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_SALT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_SALT_T": [0.0, 0.0, 1.2, 0.8],
	        "EXT_GIEX_KA_SALT_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_SALT_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_PROT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_PROT_T": [0.0, 0.0, 1.1, 0.9],
	        "EXT_GIEX_KA_PROT_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_PROT_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_T": [0.0, 0.0, 10.0, 10.0],
	        "EXT_GIEX_KD_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_LIN": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_LIN_T": [0.0, 0.0, 0.6, 0.5],
	        "EXT_GIEX_KD_LIN_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_LIN_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_QUAD": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_QUAD_T": [0.0, 0.0, 0.9, 0.7],
	        "EXT_GIEX_KD_QUAD_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_QUAD_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_SALT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_SALT_T": [0.0, 0.0, 1.0, 1.2],
	        "EXT_GIEX_KD_SALT_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_SALT_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_PROT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_PROT_T": [0.0, 0.0, 0.9, 1.1],
	        "EXT_GIEX_KD_PROT_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_PROT_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_NU": [1.2, 0.0, 0.0, 0.0],
	        "EXT_GIEX_NU_T": [0.0, 0.0, 2.0, 1.5],
	        "EXT_GIEX_NU_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_NU_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_NU_LIN": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_NU_LIN_T": [0.0, 0.0, 0.5, 1.1],
	        "EXT_GIEX_NU_LIN_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_NU_LIN_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_NU_QUAD": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_NU_QUAD_T": [0.0, 0.0, 0.8, 0.9],
	        "EXT_GIEX_NU_QUAD_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_NU_QUAD_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_SIGMA": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_SIGMA_T": [0.0, 0.0, 11.83, 10.6],
	        "EXT_GIEX_SIGMA_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_SIGMA_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_LAMBDA": 0.0,
	        "EXT_GIEX_LAMBDA_T": 100.0,
	        "EXT_GIEX_LAMBDA_TT": 0.0,
	        "EXT_GIEX_LAMBDA_TTT": 0.0
	)json", \
	R"json( "EXT_GIEX_KA": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_T": [0.0, 0.0, 3.55, 7.7, 1.59],
	        "EXT_GIEX_KA_TT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_TTT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_LIN": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_LIN_T": [0.0, 0.0, 0.8, 5.5, 0.9],
	        "EXT_GIEX_KA_LIN_TT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_LIN_TTT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_QUAD": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_QUAD_T": [0.0, 0.0, 0.5, 4.2, 0.6],
	        "EXT_GIEX_KA_QUAD_TT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_QUAD_TTT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_SALT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_SALT_T": [0.0, 0.0, 1.2, 6.3, 0.8],
	        "EXT_GIEX_KA_SALT_TT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_SALT_TTT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_PROT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_PROT_T": [0.0, 0.0, 1.1, 3.6, 0.9],
	        "EXT_GIEX_KA_PROT_TT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KA_PROT_TTT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_T": [0.0, 0.0, 10.0, 10.0, 10.0],
	        "EXT_GIEX_KD_TT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_TTT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_LIN": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_LIN_T": [0.0, 0.0, 0.6, 2.4, 0.5],
	        "EXT_GIEX_KD_LIN_TT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_LIN_TTT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_QUAD": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_QUAD_T": [0.0, 0.0, 0.9, 2.1, 0.7],
	        "EXT_GIEX_KD_QUAD_TT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_QUAD_TTT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_SALT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_SALT_T": [0.0, 0.0, 1.0, 3.4, 1.2],
	        "EXT_GIEX_KD_SALT_TT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_SALT_TTT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_PROT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_PROT_T": [0.0, 0.0, 0.9, 8.0, 1.1],
	        "EXT_GIEX_KD_PROT_TT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_KD_PROT_TTT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_NU": [1.2, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_NU_T": [0.0, 0.0, 2.0, 3.7, 1.5],
	        "EXT_GIEX_NU_TT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_NU_TTT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_NU_LIN": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_NU_LIN_T": [0.0, 0.0, 0.5, 3.0, 1.1],
	        "EXT_GIEX_NU_LIN_TT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_NU_LIN_TTT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_NU_QUAD": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_NU_QUAD_T": [0.0, 0.0, 0.8, 2.2, 0.9],
	        "EXT_GIEX_NU_QUAD_TT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_NU_QUAD_TTT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_SIGMA": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_SIGMA_T": [0.0, 0.0, 11.83, 10.0, 10.6],
	        "EXT_GIEX_SIGMA_TT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_SIGMA_TTT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_GIEX_LAMBDA": 0.0,
	        "EXT_GIEX_LAMBDA_T": 100.0,
	        "EXT_GIEX_LAMBDA_TT": 0.0,
	        "EXT_GIEX_LAMBDA_TTT": 0.0
	)json", \
	1e-10, 1e-8, CADET_NONBINDING_LIQUIDPHASE_COMP_USED, CADET_COMPARE_BINDING_VS_NONBINDING);
