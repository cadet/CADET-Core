// =============================================================================
//  CADET
//  
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
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
	1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING)

	CADET_BINDINGTEST("FREUNDLICH_LDF", "EXT_FREUNDLICH_LDF", (1, 1, 1, 1), (1, 0, 1, 1, 1), (1.0, 0.0, 2.0, -1e-4, 0.0, 0.0, 0.0, 0.0), (1.0, 3.0, 0.0, 2.0, -1e-4, 0.0, 0.0, 0.0, 0.0), \
		R"json( "FLDF_KKIN": [1.0, 1.0, 1.2, 2.0],
	        "FLDF_KF": [0.1, 0.3, 0.3, 0.2],
			"FLDF_N": [0.5, 1.0, 1.2, 0.8]
	)json", \
		R"json( "FLDF_KKIN": [1.0, 0.5, 1.0, 1.2, 2.0],
	        "FLDF_KF": [0.1, 0.3, 0.3, 0.3, 0.2],
			"FLDF_N": [0.5, 2.2, 1.0, 1.2, 0.8]
	)json", \
		R"json( "EXT_FLDF_KKIN": [0.0, 0.0, 0.0, 0.0],
	        "EXT_FLDF_KKIN_T": [1.0, 1.0, 1.2, 2.0],
	        "EXT_FLDF_KKIN_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_FLDF_KKIN_TTT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_FLDF_KF": [0.0, 0.0, 0.0, 0.0],
	        "EXT_FLDF_KF_T": [0.1, 0.3, 0.3, 0.2],
	        "EXT_FLDF_KF_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_FLDF_KF_TTT": [0.0, 0.0, 0.0, 0.0],
			"EXT_FLDF_N": [0.0, 0.0, 0.0, 0.0],
	        "EXT_FLDF_N_T": [0.5, 1.0, 1.2, 0.8],
	        "EXT_FLDF_N_TT": [0.0, 0.0, 0.0, 0.0],
	        "EXT_FLDF_N_TTT": [0.0, 0.0, 0.0, 0.0]
	)json", \
		R"json( "EXT_FLDF_KKIN": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_FLDF_KKIN_T": [1.0, 0.5, 1.0, 1.2, 2.0],
	        "EXT_FLDF_KKIN_TT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_FLDF_KKIN_TTT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_FLDF_KF": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_FLDF_KF_T": [0.1, 0.3, 0.3, 0.3, 0.2],
	        "EXT_FLDF_KF_TT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_FLDF_KF_TTT": [0.0, 0.0, 0.0, 0.0, 0.0],
			"EXT_FLDF_N": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_FLDF_N_T": [0.5, 2.0, 1.0, 1.2, 0.8],
	        "EXT_FLDF_N_TT": [0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_FLDF_N_TTT": [0.0, 0.0, 0.0, 0.0, 0.0]
	)json", \
	1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING)

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
	1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING)

  CADET_BINDINGTEST("MULTI_COMPONENT_LANGMUIR_LDF", "EXT_MULTI_COMPONENT_LANGMUIR_LDF", (1, 1), (1, 0, 1), (1.0, 2.0, 0.0, 0.0), (1.0, 3.0, 2.0, 0.0, 0.0), \
		R"json( "MCLLDF_KEQ": [1.14, 2.0],
	        "MCLLDF_KKIN": [0.004, 0.008],
	        "MCLLDF_QMAX": [4.88, 3.5]
	)json", \
		R"json( "MCLLDF_KEQ": [1.14, 1.0, 2.0],
	        "MCLLDF_KKIN": [0.004, 2.0, 0.008],
	        "MCLLDF_QMAX": [4.88, 3.0, 3.5]
	)json", \
		R"json( "EXT_MCLLDF_KEQ": [0.0, 0.0],
	        "EXT_MCLLDF_KEQ_T": [1.14, 2.0],
	        "EXT_MCLLDF_KEQ_TT": [0.0, 0.0],
	        "EXT_MCLLDF_KEQ_TTT": [0.0, 0.0],
	        "EXT_MCLLDF_KKIN": [0.0, 0.0],
	        "EXT_MCLLDF_KKIN_T": [0.004, 0.008],
	        "EXT_MCLLDF_KKIN_TT": [0.0, 0.0],
	        "EXT_MCLLDF_KKIN_TTT": [0.0, 0.0],
	        "EXT_MCLLDF_QMAX": [0.0, 0.0],
	        "EXT_MCLLDF_QMAX_T": [4.88, 3.5],
	        "EXT_MCLLDF_QMAX_TT": [0.0, 0.0],
	        "EXT_MCLLDF_QMAX_TTT": [0.0, 0.0]
	)json", \
		R"json( "EXT_MCLLDF_KEQ": [0.0, 0.0, 0.0],
	        "EXT_MCLLDF_KEQ_T": [1.14, 1.0, 2.0],
	        "EXT_MCLLDF_KEQ_TT": [0.0, 0.0, 0.0],
	        "EXT_MCLLDF_KEQ_TTT": [0.0, 0.0, 0.0],
	        "EXT_MCLLDF_KKIN": [0.0, 0.0, 0.0],
	        "EXT_MCLLDF_KKIN_T": [0.004, 2.0, 0.008],
	        "EXT_MCLLDF_KKIN_TT": [0.0, 0.0, 0.0],
	        "EXT_MCLLDF_KKIN_TTT": [0.0, 0.0, 0.0],
	        "EXT_MCLLDF_QMAX": [0.0, 0.0, 0.0],
	        "EXT_MCLLDF_QMAX_T": [4.88, 3.0, 3.5],
	        "EXT_MCLLDF_QMAX_TT": [0.0, 0.0, 0.0],
	        "EXT_MCLLDF_QMAX_TTT": [0.0, 0.0, 0.0]
	)json", \
		1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING)

  CADET_BINDINGTEST("MULTI_COMPONENT_LANGMUIR_LDF_LIQUID_PHASE", "EXT_MULTI_COMPONENT_LANGMUIR_LDF_LIQUID_PHASE", (1, 1), (1, 0, 1), (1.0, 2.0, 0.0, 0.0), (1.0, 3.0, 2.0, 0.0, 0.0), \
		R"json( "MCLLDFC_KEQ": [1.14, 2.0],
	        "MCLLDFC_KKIN": [0.004, 0.008],
	        "MCLLDFC_QMAX": [4.88, 3.5]
	)json", \
		R"json( "MCLLDFC_KEQ": [1.14, 1.0, 2.0],
	        "MCLLDFC_KKIN": [0.004, 2.0, 0.008],
	        "MCLLDFC_QMAX": [4.88, 3.0, 3.5]
	)json", \
		R"json( "EXT_MCLLDFC_KEQ": [0.0, 0.0],
	        "EXT_MCLLDFC_KEQ_T": [1.14, 2.0],
	        "EXT_MCLLDFC_KEQ_TT": [0.0, 0.0],
	        "EXT_MCLLDFC_KEQ_TTT": [0.0, 0.0],
	        "EXT_MCLLDFC_KKIN": [0.0, 0.0],
	        "EXT_MCLLDFC_KKIN_T": [0.004, 0.008],
	        "EXT_MCLLDFC_KKIN_TT": [0.0, 0.0],
	        "EXT_MCLLDFC_KKIN_TTT": [0.0, 0.0],
	        "EXT_MCLLDFC_QMAX": [0.0, 0.0],
	        "EXT_MCLLDFC_QMAX_T": [4.88, 3.5],
	        "EXT_MCLLDFC_QMAX_TT": [0.0, 0.0],
	        "EXT_MCLLDFC_QMAX_TTT": [0.0, 0.0]
	)json", \
		R"json( "EXT_MCLLDFC_KEQ": [0.0, 0.0, 0.0],
	        "EXT_MCLLDFC_KEQ_T": [1.14, 1.0, 2.0],
	        "EXT_MCLLDFC_KEQ_TT": [0.0, 0.0, 0.0],
	        "EXT_MCLLDFC_KEQ_TTT": [0.0, 0.0, 0.0],
	        "EXT_MCLLDFC_KKIN": [0.0, 0.0, 0.0],
	        "EXT_MCLLDFC_KKIN_T": [0.004, 2.0, 0.008],
	        "EXT_MCLLDFC_KKIN_TT": [0.0, 0.0, 0.0],
	        "EXT_MCLLDFC_KKIN_TTT": [0.0, 0.0, 0.0],
	        "EXT_MCLLDFC_QMAX": [0.0, 0.0, 0.0],
	        "EXT_MCLLDFC_QMAX_T": [4.88, 3.0, 3.5],
	        "EXT_MCLLDFC_QMAX_TT": [0.0, 0.0, 0.0],
	        "EXT_MCLLDFC_QMAX_TTT": [0.0, 0.0, 0.0]
	)json", \
		1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING)

CADET_BINDINGTEST_MULTI("AFFINITY_COMPLEX_TITRATION", "EXT_AFFINITY_COMPLEX_TITRATION", " with log10 values (pH)", (0, 1), (0, 0, 1), (5.5, 0.06666667, 0.0), (5.5, 0.1, 0.06666667, 0.0), \
	R"json( "ACT_KA": [1.0, 1386.9],
	        "ACT_KD": [1.0, 1.0],
	        "ACT_QMAX": [0.0000000001, 10.7439614],
	        "ACT_ETAA": [0.0, 1.81],
	        "ACT_ETAG": [0.0, 2.28],
	        "ACT_PKAA": [0.0, 2.07],
            "ACT_PKAG": [0.0, 5.29],
            "ACT_USE_ION_CONC": 0 
	)json", \
	R"json( "ACT_KA": [1.0, 1.0, 1386.9],
	        "ACT_KD": [1.0, 1.0, 1.0],
	        "ACT_QMAX": [0.0000000001, 0.0000000001, 10.7439614],
	        "ACT_ETAA": [0.0, 0.0, 1.81],
	        "ACT_ETAG": [0.0, 0.0, 2.28],
	        "ACT_PKAA": [0.0, 0.0, 2.07],
            "ACT_PKAG": [0.0, 0.0, 5.29],
            "ACT_USE_ION_CONC": 0 
	)json", \
	R"json( "EXT_ACT_KA": [0.0, 0.0],
	        "EXT_ACT_KA_T": [1.0, 1386.9],
	        "EXT_ACT_KA_TT": [0.0, 0.0],
	        "EXT_ACT_KA_TTT": [0.0, 0.0],
	        "EXT_ACT_KD": [0.0, 0.0],
	        "EXT_ACT_KD_T": [1.0, 1.0],
	        "EXT_ACT_KD_TT": [0.0, 0.0],
	        "EXT_ACT_KD_TTT": [0.0, 0.0],
	        "EXT_ACT_QMAX": [0.0, 0.0],
	        "EXT_ACT_QMAX_T": [0.0000000001, 10.7439614],
	        "EXT_ACT_QMAX_TT": [0.0, 0.0],
	        "EXT_ACT_QMAX_TTT": [0.0, 0.0],
	        "EXT_ACT_ETAA": [0.0, 0.0],
	        "EXT_ACT_ETAA_T": [0.0, 1.81],
	        "EXT_ACT_ETAA_TT": [0.0, 0.0],
	        "EXT_ACT_ETAA_TTT": [0.0, 0.0],
	        "EXT_ACT_ETAG": [0.0, 0.0],
	        "EXT_ACT_ETAG_T": [0.0, 2.28],
	        "EXT_ACT_ETAG_TT": [0.0, 0.0],
	        "EXT_ACT_ETAG_TTT": [0.0, 0.0],
	        "EXT_ACT_PKAA": [0.0, 0.0],
	        "EXT_ACT_PKAA_T": [0.0, 2.07],
	        "EXT_ACT_PKAA_TT": [0.0, 0.0],
	        "EXT_ACT_PKAA_TTT": [0.0, 0.0],
	        "EXT_ACT_PKAG": [0.0, 0.0],
	        "EXT_ACT_PKAG_T": [0.0, 5.29],
	        "EXT_ACT_PKAG_TT": [0.0, 0.0],
	        "EXT_ACT_PKAG_TTT": [0.0, 0.0],
            "EXT_ACT_USE_ION_CONC": 0 
	)json", \
	R"json( "EXT_ACT_KA": [0.0, 0.0, 0.0],
	        "EXT_ACT_KA_T": [1.0, 1.0, 1386.9],
	        "EXT_ACT_KA_TT": [0.0, 0.0, 0.0],
	        "EXT_ACT_KA_TTT": [0.0, 0.0, 0.0],
	        "EXT_ACT_KD": [0.0, 0.0, 0.0],
	        "EXT_ACT_KD_T": [1.0, 1.0, 1.0],
	        "EXT_ACT_KD_TT": [0.0, 0.0, 0.0],
	        "EXT_ACT_KD_TTT": [0.0, 0.0, 0.0],
	        "EXT_ACT_QMAX": [0.0, 0.0, 0.0],
	        "EXT_ACT_QMAX_T": [0.0000000001, 0.0000000001, 10.7439614],
	        "EXT_ACT_QMAX_TT": [0.0, 0.0, 0.0],
	        "EXT_ACT_QMAX_TTT": [0.0, 0.0, 0.0],
	        "EXT_ACT_ETAA": [0.0, 0.0, 0.0],
	        "EXT_ACT_ETAA_T": [0.0, 0.0, 1.81],
	        "EXT_ACT_ETAA_TT": [0.0, 0.0, 0.0],
	        "EXT_ACT_ETAA_TTT": [0.0, 0.0, 0.0],
	        "EXT_ACT_ETAG": [0.0, 0.0, 0.0],
	        "EXT_ACT_ETAG_T": [0.0, 0.0, 2.28],
	        "EXT_ACT_ETAG_TT": [0.0, 0.0, 0.0],
	        "EXT_ACT_ETAG_TTT": [0.0, 0.0, 0.0],
	        "EXT_ACT_PKAA": [0.0, 0.0, 0.0],
	        "EXT_ACT_PKAA_T": [0.0, 0.0, 2.07],
	        "EXT_ACT_PKAA_TT": [0.0, 0.0, 0.0],
	        "EXT_ACT_PKAA_TTT": [0.0, 0.0, 0.0],
	        "EXT_ACT_PKAG": [0.0, 0.0, 0.0],
	        "EXT_ACT_PKAG_T": [0.0, 0.0, 5.29],
	        "EXT_ACT_PKAG_TT": [0.0, 0.0, 0.0],
	        "EXT_ACT_PKAG_TTT": [0.0, 0.0, 0.0],
            "EXT_ACT_USE_ION_CONC": 0 
	)json", \
	1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_USED, CADET_COMPARE_BINDING_VS_NONBINDING)

	CADET_BINDINGTEST_MULTI("AFFINITY_COMPLEX_TITRATION", "EXT_AFFINITY_COMPLEX_TITRATION", " with raw values (ions in mol)", (0, 1), (0, 0, 1), (0.4, 0.06666667, 0.0), (0.4, 0.1, 0.06666667, 0.0), \
	R"json( "ACT_KA": [1.0, 1386.9],
	        "ACT_KD": [1.0, 1.0],
	        "ACT_QMAX": [0.0000000001, 10.7439614],
	        "ACT_ETAA": [0.0, 1.81],
	        "ACT_ETAG": [0.0, 2.28],
	        "ACT_CMID_A": [0.0, 0.380],
            "ACT_CMID_G": [0.0, 0.2],
            "ACT_USE_ION_CONC": 1
	)json", \
	R"json( "ACT_KA": [1.0, 1.0, 1386.9],
	        "ACT_KD": [1.0, 1.0, 1.0],
	        "ACT_QMAX": [0.0000000001, 0.0000000001, 10.7439614],
	        "ACT_ETAA": [0.0, 0.0, 1.81],
	        "ACT_ETAG": [0.0, 0.0, 2.28],
	        "ACT_CMID_A": [0.0, 0.0, 0.38],
            "ACT_CMID_G": [0.0, 0.0, 0.2],
            "ACT_USE_ION_CONC": 1
	)json", \
	R"json( "EXT_ACT_KA": [0.0, 0.0],
	        "EXT_ACT_KA_T": [1.0, 1386.9],
	        "EXT_ACT_KA_TT": [0.0, 0.0],
	        "EXT_ACT_KA_TTT": [0.0, 0.0],
	        "EXT_ACT_KD": [0.0, 0.0],
	        "EXT_ACT_KD_T": [1.0, 1.0],
	        "EXT_ACT_KD_TT": [0.0, 0.0],
	        "EXT_ACT_KD_TTT": [0.0, 0.0],
	        "EXT_ACT_QMAX": [0.0, 0.0],
	        "EXT_ACT_QMAX_T": [0.0000000001, 10.7439614],
	        "EXT_ACT_QMAX_TT": [0.0, 0.0],
	        "EXT_ACT_QMAX_TTT": [0.0, 0.0],
	        "EXT_ACT_ETAA": [0.0, 0.0],
	        "EXT_ACT_ETAA_T": [0.0, 1.81],
	        "EXT_ACT_ETAA_TT": [0.0, 0.0],
	        "EXT_ACT_ETAA_TTT": [0.0, 0.0],
	        "EXT_ACT_ETAG": [0.0, 0.0],
	        "EXT_ACT_ETAG_T": [0.0, 2.28],
	        "EXT_ACT_ETAG_TT": [0.0, 0.0],
	        "EXT_ACT_ETAG_TTT": [0.0, 0.0],
	        "EXT_ACT_CMID_A": [0.0, 0.0],
	        "EXT_ACT_CMID_A_T": [0.0, 0.38],
	        "EXT_ACT_CMID_A_TT": [0.0, 0.0],
	        "EXT_ACT_CMID_A_TTT": [0.0, 0.0],
	        "EXT_ACT_CMID_G": [0.0, 0.0],
	        "EXT_ACT_CMID_G_T": [0.0, 0.2],
	        "EXT_ACT_CMID_G_TT": [0.0, 0.0],
	        "EXT_ACT_CMID_G_TTT": [0.0, 0.0],
            "EXT_ACT_USE_ION_CONC": 1
	)json", \
	R"json( "EXT_ACT_KA": [0.0, 0.0, 0.0],
	        "EXT_ACT_KA_T": [1.0, 1.0, 1386.9],
	        "EXT_ACT_KA_TT": [0.0, 0.0, 0.0],
	        "EXT_ACT_KA_TTT": [0.0, 0.0, 0.0],
	        "EXT_ACT_KD": [0.0, 0.0, 0.0],
	        "EXT_ACT_KD_T": [1.0, 1.0, 1.0],
	        "EXT_ACT_KD_TT": [0.0, 0.0, 0.0],
	        "EXT_ACT_KD_TTT": [0.0, 0.0, 0.0],
	        "EXT_ACT_QMAX": [0.0, 0.0, 0.0],
	        "EXT_ACT_QMAX_T": [0.0000000001, 0.0000000001, 10.7439614],
	        "EXT_ACT_QMAX_TT": [0.0, 0.0, 0.0],
	        "EXT_ACT_QMAX_TTT": [0.0, 0.0, 0.0],
	        "EXT_ACT_ETAA": [0.0, 0.0, 0.0],
	        "EXT_ACT_ETAA_T": [0.0, 0.0, 1.81],
	        "EXT_ACT_ETAA_TT": [0.0, 0.0, 0.0],
	        "EXT_ACT_ETAA_TTT": [0.0, 0.0, 0.0],
	        "EXT_ACT_ETAG": [0.0, 0.0, 0.0],
	        "EXT_ACT_ETAG_T": [0.0, 0.0, 2.28],
	        "EXT_ACT_ETAG_TT": [0.0, 0.0, 0.0],
	        "EXT_ACT_ETAG_TTT": [0.0, 0.0, 0.0],
	        "EXT_ACT_CMID_A": [0.0, 0.0, 0.0],
	        "EXT_ACT_CMID_A_T": [0.0, 0.0, 0.38],
	        "EXT_ACT_CMID_A_TT": [0.0, 0.0, 0.0],
	        "EXT_ACT_CMID_A_TTT": [0.0, 0.0, 0.0],
	        "EXT_ACT_CMID_G": [0.0, 0.0, 0.0],
	        "EXT_ACT_CMID_G_T": [0.0, 0.0, 0.2],
	        "EXT_ACT_CMID_G_TT": [0.0, 0.0, 0.0],
	        "EXT_ACT_CMID_G_TTT": [0.0, 0.0, 0.0],
            "EXT_ACT_USE_ION_CONC": 1
	)json", \
	1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_USED, CADET_COMPARE_BINDING_VS_NONBINDING)

	// Note: the allBinding has six components to test both interpolation and extrapolation of the data, i.e. different initial values that are within or beyond the range of the provided data.
	CADET_BINDINGTEST("SPLINE_INTERPOLATION", "EXT_SPLINE_INTERPOLATION", (1, 1, 1, 1, 1, 1), (1, 0, 1, 1, 1, 1, 1), (1.5, 0.0, 0.1, 1.1, 0.3, 1.8, 0.2, 0.3, 0.8, 0.33, 0.1, 0.75), (1.5, 1.0, 0.0, 0.1, 1.1, 0.3, 1.8, 0.2, 0.3, 0.8, 0.33, 0.1, 0.75), \
	R"json( "CP_VALS_COMP_000": [
                0.0, 0.19219219, 0.37837838, 0.75075075, 1.5015015, 3.003003, 6.0
            ],
            "CP_VALS_COMP_001": [
                0.0, 0.38438438, 0.75675676, 1.5015015, 3.003003, 6.006006, 12.0
            ],
            "CP_VALS_COMP_002": [
                0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                2.5, 2.6, 2.7, 2.8, 2.9, 3.0,
                3.5, 3.6, 3.7, 3.8, 3.9, 4.0,
                4.5, 4.6, 4.7, 4.8, 4.9, 5.0
            ],
            "CP_VALS_COMP_003": [ 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ],
            "CP_VALS_COMP_004": [ 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ],
            "CP_VALS_COMP_005": [ 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ],
	        "CS_VALS_COMP_000_BND_000": [
				1.29198742, 76.84778516, 110.14325476, 141.59781958, 165.80146634, 181.33282173, 190.23846417
			],
	        "CS_VALS_COMP_001_BND_000": [
				1.93798113, 115.27167774, 165.21488214, 212.39672937, 248.70219951, 271.99923260, 285.35769626
			],
            "CS_VALS_COMP_002_BND_000": [
                0.5, 0.6, 1.7, 0.8, 2.9, 1.0,
                1.5, 1.6, 4.7, 4.8, 1.9, 2.0,
                2.5, 3.6, 2.7, 3.8, 2.9, 3.0,
                3.5, 2.6, 3.7, 2.8, 3.9, 4.0,
                4.5, 1.6, 4.7, 1.8, 4.9, 5.0
            ],
            "CS_VALS_COMP_003_BND_000": [ 0.1, 0.3, 0.35, 0.4, 0.5, 0.5 ],
            "CS_VALS_COMP_004_BND_000": [ 0.2, 0.3, 0.35, 0.4, 0.41, 0.5 ],
            "CS_VALS_COMP_005_BND_000": [ 0.3, 0.31, 0.35, 0.4, 0.45, 0.5 ],
			"SPLINE_KKIN": [2.0, 1.0, 1.2, 1.1, 3.2, 4.0]
	)json", \
	R"json( "CP_VALS_COMP_000": [
				0.0, 0.19219219, 0.37837838, 0.75075075, 1.5015015, 3.003003, 6.0
            ],
            "CP_VALS_COMP_002": [
				0.0, 0.38438438, 0.75675676, 1.5015015, 3.003003, 6.006006, 12.0
            ],
            "CP_VALS_COMP_003": [
                0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                2.5, 2.6, 2.7, 2.8, 2.9, 3.0,
                3.5, 3.6, 3.7, 3.8, 3.9, 4.0,
                4.5, 4.6, 4.7, 4.8, 4.9, 5.0
            ],
            "CP_VALS_COMP_004": [ 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ],
            "CP_VALS_COMP_005": [ 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ],
            "CP_VALS_COMP_006": [ 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ],
	        "CS_VALS_COMP_000_BND_000": [
				1.29198742, 76.84778516, 110.14325476, 141.59781958, 165.80146634, 181.33282173, 190.23846417
			],
	        "CS_VALS_COMP_002_BND_000": [
				1.93798113, 115.27167774, 165.21488214, 212.39672937, 248.70219951, 271.99923260, 285.35769626
			],
            "CS_VALS_COMP_003_BND_000": [
                0.5, 0.6, 1.7, 0.8, 2.9, 1.0,
                1.5, 1.6, 4.7, 4.8, 1.9, 2.0,
                2.5, 3.6, 2.7, 3.8, 2.9, 3.0,
                3.5, 2.6, 3.7, 2.8, 3.9, 4.0,
                4.5, 1.6, 4.7, 1.8, 4.9, 5.0
            ],
            "CS_VALS_COMP_004_BND_000": [ 0.1, 0.3, 0.35, 0.4, 0.5, 0.5 ],
            "CS_VALS_COMP_005_BND_000": [ 0.2, 0.3, 0.35, 0.4, 0.41, 0.5 ],
            "CS_VALS_COMP_006_BND_000": [ 0.3, 0.31, 0.35, 0.4, 0.45, 0.5 ],
			"SPLINE_KKIN": [2.0, 1.0, 1.2, 1.1, 3.2, 4.0]
	)json", \
	R"json( "CP_VALS_COMP_000": [
                0.0, 0.19219219, 0.37837838, 0.75075075, 1.5015015, 3.003003, 6.0
            ],
            "CP_VALS_COMP_001": [
                0.0, 0.38438438, 0.75675676, 1.5015015, 3.003003, 6.006006, 12.0
            ],
            "CP_VALS_COMP_002": [
                0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                2.5, 2.6, 2.7, 2.8, 2.9, 3.0,
                3.5, 3.6, 3.7, 3.8, 3.9, 4.0,
                4.5, 4.6, 4.7, 4.8, 4.9, 5.0
            ],
            "CP_VALS_COMP_003": [ 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ],
            "CP_VALS_COMP_004": [ 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ],
            "CP_VALS_COMP_005": [ 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ],
	        "CS_VALS_COMP_000_BND_000": [
				1.29198742, 76.84778516, 110.14325476, 141.59781958, 165.80146634, 181.33282173, 190.23846417
			],
	        "CS_VALS_COMP_001_BND_000": [
				1.93798113, 115.27167774, 165.21488214, 212.39672937, 248.70219951, 271.99923260, 285.35769626
			],
            "CS_VALS_COMP_002_BND_000": [
                0.5, 0.6, 1.7, 0.8, 2.9, 1.0,
                1.5, 1.6, 4.7, 4.8, 1.9, 2.0,
                2.5, 3.6, 2.7, 3.8, 2.9, 3.0,
                3.5, 2.6, 3.7, 2.8, 3.9, 4.0,
                4.5, 1.6, 4.7, 1.8, 4.9, 5.0
            ],
            "CS_VALS_COMP_003_BND_000": [ 0.1, 0.3, 0.35, 0.4, 0.5, 0.5 ],
            "CS_VALS_COMP_004_BND_000": [ 0.2, 0.3, 0.35, 0.4, 0.41, 0.5 ],
            "CS_VALS_COMP_005_BND_000": [ 0.3, 0.31, 0.35, 0.4, 0.45, 0.5 ],
			"EXT_SPLINE_KKIN": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
			"EXT_SPLINE_KKIN_T": [2.0, 1.0, 1.2, 1.1, 3.2, 4.0],
	        "EXT_SPLINE_KKIN_TT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_SPLINE_KKIN_TTT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	)json", \
	R"json( "CP_VALS_COMP_000": [
                0.0, 0.19219219, 0.37837838, 0.75075075, 1.5015015, 3.003003, 6.0
            ],
            "CP_VALS_COMP_002": [
                0.0, 0.38438438, 0.75675676, 1.5015015, 3.003003, 6.006006, 12.0
            ],
            "CP_VALS_COMP_003": [
                0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                2.5, 2.6, 2.7, 2.8, 2.9, 3.0,
                3.5, 3.6, 3.7, 3.8, 3.9, 4.0,
                4.5, 4.6, 4.7, 4.8, 4.9, 5.0
            ],
            "CP_VALS_COMP_004": [ 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ],
            "CP_VALS_COMP_005": [ 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ],
            "CP_VALS_COMP_006": [ 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ],
	        "CS_VALS_COMP_000_BND_000": [
				1.29198742, 76.84778516, 110.14325476, 141.59781958, 165.80146634, 181.33282173, 190.23846417
			],
	        "CS_VALS_COMP_002_BND_000": [
				1.93798113, 115.27167774, 165.21488214, 212.39672937, 248.70219951, 271.99923260, 285.35769626
			],
            "CS_VALS_COMP_003_BND_000": [
                0.5, 0.6, 1.7, 0.8, 2.9, 1.0,
                1.5, 1.6, 4.7, 4.8, 1.9, 2.0,
                2.5, 3.6, 2.7, 3.8, 2.9, 3.0,
                3.5, 2.6, 3.7, 2.8, 3.9, 4.0,
                4.5, 1.6, 4.7, 1.8, 4.9, 5.0
            ],
            "CS_VALS_COMP_004_BND_000": [ 0.1, 0.3, 0.35, 0.4, 0.5, 0.5 ],
            "CS_VALS_COMP_005_BND_000": [ 0.2, 0.3, 0.35, 0.4, 0.41, 0.5 ],
            "CS_VALS_COMP_006_BND_000": [ 0.3, 0.31, 0.35, 0.4, 0.45, 0.5 ],
			"EXT_SPLINE_KKIN": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
			"EXT_SPLINE_KKIN_T": [2.0, 1.0, 1.2, 1.1, 3.2, 4.0],
	        "EXT_SPLINE_KKIN_TT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
	        "EXT_SPLINE_KKIN_TTT": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
	)json", \
	1e-10, 1e-10, CADET_NONBINDING_LIQUIDPHASE_COMP_UNUSED, CADET_COMPARE_BINDING_VS_NONBINDING)
