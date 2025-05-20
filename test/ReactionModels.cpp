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

#include "ReactionModelTests.hpp"
#include "ColumnTests.hpp"

TEST_CASE("MassActionLaw kinetic analytic Jacobian vs AD", "[MassActionLaw],[ReactionModel],[Jacobian],[AD]")
{
	const unsigned int nBound[] = {1, 2, 1};
	const double point[] = {1.0, 2.0, 1.4, 2.1, 0.2, 1.1, 1.8};
	cadet::test::reaction::testDynamicJacobianAD("MASS_ACTION_LAW", 3, nBound,
		R"json({
			"MAL_KFWD_BULK": [1.0, 2.0, 0.4],
			"MAL_KBWD_BULK": [0.0, 0.2, 1.5],
			"MAL_STOICHIOMETRY_BULK": [ 1.0, -2.0,  3.0,
			                           -1.0,  0.0, -2.0,
			                            0.0,  1.0,  1.0],
			"MAL_EXPONENTS_BULK_FWD": [ 1.2,  0.0,  0.0,
			                            0.0,  1.3,  2.2,
			                            0.0,  1.0,  1.1],
			"MAL_EXPONENTS_BULK_BWD": [ 0.8,  2.1,  1.0,
			                            1.3,  1.0,  0.0,
			                            0.0,  0.0,  1.4],

			"MAL_KFWD_LIQUID": [1.0, 2.0, 0.8, 1.2, 2.1],
			"MAL_KBWD_LIQUID": [0.1, 0.2, 2.4, 1.9, 0.8],
			"MAL_STOICHIOMETRY_LIQUID": [ 1.0, -2.0,  3.0,  1.0, -2.0,
			                             -1.0,  0.0, -2.0, -3.0, -1.0, 
			                              0.0,  1.0,  1.0,  2.0,  3.0],
			"MAL_EXPONENTS_LIQUID_FWD": [ 0.4,  0.0,  2.0,  1.2,  0.0,
			                              1.0,  1.0,  2.0,  0.0,  1.0, 
			                              0.0,  2.0,  0.0,  2.4,  1.4],
			"MAL_EXPONENTS_LIQUID_BWD": [ 1.0,  2.0,  0.0,  0.2,  2.0,
			                              0.0,  1.0,  2.0,  3.0,  1.0, 
			                              1.0,  0.0,  0.0,  1.0,  0.0],
			"MAL_EXPONENTS_LIQUID_FWD_MODSOLID": [ 1.0,  0.0,  0.0,  1.2,  0.0,
			                                       1.0,  1.6,  2.1,  0.0,  1.0, 
			                                       0.0,  0.0,  0.0,  0.0,  1.0, 
			                                       0.0,  0.0,  0.0,  2.4,  0.0],
			"MAL_EXPONENTS_LIQUID_BWD_MODSOLID": [ 0.0,  2.0,  0.0,  0.0,  0.0,
			                                       0.0,  0.0,  0.0,  1.0,  0.0, 
			                                       0.0,  1.0,  2.0,  1.2,  0.0, 
			                                       1.0,  0.0,  0.0,  0.0,  1.4],

			"MAL_KFWD_SOLID": [1.0, 2.0],
			"MAL_KBWD_SOLID": [0.8, 1.2],
			"MAL_STOICHIOMETRY_SOLID": [ 1.0, -2.0,
			                            -1.0,  0.0,
			                            -1.0,  2.0,
			                             0.0,  1.0],
			"MAL_EXPONENTS_SOLID_FWD": [ 1.0,  0.0,
			                             0.0,  0.2,
			                             1.2,  1.8,
			                             0.0,  1.0],
			"MAL_EXPONENTS_SOLID_BWD": [ 0.0,  2.0,
			                             2.0,  0.0,
			                             0.0,  0.0,
			                             1.0,  0.8],
			"MAL_EXPONENTS_SOLID_FWD_MODLIQUID": [ 1.0, 0.0,
			                                       0.0, 2.0,
			                                       0.0, 1.5],
			"MAL_EXPONENTS_SOLID_BWD_MODLIQUID": [ 1.2, 1.6,
			                                       2.0, 0.0,
			                                       1.8, 0.0]
		})json",
		point, 1e-15, 1e-15
	);
}

TEST_CASE("MichaelisMenten kinetic and specific mass action law micro-kinetics yield same result", "[MichaelisMenten],[ReactionModel],[Simulation],[CI]")
{
	const std::string& configFilePath1 = std::string("/data/model_CSTR_MichaelisMenten_benchmark1.json");
	const std::string& configFilePath2 = std::string("/data/model_CSTR_MicroKineticsSMA_benchmark1.json");

	const double absTol = 1e-12;
	const double relTol = 5e-4;

	cadet::test::reaction::testMichaelisMentenToSMAMicroKinetic(configFilePath1, configFilePath2, absTol, relTol);
}

TEST_CASE("MichaelisMenten kinetic with two inhibitors and specific mass action law micro-kinetics yield same result", "[MichaelisMenten],[ReactionModel],[Simulation],[Reference],[CI]")
{
	const std::string& configFilePath1 = std::string("/data/model_CSTR_MichaelisMenten_twoInhib_benchmark1.json");
	const std::string& configFilePath2 = std::string("/data/model_CSTR_MicroKineticsSMA_twoInhib_benchmark1.json");

	const double absTol = 1e-3;
	const double relTol = 5e-4;

	cadet::test::reaction::testMichaelisMentenToSMAInhibitionMicroKinetic(configFilePath1, configFilePath2, absTol, relTol);
}

TEST_CASE("MichaelisMenten kinetic and numerical reference with Crank-Nicolson yield same result", "[MichaelisMenten],[ReactionModel],[Simulation],[Reference],[CI]")
{
	const std::string& configFileRelPath = std::string("/data/model_CSTR_MichaelisMenten_benchmark2.json");
	const std::string& refFileRelPath = std::string("/data/ref_CSTR_MichaelisMenten_benchmark2.h5");

	const double absTol = 1e-3;
	const double relTol = 1e-4;

	cadet::test::column::testForeignReferenceBenchmark(configFileRelPath, refFileRelPath, "000", absTol, relTol, 1);
}


TEST_CASE("MichaelisMenten kinetic analytic Jacobian vs AD without inhibition", "[MichaelisMenten],[ReactionModel],[Jacobian],[AD]")
{
	const unsigned int nBound[] = {1, 2, 1};
	const double point[] = {1.0, 2.0, 1.4, 2.1, 0.2, 1.1, 1.8};
	cadet::test::reaction::testDynamicJacobianAD("MICHAELIS_MENTEN", 3, nBound,
		R"json({
			"MM_KMM": [1.0, 2.0, 0.4],
			"MM_KI": [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
			"MM_VMAX": [1.0, 0.2, 1.5],
			"MM_STOICHIOMETRY_BULK": [ 1.0, -2.0,  3.0,
			                          -1.0,  0.0, -2.0,
			                           0.0,  1.0,  1.0]
		})json",
		point, 1e-15, 1e-15
	);
}

TEST_CASE("MichaelisMenten kinetic analytic Jacobian vs AD with inhibition", "[MichaelisMenten],[ReactionModel],[Jacobian],[AD]")
{
	const unsigned int nBound[] = {1.0, 2.0, 1.0};
	const double point[] = {1.0, 2.0, 1.4, 2.1, 0.2, 1.1, 1.8};
	cadet::test::reaction::testDynamicJacobianAD("MICHAELIS_MENTEN", 3, nBound,
		R"json({
			"MM_KMM": [1.0, 2.0, 0.4],
			"MM_KI": [-1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 3.0, 2.0, -1.0],
			"MM_VMAX": [1.0, 0.2, 1.5],
			"MM_STOICHIOMETRY_BULK": [ 1.0, -2.0,  3.0,
			                          -1.0,  0.0, -2.0,
			                           0.0,  1.0,  1.0]
		})json",
		point, 1e-15, 1e-15
	);
}
TEST_CASE("ASM3 analytic Jacobian vs AD", "[ASM3],[ReactionModel],[Jacobian],[AD],[testHere]")
{
	const unsigned int nBound[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0 };
	unsigned int ncomp = 13;
	const double point[] = { 1.0, 2.0, 1.4, 2.1, 0.2, 1.1, 1.8, 1.5, 1.0, 4.2, 1.4, 0.3, 1.4};
	cadet::test::reaction::testDynamicJacobianAD("ACTIVATED_SLUDGE_MODEL3", ncomp, nBound,
		R"json({
			"ASM3_FISS_BM_PROD": 1.0,
			"ASM3_FSI": 0.0,
			"ASM3_YH_AER": 0.8,
			"ASM3_YH_ANOX": 0.65,
			"ASM3_YSTO_AER": 0.8375,
			"ASM3_YSTO_ANOX": 0.7,
			"ASM3_FXI": 0.2,
			"ASM3_YA": 0.24,
			"ASM3_KH20": 9.0,
			"ASM3_KX": 1.0,
			"ASM3_KSTO20": 12.0,
			"ASM3_MU_H20": 3.0,
			"ASM3_BH20": 0.33,
			"ASM3_ETA_HNO3": 0.5,
			"ASM3_KHO2": 0.2,
			"ASM3_KHSS": 10.0,
			"ASM3_KHNO3": 0.5,
			"ASM3_KHNH4": 0.01,
			"ASM3_KHALK": 0.1,
			"ASM3_KHSTO": 0.1,
			"ASM3_MU_AUT20": 1.12,
			"ASM3_BAUT20": 0.18,
			"ASM3_ETAH_END": 0.5,
			"ASM3_ETAN_END": 0.5,
			"ASM3_KNO2": 0.5,
			"ASM3_KNNH4": 0.7,
			"ASM3_KNALK": 0.5,
			"ASM3_T": 12.0,
			"ASM3_V": 1000.0,
			"ASM3_IO2": 0.0,
			"ASM3_INSI": 0.01,
			"ASM3_INSS": 0.03,
			"ASM3_INXI": 0.04,
			"ASM3_INXS": 0.03,
			"ASM3_INBM": 0.07,
			"ASM3_IVSS_XI": 0.751879699,
			"ASM3_IVSS_XS": 0.555555556,
			"ASM3_IVSS_STO": 0.6,
			"ASM3_IVSS_BM": 0.704225352,
			"ASM3_ITSS_VSS_BM": 1.086956522
		})json",
		point, 1e-15, 1e-15
	);
}
