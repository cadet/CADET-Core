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
	cadet::test::reaction::testDynamicJacobianAD("MASS_ACTION_LAW_CROSS_PHASE", 3, nBound,
		R"json({
			"MAL_KFWD_LIQUID": [1.0, 2.0, 0.4],
			"MAL_KBWD_LIQUID": [0.0, 0.2, 1.5],
			"MAL_STOICHIOMETRY_LIQUID": [ 1.0, -2.0,  3.0,
			                           -1.0,  0.0, -2.0,
			                            0.0,  1.0,  1.0],
			"MAL_EXPONENTS_LIQUID_FWD": [ 1.2,  0.0,  0.0,
			                            0.0,  1.3,  2.2,
			                            0.0,  1.0,  1.1],
			"MAL_EXPONENTS_LIQUID_BWD": [ 0.8,  2.1,  1.0,
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
TEST_CASE("CSTR with MAL reaction numerical Benchmark with parameter sensitivities", "[CSTR],[MassActionLaw],[ReactionModel],[Simulation],[Reference],[Sensitivity],[CI_sens16]")
{
	std::string modelFilePath = std::string("/data/model_CSTR_reacMAL_2comp_sensbenchmark1.json");
	std::string refFilePath = std::string("/data/ref_CSTR_reacMAL_2comp_sensbenchmark1.h5");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1e-6, 1e-6, 1e-6 };

	cadet::test::column::DummyParams disc;
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("MichaelisMenten kinetic and specific mass action law micro-kinetics yield same result", "[CSTR],[MassActionLaw],[MichaelisMenten],[ReactionModel],[Simulation],[CI]")
{
	const std::string& configFilePath1 = std::string("/data/model_CSTR_MichaelisMenten_benchmark1.json");
	const std::string& configFilePath2 = std::string("/data/model_CSTR_MicroKineticsSMA_benchmark1.json");

	const double absTol = 1e-12;
	const double relTol = 5e-4;

	cadet::test::reaction::testMichaelisMentenToSMAMicroKinetic(configFilePath1, configFilePath2, absTol, relTol);
}

TEST_CASE("MichaelisMenten kinetic with two substrates and specific mass action law micro-kinetics yield same result", "[CSTR],[MassActionLaw],[MichaelisMenten],[ReactionModel],[Simulation],[Reference],[CI]")
{
	const std::string& configFilePath1 = std::string("/data/configuration_CSTR_MichaelisMenten_twoSubs_benchmark1.json");
	const std::string& configFilePath2 = std::string("/data/configuration_CSTR_MicroKineticsSMA_twoSubs_benchmark1.json");

	const double absTol = 1e-3;
	const double relTol = 5e-4;

	cadet::test::reaction::testMichaelisMentenToSMAMicroKinetic(configFilePath1, configFilePath2, absTol, relTol);
}

TEST_CASE("MichaelisMenten kinetic with two non-inhibitors and specific mass action law micro-kinetics yield same result", "[CSTR],[MassActionLaw],[MichaelisMenten],[ReactionModel],[Simulation],[Reference],[CI]")
{
	const std::string& configFilePath1 = std::string("/data/model_CSTR_MichaelisMenten_twoInhib_benchmark1.json");
	const std::string& configFilePath2 = std::string("/data/model_CSTR_MicroKineticsSMA_twoInhib_benchmark1.json");

	const double absTol = 1e-3;
	const double relTol = 5e-4;

	cadet::test::reaction::testMichaelisMentenToSMAMicroKinetic(configFilePath1, configFilePath2, absTol, relTol);
}


TEST_CASE("MichaelisMenten kinetic with two non-inhibitors and two substrates and specific mass action law micro-kinetics yield same result", "[CSTR],[MassActionLaw],[MichaelisMenten],[ReactionModel],[Simulation],[Reference],[CI]")
{
	const std::string& configFilePath1 = std::string("/data/configuration_CSTR_MichaelisMenten_twoSubs_twoInhib_benchmark1.json");
	const std::string& configFilePath2 = std::string("/data/configuration_CSTR_MichaelisMenten_twoSubs_twoInhib_benchmark1.json");

	const double absTol = 1e-3;
	const double relTol = 5e-4;

	cadet::test::reaction::testMichaelisMentenToSMAMicroKinetic(configFilePath1, configFilePath2, absTol, relTol);
}

TEST_CASE("MichaelisMenten kinetic and numerical reference with Crank-Nicolson yield same result", "[CSTR],[MassActionLaw],[MichaelisMenten],[ReactionModel],[Simulation],[Reference],[CI]")
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
			"MM_STOICHIOMETRY": [ 1.0, -2.0,  3.0,
			                          -1.0,  0.0, -2.0,
			                           0.0,  1.0,  1.0]
		})json",
		point, 1e-15, 1e-15
	);
}

TEST_CASE("MichaelisMenten kinetic analytic Jacobian vs AD with inhibition", "[MichaelisMenten],[ReactionModel],[Jacobian],[AD]")
{
	const unsigned int nBound[] = {1, 2, 1};
	const double point[] = {1.0, 2.0, 1.4, 2.1, 0.2, 1.1, 1.8};
	cadet::test::reaction::testDynamicJacobianAD("MICHAELIS_MENTEN", 3, nBound,
		R"json({
			"MM_KMM": [1.0, 2.0, 0.4],
			"MM_KI": [-1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 3.0, 2.0, -1.0],
			"MM_VMAX": [1.0, 0.2, 1.5],
			"MM_STOICHIOMETRY": [ 1.0, -2.0,  3.0,
			                          -1.0,  0.0, -2.0,
			                           0.0,  1.0,  1.0]
		})json",
		point, 1e-15, 1e-15
	);
}

TEST_CASE("MassActionLaw old interface vs. two separate reactions", "[MassActionLaw],[ReactionModel],[Simulation],[CI]")
{
	std::string modelFilePath = std::string("/data/model_CSTR_reacMAL_3comp_nreac_2.json");
	std::string refFilePath = std::string("/data/ref_CSTR_reacMAL_3comp_one_type_old_interface.h5");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1e-6, 1e-6, 1e-6 };

	cadet::test::column::DummyParams disc;
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("MassActionLaw one reaction vs. two separate reactions", "[MassActionLaw],[ReactionModel],[Simulation],[CI]")
{
	std::string modelFilePath = std::string("/data/model_CSTR_reacMAL_3comp_nreac_2.json");
	std::string refFilePath = std::string("/data/ref_CSTR_reacMAL_3comp_nreac_1.h5");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1e-6, 1e-6, 1e-6 };

	cadet::test::column::DummyParams disc;
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}
