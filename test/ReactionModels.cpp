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
#include "ReactionModelFactory.hpp"
#include "common/JsonParameterProvider.hpp"
#include "model/reaction/ConservedMoieties.hpp"

#include <cmath>
#include <memory>
#include <utility>
#include <vector>

TEST_CASE("MassActionLaw conserved moieties nullspace", "[MassActionLaw],[ReactionModel],[ConservedMoieties],[CI]")
{
	SECTION("A <-> B")
	{
		const Eigen::MatrixXd L = cadet::test::reaction::conservedMoietiesFromMassActionLaw(2,
			R"json({
				"MAL_KFWD": [1.0],
				"MAL_KBWD": [2.0],
				"MAL_STOICHIOMETRY": [-1.0,
				                       1.0],
				"MAL_IS_KINETIC": [0]
			})json");

		Eigen::MatrixXd S(2, 1);
		S << -1.0, 1.0;

		REQUIRE(L.rows() == 1);
		REQUIRE(L.cols() == 2);
		cadet::test::reaction::checkNullSpace(L, S);
	}

	SECTION("A + B <-> C")
	{
		const Eigen::MatrixXd L = cadet::test::reaction::conservedMoietiesFromMassActionLaw(3,
			R"json({
				"MAL_KFWD": [1.0],
				"MAL_KBWD": [2.0],
				"MAL_STOICHIOMETRY": [-1.0,
				                      -1.0,
				                       1.0],
				"MAL_IS_KINETIC": [0]
			})json");

		Eigen::MatrixXd S(3, 1);
		S << -1.0, -1.0, 1.0;

		REQUIRE(L.rows() == 2);
		REQUIRE(L.cols() == 3);
		cadet::test::reaction::checkNullSpace(L, S);
		
	}
	SECTION("A + B <-> C (eq) and A <-> B (kin)")
	{
		const Eigen::MatrixXd L = cadet::test::reaction::conservedMoietiesFromMassActionLaw(3,
			R"json({
				"MAL_KFWD": [1.0,1.0],
				"MAL_KBWD": [2.0,2.0],
				"MAL_STOICHIOMETRY": [-1.0,-1.0,
				                      -1.0, 1.0,
				                       1.0, 0.0],
				"MAL_IS_KINETIC": [0,1]
			})json");

		Eigen::MatrixXd S(3, 1);
		S << -1.0, -1.0, 1.0;

		REQUIRE(L.rows() == 2);
		REQUIRE(L.cols() == 3);
		cadet::test::reaction::checkNullSpace(L, S);

	}

	SECTION("No equilibrium MAL reaction")
	{
		const Eigen::MatrixXd L = cadet::test::reaction::conservedMoietiesFromMassActionLaw(3,
			R"json({
				"MAL_KFWD": [1.0],
				"MAL_KBWD": [2.0],
				"MAL_STOICHIOMETRY": [-1.0,
				                      -1.0,
				                       1.0],
				"MAL_IS_KINETIC": [1]
			})json");

 		Eigen::MatrixXd expected;
		expected = Eigen::MatrixXd::Identity(3, 3);
		
		const double tol = 1e-12;

		REQUIRE(L.rows() == expected.rows());
		REQUIRE(L.cols() == expected.cols());

		for (Eigen::Index r = 0; r < L.rows(); ++r)
		{
			for (Eigen::Index c = 0; c < L.cols(); ++c)
				CHECK(std::abs(L(r, c) - expected(r, c)) <= tol);
		}
	}
}

TEST_CASE("Conserved moiety transformations", "[ReactionModel],[ConservedMoieties],[CI]")
{
	constexpr double tol = 1e-12;

	Eigen::MatrixXd stoichiometry(3, 1);
	stoichiometry << -1.0, -1.0, 1.0;
	const std::vector<bool> equilibriumReactionFlags{true};

	cadet::model::ConservedMoieties cm;
	REQUIRE(cm.configure(3, equilibriumReactionFlags, stoichiometry, 1e-14));
	REQUIRE(cm.numMoieties() == 2);
	REQUIRE(cm.numEquilibriumReactions() == 1);

	const double source[] = {1.0, 2.0, 4.0};
	double expected[2] = {0.0, 0.0};
	const Eigen::MatrixXd& L = cm.conservedMoietyMatrix();
	for (unsigned int moiety = 0; moiety < cm.numMoieties(); ++moiety)
	{
		for (unsigned int state = 0; state < 3; ++state)
			expected[moiety] += L(moiety, state) * source[state];
	}

	SECTION("Vector")
	{
		double result[2];
		cm.applyToVector(result, source, 3);
		CHECK(std::abs(result[0] - expected[0]) <= tol);
		CHECK(std::abs(result[1] - expected[1]) <= tol);

	}

	SECTION("Derivative vector")
	{
		double result[3];
		cm.applyToDerivativeVector(result, source, 3);
		CHECK(std::abs(result[0] - expected[0]) <= tol);
		CHECK(std::abs(result[1] - expected[1]) <= tol);
		CHECK(result[2] == 0.0);
	}

	SECTION("Sparse matrix and pattern")
	{
		std::vector<Eigen::Triplet<double>> entries;
		entries.emplace_back(0, 0, 9.0);
		entries.emplace_back(1, 0, 1.0);
		entries.emplace_back(1, 3, 2.0);
		entries.emplace_back(2, 1, 3.0);
		entries.emplace_back(3, 2, 4.0);
		const Eigen::MatrixXd original = (Eigen::MatrixXd(3, 4) <<
			1.0, 0.0, 0.0, 2.0,
			0.0, 3.0, 0.0, 0.0,
			0.0, 0.0, 4.0, 0.0).finished();

		cm.addPatternToBlocks(entries, 3, 1, 1, 3, 0, 4);
		Eigen::SparseMatrix<double, Eigen::RowMajor> matrix(5, 4);
		matrix.setFromTriplets(entries.begin(), entries.end());

		const std::size_t bufferSize = cm.matrixBufferSize(matrix, 3, 1, 0, matrix.cols());
		std::vector<Eigen::Triplet<double>> buffer(bufferSize);
		CHECK(bufferSize > 0);
		cm.applyToMatrix(matrix, 3, 1, 0, matrix.cols(), buffer.data(), buffer.size());

		for (unsigned int moiety = 0; moiety < cm.numMoieties(); ++moiety)
		{
			for (unsigned int column = 0; column < 4; ++column)
			{
				double transformed = 0.0;
				for (unsigned int state = 0; state < 3; ++state)
					transformed += L(moiety, state) * original(state, column);
				CHECK(std::abs(matrix.coeff(1 + moiety, column) - transformed) <= tol);
			}
		}
		CHECK(matrix.coeff(0, 0) == 9.0);
	}

	SECTION("Repeated sparse pattern")
	{
		std::vector<Eigen::Triplet<double>> entries;
		entries.emplace_back(1, 0, 1.0);
		entries.emplace_back(5, 3, 2.0);
		cm.addPatternToBlocks(entries, 3, 1, 2, 4, 0, 4);

		for (unsigned int moiety = 0; moiety < cm.numMoieties(); ++moiety)
		{
			bool firstBlockEntryFound = false;
			bool secondBlockEntryFound = false;
			for (const auto& entry : entries)
			{
				firstBlockEntryFound = firstBlockEntryFound || ((entry.row() == 1 + moiety) && (entry.col() == 0));
				secondBlockEntryFound = secondBlockEntryFound || ((entry.row() == 5 + moiety) && (entry.col() == 3));
			}
			CHECK(firstBlockEntryFound);
			CHECK(secondBlockEntryFound);
		}
	}
}

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
