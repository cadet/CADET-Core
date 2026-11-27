// =============================================================================
//  CADET
//  
//  Copyright Â© 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include <catch.hpp>

#include "ReactionModelTests.hpp"
#include "BindingModelTests.hpp"
#include "JsonTestModels.hpp"
#include "ParticleHelper.hpp"
#include "SimHelper.hpp"
#include "UnitOperationTests.hpp"
#include "ParallelSupport.hpp"
#include "Utils.hpp"
#include "Approx.hpp"
#include "cadet/FactoryFuncs.hpp"
#include "cadet/ModelBuilder.hpp"
#include "model/UnitOperation.hpp"
#include "ColumnTests.hpp"
#include "ReactionModelFactory.hpp"
#include "SimulationTypes.hpp"
#include "cadet/Exceptions.hpp"
#include "common/JsonParameterProvider.hpp"
#include <json.hpp>
#include "linalg/DenseMatrix.hpp"
#include "model/reaction/ConservedMoieties.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
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
				"IS_KINETIC": [0]
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
				"IS_KINETIC": [0]
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
				"IS_KINETIC": [0,1]
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
				"IS_KINETIC": [1]
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
	const unsigned int nBound[] = { 1, 2, 1 };
	const double point[] = { 1.0, 2.0, 1.4, 2.1, 0.2, 1.1, 1.8 };
	cadet::test::reaction::testDynamicJacobianAD("MASS_ACTION_LAW_CROSS_PHASE", 3, nBound,
		R"json({
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
TEST_CASE("CSTR with MAL reaction numerical Benchmark with parameter sensitivities", "[CSTR],[MassActionLaw],[ReactionModel],[Simulation],[Reference],[Sensitivity],[CI_sens]")
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
	const std::string& configFilePath2 = std::string("/data/model_CSTR_MicroKineticsMAL_benchmark1.json");

	const double absTol = 1e-12;
	const double relTol = 5e-4;

	cadet::test::reaction::testMichaelisMentenToSMAMicroKinetic(configFilePath1, configFilePath2, absTol, relTol);
}

TEST_CASE("MichaelisMenten kinetic with two substrates and specific mass action law micro-kinetics yield same result", "[CSTR],[MassActionLaw],[MichaelisMenten],[ReactionModel],[Simulation],[Reference],[CI]")
{
	const std::string& configFilePath1 = std::string("/data/config_CSTR_MichaelisMenten_twoSubs_benchmark1.json");
	const std::string& configFilePath2 = std::string("/data/config_CSTR_MicroKineticsMAL_twoSubs_benchmark1.json");

	const double absTol = 1e-3;
	const double relTol = 5e-4;

	cadet::test::reaction::testMichaelisMentenToSMAMicroKinetic(configFilePath1, configFilePath2, absTol, relTol);
}

TEST_CASE("MichaelisMenten kinetic with two non-inhibitors and specific mass action law micro-kinetics yield same result", "[CSTR],[MassActionLaw],[MichaelisMenten],[ReactionModel],[Simulation],[Reference]")
{
	const std::string& configFilePath1 = std::string("/data/model_CSTR_MichaelisMenten_twoInhib_benchmark1.json");
	const std::string& configFilePath2 = std::string("/data/model_CSTR_MicroKineticsMAL_twoInhib_benchmark1.json");

	const double absTol = 1e-3;
	const double relTol = 5e-4;

	cadet::test::reaction::testMichaelisMentenToSMAMicroKinetic(configFilePath1, configFilePath2, absTol, relTol);
}


TEST_CASE("MichaelisMenten kinetic with two non-inhibitors and two substrates and specific mass action law micro-kinetics yield same result", "[CSTR],[MassActionLaw],[MichaelisMenten],[ReactionModel],[Simulation]")// fails with 0.9999999983 == RelApprox( 0.99 )
{
	const std::string& configFilePath1 = std::string("/data/config_CSTR_MichaelisMenten_twoSubs_twoInhib_benchmark1.json");
	const std::string& configFilePath2 = std::string("/data/config_CSTR_MicroKineticsMAL_twoSubs_twoInhib_benchmark1.json");

	const double absTol = 1e-3;
	const double relTol = 5e-3;

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


TEST_CASE("MichaelisMenten kinetic analytic Jacobian vs AD without inhibition", "[MichaelisMenten],[ReactionModel],[Jacobian],[AD],[tofix]") // fails with 0.0 <= 0.0
{
	const unsigned int nBound[] = { 1, 2, 1 };
	const double point[] = { 1.0, 2.0, 1.4, 2.1, 0.2, 1.1, 1.8 };
	cadet::test::reaction::testDynamicJacobianAD("MICHAELIS_MENTEN", 3, nBound,
		R"json({
			"MM_KMM": [0.0, 1.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.4, 0.0],
			"MM_KI_UC": [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
			"MM_KI_IC": [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
			"MM_VMAX": [1.0, 0.2, 1.5],
			"MM_STOICHIOMETRY": [ 1.0, -2.0,  3.0,
			                          -1.0,  0.0, -2.0,
			                           0.0,  1.0,  1.0]
		})json",
		point, 1e-15, 1e-15
	);
}

TEST_CASE("MichaelisMenten kinetic analytic Jacobian vs AD with prefactorial component", "[MichaelisMenten],[ReactionModel],[Jacobian],[AD]")
{
	const unsigned int nBound[] = {1, 2, 1};
	const double point[] = {1.0, 2.0, 1.4, 2.1, 0.2, 1.1, 1.8};
	cadet::test::reaction::testDynamicJacobianAD("MICHAELIS_MENTEN", 3, nBound,
		R"json({
			"MM_KMM": [0.0, 1.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.4, 0.0],
			"MM_KI_UC": [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
			"MM_KI_IC": [-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0],
			"MM_VMAX": [1.0, 0.2, 1.5],
			"MM_STOICHIOMETRY": [ 1.0, -2.0,  3.0,
			                          -1.0,  0.0, -2.0,
			                           0.0,  1.0,  1.0],
			"MM_PRE_K": [0.0, 0.0, 2.0, 0.0, 5.0, 0.0, 0.5, 0.0, 0.0]
		})json",
		point, 1e-15, 1e-15
	);
}

TEST_CASE("MichaelisMenten kinetic analytic Jacobian vs AD with inhibition", "[MichaelisMenten],[ReactionModel],[Jacobian],[AD]") // fails with expansion: 0.0 <= 0.0
{
	const unsigned int nBound[] = { 1, 2, 1 };
	const double point[] = { 1.0, 2.0, 1.4, 2.1, 0.2, 1.1, 1.8 };
	cadet::test::reaction::testDynamicJacobianAD("MICHAELIS_MENTEN", 3, nBound,
		R"json({
			"MM_KMM": [0.0, 1.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.4, 0.0],
			"MM_KI_UC": [-1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 3.0, 2.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 3.0, 2.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 3.0, 2.0, -1.0],
			"MM_KI_C": [-1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 3.0, 2.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 3.0, 2.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 3.0, 2.0, -1.0],
			"MM_VMAX": [1.0, 0.2, 1.5],
			"MM_STOICHIOMETRY": [ 1.0, -2.0,  3.0,
			                          -1.0,  0.0, -2.0,
			                           0.0,  1.0,  1.0]
		})json",
		point, 1e-15, 1e-15
	);
}
TEST_CASE("ASM3 analytic Jacobian vs AD with adsorbed components", "[ASM3],[ReactionModel],[Jacobian],[AD],[CI]")
{
	const unsigned int nBound[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0, 0,0,0 };
	unsigned int ncomp = 15;
	const double point[] = { 1.0, 2.0, 1.4, 2.1, 0.2, 1.1, 1.8, 1.5, 1.0, 4.2, 1.4, 0.3, 1.4, 1.0, 2.0 };
	cadet::test::reaction::testLiquidReactionJacobianAD("ACTIVATED_SLUDGE_MODEL3", ncomp, nBound,
		R"json({
			"ASM3_FISS_BM_PROD": 1.0,
			"ASM3_FSI": 0.0,
			"ASM3_FSI_AD": 0.5,
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
		point, 1e-12, 1e-12
	);
}

TEST_CASE("ASM3 analytic Jacobian vs AD without adsorbed components", "[ASM3],[ReactionModel],[Jacobian],[AD],[CI]")
{
	const unsigned int nBound[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	unsigned int ncomp = 13;
	const double point[] = { 1.0, 2.0, 1.4, 2.1, 0.2, 1.1, 1.8, 1.5, 1.0, 4.2, 1.4, 0.3, 1.4 };
	cadet::test::reaction::testLiquidReactionJacobianAD("ACTIVATED_SLUDGE_MODEL3", ncomp, nBound,
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
			"ASM3_IVSS_BM": 0.704225352,
			"ASM3_ITSS_VSS_BM": 1.086956522
		})json",
		point, 1e-12, 1e-12
	);
}

TEST_CASE("ASM3 adsorbable substrate fraction reproduces the plain setup", "[ASM3],[ReactionModel],[Simulation],[CI]")
{
	// The 15 component setup routes the complete substrate into the adsorbable fraction
	// (ASM3_FSS_AD = 1, SS starts at 0, SS_ad starts at 2), so it has to reproduce the
	// 13 component setup. SS_ad (index 14) of the first setup corresponds to SS (index 1)
	// of the second one, all remaining components correspond to themselves.
	const std::vector<int> compIdx1 = { 0, 14, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
	const std::vector<int> compIdx2 = { 0,  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

	// The tolerance is set by the two independent time integrations, which both run at
	// ABSTOL = RELTOL = 1e-6, not by the reaction model. A wrong stoichiometric coefficient
	// shows up several orders of magnitude above this.
	cadet::test::reaction::testCompareTwoSimulationReaction(
		"/data/configuration_CSTR_ASM3_with_ad.json",
		"/data/configuration_CSTR_ASM3_without_ad.json",
		1e-3, 1e-4, compIdx1, compIdx2);
}

TEST_CASE("MichaelisMenten Hill exponent changes flux and analytic Jacobian", "[MichaelisMenten],[Hill],[ReactionModel],[Jacobian],[CI]")
{
	const unsigned int nBound[] = { 0, 0 };
	const unsigned int boundOffset[] = { 0, 0 };
	cadet::JsonParameterProvider config(R"json({"MM_KMM": [2.0, 0.0], "MM_VMAX": [3.0], "MM_STOICHIOMETRY": [-1.0, 1.0], "MM_EXP": [2.0, 0.0]})json");
	cadet::ReactionModelFactory factory;
	std::unique_ptr<cadet::model::IDynamicReactionModel> model(factory.createDynamic("MICHAELIS_MENTEN"));
	REQUIRE(model != nullptr);
	REQUIRE(model->configureModelDiscretization(config, 2, nBound, boundOffset));
	REQUIRE(model->configure(config, 0, 0));

	const unsigned int workspaceSize = model->workspaceSize(2, 0, nBound);
	std::vector<char> workspaceMemory(workspaceSize);
	cadet::LinearBufferAllocator workspace(workspaceMemory.data(), workspaceMemory.data() + workspaceSize);
	const double y[] = { 1.0, 0.0 };
	double residual[] = { 0.0, 0.0 };
	model->residualFluxAdd(0.0, 0, cadet::ColumnPosition{ 0.0, 0.0, 0.0 }, 2, y, residual, 1.0, workspace);
	CHECK(residual[0] == Approx(-0.6).margin(1e-12));
	CHECK(residual[1] == Approx(0.6).margin(1e-12));

	cadet::linalg::DenseMatrix jac;
	jac.resize(2, 2);
	jac.setAll(0.0);
	model->analyticJacobianAdd(0.0, 0, cadet::ColumnPosition{ 0.0, 0.0, 0.0 }, 2, y, 1.0, jac.row(0), workspace);
	CHECK(jac.native(0, 0) == Approx(-0.96).margin(1e-12));
	CHECK(jac.native(1, 0) == Approx(0.96).margin(1e-12));
}

TEST_CASE("MichaelisMenten external variant applies inhibition like the plain model", "[MichaelisMenten],[ReactionModel],[ExternalFunction],[CI]")
{
	// One reaction 0 -> 1, component 2 inhibits substrate 0 competitively and uncompetitively
	const unsigned int nBound[] = { 0, 0, 0 };
	const nlohmann::json base = {
		{"MM_VMAX", {3.0}},
		{"MM_KMM", {2.0, 0.0, 0.0}},
		{"MM_KI_C", {-1.0, -1.0, 0.5, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0}},
		{"MM_KI_UC", {-1.0, -1.0, 0.8, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0}},
		{"MM_PRE_K", {0.0, 0.0, 0.0}},
		{"MM_EXP", {1.0, 0.0, 0.0}}
	};
	nlohmann::json plain = base;
	nlohmann::json external;
	for (auto it = base.begin(); it != base.end(); ++it)
	{
		const std::string key = "EXT_" + it.key();
		const nlohmann::json zero = std::vector<double>(it.value().size(), 0.0);
		external[key] = it.value();
		external[key + "_T"] = zero;
		external[key + "_TT"] = zero;
		external[key + "_TTT"] = zero;
	}
	plain["MM_STOICHIOMETRY"] = {-1.0, 1.0, 0.0};
	external["MM_STOICHIOMETRY"] = {-1.0, 1.0, 0.0};

	const double y[] = { 1.0, 0.0, 2.0 };
	double resPlain[] = { 0.0, 0.0, 0.0 };
	double resExternal[] = { 0.0, 0.0, 0.0 };

	auto crmPlain = cadet::test::reaction::ConfiguredDynamicReactionModel::create("MICHAELIS_MENTEN", 3, nBound, plain.dump().c_str());
	auto crmExternal = cadet::test::reaction::ConfiguredDynamicReactionModel::create("EXT_MICHAELIS_MENTEN", 3, nBound, external.dump().c_str());
	crmPlain.model().residualFluxAdd(1.0, 0u, cadet::ColumnPosition{ 0.0, 0.0, 0.0 }, 3, y, resPlain, 1.0, crmPlain.buffer());
	crmExternal.model().residualFluxAdd(1.0, 0u, cadet::ColumnPosition{ 0.0, 0.0, 0.0 }, 3, y, resExternal, 1.0, crmExternal.buffer());

	// The plain model is inhibited: its rate is below the uninhibited 3 * 1 / (2 + 1) = 1
	CHECK(std::abs(resPlain[1]) < 0.99);
	for (int i = 0; i < 3; ++i)
	{
		INFO("Component " << i);
		CHECK(resExternal[i] == Approx(resPlain[i]).margin(1e-14));
	}
}

TEST_CASE("ASM3 external variant reads its parameters like the plain model", "[ASM3],[ReactionModel],[ExternalFunction],[CI]")
{
	// Parameters that can depend on external functions are read with the EXT_ prefix, the others without
	const nlohmann::json external_capable = {
		{"ASM3_KH20", 9.0}, {"ASM3_T", 12.0}, {"ASM3_IO2", 0.5}, {"ASM3_V", 1000.0}, {"ASM3_KSTO20", 12.0},
		{"ASM3_KX", 1.0}, {"ASM3_KHO2", 0.2}, {"ASM3_KHSS", 10.0}, {"ASM3_KHNO3", 0.5}, {"ASM3_ETA_HNO3", 0.5},
		{"ASM3_KHNH4", 0.01}, {"ASM3_KHALK", 0.1}, {"ASM3_KHSTO", 0.1}, {"ASM3_MU_H20", 3.0}, {"ASM3_ETAH_END", 0.5},
		{"ASM3_BH20", 0.33}, {"ASM3_MU_AUT20", 1.12}, {"ASM3_KNO2", 0.5}, {"ASM3_KNNH4", 0.7}, {"ASM3_KNALK", 0.5},
		{"ASM3_BAUT20", 0.18}, {"ASM3_ETAN_END", 0.5}
	};
	const nlohmann::json common = {
		{"ASM3_FISS_BM_PROD", 1.0}, {"ASM3_FSI", 0.0}, {"ASM3_YH_AER", 0.8}, {"ASM3_YH_ANOX", 0.65},
		{"ASM3_YSTO_AER", 0.8375}, {"ASM3_YSTO_ANOX", 0.7}, {"ASM3_FXI", 0.2}, {"ASM3_YA", 0.24},
		{"ASM3_INSI", 0.01}, {"ASM3_INSS", 0.03}, {"ASM3_INXI", 0.04}, {"ASM3_INXS", 0.03}, {"ASM3_INBM", 0.07},
		{"ASM3_IVSS_XI", 0.751879699}, {"ASM3_IVSS_XS", 0.555555556}, {"ASM3_IVSS_STO", 0.6},
		{"ASM3_IVSS_BM", 0.704225352}, {"ASM3_ITSS_VSS_BM", 1.086956522}
	};
	nlohmann::json plain = common;
	nlohmann::json external = common;
	for (auto it = external_capable.begin(); it != external_capable.end(); ++it)
	{
		plain[it.key()] = it.value();
		const std::string key = "EXT_" + it.key();
		external[key] = it.value();
		external[key + "_T"] = 0.0;
		external[key + "_TT"] = 0.0;
		external[key + "_TTT"] = 0.0;
	}

	const unsigned int nBound[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	const double y[] = { 1.0, 2.0, 1.4, 2.1, 0.2, 1.1, 1.8, 1.5, 1.0, 4.2, 1.4, 0.3, 1.4 };
	double resPlain[13] = {};
	double resExternal[13] = {};

	auto crmPlain = cadet::test::reaction::ConfiguredDynamicReactionModel::create("ACTIVATED_SLUDGE_MODEL3", 13, nBound, plain.dump().c_str());
	auto crmExternal = cadet::test::reaction::ConfiguredDynamicReactionModel::create("EXT_ACTIVATED_SLUDGE_MODEL3", 13, nBound, external.dump().c_str());
	crmPlain.model().residualFluxAdd(1.0, 0u, cadet::ColumnPosition{ 0.0, 0.0, 0.0 }, 13, y, resPlain, 1.0, crmPlain.buffer());
	crmExternal.model().residualFluxAdd(1.0, 0u, cadet::ColumnPosition{ 0.0, 0.0, 0.0 }, 13, y, resExternal, 1.0, crmExternal.buffer());

	for (int i = 0; i < 13; ++i)
	{
		INFO("Component " << i);
		CHECK(resExternal[i] == Approx(resPlain[i]).epsilon(1e-14).margin(1e-14));
	}
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


TEST_CASE("IS_KINETIC is rejected by reaction models that do not report reactions", "[MassActionLaw],[ReactionModel],[CI]")
{
	// A cross phase mass action law couples a liquid and a solid reaction and does not report
	// individual reactions, so rapid equilibrium cannot be honoured for it. Refuse the request
	// instead of silently ignoring it.
	cadet::ReactionModelFactory rmf;
	std::unique_ptr<cadet::model::IDynamicReactionModel> model(rmf.createDynamic("MASS_ACTION_LAW_CROSS_PHASE"));
	REQUIRE(model);
	REQUIRE(model->numReactions() == 0u);

	const unsigned int nComp = 2;
	const std::vector<unsigned int> nBound(nComp, 1u);
	const std::vector<unsigned int> boundOffset{ 0u, 1u };

	cadet::JsonParameterProvider jpp(R"json({
		"MAL_KFWD_LIQUID": [1.0],
		"MAL_KBWD_LIQUID": [1.0],
		"MAL_KFWD_SOLID": [1.0],
		"MAL_KBWD_SOLID": [1.0],
		"MAL_STOICHIOMETRY_LIQUID": [-1.0, 1.0],
		"MAL_STOICHIOMETRY_SOLID": [-1.0, 1.0],
		"IS_KINETIC": [0]
	})json");

	REQUIRE(model->configureModelDiscretization(jpp, nComp, nBound.data(), boundOffset.data()));
	REQUIRE_THROWS_AS(model->configure(jpp, 0, 0), cadet::InvalidParameterException);
}

namespace
{
	/**
	 * @brief Turns the reaction group @p group of the current scope into its externally dependent variant
	 * @details Each rate constant @f$ k @f$ becomes @f$ k_0 + k_T T @f$ with @f$ k_0 = a k @f$ and
	 *          @f$ k_T = b k @f$. The external function of the test returns the time, so at @f$ t = 1 @f$
	 *          the rates equal those of the plain group if @f$ a + b = 1 @f$, provided the value of the
	 *          external function reaches the reaction.
	 * @return Number of converted groups, 0 or 1
	 */
	int makeReactionGroupExternal(cadet::JsonParameterProvider& jpp, const std::string& group, double constFraction, double slopeFraction)
	{
		if (!jpp.exists(group))
			return 0;

		jpp.pushScope(group);

		const std::string type = jpp.getString("TYPE");
		const std::vector<std::string> rates = (type == "MASS_ACTION_LAW_CROSS_PHASE")
			? std::vector<std::string>{ "MAL_KFWD_LIQUID", "MAL_KBWD_LIQUID", "MAL_KFWD_SOLID", "MAL_KBWD_SOLID" }
			: std::vector<std::string>{ "MAL_KFWD", "MAL_KBWD" };

		jpp.set("TYPE", "EXT_" + type);
		jpp.set("EXTFUN", std::vector<int>{ 0 });

		for (const std::string& name : rates)
		{
			const std::vector<double> k = jpp.getDoubleArray(name);
			std::vector<double> constant(k.size());
			std::vector<double> slope(k.size());
			for (std::size_t i = 0; i < k.size(); ++i)
			{
				constant[i] = constFraction * k[i];
				slope[i] = slopeFraction * k[i];
			}

			jpp.remove(name);
			jpp.set("EXT_" + name, constant);
			jpp.set("EXT_" + name + "_T", slope);
			jpp.set("EXT_" + name + "_TT", std::vector<double>(k.size(), 0.0));
			jpp.set("EXT_" + name + "_TTT", std::vector<double>(k.size(), 0.0));
		}

		jpp.popScope();
		return 1;
	}

	/**
	 * @brief Turns every reaction group of a unit, in the bulk and in each particle type, into its externally dependent variant
	 * @return Number of converted groups
	 */
	int makeReactionsExternal(cadet::JsonParameterProvider& jpp, double constFraction, double slopeFraction)
	{
		int nGroups = makeReactionGroupExternal(jpp, "liquid_reaction_000", constFraction, slopeFraction);

		const int nParType = jpp.exists("NPARTYPE") ? jpp.getInt("NPARTYPE") : 0;
		for (int i = 0; i < nParType; ++i)
		{
			const std::string parType = "particle_type_" + std::string(3 - std::to_string(i).length(), '0') + std::to_string(i);
			if (!jpp.exists(parType))
				continue;

			jpp.pushScope(parType);
			for (const char* group : { "liquid_reaction_000", "solid_reaction_000", "cross_phase_reaction_000" })
				nGroups += makeReactionGroupExternal(jpp, group, constFraction, slopeFraction);
			jpp.popScope();
		}

		return nGroups;
	}

	/**
	 * @brief Adds a liquid phase reaction inside each particle type, which extendModelWithDynamicReactions() does not cover
	 */
	void addParticleLiquidReactions(cadet::JsonParameterProvider& jpp)
	{
		const int nParType = jpp.getInt("NPARTYPE");
		for (int i = 0; i < nParType; ++i)
		{
			jpp.pushScope("particle_type_" + std::string(3 - std::to_string(i).length(), '0') + std::to_string(i));
			jpp.set("NREAC_LIQUID", 1);
			jpp.addScope("liquid_reaction_000");
			jpp.pushScope("liquid_reaction_000");
			jpp.set("TYPE", "MASS_ACTION_LAW");
			jpp.set("MAL_STOICHIOMETRY", std::vector<double>{ -1.0, 1.0 });
			jpp.set("MAL_KFWD", std::vector<double>{ 0.7 });
			jpp.set("MAL_KBWD", std::vector<double>{ 0.2 });
			jpp.popScope();
			jpp.popScope();
		}
	}

	/**
	 * @brief Evaluates the residual of a unit operation at @f$ t = 1 @f$, with the external functions set as a model system sets them
	 */
	std::vector<double> unitResidualAtTimeOne(cadet::JsonParameterProvider jpp)
	{
		cadet::IModelBuilder* const mb = cadet::createModelBuilder();
		REQUIRE(nullptr != mb);

		cadet::IUnitOperation* const unit = cadet::test::unitoperation::createAndConfigureUnit(jpp, *mb);

		cadet::test::binding::LinearExternalFunction extFun;
		cadet::IExternalFunction* extFuns[] = { &extFun };
		unit->setExternalFunctions(extFuns, 1);

		const unsigned int nDof = unit->numDofs();
		std::vector<double> y(nDof);
		std::vector<double> yDot(nDof, 0.0);
		std::vector<double> res(nDof, 0.0);
		cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(0.3 + idx * 0.13)) + 0.1; }, nDof);

		cadet::util::ThreadLocalStorage tls;
		tls.resize(unit->threadLocalMemorySize());

		const cadet::AdJacobianParams noAd{ nullptr, nullptr, 0u };
		unit->notifyDiscontinuousSectionTransition(1.0, 0u, { y.data(), yDot.data() }, noAd);
		unit->residual(cadet::SimulationTime{ 1.0, 0u }, cadet::ConstSimulationState{ y.data(), yDot.data() }, res.data(), tls);

		mb->destroyUnitOperation(unit);
		cadet::destroyModelBuilder(mb);
		return res;
	}

	/**
	 * @brief Checks that the externally dependent variants of all reactions of a unit give the residual of the plain reactions
	 * @param [in] jpp Unit operation with plain reactions
	 * @param [in] nGroups Number of reaction groups the unit has
	 */
	void testExternalReactionsInUnit(cadet::JsonParameterProvider& jpp, int nGroups)
	{
		// Copies are made from the JSON data, since the copy constructor keeps the original's open scopes
		cadet::JsonParameterProvider extJpp(*jpp.data());
		REQUIRE(makeReactionsExternal(extJpp, 0.5, 0.5) == nGroups);

		// Without the external function's half the rates halve. This has to change the residual, otherwise
		// the comparison below would hold whether or not the reactions and the external function act.
		cadet::JsonParameterProvider halfJpp(*jpp.data());
		REQUIRE(makeReactionsExternal(halfJpp, 0.5, 0.0) == nGroups);

		const std::vector<double> resPlain = unitResidualAtTimeOne(cadet::JsonParameterProvider(*jpp.data()));
		const std::vector<double> resExt = unitResidualAtTimeOne(extJpp);
		const std::vector<double> resHalf = unitResidualAtTimeOne(halfJpp);

		double maxChange = 0.0;
		for (std::size_t i = 0; i < resPlain.size(); ++i)
			maxChange = std::max(maxChange, std::abs(resHalf[i] - resPlain[i]));
		REQUIRE(maxChange > 1e-6);

		REQUIRE(resPlain.size() == resExt.size());
		for (std::size_t i = 0; i < resPlain.size(); ++i)
		{
			CAPTURE(i);
			CHECK(resExt[i] == cadet::test::makeApprox(resPlain[i], 1e-12, 1e-14));
		}
	}
}

TEST_CASE("External-function reactions in unit operations give the residual of the plain reactions", "[ReactionModel],[ExternalFunction],[UnitOp],[CI]")
{
	// Bulk reactions everywhere, and cross phase, solid and pore liquid reactions in two particle types where the unit has them
	const double volFrac[] = { 0.6, 0.4 };

	SECTION("CSTR")
	{
		cadet::JsonParameterProvider jpp = createCSTR(2);
		cadet::test::addBoundStates(jpp, { 1, 1 }, 1.0);
		cadet::test::addLinearBindingModel(jpp, true, { 0.1, 0.2 }, { 1.0, 0.9 });
		cadet::test::setInitialConditions(jpp, { 1.0, 2.0 }, { 3.0, 4.0 }, 6.0);
		cadet::test::setFlowRates(jpp, 0, 1.0);
		cadet::test::particle::extendModelToManyParticleTypes(jpp, 2, volFrac);
		cadet::test::reaction::extendModelWithDynamicReactions(jpp, 0, true, true, false);
		testExternalReactionsInUnit(jpp, 1 + 2 * 2);
	}

	for (const char* method : { "FV", "DG" })
	{
		SECTION(std::string("LRM ") + method)
		{
			cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("LUMPED_RATE_MODEL_WITHOUT_PORES", method);
			cadet::test::reaction::extendModelWithDynamicReactions(jpp, 0, true, false, false);
			testExternalReactionsInUnit(jpp, 1);
		}
	}

	const std::vector<std::pair<std::string, std::string>> columns = {
		{ "LUMPED_RATE_MODEL_WITH_PORES", "FV" },
		{ "GENERAL_RATE_MODEL", "FV" },
		{ "COLUMN_MODEL_1D_LRMP", "DG" },
		{ "COLUMN_MODEL_1D_GRM", "DG" },
		{ "COLUMN_MODEL_2D_GRM", "DG" },
		{ "COLUMN_MODEL_2D_GRM", "FV" },
	};
	for (const auto& column : columns)
	{
		SECTION(column.first + " " + column.second)
		{
			cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding(column.first, column.second);
			cadet::test::particle::extendModelToManyParticleTypes(jpp, 2, volFrac);
			cadet::test::reaction::extendModelWithDynamicReactions(jpp, 0, true, true, false);
			addParticleLiquidReactions(jpp);
			testExternalReactionsInUnit(jpp, 1 + 2 * 3);
		}
	}

	SECTION("MCT")
	{
		cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("MULTI_CHANNEL_TRANSPORT", "FV");
		cadet::test::reaction::extendModelWithDynamicReactions(jpp, 0, true, false, false);
		testExternalReactionsInUnit(jpp, 1);
	}
}
