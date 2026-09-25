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
#include "cadet/cadet.hpp"

#include "Approx.hpp"

#define CADET_LOGGING_DISABLE
#include "Logging.hpp"

#include "ColumnTests.hpp"
#include "JsonTestModels.hpp"
#include "ModelBuilderImpl.hpp"
#include "ParallelSupport.hpp"
#include "ReactionModelFactory.hpp"
#include "SimHelper.hpp"
#include "SimulationTypes.hpp"
#include "UnitOperationTests.hpp"
#include "Utils.hpp"
#include "common/Driver.hpp"
#include "linalg/Norms.hpp"
#include "model/StirredTankModel.hpp"
#include "model/UnitOperation.hpp"
#include "model/reaction/ConservedMoieties.hpp"

#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

cadet::model::CSTRModel* createAndConfigureCSTR(cadet::IModelBuilder& mb, cadet::JsonParameterProvider& jpp);
void checkJacobianAD(double flowRateIn, double flowRateOut, double flowRateFilter, std::function<void(cadet::JsonParameterProvider&, unsigned int)> modelRefiner);

namespace cadet
{
namespace test
{
namespace reaction
{
	Eigen::MatrixXd conservedMoietiesFromMassActionLaw(unsigned int nComp, const char* config)
	{
		std::vector<unsigned int> nBound(nComp, 0);
		std::vector<unsigned int> boundOffset(nComp, 0);

		cadet::JsonParameterProvider jpp(config);
		cadet::ReactionModelFactory rmf;
		std::unique_ptr<cadet::model::IDynamicReactionModel> model(rmf.createDynamic("MASS_ACTION_LAW"));
		REQUIRE(model);
		REQUIRE(model->configureModelDiscretization(jpp, nComp, nBound.data(), boundOffset.data()));
		REQUIRE(model->configure(jpp, 0, 0));

		const unsigned int nReactions = model->numReactions();
		Eigen::MatrixXd stoichiometry(nComp, nReactions);
		std::vector<bool> eqReactionFlags(nReactions);

		for (unsigned int r = 0; r < nReactions; ++r)
		{
			eqReactionFlags[r] = model->isInEquilibrium(r);
			for (unsigned int c = 0; c < nComp; ++c)
				stoichiometry(c, r) = model->getStoichiometry(c, r);
		}

		cadet::model::ConservedMoieties cm;
		REQUIRE(cm.configure(nComp, eqReactionFlags, stoichiometry, 1e-14));

		return cm.conservedMoietyMatrix();
	}

	void checkNullSpace(const Eigen::MatrixXd& left, const Eigen::MatrixXd& right)
	{
		const double nullspaceTol = 1e-12;
		const Eigen::MatrixXd product = left * right;

		for (Eigen::Index r = 0; r < product.rows(); ++r)
		{
			for (Eigen::Index c = 0; c < product.cols(); ++c)
				CHECK(std::abs(product(r, c)) <= nullspaceTol);
		}
	}
} // namespace reaction
} // namespace test
} // namespace cadet

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

TEST_CASE("CSTR liquid equilibrium MAL consistent initialization", "[CSTR],[MassActionLaw],[ReactionModel],[ConsistentInit],[CI],[ConservedMoieties]")
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	for (int adMode = 0; adMode < 2; ++adMode)
	{
		const bool adEnabled = (adMode > 0);
		SECTION(std::string("AD ") + (adEnabled ? "enabled" : "disabled"))
		{
			cadet::JsonParameterProvider jpp = createCSTR(2);
	
			jpp.set("NREAC_LIQUID", 1);
			jpp.addScope("liquid_reaction_000");
			jpp.pushScope("liquid_reaction_000");
			jpp.set("TYPE", "MASS_ACTION_LAW");
			jpp.set("MAL_KFWD", std::vector<double>{2.0});
			jpp.set("MAL_KBWD", std::vector<double>{1.0});
			jpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
			jpp.set("MAL_IS_KINETIC", std::vector<int>{false ? 1 : 0});
			jpp.popScope();

			cadet::model::CSTRModel* const cstr = createAndConfigureCSTR(*mb, jpp);
			cstr->setFlowRates(0.0, 0.0);
			
			std::vector<double> y(cstr->numDofs(), 0.0);
			y[2] = 2.7;
			y[3] = 0.3;
			y[4] = 1.0;

			cadet::test::unitoperation::testConsistentInitialization(cstr, adEnabled, y.data(), 1e-14, 1e-12);

			CHECK(y[2] == cadet::test::makeApprox(1.0, 1e-12, 1e-12));
			CHECK(y[3] == cadet::test::makeApprox(2.0, 1e-12, 1e-12));
			CHECK((y[2] + y[3]) == cadet::test::makeApprox(3.0, 1e-12, 1e-12));

			mb->destroyUnitOperation(cstr);
		}
	}

	destroyModelBuilder(mb);
}

TEST_CASE("CSTR liquid equilibrium MAL Jacobian vs AD", "[CSTR],[MassActionLaw],[ReactionModel],[Jacobian],[AD],[CI],[ConservedMoieties]")
{
	checkJacobianAD(1.0, 1.0, 0.0, [](cadet::JsonParameterProvider& jpp, unsigned int nComp) {
		
		jpp.set("NREAC_LIQUID", 1);
		jpp.addScope("liquid_reaction_000");
		jpp.pushScope("liquid_reaction_000");
		jpp.set("TYPE", "MASS_ACTION_LAW");
		jpp.set("MAL_KFWD", std::vector<double>{2.0});
		jpp.set("MAL_KBWD", std::vector<double>{1.0});
		jpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
		jpp.set("MAL_IS_KINETIC", std::vector<int>{false ? 1 : 0});
		jpp.popScope();

	});
}

TEST_CASE("CSTR kinetic and equilibrium MAL share liquid equilibrium", "[CSTR],[MassActionLaw],[ReactionModel],[ConsistentInit],[CI],[ConservedMoieties]")
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	cadet::JsonParameterProvider eqJpp = createCSTR(2);
	
	eqJpp.set("NREAC_LIQUID", 1);
	eqJpp.addScope("liquid_reaction_000");
	eqJpp.pushScope("liquid_reaction_000");
	eqJpp.set("TYPE", "MASS_ACTION_LAW");
	eqJpp.set("MAL_KFWD", std::vector<double>{2.0});
	eqJpp.set("MAL_KBWD", std::vector<double>{1.0});
	eqJpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
	eqJpp.set("MAL_IS_KINETIC", std::vector<int>{false ? 1 : 0});
	eqJpp.popScope();

	cadet::model::CSTRModel* const eqCstr = createAndConfigureCSTR(*mb, eqJpp);
	eqCstr->setFlowRates(0.0, 0.0);
	
	std::vector<double> yEq(eqCstr->numDofs(), 0.0);
	yEq[2] = 2.7;
	yEq[3] = 0.3;
	yEq[4] = 1.0;
	cadet::test::unitoperation::testConsistentInitialization(eqCstr, false, yEq.data(), 1e-14, 1e-12);
	
	cadet::JsonParameterProvider kinjpp = createCSTR(2);
	
	kinjpp.set("NREAC_LIQUID", 1);
	kinjpp.addScope("liquid_reaction_000");
	kinjpp.pushScope("liquid_reaction_000");
	kinjpp.set("TYPE", "MASS_ACTION_LAW");
	kinjpp.set("MAL_KFWD", std::vector<double>{2.0});
	kinjpp.set("MAL_KBWD", std::vector<double>{1.0});
	kinjpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
	kinjpp.set("MAL_IS_KINETIC", std::vector<int>{true ? 1 : 0});
	kinjpp.popScope();

	cadet::model::CSTRModel* const kinCstr = createAndConfigureCSTR(*mb, kinjpp);
	kinCstr->setFlowRates(0.0, 0.0);

	// Use the equilibrium model's consistently initialized state
	std::vector<double> yKin = yEq;
	yKin[2] = 1.0;
	yKin[3] = 2.0;
	
	std::vector<double> res(kinCstr->numDofs(), 0.0);
	cadet::util::ThreadLocalStorage tls;
	tls.resize(kinCstr->threadLocalMemorySize());

	kinCstr->notifyDiscontinuousSectionTransition(0.0, 0u, {yKin.data(), nullptr}, cadet::AdJacobianParams{nullptr, nullptr, 0u});
	kinCstr->residual(cadet::SimulationTime{0.0, 0u}, cadet::ConstSimulationState{yKin.data(), nullptr}, res.data(), tls);

	CHECK(yKin[2] == cadet::test::makeApprox(1.0, 1e-12, 1e-12));
	CHECK(yKin[3] == cadet::test::makeApprox(2.0, 1e-12, 1e-12));
	
	for (double val : res)
		CHECK(std::abs(val) <= 1e-12);
	
	mb->destroyUnitOperation(eqCstr);
	mb->destroyUnitOperation(kinCstr);
	destroyModelBuilder(mb);
}

TEST_CASE("1D column liquid equilibrium MAL consistent initialization", "[Column_1D],[MassActionLaw],[ReactionModel],[ConsistentInit],[CI],[ConservedMoieties]")
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	for (int adMode = 0; adMode < 2; ++adMode)
	{
		const bool adEnabled = (adMode > 0);
		SECTION(std::string("AD ") + (adEnabled ? "enabled" : "disabled"))
		{
			const bool kinetic = false;
			const bool withParticles = false;
			cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_LRMP", "DG");
			if (!withParticles)
			{
				// Keep the initialization tests focused on the liquid bulk phase
				jpp.set("NPARTYPE", 0);
				jpp.remove("particle_type_000");
			}
			jpp.set("INIT_C", std::vector<double>{2.7, 0.3});

			jpp.set("NREAC_LIQUID", 1);
			jpp.addScope("liquid_reaction_000");
			jpp.pushScope("liquid_reaction_000");
			jpp.set("TYPE", "MASS_ACTION_LAW");
			jpp.set("MAL_KFWD", std::vector<double>{2.0});
			jpp.set("MAL_KBWD", std::vector<double>{1.0});
			jpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
			jpp.set("MAL_IS_KINETIC", std::vector<int>{kinetic ? 1 : 0});
			jpp.popScope();

			cadet::IUnitOperation* const unit =cadet::test::unitoperation::createAndConfigureUnit(jpp, *mb);

			std::vector<double> y(unit->numDofs(), 0.0);
			y[0] = 1.0;
			y[1] = 2.0;
			// set column bulk state
			const unsigned int nComp = unit->numComponents();
			const unsigned int bulkOffset = nComp * unit->numInletPorts();
			const unsigned int inletDofs = unit->numComponents() * unit->numInletPorts();
			const unsigned columnBulkPoints =  (unit->numDofs() - inletDofs) / nComp;
			
			for (unsigned int point = 0; point < columnBulkPoints; ++point)
			{
				y[bulkOffset + point * nComp] = 2.7;
				y[bulkOffset + point * nComp + 1] = 0.3;
			}

			cadet::test::unitoperation::testConsistentInitialization(
				unit, adEnabled, y.data(), 1e-14, 1e-11);
			
			for (unsigned int point = 0; point < columnBulkPoints; ++point)
			{
				const double* const c = y.data() + bulkOffset + point * nComp;
				CHECK(c[0] == cadet::test::makeApprox(1.0, 1e-12, 1e-12));
				CHECK(c[1] == cadet::test::makeApprox(2.0, 1e-12, 1e-12));
				CHECK((c[0] + c[1]) == cadet::test::makeApprox(3.0, 1e-12, 1e-12));
			}

			mb->destroyUnitOperation(unit);
		}
	}

	destroyModelBuilder(mb);
}

TEST_CASE("1D column particle liquid and binding equilibrium consistent initialization", "[Column_1D],[MassActionLaw],[ReactionModel],[ConsistentInit],[Binding],[Particle],[CI],[ConservedMoieties]")
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_LRMP", "DG");
	jpp.pushScope("particle_type_000");
	jpp.set("INIT_CP", std::vector<double>{2.7, 0.3});
	jpp.set("INIT_CS", std::vector<double>{0.4, 0.6});
	jpp.set("PORE_ACCESSIBILITY", std::vector<double>{0.5, 0.8});
	jpp.pushScope("adsorption");
	jpp.set("IS_KINETIC", 0);
	jpp.popScope();
	jpp.set("NREAC_LIQUID", 1);
	jpp.addScope("liquid_reaction_000");
	jpp.pushScope("liquid_reaction_000");
	jpp.set("TYPE", "MASS_ACTION_LAW");
	jpp.set("MAL_KFWD", std::vector<double>{2.0});
	jpp.set("MAL_KBWD", std::vector<double>{1.0});
	jpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
	jpp.set("MAL_IS_KINETIC", std::vector<int>{0});
	jpp.popScope();
	jpp.popScope();

	cadet::IUnitOperation* const unit = cadet::test::unitoperation::createAndConfigureUnit(jpp, *mb);
	std::vector<double> y(unit->numDofs(), 0.0);
	std::vector<double> yDot(unit->numDofs(), 0.0);
	unit->readInitialCondition(jpp);
	unit->applyInitialCondition(cadet::SimulationState{y.data(), yDot.data()});

	const unsigned int nComp = unit->numComponents();
	const unsigned int inletDofs = nComp * unit->numInletPorts();
	const unsigned int nBound = 2;
	const unsigned int strideParticle = nComp + nBound;
	const unsigned int dofsPerPoint = nComp + strideParticle;
	REQUIRE((unit->numDofs() - inletDofs) % dofsPerPoint == 0);
	const unsigned int nPoints = (unit->numDofs() - inletDofs) / dofsPerPoint;
	const unsigned int particleOffset = inletDofs + nPoints * nComp;
	const double invBetaP0 = 0.25 / (0.5 * 0.75);
	const double invBetaP1 = 0.25 / (0.8 * 0.75);
	const double initialConservedMass = 2.7 + 0.3 + invBetaP0 * 0.4 + invBetaP1 * 0.6;

	cadet::test::unitoperation::testConsistentInitialization(unit, false, y.data(), 1e-14, 1e-10);

	for (unsigned int point = 0; point < nPoints; ++point)
	{
		double const* const cp = y.data() + particleOffset + point * strideParticle;
		double const* const q = cp + nComp;
		CHECK(cp[1] == cadet::test::makeApprox(2.0 * cp[0], 1e-11, 1e-11));
		CHECK(q[0] == cadet::test::makeApprox((12.3 / 45.0) * cp[0], 1e-11, 1e-11));
		CHECK(q[1] == cadet::test::makeApprox((35.5 / 20.0) * cp[1], 1e-11, 1e-11));
		CHECK((cp[0] + cp[1] + invBetaP0 * q[0] + invBetaP1 * q[1])
			== cadet::test::makeApprox(initialConservedMass, 1e-11, 1e-11));
	}

	mb->destroyUnitOperation(unit);
	destroyModelBuilder(mb);
}

TEST_CASE("1D column general rate particle liquid and binding equilibrium consistent initialization", "[Column_1D],[MassActionLaw],[ReactionModel],[ConsistentInit],[Binding],[Particle],[GRM],[CI],[ConservedMoieties]")
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_GRM", "DG");
	jpp.pushScope("particle_type_000");
	jpp.set("INIT_CP", std::vector<double>{2.7, 0.3});
	jpp.set("INIT_CS", std::vector<double>{0.4, 0.6});
	jpp.set("PORE_ACCESSIBILITY", std::vector<double>{0.5, 0.8});
	jpp.pushScope("adsorption");
	jpp.set("IS_KINETIC", 0);
	jpp.popScope();
	jpp.set("NREAC_LIQUID", 1);
	jpp.addScope("liquid_reaction_000");
	jpp.pushScope("liquid_reaction_000");
	jpp.set("TYPE", "MASS_ACTION_LAW");
	jpp.set("MAL_KFWD", std::vector<double>{2.0});
	jpp.set("MAL_KBWD", std::vector<double>{1.0});
	jpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
	jpp.set("MAL_IS_KINETIC", std::vector<int>{0});
	jpp.popScope();
	jpp.popScope();

	cadet::IUnitOperation* const unit = cadet::test::unitoperation::createAndConfigureUnit(jpp, *mb);
	std::vector<double> y(unit->numDofs(), 0.0);
	std::vector<double> yDot(unit->numDofs(), 0.0);
	unit->readInitialCondition(jpp);
	unit->applyInitialCondition(cadet::SimulationState{y.data(), yDot.data()});

	const unsigned int nComp = unit->numComponents();
	const unsigned int inletDofs = nComp * unit->numInletPorts();
	const unsigned int nBound = 2;
	const unsigned int nParticlePoints = 4;
	const unsigned int strideParticlePoint = nComp + nBound;
	const unsigned int dofsPerColumnPoint = nComp + nParticlePoints * strideParticlePoint;
	REQUIRE((unit->numDofs() - inletDofs) % dofsPerColumnPoint == 0);
	const unsigned int nColumnPoints = (unit->numDofs() - inletDofs) / dofsPerColumnPoint;
	const unsigned int particleOffset = inletDofs + nColumnPoints * nComp;
	const double invBetaP0 = 0.25 / (0.5 * 0.75);
	const double invBetaP1 = 0.25 / (0.8 * 0.75);
	const double initialConservedMass = 2.7 + 0.3 + invBetaP0 * 0.4 + invBetaP1 * 0.6;

	cadet::test::unitoperation::testConsistentInitialization(unit, false, y.data(), 1e-14, 1e-10);

	for (unsigned int columnPoint = 0; columnPoint < nColumnPoints; ++columnPoint)
	{
		for (unsigned int particlePoint = 0; particlePoint < nParticlePoints; ++particlePoint)
		{
			double const* const cp = y.data() + particleOffset
				+ (columnPoint * nParticlePoints + particlePoint) * strideParticlePoint;
			double const* const q = cp + nComp;
			CHECK(cp[1] == cadet::test::makeApprox(2.0 * cp[0], 1e-11, 1e-11));
			CHECK(q[0] == cadet::test::makeApprox((12.3 / 45.0) * cp[0], 1e-11, 1e-11));
			CHECK(q[1] == cadet::test::makeApprox((35.5 / 20.0) * cp[1], 1e-11, 1e-11));
			CHECK((cp[0] + cp[1] + invBetaP0 * q[0] + invBetaP1 * q[1])
				== cadet::test::makeApprox(initialConservedMass, 1e-11, 1e-11));
		}
	}

	mb->destroyUnitOperation(unit);
	destroyModelBuilder(mb);
}

TEST_CASE("1D column general rate particle liquid equilibrium MAL consistent initialization", "[Column_1D],[MassActionLaw],[ReactionModel],[ConsistentInit],[Particle],[GRM],[CI],[ConservedMoieties]")
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_GRM", "DG");
	jpp.pushScope("particle_type_000");
	jpp.set("INIT_CP", std::vector<double>{2.7, 0.3});
	jpp.set("NREAC_LIQUID", 1);
	jpp.addScope("liquid_reaction_000");
	jpp.pushScope("liquid_reaction_000");
	jpp.set("TYPE", "MASS_ACTION_LAW");
	jpp.set("MAL_KFWD", std::vector<double>{2.0});
	jpp.set("MAL_KBWD", std::vector<double>{1.0});
	jpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
	jpp.set("MAL_IS_KINETIC", std::vector<int>{0});
	jpp.popScope();
	jpp.popScope();

	cadet::IUnitOperation* const unit = cadet::test::unitoperation::createAndConfigureUnit(jpp, *mb);
	std::vector<double> y(unit->numDofs(), 0.0);
	std::vector<double> yDot(unit->numDofs(), 0.0);
	unit->readInitialCondition(jpp);
	unit->applyInitialCondition(cadet::SimulationState{y.data(), yDot.data()});

	const unsigned int nComp = unit->numComponents();
	const unsigned int inletDofs = nComp * unit->numInletPorts();
	const unsigned int nBound = 2;
	const unsigned int nParticlePoints = 4;
	const unsigned int strideParticlePoint = nComp + nBound;
	const unsigned int dofsPerColumnPoint = nComp + nParticlePoints * strideParticlePoint;
	REQUIRE((unit->numDofs() - inletDofs) % dofsPerColumnPoint == 0);
	const unsigned int nColumnPoints = (unit->numDofs() - inletDofs) / dofsPerColumnPoint;
	const unsigned int particleOffset = inletDofs + nColumnPoints * nComp;

	cadet::test::unitoperation::testConsistentInitialization(unit, false, y.data(), 1e-14, 1e-10);

	for (unsigned int columnPoint = 0; columnPoint < nColumnPoints; ++columnPoint)
	{
		for (unsigned int particlePoint = 0; particlePoint < nParticlePoints; ++particlePoint)
		{
			double const* const cp = y.data() + particleOffset
				+ (columnPoint * nParticlePoints + particlePoint) * strideParticlePoint;
			CHECK(cp[0] == cadet::test::makeApprox(1.0, 1e-11, 1e-11));
			CHECK(cp[1] == cadet::test::makeApprox(2.0, 1e-11, 1e-11));
			CHECK((cp[0] + cp[1]) == cadet::test::makeApprox(3.0, 1e-11, 1e-11));
		}
	}

	mb->destroyUnitOperation(unit);
	destroyModelBuilder(mb);
}

TEST_CASE("1D column particle liquid equilibrium MAL Jacobian vs AD", "[Column_1D],[MassActionLaw],[ReactionModel],[Jacobian],[AD],[Particle],[CI],[ConservedMoieties]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_LRMP", "DG");
	jpp.pushScope("particle_type_000");
	jpp.set("NREAC_LIQUID", 1);
	jpp.addScope("liquid_reaction_000");
	jpp.pushScope("liquid_reaction_000");
	jpp.set("TYPE", "MASS_ACTION_LAW");
	jpp.set("MAL_KFWD", std::vector<double>{2.0});
	jpp.set("MAL_KBWD", std::vector<double>{1.0});
	jpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
	jpp.set("MAL_IS_KINETIC", std::vector<int>{0});
	jpp.popScope();
	jpp.popScope();

	cadet::test::column::testJacobianAD(
		jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("1D column general rate particle liquid equilibrium MAL Jacobian vs AD", "[Column_1D],[MassActionLaw],[ReactionModel],[Jacobian],[AD],[Particle],[GRM],[CI],[ConservedMoieties]")
{
	for (const std::string& spatialMethod : {std::string("DG"), std::string("DGFV")})
	{
		SECTION(spatialMethod)
		{
			cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_GRM", spatialMethod);
			jpp.pushScope("particle_type_000");
			jpp.set("NREAC_LIQUID", 1);
			jpp.addScope("liquid_reaction_000");
			jpp.pushScope("liquid_reaction_000");
			jpp.set("TYPE", "MASS_ACTION_LAW");
			jpp.set("MAL_KFWD", std::vector<double>{2.0});
			jpp.set("MAL_KBWD", std::vector<double>{1.0});
			jpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
			jpp.set("MAL_IS_KINETIC", std::vector<int>{0});
			jpp.popScope();
			jpp.popScope();

			cadet::test::column::testJacobianAD(
				jpp, std::numeric_limits<float>::epsilon() * 100.0);
		}
	}
}

TEST_CASE("1D column general rate particle liquid equilibrium MAL time derivative Jacobian vs FD", "[Column_1D],[MassActionLaw],[ReactionModel],[Jacobian],[FD],[Particle],[GRM],[CI],[ConservedMoieties]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_GRM", "DG");
	jpp.pushScope("particle_type_000");
	jpp.set("NREAC_LIQUID", 1);
	jpp.addScope("liquid_reaction_000");
	jpp.pushScope("liquid_reaction_000");
	jpp.set("TYPE", "MASS_ACTION_LAW");
	jpp.set("MAL_KFWD", std::vector<double>{2.0});
	jpp.set("MAL_KBWD", std::vector<double>{1.0});
	jpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
	jpp.set("MAL_IS_KINETIC", std::vector<int>{0});
	jpp.popScope();
	jpp.popScope();

	cadet::test::unitoperation::testTimeDerivativeJacobianFD(jpp, 1e-6, 0.0, 9e-4);
}

TEST_CASE("1D column liquid equilibrium MAL initial concentration sensitivities", "[Column_1D],[MassActionLaw],[ReactionModel],[Sensitivity],[Particle],[GRM],[ConsistentInit],[CI],[ConservedMoieties]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_GRM", "DG");
	jpp.set("INIT_C", std::vector<double>{2.7, 0.3});
	jpp.set("NREAC_LIQUID", 1);
	jpp.addScope("liquid_reaction_000");
	jpp.pushScope("liquid_reaction_000");
	jpp.set("TYPE", "MASS_ACTION_LAW");
	jpp.set("MAL_KFWD", std::vector<double>{2.0});
	jpp.set("MAL_KBWD", std::vector<double>{1.0});
	jpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
	jpp.set("MAL_IS_KINETIC", std::vector<int>{0});
	jpp.popScope();
	jpp.pushScope("particle_type_000");
	jpp.set("INIT_CP", std::vector<double>{2.7, 0.3});
	jpp.set("NREAC_LIQUID", 1);
	jpp.addScope("liquid_reaction_000");
	jpp.pushScope("liquid_reaction_000");
	jpp.set("TYPE", "MASS_ACTION_LAW");
	jpp.set("MAL_KFWD", std::vector<double>{2.0});
	jpp.set("MAL_KBWD", std::vector<double>{1.0});
	jpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
	jpp.set("MAL_IS_KINETIC", std::vector<int>{0});
	jpp.popScope();
	jpp.popScope();

	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);
	cadet::ad::setDirections(cadet::ad::getMaxDirections());
	cadet::IUnitOperation* const unit = cadet::test::unitoperation::createAndConfigureUnit(jpp, *mb);
	REQUIRE(unit->setSensitiveParameter(cadet::makeParamId("INIT_C", 0, 0, cadet::ParTypeIndep,
		cadet::BoundStateIndep, cadet::ReactionIndep, cadet::SectionIndep), 0, 1.0));
	REQUIRE(unit->setSensitiveParameter(cadet::makeParamId("INIT_CP", 0, 0, 0,
		cadet::BoundStateIndep, cadet::ReactionIndep, cadet::SectionIndep), 1, 1.0));
	REQUIRE(unit->numSensParams() == 2);

	const unsigned int nDofs = unit->numDofs();
	std::vector<double> y(nDofs, 0.0);
	std::vector<double> yDot(nDofs, 0.0);
	std::vector<cadet::active> adRes(nDofs);
	const cadet::AdJacobianParams adParams{adRes.data(), nullptr, 2u};
	unit->prepareADvectors(adParams);
	unit->readInitialCondition(jpp);
	unit->applyInitialCondition(cadet::SimulationState{y.data(), yDot.data()});
	cadet::util::ThreadLocalStorage tls;
	tls.resize(unit->threadLocalMemorySize());
	const cadet::SimulationTime simTime{0.0, 0u};
	unit->notifyDiscontinuousSectionTransition(0.0, 0u, {y.data(), yDot.data()}, adParams);
	unit->consistentInitialState(simTime, y.data(), adParams, 1e-14, tls);
	unit->residualWithJacobian(simTime, {y.data(), nullptr}, yDot.data(), adParams, tls);
	unit->consistentInitialTimeDerivative(simTime, y.data(), yDot.data(), tls);
	unit->residualSensFwdWithJacobian(simTime, {y.data(), yDot.data()}, adParams, tls);

	std::vector<std::vector<double>> sensY(2, std::vector<double>(nDofs, 0.0));
	std::vector<std::vector<double>> sensYdot(2, std::vector<double>(nDofs, 0.0));
	std::vector<double*> sensYPtr{sensY[0].data(), sensY[1].data()};
	std::vector<double*> sensYdotPtr{sensYdot[0].data(), sensYdot[1].data()};
	unit->initializeSensitivityStates(sensYPtr);
	unit->consistentInitialSensitivity(simTime, {y.data(), yDot.data()}, sensYPtr, sensYdotPtr, adRes.data(), tls);

	const unsigned int nComp = unit->numComponents();
	const unsigned int inletDofs = nComp * unit->numInletPorts();
	const unsigned int nParticlePoints = 4;
	const unsigned int strideParticlePoint = 4;
	const unsigned int dofsPerColumnPoint = nComp + nParticlePoints * strideParticlePoint;
	const unsigned int nColumnPoints = (nDofs - inletDofs) / dofsPerColumnPoint;
	const unsigned int particleOffset = inletDofs + nColumnPoints * nComp;
	for (unsigned int point = 0; point < nColumnPoints; ++point)
	{
		CHECK(sensY[0][inletDofs + point * nComp] == cadet::test::makeApprox(1.0 / 3.0, 1e-11, 1e-11));
		CHECK(sensY[0][inletDofs + point * nComp + 1] == cadet::test::makeApprox(2.0 / 3.0, 1e-11, 1e-11));
		for (unsigned int shell = 0; shell < nParticlePoints; ++shell)
		{
			const unsigned int shellOffset = particleOffset + (point * nParticlePoints + shell) * strideParticlePoint;
			CHECK(sensY[1][shellOffset] == cadet::test::makeApprox(1.0 / 3.0, 1e-11, 1e-11));
			CHECK(sensY[1][shellOffset + 1] == cadet::test::makeApprox(2.0 / 3.0, 1e-11, 1e-11));
		}
	}

	std::vector<std::vector<double>> sensRes(2, std::vector<double>(nDofs, 0.0));
	std::vector<const double*> sensYConst{sensY[0].data(), sensY[1].data()};
	std::vector<const double*> sensYdotConst{sensYdot[0].data(), sensYdot[1].data()};
	std::vector<double*> sensResPtr{sensRes[0].data(), sensRes[1].data()};
	std::vector<double> tmp1(nDofs, 0.0);
	std::vector<double> tmp2(nDofs, 0.0);
	std::vector<double> tmp3(nDofs, 0.0);
	unit->residualSensFwdAdOnly(simTime, {y.data(), yDot.data()}, adRes.data(), tls);
	unit->residualSensFwdCombine(simTime, {y.data(), yDot.data()}, sensYConst, sensYdotConst,
		sensResPtr, adRes.data(), tmp1.data(), tmp2.data(), tmp3.data());
	for (unsigned int param = 0; param < 2; ++param)
		CHECK(cadet::linalg::linfNorm(sensRes[param].data() + inletDofs, nDofs - inletDofs) <= 1e-10);

	mb->destroyUnitOperation(unit);
	destroyModelBuilder(mb);
}

TEST_CASE("1D column coupled particle liquid and binding equilibrium initial sensitivities", "[Column_1D],[MassActionLaw],[ReactionModel],[Sensitivity],[Particle],[GRM],[ConsistentInit],[CI],[ConservedMoieties]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_GRM", "DG");
	cadet::test::setBindingMode(jpp, false);
	jpp.pushScope("particle_type_000");
	jpp.set("INIT_CP", std::vector<double>{2.7, 0.3});
	jpp.set("NREAC_LIQUID", 1);
	jpp.addScope("liquid_reaction_000");
	jpp.pushScope("liquid_reaction_000");
	jpp.set("TYPE", "MASS_ACTION_LAW");
	jpp.set("MAL_KFWD", std::vector<double>{2.0});
	jpp.set("MAL_KBWD", std::vector<double>{1.0});
	jpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
	jpp.set("MAL_IS_KINETIC", std::vector<int>{0});
	jpp.popScope();
	jpp.popScope();

	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);
	cadet::ad::setDirections(cadet::ad::getMaxDirections());
	cadet::IUnitOperation* const unit = cadet::test::unitoperation::createAndConfigureUnit(jpp, *mb);
	REQUIRE(unit->setSensitiveParameter(cadet::makeParamId("INIT_CP", 0, 0, 0,
		cadet::BoundStateIndep, cadet::ReactionIndep, cadet::SectionIndep), 0, 1.0));

	const unsigned int nDofs = unit->numDofs();
	std::vector<double> y(nDofs, 0.0);
	std::vector<double> yDot(nDofs, 0.0);
	std::vector<cadet::active> adRes(nDofs);
	const cadet::AdJacobianParams adParams{adRes.data(), nullptr, 1u};
	unit->prepareADvectors(adParams);
	unit->readInitialCondition(jpp);
	unit->applyInitialCondition(cadet::SimulationState{y.data(), yDot.data()});
	cadet::util::ThreadLocalStorage tls;
	tls.resize(unit->threadLocalMemorySize());
	const cadet::SimulationTime simTime{0.0, 0u};
	unit->notifyDiscontinuousSectionTransition(0.0, 0u, {y.data(), yDot.data()}, adParams);
	unit->consistentInitialState(simTime, y.data(), adParams, 1e-14, tls);
	unit->residualWithJacobian(simTime, {y.data(), nullptr}, yDot.data(), adParams, tls);
	unit->consistentInitialTimeDerivative(simTime, y.data(), yDot.data(), tls);
	unit->residualSensFwdWithJacobian(simTime, {y.data(), yDot.data()}, adParams, tls);

	std::vector<double> sensY(nDofs, 0.0);
	std::vector<double> sensYdot(nDofs, 0.0);
	std::vector<double*> sensYPtr{sensY.data()};
	std::vector<double*> sensYdotPtr{sensYdot.data()};
	unit->initializeSensitivityStates(sensYPtr);
	unit->consistentInitialSensitivity(simTime, {y.data(), yDot.data()}, sensYPtr, sensYdotPtr, adRes.data(), tls);

	const unsigned int nComp = unit->numComponents();
	const unsigned int inletDofs = nComp * unit->numInletPorts();
	const unsigned int nParticlePoints = 4;
	const unsigned int strideParticlePoint = 4;
	const unsigned int dofsPerColumnPoint = nComp + nParticlePoints * strideParticlePoint;
	const unsigned int nColumnPoints = (nDofs - inletDofs) / dofsPerColumnPoint;
	const unsigned int particleOffset = inletDofs + nColumnPoints * nComp;
	const double invBetaP = (1.0 - 0.75) / 0.75;
	for (unsigned int point = 0; point < nColumnPoints; ++point)
	{
		for (unsigned int shell = 0; shell < nParticlePoints; ++shell)
		{
			const unsigned int shellOffset = particleOffset + (point * nParticlePoints + shell) * strideParticlePoint;
			CHECK(sensY[shellOffset + 1] == cadet::test::makeApprox(2.0 * sensY[shellOffset], 1e-11, 1e-11));
			CHECK(sensY[shellOffset + 2] == cadet::test::makeApprox(12.3 / 45.0 * sensY[shellOffset], 1e-11, 1e-11));
			CHECK(sensY[shellOffset + 3] == cadet::test::makeApprox(35.5 / 20.0 * sensY[shellOffset + 1], 1e-11, 1e-11));
			const double conservedSensitivity = sensY[shellOffset] + sensY[shellOffset + 1]
				+ invBetaP * (sensY[shellOffset + 2] + sensY[shellOffset + 3]);
			CHECK(conservedSensitivity == cadet::test::makeApprox(1.0, 1e-11, 1e-11));
		}
	}

	std::vector<double> sensRes(nDofs, 0.0);
	std::vector<const double*> sensYConst{sensY.data()};
	std::vector<const double*> sensYdotConst{sensYdot.data()};
	std::vector<double*> sensResPtr{sensRes.data()};
	std::vector<double> tmp1(nDofs, 0.0);
	std::vector<double> tmp2(nDofs, 0.0);
	std::vector<double> tmp3(nDofs, 0.0);
	unit->residualSensFwdAdOnly(simTime, {y.data(), yDot.data()}, adRes.data(), tls);
	unit->residualSensFwdCombine(simTime, {y.data(), yDot.data()}, sensYConst, sensYdotConst,
		sensResPtr, adRes.data(), tmp1.data(), tmp2.data(), tmp3.data());
	CHECK(cadet::linalg::linfNorm(sensRes.data() + inletDofs, nDofs - inletDofs) <= 1e-10);

	mb->destroyUnitOperation(unit);
	destroyModelBuilder(mb);
}

TEST_CASE("1D column liquid equilibrium MAL outlet sensitivity with respect to initial concentration", "[Column_1D],[MassActionLaw],[ReactionModel],[Sensitivity],[Particle],[GRM],[Simulation],[CI],[ConservedMoieties]")
{
	cadet::JsonParameterProvider unitJpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_GRM", "DG");
	unitJpp.set("INIT_C", std::vector<double>{2.7, 0.3});
	unitJpp.set("NREAC_LIQUID", 1);
	unitJpp.addScope("liquid_reaction_000");
	unitJpp.pushScope("liquid_reaction_000");
	unitJpp.set("TYPE", "MASS_ACTION_LAW");
	unitJpp.set("MAL_KFWD", std::vector<double>{2.0});
	unitJpp.set("MAL_KBWD", std::vector<double>{1.0});
	unitJpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
	unitJpp.set("MAL_IS_KINETIC", std::vector<int>{0});
	unitJpp.popScope();
	unitJpp.pushScope("particle_type_000");
	unitJpp.set("INIT_CP", std::vector<double>{2.7, 0.3});
	unitJpp.set("NREAC_LIQUID", 1);
	unitJpp.addScope("liquid_reaction_000");
	unitJpp.pushScope("liquid_reaction_000");
	unitJpp.set("TYPE", "MASS_ACTION_LAW");
	unitJpp.set("MAL_KFWD", std::vector<double>{2.0});
	unitJpp.set("MAL_KBWD", std::vector<double>{1.0});
	unitJpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
	unitJpp.set("MAL_IS_KINETIC", std::vector<int>{0});
	unitJpp.popScope();
	unitJpp.popScope();

	cadet::JsonParameterProvider baseJpp = createLWE("COLUMN_MODEL_1D_GRM", "DG");
	nlohmann::json& root = *baseJpp.data();
	root["model"]["unit_000"] = *unitJpp.data();
	nlohmann::json& inlet = root["model"]["unit_001"];
	inlet["NCOMP"] = 2;
	for (const char* section : {"sec_000", "sec_001", "sec_002"})
	{
		inlet[section]["CONST_COEFF"] = {1.0, 2.0};
		inlet[section]["LIN_COEFF"] = {0.0, 0.0};
		inlet[section]["QUAD_COEFF"] = {0.0, 0.0};
		inlet[section]["CUBE_COEFF"] = {0.0, 0.0};
	}

	cadet::JsonParameterProvider sensitivityJpp = baseJpp;
	const cadet::ParameterId initC = cadet::makeParamId("INIT_C", 0, 0, cadet::ParTypeIndep,cadet::BoundStateIndep, cadet::ReactionIndep, cadet::SectionIndep);
	cadet::test::addSensitivity(sensitivityJpp, "INIT_C", initC, 1e-10);
	cadet::test::returnSensitivities(sensitivityJpp, 0);
	cadet::Driver sensitivityDriver;
	sensitivityDriver.configure(sensitivityJpp);
	sensitivityDriver.run();

	const double step = 1e-5;
	cadet::Driver leftDriver;
	leftDriver.configure(baseJpp);
	const double initialValue = leftDriver.simulator()->model()->getParameterDouble(initC);
	REQUIRE(!std::isnan(initialValue));
	leftDriver.simulator()->setParameterValue(initC, initialValue - step);
	leftDriver.simulator()->applyInitialCondition();
	leftDriver.run();
	cadet::Driver rightDriver;
	rightDriver.configure(baseJpp);
	rightDriver.simulator()->setParameterValue(initC, initialValue + step);
	rightDriver.simulator()->applyInitialCondition();
	rightDriver.run();

	cadet::InternalStorageUnitOpRecorder const* const sensitivityData = sensitivityDriver.solution()->unitOperation(0);
	double const* sensitivity = sensitivityData->sensOutlet(0);
	double const* left = leftDriver.solution()->unitOperation(0)->outlet();
	double const* right = rightDriver.solution()->unitOperation(0)->outlet();
	const unsigned int numValues = sensitivityData->numDataPoints() * sensitivityData->numComponents()
		* sensitivityData->numOutletPorts();
	for (unsigned int value = 0; value < numValues; ++value)
	{
		const double finiteDifference = (right[value] - left[value]) / (2.0 * step);
		CHECK(sensitivity[value] == cadet::test::makeApprox(finiteDifference, 2e-4, 2e-7));
	}
}

TEST_CASE("1D liquid equilibrium MAL Jacobian vs AD", "[Column_1D],[MassActionLaw],[ReactionModel],[Jacobian],[AD],[CI],[ConservedMoieties]")
{
	SECTION("Without particles")
	{
		const bool kinetic = false;
		const bool withParticles = false;
		cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_LRMP", "DG");
		if (!withParticles)
		{
			// Keep the initialization tests focused on the liquid bulk phase
			jpp.set("NPARTYPE", 0);
			jpp.remove("particle_type_000");
		}
		jpp.set("INIT_C", std::vector<double>{2.7, 0.3});

		jpp.set("NREAC_LIQUID", 1);
		jpp.addScope("liquid_reaction_000");
		jpp.pushScope("liquid_reaction_000");
		jpp.set("TYPE", "MASS_ACTION_LAW");
		jpp.set("MAL_KFWD", std::vector<double>{2.0});
		jpp.set("MAL_KBWD", std::vector<double>{1.0});
		jpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
		jpp.set("MAL_IS_KINETIC", std::vector<int>{kinetic ? 1 : 0});
		jpp.popScope();

		cadet::test::column::testJacobianAD(
			jpp, std::numeric_limits<float>::epsilon() * 100.0);
	}

	SECTION("With bulk-to-particle film coupling")
	{
		const bool kinetic = false;
		const bool withParticles = true;
		cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_LRMP", "DG");
		if (!withParticles)
		{
			// Keep the initialization tests focused on the liquid bulk phase
			jpp.set("NPARTYPE", 0);
			jpp.remove("particle_type_000");
		}
		jpp.set("INIT_C", std::vector<double>{2.7, 0.3});

		jpp.set("NREAC_LIQUID", 1);
		jpp.addScope("liquid_reaction_000");
		jpp.pushScope("liquid_reaction_000");
		jpp.set("TYPE", "MASS_ACTION_LAW");
		jpp.set("MAL_KFWD", std::vector<double>{2.0});
		jpp.set("MAL_KBWD", std::vector<double>{1.0});
		jpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
		jpp.set("MAL_IS_KINETIC", std::vector<int>{kinetic ? 1 : 0});
		jpp.popScope();
		cadet::test::column::testJacobianAD(
			jpp, std::numeric_limits<float>::epsilon() * 100.0);
	}
}

TEST_CASE("1D kinetic and equilibrium MAL share liquid equilibrium", "[Column_1D],[MassActionLaw],[ReactionModel],[ConsistentInit],[CI],[ConservedMoieties]")
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	bool kinetic = false;
	bool withParticles = false;
	
	cadet::JsonParameterProvider eqJpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_LRMP", "DG");
	if (!withParticles)
	{
		// Keep the initialization tests focused on the liquid bulk phase
		eqJpp.set("NPARTYPE", 0);
		eqJpp.remove("particle_type_000");
	}
	eqJpp.set("INIT_C", std::vector<double>{2.7, 0.3});

	eqJpp.set("NREAC_LIQUID", 1);
	eqJpp.addScope("liquid_reaction_000");
	eqJpp.pushScope("liquid_reaction_000");
	eqJpp.set("TYPE", "MASS_ACTION_LAW");
	eqJpp.set("MAL_KFWD", std::vector<double>{2.0});
	eqJpp.set("MAL_KBWD", std::vector<double>{1.0});
	eqJpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
	eqJpp.set("MAL_IS_KINETIC", std::vector<int>{kinetic ? 1 : 0});
	eqJpp.popScope();

	cadet::IUnitOperation* const eqUnit =
		cadet::test::unitoperation::createAndConfigureUnit(eqJpp, *mb);

	std::vector<double> yEq(eqUnit->numDofs(), 0.0);

	// set column bulk state
	const unsigned int nComp = eqUnit->numComponents();
	const unsigned int bulkOffset = nComp * eqUnit->numInletPorts();
	const unsigned int inletDofs = eqUnit->numComponents() * eqUnit->numInletPorts();
	const unsigned columnBulkPoints =  (eqUnit->numDofs() - inletDofs) / nComp;
	for (unsigned int point = 0; point < columnBulkPoints; ++point)
	{
		yEq[bulkOffset + point * nComp] = 2.7;
		yEq[bulkOffset + point * nComp + 1] = 0.3;
	}
	cadet::test::unitoperation::testConsistentInitialization(eqUnit, false, yEq.data(), 1e-14, 1e-11);

	kinetic = true;
	withParticles = false;
	cadet::JsonParameterProvider kinJpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_LRMP", "DG");
	if (!withParticles)
	{
		// Keep the initialization tests focused on the liquid bulk phase
		kinJpp.set("NPARTYPE", 0);
		kinJpp.remove("particle_type_000");
	}
	kinJpp.set("INIT_C", std::vector<double>{2.7, 0.3});

	kinJpp.set("NREAC_LIQUID", 1);
	kinJpp.addScope("liquid_reaction_000");
	kinJpp.pushScope("liquid_reaction_000");
	kinJpp.set("TYPE", "MASS_ACTION_LAW");
	kinJpp.set("MAL_KFWD", std::vector<double>{2.0});
	kinJpp.set("MAL_KBWD", std::vector<double>{1.0});
	kinJpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
	kinJpp.set("MAL_IS_KINETIC", std::vector<int>{kinetic ? 1 : 0});
	kinJpp.popScope();

	cadet::IUnitOperation* const kinUnit =
		cadet::test::unitoperation::createAndConfigureUnit(kinJpp, *mb);
	REQUIRE(kinUnit->numDofs() == eqUnit->numDofs());

	// Use the equilibrium model's consistently initialized state
	std::vector<double> yKin = yEq;
	yKin[0] = 1.0;
	yKin[1] = 2.0;

	std::vector<double> residual(kinUnit->numDofs(), 0.0);
	cadet::util::ThreadLocalStorage tls;
	tls.resize(kinUnit->threadLocalMemorySize());

	const cadet::AdJacobianParams noAdParams{nullptr, nullptr, 0u};
	kinUnit->notifyDiscontinuousSectionTransition(
		0.0, 0u, {yKin.data(), nullptr}, noAdParams);
	kinUnit->residual(
		cadet::SimulationTime{0.0, 0u},
		cadet::ConstSimulationState{yKin.data(), nullptr},
		residual.data(), tls);

	const unsigned int eqnComp = kinUnit->numComponents();
	const unsigned int eqbulkOffset = nComp * kinUnit->numInletPorts();
	const unsigned int eqinletDofs = kinUnit->numComponents() * kinUnit->numInletPorts();
	const unsigned int eqcolumnBulkPoints =  (kinUnit->numDofs() - eqinletDofs) / eqnComp;

	for (unsigned int point = 0; point < eqcolumnBulkPoints; ++point)
	{
		const double* const c = yKin.data() + eqbulkOffset + point * nComp;
		CHECK(c[0] == cadet::test::makeApprox(1.0, 1e-12, 1e-12));
		CHECK(c[1] == cadet::test::makeApprox(2.0, 1e-12, 1e-12));
	}

	for (unsigned int dof = eqbulkOffset; dof < residual.size(); ++dof)
		CHECK(std::abs(residual[dof]) <= 1e-11);

	mb->destroyUnitOperation(eqUnit);
	mb->destroyUnitOperation(kinUnit);
	destroyModelBuilder(mb);
}
