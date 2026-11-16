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

#include "ColumnTests.hpp"
#include "ParticleHelper.hpp"
#include "ReactionModelTests.hpp"
#include "JsonTestModels.hpp"
#include "Weno.hpp"
#include "Utils.hpp"
#include "UnitOperationTests.hpp"
#include "cadet/cadet.hpp"
#include "model/UnitOperation.hpp"
#include "ParallelSupport.hpp"
#include "SimulationTypes.hpp"

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

TEST_CASE("GRM LWE forward vs backward flow", "[GRM],[FV],[Simulation],[CI]")
{
	cadet::test::column::FVParams disc;

	// Test all WENO orders
	for (unsigned int i = 1; i <= cadet::Weno::maxOrder(); ++i)
	{
		disc.setBulkDiscParam("WENO_ORDER", static_cast<int>(i));
		cadet::test::column::testForwardBackward("GENERAL_RATE_MODEL", disc, 1e-9, 2e-4);
	}
}

TEST_CASE("GRM linear pulse vs analytic solution", "[GRM],[FV],[Simulation],[Analytic],[CI]")
{
	cadet::test::column::FVParams disc(512);
	cadet::test::column::testAnalyticBenchmark("GENERAL_RATE_MODEL", "/data/grm-pulseBenchmark.data", true, true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("GENERAL_RATE_MODEL", "/data/grm-pulseBenchmark.data", true, false, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("GENERAL_RATE_MODEL", "/data/grm-pulseBenchmark.data", false, true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("GENERAL_RATE_MODEL", "/data/grm-pulseBenchmark.data", false, false, disc, 6e-5, 1e-7);
}

TEST_CASE("GRM non-binding linear pulse vs analytic solution", "[GRM],[FV],[Simulation],[Analytic],[NonBinding],[CI]")
{
	cadet::test::column::FVParams disc(512);
	cadet::test::column::testAnalyticNonBindingBenchmark("GENERAL_RATE_MODEL", "/data/grm-nonBinding.data", true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticNonBindingBenchmark("GENERAL_RATE_MODEL", "/data/grm-nonBinding.data", false, disc, 6e-5, 1e-7);
}

TEST_CASE("GRM Jacobian forward vs backward flow", "[GRM],[FV],[UnitOp],[Residual],[Jacobian],[AD],[CI]")
{
	cadet::test::column::FVParams disc;

	// Test all WENO orders
	for (unsigned int i = 1; i <= cadet::Weno::maxOrder(); ++i)
	{
		disc.setBulkDiscParam("WENO_ORDER", static_cast<int>(i));
		cadet::test::column::testJacobianForwardBackward("GENERAL_RATE_MODEL", disc);
	}
}

TEST_CASE("GRM numerical Benchmark1 with parameter sensitivities for linear case", "[GRM],[FV],[Simulation],[Reference],[Sensitivity],[CI_sens]")
{
	std::string modelFilePath = std::string("/data/model_COL1D_GRM_dynLin_1comp_benchmark1.json");
	std::string refFilePath = std::string("/data/ref_GRM_dynLin_1comp_sensbenchmark1_FV_Z32parZ4.h5");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1e-4, 1e-4, 1e-4, 1e-4 };

	cadet::test::column::FVParams disc(32, 4);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("GRM numerical Benchmark2 with parameter sensitivities for linear case", "[GRM],[FV],[Simulation],[Reference],[Sensitivity],[CI_sens]")
{
	std::string modelFilePath = std::string("/data/model_COL1D_GRM_dynLin_1comp_benchmark2.json");
	std::string refFilePath = std::string("/data/ref_GRM_dynLin_1comp_sensbenchmark2_FV_Z32parZ4.h5");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1e-4, 1e-4, 1e-4, 1e-4, 1e-4 };

	cadet::test::column::FVParams disc(32, 4);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("GRM numerical Benchmark with parameter sensitivities and multiplexing for 2parType 2comp linear case", "[GRM],[FV],[Simulation],[Reference],[Sensitivity],[CI_sens]")
{
	const std::string modelFilePath = std::string("/data/model_COL1DparType2_GRM_dynLin_2comp_sensbenchmark1.json");
	const std::string refFilePath = std::string("/data/ref_GRMparType2_dynLin_2comp_sensbenchmark1_FV_Z16parZ8.h5");
	const std::vector<double> absTol = { 1e-12, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6 };
	const std::vector<double> relTol = { 1e-4, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1 };

	const cadet::test::column::FVParams disc(16, 8);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("GRM numerical Benchmark with parameter sensitivities for SMA LWE case", "[GRM],[FV],[Simulation],[Reference],[Sensitivity],[CI_sens]")
{
	const std::string& modelFilePath = std::string("/data/model_COL1D_GRM_reqSMA_4comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_GRM_reqSMA_4comp_sensbenchmark1_FV_Z16parZ2.h5");
	// The sensitivities of param_000 and param_001 peak at about 3e7 and 1e9 and deviate from the reference by
	// up to 1e-6 of their peak. The absolute tolerance only applies near their zero crossings, where 1e-6 was
	// below double precision relative to the peak (measured deviation there up to 5.4e-6).
	const std::vector<double> absTol = { 1e-4, 1e-3, 1e-3, 1e-6 };
	const std::vector<double> relTol = { 1.0, 1.0, 1.0, 1.0 };

	cadet::test::column::FVParams disc(16, 2);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "000", absTol, relTol, disc, true);
}

TEST_CASE("GRM time derivative Jacobian vs FD", "[GRM],[FV],[UnitOp],[Residual],[Jacobian],[CI],[FD]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("GENERAL_RATE_MODEL", "FV", 1e-6, 0.0, 9e-4);
}

TEST_CASE("GRM rapid-equilibrium binding flux Jacobian vs FD", "[GRM],[FV],[UnitOp],[Residual],[Jacobian],[CI],[FD]")
{
	cadet::test::column::testArrowHeadJacobianFD("GENERAL_RATE_MODEL", false, 1e-6, 2e-9);
}

TEST_CASE("GRM rapid-equilibrium binding with surf diff par dep flux Jacobian vs FD", "[GRM],[FV],[UnitOp],[Residual],[Jacobian],[ParameterDependence],[CI],[FD]")
{
	cadet::test::column::testArrowHeadJacobianFDVariableParSurfDiff("GENERAL_RATE_MODEL", 1e-6, 5e-9);
}

TEST_CASE("GRM dynamic binding with surf diff par dep Jacobian vs AD", "[GRM],[FV],[UnitOp],[Residual],[Jacobian],[ParameterDependence],[CI]")
{
	cadet::test::column::testJacobianADVariableParSurfDiff("GENERAL_RATE_MODEL", "FV", true);
}

TEST_CASE("GRM rapid-equilibrium binding with surf diff par dep Jacobian vs AD", "[GRM],[FV],[UnitOp],[Residual],[Jacobian],[ParameterDependence],[CI]")
{
	cadet::test::column::testJacobianADVariableParSurfDiff("GENERAL_RATE_MODEL", "FV", false);
}

TEST_CASE("GRM dynamic binding flux Jacobian vs FD", "[GRM],[FV],[UnitOp],[Residual],[Jacobian],[CI],[FD]")
{
	cadet::test::column::testArrowHeadJacobianFD("GENERAL_RATE_MODEL", true, 1e-6, 2e-9);
}

TEST_CASE("GRM sensitivity Jacobians", "[GRM],[FV],[UnitOp],[Sensitivity],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("GENERAL_RATE_MODEL", "FV");

	cadet::test::column::testFwdSensJacobians(jpp, 1e-4, 6e-7);
}

//TEST_CASE("GRM forward sensitivity vs FD", "[GRM],[FV],[Sensitivity],[Simulation],[failedFDtestGRM]") // todo fix. fails for SMA_KA for both binding modes
//{
//	// todo comment
//	// Relative error is checked first, we use high absolute error for letting
//	// some points that are far off pass the error test, too. This is required
//	// due to errors in finite differences.
//	const double fdStepSize[] = {5e-5, 1e-4, 1e-4, 1e-3};
//	const double absTols[] = {3e5, 2e-3, 2e-4, 5.0};
//	const double relTols[] = {5e-3, 7e-2, 8e-2, 1e-4};
//	const double passRatio[] = {0.95, 0.9, 0.91, 0.83};
//	cadet::test::column::testFwdSensSolutionFD("GENERAL_RATE_MODEL", false, fdStepSize, absTols, relTols, passRatio);
//}
//
//TEST_CASE("GRM forward sensitivity forward vs backward flow", "[GRM],[FV],[Sensitivity],[Simulation],[fixGRM]") // todo fix. fails for COL_DISPERION for both binding modes
//{
//	const double absTols[] = {4e-5, 1e-11, 1e-11, 8e-9};
//	const double relTols[] = {6e-9, 5e-8, 5e-6, 5e-10};
//	const double passRatio[] = {0.99, 0.95, 0.98, 0.98};
//	cadet::test::column::testFwdSensSolutionForwardBackward("GENERAL_RATE_MODEL", absTols, relTols, passRatio);
//}

TEST_CASE("GRM consistent initialization with linear binding", "[GRM],[FV],[ConsistentInit],[CI]")
{
	cadet::test::column::testConsistentInitializationLinearBinding("GENERAL_RATE_MODEL", "FV", 1e-12, 1e-12);
}

TEST_CASE("GRM consistent initialization with a non-binding component", "[GRM],[FV],[ConsistentInit],[CI]")
{
	cadet::test::column::testConsistentInitializationNonBindingComponent("GENERAL_RATE_MODEL", "FV", 1e-12, 1e-12);
}

//TEST_CASE("GRM consistent initialization with SMA binding", "[GRM],[FV],[ConsistentInit],[fixGRM]") // todo fix
//{
//	std::vector<double> y(4 + 4 * 16 + 16 * 4 * (4 + 4) + 4 * 16, 0.0);
//// Optimal values:
////	const double bindingCell[] = {1.2, 2.0, 1.0, 1.5, 858.034, 66.7896, 3.53273, 2.53153, 
////		1.0, 1.8, 1.5, 1.6, 856.173, 64.457, 5.73227, 2.85286};
//	const double bindingCell[] = {1.2, 2.0, 1.0, 1.5, 840.0, 63.0, 3.0, 3.0, 
//		1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0};
//	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 + 4 * 16);
//	cadet::test::util::repeat(y.data() + 4 + 4 * 16, bindingCell, 16, 4 * 16 / 2);
//	cadet::test::util::populate(y.data() + 4 + 4 * 16 + 16 * 4 * (4 + 4), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 16);
//
//	cadet::test::column::testConsistentInitializationSMABinding("GENERAL_RATE_MODEL", y.data(), 1e-14, 1e-5);
//}

TEST_CASE("GRM consistent sensitivity initialization with linear binding", "[GRM],[FV],[ConsistentInit],[Sensitivity],[CI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 4 + 4 * 16 + 16 * 4 * (4 + 4) + 4 * 16;
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, numDofs);
	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	cadet::test::column::testConsistentInitializationSensitivity("GENERAL_RATE_MODEL", "FV", y.data(), yDot.data(), true, 1e-14);
}

TEST_CASE("GRM consistent sensitivity initialization with SMA binding", "[GRM],[FV],[ConsistentInit],[Sensitivity],[CI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 4 + 4 * 16 + 16 * 4 * (4 + 4) + 4 * 16;
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);

	const double bindingCell[] = {1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0};
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 + 4 * 16);
	cadet::test::util::repeat(y.data() + 4 + 4 * 16, bindingCell, 8, 4 * 16);
	cadet::test::util::populate(y.data() + 4 + 4 * 16 + 16 * 4 * (4 + 4), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 16);

	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	cadet::test::column::testConsistentInitializationSensitivity("GENERAL_RATE_MODEL", "FV", y.data(), yDot.data(), false, 1e-9);
}

TEST_CASE("GRM inlet DOF Jacobian", "[GRM],[FV],[UnitOp],[Jacobian],[Inlet],[CI]")
{
	cadet::test::column::testInletDofJacobian("GENERAL_RATE_MODEL", "FV");
}

TEST_CASE("GRM transport Jacobian", "[GRM],[FV],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "GENERAL_RATE_MODEL", "FV");
	cadet::test::column::testJacobianAD(jpp);
}

TEST_CASE("GRM with two component linear binding Jacobian", "[GRM],[FV],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("GENERAL_RATE_MODEL", "FV");
	cadet::test::column::testJacobianAD(jpp);
}

TEST_CASE("GRM Jacobian follows a dynamic flow rate", "[GRM],[FV],[UnitOp],[Jacobian],[Inlet],[CI]")
{
	for (int dir = 0; dir < 2; ++dir)
	{
		SECTION(dir == 0 ? "Forward flow" : "Backward flow")
		{
			cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("GENERAL_RATE_MODEL", "FV");
			if (dir != 0)
				cadet::test::column::reverseFlow(jpp);

			cadet::test::column::testJacobianDynamicFlowRate(jpp, 3.0);
		}
	}
}

TEST_CASE("GRM LWE one vs two identical particle types match", "[GRM],[FV],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testOneVsTwoIdenticalParticleTypes("GENERAL_RATE_MODEL", "FV", 2e-8, 5e-5);
}

TEST_CASE("GRM LWE separate identical particle types match", "[GRM],[FV],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testSeparateIdenticalParticleTypes("GENERAL_RATE_MODEL", "FV", 1e-9, 1e-8);
}

TEST_CASE("GRM linear binding single particle matches particle distribution", "[GRM],[FV],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearMixedParticleTypes("GENERAL_RATE_MODEL", "FV", 5e-8, 5e-5);
}

TEST_CASE("GRM multiple particle types Jacobian analytic vs AD", "[GRM],[FV],[Jacobian],[AD],[ParticleType],[CI]")
{
	cadet::test::particle::testJacobianMixedParticleTypes("GENERAL_RATE_MODEL", "FV");
}

TEST_CASE("GRM multiple particle types time derivative Jacobian vs FD", "[GRM],[FV],[UnitOp],[Residual],[Jacobian],[ParticleType],[CI]")
{
	cadet::test::particle::testTimeDerivativeJacobianMixedParticleTypesFD("GENERAL_RATE_MODEL", "FV", 1e-6, 0.0, 9e-4);
}

TEST_CASE("GRM multiple spatially dependent particle types Jacobian analytic vs AD", "[GRM],[FV],[Jacobian],[AD],[ParticleType],[CI]")
{
	cadet::test::particle::testJacobianSpatiallyMixedParticleTypes("GENERAL_RATE_MODEL", "FV");
}

TEST_CASE("GRM linear binding single particle matches spatially dependent particle distribution", "[GRM],[FV],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearSpatiallyMixedParticleTypes("GENERAL_RATE_MODEL", "FV", 5e-8, 5e-5);
}

TEST_CASE("GRM multiple spatially dependent particle types flux Jacobian vs FD", "[GRM],[FV],[UnitOp],[Residual],[Jacobian],[ParticleType],[CI],[FD]")
{
	cadet::test::particle::testArrowHeadJacobianSpatiallyMixedParticleTypes("GENERAL_RATE_MODEL", 1e-6, 1e-8, 1e-5);
}

TEST_CASE("GRM dynamic reactions Jacobian vs AD bulk", "[GRM],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("GENERAL_RATE_MODEL", "FV", true, false, false);
}

TEST_CASE("GRM dynamic reactions Jacobian vs AD particle", "[GRM],[FV],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("GENERAL_RATE_MODEL", "FV", false, true, false);
}

TEST_CASE("GRM dynamic reactions Jacobian vs AD modified particle", "[GRM],[FV],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("GENERAL_RATE_MODEL", "FV", false, true, true);
}

TEST_CASE("GRM dynamic reactions Jacobian vs AD bulk and particle", "[GRM],[FV],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("GENERAL_RATE_MODEL", "FV", true, true, false);
}

TEST_CASE("GRM dynamic reactions Jacobian vs AD bulk and modified particle", "[GRM],[FV],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("GENERAL_RATE_MODEL", "FV", true, true, true);
}

TEST_CASE("GRM dynamic reactions time derivative Jacobian vs FD bulk", "[GRM],[FV],[Jacobian],[Residual],[ReactionModel],[CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("GENERAL_RATE_MODEL", "FV", true, false, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("GRM dynamic reactions time derivative Jacobian vs FD particle", "[GRM],[FV],[Jacobian],[Residual],[ReactionModel],[CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("GENERAL_RATE_MODEL", "FV", false, true, false, 1e-4, 1e-14, 9e-4);
}

TEST_CASE("GRM dynamic reactions time derivative Jacobian vs FD modified particle", "[GRM],[FV],[Jacobian],[Residual],[ReactionModel],[CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("GENERAL_RATE_MODEL", "FV", false, true, true, 1e-4, 1e-14, 9e-4);
}

TEST_CASE("GRM dynamic reactions time derivative Jacobian vs FD bulk and particle", "[GRM],[FV],[Jacobian],[Residual],[ReactionModel],[CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("GENERAL_RATE_MODEL", "FV", true, true, false, 1e-4, 1e-14, 9e-4);
}

TEST_CASE("GRM dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[GRM],[FV],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("GENERAL_RATE_MODEL", "FV", true, true, true, 1e-4, 1e-14, 9e-4);
}

inline cadet::JsonParameterProvider createColumnWithTwoCompLinearBindingThreeParticleTypesGRM()
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("GENERAL_RATE_MODEL", "FV");

	const double parVolFrac[] = {0.3, 0.6, 0.1};
	const double parFactor[] = {0.9, 0.8};
	cadet::test::particle::extendModelToManyParticleTypes(jpp, 3, parFactor, parVolFrac);

	return jpp;
}

TEST_CASE("GRM multi particle types dynamic reactions Jacobian vs AD bulk", "[GRM],[FV],[Jacobian],[AD],[ReactionModel],[ParticleType],[fixParReacTest]") // works localy in debug mode Column_1D also works
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypesGRM();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, false, false);
}

TEST_CASE("GRM multi particle types dynamic reactions Jacobian vs AD particle", "[GRM],[Jacobian],[FV],[AD],[ReactionModel],[ParticleType],[fixParReacTest]") // works localy in debug mode Column_1D also works
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypesGRM();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, false);
}

TEST_CASE("GRM multi particle types dynamic reactions Jacobian vs AD modified particle", "[GRM],[FV],[Jacobian],[AD],[ReactionModel],[ParticleType],[fixParReacTest]") // works localy in debug mode Column_1D also works
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypesGRM();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, true);
}

TEST_CASE("GRM multi particle types dynamic reactions Jacobian vs AD bulk and particle", "[GRM],[FV],[Jacobian],[AD],[ReactionModel],[ParticleType],[fixParReacTest]") // works localy in debug mode Column_1D also works
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypesGRM();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, false); 
}

TEST_CASE("GRM multi particle types dynamic reactions Jacobian vs AD bulk and modified particle", "[GRM],[FV],[Jacobian],[AD],[ReactionModel],[ParticleType],[fixParReacTest]") // works localy in debug mode Column_1D also works
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypesGRM();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, true);
}	

TEST_CASE("GRM multi particle types dynamic reactions time derivative Jacobian vs FD bulk", "[GRM],[FV],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI],[FD]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypesGRM();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, false, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("GRM multi particle types dynamic reactions time derivative Jacobian vs FD particle", "[GRM],[FV],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI],[FD]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypesGRM();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, false, 1e-4, 1e-14, 9e-4);
}

TEST_CASE("GRM multi particle types dynamic reactions time derivative Jacobian vs FD modified particle", "[GRM],[FV],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI],[FD]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypesGRM();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, true, 1e-4, 1e-14, 9e-4);
}

TEST_CASE("GRM multi particle types dynamic reactions time derivative Jacobian vs FD bulk and particle", "[GRM],[FV],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI],[FD]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypesGRM();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, false, 1e-4, 1e-14, 9e-4);
}

TEST_CASE("GRM multi particle types dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[GRM],[FV],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI],[FD]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypesGRM();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, true, 1e-4, 1e-14, 9e-4);
}

namespace
{
	/**
	 * @brief Sets geometry, core radius and radial discretization of the first particle type of a unit operation
	 * @param [in,out] jpp Parameter provider scoped to the unit operation
	 * @param [in] geom Value of PAR_GEOM
	 * @param [in] relCoreRadius Ratio of PAR_CORERADIUS and PAR_RADIUS
	 * @param [in] discType Value of PAR_DISC_TYPE
	 * @param [in] relInterfaces Relative cell interfaces (only used for USER_DEFINED)
	 */
	void setParticleGrid(cadet::JsonParameterProvider& jpp, const std::string& geom, double relCoreRadius, const std::string& discType, const std::vector<double>& relInterfaces)
	{
		jpp.pushScope("particle_type_000");
		jpp.set("PAR_GEOM", geom);
		jpp.set("PAR_CORERADIUS", relCoreRadius * jpp.getDouble("PAR_RADIUS"));
		jpp.pushScope("discretization");
		jpp.set("PAR_DISC_TYPE", discType);
		if (discType == "USER_DEFINED")
			jpp.set("PAR_DISC_VECTOR", relInterfaces);
		jpp.popScope();
		jpp.popScope();
	}

	/**
	 * @brief Returns the relative interfaces in [0, 1] of cells with equal volume
	 * @param [in] geomExponent Exponent d of the particle geometry, i.e. 3 for spheres, 2 for cylinders and 1 for slabs
	 * @param [in] relCoreRadius Ratio of core radius and particle radius
	 * @param [in] nCells Number of cells
	 */
	std::vector<double> relativeEquivolumeInterfaces(double geomExponent, double relCoreRadius, unsigned int nCells)
	{
		std::vector<double> relInterfaces(nCells + 1);
		const double coreVol = std::pow(relCoreRadius, geomExponent);
		for (unsigned int k = 0; k <= nCells; ++k)
		{
			const double r = std::pow(coreVol + static_cast<double>(k) / nCells * (1.0 - coreVol), 1.0 / geomExponent);
			relInterfaces[k] = (r - relCoreRadius) / (1.0 - relCoreRadius);
		}
		relInterfaces.front() = 0.0;
		relInterfaces.back() = 1.0;
		return relInterfaces;
	}

	/**
	 * @brief Evaluates the residual of a unit operation at a fixed sample state
	 * @param [in,out] jpp Parameter provider scoped to the unit operation
	 */
	std::vector<double> residualAtSampleState(cadet::JsonParameterProvider& jpp)
	{
		cadet::IModelBuilder* const mb = cadet::createModelBuilder();
		REQUIRE(nullptr != mb);

		cadet::IUnitOperation* const unit = cadet::test::unitoperation::createAndConfigureUnit(jpp, *mb);
		const unsigned int nDof = unit->numDofs();
		std::vector<double> y(nDof, 0.0);
		std::vector<double> res(nDof, 0.0);
		cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, nDof);

		cadet::util::ThreadLocalStorage tls;
		tls.resize(unit->threadLocalMemorySize());
		const cadet::AdJacobianParams noParams{ nullptr, nullptr, 0u };
		unit->notifyDiscontinuousSectionTransition(0.0, 0u, { y.data(), nullptr }, noParams);
		unit->residual(cadet::SimulationTime{ 0.0, 0u }, cadet::ConstSimulationState{ y.data(), nullptr }, res.data(), tls);

		mb->destroyUnitOperation(unit);
		cadet::destroyModelBuilder(mb);
		return res;
	}

	/**
	 * @brief Returns the maximum deviation of two residuals relative to the largest absolute entry of the first one
	 */
	double relativeResidualDeviation(const std::vector<double>& a, const std::vector<double>& b)
	{
		REQUIRE(a.size() == b.size());
		double maxDiff = 0.0;
		double maxAbs = 0.0;
		for (std::size_t i = 0; i < a.size(); ++i)
		{
			maxDiff = std::max(maxDiff, std::abs(a[i] - b[i]));
			maxAbs = std::max(maxAbs, std::abs(a[i]));
		}
		return maxDiff / maxAbs;
	}
}

TEST_CASE("Column units with FV particles map USER_DEFINED particle cells to the particle radius", "[GRM],[FV],[DGFV],[GRM2D],[UnitOp],[Residual],[Jacobian],[AD],[CI]")
{
	// A user defined grid with the interfaces of the equidistant (equivolume) grid has to reproduce the
	// EQUIDISTANT (EQUIVOLUME) discretization, which requires the relative interfaces in PAR_DISC_VECTOR to be
	// mapped onto [PAR_CORERADIUS, PAR_RADIUS]. This is checked for the arrow head FV GRM, the 1D column model
	// with FV particles and the 2D FV GRM. Additionally, the analytic Jacobian of a non-uniform user defined grid
	// has to match AD (not available for the 2D FV GRM). The 2D FV GRM is only checked with spherical particles,
	// as its residual does not depend on PAR_GEOM.
	const std::vector<std::pair<std::string, std::string>> units = { { "GENERAL_RATE_MODEL", "FV" }, { "COLUMN_MODEL_1D_GRM", "DGFV" }, { "COLUMN_MODEL_2D_GRM", "FV" } };
	for (const auto& unit : units)
	{
		for (const std::string geom : { "SPHERE", "CYLINDER", "SLAB" })
		{
			if ((unit.first == "COLUMN_MODEL_2D_GRM") && (geom != "SPHERE"))
				continue;

			const double geomExponent = (geom == "SPHERE") ? 3.0 : ((geom == "CYLINDER") ? 2.0 : 1.0);
			for (const double relCoreRadius : { 0.0, 0.3 })
			{
				INFO("Unit " << unit.first << " (" << unit.second << "), geometry " << geom << ", PAR_CORERADIUS / PAR_RADIUS " << relCoreRadius);

				cadet::JsonParameterProvider jppRef = createColumnWithTwoCompLinearBinding(unit.first, unit.second);
				jppRef.pushScope("particle_type_000");
				jppRef.pushScope("discretization");
				const unsigned int nCells = jppRef.getInt("NCELLS");
				jppRef.popScope();
				jppRef.popScope();

				std::vector<double> relEquidistant(nCells + 1);
				std::vector<double> relNonUniform(nCells + 1);
				for (unsigned int k = 0; k <= nCells; ++k)
				{
					relEquidistant[k] = static_cast<double>(k) / nCells;
					relNonUniform[k] = std::sqrt(static_cast<double>(k) / nCells);
				}

				for (const std::string discType : { "EQUIDISTANT", "EQUIVOLUME" })
				{
					cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding(unit.first, unit.second);
					setParticleGrid(jpp, geom, relCoreRadius, discType, {});
					const std::vector<double> resRef = residualAtSampleState(jpp);

					cadet::JsonParameterProvider jppUser = createColumnWithTwoCompLinearBinding(unit.first, unit.second);
					setParticleGrid(jppUser, geom, relCoreRadius, "USER_DEFINED", (discType == "EQUIDISTANT") ? relEquidistant : relativeEquivolumeInterfaces(geomExponent, relCoreRadius, nCells));
					const std::vector<double> resUser = residualAtSampleState(jppUser);

					const double deviation = relativeResidualDeviation(resRef, resUser);
					CAPTURE(discType, deviation);
					CHECK(deviation <= 1e-10);
				}

				if (unit.first != "COLUMN_MODEL_2D_GRM")
				{
					cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding(unit.first, unit.second);
					setParticleGrid(jpp, geom, relCoreRadius, "USER_DEFINED", relNonUniform);
					cadet::test::column::testJacobianAD(jpp);
				}
			}
		}
	}
}
