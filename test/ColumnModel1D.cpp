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
#include "Approx.hpp"

#define CADET_LOGGING_DISABLE
#include "Logging.hpp"

#include "ColumnTests.hpp"
#include "ParticleHelper.hpp"
#include "ReactionModelTests.hpp"
#include "JsonTestModels.hpp"
#include "Utils.hpp"
#include "common/Driver.hpp"

TEST_CASE("Column_Model_1D LWE forward vs backward flow", "[Column_1D],[DG],[DG1D],[Simulation],[CI]")
{
	cadet::test::column::DGparams disc;

	// Test all integration modes
	for (int i = 0; i <= 1; i++)
	{
		disc.setIntegrationMode(i);
		cadet::test::column::testForwardBackward("COLUMN_MODEL_1D", disc, 1e-9, 2e-4);
	}
}

TEST_CASE("Column_Model_1D linear pulse vs analytic solution", "[Column_1D],[DG],[DG1D],[Simulation],[Analytic],[todo_CI]")
{
	cadet::test::column::DGparams disc;
	cadet::test::column::testAnalyticBenchmark("GENERAL_RATE_MODEL", "/data/grm-pulseBenchmark.data", true, true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("GENERAL_RATE_MODEL", "/data/grm-pulseBenchmark.data", true, false, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("GENERAL_RATE_MODEL", "/data/grm-pulseBenchmark.data", false, true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("GENERAL_RATE_MODEL", "/data/grm-pulseBenchmark.data", false, false, disc, 6e-5, 1e-7);
}

TEST_CASE("Column_Model_1D non-binding linear pulse vs analytic solution", "[Column_1D],[DG],[DG1D],[Simulation],[Analytic],[NonBinding],[todo_CI]")
{
	cadet::test::column::DGparams disc;
	cadet::test::column::testAnalyticNonBindingBenchmark("GENERAL_RATE_MODEL", "/data/grm-nonBinding.data", true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticNonBindingBenchmark("GENERAL_RATE_MODEL", "/data/grm-nonBinding.data", false, disc, 6e-5, 1e-7);
}

// todo FIX (scheitert bei backward flow jacobian vs AD)
//TEST_CASE("Column_Model_1D Jacobian forward vs backward flow", "[Column_1D],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[AD],[todo]")
//{
//	cadet::test::column::DGparams disc;
//
//	// Test all integration modes
//	for (int i = 0; i <= 1; i++)
//	{
//		disc.setIntegrationMode(i);
//		cadet::test::column::testJacobianForwardBackward("GENERAL_RATE_MODEL", disc, std::numeric_limits<float>::epsilon() * 100.0);
//	}
//}

TEST_CASE("Column_Model_1D numerical Benchmark with parameter sensitivities for linear case", "[Column_1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[todo_CI_sens3]")
{
	std::string modelFilePath = std::string("/data/model_GRM_dynLin_1comp_benchmark1.json");
	std::string refFilePath = std::string("/data/ref_GRM_dynLin_1comp_sensbenchmark1_cDG_P3Z8_GSM_parP3parZ1.h5");
	const std::vector<double> absTol = { 1e-12, 1e-6, 1e-6, 1e-12 };
	const std::vector<double> relTol = { 1e-4, 1e-3, 1e-4, 1e-4 };

	cadet::test::column::DGparams disc(0, 3, 8, 3, 1);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);

	modelFilePath = std::string("/data/model_GRM_dynLin_1comp_sensbenchmark2.json");
	refFilePath = std::string("/data/ref_GRM_dynLin_1comp_sensbenchmark2_cDG_P3Z8_GSM_parP3parZ1.h5");

	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("Column_Model_1D numerical Benchmark with parameter sensitivities and multiplexing for 2parType 2comp linear case", "[Column_1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[todo_CI_sens15]")
{
	const std::string modelFilePath = std::string("/data/model_GRMparType2_dynLin_2comp_sensbenchmark1.json");
	const std::string refFilePath = std::string("/data/ref_GRMparType2_dynLin_2comp_sensbenchmark1_cDG_P3Z4_DGexInt_parP3parZ2.h5");
	const std::vector<double> absTol = { 1e-12, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6 };
	const std::vector<double> relTol = { 1e-4, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1 };

	const cadet::test::column::DGparams disc(0, 3, 4, 3, 2);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("Column_Model_1D numerical Benchmark with parameter sensitivities for SMA LWE case", "[Column_1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[todo_CI_sens4]")
{
	const std::string modelFilePath = std::string("/data/model_GRM_reqSMA_4comp_benchmark1.json");
	const std::string refFilePath = std::string("/data/ref_GRM_reqSMA_4comp_sensbenchmark1_cDG_P3Z8_GSM_parP3parZ1.h5");
	const std::vector<double> absTol = { 1e-12, 1e-6, 1e-6, 1e-12 };
	const std::vector<double> relTol = { 1.0, 1.0, 1.0, 1.0 };

	cadet::test::column::DGparams disc(0, 3, 8, 3, 1);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "000", absTol, relTol, disc, true);
}

TEST_CASE("Column_Model_1D LWE DGSEM and GSM particle discretization yields similar accuracy", "[Column_1D],[DG],[DG1D],[Simulation],[todo_CI]")
{
	cadet::JsonParameterProvider jpp = createLWE("GENERAL_RATE_MODEL", "DG");
	cadet::test::column::DGparams disc(0, 4, 2, 3, 1); // Note that we want to employ only a single particle element
	disc.setDisc(jpp);

	const double absTol = 1e-9;
	const double relTol = 2e-4;

	// GSM discretization
	cadet::Driver drvGSM;
	drvGSM.configure(jpp);
	drvGSM.run();

	// Force single element DGSEM discretization (GSM is default for single particle element discretization)
	jpp.pushScope("model");
	jpp.pushScope("unit_000");
	jpp.pushScope("discretization");
	jpp.set("PAR_GSM", false);
	jpp.popScope();
	jpp.popScope();
	jpp.popScope();

	cadet::Driver drvDGSEM;
	drvDGSEM.configure(jpp);
	drvDGSEM.run();

	cadet::InternalStorageUnitOpRecorder const* const GSMData = drvGSM.solution()->unitOperation(0);
	cadet::InternalStorageUnitOpRecorder const* const DGSEMData = drvDGSEM.solution()->unitOperation(0);

	double const* GSMOutlet = GSMData->outlet();
	double const* DGSEMOutlet = DGSEMData->outlet();

	const unsigned int nComp = GSMData->numComponents();
	for (unsigned int i = 0; i < GSMData->numDataPoints() * GSMData->numInletPorts() * nComp; ++i, ++GSMOutlet, ++DGSEMOutlet)
	{
		// Forward flow inlet = backward flow outlet
		CAPTURE(i);
		CHECK((*GSMOutlet) == cadet::test::makeApprox(*DGSEMOutlet, relTol, absTol));
	}
}

TEST_CASE("Column_Model_1D time derivative Jacobian vs FD", "[Column_1D],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[todo_CI],[FD]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("GENERAL_RATE_MODEL", "DG", 1e-6, 0.0, 9e-4);
}

//TEST_CASE("Column_Model_1D dynamic binding with surf diff par dep Jacobian vs AD", "[Column_1D],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[ParameterDependence],[fix]")
//{
//	cadet::test::column::testJacobianADVariableParSurfDiff("GENERAL_RATE_MODEL", "DG", true);
//}

TEST_CASE("Column_Model_1D rapid-equilibrium binding with surf diff par dep Jacobian vs AD", "[Column_1D],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[ParameterDependence]") // todo include in CI (runs locally but fails on server with linux)
{
	cadet::test::column::testJacobianADVariableParSurfDiff("GENERAL_RATE_MODEL", "DG", false);
}

TEST_CASE("Column_Model_1D sensitivity Jacobians", "[Column_1D],[DG],[DG1D],[UnitOp],[Sensitivity]") // todo does not run on CI but locally on windows
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("GENERAL_RATE_MODEL", "DG");

	cadet::test::column::testFwdSensJacobians(jpp, 1e-4, 6e-7);
}

//// todo fix: not just adjust tolerances as in FV but theres an actual error here: access violation in densematrix
//TEST_CASE("Column_Model_1D forward sensitivity vs FD", "[Column_1D],[DG],[DG1D],[Sensitivity],[Simulation],[todo]")
//{
//	// Relative error is checked first, we use high absolute error for letting
//	// some points that are far off pass the error test, too. This is required
//	// due to errors in finite differences.
//	const double fdStepSize[] = { 5e-5, 1e-4, 1e-4, 1e-3 };
//	const double absTols[] = { 3e5, 2e-3, 2e-4, 5.0 };
//	const double relTols[] = { 5e-3, 7e-2, 8e-2, 1e-4 };
//	const double passRatio[] = { 0.95, 0.9, 0.91, 0.83 };
//	cadet::test::column::testFwdSensSolutionFD("GENERAL_RATE_MODEL", "DG", false, fdStepSize, absTols, relTols, passRatio);
//}

//// todo fix: not just adjust tolerances as in FV but theres an actual error here: access violation in densematrix
//TEST_CASE("Column_Model_1D forward sensitivity forward vs backward flow", "[Column_1D],[DG],[DG1D],[Sensitivity],[Simulation],[todo]")
//{
//	const double absTols[] = { 4e-5, 1e-11, 1e-11, 8e-9 };
//	const double relTols[] = { 6e-9, 5e-8, 5e-6, 5e-10 };
//	const double passRatio[] = { 0.99, 0.95, 0.98, 0.98 };
//	cadet::test::column::testFwdSensSolutionForwardBackward("GENERAL_RATE_MODEL", "DG", absTols, relTols, passRatio);
//}

// todo fix consistent initialization for AD with req binding
TEST_CASE("Column_Model_1D consistent initialization with linear binding", "[Column_1D],[DG],[DG1D],[ConsistentInit],[todo_CI]")
{
	cadet::test::column::testConsistentInitializationLinearBinding("GENERAL_RATE_MODEL", "DG", 1e-12, 1e-14, 0, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("GENERAL_RATE_MODEL", "DG", 1e-12, 1e-12, 1, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("GENERAL_RATE_MODEL", "DG", 1e-12, 1e-14, 0, 1);
	//cadet::test::column::testConsistentInitializationLinearBinding("GENERAL_RATE_MODEL", "DG", 1e-12, 1e-14, 1, 1);
}

//// todo fix consistent initialization for SMA (initialization not completely correct; AD gives assertion error)
//TEST_CASE("Column_Model_1D consistent initialization with SMA binding", "[Column_1D],[DG],[DG1D],[ConsistentInit],[todo]")
//{
//	std::vector<double> y(4 + 4 * 16 + 16 * 4 * (4 + 4) + 4 * 16, 0.0);
//	// Optimal values:
//	//	const double bindingCell[] = {1.2, 2.0, 1.0, 1.5, 858.034, 66.7896, 3.53273, 2.53153, 
//	//		1.0, 1.8, 1.5, 1.6, 856.173, 64.457, 5.73227, 2.85286};
//	const double bindingCell[] = { 1.2, 2.0, 1.0, 1.5, 840.0, 63.0, 3.0, 3.0,
//		1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0 };
//	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 + 4 * 16);
//	cadet::test::util::repeat(y.data() + 4 + 4 * 16, bindingCell, 16, 4 * 16 / 2);
//	cadet::test::util::populate(y.data() + 4 + 4 * 16 + 16 * 4 * (4 + 4), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 16);
//
//	cadet::test::column::testConsistentInitializationSMABinding("GENERAL_RATE_MODEL", "DG", y.data(), 1e-14, 1e-5);
//}

// todo fix kinetic binding sensitivity init
TEST_CASE("Column_Model_1D consistent sensitivity initialization with linear binding", "[Column_1D],[DG],[DG1D],[ConsistentInit],[Sensitivity],[todo_CI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 4 + 4 * 16 + 16 * 4 * (4 + 4) + 4 * 16;
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, numDofs);
	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	//cadet::test::column::testConsistentInitializationSensitivity("GENERAL_RATE_MODEL", "DG", y.data(), yDot.data(), true, 1e-14, 0, 0);
	cadet::test::column::testConsistentInitializationSensitivity("GENERAL_RATE_MODEL", "DG", y.data(), yDot.data(), true, 1e-14, 1, 0);
	//cadet::test::column::testConsistentInitializationSensitivity("GENERAL_RATE_MODEL", "DG", y.data(), yDot.data(), true, 1e-14, 0, 1);
	cadet::test::column::testConsistentInitializationSensitivity("GENERAL_RATE_MODEL", "DG", y.data(), yDot.data(), true, 1e-14, 1, 1);
}

//// todo fix memory stuff (works for FV) 
//TEST_CASE("Column_Model_1D consistent sensitivity initialization with SMA binding", "[Column_1D],[DG],[DG1D],[ConsistentInit],[Sensitivity],[fffffffiujbnlk]")
//{
//	// Fill state vector with given initial values
//	const unsigned int numDofs = 4 + 4 * 16 + 16 * 4 * (4 + 4) + 4 * 16;
//	std::vector<double> y(numDofs, 0.0);
//	std::vector<double> yDot(numDofs, 0.0);
//
//	const double bindingCell[] = { 1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0 };
//	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 + 4 * 16);
//	cadet::test::util::repeat(y.data() + 4 + 4 * 16, bindingCell, 8, 4 * 16);
//	cadet::test::util::populate(y.data() + 4 + 4 * 16 + 16 * 4 * (4 + 4), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 16);
//
//	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);
//
//	cadet::test::column::testConsistentInitializationSensitivity("GENERAL_RATE_MODEL", "DG", y.data(), yDot.data(), false, 1e-9);
//}

TEST_CASE("Column_Model_1D inlet DOF Jacobian", "[Column_1D],[DG],[DG1D],[UnitOp],[Jacobian],[Inlet],[todo_CI]")
{
	cadet::test::column::testInletDofJacobian("GENERAL_RATE_MODEL", "DG");
}

TEST_CASE("Column_Model_1D transport Jacobian", "[Column_1D],[DG],[DG1D],[UnitOp],[Jacobian],[todo_CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "GENERAL_RATE_MODEL", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_Model_1D with two component linear binding Jacobian", "[Column_1D],[DG],[DG1D],[UnitOp],[Jacobian],[todo_CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("GENERAL_RATE_MODEL", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_Model_1D LWE one vs two identical particle types match", "[Column_1D],[DG],[DG1D],[Simulation],[ParticleType],[todo_CI]")
{
	cadet::test::particle::testOneVsTwoIdenticalParticleTypes("GENERAL_RATE_MODEL", "DG", 2e-8, 5e-5);
}

TEST_CASE("Column_Model_1D LWE separate identical particle types match", "[Column_1D],[DG],[DG1D],[Simulation],[ParticleType],[todo_CI]")
{
	cadet::test::particle::testSeparateIdenticalParticleTypes("GENERAL_RATE_MODEL", "DG", 2e-8, 5e-5);
}

TEST_CASE("Column_Model_1D linear binding single particle matches particle distribution", "[Column_1D],[DG],[DG1D],[Simulation],[ParticleType],[todo_CI]")
{
	cadet::test::particle::testLinearMixedParticleTypes("GENERAL_RATE_MODEL", "DG", 5e-8, 5e-5);
}

TEST_CASE("Column_Model_1D multiple particle types Jacobian analytic vs AD", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ParticleType],[todo_CI]")
{
	cadet::test::particle::testJacobianMixedParticleTypes("GENERAL_RATE_MODEL", "DG", 1e-14);
}

TEST_CASE("Column_Model_1D multiple particle types time derivative Jacobian vs FD", "[Column_1D],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[ParticleType],[todo_CI]")
{
	cadet::test::particle::testTimeDerivativeJacobianMixedParticleTypesFD("GENERAL_RATE_MODEL", "DG", 1e-6, 0.0, 9e-4);
}

TEST_CASE("Column_Model_1D multiple spatially dependent particle types Jacobian analytic vs AD", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ParticleType],[todo_CI]")
{
	cadet::test::particle::testJacobianSpatiallyMixedParticleTypes("GENERAL_RATE_MODEL", "DG", 1e-11);
}

TEST_CASE("Column_Model_1D linear binding single particle matches spatially dependent particle distribution", "[Column_1D],[DG],[DG1D],[Simulation],[ParticleType],[todo_CI]")
{
	cadet::test::particle::testLinearSpatiallyMixedParticleTypes("GENERAL_RATE_MODEL", "DG", 5e-8, 5e-5);
}

TEST_CASE("Column_Model_1D dynamic reactions Jacobian vs AD bulk", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[todo_CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("GENERAL_RATE_MODEL", "DG", true, false, false, 1e-14);
}

TEST_CASE("Column_Model_1D dynamic reactions Jacobian vs AD particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[todo_CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("GENERAL_RATE_MODEL", "DG", false, true, false, 1e-14);
}

TEST_CASE("Column_Model_1D dynamic reactions Jacobian vs AD modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[todo_CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("GENERAL_RATE_MODEL", "DG", false, true, true, 1e-14);
}

TEST_CASE("Column_Model_1D dynamic reactions Jacobian vs AD bulk and particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[todo_CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("GENERAL_RATE_MODEL", "DG", true, true, false, 1e-14);
}

TEST_CASE("Column_Model_1D dynamic reactions Jacobian vs AD bulk and modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[todo_CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("GENERAL_RATE_MODEL", "DG", true, true, true, 1e-14);
}

TEST_CASE("Column_Model_1D dynamic reactions time derivative Jacobian vs FD bulk", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[todo_CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("GENERAL_RATE_MODEL", "DG", true, false, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_Model_1D dynamic reactions time derivative Jacobian vs FD particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[todo_CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("GENERAL_RATE_MODEL", "DG", false, true, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_Model_1D dynamic reactions time derivative Jacobian vs FD modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[todo_CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("GENERAL_RATE_MODEL", "DG", false, true, true, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_Model_1D dynamic reactions time derivative Jacobian vs FD bulk and particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[todo_CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("GENERAL_RATE_MODEL", "DG", true, true, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_Model_1D dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[todo_CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("GENERAL_RATE_MODEL", "DG", true, true, true, 1e-6, 1e-14, 9e-4);
}

inline cadet::JsonParameterProvider createColumnWithTwoCompLinearBindingThreeParticleTypes()
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("GENERAL_RATE_MODEL", "DG");

	const double parVolFrac[] = { 0.3, 0.6, 0.1 };
	const double parFactor[] = { 0.9, 0.8 };
	cadet::test::particle::extendModelToManyParticleTypes(jpp, 3, parFactor, parVolFrac);

	return jpp;
}

TEST_CASE("Column_Model_1D multi particle types dynamic reactions Jacobian vs AD bulk", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[todo_CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, false, false);
}

TEST_CASE("Column_Model_1D multi particle types dynamic reactions Jacobian vs AD particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[todo_CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, false);
}

TEST_CASE("Column_Model_1D multi particle types dynamic reactions Jacobian vs AD modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[todo_CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, true);
}

TEST_CASE("Column_Model_1D multi particle types dynamic reactions Jacobian vs AD bulk and particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[todo_CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, false);
}

TEST_CASE("Column_Model_1D multi particle types dynamic reactions Jacobian vs AD bulk and modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[todo_CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, true);
}

TEST_CASE("Column_Model_1D multi particle types dynamic reactions time derivative Jacobian vs FD bulk", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[todo_CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, false, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_Model_1D multi particle types dynamic reactions time derivative Jacobian vs FD particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[todo_CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_Model_1D multi particle types dynamic reactions time derivative Jacobian vs FD modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[todo_CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, true, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_Model_1D multi particle types dynamic reactions time derivative Jacobian vs FD bulk and particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[todo_CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_Model_1D multi particle types dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[todo_CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, true, 1e-6, 1e-14, 9e-4);
}
