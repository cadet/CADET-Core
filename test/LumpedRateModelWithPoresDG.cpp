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
#include "Utils.hpp"

TEST_CASE("LRMP_DG LWE forward vs backward flow", "[LRMP],[DG],[DG1D],[Simulation],[CI]")
{
	cadet::test::column::DGparams disc;

	// Test all integration modes
	for (int i = 0; i <= 1; i++)
	{
		disc.setIntegrationMode(i);
		cadet::test::column::testForwardBackward("LUMPED_RATE_MODEL_WITH_PORES", disc, 6e-8, 4e-6);
	}
}

TEST_CASE("LRMP_DG linear pulse vs analytic solution", "[LRMP],[DG],[DG1D],[Simulation],[Analytic],[CI]")
{
	cadet::test::column::DGparams disc;
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITH_PORES", "/data/lrmp-pulseBenchmark.data", true, true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITH_PORES", "/data/lrmp-pulseBenchmark.data", true, false, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITH_PORES", "/data/lrmp-pulseBenchmark.data", false, true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITH_PORES", "/data/lrmp-pulseBenchmark.data", false, false, disc, 6e-5, 1e-7);
}

TEST_CASE("LRMP_DG non-binding linear pulse vs analytic solution", "[LRMP],[DG],[DG1D],[Simulation],[Analytic],[NonBinding],[CI]")
{
	cadet::test::column::DGparams disc;
	cadet::test::column::testAnalyticNonBindingBenchmark("LUMPED_RATE_MODEL_WITH_PORES", "/data/lrmp-nonBinding.data", true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticNonBindingBenchmark("LUMPED_RATE_MODEL_WITH_PORES", "/data/lrmp-nonBinding.data", false, disc, 6e-5, 1e-7);
}

//TEST_CASE("LRMP_DG Jacobian forward vs backward flow", "[LRMP],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[AD],[fix]")
//{
//	cadet::test::column::DGparams disc;
//
//	// Test all integration modes
//	for (int i = 0; i <= 1; i++)
//	{
//		disc.setIntegrationMode(i);
//		cadet::test::column::testJacobianForwardBackward("LUMPED_RATE_MODEL_WITH_PORES", disc, std::numeric_limits<float>::epsilon() * 100.0);
//	}
//}

TEST_CASE("LRMP_DG numerical Benchmark with parameter sensitivities for linear case", "[LRMP],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_LRMP_dynLin_1comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_LRMP_dynLin_1comp_sensbenchmark1_DG_P3Z8.h5");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1.0, 1.0, 1.0, 1.0 };

	cadet::test::column::DGparams disc(0, 3, 8);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("LRMP_DG numerical Benchmark with parameter sensitivities for SMA LWE case", "[LRMP],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_LRMP_reqSMA_4comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_LRMP_reqSMA_4comp_sensbenchmark1_DG_P3Z8.h5");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1.0, 1.0, 1.0, 1.0 };

	cadet::test::column::DGparams disc(0, 3, 8);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "000", absTol, relTol, disc, true);
}

TEST_CASE("LRMP_DG time derivative Jacobian vs FD", "[LRMP],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[CI]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("LUMPED_RATE_MODEL_WITH_PORES", "DG");
}

TEST_CASE("LRMP_DG sensitivity Jacobians", "[LRMP],[DG],[DG1D],[UnitOp],[Sensitivity],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("LUMPED_RATE_MODEL_WITH_PORES", "DG");

	cadet::test::column::testFwdSensJacobians(jpp, 1e-4, 6e-7);
}
 
//// todo fix: not just adjust tolerances as in FV but theres an actual error here: access violation in densematrix
//TEST_CASE("LRMP_DG forward sensitivity vs FD", "[LRMP],[DG],[DG1D],[Sensitivity],[Simulation],[todo]")
//{
//	// Relative error is checked first, we use high absolute error for letting
//	// some points that are far off pass the error test, too. This is required
//	// due to errors in finite differences.
//	const double fdStepSize[] = { 1e-5, 1e-6, 1e-3, 1e-5 };
//	const double absTols[] = { 6e5, 2e-2, 2e-2, 1.0 };
//	const double relTols[] = { 5e-3, 1e-1, 5e-1, 6e-3 };
//	const double passRatio[] = { 0.87, 0.84, 0.88, 0.95 };
//	cadet::test::column::testFwdSensSolutionFD("LUMPED_RATE_MODEL_WITH_PORES", "DG", true, fdStepSize, absTols, relTols, passRatio);
//}
 
// todo fix: not just adjust tolerances as in FV but theres an actual error here: access violation in densematrix
//TEST_CASE("LRMP_DG forward sensitivity forward vs backward flow", "[LRMP],[DG],[DG1D],[Sensitivity],[Simulation],[todo]")
//{
//	const double absTols[] = { 50.0, 2e-10, 1.0, 5e-7 };
//	const double relTols[] = { 2e-4, 9e-6, 5e-7, 1e-7 };
//	const double passRatio[] = { 1.0, 0.99, 0.98, 0.99 };
//	cadet::test::column::testFwdSensSolutionForwardBackward("LUMPED_RATE_MODEL_WITH_PORES", "DG", absTols, relTols, passRatio);
//}

// todo fix consistent initialization for AD with req binding
TEST_CASE("LRMP_DG consistent initialization with linear binding", "[LRMP],[DG],[DG1D],[ConsistentInit],[CI]")
{
	cadet::test::column::testConsistentInitializationLinearBinding("LUMPED_RATE_MODEL_WITH_PORES", "DG", 1e-12, 1e-12, 0, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("LUMPED_RATE_MODEL_WITH_PORES", "DG", 1e-12, 1e-12, 1, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("LUMPED_RATE_MODEL_WITH_PORES", "DG", 1e-12, 1e-12, 0, 1);
	//cadet::test::column::testConsistentInitializationLinearBinding("LUMPED_RATE_MODEL_WITH_PORES", "DG", 1e-12, 1e-12, 1, 1);
}

////// todo fix consistent initialization for SMA (initialization not completely correct; AD gives assertion error)
//TEST_CASE("LRMP_DG consistent initialization with SMA binding", "[LRMP],[DG],[DG1D],[ConsistentInit],[todo]")
//{
//	std::vector<double> y(4 + 4 * 16 + 16 * (4 + 4) + 4 * 16, 0.0);
//	const double bindingCell[] = { 1.2, 2.0, 1.0, 1.5, 840.0, 63.0, 3.0, 3.0,
//		1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0 };
//	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 + 4 * 16);
//	cadet::test::util::repeat(y.data() + 4 + 4 * 16, bindingCell, 16, 8);
//	cadet::test::util::populate(y.data() + 4 + 4 * 16 + 16 * (4 + 4), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 16);
//
//	//cadet::test::column::testConsistentInitializationSMABinding("LUMPED_RATE_MODEL_WITH_PORES", "DG", y.data(), 1e-14, 1e-5, 0, 0);
//	//cadet::test::column::testConsistentInitializationSMABinding("LUMPED_RATE_MODEL_WITH_PORES", "DG", y.data(), 1e-14, 1e-5, 1, 0);
//	//cadet::test::column::testConsistentInitializationSMABinding("LUMPED_RATE_MODEL_WITH_PORES", "DG", y.data(), 1e-14, 1e-5, 0, 1);
//	//cadet::test::column::testConsistentInitializationSMABinding("LUMPED_RATE_MODEL_WITH_PORES", "DG", y.data(), 1e-14, 1e-5, 1, 1);
//}

// todo fix kinetic binding sensitivity init
TEST_CASE("LRMP_DG consistent sensitivity initialization with linear binding", "[LRMP],[DG],[DG1D],[ConsistentInit],[Sensitivity],[CI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 2 + 2 * 10 + 10 * (2 + 2);
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, numDofs);
	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	//cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITH_PORES", "DG", y.data(), yDot.data(), true, 1e-14, 0, 0); // doesnt work
	cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITH_PORES", "DG", y.data(), yDot.data(), true, 1e-14, 1, 0); // works
	//cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITH_PORES", "DG", y.data(), yDot.data(), true, 1e-14, 0, 1); // doesnt work
	cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITH_PORES", "DG", y.data(), yDot.data(), true, 1e-14, 1, 1); // works
}

//// todo fix memory stuff (works for FV) 
//TEST_CASE("LRMP_DG consistent sensitivity initialization with SMA binding", "[LRMP],[DG],[DG1D],[ConsistentInit],[Sensitivity],[todo]")
//{
//	// Fill state vector with given initial values
//  const unsigned int numDofs = 4 + 4 * 10 + 10 * (4 + 4);
//	std::vector<double> y(numDofs, 0.0);
//	std::vector<double> yDot(numDofs, 0.0);
//
//	const double bindingCell[] = { 1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0 };
//	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 + 4 * 10);
//	cadet::test::util::repeat(y.data() + 4 + 4 * 10, bindingCell, 8, 10);
//	cadet::test::util::populate(y.data() + 4 + 4 * 10 + 10 * (4 + 4), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 10);
//
//	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);
//
//	//cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITH_PORES", "DG", y.data(), yDot.data(), false, 1e-10, 0, 0);
//	//cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITH_PORES", "DG", y.data(), yDot.data(), false, 1e-10, 1, 0);
//	//cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITH_PORES", "DG", y.data(), yDot.data(), false, 1e-10, 0, 1);
//	//cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITH_PORES", "DG", y.data(), yDot.data(), false, 1e-10, 1, 1);
//}

TEST_CASE("LRMP_DG inlet DOF Jacobian", "[LRMP],[DG],[DG1D],[UnitOp],[Jacobian],[Inlet],[CI]")
{
	cadet::test::column::testInletDofJacobian("LUMPED_RATE_MODEL_WITH_PORES", "DG");
}

TEST_CASE("LRMP_DG transport Jacobian", "[LRMP],[DG],[DG1D],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "LUMPED_RATE_MODEL_WITH_PORES", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("LRMP_DG with two component linear binding Jacobian", "[LRMP],[DG],[DG1D],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("LUMPED_RATE_MODEL_WITH_PORES", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("LRMP_DG LWE one vs two identical particle types match", "[LRMP],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testOneVsTwoIdenticalParticleTypes("LUMPED_RATE_MODEL_WITH_PORES", "DG", 2.2e-8, 6e-5);
}

TEST_CASE("LRMP_DG LWE separate identical particle types match", "[LRMP],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testSeparateIdenticalParticleTypes("LUMPED_RATE_MODEL_WITH_PORES", "DG", 1e-12, 1e-12);
}

TEST_CASE("LRMP_DG linear binding single particle matches particle distribution", "[LRMP],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearMixedParticleTypes("LUMPED_RATE_MODEL_WITH_PORES", "DG", 5e-8, 5e-5);
}

TEST_CASE("LRMP_DG multiple particle types Jacobian analytic vs AD", "[LRMP],[DG],[DG1D],[Jacobian],[AD],[ParticleType],[CI]")
{
	cadet::test::particle::testJacobianMixedParticleTypes("LUMPED_RATE_MODEL_WITH_PORES", "DG", 1e-11);
}

TEST_CASE("LRMP_DG multiple particle types time derivative Jacobian vs FD", "[LRMP],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[ParticleType],[CI]")
{
	cadet::test::particle::testTimeDerivativeJacobianMixedParticleTypesFD("LUMPED_RATE_MODEL_WITH_PORES", "DG", 1e-6, 0.0, 9e-4);
}

TEST_CASE("LRMP_DG multiple spatially dependent particle types Jacobian analytic vs AD", "[LRMP],[DG],[DG1D],[Jacobian],[AD],[ParticleType],[CI]")
{
	cadet::test::particle::testJacobianSpatiallyMixedParticleTypes("LUMPED_RATE_MODEL_WITH_PORES", "DG", std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("LRMP_DG linear binding single particle matches spatially dependent particle distribution", "[LRMP],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearSpatiallyMixedParticleTypes("LUMPED_RATE_MODEL_WITH_PORES", "DG", 5e-8, 5e-5);
}

TEST_CASE("LRMP_DG dynamic reactions Jacobian vs AD bulk", "[LRMP],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("LUMPED_RATE_MODEL_WITH_PORES", "DG", true, false, false, std::numeric_limits<float>::epsilon() * 100.0);
}
TEST_CASE("LRMP_DG dynamic reactions Jacobian vs AD particle", "[LRMP],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("LUMPED_RATE_MODEL_WITH_PORES", "DG", false, true, false, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("LRMP_DG dynamic reactions Jacobian vs AD modified particle", "[LRMP],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("LUMPED_RATE_MODEL_WITH_PORES", "DG", false, true, true, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("LRMP_DG dynamic reactions Jacobian vs AD bulk and particle", "[LRMP],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("LUMPED_RATE_MODEL_WITH_PORES", "DG", true, true, false, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("LRMP_DG dynamic reactions Jacobian vs AD bulk and modified particle", "[LRMP],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("LUMPED_RATE_MODEL_WITH_PORES", "DG", true, true, true, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("LRMP_DG dynamic reactions time derivative Jacobian vs FD bulk", "[LRMP],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("LUMPED_RATE_MODEL_WITH_PORES", "DG", true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("LRMP_DG dynamic reactions time derivative Jacobian vs FD particle", "[LRMP],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("LUMPED_RATE_MODEL_WITH_PORES", "DG", false, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("LRMP_DG dynamic reactions time derivative Jacobian vs FD modified particle", "[LRMP],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("LUMPED_RATE_MODEL_WITH_PORES", "DG", false, true, true, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("LRMP_DG dynamic reactions time derivative Jacobian vs FD bulk and particle", "[LRMP],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("LUMPED_RATE_MODEL_WITH_PORES", "DG", true, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("LRMP_DG dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[LRMP],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("LUMPED_RATE_MODEL_WITH_PORES", "DG", true, true, true, 1e-6, 1e-14, 8e-4);
}

inline cadet::JsonParameterProvider createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes()
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("LUMPED_RATE_MODEL_WITH_PORES", "DG");

	const double parVolFrac[] = { 0.3, 0.6, 0.1 };
	const double parFactor[] = { 0.9, 0.8 };
	cadet::test::particle::extendModelToManyParticleTypes(jpp, 3, parFactor, parVolFrac);

	return jpp;
}

TEST_CASE("LRMP_DG multi particle types dynamic reactions Jacobian vs AD bulk", "[LRMP],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, false, false, 1e-11);
}

TEST_CASE("LRMP_DG multi particle types dynamic reactions Jacobian vs AD particle", "[LRMP],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, false, 1e-11);
}

TEST_CASE("LRMP_DG multi particle types dynamic reactions Jacobian vs AD modified particle", "[LRMP],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, true, 1e-11);
}

TEST_CASE("LRMP_DG multi particle types dynamic reactions Jacobian vs AD bulk and particle", "[LRMP],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, false, 1e-11);
}

TEST_CASE("LRMP_DG multi particle types dynamic reactions Jacobian vs AD bulk and modified particle", "[LRMP],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, true, 1e-11);
}

TEST_CASE("LRMP_DG multi particle types dynamic reactions time derivative Jacobian vs FD bulk", "[LRMP],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("LRMP_DG multi particle types dynamic reactions time derivative Jacobian vs FD particle", "[LRMP],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("LRMP_DG multi particle types dynamic reactions time derivative Jacobian vs FD modified particle", "[LRMP],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, true, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("LRMP_DG multi particle types dynamic reactions time derivative Jacobian vs FD bulk and particle", "[LRMP],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("LRMP_DG multi particle types dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[LRMP],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, true, 1e-6, 1e-14, 8e-4);
}
