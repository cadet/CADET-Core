// =============================================================================
//  CADET
//  
//  Copyright © 2008-2024: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
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

// todo add a meaningful backward flow test

// todo find analytical solution with linear binding

// todo find analytical solution without binding

// todo add (more) numerical reference (and EOC) tests
 
TEST_CASE("Radial LRMP numerical Benchmark with parameter sensitivities for linear case", "[RadLRMP],[Simulation],[Reference],[Sensitivity]") // todo CI flag: currently only runs locally but fails on server
{
	const std::string& modelFilePath = std::string("/data/model_radLRMP_dynLin_1comp_sensbenchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_radLRMP_dynLin_1comp_sensbenchmark1_FV_Z32.h5");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1e-6, 1e-6, 1e-6, 1e-6 };
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, 32, 0, true);
}

TEST_CASE("Radial LRMP transport Jacobian", "[RadLRMP],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "RADIAL_LUMPED_RATE_MODEL_WITH_PORES");
	cadet::test::column::testJacobianAD(jpp);
}

// NOTE: THE FOLLOWING TESTS ARE ONLY INCLUDED IN THE RELEASE CI, NOT THE STANDARD CI SINCE THEY ARE (TO A HIGH DEGREE) REDUNDANT WITH THE AXIAL FLOW TESTS

TEST_CASE("Radial LRMP Jacobian forward vs backward flow", "[RadLRMP],[UnitOp],[Residual],[Jacobian],[AD],[ReleaseCI]")
{
	cadet::test::column::testJacobianWenoForwardBackward("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", 0); // note: no weno for radial flow models atm
}

TEST_CASE("Radial LRMP time derivative Jacobian vs FD", "[RadLRMP],[UnitOp],[Residual],[Jacobian],[ReleaseCI]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("RADIAL_LUMPED_RATE_MODEL_WITH_PORES");
}

TEST_CASE("Radial LRMP flux Jacobian vs FD", "[RadLRMP],[UnitOp],[Residual],[Jacobian],[ReleaseCI]")
{
	cadet::test::column::testArrowHeadJacobianFD("RADIAL_LUMPED_RATE_MODEL_WITH_PORES");
}

TEST_CASE("Radial LRMP sensitivity Jacobians", "[RadLRMP],[UnitOp],[Sensitivity],[ReleaseCI]")
{
	cadet::test::column::testFwdSensJacobians("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", 1e-4, 6e-7);
}

//TEST_CASE("Radial LRMP forward sensitivity vs FD", "[RadLRMP],[Sensitivity],[Simulation],[failedFDtestLRMP]") // todo fix
//{
//	// todo comment on FD issue
//	// Relative error is checked first, we use high absolute error for letting
//	// some points that are far off pass the error test, too. This is required
//	// due to errors in finite differences.
//	const double fdStepSize[] = {1e-5, 1e-6, 1e-3, 1e-5};
//	const double absTols[] = {6e5, 2e-2, 2e-2, 1.0};
//	const double relTols[] = {5e-3, 1e-1, 5e-1, 6e-3};
//	const double passRatio[] = {0.87, 0.84, 0.88, 0.95};
//	cadet::test::column::testFwdSensSolutionFD("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", true, fdStepSize, absTols, relTols, passRatio);
//}

//TEST_CASE("Radial LRMP forward sensitivity forward vs backward flow", "[RadLRMP],[Sensitivity],[Simulation],[fixLRMP]") // todo fix
//{
//	// todo why is there a pass ratio when we compare fwd and bwd flow?
//	const double absTols[] = {50.0, 2e-10, 1.0, 5e-7};
//	const double relTols[] = {2e-4, 9e-6, 5e-7, 1e-7};
//	const double passRatio[] = {1.0, 0.99, 0.98, 0.99}; // todo? works for 0.79, 0.99, 0.98, 0.99
//	cadet::test::column::testFwdSensSolutionForwardBackward("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", absTols, relTols, passRatio);
//}

TEST_CASE("Radial LRMP consistent initialization with linear binding", "[RadLRMP],[ConsistentInit],[ReleaseCI]")
{
	cadet::test::column::testConsistentInitializationLinearBinding("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", 1e-12, 1e-12);
}

//TEST_CASE("Radial LRMP consistent initialization with SMA binding", "[RadLRMP],[ConsistentInit],[fixLRMP]")
//{
//	std::vector<double> y(4 + 4 * 16 + 16 * (4 + 4) + 4 * 16, 0.0);
//// Optimal values:
////	const double bindingCell[] = {1.2, 2.0, 1.0, 1.5, 858.034, 66.7896, 3.53273, 2.53153, 
////		1.0, 1.8, 1.5, 1.6, 856.173, 64.457, 5.73227, 2.85286};
//	const double bindingCell[] = {1.2, 2.0, 1.0, 1.5, 840.0, 63.0, 3.0, 3.0, 
//		1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0};
//	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 + 4 * 16);
//	cadet::test::util::repeat(y.data() + 4 + 4 * 16, bindingCell, 16, 8);
//	cadet::test::util::populate(y.data() + 4 + 4 * 16 + 16 * (4 + 4), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 16);
//
//	cadet::test::column::testConsistentInitializationSMABinding("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", y.data(), 1e-14, 1e-5);
//}

TEST_CASE("Radial LRMP consistent sensitivity initialization with linear binding", "[RadLRMP],[ConsistentInit],[Sensitivity],[ReleaseCI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 4 + 4 * 16 + 16 * (4 + 4) + 4 * 16;
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, numDofs);
	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	cadet::test::column::testConsistentInitializationSensitivity("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", y.data(), yDot.data(), true, 1e-14);
}

TEST_CASE("Radial LRMP consistent sensitivity initialization with SMA binding", "[RadLRMP],[ConsistentInit],[Sensitivity],[ReleaseCI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 4 + 4 * 16 + 16 * (4 + 4) + 4 * 16;
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);

	const double bindingCell[] = {1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0};
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 + 4 * 16);
	cadet::test::util::repeat(y.data() + 4 + 4 * 16, bindingCell, 8, 16);
	cadet::test::util::populate(y.data() + 4 + 4 * 16 + 16 * (4 + 4), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 16);

	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	cadet::test::column::testConsistentInitializationSensitivity("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", y.data(), yDot.data(), false, 1e-10);
}

TEST_CASE("Radial LRMP inlet DOF Jacobian", "[RadLRMP],[UnitOp],[Jacobian],[Inlet],[ReleaseCI]")
{
	cadet::test::column::testInletDofJacobian("RADIAL_LUMPED_RATE_MODEL_WITH_PORES");
}

TEST_CASE("Radial LRMP with two component linear binding Jacobian", "[RadLRMP],[UnitOp],[Jacobian],[ReleaseCI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("RADIAL_LUMPED_RATE_MODEL_WITH_PORES");
	cadet::test::column::testJacobianAD(jpp);
}

TEST_CASE("Radial LRMP LWE one vs two identical particle types match", "[RadLRMP],[Simulation],[ParticleType],[ReleaseCI]")
{
	cadet::test::particle::testOneVsTwoIdenticalParticleTypes("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", 2.2e-8, 6e-5);
}

TEST_CASE("Radial LRMP LWE separate identical particle types match", "[RadLRMP],[Simulation],[ParticleType],[ReleaseCI]")
{
	cadet::test::particle::testSeparateIdenticalParticleTypes("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", 1e-15, 1e-15);
}

TEST_CASE("Radial LRMP linear binding single particle matches particle distribution", "[RadLRMP],[Simulation],[ParticleType],[ReleaseCI]")
{
	cadet::test::particle::testLinearMixedParticleTypes("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", 5e-8, 5e-5);
}

TEST_CASE("Radial LRMP multiple particle types Jacobian analytic vs AD", "[RadLRMP],[Jacobian],[AD],[ParticleType],[ReleaseCI]")
{
	cadet::test::particle::testJacobianMixedParticleTypes("RADIAL_LUMPED_RATE_MODEL_WITH_PORES");
}

TEST_CASE("Radial LRMP multiple particle types time derivative Jacobian vs FD", "[RadLRMP],[UnitOp],[Residual],[Jacobian],[ParticleType],[ReleaseCI],[FD]")
{
	cadet::test::particle::testTimeDerivativeJacobianMixedParticleTypesFD("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", 1e-6, 0.0, 9e-4);
}

TEST_CASE("Radial LRMP multiple spatially dependent particle types Jacobian analytic vs AD", "[RadLRMP],[Jacobian],[AD],[ParticleType],[ReleaseCI]")
{
	cadet::test::particle::testJacobianSpatiallyMixedParticleTypes("RADIAL_LUMPED_RATE_MODEL_WITH_PORES");
}

TEST_CASE("Radial LRMP linear binding single particle matches spatially dependent particle distribution", "[RadLRMP],[Simulation],[ParticleType],[ReleaseCI]")
{
	cadet::test::particle::testLinearSpatiallyMixedParticleTypes("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", 5e-8, 5e-5);
}

TEST_CASE("Radial LRMP multiple spatially dependent particle types flux Jacobian vs FD", "[RadLRMP],[UnitOp],[Residual],[Jacobian],[ParticleType],[ReleaseCI],[FD]")
{
	cadet::test::particle::testArrowHeadJacobianSpatiallyMixedParticleTypes("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", 1e-6, 1e-8, 1e-5);
}

TEST_CASE("Radial LRMP dynamic reactions Jacobian vs AD bulk", "[RadLRMP],[Jacobian],[AD],[ReactionModel],[ReleaseCI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", true, false, false);
}

TEST_CASE("Radial LRMP dynamic reactions Jacobian vs AD particle", "[RadLRMP],[Jacobian],[AD],[ReactionModel],[ReleaseCI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", false, true, false);
}

TEST_CASE("Radial LRMP dynamic reactions Jacobian vs AD modified particle", "[RadLRMP],[Jacobian],[AD],[ReactionModel],[ReleaseCI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", false, true, true);
}

TEST_CASE("Radial LRMP dynamic reactions Jacobian vs AD bulk and particle", "[RadLRMP],[Jacobian],[AD],[ReactionModel],[ReleaseCI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", true, true, false);
}

TEST_CASE("Radial LRMP dynamic reactions Jacobian vs AD bulk and modified particle", "[RadLRMP],[Jacobian],[AD],[ReactionModel],[ReleaseCI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", true, true, true);
}

TEST_CASE("Radial LRMP dynamic reactions time derivative Jacobian vs FD bulk", "[RadLRMP],[Jacobian],[Residual],[ReactionModel],[ReleaseCI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Radial LRMP dynamic reactions time derivative Jacobian vs FD particle", "[RadLRMP],[Jacobian],[Residual],[ReactionModel],[ReleaseCI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", false, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Radial LRMP dynamic reactions time derivative Jacobian vs FD modified particle", "[RadLRMP],[Jacobian],[Residual],[ReactionModel],[ReleaseCI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", false, true, true, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Radial LRMP dynamic reactions time derivative Jacobian vs FD bulk and particle", "[RadLRMP],[Jacobian],[Residual],[ReactionModel],[ReleaseCI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", true, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Radial LRMP dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[RadLRMP],[Jacobian],[Residual],[ReactionModel],[ReleaseCI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("RADIAL_LUMPED_RATE_MODEL_WITH_PORES", true, true, true, 1e-6, 1e-14, 8e-4);
}

inline cadet::JsonParameterProvider createColumnWithTwoCompLinearBindingThreeParticleTypes()
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("RADIAL_LUMPED_RATE_MODEL_WITH_PORES");

	const double parVolFrac[] = {0.3, 0.6, 0.1};
	const double parFactor[] = {0.9, 0.8};
	cadet::test::particle::extendModelToManyParticleTypes(jpp, 3, parFactor, parVolFrac);

	return jpp;
}

TEST_CASE("Radial LRMP multi particle types dynamic reactions Jacobian vs AD bulk", "[RadLRMP],[Jacobian],[AD],[ReactionModel],[ParticleType],[ReleaseCI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, false, false);
}

TEST_CASE("Radial LRMP multi particle types dynamic reactions Jacobian vs AD particle", "[RadLRMP],[Jacobian],[AD],[ReactionModel],[ParticleType],[ReleaseCI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, false);
}

TEST_CASE("Radial LRMP multi particle types dynamic reactions Jacobian vs AD modified particle", "[RadLRMP],[Jacobian],[AD],[ReactionModel],[ParticleType],[ReleaseCI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, true);
}

TEST_CASE("Radial LRMP multi particle types dynamic reactions Jacobian vs AD bulk and particle", "[RadLRMP],[Jacobian],[AD],[ReactionModel],[ParticleType],[ReleaseCI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, false);
}

TEST_CASE("Radial LRMP multi particle types dynamic reactions Jacobian vs AD bulk and modified particle", "[RadLRMP],[Jacobian],[AD],[ReactionModel],[ParticleType],[ReleaseCI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, true);
}

TEST_CASE("Radial LRMP multi particle types dynamic reactions time derivative Jacobian vs FD bulk", "[RadLRMP],[Jacobian],[Residual],[ReactionModel],[ParticleType],[ReleaseCI],[FD]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Radial LRMP multi particle types dynamic reactions time derivative Jacobian vs FD particle", "[RadLRMP],[Jacobian],[Residual],[ReactionModel],[ParticleType],[ReleaseCI],[FD]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Radial LRMP multi particle types dynamic reactions time derivative Jacobian vs FD modified particle", "[RadLRMP],[Jacobian],[Residual],[ReactionModel],[ParticleType],[ReleaseCI],[FD]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, true, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Radial LRMP multi particle types dynamic reactions time derivative Jacobian vs FD bulk and particle", "[RadLRMP],[Jacobian],[Residual],[ReactionModel],[ParticleType],[ReleaseCI],[FD]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Radial LRMP multi particle types dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[RadLRMP],[Jacobian],[Residual],[ReactionModel],[ParticleType],[ReleaseCI],[FD]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, true, 1e-6, 1e-14, 8e-4);
}
