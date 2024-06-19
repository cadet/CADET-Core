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

TEST_CASE("Radial GRM numerical Benchmark with parameter sensitivities for linear case", "[RadGRM],[Simulation],[Reference],[Sensitivity]") // todo CI flag: currently only runs locally but fails on server
{
	const std::string& modelFilePath = std::string("/data/model_radGRM_dynLin_1comp_sensbenchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_radGRM_dynLin_1comp_sensbenchmark1_FV_Z32parZ4.h5");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1e-6, 1e-6, 1e-6, 1e-6 };
	cadet::test::column::FVparams disc(32, 4);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("Radial GRM transport Jacobian", "[RadGRM],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "RADIAL_GENERAL_RATE_MODEL", "FV");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

// NOTE: THE FOLLOWING TESTS ARE ONLY INCLUDED IN THE RELEASE CI, NOT THE STANDARD CI SINCE THEY ARE (TO A HIGH DEGREE) REDUNDANT WITH THE AXIAL FLOW TESTS

TEST_CASE("Radial GRM Jacobian forward vs backward flow", "[RadGRM],[UnitOp],[Residual],[Jacobian],[AD],[ReleaseCI]")
{
	cadet::test::column::FVparams disc(16);
	cadet::test::column::testJacobianForwardBackward("RADIAL_GENERAL_RATE_MODEL", disc);
}

TEST_CASE("Radial GRM time derivative Jacobian vs FD", "[RadGRM],[UnitOp],[Residual],[Jacobian],[ReleaseCI],[FD]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("RADIAL_GENERAL_RATE_MODEL", "FV", 1e-6, 0.0, 9e-4);
}

TEST_CASE("Radial GRM rapid-equilibrium binding flux Jacobian vs FD", "[RadGRM],[UnitOp],[Residual],[Jacobian],[ReleaseCI],[FD]")
{
	cadet::test::column::testArrowHeadJacobianFD("RADIAL_GENERAL_RATE_MODEL", false, 1e-6, 2e-9);
}

TEST_CASE("Radial GRM rapid-equilibrium binding with surf diff par dep flux Jacobian vs FD", "[RadGRM],[UnitOp],[Residual],[Jacobian],[ParameterDependence],[ReleaseCI],[FD]")
{
	cadet::test::column::testArrowHeadJacobianFDVariableParSurfDiff("RADIAL_GENERAL_RATE_MODEL", 1e-6, 5e-9);
}

TEST_CASE("Radial GRM dynamic binding with surf diff par dep Jacobian vs AD", "[RadGRM],[UnitOp],[Residual],[Jacobian],[ParameterDependence],[ReleaseCI]")
{
	cadet::test::column::testJacobianADVariableParSurfDiff("RADIAL_GENERAL_RATE_MODEL", "FV", true);
}

TEST_CASE("Radial GRM rapid-equilibrium binding with surf diff par dep Jacobian vs AD", "[RadGRM],[UnitOp],[Residual],[Jacobian],[ParameterDependence],[ReleaseCI]")
{
	cadet::test::column::testJacobianADVariableParSurfDiff("RADIAL_GENERAL_RATE_MODEL", "FV", false);
}

TEST_CASE("Radial GRM dynamic binding flux Jacobian vs FD", "[RadGRM],[UnitOp],[Residual],[Jacobian],[ReleaseCI],[FD]")
{
	cadet::test::column::testArrowHeadJacobianFD("RADIAL_GENERAL_RATE_MODEL", true, 1e-6, 2e-9);
}

TEST_CASE("Radial GRM sensitivity Jacobians", "[RadGRM],[UnitOp],[Sensitivity],[ReleaseCI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("RADIAL_GENERAL_RATE_MODEL", "FV");

	cadet::test::column::testFwdSensJacobians(jpp, 1e-4, 6e-7);
}

//TEST_CASE("Radial GRM forward sensitivity vs FD", "[RadGRM],[Sensitivity],[Simulation],[failedFDtestGRM]") // todo fix. fails for SMA_KA for both binding modes
//{
//	// todo comment
//	// Relative error is checked first, we use high absolute error for letting
//	// some points that are far off pass the error test, too. This is required
//	// due to errors in finite differences.
//	const double fdStepSize[] = {5e-5, 1e-4, 1e-4, 1e-3};
//	const double absTols[] = {3e5, 2e-3, 2e-4, 5.0};
//	const double relTols[] = {5e-3, 7e-2, 8e-2, 1e-4};
//	const double passRatio[] = {0.95, 0.9, 0.91, 0.83};
//	cadet::test::column::testFwdSensSolutionFD("RADIAL_GENERAL_RATE_MODEL", false, fdStepSize, absTols, relTols, passRatio);
//}
//
//TEST_CASE("Radial GRM forward sensitivity forward vs backward flow", "[RadGRM],[Sensitivity],[Simulation],[fixGRM]") // todo fix. fails for COL_DISPERION for both binding modes
//{
//	const double absTols[] = {4e-5, 1e-11, 1e-11, 8e-9};
//	const double relTols[] = {6e-9, 5e-8, 5e-6, 5e-10};
//	const double passRatio[] = {0.99, 0.95, 0.98, 0.98};
//	cadet::test::column::testFwdSensSolutionForwardBackward("RADIAL_GENERAL_RATE_MODEL", absTols, relTols, passRatio);
//}

TEST_CASE("Radial GRM consistent initialization with linear binding", "[RadGRM],[ConsistentInit],[ReleaseCI]")
{
	cadet::test::column::testConsistentInitializationLinearBinding("RADIAL_GENERAL_RATE_MODEL", "FV", 1e-12, 1e-12);
}

//TEST_CASE("Radial GRM consistent initialization with SMA binding", "[RadGRM],[ConsistentInit],[fixGRM]") // todo fix
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
//	cadet::test::column::testConsistentInitializationSMABinding("RADIAL_GENERAL_RATE_MODEL", y.data(), 1e-14, 1e-5);
//}

TEST_CASE("Radial GRM consistent sensitivity initialization with linear binding", "[RadGRM],[ConsistentInit],[Sensitivity],[ReleaseCI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 4 + 4 * 16 + 16 * 4 * (4 + 4) + 4 * 16;
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, numDofs);
	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	cadet::test::column::testConsistentInitializationSensitivity("RADIAL_GENERAL_RATE_MODEL", "FV", y.data(), yDot.data(), true, 1e-14);
}

TEST_CASE("Radial GRM consistent sensitivity initialization with SMA binding", "[RadGRM],[ConsistentInit],[Sensitivity],[ReleaseCI]")
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

	cadet::test::column::testConsistentInitializationSensitivity("RADIAL_GENERAL_RATE_MODEL", "FV", y.data(), yDot.data(), false, 1e-9);
}

TEST_CASE("Radial GRM inlet DOF Jacobian", "[RadGRM],[UnitOp],[Jacobian],[Inlet],[ReleaseCI]")
{
	cadet::test::column::testInletDofJacobian("RADIAL_GENERAL_RATE_MODEL", "FV");
}

TEST_CASE("Radial GRM with two component linear binding Jacobian", "[RadGRM],[UnitOp],[Jacobian],[ReleaseCI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("RADIAL_GENERAL_RATE_MODEL", "FV");
	cadet::test::column::testJacobianAD(jpp);
}

TEST_CASE("Radial GRM LWE one vs two identical particle types match", "[RadGRM],[Simulation],[ParticleType],[ReleaseCI]")
{
	cadet::test::particle::testOneVsTwoIdenticalParticleTypes("RADIAL_GENERAL_RATE_MODEL", "FV", 2e-8, 5e-5);
}

TEST_CASE("Radial GRM LWE separate identical particle types match", "[RadGRM],[Simulation],[ParticleType],[ReleaseCI]")
{
	cadet::test::particle::testSeparateIdenticalParticleTypes("RADIAL_GENERAL_RATE_MODEL", "FV", 1e-15, 1e-15);
}

TEST_CASE("Radial GRM linear binding single particle matches particle distribution", "[RadGRM],[Simulation],[ParticleType],[ReleaseCI]")
{
	cadet::test::particle::testLinearMixedParticleTypes("RADIAL_GENERAL_RATE_MODEL", "FV", 5e-8, 5e-5);
}

TEST_CASE("Radial GRM multiple particle types Jacobian analytic vs AD", "[RadGRM],[Jacobian],[AD],[ParticleType],[ReleaseCI]")
{
	cadet::test::particle::testJacobianMixedParticleTypes("RADIAL_GENERAL_RATE_MODEL", "FV");
}

TEST_CASE("Radial GRM multiple particle types time derivative Jacobian vs FD", "[RadGRM],[UnitOp],[Residual],[Jacobian],[ParticleType],[ReleaseCI]")
{
	cadet::test::particle::testTimeDerivativeJacobianMixedParticleTypesFD("RADIAL_GENERAL_RATE_MODEL", "FV", 1e-6, 0.0, 9e-4);
}

TEST_CASE("Radial GRM multiple spatially dependent particle types Jacobian analytic vs AD", "[RadGRM],[Jacobian],[AD],[ParticleType],[ReleaseCI]")
{
	cadet::test::particle::testJacobianSpatiallyMixedParticleTypes("RADIAL_GENERAL_RATE_MODEL", "FV");
}

TEST_CASE("Radial GRM linear binding single particle matches spatially dependent particle distribution", "[RadGRM],[Simulation],[ParticleType],[ReleaseCI]")
{
	cadet::test::particle::testLinearSpatiallyMixedParticleTypes("RADIAL_GENERAL_RATE_MODEL", "FV", 5e-8, 5e-5);
}

TEST_CASE("Radial GRM multiple spatially dependent particle types flux Jacobian vs FD", "[RadGRM],[UnitOp],[Residual],[Jacobian],[ParticleType],[ReleaseCI],[FD]")
{
	cadet::test::particle::testArrowHeadJacobianSpatiallyMixedParticleTypes("RADIAL_GENERAL_RATE_MODEL", 1e-6, 1e-8, 1e-5);
}

TEST_CASE("Radial GRM dynamic reactions Jacobian vs AD bulk", "[RadGRM],[Jacobian],[AD],[ReactionModel],[ReleaseCI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_GENERAL_RATE_MODEL", "FV", true, false, false);
}

TEST_CASE("Radial GRM dynamic reactions Jacobian vs AD particle", "[RadGRM],[Jacobian],[AD],[ReactionModel],[ReleaseCI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_GENERAL_RATE_MODEL", "FV", false, true, false);
}

TEST_CASE("Radial GRM dynamic reactions Jacobian vs AD modified particle", "[RadGRM],[Jacobian],[AD],[ReactionModel],[ReleaseCI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_GENERAL_RATE_MODEL", "FV", false, true, true);
}

TEST_CASE("Radial GRM dynamic reactions Jacobian vs AD bulk and particle", "[RadGRM],[Jacobian],[AD],[ReactionModel],[ReleaseCI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_GENERAL_RATE_MODEL", "FV", true, true, false);
}

TEST_CASE("Radial GRM dynamic reactions Jacobian vs AD bulk and modified particle", "[RadGRM],[Jacobian],[AD],[ReactionModel],[ReleaseCI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_GENERAL_RATE_MODEL", "FV", true, true, true);
}

TEST_CASE("Radial GRM dynamic reactions time derivative Jacobian vs FD bulk", "[RadGRM],[Jacobian],[Residual],[ReactionModel],[ReleaseCI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("RADIAL_GENERAL_RATE_MODEL", "FV", true, false, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Radial GRM dynamic reactions time derivative Jacobian vs FD particle", "[RadGRM],[Jacobian],[Residual],[ReactionModel],[ReleaseCI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("RADIAL_GENERAL_RATE_MODEL", "FV", false, true, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Radial GRM dynamic reactions time derivative Jacobian vs FD modified particle", "[RadGRM],[Jacobian],[Residual],[ReactionModel],[ReleaseCI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("RADIAL_GENERAL_RATE_MODEL", "FV", false, true, true, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Radial GRM dynamic reactions time derivative Jacobian vs FD bulk and particle", "[RadGRM],[Jacobian],[Residual],[ReactionModel],[ReleaseCI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("RADIAL_GENERAL_RATE_MODEL", "FV", true, true, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Radial GRM dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[RadGRM],[Jacobian],[Residual],[ReactionModel],[ReleaseCI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("RADIAL_GENERAL_RATE_MODEL", "FV", true, true, true, 1e-6, 1e-14, 9e-4);
}

inline cadet::JsonParameterProvider createColumnWithTwoCompLinearBindingThreeParticleTypes()
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("RADIAL_GENERAL_RATE_MODEL", "FV");

	const double parVolFrac[] = {0.3, 0.6, 0.1};
	const double parFactor[] = {0.9, 0.8};
	cadet::test::particle::extendModelToManyParticleTypes(jpp, 3, parFactor, parVolFrac);

	return jpp;
}

TEST_CASE("Radial GRM multi particle types dynamic reactions Jacobian vs AD bulk", "[RadGRM],[Jacobian],[AD],[ReactionModel],[ParticleType],[ReleaseCI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, false, false);
}

TEST_CASE("Radial GRM multi particle types dynamic reactions Jacobian vs AD particle", "[RadGRM],[Jacobian],[AD],[ReactionModel],[ParticleType],[ReleaseCI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, false);
}

TEST_CASE("Radial GRM multi particle types dynamic reactions Jacobian vs AD modified particle", "[RadGRM],[Jacobian],[AD],[ReactionModel],[ParticleType],[ReleaseCI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, true);
}

TEST_CASE("Radial GRM multi particle types dynamic reactions Jacobian vs AD bulk and particle", "[RadGRM],[Jacobian],[AD],[ReactionModel],[ParticleType],[ReleaseCI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, false);
}

TEST_CASE("Radial GRM multi particle types dynamic reactions Jacobian vs AD bulk and modified particle", "[RadGRM],[Jacobian],[AD],[ReactionModel],[ParticleType],[ReleaseCI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, true);
}

TEST_CASE("Radial GRM multi particle types dynamic reactions time derivative Jacobian vs FD bulk", "[RadGRM],[Jacobian],[Residual],[ReactionModel],[ParticleType],[ReleaseCI],[FD]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, false, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Radial GRM multi particle types dynamic reactions time derivative Jacobian vs FD particle", "[RadGRM],[Jacobian],[Residual],[ReactionModel],[ParticleType],[ReleaseCI],[FD]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Radial GRM multi particle types dynamic reactions time derivative Jacobian vs FD modified particle", "[RadGRM],[Jacobian],[Residual],[ReactionModel],[ParticleType],[ReleaseCI],[FD]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, true, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Radial GRM multi particle types dynamic reactions time derivative Jacobian vs FD bulk and particle", "[RadGRM],[Jacobian],[Residual],[ReactionModel],[ParticleType],[ReleaseCI],[FD]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Radial GRM multi particle types dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[RadGRM],[Jacobian],[Residual],[ReactionModel],[ParticleType],[ReleaseCI],[FD]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, true, 1e-6, 1e-14, 9e-4);
}
