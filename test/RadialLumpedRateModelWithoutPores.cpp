// =============================================================================
//  CADET
//  
//  Copyright © The CADET Authors
//            Please see the CONTRIBUTING.md file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include <catch.hpp>

#include "ColumnTests.hpp"
#include "ReactionModelTests.hpp"
#include "Utils.hpp"
#include "JsonTestModels.hpp"

// todo add a meaningful backward flow test

// todo find analytical solution with linear binding

// todo find analytical solution without binding

// todo add (more) numerical reference (and EOC) tests

TEST_CASE("Radial LRM numerical Benchmark with parameter sensitivities for linear case", "[RadLRM],[Simulation],[Reference],[Sensitivity]") // todo CI flag: currently only runs locally but fails on server
{
	const std::string& modelFilePath = std::string("/data/model_radLRM_dynLin_1comp_sensbenchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_radLRM_dynLin_1comp_sensbenchmark1_FV_Z32.h5");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1e-5, 1e-6, 1e-6, 1e-6 };
	cadet::test::column::FVparams disc(32);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Radial LRM transport Jacobian", "[RadLRM],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV");
	cadet::test::column::testJacobianAD(jpp);
}

// NOTE: THE FOLLOWING TESTS ARE ONLY INCLUDED IN THE RELEASE CI, NOT THE STANDARD CI SINCE THEY ARE (TO A HIGH DEGREE) REDUNDANT WITH THE AXIAL FLOW TESTS

TEST_CASE("Radial LRM Jacobian forward vs backward flow", "[RadLRM],[UnitOp],[Residual],[Jacobian],[AD],[ReleaseCI]")
{
	cadet::test::column::FVparams disc(16);
	cadet::test::column::testJacobianForwardBackward("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", 16);
}

TEST_CASE("Radial LRM time derivative Jacobian vs FD", "[RadLRM],[UnitOp],[Residual],[Jacobian],[ReleaseCI],[FD]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV");
}

TEST_CASE("Radial LRM sensitivity Jacobians", "[RadLRM],[UnitOp],[Sensitivity],[ReleaseCI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV");

	cadet::test::column::testFwdSensJacobians(jpp, 1e-4, 3e-7, 5e-5);
}

//TEST_CASE("Radial LRM forward sensitivity vs FD", "[RadLRM],[Sensitivity],[Simulation],[failedFDtestLRM],[FD]") // todo (for all models) find tolerances
//{
//	// Relative error is checked first, we use high absolute error for letting
//	// some points that are far off pass the error test, too. This is required
//	// due to errors in finite differences.
//	const double fdStepSize[] = {5e-3, 5e-3, 5e-3, 1e-3};
//	const double absTols[] = {2e8, 8e-3, 2e-2, 3e-1};
//	const double relTols[] = {1e-1, 5e-1, 5e-2, 1e-2};
//	const double passRatio[] = {0.88, 0.84, 0.73, 0.87};
//	cadet::test::column::testFwdSensSolutionFD("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", false, fdStepSize, absTols, relTols, passRatio);
//}
//
//TEST_CASE("Radial LRM forward sensitivity forward vs backward flow", "[RadLRM],[Sensitivity],[Simulation],[fixLRM]") // todo (for all models) find tolerances? why is there a pass ratio here, shouldnt this be precise?
//{
//	const double absTols[] = {500.0, 8e-7, 9e-7, 2e-3};
//	const double relTols[] = {7e-3, 5e-5, 5e-5, 9e-4};
//	const double passRatio[] = {0.99, 0.97, 0.97, 0.98};
//	cadet::test::column::testFwdSensSolutionForwardBackward("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", absTols, relTols, passRatio);
//}

TEST_CASE("Radial LRM consistent initialization with linear binding", "[RadLRM],[ConsistentInit],[ReleaseCI]")
{
	cadet::test::column::testConsistentInitializationLinearBinding("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV", 1e-12, 1e-12);
}

//TEST_CASE("Radial LRM consistent initialization with SMA binding", "[RadLRM],[ConsistentInit],[fixLRM]") // todo (for all models) fix
//{
//	std::vector<double> y(4 + 16 * (4 + 4), 0.0);
//	// Optimal values:
//	//	const double bindingCell[] = {1.2, 2.0, 1.0, 1.5, 858.034, 66.7896, 3.53273, 2.53153, 
//	//		1.0, 1.8, 1.5, 1.6, 856.173, 64.457, 5.73227, 2.85286};
//	const double bindingCell[] = {1.2, 2.0, 1.0, 1.5, 840.0, 63.0, 3.0, 3.0, 
//		1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0};
//	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4);
//	cadet::test::util::repeat(y.data() + 4, bindingCell, 16, 8);
//
//	cadet::test::column::testConsistentInitializationSMABinding("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", y.data(), 1e-14, 1e-5);
//}

TEST_CASE("Radial LRM consistent sensitivity initialization with linear binding", "[RadLRM],[ConsistentInit],[Sensitivity],[ReleaseCI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 4 + 16 * (4 + 4);
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, numDofs);
	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	cadet::test::column::testConsistentInitializationSensitivity("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV", y.data(), yDot.data(), true, 1e-12);
}

TEST_CASE("Radial LRM consistent sensitivity initialization with SMA binding", "[RadLRM],[ConsistentInit],[Sensitivity],[ReleaseCI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 4 + 16 * (4 + 4);
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);

	const double bindingCell[] = {1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0};
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4);
	cadet::test::util::repeat(y.data() + 4, bindingCell, 8, 16);

	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	cadet::test::column::testConsistentInitializationSensitivity("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV", y.data(), yDot.data(), false, 1e-9);
}

TEST_CASE("Radial LRM inlet DOF Jacobian", "[RadLRM],[UnitOp],[Jacobian],[Inlet],[ReleaseCI]")
{
	cadet::test::column::testInletDofJacobian("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV");
}

TEST_CASE("Radial LRM with two component linear binding Jacobian", "[RadLRM],[UnitOp],[Jacobian],[ReleaseCI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV");
	cadet::test::column::testJacobianAD(jpp);
}

TEST_CASE("Radial LRM dynamic reactions Jacobian vs AD bulk", "[RadLRM],[Jacobian],[AD],[ReactionModel],[ReleaseCI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV", true, false, false);
}

TEST_CASE("Radial LRM dynamic reactions Jacobian vs AD modified bulk", "[RadLRM],[Jacobian],[AD],[ReactionModel],[ReleaseCI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV", true, false, true);
}

TEST_CASE("Radial LRM dynamic reactions time derivative Jacobian vs FD bulk", "[RadLRM],[Jacobian],[Residual],[ReactionModel],[ReleaseCI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV", true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Radial LRM dynamic reactions time derivative Jacobian vs FD modified bulk", "[RadLRM],[Jacobian],[Residual],[ReactionModel],[ReleaseCI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV", true, false, true, 1e-6, 1e-14, 8e-4);
}
