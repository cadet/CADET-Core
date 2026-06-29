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
#include "ReactionModelTests.hpp"
#include "Utils.hpp"
#include "JsonTestModels.hpp"

TEST_CASE("Radial LRM numerical Benchmark with parameter sensitivities for linear case", "[RadLRM],[Simulation],[Reference],[Sensitivity],[CI_sens12]")
{
	const std::string& modelFilePath = std::string("/data/model_radLRM_dynLin_1comp_sensbenchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_radLRM_dynLin_1comp_sensbenchmark1_FV_Z32.h5");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1e-5, 1e-6, 1e-6, 1e-6 };
	cadet::test::column::FVParams disc(32);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("Radial LRM transport Jacobian", "[RadLRM],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV");
	cadet::test::column::testJacobianAD(jpp);
}

TEST_CASE("Radial LRM Jacobian forward vs backward flow", "[RadLRM],[UnitOp],[Residual],[Jacobian],[AD],[CI]")
{
	cadet::test::column::FVParams disc(16);
	cadet::test::column::testJacobianForwardBackward("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", disc);
}

TEST_CASE("Radial LRM time derivative Jacobian vs FD", "[RadLRM],[UnitOp],[Residual],[Jacobian],[FD],[CI]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV");
}

TEST_CASE("Radial LRM sensitivity Jacobians", "[RadLRM],[UnitOp],[Sensitivity],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV");

	cadet::test::column::testFwdSensJacobians(jpp, 1e-4, 3e-7, 5e-5);
}

// NOTE: Consistent initialization tests ARE NOT INCLUDED IN THE CI SINCE THEY ARE REDUNDANT WITH THE AXIAL FLOW TESTS

TEST_CASE("Radial LRM inlet DOF Jacobian", "[RadLRM],[UnitOp],[Jacobian],[Inlet],[CI]")
{
	cadet::test::column::testInletDofJacobian("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV");
}

TEST_CASE("Radial LRM with two component linear binding Jacobian", "[RadLRM],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV");
	cadet::test::column::testJacobianAD(jpp);
}

TEST_CASE("Radial LRM dynamic reactions Jacobian vs AD bulk", "[RadLRM],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV", true, false, false);
}

TEST_CASE("Radial LRM dynamic reactions Jacobian vs AD modified bulk", "[RadLRM],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV", true, false, true);
}

TEST_CASE("Radial LRM dynamic reactions time derivative Jacobian vs FD bulk", "[RadLRM],[Jacobian],[Residual],[ReactionModel],[FD],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV", true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Radial LRM dynamic reactions time derivative Jacobian vs FD modified bulk", "[RadLRM],[Jacobian],[Residual],[ReactionModel],[FD],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "FV", true, false, true, 1e-6, 1e-14, 8e-4);
}
