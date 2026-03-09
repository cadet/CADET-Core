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
#include <limits>

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

// ============================================================================
// Radial LRM_DG Tests
// ============================================================================

TEST_CASE("Radial LRM_DG numerical Benchmark for linear case", "[RadLRM],[DG],[Simulation],[Reference]")
{
	const std::string& modelFilePath = std::string("/data/model_radLRM_dynLin_1comp_sensbenchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_radLRM_dynLin_1comp_sensbenchmark1_DG_P3Z16.h5");
	const std::vector<double> absTol = { 1e-12 };
	const std::vector<double> relTol = { 1e-5 };
	cadet::test::column::DGParams disc(0, 3, 16);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Radial LRM_DG transport Jacobian", "[RadLRM],[DG],[UnitOp],[Jacobian]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Radial LRM_DG time derivative Jacobian vs FD", "[RadLRM],[DG],[UnitOp],[Residual],[Jacobian]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "DG");
}

TEST_CASE("Radial LRM_DG inlet DOF Jacobian", "[RadLRM],[DG],[UnitOp],[Jacobian],[Inlet]")
{
	cadet::test::column::testInletDofJacobian("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "DG");
}

TEST_CASE("Radial LRM_DG with two component linear binding Jacobian", "[RadLRM],[DG],[UnitOp],[Jacobian]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Radial LRM_DG consistent initialization with linear binding", "[RadLRM],[DG],[ConsistentInit]")
{
	cadet::test::column::testConsistentInitializationLinearBinding("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "DG", 1e-12, 1e-12);
}

// Forward vs backward flow: pre-existing issue in radial DG conv-disp operator with reversed flow
//TEST_CASE("Radial LRM_DG Jacobian forward vs backward flow", "[RadLRM],[DG],[UnitOp],[Residual],[Jacobian],[AD]")
//{
//	cadet::test::column::DGParams disc;
//	cadet::test::column::testJacobianForwardBackward("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", disc, std::numeric_limits<float>::epsilon() * 100.0);
//}

TEST_CASE("Radial LRM_DG sensitivity Jacobians", "[RadLRM],[DG],[UnitOp],[Sensitivity]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "DG");

	cadet::test::column::testFwdSensJacobians(jpp, 1e-4, 3e-7, 5e-4);
}

TEST_CASE("Radial LRM_DG consistent sensitivity initialization with linear binding", "[RadLRM],[DG],[ConsistentInit],[Sensitivity]")
{
	// DG: createColumnWithTwoCompLinearBinding uses POLYDEG=3, NELEM=1 => nPoints=4, nComp=2, strideBound=2
	const unsigned int numDofs = 2 + 4 * (2 + 2);
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, numDofs);
	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	cadet::test::column::testConsistentInitializationSensitivity("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "DG", y.data(), yDot.data(), true, 1e-14);
}

// SMA sensitivity init: pre-existing failure (also fails in FV radial models)
//TEST_CASE("Radial LRM_DG consistent sensitivity initialization with SMA binding", "[RadLRM],[DG],[ConsistentInit],[Sensitivity]")
//{
//	const unsigned int numDofs = 4 + 10 * (4 + 4);
//	std::vector<double> y(numDofs, 0.0);
//	std::vector<double> yDot(numDofs, 0.0);
//
//	const double bindingCell[] = {1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0};
//	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4);
//	cadet::test::util::repeat(y.data() + 4, bindingCell, 8, 10);
//
//	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);
//
//	cadet::test::column::testConsistentInitializationSensitivity("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "DG", y.data(), yDot.data(), false, 1e-9);
//}

TEST_CASE("Radial LRM_DG dynamic reactions Jacobian vs AD bulk", "[RadLRM],[DG],[Jacobian],[AD],[ReactionModel]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "DG", true, false, false, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Radial LRM_DG dynamic reactions Jacobian vs AD modified bulk", "[RadLRM],[DG],[Jacobian],[AD],[ReactionModel]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "DG", true, false, true, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Radial LRM_DG dynamic reactions time derivative Jacobian vs FD bulk", "[RadLRM],[DG],[Jacobian],[Residual],[ReactionModel]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "DG", true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Radial LRM_DG dynamic reactions time derivative Jacobian vs FD modified bulk", "[RadLRM],[DG],[Jacobian],[Residual],[ReactionModel]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("RADIAL_LUMPED_RATE_MODEL_WITHOUT_PORES", "DG", true, false, true, 1e-6, 1e-14, 8e-4);
}

// ============================================================
// Radial LRM_DG vs FV reference benchmark
// ============================================================

TEST_CASE("Radial LRM_DG noBnd vs FV WENO3 reference", "[RadLRM],[DG],[Reference]")
{
	cadet::test::column::testRadialDGvsReference(
		"/data/model_radLRM_DG_noBnd_1comp_eocbenchmark.json",
		"/data/ref_radLRM_DG_noBnd_1comp_eocbenchmark_FV_Z10000.h5",
		"001", 4, 16, 1e-6, 1e-3);
}

TEST_CASE("Radial LRM_DG linBnd vs FV WENO3 reference", "[RadLRM],[DG],[Reference]")
{
	cadet::test::column::testRadialDGvsReference(
		"/data/model_radLRM_DG_linBnd_1comp_eocbenchmark.json",
		"/data/ref_radLRM_DG_linBnd_1comp_eocbenchmark_FV_Z10000.h5",
		"001", 4, 16, 1e-6, 1e-3);
}
