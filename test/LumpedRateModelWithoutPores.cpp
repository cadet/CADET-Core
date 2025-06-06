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
#include "Weno.hpp"
#include "Utils.hpp"
#include "JsonTestModels.hpp"

TEST_CASE("LRM LWE forward vs backward flow", "[LRM],[FV],[Simulation],[CI]")
{
	cadet::test::column::FVparams disc;

	// Test all WENO orders
	for (unsigned int i = 1; i <= cadet::Weno::maxOrder(); ++i)
	{
		disc.setWenoOrder(i);
		cadet::test::column::testForwardBackward("LUMPED_RATE_MODEL_WITHOUT_PORES", disc, 6e-9, 6e-4);
	}
}

TEST_CASE("LRM linear pulse vs analytic solution", "[LRM],[FV],[Simulation],[Reference],[Analytic],[CI]")
{
	cadet::test::column::FVparams disc(1024);
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITHOUT_PORES", "/data/lrm-pulseBenchmark.data", true, true, disc, 2e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITHOUT_PORES", "/data/lrm-pulseBenchmark.data", true, false, disc, 2e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITHOUT_PORES", "/data/lrm-pulseBenchmark.data", false, true, disc, 2e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITHOUT_PORES", "/data/lrm-pulseBenchmark.data", false, false, disc, 2e-5, 1e-7);
}

TEST_CASE("LRM non-binding linear pulse vs analytic solution", "[LRM],[FV],[Simulation],[Reference],[Analytic],[NonBinding],[CI]")
{
	cadet::test::column::FVparams disc(1024);
	cadet::test::column::testAnalyticNonBindingBenchmark("LUMPED_RATE_MODEL_WITHOUT_PORES", "/data/lrm-nonBinding.data", true, disc, 2e-5, 1e-7);
	cadet::test::column::testAnalyticNonBindingBenchmark("LUMPED_RATE_MODEL_WITHOUT_PORES", "/data/lrm-nonBinding.data", false, disc, 2e-5, 1e-7);
}

TEST_CASE("LRM Jacobian forward vs backward flow", "[LRM],[FV],[UnitOp],[Residual],[Jacobian],[AD],[CI]")
{
	cadet::test::column::FVparams disc;

	// Test all WENO orders
	for (unsigned int i = 1; i <= cadet::Weno::maxOrder(); ++i)
	{
		disc.setWenoOrder(i);
		cadet::test::column::testJacobianForwardBackward("LUMPED_RATE_MODEL_WITHOUT_PORES", disc);
	}
}

TEST_CASE("LRM numerical Benchmark with parameter sensitivities for linear case", "[LRM],[FV],[Simulation],[Reference],[Sensitivity],[CI_sens5]")
{
	const std::string& modelFilePath = std::string("/data/model_LRM_dynLin_1comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_LRM_dynLin_1comp_sensbenchmark1_FV_Z32.h5");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1e-4, 1e-4, 1e-4, 1e-4 };

	cadet::test::column::FVparams disc(32);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("LRM numerical Benchmark with parameter sensitivities for SMA LWE case", "[LRM],[FV],[Simulation],[Reference],[Sensitivity],[CI_sens6]")
{
	const std::string& modelFilePath = std::string("/data/model_LRM_reqSMA_4comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_LRM_reqSMA_4comp_sensbenchmark1_FV_Z32.h5");
	const std::vector<double> absTol = { 1e-12, 1e-3, 1e-3, 1e-6 };
	const std::vector<double> relTol = { 1e-4, 1.0, 1.0, 1.0 };

	cadet::test::column::FVparams disc(32);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "000", absTol, relTol, disc, true);
}

TEST_CASE("LRM numerical EOC Benchmark with parameter sensitivities for linear case", "[LRM],[FV],[EOC],[EOC_LRM_FV]")
{
	const std::string& modelFilePath = std::string("/data/model_LRM_dynLin_1comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_LRM_dynLin_1comp_sensbenchmark1_FV_Z131072.h5");
	const std::string& convFilePath = std::string("/data/convergence_LRM_dynLin_1comp_sensbenchmark1.json");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1e-2, 1e-3, 1e-1, 1e-2 };

	cadet::test::column::FVparams disc(16);
	cadet::test::column::testEOCReferenceBenchmark(modelFilePath, refFilePath, convFilePath, "001", absTol, relTol, 4, disc, true);
}

TEST_CASE("LRM numerical EOC Benchmark with parameter sensitivities for SMA LWE case", "[LRM],[FV],[EOC],[EOC_LRM_FV]")
{
	const std::string& modelFilePath = std::string("/data/model_LRM_reqSMA_4comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_LRM_reqSMA_4comp_sensbenchmark1_FV_Z4096.h5");
	const std::string& convFilePath = std::string("/data/convergence_LRM_reqSMA_4comp_sensbenchmark1.json");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1e-4, 1e-4, 1e-4, 1e-4 };

	cadet::test::column::FVparams disc(8);
	cadet::test::column::testEOCReferenceBenchmark(modelFilePath, refFilePath, convFilePath, "000", absTol, relTol, 2, disc, true);
}

TEST_CASE("LRM time derivative Jacobian vs FD", "[LRM],[FV],[UnitOp],[Residual],[Jacobian],[CI],[FD]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("LUMPED_RATE_MODEL_WITHOUT_PORES", "FV");
}

TEST_CASE("LRM sensitivity Jacobians", "[LRM],[FV],[UnitOp],[Sensitivity],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("LUMPED_RATE_MODEL_WITHOUT_PORES", "FV");

	cadet::test::column::testFwdSensJacobians(jpp, 1e-4, 3e-7, 5e-5);
}

//TEST_CASE("LRM forward sensitivity vs FD", "[LRM],[FV],[Sensitivity],[Simulation],[failedFDtestLRM],[FD]") // todo fix tolerances
//{
//	// Relative error is checked first, we use high absolute error for letting
//	// some points that are far off pass the error test, too. This is required
//	// due to errors in finite differences.
//	const double fdStepSize[] = {5e-3, 5e-3, 5e-3, 1e-3};
//	const double absTols[] = {2e8, 8e-3, 2e-2, 3e-1};
//	const double relTols[] = {1e-1, 5e-1, 5e-2, 1e-2};
//	const double passRatio[] = {0.88, 0.84, 0.73, 0.87};
//	cadet::test::column::testFwdSensSolutionFD("LUMPED_RATE_MODEL_WITHOUT_PORES", false, fdStepSize, absTols, relTols, passRatio);
//}
//
//TEST_CASE("LRM forward sensitivity forward vs backward flow", "[LRM],[FV],[Sensitivity],[Simulation],[fixLRM]") // todo fix tolerances? why is there a pass ratio here, shouldnt this be precise?
//{
//	const double absTols[] = {500.0, 8e-7, 9e-7, 2e-3};
//	const double relTols[] = {7e-3, 5e-5, 5e-5, 9e-4};
//	const double passRatio[] = {0.99, 0.97, 0.97, 0.98};
//	cadet::test::column::testFwdSensSolutionForwardBackward("LUMPED_RATE_MODEL_WITHOUT_PORES", absTols, relTols, passRatio);
//}

TEST_CASE("LRM consistent initialization with linear binding", "[LRM],[FV],[ConsistentInit],[CI]")
{
	cadet::test::column::testConsistentInitializationLinearBinding("LUMPED_RATE_MODEL_WITHOUT_PORES", "FV", 1e-12, 1e-12);
}

//TEST_CASE("LRM consistent initialization with SMA binding", "[LRM],[FV],[ConsistentInit],[fixLRM]") // todo fix
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
//	cadet::test::column::testConsistentInitializationSMABinding("LUMPED_RATE_MODEL_WITHOUT_PORES", y.data(), 1e-14, 1e-5);
//}

TEST_CASE("LRM consistent sensitivity initialization with linear binding", "[LRM],[FV],[ConsistentInit],[Sensitivity],[CI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 4 + 16 * (4 + 4);
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, numDofs);
	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITHOUT_PORES", "FV", y.data(), yDot.data(), true, 1e-14);
}

TEST_CASE("LRM consistent sensitivity initialization with SMA binding", "[LRM],[FV],[ConsistentInit],[Sensitivity],[CI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 4 + 16 * (4 + 4);
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);

	const double bindingCell[] = {1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0};
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4);
	cadet::test::util::repeat(y.data() + 4, bindingCell, 8, 16);

	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITHOUT_PORES", "FV", y.data(), yDot.data(), false, 1e-9);
}

TEST_CASE("LRM inlet DOF Jacobian", "[LRM],[FV],[UnitOp],[Jacobian],[Inlet],[CI]")
{
	cadet::test::column::testInletDofJacobian("LUMPED_RATE_MODEL_WITHOUT_PORES", "FV");
}

TEST_CASE("LRM transport Jacobian", "[LRM],[FV],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "LUMPED_RATE_MODEL_WITHOUT_PORES", "FV");
	cadet::test::column::testJacobianAD(jpp);
}

TEST_CASE("LRM with two component linear binding Jacobian", "[LRM],[FV],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("LUMPED_RATE_MODEL_WITHOUT_PORES", "FV");
	cadet::test::column::testJacobianAD(jpp);
}

TEST_CASE("LRM dynamic reactions Jacobian vs AD bulk", "[LRM],[FV],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("LUMPED_RATE_MODEL_WITHOUT_PORES", "FV", true, false, false);
}

TEST_CASE("LRM dynamic reactions Jacobian vs AD modified bulk", "[LRM],[FV],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("LUMPED_RATE_MODEL_WITHOUT_PORES", "FV", true, false, true);
}

TEST_CASE("LRM dynamic reactions time derivative Jacobian vs FD bulk", "[LRM],[FV],[Jacobian],[Residual],[ReactionModel],[CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("LUMPED_RATE_MODEL_WITHOUT_PORES", "FV", true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("LRM dynamic reactions time derivative Jacobian vs FD modified bulk", "[LRM],[FV],[Jacobian],[Residual],[ReactionModel],[CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("LUMPED_RATE_MODEL_WITHOUT_PORES", "FV", true, false, true, 1e-6, 1e-14, 8e-4);
}
