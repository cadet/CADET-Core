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

TEST_CASE("Column_1D as GRM LWE forward vs backward flow", "[Column_1D],[DG],[DG1D],[Simulation],[CI]")
{
	cadet::test::column::DGParams disc;

	// Test all integration modes
	for (int i = 0; i <= 1; i++)
	{
		disc.setBulkDiscParam("EXACT_INTEGRATION", i);
		cadet::test::column::testForwardBackward("COLUMN_1D_GRM", disc, 1e-9, 2e-4);
	}
}

TEST_CASE("Column_1D as GRM linear pulse vs analytic solution", "[Column_1D],[DG],[DG1D],[Simulation],[Analytic],[CI]")
{
	cadet::test::column::DGParams disc;
	cadet::test::column::testAnalyticBenchmark("COLUMN_1D_GRM", "/data/grm-pulseBenchmark.data", true, true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("COLUMN_1D_GRM", "/data/grm-pulseBenchmark.data", true, false, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("COLUMN_1D_GRM", "/data/grm-pulseBenchmark.data", false, true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("COLUMN_1D_GRM", "/data/grm-pulseBenchmark.data", false, false, disc, 6e-5, 1e-7);
}

TEST_CASE("Column_1D as LRMP linear pulse vs analytic solution", "[Column_1D],[DG],[DG1D],[Simulation],[Analytic],[CI]")
{
	cadet::test::column::DGParams disc;
	cadet::test::column::testAnalyticBenchmark("COLUMN_1D_LRMP", "/data/lrmp-pulseBenchmark.data", true, true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("COLUMN_1D_LRMP", "/data/lrmp-pulseBenchmark.data", true, false, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("COLUMN_1D_LRMP", "/data/lrmp-pulseBenchmark.data", false, true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("COLUMN_1D_LRMP", "/data/lrmp-pulseBenchmark.data", false, false, disc, 6e-5, 1e-7);
}

TEST_CASE("Column_1D as GRM non-binding linear pulse vs analytic solution", "[Column_1D],[DG],[DG1D],[Simulation],[Analytic],[NonBinding],[CI]")
{
	cadet::test::column::DGParams disc;
	cadet::test::column::testAnalyticNonBindingBenchmark("COLUMN_1D_GRM", "/data/grm-nonBinding.data", true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticNonBindingBenchmark("COLUMN_1D_GRM", "/data/grm-nonBinding.data", false, disc, 6e-5, 1e-7);
}

TEST_CASE("Column_1D as LRMP non-binding linear pulse vs analytic solution", "[Column_1D],[DG],[DG1D],[Simulation],[Analytic],[NonBinding],[CI]")
{
	cadet::test::column::DGParams disc;
	cadet::test::column::testAnalyticNonBindingBenchmark("COLUMN_1D_LRMP", "/data/lrmp-nonBinding.data", true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticNonBindingBenchmark("COLUMN_1D_LRMP", "/data/lrmp-nonBinding.data", false, disc, 6e-5, 1e-7);
}

TEST_CASE("Column_1D as GRM Jacobian forward vs backward flow", "[Column_1D],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[AD],[CI]")
{
	cadet::test::column::DGParams disc;
	disc.setBulkDiscParam("EXACT_INTEGRATION", 1);
	cadet::test::column::testJacobianForwardBackward("COLUMN_1D_GRM", disc, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as DPF numerical Benchmark1", "[Column_1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI]")
{
	std::string modelFilePath = std::string("/data/model_COL1D_DPFR_1comp_benchmark1.json");
	std::string refFilePath = std::string("/data/ref_COL1D_DPFR_1comp_benchmark1_DG_P3Z8.h5");
	const std::vector<double> absTol = { 1e-12 };
	const std::vector<double> relTol = { 1e-4 };

	cadet::test::column::DGParams disc(0, 3, 8, 3, 1);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("Column_1D as GRM numerical Benchmark1 with parameter sensitivities for linear case", "[Column_1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI_sens18]")
{
	std::string modelFilePath = std::string("/data/model_COL1D_GRM_dynLin_1comp_benchmark1.json");
	std::string refFilePath = std::string("/data/ref_GRM_dynLin_1comp_sensbenchmark1_cDG_P3Z8_GSM_parP3parZ1.h5");
	const std::vector<double> absTol = { 1e-12, 1e-6, 1e-6, 1e-12 };
	const std::vector<double> relTol = { 1e-4, 1e-3, 1e-4, 1e-4 };

	cadet::test::column::DGParams disc(0, 3, 8, 3, 1);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("Column_1D as LRMP numerical Benchmark with parameter sensitivities for linear case", "[Column_1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI_sens22]")
{
	const std::string& modelFilePath = std::string("/data/model_COL1D_LRMP_dynLin_1comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_LRMP_dynLin_1comp_sensbenchmark1_DG_P3Z8.h5");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1.0, 1.0, 1.0, 1.0 };

	cadet::test::column::DGParams disc(0, 3, 8);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("Column_1D as GRM numerical Benchmark2 with parameter sensitivities for linear case", "[Column_1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI_sens18]")
{
	std::string modelFilePath = std::string("/data/model_COL1D_GRM_dynLin_1comp_sensbenchmark2.json");
	std::string refFilePath = std::string("/data/ref_GRM_dynLin_1comp_sensbenchmark2_cDG_P3Z8_GSM_parP3parZ1.h5");
	const std::vector<double> absTol = { 1e-12, 1e-6, 1e-6, 1e-12 };
	const std::vector<double> relTol = { 1e-4, 1e-3, 1e-4, 1e-4 };

	cadet::test::column::DGParams disc(0, 3, 8, 3, 1);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("Column_1D as GRM numerical Benchmark with parameter sensitivities and multiplexing for 2parType 2comp linear case", "[Column_1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI_sens19]")
{
	const std::string modelFilePath = std::string("/data/model_COL1DparType2_GRM_dynLin_2comp_sensbenchmark1.json");
	const std::string refFilePath = std::string("/data/ref_GRMparType2_dynLin_2comp_sensbenchmark1_cDG_P3Z4_DGexInt_parP3parZ2.h5");
	const std::vector<double> absTol = { 1e-12, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6 };
	const std::vector<double> relTol = { 1e-4, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1 };

	const cadet::test::column::DGParams disc(0, 3, 4, 3, 2);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("Column_1D as GRM numerical Benchmark with parameter sensitivities for SMA LWE case", "[Column_1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI_sens20]")
{
	const std::string modelFilePath = std::string("/data/model_COL1D_GRM_reqSMA_4comp_benchmark1.json");
	const std::string refFilePath = std::string("/data/ref_GRM_reqSMA_4comp_sensbenchmark1_cDG_P3Z8_GSM_parP3parZ1.h5");
	const std::vector<double> absTol = { 1e-12, 1e-6, 1e-6, 1e-12 };
	const std::vector<double> relTol = { 1.0, 1.0, 1.0, 1.0 };

	cadet::test::column::DGParams disc(0, 3, 8, 3, 1);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "000", absTol, relTol, disc, true);
}

TEST_CASE("Column_1D as LRMP numerical Benchmark with parameter sensitivities for SMA LWE case", "[Column_1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI_sens23]")
{
	const std::string& modelFilePath = std::string("/data/model_COL1D_LRMP_reqSMA_4comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_LRMP_reqSMA_4comp_sensbenchmark1_DG_P3Z8.h5");
	const std::vector<double> absTol = { 1e-12, 1e-12, 5e-10, 1e-12 };
	const std::vector<double> relTol = { 1.0, 1.0, 1.0, 1.0 };

	cadet::test::column::DGParams disc(0, 3, 8);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "000", absTol, relTol, disc, true);
}

TEST_CASE("Column_1D with mixed general rate and homogeneous particle types and linear binding numerical Benchmark with parameter sensitivities", "[Column_1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI_sens24]")
{
	const std::string& modelFilePath = std::string("/data/model_COL1D_2parTypeMixed_2comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_COL1D_2parTypeMixed_2comp_benchmark1.h5");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

	cadet::test::column::DGParams disc(0, 3, 5, 3, 1);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("Column_1D as GRM LWE DGSEM and GSM particle discretization yields similar accuracy", "[Column_1D],[DG],[DG1D],[Simulation],[CI]")
{
	cadet::JsonParameterProvider jpp = createLWE("COLUMN_1D_GRM", "DG");
	cadet::test::column::DGParams disc(0, 4, 2, 3, 1); // Note that we want to employ only a single particle element
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

TEST_CASE("Column_1D as GRM time derivative Jacobian vs FD", "[Column_1D],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[CI],[FD]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("COLUMN_1D_GRM", "DG", 1e-6, 0.0, 9e-4);
}

TEST_CASE("Column_1D as LRMP time derivative Jacobian vs FD", "[Column_1D],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[CI]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("COLUMN_1D_LRMP", "DG");
}

TEST_CASE("Column_1D as GRM sensitivity Jacobians", "[Column_1D],[DG],[DG1D],[UnitOp],[Sensitivity],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_1D_GRM", "DG");

	cadet::test::column::testFwdSensJacobians(jpp, 1e-4, 6e-7, 1e-4);
}

TEST_CASE("Column_1D as LRMP sensitivity Jacobians", "[Column_1D],[DG],[DG1D],[UnitOp],[Sensitivity],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_1D_LRMP", "DG");

	cadet::test::column::testFwdSensJacobians(jpp, 1e-4, 6e-7);
}

//// todo fix: not just adjust tolerances as in FV but theres an actual error here: access violation in densematrix
//TEST_CASE("Column_1D as GRM forward sensitivity vs FD", "[Column_1D],[DG],[DG1D],[Sensitivity],[Simulation],[todo]")
//{
//	// Relative error is checked first, we use high absolute error for letting
//	// some points that are far off pass the error test, too. This is required
//	// due to errors in finite differences.
//	const double fdStepSize[] = { 5e-5, 1e-4, 1e-4, 1e-3 };
//	const double absTols[] = { 3e5, 2e-3, 2e-4, 5.0 };
//	const double relTols[] = { 5e-3, 7e-2, 8e-2, 1e-4 };
//	const double passRatio[] = { 0.95, 0.9, 0.91, 0.83 };
//	cadet::test::column::testFwdSensSolutionFD("COLUMN_1D_GRM", "DG", false, fdStepSize, absTols, relTols, passRatio);
//}

//// todo fix: not just adjust tolerances as in FV but theres an actual error here: access violation in densematrix
//TEST_CASE("Column_1D as GRM forward sensitivity forward vs backward flow", "[Column_1D],[DG],[DG1D],[Sensitivity],[Simulation],[todo]")
//{
//	const double absTols[] = { 4e-5, 1e-11, 1e-11, 8e-9 };
//	const double relTols[] = { 6e-9, 5e-8, 5e-6, 5e-10 };
//	const double passRatio[] = { 0.99, 0.95, 0.98, 0.98 };
//	cadet::test::column::testFwdSensSolutionForwardBackward("COLUMN_1D_GRM", "DG", absTols, relTols, passRatio);
//}

// todo fix consistent initialization for AD with req binding
TEST_CASE("Column_1D as GRM consistent initialization with linear binding", "[Column_1D],[DG],[DG1D],[ConsistentInit],[CI]")
{
	cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_1D_GRM", "DG", 1e-12, 1e-14, 0, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_1D_GRM", "DG", 1e-12, 1e-12, 1, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_1D_GRM", "DG", 1e-12, 1e-14, 0, 1);
	//cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_1D_GRM", "DG", 1e-12, 1e-14, 1, 1);
}

// todo fix consistent initialization for AD with req binding
TEST_CASE("Column_1D as LRMP consistent initialization with linear binding", "[Column_1D],[DG],[DG1D],[ConsistentInit],[CI]")
{
	cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_1D_LRMP", "DG", 1e-12, 1e-12, 0, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_1D_LRMP", "DG", 1e-12, 1e-12, 1, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_1D_LRMP", "DG", 1e-12, 1e-12, 0, 1);
	//cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_1D_LRMP", "DG", 1e-12, 1e-12, 1, 1);
}

//// todo fix consistent initialization for SMA (initialization not completely correct; AD gives assertion error)
//TEST_CASE("Column_1D as GRM consistent initialization with SMA binding", "[Column_1D],[DG],[DG1D],[ConsistentInit],[todo]")
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
//	cadet::test::column::testConsistentInitializationSMABinding("COLUMN_1D_GRM", "DG", y.data(), 1e-14, 1e-5);
//}

// todo fix kinetic binding sensitivity init
TEST_CASE("Column_1D as GRM consistent sensitivity initialization with linear binding", "[Column_1D],[DG],[DG1D],[ConsistentInit],[Sensitivity],[CI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 4 + 4 * 16 + 16 * 4 * (4 + 4) + 4 * 16;
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, numDofs);
	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	//cadet::test::column::testConsistentInitializationSensitivity("COLUMN_1D_GRM", "DG", y.data(), yDot.data(), true, 1e-14, 0, 0);
	cadet::test::column::testConsistentInitializationSensitivity("COLUMN_1D_GRM", "DG", y.data(), yDot.data(), true, 1e-14, 1, 0);
	//cadet::test::column::testConsistentInitializationSensitivity("COLUMN_1D_GRM", "DG", y.data(), yDot.data(), true, 1e-14, 0, 1);
	cadet::test::column::testConsistentInitializationSensitivity("COLUMN_1D_GRM", "DG", y.data(), yDot.data(), true, 1e-14, 1, 1);
}

// todo fix kinetic binding sensitivity init
TEST_CASE("Column_1D LRMP consistent sensitivity initialization with linear binding", "[Column_1D],[DG],[DG1D],[ConsistentInit],[Sensitivity],[CI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 2 + 2 * 10 + 10 * (2 + 2);
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, numDofs);
	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	//cadet::test::column::testConsistentInitializationSensitivity("COLUMN_1D_LRMP", "DG", y.data(), yDot.data(), true, 1e-14, 0, 0); // doesnt work
	cadet::test::column::testConsistentInitializationSensitivity("COLUMN_1D_LRMP", "DG", y.data(), yDot.data(), true, 1e-14, 1, 0); // works
	//cadet::test::column::testConsistentInitializationSensitivity("COLUMN_1D_LRMP", "DG", y.data(), yDot.data(), true, 1e-14, 0, 1); // doesnt work
	cadet::test::column::testConsistentInitializationSensitivity("COLUMN_1D_LRMP", "DG", y.data(), yDot.data(), true, 1e-14, 1, 1); // works
}

//// todo fix memory stuff (works for FV) 
//TEST_CASE("Column_1D as GRM consistent sensitivity initialization with SMA binding", "[Column_1D],[DG],[DG1D],[ConsistentInit],[Sensitivity],[fffffffiujbnlk]")
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
//	cadet::test::column::testConsistentInitializationSensitivity("COLUMN_1D_GRM", "DG", y.data(), yDot.data(), false, 1e-9);
//}

TEST_CASE("Column_1D Jacobian for 2ParType with general rate and homogeneous particle with two component linear binding", "[Column_1D],[DG],[DG1D],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumn2ParType1GeneralRate1HomoParticleBothWithTwoCompLinearJson("COLUMN_MODEL_1D", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as GRM inlet DOF Jacobian", "[Column_1D],[DG],[DG1D],[UnitOp],[Jacobian],[Inlet],[CI]")
{
	cadet::test::column::testInletDofJacobian("COLUMN_1D_GRM", "DG");
}

TEST_CASE("Column_1D as GRM transport Jacobian", "[Column_1D],[DG],[DG1D],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "COLUMN_1D_GRM", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as GRM with two component linear binding Jacobian", "[Column_1D],[DG],[DG1D],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_1D_GRM", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as GRM LWE one vs two identical particle types match", "[Column_1D],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testOneVsTwoIdenticalParticleTypes("COLUMN_1D_GRM", "DG", 2e-8, 5e-5);
}

TEST_CASE("Column_1D as GRM LWE separate identical particle types match", "[Column_1D],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testSeparateIdenticalParticleTypes("COLUMN_1D_GRM", "DG", 2e-8, 5e-5);
}

TEST_CASE("Column_1D as GRM linear binding single particle matches particle distribution", "[Column_1D],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearMixedParticleTypes("COLUMN_1D_GRM", "DG", 5e-8, 5e-5);
}

TEST_CASE("Column_1D as GRM multiple particle types Jacobian analytic vs AD", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ParticleType],[CI]")
{
	cadet::test::particle::testJacobianMixedParticleTypes("COLUMN_1D_GRM", "DG", 5e-12, true);
}

TEST_CASE("Column_1D as GRM multiple particle types time derivative Jacobian vs FD", "[Column_1D],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[ParticleType],[CI]")
{
	cadet::test::particle::testTimeDerivativeJacobianMixedParticleTypesFD("COLUMN_1D_GRM", "DG", 1e-6, 0.0, 9e-4, true);
}

TEST_CASE("Column_1D as GRM multiple spatially dependent particle types Jacobian analytic vs AD", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ParticleType],[CI]")
{
	cadet::test::particle::testJacobianSpatiallyMixedParticleTypes("COLUMN_1D_GRM", "DG", 1e-11, true);
}

TEST_CASE("Column_1D as GRM linear binding single particle matches spatially dependent particle distribution", "[Column_1D],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearSpatiallyMixedParticleTypes("COLUMN_1D_GRM", "DG", 5e-8, 5e-5);
}

TEST_CASE("Column_1D as GRM dynamic reactions Jacobian vs AD bulk", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_1D_GRM", "DG", true, false, false, 1e-10);
}

TEST_CASE("Column_1D as GRM dynamic reactions Jacobian vs AD particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_1D_GRM", "DG", false, true, false, 1e-10);
}

TEST_CASE("Column_1D as GRM dynamic reactions Jacobian vs AD modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_1D_GRM", "DG", false, true, true, 1e-10);
}

TEST_CASE("Column_1D as GRM dynamic reactions Jacobian vs AD bulk and particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_1D_GRM", "DG", true, true, false, 1e-10);
}

TEST_CASE("Column_1D as GRM dynamic reactions Jacobian vs AD bulk and modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_1D_GRM", "DG", true, true, true, 1e10);
}

TEST_CASE("Column_1D as GRM dynamic reactions time derivative Jacobian vs FD bulk", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_1D_GRM", "DG", true, false, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_1D as GRM dynamic reactions time derivative Jacobian vs FD particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_1D_GRM", "DG", false, true, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_1D as GRM dynamic reactions time derivative Jacobian vs FD modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_1D_GRM", "DG", false, true, true, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_1D as GRM dynamic reactions time derivative Jacobian vs FD bulk and particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_1D_GRM", "DG", true, true, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_1D as GRM dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_1D_GRM", "DG", true, true, true, 1e-6, 1e-14, 9e-4);
}

inline cadet::JsonParameterProvider createColumnWithTwoCompLinearBindingThreeParticleTypes()
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_1D_GRM", "DG");

	const double parVolFrac[] = { 0.3, 0.6, 0.1 };
	const double parFactor[] = { 0.9, 0.8 };
	cadet::test::particle::extendModelToManyParticleTypes(jpp, 3, parFactor, parVolFrac, true);

	return jpp;
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions Jacobian vs AD bulk", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, false, false);
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions Jacobian vs AD particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, false);
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions Jacobian vs AD modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, true);
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions Jacobian vs AD bulk and particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, false);
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions Jacobian vs AD bulk and modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, true);
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD bulk", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, false, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, true, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD bulk and particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, true, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_1D as LRMP inlet DOF Jacobian", "[Column_1D],[DG],[DG1D],[UnitOp],[Jacobian],[Inlet],[CI]")
{
	cadet::test::column::testInletDofJacobian("COLUMN_1D_LRMP", "DG");
}

TEST_CASE("Column_1D as LRMP transport Jacobian", "[Column_1D],[DG],[DG1D],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "COLUMN_1D_LRMP", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as LRMP with two component linear binding Jacobian", "[Column_1D],[DG],[DG1D],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_1D_LRMP", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as LRMP LWE one vs two identical particle types match", "[Column_1D],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testOneVsTwoIdenticalParticleTypes("COLUMN_1D_LRMP", "DG", 2.2e-8, 6e-5);
}

TEST_CASE("Column_1D as LRMP LWE separate identical particle types match", "[Column_1D],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testSeparateIdenticalParticleTypes("COLUMN_1D_LRMP", "DG", 1e-12, 1e-12);
}

TEST_CASE("Column_1D as LRMP linear binding single particle matches particle distribution", "[Column_1D],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearMixedParticleTypes("COLUMN_1D_LRMP", "DG", 5e-8, 5e-5);
}

TEST_CASE("Column_1D as LRMP multiple particle types Jacobian analytic vs AD", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ParticleType],[CI]")
{
	cadet::test::particle::testJacobianMixedParticleTypes("COLUMN_1D_LRMP", "DG", 1e-11, true);
}

TEST_CASE("Column_1D as LRMP multiple particle types time derivative Jacobian vs FD", "[Column_1D],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[ParticleType],[CI]")
{
	cadet::test::particle::testTimeDerivativeJacobianMixedParticleTypesFD("COLUMN_1D_LRMP", "DG", 1e-6, 0.0, 9e-4, true);
}

TEST_CASE("Column_1D as LRMP multiple spatially dependent particle types Jacobian analytic vs AD", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ParticleType],[CI]")
{
	cadet::test::particle::testJacobianSpatiallyMixedParticleTypes("COLUMN_1D_LRMP", "DG", std::numeric_limits<float>::epsilon() * 100.0, true);
}

TEST_CASE("Column_1D as LRMP linear binding single particle matches spatially dependent particle distribution", "[Column_1D],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearSpatiallyMixedParticleTypes("COLUMN_1D_LRMP", "DG", 5e-8, 5e-5);
}

TEST_CASE("Column_1D as LRMP dynamic reactions Jacobian vs AD bulk", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_1D_LRMP", "DG", true, false, false, std::numeric_limits<float>::epsilon() * 100.0);
}
TEST_CASE("Column_1D as LRMP dynamic reactions Jacobian vs AD particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_1D_LRMP", "DG", false, true, false, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as LRMP dynamic reactions Jacobian vs AD modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_1D_LRMP", "DG", false, true, true, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as LRMP dynamic reactions Jacobian vs AD bulk and particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_1D_LRMP", "DG", true, true, false, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as LRMP dynamic reactions Jacobian vs AD bulk and modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_1D_LRMP", "DG", true, true, true, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as LRMP dynamic reactions time derivative Jacobian vs FD bulk", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_1D_LRMP", "DG", true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_1D as LRMP dynamic reactions time derivative Jacobian vs FD particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_1D_LRMP", "DG", false, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_1D as LRMP dynamic reactions time derivative Jacobian vs FD modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_1D_LRMP", "DG", false, true, true, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_1D as LRMP dynamic reactions time derivative Jacobian vs FD bulk and particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_1D_LRMP", "DG", true, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_1D as LRMP dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_1D_LRMP", "DG", true, true, true, 1e-6, 1e-14, 8e-4);
}

inline cadet::JsonParameterProvider createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes()
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_1D_LRMP", "DG");

	const double parVolFrac[] = { 0.3, 0.6, 0.1 };
	const double parFactor[] = { 0.9, 0.8 };
	cadet::test::particle::extendModelToManyParticleTypes(jpp, 3, parFactor, parVolFrac, true);

	return jpp;
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions Jacobian vs AD bulk", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, false, false, 1e-8);
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions Jacobian vs AD particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, false, 1e-8);
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions Jacobian vs AD modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, true, 1e-8);
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions Jacobian vs AD bulk and particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, false, 1e-8);
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions Jacobian vs AD bulk and modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, true, 1e-8);
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD bulk", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, true, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD bulk and particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[Column_1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, true, 1e-6, 1e-14, 8e-4);
}