// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2018: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include <catch.hpp>

#include "ColumnTests.hpp"
#include "Weno.hpp"
#include "Utils.hpp"

TEST_CASE("LRMP LWE forward vs backward flow", "[LRMP],[Simulation]")
{
	// Test all WENO orders
	for (unsigned int i = 1; i <= cadet::Weno::maxOrder(); ++i)
		cadet::test::column::testWenoForwardBackward("LUMPED_RATE_MODEL_WITH_PORES", i, 6e-8, 4e-6);
}

TEST_CASE("LRMP linear pulse vs analytic solution", "[LRMP],[Simulation],[Analytic]")
{
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITH_PORES", "/data/lrmp-pulseBenchmark.data", true, true, 512, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITH_PORES", "/data/lrmp-pulseBenchmark.data", true, false, 512, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITH_PORES", "/data/lrmp-pulseBenchmark.data", false, true, 512, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITH_PORES", "/data/lrmp-pulseBenchmark.data", false, false, 512, 6e-5, 1e-7);
}

TEST_CASE("LRMP non-binding linear pulse vs analytic solution", "[LRMP],[Simulation],[Analytic],[NonBinding]")
{
	cadet::test::column::testAnalyticNonBindingBenchmark("LUMPED_RATE_MODEL_WITH_PORES", "/data/lrmp-nonBinding.data", true, 512, 6e-5, 1e-7);
	cadet::test::column::testAnalyticNonBindingBenchmark("LUMPED_RATE_MODEL_WITH_PORES", "/data/lrmp-nonBinding.data", false, 512, 6e-5, 1e-7);
}

TEST_CASE("LRMP Jacobian forward vs backward flow", "[LRMP],[UnitOp],[Residual],[Jacobian],[AD]")
{
	// Test all WENO orders
	for (unsigned int i = 1; i <= cadet::Weno::maxOrder(); ++i)
		cadet::test::column::testJacobianWenoForwardBackward("LUMPED_RATE_MODEL_WITH_PORES", i);
}

TEST_CASE("LRMP time derivative Jacobian vs FD", "[LRMP],[UnitOp],[Residual],[Jacobian]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("LUMPED_RATE_MODEL_WITH_PORES");
}

TEST_CASE("LRMP sensitivity Jacobians", "[LRMP],[UnitOp],[Sensitivity]")
{
	cadet::test::column::testFwdSensJacobians("LUMPED_RATE_MODEL_WITH_PORES", 1e-4, 6e-7);
}

TEST_CASE("LRMP forward sensitivity vs FD", "[LRMP],[Sensitivity],[Simulation]")
{
	// Relative error is checked first, we use high absolute error for letting
	// some points that are far off pass the error test, too. This is required
	// due to errors in finite differences.
	const double fdStepSize[] = {1e-5, 1e-6, 1e-3, 1e-5};
	const double absTols[] = {6e5, 2e-2, 2e-2, 1.0};
	const double relTols[] = {5e-3, 1e-1, 5e-1, 6e-3};
	const double passRatio[] = {0.87, 0.84, 0.88, 0.95};
	cadet::test::column::testFwdSensSolutionFD("LUMPED_RATE_MODEL_WITH_PORES", true, fdStepSize, absTols, relTols, passRatio);
}

TEST_CASE("LRMP forward sensitivity forward vs backward flow", "[LRMP],[Sensitivity],[Simulation]")
{
	const double absTols[] = {50.0, 2e-10, 1.0, 5e-7};
	const double relTols[] = {2e-4, 9e-6, 5e-7, 1e-7};
	const double passRatio[] = {1.0, 0.99, 0.98, 0.99};
	cadet::test::column::testFwdSensSolutionForwardBackward("LUMPED_RATE_MODEL_WITH_PORES", absTols, relTols, passRatio);
}

TEST_CASE("LRMP consistent initialization with linear binding", "[LRMP],[ConsistentInit]")
{
	cadet::test::column::testConsistentInitializationLinearBinding("LUMPED_RATE_MODEL_WITH_PORES", 1e-12, 1e-14);
}

TEST_CASE("LRMP consistent initialization with SMA binding", "[LRMP],[ConsistentInit]")
{
	std::vector<double> y(4 + 4 * 16 + 16 * (4 + 4) + 4 * 16, 0.0);
// Optimal values:
//	const double bindingCell[] = {1.2, 2.0, 1.0, 1.5, 858.034, 66.7896, 3.53273, 2.53153, 
//		1.0, 1.8, 1.5, 1.6, 856.173, 64.457, 5.73227, 2.85286};
	const double bindingCell[] = {1.2, 2.0, 1.0, 1.5, 840.0, 63.0, 3.0, 3.0, 
		1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0};
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 + 4 * 16);
	cadet::test::util::repeat(y.data() + 4 + 4 * 16, bindingCell, 16, 8);
	cadet::test::util::populate(y.data() + 4 + 4 * 16 + 16 * (4 + 4), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 16);

	cadet::test::column::testConsistentInitializationSMABinding("LUMPED_RATE_MODEL_WITH_PORES", y.data(), 1e-14, 1e-5);
}

TEST_CASE("LRMP consistent sensitivity initialization with linear binding", "[LRMP],[ConsistentInit],[Sensitivity]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 4 + 4 * 16 + 16 * (4 + 4) + 4 * 16;
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, numDofs);
	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITH_PORES", y.data(), yDot.data(), true, 1e-14);
}

TEST_CASE("LRMP consistent sensitivity initialization with SMA binding", "[LRMP],[ConsistentInit],[Sensitivity]")
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

	cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITH_PORES", y.data(), yDot.data(), false, 1e-10);
}
