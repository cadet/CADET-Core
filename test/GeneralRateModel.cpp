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

TEST_CASE("GRM LWE forward vs backward flow", "[GRM],[Simulation]")
{
	// Test all WENO orders
	for (unsigned int i = 1; i <= cadet::Weno::maxOrder(); ++i)
		cadet::test::column::testWenoForwardBackward("GENERAL_RATE_MODEL", i, 1e-9, 2e-4);
}

TEST_CASE("GRM linear pulse vs analytic solution", "[GRM],[Simulation],[Analytic]")
{
	cadet::test::column::testAnalyticBenchmark("GENERAL_RATE_MODEL", "/data/grm-pulseBenchmark.data", true, true, 512, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("GENERAL_RATE_MODEL", "/data/grm-pulseBenchmark.data", true, false, 512, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("GENERAL_RATE_MODEL", "/data/grm-pulseBenchmark.data", false, true, 512, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("GENERAL_RATE_MODEL", "/data/grm-pulseBenchmark.data", false, false, 512, 6e-5, 1e-7);
}

TEST_CASE("GRM non-binding linear pulse vs analytic solution", "[GRM],[Simulation],[Analytic],[NonBinding]")
{
	cadet::test::column::testAnalyticNonBindingBenchmark("GENERAL_RATE_MODEL", "/data/grm-nonBinding.data", true, 512, 6e-5, 1e-7);
	cadet::test::column::testAnalyticNonBindingBenchmark("GENERAL_RATE_MODEL", "/data/grm-nonBinding.data", false, 512, 6e-5, 1e-7);
}

TEST_CASE("GRM Jacobian forward vs backward flow", "[GRM],[UnitOp],[Residual],[Jacobian],[AD]")
{
	// Test all WENO orders
	for (unsigned int i = 1; i <= cadet::Weno::maxOrder(); ++i)
		cadet::test::column::testJacobianWenoForwardBackward("GENERAL_RATE_MODEL", i);
}

TEST_CASE("GRM time derivative Jacobian vs FD", "[GRM],[UnitOp],[Residual],[Jacobian]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("GENERAL_RATE_MODEL", 1e-6, 0.0, 9e-4);
}

TEST_CASE("GRM sensitivity Jacobians", "[GRM],[UnitOp],[Sensitivity]")
{
	cadet::test::column::testFwdSensJacobians("GENERAL_RATE_MODEL", 1e-5);
}

TEST_CASE("GRM forward sensitivity vs FD", "[GRM],[Sensitivity],[Simulation]")
{
	// Relative error (5e-3, 8e-4) is checkd first, we use high absolute error (3e5, 2e-3)
	// for letting some very few points pass the error test, too. This is required due to
	// errors in finite differences.
	const double absTols[] = {3e5, 2e-3, 1.0, 5.0};
	const double relTols[] = {5e-3, 8e-4, 1e-5, 1e-4};
	const double passRatio[] = {0.97, 0.99, 0.999, 0.83};
	cadet::test::column::testFwdSensSolutionFD("GENERAL_RATE_MODEL", 5e-5, absTols, relTols, passRatio);
}

TEST_CASE("GRM forward sensitivity forward vs backward flow", "[GRM],[Sensitivity],[Simulation]")
{
	const double absTols[] = {4e-5, 1e-11, 1e-11, 2e-9};
	const double relTols[] = {6e-9, 1e-13, 1e-13, 3e-10};
	const double passRatio[] = {1.0, 0.99, 0.96, 1.0};
	cadet::test::column::testFwdSensSolutionForwardBackward("GENERAL_RATE_MODEL", absTols, relTols, passRatio);
}
