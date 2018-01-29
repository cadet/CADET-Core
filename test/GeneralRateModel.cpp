// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2017: The CADET Authors
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
		cadet::test::column::testWenoForwardBackward("GENERAL_RATE_MODEL", i, 1e-9, 1e-6);
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
