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

TEST_CASE("LRM linear pulse vs analytic solution", "[LRM],[Simulation],[Analytic]")
{
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITHOUT_PORES", "/data/lrm-pulseBenchmark.data", true, true, 1024, 2e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITHOUT_PORES", "/data/lrm-pulseBenchmark.data", true, false, 1024, 2e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITHOUT_PORES", "/data/lrm-pulseBenchmark.data", false, true, 1024, 2e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITHOUT_PORES", "/data/lrm-pulseBenchmark.data", false, false, 1024, 2e-5, 1e-7);
}

TEST_CASE("LRM non-binding linear pulse vs analytic solution", "[LRM],[Simulation],[Analytic],[NonBinding]")
{
	cadet::test::column::testAnalyticNonBindingBenchmark("LUMPED_RATE_MODEL_WITHOUT_PORES", "/data/lrm-nonBinding.data", true, 1024, 2e-5, 1e-7);
	cadet::test::column::testAnalyticNonBindingBenchmark("LUMPED_RATE_MODEL_WITHOUT_PORES", "/data/lrm-nonBinding.data", false, 1024, 2e-5, 1e-7);
}

TEST_CASE("LRM Jacobian forward vs backward flow", "[LRM],[UnitOp],[Residual],[Jacobian],[AD]")
{
	// Test all WENO orders
	for (unsigned int i = 1; i <= cadet::Weno::maxOrder(); ++i)
		cadet::test::column::testJacobianWenoForwardBackward("LUMPED_RATE_MODEL_WITHOUT_PORES", i);
}

TEST_CASE("LRM time derivative Jacobian vs FD", "[LRM],[UnitOp],[Residual],[Jacobian]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("LUMPED_RATE_MODEL_WITHOUT_PORES");
}

TEST_CASE("LRM sensitivity Jacobians", "[LRM],[UnitOp],[Sensitivity]")
{
	cadet::test::column::testFwdSensJacobians("LUMPED_RATE_MODEL_WITHOUT_PORES", 1e-5);
}

