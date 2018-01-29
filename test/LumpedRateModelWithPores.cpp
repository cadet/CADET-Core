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
	cadet::test::column::testFwdSensJacobians("LUMPED_RATE_MODEL_WITH_PORES", 1e-5);
}
