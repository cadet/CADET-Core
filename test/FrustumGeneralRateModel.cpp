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
#include "ParticleHelper.hpp"
#include "ReactionModelTests.hpp"
#include "JsonTestModels.hpp"
#include "Utils.hpp"


TEST_CASE("Frustum GRM with constant radius equals axial model numerical Benchmark with parameter sensitivities for linear case", "[FrustumGRM],[Simulation],[Reference],[Sensitivity],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_frustGRM_dynLin_1comp_sensbenchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_GRM_dynLin_1comp_sensbenchmark3_FVupwind_Z32parZ4.h5");
	const std::vector<double> absTol = { 1e-10, 1e-8, 1e-4, 1e-4 };
	const std::vector<double> relTol = { 1.0, 1.0, 1.0, 1.0 };
	cadet::test::column::FVParams disc(32, 4);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Frustum GRM transport Jacobian", "[FrustumGRM],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "FRUSTUM_GENERAL_RATE_MODEL", "FV");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Frustum GRM Jacobian forward vs backward flow", "[FrustumGRM],[UnitOp],[Residual],[Jacobian],[AD],[CI]")
{
	cadet::test::column::FVParams disc(16);
	cadet::test::column::testJacobianForwardBackward("FRUSTUM_GENERAL_RATE_MODEL", disc);
}

