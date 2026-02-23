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

TEST_CASE("Reference test: spline binding single component", "[Column_1D],[SplineBinding],[Simulation],[CI]")
{
	std::string modelFilePath = std::string("/data/config_COL1D_GRM_spline_knots_shallow_7.json");
	std::string refFilePath = std::string("/data/ref_COL1D_GRM_spline_knots_shallow_7.h5");
	const std::vector<double> absTol = { 5e-7 };
	const std::vector<double> relTol = { 1e-5};

	cadet::test::column::FVParams disc(100, 15, 3);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}