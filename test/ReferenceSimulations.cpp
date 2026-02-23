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

TEST_CASE("Reference test: FV KOREN", "[Column_1D],[FV_KOREN],[Simulation],[CI]")
{
	std::string modelFilePath = std::string("/data/config_axKOREN_eq_N128.json");
	std::string refFilePath = std::string("/data/ref_axKOREN_eq_N128.h5");
	const std::vector<double> absTol = { 5e-12 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: FV KOREN non-equidistant grid", "[Column_1D],[FV_KOREN],[Simulation],[CI]")
{
	std::string modelFilePath = std::string("/data/config_axKOREN_noneq_N128.json");
	std::string refFilePath = std::string("/data/ref_axKOREN_noneq_N128.h5");
	const std::vector<double> absTol = { 5e-12 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: FV WENO2", "[Column_1D],[FV_WENO2],[Simulation],[CI]")
{
	std::string modelFilePath = std::string("/data/config_axWENO2_eq_N128.json");
	std::string refFilePath = std::string("/data/ref_axWENO2_eq_N128.h5");
	const std::vector<double> absTol = { 5e-12 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: FV WENO2 non-equidistant grid", "[Column_1D],[FV_WENO2],[Simulation],[CI]")
{
	std::string modelFilePath = std::string("/data/config_axWENO2_noneq_N128.json");
	std::string refFilePath = std::string("/data/ref_axWENO2_noneq_N128.h5");
	const std::vector<double> absTol = { 5e-12 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: FV WENO3", "[Column_1D],[FV_WENO3],[Simulation],[CI]")
{
	std::string modelFilePath = std::string("/data/config_axWENO3_eq_N128.json");
	std::string refFilePath = std::string("/data/ref_axWENO3_eq_N128.h5");
	const std::vector<double> absTol = { 5e-12 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: FV WENO3 non-equidistant grid", "[Column_1D],[FV_WENO3],[Simulation],[CI]")
{
	std::string modelFilePath = std::string("/data/config_axWENO3_noneq_N128.json");
	std::string refFilePath = std::string("/data/ref_axWENO3_noneq_N128.h5");
	const std::vector<double> absTol = { 5e-12 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}