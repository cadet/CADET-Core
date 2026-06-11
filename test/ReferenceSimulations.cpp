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

TEST_CASE("Reference test: Frustum FV KOREN", "[Column_1D],[frustumFV_KOREN],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_radKOREN_eq_Z128.json");
	std::string refFilePath = std::string("/data/ref_radKOREN_eq_Z128.h5");
	const std::vector<double> absTol = { 5e-7 };
	const std::vector<double> relTol = { 1e-5 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: Frustum FV KOREN non-equidistant grid", "[Column_1D],[frustumFV_KOREN],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_radKOREN_nonEq_Z128.json");
	std::string refFilePath = std::string("/data/ref_radKOREN_nonEq_Z128.h5");
	const std::vector<double> absTol = { 5e-10 };
	const std::vector<double> relTol = { 1e-5 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: Frustum FV WENO2", "[Column_1D],[frustumFV_WENO2],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_radWENO2_eq_Z128.json");
	std::string refFilePath = std::string("/data/ref_radWENO2_eq_Z128.h5");
	const std::vector<double> absTol = { 5e-7 };
	const std::vector<double> relTol = { 1e-5 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: Frustum FV WENO2 non-equidistant grid", "[Column_1D],[frustumFV_WENO2],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_radWENO2_nonEq_Z128.json");
	std::string refFilePath = std::string("/data/ref_radWENO2_nonEq_Z128.h5");
	const std::vector<double> absTol = { 5e-7 };
	const std::vector<double> relTol = { 1e-5 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: Frustum FV WENO3", "[Column_1D],[frustumFV_WENO3],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_radWENO3_eq_Z128.json");
	std::string refFilePath = std::string("/data/ref_radWENO3_eq_Z128.h5");
	const std::vector<double> absTol = { 5e-10 };
	const std::vector<double> relTol = { 1e-5 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: Frustum FV WENO3 non-equidistant grid", "[Column_1D],[frustumFV_WENO3],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_radWENO3_nonEq_Z128.json");
	std::string refFilePath = std::string("/data/ref_radWENO3_nonEq_Z128.h5");
	const std::vector<double> absTol = { 5e-10 };
	const std::vector<double> relTol = { 1e-5 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: Radial FV KOREN", "[Column_1D],[radFV_KOREN],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_radKOREN_eq_Z128.json");
	std::string refFilePath = std::string("/data/ref_radKOREN_eq_Z128.h5");
	const std::vector<double> absTol = { 5e-10 };
	const std::vector<double> relTol = { 1e-4 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: Radial FV KOREN non-equidistant grid", "[Column_1D],[radFV_KOREN],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_radKOREN_nonEq_Z128.json");
	std::string refFilePath = std::string("/data/ref_radKOREN_nonEq_Z128.h5");
	const std::vector<double> absTol = { 5e-10 };
	const std::vector<double> relTol = { 1e-4 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: Radial FV WENO2", "[Column_1D],[radFV_WENO2],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_radWENO2_eq_Z128.json");
	std::string refFilePath = std::string("/data/ref_radWENO2_eq_Z128.h5");
	const std::vector<double> absTol = { 5e-10 };
	const std::vector<double> relTol = { 1e-4 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: Radial FV WENO2 non-equidistant grid", "[Column_1D],[radFV_WENO2],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_radWENO2_nonEq_Z128.json");
	std::string refFilePath = std::string("/data/ref_radWENO2_nonEq_Z128.h5");
	const std::vector<double> absTol = { 1e-7 };
	const std::vector<double> relTol = { 1e-4 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: Radial FV WENO3", "[Column_1D],[radFV_WENO3],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_radWENO3_eq_Z128.json");
	std::string refFilePath = std::string("/data/ref_radWENO3_eq_Z128.h5");
	const std::vector<double> absTol = { 5e-12 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: Radial FV WENO3 non-equidistant grid", "[Column_1D],[radFV_WENO3],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_radWENO3_nonEq_Z128.json");
	std::string refFilePath = std::string("/data/ref_radWENO3_nonEq_Z128.h5");
	const std::vector<double> absTol = { 5e-12 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: FV KOREN", "[Column_1D],[FV_KOREN],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_axKOREN_eq_Z128.json");
	std::string refFilePath = std::string("/data/ref_axKOREN_eq_Z128.h5");
	const std::vector<double> absTol = { 5e-10 };
	const std::vector<double> relTol = { 1e-4 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: FV KOREN non-equidistant grid", "[Column_1D],[FV_KOREN],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_axKOREN_nonEq_Z128.json");
	std::string refFilePath = std::string("/data/ref_axKOREN_nonEq_Z128.h5");
	const std::vector<double> absTol = { 5e-10 };
	const std::vector<double> relTol = { 1e-4 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: FV WENO2", "[Column_1D],[FV_WENO2],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_axWENO2_eq_Z128.json");
	std::string refFilePath = std::string("/data/ref_axWENO2_eq_Z128.h5");
	const std::vector<double> absTol = { 5e-10 };
	const std::vector<double> relTol = { 1e-4 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: FV WENO2 non-equidistant grid", "[Column_1D],[FV_WENO2],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_axWENO2_nonEq_Z128.json");
	std::string refFilePath = std::string("/data/ref_axWENO2_nonEq_Z128.h5");
	const std::vector<double> absTol = { 5e-12 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: FV WENO3", "[Column_1D],[FV_WENO3],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_axWENO3_eq_Z128.json");
	std::string refFilePath = std::string("/data/ref_axWENO3_eq_Z128.h5");
	const std::vector<double> absTol = { 5e-12 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: FV WENO3 non-equidistant grid", "[Column_1D],[FV_WENO3],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_axWENO3_nonEq_Z128.json");
	std::string refFilePath = std::string("/data/ref_axWENO3_nonEq_Z128.h5");
	const std::vector<double> absTol = { 5e-12 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

// MCT FV variant reference tests
TEST_CASE("Reference test: MCT FV KOREN", "[MCT],[FV_KOREN],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_MCT_KOREN_eq_Z128.json");
	std::string refFilePath = std::string("/data/ref_MCT_KOREN_eq_Z128.h5");
	const std::vector<double> absTol = { 5e-10 };
	const std::vector<double> relTol = { 1e-4 };

	cadet::test::column::FVParams disc(128);
	disc.setBulkDiscParam("NCHANNEL", 1);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: MCT FV KOREN non-equidistant grid", "[MCT],[FV_KOREN],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_MCT_KOREN_nonEq_Z128.json");
	std::string refFilePath = std::string("/data/ref_MCT_KOREN_nonEq_Z128.h5");
	const std::vector<double> absTol = { 5e-10 };
	const std::vector<double> relTol = { 1e-4 };

	cadet::test::column::FVParams disc(128);
	disc.setBulkDiscParam("NCHANNEL", 1);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: MCT FV WENO3 non-equidistant grid", "[MCT],[FV_WENO3],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_MCT_WENO3_nonEq_Z128.json");
	std::string refFilePath = std::string("/data/ref_MCT_WENO3_nonEq_Z128.h5");
	const std::vector<double> absTol = { 5e-12 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::FVParams disc(128);
	disc.setBulkDiscParam("NCHANNEL", 1);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

// 2DGRM FV variant reference tests
TEST_CASE("Reference test: 2DGRM FV KOREN", "[GRM2D],[FV_KOREN],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_2DGRM_KOREN_eq_Z128.json");
	std::string refFilePath = std::string("/data/ref_2DGRM_KOREN_eq_Z128.h5");
	const std::vector<double> absTol = { 5e-10 };
	const std::vector<double> relTol = { 1e-4 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: 2DGRM FV KOREN non-equidistant grid", "[GRM2D],[FV_KOREN],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_2DGRM_KOREN_nonEq_Z128.json");
	std::string refFilePath = std::string("/data/ref_2DGRM_KOREN_nonEq_Z128.h5");
	const std::vector<double> absTol = { 5e-10 };
	const std::vector<double> relTol = { 1e-4 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: 2DGRM FV WENO3 non-equidistant grid", "[GRM2D],[FV_WENO3],[Simulation],[CI],[numRef]")
{
	std::string modelFilePath = std::string("/data/config_2DGRM_WENO3_nonEq_Z128.json");
	std::string refFilePath = std::string("/data/ref_2DGRM_WENO3_nonEq_Z128.h5");
	const std::vector<double> absTol = { 5e-12 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::FVParams disc(128);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: spline binding single component", "[Column_1D],[SplineBinding],[Simulation],[CI]")
{
	std::string modelFilePath = std::string("/data/config_COL1D_GRM_spline_knots_shallow_7.json");
	std::string refFilePath = std::string("/data/ref_COL1D_GRM_spline_knots_shallow_7.h5");
	const std::vector<double> absTol = { 5e-7 };
	const std::vector<double> relTol = { 1e-5};

	cadet::test::column::FVParams disc(100, 15, 3);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: GPR binding with MLP kernel single component", "[Column_1D],[GPRBinding],[Simulation],[CI]")
{
	std::string modelFilePath = std::string("/data/config_COL1D_GRM_GPR_Shallow_MLP_7.json");
	std::string refFilePath = std::string("/data/ref_COL1D_GRM_GPR_Shallow_MLP_7.h5");
	const std::vector<double> absTol = { 5e-4 };
	const std::vector<double> relTol = { 1e-5 };

	cadet::test::column::FVParams disc(32, 8, 3);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: GPR binding with RBF kernel single component", "[Column_1D],[GPRBinding],[Simulation],[CI]")
{
	std::string modelFilePath = std::string("/data/config_COL1D_GRM_GPR_Shallow_RBF_15.json");
	std::string refFilePath = std::string("/data/ref_COL1D_GRM_GPR_Shallow_RBF_15.h5");
	const std::vector<double> absTol = { 1e-8 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::FVParams disc(32, 8, 3);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: GPR binding with MLP kernel for competitive Langmuir 2 component", "[Column_1D],[GPRBinding],[Simulation],[CI]")
{
	std::string modelFilePath = std::string("/data/config_Col1D_LRM_langmuirGPR_2comp_benchmark1.json");
	std::string refFilePath = std::string("/data/ref_Col1D_LRM_langmuirGPR_2comp_benchmark1.h5");
	const std::vector<double> absTol = { 1e-8 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::FVParams disc(32);
	disc.setBulkDiscParam("WENO_ORDER", static_cast<int>(2));
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: ANN binding for Langmuir single component", "[Column_1D],[ANNBinding],[Simulation],[CI]")
{
	std::string modelFilePath = std::string("/data/config_Col1D_GRM_langANN_1comp_benchmark1.json");
	std::string refFilePath = std::string("/data/ref_Col1D_GRM_langANN_1comp_benchmark1.h5");
	const std::vector<double> absTol = { 1e-8 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::FVParams disc(32, 8, 3);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: ANN binding for competitive Langmuir 2 component", "[Column_1D],[ANNBinding],[Simulation],[CI]")
{
	std::string modelFilePath = std::string("/data/config_Col1D_LRM_langANN_2comp_benchmark1.json");
	std::string refFilePath = std::string("/data/ref_Col1D_LRM_langANN_2comp_benchmark1.h5");
	const std::vector<double> absTol = { 1e-8 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::FVParams disc(32);
	disc.setBulkDiscParam("WENO_ORDER", static_cast<int>(2));
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: cyclic system with LRMP", "[Column_1D],[LRMP],[Simulation],[System],[CI]")
{
	std::string modelFilePath = std::string("/data/config_cyclicSystem1_LRMP_linBnd_1comp.json");
	std::string refFilePath = std::string("/data/ref_cyclicSystem1_LRMP_linBnd_1comp.h5");
	const std::vector<double> absTol = { 1e-10 };
	const std::vector<double> relTol = { 1e-8 };

	cadet::test::column::DGParams disc;
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "003", absTol, relTol, disc, false);
}

TEST_CASE("Reference test: acyclic system with LRMP", "[Column_1D],[LRMP],[Simulation],[System],[CI]")
{
	std::string modelFilePath = std::string("/data/config_acyclicSystem1_LRMP_linBnd_1comp.json");
	std::string refFilePath = std::string("/data/ref_acyclicSystem1_LRMP_linBnd_1comp.h5");
	const std::vector<double> absTol = { 1e-10 };
	const std::vector<double> relTol = { 1e-8 };

	cadet::test::column::DGParams disc;
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "006", absTol, relTol, disc, false);
}
