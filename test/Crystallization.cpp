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
#include "ReactionModelTests.hpp"
#include "Weno.hpp"
#include "Utils.hpp"
#include "JsonTestModels.hpp"

/**
 * @brief Returns the absolute path to the test/ folder of the project
 * @details Absolute path to the test/ folder of the project without trailing slash
 * @return Absolute path to the test/ folder
 */
const char* getTestDirectory();


/*
 * Pure PBM tests
*/
TEST_CASE("Crystallization in a CSTR with initial distribution and growth", "[Crystallization],[PBM],[Simulation],[Reference],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_cry_CSTR_PBM_growth_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_cry_CSTR_PBM_growth_benchmark1.h5");
	const std::vector<double> absTol = { 1e-10 };
	const std::vector<double> relTol = { 1e-10 };

	cadet::test::column::Dummyparams disc; // CSTR, so no spatial resolution
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Crystallization in a CSTR with initial distribution and size-dependent growth", "[Crystallization],[PBM],[Simulation],[Reference],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_cry_CSTR_PBM_growthSizeDep_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_cry_CSTR_PBM_growthSizeDep_benchmark1.h5");
	const std::vector<double> absTol = { 1e-10 };
	const std::vector<double> relTol = { 1e-10 };

	cadet::test::column::Dummyparams disc; // CSTR, so no spatial resolution
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Crystallization in a CSTR with primary nucleation and growth", "[Crystallization],[PBM],[Simulation],[Reference],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_cry_CSTR_PBM_primaryNucleationAndGrowth_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_cry_CSTR_PBM_primaryNucleationAndGrowth_benchmark1.h5");
	const std::vector<double> absTol = { 1e-10 };
	const std::vector<double> relTol = { 1e-10 };

	cadet::test::column::Dummyparams disc; // CSTR, so no spatial resolution
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Crystallization in a CSTR with primary nucleation, growth and growth rate dispersion", "[Crystallization],[PBM],[Simulation],[Reference],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_cry_CSTR_PBM_primaryNucleationGrowthGrowthRateDispersion_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_cry_CSTR_PBM_primaryNucleationGrowthGrowthRateDispersion_benchmark1.h5");
	const std::vector<double> absTol = { 1e-10 };
	const std::vector<double> relTol = { 1e-10 };

	cadet::test::column::Dummyparams disc; // CSTR, so no spatial resolution
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Crystallization in a CSTR with primary and secondary nucleation and growth", "[Crystallization],[PBM],[Simulation],[Reference],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_cry_CSTR_PBM_primarySecondaryNucleationAndGrowth_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_cry_CSTR_PBM_primarySecondaryNucleationAndGrowth_benchmark1.h5");
	const std::vector<double> absTol = { 5e+6 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::Dummyparams disc; // CSTR, so no spatial resolution
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Crystallization in a DPFR/LRM with primary and secondary nucleation and growth", "[CrysToFix0],[Crystallization],[PBM],[Simulation],[Reference],[CI]") 
{
	const std::string& modelFilePath = std::string("/data/model_cry_DPFR_PBM_primarySecondaryNucleationGrowth_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_cry_DPFR_PBM_primarySecondaryNucleationGrowth_benchmark1.h5");
	const std::vector<double> absTol = { 2e+8 };
	const std::vector<double> relTol = { 5e-6 };

	cadet::test::column::FVparams disc(25);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Crystallization Jacobian verification for a CSTR with initial distribution and growth", "[Crystallization],[PBM],[UnitOp],[Jacobian],[CI]")
{
	// read json model setup file
	const std::string& modelFileRelPath = std::string("/data/model_cry_CSTR_PBM_growth_benchmark1.json");
	const std::string setupFile = std::string(getTestDirectory()) + modelFileRelPath;
	cadet::JsonParameterProvider pp_setup(cadet::JsonParameterProvider::fromFile(setupFile));

	pp_setup.pushScope("model");
	pp_setup.pushScope("unit_001");

	// set the discretization scheme
	pp_setup.pushScope("reaction_bulk");
	pp_setup.set("CRY_GROWTH_SCHEME_ORDER", 4);
	pp_setup.popScope();

	cadet::test::column::testJacobianAD(pp_setup, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Crystallization Jacobian verification for a CSTR with initial distribution and size-dependent growth", "[Crystallization],[PBM],[UnitOp],[Jacobian],[CI]")
{
	// read json model setup file
	const std::string& modelFileRelPath = std::string("/data/model_cry_CSTR_PBM_growthSizeDep_benchmark1.json");
	const std::string setupFile = std::string(getTestDirectory()) + modelFileRelPath;
	cadet::JsonParameterProvider pp_setup(cadet::JsonParameterProvider::fromFile(setupFile));

	pp_setup.pushScope("model");
	pp_setup.pushScope("unit_001");

	// set the discretization scheme
	pp_setup.pushScope("reaction_bulk");
	pp_setup.set("CRY_GROWTH_SCHEME_ORDER", 3);
	pp_setup.popScope();

	cadet::test::column::testJacobianAD(pp_setup, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Crystallization Jacobian verification for a CSTR with primary nucleation and growth", "[Crystallization],[PBM],[UnitOp],[Jacobian],[CI]")
{
	// read json model setup file
	const std::string& modelFileRelPath = std::string("/data/model_cry_CSTR_PBM_primaryNucleationAndGrowth_benchmark1.json");
	const std::string setupFile = std::string(getTestDirectory()) + modelFileRelPath;
	cadet::JsonParameterProvider pp_setup(cadet::JsonParameterProvider::fromFile(setupFile));

	pp_setup.pushScope("model");
	pp_setup.pushScope("unit_001");

	// set the discretization scheme
	pp_setup.pushScope("reaction_bulk");
	pp_setup.set("CRY_GROWTH_SCHEME_ORDER", 1);
	pp_setup.popScope();

	cadet::test::column::testJacobianAD(pp_setup, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Crystallization Jacobian verification for a CSTR with primary nucleation, growth and growth rate dispersion", "[Crystallization],[PBM],[UnitOp],[Jacobian],[CI]")
{
	// read json model setup file
	const std::string& modelFileRelPath = std::string("/data/model_cry_CSTR_PBM_primaryNucleationGrowthGrowthRateDispersion_benchmark1.json");
	const std::string setupFile = std::string(getTestDirectory()) + modelFileRelPath;
	cadet::JsonParameterProvider pp_setup(cadet::JsonParameterProvider::fromFile(setupFile));

	pp_setup.pushScope("model");
	pp_setup.pushScope("unit_001");

	// set the discretization scheme
	pp_setup.pushScope("reaction_bulk");
	pp_setup.set("CRY_GROWTH_SCHEME_ORDER", 1);
	pp_setup.popScope();

	cadet::test::column::testJacobianAD(pp_setup, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Crystallization Jacobian verification for a CSTR with primary and secondary nucleation and growth", "[Crystallization],[PBM],[UnitOp],[Jacobian],[CI]")
{
	// read json model setup file
	const std::string& modelFileRelPath = std::string("/data/model_cry_CSTR_PBM_primarySecondaryNucleationAndGrowth_benchmark1.json");
	const std::string setupFile = std::string(getTestDirectory()) + modelFileRelPath;
	cadet::JsonParameterProvider pp_setup(cadet::JsonParameterProvider::fromFile(setupFile));

	pp_setup.pushScope("model");
	pp_setup.pushScope("unit_001");

	// set the discretization scheme
	pp_setup.pushScope("reaction_bulk");
	pp_setup.set("CRY_GROWTH_SCHEME_ORDER", 3);
	pp_setup.popScope();

	cadet::test::column::testJacobianAD(pp_setup, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Crystallization Jacobian verification for a DPFR/LRM with primary and secondary nucleation and growth", "[Crystallization1],[UnitOp],[Jacobian],[CI]")
{
	// read json model setup file
	const std::string& modelFileRelPath = std::string("/data/model_cry_DPFR_PBM_primarySecondaryNucleationGrowth_benchmark1.json");
	const std::string setupFile = std::string(getTestDirectory()) + modelFileRelPath;
	cadet::JsonParameterProvider pp_setup(cadet::JsonParameterProvider::fromFile(setupFile));

	pp_setup.pushScope("model");
	pp_setup.pushScope("unit_001");

	// for this specific test, we need to define a (high) tolerances as the values in this test are numerically very challenging (values of ca. 1E+24)
	const double ADabsTol = 5e+8;
	const double FDabsTol = 1e+10;

	cadet::test::column::testJacobianAD(pp_setup, FDabsTol, ADabsTol);
}

/*
 * Pure aggregation tests
*/
TEST_CASE("Crystallization Aggregation in a CSTR", "[Crystallization],[Aggregation],[Simulation],[Reference],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_cry_CSTR_aggregation_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_cry_CSTR_aggregation_benchmark1.h5");
	const std::vector<double> absTol = { 1e-10 };
	const std::vector<double> relTol = { 1e-10 };

	cadet::test::column::Dummyparams disc; // CSTR, so no spatial resolution
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Crystallization Jacobian verification for Aggregation in a CSTR", "[Crystallization],[Aggregation],[UnitOp],[Jacobian],[CI]")
{
	// read json model setup file
	const std::string& modelFileRelPath = std::string("/data/model_cry_CSTR_aggregation_benchmark1.json");
	const std::string setupFile = std::string(getTestDirectory()) + modelFileRelPath;
	cadet::JsonParameterProvider pp_setup(cadet::JsonParameterProvider::fromFile(setupFile));

	pp_setup.pushScope("model");
	pp_setup.pushScope("unit_001");

	const double ADabsTol = 1E-4;
	const double FDabsTol = 1E+2;

	cadet::test::column::testJacobianAD(pp_setup, FDabsTol, ADabsTol);
}

/*
 * Pure fragmentation/breakage tests
*/
TEST_CASE("Crystallization Fragmentation in a CSTR", "[Crystallization],[Fragmentation],[Simulation],[Reference],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_cry_CSTR_fragmentation_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_cry_CSTR_fragmentation_benchmark1.h5");
	const std::vector<double> absTol = { 1e-7 };
	const std::vector<double> relTol = { 1e-7 };

	cadet::test::column::Dummyparams disc; // CSTR, so no spatial resolution
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Crystallization Jacobian verification for Fragmentation in a CSTR", "[Crystallization],[Fragmentation],[UnitOp],[Jacobian],[CI]")
{
	// read json model setup file
	const std::string& modelFileRelPath = std::string("/data/model_cry_CSTR_fragmentation_benchmark1.json");
	const std::string setupFile = std::string(getTestDirectory()) + modelFileRelPath;
	cadet::JsonParameterProvider pp_setup(cadet::JsonParameterProvider::fromFile(setupFile));

	pp_setup.pushScope("model");
	pp_setup.pushScope("unit_001");

	const double ADabsTol = 1E-4;
	const double FDabsTol = 1E+2;

	cadet::test::column::testJacobianAD(pp_setup, FDabsTol, ADabsTol);
}

/*
 * Combined fragmentation, breakage tests
*/
TEST_CASE("Crystallization combined Aggregation and Fragmentation in a CSTR", "[Crystallization],[AggFrag],[Simulation],[Reference],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_cry_CSTR_aggFrag_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_cry_CSTR_aggFrag_benchmark1.h5");
	const std::vector<double> absTol = { 1e-8 };
	const std::vector<double> relTol = { 1e-8 };

	cadet::test::column::Dummyparams disc; // CSTR, so no spatial resolution
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Crystallization Jacobian verification for combined Aggregation and Fragmentation in a CSTR", "[Crystallization],[AggFrag],[UnitOp],[Jacobian],[CI]")
{
	// read json model setup file
	const std::string& modelFileRelPath = std::string("/data/model_cry_CSTR_aggFrag_benchmark1.json");
	const std::string setupFile = std::string(getTestDirectory()) + modelFileRelPath;
	cadet::JsonParameterProvider pp_setup(cadet::JsonParameterProvider::fromFile(setupFile));

	pp_setup.pushScope("model");
	pp_setup.pushScope("unit_001");

	const double ADabsTol = 1E-4;
	const double FDabsTol = 1E+2;

	cadet::test::column::testJacobianAD(pp_setup, FDabsTol, ADabsTol);
}

/*
 * Combined PBM, fragmentation, breakage tests
*/
TEST_CASE("Crystallization combined PBM, Aggregation and Fragmentation in a CSTR", "[Crystallization],[PBMAggFrag],[Simulation],[Reference],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_cry_CSTR_PBM_Agg_Frag_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_cry_CSTR_PBM_Agg_Frag_benchmark1.h5");
	const std::vector<double> absTol = { 2e-7 };
	const std::vector<double> relTol = { 8e-1 };

	cadet::test::column::Dummyparams disc; // CSTR, so no spatial resolution
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Crystallization Jacobian verification for combined PBM, Aggregation and Fragmentation in a CSTR", "[Crystallization],[PBMAggFrag],[UnitOp],[Jacobian],[CI]")
{
	// read json model setup file
	const std::string& modelFileRelPath = std::string("/data/model_cry_CSTR_PBM_Agg_Frag_benchmark1.json");
	const std::string setupFile = std::string(getTestDirectory()) + modelFileRelPath;
	cadet::JsonParameterProvider pp_setup(cadet::JsonParameterProvider::fromFile(setupFile));

	pp_setup.pushScope("model");
	pp_setup.pushScope("unit_001");

	const double ADabsTol = 1E-4;
	const double FDabsTol = 1E+2;

	cadet::test::column::testJacobianAD(pp_setup, FDabsTol, ADabsTol);
}

/*
 * Combined PBM, aggregation test in a DPFR
*/
TEST_CASE("Crystallization combined PBM and Aggregation in a DPFR", "[Crystallization],[DPFR_PBMAgg],[Simulation],[Reference],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_cry_DPFR_PBM_aggregation_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_cry_DPFR_PBM_aggregation_benchmark1.h5");
	const std::vector<double> absTol = { 1E+10 }; // we need to define a (high) tolerances as the numerical values in this test are extremely high, values of up to xE+27
	const std::vector<double> relTol = { 1E-12 };

	cadet::test::column::FVparams disc(16);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Crystallization Jacobian verification for combined PBM and Aggregation in a DPFR", "[Crystallization],[DPFR_PBMAgg],[UnitOp],[Jacobian],[CI]")
{
	// read json model setup file
	const std::string& modelFileRelPath = std::string("/data/model_cry_DPFR_PBM_aggregation_benchmark1.json");
	const std::string setupFile = std::string(getTestDirectory()) + modelFileRelPath;
	cadet::JsonParameterProvider pp_setup(cadet::JsonParameterProvider::fromFile(setupFile));

	pp_setup.pushScope("model");
	pp_setup.pushScope("unit_001");

	// reduce axial discretization to stay within the allowed number of AD directions
	pp_setup.pushScope("discretization");
	pp_setup.set("NCOL", 4); // 4 * 34 < 157
	pp_setup.popScope();

	// for this specific test, we need to define a (high) tolerances as the numerical values in this test are extremely high, values of up to xE+27
	const double ADabsTol = 1E+11;
	const double FDabsTol = 1E+10;

	cadet::test::column::testJacobianAD(pp_setup, FDabsTol, ADabsTol);
}
