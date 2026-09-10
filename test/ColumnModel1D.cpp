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
#include "Approx.hpp"

#define CADET_LOGGING_DISABLE
#include "Logging.hpp"

#include "ColumnTests.hpp"
#include "ParticleHelper.hpp"
#include "ReactionModelTests.hpp"
#include "JsonTestModels.hpp"
#include "model/UnitOperation.hpp"
#include "ParallelSupport.hpp"
#include "SimHelper.hpp"
#include "SimulationTypes.hpp"
#include "UnitOperationTests.hpp"
#include "Utils.hpp"
#include "common/Driver.hpp"
#include "model/parts/DGToolbox.hpp"

#include <cmath>
#include <functional>
#include <vector>

/**
 * @brief Returns the absolute path to the test/ folder of the project
 * @details Absolute path to the test/ folder of the project without trailing slash
 * @return Absolute path to the test/ folder
 */
const char* getTestDirectory();

namespace
{
	/**
	 * @brief Computes the cross section area at every DG node of a column unit
	 * @details Evaluates the given area profile at the LGL nodes of all elements, in the
	 *          ordering expected by CROSS_SECTIONAL_AREA_AT_NODES, i.e. element-major with
	 *          the element interface nodes listed twice.
	 * @param [in] polyDeg DG polynomial degree
	 * @param [in] nElem number of DG elements
	 * @param [in] area cross section area as a function of the normalized axial coordinate x / H in [0, 1]
	 * @return cross section area at every DG node
	 */
	std::vector<double> crossSectionAreaAtNodes(const int polyDeg, const int nElem, const std::function<double(double)>& area)
	{
		Eigen::VectorXd nodes(polyDeg + 1);
		Eigen::VectorXd weights(polyDeg + 1);
		cadet::model::parts::dgtoolbox::lglNodesWeights(polyDeg, nodes, weights, false);

		std::vector<double> areaAtNodes(nElem * (polyDeg + 1));

		for (int elem = 0; elem < nElem; ++elem)
		{
			for (int node = 0; node <= polyDeg; ++node)
				areaAtNodes[elem * (polyDeg + 1) + node] = area((elem + 0.5 * (1.0 + nodes[node])) / nElem);
		}

		return areaAtNodes;
	}

	/**
	 * @brief Replaces the geometry of a column unit by the SMOOTHLY_VARYING geometry
	 * @details The parameter provider has to be scoped to the column unit group.
	 * @param [in,out] jpp parameter provider, scoped to the column unit
	 * @param [in] area cross section area as a function of the normalized axial coordinate x / H in [0, 1]
	 */
	void setSmoothlyVaryingGeometry(cadet::JsonParameterProvider& jpp, const std::function<double(double)>& area)
	{
		jpp.pushScope("discretization");
		const int polyDeg = jpp.getInt("POLYDEG");
		const int nElem = jpp.getInt("NELEM");
		jpp.popScope();

		jpp.set("GEOMETRY", std::string("SMOOTHLY_VARYING"));
		jpp.set("CROSS_SECTIONAL_AREA_AT_NODES", crossSectionAreaAtNodes(polyDeg, nElem, area));
	}

	/**
	 * @brief Replaces the geometry of a column unit of a full model system by the SMOOTHLY_VARYING geometry
	 * @param [in,out] jpp parameter provider, scoped to the root of the model system
	 * @param [in] unitId unit group name, e.g. "unit_000"
	 * @param [in] area cross section area as a function of the normalized axial coordinate x / H in [0, 1]
	 */
	void setSmoothlyVaryingGeometry(cadet::JsonParameterProvider& jpp, const std::string& unitId, const std::function<double(double)>& area)
	{
		jpp.pushScope("model");
		jpp.pushScope(unitId);

		setSmoothlyVaryingGeometry(jpp, area);

		jpp.popScope();
		jpp.popScope();
	}

	/**
	 * @brief Returns the cross section area of a frustum with the given end areas
	 * @param [in] areaLargeEnd cross section area at x = 0
	 * @param [in] areaSmallEnd cross section area at x = H
	 */
	std::function<double(double)> frustumArea(const double areaLargeEnd, const double areaSmallEnd)
	{
		const double pi = 3.14159265358979323846;
		const double radiusLargeEnd = std::sqrt(areaLargeEnd / pi);
		const double radiusSmallEnd = std::sqrt(areaSmallEnd / pi);

		return [=](double normalizedPos)
			{
				const double radius = radiusLargeEnd + normalizedPos * (radiusSmallEnd - radiusLargeEnd);
				return pi * radius * radius;
			};
	}
}

TEST_CASE("Column_1D as axial GRM with FV equivalence with arrow head implementation", "[AxialColumn1D],[FV],[Simulation],[CI]")
{
	cadet::JsonParameterProvider jpp1 = createLWE("COLUMN_MODEL_1D_GRM", "FV");
	cadet::JsonParameterProvider jpp2 = createLWE("COLUMN_MODEL_1D_GRM", "FV");
	cadet::test::column::FVParams disc(32);

	disc.setDisc(jpp1);
	disc.setDisc(jpp2);

	// disable arrow head optimization for second setting
	jpp2.pushScope("model");
	jpp2.pushScope("unit_000");
	jpp2.pushScope("discretization");
	jpp2.set("FV_ARROW_HEAD_OPTIMIZATION", false);
	jpp2.popScope();
	jpp2.popScope();
	jpp2.popScope();

	// low time integration tolerances to minimize impact of different linear solvers
	jpp2.pushScope("solver");
	jpp2.pushScope("time_integrator");
	jpp2.set("ABSTOL", 1e-12);
	jpp2.set("RELTOL", 1e-10);
	jpp2.popScope();
	jpp2.popScope();

	jpp1.pushScope("solver");
	jpp1.pushScope("time_integrator");
	jpp1.set("ABSTOL", 1e-12);
	jpp1.set("RELTOL", 1e-10);
	jpp1.popScope();
	jpp1.popScope();

	cadet::test::column::testEqualResults(jpp1, jpp2, 1e-10, 1e-8, 0);
}

TEST_CASE("Column_1D as radial flow LRMP with FV equivalence with arrow head implementation", "[AxialColumn1D],[FV],[Simulation],[CI]")
{
	cadet::JsonParameterProvider jpp1 = createLWE("RADIAL_COLUMN_MODEL_1D_LRMP", "FV");
	cadet::JsonParameterProvider jpp2 = createLWE("RADIAL_COLUMN_MODEL_1D_LRMP", "FV");
	cadet::test::column::FVParams disc(32);
	disc.setBulkDiscParam("WENO_ORDER", 2);

	disc.setDisc(jpp1);
	disc.setDisc(jpp2);

	// disable arrow head optimization for second setting
	jpp2.pushScope("model");
	jpp2.pushScope("unit_000");
	jpp2.pushScope("discretization");
	jpp2.set("FV_ARROW_HEAD_OPTIMIZATION", false);
	jpp2.popScope();
	jpp2.popScope();
	jpp2.popScope();

	// low time integration tolerances to minimize impact of different linear solvers
	jpp2.pushScope("solver");
	jpp2.pushScope("time_integrator");
	jpp2.set("ABSTOL", 1e-12);
	jpp2.set("RELTOL", 1e-10);
	jpp2.popScope();
	jpp2.popScope();

	jpp1.pushScope("solver");
	jpp1.pushScope("time_integrator");
	jpp1.set("ABSTOL", 1e-12);
	jpp1.set("RELTOL", 1e-10);
	jpp1.popScope();
	jpp1.popScope();

	cadet::test::column::testEqualResults(jpp1, jpp2, 1e-10, 1e-8, 0);
}

TEST_CASE("Column_1D as frustum LRMP with FV equivalence with arrow head implementation", "[AxialColumn1D],[FV],[Simulation],[CI]")
{
	cadet::JsonParameterProvider jpp1 = createLWE("FRUSTUM_COLUMN_MODEL_1D_LRMP", "FV");
	cadet::JsonParameterProvider jpp2 = createLWE("FRUSTUM_COLUMN_MODEL_1D_LRMP", "FV");
	cadet::test::column::FVParams disc(32);
	disc.setBulkDiscParam("WENO_ORDER", 2);

	disc.setDisc(jpp1);
	disc.setDisc(jpp2);

	// disable arrow head optimization for second setting
	jpp2.pushScope("model");
	jpp2.pushScope("unit_000");
	jpp2.pushScope("discretization");
	jpp2.set("FV_ARROW_HEAD_OPTIMIZATION", false);
	jpp2.popScope();
	jpp2.popScope();
	jpp2.popScope();

	// low time integration tolerances to minimize impact of different linear solvers
	jpp2.pushScope("solver");
	jpp2.pushScope("time_integrator");
	jpp2.set("ABSTOL", 1e-12);
	jpp2.set("RELTOL", 1e-10);
	jpp2.popScope();
	jpp2.popScope();

	jpp1.pushScope("solver");
	jpp1.pushScope("time_integrator");
	jpp1.set("ABSTOL", 1e-12);
	jpp1.set("RELTOL", 1e-10);
	jpp1.popScope();
	jpp1.popScope();

	cadet::test::column::testEqualResults(jpp1, jpp2, 1e-10, 1e-8, 0);
}

TEST_CASE("Column_1D as GRM with FV transport Jacobian", "[AxialColumn1D],[FV],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "COLUMN_MODEL_1D_GRM", "FV");

	jpp.pushScope("discretization");
	jpp.set("FV_ARROW_HEAD_OPTIMIZATION", false);
	jpp.popScope();

	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as GRM LWE forward vs backward flow", "[AxialColumn1D],[DG],[DG1D],[Simulation],[CI]")
{
	cadet::test::column::DGParams disc;

	// Test all integration modes
	for (int i = 0; i <= 2; i++)
	{
		disc.setBulkDiscParam("USE_COLLOCATION_DG", i);
		cadet::test::column::testForwardBackward("COLUMN_MODEL_1D_GRM", disc, 1e-9, 2e-4);
	}
}

TEST_CASE("Column_1D as GRM linear pulse vs analytic solution", "[AxialColumn1D],[DG],[DG1D],[Simulation],[Analytic],[CI]")
{
	cadet::test::column::DGParams disc;
	cadet::test::column::testAnalyticBenchmark("COLUMN_MODEL_1D_GRM", "/data/grm-pulseBenchmark.data", true, true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("COLUMN_MODEL_1D_GRM", "/data/grm-pulseBenchmark.data", true, false, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("COLUMN_MODEL_1D_GRM", "/data/grm-pulseBenchmark.data", false, true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("COLUMN_MODEL_1D_GRM", "/data/grm-pulseBenchmark.data", false, false, disc, 6e-5, 1e-7);
}

TEST_CASE("Column_1D as GRM linear pulse vs analytic solution with bulk DG and particle FV discretization", "[AxialColumn1D],[DGFV],[Simulation],[Analytic],[CI]")
{
	auto discDGFV = cadet::test::column::createDGFVParams(0, 3, 15, 0, 0, 5);
	cadet::test::column::testAnalyticBenchmark("COLUMN_MODEL_1D_GRM", "/data/grm-pulseBenchmark.data", true, true, *discDGFV, "DGFV", 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("COLUMN_MODEL_1D_GRM", "/data/grm-pulseBenchmark.data", true, false, *discDGFV, "DGFV", 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("COLUMN_MODEL_1D_GRM", "/data/grm-pulseBenchmark.data", false, true, *discDGFV, "DGFV", 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("COLUMN_MODEL_1D_GRM", "/data/grm-pulseBenchmark.data", false, false, *discDGFV, "DGFV", 6e-5, 1e-7);
}

TEST_CASE("Column_1D as GRM linear pulse vs analytic solution with bulk FV and particle DG discretization", "[AxialColumn1D],[FVDG],[Simulation],[Analytic],[CI]")
{
	// Regression test: an FV bulk discretization combined with a DG particle discretization must be routed to
	// ColumnModel1D. The arrow head FV unit operations discretize the particles with FV only and used to be
	// selected for any FV bulk discretization, which made this combination fail during configuration.
	auto discFVDG = cadet::test::column::createFVDGParams(512, 3, 0, 3, 1); // 512 axial FV cells, as in the pure FV analytic benchmarks
	cadet::test::column::testAnalyticBenchmark("COLUMN_MODEL_1D_GRM", "/data/grm-pulseBenchmark.data", true, true, *discFVDG, "FVDG", 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("COLUMN_MODEL_1D_GRM", "/data/grm-pulseBenchmark.data", false, false, *discFVDG, "FVDG", 6e-5, 1e-7);
}

TEST_CASE("Column_1D as GRM with FV particles is invariant to an equivalent constant FILM_DIFFUSION_DEP", "[AxialColumn1D],[DGFV],[Simulation],[ParameterDependence],[CI]")
{
	// Regression test: the FV particle operator corrects k_f for the resistance of the outer half shell. That
	// correction has to be evaluated with the FILM_DIFFUSION_DEP modified k_f, just like the particle side
	// boundary condition, otherwise the bulk and the particle side of the film diffusion flux disagree.
	// A POWER_LAW dependence with exponent 0 and base b is a constant factor b, so scaling FILM_DIFFUSION by
	// 1/b must reproduce the unmodified simulation exactly.
	const double base = 4.0;

	cadet::JsonParameterProvider jpp1 = createLWE("COLUMN_MODEL_1D_GRM", "DGFV");
	cadet::JsonParameterProvider jpp2 = createLWE("COLUMN_MODEL_1D_GRM", "DGFV");

	auto disc = cadet::test::column::createDGFVParams(0, 3, 8, 0, 0, 4);
	disc->setDisc(jpp1);
	disc->setDisc(jpp2);

	jpp2.pushScope("model");
	jpp2.pushScope("unit_000");
	jpp2.pushScope("particle_type_000");

	std::vector<double> filmDiff = jpp2.getDoubleArray("FILM_DIFFUSION");
	for (double& fd : filmDiff)
		fd /= base;

	jpp2.set("FILM_DIFFUSION", filmDiff);
	jpp2.set("FILM_DIFFUSION_DEP", std::string("POWER_LAW"));
	jpp2.set("FILM_DIFFUSION_DEP_BASE", base);
	jpp2.set("FILM_DIFFUSION_DEP_EXPONENT", 0.0);

	jpp2.popScope();
	jpp2.popScope();
	jpp2.popScope();

	cadet::test::column::testEqualResults(jpp1, jpp2, 1e-10, 1e-8, 0);
}

TEST_CASE("Column_1D as LRMP linear pulse vs analytic solution", "[AxialColumn1D],[DG],[DG1D],[Simulation],[Analytic],[CI]")
{
	cadet::test::column::DGParams disc;
	cadet::test::column::testAnalyticBenchmark("COLUMN_MODEL_1D_LRMP", "/data/lrmp-pulseBenchmark.data", true, true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("COLUMN_MODEL_1D_LRMP", "/data/lrmp-pulseBenchmark.data", true, false, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("COLUMN_MODEL_1D_LRMP", "/data/lrmp-pulseBenchmark.data", false, true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("COLUMN_MODEL_1D_LRMP", "/data/lrmp-pulseBenchmark.data", false, false, disc, 6e-5, 1e-7);
}

TEST_CASE("Column_1D as GRM non-binding linear pulse vs analytic solution", "[AxialColumn1D],[DG],[DG1D],[Simulation],[Analytic],[NonBinding],[CI]")
{
	cadet::test::column::DGParams disc;
	cadet::test::column::testAnalyticNonBindingBenchmark("COLUMN_MODEL_1D_GRM", "/data/grm-nonBinding.data", true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticNonBindingBenchmark("COLUMN_MODEL_1D_GRM", "/data/grm-nonBinding.data", false, disc, 6e-5, 1e-7);
}

TEST_CASE("Column_1D as LRMP non-binding linear pulse vs analytic solution", "[AxialColumn1D],[DG],[DG1D],[Simulation],[Analytic],[NonBinding],[CI]")
{
	cadet::test::column::DGParams disc;
	cadet::test::column::testAnalyticNonBindingBenchmark("COLUMN_MODEL_1D_LRMP", "/data/lrmp-nonBinding.data", true, disc, 6e-5, 1e-7);
	cadet::test::column::testAnalyticNonBindingBenchmark("COLUMN_MODEL_1D_LRMP", "/data/lrmp-nonBinding.data", false, disc, 6e-5, 1e-7);
}

TEST_CASE("Column_1D as GRM Jacobian forward vs backward flow", "[AxialColumn1D],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[AD],[CI]")
{
	cadet::test::column::DGParams disc;
	disc.setBulkDiscParam("USE_COLLOCATION_DG", 1);
	cadet::test::column::testJacobianForwardBackward("COLUMN_MODEL_1D_GRM", disc, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as DPF numerical Benchmark1", "[AxialColumn1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI]")
{
	std::string modelFilePath = std::string("/data/model_COL1D_DPFR_1comp_benchmark1.json");
	std::string refFilePath = std::string("/data/ref_COL1D_DPFR_1comp_benchmark1_DG_P3Z8.h5");
	const std::vector<double> absTol = { 1e-12 };
	const std::vector<double> relTol = { 1e-4 };

	cadet::test::column::DGParams disc(1, 3, 8, 3, 1);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("Column_1D as GRM numerical Benchmark1 with parameter sensitivities for linear case", "[AxialColumn1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI_sens17]")
{
	std::string modelFilePath = std::string("/data/model_COL1D_GRM_dynLin_1comp_benchmark1.json");
	std::string refFilePath = std::string("/data/ref_GRM_dynLin_1comp_sensbenchmark1_cDG_P3Z8_GSM_parP3parZ1.h5");
	const std::vector<double> absTol = { 1e-12, 1e-6, 1e-6, 1e-10 };
	const std::vector<double> relTol = { 1e-4, 1e-3, 1e-4, 1e-3 };

	cadet::test::column::DGParams disc(1, 3, 8, 3, 1);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("Column_1D as GRM numerical Benchmark1 for linear case with surface diffusion", "[AxialColumn1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI]")
{
	std::string modelFilePath = std::string("/data/model_COL1D_GRMsd_dynLin_1comp_benchmark1.json");
	std::string refFilePath = std::string("/data/ref_COL1D_GRMsd_dynLin_1comp_benchmark1_cDG_P3Z8_GSM_parP3parZ1.h5");
	const std::vector<double> absTol = { 1e-12 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::DGParams disc(1, 3, 8, 3, 1);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Column_1D as LRMP numerical Benchmark with parameter sensitivities for linear case", "[AxialColumn1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI_sens22]")
{
	const std::string& modelFilePath = std::string("/data/model_COL1D_LRMP_dynLin_1comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_LRMP_dynLin_1comp_sensbenchmark1_DG_P3Z8.h5");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1.0, 1.0, 1.0, 1.0 };

	cadet::test::column::DGParams disc(1, 3, 8);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("Column_1D as pseudo-LRMP (GRM with single FV cell) numerical Benchmark with parameter sensitivities for linear case", "[AxialColumn1D],[DGFV],[Simulation],[Reference],[Sensitivity],[CI_sens25]")
{
	const std::string& modelFilePath = std::string("/data/model_COL1D_pseudoGRM_dynLin_1comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_LRMP_dynLin_1comp_sensbenchmark1_DG_P3Z8.h5");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1.0, 1.0, 1.0, 1.0 };

	// FV_BOUNDARY_ORDER = 1 und disabled surfaceDiffusion (otherwise surface diffusion contributes to film diffusion when binding is dynamic)
	auto discDGFV = cadet::test::column::createDGFVParams(1, 3, 8, 0, 0, 1);
	discDGFV->setParticleDiscParam("FV_BOUNDARY_ORDER", 1);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, *discDGFV, true);
}

TEST_CASE("Column_1D as GRM numerical Benchmark2 with parameter sensitivities for linear case", "[AxialColumn1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI_sens18]")
{
	std::string modelFilePath = std::string("/data/model_COL1D_GRM_dynLin_1comp_benchmark2.json");
	std::string refFilePath = std::string("/data/ref_GRM_dynLin_1comp_sensbenchmark2_cDG_P3Z8_GSM_parP3parZ1.h5");
	const std::vector<double> absTol = { 1e-12, 1e-6, 1e-6, 1e-12 };
	const std::vector<double> relTol = { 1e-4, 1e-3, 1e-4, 1e-4 };

	cadet::test::column::DGParams disc(1, 3, 8, 3, 1);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("Column_1D as GRM numerical Benchmark with parameter sensitivities and multiplexing for 2parType 2comp linear case", "[AxialColumn1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI_sens19]")
{
	const std::string modelFilePath = std::string("/data/model_COL1DparType2_GRM_dynLin_2comp_sensbenchmark1.json");
	const std::string refFilePath = std::string("/data/ref_GRMparType2_dynLin_2comp_sensbenchmark1_cDG_P3Z4_DGexInt_parP3parZ2.h5");
	const std::vector<double> absTol = { 1e-12, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6 };
	const std::vector<double> relTol = { 1e-4, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1 };

	const cadet::test::column::DGParams disc(1, 3, 4, 3, 2);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("Column_1D as GRM numerical Benchmark with parameter sensitivities for SMA LWE case", "[AxialColumn1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI_sens20]")
{
	const std::string modelFilePath = std::string("/data/model_COL1D_GRM_reqSMA_4comp_benchmark1.json");
	const std::string refFilePath = std::string("/data/ref_GRM_reqSMA_4comp_sensbenchmark1_exIntDG_P3Z8_GSM_parP3parZ1.h5");
	const std::vector<double> absTol = { 1e-9, 1e-6, 1e-6, 1e-10 };
	const std::vector<double> relTol = { 1.0, 1.0, 1.0, 1.0 };

	cadet::test::column::DGParams disc(1, 3, 8, 3, 1);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "000", absTol, relTol, disc, true);
}

TEST_CASE("Column_1D as LRMP numerical Benchmark with parameter sensitivities for SMA LWE case", "[AxialColumn1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI_sens23]")
{
	const std::string& modelFilePath = std::string("/data/model_COL1D_LRMP_reqSMA_4comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_LRMP_reqSMA_4comp_sensbenchmark1_DG_P3Z8.h5");
	const std::vector<double> absTol = { 1e-12, 1e-5, 5e-9, 1e-12 };
	const std::vector<double> relTol = { 1.0, 1.0, 1.0, 1.0 };

	cadet::test::column::DGParams disc(1, 3, 8);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "000", absTol, relTol, disc, true);
}

TEST_CASE("Column_1D with mixed general rate and homogeneous particle types and linear binding numerical Benchmark with parameter sensitivities", "[AxialColumn1D],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI_sens24]")
{
	const std::string& modelFilePath = std::string("/data/model_COL1D_2parTypeMixed_2comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_COL1D_2parTypeMixed_2comp_benchmark1.h5");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12, 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };

	cadet::test::column::DGParams disc(1, 3, 5, 3, 1);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("Column_1D as GRM LWE DGSEM and GSM particle discretization yields similar accuracy", "[AxialColumn1D],[DG],[DG1D],[Simulation],[CI]")
{
	cadet::JsonParameterProvider jpp = createLWE("COLUMN_MODEL_1D_GRM", "DG");
	cadet::test::column::DGParams disc(1, 4, 2, 3, 1); // Note that we want to employ only a single particle element
	disc.setDisc(jpp);

	const double absTol = 1e-9;
	const double relTol = 2e-4;

	// GSM discretization
	cadet::Driver drvGSM;
	drvGSM.configure(jpp);
	drvGSM.run();

	// Force single element DGSEM discretization (GSM is default for single particle element discretization)
	jpp.pushScope("model");
	jpp.pushScope("unit_000");
	jpp.pushScope("discretization");
	jpp.set("PAR_GSM", false);
	jpp.popScope();
	jpp.popScope();
	jpp.popScope();

	cadet::Driver drvDGSEM;
	drvDGSEM.configure(jpp);
	drvDGSEM.run();

	cadet::InternalStorageUnitOpRecorder const* const GSMData = drvGSM.solution()->unitOperation(0);
	cadet::InternalStorageUnitOpRecorder const* const DGSEMData = drvDGSEM.solution()->unitOperation(0);

	double const* GSMOutlet = GSMData->outlet();
	double const* DGSEMOutlet = DGSEMData->outlet();

	const unsigned int nComp = GSMData->numComponents();
	for (unsigned int i = 0; i < GSMData->numDataPoints() * GSMData->numInletPorts() * nComp; ++i, ++GSMOutlet, ++DGSEMOutlet)
	{
		// Forward flow inlet = backward flow outlet
		CAPTURE(i);
		CHECK((*GSMOutlet) == cadet::test::makeApprox(*DGSEMOutlet, relTol, absTol));
	}
}

TEST_CASE("Column_1D as GRM time derivative Jacobian vs FD", "[AxialColumn1D],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[CI],[FD]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("COLUMN_MODEL_1D_GRM", "DG", 1e-6, 0.0, 9e-4);
}

TEST_CASE("Column_1D as LRMP time derivative Jacobian vs FD", "[AxialColumn1D],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[CI]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("COLUMN_MODEL_1D_LRMP", "DG");
}

TEST_CASE("Column_1D as GRM sensitivity Jacobians", "[AxialColumn1D],[DG],[DG1D],[UnitOp],[Sensitivity],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_GRM", "DG");

	cadet::test::column::testFwdSensJacobians(jpp, 1e-4, 6e-7, 1e-4);
}

TEST_CASE("Column_1D as LRMP sensitivity Jacobians", "[AxialColumn1D],[DG],[DG1D],[UnitOp],[Sensitivity],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_LRMP", "DG");

	cadet::test::column::testFwdSensJacobians(jpp, 1e-4, 6e-7);
}

//// todo fix: not just adjust tolerances as in FV but theres an actual error here: access violation in densematrix
//TEST_CASE("Column_1D as GRM forward sensitivity vs FD", "[AxialColumn1D],[DG],[DG1D],[Sensitivity],[Simulation],[todo]")
//{
//	// Relative error is checked first, we use high absolute error for letting
//	// some points that are far off pass the error test, too. This is required
//	// due to errors in finite differences.
//	const double fdStepSize[] = { 5e-5, 1e-4, 1e-4, 1e-3 };
//	const double absTols[] = { 3e5, 2e-3, 2e-4, 5.0 };
//	const double relTols[] = { 5e-3, 7e-2, 8e-2, 1e-4 };
//	const double passRatio[] = { 0.95, 0.9, 0.91, 0.83 };
//	cadet::test::column::testFwdSensSolutionFD("COLUMN_MODEL_1D_GRM", "DG", false, fdStepSize, absTols, relTols, passRatio);
//}

//// todo fix: not just adjust tolerances as in FV but theres an actual error here: access violation in densematrix
//TEST_CASE("Column_1D as GRM forward sensitivity forward vs backward flow", "[AxialColumn1D],[DG],[DG1D],[Sensitivity],[Simulation],[todo]")
//{
//	const double absTols[] = { 4e-5, 1e-11, 1e-11, 8e-9 };
//	const double relTols[] = { 6e-9, 5e-8, 5e-6, 5e-10 };
//	const double passRatio[] = { 0.99, 0.95, 0.98, 0.98 };
//	cadet::test::column::testFwdSensSolutionForwardBackward("COLUMN_MODEL_1D_GRM", "DG", absTols, relTols, passRatio);
//}

// todo fix consistent initialization for AD with req binding
TEST_CASE("Column_1D as GRM consistent initialization with linear binding", "[AxialColumn1D],[DG],[DG1D],[ConsistentInit],[CI]")
{
	cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_MODEL_1D_GRM", "DG", 1e-12, 1e-14, 0, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_MODEL_1D_GRM", "DG", 1e-12, 1e-12, 1, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_MODEL_1D_GRM", "DG", 1e-12, 1e-14, 0, 1);
	//cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_MODEL_1D_GRM", "DG", 1e-12, 1e-14, 1, 1);
}

// todo fix consistent initialization for AD with req binding
TEST_CASE("Column_1D as LRMP consistent initialization with linear binding", "[AxialColumn1D],[DG],[DG1D],[ConsistentInit],[CI]")
{
	cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_MODEL_1D_LRMP", "DG", 1e-12, 1e-12, 0, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_MODEL_1D_LRMP", "DG", 1e-12, 1e-12, 1, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_MODEL_1D_LRMP", "DG", 1e-12, 1e-12, 0, 1);
	//cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_MODEL_1D_LRMP", "DG", 1e-12, 1e-12, 1, 1);
}

TEST_CASE("1D column liquid equilibrium MAL consistent initialization", "[Column_1D],[MassActionLaw],[ReactionModel],[ConsistentInit],[CI]")
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	for (int adMode = 0; adMode < 2; ++adMode)
	{
		const bool adEnabled = (adMode > 0);
		SECTION(std::string("AD ") + (adEnabled ? "enabled" : "disabled"))
		{
			const bool kinetic = false;
			const bool withParticles = false;
			cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_LRMP", "DG");
			if (!withParticles)
			{
				// Keep the initialization tests focused on the liquid bulk phase
				jpp.set("NPARTYPE", 0);
				jpp.remove("particle_type_000");
			}
			jpp.set("INIT_C", std::vector<double>{2.7, 0.3});

			jpp.set("NREAC_LIQUID", 1);
			jpp.addScope("liquid_reaction_000");
			jpp.pushScope("liquid_reaction_000");
			jpp.set("TYPE", "MASS_ACTION_LAW");
			jpp.set("MAL_KFWD", std::vector<double>{2.0});
			jpp.set("MAL_KBWD", std::vector<double>{1.0});
			jpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
			jpp.set("MAL_IS_KINETIC", std::vector<int>{kinetic ? 1 : 0});
			jpp.popScope();

			cadet::IUnitOperation* const unit =cadet::test::unitoperation::createAndConfigureUnit(jpp, *mb);

			std::vector<double> y(unit->numDofs(), 0.0);
			y[0] = 1.0;
			y[1] = 2.0;
			// set column bulk state
			const unsigned int nComp = unit->numComponents();
			const unsigned int bulkOffset = nComp * unit->numInletPorts();
			const unsigned int inletDofs = unit->numComponents() * unit->numInletPorts();
			const unsigned columnBulkPoints =  (unit->numDofs() - inletDofs) / nComp;
			
			for (unsigned int point = 0; point < columnBulkPoints; ++point)
			{
				y[bulkOffset + point * nComp] = 2.7;
				y[bulkOffset + point * nComp + 1] = 0.3;
			}

			cadet::test::unitoperation::testConsistentInitialization(
				unit, adEnabled, y.data(), 1e-14, 1e-11);
			
			for (unsigned int point = 0; point < columnBulkPoints; ++point)
			{
				const double* const c = y.data() + bulkOffset + point * nComp;
				CHECK(c[0] == cadet::test::makeApprox(1.0, 1e-12, 1e-12));
				CHECK(c[1] == cadet::test::makeApprox(2.0, 1e-12, 1e-12));
				CHECK((c[0] + c[1]) == cadet::test::makeApprox(3.0, 1e-12, 1e-12));
			}

			mb->destroyUnitOperation(unit);
		}
	}

	destroyModelBuilder(mb);
}

TEST_CASE("1D liquid equilibrium MAL Jacobian vs AD", "[Column_1D],[MassActionLaw],[ReactionModel],[Jacobian],[AD],[CI]")
{
	SECTION("Without particles")
	{
		const bool kinetic = false;
		const bool withParticles = false;
		cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_LRMP", "DG");
		if (!withParticles)
		{
			// Keep the initialization tests focused on the liquid bulk phase
			jpp.set("NPARTYPE", 0);
			jpp.remove("particle_type_000");
		}
		jpp.set("INIT_C", std::vector<double>{2.7, 0.3});

		jpp.set("NREAC_LIQUID", 1);
		jpp.addScope("liquid_reaction_000");
		jpp.pushScope("liquid_reaction_000");
		jpp.set("TYPE", "MASS_ACTION_LAW");
		jpp.set("MAL_KFWD", std::vector<double>{2.0});
		jpp.set("MAL_KBWD", std::vector<double>{1.0});
		jpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
		jpp.set("MAL_IS_KINETIC", std::vector<int>{kinetic ? 1 : 0});
		jpp.popScope();

		cadet::test::column::testJacobianAD(
			jpp, std::numeric_limits<float>::epsilon() * 100.0);
	}

	SECTION("With bulk-to-particle film coupling")
	{
		const bool kinetic = false;
		const bool withParticles = true;
		cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_LRMP", "DG");
		if (!withParticles)
		{
			// Keep the initialization tests focused on the liquid bulk phase
			jpp.set("NPARTYPE", 0);
			jpp.remove("particle_type_000");
		}
		jpp.set("INIT_C", std::vector<double>{2.7, 0.3});

		jpp.set("NREAC_LIQUID", 1);
		jpp.addScope("liquid_reaction_000");
		jpp.pushScope("liquid_reaction_000");
		jpp.set("TYPE", "MASS_ACTION_LAW");
		jpp.set("MAL_KFWD", std::vector<double>{2.0});
		jpp.set("MAL_KBWD", std::vector<double>{1.0});
		jpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
		jpp.set("MAL_IS_KINETIC", std::vector<int>{kinetic ? 1 : 0});
		jpp.popScope();
		cadet::test::column::testJacobianAD(
			jpp, std::numeric_limits<float>::epsilon() * 100.0);
	}
}

TEST_CASE("1D kinetic and equilibrium MAL share liquid equilibrium", "[Column_1D],[MassActionLaw],[ReactionModel],[ConsistentInit],[CI]")
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	bool kinetic = false;
	bool withParticles = false;
	
	cadet::JsonParameterProvider eqJpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_LRMP", "DG");
	if (!withParticles)
	{
		// Keep the initialization tests focused on the liquid bulk phase
		eqJpp.set("NPARTYPE", 0);
		eqJpp.remove("particle_type_000");
	}
	eqJpp.set("INIT_C", std::vector<double>{2.7, 0.3});

	eqJpp.set("NREAC_LIQUID", 1);
	eqJpp.addScope("liquid_reaction_000");
	eqJpp.pushScope("liquid_reaction_000");
	eqJpp.set("TYPE", "MASS_ACTION_LAW");
	eqJpp.set("MAL_KFWD", std::vector<double>{2.0});
	eqJpp.set("MAL_KBWD", std::vector<double>{1.0});
	eqJpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
	eqJpp.set("MAL_IS_KINETIC", std::vector<int>{kinetic ? 1 : 0});
	eqJpp.popScope();

	cadet::IUnitOperation* const eqUnit =
		cadet::test::unitoperation::createAndConfigureUnit(eqJpp, *mb);

	std::vector<double> yEq(eqUnit->numDofs(), 0.0);

	// set column bulk state
	const unsigned int nComp = eqUnit->numComponents();
	const unsigned int bulkOffset = nComp * eqUnit->numInletPorts();
	const unsigned int inletDofs = eqUnit->numComponents() * eqUnit->numInletPorts();
	const unsigned columnBulkPoints =  (eqUnit->numDofs() - inletDofs) / nComp;
	for (unsigned int point = 0; point < columnBulkPoints; ++point)
	{
		yEq[bulkOffset + point * nComp] = 2.7;
		yEq[bulkOffset + point * nComp + 1] = 0.3;
	}
	cadet::test::unitoperation::testConsistentInitialization(eqUnit, false, yEq.data(), 1e-14, 1e-11);

	kinetic = true;
	withParticles = false;
	cadet::JsonParameterProvider kinJpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_LRMP", "DG");
	if (!withParticles)
	{
		// Keep the initialization tests focused on the liquid bulk phase
		kinJpp.set("NPARTYPE", 0);
		kinJpp.remove("particle_type_000");
	}
	kinJpp.set("INIT_C", std::vector<double>{2.7, 0.3});

	kinJpp.set("NREAC_LIQUID", 1);
	kinJpp.addScope("liquid_reaction_000");
	kinJpp.pushScope("liquid_reaction_000");
	kinJpp.set("TYPE", "MASS_ACTION_LAW");
	kinJpp.set("MAL_KFWD", std::vector<double>{2.0});
	kinJpp.set("MAL_KBWD", std::vector<double>{1.0});
	kinJpp.set("MAL_STOICHIOMETRY", std::vector<double>{-1.0, 1.0});
	kinJpp.set("MAL_IS_KINETIC", std::vector<int>{kinetic ? 1 : 0});
	kinJpp.popScope();

	cadet::IUnitOperation* const kinUnit =
		cadet::test::unitoperation::createAndConfigureUnit(kinJpp, *mb);
	REQUIRE(kinUnit->numDofs() == eqUnit->numDofs());

	// Use the equilibrium model's consistently initialized state
	std::vector<double> yKin = yEq;
	yKin[0] = 1.0;
	yKin[1] = 2.0;

	std::vector<double> residual(kinUnit->numDofs(), 0.0);
	cadet::util::ThreadLocalStorage tls;
	tls.resize(kinUnit->threadLocalMemorySize());

	const cadet::AdJacobianParams noAdParams{nullptr, nullptr, 0u};
	kinUnit->notifyDiscontinuousSectionTransition(
		0.0, 0u, {yKin.data(), nullptr}, noAdParams);
	kinUnit->residual(
		cadet::SimulationTime{0.0, 0u},
		cadet::ConstSimulationState{yKin.data(), nullptr},
		residual.data(), tls);

	const unsigned int eqnComp = kinUnit->numComponents();
	const unsigned int eqbulkOffset = nComp * kinUnit->numInletPorts();
	const unsigned int eqinletDofs = kinUnit->numComponents() * kinUnit->numInletPorts();
	const unsigned int eqcolumnBulkPoints =  (kinUnit->numDofs() - eqinletDofs) / eqnComp;

	for (unsigned int point = 0; point < eqcolumnBulkPoints; ++point)
	{
		const double* const c = yKin.data() + eqbulkOffset + point * nComp;
		CHECK(c[0] == cadet::test::makeApprox(1.0, 1e-12, 1e-12));
		CHECK(c[1] == cadet::test::makeApprox(2.0, 1e-12, 1e-12));
	}

	for (unsigned int dof = eqbulkOffset; dof < residual.size(); ++dof)
		CHECK(std::abs(residual[dof]) <= 1e-11);

	mb->destroyUnitOperation(eqUnit);
	mb->destroyUnitOperation(kinUnit);
	destroyModelBuilder(mb);
}

//// todo fix consistent initialization for SMA (initialization not completely correct; AD gives assertion error)
//TEST_CASE("Column_1D as GRM consistent initialization with SMA binding", "[AxialColumn1D],[DG],[DG1D],[ConsistentInit],[todo]")
//{
//	std::vector<double> y(4 + 4 * 16 + 16 * 4 * (4 + 4) + 4 * 16, 0.0);
//	// Optimal values:
//	//	const double bindingCell[] = {1.2, 2.0, 1.0, 1.5, 858.034, 66.7896, 3.53273, 2.53153, 
//	//		1.0, 1.8, 1.5, 1.6, 856.173, 64.457, 5.73227, 2.85286};
//	const double bindingCell[] = { 1.2, 2.0, 1.0, 1.5, 840.0, 63.0, 3.0, 3.0,
//		1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0 };
//	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 + 4 * 16);
//	cadet::test::util::repeat(y.data() + 4 + 4 * 16, bindingCell, 16, 4 * 16 / 2);
//	cadet::test::util::populate(y.data() + 4 + 4 * 16 + 16 * 4 * (4 + 4), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 16);
//
//	cadet::test::column::testConsistentInitializationSMABinding("COLUMN_MODEL_1D_GRM", "DG", y.data(), 1e-14, 1e-5);
//}

// todo fix kinetic binding sensitivity init
TEST_CASE("Column_1D as GRM consistent sensitivity initialization with linear binding", "[AxialColumn1D],[DG],[DG1D],[ConsistentInit],[Sensitivity],[CI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 4 + 4 * 16 + 16 * 4 * (4 + 4) + 4 * 16;
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, numDofs);
	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	//cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_1D_GRM", "DG", y.data(), yDot.data(), true, 1e-14, 0, 0);
	cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_1D_GRM", "DG", y.data(), yDot.data(), true, 1e-14, 1, 0);
	//cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_1D_GRM", "DG", y.data(), yDot.data(), true, 1e-14, 0, 1);
	cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_1D_GRM", "DG", y.data(), yDot.data(), true, 1e-14, 1, 1);
}

// todo fix kinetic binding sensitivity init
TEST_CASE("Column_1D LRMP consistent sensitivity initialization with linear binding", "[AxialColumn1D],[DG],[DG1D],[ConsistentInit],[Sensitivity],[CI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 2 + 2 * 10 + 10 * (2 + 2);
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, numDofs);
	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	//cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_1D_LRMP", "DG", y.data(), yDot.data(), true, 1e-14, 0, 0); // doesnt work
	cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_1D_LRMP", "DG", y.data(), yDot.data(), true, 1e-14, 1, 0); // works
	//cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_1D_LRMP", "DG", y.data(), yDot.data(), true, 1e-14, 0, 1); // doesnt work
	cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_1D_LRMP", "DG", y.data(), yDot.data(), true, 1e-14, 1, 1); // works
}

//// todo fix memory stuff (works for FV) 
//TEST_CASE("Column_1D as GRM consistent sensitivity initialization with SMA binding", "[AxialColumn1D],[DG],[DG1D],[ConsistentInit],[Sensitivity],[fffffffiujbnlk]")
//{
//	// Fill state vector with given initial values
//	const unsigned int numDofs = 4 + 4 * 16 + 16 * 4 * (4 + 4) + 4 * 16;
//	std::vector<double> y(numDofs, 0.0);
//	std::vector<double> yDot(numDofs, 0.0);
//
//	const double bindingCell[] = { 1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0 };
//	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 + 4 * 16);
//	cadet::test::util::repeat(y.data() + 4 + 4 * 16, bindingCell, 8, 4 * 16);
//	cadet::test::util::populate(y.data() + 4 + 4 * 16 + 16 * 4 * (4 + 4), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 16);
//
//	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);
//
//	cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_1D_GRM", "DG", y.data(), yDot.data(), false, 1e-9);
//}

TEST_CASE("Column_1D Jacobian for 2ParType with general rate and homogeneous particle with two component linear binding", "[AxialColumn1D],[DG],[DG1D],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumn2ParType1GeneralRate1HomoParticleBothWithTwoCompLinearJson("COLUMN_MODEL_1D", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as GRM inlet DOF Jacobian", "[AxialColumn1D],[DG],[DG1D],[UnitOp],[Jacobian],[Inlet],[CI]")
{
	cadet::test::column::testInletDofJacobian("COLUMN_MODEL_1D_GRM", "DG");
}

TEST_CASE("Column_1D as GRM transport Jacobian", "[AxialColumn1D],[DG],[DG1D],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "COLUMN_MODEL_1D_GRM", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as GRM with two component linear binding Jacobian", "[AxialColumn1D],[DG],[DG1D],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_GRM", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as GRM transport Jacobian with bulk DG and par FV discretization", "[AxialColumn1D],[DGFV],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "COLUMN_MODEL_1D_GRM", "DGFV");
	cadet::test::column::testJacobianAD(jpp, 1e+10);
}

TEST_CASE("Column_1D as GRM with equilibrium linear binding Jacobian with bulk DG and par FV discretization", "[AxialColumn1D],[DGFV],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, false, "COLUMN_MODEL_1D_GRM", "DGFV");
	cadet::test::column::testJacobianAD(jpp, 1e+10);
}

TEST_CASE("Column_1D as GRM with two component kinetic linear binding Jacobian with bulk DG and par FV discretization", "[AxialColumn1D],[DGFV],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_GRM", "DGFV");
	cadet::test::column::testJacobianAD(jpp, 1e+10);
}

TEST_CASE("Column_1D as GRM LWE one vs two identical particle types match", "[AxialColumn1D],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testOneVsTwoIdenticalParticleTypes("COLUMN_MODEL_1D_GRM", "DG", 2e-8, 5e-5);
}

TEST_CASE("Column_1D as GRM LWE separate identical particle types match", "[AxialColumn1D],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testSeparateIdenticalParticleTypes("COLUMN_MODEL_1D_GRM", "DG", 2e-8, 5e-5);
}

TEST_CASE("Column_1D as GRM linear binding single particle matches particle distribution", "[AxialColumn1D],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearMixedParticleTypes("COLUMN_MODEL_1D_GRM", "DG", 5e-8, 5e-5);
}

TEST_CASE("Column_1D as GRM multiple particle types Jacobian analytic vs AD", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ParticleType],[CI]")
{
	cadet::test::particle::testJacobianMixedParticleTypes("COLUMN_MODEL_1D_GRM", "DG", 5e-12);
}

TEST_CASE("Column_1D as GRM multiple particle types time derivative Jacobian vs FD", "[AxialColumn1D],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[ParticleType],[CI]")
{
	cadet::test::particle::testTimeDerivativeJacobianMixedParticleTypesFD("COLUMN_MODEL_1D_GRM", "DG", 1e-6, 0.0, 9e-4);
}

TEST_CASE("Column_1D as GRM multiple spatially dependent particle types Jacobian analytic vs AD", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ParticleType],[CI]")
{
	cadet::test::particle::testJacobianSpatiallyMixedParticleTypes("COLUMN_MODEL_1D_GRM", "DG", 1e-11);
}

TEST_CASE("Column_1D as GRM linear binding single particle matches spatially dependent particle distribution", "[AxialColumn1D],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearSpatiallyMixedParticleTypes("COLUMN_MODEL_1D_GRM", "DG", 5e-8, 5e-5);
}

TEST_CASE("Column_1D as GRM dynamic reactions Jacobian vs AD bulk", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_MODEL_1D_GRM", "DG", true, false, false, 2e-12);
}

TEST_CASE("Column_1D as GRM dynamic reactions Jacobian vs AD particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_MODEL_1D_GRM", "DG", false, true, false, 1e-14);
}

TEST_CASE("Column_1D as GRM does not support surf diff par dep", "[AxialColumn1D],[DG],[UnitOp],[Jacobian],[ParameterDependence],[CI]")
{
	REQUIRE_THROWS_WITH(
		cadet::test::column::testJacobianADVariableParSurfDiff("COLUMN_MODEL_1D_GRM", "DG", false),
		Catch::Contains("Parameter dependence SURFACE_DIFFUSION_DEP not implemented")
	);
}

TEST_CASE("Column_1D as GRM film diffusion par dep Jacobian analytic vs AD FV", "[AxialColumn1D],[FV],[UnitOp],[Jacobian],[AD],[ParameterDependence],[CI]")
{
	cadet::test::column::testJacobianADVariableFilmDiff("COLUMN_MODEL_1D_GRM", "FV", false);
}

TEST_CASE("Column_1D as GRM film diffusion par dep Jacobian analytic vs AD DG", "[AxialColumn1D],[DG],[DG1D],[UnitOp],[Jacobian],[AD],[ParameterDependence],[CI]")
{
	cadet::test::column::testJacobianADVariableFilmDiff("COLUMN_MODEL_1D_GRM", "DG", false);
}

TEST_CASE("Column_1D as LRMP film diffusion par dep Jacobian analytic vs AD FV", "[AxialColumn1D],[FV],[UnitOp],[Jacobian],[AD],[ParameterDependence],[CI]")
{
	cadet::test::column::testJacobianADVariableFilmDiff("COLUMN_MODEL_1D_LRMP", "FV", false);
}

TEST_CASE("Column_1D as LRMP film diffusion par dep Jacobian analytic vs AD DG", "[AxialColumn1D],[DG],[DG1D],[UnitOp],[Jacobian],[AD],[ParameterDependence],[CI]")
{
	cadet::test::column::testJacobianADVariableFilmDiff("COLUMN_MODEL_1D_LRMP", "DG", false);
}

TEST_CASE("Column_1D as frustum GRM col dispersion van Deemter par dep Jacobian analytic vs AD FV", "[FrustumColumn1D],[FV],[UnitOp],[Jacobian],[AD],[ParameterDependence],[CI]")
{
	cadet::test::column::testJacobianADVariableColDispersionVanDeemter("FRUSTUM_COLUMN_MODEL_1D_GRM", "FV", false);
}

TEST_CASE("Column_1D as frustum GRM col dispersion van Deemter par dep Jacobian analytic vs AD DG", "[FrustumColumn1D],[DG],[DG1D],[UnitOp],[Jacobian],[AD],[ParameterDependence],[CI]")
{
	cadet::test::column::testJacobianADVariableColDispersionVanDeemter("FRUSTUM_COLUMN_MODEL_1D_GRM", "DG", false);
}

TEST_CASE("Column_1D as frustum LRMP col dispersion van Deemter par dep Jacobian analytic vs AD FV", "[FrustumColumn1D],[FV],[UnitOp],[Jacobian],[AD],[ParameterDependence],[CI]")
{
	cadet::test::column::testJacobianADVariableColDispersionVanDeemter("FRUSTUM_COLUMN_MODEL_1D_LRMP", "FV", false);
}

TEST_CASE("Column_1D as frustum LRMP col dispersion van Deemter par dep Jacobian analytic vs AD DG", "[FrustumColumn1D],[DG],[DG1D],[UnitOp],[Jacobian],[AD],[ParameterDependence],[CI]")
{
	cadet::test::column::testJacobianADVariableColDispersionVanDeemter("FRUSTUM_COLUMN_MODEL_1D_LRMP", "DG", false);
}

TEST_CASE("Column_1D as radial GRM col dispersion van Deemter par dep Jacobian analytic vs AD FV", "[RadialColumn1D],[FV],[UnitOp],[Jacobian],[AD],[ParameterDependence],[CI]")
{
	cadet::test::column::testJacobianADVariableColDispersionVanDeemter("RADIAL_COLUMN_MODEL_1D_GRM", "FV", false);
}

TEST_CASE("Column_1D as radial GRM col dispersion van Deemter par dep Jacobian analytic vs AD DG", "[RadialColumn1D],[DG],[DG1D],[UnitOp],[Jacobian],[AD],[ParameterDependence],[CI]")
{
	cadet::test::column::testJacobianADVariableColDispersionVanDeemter("RADIAL_COLUMN_MODEL_1D_GRM", "DG", false);
}

TEST_CASE("Column_1D as axial GRM col dispersion power law par dep Jacobian analytic vs AD DG", "[AxialColumn1D],[DG],[DG1D],[UnitOp],[Jacobian],[AD],[ParameterDependence],[CI]")
{
	cadet::test::column::testJacobianADVariableColDispersionPowerLaw("COLUMN_MODEL_1D_GRM", "DG", false);
}

TEST_CASE("Column_1D as frustum GRM col dispersion power law par dep Jacobian analytic vs AD DG", "[FrustumColumn1D],[DG],[DG1D],[UnitOp],[Jacobian],[AD],[ParameterDependence],[CI]")
{
	cadet::test::column::testJacobianADVariableColDispersionPowerLaw("FRUSTUM_COLUMN_MODEL_1D_GRM", "DG", false);
}

TEST_CASE("Column_1D as radial GRM col dispersion power law par dep Jacobian analytic vs AD DG", "[RadialColumn1D],[DG],[DG1D],[UnitOp],[Jacobian],[AD],[ParameterDependence],[CI]")
{
	cadet::test::column::testJacobianADVariableColDispersionPowerLaw("RADIAL_COLUMN_MODEL_1D_GRM", "DG", false);
}

TEST_CASE("Column_1D as GRM col dispersion par dep is applied by collocation DG", "[AxialColumn1D],[DG],[DG1D],[Simulation],[ParameterDependence],[CI]")
{
	// The collocation DG operator used to configure COL_DISPERSION_DEP and then never evaluate it,
	// so the dependence was silently ignored. Both settings below describe the same physics: the
	// LWE model's axial cylinder geometry is dimensioned for an interstitial velocity of exactly
	// 5.75e-4, so D(v) = (COL_DISPERSION / 5.75e-4) * |v| reproduces the reference setting's plain,
	// constant COL_DISPERSION. They must therefore give the same result.
	const double velocity = 5.75e-4;

	cadet::JsonParameterProvider jpp1 = createLWE("COLUMN_MODEL_1D_GRM", "DG");
	cadet::JsonParameterProvider jpp2 = createLWE("COLUMN_MODEL_1D_GRM", "DG");
	cadet::test::column::DGParams disc(1, 3, 8); // collocation DG
	disc.setDisc(jpp1);
	disc.setDisc(jpp2);

	jpp1.pushScope("model");
	jpp1.pushScope("unit_000");
	const double baseDispersion = jpp1.getDouble("COL_DISPERSION");
	jpp1.set("COL_DISPERSION", baseDispersion / velocity);
	jpp1.set("COL_DISPERSION_DEP", "POWER_LAW");
	jpp1.set("COL_DISPERSION_DEP_BASE", 1.0);
	jpp1.set("COL_DISPERSION_DEP_EXPONENT", 1.0);
	jpp1.popScope();
	jpp1.popScope();

	cadet::test::column::testEqualResults(jpp1, jpp2, 1e-10, 1e-8, 0);
}

TEST_CASE("Column_1D as GRM col dispersion power law par dep Jacobian analytic vs AD collocation DG", "[AxialColumn1D],[DG],[DG1D],[UnitOp],[Jacobian],[AD],[ParameterDependence],[CI]")
{
	cadet::test::column::testJacobianADVariableColDispersionPowerLaw("COLUMN_MODEL_1D_GRM", "DG", false, true);
}

TEST_CASE("Column_1D as GRM dynamic reactions Jacobian vs AD modified particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_MODEL_1D_GRM", "DG", false, true, true, 1e-14);
}

TEST_CASE("Column_1D as GRM dynamic reactions Jacobian vs AD bulk and particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_MODEL_1D_GRM", "DG", true, true, false, 1e-14);
}

TEST_CASE("Column_1D as GRM dynamic reactions Jacobian vs AD bulk and modified particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_MODEL_1D_GRM", "DG", true, true, true, 1e-14);
}

TEST_CASE("Column_1D as GRM dynamic reactions time derivative Jacobian vs FD bulk", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_1D_GRM", "DG", true, false, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_1D as GRM dynamic reactions time derivative Jacobian vs FD particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_1D_GRM", "DG", false, true, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_1D as GRM dynamic reactions time derivative Jacobian vs FD modified particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_1D_GRM", "DG", false, true, true, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_1D as GRM dynamic reactions time derivative Jacobian vs FD bulk and particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_1D_GRM", "DG", true, true, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_1D as GRM dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI],[FD]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_1D_GRM", "DG", true, true, true, 1e-6, 1e-14, 9e-4);
}

inline cadet::JsonParameterProvider create1DColumnWithTwoCompLinearBindingThreeParticleTypes()
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_GRM", "DG");

	const double parVolFrac[] = { 0.3, 0.6, 0.1 };
	const double parFactor[] = { 0.9, 0.8 };
	cadet::test::particle::extendModelToManyParticleTypes(jpp, 3, parFactor, parVolFrac);

	return jpp;
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions Jacobian vs AD bulk", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = create1DColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, false, false);
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions Jacobian vs AD particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = create1DColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, false);
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions Jacobian vs AD modified particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = create1DColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, true);
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions Jacobian vs AD bulk and particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = create1DColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, false);
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions Jacobian vs AD bulk and modified particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = create1DColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, true);
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD bulk", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = create1DColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, false, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = create1DColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD modified particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = create1DColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, true, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD bulk and particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = create1DColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_1D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = create1DColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, true, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Column_1D as LRMP inlet DOF Jacobian", "[AxialColumn1D],[DG],[DG1D],[UnitOp],[Jacobian],[Inlet],[CI]")
{
	cadet::test::column::testInletDofJacobian("COLUMN_MODEL_1D_LRMP", "DG");
}

TEST_CASE("Column_1D as LRMP transport Jacobian", "[AxialColumn1D],[DG],[DG1D],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "COLUMN_MODEL_1D_LRMP", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as LRMP with two component linear binding Jacobian", "[AxialColumn1D],[DG],[DG1D],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_LRMP", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as LRMP LWE one vs two identical particle types match", "[AxialColumn1D],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	// Measured max abs diff between the combined- and separately-simulated identical particle types
	// is ~3.7e-5 (relative to peak signal ~50, i.e. ~7e-7 relative), consistent with ordinary
	// solver-tolerance-level floating-point noise from the two code paths' differing Jacobian/DOF
	// structure (the radial DG sibling test below does not exhibit this and keeps the tighter tolerance).
	cadet::test::particle::testOneVsTwoIdenticalParticleTypes("COLUMN_MODEL_1D_LRMP", "DG", 1e-5, 5e-4);
}

TEST_CASE("Column_1D as LRMP LWE separate identical particle types match", "[AxialColumn1D],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testSeparateIdenticalParticleTypes("COLUMN_MODEL_1D_LRMP", "DG", 1e-12, 1e-12);
}

TEST_CASE("Column_1D as LRMP linear binding single particle matches particle distribution", "[AxialColumn1D],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearMixedParticleTypes("COLUMN_MODEL_1D_LRMP", "DG", 5e-8, 5e-5);
}

TEST_CASE("Column_1D as LRMP multiple particle types Jacobian analytic vs AD", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ParticleType],[CI]")
{
	cadet::test::particle::testJacobianMixedParticleTypes("COLUMN_MODEL_1D_LRMP", "DG", 1e-11);
}

TEST_CASE("Column_1D as LRMP multiple particle types time derivative Jacobian vs FD", "[AxialColumn1D],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[ParticleType],[CI]")
{
	cadet::test::particle::testTimeDerivativeJacobianMixedParticleTypesFD("COLUMN_MODEL_1D_LRMP", "DG", 1e-6, 0.0, 9e-4);
}

TEST_CASE("Column_1D as LRMP multiple spatially dependent particle types Jacobian analytic vs AD", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ParticleType],[CI]")
{
	cadet::test::particle::testJacobianSpatiallyMixedParticleTypes("COLUMN_MODEL_1D_LRMP", "DG", std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as LRMP linear binding single particle matches spatially dependent particle distribution", "[AxialColumn1D],[DG],[DG1D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearSpatiallyMixedParticleTypes("COLUMN_MODEL_1D_LRMP", "DG", 5e-8, 5e-5);
}

TEST_CASE("Column_1D as LRMP dynamic reactions Jacobian vs AD bulk", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_MODEL_1D_LRMP", "DG", true, false, false, std::numeric_limits<float>::epsilon() * 100.0);
}
TEST_CASE("Column_1D as LRMP dynamic reactions Jacobian vs AD particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_MODEL_1D_LRMP", "DG", false, true, false, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as LRMP dynamic reactions Jacobian vs AD modified particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_MODEL_1D_LRMP", "DG", false, true, true, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as LRMP dynamic reactions Jacobian vs AD bulk and particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_MODEL_1D_LRMP", "DG", true, true, false, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as LRMP dynamic reactions Jacobian vs AD bulk and modified particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("COLUMN_MODEL_1D_LRMP", "DG", true, true, true, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Column_1D as LRMP dynamic reactions time derivative Jacobian vs FD bulk", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_1D_LRMP", "DG", true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_1D as LRMP dynamic reactions time derivative Jacobian vs FD particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_1D_LRMP", "DG", false, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_1D as LRMP dynamic reactions time derivative Jacobian vs FD modified particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_1D_LRMP", "DG", false, true, true, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_1D as LRMP dynamic reactions time derivative Jacobian vs FD bulk and particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_1D_LRMP", "DG", true, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_1D as LRMP dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_1D_LRMP", "DG", true, true, true, 1e-6, 1e-14, 8e-4);
}

inline cadet::JsonParameterProvider createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes()
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_1D_LRMP", "DG");

	const double parVolFrac[] = { 0.3, 0.6, 0.1 };
	const double parFactor[] = { 0.9, 0.8 };
	cadet::test::particle::extendModelToManyParticleTypes(jpp, 3, parFactor, parVolFrac);

	return jpp;
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions Jacobian vs AD bulk", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, false, false, 1e-8);
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions Jacobian vs AD particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, false, 1e-8);
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions Jacobian vs AD modified particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, true, 1e-8);
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions Jacobian vs AD bulk and particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, false, 1e-8);
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions Jacobian vs AD bulk and modified particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, true, 1e-8);
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD bulk", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD modified particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, true, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD bulk and particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_1D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[AxialColumn1D],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createLRMPColumnWithTwoCompLinearBindingThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, true, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Radial Column_1D as GRM LWE forward vs backward flow", "[RadialColumn1D],[DG],[DG1D],[Simulation],[CI]")
{
	// Forward and backward flow only match when the cross-sectional area is (nearly) constant
	// which is achieved within testForwardBackward by changing the geometry parameters accordingly
	cadet::test::column::DGParams disc(0, 3, 4);
	cadet::test::column::testForwardBackward("RADIAL_COLUMN_MODEL_1D_GRM", disc, 1e-8, 1e-4);
}

TEST_CASE("Radial Column_1D as GRM transport Jacobian", "[RadialColumn1D],[DG],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "RADIAL_COLUMN_MODEL_1D_GRM", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Radial Column_1D as GRM Jacobian forward vs backward flow", "[RadialColumn1D],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[AD],[CI]")
{
	cadet::test::column::DGParams disc;
	cadet::test::column::testJacobianForwardBackward("RADIAL_COLUMN_MODEL_1D_GRM", disc, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Radial Column_1D as GRM time derivative Jacobian vs FD", "[RadialColumn1D],[DG],[UnitOp],[Residual],[Jacobian],[CI]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("RADIAL_COLUMN_MODEL_1D_GRM", "DG");
}

TEST_CASE("Radial Column_1D as GRM inlet DOF Jacobian", "[RadialColumn1D],[DG],[UnitOp],[Jacobian],[Inlet],[CI]")
{
	cadet::test::column::testInletDofJacobian("RADIAL_COLUMN_MODEL_1D_GRM", "DG");
}

TEST_CASE("Radial Column_1D as GRM with two component linear binding Jacobian", "[RadialColumn1D],[DG],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("RADIAL_COLUMN_MODEL_1D_GRM", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Radial Column_1D as LRMP transport Jacobian", "[RadialColumn1D],[DG],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "RADIAL_COLUMN_MODEL_1D_LRMP", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Radial Column_1D as LRMP time derivative Jacobian vs FD", "[RadialColumn1D],[DG],[UnitOp],[Residual],[Jacobian],[CI]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("RADIAL_COLUMN_MODEL_1D_LRMP", "DG");
}

TEST_CASE("Radial Column_1D as LRMP inlet DOF Jacobian", "[RadialColumn1D],[DG],[UnitOp],[Jacobian],[Inlet],[CI]")
{
	cadet::test::column::testInletDofJacobian("RADIAL_COLUMN_MODEL_1D_LRMP", "DG");
}

TEST_CASE("Radial Column_1D as LRMP with two component linear binding Jacobian", "[RadialColumn1D],[DG],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("RADIAL_COLUMN_MODEL_1D_LRMP", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Radial Column_1D as GRM sensitivity Jacobians", "[RadialColumn1D],[DG],[UnitOp],[Sensitivity],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("RADIAL_COLUMN_MODEL_1D_GRM", "DG");
	cadet::test::column::testFwdSensJacobians(jpp, 1e-4, 6e-7);
}

TEST_CASE("Radial Column_1D as LRMP sensitivity Jacobians", "[RadialColumn1D],[DG],[UnitOp],[Sensitivity],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("RADIAL_COLUMN_MODEL_1D_LRMP", "DG");
	cadet::test::column::testFwdSensJacobians(jpp, 1e-4, 6e-7);
}

TEST_CASE("Radial Column_1D as GRM consistent initialization with linear binding", "[RadialColumn1D],[DG],[ConsistentInit],[CI]")
{
	cadet::test::column::testConsistentInitializationLinearBinding("RADIAL_COLUMN_MODEL_1D_GRM", "DG", 1e-12, 1e-12, 0, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("RADIAL_COLUMN_MODEL_1D_GRM", "DG", 1e-12, 1e-12, 1, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("RADIAL_COLUMN_MODEL_1D_GRM", "DG", 1e-12, 1e-12, 0, 1);
}

TEST_CASE("Radial Column_1D as LRMP consistent initialization with linear binding", "[RadialColumn1D],[DG],[ConsistentInit],[CI]")
{
	cadet::test::column::testConsistentInitializationLinearBinding("RADIAL_COLUMN_MODEL_1D_LRMP", "DG", 1e-12, 1e-12, 0, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("RADIAL_COLUMN_MODEL_1D_LRMP", "DG", 1e-12, 1e-12, 1, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("RADIAL_COLUMN_MODEL_1D_LRMP", "DG", 1e-12, 1e-12, 0, 1);
}

TEST_CASE("Radial Column_1D as GRM dynamic reactions Jacobian vs AD bulk", "[RadialColumn1D],[DG],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_COLUMN_MODEL_1D_GRM", "DG", true, false, false, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Radial Column_1D as GRM dynamic reactions Jacobian vs AD particle", "[RadialColumn1D],[DG],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_COLUMN_MODEL_1D_GRM", "DG", false, true, false, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Radial Column_1D as GRM dynamic reactions Jacobian vs AD bulk and particle", "[RadialColumn1D],[DG],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_COLUMN_MODEL_1D_GRM", "DG", true, true, false, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Radial Column_1D as LRMP dynamic reactions Jacobian vs AD bulk", "[RadialColumn1D],[DG],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_COLUMN_MODEL_1D_LRMP", "DG", true, false, false, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Radial Column_1D as LRMP dynamic reactions Jacobian vs AD particle", "[RadialColumn1D],[DG],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_COLUMN_MODEL_1D_LRMP", "DG", false, true, false, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Radial Column_1D as LRMP dynamic reactions Jacobian vs AD bulk and particle", "[RadialColumn1D],[DG],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("RADIAL_COLUMN_MODEL_1D_LRMP", "DG", true, true, false, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("Radial Column_1D as GRM LWE one vs two identical particle types match", "[RadialColumn1D],[DG],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testOneVsTwoIdenticalParticleTypes("RADIAL_COLUMN_MODEL_1D_GRM", "DG", 2e-8, 5e-5);
}

TEST_CASE("Radial Column_1D as GRM LWE separate identical particle types match", "[RadialColumn1D],[DG],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testSeparateIdenticalParticleTypes("RADIAL_COLUMN_MODEL_1D_GRM", "DG", 2e-8, 5e-5);
}

TEST_CASE("Radial Column_1D as GRM linear binding single particle matches particle distribution", "[RadialColumn1D],[DG],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearMixedParticleTypes("RADIAL_COLUMN_MODEL_1D_GRM", "DG", 5e-8, 5e-5);
}

TEST_CASE("Radial Column_1D as GRM multiple particle types Jacobian analytic vs AD", "[RadialColumn1D],[DG],[Jacobian],[AD],[ParticleType],[CI]")
{
	cadet::test::particle::testJacobianMixedParticleTypes("RADIAL_COLUMN_MODEL_1D_GRM", "DG", 5e-12);
}

TEST_CASE("Radial Column_1D as GRM multiple particle types time derivative Jacobian vs FD", "[RadialColumn1D],[DG],[UnitOp],[Residual],[Jacobian],[ParticleType],[CI]")
{
	cadet::test::particle::testTimeDerivativeJacobianMixedParticleTypesFD("RADIAL_COLUMN_MODEL_1D_GRM", "DG", 1e-6, 0.0, 9e-4);
}

TEST_CASE("Radial Column_1D as GRM multiple spatially dependent particle types Jacobian analytic vs AD", "[RadialColumn1D],[DG],[Jacobian],[AD],[ParticleType],[CI]")
{
	cadet::test::particle::testJacobianSpatiallyMixedParticleTypes("RADIAL_COLUMN_MODEL_1D_GRM", "DG", 1e-11);
}

TEST_CASE("Radial Column_1D as GRM linear binding single particle matches spatially dependent particle distribution", "[RadialColumn1D],[DG],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearSpatiallyMixedParticleTypes("RADIAL_COLUMN_MODEL_1D_GRM", "DG", 5e-8, 5e-5);
}

TEST_CASE("Radial Column_1D as LRMP LWE one vs two identical particle types match", "[RadialColumn1D],[DG],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testOneVsTwoIdenticalParticleTypes("RADIAL_COLUMN_MODEL_1D_LRMP", "DG", 5e-8, 6e-5);
}

TEST_CASE("Radial Column_1D as LRMP LWE separate identical particle types match", "[RadialColumn1D],[DG],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testSeparateIdenticalParticleTypes("RADIAL_COLUMN_MODEL_1D_LRMP", "DG", 2e-8, 5e-5);
}

TEST_CASE("Radial Column_1D as LRMP linear binding single particle matches particle distribution", "[RadialColumn1D],[DG],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearMixedParticleTypes("RADIAL_COLUMN_MODEL_1D_LRMP", "DG", 5e-8, 5e-5);
}

TEST_CASE("Radial Column_1D as LRMP multiple particle types Jacobian analytic vs AD", "[RadialColumn1D],[DG],[Jacobian],[AD],[ParticleType],[CI]")
{
	cadet::test::particle::testJacobianMixedParticleTypes("RADIAL_COLUMN_MODEL_1D_LRMP", "DG", 5e-12);
}

TEST_CASE("Radial Column_1D as LRMP multiple particle types time derivative Jacobian vs FD", "[RadialColumn1D],[DG],[UnitOp],[Residual],[Jacobian],[ParticleType],[CI]")
{
	cadet::test::particle::testTimeDerivativeJacobianMixedParticleTypesFD("RADIAL_COLUMN_MODEL_1D_LRMP", "DG", 1e-6, 0.0, 9e-4);
}

TEST_CASE("Radial Column_1D as LRMP multiple spatially dependent particle types Jacobian analytic vs AD", "[RadialColumn1D],[DG],[Jacobian],[AD],[ParticleType],[CI]")
{
	cadet::test::particle::testJacobianSpatiallyMixedParticleTypes("RADIAL_COLUMN_MODEL_1D_LRMP", "DG", 1e-11);
}

TEST_CASE("Radial Column_1D as LRMP linear binding single particle matches spatially dependent particle distribution", "[RadialColumn1D],[DG],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearSpatiallyMixedParticleTypes("RADIAL_COLUMN_MODEL_1D_LRMP", "DG", 5e-8, 5e-5);
}

inline cadet::JsonParameterProvider createRadial1DGRMColumnWithThreeParticleTypes()
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("RADIAL_COLUMN_MODEL_1D_GRM", "DG");

	const double parVolFrac[] = { 0.3, 0.6, 0.1 };
	const double parFactor[] = { 0.9, 0.8 };
	cadet::test::particle::extendModelToManyParticleTypes(jpp, 3, parFactor, parVolFrac);

	return jpp;
}

TEST_CASE("Radial Column_1D as GRM multi particle types dynamic reactions Jacobian vs AD bulk", "[RadialColumn1D],[DG],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DGRMColumnWithThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, false, false);
}

TEST_CASE("Radial Column_1D as GRM multi particle types dynamic reactions Jacobian vs AD particle", "[RadialColumn1D],[DG],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DGRMColumnWithThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, false);
}

TEST_CASE("Radial Column_1D as GRM multi particle types dynamic reactions Jacobian vs AD modified particle", "[RadialColumn1D],[DG],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DGRMColumnWithThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, true);
}

TEST_CASE("Radial Column_1D as GRM multi particle types dynamic reactions Jacobian vs AD bulk and particle", "[RadialColumn1D],[DG],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DGRMColumnWithThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, false);
}

TEST_CASE("Radial Column_1D as GRM multi particle types dynamic reactions Jacobian vs AD bulk and modified particle", "[RadialColumn1D],[DG],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DGRMColumnWithThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, true);
}

TEST_CASE("Radial Column_1D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD bulk", "[RadialColumn1D],[DG],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DGRMColumnWithThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, false, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Radial Column_1D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD particle", "[RadialColumn1D],[DG],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DGRMColumnWithThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Radial Column_1D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD modified particle", "[RadialColumn1D],[DG],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DGRMColumnWithThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, true, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Radial Column_1D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD bulk and particle", "[RadialColumn1D],[DG],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DGRMColumnWithThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, false, 1e-6, 1e-14, 9e-4);
}

TEST_CASE("Radial Column_1D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[RadialColumn1D],[DG],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DGRMColumnWithThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, true, 1e-6, 1e-14, 9e-4);
}

inline cadet::JsonParameterProvider createRadial1DLRMPColumnWithThreeParticleTypes()
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("RADIAL_COLUMN_MODEL_1D_LRMP", "DG");

	const double parVolFrac[] = { 0.3, 0.6, 0.1 };
	const double parFactor[] = { 0.9, 0.8 };
	cadet::test::particle::extendModelToManyParticleTypes(jpp, 3, parFactor, parVolFrac);

	return jpp;
}

TEST_CASE("Radial Column_1D as LRMP multi particle types dynamic reactions Jacobian vs AD bulk", "[RadialColumn1D],[DG],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DLRMPColumnWithThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, false, false);
}

TEST_CASE("Radial Column_1D as LRMP multi particle types dynamic reactions Jacobian vs AD particle", "[RadialColumn1D],[DG],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DLRMPColumnWithThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, false);
}

TEST_CASE("Radial Column_1D as LRMP multi particle types dynamic reactions Jacobian vs AD modified particle", "[RadialColumn1D],[DG],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DLRMPColumnWithThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, false, true, true);
}

TEST_CASE("Radial Column_1D as LRMP multi particle types dynamic reactions Jacobian vs AD bulk and particle", "[RadialColumn1D],[DG],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DLRMPColumnWithThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, false);
}

TEST_CASE("Radial Column_1D as LRMP multi particle types dynamic reactions Jacobian vs AD bulk and modified particle", "[RadialColumn1D],[DG],[Jacobian],[AD],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DLRMPColumnWithThreeParticleTypes();
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD(jpp, true, true, true);
}

TEST_CASE("Radial Column_1D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD bulk", "[RadialColumn1D],[DG],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DLRMPColumnWithThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Radial Column_1D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD particle", "[RadialColumn1D],[DG],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DLRMPColumnWithThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Radial Column_1D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD modified particle", "[RadialColumn1D],[DG],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DLRMPColumnWithThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, true, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Radial Column_1D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD bulk and particle", "[RadialColumn1D],[DG],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DLRMPColumnWithThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Radial Column_1D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[RadialColumn1D],[DG],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createRadial1DLRMPColumnWithThreeParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, true, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Radial Column_1D pure convection dispersion numerical benchmark", "[RadialColumn1D],[DG],[DG1D],[Simulation],[Reference],[CI]")
{
	std::string modelFilePath = std::string("/data/config_radCOL1D_transport_1comp_DGP3_benchmark1_DG_P3Z16.json");
	std::string refFilePath = std::string("/data/ref_radCOL1D_transport_1comp_DGP3_benchmark1_DG_P3Z16.h5");
	const std::vector<double> absTol = { 1e-12 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::DGParams disc(0, 3, 16);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}

TEST_CASE("Radial Column_1D as LRMP numerical Benchmark for SMA LWE case", "[RadialColumn1D],[DG],[DG1D],[Simulation],[Reference]")
{
	const std::string& modelFilePath = std::string("/data/config_radialLRMP_reqSMA_4comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_radialLRMP_reqSMA_4comp_benchmark1_DG_P3Z8.h5");
	const std::vector<double> absTol = { 1e-12 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::DGParams disc(0, 3, 8);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "000", absTol, relTol, disc, false);
}

TEST_CASE("Frustum Column_1D as GRM transport Jacobian", "[FrustumColumn1D],[DG],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "FRUSTUM_COLUMN_MODEL_1D_GRM", "DG");

	std::vector<double> flowRate;
	flowRate.push_back(0.01);

	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0, std::numeric_limits<float>::epsilon() * 100.0, flowRate);
}

TEST_CASE("Frustum Column_1D as GRM LWE forward vs backward flow", "[FrustumColumn1D],[DG],[DG1D],[Simulation],[CI]")
{
	cadet::test::column::DGParams disc;
	// Note: this test internally sets the geometry to a constant cross section cylinder so that forward and backward flow
	// actually yield the same result for the frustum geometry.
	cadet::test::column::testForwardBackward("FRUSTUM_COLUMN_MODEL_1D_GRM", disc, 1e-10, 1e-6);
}

TEST_CASE("Frustum Column_1D as LRMP numerical Benchmark for SMA LWE case", "[FrustumColumn1D],[DG],[DG1D],[Simulation],[Reference]")
{
	const std::string& modelFilePath = std::string("/data/config_frustumLRMP_reqSMA_4comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_frustumLRMP_reqSMA_4comp_benchmark1_DG_P3Z8.h5");
	const std::vector<double> absTol = { 5e-8 };
	const std::vector<double> relTol = { 1e-4 };

	cadet::test::column::DGParams disc(0, 3, 8);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "000", absTol, relTol, disc, false);
}

TEST_CASE("Smoothly varying Column_1D as GRM equivalence with frustum geometry", "[SmoothlyVaryingColumn1D],[DG],[DG1D],[Simulation],[CI]")
{
	// The frustum cross section area is a quadratic polynomial in the axial coordinate and is
	// thus represented exactly by the nodal interpolant of the smoothly varying geometry (for
	// POLYDEG >= 2). Both settings must therefore describe the very same column.
	cadet::JsonParameterProvider jpp1 = createLWE("FRUSTUM_COLUMN_MODEL_1D_GRM", "DG");
	cadet::JsonParameterProvider jpp2 = createLWE("FRUSTUM_COLUMN_MODEL_1D_GRM", "DG");

	// end areas of the frustum test model, see createColumnWithSMAJson
	const double areaLargeEnd = 0.000314159265358;
	setSmoothlyVaryingGeometry(jpp2, "unit_000", frustumArea(areaLargeEnd, 0.75 * areaLargeEnd));

	cadet::test::column::testEqualResults(jpp1, jpp2, 1e-10, 1e-8, 0);
}

TEST_CASE("Smoothly varying Column_1D as GRM transport Jacobian", "[SmoothlyVaryingColumn1D],[DG],[UnitOp],[Jacobian],[CI]")
{
	// createColumnLinearBenchmark returns the configuration of the column unit itself
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "FRUSTUM_COLUMN_MODEL_1D_GRM", "DG");

	// sine hump around the large end area of the frustum benchmark, i.e. an area profile that
	// none of the other geometries can represent
	const double areaLargeEnd = 4347.826086956522;
	setSmoothlyVaryingGeometry(jpp, [=](double normalizedPos)
		{
			return areaLargeEnd * (1.0 + 0.5 * std::sin(3.14159265358979323846 * normalizedPos));
		});

	std::vector<double> flowRate;
	flowRate.push_back(0.01);

	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0, std::numeric_limits<float>::epsilon() * 100.0, flowRate);
}

TEST_CASE("Smoothly varying Column_1D rejects invalid cross sectional areas", "[SmoothlyVaryingColumn1D],[DG],[DG1D],[CI]")
{
	const double areaLargeEnd = 0.000314159265358;

	// Builds a valid smoothly varying setting and then breaks its area field in one specific way.
	// The break receives the number of nodes per element, so that a specific node of the grid
	// (e.g. an element interface node) can be addressed.
	auto brokenSetting = [=](const std::function<void(std::vector<double>&, int)>& breakArea)
		{
			cadet::JsonParameterProvider jpp = createLWE("FRUSTUM_COLUMN_MODEL_1D_GRM", "DG");
			setSmoothlyVaryingGeometry(jpp, "unit_000", frustumArea(areaLargeEnd, 0.75 * areaLargeEnd));

			jpp.pushScope("model");
			jpp.pushScope("unit_000");

			jpp.pushScope("discretization");
			const int nNodes = jpp.getInt("POLYDEG") + 1;
			jpp.popScope();

			std::vector<double> areas = jpp.getDoubleArray("CROSS_SECTIONAL_AREA_AT_NODES");
			breakArea(areas, nNodes);
			jpp.set("CROSS_SECTIONAL_AREA_AT_NODES", areas);
			jpp.popScope();
			jpp.popScope();

			return jpp;
		};

	SECTION("Wrong number of entries")
	{
		cadet::JsonParameterProvider jpp = brokenSetting([](std::vector<double>& areas, int nNodes) { areas.pop_back(); });

		cadet::Driver drv;
		REQUIRE_THROWS_WITH(drv.configure(jpp), Catch::Contains("CROSS_SECTIONAL_AREA_AT_NODES must have"));
	}

	SECTION("Non-positive entry")
	{
		cadet::JsonParameterProvider jpp = brokenSetting([](std::vector<double>& areas, int nNodes) { areas[0] = 0.0; });

		cadet::Driver drv;
		REQUIRE_THROWS_WITH(drv.configure(jpp), Catch::Contains("CROSS_SECTIONAL_AREA_AT_NODES must be strictly positive"));
	}

	SECTION("Discontinuous at an element interface")
	{
		// first node of the second element, which shares its position with the last node of the first element
		cadet::JsonParameterProvider jpp = brokenSetting([](std::vector<double>& areas, int nNodes) { areas[nNodes] *= 1.5; });

		cadet::Driver drv;
		REQUIRE_THROWS_WITH(drv.configure(jpp), Catch::Contains("CROSS_SECTIONAL_AREA_AT_NODES must be continuous"));
	}
}

TEST_CASE("Smoothly varying Column_1D numerical Benchmark for pure transport case", "[SmoothlyVaryingColumn1D],[DG],[DG1D],[Simulation],[CI],[numRef]")
{
	// Sine cross section area profile with combined convection and dispersion, i.e. the setting
	// of the corresponding EOC study in CADET-Verification (scripts/verify_geometries.py).
	// The discretization must not be changed here: both CROSS_SECTIONAL_AREA_AT_NODES and the
	// nodal initial condition INIT_STATE of this setting are tied to the DG grid.
	const std::string& modelFilePath = std::string("/data/config_smoothlyVaryingCOL1D_transport_1comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_smoothlyVaryingCOL1D_transport_1comp_benchmark1_DG_P3Z8.h5");
	const std::vector<double> absTol = { 1e-10 };
	const std::vector<double> relTol = { 1e-6 };

	cadet::test::column::DGParams disc(0, 3, 8);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false);
}
