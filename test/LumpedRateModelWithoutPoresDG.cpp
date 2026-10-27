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
#include "Utils.hpp"
#include "JsonTestModels.hpp"
#include "model/parts/DGToolbox.hpp"

#include <cmath>
#include <functional>
#include <vector>

namespace
{
	/**
	 * @brief Sets an axially varying total porosity on a column unit
	 * @details Evaluates the given profile at the LGL nodes of all elements, in the ordering TOTAL_POROSITY
	 *          expects, i.e. element-major with the element interface nodes listed twice. CADET-Core
	 *          multiplies the profile with the geometric cross section area of the unit, so the geometry
	 *          is left untouched.
	 * @param [in,out] jpp parameter provider, scoped to the column unit
	 * @param [in] porosity total porosity as a function of the normalized axial coordinate x / H in [0, 1]
	 */
	void setAxiallyVaryingTotalPorosity(cadet::JsonParameterProvider& jpp, const std::function<double(double)>& porosity)
	{
		jpp.pushScope("discretization");
		const int polyDeg = jpp.getInt("POLYDEG");
		const int nElem = jpp.getInt("NELEM");
		jpp.popScope();

		Eigen::VectorXd nodes(polyDeg + 1);
		Eigen::VectorXd weights(polyDeg + 1);
		cadet::model::parts::dgtoolbox::lglNodesWeights(polyDeg, nodes, weights, false);

		std::vector<double> porosityAtNodes(nElem * (polyDeg + 1));
		for (int elem = 0; elem < nElem; ++elem)
		{
			for (int node = 0; node <= polyDeg; ++node)
				porosityAtNodes[elem * (polyDeg + 1) + node] = porosity((elem + 0.5 * (1.0 + nodes[node])) / nElem);
		}

		jpp.set("TOTAL_POROSITY", porosityAtNodes);
	}
}

TEST_CASE("LRM_DG LWE forward vs backward flow", "[LRM],[DG],[DG1D],[Simulation],[CI]")
{
	cadet::test::column::DGParams disc;

	// Test all integration modes
	for (int i = 0; i <= 2; i++)
	{
		disc.setBulkDiscParam("USE_COLLOCATION_DG", i);
		cadet::test::column::testForwardBackward("LUMPED_RATE_MODEL_WITHOUT_PORES", disc, 6e-9, 6e-2);
	}
}

TEST_CASE("LRM_DG linear pulse vs analytic solution", "[LRM],[DG],[DG1D],[Simulation],[Analytic],[CI]")
{
	cadet::test::column::DGParams disc;
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITHOUT_PORES", "/data/lrm-pulseBenchmark.data", true, true, disc, 2e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITHOUT_PORES", "/data/lrm-pulseBenchmark.data", true, false, disc, 2e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITHOUT_PORES", "/data/lrm-pulseBenchmark.data", false, true, disc, 2e-5, 1e-7);
	cadet::test::column::testAnalyticBenchmark("LUMPED_RATE_MODEL_WITHOUT_PORES", "/data/lrm-pulseBenchmark.data", false, false, disc, 2e-5, 1e-7);
}

TEST_CASE("LRM_DG non-binding linear pulse vs analytic solution", "[LRM],[DG],[DG1D],[Simulation],[Analytic],[NonBinding],[CI]")
{
	cadet::test::column::DGParams disc;
	cadet::test::column::testAnalyticNonBindingBenchmark("LUMPED_RATE_MODEL_WITHOUT_PORES", "/data/lrm-nonBinding.data", true, disc, 2e-5, 1e-7);
	cadet::test::column::testAnalyticNonBindingBenchmark("LUMPED_RATE_MODEL_WITHOUT_PORES", "/data/lrm-nonBinding.data", false, disc, 2e-5, 1e-7);
}

//TEST_CASE("LRM_DG Jacobian forward vs backward flow", "[LRM],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[AD],[fix]")
//{
//	cadet::test::column::DGParams disc;
//
//	// Test all integration modes
//	for (int i = 0; i < 2; i++)
//	{
//		disc.setBulkDiscParam("USE_COLLOCATION_DG", i);
//		cadet::test::column::testJacobianForwardBackward("LUMPED_RATE_MODEL_WITHOUT_PORES", disc, std::numeric_limits<float>::epsilon() * 100.0);
//	}
//}

TEST_CASE("LRM_DG numerical Benchmark with parameter sensitivities for linear case", "[LRM],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI_sens7]")
{
	const std::string& modelFilePath = std::string("/data/model_LRM_dynLin_1comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_LRM_dynLin_1comp_sensbenchmark1_DG_P3Z8.h5");
	const std::vector<double> absTol = { 1e-12, 1e-12, 1e-12, 1e-12 };
	const std::vector<double> relTol = { 1.0, 1.0, 1.0, 1.0 };

	cadet::test::column::DGParams disc(1, 3, 8);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, true);
}

TEST_CASE("LRM_DG numerical Benchmark with parameter sensitivities for SMA LWE case", "[LRM],[DG],[DG1D],[Simulation],[Reference],[Sensitivity],[CI_sens8]")
{
	const std::string& modelFilePath = std::string("/data/model_LRM_reqSMA_4comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_LRM_reqSMA_4comp_sensbenchmark1_DG_P3Z8.h5");
	const std::vector<double> absTol = { 1e-8, 1e-6, 1e-6, 1e-12 };
	const std::vector<double> relTol = { 1.0, 1.0, 1.0, 1.0 };

	cadet::test::column::DGParams disc(1, 3, 8);
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "000", absTol, relTol, disc, true);
}

TEST_CASE("LRM_DG time derivative Jacobian vs FD", "[LRM],[DG],[DG1D],[UnitOp],[Residual],[Jacobian],[CI]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("LUMPED_RATE_MODEL_WITHOUT_PORES", "DG");
}

TEST_CASE("LRM_DG sensitivity Jacobians", "[LRM],[DG],[DG1D],[UnitOp],[Sensitivity],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("LUMPED_RATE_MODEL_WITHOUT_PORES", "DG");

	cadet::test::column::testFwdSensJacobians(jpp, 1e-4, 3e-7, 5e-4);
}

//TEST_CASE("LRM_DG forward sensitivity vs FD", "[LRM],[DG],[DG1D],[Sensitivity],[Simulation]") // todo fix tolerances (also for FV)
//{
//	// Relative error is checked first, we use high absolute error for letting
//	// some points that are far off pass the error test, too. This is required
//	// due to errors in finite differences.
//	const double fdStepSize[] = { 5e-3, 5e-3, 5e-3, 1e-3 };
//	const double absTols[] = { 2e8, 8e-3, 2e-2, 3e-1 };
//	const double relTols[] = { 1e-1, 5e-1, 5e-2, 1e-2 };
//	const double passRatio[] = { 0.88, 0.84, 0.73, 0.87 };
//	cadet::test::column::testFwdSensSolutionFD("LUMPED_RATE_MODEL_WITHOUT_PORES", "DG", false, fdStepSize, absTols, relTols, passRatio);
//}

//TEST_CASE("LRM_DG forward sensitivity forward vs backward flow", "[LRM],[DG],[DG1D],[Sensitivity],[Simulation]") // todo fix  (also for FV) tolerances? why is there a pass ratio here, shouldnt this be precise?
//{
//	const double absTols[] = { 500.0, 8e-7, 9e-7, 2e-3 };
//	const double relTols[] = { 7e-3, 5e-5, 5e-5, 9e-4 };
//	const double passRatio[] = { 0.99, 0.97, 0.97, 0.98 };
//	cadet::test::column::testFwdSensSolutionForwardBackward("LUMPED_RATE_MODEL_WITHOUT_PORES", "DG", absTols, relTols, passRatio);
//}

TEST_CASE("LRM_DG consistent initialization with linear binding", "[LRM],[DG],[DG1D],[ConsistentInit],[CI]")
{
	cadet::test::column::testConsistentInitializationLinearBinding("LUMPED_RATE_MODEL_WITHOUT_PORES", "DG", 1e-12, 1e-12);
}

//TEST_CASE("LRM_DG consistent initialization with SMA binding", "[LRM],[DG],[DG1D],[ConsistentInit]") // todo fix (also for FV)
//{
//	std::vector<double> y(4 + 16 * (4 + 4), 0.0);
//	// Optimal values:
//	//	const double bindingCell[] = {1.2, 2.0, 1.0, 1.5, 858.034, 66.7896, 3.53273, 2.53153, 
//	//		1.0, 1.8, 1.5, 1.6, 856.173, 64.457, 5.73227, 2.85286};
//	const double bindingCell[] = { 1.2, 2.0, 1.0, 1.5, 840.0, 63.0, 3.0, 3.0,
//		1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0 };
//	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4);
//	cadet::test::util::repeat(y.data() + 4, bindingCell, 16, 8);
//
//	cadet::test::column::testConsistentInitializationSMABinding("LUMPED_RATE_MODEL_WITHOUT_PORES", "DG", y.data(), 1e-14, 1e-5);
//}

TEST_CASE("LRM_DG consistent sensitivity initialization with linear binding", "[LRM],[DG],[DG1D],[ConsistentInit],[Sensitivity],[CI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 4 + 10 * (4 + 4);
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, numDofs);
	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITHOUT_PORES", "DG", y.data(), yDot.data(), true, 1e-14);
}

//// todo fix: sigsev read access violation when allocating _tempState = new double[numDofs()] in configure model discretization
//TEST_CASE("LRM_DG consistent sensitivity initialization with SMA binding", "[LRM],[DG],[DG1D],[ConsistentInit],[Sensitivity],[todo]")
//{
//	// Fill state vector with given initial values
//	const unsigned int numDofs = 4 + 10 * (4 + 4);
//	std::vector<double> y(numDofs, 0.0);
//	std::vector<double> yDot(numDofs, 0.0);
//
//	const double bindingCell[] = { 1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0 };
//	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4);
//	cadet::test::util::repeat(y.data() + 4, bindingCell, 8, 16);
//
//	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);
//
//	cadet::test::column::testConsistentInitializationSensitivity("LUMPED_RATE_MODEL_WITHOUT_PORES", "DG", y.data(), yDot.data(), false, 1e-10);
//}

TEST_CASE("LRM_DG inlet DOF Jacobian", "[LRM],[DG],[DG1D],[UnitOp],[Jacobian],[Inlet],[CI]")
{
	cadet::test::column::testInletDofJacobian("LUMPED_RATE_MODEL_WITHOUT_PORES", "DG");
}

TEST_CASE("LRM_DG transport Jacobian", "[LRM],[DG],[DG1D],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, true, "LUMPED_RATE_MODEL_WITHOUT_PORES", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("LRM_DG with two component linear binding Jacobian", "[LRM],[DG],[DG1D],[UnitOp],[Jacobian],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("LUMPED_RATE_MODEL_WITHOUT_PORES", "DG");
	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("LRM_DG dynamic reactions Jacobian vs AD bulk", "[LRM],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("LUMPED_RATE_MODEL_WITHOUT_PORES", "DG", true, false, false, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("LRM_DG dynamic reactions Jacobian vs AD modified bulk", "[LRM],[DG],[DG1D],[Jacobian],[AD],[ReactionModel],[CI]")
{
	cadet::test::reaction::testUnitJacobianDynamicReactionsAD("LUMPED_RATE_MODEL_WITHOUT_PORES", "DG", true, false, true, std::numeric_limits<float>::epsilon() * 100.0);
}

TEST_CASE("LRM_DG dynamic reactions time derivative Jacobian vs FD bulk", "[LRM],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("LUMPED_RATE_MODEL_WITHOUT_PORES", "DG", true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("LRM_DG dynamic reactions time derivative Jacobian vs FD modified bulk", "[LRM],[DG],[DG1D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("LUMPED_RATE_MODEL_WITHOUT_PORES", "DG", true, false, true, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("LRM_DG with a constant nodal TOTAL_POROSITY matches the scalar one", "[LRM],[DG],[DG1D],[Simulation],[CI]")
{
	// A nodal TOTAL_POROSITY whose entries all agree with the scalar one describes the very same
	// column. It takes the other route through the code: the porosity is folded into the area
	// weight of the transport operator instead of dividing the flow rate, the operators are
	// assembled by nodal quadrature rather than analytically, and the solid prefactor is evaluated
	// per discrete point. Were the porosity applied twice, or not at all, the interstitial
	// velocities of the two settings would differ by the factor eps.
	//
	// Both settings are run on the variable coefficient operator, since a nodal porosity routes to
	// it: comparing it against the collocation operator that an axial flow cylinder uses by default
	// would measure the difference between two discretizations rather than the porosity.
	const auto setting = [](const std::function<void(cadet::JsonParameterProvider&)>& modify)
		{
			cadet::JsonParameterProvider jpp = createLinearBenchmark(false, false, "LUMPED_RATE_MODEL_WITHOUT_PORES", "DG");

			jpp.pushScope("model");
			jpp.pushScope("unit_000");
			jpp.pushScope("discretization");
			jpp.set("USE_COLLOCATION_DG", false);
			jpp.popScope();

			modify(jpp);

			jpp.popScope();
			jpp.popScope();

			// The two settings assemble the same operators by different arithmetic, so their
			// round-off differs and the adaptive time integration amplifies it to its own
			// tolerance. Tighten it so that the comparison measures the discretization.
			jpp.pushScope("solver");
			jpp.pushScope("time_integrator");
			jpp.set("ABSTOL", 1e-12);
			jpp.set("RELTOL", 1e-10);
			jpp.popScope();
			jpp.popScope();

			return jpp;
		};

	double totalPorosity = 0.0;
	cadet::JsonParameterProvider jpp1 = setting([&](cadet::JsonParameterProvider& jpp) { totalPorosity = jpp.getDouble("TOTAL_POROSITY"); });
	cadet::JsonParameterProvider jpp2 = setting([&](cadet::JsonParameterProvider& jpp) { setAxiallyVaryingTotalPorosity(jpp, [=](double normalizedPos) { return totalPorosity; }); });

	cadet::test::column::testEqualResults(jpp1, jpp2, 1e-10, 1e-8, 0);
}

TEST_CASE("LRM_DG Jacobian analytic vs AD with axially varying TOTAL_POROSITY", "[LRM],[DG],[DG1D],[UnitOp],[Jacobian],[AD],[CI]")
{
	// The solid phase prefactor (1 - eps_t) / eps_t is the only part of the model that does not
	// follow from the flow through area, and it is the only part of the Jacobian that an axially
	// varying porosity changes. AD differentiates the residual independently of the analytic
	// Jacobian assembly, so this pins the per node evaluation of the prefactor in both. The model
	// therefore has to carry a bound state: without one the prefactor multiplies nothing and the
	// test would only cover the transport operator.
	cadet::JsonParameterProvider jpp = createColumnLinearBenchmark(false, false, "LUMPED_RATE_MODEL_WITHOUT_PORES", "DG");

	const double pi = 3.14159265358979323846;
	setAxiallyVaryingTotalPorosity(jpp, [=](double normalizedPos) { return 0.3 + 0.1 * std::sin(pi * normalizedPos); });

	std::vector<double> flowRate;
	flowRate.push_back(0.01);

	cadet::test::column::testJacobianAD(jpp, std::numeric_limits<float>::epsilon() * 100.0, std::numeric_limits<float>::epsilon() * 100.0, flowRate);
}
