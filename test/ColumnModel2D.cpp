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
#include "Utils.hpp"
#include "JsonTestModels.hpp"

void test2DColumnJacobian(const std::string relModelFilePath, const int maxAxElem, const int maxRadElem, const int axPolyDeg = 0, const int radPolyDeg = 0, const int minAxElem = 1, const int minRadElem = 1)
{
	if (maxAxElem < minAxElem || maxRadElem < minRadElem)
		REQUIRE(false);

	cadet::JsonParameterProvider jpp = cadet::test::column::getReferenceFile(relModelFilePath);

	// get the number of radial ports
	jpp.pushScope("model");
	const int nUnits = jpp.getInt("NUNITS"); // there is one column and (nUnits-1) inlets, one per radial port
	const int columnIdx = 0;
	const std::string unitID = "000";

	// we need to set flowRates for 2D models for the JacobianAD test, since velocity is only set with a call to the setFlowRate function,
	// which in turn is only called if we specify flow rates. We get the flow rates from the connections matrix
	jpp.pushScope("connections");
	jpp.pushScope("switch_000");

	const std::vector<double> connections = jpp.getDoubleArray("CONNECTIONS");

	jpp.popScope();
	jpp.popScope();
	jpp.pushScope("unit_" + unitID);
	jpp.pushScope("discretization");

	const int nRad = (jpp.getInt("RAD_POLYDEG") + 1) * jpp.getInt("RAD_NELEM"); // original number of radial points from file

	// make sure we have enough flow rate entries
	if (radPolyDeg == 0)
		REQUIRE(maxRadElem <= jpp.getInt("RAD_NELEM"));
	else
		REQUIRE((radPolyDeg + 1) * maxRadElem <= nRad);

	std::vector<cadet::active> flowRate;
	for (int i = 1; i <= nRad; i++)
		flowRate.push_back(static_cast<cadet::active>(connections[i * 7 - 1]));

	if (axPolyDeg > 0)
		jpp.set("AX_POLYDEG", axPolyDeg);
	if (radPolyDeg > 0)
		jpp.set("RAD_POLYDEG", radPolyDeg);
	jpp.popScope();

	// This test might run out of memory due to the required AD directions:
	// (axPolyDeg + 1) * axNElem * (radPolyDeg + 1) * radNElem * (nComp + nParType * (nComp + nBound))
	for (int zElem = minAxElem; zElem <= maxAxElem; zElem++)
	{
		for (int rElem = minRadElem; rElem <= maxRadElem; rElem++)
		{
			jpp.pushScope("discretization");
			jpp.set("AX_NELEM", zElem);
			jpp.set("RAD_NELEM", rElem);
			jpp.popScope();

			cadet::test::column::testJacobianAD(jpp, 1e10, std::numeric_limits<float>::epsilon() * 100.0, &flowRate[0]);
		}
	}
}

/* 2D LRMP test cases */

TEST_CASE("Column_2D as LRMP inlet DOF Jacobian", "[Column_2D],[DG],[DG2D],[UnitOp],[Jacobian],[Inlet],[CI]")
{
	cadet::test::column::testInletDofJacobian("COLUMN_MODEL_2D_LRMP", "DG");
}

TEST_CASE("Column_2D as LRMP time derivative Jacobian vs FD", "[Column_2D],[DG],[DG2D],[UnitOp],[Residual],[Jacobian],[CI],[FD]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("COLUMN_MODEL_2D_LRMP", "DG", 1e-6, 0.0, 9e-4);
}

TEST_CASE("Column_2D as LRMP 1comp lin. binding Jacobian", "[Column_2D],[DG],[DG2D],[UnitOp],[Jacobian],[CI]")
{
	const std::string relModelFilePath = std::string("/data/model_COL2D_LRMP_linBnd_1comp.json");

	// Required AD directions:
	// inletDof + (axPolyDeg + 1) * axNElem * (radPolyDeg + 1) * radNElem * (nComp + nParType * (nComp + nBound))
	// req. AD dirs: 152 = 8radPoints + 144 pure dofs (6axPoints*8radPoints) * (1 + 2)
	test2DColumnJacobian(relModelFilePath, 3, 4, 1, 1, 3, 4);
}

TEST_CASE("Column_2D as LRMP with 2parType no binding Jacobian", "[Column_2D],[DG],[DG2D],[UnitOp],[Jacobian],[CI]")
{
	const std::string relModelFilePath = std::string("/data/model_COL2D_LRMP_2parType_1comp.json");

	// Required AD directions:
	// inletDof + (axPolyDeg + 1) * axNElem * (radPolyDeg + 1) * radNElem * (nComp + nParType * (nComp + nBound))
	// req. AD dirs: 152 = 8radPoints + 144 pure dofs (6axPoints*8radPoints) * (1 + 2 * 1)
	test2DColumnJacobian(relModelFilePath, 3, 4, 1, 1, 3, 4);
}

TEST_CASE("Column_2D as LRMP with two component linear binding Jacobian", "[Column_2D],[DG],[DG2D],[UnitOp],[Jacobian],[CI]")
{
	const std::string relModelFilePath = std::string("/data/model_COL2D_LRMP_dynLin_2comp.json");

	// This test might run out of memory due to the required AD directions:
	// inletDof + (axPolyDeg + 1) * axNElem * (radPolyDeg + 1) * radNElem * (nComp + nParType * (nComp + nBound))
	// req. AD dirs: 156 = 2*6radPoints + 144 pure dofs (4axPoints*6radPoints) * (2 + 4)
	test2DColumnJacobian(relModelFilePath, 2, 3, 1, 1);
}

TEST_CASE("Column_2D as LRMP2D analytical reference test for a three zone linear binding case", "[Column_2D],[DG],[DG2D],[Simulation],[Reference],[Analytical],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_COL2D_LRMP3Zone_dynLin_1Comp_benchmark1.json");
	// Note that the analytical reference is actually a 2DGRM but with D^p -> \infty
	const std::string& refFilePath = std::string("/data/refAna_2DLRMP3Zone_dynLin_1Comp_radZ3_benchmark1.h5");
	const std::vector<double> absTol = { 1E-2 }; // relatively high tolerance needed here, since the analyrtical solution computes cross sectional averages
	const std::vector<double> relTol = { 5E-2 };

	//(int exact, int polyDeg, int elem, int parPolyDeg, int parNelem, int radPolyDeg, int radNelem)
	cadet::test::column::DGParams disc(1, 3, 8, 3, 1, 3, 3);
	const int simDataStride = 12; // number of radial ports
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "000", absTol, relTol, disc, true, simDataStride);
}

TEST_CASE("Column_2D as LRMP sensitivity Jacobians", "[Column_2D],[UnitOp],[Sensitivity],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_2D_LRMP", "DG");

	cadet::test::column::testFwdSensJacobians(jpp, 1e-6, 5e-4, 1e-3);
}

TEST_CASE("Column_2D as LRMP consistent initialization with linear binding", "[Column_2D],[ConsistentInit],[CI]")
{
	cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_MODEL_2D_LRMP", "DG", 1e-12, 1e-12, 0, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_MODEL_2D_LRMP", "DG", 1e-12, 1e-12, 1, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_MODEL_2D_LRMP", "DG", 1e-12, 1e-12, 0, 1);
	//cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_MODEL_2D_LRMP", "DG", 1e-12, 1e-12, 1, 1); // @todo AD with req binding does not work
}

//TEST_CASE("Column_2D as LRMP consistent initialization with SMA binding", "[Column_2D],[ConsistentInit],[todo]")  // todo fix (also doesnt work for other models)
//{
//	std::vector<double> y(4 * 3 + 4 * 8 * 3 + 8 * 3 * 3 * (4 + 4) + 4 * 8 * 3, 0.0);
//	// Optimal values:
//	//	const double bindingCell[] = {1.2, 2.0, 1.0, 1.5, 858.034, 66.7896, 3.53273, 2.53153, 
//	//		1.0, 1.8, 1.5, 1.6, 856.173, 64.457, 5.73227, 2.85286};
//	const double bindingCell[] = { 1.2, 2.0, 1.0, 1.5, 840.0, 63.0, 3.0, 3.0,
//		1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0 };
//	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 3 + 4 * 8 * 3);
//	cadet::test::util::repeat(y.data() + 4 * 3 + 4 * 8 * 3, bindingCell, 16, 3 * 8 * 3 / 2);
//	cadet::test::util::populate(y.data() + 4 * 3 + 4 * 8 * 3 + 8 * 3 * 3 * (4 + 4), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 8 * 3);
//
//	cadet::test::column::testConsistentInitializationSMABinding("COLUMN_MODEL_2D_LRMP", "DG", y.data(), 1e-14, 1e-5, 0, 0);
//	cadet::test::column::testConsistentInitializationSMABinding("COLUMN_MODEL_2D_LRMP", "DG", y.data(), 1e-14, 1e-5, 1, 0);
//	//cadet::test::column::testConsistentInitializationSMABinding("COLUMN_MODEL_2D_LRMP", "DG", y.data(), 1e-14, 1e-5, 0, 1);
//	//cadet::test::column::testConsistentInitializationSMABinding("COLUMN_MODEL_2D_LRMP", "DG", y.data(), 1e-14, 1e-5, 1, 1);
//}

TEST_CASE("Column_2D as LRMP consistent sensitivity initialization with linear binding", "[Column_2D],[ConsistentInit],[Sensitivity],[CI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 2 * 3 + 2 * 8 * 3 + 8 * 3 * 3 * (2 + 2) + 2 * 8 * 3;
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, numDofs);
	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	// todo: kinetic binding doesnt work here
	//cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_2D_LRMP", "DG", y.data(), yDot.data(), true, 1e-14, 0, 0);
	cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_2D_LRMP", "DG", y.data(), yDot.data(), true, 1e-14, 1, 0);
	//cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_2D_LRMP", "DG", y.data(), yDot.data(), true, 1e-14, 0, 1);
	cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_2D_LRMP", "DG", y.data(), yDot.data(), true, 1e-14, 1, 1);
}

TEST_CASE("Column_2D as LRMP consistent sensitivity initialization with SMA binding", "[Column_2D],[ConsistentInit],[Sensitivity],[CI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 4 * 3 + 4 * 8 * 3 + 8 * 3 * 3 * (4 + 4) + 4 * 8 * 3;
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);

	const double bindingCell[] = { 1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0 };
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 3 + 4 * 8 * 3);
	cadet::test::util::repeat(y.data() + 4 * 3 + 4 * 8 * 3, bindingCell, 8, 3 * 8 * 3);
	cadet::test::util::populate(y.data() + 4 * 3 + 4 * 8 * 3 + 8 * 3 * 3 * (4 + 4), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 8 * 3);

	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	// todo: kinetic binding doesnt work here
	//cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_2D_LRMP", "DG", y.data(), yDot.data(), false, 1e-9, 0, 0);
	cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_2D_LRMP", "DG", y.data(), yDot.data(), false, 1e-9, 1, 0);
	//cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_2D_LRMP", "DG", y.data(), yDot.data(), false, 1e-9, 0, 1);
	cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_2D_LRMP", "DG", y.data(), yDot.data(), false, 1e-9, 1, 1);
}

TEST_CASE("Column_2D as LRMP LWE one vs two identical particle types match", "[Column_2D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testOneVsTwoIdenticalParticleTypes("COLUMN_MODEL_2D_LRMP", "DG", 1e-7, 5e-5);
}

TEST_CASE("Column_2D as LRMP LWE separate identical particle types match", "[Column_2D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testSeparateIdenticalParticleTypes("COLUMN_MODEL_2D_LRMP", "DG", 1e-7, 5e-5);
}

TEST_CASE("Column_2D as LRMP linear binding single particle matches particle distribution", "[Column_2D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearMixedParticleTypes("COLUMN_MODEL_2D_LRMP", "DG", 5e-8, 5e-5);
}

TEST_CASE("Column_2D as LRMP multiple particle types Jacobian analytic vs AD", "[Column_2D],[Jacobian],[AD],[ParticleType]")
{
	cadet::test::particle::testJacobianMixedParticleTypes("COLUMN_MODEL_2D_LRMP", "DG", 1e10); // @todo figure out why FD Jacobian pattern comparison doesnt work but AD Jacobian comparison does
}

TEST_CASE("Column_2D as LRMP multiple particle types time derivative Jacobian vs FD", "[Column_2D],[UnitOp],[Residual],[Jacobian],[ParticleType],[CI]")
{
	cadet::test::particle::testTimeDerivativeJacobianMixedParticleTypesFD("COLUMN_MODEL_2D_LRMP", "DG", 1e-6, 0.0, 5e-3, true);
}

TEST_CASE("Column_2D as LRMP linear binding single particle matches spatially dependent particle distribution", "[Column_2D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearSpatiallyMixedParticleTypes("COLUMN_MODEL_2D_LRMP", "DG", 5e-8, 5e-5);
}

TEST_CASE("Column_2D as LRMP dynamic reactions time derivative Jacobian vs FD bulk", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_2D_LRMP", "DG", true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_2D as LRMP dynamic reactions time derivative Jacobian vs FD particle", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_2D_LRMP", "DG", false, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_2D as LRMP dynamic reactions time derivative Jacobian vs FD modified particle", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_2D_LRMP", "DG", false, true, true, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_2D as LRMP dynamic reactions time derivative Jacobian vs FD bulk and particle", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_2D_LRMP", "DG", true, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_2D as LRMP dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_2D_LRMP", "DG", true, true, true, 1e-6, 1e-14, 8e-4);
}

inline cadet::JsonParameterProvider createColumnWithTwoCompLinearBindingThreeHomParticleTypes()
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_2D_LRMP", "DG");

	const double parVolFrac[] = { 0.3, 0.6, 0.1 };
	const double parFactor[] = { 0.9, 0.8 };
	cadet::test::particle::extendModelToManyParticleTypes(jpp, 3, parFactor, parVolFrac, true);

	return jpp;
}

TEST_CASE("Column_2D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD bulk", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeHomParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_2D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD particle", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeHomParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_2D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD modified particle", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeHomParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, true, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_2D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD bulk and particle", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeHomParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_2D as LRMP multi particle types dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeHomParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, true, 1e-6, 1e-14, 8e-4);
}

/* 2D DPF test cases */

TEST_CASE("Column_2D as DPF for pure bulk transport Jacobian with radially variable parameters", "[Column_2D],[DG],[DG2D],[UnitOp],[Jacobian],[CI]")
{
	const std::string relModelFilePath = std::string("/data/model_COL2D_DPF_1comp.json");

	// Required AD directions:
	// inletDof + (axPolyDeg + 1) * axNElem * (radPolyDeg + 1) * radNElem * (nComp + nParType * (nComp + nBound))
	// req. AD dirs: 156 = 12radPoints + 144 pure dofs (12axPoints*12radPoints) * (1 + 0)
	test2DColumnJacobian(relModelFilePath, 6, 6, 1, 1, 6, 6);
}

TEST_CASE("Column_2D as DPF numerical Benchmark for pure bulk transport case with three radial zones", "[Column_2D],[DG],[DG2D],[Simulation],[Reference],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_COL2D_DPF3Zone_noFilmDiff_1Comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/ref_2DLRMP3Zone_noFilmDiff_1Comp_benchmark1.h5");
	const std::vector<double> absTol = { 5E-9 };
	const std::vector<double> relTol = { 5E-4 };

	cadet::test::column::DGParams disc;
	const int simDataStride = (3 + 1) * 6; // number of radial ports
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "000", absTol, relTol, disc, true, simDataStride);
}

/* 2D GRM test cases */

TEST_CASE("Column_2D as GRM without radial variance equals Column_1D as GRM numerical Benchmark", "[Column_2D],[DG],[DG2D],[Simulation],[Reference],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_COL2D_GRM_dynLin_1comp_benchmark1.json");
	std::string refFilePath = std::string("/data/ref_GRM_dynLin_1comp_sensbenchmark1_cDG_P3Z8_GSM_parP3parZ1.h5");
	const std::vector<double> absTol = { 5E-7 };
	const std::vector<double> relTol = { 1E-5 };

	cadet::test::column::DGParams disc(0, 3, 8, 3, 1, 3, 3);
	const int simDataStride = (3 + 1) * 3; // number of radial ports
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "001", absTol, relTol, disc, false, simDataStride);
}

TEST_CASE("Column_2D as GRM inlet DOF Jacobian", "[Column_2D],[DG],[DG2D],[UnitOp],[Jacobian],[Inlet],[CI]")
{
	cadet::test::column::testInletDofJacobian("COLUMN_MODEL_2D_GRM", "DG");
}

TEST_CASE("Column_2D as GRM time derivative Jacobian vs FD", "[Column_2D],[DG],[DG2D],[UnitOp],[Residual],[Jacobian],[CI],[FD]")
{
	cadet::test::column::testTimeDerivativeJacobianFD("COLUMN_MODEL_2D_GRM", "DG", 1e-6, 0.0, 9e-4);
}

TEST_CASE("Column_2D as GRM 1comp lin. binding Jacobian", "[Column_2D],[DG],[DG2D],[UnitOp],[Jacobian],[CI]")
{
	const std::string relModelFilePath = std::string("/data/model_COL2D_GRM_linBnd_1comp.json");

	// Required AD directions:
	// inletDof + (axPolyDeg + 1) * axNElem * (radPolyDeg + 1) * radNElem * (nComp + ([parPolyDeg + 1] * parNElem) * (nComp + nBound))
	// req. AD dirs: 116 = 6radPoints + 108 pure dofs (2axPoints*6radPoints) * (1 + 4parPoints * 2)
	test2DColumnJacobian(relModelFilePath, 1, 3, 1, 1, 1, 3);
}

TEST_CASE("Column_2D as GRM with 2parType no binding Jacobian", "[Column_2D],[DG],[DG2D],[UnitOp],[Jacobian],[CI]")
{
	const std::string relModelFilePath = std::string("/data/model_COL2D_GRM_2parType_1comp.json");

	// Required AD directions:
	// inletDof + (axPolyDeg + 1) * axNElem * (radPolyDeg + 1) * radNElem * (nComp + nParTypesPoints * (nComp + nBound))
	// req. AD dirs: 136 = 8radPoints + 128 pure dofs (2axPoints*8radPoints) * (1 + 7 * 1)
	test2DColumnJacobian(relModelFilePath, 1, 4, 1, 1, 1, 4);
}

TEST_CASE("Column_2D as GRM with two component linear binding Jacobian", "[Column_2D],[DG],[DG2D],[UnitOp],[Jacobian],[CI]")
{
	const std::string relModelFilePath = std::string("/data/model_COL2D_GRM_dynLin_2comp.json");

	// This test might run out of memory due to the required AD directions:
	// inletDof + (axPolyDeg + 1) * axNElem * (radPolyDeg + 1) * radNElem * (nComp + nParPoints * (nComp + nBound))
	// req. AD dirs: 142 = 2*6radPoints + 120 pure dofs (2axPoints*6radPoints) * (2 + 2 * 4)
	test2DColumnJacobian(relModelFilePath, 1, 3, 1, 1);
}

TEST_CASE("Column_2D as GRM non limiting particle diffusion analytical reference test for a three zone linear binding case", "[Column_2D],[DG],[DG2D],[Simulation],[Reference],[Analytical],[CI]")
{
	// Note that for this setting we have D^p -> \infty (D^p = 1e-6 suffices here)
	// Hence, this test in combination with the similar test for the LRMP above shows equivalence of the LRMP and GRM with non limiting particle diffusion

	const std::string& modelFilePath = std::string("/data/model_COL2D_GRM3Zone_dynLin_1Comp_benchmark1.json");
	const std::string& refFilePath = std::string("/data/refAna_2DLRMP3Zone_dynLin_1Comp_radZ3_benchmark1.h5");
	const std::vector<double> absTol = { 1E-2 }; // relatively high tolerance needed here, since the analyrtical solution computes cross sectional averages
	const std::vector<double> relTol = { 5E-2 };

	//(int exact, int polyDeg, int elem, int parPolyDeg, int parNelem, int radPolyDeg, int radNelem)
	cadet::test::column::DGParams disc(1, 3, 8, 3, 1, 3, 3);
	const int simDataStride = 12; // number of radial ports
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "000", absTol, relTol, disc, true, simDataStride);
}

TEST_CASE("Column_2D as GRM analytical reference test for a three zone linear binding GRM with surface diffusion", "[GRM2D],[DG],[DG2D],[Simulation],[Reference],[Analytical],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_COL2D_GRMsd3Zone_dynLin_1Comp_benchmark1.json");
	const std::string refFilePath = std::string("/data/refAna_2DGRMsd3Zone_dynLin_1Comp_radZ3_benchmark1.h5");
	const std::vector<double> absTol = { 5E-4 }; // relatively high tolerance needed here, since the analytical solution is given as cross sectional averages
	const std::vector<double> relTol = { 2E-1 };

	cadet::test::column::DGParams disc;
	const int simDataStride = 12; // number of radial ports

	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "000", absTol, relTol, disc, false, simDataStride);
}

TEST_CASE("Column_2D as GRM numerical reference test for a three zone linear binding GRM with surface diffusion", "[GRM2D],[DG],[DG2D],[Simulation],[Reference],[Analytical],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_COL2D_GRMsd3Zone_dynLin_1Comp_benchmark2.json");
	const std::string refFilePath = std::string("/data/ref_COL2D_GRMsd3Zone_dynLin_1Comp_benchmark1_DG_axP3Z8_radP3Z3_parP3Z1.h5");
	const std::vector<double> absTol = { 1E-10 };
	const std::vector<double> relTol = { 1E-8 };

	cadet::test::column::DGParams disc;
	const int simDataStride = 12; // number of radial ports
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "000", absTol, relTol, disc, false, simDataStride);
}

TEST_CASE("Column_2D as GRM without radial variation SMA LWE numerical 1D reference test", "[GRM2D],[DG],[DG2D],[Simulation],[Reference],[Analytical],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_COL2D_GRM2Zone_noRadVar_SMA_LWE.json");
	const std::string refFilePath = std::string("/data/ref_GRM_reqSMA_4comp_sensbenchmark1_exIntDG_P3Z8_GSM_parP3parZ1.h5");
	const std::vector<double> absTol = { 5E-6 };
	const std::vector<double> relTol = { 5E-2 };

	cadet::test::column::DGParams disc(1, 3, 8, 3, 1, 3, 2); // (int exact, int polyDeg, int elem, int parPolyDeg, int parNelem, int radPolyDeg, int radNelem)
	const int simDataStride = 8; // number of radial ports
	const int outletDataStride = 4; // number of components
	const int outletDataOffset = 1; // offset to component to be compared
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "000", absTol, relTol, disc, false, simDataStride, outletDataStride, outletDataOffset);
}

TEST_CASE("Column_2D as GRM sensitivity Jacobians", "[Column_2D],[UnitOp],[Sensitivity],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_2D_GRM", "DG");

	cadet::test::column::testFwdSensJacobians(jpp, 1e-6, 5e-4, 1e-3);
}

TEST_CASE("Column_2D as GRM consistent initialization with linear binding", "[Column_2D],[ConsistentInit],[CI]")
{
	cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_MODEL_2D_GRM", "DG", 1e-12, 1e-12, 0, 0);
	cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_MODEL_2D_GRM", "DG", 1e-12, 1e-12, 1, 0);

	// The following tests are disabled since they require more than the maximum number of AD directions for the CI
	//cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_MODEL_2D_GRM", "DG", 1e-12, 1e-12, 0, 1);
	//cadet::test::column::testConsistentInitializationLinearBinding("COLUMN_MODEL_2D_GRM", "DG", 1e-12, 1e-12, 1, 1);
}

//TEST_CASE("Column_2D as GRM consistent initialization with SMA binding", "[Column_2D],[ConsistentInit]")  // todo fix (also doesnt work for other models)
//{
//	std::vector<double> y(4 * 3 + 4 * 8 * 3 + 8 * 3 * 3 * (4 + 4) + 4 * 8 * 3, 0.0);
//	// Optimal values:
//	//	const double bindingCell[] = {1.2, 2.0, 1.0, 1.5, 858.034, 66.7896, 3.53273, 2.53153, 
//	//		1.0, 1.8, 1.5, 1.6, 856.173, 64.457, 5.73227, 2.85286};
//	const double bindingCell[] = { 1.2, 2.0, 1.0, 1.5, 840.0, 63.0, 3.0, 3.0,
//		1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0 };
//	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 3 + 4 * 8 * 3);
//	cadet::test::util::repeat(y.data() + 4 * 3 + 4 * 8 * 3, bindingCell, 16, 3 * 8 * 3 / 2);
//	cadet::test::util::populate(y.data() + 4 * 3 + 4 * 8 * 3 + 8 * 3 * 3 * (4 + 4), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 8 * 3);
//
//	cadet::test::column::testConsistentInitializationSMABinding("COLUMN_MODEL_2D_GRM", "DG", y.data(), 1e-14, 1e-5, 0, 0);
//	cadet::test::column::testConsistentInitializationSMABinding("COLUMN_MODEL_2D_GRM", "DG", y.data(), 1e-14, 1e-5, 1, 0);
//	//cadet::test::column::testConsistentInitializationSMABinding("COLUMN_MODEL_2D_GRM", "DG", y.data(), 1e-14, 1e-5, 0, 1);
//	//cadet::test::column::testConsistentInitializationSMABinding("COLUMN_MODEL_2D_GRM", "DG", y.data(), 1e-14, 1e-5, 1, 1);
//}

TEST_CASE("Column_2D as GRM consistent sensitivity initialization with linear binding", "[Column_2D],[ConsistentInit],[Sensitivity],[CI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 2 * 3 + 2 * 8 * 3 + 8 * 3 * 3 * (2 + 2) + 2 * 8 * 3;
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, numDofs);
	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	// todo: kinetic binding doesnt work here, same as in other units
	//cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_2D_GRM", "DG", y.data(), yDot.data(), true, 1e-14, 0, 0);
	cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_2D_GRM", "DG", y.data(), yDot.data(), true, 1e-14, 1, 0);
	//cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_2D_GRM", "DG", y.data(), yDot.data(), true, 1e-14, 0, 1);
	cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_2D_GRM", "DG", y.data(), yDot.data(), true, 1e-14, 1, 1);
}

TEST_CASE("Column_2D as GRM consistent sensitivity initialization with SMA binding", "[Column_2D],[ConsistentInit],[Sensitivity],[CI]")
{
	// Fill state vector with given initial values
	const unsigned int numDofs = 4 * 3 + 4 * 8 * 3 + 8 * 3 * 3 * (4 + 4) + 4 * 8 * 3;
	std::vector<double> y(numDofs, 0.0);
	std::vector<double> yDot(numDofs, 0.0);

	const double bindingCell[] = { 1.0, 1.8, 1.5, 1.6, 840.0, 63.0, 6.0, 3.0 };
	cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 3 + 4 * 8 * 3);
	cadet::test::util::repeat(y.data() + 4 * 3 + 4 * 8 * 3, bindingCell, 8, 3 * 8 * 3);
	cadet::test::util::populate(y.data() + 4 * 3 + 4 * 8 * 3 + 8 * 3 * 3 * (4 + 4), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, 4 * 8 * 3);

	cadet::test::util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, numDofs);

	// todo: kinetic binding doesnt work here, same as in other units
	//cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_2D_GRM", "DG", y.data(), yDot.data(), false, 1e-9, 0, 0);
	cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_2D_GRM", "DG", y.data(), yDot.data(), false, 1e-9, 1, 0);
	//cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_2D_GRM", "DG", y.data(), yDot.data(), false, 1e-9, 0, 1);
	cadet::test::column::testConsistentInitializationSensitivity("COLUMN_MODEL_2D_GRM", "DG", y.data(), yDot.data(), false, 1e-9, 1, 1);
}

// todo modify SMA test case so that 2DDG does not produce negative values
TEST_CASE("Column_2D as GRM LWE one vs two identical particle types match", "[Column_2D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testOneVsTwoIdenticalParticleTypes("COLUMN_MODEL_2D_GRM", "DG", 1e-7, 5e-5);
}

// todo modify SMA test case so that 2DDG does not produce negative values
TEST_CASE("Column_2D as GRM LWE separate identical particle types match", "[Column_2D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testSeparateIdenticalParticleTypes("COLUMN_MODEL_2D_GRM", "DG", 1e-7, 5e-5);
}

TEST_CASE("Column_2D as GRM linear binding single particle matches particle distribution", "[Column_2D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearMixedParticleTypes("COLUMN_MODEL_2D_GRM", "DG", 5e-8, 5e-5);
}

TEST_CASE("Column_2D as GRM multiple particle types time derivative Jacobian vs FD", "[Column_2D],[UnitOp],[Residual],[Jacobian],[ParticleType],[CI]")
{
	cadet::test::particle::testTimeDerivativeJacobianMixedParticleTypesFD("COLUMN_MODEL_2D_GRM", "DG", 1e-6, 0.0, 5e-3, true);
}

TEST_CASE("Column_2D as GRM linear binding single particle matches spatially dependent particle distribution", "[Column_2D],[Simulation],[ParticleType],[CI]")
{
	cadet::test::particle::testLinearSpatiallyMixedParticleTypes("COLUMN_MODEL_2D_GRM", "DG", 5e-8, 5e-5);
}

TEST_CASE("Column_2D as GRM dynamic reactions time derivative Jacobian vs FD bulk", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_2D_GRM", "DG", true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_2D as GRM dynamic reactions time derivative Jacobian vs FD particle", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_2D_GRM", "DG", false, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_2D as GRM dynamic reactions time derivative Jacobian vs FD modified particle", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_2D_GRM", "DG", false, true, true, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_2D as GRM dynamic reactions time derivative Jacobian vs FD bulk and particle", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_2D_GRM", "DG", true, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_2D as GRM dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[CI]")
{
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD("COLUMN_MODEL_2D_GRM", "DG", true, true, true, 1e-6, 1e-14, 8e-4);
}

inline cadet::JsonParameterProvider createColumnWithTwoCompLinearBindingThreeGRMParticleTypes()
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding("COLUMN_MODEL_2D_GRM", "DG");

	const double parVolFrac[] = { 0.3, 0.6, 0.1 };
	const double parFactor[] = { 0.9, 0.8 };
	cadet::test::particle::extendModelToManyParticleTypes(jpp, 3, parFactor, parVolFrac, true);

	return jpp;
}

TEST_CASE("Column_2D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD bulk", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeGRMParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, false, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_2D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD particle", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeGRMParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_2D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD modified particle", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeGRMParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, false, true, true, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_2D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD bulk and particle", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeGRMParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, false, 1e-6, 1e-14, 8e-4);
}

TEST_CASE("Column_2D as GRM multi particle types dynamic reactions time derivative Jacobian vs FD bulk and modified particle", "[Column_2D],[Jacobian],[Residual],[ReactionModel],[ParticleType],[CI]")
{
	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBindingThreeGRMParticleTypes();
	cadet::test::reaction::testTimeDerivativeJacobianDynamicReactionsFD(jpp, true, true, true, 1e-6, 1e-14, 8e-4);
}

/* Mixed homogeneous and 1D particles test cases */

TEST_CASE("Column_2D with 3 zones and mixed homogeneous and 1D particles numerical reference test", "[GRM2D],[DG],[DG2D],[Simulation],[Reference],[Analytical],[CI]")
{
	const std::string& modelFilePath = std::string("/data/model_COL2D_GRMsd3Zone2ParType_dynLin_1Comp_benchmark1.json");
	const std::string refFilePath = std::string("/data/ref_COL2D_GRMsd3Zone2ParType_dynLin_1Comp_benchmark1_DG_axP3Z8_radP3Z3_parP3Z1.h5");
	const std::vector<double> absTol = { 1E-12 };
	const std::vector<double> relTol = { 1E-12 };

	cadet::test::column::DGParams disc;
	const int simDataStride = 12; // number of radial ports
	cadet::test::column::testReferenceBenchmark(modelFilePath, refFilePath, "000", absTol, relTol, disc, false, simDataStride);
}