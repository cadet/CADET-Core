// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2018: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include <catch.hpp>
#include "cadet/cadet.hpp"

#define CADET_LOGGING_DISABLE
#include "Logging.hpp"

#include "ModelBuilderImpl.hpp"
#include "model/ModelSystemImpl.hpp"
#include "JsonTestModels.hpp"
#include "JacobianHelper.hpp"
#include "ColumnTests.hpp"
#include "Utils.hpp"

#include <limits>
#include <vector>

TEST_CASE("ModelSystem Jacobian AD vs analytic", "[ModelSystem],[Jacobian],[AD]")
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	for (int bindMode = 0; bindMode < 2; ++bindMode)
	{
		const bool isKinetic = bindMode;
		SECTION(isKinetic ? "Kinetic binding" : "Quasi-stationary binding")
		{
			// Use some test case parameters
			cadet::JsonParameterProvider jpp = createLinearBenchmark(isKinetic, false, "GENERAL_RATE_MODEL");

			// Extract section times
			jpp.pushScope("solver");
			jpp.pushScope("sections");

			const std::vector<double> secTimes = jpp.getDoubleArray("SECTION_TIMES");
			std::vector<bool> secCont(secTimes.size() - 2, false);
			if (jpp.exists("SECTION_CONTINUITY")) 
				secCont = jpp.getBoolArray("SECTION_CONTINUITY");

			jpp.popScope();
			jpp.popScope();
			
			// Create and configure ModelSystem
			jpp.pushScope("model");
			cadet::test::column::setNumAxialCells(jpp, 10);
			cadet::IModelSystem* const cadSysAna = mb->createSystem(jpp);
			REQUIRE(cadSysAna);
			cadet::model::ModelSystem* const sysAna = reinterpret_cast<cadet::model::ModelSystem*>(cadSysAna);

			cadet::IModelSystem* const cadSysAD = mb->createSystem(jpp);
			REQUIRE(cadSysAD);
			cadet::model::ModelSystem* const sysAD = reinterpret_cast<cadet::model::ModelSystem*>(cadSysAD);

			bool* const secContArray = new bool[secCont.size()];
			std::copy(secCont.begin(), secCont.end(), secContArray);
			sysAna->setSectionTimes(secTimes.data(), secContArray, secTimes.size() - 1);
			sysAD->setSectionTimes(secTimes.data(), secContArray, secTimes.size() - 1);
			delete[] secContArray;

			// Enable AD
			cadet::ad::setDirections(cadet::ad::getMaxDirections());
			for (unsigned int i = 0; i < sysAD->numModels(); ++i)
				sysAD->getUnitOperationModel(i)->useAnalyticJacobian(false);

			cadet::active* adRes = new cadet::active[sysAD->numDofs()];
			cadet::active* adY = new cadet::active[sysAD->numDofs()];

			sysAD->prepareADvectors(adRes, adY, 0);

			// Setup matrices
			sysAna->notifyDiscontinuousSectionTransition(0.0, 0u, nullptr, nullptr, 0u);
			sysAD->notifyDiscontinuousSectionTransition(0.0, 0u, adRes, adY, 0u);

			// Obtain memory for state, Jacobian multiply direction, Jacobian column
			const unsigned int nDof = sysAD->numDofs();
			std::vector<double> y(nDof, 0.0);
			std::vector<double> yDot(nDof, 0.0);
			std::vector<double> jacDir(nDof, 0.0);
			std::vector<double> jacCol1(nDof, 0.0);
			std::vector<double> jacCol2(nDof, 0.0);

			// Fill state vector with some values
			cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, nDof);
			cadet::test::util::populate(yDot.data(), [=](unsigned int idx) { return std::abs(std::sin((idx + nDof) * 0.13)) + 1e-4; }, nDof);

			// Compute state Jacobian
			sysAna->residualWithJacobian(0.0, 0u, 1.0, y.data(), yDot.data(), jacDir.data(), nullptr, nullptr, 0u);
			sysAD->residualWithJacobian(0.0, 0u, 1.0, y.data(), yDot.data(), jacDir.data(), adRes, adY, 0u);
			std::fill(jacDir.begin(), jacDir.end(), 0.0);

			// Compare Jacobians
			cadet::test::checkJacobianPatternFD(
				[sysAna, &yDot](double const* lDir, double* res) -> void { sysAna->residual(0.0, 0u, 1.0, lDir, yDot.data(), res); },
				[sysAD, &y, &yDot](double const* lDir, double* res) -> void { sysAD->multiplyWithJacobian(0.0, 0u, 1.0, y.data(), yDot.data(), lDir, 1.0, 0.0, res); },
				y.data(), jacDir.data(), jacCol1.data(), jacCol2.data(), nDof);

			cadet::test::checkJacobianPatternFD(
				[sysAna, &yDot](double const* lDir, double* res) -> void { sysAna->residual(0.0, 0u, 1.0, lDir, yDot.data(), res); },
				[sysAna, &y, &yDot](double const* lDir, double* res) -> void { sysAna->multiplyWithJacobian(0.0, 0u, 1.0, y.data(), yDot.data(), lDir, 1.0, 0.0, res); },
				y.data(), jacDir.data(), jacCol1.data(), jacCol2.data(), nDof);

			cadet::test::compareJacobian(
				[sysAna, &y, &yDot](double const* lDir, double* res) -> void { sysAna->multiplyWithJacobian(0.0, 0u, 1.0, y.data(), yDot.data(), lDir, 1.0, 0.0, res); },
				[sysAD, &y, &yDot](double const* lDir, double* res) -> void { sysAD->multiplyWithJacobian(0.0, 0u, 1.0, y.data(), yDot.data(), lDir, 1.0, 0.0, res); },
				jacDir.data(), jacCol1.data(), jacCol2.data(), nDof);

			delete[] adRes;
			delete[] adY;
		}
	}
	destroyModelBuilder(mb);
}

TEST_CASE("ModelSystem time derivative Jacobian FD vs analytic", "[ModelSystem],[Jacobian],[AD]")
{
	const double h = 5e-4;
	const double absTol = 0.0;
	const double relTol = std::numeric_limits<float>::epsilon() * 100.0;

	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	for (int bindMode = 0; bindMode < 2; ++bindMode)
	{
		const bool isKinetic = bindMode;
		SECTION(isKinetic ? "Kinetic binding" : "Quasi-stationary binding")
		{
			// Use some test case parameters
			cadet::JsonParameterProvider jpp = createLinearBenchmark(isKinetic, false, "GENERAL_RATE_MODEL");

			// Extract section times
			jpp.pushScope("solver");
			jpp.pushScope("sections");

			const std::vector<double> secTimes = jpp.getDoubleArray("SECTION_TIMES");
			std::vector<bool> secCont(secTimes.size() - 2, false);
			if (jpp.exists("SECTION_CONTINUITY")) 
				secCont = jpp.getBoolArray("SECTION_CONTINUITY");

			jpp.popScope();
			jpp.popScope();
			
			// Create and configure ModelSystem
			jpp.pushScope("model");
			cadet::test::column::setNumAxialCells(jpp, 10);
			cadet::IModelSystem* const cadSys = mb->createSystem(jpp);
			REQUIRE(cadSys);
			cadet::model::ModelSystem* const sys = reinterpret_cast<cadet::model::ModelSystem*>(cadSys);

			bool* const secContArray = new bool[secCont.size()];
			std::copy(secCont.begin(), secCont.end(), secContArray);
			sys->setSectionTimes(secTimes.data(), secContArray, secTimes.size() - 1);
			delete[] secContArray;

			// Setup matrices
			sys->notifyDiscontinuousSectionTransition(0.0, 0u, nullptr, nullptr, 0u);

			// Obtain memory for state, Jacobian multiply direction, Jacobian column
			const unsigned int nDof = sys->numDofs();
			std::vector<double> y(nDof, 0.0);
			std::vector<double> yDot(nDof, 0.0);
			std::vector<double> jacDir(nDof, 0.0);
			std::vector<double> jacCol1(nDof, 0.0);
			std::vector<double> jacCol2(nDof, 0.0);

			// Fill state vector with some values
			cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, nDof);
			cadet::test::util::populate(yDot.data(), [=](unsigned int idx) { return std::abs(std::sin((idx + nDof) * 0.13)) + 1e-4; }, nDof);

			// Compute state Jacobian
			sys->residualWithJacobian(0.0, 0u, 1.0, y.data(), yDot.data(), jacDir.data(), nullptr, nullptr, 0u);
			std::fill(jacDir.begin(), jacDir.end(), 0.0);

			// Compare Jacobians
			cadet::test::compareJacobianFD(
				[sys, &y](double const* lDir, double* res) -> void { sys->residual(0.0, 0u, 1.0, y.data(), lDir, res); },
				[sys, &y, &yDot](double const* lDir, double* res) -> void { sys->multiplyWithDerivativeJacobian(0.0, 0u, 1.0, y.data(), yDot.data(), lDir, res); },
				yDot.data(), jacDir.data(), jacCol1.data(), jacCol2.data(), nDof, h, absTol, relTol);
		}
	}
	destroyModelBuilder(mb);
}

TEST_CASE("ModelSystem sensitivity Jacobians", "[ModelSystem],[Sensitivity]")
{
	const double h = 5e-5;
	const double absTol = 5e-8;
	const double relTol = 5e-6; // std::numeric_limits<float>::epsilon() * 100.0;

	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	for (int bindMode = 0; bindMode < 2; ++bindMode)
	{
		const bool isKinetic = bindMode;
		SECTION(isKinetic ? "Kinetic binding" : "Quasi-stationary binding")
		{
			// Use some test case parameters
			cadet::JsonParameterProvider jpp = createLinearBenchmark(isKinetic, false, "GENERAL_RATE_MODEL");

			// Extract section times
			jpp.pushScope("solver");
			jpp.pushScope("sections");

			const std::vector<double> secTimes = jpp.getDoubleArray("SECTION_TIMES");
			std::vector<bool> secCont(secTimes.size() - 2, false);
			if (jpp.exists("SECTION_CONTINUITY")) 
				secCont = jpp.getBoolArray("SECTION_CONTINUITY");

			jpp.popScope();
			jpp.popScope();
			
			// Create and configure ModelSystem
			jpp.pushScope("model");
			cadet::test::column::setNumAxialCells(jpp, 10);
			cadet::IModelSystem* const cadSys = mb->createSystem(jpp);
			REQUIRE(cadSys);
			cadet::model::ModelSystem* const sys = reinterpret_cast<cadet::model::ModelSystem*>(cadSys);

			bool* const secContArray = new bool[secCont.size()];
			std::copy(secCont.begin(), secCont.end(), secContArray);
			sys->setSectionTimes(secTimes.data(), secContArray, secTimes.size() - 1);
			delete[] secContArray;

			// Enable AD
			cadet::ad::setDirections(cadet::ad::getMaxDirections());
			cadet::active* adRes = new cadet::active[sys->numDofs()];
			sys->prepareADvectors(adRes, nullptr, 0);

			// Add dispersion parameter sensitivity
			REQUIRE(sys->setSensitiveParameter(cadet::makeParamId(cadet::hashString("COL_DISPERSION"), 0, cadet::CompIndep, cadet::BoundPhaseIndep, cadet::ReactionIndep, cadet::SectionIndep), 0, 1.0));

			// Setup matrices
			sys->notifyDiscontinuousSectionTransition(0.0, 0u, adRes, nullptr, 0u);

			// Obtain memory for state, Jacobian multiply direction, Jacobian column
			const unsigned int nDof = sys->numDofs();
			const std::vector<double> zeros(nDof, 0.0);
			const std::vector<double> ones(nDof, 1.0);
			std::vector<double> y(nDof, 0.0);
			std::vector<double> yDot(nDof, 0.0);
			std::vector<double> jacDir(nDof, 0.0);
			std::vector<double> jacCol1(nDof, 0.0);
			std::vector<double> jacCol2(nDof, 0.0);
			std::vector<double> temp1(nDof, 0.0);
			std::vector<double> temp2(nDof, 0.0);
			std::vector<double> temp3(nDof, 0.0);

			std::vector<const double*> yS(1, zeros.data());
			std::vector<const double*> ySdot(1, zeros.data());
			std::vector<double*> resS(1, nullptr);

			// Fill state vector with some values
			cadet::test::util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, nDof);
			cadet::test::util::populate(yDot.data(), [=](unsigned int idx) { return std::abs(std::sin((idx + nDof) * 0.13)) + 1e-4; }, nDof);

			// Calculate Jacobian
			sys->residualWithJacobian(0.0, 0u, 1.0, y.data(), yDot.data(), jacDir.data(), adRes, nullptr, 0u);

			// Check state Jacobian
			cadet::test::compareJacobianFD(
				[&](double const* lDir, double* res) -> void {
					yS[0] = lDir;
					resS[0] = res;
					sys->residualSensFwd(1, 0.0, 0u, 1.0, y.data(), yDot.data(), nullptr, yS, ySdot, resS, adRes, temp1.data(), temp2.data(), temp3.data());
				}, 
				[&](double const* lDir, double* res) -> void { sys->multiplyWithJacobian(0.0, 0u, 1.0, y.data(), yDot.data(), lDir, 1.0, 0.0, res); }, 
				zeros.data(), jacDir.data(), jacCol1.data(), jacCol2.data(), nDof, h, absTol, relTol);

			// Reset evaluation point
			yS[0] = zeros.data();
			ySdot[0] = zeros.data();

			// Check time derivative Jacobian
			cadet::test::compareJacobianFD(
				[&](double const* lDir, double* res) -> void {
					ySdot[0] = lDir;
					resS[0] = res;
					sys->residualSensFwd(1, 0.0, 0u, 1.0, y.data(), yDot.data(), nullptr, yS, ySdot, resS, adRes, temp1.data(), temp2.data(), temp3.data());
				}, 
				[&](double const* lDir, double* res) -> void { sys->multiplyWithDerivativeJacobian(0.0, 0u, 1.0, y.data(), yDot.data(), lDir, res); }, 
				zeros.data(), jacDir.data(), jacCol1.data(), jacCol2.data(), nDof, h, absTol, relTol);

			delete[] adRes;
		}
	}

	destroyModelBuilder(mb);
}
