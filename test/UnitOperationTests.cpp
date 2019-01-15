// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2019: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include <catch.hpp>

#define CADET_LOGGING_DISABLE
#include "Logging.hpp"

#include "AutoDiff.hpp"
#include "linalg/Norms.hpp"
#include "UnitOperation.hpp"
#include "UnitOperationTests.hpp"

#include "Utils.hpp"

#include <vector>

namespace cadet
{

namespace test
{

namespace unitoperation
{

	void testConsistentInitialization(cadet::IUnitOperation* const unit, bool adEnabled, double* const y, double consTol, double absTol)
	{
		cadet::active* adRes = nullptr;
		cadet::active* adY = nullptr;

		// Enable AD
		if (adEnabled)
		{
			cadet::ad::setDirections(cadet::ad::getMaxDirections());
			unit->useAnalyticJacobian(false);

			adRes = new cadet::active[unit->numDofs()];
			adY = new cadet::active[unit->numDofs()];
			unit->prepareADvectors(adRes, adY, 0);
		}
		else
		{
			unit->useAnalyticJacobian(true);
		}

		// Setup matrices
		unit->notifyDiscontinuousSectionTransition(0.0, 0u, adRes, adY, 0u);

		// Obtain memory for state, Jacobian multiply direction, Jacobian column
		std::vector<double> yDot(unit->numDofs(), 0.0);
		std::vector<double> res(unit->numDofs(), 0.0);

		// Initialize algebraic variables in state vector
		unit->consistentInitialState(0.0, 0u, 1.0, y, adRes, adY, 0u, consTol);

		// Evaluate residual without yDot and update Jacobians
		// Store residual in yDot
		unit->residualWithJacobian(0.0, 0u, 1.0, y, nullptr, yDot.data(), adRes, adY, 0u);
		CAPTURE(yDot);

		// Initialize time derivative of state
		unit->consistentInitialTimeDerivative(0.0, 0u, 1.0, y, yDot.data());

		// Check norm of residual (but skip over connection DOFs at the beginning of the local state vector slice)
		unit->residual(0.0, 0u, 1.0, y, yDot.data(), res.data());
		INFO("Residual " << linalg::linfNorm(res.data() + unit->numComponents(), unit->numDofs() - unit->numComponents()));
		CAPTURE(res);
		CHECK(linalg::linfNorm(res.data() + unit->numComponents(), unit->numDofs() - unit->numComponents()) <= absTol);

		// Clean up
		delete[] adRes;
		delete[] adY;
	}

	void testConsistentInitializationSensitivity(cadet::IUnitOperation* const unit, bool adEnabled, double const* const y, double const* const yDot, double absTol)
	{
		const unsigned int nSens = unit->numSensParams();
		REQUIRE(nSens > 0);

		cadet::active* adRes = new cadet::active[unit->numDofs()];
		cadet::active* adY = nullptr;

		// Enable AD
		if (adEnabled)
		{
			cadet::ad::setDirections(cadet::ad::getMaxDirections());
			unit->useAnalyticJacobian(false);

			adY = new cadet::active[unit->numDofs()];
			unit->prepareADvectors(adRes, adY, nSens);
		}
		else
		{
			unit->useAnalyticJacobian(true);
		}

		// Setup matrices
		unit->notifyDiscontinuousSectionTransition(0.0, 0u, adRes, adY, nSens);

		// Calculate dres / dp and Jacobians
		unit->residualSensFwdWithJacobian(0.0, 0u, 1.0, y, yDot, adRes, adY, nSens);

		// Obtain memory for state, Jacobian multiply direction, Jacobian column
		std::vector<double> sensY(unit->numDofs() * nSens, 0.0);
		std::vector<double> sensYdot(unit->numDofs() * nSens, 0.0);
		std::vector<double> res(unit->numDofs() * nSens, 0.0);

		std::vector<double*> vecSensY(nSens, nullptr);
		std::vector<double*> vecSensYdot(nSens, nullptr);
		std::vector<double const*> cVecSensY(nSens, nullptr);
		std::vector<double const*> cVecSensYdot(nSens, nullptr);
		std::vector<double*> vecRes(nSens, nullptr);

		for (unsigned int i = 0; i < nSens; ++i)
		{
			vecSensY[i] = sensY.data() + i * unit->numDofs();
			vecSensYdot[i] = sensYdot.data() + i * unit->numDofs();
			cVecSensY[i] = sensY.data() + i * unit->numDofs();
			cVecSensYdot[i] = sensYdot.data() + i * unit->numDofs();
			vecRes[i] = res.data() + i * unit->numDofs();
		}

		std::vector<double> tmp1(unit->numDofs(), 0.0);
		std::vector<double> tmp2(unit->numDofs(), 0.0);
		std::vector<double> tmp3(unit->numDofs(), 0.0);

		// Fill sensY and sensYdot with some values
		util::populate(sensY.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, sensY.size());
		util::populate(sensYdot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.9)) + 1e-4; }, sensYdot.size());

		// Initialize sensitivity state vectors
		unit->consistentInitialSensitivity(0.0, 0u, 1.0, y, yDot, vecSensY, vecSensYdot, adRes);

		// Check norm of residual (but skip over connection DOFs at the beginning of the local state vector slice)
		unit->residualSensFwdAdOnly(0.0, 0u, 1.0, y, yDot, adRes);
		unit->residualSensFwdCombine(0.0, 0u, 1.0, y, yDot, cVecSensY, cVecSensYdot, vecRes, adRes, tmp1.data(), tmp2.data(), tmp3.data());

		for (unsigned int i = 0; i < nSens; ++i)
		{
			INFO("Residual " << i << ": " << linalg::linfNorm(vecRes[i] + unit->numComponents(), unit->numDofs() - unit->numComponents()));
			CHECK(linalg::linfNorm(vecRes[i] + unit->numComponents(), unit->numDofs() - unit->numComponents()) <= absTol);
		}

		// Clean up
		delete[] adRes;
		delete[] adY;
	}

} // namespace unitoperation
} // namespace test
} // namespace cadet
