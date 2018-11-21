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

#define CADET_LOGGING_DISABLE
#include "Logging.hpp"

#include "AutoDiff.hpp"
#include "linalg/Norms.hpp"
#include "UnitOperation.hpp"
#include "UnitOperationTests.hpp"

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

} // namespace unitoperation
} // namespace test
} // namespace cadet
