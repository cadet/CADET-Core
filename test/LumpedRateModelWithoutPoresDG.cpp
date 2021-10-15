// =============================================================================
//  CADET
//  
//  Copyright © 2008-2021: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include <catch.hpp>
#include "Approx.hpp"

#include "cadet/cadet.hpp"
#include "UnitOperationTests.hpp"
#include "JsonTestModels.hpp"
#include "Utils.hpp"
#include "linalg/Norms.hpp"
#include "SimulationTypes.hpp"
#include "ParallelSupport.hpp"

#include <cmath>
#include <functional>
#include <cstdint>


namespace
{
	void setParameters(cadet::JsonParameterProvider& jpp, unsigned int nCol, unsigned int nNodes)
	{
		int level = 0;

		if (jpp.exists("model"))
		{
			jpp.pushScope("model");
			++level;
		}
		if (jpp.exists("unit_000"))
		{
			jpp.pushScope("unit_000");
			++level;
		}

		// Set model parameters
		// jpp.set("PARAM", value);

		jpp.pushScope("discretization");

		// Set discretization parameters
		jpp.set("NCOL", static_cast<int>(nCol));
		jpp.set("NNODES", static_cast<int>(nNodes));

		jpp.popScope();
	
		for (int l = 0; l < level; ++l)
			jpp.popScope();
	}
}

TEST_CASE("LRM_DG test", "[LRM_DG],[UnitOp]")
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	cadet::JsonParameterProvider jpp = createColumnWithTwoCompLinearBinding(uoType);
	setParameters(jpp, 10, 2);
	cadet::IUnitOperation* const unit = unitoperation::createAndConfigureUnit(jpp, *mb);

	// Disable AD
	unit->useAnalyticJacobian(true);

	// Obtain memory for state, etc.
	std::vector<double> y(unit->numDofs(), 0.0);
	std::vector<double> yDot(unit->numDofs(), 0.0);
	std::vector<double> res(unit->numDofs(), 0.0);

	// Fill state vector with some values
	util::populate(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, unit->numDofs());
	util::populate(yDot.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.11)) + 2e-4; }, unit->numDofs());

	const AdJacobianParams noAdParams{nullptr, nullptr, 0u};
	const ConstSimulationState simState{y.data(), yDot.data()};
	unit->notifyDiscontinuousSectionTransition(0.0, 0u, simState, noAdParams);

	// Evaluate residual
	const SimulationTime simTime{0.0, 0u};
	cadet::util::ThreadLocalStorage tls;
	tls.resize(unit->threadLocalMemorySize());

	unit->residualWithJacobian(simTime, simState, res.data(), noAdParams, tls);

	mb->destroyUnitOperation(unit);
	destroyModelBuilder(mb);
}

