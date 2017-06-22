// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2017: The CADET Authors
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
#include "JsonParameterProvider.hpp"
#include "common/Driver.hpp"
#include "Weno.hpp"

#include <cmath>
#include <functional>

namespace
{
	/**
	 * @brief Sets the WENO order in a configuration
	 * @details Overwrites the WENO_ORDER field in the weno group of the given ParameterProvider.
	 * @param [in,out] jpp ParameterProvider to change the WENO order in
	 * @param [in] order Target order
	 */
	inline void setWenoOrder(cadet::JsonParameterProvider& jpp, int order)
	{
		jpp.pushScope("model");
		jpp.pushScope("unit_000");
		jpp.pushScope("discretization");
		jpp.pushScope("weno");

		jpp.set("WENO_ORDER", order);

		jpp.popScope();
		jpp.popScope();
		jpp.popScope();
		jpp.popScope();
	}
}

/**
 * @brief Reverses the flow of the GRM unit operation
 * @param [in,out] jpp ParameterProvider to change the flow direction in
 */
void reverseFlow(cadet::JsonParameterProvider& jpp)
{
	jpp.pushScope("model");
	jpp.pushScope("unit_000");

	jpp.set("VELOCITY", -jpp.getDouble("VELOCITY"));

	jpp.popScope();
	jpp.popScope();
}

inline Approx makeApprox(double val, double relTol, double absTol)
{
	return Approx(val).epsilon(relTol).margin(absTol);
}

void testWenoForwardBackward(int wenoOrder)
{
	SECTION("Forward vs backward flow (WENO=" + std::to_string(wenoOrder) + ")")
	{
		// Use Load-Wash-Elution test case
		cadet::JsonParameterProvider jpp = createLWE();
		setWenoOrder(jpp, wenoOrder);

		// Forward flow
		cadet::Driver drvFwd;
		drvFwd.configure(jpp);
		drvFwd.run();

		// Backward flow
		reverseFlow(jpp);	
		cadet::Driver drvBwd;
		drvBwd.configure(jpp);
		drvBwd.run();

		cadet::InternalStorageUnitOpRecorder const* const fwdData = drvFwd.solution()->unitOperation(0);
		cadet::InternalStorageUnitOpRecorder const* const bwdData = drvBwd.solution()->unitOperation(0);

		double const* fwdInlet = fwdData->inlet();
		double const* fwdOutlet = fwdData->outlet();
		double const* bwdInlet = bwdData->inlet();
		double const* bwdOutlet = bwdData->outlet();

		const unsigned int nComp = fwdData->numComponents();
		for (unsigned int i = 0; i < fwdData->numDataPoints() * nComp; ++i, ++fwdInlet, ++fwdOutlet, ++bwdInlet, ++bwdOutlet)
		{
			// Compare with relative error 1e-6 and absolute error 1e-9

			// Forward flow inlet = backward flow outlet
			CHECK((*fwdInlet) == makeApprox(*bwdOutlet, 1e-6, 1e-9));
			// Forward flow outlet = backward flow inlet
			CHECK((*fwdOutlet) == makeApprox(*bwdInlet, 1e-6, 1e-9));
		}
	}
}

TEST_CASE("LWE forward vs backward flow", "[GRM],[Simulation]")
{
	// Test all WENO orders
	for (unsigned int i = 1; i < cadet::Weno::maxOrder(); ++i)
		testWenoForwardBackward(i);
}
