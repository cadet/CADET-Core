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

/**
 * @file 
 * Tests for Driver timeout functionality
 */

#include <catch.hpp>
#include <chrono>
#include <json.hpp>
using json = nlohmann::json;

#include "Logging.hpp"
#include "common/Driver.hpp"
#include "common/JsonParameterProvider.hpp"
#include "JsonTestModels.hpp"

namespace cadet
{
namespace test
{

namespace
{
	/**
	 * @brief Creates a LWE configuration with DG discretization for timeout testing
	 * @details Uses increased discretization to make simulation take longer, allowing
	 *          timeout to be tested reliably
	 * @param [in] uoType Unit operation type
	 * @param [in] timeout timeout
	 * @return Configuration with timeout and DG discretization
	 */
	cadet::JsonParameterProvider createLWEConfigForTimeoutTest(const std::string& uoType, const double timeout, const int polyDeg, const int nElem)
	{
		// Start with standard LWE configuration
		cadet::JsonParameterProvider paramProvider = createLWE(uoType, "DG");

		// Set timeout
		paramProvider.pushScope("solver");
		paramProvider.set("TIMEOUT", timeout);
		paramProvider.popScope();

		// Increase discretization and thus required compute time
		paramProvider.pushScope("model");
		paramProvider.pushScope("unit_000");
		paramProvider.pushScope("discretization");
		paramProvider.set("POLYDEG", polyDeg);
		paramProvider.set("NELEM", nElem);
		paramProvider.popScope();
		paramProvider.popScope();
		paramProvider.popScope();

		return paramProvider;
	}

} // namespace

TEST_CASE("Test Callback: timeout interrupts simulation but data is saved", "[Callback],[Timeout],[CI_callback]")
{
		const double timeout = 2.5;
		const int polyDeg = 3;
		const int nElem = 16;

		Driver driver;
		cadet::JsonParameterProvider pp = createLWEConfigForTimeoutTest("COLUMN_MODEL_1D_GRM", timeout, polyDeg, nElem);

		driver.configure(pp);
		driver.run();

		REQUIRE(driver.simulator()->stoppedByNotificationCallback());
		// factor applied since callback is invoked at the end of a time step, which might be a little later than the specified timeout
		REQUIRE(driver.simulator()->lastSimulationDuration() <= 1.2 * timeout);

		// test that part of the solution is written
		cadet::InternalStorageUnitOpRecorder const* const simData = driver.solution()->unitOperation(0);
		REQUIRE(simData->numDataPoints() > 3);
		REQUIRE(simData->numDataPoints() < 1500);
}

TEST_CASE("Test Callback: timeout is ignored when set to zero", "[Callback],[Timeout],[CI_callback]")
{
	const double timeout = 0.0;
	const int polyDeg = 2;
	const int nElem = 4;

	Driver driver;
	cadet::JsonParameterProvider pp = createLWEConfigForTimeoutTest("COLUMN_MODEL_1D_GRM", timeout, polyDeg, nElem);

	driver.configure(pp);
	driver.run();

	// test that whole solution is written
	cadet::InternalStorageUnitOpRecorder const* const simData = driver.solution()->unitOperation(0);

	REQUIRE(simData->numDataPoints() == 1501);
}

TEST_CASE("Test Callback: timeout is ignored when set to < 0", "[Callback],[Timeout],[CI_callback]")
{
	const double timeout = -5.0;
	const int polyDeg = 2;
	const int nElem = 4;

	Driver driver;
	cadet::JsonParameterProvider pp = createLWEConfigForTimeoutTest("COLUMN_MODEL_1D_GRM", timeout, polyDeg, nElem);

	driver.configure(pp);
	driver.run();

	// test that whole solution is written
	cadet::InternalStorageUnitOpRecorder const* const simData = driver.solution()->unitOperation(0);
	REQUIRE(simData->numDataPoints() == 1501);
}

} // namespace test
} // namespace cadet
