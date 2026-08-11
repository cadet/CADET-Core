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
#include <thread>
#include <json.hpp>
using json = nlohmann::json;

#include "Logging.hpp"
#include "cadet/Notification.hpp"
#include "common/Driver.hpp"
#include "common/JsonParameterProvider.hpp"
#include "common/TimeoutCallback.hpp"
#include "JsonTestModels.hpp"

namespace cadet
{
namespace test
{

namespace
{
	/**
	 * @brief Creates a LWE configuration with DG discretization for timeout testing
	 * @param [in] uoType Unit operation type
	 * @param [in] timeout timeout
	 * @param [in] polyDeg Polynomial degree of the DG discretization
	 * @param [in] nElem Number of DG elements
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

		// Set discretization to which the required compute time relates
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

TEST_CASE("Test Callback: timeout stops time integration at every checkpoint", "[Callback],[Timeout],[CI_callback]")
{
	cadet::TimeoutCallback callback;
	callback.setTimeout(0.01);
	callback.timeIntegrationStart();

	// Sleeping is guaranteed to take at least the requested time, so the timeout has surely elapsed
	std::this_thread::sleep_for(std::chrono::milliseconds(100));

	REQUIRE(!callback.timeIntegrationSection(0, 0.0, nullptr, nullptr, 0.0));
	REQUIRE(!callback.timeIntegrationStep(0, 0.0, nullptr, nullptr, 0.0));
	REQUIRE(!callback.timeIntegrationLinearSolve(0, 0.0, nullptr, nullptr));
}

TEST_CASE("Test Callback: timeout lets time integration continue before it has elapsed", "[Callback],[Timeout],[CI_callback]")
{
	cadet::TimeoutCallback callback;
	callback.setTimeout(3600.0);
	callback.timeIntegrationStart();

	REQUIRE(callback.timeIntegrationSection(0, 0.0, nullptr, nullptr, 0.0));
	REQUIRE(callback.timeIntegrationStep(0, 0.0, nullptr, nullptr, 0.0));
	REQUIRE(callback.timeIntegrationLinearSolve(0, 0.0, nullptr, nullptr));
}

TEST_CASE("Test Callback: timeout grants every time integration the full time budget", "[Callback],[Timeout],[CI_callback]")
{
	cadet::TimeoutCallback callback;
	callback.setTimeout(0.01);

	callback.timeIntegrationStart();
	std::this_thread::sleep_for(std::chrono::milliseconds(100));
	REQUIRE(!callback.timeIntegrationStep(0, 0.0, nullptr, nullptr, 0.0));

	// A time integration that was stopped by the timeout never reaches timeIntegrationEnd(), so the
	// time spent on it must not be charged against the time budget of the following time integration
	callback.timeIntegrationStart();
	REQUIRE(callback.timeIntegrationStep(0, 0.0, nullptr, nullptr, 0.0));
}

TEST_CASE("Test Callback: timeout interrupts simulation but data is saved", "[Callback],[Timeout],[CI_callback]")
{
	const int polyDeg = 3;
	const int nElem = 16;

	// Reference run without timeout. Deriving the timeout from the measured duration keeps the test
	// independent of the speed of the machine it runs on
	Driver refDriver;
	cadet::JsonParameterProvider refPP = createLWEConfigForTimeoutTest("COLUMN_MODEL_1D_GRM", 0.0, polyDeg, nElem);

	refDriver.configure(refPP);
	refDriver.run();

	const double refDuration = refDriver.simulator()->lastSimulationDuration();
	REQUIRE(refDriver.solution()->unitOperation(0)->numDataPoints() == 1501);

	// The timeout is hit long before the simulation is finished, even if this run happens to be
	// substantially faster than the reference run
	const double timeout = 0.1 * refDuration;

	Driver driver;
	cadet::JsonParameterProvider pp = createLWEConfigForTimeoutTest("COLUMN_MODEL_1D_GRM", timeout, polyDeg, nElem);

	driver.configure(pp);
	driver.run();

	REQUIRE(driver.simulator()->stoppedByNotificationCallback());
	// The timeout is only checked in between linear solves, so the simulation may overshoot it by
	// the duration of a single linear solve
	REQUIRE(driver.simulator()->lastSimulationDuration() <= timeout + 0.1 * refDuration);
	// test that part of the solution is written
	REQUIRE(driver.solution()->unitOperation(0)->numDataPoints() > 1);
	REQUIRE(driver.solution()->unitOperation(0)->numDataPoints() < 1501);
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
