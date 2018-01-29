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
#include "JsonTestModels.hpp"
#include "SimHelper.hpp"
#include "common/Driver.hpp"

#include <cmath>
#include <functional>
#include <vector>
#include <algorithm>
#include <iterator>

inline Approx makeApprox(double val, double relTol, double absTol)
{
	return Approx(val).epsilon(relTol).margin(absTol);
}

inline void setFlowRateFilter(cadet::JsonParameterProvider& jpp, double filter)
{
	jpp.pushScope("model");
	jpp.pushScope("unit_000");

	jpp.set("FLOWRATE_FILTER", filter);

	jpp.popScope();
	jpp.popScope();
}

	inline void runSim(cadet::JsonParameterProvider& jpp, std::function<double(double)> solC, std::function<double(double)> solV)
{
	// Run simulation
	cadet::Driver drv;
	drv.configure(jpp);
	drv.run();

	// Get data from simulation
	cadet::InternalStorageUnitOpRecorder const* const simData = drv.solution()->unitOperation(0);
	double const* outlet = simData->outlet();
	double const* volume = simData->volume();
	double const* time = drv.solution()->time();

	// Compare
	for (unsigned int i = 0; i < simData->numDataPoints(); ++i, ++outlet, ++volume, ++time)
	{
		// Compare with relative error 1e-6 and absolute error 4e-5
		CAPTURE(*time);
		CHECK((*outlet) == makeApprox(solC(*time), 1e-6, 4e-5));
		CHECK((*volume) == makeApprox(solV(*time), 1e-6, 4e-5));
	}
}

inline void runSim(cadet::JsonParameterProvider& jpp, std::function<double(double)> solC, std::function<double(double)> solQ, std::function<double(double)> solV, double absTol = 4e-5, double relTol = 1e-6)
{
	// Run simulation
	cadet::Driver drv;
	drv.configure(jpp);
	drv.run();

	// Get data from simulation
	cadet::InternalStorageUnitOpRecorder const* const simData = drv.solution()->unitOperation(0);
	double const* outlet = simData->outlet();
	double const* volume = simData->volume();
	double const* solid = simData->solid();
	double const* time = drv.solution()->time();

	// Compare
	for (unsigned int i = 0; i < simData->numDataPoints(); ++i, ++outlet, ++volume, ++solid, ++time)
	{
		CAPTURE(*time);
		CHECK((*outlet) == makeApprox(solC(*time), relTol, absTol));
		CHECK((*volume) == makeApprox(solV(*time), relTol, absTol));
		CHECK((*solid) == makeApprox(solQ(*time), relTol, absTol));
	}
}

inline void runSensSim(cadet::JsonParameterProvider& jpp, std::function<double(double)> solC, std::function<double(double)> solV, double absTol = 1e-5, double relTol = 1e-8)
{
	// Run simulation
	cadet::Driver drv;
	drv.configure(jpp);
	drv.run();

	// Get data from simulation
	cadet::InternalStorageUnitOpRecorder const* const simData = drv.solution()->unitOperation(0);
	double const* outlet = simData->sensOutlet(0);
	double const* volume = simData->sensVolume(0);
	double const* time = drv.solution()->time();

	// Compare
	for (unsigned int i = 0; i < simData->numDataPoints(); ++i, ++outlet, ++volume, ++time)
	{
		CAPTURE(*time);
		CHECK((*outlet) == makeApprox(solC(*time), relTol, absTol));
		CHECK((*volume) == makeApprox(solV(*time), relTol, absTol));
	}
}

TEST_CASE("CSTR vs analytic solution (V constant) w/o binding model", "[CSTR],[Simulation]")
{
	cadet::JsonParameterProvider jpp = createCSTRBenchmark(3, 119.0, 1.0);
	cadet::test::setSectionTimes(jpp, {0.0, 10.0, 100.0, 119.0});
	cadet::test::setInitialConditions(jpp, {0.0}, {}, 10.0);
	cadet::test::setInletProfile(jpp, 0, 0, 1.0, 0.0, 0.0, 0.0);
	cadet::test::setInletProfile(jpp, 1, 0, 1.0, -1.0 / 90.0, 0.0, 0.0);
	cadet::test::setInletProfile(jpp, 2, 0, 0.0, 0.0, 0.0, 0.0);
	cadet::test::setFlowRates(jpp, 0, 1.0, 0.5, 0.5);
	cadet::test::setFlowRates(jpp, 1, 1.0, 0.5, 0.5);
	cadet::test::setFlowRates(jpp, 2, 1.0, 0.5, 0.5);

	const double temp = 10.0 * (9.0 + 2.0 * std::sqrt(std::exp(1.0)));
	const double temp2 = 2.0 / 9.0 * (-9.0 - 2.0 * std::sqrt(std::exp(1.0)) + 2 * std::exp(5));
	runSim(jpp, [=](double t) {
			if (t <= 10.0)
				return -2.0 * std::expm1(-t / 20.0);
			else if (t <= 100.0)
				return (120.0 - temp * std::exp(-t / 20.0) - t)  / 45.0;
			else
				return std::exp(-5.0 - (t - 100.0) / 20.0) * temp2;
		}, 
		[](double t) {
			return 10.0;
	});
}

TEST_CASE("CSTR vs analytic solution (V increasing) w/o binding model", "[CSTR],[Simulation]")
{
	cadet::JsonParameterProvider jpp = createCSTRBenchmark(1, 100.0, 1.0);
	cadet::test::setSectionTimes(jpp, {0.0, 100.0});
	cadet::test::setInitialConditions(jpp, {1.0}, {}, 10.0);
	cadet::test::setInletProfile(jpp, 0, 0, 1.0, 0.0, 0.0, 0.0);
	cadet::test::setFlowRates(jpp, 0, 2.0, 1.0, 0.5);

	runSim(jpp, [=](double t) {
			return 4.0 * (6000.0 + t * (1200.0 + t * (60.0 + t))) / (3.0 * std::pow(20.0 + t, 3.0));
		}, 
		[](double t) {
			return 10.0 + 0.5 * t;
	});
}

TEST_CASE("CSTR vs analytic solution (V decreasing) w/o binding model", "[CSTR],[Simulation]")
{
	cadet::JsonParameterProvider jpp = createCSTRBenchmark(1, 100.0, 1.0);
	cadet::test::setSectionTimes(jpp, {0.0, 100.0});
	cadet::test::setInitialConditions(jpp, {1.0}, {}, 10.0);
	cadet::test::setInletProfile(jpp, 0, 0, 1.0, 0.0, 0.0, 0.0);
	cadet::test::setFlowRates(jpp, 0, 1.5, 1.5, 0.5);

	runSim(jpp, [=](double t) {
			return 1.0 + t * (1.0 / 20.0 - t / 800.0);
		}, 
		[](double t) {
			return 10.0 - 0.5 * t;
	});
}

TEST_CASE("CSTR vs analytic solution (V constant) with dynamic linear binding", "[CSTR],[Simulation]")
{
	cadet::JsonParameterProvider jpp = createCSTRBenchmark(1, 100.0, 1.0);
	cadet::test::setSectionTimes(jpp, {0.0, 100.0});
	cadet::test::addBoundStates(jpp, {1}, 0.5);
	cadet::test::setInitialConditions(jpp, {0.0}, {0.0}, 1.0);
	cadet::test::setInletProfile(jpp, 0, 0, 1.0, 0.0, 0.0, 0.0);
	cadet::test::setFlowRates(jpp, 0, 0.1, 0.1, 0.0);
	cadet::test::addLinearBindingModel(jpp, true, {0.1}, {10.0});

	const double sqrt2501 = std::sqrt(2501.0);
	runSim(jpp, [=](double t) {
			return 1.0 - std::exp(-5.1 * t) * (2501.0 * std::cosh(sqrt2501 * t / 10.0) + 50.0 * sqrt2501 * std::sinh(sqrt2501 * t / 10.0)) / 2501.0;
		}, 
		[=](double t) {
			return 0.01 - std::exp(-5.1 * t) * (2501.0 * std::cosh(sqrt2501 * t / 10.0) + 51.0 * sqrt2501 * std::sinh(sqrt2501 * t / 10.0)) / 250100.0;
		}, 
		[](double t) {
			return 1.0;
	});
}

TEST_CASE("CSTR vs analytic solution (V constant) with quasi-stationary linear binding", "[CSTR],[Simulation]")
{
	cadet::JsonParameterProvider jpp = createCSTRBenchmark(1, 100.0, 1.0);
	cadet::test::setSectionTimes(jpp, {0.0, 100.0});
	cadet::test::addBoundStates(jpp, {1}, 0.5);
	cadet::test::setInitialConditions(jpp, {0.0}, {0.0}, 1.0);
	cadet::test::setInletProfile(jpp, 0, 0, 1.0, 0.0, 0.0, 0.0);
	cadet::test::setFlowRates(jpp, 0, 0.1, 0.1, 0.0);
	cadet::test::addLinearBindingModel(jpp, false, {0.1}, {10.0});

	const double sqrt2501 = std::sqrt(2501.0);
	runSim(jpp, [=](double t) {
			return -std::expm1(-10.0 / 101.0 * t);
		}, 
		[=](double t) {
			return -std::expm1(-10.0 / 101.0 * t) * 0.01;
		}, 
		[](double t) {
			return 1.0;
	});
}

TEST_CASE("CSTR filter flowrate sensitivity vs analytic solution (V constant) w/o binding model", "[CSTR],[Simulation],[AD],[Sensitivity]")
{
	cadet::JsonParameterProvider jpp = createCSTRBenchmark(3, 119.0, 1.0);
	cadet::test::setSectionTimes(jpp, {0.0, 10.0, 100.0, 119.0});
	cadet::test::setInitialConditions(jpp, {0.0}, {}, 10.0);
	cadet::test::setInletProfile(jpp, 0, 0, 1.0, 0.0, 0.0, 0.0);
	cadet::test::setInletProfile(jpp, 1, 0, 1.0, -1.0 / 90.0, 0.0, 0.0);
	cadet::test::setInletProfile(jpp, 2, 0, 0.0, 0.0, 0.0, 0.0);
	cadet::test::setFlowRates(jpp, 0, 1.0, 0.5, 0.5);
	cadet::test::setFlowRates(jpp, 1, 1.0, 0.5, 0.5);
	cadet::test::setFlowRates(jpp, 2, 1.0, 0.5, 0.5);
	setFlowRateFilter(jpp, 0.5);
	cadet::test::addSensitivity(jpp, "FLOWRATE_FILTER", cadet::makeParamId("FLOWRATE_FILTER", 0, cadet::CompIndep, cadet::BoundPhaseIndep, cadet::ReactionIndep, cadet::SectionIndep), 1e-6);
	cadet::test::returnSensitivities(jpp, 0);

	const double sqrtE = std::sqrt(std::exp(1.0));
	runSensSim(jpp, [=](double t) {
			if (t <= 10.0)
				return 4.0 + 1.0/200.0 * std::exp(-t / 20.0) * (-800.0 + (-40.0 + t) * t);
			else if (t <= 100.0)
				return -160.0 * (-80.0 + t) / 1800.0 + std::exp(-t / 20.0) * (-4.0 + (9.0 * (-40.0 + t) * t + 2 * sqrtE * (-1700.0 + (-40.0 + t) * t)) / 1800.0);
			else
				return std::exp(-t / 20.0) * (-4.0 + (9.0 * (-40.0 + t) * t + std::exp(5.0) * (8800.0 - 2.0 * (-40.0 + t) * t) + 2.0 * sqrtE * (-1700.0 + (-40.0 + t) * t)) / 1800.0);
		}, 
		[](double t) {
			return -t;
	});
}

TEST_CASE("CSTR LIN_COEFF sensitivity vs analytic solution (V constant) w/o binding model", "[CSTR],[Simulation],[AD],[Sensitivity]")
{
	cadet::JsonParameterProvider jpp = createCSTRBenchmark(1, 100.0, 1.0);
	cadet::test::setSectionTimes(jpp, {0.0, 100.0});
	cadet::test::setInitialConditions(jpp, {0.0}, {}, 10.0);
	cadet::test::setInletProfile(jpp, 0, 0, 1.0, 1.0, 0.0, 0.0);
	cadet::test::setFlowRates(jpp, 0, 1.0, 0.5, 0.5);
	cadet::test::addSensitivity(jpp, "LIN_COEFF", cadet::makeParamId("LIN_COEFF", 1, 0, cadet::BoundPhaseIndep, cadet::ReactionIndep, 0), 1e-6);
	cadet::test::returnSensitivities(jpp, 0);

	runSensSim(jpp, [=](double t) { return 2.0 * (20.0 * std::expm1(-t / 20.0) + t); }, [](double t) { return 0.0; }, 1e-5, 7e-8);
}
