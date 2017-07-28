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

#include <cmath>
#include <functional>
#include <vector>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <iomanip>

namespace
{
	/**
	 * @brief Sets the initial conditions of a CSTR
	 * @details Overwrites the INIT_C and INIT_VOLUME fields of the given ParameterProvider.
	 * @param [in,out] jpp ParameterProvider
	 * @param [in] c Initial concentration
	 * @param [in] v Initial volume
	 */
	inline void setInitialConditions(cadet::JsonParameterProvider& jpp, const std::vector<double>& c, double v)
	{
		std::ostringstream ss;

		jpp.pushScope("model");
		jpp.pushScope("unit_000");

		jpp.set("INIT_C", c);
		jpp.set("INIT_VOLUME", v);

		jpp.popScope();
		jpp.popScope();
	}

	/**
	 * @brief Sets the flow rates of a section
	 * @details Sets inflow and outflow as well as filter flow rate.
	 * @param [in,out] jpp ParameterProvider
	 * @param [in] secIdx Section index
	 * @param [in] in Inflow rate
	 * @param [in] out Outflow rate
	 * @param [in] filter Filter flow rate
	 */
	inline void setFlowRates(cadet::JsonParameterProvider& jpp, unsigned int secIdx, double in, double out, double filter)
	{
		std::ostringstream ss;

		jpp.pushScope("model");
		jpp.pushScope("unit_000");

		std::vector<double> frf = jpp.getDoubleArray("FLOWRATE_FILTER");
		if (frf.size() <= secIdx)
			std::fill_n(std::back_inserter(frf), frf.size() - secIdx + 1, 0.0);

		frf[secIdx] = filter;
		jpp.set("FLOWRATE_FILTER", frf);

		jpp.popScope();
		jpp.pushScope("connections");

		ss << "switch_" << std::setfill('0') << std::setw(3) << secIdx;
		jpp.pushScope(ss.str());

		std::vector<double> con = jpp.getDoubleArray("CONNECTIONS");
		con[4] = in;
		con[9] = out;
		jpp.set("CONNECTIONS", con);

		jpp.popScope();
		jpp.popScope();
		jpp.popScope();
	}

	/**
	 * @brief Sets the inlet profile of a section and component
	 * @details Overwrites the spline piece of the given section and component.
	 * @param [in,out] jpp ParameterProvider
	 * @param [in] secIdx Section index
	 * @param [in] comp Component index
	 * @param [in] con Constant polynomial coefficient
	 * @param [in] lin Linear polynomial coefficient
	 * @param [in] quad Quadratic polynomial coefficient
	 * @param [in] cub Cubic polynomial coefficient
	 */
	inline void setInletProfile(cadet::JsonParameterProvider& jpp, unsigned int secIdx, unsigned int comp, double con, double lin, double quad, double cub)
	{
		std::ostringstream ss;

		jpp.pushScope("model");
		jpp.pushScope("unit_001");

		ss << "sec_" << std::setfill('0') << std::setw(3) << secIdx;
		jpp.pushScope(ss.str());

		std::vector<double> cCon = jpp.getDoubleArray("CONST_COEFF");
		cCon[comp] = con;
		jpp.set("CONST_COEFF", cCon);

		std::vector<double> cLin = jpp.getDoubleArray("LIN_COEFF");
		cLin[comp] = lin;
		jpp.set("LIN_COEFF", cLin);

		std::vector<double> cQuad = jpp.getDoubleArray("QUAD_COEFF");
		cQuad[comp] = quad;
		jpp.set("QUAD_COEFF", cQuad);

		std::vector<double> cCube = jpp.getDoubleArray("CUBE_COEFF");
		cCube[comp] = cub;
		jpp.set("CUBE_COEFF", cCube);

		jpp.popScope();
		jpp.popScope();
		jpp.popScope();
	}

	/**
	 * @brief Sets the section times
	 * @details Overwrites the SECTION_TIMES field of the given ParameterProvider.
	 * @param [in,out] jpp ParameterProvider
	 * @param [in] secTimes Section times vector
	 */
	inline void setSectionTimes(cadet::JsonParameterProvider& jpp, const std::vector<double>& secTimes)
	{
		jpp.pushScope("solver");
		jpp.pushScope("sections");

		jpp.set("SECTION_TIMES", secTimes);

		jpp.popScope();
		jpp.popScope();
	}
}

inline Approx makeApprox(double val, double relTol, double absTol)
{
	return Approx(val).epsilon(relTol).margin(absTol);
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

TEST_CASE("CSTR vs analytic solution (V = const) w/o binding model", "[CSTR],[Simulation]")
{
	cadet::JsonParameterProvider jpp = createCSTRBenchmark(3, 119.0, 1.0);
	setSectionTimes(jpp, {0.0, 10.0, 100.0, 119.0});
	setInitialConditions(jpp, {0.0}, 10.0);
	setInletProfile(jpp, 0, 0, 1.0, 0.0, 0.0, 0.0);
	setInletProfile(jpp, 1, 0, 1.0, -1.0 / 90.0, 0.0, 0.0);
	setInletProfile(jpp, 2, 0, 0.0, 0.0, 0.0, 0.0);
	setFlowRates(jpp, 0, 1.0, 0.5, 0.5);
	setFlowRates(jpp, 1, 1.0, 0.5, 0.5);
	setFlowRates(jpp, 2, 1.0, 0.5, 0.5);

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
	setSectionTimes(jpp, {0.0, 100.0});
	setInitialConditions(jpp, {1.0}, 10.0);
	setInletProfile(jpp, 0, 0, 1.0, 0.0, 0.0, 0.0);
	setFlowRates(jpp, 0, 2.0, 1.0, 0.5);

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
	setSectionTimes(jpp, {0.0, 100.0});
	setInitialConditions(jpp, {1.0}, 10.0);
	setInletProfile(jpp, 0, 0, 1.0, 0.0, 0.0, 0.0);
	setFlowRates(jpp, 0, 1.5, 1.5, 0.5);

	runSim(jpp, [=](double t) {
			return 1.0 + t * (1.0 / 20.0 - t / 800.0);
		}, 
		[](double t) {
			return 10.0 - 0.5 * t;
	});
}
