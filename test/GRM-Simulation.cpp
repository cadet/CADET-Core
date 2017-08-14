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
#include <cstdio>

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

	/**
	 * @brief Reads reference chromatograms from a test data file
	 * @details The file format is as follows:
	 *          Number of data points (uint32)
	 *          Time points (array of doubles)
	 *          Chromatogram for dynamic binding (array of doubles)
	 *          Chromatogram for quasi-stationary binding (array of doubles)
	 */
	class ReferenceDataReader
	{
	public:
		ReferenceDataReader(const char* fileName) : _f(nullptr)
		{
			_f = std::fopen(fileName, "rb");
			std::fread(&_numElements, 4, 1, _f);
		}

		std::vector<double> time()
		{
			std::vector<double> v(_numElements, 0.0);
			std::fseek(_f, 4, SEEK_SET);
			std::fread(v.data(), 8, _numElements, _f);
			return v;
		}

		std::vector<double> analyticDynamic()
		{
			std::vector<double> v(_numElements, 0.0);
			std::fseek(_f, 4 + _numElements * 8, SEEK_SET);
			std::fread(v.data(), 8, _numElements, _f);
			return v;
		}

		std::vector<double> analyticQuasiStationary()
		{
			std::vector<double> v(_numElements, 0.0);
			std::fseek(_f, 4 + 2 * _numElements * 8, SEEK_SET);
			std::fread(v.data(), 8, _numElements, _f);
			return v;
		}

	private:
		std::FILE* _f;
		uint32_t _numElements;
	};
}

/**
 * @brief Returns the absolute path to the test/ folder of the project
 * @details Absolute path to the test/ folder of the project without trailing slash
 * @return Absolute path to the test/ folder
 */
const char* getTestDirectory();

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

void testAnalyticBenchmark(bool forwardFlow, bool dynamicBinding)
{
	const std::string fwdStr = (forwardFlow ? "forward" : "backward");
	SECTION("Analytic" + fwdStr + " flow with " + (dynamicBinding ? "dynamic" : "quasi-stationary") + " binding")
	{
		// Setup simulation
		cadet::JsonParameterProvider jpp = createLinearBenchmark(dynamicBinding, false);
		if (!forwardFlow)
			reverseFlow(jpp);

		// Run simulation
		cadet::Driver drv;
		drv.configure(jpp);
		drv.run();

		// Read reference data from test file
		const std::string refFile = std::string(getTestDirectory()) + "/data/grm-pulseBenchmark.data";
		ReferenceDataReader rd(refFile.c_str());
		const std::vector<double> time = rd.time();
		const std::vector<double> ref = (dynamicBinding ? rd.analyticDynamic() : rd.analyticQuasiStationary());

		// Get data from simulation
		cadet::InternalStorageUnitOpRecorder const* const simData = drv.solution()->unitOperation(0);
		double const* outlet = (forwardFlow ? simData->outlet() : simData->inlet());

		// Compare
		for (unsigned int i = 0; i < simData->numDataPoints(); ++i, ++outlet)
		{
			// Compare with relative error 1e-6 and absolute error 4e-5

			// Note that the simulation only saves the chromatogram at multiples of 2 (i.e., 0s, 2s, 4s, ...)
			// whereas the reference solution is given at every second (0s, 1s, 2s, 3s, ...)
			// Thus, we only take the even indices of the reference array
			CAPTURE(time[2 * i]);
			CHECK((*outlet) == makeApprox(ref[2 * i], 1e-6, 4e-5));
		}
	}
}

inline void testAnalyticNonBindingBenchmark(bool forwardFlow)
{
	const std::string fwdStr = (forwardFlow ? "forward" : "backward");
	SECTION("Analytic" + fwdStr + " flow")
	{
		// Setup simulation
		cadet::JsonParameterProvider jpp = createLinearBenchmark(true, true);
		if (!forwardFlow)
			reverseFlow(jpp);

		// Run simulation
		cadet::Driver drv;
		drv.configure(jpp);
		drv.run();

		// Read reference data from test file
		const std::string refFile = std::string(getTestDirectory()) + "/data/grm-nonBinding.data";
		ReferenceDataReader rd(refFile.c_str());
		const std::vector<double> time = rd.time();
		const std::vector<double> ref = rd.analyticDynamic();

		// Get data from simulation
		cadet::InternalStorageUnitOpRecorder const* const simData = drv.solution()->unitOperation(0);
		double const* outlet = (forwardFlow ? simData->outlet() : simData->inlet());

		// Compare
		for (unsigned int i = 0; i < simData->numDataPoints(); ++i, ++outlet)
		{
			// Compare with relative error 1e-6 and absolute error 6e-5
			CAPTURE(time[i]);
			CHECK((*outlet) == makeApprox(ref[i], 1e-6, 6e-5));
		}
	}
}

TEST_CASE("LWE forward vs backward flow", "[GRM],[Simulation]")
{
	// Test all WENO orders
	for (unsigned int i = 1; i <= cadet::Weno::maxOrder(); ++i)
		testWenoForwardBackward(i);
}

TEST_CASE("GRM linear pulse vs analytic solution", "[GRM],[Simulation],[Analytic]")
{
	testAnalyticBenchmark(true, true);
	testAnalyticBenchmark(true, false);
	testAnalyticBenchmark(false, true);
	testAnalyticBenchmark(false, false);
}

TEST_CASE("GRM non-binding linear pulse vs analytic solution", "[GRM],[Simulation],[Analytic],[NonBinding]")
{
	testAnalyticNonBindingBenchmark(true);
	testAnalyticNonBindingBenchmark(false);
}
