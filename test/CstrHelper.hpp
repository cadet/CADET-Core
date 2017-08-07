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

/**
 * @file 
 * Defines helper functions for building CSTR models.
 */

#ifndef CADETTEST_CSTRHELPER_HPP_
#define CADETTEST_CSTRHELPER_HPP_

#include "JsonParameterProvider.hpp"

#include <sstream>
#include <iomanip>

namespace cadet
{

namespace test
{

	/**
	 * @brief Sets the initial conditions of a CSTR
	 * @details Overwrites the INIT_C and INIT_VOLUME fields of the given ParameterProvider.
	 * @param [in,out] jpp ParameterProvider
	 * @param [in] c Initial liquid phase concentration
	 * @param [in] q Initial bound phase concentration
	 * @param [in] v Initial volume
	 */
	inline void setInitialConditions(cadet::JsonParameterProvider& jpp, const std::vector<double>& c, const std::vector<double>& q, double v)
	{
		std::ostringstream ss;

		const bool inSystem = jpp.exists("model");
		if (inSystem)
		{
			jpp.pushScope("model");
			jpp.pushScope("unit_000");
		}

		jpp.set("INIT_C", c);
		jpp.set("INIT_VOLUME", v);
		if (!q.empty())
			jpp.set("INIT_Q", q);

		if (inSystem)
		{
			jpp.popScope();
			jpp.popScope();
		}
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

	/**
	 * @brief Adds bound states to a CSTR model
	 * @param [in] jpp ParameterProvider
	 * @param [in] nBound Array with number of bound states for each component
	 * @param [in] porosity Porosity
	 */
	inline void addBoundStates(cadet::JsonParameterProvider& jpp, const std::vector<int>& nBound, double porosity)
	{
		const bool inSystem = jpp.exists("model");
		if (inSystem)
		{
			jpp.pushScope("model");
			jpp.pushScope("unit_000");
		}

		jpp.set("NBOUND", nBound);
		jpp.set("POROSITY", porosity);

		if (inSystem)
		{
			jpp.popScope();
			jpp.popScope();
		}
	}

	/**
	 * @brief Adds linear binding model to a CSTR model
	 * @param [in] jpp ParameterProvider
	 * @param [in] kinetic Determines whether kinetic or quasi-stationary binding is used
	 * @param [in] kA Vector with kA rates
	 * @param [in] kD Vector with kD rates
	 */
	inline void addLinearBindingModel(cadet::JsonParameterProvider& jpp, bool kinetic, const std::vector<double>& kA, const std::vector<double>& kD)
	{
		const bool inSystem = jpp.exists("model");
		if (inSystem)
		{
			jpp.pushScope("model");
			jpp.pushScope("unit_000");
		}

		jpp.set("ADSORPTION_MODEL", "LINEAR");

		jpp.addScope("adsorption");
		jpp.pushScope("adsorption");

		jpp.set("IS_KINETIC", kinetic);
		jpp.set("LIN_KA", kA);
		jpp.set("LIN_KD", kD);

		jpp.popScope();
		if (inSystem)
		{
			jpp.popScope();
			jpp.popScope();
		}
	}

	/**
	 * @brief Adds Langmuir binding model to a CSTR model
	 * @param [in] jpp ParameterProvider
	 * @param [in] kinetic Determines whether kinetic or quasi-stationary binding is used
	 * @param [in] kA Vector with kA rates
	 * @param [in] kD Vector with kD rates
	 * @param [in] qMax Vector with maximum capacities qMax
	 */
	inline void addLangmuirBindingModel(cadet::JsonParameterProvider& jpp, bool kinetic, const std::vector<double>& kA, const std::vector<double>& kD, const std::vector<double>& qMax)
	{
		const bool inSystem = jpp.exists("model");
		if (inSystem)
		{
			jpp.pushScope("model");
			jpp.pushScope("unit_000");
		}

		jpp.set("ADSORPTION_MODEL", "MULTI_COMPONENT_LANGMUIR");

		jpp.addScope("binding");
		jpp.pushScope("binding");

		jpp.set("IS_KINETIC", kinetic);
		jpp.set("MCL_KA", kA);
		jpp.set("MCL_KD", kD);
		jpp.set("MCL_QMAX", qMax);

		jpp.popScope();
		if (inSystem)
		{
			jpp.popScope();
			jpp.popScope();
		}
	}

} // namespace test
} // namespace cadet

#endif  // CADETTEST_CSTRHELPER_HPP_
