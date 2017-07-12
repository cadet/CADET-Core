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

#include "model/StirredTankModel.hpp"
#include "ModelBuilderImpl.hpp"

#include "JsonParameterProvider.hpp"
#include "JacobianHelper.hpp"

#include <cmath>
#include <functional>

/**
 * @brief Fills the state vector with a given function
 * @details The function @p f uses the current index to assign a value.
 * @param [out] y Filled state vector
 * @param [in] f Function for computing the content of the state vector
 * @param [in] numDofs Size of the state vector
 */
inline void fillState(double* y, std::function<double(unsigned int)> f, unsigned int numDofs)
{
	for (unsigned int i = 0; i < numDofs; ++i)
		y[i] = f(i);
}

/**
 * @brief Creates a runnable CSTR model with given WENO order
 * @details Creates a CSTR model and configures it using the given IParameterProvider @p jpp.
 * @param [in] mb ModelBuilder
 * @param [in] jpp Configuration of the CSTR
 * @return Runnable CSTR
 */
cadet::model::CSTRModel* createAndConfigureCSTR(cadet::IModelBuilder& mb, cadet::JsonParameterProvider& jpp)
{
	// Create a CSTR
	cadet::IModel* const iCstr = mb.createUnitOperation("CSTR", 0);
	REQUIRE(nullptr != iCstr);

	cadet::model::CSTRModel* const cstr = reinterpret_cast<cadet::model::CSTRModel*>(iCstr);

	// Configure
	cadet::ModelBuilder& temp = *reinterpret_cast<cadet::ModelBuilder*>(&mb);
	REQUIRE(cstr->configure(jpp, temp));

	return cstr;
}

inline void checkJacobian(double flowRateIn, double flowRateOut, double flowRateFilter)
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	// CSTR with 2 components
	const unsigned int nComp = 2;
	cadet::JsonParameterProvider jpp = createCSTR(nComp);

	cadet::model::CSTRModel* const cstr = createAndConfigureCSTR(*mb, jpp);
	if (flowRateFilter > 0.0)
		cstr->flowRateFilter().push_back(flowRateFilter);

	cstr->setFlowRates(flowRateIn, flowRateOut);

	// Setup matrices
	cstr->notifyDiscontinuousSectionTransition(0.0, 0u, nullptr, nullptr, 0u);

	// Obtain memory for state, Jacobian multiply direction, Jacobian column
	const unsigned int nDof = cstr->numDofs();
	std::vector<double> y(nDof, 0.0);
	std::vector<double> yDot(nDof, 0.0);
	std::vector<double> jacDir(nDof, 0.0);
	std::vector<double> jacCol1(nDof, 0.0);
	std::vector<double> jacCol2(nDof, 0.0);

	// Fill state vectors with some values
	fillState(y.data(), [=](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, nDof);
	fillState(yDot.data(), [=](unsigned int idx) { return std::abs(std::sin((idx + nDof) * 0.13)) + 1e-4; }, nDof);

	// Compare Jacobians
//	cadet::test::checkJacobianPatternFD(cstr, cstr, y.data(), yDot.data(), jacDir.data(), jacCol1.data(), jacCol2.data());
	cadet::test::compareJacobianFD(cstr, cstr, y.data(), yDot.data(), jacDir.data(), jacCol1.data(), jacCol2.data());

	mb->destroyUnitOperation(cstr);
	destroyModelBuilder(mb);	
}

inline void checkTimeDerivativeJacobian(double flowRateIn, double flowRateOut, double flowRateFilter)
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	// CSTR with 2 components
	const unsigned int nComp = 2;
	cadet::JsonParameterProvider jpp = createCSTR(nComp);

	cadet::model::CSTRModel* const cstr = createAndConfigureCSTR(*mb, jpp);
	if (flowRateFilter > 0.0)
		cstr->flowRateFilter().push_back(flowRateFilter);

	cstr->setFlowRates(flowRateIn, flowRateOut);

	// Setup matrices
	cstr->notifyDiscontinuousSectionTransition(0.0, 0u, nullptr, nullptr, 0u);

	// Obtain memory for state, Jacobian multiply direction, Jacobian column
	const unsigned int nDof = cstr->numDofs();
	std::vector<double> y(nDof, 0.0);
	std::vector<double> yDot(nDof, 0.0);
	std::vector<double> jacDir(nDof, 0.0);
	std::vector<double> jacCol1(nDof, 0.0);
	std::vector<double> jacCol2(nDof, 0.0);

	// Fill state vectors with some values
	fillState(y.data(), [=](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, nDof);
	fillState(yDot.data(), [=](unsigned int idx) { return std::abs(std::sin((idx + nDof) * 0.13)) + 1e-4; }, nDof);

	// Compare Jacobians
	cadet::test::compareTimeDerivativeJacobianFD(cstr, cstr, y.data(), yDot.data(), jacDir.data(), jacCol1.data(), jacCol2.data());

	mb->destroyUnitOperation(cstr);
	destroyModelBuilder(mb);
}

TEST_CASE("StirredTankModel Jacobian", "[CSTR],[UnitOp],[Residual],[Jacobian]")
{
	const double rateList[] = {0.0, 0.0, 0.0,
	                           1.0, 0.0, 0.0,
	                           0.0, 1.0, 0.0,
	                           0.0, 0.0, 1.0,
	                           1.0, 1.0, 0.0,
	                           1.0, 0.0, 1.0,
	                           0.0, 1.0, 1.0,
	                           1.0, 1.0, 1.0,
	                           1.0, 2.0, 3.0,
	                           2.0, 3.0, 1.0,
	                           3.0, 1.0, 2.0,
	                           3.0, 2.0, 1.0,
	                           2.0, 1.0, 3.0,
	                           1.0, 3.0, 2.0,
	                           1.0, 2.0, 0.0,
	                           2.0, 1.0, 0.0};
	double const* row = &rateList[0];
	for (unsigned int i = 0; i < sizeof(rateList) / sizeof(double) / 3; ++i, row += 3)
	{
		SECTION("Fin = " + std::to_string(row[0]) + " Fout = " + std::to_string(row[1]) + " Ffilter = " + std::to_string(row[2]))
		{
			checkJacobian(row[0], row[1], row[2]);
		}
	}
}

TEST_CASE("StirredTankModel time derivative Jacobian", "[CSTR],[UnitOp],[Residual],[Jacobian]")
{
	const double rateList[] = {0.0, 0.0, 0.0,
	                           1.0, 0.0, 0.0,
	                           0.0, 1.0, 0.0,
	                           0.0, 0.0, 1.0,
	                           1.0, 1.0, 0.0,
	                           1.0, 0.0, 1.0,
	                           0.0, 1.0, 1.0,
	                           1.0, 1.0, 1.0,
	                           1.0, 2.0, 3.0,
	                           2.0, 3.0, 1.0,
	                           3.0, 1.0, 2.0,
	                           3.0, 2.0, 1.0,
	                           2.0, 1.0, 3.0,
	                           1.0, 3.0, 2.0,
	                           1.0, 2.0, 0.0,
	                           2.0, 1.0, 0.0};
	double const* row = &rateList[0];
	for (unsigned int i = 0; i < sizeof(rateList) / sizeof(double) / 3; ++i, row += 3)
	{
		SECTION("Fin = " + std::to_string(row[0]) + " Fout = " + std::to_string(row[1]) + " Ffilter = " + std::to_string(row[2]))
		{
			checkTimeDerivativeJacobian(row[0], row[1], row[2]);
		}
	}
}
