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
#include "CstrHelper.hpp"

#include <cmath>
#include <functional>

namespace
{
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

inline void checkJacobianAD(double flowRateIn, double flowRateOut, double flowRateFilter, std::function<void(cadet::JsonParameterProvider&, unsigned int)> modelRefiner)
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	// CSTR with 2 components
	const unsigned int nComp = 2;
	cadet::JsonParameterProvider jpp = createCSTR(nComp);
	modelRefiner(jpp, nComp);

	cadet::model::CSTRModel* const cstrAna = createAndConfigureCSTR(*mb, jpp);
	cadet::model::CSTRModel* const cstrAD = createAndConfigureCSTR(*mb, jpp);
	if (flowRateFilter > 0.0)
	{
		cstrAna->flowRateFilter().push_back(flowRateFilter);
		cstrAD->flowRateFilter().push_back(flowRateFilter);
	}

	cstrAna->setFlowRates(flowRateIn, flowRateOut);
	cstrAD->setFlowRates(flowRateIn, flowRateOut);

	// Enable AD
	cadet::ad::setDirections(cadet::ad::getMaxDirections());
	cstrAD->useAnalyticJacobian(false);

	cadet::active* adRes = new cadet::active[cstrAD->numDofs()];
	cadet::active* adY = new cadet::active[cstrAD->numDofs()];

	cstrAD->prepareADvectors(adRes, adY, 0);

	// Setup matrices
	cstrAna->notifyDiscontinuousSectionTransition(0.0, 0u, nullptr, nullptr, 0u);
	cstrAD->notifyDiscontinuousSectionTransition(0.0, 0u, adRes, adY, 0u);

	// Obtain memory for state, Jacobian multiply direction, Jacobian column
	const unsigned int nDof = cstrAna->numDofs();
	std::vector<double> y(nDof, 0.0);
	std::vector<double> yDot(nDof, 0.0);
	std::vector<double> jacDir(nDof, 0.0);
	std::vector<double> jacCol1(nDof, 0.0);
	std::vector<double> jacCol2(nDof, 0.0);

	// Fill state vectors with some values
	fillState(y.data(), [=](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, nDof);
	fillState(yDot.data(), [=](unsigned int idx) { return std::abs(std::sin((idx + nDof) * 0.13)) + 1e-4; }, nDof);

	// Compute state Jacobian
	cstrAna->residualWithJacobian(0.0, 0u, 1.0, y.data(), nullptr, jacDir.data(), nullptr, nullptr, 0u);
	cstrAD->residualWithJacobian(0.0, 0u, 1.0, y.data(), nullptr, jacDir.data(), adRes, adY, 0u);
	std::fill(jacDir.begin(), jacDir.end(), 0.0);

	// Compare Jacobians
	cadet::test::compareJacobian(cstrAna, cstrAD, y.data(), yDot.data(), jacDir.data(), jacCol1.data(), jacCol2.data());

	delete[] adRes;
	delete[] adY;
	mb->destroyUnitOperation(cstrAna);
	mb->destroyUnitOperation(cstrAD);
	destroyModelBuilder(mb);	
}

inline void checkJacobianFD(double flowRateIn, double flowRateOut, double flowRateFilter, std::function<void(cadet::JsonParameterProvider&, unsigned int)> modelRefiner)
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	// CSTR with 2 components
	const unsigned int nComp = 2;
	cadet::JsonParameterProvider jpp = createCSTR(nComp);
	modelRefiner(jpp, nComp);

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

	// Compute state Jacobian
	cstr->residualWithJacobian(0.0, 0u, 1.0, y.data(), yDot.data(), jacDir.data(), nullptr, nullptr, 0u);
	std::fill(jacDir.begin(), jacDir.end(), 0.0);

	// Compare Jacobians
//	cadet::test::checkJacobianPatternFD(cstr, cstr, y.data(), yDot.data(), jacDir.data(), jacCol1.data(), jacCol2.data());
	cadet::test::compareJacobianFD(cstr, cstr, y.data(), yDot.data(), jacDir.data(), jacCol1.data(), jacCol2.data());

	mb->destroyUnitOperation(cstr);
	destroyModelBuilder(mb);	
}

inline void checkTimeDerivativeJacobianFD(double flowRateIn, double flowRateOut, double flowRateFilter, std::function<void(cadet::JsonParameterProvider&, unsigned int)> modelRefiner)
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	// CSTR with 2 components
	const unsigned int nComp = 2;
	cadet::JsonParameterProvider jpp = createCSTR(nComp);
	modelRefiner(jpp, nComp);

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

TEST_CASE("StirredTankModel Jacobian vs FD w/o binding model", "[CSTR],[UnitOp],[Residual],[Jacobian]")
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
			checkJacobianFD(row[0], row[1], row[2], [](cadet::JsonParameterProvider& jpp, unsigned int nComp) { });
		}
	}
}

TEST_CASE("StirredTankModel Jacobian vs AD w/o binding model", "[CSTR],[UnitOp],[Residual],[Jacobian],[AD]")
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
			checkJacobianAD(row[0], row[1], row[2], [](cadet::JsonParameterProvider& jpp, unsigned int nComp) { });
		}
	}
}

TEST_CASE("StirredTankModel time derivative Jacobian vs FD w/o binding model", "[CSTR],[UnitOp],[Residual],[Jacobian]")
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
			checkTimeDerivativeJacobianFD(row[0], row[1], row[2], [](cadet::JsonParameterProvider& jpp, unsigned int nComp) { });
		}
	}
}


TEST_CASE("StirredTankModel Jacobian vs FD with linear binding", "[CSTR],[UnitOp],[Residual],[Jacobian]")
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
	for (int j = 0; j < 2; ++j)
	{
		const bool dynamic = (j == 0);
		const std::string bndMode = dynamic ? " dynamic" : " quasi-stationary";
		double const* row = &rateList[0];
		for (unsigned int i = 0; i < sizeof(rateList) / sizeof(double) / 3; ++i, row += 3)
		{
			SECTION("Fin = " + std::to_string(row[0]) + " Fout = " + std::to_string(row[1]) + " Ffilter = " + std::to_string(row[2]) + bndMode)
			{
				checkJacobianFD(row[0], row[1], row[2], [=](cadet::JsonParameterProvider& jpp, unsigned int nComp) {
					cadet::test::addBoundStates(jpp, {1, 1}, 0.5);
					cadet::test::addLinearBindingModel(jpp, dynamic, {5.0, 4.0}, {2.0, 3.0});
				});
			}
		}
	}
}

TEST_CASE("StirredTankModel Jacobian vs AD with linear binding", "[CSTR],[UnitOp],[Residual],[Jacobian],[AD]")
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
	for (int j = 0; j < 2; ++j)
	{
		const bool dynamic = (j == 0);
		const std::string bndMode = dynamic ? " dynamic" : " quasi-stationary";
		double const* row = &rateList[0];
		for (unsigned int i = 0; i < sizeof(rateList) / sizeof(double) / 3; ++i, row += 3)
		{
			SECTION("Fin = " + std::to_string(row[0]) + " Fout = " + std::to_string(row[1]) + " Ffilter = " + std::to_string(row[2]) + bndMode)
			{
				checkJacobianAD(row[0], row[1], row[2], [=](cadet::JsonParameterProvider& jpp, unsigned int nComp) {
					cadet::test::addBoundStates(jpp, {1, 1}, 0.5);
					cadet::test::addLinearBindingModel(jpp, dynamic, {5.0, 4.0}, {2.0, 3.0});
				});
			}
		}
	}
}

TEST_CASE("StirredTankModel time derivative Jacobian vs FD with linear binding", "[CSTR],[UnitOp],[Residual],[Jacobian]")
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
	for (int j = 0; j < 2; ++j)
	{
		const bool dynamic = (j == 0);
		const std::string bndMode = dynamic ? " dynamic" : " quasi-stationary";
		double const* row = &rateList[0];
		for (unsigned int i = 0; i < sizeof(rateList) / sizeof(double) / 3; ++i, row += 3)
		{
			SECTION("Fin = " + std::to_string(row[0]) + " Fout = " + std::to_string(row[1]) + " Ffilter = " + std::to_string(row[2]) + bndMode)
			{
				checkTimeDerivativeJacobianFD(row[0], row[1], row[2], [=](cadet::JsonParameterProvider& jpp, unsigned int nComp) {
					cadet::test::addBoundStates(jpp, {1, 1}, 0.5);
					cadet::test::addLinearBindingModel(jpp, dynamic, {5.0, 4.0}, {2.0, 3.0});
				});
			}
		}
	}
}
