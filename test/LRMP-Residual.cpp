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

#include "model/LumpedRateModelWithPores.hpp"
#include "Weno.hpp"
#include "ModelBuilderImpl.hpp"

#include "JsonParameterProvider.hpp"
#include "JacobianHelper.hpp"

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
		jpp.pushScope("discretization");
		jpp.pushScope("weno");

		jpp.set("WENO_ORDER", order);

		jpp.popScope();
		jpp.popScope();
	}
}

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
 * @brief Fills the bulk part of the state vector with a given function
 * @details The function @p f uses the component index, the column cell index, and the 
 *          current (forward flow) index to assign a value.
 * @param [out] y Filled state vector
 * @param [in] f Function for computing the content of the bulk part of the state vector
 * @param [in] nComp Number of components
 * @param [in] nCol Number of bulk cells per component
 */
inline void fillStateBulkFwd(double* y, std::function<double(unsigned int, unsigned int, unsigned int)> f, unsigned int nComp, unsigned int nCol)
{
	for (unsigned int col = 0; col < nCol; ++col)
	{
		double* const r = y + nComp + nComp * col;
		for (unsigned int comp = 0; comp < nComp; ++comp)
		{
			r[comp] = f(comp, col, col * nComp + comp);
		}
	}
}

/**
 * @brief Fills the bulk part of the state vector with a given function assuming backwards flow
 * @details The function @p f uses the component index, the column cell index, and the 
 *          current (forward flow) index to assign a value.
 * @param [out] y Filled state vector
 * @param [in] f Function for computing the content of the bulk part of the state vector
 * @param [in] nComp Number of components
 * @param [in] nCol Number of bulk cells per component
 */
inline void fillStateBulkBwd(double* y, std::function<double(unsigned int, unsigned int, unsigned int)> f, unsigned int nComp, unsigned int nCol)
{
	for (unsigned int col = nCol-1; col < nCol; --col)
	{
		double* const r = y + nComp + nComp * col;
		for (unsigned int comp = 0; comp < nComp; ++comp)
		{
			r[comp] = f(comp, nCol - col - 1, (nCol - col - 1) * nComp + comp);
		}
	}
}

/**
 * @brief Compares two residual bulk parts where one is created by forward and the other by backwards flow
 * @param [in] r1 Residual of one mode (forward / backward)
 * @param [in] r2 Residual of the other mode (backward / forward)
 * @param [in] nComp Number of components
 * @param [in] nCol Number of bulk cells per component
 */
inline void compareResidualBulkFwdBwd(double const* r1, double const* r2, unsigned int nComp, unsigned int nCol)
{
	for (unsigned int col = 0; col < nCol; ++col)
	{
		double const* const a = r1 + nComp + nComp * col;
		double const* const b = r2 + nComp + nComp * (nCol - col - 1);
		for (unsigned int comp = 0; comp < nComp; ++comp)
		{
			CHECK(a[comp] == Approx(b[comp]));
		}
	}	
}

/**
 * @brief Compares two residual bulk parts that have the same flow direction
 * @param [in] r1 Residual
 * @param [in] r2 Residual
 * @param [in] nComp Number of components
 * @param [in] nCol Number of bulk cells per component
 */
inline void compareResidualBulkFwdFwd(double const* r1, double const* r2, unsigned int nComp, unsigned int nCol)
{
	for (unsigned int col = 0; col < nCol; ++col)
	{
		double const* const a = r1 + nComp + nComp * col;
		double const* const b = r2 + nComp + nComp * col;
		for (unsigned int comp = 0; comp < nComp; ++comp)
		{
			CHECK(a[comp] == Approx(b[comp]));
		}
	}	
}

/**
 * @brief Creates a runnable LRMP model with given WENO order
 * @details Creates a LRMP model and configures it using the given IParameterProvider @p jpp.
 * @param [in] mb ModelBuilder
 * @param [in] jpp Configuration of the LRMP
 * @param [in] wenoOrder WENO order
 * @return Runnable LRMP
 */
inline cadet::model::LumpedRateModelWithPores* createAndConfigureLRMP(cadet::IModelBuilder& mb, cadet::JsonParameterProvider& jpp, int wenoOrder)
{
	// Create a LRMP
	cadet::IModel* const iLrmp = mb.createUnitOperation("LUMPED_RATE_MODEL_WITH_PORES", 0);
	REQUIRE(nullptr != iLrmp);

	cadet::model::LumpedRateModelWithPores* const lrmp = reinterpret_cast<cadet::model::LumpedRateModelWithPores*>(iLrmp);

	// Set WENO order
	setWenoOrder(jpp, wenoOrder);
	// Configure
	cadet::ModelBuilder& temp = *reinterpret_cast<cadet::ModelBuilder*>(&mb);
	REQUIRE(lrmp->configure(jpp, temp));

	// Do some checks
	const unsigned int nComp = jpp.getInt("NCOMP");
	REQUIRE(lrmp->numComponents() == nComp);

	return lrmp;
}

inline void testJacobianWenoForwardBackward(int wenoOrder)
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	SECTION("Forward vs backward flow Jacobian (WENO=" + std::to_string(wenoOrder) + ")")
	{
		// Use LRMP with SMA binding as test case
		cadet::JsonParameterProvider jpp = createGRMwithLinear();
		const unsigned int nComp = jpp.getInt("NCOMP");

		cadet::model::LumpedRateModelWithPores* const lrmpAna = createAndConfigureLRMP(*mb, jpp, wenoOrder);
		cadet::model::LumpedRateModelWithPores* const lrmpAD = createAndConfigureLRMP(*mb, jpp, wenoOrder);

		// Enable AD
		cadet::ad::setDirections(cadet::ad::getMaxDirections());
		lrmpAD->useAnalyticJacobian(false);

		cadet::active* adRes = new cadet::active[lrmpAD->numDofs()];
		cadet::active* adY = new cadet::active[lrmpAD->numDofs()];

		lrmpAD->prepareADvectors(adRes, adY, 0);

		// Setup matrices
		lrmpAna->notifyDiscontinuousSectionTransition(0.0, 0u, nullptr, nullptr, 0u);
		lrmpAD->notifyDiscontinuousSectionTransition(0.0, 0u, adRes, adY, 0u);

		// Obtain memory for state, Jacobian multiply direction, Jacobian column
		std::vector<double> y(lrmpAD->numDofs(), 0.0);
		std::vector<double> jacDir(lrmpAD->numDofs(), 0.0);
		std::vector<double> jacCol1(lrmpAD->numDofs(), 0.0);
		std::vector<double> jacCol2(lrmpAD->numDofs(), 0.0);

		// Fill state vector with some values
		fillState(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, lrmpAna->numDofs());
//		fillState(y.data(), [](unsigned int idx) { return 1.0; }, lrmpAna->numDofs());

		// Obtain number of column cells
		jpp.pushScope("discretization");
		const unsigned int nCol = jpp.getInt("NCOL");
		REQUIRE(nCol == 15u);
		jpp.popScope();

		SECTION("Forward then backward flow (nonzero state)")
		{
			// Compute state Jacobian
			lrmpAna->residualWithJacobian(0.0, 0u, 1.0, y.data(), nullptr, jacDir.data(), nullptr, nullptr, 0u);
			lrmpAD->residualWithJacobian(0.0, 0u, 1.0, y.data(), nullptr, jacDir.data(), adRes, adY, 0u);
			std::fill(jacDir.begin(), jacDir.end(), 0.0);

			// Compare Jacobians
			cadet::test::checkJacobianPatternFD(lrmpAna, lrmpAD, y.data(), nullptr, jacDir.data(), jacCol1.data(), jacCol2.data());
			cadet::test::checkJacobianPatternFD(lrmpAna, lrmpAna, y.data(), nullptr, jacDir.data(), jacCol1.data(), jacCol2.data());
			cadet::test::compareJacobian(lrmpAna, lrmpAD, nullptr, nullptr, jacDir.data(), jacCol1.data(), jacCol2.data());
//			cadet::test::compareJacobianFD(lrmpAna, lrmpAD, y.data(), nullptr, jacDir.data(), jacCol1.data(), jacCol2.data());

			// Reverse flow
			const bool paramSet = lrmpAna->setParameter(cadet::makeParamId(cadet::hashString("VELOCITY"), 0, cadet::CompIndep, cadet::BoundPhaseIndep, cadet::ReactionIndep, cadet::SectionIndep), -jpp.getDouble("VELOCITY"));
			REQUIRE(paramSet);
			// Reverse flow
			const bool paramSet2 = lrmpAD->setParameter(cadet::makeParamId(cadet::hashString("VELOCITY"), 0, cadet::CompIndep, cadet::BoundPhaseIndep, cadet::ReactionIndep, cadet::SectionIndep), -jpp.getDouble("VELOCITY"));
			REQUIRE(paramSet2);

			// Setup
			lrmpAna->notifyDiscontinuousSectionTransition(0.0, 0u, nullptr, nullptr, 0u);
			lrmpAD->notifyDiscontinuousSectionTransition(0.0, 0u, adRes, adY, 0u);

			// Compute state Jacobian
			lrmpAna->residualWithJacobian(0.0, 0u, 1.0, y.data(), nullptr, jacDir.data(), nullptr, nullptr, 0u);
			lrmpAD->residualWithJacobian(0.0, 0u, 1.0, y.data(), nullptr, jacDir.data(), adRes, adY, 0u);
			std::fill(jacDir.begin(), jacDir.end(), 0.0);

			// Compare Jacobians
			cadet::test::checkJacobianPatternFD(lrmpAna, lrmpAD, y.data(), nullptr, jacDir.data(), jacCol1.data(), jacCol2.data());
			cadet::test::checkJacobianPatternFD(lrmpAna, lrmpAna, y.data(), nullptr, jacDir.data(), jacCol1.data(), jacCol2.data());
//			cadet::test::compareJacobianFD(lrmpAD, lrmpAna, y.data(), jacDir.data(), nullptr, jacCol1.data(), jacCol2.data());
//			cadet::test::compareJacobianFD(lrmpAna, lrmpAD, y.data(), jacDir.data(), nullptr, jacCol1.data(), jacCol2.data());
			cadet::test::compareJacobian(lrmpAna, lrmpAD, nullptr, nullptr, jacDir.data(), jacCol1.data(), jacCol2.data());
		}

		delete[] adRes;
		delete[] adY;
		mb->destroyUnitOperation(lrmpAna);
		mb->destroyUnitOperation(lrmpAD);
	}
	destroyModelBuilder(mb);
}

TEST_CASE("LumpedRateModelWithPores Jacobian forward vs backward flow", "[LRMP],[UnitOp],[Residual],[Jacobian],[AD]")
{
	// Test all WENO orders
	for (unsigned int i = 1; i <= cadet::Weno::maxOrder(); ++i)
		testJacobianWenoForwardBackward(i);
}
