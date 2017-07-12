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

#include "model/GeneralRateModel.hpp"
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
	for (unsigned int comp = 0; comp < nComp; ++comp)
	{
		double* const r = y + nComp + nCol * comp;
		for (unsigned int i = 0; i < nCol; ++i)
		{
			r[i] = f(comp, i, i + comp * nCol);
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
	for (unsigned int comp = 0; comp < nComp; ++comp)
	{
		double* const r = y + nComp + nCol * comp;
		for (unsigned int i = 0; i < nCol; ++i)
		{
			r[nCol - i - 1] = f(comp, i, i + comp * nCol);
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
	for (unsigned int comp = 0; comp < nComp; ++comp)
	{
		double const* const a = r1 + nComp + nCol * comp;
		double const* const b = r2 + nComp + nCol * comp;
		for (unsigned int i = 0; i < nCol; ++i)
		{
			CHECK(a[i] == Approx(b[nCol - i - 1]));
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
	for (unsigned int comp = 0; comp < nComp; ++comp)
	{
		double const* const a = r1 + nComp + nCol * comp;
		double const* const b = r2 + nComp + nCol * comp;
		for (unsigned int i = 0; i < nCol; ++i)
		{
			CHECK(a[i] == Approx(b[i]));
		}
	}	
}

/**
 * @brief Creates a runnable GRM model with given WENO order
 * @details Creates a GRM model and configures it using the given IParameterProvider @p jpp.
 * @param [in] mb ModelBuilder
 * @param [in] jpp Configuration of the GRM
 * @param [in] wenoOrder WENO order
 * @return Runnable GRM
 */
cadet::model::GeneralRateModel* createAndConfigureGRM(cadet::IModelBuilder& mb, cadet::JsonParameterProvider& jpp, int wenoOrder)
{
	// Create a GRM
	cadet::IModel* const iGrm = mb.createUnitOperation("GENERAL_RATE_MODEL", 0);
	REQUIRE(nullptr != iGrm);

	cadet::model::GeneralRateModel* const grm = reinterpret_cast<cadet::model::GeneralRateModel*>(iGrm);

	// Set WENO order
	setWenoOrder(jpp, wenoOrder);
	// Configure
	cadet::ModelBuilder& temp = *reinterpret_cast<cadet::ModelBuilder*>(&mb);
	REQUIRE(grm->configure(jpp, temp));

	// Do some checks
	const unsigned int nComp = jpp.getInt("NCOMP");
	REQUIRE(grm->numComponents() == nComp);

	return grm;
}

void testJacobianWenoForwardBackward(int wenoOrder)
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	SECTION("Forward vs backward flow Jacobian (WENO=" + std::to_string(wenoOrder) + ")")
	{
		// Use GRM with SMA binding as test case
		cadet::JsonParameterProvider jpp = createGRMwithLinear();
		const unsigned int nComp = jpp.getInt("NCOMP");

		cadet::model::GeneralRateModel* const grmAna = createAndConfigureGRM(*mb, jpp, wenoOrder);
		cadet::model::GeneralRateModel* const grmAD = createAndConfigureGRM(*mb, jpp, wenoOrder);

		// Enable AD
		cadet::ad::setDirections(cadet::ad::getMaxDirections());
		grmAD->useAnalyticJacobian(false);

		cadet::active* adRes = new cadet::active[grmAD->numDofs()];
		cadet::active* adY = new cadet::active[grmAD->numDofs()];

		grmAD->prepareADvectors(adRes, adY, 0);

		// Setup matrices
		grmAna->notifyDiscontinuousSectionTransition(0.0, 0u, nullptr, nullptr, 0u);
		grmAD->notifyDiscontinuousSectionTransition(0.0, 0u, adRes, adY, 0u);

		// Obtain memory for state, Jacobian multiply direction, Jacobian column
		std::vector<double> y(grmAD->numDofs(), 0.0);
		std::vector<double> jacDir(grmAD->numDofs(), 0.0);
		std::vector<double> jacCol1(grmAD->numDofs(), 0.0);
		std::vector<double> jacCol2(grmAD->numDofs(), 0.0);

		// Fill state vector with some values
		fillState(y.data(), [](unsigned int idx) { return std::abs(std::sin(idx * 0.13)) + 1e-4; }, grmAna->numDofs());
//		fillState(y.data(), [](unsigned int idx) { return 1.0; }, grmAna->numDofs());

		// Obtain number of column cells
		jpp.pushScope("discretization");
		const unsigned int nCol = jpp.getInt("NCOL");
		REQUIRE(nCol == 15u);
		jpp.popScope();

		SECTION("Forward then backward flow (nonzero state)")
		{
			// Compute state Jacobian
			grmAna->residualWithJacobian(0.0, 0u, 1.0, y.data(), nullptr, jacDir.data(), nullptr, nullptr, 0u);
			grmAD->residualWithJacobian(0.0, 0u, 1.0, y.data(), nullptr, jacDir.data(), adRes, adY, 0u);
			std::fill(jacDir.begin(), jacDir.end(), 0.0);

			// Compare Jacobians
			cadet::test::checkJacobianPatternFD(grmAna, grmAD, y.data(), nullptr, jacDir.data(), jacCol1.data(), jacCol2.data());
			cadet::test::checkJacobianPatternFD(grmAna, grmAna, y.data(), nullptr, jacDir.data(), jacCol1.data(), jacCol2.data());
			cadet::test::compareJacobian(grmAna, grmAD, nullptr, nullptr, jacDir.data(), jacCol1.data(), jacCol2.data());
//			cadet::test::compareJacobianFD(grmAna, grmAD, y.data(), nullptr, jacDir.data(), jacCol1.data(), jacCol2.data());

			// Reverse flow
			const bool paramSet = grmAna->setParameter(cadet::makeParamId(cadet::hashString("VELOCITY"), 0, cadet::CompIndep, cadet::BoundPhaseIndep, cadet::ReactionIndep, cadet::SectionIndep), -jpp.getDouble("VELOCITY"));
			REQUIRE(paramSet);
			// Reverse flow
			const bool paramSet2 = grmAD->setParameter(cadet::makeParamId(cadet::hashString("VELOCITY"), 0, cadet::CompIndep, cadet::BoundPhaseIndep, cadet::ReactionIndep, cadet::SectionIndep), -jpp.getDouble("VELOCITY"));
			REQUIRE(paramSet2);

			// Setup
			grmAna->notifyDiscontinuousSectionTransition(0.0, 0u, nullptr, nullptr, 0u);
			grmAD->notifyDiscontinuousSectionTransition(0.0, 0u, adRes, adY, 0u);

			// Compute state Jacobian
			grmAna->residualWithJacobian(0.0, 0u, 1.0, y.data(), nullptr, jacDir.data(), nullptr, nullptr, 0u);
			grmAD->residualWithJacobian(0.0, 0u, 1.0, y.data(), nullptr, jacDir.data(), adRes, adY, 0u);
			std::fill(jacDir.begin(), jacDir.end(), 0.0);

			// Compare Jacobians
			cadet::test::checkJacobianPatternFD(grmAna, grmAD, y.data(), nullptr, jacDir.data(), jacCol1.data(), jacCol2.data());
			cadet::test::checkJacobianPatternFD(grmAna, grmAna, y.data(), nullptr, jacDir.data(), jacCol1.data(), jacCol2.data());
//			cadet::test::compareJacobianFD(grmAD, grmAna, y.data(), jacDir.data(), nullptr, jacCol1.data(), jacCol2.data());
//			cadet::test::compareJacobianFD(grmAna, grmAD, y.data(), jacDir.data(), nullptr, jacCol1.data(), jacCol2.data());
			cadet::test::compareJacobian(grmAna, grmAD, nullptr, nullptr, jacDir.data(), jacCol1.data(), jacCol2.data());
		}
		mb->destroyUnitOperation(grmAna);
		mb->destroyUnitOperation(grmAD);
	}
	destroyModelBuilder(mb);
}

void testResidualBulkWenoForwardBackward(int wenoOrder)
{
	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	SECTION("Forward vs backward flow residual bulk (WENO=" + std::to_string(wenoOrder) + ")")
	{
		// Use GRM with SMA binding as test case
		cadet::JsonParameterProvider jpp = createGRMwithSMA();
		const unsigned int nComp = jpp.getInt("NCOMP");

		cadet::model::GeneralRateModel* const grm = createAndConfigureGRM(*mb, jpp, wenoOrder);

		// Setup matrices
		grm->notifyDiscontinuousSectionTransition(0.0, 0u, nullptr, nullptr, 0u);

		// Obtain memory for state and residual
		std::vector<double> y(grm->numDofs(), 0.0);
		std::vector<double> res(grm->numDofs(), 0.0);

		// The first nComp DOFs are inlet DOFs
		for (unsigned int i = 0; i < nComp; ++i)
			y[i] = i + 1.0;

		// Obtain number of column cells
		jpp.pushScope("discretization");
		const unsigned int nCol = jpp.getInt("NCOL");
		REQUIRE(nCol == 16u);
		jpp.popScope();

		// Check forward against backward residual with zero state
		SECTION("Forward flow yields backwards flow residual (zero state)")
		{
			// Forward flow residual
			grm->residual(0.0, 0u, 1.0, y.data(), nullptr, res.data());

			// Reverse flow
			const bool paramSet = grm->setParameter(cadet::makeParamId(cadet::hashString("VELOCITY"), 0, cadet::CompIndep, cadet::BoundPhaseIndep, cadet::ReactionIndep, cadet::SectionIndep), -jpp.getDouble("VELOCITY"));
			REQUIRE(paramSet);

			// Setup
			grm->notifyDiscontinuousSectionTransition(0.0, 0u, nullptr, nullptr, 0u);
			std::vector<double> resRev(grm->numDofs(), 0.0);

			// Backward flow residual
			grm->residual(0.0, 0u, 1.0, y.data(), nullptr, resRev.data());

			// Compare
			compareResidualBulkFwdBwd(res.data(), resRev.data(), nComp, nCol);

			// Reverse flow again to go back to forward
			const bool paramSet2 = grm->setParameter(cadet::makeParamId(cadet::hashString("VELOCITY"), 0, cadet::CompIndep, cadet::BoundPhaseIndep, cadet::ReactionIndep, cadet::SectionIndep), jpp.getDouble("VELOCITY"));
			REQUIRE(paramSet2);

			// Setup again to update matrices
			grm->notifyDiscontinuousSectionTransition(0.0, 0u, nullptr, nullptr, 0u);
			std::vector<double> resFwd2(grm->numDofs(), 0.0);

			// Forward flow residual
			grm->residual(0.0, 0u, 1.0, y.data(), nullptr, resFwd2.data());

			// Compare against first forward flow residual
			compareResidualBulkFwdFwd(res.data(), resFwd2.data(), nComp, nCol);
		}

		// Check forward against backward residual with some non-zero state
		SECTION("Forward flow yields backwards flow residual (nonzero state)")
		{
			// Fill state vector with some values
			fillStateBulkFwd(y.data(), [](unsigned int comp, unsigned int col, unsigned int idx) { return std::abs(std::sin(idx * 0.13)); }, nComp, nCol);

			grm->residual(0.0, 0u, 1.0, y.data(), nullptr, res.data());

			// Reverse state for backwards flow
			fillStateBulkBwd(y.data(), [](unsigned int comp, unsigned int col, unsigned int idx) { return std::abs(std::sin(idx * 0.13)); }, nComp, nCol);

			// Reverse flow
			const bool paramSet = grm->setParameter(cadet::makeParamId(cadet::hashString("VELOCITY"), 0, cadet::CompIndep, cadet::BoundPhaseIndep, cadet::ReactionIndep, cadet::SectionIndep), -jpp.getDouble("VELOCITY"));
			REQUIRE(paramSet);

			// Setup
			grm->notifyDiscontinuousSectionTransition(0.0, 0u, nullptr, nullptr, 0u);
			std::vector<double> resRev(grm->numDofs(), 0.0);

			// Backward flow residual
			grm->residual(0.0, 0u, 1.0, y.data(), nullptr, resRev.data());

			// Compare
			compareResidualBulkFwdBwd(res.data(), resRev.data(), nComp, nCol);

			// Reverse flow again to go back to forward
			const bool paramSet2 = grm->setParameter(cadet::makeParamId(cadet::hashString("VELOCITY"), 0, cadet::CompIndep, cadet::BoundPhaseIndep, cadet::ReactionIndep, cadet::SectionIndep), jpp.getDouble("VELOCITY"));
			REQUIRE(paramSet2);

			// Setup again to update matrices
			grm->notifyDiscontinuousSectionTransition(0.0, 0u, nullptr, nullptr, 0u);
			std::vector<double> resFwd2(grm->numDofs(), 0.0);

			// Fill state vector with forward values again
			fillStateBulkFwd(y.data(), [](unsigned int comp, unsigned int col, unsigned int idx) { return std::abs(std::sin(idx * 0.13)); }, nComp, nCol);

			// Forward flow residual
			grm->residual(0.0, 0u, 1.0, y.data(), nullptr, resFwd2.data());

			// Compare against first forward flow residual
			compareResidualBulkFwdFwd(res.data(), resFwd2.data(), nComp, nCol);
		}
		mb->destroyUnitOperation(grm);
	}
	destroyModelBuilder(mb);
}

TEST_CASE("GeneralRateModel bulk residual forward vs backward flow", "[GRM],[UnitOp],[Residual]")
{
	// Test all WENO orders
	for (unsigned int i = 1; i < cadet::Weno::maxOrder(); ++i)
		testResidualBulkWenoForwardBackward(i);
}

TEST_CASE("GeneralRateModel Jacobian forward vs backward flow", "[GRM],[UnitOp],[Residual],[Jacobian]")
{
	// Test all WENO orders
	for (unsigned int i = 1; i < cadet::Weno::maxOrder(); ++i)
		testJacobianWenoForwardBackward(i);
}
