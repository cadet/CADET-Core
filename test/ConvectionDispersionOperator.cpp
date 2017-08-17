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

#include "model/operator/ConvectionDispersionOperator.hpp"
#include "Weno.hpp"

#include "JsonParameterProvider.hpp"
#include "ColumnTests.hpp"
#include "JacobianHelper.hpp"

#include <cmath>
#include <functional>

namespace
{

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
	inline void compareResidualBulkFwdBwd(double const* r1, double const* r2, int nComp, int nCol)
	{
		for (int col = 0; col < nCol; ++col)
		{
			double const* const a = r1 + nComp + nComp * col;
			double const* const b = r2 + nComp + nComp * (nCol - col - 1);
			for (int comp = 0; comp < nComp; ++comp)
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
	inline void compareResidualBulkFwdFwd(double const* r1, double const* r2, int nComp, int nCol)
	{
		for (int col = 0; col < nCol; ++col)
		{
			double const* const a = r1 + nComp + nComp * col;
			double const* const b = r2 + nComp + nComp * col;
			for (int comp = 0; comp < nComp; ++comp)
			{
				CHECK(a[comp] == Approx(b[comp]));
			}
		}	
	}

} // namespace

void testResidualBulkWenoForwardBackward(int wenoOrder)
{
	SECTION("Forward vs backward flow residual bulk (WENO=" + std::to_string(wenoOrder) + ")")
	{
		// Obtain parameters from some test case
		cadet::JsonParameterProvider jpp = createColumnWithSMA("GENERAL_RATE_MODEL");
		cadet::test::column::setWenoOrder(jpp, wenoOrder);

		const int nComp = jpp.getInt("NCOMP");
		jpp.pushScope("discretization");
		const int nCol = jpp.getInt("NCOL");
		jpp.popScope();

		cadet::model::operators::ConvectionDispersionOperator convDispOp;

		// Configure the operator
		typedef std::unordered_map<cadet::ParameterId, cadet::active*> ParameterMap;
		ParameterMap parameters;
		REQUIRE(convDispOp.configure(0, jpp, parameters, nComp, nCol));

		// Make sure that VELOCITY parameter is present
		const cadet::ParameterId paramVelocity = cadet::makeParamId(cadet::hashString("VELOCITY"), 0, cadet::CompIndep, cadet::BoundPhaseIndep, cadet::ReactionIndep, cadet::SectionIndep);
		const typename ParameterMap::iterator itVelocity = parameters.find(paramVelocity);
		REQUIRE(itVelocity != parameters.end());

		// Setup matrices
		convDispOp.notifyDiscontinuousSectionTransition(0.0, 0u, nullptr, nullptr, 0u);

		// Obtain memory for state and residual
		const int nDof = (nCol + 1) * nComp;
		std::vector<double> y(nDof, 0.0);
		std::vector<double> res(nDof, 0.0);

		// The first nComp DOFs are inlet DOFs
		for (unsigned int i = 0; i < nComp; ++i)
			y[i] = i + 1.0;

		// Check forward against backward residual with zero state
		SECTION("Forward flow yields backwards flow residual (zero state)")
		{
			// Forward flow residual
			convDispOp.residual(0.0, 0u, 1.0, y.data(), nullptr, res.data(), false);

			// Reverse flow
			itVelocity->second->setValue(-jpp.getDouble("VELOCITY"));
			convDispOp.notifyDiscontinuousSectionTransition(0.0, 0u, nullptr, nullptr, 0u);
			std::vector<double> resRev(nDof, 0.0);

			// Backward flow residual
			convDispOp.residual(0.0, 0u, 1.0, y.data(), nullptr, resRev.data(), false);

			// Compare
			compareResidualBulkFwdBwd(res.data(), resRev.data(), nComp, nCol);

			// Reverse flow again to go back to forward
			itVelocity->second->setValue(jpp.getDouble("VELOCITY"));
			convDispOp.notifyDiscontinuousSectionTransition(0.0, 0u, nullptr, nullptr, 0u);
			std::vector<double> resFwd2(nDof, 0.0);

			// Forward flow residual
			convDispOp.residual(0.0, 0u, 1.0, y.data(), nullptr, resFwd2.data(), false);

			// Compare against first forward flow residual
			compareResidualBulkFwdFwd(res.data(), resFwd2.data(), nComp, nCol);
		}

		// Check forward against backward residual with some non-zero state
		SECTION("Forward flow yields backwards flow residual (nonzero state)")
		{
			// Fill state vector with some values
			fillStateBulkFwd(y.data(), [](unsigned int comp, unsigned int col, unsigned int idx) { return std::abs(std::sin(idx * 0.13)); }, nComp, nCol);

			convDispOp.residual(0.0, 0u, 1.0, y.data(), nullptr, res.data(), false);

			// Reverse state for backwards flow
			fillStateBulkBwd(y.data(), [](unsigned int comp, unsigned int col, unsigned int idx) { return std::abs(std::sin(idx * 0.13)); }, nComp, nCol);

			// Reverse flow
			itVelocity->second->setValue(-jpp.getDouble("VELOCITY"));
			convDispOp.notifyDiscontinuousSectionTransition(0.0, 0u, nullptr, nullptr, 0u);
			std::vector<double> resRev(nDof, 0.0);

			// Backward flow residual
			convDispOp.residual(0.0, 0u, 1.0, y.data(), nullptr, resRev.data(), false);

			// Compare
			compareResidualBulkFwdBwd(res.data(), resRev.data(), nComp, nCol);

			// Reverse flow again to go back to forward
			itVelocity->second->setValue(jpp.getDouble("VELOCITY"));
			convDispOp.notifyDiscontinuousSectionTransition(0.0, 0u, nullptr, nullptr, 0u);
			std::vector<double> resFwd2(nDof, 0.0);

			// Fill state vector with forward values again
			fillStateBulkFwd(y.data(), [](unsigned int comp, unsigned int col, unsigned int idx) { return std::abs(std::sin(idx * 0.13)); }, nComp, nCol);

			// Forward flow residual
			convDispOp.residual(0.0, 0u, 1.0, y.data(), nullptr, resFwd2.data(), false);

			// Compare against first forward flow residual
			compareResidualBulkFwdFwd(res.data(), resFwd2.data(), nComp, nCol);
		}
	}
}

TEST_CASE("ConvectionDispersionOperator residual forward vs backward flow", "[Operator],[Residual]")
{
	// Test all WENO orders
	for (unsigned int i = 1; i <= cadet::Weno::maxOrder(); ++i)
		testResidualBulkWenoForwardBackward(i);
}
