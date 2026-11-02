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

#include <catch.hpp>
#include "cadet/cadet.hpp"
#include "cadet/Notification.hpp"

#define CADET_LOGGING_DISABLE
#include "Logging.hpp"

#include <vector>
#include <limits>
#include <random>
#include <chrono>
#include <algorithm>

#include "linalg/DenseMatrix.hpp"
#include "linalg/BandMatrix.hpp"
#include "AdUtils.hpp"
#include "AutoDiff.hpp"
#include "ModelBuilderImpl.hpp"
#include "model/UnitOperation.hpp"
#include "common/Driver.hpp"

#include "MatrixHelper.hpp"
#include "JacobianHelper.hpp"
#include "JsonTestModels.hpp"
#include "SimHelper.hpp"
#include "UnitOperationTests.hpp"

/**
 * @brief Creates a residual that has a banded Jacobian matrix as created by createBandMatrix()
 * @param [in] x Residual argument
 * @param [out] out Vector that holds the residual
 * @param [in] rows Number of rows
 * @param [in] lowerBand Lower bandwidth (excluding the main diagonal)
 * @param [in] upperBand Upper bandwidth (excluding the main diagonal)
 */
template <typename T>
void bandMatrixJacobian(T const* x, T* out, unsigned int rows, unsigned int lowerBand, unsigned int upperBand)
{
	double counter = 1.0;
	for (unsigned int r = 0; r < rows; ++r)
	{
		for (unsigned int c = 0; c < rows; ++c)
		{
			// Compute diagonal index of current position
			const int curDiag = static_cast<int>(c) - static_cast<int>(r);

			if ((curDiag >= -static_cast<int>(lowerBand)) && (curDiag <= static_cast<int>(upperBand)))
			{
				out[r] += counter * x[c];
				counter += 1.0;
			}
		}
	}
}

TEST_CASE("Extract banded Jacobian via AD", "[AD],[BandMatrix]")
{
	// Matrix size
	const unsigned int matSize = 10;
	const unsigned int lowerBand = 2;
	const unsigned int upperBand = 3;

	// Initialize AD and allocate AD vectors
	cadet::ad::setDirections(lowerBand + 1 + upperBand);

	cadet::active* res = new cadet::active[matSize];
	cadet::active* x = new cadet::active[matSize];

	// Set seed vectors
	cadet::ad::prepareAdVectorSeedsForBandMatrix(x, 0, matSize, lowerBand, upperBand, lowerBand);
	cadet::ad::fillAd(x, matSize, 0.0);

	// Compute residual with banded Jacobian
	bandMatrixJacobian(x, res, matSize, lowerBand, upperBand);

	// Extract the banded Jacobian matrix from AD
	cadet::linalg::BandMatrix bm;
	bm.resize(matSize, lowerBand, upperBand);

	cadet::ad::extractBandedJacobianFromAd(res, 0, lowerBand, bm);

	// Get reference matrix
	const cadet::linalg::BandMatrix ref = cadet::test::createBandMatrix<cadet::linalg::BandMatrix>(matSize, lowerBand, upperBand);

	// Compare matrices
	const unsigned int n = ref.rows() * ref.stride();
	double const* const adMat = bm.data();
	double const* const refMat = ref.data();
	for (unsigned int i = 0; i < n; ++i)
		CHECK(refMat[i] == adMat[i]);

	delete[] x;
	delete[] res;
}

TEST_CASE("Extract dense submatrix from banded Jacobian via AD", "[AD],[DenseMatrix]")
{
	// Matrix size
	const unsigned int matSize = 10;
	const unsigned int lowerBand = 2;
	const unsigned int upperBand = 3;

	// Initialize AD and allocate AD vectors
	cadet::ad::setDirections(lowerBand + 1 + upperBand);

	cadet::active* res = new cadet::active[matSize];
	cadet::active* x = new cadet::active[matSize];

	// Set seed vectors
	cadet::ad::prepareAdVectorSeedsForBandMatrix(x, 0, matSize, lowerBand, upperBand, lowerBand);
	cadet::ad::fillAd(x, matSize, 0.0);

	// Compute residual with banded Jacobian
	bandMatrixJacobian(x, res, matSize, lowerBand, upperBand);

	// Extract the banded Jacobian matrix from AD
	cadet::linalg::BandMatrix bm;
	bm.resize(matSize, lowerBand, upperBand);

	cadet::ad::extractBandedJacobianFromAd(res, 0, lowerBand, bm);

	cadet::linalg::DenseMatrix dm;
	dm.resize(3, 6);

	dm.setAll(0.0);
	cadet::ad::extractDenseJacobianFromBandedAd(res, 0, 0, lowerBand, lowerBand, upperBand, dm);
	cadet::test::checkMatrixAgainstLinearArray(dm.data(), {1, 2, 3, 4, 0, 0, 5, 6, 7, 8, 9, 0, 10, 11, 12, 13, 14, 15});

	dm.setAll(0.0);
	cadet::ad::extractDenseJacobianFromBandedAd(res, 1, 0, lowerBand, lowerBand, upperBand, dm);
	cadet::test::checkMatrixAgainstLinearArray(dm.data(), {6, 7, 8, 9, 0, 0, 11, 12, 13, 14, 15, 0, 16, 17, 18, 19, 20, 21});

	dm.setAll(0.0);
	cadet::ad::extractDenseJacobianFromBandedAd(res, 2, 0, lowerBand, lowerBand, upperBand, dm);
	cadet::test::checkMatrixAgainstLinearArray(dm.data(), {12, 13, 14, 15, 0, 0, 17, 18, 19, 20, 21, 0, 22, 23, 24, 25, 26, 27});

	dm.setAll(0.0);
	cadet::ad::extractDenseJacobianFromBandedAd(res, 5, 0, lowerBand, lowerBand, upperBand, dm);
	cadet::test::checkMatrixAgainstLinearArray(dm.data(), {30, 31, 32, 33, 0, 0, 35, 36, 37, 38, 39, 0, 40, 41, 42, 43, 44, 0});

	dm.setAll(0.0);
	cadet::ad::extractDenseJacobianFromBandedAd(res, 6, 0, lowerBand, lowerBand, upperBand, dm);
	cadet::test::checkMatrixAgainstLinearArray(dm.data(), {36, 37, 38, 39, 0, 0, 41, 42, 43, 44, 0, 0, 45, 46, 47, 48, 0, 0});

	dm.setAll(0.0);
	cadet::ad::extractDenseJacobianFromBandedAd(res, 7, 0, lowerBand, lowerBand, upperBand, dm);
	cadet::test::checkMatrixAgainstLinearArray(dm.data(), {42, 43, 44, 0, 0, 0, 46, 47, 48, 0, 0, 0, 49, 50, 51, 0, 0, 0});

	delete[] x;
	delete[] res;
}

TEST_CASE("Extract square dense submatrix from banded Jacobian via AD", "[AD],[DenseMatrix]")
{
	// Matrix size
	const unsigned int matSize = 32;
	const unsigned int lowerBand = 8;
	const unsigned int upperBand = 12;

	// Initialize AD and allocate AD vectors
	cadet::ad::setDirections(lowerBand + 1 + upperBand);

	cadet::active* res = new cadet::active[matSize];
	cadet::active* x = new cadet::active[matSize];

	// Set seed vectors
	cadet::ad::prepareAdVectorSeedsForBandMatrix(x, 0, matSize, lowerBand, upperBand, lowerBand);
	cadet::ad::fillAd(x, matSize, 0.0);

	// Compute residual with banded Jacobian
	bandMatrixJacobian(x, res, matSize, lowerBand, upperBand);

	cadet::linalg::DenseMatrix dm;
	dm.resize(8, 8);

	dm.setAll(0.0);
	cadet::ad::extractDenseJacobianFromBandedAd(res, 0, 0, lowerBand, lowerBand, upperBand, dm);
	cadet::test::checkMatrixAgainstLinearArray(dm.data(), {1, 2, 3, 4, 5, 6, 7, 8,
	                                                       14, 15, 16, 17, 18, 19, 20, 21,
	                                                       28, 29, 30, 31, 32, 33, 34, 35,
	                                                       43, 44, 45, 46, 47, 48, 49, 50,
	                                                       59, 60, 61, 62, 63, 64, 65, 66,
	                                                       76, 77, 78, 79, 80, 81, 82, 83,
	                                                       94, 95, 96, 97, 98, 99, 100, 101,
	                                                       113, 114, 115, 116, 117, 118, 119, 120});

	dm.setAll(0.0);
	cadet::ad::extractDenseJacobianFromBandedAd(res, 8, 0, lowerBand, lowerBand, upperBand, dm);
	cadet::test::checkMatrixAgainstLinearArray(dm.data(), {141, 142, 143, 144, 145, 146, 147, 148,
	                                                       161, 162, 163, 164, 165, 166, 167, 168,
	                                                       181, 182, 183, 184, 185, 186, 187, 188,
	                                                       201, 202, 203, 204, 205, 206, 207, 208,
	                                                       221, 222, 223, 224, 225, 226, 227, 228,
	                                                       241, 242, 243, 244, 245, 246, 247, 248,
	                                                       261, 262, 263, 264, 265, 266, 267, 268,
	                                                       281, 282, 283, 284, 285, 286, 287, 288});

	delete[] x;
	delete[] res;
}

TEST_CASE("Banded AD Jacobian vs FD", "[AD],[BandMatrix]")
{
	// Matrix size
	const unsigned int matSize = 32;
	const unsigned int lowerBand = 8;
	const unsigned int upperBand = 12;

	// Initialize AD and allocate AD vectors
	cadet::ad::setDirections(lowerBand + 1 + upperBand);

	cadet::active* res = new cadet::active[matSize];
	cadet::active* x = new cadet::active[matSize];

	// Set seed vectors
	cadet::ad::prepareAdVectorSeedsForBandMatrix(x, 0, matSize, lowerBand, upperBand, lowerBand);
	cadet::ad::fillAd(x, matSize, 0.0);

	// Compute residual with banded Jacobian
	bandMatrixJacobian(x, res, matSize, lowerBand, upperBand);

	// Extract the banded Jacobian matrix from AD
	cadet::linalg::BandMatrix bm;
	bm.resize(matSize, lowerBand, upperBand);

	cadet::ad::extractBandedJacobianFromAd(res, 0, lowerBand, bm);

	delete[] x;
	delete[] res;

	std::vector<double> y(matSize, 0.0);
	std::vector<double> dir(matSize, 0.0);
	std::vector<double> colA(matSize, 0.0);
	std::vector<double> colB(matSize, 0.0);

	cadet::test::compareJacobianFD(
		[=](double const* y, double* r) { std::fill_n(r, matSize, 0.0); bandMatrixJacobian(y, r, matSize, lowerBand, upperBand); },
		[&](double const* y, double* r) { bm.multiplyVector(y, r); },
		y.data(), dir.data(), colA.data(), colB.data(), matSize, matSize, 1e-7, 0.0, 1e-15
	);
}

TEST_CASE("Models are built with the maximum number of AD directions", "[AD],[UnitOp],[CI]")
{
	// An active zeroes its gradient only over the global direction count, so a model must not be built
	// while that count is lowered: its parameters would carry uninitialized gradients into the next
	// simulation, which then seeds and reads more directions than were ever initialized. Time
	// integration restores the count when it returns, but the stepwise interface leaves it lowered
	// between steps, so building a model raises it in its own right.
	cadet::ad::setDirections(1);

	cadet::IModelBuilder* const mb = cadet::createModelBuilder();
	REQUIRE(nullptr != mb);

	cadet::JsonParameterProvider jpp = createCSTR(2);
	cadet::IUnitOperation* const unit = cadet::test::unitoperation::createAndConfigureUnit(jpp, *mb);

	CHECK(cadet::ad::getDirections() == cadet::ad::getMaxDirections());

	mb->destroyUnitOperation(unit);
	destroyModelBuilder(mb);
}

namespace
{
	/**
	 * @brief Records the number of AD directions that is in effect while a simulation runs
	 * @details The direction count is lowered before the integration starts and restored when it returns,
	 *          so it cannot be observed from the outside. This callback reads it from within the run.
	 */
	class AdDirectionRecorder : public cadet::INotificationCallback
	{
	public:
		AdDirectionRecorder() : _directions(-1) { }

		void timeIntegrationStart() override { _directions = static_cast<int>(cadet::ad::getDirections()); }
		void timeIntegrationEnd() override { }
		void timeIntegrationError(char const* message, unsigned int section, double time, double progress) override { }
		bool timeIntegrationSection(unsigned int section, double time, double const* state, double const* stateDot, double progress) override { return true; }
		bool timeIntegrationStep(unsigned int section, double time, double const* state, double const* stateDot, double progress) override { return true; }
		bool timeIntegrationLinearSolve(unsigned int section, double time, double const* state, double const* stateDot) override { return true; }

		/**
		 * @brief Returns the direction count seen during the run, or @c -1 if the run never started
		 */
		int directions() const { return _directions; }

	private:
		int _directions;
	};
}

TEST_CASE("Time integration lowers the global direction count to what its model needs and restores it", "[AD],[Simulation],[CI]")
{
	// A simulation lowers the global direction count to one direction per sensitivity plus the
	// requirement of its model, so that the actives of its residual evaluations zero no more of their
	// gradient than is read back. The lowered count is only valid while the integration runs: an active
	// created afterwards would zero too few directions and hand uninitialized memory to the next run.
	// So integrate() puts the count back to the maximum on return.

	// Requirement of the model, determined without running a simulation
	const auto modelRequirement = [](bool analyticJac)
	{
		cadet::ad::setDirections(cadet::ad::getMaxDirections());

		cadet::IModelBuilder* const mb = cadet::createModelBuilder();
		REQUIRE(nullptr != mb);

		cadet::JsonParameterProvider unitJpp = createCSTR(1);
		unitJpp.set("USE_ANALYTIC_JACOBIAN", analyticJac);

		cadet::IUnitOperation* const unit = cadet::test::unitoperation::createAndConfigureUnit(unitJpp, *mb);
		const unsigned int dirs = unit->requiredADdirs();

		mb->destroyUnitOperation(unit);
		destroyModelBuilder(mb);

		return dirs;
	};

	const auto useAdJacobian = [](cadet::JsonParameterProvider& p)
	{
		p.pushScope("model");
		p.pushScope("unit_000");
		p.set("USE_ANALYTIC_JACOBIAN", false);
		p.popScope();
		p.popScope();
	};

	// Unit 000 of this configuration is the model that modelRequirement() builds; the inlet and the
	// outlet unit require no AD directions, so the system requires what unit 000 requires
	cadet::JsonParameterProvider jpp = createCSTRBenchmark(1, 10.0, 1.0);
	cadet::test::setSectionTimes(jpp, {0.0, 10.0});
	cadet::test::setInitialConditions(jpp, {0.0}, {}, 10.0);
	cadet::test::setInletProfile(jpp, 0, 0, 1.0, 0.0, 0.0, 0.0);
	cadet::test::setFlowRates(jpp, 0, 1.0, 0.5, 0.5);

	SECTION("Analytic Jacobian without sensitivities")
	{
		const unsigned int required = modelRequirement(true);

		cadet::ad::setDirections(cadet::ad::getMaxDirections());

		AdDirectionRecorder rec;
		cadet::Driver drv;
		drv.configure(jpp);
		drv.simulator()->setNotificationCallback(&rec);
		drv.run();

		CHECK(rec.directions() == static_cast<int>(required));
		CHECK(cadet::ad::getDirections() == cadet::ad::getMaxDirections());
	}

	SECTION("AD Jacobian without sensitivities")
	{
		const unsigned int required = modelRequirement(false);
		REQUIRE(required > 0);

		useAdJacobian(jpp);

		cadet::ad::setDirections(cadet::ad::getMaxDirections());

		AdDirectionRecorder rec;
		cadet::Driver drv;
		drv.configure(jpp);
		drv.simulator()->setNotificationCallback(&rec);
		drv.run();

		CHECK(rec.directions() == static_cast<int>(required));
		CHECK(cadet::ad::getDirections() == cadet::ad::getMaxDirections());
	}

	SECTION("AD Jacobian with a parameter sensitivity")
	{
		const unsigned int required = modelRequirement(false);

		useAdJacobian(jpp);
		cadet::test::addSensitivity(jpp, "FLOWRATE_FILTER", cadet::makeParamId("FLOWRATE_FILTER", 0, cadet::CompIndep, cadet::ParTypeIndep, cadet::BoundStateIndep, cadet::ReactionIndep, cadet::SectionIndep), 1e-6);
		cadet::test::returnSensitivities(jpp, 0);

		cadet::ad::setDirections(cadet::ad::getMaxDirections());

		AdDirectionRecorder rec;
		cadet::Driver drv;
		drv.configure(jpp);
		drv.simulator()->setNotificationCallback(&rec);
		drv.run();

		REQUIRE(drv.simulator()->numSensParams() == 1);
		CHECK(rec.directions() == static_cast<int>(1 + required));
		CHECK(cadet::ad::getDirections() == cadet::ad::getMaxDirections());
	}
}
