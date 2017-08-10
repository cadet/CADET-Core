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

#include <vector>
#include <limits>
#include <random>
#include <chrono>

#include "linalg/DenseMatrix.hpp"
#include "linalg/BandMatrix.hpp"
#include "AdUtils.hpp"
#include "AutoDiff.hpp"

#include "MatrixHelper.hpp"

/**
 * @brief Create a BandMatrix of given type and fill entries with their linear array index (1-based)
 * @param [in] Number of rows
 * @param [in] Lower bandwidth (excluding main diagonal)
 * @param [in] Upper bandwidth (excluding main diagonal)
 * @tparam Matrix_t Type of banded matrix to create
 * @return Banded matrix of given shape
 */
template <typename Matrix_t>
Matrix_t createBandMatrix(unsigned int rows, unsigned int lower, unsigned int upper)
{
	Matrix_t bm;
	bm.resize(rows, lower, upper);

	double val = 1.0;
	for (unsigned int row = 0; row < bm.rows(); ++row)
	{
		const int lower = std::max(-static_cast<int>(bm.lowerBandwidth()), -static_cast<int>(row));
		const int upper = std::min(static_cast<int>(bm.upperBandwidth()), static_cast<int>(bm.rows() - row) - 1);
		for (int col = lower; col <= upper; ++col)
		{
			bm.centered(row, col) = val;
			val += 1.0;
		}
	}
	return bm;
}

/**
 * @brief Creates a residual that has a banded Jacobian matrix as created by createBandMatrix()
 * @param [in] x AD enabled residual argument
 * @param [out] out Vector that holds the residual
 * @param [in] rows Number of rows
 * @param [in] lowerBand Lower bandwidth (excluding the main diagonal)
 * @param [in] upperBand Upper bandwidth (excluding the main diagonal)
 */
void bandMatrixJacobian(cadet::active const* x, cadet::active* out, unsigned int rows, unsigned int lowerBand, unsigned int upperBand)
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

	// Compute residual with banded Jacobian
	bandMatrixJacobian(x, res, matSize, lowerBand, upperBand);

	// Extract the banded Jacobian matrix from AD
	cadet::linalg::BandMatrix bm;
	bm.resize(matSize, lowerBand, upperBand);

	cadet::ad::extractBandedJacobianFromAd(res, 0, lowerBand, bm);

	// Get reference matrix
	const cadet::linalg::BandMatrix ref = createBandMatrix<cadet::linalg::BandMatrix>(matSize, lowerBand, upperBand);

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

	// Compute residual with banded Jacobian
	bandMatrixJacobian(x, res, matSize, lowerBand, upperBand);

	// Extract the banded Jacobian matrix from AD
	cadet::linalg::BandMatrix bm;
	bm.resize(matSize, lowerBand, upperBand);

	cadet::ad::extractBandedJacobianFromAd(res, 0, lowerBand, bm);

	cadet::linalg::DenseMatrix dm;
	dm.resize(3, 6);

	cadet::ad::extractDenseJacobianFromBandedAd(res, 0, 0, lowerBand, lowerBand, upperBand, dm);
	cadet::test::checkMatrixAgainstLinearArray(dm.data(), {1, 2, 3, 4, 0, 0, 5, 6, 7, 8, 9, 0, 10, 11, 12, 13, 14, 15});

	cadet::ad::extractDenseJacobianFromBandedAd(res, 1, 0, lowerBand, lowerBand, upperBand, dm);
	cadet::test::checkMatrixAgainstLinearArray(dm.data(), {6, 7, 8, 9, 0, 0, 11, 12, 13, 14, 15, 0, 16, 17, 18, 19, 20, 21});

	cadet::ad::extractDenseJacobianFromBandedAd(res, 2, 0, lowerBand, lowerBand, upperBand, dm);
	cadet::test::checkMatrixAgainstLinearArray(dm.data(), {12, 13, 14, 15, 0, 0, 17, 18, 19, 20, 21, 0, 22, 23, 24, 25, 26, 27});

	cadet::ad::extractDenseJacobianFromBandedAd(res, 5, 0, lowerBand, lowerBand, upperBand, dm);
	cadet::test::checkMatrixAgainstLinearArray(dm.data(), {30, 31, 32, 33, 0, 0, 35, 36, 37, 38, 39, 0, 40, 41, 42, 43, 44, 0});

	cadet::ad::extractDenseJacobianFromBandedAd(res, 6, 0, lowerBand, lowerBand, upperBand, dm);
	cadet::test::checkMatrixAgainstLinearArray(dm.data(), {36, 37, 38, 39, 0, 0, 41, 42, 43, 44, 0, 0, 45, 46, 47, 48, 0, 0});

	cadet::ad::extractDenseJacobianFromBandedAd(res, 7, 0, lowerBand, lowerBand, upperBand, dm);
	cadet::test::checkMatrixAgainstLinearArray(dm.data(), {42, 43, 44, 0, 0, 0, 46, 47, 48, 0, 0, 0, 49, 50, 51, 0, 0, 0});

	delete[] x;
	delete[] res;
}
