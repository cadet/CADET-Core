// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//
//  Copyright © 2008-2016: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "linalg/BandMatrix.hpp"
#include "linalg/DenseMatrix.hpp"
#include "AutoDiff.hpp"

#include <limits>
#include <algorithm>

namespace cadet
{

namespace ad
{

void prepareAdVectorSeedsForBandMatrix(active* const adVec, unsigned int adDirOffset, unsigned int rows, 
	unsigned int lowerBandwidth, unsigned int upperBandwidth, unsigned int diagDir)
{
	// Start with diagonal Jacobian element
	unsigned int dir = diagDir;
	for (unsigned int eq = 0; eq < rows; ++eq)
	{
		adVec[eq].setADValue(adDirOffset + dir, 1.0);

		// Wrap around at end of row and jump to lowest subdiagonal
		if (dir == diagDir + upperBandwidth)
			dir = diagDir - lowerBandwidth;
		else
			++dir;
	}
}

void extractBandedJacobianFromAd(active const* const adVec, unsigned int adDirOffset, unsigned int diagDir, linalg::BandMatrix& mat)
{
	const unsigned int lowerBandwidth = mat.lowerBandwidth();
	const unsigned int upperBandwidth = mat.upperBandwidth();
	const unsigned int stride = lowerBandwidth + 1 + upperBandwidth;
	for (unsigned int eq = 0; eq < mat.rows(); ++eq)
	{
		// Start with lowest subdiagonal and stay in the range of the columns:
		// diagDir might be too big for the matrix and, hence, dir ranges between
		// diagDir - lowerBandwidth <= dir <= diagDir + upperBandwidth
		unsigned int dir = diagDir - lowerBandwidth + eq % stride;

		// Loop over diagonals
		for (unsigned int diag = 0; diag < stride; ++diag)
		{
			mat.native(eq, diag) = adVec[eq].getADValue(adDirOffset + dir);

			// Wrap around at end of row and jump to lowest subdiagonal
			if (dir == diagDir + upperBandwidth)
				dir = diagDir - lowerBandwidth;
			else
				++dir;
		}
	}
}

void extractDenseJacobianFromBandedAd(active const* const adVec, unsigned int row, unsigned int adDirOffset, unsigned int diagDir, 
	unsigned int lowerBandwidth, unsigned int upperBandwidth, linalg::detail::DenseMatrixBase& mat)
{
	const unsigned int stride = lowerBandwidth + 1 + upperBandwidth;
	for (unsigned int eq = 0; eq < mat.rows(); ++eq)
	{
		// Start with lowest subdiagonal and stay in the range of the columns:
		// diagDir might be too big for the matrix and, hence, dir ranges between
		// diagDir - lowerBandwidth <= dir <= diagDir + upperBandwidth
		unsigned int dir = diagDir - lowerBandwidth + ((eq + row) % stride);

		// Loop over diagonals
		for (unsigned int diag = 0; diag < stride; ++diag)
		{
			// Calculate column index from subdiagonal
			const int diagCol = static_cast<int>(diag) - static_cast<int>(lowerBandwidth) + static_cast<int>(eq);
			if ((diagCol >= 0) && (diagCol < static_cast<int>(mat.columns())))
				mat.native(eq, diagCol) = adVec[row + eq].getADValue(adDirOffset + dir);

			// Wrap around at end of row and jump to lowest subdiagonal
			if (dir == diagDir + upperBandwidth)
				dir = diagDir - lowerBandwidth;
			else
				++dir;
		}
	}
}

double compareBandedJacobianWithAd(active const* const adVec, unsigned int adDirOffset, unsigned int diagDir, const linalg::BandMatrix& mat)
{
	const unsigned int lowerBandwidth = mat.lowerBandwidth();
	const unsigned int upperBandwidth = mat.upperBandwidth();
	const unsigned int stride = lowerBandwidth + 1 + upperBandwidth;

	double maxDiff = 0.0;
	for (unsigned int eq = 0; eq < mat.rows(); ++eq)
	{
		// Start with lowest subdiagonal and stay in the range of the columns:
		// diagDir might be too big for the matrix and, hence, dir ranges between
		// diagDir - lowerBandwidth <= dir <= diagDir + upperBandwidth
		unsigned int dir = diagDir - lowerBandwidth + eq % stride;

		// Loop over diagonals
		for (unsigned int diag = 0; diag < stride; ++diag)
		{
			double baseVal = adVec[eq].getADValue(adDirOffset + dir);
			if (std::isnan(mat.native(eq, diag)) || std::isnan(baseVal))
				return std::numeric_limits<double>::quiet_NaN();
			const double diff = std::abs(mat.native(eq, diag) - baseVal);

			baseVal = std::abs(baseVal);
			if (baseVal > 0.0)
				maxDiff = std::max(maxDiff, diff / baseVal);
			else
				maxDiff = std::max(maxDiff, diff);

			// Wrap around at end of row and jump to lowest subdiagonal
			if (dir == diagDir + upperBandwidth)
				dir = diagDir - lowerBandwidth;
			else
				++dir;
		}
	}
	return maxDiff;
}

double compareDenseJacobianWithBandedAd(active const* const adVec, unsigned int row, unsigned int adDirOffset, unsigned int diagDir, 
	unsigned int lowerBandwidth, unsigned int upperBandwidth, const linalg::detail::DenseMatrixBase& mat)
{
	double maxDiff = 0.0;
	const unsigned int stride = lowerBandwidth + 1 + upperBandwidth;
	for (unsigned int eq = 0; eq < mat.rows(); ++eq)
	{
		// Start with lowest subdiagonal and stay in the range of the columns:
		// diagDir might be too big for the matrix and, hence, dir ranges between
		// diagDir - lowerBandwidth <= dir <= diagDir + upperBandwidth
		unsigned int dir = diagDir - lowerBandwidth + ((eq + row) % stride);

		// Loop over diagonals
		for (unsigned int diag = 0; diag < stride; ++diag)
		{
			// Calculate column index from subdiagonal
			const int diagCol = static_cast<int>(diag) - static_cast<int>(lowerBandwidth) + static_cast<int>(eq);
			if ((diagCol >= 0) && (diagCol < static_cast<int>(mat.columns())))
			{
				double baseVal = adVec[row + eq].getADValue(adDirOffset + dir);
				if (std::isnan(mat.native(eq, diagCol)) || std::isnan(baseVal))
					return std::numeric_limits<double>::quiet_NaN();
				const double diff = std::abs(mat.native(eq, diagCol) - baseVal);
				
				baseVal = std::abs(baseVal);
				if (baseVal > 0.0)
					maxDiff = std::max(maxDiff, diff / baseVal);
				else
					maxDiff = std::max(maxDiff, diff);
			}

			// Wrap around at end of row and jump to lowest subdiagonal
			if (dir == diagDir + upperBandwidth)
				dir = diagDir - lowerBandwidth;
			else
				++dir;
		}
	}
	return maxDiff;
}

}  // namespace ad

}  // namespace cadet
