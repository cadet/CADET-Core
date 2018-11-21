// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//
//  Copyright © 2008-2018: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "linalg/BandMatrix.hpp"
#include "linalg/DenseMatrix.hpp"
#include "linalg/SparseMatrix.hpp"
#include "AdUtils.hpp"

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
		// Clear previously set directions
		adVec[eq].fillADValue(adDirOffset, 0.0);
		// Set direction
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

void prepareAdVectorSeedsForDenseMatrix(active* const adVec, unsigned int adDirOffset, unsigned int cols)
{
	for (unsigned int col = 0; col < cols; ++col)
	{
		// Clear previously set directions
		adVec[col].fillADValue(adDirOffset, 0.0);
		// Set direction
		adVec[col].setADValue(adDirOffset + col, 1.0);
	}
}

void extractDenseJacobianFromAd(active const* const adVec, unsigned int adDirOffset, linalg::detail::DenseMatrixBase& mat)
{
	for (unsigned int eq = 0; eq < mat.rows(); ++eq)
	{
		// Loop over columns
		for (unsigned int col = 0; col < mat.columns(); ++col)
		{
			mat.native(eq, col) = adVec[eq].getADValue(adDirOffset + col);
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

double compareDenseJacobianWithAd(active const* const adVec, unsigned int adDirOffset, const linalg::detail::DenseMatrixBase& mat)
{
	double maxDiff = 0.0;
	for (unsigned int eq = 0; eq < mat.rows(); ++eq)
	{
		// Loop over columns
		for (unsigned int col = 0; col < mat.columns(); ++col)
		{
			double baseVal = adVec[eq].getADValue(adDirOffset + col);
			if (std::isnan(mat.native(eq, col)) || std::isnan(baseVal))
				return std::numeric_limits<double>::quiet_NaN();
			const double diff = std::abs(mat.native(eq, col) - baseVal);

			baseVal = std::abs(baseVal);
			if (baseVal > 0.0)
				maxDiff = std::max(maxDiff, diff / baseVal);
			else
				maxDiff = std::max(maxDiff, diff);
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

void adMatrixVectorMultiply(const linalg::SparseMatrix<active>& mat, double const* x, double* y, double alpha, double beta, unsigned int adDir)
{
	const std::vector<unsigned int>& rows = mat.rows();
	const std::vector<unsigned int>& cols = mat.cols();
	const std::vector<active>& values = mat.values();
	const unsigned int numNonZero = mat.numNonZero();

	for (unsigned int i = 0; i < numNonZero; ++i)
	{
		y[rows[i]] = alpha * values[i].getADValue(adDir) * x[cols[i]] + beta * y[rows[i]];
	}
}

DenseJacobianExtractor::DenseJacobianExtractor() { }

void DenseJacobianExtractor::extractJacobian(active const* adRes, unsigned int row, unsigned int adDirOffset, linalg::detail::DenseMatrixBase& mat) const
{
	extractDenseJacobianFromAd(adRes + row, adDirOffset + row, mat);
}

double DenseJacobianExtractor::compareWithJacobian(active const* adRes, unsigned int row, unsigned int adDirOffset, linalg::detail::DenseMatrixBase& mat) const
{
	return compareDenseJacobianWithAd(adRes + row, adDirOffset + row, mat);
}


BandedJacobianExtractor::BandedJacobianExtractor(unsigned int diagDir, unsigned int lowerBandwidth, unsigned int upperBandwidth)
	: _diagDir(diagDir), _lowerBandwidth(lowerBandwidth), _upperBandwidth(upperBandwidth)
	{ }

void BandedJacobianExtractor::extractJacobian(active const* adRes, unsigned int row, unsigned int adDirOffset, linalg::detail::DenseMatrixBase& mat) const
{
	extractDenseJacobianFromBandedAd(adRes, row, adDirOffset, _diagDir, _lowerBandwidth, _upperBandwidth, mat);
}

double BandedJacobianExtractor::compareWithJacobian(active const* adRes, unsigned int row, unsigned int adDirOffset, linalg::detail::DenseMatrixBase& mat) const
{
	return compareDenseJacobianWithBandedAd(adRes, row, adDirOffset, _diagDir, _lowerBandwidth, _upperBandwidth, mat);
}

}  // namespace ad

}  // namespace cadet
