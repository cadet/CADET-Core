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

#define ACTIVE_SFAD

#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>

#include "linalg/DenseMatrix.hpp"
#include "linalg/BandMatrix.hpp"
#include "AdUtils.hpp"
#include "AutoDiff.hpp"

ACTIVE_INIT;

void bandMatrixJacobian(cadet::active const* x, cadet::active* out, unsigned int lowerBand, unsigned int upperBand, unsigned int rows)
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

int main(int argc, char** argv)
{
	using cadet::linalg::BandMatrix;

	const unsigned int matSize = 10;
	const unsigned int lowerBand = 2;
	const unsigned int upperBand = 3;

	cadet::ad::setDirections(lowerBand + 1 + upperBand);

	cadet::active* res = new cadet::active[matSize];
	cadet::active* x = new cadet::active[matSize];

	cadet::ad::prepareAdVectorSeedsForBandMatrix(x, 0, matSize, lowerBand, upperBand, lowerBand);

	bandMatrixJacobian(x, res, lowerBand, upperBand, matSize);

	BandMatrix bm;
	bm.resize(matSize, lowerBand, upperBand);

	cadet::ad::extractBandedJacobianFromAd(res, 0, lowerBand, bm);
	std::cout << "Bandmatrix from AD:\n" << bm << std::endl;

	cadet::linalg::DenseMatrix dm;
	dm.resize(3,6);

	cadet::ad::extractDenseJacobianFromBandedAd(res, 0, 0, lowerBand, lowerBand, upperBand, dm);
	std::cout << "Dense submatrix @0:\n" << dm << std::endl;

	cadet::ad::extractDenseJacobianFromBandedAd(res, 1, 0, lowerBand, lowerBand, upperBand, dm);
	std::cout << "Dense submatrix @1:\n" << dm << std::endl;

	cadet::ad::extractDenseJacobianFromBandedAd(res, 2, 0, lowerBand, lowerBand, upperBand, dm);
	std::cout << "Dense submatrix @2:\n" << dm << std::endl;

	cadet::ad::extractDenseJacobianFromBandedAd(res, 5, 0, lowerBand, lowerBand, upperBand, dm);
	std::cout << "Dense submatrix @5:\n" << dm << std::endl;

	cadet::ad::extractDenseJacobianFromBandedAd(res, 6, 0, lowerBand, lowerBand, upperBand, dm);
	std::cout << "Dense submatrix @6:\n" << dm << std::endl;

	cadet::ad::extractDenseJacobianFromBandedAd(res, 7, 0, lowerBand, lowerBand, upperBand, dm);
	std::cout << "Dense submatrix @7:\n" << dm << std::endl;

	delete[] x;
	delete[] res;

	return 0;
}
