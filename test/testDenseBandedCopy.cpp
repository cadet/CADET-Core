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

#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>

#include "linalg/DenseMatrix.hpp"
#include "linalg/BandMatrix.hpp"

template <typename MatType_T>
void runTest()
{
	using cadet::linalg::DenseMatrix;

	std::cout << "Build matrix" << std::endl;

	MatType_T bm;
	bm.resize(8, 2, 3);
	std::cout << bm.rows() << " rows, Bandwidth " << bm.lowerBandwidth() << " + 1 + " << bm.upperBandwidth() << std::endl;  

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

	std::cout << bm << std::endl;

	// Copy with bottom left 0
	DenseMatrix dm1;
	dm1.resize(4, 4);
	dm1.copySubmatrixFromBanded(bm, 0, 0, 4, 4);
	std::cout << "Submatrix 0,0:\n" << dm1 << std::endl;

	// Copy with top right 0
	DenseMatrix dm2;
	dm2.resize(4, 4);
	dm2.copySubmatrixFromBanded(bm, 4, 0, 4, 4);
	std::cout << "Submatrix 4,0:\n" << dm2 << std::endl;

	// Copy from middle
	DenseMatrix dm3;
	dm3.resize(4, 4);
	dm3.copySubmatrixFromBanded(bm, 2, 2, 4, 4);
	std::cout << "Submatrix 2,2:\n" << dm3 << std::endl;

	// Copy from middle
	DenseMatrix dm4;
	dm4.resize(4, 4);
	dm4.copySubmatrixFromBanded(bm, 2, -1, 4, 4);
	std::cout << "Submatrix 2,-1:\n" << dm4 << std::endl;

	// Copy from middle
	DenseMatrix dm5;
	dm5.resize(4, 4);
	dm5.copySubmatrixFromBanded(bm, 1, 0, 4, 4);
	std::cout << "Submatrix 1,0:\n" << dm5 << std::endl;

	// Copy from middle
	DenseMatrix dm6;
	dm6.resize(4, 4);
	dm6.copySubmatrixFromBanded(bm, 0, 4, 4, 4);
	std::cout << "Submatrix 0,4:\n" << dm6 << std::endl;

	// Copy one row
	DenseMatrix dm7;
	dm7.resize(1, 2);
	dm7.copySubmatrixFromBanded(bm, 2, 1, 1, 2);
	std::cout << "Submatrix 2,1:\n" << dm7 << std::endl;

	// Copy one row
	DenseMatrix dm8;
	dm8.resize(1, 3);
	dm8.copySubmatrixFromBanded(bm, 2, -1, 1, 3);
	std::cout << "Submatrix 2,-1:\n" << dm8 << std::endl;
}

int main(int argc, char** argv)
{
	using cadet::linalg::BandMatrix;
	using cadet::linalg::FactorizableBandMatrix;

	std::cout << "BandMatrix ====================" << std::endl;
	runTest<BandMatrix>();

	std::cout << "\nFactorizableBandMatrix ====================" << std::endl;
	runTest<FactorizableBandMatrix>();

	return 0;
}
