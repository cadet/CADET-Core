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

#include "linalg/BandMatrix.hpp"

int main(int argc, char** argv)
{
	using cadet::linalg::BandMatrix;
	using cadet::linalg::FactorizableBandMatrix;

	std::cout << "Build matrix" << std::endl;

	BandMatrix bm;
	bm.resize(10, 2, 3);
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
	
	FactorizableBandMatrix fbm;
	fbm.resize(bm.rows(), bm.lowerBandwidth(), bm.upperBandwidth());
	fbm.copyOver(bm);
	std::cout << "\n" << fbm << "\n" << std::endl;

	FactorizableBandMatrix::RowIterator it = fbm.row(0);
	for (unsigned int row = 0; row < bm.rows(); ++row, ++it)
	{
		const int lower = std::max(-static_cast<int>(bm.lowerBandwidth()), -static_cast<int>(row));
		const int upper = std::min(static_cast<int>(bm.upperBandwidth()), static_cast<int>(bm.rows() - row) - 1);
		for (int col = lower; col <= upper; ++col)
		{
			if (fbm.centered(row, col) != it[col])
				std::cout << "Error in row " << row << " col " << col << ": " << fbm.centered(row, col) << " vs " << it[col] << std::endl;
		}
	}


	if (!fbm.factorize())
		std::cout << "Failed to factorize matrix" << std::endl;
	std::cout << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1);
	std::vector<double> vec(bm.rows(), 0.0);
	std::cout << "invMat = [";
	for (unsigned int i = 0; i < bm.rows(); ++i)
	{
		vec[i] = 1.0;
		fbm.solve(vec.data());

		for (unsigned int j = 0; j < bm.rows(); ++j)
		{
			std::cout << vec[j] << ", ";
			vec[j] = 0.0;
		}
		std::cout << "; ";
	}
	std::cout << "];" << std::endl;

	return 0;
}
