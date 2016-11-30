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
#include <random>
#include <chrono>

#include "linalg/DenseMatrix.hpp"

int main(int argc, char** argv)
{
	using cadet::linalg::DenseMatrix;

	// Initialize standard normal RNG
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::normal_distribution<double> distribution(0.0, 1);

	std::cout << "Build matrix" << std::endl;

	DenseMatrix dm;
	dm.resize(4, 4);
	std::cout << dm.rows() << " size" << std::endl;  

	double val = 1.0;
	for (unsigned int row = 0; row < dm.rows(); ++row)
	{
		for (unsigned int col = 0; col < dm.columns(); ++col)
		{
			dm.native(row, col) = distribution(generator);
		}
	}

	std::cout << dm << std::endl;

	// Copy for QR
	DenseMatrix dm2 = dm;

	// Perform LU factorization and compute inverse matrix
	if (!dm.factorize())
		std::cout << "Failed to factorize matrix" << std::endl;

	std::cout << std::scientific << std::setprecision(std::numeric_limits<double>::digits10 + 1);
	std::vector<double> vecLU(dm.rows(), 0.0);
	std::cout << "invMatLU = [";
	for (unsigned int i = 0; i < dm.rows(); ++i)
	{
		vecLU[i] = 1.0;
		dm.solve(vecLU.data());

		for (unsigned int j = 0; j < dm.rows(); ++j)
		{
			std::cout << vecLU[j] << ", ";
			vecLU[j] = 0.0;
		}
		std::cout << "; ";
	}
	std::cout << "]';" << std::endl;

	// Perform QR factorization and compute inverse matrix
	std::vector<double> workingMemory(2 * dm2.rows(), 0.0);
	if (!dm2.robustFactorize(workingMemory.data()))
		std::cout << "Failed to factorize matrix" << std::endl;

	std::vector<double> vecQR(dm2.rows(), 0.0);
	std::cout << "invMatQR = [";
	for (unsigned int i = 0; i < dm2.rows(); ++i)
	{
		vecQR[i] = 1.0;
		dm2.robustSolve(vecQR.data(), workingMemory.data());

		for (unsigned int j = 0; j < dm2.rows(); ++j)
		{
			std::cout << vecQR[j] << ", ";
			vecQR[j] = 0.0;
		}
		std::cout << "; ";
	}
	std::cout << "]';" << std::endl;

	return 0;
}
