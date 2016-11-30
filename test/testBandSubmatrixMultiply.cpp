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
#include <cmath>

#include "linalg/DenseMatrix.hpp"
#include "linalg/BandMatrix.hpp"

void printMatrix(double const* mat, unsigned int rows, unsigned int cols, bool candidate)
{
	if (candidate)
		std::cout << "out = [" << mat[0];
	else
		std::cout << "ref = [" << mat[0];

	for (unsigned int i = 1; i < cols; ++i)
	{
		std::cout << ", " << mat[i];
	}
	for (unsigned int i = 1; i < rows; ++i)
	{
		std::cout << "; " << mat[i * cols];
		for (unsigned int j = 1; j < cols; ++j)
		{
			std::cout << ", " << mat[i * cols + j];
		}
	}
	std::cout << "];\n";
}

bool checkMatrix(double const* mat, const std::vector<double>& matRef)
{
	for (unsigned int i = 0; i < matRef.size(); ++i)
	{
		if (std::abs(matRef[i] - mat[i]) >= 1e-14)
		{
			std::cout << " => FAILED at element " << i << "\n";
			return false;
		}
	}
	std::cout << " => PASSED\n";
	return true;
}

void recoverMatrix(const cadet::linalg::BandMatrix& mat, unsigned int startRow, int startDiag, unsigned int numRows, unsigned int numCols, double* const out = nullptr)
{
	std::vector<double> x(numCols, 0.0);
	std::vector<double> y(numRows, 0.0);

	if (!out)
		std::cout << "B = [";

	for (unsigned int i = 0; i < numCols; ++i)
	{
		// Multiply with basis vectors
		x[i] = 1.0;
		mat.submatrixMultiplyVector(x.data(), startRow, startDiag, numRows, numCols, y.data());
		x[i] = 0.0;

		if (out)
		{
			for (unsigned int j = 0; j < numRows; ++j)
				out[i + j * numCols] = y[j];
		}
		else
		{
			std::cout << y[0];
			for (unsigned int j = 1; j < numRows; ++j)
				std::cout << "," << y[j];
			std::cout << ";";
		}
	}
	if (!out)
		std::cout << "]';" << std::endl;
}

int main(int argc, char** argv)
{
	using cadet::linalg::BandMatrix;

	{
		std::cout << "Build matrix" << std::endl;

		BandMatrix bm;
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

		// Perform some multiplications
		std::vector<double> out(8 * 6, 0.0);

		{
			std::cout << "0,0 - 5,5\n";
			recoverMatrix(bm, 0, 0, 5, 5, out.data());
			std::vector<double> cmp {1, 2, 3, 4, 0, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 0, 16, 17, 18, 19, 0, 0, 22, 23, 24};
			if (!checkMatrix(out.data(), cmp))
			{
				printMatrix(out.data(), 5, 5, true);
				printMatrix(cmp.data(), 5, 5, false);
			}
		}

		{
			std::cout << "2,0 - 5,5\n";
			recoverMatrix(bm, 2, 0, 5, 5, out.data());
			std::vector<double> cmp {12, 13, 14, 15, 0, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 0, 28, 29, 30, 31, 0, 0, 33, 34, 35};
			if (!checkMatrix(out.data(), cmp))
			{
				printMatrix(out.data(), 5, 5, true);
				printMatrix(cmp.data(), 5, 5, false);
			}
		}

		{
			std::cout << "2,-2 -, 5,5\n";
			recoverMatrix(bm, 2, -2, 5, 5, out.data());
			std::vector<double> cmp {10, 11, 12, 13, 14, 0, 16, 17, 18, 19, 0, 0, 22, 23, 24, 0, 0, 0, 28, 29, 0, 0, 0, 0, 33};
			if (!checkMatrix(out.data(), cmp))
			{
				printMatrix(out.data(), 5, 5, true);
				printMatrix(cmp.data(), 5, 5, false);
			}
		}

		{
			std::cout << "2,-1 -, 5,5\n";
			recoverMatrix(bm, 2, -1, 5, 5, out.data());
			std::vector<double> cmp {11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 0, 22, 23, 24, 25, 0, 0, 28, 29, 30, 0, 0, 0, 33, 34};
			if (!checkMatrix(out.data(), cmp))
			{
				printMatrix(out.data(), 5, 5, true);
				printMatrix(cmp.data(), 5, 5, false);
			}
		}

		{
			std::cout << "2,1 - 5,5\n";
			recoverMatrix(bm, 2, 1, 5, 5, out.data());
			std::vector<double> cmp {13, 14, 15, 0, 0, 18, 19, 20, 21, 0, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 0, 33, 34, 35, 36};
			if (!checkMatrix(out.data(), cmp))
			{
				printMatrix(out.data(), 5, 5, true);
				printMatrix(cmp.data(), 5, 5, false);
			}
		}

		{
			std::cout << "2,3 - 5,5\n";
			recoverMatrix(bm, 2, 3, 5, 5, out.data());
			std::vector<double> cmp {15, 0, 0, 0, 0, 20, 21, 0, 0, 0, 25, 26, 27, 0, 0, 30, 31, 32, 0, 0, 34, 35, 36, 0, 0};
			if (!checkMatrix(out.data(), cmp))
			{
				printMatrix(out.data(), 5, 5, true);
				printMatrix(cmp.data(), 5, 5, false);
			}
		}

		{
			std::cout << "2,0 - 5,7\n";
			recoverMatrix(bm, 2, 0, 5, 7, out.data());
			std::vector<double> cmp {12, 13, 14, 15, 0, 0, 0, 17, 18, 19, 20, 21, 0, 0, 22, 23, 24, 25, 26, 27, 0, 0, 28, 29, 30, 31, 32, 0, 0, 0, 33, 34, 35, 36, 0};
			if (!checkMatrix(out.data(), cmp))
			{
				printMatrix(out.data(), 5, 7, true);
				printMatrix(cmp.data(), 5, 7, false);
			}
		}

		{
			std::cout << "2,-1 -, 5,3\n";
			recoverMatrix(bm, 2, -1, 5, 3, out.data());
			std::vector<double> cmp {11, 12, 13, 16, 17, 18, 0, 22, 23, 0, 0, 28, 0, 0, 0};
			if (!checkMatrix(out.data(), cmp))
			{
				printMatrix(out.data(), 5, 3, true);
				printMatrix(cmp.data(), 5, 3, false);
			}
		}

		{
			std::cout << "2,1 - 5,3\n";
			recoverMatrix(bm, 2, 1, 5, 3, out.data());
			std::vector<double> cmp {13, 14, 15, 18, 19, 20, 23, 24, 25, 28, 29, 30, 0, 33, 34};
			if (!checkMatrix(out.data(), cmp))
			{
				printMatrix(out.data(), 5, 3, true);
				printMatrix(cmp.data(), 5, 3, false);
			}
		}

		{
			std::cout << "2,2 - 5,3\n";
			recoverMatrix(bm, 2, 2, 5, 3, out.data());
			std::vector<double> cmp {14, 15, 0, 19, 20, 21, 24, 25, 26, 29, 30, 31, 33, 34, 35};
			if (!checkMatrix(out.data(), cmp))
			{
				printMatrix(out.data(), 5, 3, true);
				printMatrix(cmp.data(), 5, 3, false);
			}
		}

		{
			std::cout << "2,3 - 5,3\n";
			recoverMatrix(bm, 2, 3, 5, 3, out.data());
			std::vector<double> cmp {15, 0, 0, 20, 21, 0, 25, 26, 27, 30, 31, 32, 34, 35, 36};
			if (!checkMatrix(out.data(), cmp))
			{
				printMatrix(out.data(), 5, 3, true);
				printMatrix(cmp.data(), 5, 3, false);
			}
		}

		{
			std::cout << "2,-2 - 3,3\n";
			recoverMatrix(bm, 2, -2, 3, 3, out.data());
			std::vector<double> cmp {10, 11, 12, 0, 16, 17, 0, 0, 22};
			if (!checkMatrix(out.data(), cmp))
			{
				printMatrix(out.data(), 3, 3, true);
				printMatrix(cmp.data(), 3, 3, false);
			}
		}

		{
			std::cout << "2,1 - 3,3\n";
			recoverMatrix(bm, 2, 1, 3, 3, out.data());
			std::vector<double> cmp {13, 14, 15, 18, 19, 20, 23, 24, 25};
			if (!checkMatrix(out.data(), cmp))
			{
				printMatrix(out.data(), 3, 3, true);
				printMatrix(cmp.data(), 3, 3, false);
			}
		}

		{
			std::cout << "2,1 - 3,5\n";
			recoverMatrix(bm, 2, 1, 3, 5, out.data());
			std::vector<double> cmp {13, 14, 15, 0, 0, 18, 19, 20, 21, 0, 23, 24, 25, 26, 27};
			if (!checkMatrix(out.data(), cmp))
			{
				printMatrix(out.data(), 3, 5, true);
				printMatrix(cmp.data(), 3, 5, false);
			}
		}

		{
			std::cout << "3,0 - 1,2\n";
			recoverMatrix(bm, 3, 0, 1, 2, out.data());
			std::vector<double> cmp {18, 19};
			if (!checkMatrix(out.data(), cmp))
			{
				printMatrix(out.data(), 1, 2, true);
				printMatrix(cmp.data(), 1, 2, false);
			}
		}

		{
			std::cout << "3,1 - 1,2\n";
			recoverMatrix(bm, 3, 1, 1, 2, out.data());
			std::vector<double> cmp {19, 20};
			if (!checkMatrix(out.data(), cmp))
			{
				printMatrix(out.data(), 1, 2, true);
				printMatrix(cmp.data(), 1, 2, false);
			}
		}
	}

	{
		std::cout << "Build matrix" << std::endl;

		BandMatrix bm;
		bm.resize(24, 6, 9);
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
		
		// Perform some multiplications
		std::vector<double> out(24*16, 0.0);

		{
			std::cout << "3,1 - 1,2\n";
			recoverMatrix(bm, 3, 1, 1, 2, out.data());
			std::vector<double> cmp {38, 39};
			if (!checkMatrix(out.data(), cmp))
			{
				printMatrix(out.data(), 1, 2, true);
				printMatrix(cmp.data(), 1, 2, false);
			}
		}

		{
			std::cout << "3,-1 - 1,3\n";
			recoverMatrix(bm, 3, -1, 1, 3, out.data());
			std::vector<double> cmp {36, 37, 38};
			if (!checkMatrix(out.data(), cmp))
			{
				printMatrix(out.data(), 1, 3, true);
				printMatrix(cmp.data(), 1, 3, false);
			}
		}
	}

	return 0;
}
