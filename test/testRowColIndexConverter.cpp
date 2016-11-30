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

#include "common/OrderingConverter.hpp"

template <typename T>
std::ostream& printVector(const std::vector<T>& vec)
{
	std::cout << "[" << vec[0];
	for (unsigned int i = 1; i < vec.size(); ++i)
		std::cout << ", " << vec[i];
	std::cout << "];";
	return std::cout;
}

int main(int argc, char** argv)
{
	std::cout << "Test 3x2x4" << std::endl;
	std::vector<double> dataRow(3*2*4, 0.0);
	std::vector<double> dataCol(3*2*4, 0.0);

	std::vector<unsigned int> dims({3,2,4});
	std::cout << "dims = "; printVector(dims) << std::endl;

	cadet::OrderingConverter conv(dims);

	// Create data
	std::vector<unsigned int> subscript = dims;
	for (unsigned int i = 0; i < dataRow.size(); ++i)
	{
		const unsigned int linearCol = conv.rowToCol(i);

		dataRow[i] = static_cast<double>(i);
		dataCol[linearCol] = static_cast<double>(i);
	}
	std::cout << "dataRow = "; printVector(dataRow) << std::endl;
	std::cout << "dataCol = "; printVector(dataCol) << std::endl;

	// Compare
	for (unsigned int i = 0; i < dataRow.size(); ++i)
	{
		const unsigned int linearRow = conv.colToRow(i);

		if (dataCol[i] != dataRow[linearRow])
		{
			std::cout << "==== Error in col-major index " << i << " "; printVector(conv.subscriptIndex()) << std::endl;
			std::cout << "  => row-major index " << linearRow << std::endl;
			std::cout << "  => " << dataCol[i] << " (col) vs " << dataRow[linearRow] << " (row)" << std::endl;
		}
	}

	return 0;
}
