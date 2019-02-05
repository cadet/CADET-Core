// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2019: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "umfpack.h"

#include "linalg/UMFPackSparseMatrix.hpp"
#include <algorithm>

namespace cadet
{

namespace linalg
{

UMFPackSparseMatrix::UMFPackSparseMatrix() : CompressedSparseMatrix(), 
	_symbolic(nullptr), _numeric(nullptr), _info(UMFPACK_INFO, 0.0), _options(UMFPACK_CONTROL, 0.0), _result(0), _workSpaceIdx(0, 0), _workSpace(0, 0.0)
{
	umfpack_di_defaults(_options.data());
	_options[UMFPACK_ORDERING] = UMFPACK_ORDERING_AMD;
}

UMFPackSparseMatrix::UMFPackSparseMatrix(unsigned int numRows, unsigned int numNonZeros) : CompressedSparseMatrix(numRows, numNonZeros), 
	_symbolic(nullptr), _numeric(nullptr), _info(UMFPACK_INFO, 0.0), _options(UMFPACK_CONTROL, 0.0), _result(numRows, 0.0), _workSpaceIdx(numRows, 0), _workSpace(numRows, 0.0)
{
	umfpack_di_defaults(_options.data());
	_options[UMFPACK_ORDERING] = UMFPACK_ORDERING_AMD;
}

UMFPackSparseMatrix::~UMFPackSparseMatrix() CADET_NOEXCEPT
{
	if (_symbolic)
		umfpack_di_free_symbolic(&_symbolic);

	if (_numeric)
		umfpack_di_free_numeric(&_numeric);
}

void UMFPackSparseMatrix::prepare()
{
	if (_result.size() == 0)
	{
		_result.resize(rows(), 0.0);
		_workSpaceIdx.resize(rows(), 0);
		if (_options[UMFPACK_IRSTEP] > 0)
			_workSpace.resize(5 * rows(), 0.0);
		else
			_workSpace.resize(rows(), 0.0);
	}

	if (_symbolic)
	{
		// This call also sets _symbolic to nullptr
		umfpack_di_free_symbolic(&_symbolic);
	}

	const int status = umfpack_di_symbolic(rows(), rows(), _rowStart.data(), _colIdx.data(), _values.data(), &_symbolic, _options.data(), _info.data());
	if (status != UMFPACK_OK)
	{
//		umfpack_di_report_info(_options.data(), _info.data());
//		umfpack_di_report_status(_options.data(), status);
	}
}

bool UMFPackSparseMatrix::factorize()
{
	if (_numeric)
	{
		// This call also sets _numeric to nullptr
		umfpack_di_free_numeric(&_numeric);
	}

	const int status = umfpack_di_numeric(_rowStart.data(), _colIdx.data(), _values.data(), _symbolic, &_numeric, _options.data(), _info.data());
	if (status != UMFPACK_OK)
	{
//		umfpack_di_report_info(_options.data(), _info.data());
//		umfpack_di_report_status(_options.data(), status);
	}

	return status == UMFPACK_OK;
}

bool UMFPackSparseMatrix::solve(double* rhs) const
{
	const int status = umfpack_di_wsolve(UMFPACK_Aat, _rowStart.data(), _colIdx.data(), _values.data(), _result.data(), rhs, _numeric, _options.data(), _info.data(), _workSpaceIdx.data(), _workSpace.data());
	std::copy(_result.begin(), _result.end(), rhs);

	return status == UMFPACK_OK;
}

}  // namespace linalg

}  // namespace cadet
