// =============================================================================
//  CADET
//
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file
 * Defines the ExchangeModel interface.
 */

#ifndef LIBCADET_EXCHANGEMODELINTERFACE_HPP_
#define LIBCADET_EXCHANGEMODELINTERFACE_HPP_

#include <unordered_map>

#include "CompileTimeConfig.hpp"
#include "cadet/ParameterProvider.hpp"
#include "cadet/ParameterId.hpp"
#include "linalg/DenseMatrix.hpp"
#include "linalg/BandMatrix.hpp"

#ifdef ENABLE_DG
	#include "linalg/BandedEigenSparseRowIterator.hpp"
#endif

#include "AutoDiff.hpp"
#include "SimulationTypes.hpp"
#include "Memory.hpp"

namespace cadet
{

class IParameterProvider;
class IExternalFunction;

struct ColumnPosition;

namespace model
{

class IExchangeModel
{
public:

	virtual ~IExchangeModel() CADET_NOEXCEPT { }

	virtual const char* name() const CADET_NOEXCEPT = 0;

	//virtual bool requiresConfiguration() const CADET_NOEXCEPT = 0;

	virtual bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int nChannel, unsigned int nCol) = 0;

	virtual bool configure(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx) = 0;

	virtual int residual(active const* y, active* res, WithParamSensitivity, bool wantJac, linalg::BandedSparseRowIterator jacBegin) const = 0;
	virtual int residual(active const* y, active* res, WithoutParamSensitivity, bool wantJac, linalg::BandedSparseRowIterator jacBegin) const = 0;
	virtual int residual(double const* y, active* res, WithParamSensitivity, bool wantJac, linalg::BandedSparseRowIterator jacBegin) const = 0;
	virtual int residual(double const* y, double* res, WithoutParamSensitivity, bool wantJac, linalg::BandedSparseRowIterator jacBegin) const = 0;

protected:
};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_EXCHANGEMODELINTERFACE_HPP_
