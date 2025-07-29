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

	virtual bool hasQuasiStationary(int comp) const{ return false;}
	
	virtual void quasiStationarityMap(int comp, std::vector<std::pair<unsigned int, unsigned int>>& qsOrgDestMask) const {}

	virtual int residual(active const* y, active* res, WithParamSensitivity, bool wantJac, linalg::BandedSparseRowIterator jacBegin) const = 0;
	virtual int residual(active const* y, active* res, WithoutParamSensitivity, bool wantJac, linalg::BandedSparseRowIterator jacBegin) const = 0;
	virtual int residual(double const* y, active* res, WithParamSensitivity, bool wantJac, linalg::BandedSparseRowIterator jacBegin) const = 0;
	virtual int residual(double const* y, double* res, WithoutParamSensitivity, bool wantJac, linalg::BandedSparseRowIterator jacBegin) const = 0;

	virtual int numComp() const { return 0; }
	virtual int numChannel() const  { return 0; }
	virtual int numColums() const { return 0; }

	virtual double getCrossSectionRation(int idxOrig, int idxDest) const { return 0.0; }

	virtual void timeDerivativeQuasiStationaryExchange(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* y, double* dResDt, LinearBufferAllocator workSpace) const {}

	virtual bool dependsOnTime() const { return false; }

protected:
};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_EXCHANGEMODELINTERFACE_HPP_
