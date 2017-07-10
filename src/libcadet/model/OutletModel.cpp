// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//
//  Copyright © 2008-2017: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "model/OutletModel.hpp"
#include "ParamReaderHelper.hpp"
#include "cadet/Exceptions.hpp"
#include "cadet/SolutionRecorder.hpp"

#include "ConfigurationHelper.hpp"
#include "ParamIdUtil.hpp"

#include <algorithm>
#include <functional>

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include "AdUtils.hpp"


inline void residual(double const* y, unsigned int nComp, double* res)
{
	std::copy(y, y + nComp, res);
}


namespace cadet
{

namespace model
{

OutletModel::OutletModel(UnitOpIdx unitOpIdx) : _unitOpIdx(unitOpIdx)
{
}

OutletModel::~OutletModel() CADET_NOEXCEPT
{
}

unsigned int OutletModel::numDofs() const CADET_NOEXCEPT
{
	return _nComp;
}

unsigned int OutletModel::numPureDofs() const CADET_NOEXCEPT
{
	return 0;
}

bool OutletModel::usesAD() const CADET_NOEXCEPT
{
	return false;
}

bool OutletModel::configure(IParameterProvider& paramProvider, IConfigHelper& helper)
{
	_nComp = paramProvider.getInt("NCOMP");
	return true;
}

bool OutletModel::reconfigure(IParameterProvider& paramProvider)
{
	return true;
}

void OutletModel::setSectionTimes(double const* secTimes, bool const* secContinuity, unsigned int nSections) { }

std::unordered_map<ParameterId, double> OutletModel::getAllParameterValues() const
{
	return std::unordered_map<ParameterId, double>();
}

bool OutletModel::hasParameter(const ParameterId& pId) const
{
	return false;
}

bool OutletModel::setParameter(const ParameterId& pId, int value)
{
	return false;
}

bool OutletModel::setParameter(const ParameterId& pId, double value)
{
	return false;
}

bool OutletModel::setParameter(const ParameterId& pId, bool value)
{
	return false;
}

void OutletModel::setSensitiveParameterValue(const ParameterId& pId, double value)
{
}

bool OutletModel::setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue)
{
	return false;
}

void OutletModel::clearSensParams()
{
}

void OutletModel::useAnalyticJacobian(const bool analyticJac) { }
void OutletModel::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active* const adRes, active* const adY, unsigned int adDirOffset) { }

void OutletModel::reportSolution(ISolutionRecorder& recorder, double const* const solution) const
{
	Exporter expr(_nComp, solution);
	recorder.beginUnitOperation(_unitOpIdx, *this, expr);
	recorder.endUnitOperation();
}

void OutletModel::reportSolutionStructure(ISolutionRecorder& recorder) const
{
	Exporter expr(_nComp, nullptr);
	recorder.unitOperationStructure(_unitOpIdx, *this, expr);
}

unsigned int OutletModel::requiredADdirs() const CADET_NOEXCEPT
{
	return 0;
}

void OutletModel::prepareADvectors(active* const adRes, active* const adY, unsigned int adDirOffset) const { }

void OutletModel::applyInitialCondition(double* const vecStateY, double* const vecStateYdot)
{
	std::fill(vecStateY, vecStateY + _nComp, 0.0);
	std::fill(vecStateYdot, vecStateYdot + _nComp, 0.0);
}

void OutletModel::applyInitialCondition(IParameterProvider& paramProvider, double* const vecStateY, double* const vecStateYdot)
{
	std::fill(vecStateY, vecStateY + _nComp, 0.0);
	std::fill(vecStateYdot, vecStateYdot + _nComp, 0.0);
}

int OutletModel::residual(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double* const res)
{
	::residual(y, _nComp, res);
	return 0;
}

int OutletModel::residualWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, 
	double* const res, active* const adRes, active* const adY, unsigned int adDirOffset)
{
	// Jacobian is always identity
	::residual(y, _nComp, res);
	return 0;
}

int OutletModel::residualSensFwdAdOnly(const active& t, unsigned int secIdx, const active& timeFactor,
	double const* const y, double const* const yDot, active* const adRes)
{
	for (unsigned int i = 0; i < _nComp; ++i)
		adRes[i] = y[i];

	return 0;
}

int OutletModel::residualSensFwdCombine(const active& timeFactor, const std::vector<const double*>& yS, const std::vector<const double*>& ySdot,
	const std::vector<double*>& resS, active const* adRes, double* const tmp1, double* const tmp2, double* const tmp3)
{
	// Directional derivative (dF / dy) * s does nothing since dF / dy = I (identity)
	for (unsigned int param = 0; param < resS.size(); ++param)
	{
		double const* const y = yS[param];
		double* const res = resS[param];
		std::copy(y, y + _nComp, res);
	}
	return 0;
}

int OutletModel::residualSensFwdWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, 
	double const* const yDot, active* const adRes, active* const adY, unsigned int adDirOffset)
{
	for (unsigned int i = 0; i < _nComp; ++i)
		adRes[i] = y[i];

	return 0;
}

void OutletModel::consistentInitialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes)
{
	// Nothing to do here as inlet DOFs are initialized by ModelSystem
}

void OutletModel::leanConsistentInitialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes)
{
	// Nothing to do here as inlet DOFs are initialized by ModelSystem
}

void OutletModel::multiplyWithJacobian(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double const* yS, double alpha, double beta, double* ret)
{
	// dF / dy = I (identity matrix)
	for (unsigned int i = 0; i < _nComp; ++i)
	{
		ret[i] = alpha * yS[i] + beta * ret[i];
	}
}

void OutletModel::multiplyWithDerivativeJacobian(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double const* sDot, double* ret)
{
	std::fill_n(ret, numDofs(), 0.0);
}

}  // namespace model

}  // namespace cadet
