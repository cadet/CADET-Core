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

#include "model/OutletModel.hpp"
#include "ParamReaderHelper.hpp"
#include "cadet/Exceptions.hpp"
#include "cadet/SolutionRecorder.hpp"

#include "ConfigurationHelper.hpp"
#include "ParamIdUtil.hpp"
#include "linalg/Norms.hpp"

#include <algorithm>
#include <functional>

#include "LoggingUtils.hpp"
#include "Logging.hpp"


#include "GeneralRateModel.hpp"

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
//	return _nComp;
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

void OutletModel::setSensitiveParameterValue(const ParameterId& pId, double value) { }

bool OutletModel::setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue)
{
	return false;
}

void OutletModel::clearSensParams() { }
void OutletModel::useAnalyticJacobian(const bool analyticJac) { }
void OutletModel::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx) { }

void OutletModel::reportSolution(ISolutionRecorder& recorder, double const* const solution, const GeneralRateModel& grm) const
{
	Exporter expr(_nComp, solution + grm.localOutletComponentIndex(), grm.localOutletComponentStride());
	recorder.beginUnitOperation(_unitOpIdx, *this, expr);
	recorder.endUnitOperation();
}

void OutletModel::reportSolution(ISolutionRecorder& recorder, double const* const solution) const
{
	Exporter expr(_nComp, solution, 1);
	recorder.beginUnitOperation(_unitOpIdx, *this, expr);
	recorder.endUnitOperation();
}

void OutletModel::reportSolutionStructure(ISolutionRecorder& recorder) const
{
	Exporter expr(_nComp, nullptr, 0);
	recorder.unitOperationStructure(_unitOpIdx, *this, expr);
}

unsigned int OutletModel::requiredADdirs() const CADET_NOEXCEPT
{
	return 0;
}

void OutletModel::prepareADvectors(active* const adRes, active* const adY, unsigned int numSensAdDirs) const { }

void OutletModel::applyInitialCondition(double* const vecStateY, double* const vecStateYdot)
{
//	std::fill(vecStateY, vecStateY + _nComp, 0.0);
//	std::fill(vecStateYdot, vecStateYdot + _nComp, 0.0);
}

void OutletModel::applyInitialCondition(IParameterProvider& paramProvider, double* const vecStateY, double* const vecStateYdot)
{
//	std::fill(vecStateY, vecStateY + _nComp, 0.0);
//	std::fill(vecStateYdot, vecStateYdot + _nComp, 0.0);
}

int OutletModel::residual(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double* const res)
{
//	std::copy(y, y + _nComp, res);
	return 0;
}

int OutletModel::residualWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, 
	double* const res, active* const adRes, active* const adY, unsigned int numSensAdDirs)
{
//	std::copy(y, y + _nComp, res);
	return 0;
}

double OutletModel::residualNorm(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot)
{
//	return linalg::linfNorm(y, _nComp);
	return 0.0;
}

int OutletModel::residualSensFwdAdOnly(const active& t, unsigned int secIdx, const active& timeFactor,
	double const* const y, double const* const yDot, active* const adRes)
{
/*
	for (unsigned int i = 0; i < _nComp; ++i)
		adRes[i] = y[i];
*/
	return 0;
}

int OutletModel::residualSensFwdCombine(const active& timeFactor, const std::vector<const double*>& yS, const std::vector<const double*>& ySdot,
	const std::vector<double*>& resS, active const* adRes, double* const tmp1, double* const tmp2, double* const tmp3)
{
/*
	for (unsigned int i = 0; i < resS.size(); ++i)
		std::fill(resS[i], resS[i] + _nComp, 0.0);
*/
	return 0;
}

int OutletModel::residualSensFwd(unsigned int nSens, const active& t, unsigned int secIdx,
	const active& timeFactor, double const* const y, double const* const yDot, double const* const res,
	const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS,
	active* const adRes, double* const tmp1, double* const tmp2, double* const tmp3)
{
/*
	if (adRes)
	{
		for (unsigned int i = 0; i < _nComp; ++i)
		{
			adRes[i] = 0.0;
		}
	}

	for (unsigned int i = 0; i < resS.size(); ++i)
		std::fill(resS[i], resS[i] + _nComp, 0.0);
*/
	return 0;
}

int OutletModel::residualSensFwdWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, 
	double const* const yDot, active* const adRes, active* const adY, unsigned int numSensAdDirs)
{
/*
	if (adRes)
	{
		for (unsigned int i = 0; i < _nComp; ++i)
		{
			adRes[i] = 0.0;
		}
	}
*/
	return 0;
}

void OutletModel::residualSensFwdNorm(unsigned int nSens, const active& t, unsigned int secIdx,
		const active& timeFactor, double const* const y, double const* const yDot,
		const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, double* const norms,
		active* const adRes, double* const tmp)
{
//	std::fill(norms, norms + nSens, 0.0);
}

int OutletModel::linearSolve(double t, double timeFactor, double alpha, double tol, double* const rhs, double const* const weight,
		double const* const y, double const* const yDot, double const* const res)
{
//	std::fill(rhs, rhs + _nComp, 0.0);
	return 0;
}

void OutletModel::consistentIntialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active* const adRes, active* const adY)
{
/*
	for (unsigned int i = 0; i < vecSensY.size(); ++i)
	{
		std::fill(vecSensY[i], vecSensY[i] + _nComp, 0.0);
		std::fill(vecSensYdot[i], vecSensYdot[i] + _nComp, 0.0);	
	}
*/
}

void OutletModel::consistentIntialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes)
{
/*
	for (unsigned int i = 0; i < vecSensY.size(); ++i)
	{
		std::fill(vecSensY[i], vecSensY[i] + _nComp, 0.0);
		std::fill(vecSensYdot[i], vecSensYdot[i] + _nComp, 0.0);	
	}
*/
}

active OutletModel::inletConnectionFactorActive(unsigned int compIdx, unsigned int secIdx) const CADET_NOEXCEPT
{
	return active(1.0);
}

double OutletModel::inletConnectionFactor(unsigned int compIdx, unsigned int secIdx) const CADET_NOEXCEPT
{
	return 1.0;
}

}  // namespace model

}  // namespace cadet
