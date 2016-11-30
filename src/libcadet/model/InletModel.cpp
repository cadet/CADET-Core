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

#include "model/InletModel.hpp"
#include "ParamReaderHelper.hpp"
#include "cadet/Exceptions.hpp"
#include "cadet/InletProfile.hpp"
#include "cadet/SolutionRecorder.hpp"

#include "ConfigurationHelper.hpp"
#include "ParamIdUtil.hpp"

#include <algorithm>
#include <functional>

#include "LoggingUtils.hpp"
#include "Logging.hpp"

namespace cadet
{

namespace model
{

InletModel::InletModel(UnitOpIdx unitOpIdx) : _unitOpIdx(unitOpIdx), _inlet(nullptr),
	_inletConcentrationsRaw(nullptr), _inletDerivatives(nullptr), _inletConcentrations(nullptr)
{
}

InletModel::~InletModel() CADET_NOEXCEPT
{
	delete[] _inletConcentrationsRaw;
	// Do not delete _inletDerivatives since its memory is owned by _inletConcentrationsRaw
	delete[] _inletConcentrations;

	delete _inlet;
}

unsigned int InletModel::numDofs() const CADET_NOEXCEPT
{
	return 0;
}

bool InletModel::usesAD() const CADET_NOEXCEPT
{
	return false;
}

bool InletModel::configure(IParameterProvider& paramProvider, IConfigHelper& helper)
{
	_nComp = paramProvider.getInt("NCOMP");
	const std::string inType = paramProvider.getString("INLET_TYPE");

	_inlet = helper.createInletProfile(inType);
	if (!_inlet)
		throw InvalidParameterException("INLET_TYPE " + inType + " is unknown");

	_inlet->configure(&paramProvider, _nComp);

	// Allocate memory
	_inletConcentrationsRaw = new double[3 * _nComp]; // Thrice the size because _inletDerivatives uses two thirds
	_inletDerivatives = _inletConcentrationsRaw + _nComp;
	_inletConcentrations = new active[_nComp];

	return true;
}

bool InletModel::reconfigure(IParameterProvider& paramProvider)
{
	return _inlet->configure(&paramProvider, _nComp);
}

void InletModel::setSectionTimes(double const* secTimes, bool const* secContinuity, unsigned int nSections)
{
	if (_inlet)
		_inlet->setSectionTimes(secTimes, secContinuity, nSections);
}

std::unordered_map<ParameterId, double> InletModel::getAllParameterValues() const
{
	std::unordered_map<ParameterId, double> data;
	if (!_inlet)
		return data;

	const std::vector<ParameterId> inletParams = _inlet->availableParameters(_unitOpIdx);
	for (const ParameterId& pId : inletParams)
	{
		// TODO: Query value here
		data[pId] = 0.0;
	}

	return data;
}

bool InletModel::hasParameter(const ParameterId& pId) const
{
	if (!_inlet || ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep)))
		return false;
	
	const std::vector<ParameterId> inletParams = _inlet->availableParameters(_unitOpIdx);

	// Search for parameter
	for (const ParameterId& ipid : inletParams)
	{
		if (ipid == pId)
			return true;
	}
	return false;
}

bool InletModel::setParameter(const ParameterId& pId, int value)
{
	return false;
}

bool InletModel::setParameter(const ParameterId& pId, double value)
{
	// Check inlet and filter parameters
	if (_inlet && ((pId.unitOperation == _unitOpIdx) || (pId.unitOperation == UnitOpIndep)))
	{
		_inlet->setParameterValue(pId, value);
		return true;
	}

	return false;
}

bool InletModel::setParameter(const ParameterId& pId, bool value)
{
	return false;
}

void InletModel::setSensitiveParameterValue(const ParameterId& pId, double value)
{
	// Check inlet and filter parameters
	if (_inlet && ((pId.unitOperation == _unitOpIdx) || (pId.unitOperation == UnitOpIndep)))
		_inlet->setParameterValue(pId, value);
}

bool InletModel::setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue)
{
	if (_inlet && ((pId.unitOperation == _unitOpIdx) || (pId.unitOperation == UnitOpIndep)))
	{
		const std::vector<ParameterId> inletParams = _inlet->availableParameters(_unitOpIdx);

		// Search for parameter
		bool found = false;
		for (ParameterId ipid : inletParams)
		{
			if (ipid == pId)
			{
				found = true;
				break;
			}
		}

		if (found)
		{
			if (_sensParamsInlet.find(pId) == _sensParamsInlet.end())
			{
				LOG(Debug) << "Found parameter " << pId << " in Inlet: Dir " << adDirection << " is set to " << adValue;
			}
			
			// Register parameter and reserve AD direction
			_sensParamsInlet[pId] = std::make_tuple(adDirection, adValue);
			return true;
		}
	}

	return false;
}

void InletModel::clearSensParams()
{
	_sensParamsInlet.clear();
}

template<> double const* InletModel::moveInletValues(double const* const rawValues, const active& t, unsigned int secIdx) const
{
	// No parameter derivatives required, return raw values
	return rawValues;
}

template<> active const* InletModel::moveInletValues(double const* const rawValues, const active& t, unsigned int secIdx) const
{
	// Convert to active
	for (unsigned int i = 0; i < _nComp; ++i)
		_inletConcentrations[i] = rawValues[i];

	// Get time derivative
	_inlet->timeDerivative(static_cast<double>(t), secIdx, _inletDerivatives + _nComp);
	LOG(Debug) << "dInlet / dt = " << cadet::log::VectorPtr<double>(_inletDerivatives + _nComp, _nComp); 

	// Retrieve parameter derivatives
	for (auto sp : _sensParamsInlet)
	{
		const unsigned int adDir = std::get<0>(sp.second);
		const double adValue = std::get<1>(sp.second);
		const double tAdVal = t.getADValue(adDir);

		_inlet->parameterDerivative(static_cast<double>(t), secIdx, sp.first, _inletDerivatives);
		LOG(Debug) << "dInlet / dp = " << cadet::log::VectorPtr<double>(_inletDerivatives, _nComp); 

		// Copy derivatives into AD datatypes
		for (unsigned int i = 0; i < _nComp; ++i)
			_inletConcentrations[i].setADValue(adDir, _inletDerivatives[i] * adValue + _inletDerivatives[_nComp + i] * tAdVal);

		LOG(Debug) << "totalInlet " << adDir << " " << cadet::log::VectorPtr<cadet::active>(_inletConcentrations, _nComp);
	}

	return _inletConcentrations;
}

void InletModel::useAnalyticJacobian(const bool analyticJac) { }
void InletModel::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx) { }

void InletModel::reportSolution(ISolutionRecorder& recorder, double const* const solution) const
{
	Exporter expr(_nComp, _inletConcentrationsRaw);
	recorder.beginUnitOperation(_unitOpIdx, *this, expr);
	recorder.endUnitOperation();
}

void InletModel::reportSolutionStructure(ISolutionRecorder& recorder) const
{
	Exporter expr(_nComp, nullptr);
	recorder.unitOperationStructure(_unitOpIdx, *this, expr);
}

unsigned int InletModel::requiredADdirs() const CADET_NOEXCEPT
{
	return 0;
}

void InletModel::prepareADvectors(active* const adRes, active* const adY, unsigned int numSensAdDirs) const { }

int InletModel::residual(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double* const res)
{
//	std::cout << "======= INLET RESIDUAL ========== t = " << static_cast<double>(t) << " sec = " << secIdx << " dt = " << static_cast<double>(timeFactor) << std::endl;

	// Evaluate the user-specified function for the inlet concentration
	_inlet->inletConcentration(static_cast<double>(t), secIdx, _inletConcentrationsRaw);

	// Inlet concentrations are now in _inletConcentrationsRaw

	return 0;
}

int InletModel::residualWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, double* const res, active* const adRes, active* const adY, unsigned int numSensAdDirs)
{
	return residual(static_cast<double>(t), secIdx, static_cast<double>(timeFactor), y, yDot, res);
}

double InletModel::residualNorm(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot)
{
	return 0.0;
}

int InletModel::residualSensFwdAdOnly(const active& t, unsigned int secIdx, const active& timeFactor,
	double const* const y, double const* const yDot, active* const adRes)
{
	// Evaluate the user-specified function for the inlet concentration
	_inlet->inletConcentration(static_cast<double>(t), secIdx, _inletConcentrationsRaw);

	// Copy inlet concentrations over to active types if necessary
	moveInletValues<active>(_inletConcentrationsRaw, t, secIdx);

	// Active inlet concentrations are now in _inletConcentrations

	return 0;
}

int InletModel::residualSensFwdCombine(const active& timeFactor, const std::vector<const double*>& yS, const std::vector<const double*>& ySdot,
	const std::vector<double*>& resS, active const* adRes, double* const tmp1, double* const tmp2, double* const tmp3)
{
	return 0;
}

int InletModel::residualSensFwd(unsigned int nSens, const active& t, unsigned int secIdx,
	const active& timeFactor, double const* const y, double const* const yDot, double const* const res,
	const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS,
	active* const adRes, double* const tmp1, double* const tmp2, double* const tmp3)
{
	return 0;
}

int InletModel::residualSensFwdWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, active* const adRes, active* const adY, unsigned int numSensAdDirs)
{
	return residualSensFwdAdOnly(t, secIdx, timeFactor, y, yDot, adRes);
}

void InletModel::residualSensFwdNorm(unsigned int nSens, const active& t, unsigned int secIdx,
		const active& timeFactor, double const* const y, double const* const yDot,
		const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, double* const norms,
		active* const adRes, double* const tmp)
{
	std::fill(norms, norms + nSens, 0.0);
}

active InletModel::inletConnectionFactorActive(unsigned int compIdx, unsigned int secIdx) const CADET_NOEXCEPT
{
	return 0.0;
}

double InletModel::inletConnectionFactor(unsigned int compIdx, unsigned int secIdx) const CADET_NOEXCEPT
{
	return 0.0;
}

}  // namespace model

}  // namespace cadet
