// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//
//  Copyright © 2008-2018: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "model/UnitOperationBase.hpp"
#include "model/BindingModel.hpp"

#include "SensParamUtil.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <iterator>
#include <limits>

namespace cadet
{

namespace model
{

UnitOperationBase::UnitOperationBase(UnitOpIdx unitOpIdx) : _unitOpIdx(unitOpIdx), _binding(nullptr)
{
}

UnitOperationBase::~UnitOperationBase() CADET_NOEXCEPT
{
	delete _binding;
}

std::unordered_map<ParameterId, double> UnitOperationBase::getAllParameterValues() const
{
	std::unordered_map<ParameterId, double> data;
	std::transform(_parameters.begin(), _parameters.end(), std::inserter(data, data.end()),
	               [](const std::pair<const ParameterId, active*>& p) { return std::make_pair(p.first, static_cast<double>(*p.second)); });

	if (!_binding)
		return data;

	const std::unordered_map<ParameterId, double> localData = _binding->getAllParameterValues();
	for (const std::pair<ParameterId, double>& val : localData)
		data[val.first] = val.second;

	return data;
}

double UnitOperationBase::getParameterDouble(const ParameterId& pId) const
{
	// Check our own parameters
	const paramMap_t::const_iterator paramHandle = _parameters.find(pId);
	if (paramHandle != _parameters.end())
		return static_cast<double>(*paramHandle->second);

	// Check binding model parameters
	if (_binding)
	{
		active const* const val = _binding->getParameter(pId);
		if (val)
			return static_cast<double>(*val);
	}

	// Not found
	return std::numeric_limits<double>::quiet_NaN();
}

bool UnitOperationBase::hasParameter(const ParameterId& pId) const
{
	const bool hasParam = _parameters.find(pId) != _parameters.end();
	if (_binding)
		return hasParam || _binding->hasParameter(pId);
	return hasParam;
}

bool UnitOperationBase::setParameter(const ParameterId& pId, int value)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return false;

	if (_binding)
		return _binding->setParameter(pId, value);
	return false;
}

bool UnitOperationBase::setParameter(const ParameterId& pId, double value)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return false;

	paramMap_t::iterator paramHandle = _parameters.find(pId);
	if (paramHandle != _parameters.end())
	{
		paramHandle->second->setValue(value);
		return true;
	}
	else if (_binding)
		return _binding->setParameter(pId, value);

	return false;
}

bool UnitOperationBase::setParameter(const ParameterId& pId, bool value)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return false;

	if (_binding)
		return _binding->setParameter(pId, value);
	return false;
}

void UnitOperationBase::setSensitiveParameterValue(const ParameterId& pId, double value)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return;

	// Check our own parameters
	paramMap_t::iterator paramHandle = _parameters.find(pId);
	if ((paramHandle != _parameters.end()) && contains(_sensParams, paramHandle->second))
	{
		paramHandle->second->setValue(value);
		return;
	}

	// Check binding model parameters
	if (_binding)
	{
		active* const val = _binding->getParameter(pId);
		if (val && contains(_sensParams, val))
		{
			val->setValue(value);
			return;
		}
	}
}

bool UnitOperationBase::setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return false;

	// Check own parameters
	paramMap_t::iterator paramHandle = _parameters.find(pId);
	if (paramHandle != _parameters.end())
	{
		LOG(Debug) << "Found parameter " << pId << " in GRM: Dir " << adDirection << " is set to " << adValue;

		// Register parameter and set AD seed / direction
		_sensParams.insert(paramHandle->second);
		paramHandle->second->setADValue(adDirection, adValue);
		return true;
	}

	// Check binding model parameters
	if (_binding)
	{
		active* const paramBinding = _binding->getParameter(pId);
		if (paramBinding)
		{
			LOG(Debug) << "Found parameter " << pId << " in AdsorptionModel: Dir " << adDirection << " is set to " << adValue;

			// Register parameter and set AD seed / direction
			_sensParams.insert(paramBinding);
			paramBinding->setADValue(adDirection, adValue);
			return true;
		}
	}

	return false;
}

void UnitOperationBase::clearSensParams()
{
	// Remove AD directions from parameters
	for (auto sp : _sensParams)
		sp->setADValue(0.0);

	_sensParams.clear();
}

int UnitOperationBase::residualSensFwdCombine(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, 
	const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS, active const* adRes, 
	double* const tmp1, double* const tmp2, double* const tmp3)
{
	for (unsigned int param = 0; param < yS.size(); param++)
	{
		// tmp1 stores result of (dF / dy) * s
		// tmp2 stores result of (dF / dyDot) * sDot

		// Directional derivative (dF / dy) * s
		multiplyWithJacobian(static_cast<double>(t), secIdx, static_cast<double>(timeFactor), y, yDot, yS[param], 1.0, 0.0, tmp1);

		// Directional derivative (dF / dyDot) * sDot
		multiplyWithDerivativeJacobian(static_cast<double>(t), secIdx, static_cast<double>(timeFactor), y, yDot, ySdot[param], tmp2);

		// Complete sens residual is the sum:
		double* const ptrResS = resS[param];
		for (unsigned int i = 0; i < numDofs(); ++i)
		{
			ptrResS[i] = tmp1[i] + tmp2[i] + adRes[i].getADValue(param);
		}
	}
	return 0;
}

}  // namespace model

}  // namespace cadet
