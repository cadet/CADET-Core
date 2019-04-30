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

#include "model/UnitOperationBase.hpp"
#include "model/BindingModel.hpp"
#include "model/ReactionModel.hpp"
#include "SimulationTypes.hpp"

#include "SensParamUtil.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <iterator>
#include <limits>

namespace cadet
{

namespace model
{

namespace
{
	template <typename T>
	void getAllParameterValuesImpl(std::unordered_map<ParameterId, double>& data, const std::vector<T*>& items, bool singleItem)
	{
		if (!items.empty())
		{
			if (singleItem)
			{
				const std::unordered_map<ParameterId, double> localData = items[0]->getAllParameterValues();
				for (const std::pair<ParameterId, double>& val : localData)
					data[val.first] = val.second;
			}
			else
			{
				for (T const* bm : items)
				{
					const std::unordered_map<ParameterId, double> localData = bm->getAllParameterValues();
					for (const std::pair<ParameterId, double>& val : localData)
						data[val.first] = val.second;
				}
			}
		}		
	}

	template <typename T>
	bool getParameterDoubleImpl(const ParameterId& pId, const std::vector<T*>& items, bool singleItem, double& out)
	{
		// Check binding model parameters
		if (!items.empty())
		{
			if (singleItem)
			{
				active const* const val = items[0]->getParameter(pId);
				if (val)
				{
					out = static_cast<double>(*val);
					return true;
				}
			}
			else
			{
				for (T* bm : items)
				{
					active const* const val = bm->getParameter(pId);
					if (val)
					{
						out = static_cast<double>(*val);
						return true;
					}
				}
			}
		}

		// Not found
		return false;
	}

	template <typename T>
	bool hasParameterImpl(const ParameterId& pId, const std::vector<T*>& items, bool singleItem)
	{
		if (!items.empty())
		{
			if (singleItem)
			{
				if (items[0]->hasParameter(pId))
					return true;
			}
			else
			{
				for (T const* bm : items)
				{
					if (bm->hasParameter(pId))
						return true;
				}
			}
		}

		return false;
	}

	template <typename T, typename param_t>
	bool setParameterImpl(const ParameterId& pId, param_t value, const std::vector<T*>& items, bool singleItem)
	{
		if (!items.empty())
		{
			if (singleItem)
			{
				if (items[0]->setParameter(pId, value))
					return true;
			}
			else
			{
				for (T* bm : items)
				{
					if (bm->setParameter(pId, value))
						return true;
				}
			}
		}

		return false;
	}

	template <typename T>
	bool setSensitiveParameterValueImpl(const ParameterId& pId, double value, const std::unordered_set<active*>& sensParams, const std::vector<T*>& items, bool singleItem)
	{
		if (!items.empty())
		{
			if (singleItem)
			{
				active* const val = items[0]->getParameter(pId);
				if (val && contains(sensParams, val))
				{
					val->setValue(value);
					return true;
				}
			}
			else
			{
				for (T* bm : items)
				{
					active* const val = bm->getParameter(pId);
					if (val && contains(sensParams, val))
					{
						val->setValue(value);
						return true;
					}
				}
			}
		}

		return false;
	}

	template <typename T>
	bool setSensitiveParameterImpl(const ParameterId& pId, unsigned int adDirection, double adValue, std::unordered_set<active*>& sensParams, const std::vector<T*>& items, bool singleItem)
	{
		if (!items.empty())
		{
			if (singleItem)
			{
				active* const paramBinding = items[0]->getParameter(pId);
				if (paramBinding)
				{
					// Register parameter and set AD seed / direction
					sensParams.insert(paramBinding);
					paramBinding->setADValue(adDirection, adValue);
					return true;
				}
			}
			else
			{
				for (T* bm : items)
				{
					active* const paramBinding = bm->getParameter(pId);
					if (paramBinding)
					{
						// Register parameter and set AD seed / direction
						sensParams.insert(paramBinding);
						paramBinding->setADValue(adDirection, adValue);
						return true;
					}
				}
			}
		}

		return false;
	}
}

UnitOperationBase::UnitOperationBase(UnitOpIdx unitOpIdx) : _unitOpIdx(unitOpIdx), _binding(0, nullptr), _singleBinding(false),
	_dynReaction(0, nullptr), _singleDynReaction(false)
{
}

UnitOperationBase::~UnitOperationBase() CADET_NOEXCEPT
{
	clearBindingModels();
	clearDynamicReactionModels();
}

void UnitOperationBase::clearBindingModels() CADET_NOEXCEPT
{
	if (_singleBinding)
	{
		if (!_binding.empty())
			delete _binding[0];
	}
	else
	{
		for (IBindingModel* bm : _binding)
			delete bm;
	}

	_binding.clear();
}

void UnitOperationBase::clearDynamicReactionModels() CADET_NOEXCEPT
{
	if (_singleDynReaction)
	{
		if (!_dynReaction.empty())
			delete _dynReaction[0];
	}
	else
	{
		for (IDynamicReactionModel* drm : _dynReaction)
			delete drm;
	}

	_dynReaction.clear();
}

std::unordered_map<ParameterId, double> UnitOperationBase::getAllParameterValues() const
{
	std::unordered_map<ParameterId, double> data;
	std::transform(_parameters.begin(), _parameters.end(), std::inserter(data, data.end()),
	               [](const std::pair<const ParameterId, active*>& p) { return std::make_pair(p.first, static_cast<double>(*p.second)); });

	getAllParameterValuesImpl(data, _binding, _singleBinding);
	getAllParameterValuesImpl(data, _dynReaction, _singleDynReaction);

	return data;
}

double UnitOperationBase::getParameterDouble(const ParameterId& pId) const
{
	// Check our own parameters
	const paramMap_t::const_iterator paramHandle = _parameters.find(pId);
	if (paramHandle != _parameters.end())
		return static_cast<double>(*paramHandle->second);

	double val = std::numeric_limits<double>::quiet_NaN();
	if (getParameterDoubleImpl(pId, _binding, _singleBinding, val))
		return val;
	if (getParameterDoubleImpl(pId, _dynReaction, _singleDynReaction, val))
		return val;

	// Not found
	return std::numeric_limits<double>::quiet_NaN();
}

bool UnitOperationBase::hasParameter(const ParameterId& pId) const
{
	if (_parameters.find(pId) != _parameters.end())
		return true;

	if (hasParameterImpl(pId, _binding, _singleBinding))
		return true;
	if (hasParameterImpl(pId, _dynReaction, _singleDynReaction))
		return true;
	
	return false;
}

bool UnitOperationBase::setParameter(const ParameterId& pId, int value)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return false;

	if (setParameterImpl(pId, value, _binding, _singleBinding))
		return true;
	if (setParameterImpl(pId, value, _dynReaction, _singleDynReaction))
		return true;

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

	if (setParameterImpl(pId, value, _binding, _singleBinding))
		return true;
	if (setParameterImpl(pId, value, _dynReaction, _singleDynReaction))
		return true;

	return false;
}

bool UnitOperationBase::setParameter(const ParameterId& pId, bool value)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return false;

	if (setParameterImpl(pId, value, _binding, _singleBinding))
		return true;
	if (setParameterImpl(pId, value, _dynReaction, _singleDynReaction))
		return true;

	return false;
}

void UnitOperationBase::setSensitiveParameterValue(const ParameterId& pId, double value)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return;

	// Check our own model parameters
	paramMap_t::iterator paramHandle = _parameters.find(pId);
	if ((paramHandle != _parameters.end()) && contains(_sensParams, paramHandle->second))
	{
		paramHandle->second->setValue(value);
		return;
	}

	if (setSensitiveParameterValueImpl(pId, value, _sensParams, _binding, _singleBinding))
		return;
	if (setSensitiveParameterValueImpl(pId, value, _sensParams, _dynReaction, _singleDynReaction))
		return;
}

bool UnitOperationBase::setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return false;

	// Check own model parameters
	paramMap_t::iterator paramHandle = _parameters.find(pId);
	if (paramHandle != _parameters.end())
	{
		LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;

		// Register parameter and set AD seed / direction
		_sensParams.insert(paramHandle->second);
		paramHandle->second->setADValue(adDirection, adValue);
		return true;
	}

	if (setSensitiveParameterImpl(pId, adDirection, adValue, _sensParams, _binding, _singleBinding))
	{
		LOG(Debug) << "Found parameter " << pId << " in AdsorptionModel: Dir " << adDirection << " is set to " << adValue;
		return true;
	}
	if (setSensitiveParameterImpl(pId, adDirection, adValue, _sensParams, _dynReaction, _singleDynReaction))
	{
		LOG(Debug) << "Found parameter " << pId << " in DynamicReactionModel: Dir " << adDirection << " is set to " << adValue;
		return true;
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

unsigned int UnitOperationBase::numSensParams() const
{
	return _sensParams.size();
}

int UnitOperationBase::residualSensFwdCombine(const ActiveSimulationTime& simTime, const ConstSimulationState& simState,
	const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS, active const* adRes,
	double* const tmp1, double* const tmp2, double* const tmp3)
{
	for (unsigned int param = 0; param < yS.size(); ++param)
	{
		// tmp1 stores result of (dF / dy) * s
		// tmp2 stores result of (dF / dyDot) * sDot

		// Directional derivative (dF / dy) * s
		multiplyWithJacobian(toSimple(simTime), simState, yS[param], 1.0, 0.0, tmp1);

		// Directional derivative (dF / dyDot) * sDot
		multiplyWithDerivativeJacobian(toSimple(simTime), simState, ySdot[param], tmp2);

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
