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

#include "model/ModelSystemImpl.hpp"
#include "ParamIdUtil.hpp"
#include "cadet/ParameterProvider.hpp"
#include "cadet/Exceptions.hpp"
#include "cadet/ExternalFunction.hpp"
#include "ConfigurationHelper.hpp"
#include "linalg/Norms.hpp"

#include <algorithm>
#include <string>
#include <iomanip>
#include <sstream>

#include "LoggingUtils.hpp"
#include "Logging.hpp"


#include "OutletModel.hpp"
#include "GeneralRateModel.hpp"

namespace
{
	template <typename ParamType>
	struct InOutFactorProxy { };

	template <>
	struct InOutFactorProxy<cadet::active>
	{
		static const cadet::active inletFactor(cadet::IUnitOperation const* const model, unsigned int compIdx, unsigned int secIdx)
		{
			return model->inletConnectionFactorActive(compIdx, secIdx);
		}

		static cadet::active const* const data(cadet::IUnitOperation const* const model) { return model->getDataActive(); }
	};

	template <>
	struct InOutFactorProxy<double>
	{
		static const double inletFactor(cadet::IUnitOperation const* const model, unsigned int compIdx, unsigned int secIdx)
		{
			return model->inletConnectionFactor(compIdx, secIdx);
		}

		static double const* const data(cadet::IUnitOperation const* const model) { return model->getData(); }
	};

	unsigned int indexOfUnitOp(const std::vector<cadet::IUnitOperation*>& models, unsigned int unitOpIdx)
	{
		for (unsigned int i = 0; i < models.size(); ++i)
		{
			if (models[i]->unitOperationId() == unitOpIdx)
				return i;
		}
		return models.size();
	}
}

namespace cadet
{

namespace model
{

ModelSystem::ModelSystem() : _curSwitchIndex(0)
{
}

ModelSystem::~ModelSystem() CADET_NOEXCEPT
{
	for (IUnitOperation* model : _models)
		delete model;

	for (IExternalFunction* extFun : _extFunctions)
		delete extFun;
}

void ModelSystem::addModel(IModel* unitOp)
{
	// Check for unique unit operation id
	if (indexOfUnitOp(_models, unitOp->unitOperationId()) < _models.size())
		throw InvalidParameterException("Cannot add model because of already existing unit operation id " + std::to_string(unitOp->unitOperationId()));

	_models.push_back(static_cast<IUnitOperation*>(unitOp));

	// Propagate external functions to submodel
	_models.back()->setExternalFunctions(_extFunctions.data(), _extFunctions.size());
}

IModel* ModelSystem::getModel(unsigned int index)
{
	if (index < _models.size())
		return _models[index];
	else
		return nullptr;
}

IModel const* ModelSystem::getModel(unsigned int index) const
{
	if (index < _models.size())
		return _models[index];
	else
		return nullptr;
}

IModel* ModelSystem::getUnitOperationModel(unsigned int unitOpIdx)
{
	for (IUnitOperation* m : _models)
	{
		if (m->unitOperationId() == unitOpIdx)
			return m;
	}
	return nullptr;
}

IModel const* ModelSystem::getUnitOperationModel(unsigned int unitOpIdx) const
{
	for (IUnitOperation* m : _models)
	{
		if (m->unitOperationId() == unitOpIdx)
			return m;
	}
	return nullptr;
}

unsigned int ModelSystem::numModels() const CADET_NOEXCEPT
{
	return _models.size();
}

void ModelSystem::removeModel(IModel const* unitOp)
{
	std::vector<IUnitOperation*>::iterator it = std::find(_models.begin(), _models.end(), static_cast<IUnitOperation const*>(unitOp));
	if (it != _models.end())
		_models.erase(it);
}

IModel* ModelSystem::removeModel(UnitOpIdx unitOp)
{
	for (std::vector<IUnitOperation*>::iterator it = _models.begin(); it != _models.end(); ++it)
	{
		if ((*it)->unitOperationId() == unitOp)
		{
			IModel* const mod = static_cast<IModel*>(*it);
			_models.erase(it);
			return mod;
		}
	}
	return nullptr;
}

UnitOpIdx ModelSystem::maxUnitOperationId() const CADET_NOEXCEPT
{
	if (_models.empty())
		return UnitOpIndep;

	UnitOpIdx id = 0;
	for (IUnitOperation* m : _models)
		id = std::max(id, m->unitOperationId());
	return id;
}

unsigned int ModelSystem::addExternalFunction(IExternalFunction& extFun)
{
	_extFunctions.push_back(&extFun);

	// Propagate external functions to submodels
	for (IUnitOperation* m : _models)
		m->setExternalFunctions(_extFunctions.data(), _extFunctions.size());

	return _extFunctions.size() - 1;
}

IExternalFunction* ModelSystem::getExternalFunction(unsigned int index)
{
	if (index < _extFunctions.size())
		return _extFunctions[index];
	else
		return nullptr;
}

IExternalFunction const* ModelSystem::getExternalFunction(unsigned int index) const
{
	if (index < _extFunctions.size())
		return _extFunctions[index];
	else
		return nullptr;
}

unsigned int ModelSystem::numExternalFunctions() const CADET_NOEXCEPT
{
	return _extFunctions.size();
}

void ModelSystem::removeExternalFunction(IExternalFunction const* extFun)
{
	for (std::vector<IExternalFunction*>::iterator it = _extFunctions.begin(); it != _extFunctions.end(); ++it)
	{
		if (*it == extFun)
		{
			_extFunctions.erase(it);
			break;
		}
	}

	// Update external functions in submodels
	for (IUnitOperation* m : _models)
		m->setExternalFunctions(_extFunctions.data(), _extFunctions.size());
}

unsigned int ModelSystem::numDofs() const CADET_NOEXCEPT
{
	unsigned int dofs = numPureDofs();
	return dofs;
}

unsigned int ModelSystem::numPureDofs() const CADET_NOEXCEPT
{
	unsigned int dofs = 0;
	for (IUnitOperation* m : _models)
		dofs += m->numDofs();
	
	return dofs;
}

bool ModelSystem::usesAD() const CADET_NOEXCEPT
{
	for (IUnitOperation* m : _models)
	{
		if (m->usesAD())
			return true;
	}
	return false;
}

void ModelSystem::rebuildInternalDataStructures()
{
	// Sort models by unit operation Id

	// Calculate array with DOF offsets
	_dofOffset.clear();
	_dofOffset.reserve(_models.size());

	unsigned int totalDof = 0;
	for (IUnitOperation* m : _models)
	{
		_dofOffset.push_back(totalDof);
		totalDof += m->numDofs();
	}

	LOG(Debug) << "DOF offsets: " << _dofOffset;
}

bool ModelSystem::configure(IParameterProvider& paramProvider, IConfigHelper& helper)
{
	// Unit operation models are already configured
	rebuildInternalDataStructures();

	configureSwitches(paramProvider);
	_curSwitchIndex = 0;

	// Create and configure all external functions
	bool success = true;
	if (paramProvider.exists("external"))
	{
		paramProvider.pushScope("external");

		unsigned int i = 0;
		std::ostringstream oss;
		oss << "source_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << i;
		while (paramProvider.exists(oss.str()))
		{
			// Create and configure unit operation
			paramProvider.pushScope(oss.str());

			const std::string extType = paramProvider.getString("EXTFUN_TYPE");
			IExternalFunction* const func = helper.createExternalFunction(extType);
			if (func)
			{
				const bool confSuccess = func->configure(&paramProvider);

				if (confSuccess)
					_extFunctions.push_back(func);
				else
				{
					// Ignore the external function and delete this instance
					_extFunctions.push_back(nullptr);
					delete func;
					success = false;

					LOG(Error) << "Failed to configure external source " << i << " (" << extType << "), source is ignored";
				}
			}
			else
			{
				// Unknown type of external function
				_extFunctions.push_back(nullptr);
				success = false;

				LOG(Error) << "Failed to create external source " << i << " as type " << extType << " is unknown, source is ignored";
			}

			paramProvider.popScope();

			// Next group in file format
			++i;
			oss.str("");
			oss << "source_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << i;
		}

		paramProvider.popScope();
	}

	// Propagate external functions to submodels
	for (IUnitOperation* m : _models)
		m->setExternalFunctions(_extFunctions.data(), _extFunctions.size());

	return success;
}

bool ModelSystem::reconfigure(IParameterProvider& paramProvider)
{
	configureSwitches(paramProvider);

	// Reconfigure all external functions
	bool success = true;
	if (paramProvider.exists("external"))
	{
		paramProvider.pushScope("external");

		std::ostringstream oss;
		for (unsigned int i = 0; i < _extFunctions.size(); ++i)
		{
			IExternalFunction* const func = _extFunctions[i];
			if (!func)
				continue;

			oss.str("");
			oss << "source_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << i;
			if (!paramProvider.exists(oss.str()))
				continue;

			paramProvider.pushScope(oss.str());
			const bool localSuccess = func->configure(&paramProvider);
			paramProvider.popScope();

			if (!localSuccess)
			{
				LOG(Error) << "Failed to reconfigure external source " << i;
			}

			success = localSuccess && success;
		}

		paramProvider.popScope();
	}

	return success;
}

/**
 * @brief Reads valve switches from the given parameter provider
 * @param [in] paramProvider Parameter provider
 */
void ModelSystem::configureSwitches(IParameterProvider& paramProvider)
{
	// Read connections of unit operations
	paramProvider.pushScope("connections");

	const unsigned int numSwitches = paramProvider.getInt("NSWITCHES");
	_switchSectionIndex.clear();
	_switchSectionIndex.reserve(numSwitches);
	_connections.clear();
	_connections.reserve(numSwitches * 4 * _models.size() * _models.size(), numSwitches);

	std::ostringstream oss;
	for (unsigned int i = 0; i < numSwitches; ++i)
	{
		oss.str("");
		oss << "switch_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << i;

		paramProvider.pushScope(oss.str());

		_switchSectionIndex.push_back(paramProvider.getInt("SECTION"));
		std::vector<int> conn = paramProvider.getIntArray("CONNECTIONS");
		if ((conn.size() % 4) != 0)
			throw InvalidParameterException("CONNECTIONS matrix has to have 4 columns");

		checkConnectionList(conn);

		_connections.pushBackSlice(conn);
		paramProvider.popScope();
	}

	paramProvider.popScope();

	if (_switchSectionIndex[0] != 0)
		throw InvalidParameterException("First element of SECTION in connections group has to be 0");

	// TODO: Sanity check all connections (existence of inlets and outlets)
}

bool ModelSystem::reconfigureModel(IParameterProvider& paramProvider, unsigned int unitOpIdx)
{
	IUnitOperation* const model = static_cast<IUnitOperation*>(getUnitOperationModel(unitOpIdx));
	if (!model)
		return false;

	return model->reconfigure(paramProvider);
}

/**
 * @brief Checks the given unit operation connection list and reformats it
 * @details Throws an exception if something is incorrect. Reformats the connection list by
 *          substituting unit operation IDs with local indices.
 * @param [in,out] conn Matrix with 4 columns holding all connections. The matrix is expected
 *                      to be in row-major storage format. On exit, the unit operation IDs are
 *                      substituted by the corresponding indices of the unit operations in the
 *                      local _models vector.
 */
void ModelSystem::checkConnectionList(std::vector<int>& conn) const
{
	for (unsigned int i = 0; i < conn.size() / 4; ++i)
	{
		// Extract current connection
		int& uoSource = conn[4*i];
		int& uoDest = conn[4*i+1];
		const int compSource = conn[4*i+2];
		const int compDest = conn[4*i+3];

		if (uoSource < 0)
			throw InvalidParameterException("In CONNECTIONS matrix (row " + std::to_string(i) + "): Source unit operation id has to be at least 0 in connection");
		if (uoDest < 0)
			throw InvalidParameterException("In CONNECTIONS matrix (row " + std::to_string(i) + "): Destination unit operation id has to be at least 0 in connection");

		// Convert to index
		uoSource = indexOfUnitOp(_models, uoSource);
		uoDest = indexOfUnitOp(_models, uoDest);

		if (static_cast<unsigned int>(uoSource) >= _models.size())
			throw InvalidParameterException("In CONNECTIONS matrix (row " + std::to_string(i) + "): Source unit operation id not found in connection");
		if (static_cast<unsigned int>(uoDest) >= _models.size())
			throw InvalidParameterException("In CONNECTIONS matrix (row " + std::to_string(i) + "): Destination unit operation id not found in connection");

		// Check component indices
		if ((compSource >= 0) && (static_cast<unsigned int>(compSource) >= _models[uoSource]->numComponents()))
			throw InvalidParameterException("In CONNECTIONS matrix (row " + std::to_string(i) + "): Source component index exceeds number of components " + std::to_string(_models[uoSource]->numComponents()));
		if ((compDest >= 0) && (static_cast<unsigned int>(compDest) >= _models[uoDest]->numComponents()))
			throw InvalidParameterException("In CONNECTIONS matrix (row " + std::to_string(i) + "): Destination component index exceeds number of components " + std::to_string(_models[uoDest]->numComponents()));

		if (((compSource < 0) && (compDest >= 0)) || ((compSource >= 0) && (compDest < 0)))
			throw InvalidParameterException("In CONNECTIONS matrix (row " + std::to_string(i) + "): Only source or destination (not both) are set to connect all components in connection from unit " 
				+ std::to_string(_models[uoSource]->unitOperationId()) + " to " + std::to_string(_models[uoDest]->unitOperationId()));

		if ((compSource < 0) && (compDest < 0) && (_models[uoSource]->numComponents() != _models[uoDest]->numComponents()))
			throw InvalidParameterException("In CONNECTIONS matrix (row " + std::to_string(i) + "): Number of components not equal when connecting all components from unit " 
				+ std::to_string(_models[uoSource]->unitOperationId()) + " to " + std::to_string(_models[uoDest]->unitOperationId()));
	}

	// TODO: Check for conflicting entries
	// TODO: Plausibility check of total connections
}

std::unordered_map<ParameterId, double> ModelSystem::getAllParameterValues() const
{
	std::unordered_map<ParameterId, double> data;

	for (IUnitOperation* m : _models)
	{
		const std::unordered_map<ParameterId, double> localData = m->getAllParameterValues();
		for (const std::pair<ParameterId, double>& val : localData)
			data[val.first] = val.second;
	}

	return data;
}

bool ModelSystem::hasParameter(const ParameterId& pId) const
{
	for (IUnitOperation* m : _models)
	{
		if ((m->unitOperationId() == pId.unitOperation) || (pId.unitOperation == UnitOpIndep))
		{
			if (m->hasParameter(pId))
				return true;
		}
	}
	return false;
}

template <typename ParamType>
bool ModelSystem::setParameterImpl(const ParameterId& pId, const ParamType value)
{
	// Filter by unit operation ID
	for (IUnitOperation* m : _models)
	{
		if ((m->unitOperationId() == pId.unitOperation) || (pId.unitOperation == UnitOpIndep))
			return m->setParameter(pId, value);
	}
	return false;
}

bool ModelSystem::setParameter(const ParameterId& pId, int value)
{
	return setParameterImpl(pId, value);
}

bool ModelSystem::setParameter(const ParameterId& pId, double value)
{
	return setParameterImpl(pId, value);
}

bool ModelSystem::setParameter(const ParameterId& pId, bool value)
{
	return setParameterImpl(pId, value);
}

void ModelSystem::setSensitiveParameterValue(const ParameterId& pId, double value)
{
	// Filter by unit operation ID
	for (IUnitOperation* m : _models)
	{
		if ((m->unitOperationId() == pId.unitOperation) || (pId.unitOperation == UnitOpIndep))
		{
			m->setSensitiveParameterValue(pId, value);
		}
	}
}

bool ModelSystem::setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue)
{
	bool found = false;

	// Filter by unit operation ID
	for (IUnitOperation* m : _models)
	{
		if ((m->unitOperationId() == pId.unitOperation) || (pId.unitOperation == UnitOpIndep))
		{
			found = m->setSensitiveParameter(pId, adDirection, adValue) || found;
		}
	}
	return found;
}

void ModelSystem::clearSensParams()
{
	for (IUnitOperation* m : _models)
		m->clearSensParams();
}

void ModelSystem::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx)
{
	for (IUnitOperation* m : _models)
		m->notifyDiscontinuousSectionTransition(t, secIdx);

	// Check if simulation is (re-)starting from the very beginning
	if (secIdx == 0)
		_curSwitchIndex = 0;

	const unsigned int prevSwitch = _curSwitchIndex;

	// If there are still some switches left and the next switch occurrs in this section, advance index
	if ((_curSwitchIndex < _switchSectionIndex.size() - 1) && (_switchSectionIndex[_curSwitchIndex + 1] >= secIdx))
		++_curSwitchIndex;


	LOG(Debug) << "Valve switched from connection " << prevSwitch << " to " << _curSwitchIndex;
	int const* ptrConn = _connections[_curSwitchIndex];
	for (unsigned int i = 0; i < _connections.sliceSize(_curSwitchIndex) / 4; ++i, ptrConn += 4)
	{
		LOG(Debug) << "Unit op " << ptrConn[0] << " (" << _models[ptrConn[0]]->unitOperationName() << ") comp " << ptrConn[2] << " => " 
		           << ptrConn[1] << " (" << _models[ptrConn[1]]->unitOperationName() << ") comp " << ptrConn[3];
	}

	// TODO: Assemble outer Jacobian submatrices
}

void ModelSystem::reportSolution(ISolutionRecorder& recorder, double const* const solution) const
{
	// TODO: Adjust indexing / offset of solution vector
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		if (m->hasInlet() && !m->hasOutlet())
		{
			// Detected outlet model
			// Seach for GRM and hand it over

			for (unsigned int j = 0; j < _models.size(); ++j)
			{
				IUnitOperation* const m2 = _models[j];
				if (m2->hasInlet() && m2->hasOutlet())
				{
					static_cast<OutletModel* const>(m)->reportSolution(recorder, solution + _dofOffset[j], *static_cast<GeneralRateModel* const>(m2));
					break;
				}
			}
		}
		else
		{
			m->reportSolution(recorder, solution + _dofOffset[i]);
		}	
	}
}

void ModelSystem::reportSolutionStructure(ISolutionRecorder& recorder) const
{
	for (IUnitOperation* m : _models)
		m->reportSolutionStructure(recorder);
}

unsigned int ModelSystem::requiredADdirs() const CADET_NOEXCEPT
{
	// TODO: Apply compression here
	unsigned int dirs = 0;
	for (IUnitOperation* m : _models)
		dirs += m->requiredADdirs();
	
	return dirs;
}

void ModelSystem::prepareADvectors(active* const adRes, active* const adY, unsigned int numSensAdDirs) const
{
	// TODO: Adjust indexing / offset of AD vectors
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		const unsigned int offset = _dofOffset[i];
		_models[i]->prepareADvectors(adRes + offset, adY + offset, numSensAdDirs);
	}
}

double ModelSystem::residualNorm(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot)
{
	std::vector<double> temp(numDofs(), 0.0);
	residual(t, secIdx, timeFactor, y, yDot, temp.data());
	return linalg::linfNorm(temp.data(), temp.size());
}

int ModelSystem::residual(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double* const res)
{
	BENCH_START(_timerResidual);

	int result = 0;

	// TODO: Adjust indexing / offset of vectors
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];
		const int intermediateRes = m->residual(t, secIdx, timeFactor, y + offset, yDot + offset, res + offset);

		// If result is already -1 (non-recoverable error), then we stick to it
		// If result is ok or recoverable and intermediate result is recoverable, then we take intermediate result
		if ((result >= 0) && (intermediateRes > 0))
		{
			result = intermediateRes;
		}
	}

	// Handle connections
	residualConnectUnitOps<double, double, double>(secIdx, y, yDot, res);

	BENCH_STOP(_timerResidual);
	return result;
}

int ModelSystem::residualWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, double* const res, 
	active* const adRes, active* const adY, unsigned int numSensAdDirs)
{
	BENCH_START(_timerResidual);

	int result = 0;

	// TODO: Adjust indexing / offset of vectors
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];
		const int intermediateRes = m->residualWithJacobian(t, secIdx, timeFactor, y + offset, yDot + offset, res + offset, adRes + offset, adY + offset, numSensAdDirs);

		// If result is already -1 (non-recoverable error), then we stick to it
		// If result is ok or recoverable and intermediate result is recoverable, then we take intermediate result
		if ((result >= 0) && (intermediateRes > 0))
		{
			result = intermediateRes;
		}
	}

	// Handle connections
	residualConnectUnitOps<double, double, double>(secIdx, y, yDot, res);

	BENCH_STOP(_timerResidual);
	return result;
}

template <typename StateType, typename ResidualType, typename ParamType>
void ModelSystem::residualConnectUnitOps(unsigned int secIdx, StateType const* const y, double const* const yDot, ResidualType* const res)
{
	// Current connection list
	int const* conList = _connections[_curSwitchIndex];
	const unsigned int numCon = _connections.sliceSize(_curSwitchIndex) / 4;
	for (unsigned int i = 0; i < numCon; ++i, conList += 4)
	{
		IUnitOperation const* const fromModel = _models[conList[0]];
		const unsigned int fromStride = fromModel->localOutletComponentStride();
		const unsigned int fromBegin = fromModel->localOutletComponentIndex();
		ResidualType const* const fromData = InOutFactorProxy<ResidualType>::data(fromModel);

		IUnitOperation const* const toModel = _models[conList[1]];
		const unsigned int toStride = toModel->localInletComponentStride();
		const unsigned int toBegin = toModel->localInletComponentIndex();

		if (conList[2] < 0)
		{
			// Connect all components
			for (unsigned int comp = 0; comp < fromModel->numComponents(); ++comp)
			{
				const ParamType inFactor = InOutFactorProxy<ParamType>::inletFactor(toModel, comp, secIdx);
				if (fromModel->numDofs() == 0)
				{
					res[toBegin + comp * toStride] += inFactor * fromData[comp];
				}
			}
		}
		else
		{
			// Connect only specific components
			const ParamType inFactor = InOutFactorProxy<ParamType>::inletFactor(toModel, static_cast<unsigned int>(conList[3]), secIdx);
			if (fromModel->numDofs() == 0)
			{
				res[toBegin + static_cast<unsigned int>(conList[3]) * toStride] += inFactor * fromData[conList[2]];
			}
		}
	}
}

int ModelSystem::residualSensFwd(unsigned int nSens, const active& t, unsigned int secIdx,
	const active& timeFactor, double const* const y, double const* const yDot, double const* const res,
	const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS,
	active* const adRes, double* const tmp1, double* const tmp2, double* const tmp3)
{
	BENCH_START(_timerResidualSens);

	int result = 0;

	// Step 1: Calculate sensitivities using AD in vector mode

	// TODO: Adjust indexing / offset of vectors
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];
		const int intermediateRes = m->residualSensFwdAdOnly(t, secIdx, timeFactor, y + offset, yDot + offset, adRes + offset);

		// If result is already -1 (non-recoverable error), then we stick to it
		// If result is ok or recoverable and intermediate result is recoverable, then we take intermediate result
		if ((result >= 0) && (intermediateRes > 0))
		{
			result = intermediateRes;
		}
	}

	// Connect units
	residualConnectUnitOps<double, active, active>(secIdx, y, yDot, adRes);

	// Step 2: Compute forward sensitivity residuals by multiplying with system Jacobians
	std::vector<const double*> ySlocal(yS.size(), nullptr);
	std::vector<const double*> ySdotLocal(ySdot.size(), nullptr);
	std::vector<double*> resSlocal(resS.size(), nullptr);
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];

		// Use correct offset in sensitivity state vectors
		for (unsigned int j = 0; j < yS.size(); ++j)
		{
			ySlocal[j] = yS[j] + offset;
			ySdotLocal[j] = ySdot[j] + offset;
			resSlocal[j] = resS[j] + offset;
		}

		const int intermediateRes = m->residualSensFwdCombine(timeFactor, ySlocal, ySdotLocal, resSlocal, adRes + offset, tmp1 + offset, tmp2 + offset, tmp3 + offset);

		// If result is already -1 (non-recoverable error), then we stick to it
		// If result is ok or recoverable and intermediate result is recoverable, then we take intermediate result
		if ((result >= 0) && (intermediateRes > 0))
		{
			result = intermediateRes;
		}
	}

	BENCH_STOP(_timerResidualSens);
	return result;
}

void ModelSystem::residualSensFwdNorm(unsigned int nSens, const active& t, unsigned int secIdx,
		const active& timeFactor, double const* const y, double const* const yDot,
		const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, double* const norms,
		active* const adRes, double* const tmp)
{
	const unsigned int nDOFs = numDofs();

	// Reserve memory for nSens residual vectors
	util::SlicedVector<double> tempRes;
	tempRes.reserve(nSens * nDOFs, nSens);

	std::vector<double*> resPtr(nSens, nullptr);
	for (unsigned int i = 0; i < resPtr.size(); ++i)
	{
		tempRes.pushBackSlice(nDOFs);
		resPtr[i] = tempRes[i];
	}

	// Reserve some more temporary memory
	std::vector<double> tempMem(nDOFs * 2, 0.0);

	// Evaluate all the sensitivity system residuals at once
	residualSensFwd(nSens, t, secIdx, timeFactor, y, yDot, nullptr, yS, ySdot, resPtr, adRes, tmp, tempMem.data(), tempMem.data() + nDOFs);

	// Calculate norms
	for (unsigned int i = 0; i < nSens; ++i)
		norms[i] = linalg::linfNorm(tempRes[i], nDOFs);
}

void ModelSystem::applyInitialCondition(double* const vecStateY, double* const vecStateYdot)
{
	// TODO: Adjust indexing / offset of vectors
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];
		m->applyInitialCondition(vecStateY + offset, vecStateYdot + offset);
	}
}

void ModelSystem::applyInitialCondition(IParameterProvider& paramProvider, double* const vecStateY, double* const vecStateYdot)
{
	bool skipModels = false;

	// Check if INIT_STATE_Y is present
	if (paramProvider.exists("INIT_STATE_Y"))
	{
		const std::vector<double> initState = paramProvider.getDoubleArray("INIT_STATE_Y");
		if (initState.size() >= numDofs())
		{
			std::copy(initState.data(), initState.data() + numDofs(), vecStateY);
			skipModels = true;
		}
	}

	// Check if INIT_STATE_YDOT is present
	if (paramProvider.exists("INIT_STATE_YDOT"))
	{
		const std::vector<double> initState = paramProvider.getDoubleArray("INIT_STATE_YDOT");
		if (initState.size() >= numDofs())
		{
			std::copy(initState.data(), initState.data() + numDofs(), vecStateYdot);
		}
	}

	if (skipModels)
		return;

	// TODO: Adjust indexing / offset of vectors
	std::ostringstream oss;
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];

		oss.str("");
		oss << "unit_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << static_cast<int>(m->unitOperationId());

		paramProvider.pushScope(oss.str());
		m->applyInitialCondition(paramProvider, vecStateY + offset, vecStateYdot + offset);
		paramProvider.popScope();
	}
}

void ModelSystem::consistentInitialConditions(double t, unsigned int secIdx, double timeFactor, double* const vecStateY, double* const vecStateYdot, 
	active* const adRes, active* const adY, unsigned int numSensAdDirs, double errorTol)
{
	BENCH_START(_timerConsistentInit);

	// TODO: Adjust indexing / offset of vectors

	// Phase 1: Solve algebraic equations and update state
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];
		m->consistentInitialState(t, secIdx, timeFactor, vecStateY + offset, adRes + offset, adY + offset, numSensAdDirs, errorTol);
	}

	// Phase 2: Calculate residual with current state

	// Evaluate residual for right hand side without time derivatives \dot{y} and store it in vecStateYdot
	// Also evaluate the Jacobian at the new position
	residualWithJacobian(active(t), secIdx, active(timeFactor), vecStateY, nullptr, vecStateYdot, adRes, adY, numSensAdDirs);

	// Note that we have omitted negating the residual as required. We will fix that later.

	// Phase 3: Solve for time derivatives
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];
		m->consistentInitialTimeDerivative(t, timeFactor, vecStateYdot + offset);
	}

	BENCH_STOP(_timerConsistentInit);
}

void ModelSystem::consistentIntialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active* const adRes, active* const adY)
{
	BENCH_START(_timerConsistentInit);

	// Compute parameter sensitivities and update the Jacobian
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];
		m->residualSensFwdWithJacobian(t, secIdx, timeFactor, vecStateY + offset, vecStateYdot + offset, adRes + offset, adY + offset, vecSensY.size());
	}

	// Connect units
	residualConnectUnitOps<double, active, active>(secIdx, vecStateY, vecStateYdot, adRes);

	// TODO: Adjust indexing / offset of vectors
	std::vector<double*> vecSensYlocal(vecSensY.size(), nullptr);
	std::vector<double*> vecSensYdotLocal(vecSensYdot.size(), nullptr);
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];

		// Use correct offset in sensitivity state vectors
		for (unsigned int j = 0; j < vecSensY.size(); ++j)
		{
			vecSensYlocal[j] = vecSensY[j] + offset;
			vecSensYdotLocal[j] = vecSensYdot[j] + offset;
		}

		m->consistentIntialSensitivity(t, secIdx, timeFactor, vecStateY + offset, vecStateYdot + offset, vecSensYlocal, vecSensYdotLocal, adRes + offset);
	}

	BENCH_STOP(_timerConsistentInit);
}

void ModelSystem::leanConsistentInitialConditions(double t, unsigned int secIdx, double timeFactor, double* const vecStateY, double* const vecStateYdot, 
	active* const adRes, active* const adY, unsigned int numSensAdDirs, double errorTol)
{
	BENCH_START(_timerConsistentInit);

	// Phase 1: Solve algebraic equations and update state
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];
		m->leanConsistentInitialState(t, secIdx, timeFactor, vecStateY + offset, adRes + offset, adY + offset, numSensAdDirs, errorTol);
	}

	// Phase 2: Calculate residual with current state

	std::vector<double> tempRes(numDofs(), 0.0);

	// Evaluate residual for right hand side without time derivatives \dot{y} and store it in tempRes
	// Also evaluate the Jacobian at the new position
	residualWithJacobian(active(t), secIdx, active(timeFactor), vecStateY, nullptr, tempRes.data(), adRes, adY, numSensAdDirs);

	// Note that we have omitted negating the residual as required. We will fix that later.

	// Phase 3: Solve for time derivatives
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];
		m->leanConsistentInitialTimeDerivative(t, timeFactor, vecStateYdot + offset, tempRes.data() + offset);
	}

	BENCH_STOP(_timerConsistentInit);
}

void ModelSystem::leanConsistentIntialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active* const adRes, active* const adY)
{
	BENCH_START(_timerConsistentInit);

	// Compute parameter sensitivities and update the Jacobian
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];
		m->residualSensFwdWithJacobian(t, secIdx, timeFactor, vecStateY + offset, vecStateYdot + offset, adRes + offset, adY + offset, vecSensY.size());
	}

	// Connect units
	residualConnectUnitOps<double, active, active>(secIdx, vecStateY, vecStateYdot, adRes);

	// TODO: Adjust indexing / offset of vectors
	std::vector<double*> vecSensYlocal(vecSensY.size(), nullptr);
	std::vector<double*> vecSensYdotLocal(vecSensYdot.size(), nullptr);
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];

		// Use correct offset in sensitivity state vectors
		for (unsigned int j = 0; j < vecSensY.size(); ++j)
		{
			vecSensYlocal[j] = vecSensY[j] + offset;
			vecSensYdotLocal[j] = vecSensYdot[j] + offset;
		}

		m->leanConsistentIntialSensitivity(t, secIdx, timeFactor, vecStateY + offset, vecStateYdot + offset, vecSensYlocal, vecSensYdotLocal, adRes + offset);
	}

	BENCH_STOP(_timerConsistentInit);
}

int ModelSystem::linearSolve(double t, double timeFactor, double alpha, double outerTol, double* const rhs, double const* const weight,
	double const* const y, double const* const yDot, double const* const res)
{
	BENCH_START(_timerLinearSolve);

	int result = 0;

	// TODO: Actually solve a system
	// TODO: Adjust indexing / offset of vectors
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];
		const int intermediateRes = m->linearSolve(t, timeFactor, alpha, outerTol, rhs + offset, weight + offset, y + offset, yDot + offset, res + offset);

		// If result is already -1 (non-recoverable error), then we stick to it
		// If result is ok or recoverable and intermediate result is recoverable, then we take intermediate result
		if ((result >= 0) && (intermediateRes > 0))
		{
			result = intermediateRes;
		}
	}

	BENCH_STOP(_timerLinearSolve);
	return result;
}

void ModelSystem::setSectionTimes(double const* secTimes, bool const* secContinuity, unsigned int nSections)
{
	for (IUnitOperation* m : _models)
		m->setSectionTimes(secTimes, secContinuity, nSections);	

	for (IExternalFunction* extFun : _extFunctions)
		extFun->setSectionTimes(secTimes, secContinuity, nSections);
}

/**
 * @brief Calculates error tolerances for additional coupling DOFs
 * @details ModelSystem uses additional DOFs to decouple a system of unit operations for parallelization.
 *          These additional DOFs don't get an error tolerance from the user because he shouldn't be
 *          aware of those (implementation detail). This function is responsible for calculating error
 *          tolerances for these additional coupling DOFs.
 * 
 * @param [in] errorTol Pointer to array of error tolerances for system without coupling DOFs
 * @param [in] errorTolLength Length of @p errorTol array
 * @return Vector with error tolerances for additional coupling DOFs
 */
std::vector<double> ModelSystem::calculateErrorTolsForAdditionalDofs(double const* errorTol, unsigned int errorTolLength)
{
	// Return empty vector since we don't have coupling DOFs, yet
	return std::vector<double>(0, 0.0);
}

void ModelSystem::expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut)
{
	// TODO: Adjust indexing / offset of vectors
	// TODO: Handle connection of unit operations
	for (IUnitOperation* m : _models)
	{
		m->expandErrorTol(errorSpec, errorSpecSize, expandOut);
	}
}

}  // namespace model

}  // namespace cadet
