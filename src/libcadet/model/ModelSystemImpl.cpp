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

#include "model/ModelSystemImpl.hpp"
#include "cadet/ParameterProvider.hpp"
#include "cadet/Exceptions.hpp"
#include "cadet/ExternalFunction.hpp"

#include "ConfigurationHelper.hpp"
#include "linalg/SparseMatrix.hpp"
#include "linalg/Norms.hpp"
#include "SensParamUtil.hpp"
#include "AdUtils.hpp"
#include "SimulationTypes.hpp"

#include <algorithm>
#include <string>
#include <iomanip>
#include <sstream>
#include <functional>
#include <iterator>
#include <limits>

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include "ParallelSupport.hpp"
#ifdef CADET_PARALLELIZE
	#include <tbb/tbb.h>
#endif

namespace
{
	/**
	 * @brief Find array index of unit operation with given id
	 * @details Unit operation index does not need to match index of unit operation in an array.
	 * @param [in] models List of models
	 * @param [in] unitOpIdx Unit operation index to look for
	 * @return Array index of the unit operation identified by unit operation index or invalid array index if the unit operation has not been found
	 */
	inline unsigned int indexOfUnitOp(const std::vector<cadet::IUnitOperation*>& models, unsigned int unitOpIdx)
	{
		for (unsigned int i = 0; i < models.size(); ++i)
		{
			if (models[i]->unitOperationId() == unitOpIdx)
				return i;
		}
		return models.size();
	}

	/**
	 * @brief Returns whether a given unit operation is a terminal node in the network
	 * @param [in] conn Connection list
	 * @param [in] size Number of rows in the connection list
	 * @param [in] unitOpIdx Index of the unit operation to check
	 * @return @c true if the unit operation is a terminal node in the network, @c false otherwise
	 */
	inline bool isTerminal(int const* const conn, unsigned int size, int unitOpIdx)
	{
		for (unsigned int i = 0; i < size; ++i)
		{
			if (conn[6*i] == unitOpIdx)
				return false;
		}
		return true;
	}

	/**
	 * @brief Selects either double or active SparseMatrix based on template argument
	 * @details Helper function that returns either @p a or @p b depending on the template argument.
	 * @param [in] a SparseMatrix of double elements
	 * @param [in] b SparseMatrix of active elements
	 * @tparam selector_t One of @c double or @c active
	 * @return Either @p a or @p b depending on the template argument
	 */
	template <class selector_t>
	const cadet::linalg::SparseMatrix<selector_t>& select(const cadet::linalg::SparseMatrix<double>& a, const cadet::linalg::SparseMatrix<cadet::active>& b)
	{
		cadet_assert(false);
	}

	template <>
	const cadet::linalg::SparseMatrix<double>& select<double>(const cadet::linalg::SparseMatrix<double>& a, const cadet::linalg::SparseMatrix<cadet::active>& b)
	{
		return a;
	}

	template <>
	const cadet::linalg::SparseMatrix<cadet::active>& select<cadet::active>(const cadet::linalg::SparseMatrix<double>& a, const cadet::linalg::SparseMatrix<cadet::active>& b)
	{
		return b;
	}

	/**
	 * @brief Computes a total return code from a list of separate return codes
	 * @details A negative return code indicates a non-recoverable error. Positive
	 *          values indicate recoverable errors and a value of @c 0 indicates
	 *          no error.
	 * @param [in] err List of error codes to be fused into one
	 * @return Total error code summarizing all codes in the list
	 */
	inline int totalErrorIndicatorFromLocal(const std::vector<int>& err)
	{
		int totalError = 0;
		for (unsigned int i = 0; i < err.size(); ++i)
		{
			// Negative values are non-recoverable errors
			if (err[i] < 0)
				return err[i];

			// 0 = okay, positive values = recoverable error
			totalError = std::max(totalError, err[i]);
		}
		return totalError;
	}

	/**
	 * @brief Fuses two error codes into one
	 * @details A negative return code indicates a non-recoverable error. Positive
	 *          values indicate recoverable errors and a value of @c 0 indicates
	 *          no error.
	 * @param [in] curCode Current error code
	 * @param [in] nextCode Next error code that is fused into the current one
	 * @return Fused error code summarizing both inputs
	 */
	inline int updateErrorIndicator(int curCode, int nextCode)
	{
		if ((curCode < 0) || (nextCode < 0))
			return std::min(curCode, nextCode);
		return std::max(curCode, nextCode);
	}

	struct FullTag {};
	struct LeanTag {};

	template <class tag_t>
	struct ConsistentInit {};

	template <>
	struct ConsistentInit<FullTag>
	{
		static inline void state(cadet::IUnitOperation* model, const cadet::SimulationTime& simTime, double* const vecStateY, const cadet::AdJacobianParams& adJac, double errorTol, cadet::util::ThreadLocalStorage& threadLocalMem)
		{
			model->consistentInitialState(simTime, vecStateY, adJac, errorTol, threadLocalMem);
		}

		static inline void timeDerivative(cadet::IUnitOperation* model, const cadet::SimulationTime& simTime, double const* vecStateY, double* const vecStateYdot, double* const res, cadet::util::ThreadLocalStorage& threadLocalMem)
		{
			model->consistentInitialTimeDerivative(simTime, vecStateY, vecStateYdot, threadLocalMem);
		}

		static inline int residualWithJacobian(cadet::model::ModelSystem& ms, const cadet::SimulationTime& simTime, const cadet::ConstSimulationState& simState, double* const res, double* const temp,
			const cadet::AdJacobianParams& adJac)
		{
			return ms.residualWithJacobian(simTime, simState, res, adJac);
		}

		static inline void parameterSensitivity(cadet::IUnitOperation* model, const cadet::SimulationTime& simTime, const cadet::ConstSimulationState& simState,
			std::vector<double*>& vecSensYlocal, std::vector<double*>& vecSensYdotLocal, cadet::active const* const adRes, cadet::util::ThreadLocalStorage& threadLocalMem)
		{
			model->consistentInitialSensitivity(simTime, simState, vecSensYlocal, vecSensYdotLocal, adRes, threadLocalMem);
		}
	};

	template <>
	struct ConsistentInit<LeanTag>
	{
		static inline void state(cadet::IUnitOperation* model, const cadet::SimulationTime& simTime, double* const vecStateY, const cadet::AdJacobianParams& adJac, double errorTol, cadet::util::ThreadLocalStorage& threadLocalMem)
		{
			model->leanConsistentInitialState(simTime, vecStateY, adJac, errorTol, threadLocalMem);
		}

		static inline void timeDerivative(cadet::IUnitOperation* model, const cadet::SimulationTime& simTime, double const* vecStateY, double* const vecStateYdot, double* const res, cadet::util::ThreadLocalStorage& threadLocalMem)
		{
			model->leanConsistentInitialTimeDerivative(simTime.t, vecStateY, vecStateYdot, res, threadLocalMem);
		}

		static inline int residualWithJacobian(cadet::model::ModelSystem& ms, const cadet::SimulationTime& simTime, const cadet::ConstSimulationState& simState, double* const res, double* const temp,
			const cadet::AdJacobianParams& adJac)
		{
			return ms.residualWithJacobian(simTime, simState, temp, adJac);
		}

		static inline void parameterSensitivity(cadet::IUnitOperation* model, const cadet::SimulationTime& simTime, const cadet::ConstSimulationState& simState,
			std::vector<double*>& vecSensYlocal, std::vector<double*>& vecSensYdotLocal, cadet::active const* const adRes, cadet::util::ThreadLocalStorage& threadLocalMem)
		{
			model->leanConsistentInitialSensitivity(simTime, simState, vecSensYlocal, vecSensYdotLocal, adRes, threadLocalMem);
		}
	};


	template <bool evalJacobian>
	struct ResidualSensCaller {};

	template <>
	struct ResidualSensCaller<true>
	{
		static inline int call(cadet::IUnitOperation* model, const cadet::SimulationTime& simTime, 
			const cadet::ConstSimulationState& simState, const cadet::AdJacobianParams& adJac, cadet::util::ThreadLocalStorage& threadLocalMem)
		{
			return model->residualSensFwdWithJacobian(simTime, simState, adJac, threadLocalMem);
		}
	};

	template <>
	struct ResidualSensCaller<false>
	{
		static inline int call(cadet::IUnitOperation* model, const cadet::SimulationTime& simTime, 
			const cadet::ConstSimulationState& simState, const cadet::AdJacobianParams& adJac, cadet::util::ThreadLocalStorage& threadLocalMem)
		{
			return model->residualSensFwdAdOnly(simTime, simState, adJac.adRes, threadLocalMem);
		}
	};

	template <typename state_t>
	inline state_t applyOffset(const state_t& state, unsigned int offset)
	{
		return state_t{
			state.vecStateY + offset,
			(state.vecStateYdot) ? (state.vecStateYdot + offset) : nullptr
		};
	}

	template <>
	inline cadet::AdJacobianParams applyOffset(const cadet::AdJacobianParams& adJac, unsigned int offset)
	{
		return cadet::AdJacobianParams{
			(adJac.adRes) ? (adJac.adRes + offset) : nullptr,
			(adJac.adY) ? (adJac.adY + offset) : nullptr,
			adJac.adDirOffset
		};
	}
}

namespace cadet
{

namespace model
{

ModelSystem::ModelSystem() : _jacNF(nullptr), _jacFN(nullptr), _jacActiveFN(nullptr), _curSwitchIndex(0), _tempState(nullptr), _initState(0, 0.0), _initStateDot(0, 0.0)
{
}

ModelSystem::~ModelSystem() CADET_NOEXCEPT
{
	for (IUnitOperation* model : _models)
		delete model;

	for (IExternalFunction* extFun : _extFunctions)
		delete extFun;

	delete[] _tempState;
	delete[] _jacNF;
	delete[] _jacFN;
	delete[] _jacActiveFN;
}

void ModelSystem::addModel(IModel* unitOp)
{
	// Check for unique unit operation id
	if (indexOfUnitOp(_models, unitOp->unitOperationId()) < _models.size())
		throw InvalidParameterException("Cannot add model because of already existing unit operation id " + std::to_string(unitOp->unitOperationId()));

	IUnitOperation* const uo = static_cast<IUnitOperation*>(unitOp);
	_models.push_back(uo);

	if (uo->hasInlet() && uo->hasOutlet())
		_inOutModels.push_back(_models.size() - 1);

	// Propagate external functions to submodel
	uo->setExternalFunctions(_extFunctions.data(), _extFunctions.size());
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

IUnitOperation* ModelSystem::getUnitOperationModel(unsigned int unitOpIdx)
{
	for (IUnitOperation* m : _models)
	{
		if (m->unitOperationId() == unitOpIdx)
			return m;
	}
	return nullptr;
}

IUnitOperation const* ModelSystem::getUnitOperationModel(unsigned int unitOpIdx) const
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

unsigned int ModelSystem::totalNumInletPorts() const CADET_NOEXCEPT
{
	unsigned int nPorts = 0;
	for (IUnitOperation const* m : _models)
		nPorts += m->numInletPorts();

	return nPorts;
}

unsigned int ModelSystem::totalNumOutletPorts() const CADET_NOEXCEPT
{
	unsigned int nPorts = 0;
	for (IUnitOperation const* m : _models)
		nPorts += m->numOutletPorts();

	return nPorts;
}

unsigned int ModelSystem::maxUnitInletPorts() const CADET_NOEXCEPT
{
	unsigned int nPorts = 0;
	for (IUnitOperation const* m : _models)
		nPorts = std::max(nPorts, m->numInletPorts());

	return nPorts;
}

unsigned int ModelSystem::maxUnitOutletPorts() const CADET_NOEXCEPT
{
	unsigned int nPorts = 0;
	for (IUnitOperation const* m : _models)
		nPorts = std::max(nPorts, m->numOutletPorts());

	return nPorts;
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
	return _dofOffset.back() + numCouplingDOF();
}

unsigned int ModelSystem::numPureDofs() const CADET_NOEXCEPT
{
	unsigned int dofs = 0;
	for (IUnitOperation* m : _models)
		dofs += m->numPureDofs();
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

/**
* @brief Create data structures to keep track of entries and locations in the state vector
* @details Three data structures are created. One keeps track of the offsets to each unit operation.
*          The second keeps track of the size of each unit operation. The third is mapping from
*          unit operation and component index to unique inlet DOF index.
*/
void ModelSystem::rebuildInternalDataStructures()
{
	// Calculate array with DOF offsets
	_dofOffset.clear();
	_dofs.clear();

	// The additional entry holds the offset for the superstructure
	_dofOffset.reserve(_models.size()+1);
	_dofs.reserve(_models.size() + 1);

	// Process DOF from models
	unsigned int totalDof = 0;
	for (IUnitOperation const* m : _models)
	{
		_dofOffset.push_back(totalDof);
		totalDof += m->numDofs();
		_dofs.push_back(m->numDofs());
	}

	// Process DOF from superstructure
	_dofOffset.push_back(totalDof);

	/*
		A mapping is required that turns a local model, port index, and component index into the location of the inlet DOF in
		the global state vector. Some unit operations do not have inlet DOFs (e.g., inlet unit operation). Hence,
		a map is constructed which converts local indices into global ones.
	*/

	// Build a mapping (unitOpIdx, PortIdx, CompIdx) -> local coupling DOF index
	_couplingIdxMap.clear();
	unsigned int counter = 0;
	for (unsigned int i = 0; i < numModels(); ++i)
	{
		IUnitOperation const* const model = _models[i];
		
		//Only unit operations with an inlet have dedicated inlet DOFs
		if (model->hasInlet())
		{
			for (unsigned int port = 0; port < model->numInletPorts(); ++port)
			{
				for (unsigned int comp = 0; comp < model->numComponents(); ++comp)
				{
					_couplingIdxMap.insert({ std::make_tuple(i, port, comp), counter });
					++counter;
				}
			}
		}
	}

	_dofs.push_back(numCouplingDOF());

	// Allocate error indicator vector
	_errorIndicator.resize(_models.size(), 0);

	LOG(Debug) << "DOF offsets: " << _dofOffset;
}

/**
* @brief Allocates memory for the superstructure matrices
* @details How many connections each unit has determines how much memory has to be allocated for the coupling matrices.
*          This function walks the connections over the entire simulation in order to determine the maximum number of 
*          connections over the whole simulation which governs the number of entries in the sparse superstructure coupling
*          matrices. Finally, the required memory is allocated in the matrices.
*/
void ModelSystem::allocateSuperStructMatrices()
{
	// Step 1: Calculate number of connections per unit operation per valve switch
	// We record the number of outgoing connections (i.e., components) which act as 
	// sources in the bottom macro-row of the superstructure
	const unsigned int nModels = numModels();
	std::vector<unsigned int> sourcesPerUnitOpPerConfig(nModels * _switchSectionIndex.size(), 0u);

	for (unsigned int idx = 0; idx < _switchSectionIndex.size(); ++idx)
	{
		int const* ptrConn = _connections[idx];

		for (unsigned int i = 0; i < _connections.sliceSize(idx) / 6; ++i, ptrConn += 6)
		{
			// Extract current connection
			const int uoSource = ptrConn[0];
			const int portSource = ptrConn[2];
			const int compSource = ptrConn[4];

			IUnitOperation const* const model = _models[uoSource];
			if (portSource == -1)
			{
				if (compSource == -1)
					sourcesPerUnitOpPerConfig[nModels * idx + uoSource] += model->numOutletPorts() * model->numComponents();
				else
					sourcesPerUnitOpPerConfig[nModels * idx + uoSource] += model->numOutletPorts();
			}
			else
			{
				if (compSource == -1)
					sourcesPerUnitOpPerConfig[nModels * idx + uoSource] += model->numComponents();
				else
					sourcesPerUnitOpPerConfig[nModels * idx + uoSource] += 1;
			}
		}
	}

	// Step 2: Take maximum over valve switches to obtain maximum number of matrix entries per unit operation
	std::vector<unsigned int> numOutgoing(nModels, 0u);

	// Assign the first row of sourcesPerUnitOpPerConfig to the maximum number of outputs
	std::copy(sourcesPerUnitOpPerConfig.begin(), sourcesPerUnitOpPerConfig.begin() + nModels, numOutgoing.begin());

	// Loop over remaining lines and take maximum per element
	for (unsigned int sectionIdx = 1; sectionIdx < _switchSectionIndex.size(); ++sectionIdx)
	{
		for (unsigned int i = 0; i < nModels; ++i)
		{
			numOutgoing[i] = std::max(numOutgoing[i], sourcesPerUnitOpPerConfig[nModels * sectionIdx + i]);
		}
	}

	// Step 3: Allocate memory based on maximum number of connections for each unit operation
	for (unsigned int i = 0; i < numModels(); ++i)
	{
		// Bottom macro-row
		_jacActiveFN[i].resize(numOutgoing[i]);

		// Right macro-column
		// Each unit operation has inlets equal to numComponents * numInletPorts so long as the unit operation has an inlet
		IUnitOperation const* const model = _models[i];
		if (model->hasInlet())
			_jacNF[i].resize(model->numComponents() * model->numInletPorts());
	}
}

bool ModelSystem::configureModelDiscretization(IParameterProvider& paramProvider, IConfigHelper& helper)
{
	// Discretizations of unit operation models are already configured
	rebuildInternalDataStructures();

	_parameters.clear();
	configureSwitches(paramProvider);
	_curSwitchIndex = 0;

	// Allocate memory to coupling matrices
	_jacActiveFN = new linalg::SparseMatrix<active>[numModels()];
	_jacNF = new linalg::SparseMatrix<double>[numModels()];
	_jacFN = new linalg::SparseMatrix<double>[numModels()];

	// Calculate the sizes that need to be allocated for all the inlets and outlets
	allocateSuperStructMatrices();

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

	// Read solver settings
	paramProvider.pushScope("solver");
	const int maxKrylov = paramProvider.getInt("MAX_KRYLOV");
	const int gsType = paramProvider.getInt("GS_TYPE");
	const int maxRestarts = paramProvider.getInt("MAX_RESTARTS");
	_schurSafety = paramProvider.getDouble("SCHUR_SAFETY");
	paramProvider.popScope();

	// Initialize and configure GMRES for solving the Schur-complement
    _gmres.initialize(numCouplingDOF(), maxKrylov, linalg::toOrthogonalization(gsType), maxRestarts);

	// Allocate tempState vector
	delete[] _tempState;
	_tempState = new double[numDofs()];

//	_tempSchur = new double[*std::max_element(_dofs.begin(), _dofs.end())];
	_flowRateIn.resize(maxUnitInletPorts(), 0.0);
	_flowRateOut.resize(maxUnitOutletPorts(), 0.0);

	_totalInletFlow.reserve(totalNumInletPorts(), _models.size());
	for (IUnitOperation const* m : _models)
		_totalInletFlow.pushBackSlice(m->numInletPorts());

#ifdef CADET_DEBUG

	LOG(Debug) << "Num units " << _models.size() << " max ID " << maxUnitOperationId() << " total DOF " << numDofs();
	LOG(Debug) << "Unit op state vector size: uoDOFs = " << log::VectorPtr<unsigned int>(_dofs.data(), _dofs.size());
	LOG(Debug) << "Unit op state vector offsets: uoOffset = " << log::VectorPtr<unsigned int>(_dofOffset.data(), _dofOffset.size());

#endif

	return success;
}

bool ModelSystem::configure(IParameterProvider& paramProvider)
{
	_parameters.clear();

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

	// Reconfigure solver settings
	paramProvider.pushScope("solver");
	const int gsType = paramProvider.getInt("GS_TYPE");
	const int maxRestarts = paramProvider.getInt("MAX_RESTARTS");
	_schurSafety = paramProvider.getDouble("SCHUR_SAFETY");
	paramProvider.popScope();

	_gmres.orthoMethod(linalg::toOrthogonalization(gsType));
	_gmres.maxRestarts(maxRestarts);

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
	
	// TODO Improve those very conservative size estimates
	_switchSectionIndex.clear();
	_switchSectionIndex.reserve(numSwitches);
	_connections.clear();
	_connections.reserve(numSwitches * 6 * _models.size() * _models.size(), numSwitches);
	_flowRates.clear();
	_flowRates.reserve(numSwitches * _models.size() * _models.size(), numSwitches);

#if CADET_COMPILER_CXX_CONSTEXPR
	constexpr StringHash flowHash = hashString("CONNECTION");
#else
	const StringHash flowHash = hashString("CONNECTION");
#endif

	std::ostringstream oss;
	for (unsigned int i = 0; i < numSwitches; ++i)
	{
		oss.str("");
		oss << "switch_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << i;

		paramProvider.pushScope(oss.str());

		_switchSectionIndex.push_back(paramProvider.getInt("SECTION"));

		if ((i > 0) && (_switchSectionIndex.back() <= _switchSectionIndex[_switchSectionIndex.size() - 2]))
			throw InvalidParameterException("SECTION index has to be monotonically increasing (switch " + std::to_string(i) + ")");

		std::vector<double> connFlow = paramProvider.getDoubleArray("CONNECTIONS");
		if ((connFlow.size() % 7) != 0)
			throw InvalidParameterException("CONNECTIONS matrix has to have 7 columns");

		std::vector<int> conn(connFlow.size() / 7 * 6, 0);
		std::vector<double> fr(connFlow.size() / 7, 0.0);
		
		checkConnectionList(connFlow, conn, fr, i);

		_connections.pushBackSlice(conn);

		// Convert double to active while pushing into the SlicedVector
		// also register parameter to enable sensitivities
		if (fr.size() > 0)
		{
			_flowRates.pushBack(fr[0]);
			_parameters[makeParamId(flowHash, UnitOpIndep, conn[2], conn[3], conn[0], conn[1], _switchSectionIndex.back())] = _flowRates.back();
			for (unsigned int j = 1; j < fr.size(); ++j)
			{
				_flowRates.pushBackInLastSlice(fr[j]);

				// Check if a previous identical connection (except for component indices) exists
				bool found = false;
				for (unsigned int k = 0; k < j; ++k)
				{
					if ((conn[6*k] == conn[6*j]) && (conn[6*k+1] == conn[6*j+1]) && (conn[6*k+2] == conn[6*j+2]) && (conn[6*k+3] == conn[6*j+3]))
					{
						found = true;
						break;
					}
				}

				// Only register the first occurrence of a flow parameter
				if (!found)
					_parameters[makeParamId(flowHash, UnitOpIndep, conn[6*j+2], conn[6*j+3], conn[6*j], conn[6*j+1], _switchSectionIndex.back())] = (_flowRates.back() + j);
			}
		}
		else
			// Add empty slice
			_flowRates.pushBackSlice(nullptr, 0);

		paramProvider.popScope();
	}

	paramProvider.popScope();

	if (_switchSectionIndex[0] != 0)
		throw InvalidParameterException("First element of SECTION in connections group has to be 0");
}

bool ModelSystem::configureModel(IParameterProvider& paramProvider, unsigned int unitOpIdx)
{
	IUnitOperation* const model = getUnitOperationModel(unitOpIdx);
	if (!model)
		return false;

	return model->configure(paramProvider);
}

/**
 * @brief Checks the given unit operation connection list and reformats it
 * @details Throws an exception if something is incorrect. Reformats the connection list by
 *          substituting unit operation IDs with local indices.
 * @param [in] conn Matrix with 7 columns holding all connections. The matrix is expected
 *             to be in row-major storage format.
 * @param [out] connOnly Matrix with 6 columns holding all connections. While @p conn contains
 *              connection indices and flow rate, this matrix only holds connection indices.
 *              It is assumed to be pre-allocated (same number of rows as @p conn). The unit
 *              operation IDs are substituted by the corresponding indices of the unit operations
 *              in the local _models vector.
 * @param [out] flowRates Vector with flow rates for each connection in the list. It is assumed
 *              to be pre-allocated (same number of rows as @p conn).
 * @param [in] idxSwitch Index of the valve switch that corresponds to this connection list
 */
void ModelSystem::checkConnectionList(const std::vector<double>& conn, std::vector<int>& connOnly, std::vector<double>& flowRates, unsigned int idxSwitch) const
{
	std::vector<double> totalInflow(_models.size(), 0.0);
	std::vector<double> totalOutflow(_models.size(), 0.0);
	for (unsigned int i = 0; i < conn.size() / 7; ++i)
	{
		// Extract current connection
		int uoSource = static_cast<int>(conn[7*i]);
		int uoDest = static_cast<int>(conn[7*i+1]);
		const int portSource = static_cast<int>(conn[7*i+2]);
		const int portDest = static_cast<int>(conn[7*i+3]);
		const int compSource = static_cast<int>(conn[7*i+4]);
		const int compDest = static_cast<int>(conn[7*i+5]);
		double fr = conn[7*i + 6];

		if (uoSource < 0)
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Source unit operation id has to be at least 0 in connection");
		if (uoDest < 0)
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Destination unit operation id has to be at least 0 in connection");

		// Convert to index
		uoSource = indexOfUnitOp(_models, uoSource);
		uoDest = indexOfUnitOp(_models, uoDest);

		if (static_cast<unsigned int>(uoSource) >= _models.size())
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Source unit operation id " + std::to_string(uoSource) + " not found in connection");
		if (static_cast<unsigned int>(uoDest) >= _models.size())
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Destination unit operation id " + std::to_string(uoDest) + " not found in connection");

		// Check if unit operations have inlets and outlets
		if (!_models[uoSource]->hasOutlet())
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Source unit operation " + std::to_string(_models[uoSource]->unitOperationId()) + " does not have an outlet");
		if (!_models[uoDest]->hasInlet())
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Source unit operation " + std::to_string(_models[uoDest]->unitOperationId()) + " does not have an inlet");

		// Check port indices
		if ((portSource >= 0) && (static_cast<unsigned int>(portSource) >= _models[uoSource]->numOutletPorts()))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Source port index (" + std::to_string(portSource) + ") exceeds number of outlet ports (" + std::to_string(_models[uoSource]->numOutletPorts()) + ")");
		if ((portDest >= 0) && (static_cast<unsigned int>(portDest) >= _models[uoDest]->numInletPorts()))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Destination port index (" + std::to_string(portDest) + ") exceeds number of inlet ports (" + std::to_string(_models[uoDest]->numInletPorts()) + ")");

		if (((portSource < 0) && (portDest >= 0)) || ((portSource >= 0) && (portDest < 0)))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Only source or destination (not both) are set to connect all ports in connection from unit " 
				+ std::to_string(_models[uoSource]->unitOperationId()) + " to " + std::to_string(_models[uoDest]->unitOperationId()));

		if ((portSource < 0) && (portDest < 0) && (_models[uoSource]->numOutletPorts() != _models[uoDest]->numInletPorts()))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Number of ports not equal when connecting all ports from unit " 
				+ std::to_string(_models[uoSource]->unitOperationId()) + " (" + std::to_string(_models[uoSource]->numOutletPorts()) + " ports) to " + std::to_string(_models[uoDest]->unitOperationId()) + " (" + std::to_string(_models[uoDest]->numInletPorts()) + " ports)");

		// Check component indices
		if ((compSource >= 0) && (static_cast<unsigned int>(compSource) >= _models[uoSource]->numComponents()))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Source component index (" + std::to_string(compSource) + ") exceeds number of components (" + std::to_string(_models[uoSource]->numComponents()) + ")");
		if ((compDest >= 0) && (static_cast<unsigned int>(compDest) >= _models[uoDest]->numComponents()))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Destination component index (" + std::to_string(compDest) + ") exceeds number of components (" + std::to_string(_models[uoDest]->numComponents()) + ")");

		if (((compSource < 0) && (compDest >= 0)) || ((compSource >= 0) && (compDest < 0)))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Only source or destination (not both) are set to connect all components in connection from unit " 
				+ std::to_string(_models[uoSource]->unitOperationId()) + " to " + std::to_string(_models[uoDest]->unitOperationId()));

		if ((compSource < 0) && (compDest < 0) && (_models[uoSource]->numComponents() != _models[uoDest]->numComponents()))
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + " row " + std::to_string(i) + "): Number of components not equal when connecting all components from unit " 
				+ std::to_string(_models[uoSource]->unitOperationId()) + " to " + std::to_string(_models[uoDest]->unitOperationId()));
	
		// Add connection to index matrix
		connOnly[6*i] = uoSource;
		connOnly[6*i+1] = uoDest;
		connOnly[6*i+2] = portSource;
		connOnly[6*i+3] = portDest;
		connOnly[6*i+4] = compSource;
		connOnly[6*i+5] = compDest;

		// Add flow rate of connection to balance

		// Check if such a connection has occurred before (for a different component)
		bool found = false;
		for (unsigned int j = 0; j < i; ++j)
		{
			if ((conn[7*j] == uoSource) && (uoDest == conn[7*j+1]) && (conn[7*j+2] == portSource) && (portDest == conn[7*j+3]))
			{
				// Take flow rate that appears first
				fr = conn[7*j+6];
				found = true;
				break;
			}
		}

		// Total flow rates: Only add flow rate once (not for each component)
		if (!found)
		{
			// Add flow rates to balance
			totalInflow[uoDest] += fr;
			totalOutflow[uoSource] += fr;
		}

		// Add flow rate to list
		flowRates[i] = fr;
	}

	// Check flow rate balance
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		// Unit operations with only one port (inlet or outlet) do not need to balance their flows
		if ((totalInflow[i] >= 0.0) && (totalOutflow[i] == 0.0) && _models[i]->hasInlet() && !_models[i]->hasOutlet())
			continue;
		if ((totalInflow[i] == 0.0) && (totalOutflow[i] >= 0.0) && !_models[i]->hasInlet() && _models[i]->hasOutlet())
			continue;

		// Terminal unit operations do not need to balance their flows
		if ((totalOutflow[i] >= 0.0) && isTerminal(connOnly.data(), connOnly.size() / 6, i))
			continue;

		// Check balance and account for whether accumulation is allowed
		const double diff = std::abs(totalInflow[i] - totalOutflow[i]);
		if (((diff >= 1e-15) || (diff >= 1e-15 * std::abs(totalOutflow[i]))) && !_models[i]->canAccumulate())
			throw InvalidParameterException("In CONNECTIONS matrix (switch " + std::to_string(idxSwitch) + "): Flow rate balance is not closed for unit operation " + std::to_string(i) + ", imbalanced by " + std::to_string(diff));
	}

	// TODO: Check for conflicting entries
	// TODO: Plausibility check of total connections
}

std::unordered_map<ParameterId, double> ModelSystem::getAllParameterValues() const
{
	std::unordered_map<ParameterId, double> data;
	std::transform(_parameters.begin(), _parameters.end(), std::inserter(data, data.end()),
	               [](const std::pair<const ParameterId, active*>& p) { return std::make_pair(p.first, static_cast<double>(*p.second)); });

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
	return _parameters.find(pId) != _parameters.end();
}

double ModelSystem::getParameterDouble(const ParameterId& pId) const
{
	for (IUnitOperation* m : _models)
	{
		if ((m->unitOperationId() == pId.unitOperation) || (pId.unitOperation == UnitOpIndep))
		{
			if (m->hasParameter(pId))
			{
				const double val = m->getParameterDouble(pId);
				if (!std::isnan(val))
					return val;
			}
		}
	}
	const std::unordered_map<ParameterId, active*>::const_iterator it = _parameters.find(pId);
	if (it == _parameters.end())
		return std::numeric_limits<double>::quiet_NaN();

	return static_cast<double>(*it->second);
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
	bool found = false;
	auto paramHandle = _parameters.find(pId);
	if (paramHandle != _parameters.end())
	{
		paramHandle->second->setValue(value);
		found = true;

		// Multiplex flow rate parameters
#if CADET_COMPILER_CXX_CONSTEXPR
		constexpr StringHash flowHash = hashString("CONNECTION");
#else
		const StringHash flowHash = hashString("CONNECTION");
#endif
		if (flowHash == pId.name)
		{
			// Find the index of the valve switch
			const auto it = std::find(_switchSectionIndex.begin(), _switchSectionIndex.end(), pId.section);
			if (it != _switchSectionIndex.end())
			{
				const unsigned int idxSwitch = std::distance(_switchSectionIndex.begin(), it);
				int const* ptrConn = _connections[idxSwitch];
				active* conRates = _flowRates[idxSwitch];

				// Set value for all flow rates of the same connection (except for components)
				for (unsigned int i = 0; i < _connections.sliceSize(idxSwitch) / 6; ++i, ptrConn += 6, ++conRates)
				{
					if ((ptrConn[2] != pId.component) || (ptrConn[3] != pId.particleType) || (ptrConn[0] != pId.boundState) || (ptrConn[1] != pId.reaction))
						continue;

					conRates->setValue(value);
				}
			}
		}
	}

	return setParameterImpl(pId, value) || found;
}

bool ModelSystem::setParameter(const ParameterId& pId, bool value)
{
	return setParameterImpl(pId, value);
}

void ModelSystem::setSensitiveParameterValue(const ParameterId& pId, double value)
{
	if (pId.unitOperation == UnitOpIndep)
	{
		// Handle flow rates
		auto paramHandle = _parameters.find(pId);
		if ((paramHandle != _parameters.end()) && contains(_sensParams, paramHandle->second))
		{
			paramHandle->second->setValue(value);

			// Multiplex flow rate parameters
#if CADET_COMPILER_CXX_CONSTEXPR
			constexpr StringHash flowHash = hashString("CONNECTION");
#else
			const StringHash flowHash = hashString("CONNECTION");
#endif
			if (flowHash == pId.name)
			{
				// Find the index of the valve switch
				const auto it = std::find(_switchSectionIndex.begin(), _switchSectionIndex.end(), pId.section);
				if (it != _switchSectionIndex.end())
				{
					const unsigned int idxSwitch = std::distance(_switchSectionIndex.begin(), it);
					int const* ptrConn = _connections[idxSwitch];
					active* conRates = _flowRates[idxSwitch];

					// Set value for all flow rates of the same connection (except for components)
					for (unsigned int i = 0; i < _connections.sliceSize(idxSwitch) / 6; ++i, ptrConn += 6, ++conRates)
					{
						if ((ptrConn[2] != pId.component) || (ptrConn[3] != pId.particleType) || (ptrConn[0] != pId.boundState) || (ptrConn[1] != pId.reaction))
							continue;

						conRates->setValue(value);
					}
				}
			}
		}
	}

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

	// Check own parameters
	if (pId.unitOperation == UnitOpIndep)
	{
		auto paramHandle = _parameters.find(pId);
		if (paramHandle != _parameters.end())
		{
			LOG(Debug) << "Found parameter " << pId << " in ModelSystem: Dir " << adDirection << " is set to " << adValue;

			// Register parameter and set AD seed / direction
			_sensParams.insert(paramHandle->second);
			paramHandle->second->setADValue(adDirection, adValue);

			// Multiplex flow rate parameters
#if CADET_COMPILER_CXX_CONSTEXPR
			constexpr StringHash flowHash = hashString("CONNECTION");
#else
			const StringHash flowHash = hashString("CONNECTION");
#endif
			if (flowHash == pId.name)
			{
				// Find the index of the valve switch
				const auto it = std::find(_switchSectionIndex.begin(), _switchSectionIndex.end(), pId.section);
				if (it != _switchSectionIndex.end())
				{
					const unsigned int idxSwitch = std::distance(_switchSectionIndex.begin(), it);
					int const* ptrConn = _connections[idxSwitch];
					active* conRates = _flowRates[idxSwitch];

					// Set AD direction for all flow rates of the same connection (except for components)
					for (unsigned int i = 0; i < _connections.sliceSize(idxSwitch) / 6; ++i, ptrConn += 6, ++conRates)
					{
						if ((ptrConn[2] != pId.component) || (ptrConn[3] != pId.particleType) || (ptrConn[0] != pId.boundState) || (ptrConn[1] != pId.reaction))
							continue;

						conRates->setADValue(adDirection, adValue);
					}
				}
			}

			found = true;
		}
	}

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
	// Remove AD directions from parameters
	for (auto sp : _sensParams)
		sp->setADValue(0.0);

	_sensParams.clear();

	// Propagate call to models
	for (IUnitOperation* m : _models)
		m->clearSensParams();
}

void ModelSystem::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, const AdJacobianParams& adJac)
{
	// Check if simulation is (re-)starting from the very beginning
	if (secIdx == 0)
		_curSwitchIndex = 0;

	const unsigned int wrapSec = secIdx % _switchSectionIndex.size();
	const unsigned int prevSwitch = _curSwitchIndex;

	// If there are still some switches left and the next switch occurs in this section, advance index
	if ((_curSwitchIndex < _switchSectionIndex.size() - 1) && (_switchSectionIndex[_curSwitchIndex + 1] <= wrapSec))
	{
		++_curSwitchIndex;
	}
	else if (_curSwitchIndex == _switchSectionIndex.size() - 1)
	{
		// We're in the last valve configuration, let's check if we should cycle back to the first one
		if (_switchSectionIndex[0] == wrapSec)
			_curSwitchIndex = 0;
	}

	// Notify models that a discontinuous section transition has happened
	int const* ptrConn = _connections[_curSwitchIndex];
	active const* const conRates = _flowRates[_curSwitchIndex];
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		const unsigned int offset = _dofOffset[i];
		std::fill(_flowRateIn.begin(), _flowRateIn.end(), 0.0);
		std::fill(_flowRateOut.begin(), _flowRateOut.end(), 0.0);

		// Compute total inlet and outlet flow rate for this unit operation by traversing connection list
		for (unsigned int j = 0; j < _connections.sliceSize(_curSwitchIndex) / 6; ++j)
		{
			const int uoSource = ptrConn[6*j];
			const int uoDest = ptrConn[6*j+1];
			const int portSource = ptrConn[6*j+2];
			const int portDest = ptrConn[6*j+3];

			// Make sure this is the first connection (there may be several with different components)
			bool skip = false;
			for (unsigned int k = 0; k < j; ++k)
			{
				if ((ptrConn[6*k] == uoSource) && (ptrConn[6*k+1] == uoDest) && (ptrConn[6*k+2] == portSource) && (ptrConn[6*k+3] == portDest))
				{
					skip = true;
					break;
				}
			}

			// Skip this row in connection list if there was a previous connection
			if (skip)
				continue;

			if (uoSource == i)
			{
				if (portSource >= 0)
					_flowRateOut[portSource] += conRates[j];
				else
				{
					for (unsigned int k = 0; k < _models[i]->numOutletPorts(); ++k)
						_flowRateOut[k] += conRates[j];
				}
			}
			if (uoDest == i)
			{
				if (portDest >= 0)
					_flowRateIn[portDest] += conRates[j];
				else
				{
					for (unsigned int k = 0; k < _models[i]->numInletPorts(); ++k)
						_flowRateIn[k] += conRates[j];
				}
			}
		}

		_models[i]->setFlowRates(_flowRateIn.data(), _flowRateOut.data());
		_models[i]->notifyDiscontinuousSectionTransition(t, secIdx, applyOffset(adJac, offset));
	}

#ifdef CADET_DEBUG
	LOG(Debug) << "Switching from valve configuration " << prevSwitch << " to " << _curSwitchIndex << " (sec = " << secIdx << " wrapSec = " << wrapSec << ")";
	for (unsigned int i = 0; i < _connections.sliceSize(_curSwitchIndex) / 6; ++i, ptrConn += 6)
	{
		// Extract current connection
		const int uoSource = ptrConn[0];
		const int uoDest = ptrConn[1];
		const int portSource = ptrConn[2];
		const int portDest = ptrConn[3];
		const int compSource = ptrConn[4];
		const int compDest = ptrConn[5];

		//Number of components was already verified so assume they are all correct

		LOG(Debug) << "Unit op " << uoSource << " (" << _models[uoSource]->unitOperationName() << ") port " << portSource << " comp " << compSource << " => "
		           << uoDest << " (" << _models[uoDest]->unitOperationName() << ") port " << portDest << " comp " << compDest;
	}
#endif

	if ((0 == secIdx) || (prevSwitch != _curSwitchIndex))
		assembleSuperStructMatrices(secIdx);		
}

/**
* @brief Rebuild the outer network connection matrices in the super structure
* @details Rebuild NF and FN matrices. This should only be called if the connections have changed. 
* @param [in] secIdx Section index
*/
void ModelSystem::assembleSuperStructMatrices(unsigned int secIdx)
{
	// Clear the matrices before we set new entries
	for (unsigned int i = 0; i < numModels(); ++i)
	{
		_jacNF[i].clear();
		_jacActiveFN[i].clear();
	}

	// Assemble Jacobian submatrices

	// Right macro-column
	// NF
	unsigned int couplingIdx = 0;
	for (unsigned int i = 0; i < numModels(); ++i)
	{
		IUnitOperation const* const model = _models[i];
		
		// Only items with an inlet have non-zero entries in the NF matrices
		if (model->hasInlet())
		{
			for (unsigned int port = 0; port < model->numInletPorts(); ++port)
			{
				// Each component generates a -1 for its inlet in the NF[i] matrix and increases the couplingIdx by 1
				const unsigned int localInletComponentIndex = model->localInletComponentIndex(port);
				const unsigned int localInletComponentStride = model->localInletComponentStride(port);
				for (unsigned int comp = 0; comp < model->numComponents(); ++comp)
				{
					_jacNF[i].addElement(localInletComponentIndex + comp * localInletComponentStride, couplingIdx, -1.0);
					++couplingIdx;
				}
			}
		}
	}

	// Calculate total flow rate for each inlet
	int const* const ptrConn = _connections[_curSwitchIndex];
	active const* const ptrRate = _flowRates[_curSwitchIndex];

	// Reset _totalInletFlow back to zero
	_totalInletFlow.fill(0.0);

	// Compute total volumetric inflow for each unit operation port
	for (unsigned int i = 0; i < _connections.sliceSize(_curSwitchIndex) / 6; ++i)
	{
		// Extract current connection
		const int uoSource = ptrConn[6*i];
		const int uoDest = ptrConn[6*i + 1];
		const int portSource = ptrConn[6*i + 2];
		const int portDest = ptrConn[6*i + 3];

		// Check if the same connection has appeared before (with different components)
		bool skip = false;
		for (unsigned int j = 0; j < i; ++j)
		{
			if ((ptrConn[6*j] == uoSource) && (ptrConn[6*j + 1] == uoDest) && (ptrConn[6*j + 2] == portSource) && (ptrConn[6*j + 3] == portDest))
			{
				skip = true;
				break;
			}
		}

		// Skip this row in connection list if there was an identical previous connection (except for component indices)
		if (skip)
			continue;

		// Use the first flow rate from uoSource to uoDest
		if (portDest < 0)
		{
			for (unsigned int j = 0; j < _models[uoDest]->numInletPorts(); ++j)
				_totalInletFlow(uoDest, j) += ptrRate[i];
		}
		else
			_totalInletFlow(uoDest, portDest) += ptrRate[i];
	}

	// Bottom macro-row
	// FN
	for (unsigned int i = 0; i < _connections.sliceSize(_curSwitchIndex) / 6; ++i)
	{
		// Extract current connection
		const int uoSource = ptrConn[6*i];
		const int uoDest = ptrConn[6*i + 1];
		const int portSource = ptrConn[6*i + 2];
		const int portDest = ptrConn[6*i + 3];
		const int compSource = ptrConn[6*i + 4];
		const int compDest = ptrConn[6*i + 5];

		// Obtain index of first connection from uoSource to uoDest
		unsigned int idx = i;
		for (unsigned int j = 0; j < i; ++j)
		{
			if ((ptrConn[6*j] == uoSource) && (ptrConn[6*j + 1] == uoDest) && (ptrConn[6*j + 2] == portSource) && (ptrConn[6*j + 3] == portDest))
			{
				idx = j;
				break;
			}
		}

		// idx contains the index of the first connection from uoSource to uoDest
		// Hence, ptrRate[idx] is the flow rate to use for this connection

		IUnitOperation const* const modelSource = _models[uoSource];

		// The outlet column is the outlet index + component number * outlet stride

		if (portSource == -1)
		{
			for (unsigned int j = 0; j < modelSource->numOutletPorts(); ++j)
			{
				const unsigned int outletIndex = modelSource->localOutletComponentIndex(j);
				const unsigned int outletStride = modelSource->localOutletComponentStride(j);

				if (compSource == -1)
				{
					// Connect all components with the same flow rate
					for (unsigned int comp = 0; comp < modelSource->numComponents(); ++comp)
					{
						const unsigned int row = _couplingIdxMap[std::make_tuple(uoDest, j, comp)];  // destination coupling DOF
						const unsigned int col = outletIndex + outletStride * comp;
						_jacActiveFN[uoSource].addElement(row, col, -ptrRate[idx] / _totalInletFlow(uoDest, j));
					}
				}
				else
				{
					const unsigned int row = _couplingIdxMap[std::make_tuple(uoDest, j, compDest)];  // destination coupling DOF
					const unsigned int col = outletIndex + outletStride * compSource;
					_jacActiveFN[uoSource].addElement(row, col, -ptrRate[idx] / _totalInletFlow(uoDest,j));
				}
			}
		}
		else
		{
			const unsigned int outletIndex = modelSource->localOutletComponentIndex(portSource);
			const unsigned int outletStride = modelSource->localOutletComponentStride(portSource);

			if (compSource == -1)
			{
				// Connect all components with the same flow rate
				for (unsigned int comp = 0; comp < modelSource->numComponents(); ++comp)
				{
					const unsigned int row = _couplingIdxMap[std::make_tuple(uoDest, portDest, comp)];  // destination coupling DOF
					const unsigned int col = outletIndex + outletStride * comp;
					_jacActiveFN[uoSource].addElement(row, col, -ptrRate[idx] / _totalInletFlow(uoDest, portDest));
				}
			}
			else
			{
				const unsigned int row = _couplingIdxMap[std::make_tuple(uoDest, portDest, compDest)];  // destination coupling DOF
				const unsigned int col = outletIndex + outletStride * compSource;
				_jacActiveFN[uoSource].addElement(row, col, -ptrRate[idx] / _totalInletFlow(uoDest, portDest));
			}
		}
	}

	// Copy active sparse matrices to their double pendants
	for (unsigned int i = 0; i < numModels(); ++i)
		_jacFN[i].copyFrom(_jacActiveFN[i]);
}

void ModelSystem::reportSolution(ISolutionRecorder& recorder, double const* const solution) const
{
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		m->reportSolution(recorder, solution + _dofOffset[i]);
	}
}

void ModelSystem::reportSolutionStructure(ISolutionRecorder& recorder) const
{
	for (IUnitOperation* m : _models)
		m->reportSolutionStructure(recorder);
}

unsigned int ModelSystem::requiredADdirs() const CADET_NOEXCEPT
{
	// Take maximum of required AD directions since each unit operation is (locally) independent from the rest
	unsigned int dirs = 0;
	for (IUnitOperation* m : _models)
		dirs = std::max(dirs, m->requiredADdirs());
	return dirs;
}

void ModelSystem::prepareADvectors(const AdJacobianParams& adJac) const
{
	// Early out if AD is disabled
	if (!adJac.adY)
		return;

	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		if (_models[i]->usesAD())
		{
			const unsigned int offset = _dofOffset[i];
			_models[i]->prepareADvectors(applyOffset(adJac, offset));
		}
	}
}

double ModelSystem::residualNorm(const SimulationTime& simTime, const ConstSimulationState& simState)
{
	residual(simTime, simState, _tempState);
	LOG(Debug) << "Residual: " << log::VectorPtr<double>(_tempState, numDofs());
	return linalg::linfNorm(_tempState, numDofs());
}

int ModelSystem::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res)
{
	BENCH_START(_timerResidual);

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(size_t(0), _models.size(), [&](size_t i)
#else
	for (unsigned int i = 0; i < _models.size(); ++i)
#endif
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];
		_errorIndicator[i] = m->residual(simTime, applyOffset(simState, offset), res + offset, _threadLocalStorage);
	} CADET_PARFOR_END;

	// Handle connections
	residualConnectUnitOps<double, double, double>(simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res);

	BENCH_STOP(_timerResidual);
	return totalErrorIndicatorFromLocal(_errorIndicator);
}

int ModelSystem::residualWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState,
	double* const res, const AdJacobianParams& adJac)
{
	BENCH_START(_timerResidual);

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(size_t(0), _models.size(), [&](size_t i)
#else
	for (unsigned int i = 0; i < _models.size(); ++i)
#endif
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];

		_errorIndicator[i] = m->residualWithJacobian(simTime, applyOffset(simState, offset),
			res + offset, applyOffset(adJac, offset), _threadLocalStorage);

	} CADET_PARFOR_END;

	// Handle connections
	residualConnectUnitOps<double, double, double>(simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res);

	BENCH_STOP(_timerResidual);
	return totalErrorIndicatorFromLocal(_errorIndicator);
}

/**
* @brief Calculate coupling DOF residual
* @param [in] secIdx  Section ID
* @param [in] y State vector
* @param [in] yDot Derivative state vector
* @param [in,out] res Residual vector
* @tparam StateType Type of the state vector
* @tparam ResidualType Type of the residual vector
* @tparam ParamType Type of the parameters
*/
template <typename StateType, typename ResidualType, typename ParamType>
void ModelSystem::residualConnectUnitOps(unsigned int secIdx, StateType const* const y, double const* const yDot, ResidualType* const res) CADET_NOEXCEPT
{
	// Use connection matrices for the residual
	const unsigned int finalOffset = _dofOffset.back();

	// N_f (Inlets to Inlets) Lower Right diagonal (Identity matrix)
	// The lower right matrix is Identity so residual equals y value
	for (unsigned int i = finalOffset; i < numDofs(); ++i)
		res[i] = y[i];

	// These could technically be done in parallel but from profiling no time is spent here 
	// and the parallelization has more overhead than can be gained.

	// N_{x,f} Inlets (Right) matrices; Right macro-column
	unsigned int offset;
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		offset = _dofOffset[i];
		_jacNF[i].multiplyAdd(y + finalOffset, res + offset);
	}

	// N_{f,x} Outlet (Lower) matrices; Bottom macro-row
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		offset = _dofOffset[i];
		select<ParamType>(_jacFN[i], _jacActiveFN[i]).multiplyAdd(y + offset, res + finalOffset);
	}
}

int ModelSystem::residualSensFwd(unsigned int nSens, const SimulationTime& simTime,
	const ConstSimulationState& simState, double const* const res,
	const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS,
	active* const adRes, double* const tmp1, double* const tmp2, double* const tmp3)
{
	BENCH_SCOPE(_timerResidualSens);
	return residualSensFwdWithJacobianAlgorithm<false>(nSens, simTime, simState, res, yS, ySdot, resS, AdJacobianParams{adRes, nullptr, 0}, tmp1, tmp2, tmp3);
}

void ModelSystem::multiplyWithMacroJacobian(double const* yS, double alpha, double beta, double* ret)
{
	const unsigned int finalOffset = _dofOffset.back();
	
	// Set ret_con = yS_con
	// This applies the identity matrix in the bottom right corner of the Jaocbian (network coupling equation)

	for (unsigned int i = finalOffset; i < numDofs(); ++i)
	{
		ret[i] = alpha * yS[i] + beta * ret[i];
	}

	// N_{x,f} Inlets (Right) matrices
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		const unsigned int offset = _dofOffset[i];
		_jacNF[i].multiplyAdd(yS + finalOffset, ret + offset, alpha);
	}

	// N_{f,x} Outlet (Lower) matrices
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		const unsigned int offset = _dofOffset[i];
		_jacFN[i].multiplyAdd(yS + offset, ret + finalOffset, alpha);
	}
}

void ModelSystem::residualSensFwdNorm(unsigned int nSens, const SimulationTime& simTime,
		const ConstSimulationState& simState,
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
	residualSensFwd(nSens, simTime, simState, nullptr, yS, ySdot, resPtr, adRes, tmp, tempMem.data(), tempMem.data() + nDOFs);

	// Calculate norms
	for (unsigned int i = 0; i < nSens; ++i)
		norms[i] = linalg::linfNorm(tempRes[i], nDOFs);
}

int ModelSystem::residualSensFwdWithJacobian(unsigned int nSens, const SimulationTime& simTime,
		const ConstSimulationState& simState, double const* const res,
		const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS,
		const AdJacobianParams& adJac, double* const tmp1, double* const tmp2, double* const tmp3)
{
	return residualSensFwdWithJacobianAlgorithm<true>(nSens, simTime, simState, res, yS, ySdot, resS, adJac, tmp1, tmp2, tmp3);
}

template <bool evalJacobian>
int ModelSystem::residualSensFwdWithJacobianAlgorithm(unsigned int nSens, const SimulationTime& simTime,
	const ConstSimulationState& simState, double const* const res,
	const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS,
	const AdJacobianParams& adJac, double* const tmp1, double* const tmp2, double* const tmp3)
{
	BENCH_START(_timerResidualSens);

	const unsigned int nModels = _models.size();

	//Resize yStemp and yStempDot (this should be a noop except for the first time)
	_yStemp.resize(nModels);
	_yStempDot.resize(nModels);
	_resSTemp.resize(nModels);

	for (unsigned int i = 0; i < nModels; ++i)
	{
		_yStemp[i].resize(yS.size());
		_yStempDot[i].resize(ySdot.size());
		_resSTemp[i].resize(resS.size());
	}

	// Step 1: Calculate sensitivities using AD in vector mode

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(size_t(0), size_t(nModels), [&](size_t i)
#else
	for (unsigned int i = 0; i < nModels; ++i)
#endif
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];

		_errorIndicator[i] = ResidualSensCaller<evalJacobian>::call(m, simTime, applyOffset(simState, offset), applyOffset(adJac, offset), _threadLocalStorage);
	} CADET_PARFOR_END;

	// Connect units
	residualConnectUnitOps<double, active, active>(simTime.secIdx, simState.vecStateY, simState.vecStateYdot, adJac.adRes);

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(size_t(0), size_t(nModels), [&](size_t i)
#else
	for (unsigned int i = 0; i < nModels; ++i)
#endif
	{
		// Step 2: Compute forward sensitivity residuals by multiplying with system Jacobians
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];

		// Move this outside the loop, these are memory addresses and should never change
		// Use correct offset in sensitivity state vectors
		for (unsigned int j = 0; j < yS.size(); ++j)
		{
			_yStemp[i][j] = yS[j] + offset;
			_yStempDot[i][j] = ySdot[j] + offset;
			_resSTemp[i][j] = resS[j] + offset;
		}

		const int intermediateRes = m->residualSensFwdCombine(simTime, applyOffset(simState, offset), _yStemp[i], _yStempDot[i], _resSTemp[i], adJac.adRes + offset, tmp1 + offset, tmp2 + offset, tmp3 + offset);
		_errorIndicator[i] = updateErrorIndicator(_errorIndicator[i], intermediateRes);
	} CADET_PARFOR_END;

	// tmp1 stores result of (dF / dy) * s
	// tmp2 stores result of (dF / dyDot) * sDot

	const unsigned int finalOffset = _dofOffset.back();

	// Handle super structure (i.e., right macro column and lower macro row)

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(size_t(0), yS.size(), [&](size_t param)
#else
	for (unsigned int param = 0; param < yS.size(); ++param)
#endif
	{
		double* const ptrResS = resS[param];

		// Directional derivative: res_{con} = (dF / dy) * s
		// Also adds contribution of the right macro column blocks
		multiplyWithMacroJacobian(yS[param], ptrResS);

		// Directional derivative (dF / dyDot) * sDot  (always zero so ignore it)

		//The other adRes values have already been taken care of in the unit operations
		for (unsigned int i = finalOffset; i < numDofs(); ++i)
		{
			ptrResS[i] += adJac.adRes[i].getADValue(param);
		}
	} CADET_PARFOR_END;

	BENCH_STOP(_timerResidualSens);
	return totalErrorIndicatorFromLocal(_errorIndicator);
}

int ModelSystem::dResDpFwdWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState,
	const AdJacobianParams& adJac)
{
	BENCH_SCOPE(_timerResidualSens);

	// Evaluate residual for all parameters using AD in vector mode and at the same time update the 
	// Jacobian (in one AD run, if analytic Jacobians are disabled)
#ifdef CADET_PARALLELIZE
	tbb::parallel_for(size_t(0), _models.size(), [&](size_t i)
#else
	for (unsigned int i = 0; i < _models.size(); ++i)
#endif
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];

		_errorIndicator[i] = m->residualSensFwdWithJacobian(simTime, applyOffset(simState, offset),
			applyOffset(adJac, offset), _threadLocalStorage);

	} CADET_PARFOR_END;

	// Handle connections
	residualConnectUnitOps<double, active, active>(simTime.secIdx, simState.vecStateY, simState.vecStateYdot, adJac.adRes);

	return totalErrorIndicatorFromLocal(_errorIndicator);
}

void ModelSystem::applyInitialCondition(const SimulationState& simState) const
{
	// If we have the full state vector available, use that and skip unit operations
	if (_initState.size() >= numDofs())
	{
		std::copy(_initState.data(), _initState.data() + numDofs(), simState.vecStateY);

		if (_initStateDot.size() >= numDofs())
			std::copy(_initStateDot.data(), _initStateDot.data() + numDofs(), simState.vecStateYdot);

		return;
	}

	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation const* const m = _models[i];
		const unsigned int offset = _dofOffset[i];
		m->applyInitialCondition(applyOffset(simState, offset));
	}
}

void ModelSystem::readInitialCondition(IParameterProvider& paramProvider)
{
	// Check if INIT_STATE_Y is present
	if (paramProvider.exists("INIT_STATE_Y"))
		_initState = paramProvider.getDoubleArray("INIT_STATE_Y");

	// Check if INIT_STATE_YDOT is present
	if (paramProvider.exists("INIT_STATE_YDOT"))
		_initStateDot = paramProvider.getDoubleArray("INIT_STATE_YDOT");

	std::ostringstream oss;
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];

		oss.str("");
		oss << "unit_" << std::setfill('0') << std::setw(3) << std::setprecision(0) << static_cast<int>(m->unitOperationId());

		const std::string subScope = oss.str();
		if (paramProvider.exists(subScope))
		{
			paramProvider.pushScope(subScope);
			m->readInitialCondition(paramProvider);
			paramProvider.popScope();
		}
	}
}

void ModelSystem::initializeSensitivityStates(const std::vector<double*>& vecSensY) const
{
	std::vector<double*> vecSensYlocal(vecSensY.size(), nullptr);
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];

		// Use correct offset in sensitivity state vectors
		for (unsigned int j = 0; j < vecSensY.size(); ++j)
			vecSensYlocal[j] = vecSensY[j] + offset;

		m->initializeSensitivityStates(vecSensYlocal);
	}
}

void ModelSystem::solveCouplingDOF(double* const vec)
{
	const unsigned int finalOffset = _dofOffset.back();

	// N_{f,x} Outlet (lower) matrices; Bottom macro-row
	// N_{f,x,1} * y_1 + ... + N_{f,x,nModels} * y_{nModels} + y_{coupling} = f
	// y_{coupling} = f - N_{f,x,1} * y_1 - ... - N_{f,x,nModels} * y_{nModels}
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		const unsigned int offset = _dofOffset[i];
		_jacFN[i].multiplySubtract(vec + offset, vec + finalOffset);
	}

	// Calculate inlet DOF for unit operations based on the coupling conditions. Depends on coupling conditions.
	// y_{unit op inlet} - y_{coupling} = 0
	// y_{unit op inlet} = y_{coupling}
	unsigned int idxCoupling = finalOffset;
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];
		if (!m->hasInlet())
			continue;
		
		for (unsigned int port = 0; port < m->numInletPorts(); ++port)
		{
			const unsigned int localIndex = m->localInletComponentIndex(port);
			const unsigned int localStride = m->localInletComponentStride(port);
			for (unsigned int comp = 0; comp < m->numComponents(); ++comp)
			{
				vec[offset + localIndex + comp*localStride] = vec[idxCoupling];
				++idxCoupling;
			}
		}
	}
}

template <typename tag_t>
void ModelSystem::consistentInitialConditionAlgorithm(const SimulationTime& simTime, const SimulationState& simState,
	const AdJacobianParams& adJac, double errorTol)
{
	BENCH_SCOPE(_timerConsistentInit);

	// Phase 1: Compute algebraic state variables

	// Consistent initial state for unit operations that only have outlets (system input, Inlet unit operation)
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];
		if (!m->hasInlet())
		{
			ConsistentInit<tag_t>::state(m, simTime, simState.vecStateY + offset, applyOffset(adJac, offset), errorTol, _threadLocalStorage);
		}
	}

	// Calculate coupling DOFs
	// These operations only requires correct unit operation outlet DOFs.
	// The outlets of the inlet unit operations have already been set above.
	// All other units are assumed to have correct outputs since their outlet DOFs are dynamic.
	const unsigned int finalOffset = _dofOffset.back();

	// Zero out the coupling DOFs (provides right hand side of 0 for solveCouplingDOF())
	std::fill(simState.vecStateY + finalOffset, simState.vecStateY + numDofs(), 0.0);

	solveCouplingDOF(simState.vecStateY);

	// Consistent initial state for all other unit operations (unit operations that have inlets)
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];
		if (m->hasInlet())
		{
			ConsistentInit<tag_t>::state(m, simTime, simState.vecStateY + offset, applyOffset(adJac, offset), errorTol, _threadLocalStorage);
		}
	}


	// Phase 2: Calculate residual with current state

	// Evaluate residual for right hand side without time derivatives \dot{y} and store it in vecStateYdot (or _tempState in case of lean initialization)
	// Also evaluate the Jacobian at the current position
	ConsistentInit<tag_t>::residualWithJacobian(*this, simTime, ConstSimulationState{simState.vecStateY, nullptr}, simState.vecStateYdot, _tempState, adJac);

	LOG(Debug) << "Residual post state: " << log::VectorPtr<double>(simState.vecStateYdot, numDofs());

	// Phase3 3: Calculate dynamic state variables yDot

	// Calculate all local yDot state variables
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];
		ConsistentInit<tag_t>::timeDerivative(m, simTime, simState.vecStateY + offset, simState.vecStateYdot + offset, _tempState + offset, _threadLocalStorage);
	}

	// Zero out the coupling DOFs (provides right hand side of 0 for solveCouplingDOF())
	std::fill(simState.vecStateYdot + finalOffset, simState.vecStateYdot + numDofs(), 0.0);
	// Calculate coupling DOFs
	solveCouplingDOF(simState.vecStateYdot);

	// Only enable this when you need to see the full jacobian for the system.
	// genJacobian(simTime, vecStateY, vecStateYdot);
}

void ModelSystem::consistentInitialConditions(const SimulationTime& simTime, const SimulationState& simState, 
	const AdJacobianParams& adJac, double errorTol)
{
	consistentInitialConditionAlgorithm<FullTag>(simTime, simState, adJac, errorTol);
}

void ModelSystem::consistentInitialSensitivity(const SimulationTime& simTime, 
	const ConstSimulationState& simState, std::vector<double*>& vecSensY, 
	std::vector<double*>& vecSensYdot, active* const adRes, active* const adY)
{
	consistentInitialSensitivityAlgorithm<FullTag>(simTime, simState, vecSensY, vecSensYdot, adRes, adY);
}

template <typename tag_t>
void ModelSystem::consistentInitialSensitivityAlgorithm(const SimulationTime& simTime, 
	const ConstSimulationState& simState, std::vector<double*>& vecSensY, 
	std::vector<double*>& vecSensYdot, active* const adRes, active* const adY)
{
	BENCH_SCOPE(_timerConsistentInit);

	// Compute parameter sensitivities and update the Jacobian
	dResDpFwdWithJacobian(simTime, simState, AdJacobianParams{adRes, adY, static_cast<unsigned int>(vecSensY.size())});

	std::vector<double*> vecSensYlocal(vecSensY.size(), nullptr);
	std::vector<double*> vecSensYdotLocal(vecSensYdot.size(), nullptr);
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];

		if (!m->hasInlet())
		{
			const unsigned int offset = _dofOffset[i];

			// Use correct offset in sensitivity state vectors
			for (unsigned int j = 0; j < vecSensY.size(); ++j)
			{
				vecSensYlocal[j] = vecSensY[j] + offset;
				vecSensYdotLocal[j] = vecSensYdot[j] + offset;
			}

			ConsistentInit<tag_t>::parameterSensitivity(m, simTime, applyOffset(simState, offset), vecSensYlocal, vecSensYdotLocal, adRes + offset, _threadLocalStorage);
		}
	}

	const unsigned int finalOffset = _dofOffset.back();

	for (unsigned int param = 0; param < vecSensY.size(); ++param)
	{
		double* const vsy = vecSensY[param];
		for (unsigned int i = finalOffset; i < numDofs(); ++i)
			vsy[i] = -adRes[i].getADValue(param);

		solveCouplingDOF(vsy);
	}

	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		IUnitOperation* const m = _models[i];

		if (m->hasInlet())
		{
			const unsigned int offset = _dofOffset[i];

			// Use correct offset in sensitivity state vectors
			for (unsigned int j = 0; j < vecSensY.size(); ++j)
			{
				vecSensYlocal[j] = vecSensY[j] + offset;
				vecSensYdotLocal[j] = vecSensYdot[j] + offset;
			}

			ConsistentInit<tag_t>::parameterSensitivity(m, simTime, applyOffset(simState, offset), vecSensYlocal, vecSensYdotLocal, adRes + offset, _threadLocalStorage);
		}
	}
		
	for (unsigned int i = 0; i < vecSensY.size(); ++i)
	{
		double* const vsyd = vecSensYdot[i];

		// Calculate -(d^2 res_con / (dy dp)) * \dot{y}
		if (_models.empty())
		{
			std::fill(vsyd + finalOffset, vsyd + numDofs(), 0.0);
		}
		else
		{
			ad::adMatrixVectorMultiply(_jacActiveFN[0], simState.vecStateYdot + _dofOffset[0], vsyd + finalOffset, -1.0, 0.0, i);
			for (unsigned int j = 1; j < _models.size(); ++j)
			{
				const unsigned int offset = _dofOffset[j];
				ad::adMatrixVectorMultiply(_jacActiveFN[j], simState.vecStateYdot + offset, vsyd + finalOffset, -1.0, 1.0, i);
			}
		}
		solveCouplingDOF(vsyd);
	}
}

void ModelSystem::leanConsistentInitialConditions(const SimulationTime& simTime, const SimulationState& simState, 
	const AdJacobianParams& adJac, double errorTol)
{
	consistentInitialConditionAlgorithm<LeanTag>(simTime, simState, adJac, errorTol);
}

void ModelSystem::leanConsistentInitialSensitivity(const SimulationTime& simTime, const ConstSimulationState& simState,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active* const adRes, active* const adY)
{
	consistentInitialSensitivityAlgorithm<LeanTag>(simTime, simState, vecSensY, vecSensYdot, adRes, adY);
}

/**
 * @brief Multiplies a vector with the full Jacobian of the entire system (i.e., @f$ \frac{\partial F}{\partial y}\left(t, y, \dot{y}\right) @f$)
 * @details Actually, the operation @f$ z = \alpha \frac{\partial F}{\partial y} x + \beta z @f$ is performed. 
 * @param [in] simTime Current simulation time point
 * @param [in] simState Simulation state vectors
 * @param [in] yS Vector @f$ x @f$ that is transformed by the Jacobian @f$ \frac{\partial F}{\partial y} @f$
 * @param [in] alpha Factor @f$ \alpha @f$ in front of @f$ \frac{\partial F}{\partial y} @f$
 * @param [in] beta Factor @f$ \beta @f$ in front of @f$ z @f$
 * @param [in,out] ret Vector @f$ z @f$ which stores the result of the operation
 */
void ModelSystem::multiplyWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* yS, double alpha, double beta, double* ret)
{
	for (unsigned int idxModel = 0; idxModel < _models.size(); ++idxModel)
	{
		IUnitOperation* const m = _models[idxModel];
		const unsigned int offset = _dofOffset[idxModel];
		m->multiplyWithJacobian(simTime, simState, yS + offset, alpha, beta, ret + offset);
	}
	multiplyWithMacroJacobian(yS, alpha, beta, ret);
}

/**
 * @brief Multiplies a vector with the full time derivative Jacobian of the entire system (i.e., @f$ \frac{\partial F}{\partial \dot{y}}\left(t, y, \dot{y}\right) @f$)
 * @details The operation @f$ z = \frac{\partial F}{\partial \dot{y}} x @f$ is performed. 
 * @param [in] simTime Current simulation time point
 * @param [in] simState Simulation state vectors
 * @param [in] yS Vector @f$ x @f$ that is transformed by the Jacobian @f$ \frac{\partial F}{\partial \dot{y}} @f$
 * @param [in,out] ret Vector @f$ z @f$ which stores the result of the operation
 */
void ModelSystem::multiplyWithDerivativeJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* yS, double* ret)
{
	for (unsigned int idxModel = 0; idxModel < _models.size(); ++idxModel)
	{
		IUnitOperation* const m = _models[idxModel];
		const unsigned int offset = _dofOffset[idxModel];
		m->multiplyWithDerivativeJacobian(simTime, simState, yS + offset, ret + offset);
	}
	std::fill(ret + _dofOffset.back(), ret + numDofs(), 0.0);
}

/**
 * @brief Generate full system Jacobian FD and multiplyWithJacobian
 * @details During debugging this allows you to generate the full jacobian and verify the jacobian structure
 *          is what it should be. The system uses FD and multiplyWithJacobian to create the full jacobian.
 *          Use this function with a debugger and pull the values out of memory to visualize it.
 *
 * @param [in] simTime Current simulation time point
 * @param [in] simState Simulation state vectors
 */
void ModelSystem::genJacobian(const SimulationTime& simTime, const ConstSimulationState& simState)
{
	// This method is only for debugging. No point in optimizing it
	const unsigned int size = numDofs();

	// Jacobians are saved in column-major ordering (i.e., each column is added to the array sequentially / columns are stacked together)
	std::vector<double> jacobian(size*size, 0.0);
	std::vector<double> jacobianDot(size*size, 0.0);

	std::vector<double> jacobianFD(size*size, 0.0);
	std::vector<double> jacobianFDDot(size*size, 0.0);

	const double h = 1e-5;

	std::vector<double> f(size, 0.0);
	std::vector<double> fdot(size, 0.0);
	std::vector<double> fh(size, 0.0);
	std::vector<double> fhdot(size, 0.0);

	std::vector<double> res(size, 0.0);
	std::vector<double> resh(size, 0.0);

	// create Jacobian
	for (unsigned int i = 0; i < size; ++i)
	{
		// Clear res and resh
		std::fill(res.begin(), res.end(), 0.0);
		std::fill(resh.begin(), resh.end(), 0.0);
		
		// Copy y and yDot
		std::copy_n(simState.vecStateY, size, &f[0]);
		std::copy_n(simState.vecStateY, size, &fh[0]);

		std::copy_n(simState.vecStateYdot, size, &fdot[0]);
		std::copy_n(simState.vecStateYdot, size, &fhdot[0]);

		// Change ith entry
		double stepSize = h;
		if (f[i] != 0.0)
			stepSize = f[i] * h;

		f[i] -= stepSize / 2;
		fh[i] += stepSize / 2;

		residual(simTime, ConstSimulationState{&f[0], &fdot[0]}, &res[0]);
		residual(simTime, ConstSimulationState{&fh[0], &fhdot[0]}, &resh[0]);

		for (unsigned int j = 0; j < size; ++j)
		{
			jacobianFD[i*size + j] = (resh[j] - res[j]) / stepSize;
		}
	}
	
	// create JacobianDot
	for (unsigned int i = 0; i < size; ++i)
	{
		// Clear res and resh
		std::fill(res.begin(), res.end(), 0.0);
		std::fill(resh.begin(), resh.end(), 0.0);

		// Copy y and yDot
		std::copy_n(simState.vecStateY, size, &f[0]);
		std::copy_n(simState.vecStateY, size, &fh[0]);

		std::copy_n(simState.vecStateYdot, size, &fdot[0]);
		std::copy_n(simState.vecStateYdot, size, &fhdot[0]);

		// Change ith entry
		double stepSize = h;
		if (fdot[i] != 0.0)
			stepSize = fdot[i] * h;

		fdot[i] -= stepSize / 2;
		fhdot[i] += stepSize / 2;

		residual(simTime, ConstSimulationState{&f[0], &fdot[0]}, &res[0]);
		residual(simTime, ConstSimulationState{&fh[0], &fhdot[0]}, &resh[0]);

		for (unsigned int j = 0; j < size; ++j)
		{
			jacobianFDDot[i*size + j] = (resh[j] - res[j]) / stepSize;
		}
	}

	std::vector<double> unit(size, 0.0);

	for (unsigned int i = 0; i < size; ++i)
	{
		std::fill(res.begin(), res.end(), 0.0);
		// Clear res and resh
		unit[i] = 1.0;

		multiplyWithJacobian(simTime, simState, unit.data(), 1.0, 0.0, res.data());
		std::copy(res.begin(), res.end(), jacobian.begin() + i * size);

		unit[i] = 0.0;
	}

	for (unsigned int i = 0; i < size; ++i)
	{
		std::fill(res.begin(), res.end(), 0.0);
		// Clear res and resh
		unit[i] = 1.0;

		multiplyWithDerivativeJacobian(simTime, simState, unit.data(), res.data());
		std::copy(res.begin(), res.end(), jacobianDot.begin() + i * size);

		unit[i] = 0.0;
	}

	LOG(Debug) << "jacFD = " << log::MatrixPtr<double>(jacobianFD.data(), size, size, true);
	LOG(Debug) << "jacFDDot = " << log::MatrixPtr<double>(jacobianFDDot.data(), size, size, true);
	LOG(Debug) << "jac = " << log::MatrixPtr<double>(jacobian.data(), size, size, true);
	LOG(Debug) << "jacDot = " << log::MatrixPtr<double>(jacobianDot.data(), size, size, true);
}

/**
 * @brief Generate full system Jacobian with Sensitivities using FD and multiplyWithJacobian
 * @details During debugging this allows you to generate the full sensitivity jacobian and verify the jacobian structure
 *          is what it should be. The system uses FD and multiplyWithJacobian to create the full jacobian.
 *          Use this function with a debugger and pull the values out of memory to visualize it.
 *
 * @param [in] t Current time point
 * @param [in] simTime Current simulation time point
 * @param [in] simState Simulation state vectors
 * @param [in] residual vector
 * @param [in] yS Sensitivity State Vector
 * @param [in] ySdot Sensitivity State Vector
 * @param [in] resS Sensitivity residual vector
 * @param [in] adRes
 * @param [in] tmp1
 * @param [in] tmp2
 * @param [in] tmp3
 */
void ModelSystem::genJacobian(unsigned int nSens, const SimulationTime& simTime,
	const ConstSimulationState& simState, double const* const res,
	const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS,
	active* const adRes, double* const tmp1, double* const tmp2, double* const tmp3)
{
	// This method is only for debugging. Don't bother optimizing it
	const unsigned int size = numDofs();

	// Jacobians are saved in column-major ordering (i.e., each column is added to the array sequentially / columns are stacked together)
	std::vector<std::vector<double>> jacobianFD(nSens, std::vector<double>(size*size));
	std::vector<std::vector<double>> jacobianFDDot(nSens, std::vector<double>(size*size));

	const double h = 1e-5;

	// -h/2
	std::vector<double> tmp1mh(size, 0.0);
	std::vector<double> tmp2mh(size, 0.0);
	std::vector<double> tmp3mh(size, 0.0);

	//  h/2
	std::vector<double> tmp1ph(size, 0.0);
	std::vector<double> tmp2ph(size, 0.0);
	std::vector<double> tmp3ph(size, 0.0);

	std::vector<active> adResmh(size, 0.0);
	std::vector<active> adResph(size, 0.0);

	std::vector<double *> ySmh(nSens);
	std::vector<double *> ySdotmh(nSens);
	std::vector<double *> resSmh(nSens);

	std::vector<double *> ySph(nSens);
	std::vector<double *> ySdotph(nSens);
	std::vector<double *> resSph(nSens);

	std::vector<const double *> CySmh(nSens);
	std::vector<const double *> CySdotmh(nSens);

	std::vector<const double *> CySph(nSens);
	std::vector<const double *> CySdotph(nSens);

	// Allocate memory
	for (unsigned int j = 0; j < nSens; ++j)
	{
		ySmh[j] = new double[size];
		ySdotmh[j] = new double[size];
		resSmh[j] = new double[size];

		ySph[j] = new double[size];
		ySdotph[j] = new double[size];
		resSph[j] = new double[size];
	}


	for (unsigned int j = 0; j < nSens; ++j)
	{
		CySmh[j] = ySmh[j];
		CySdotmh[j] = ySdotmh[j];
		CySph[j] = ySph[j];
		CySdotph[j] = ySdotph[j];
	}
	
	// create Jacobian
	for (unsigned int i = 0; i < size; ++i)
	{
		// need to make copies of yS, ySdot, resS, adRes, tmp1, tmp2, tmp3

		// adRes
		std::copy_n(adRes, size, &adResmh[0]);
		std::copy_n(adRes, size, &adResph[0]);

		// tmp1
		std::copy_n(tmp1, size, &tmp1mh[0]);
		std::copy_n(tmp1, size, &tmp1ph[0]);

		// tmp2
		std::copy_n(tmp2, size, &tmp2mh[0]);
		std::copy_n(tmp2, size, &tmp2ph[0]);

		// tmp3
		std::copy_n(tmp3, size, &tmp3mh[0]);
		std::copy_n(tmp3, size, &tmp3ph[0]);
		
		// Clear sync up
		for (unsigned int j = 0; j < nSens; ++j)
		{
			std::copy_n(yS[j], size, ySmh[j]);
			std::copy_n(yS[j], size, ySph[j]);

			std::copy_n(ySdot[j], size, ySdotmh[j]);
			std::copy_n(ySdot[j], size, ySdotph[j]);

			std::copy_n(resS[j], size, resSmh[j]);
			std::copy_n(resS[j], size, resSph[j]);
		}

		std::vector<double> stepSize(nSens, false);

		// Change ith entry
		for (unsigned int j = 0; j < nSens; ++j)
		{
			const double val = ySmh[j][i];
			if (val == 0.0)
			{
				ySmh[j][i] -= h / 2;
				ySph[j][i] += h / 2;
				stepSize[j] = h;
			}
			else
			{
				ySmh[j][i] -= val * h / 2;
				ySph[j][i] += val * h / 2;
				stepSize[j] = val * h;
			}
		}

		// clear jacobian

		// -h/2
		residualSensFwd(nSens, simTime, simState, res, CySmh, CySdotmh, resSmh, &adResmh[0], &tmp1mh[0], &tmp2mh[0], &tmp3mh[0]);

		// +h/2
		residualSensFwd(nSens, simTime, simState, res, CySph, CySdotph, resSph, &adResph[0], &tmp1ph[0], &tmp2ph[0], &tmp3ph[0]);

		for (unsigned int sens = 0; sens < nSens; ++sens)
		{
			for (unsigned int j = 0; j < size; ++j)
			{
				// Residual is negative so it has to be negated to get the correct jacobian
				jacobianFD[sens][i*size + j] = (resSph[sens][j] - resSmh[sens][j]) / stepSize[sens];
			}
		}
	}

	//create Jacobian
	for (unsigned int i = 0; i < size; ++i)
	{
		// need to make copies of yS, ySdot, resS, adRes, tmp1, tmp2, tmp3

		// adRes
		std::copy_n(adRes, size, &adResmh[0]);
		std::copy_n(adRes, size, &adResph[0]);

		// tmp1
		std::copy_n(tmp1, size, &tmp1mh[0]);
		std::copy_n(tmp1, size, &tmp1ph[0]);

		// tmp2
		std::copy_n(tmp2, size, &tmp2mh[0]);
		std::copy_n(tmp2, size, &tmp2ph[0]);

		// tmp3
		std::copy_n(tmp3, size, &tmp3mh[0]);
		std::copy_n(tmp3, size, &tmp3ph[0]);

		// Clear sync up
		for (unsigned int j = 0; j < nSens; ++j)
		{
			std::copy_n(yS[j], size, ySmh[j]);
			std::copy_n(yS[j], size, ySph[j]);

			std::copy_n(ySdot[j], size, ySdotmh[j]);
			std::copy_n(ySdot[j], size, ySdotph[j]);

			std::copy_n(resS[j], size, resSmh[j]);
			std::copy_n(resS[j], size, resSph[j]);
		}

		std::vector<double> stepSize(nSens, false);

		// Change ith entry
		for (unsigned int j = 0; j < nSens; ++j)
		{
			const double val = ySdotmh[j][i];
			if (val == 0.0)
			{
				ySdotmh[j][i] -= h / 2;
				ySdotph[j][i] += h / 2;
				stepSize[j] = h;
			}
			else
			{
				ySdotmh[j][i] -= val * h / 2;
				ySdotph[j][i] += val * h / 2;
				stepSize[j] = val * h;
			}
		}

		// clear jacobian

		// -h/2
		residualSensFwd(nSens, simTime, simState, res, CySmh, CySdotmh, resSmh, &adResmh[0], &tmp1mh[0], &tmp2mh[0], &tmp3mh[0]);

		// +h/2
		residualSensFwd(nSens, simTime, simState, res, CySph, CySdotph, resSph, &adResph[0], &tmp1ph[0], &tmp2ph[0], &tmp3ph[0]);

		for (unsigned int sens = 0; sens < nSens; ++sens)
		{
			for (unsigned int j = 0; j < size; ++j)
			{
				//Residual is negative so it has to be negated to get the correct jacobian
				jacobianFDDot[sens][i*size + j] = (resSph[sens][j] - resSmh[sens][j]) / stepSize[sens];
			}
		}
	}

	// Free memory 

	for (unsigned int j = 0; j < nSens; ++j)
	{
		delete ySmh[j];
		delete ySdotmh[j];
		delete resSmh[j];

		delete ySph[j];
		delete ySdotph[j];
		delete resSph[j];
	}

	for (unsigned int i = 0; i < nSens; ++i)
	{
		LOG(Debug) << "jacSens" << i << " = " << log::MatrixPtr<double>(jacobianFD[i].data(), size, size, true);
		LOG(Debug) << "jacSensDot" << i << " = " << log::MatrixPtr<double>(jacobianFDDot[i].data(), size, size, true);
	}
}

int ModelSystem::linearSolve(double t, double alpha, double outerTol, double* const rhs, double const* const weight,
	const ConstSimulationState& simState)
{
	// TODO: Add early out error checks

	BENCH_SCOPE(_timerLinearSolve);

	const unsigned int finalOffset = _dofOffset[_models.size()];

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(size_t(0), _models.size(), [=](size_t i)
#else
	for (unsigned int i = 0; i < _models.size(); ++i)
#endif
	{
		IUnitOperation* const m = _models[i];
		const unsigned int offset = _dofOffset[i];
		_errorIndicator[i] = m->linearSolve(t, alpha, outerTol, rhs + offset, weight + offset, applyOffset(simState, offset));
	} CADET_PARFOR_END;

	// Solve last row of L with backwards substitution: y_f = b_f - \sum_{i=0}^{N_z} J_{f,i} y_i
	// Note that we cannot easily parallelize this loop since the results of the sparse
	// matrix-vector multiplications are added in-place to rhs. We would need one copy of rhs
	// for each thread and later fuse them together (reduction statement).
	for (unsigned int i = 0; i < _models.size(); ++i)
	{
		const unsigned int offset = _dofOffset[i];
		_jacFN[i].multiplySubtract(rhs + offset, rhs + finalOffset);
	}


	// Now, rhs contains the full intermediate solution y = L^{-1} b

	// Initialize temporary storage by copying over the fluxes
	std::fill_n(_tempState, finalOffset, 0.0);
	std::copy_n(rhs + finalOffset, numCouplingDOF(), _tempState + finalOffset);


	// ==== Step 3: Solve Schur-complement to get x_f = S^{-1} y_f
	// Column and particle parts remain unchanged.
	// The only thing to be done is the iterative (and approximate)
	// solution of the Schur complement system:
	//     S * x_f = y_f

	// Note that rhs is updated in-place with the solution of the Schur-complement
	// The temporary storage is only needed to hold the right hand side of the Schur-complement
	const double tolerance = std::sqrt(static_cast<double>(numDofs())) * outerTol * _schurSafety;

	// The network version of the schurCompletmentMatrixVector function need access to more information than the current interface
	// Instead of changing the interface a lambda function is used and closed over the additional variables
	auto schurComplementMatrixVectorPartial = [&, this](void* userData, double const* x, double* z) -> int 
	{
		return ModelSystem::schurComplementMatrixVector(x, z, t, alpha, outerTol, weight, simState);
	};

	_gmres.matrixVectorMultiplier(schurComplementMatrixVectorPartial);
	
	// Reset error indicator as it is used in schurComplementMatrixVector()
	const int curError = totalErrorIndicatorFromLocal(_errorIndicator);
	std::fill(_errorIndicator.begin(), _errorIndicator.end(), 0);

	const int gmresResult = _gmres.solve(tolerance, weight + finalOffset, _tempState + finalOffset, rhs + finalOffset);

	// Set last cumulative error to all elements to restore state (in the end only total error matters)
	std::fill(_errorIndicator.begin(), _errorIndicator.end(), updateErrorIndicator(curError, gmresResult));

	// Reset temporary memory
	std::fill_n(_tempState, finalOffset, 0.0);

	// At this point, rhs contains the intermediate solution [y_0, ..., y_{N_z}, x_f]

	// ==== Step 4: Solve U * x = y by backward substitution
	// The fluxes are already solved and remain unchanged
#ifdef CADET_PARALLELIZE
	tbb::parallel_for(size_t(0), _models.size(), [=](size_t idxModel)
#else
	for (unsigned int idxModel = 0; idxModel < _models.size(); ++idxModel)
#endif
	{
		IUnitOperation* const m = _models[idxModel];
		const unsigned int offset = _dofOffset[idxModel];

		// Compute tempState_i = N_{i,f} * y_f
		_jacNF[idxModel].multiplyVector(rhs + finalOffset, _tempState + offset);

		// Apply N_i^{-1} to tempState_i
		const int linSolve = m->linearSolve(t, alpha, outerTol, _tempState + offset, weight + offset, applyOffset(simState, offset));
		_errorIndicator[idxModel] = updateErrorIndicator(_errorIndicator[idxModel], linSolve);

		// Compute rhs_i = y_i - N_i^{-1} * N_{i,f} * y_f = y_i - tempState_i
		const unsigned int offsetNext = _dofOffset[idxModel + 1];
		for (unsigned int i = offset; i < offsetNext; ++i)
		{
			rhs[i] -= _tempState[i];
		}
	} CADET_PARFOR_END;

	return totalErrorIndicatorFromLocal(_errorIndicator);
}

/**
* @brief Performs the matrix-vector product @f$ z = Sx @f$ with the Schur-complement @f$ S @f$ from the Jacobian
* @details The Schur-complement @f$ S @f$ is given by
*          @f[ \begin{align}
S &= J_f - J_{f,0} \, J_0^{-1} \, J_{0,f} - \sum_{p=1}^{N_z}{J_{f,p} \, J_p^{-1} \, J_{p,f}} \\
&= I - \sum_{p=0}^{N_z}{J_{f,p} \, J_p^{-1} \, J_{p,f}},
\end{align} @f]
*          where @f$ J_f = I @f$ is the identity matrix and the off-diagonal blocks @f$ J_{i,f} @f$
*          and @f$ J_{f,i} @f$ for @f$ i = 0, \dots, N_{z} @f$ are sparse.
*
*          The matrix-vector multiplication is executed in parallel as follows:
*              -# Compute @f$ J_{f,i} \, J_i^{-1} \, J_{i,f} @f$ independently (in parallel with respect to index @f$ i @f$)
*              -# Subtract the result from @f$ z @f$ in a critical section to avoid race conditions
*
* @param [in] x Vector @f$ x @f$ the matrix @f$ S @f$ is multiplied with
* @param [out] z Result of the matrix-vector multiplication
* @return @c 0 if successful, any other value in case of failure
*/
int ModelSystem::schurComplementMatrixVector(double const* x, double* z, double t, double alpha, double outerTol, double const* const weight,
	const ConstSimulationState& simState) const
{
	BENCH_SCOPE(_timerMatVec);

	// Copy x over to result z, which corresponds to the application of the identity matrix
	std::copy(x, x + numCouplingDOF(), z);

	// Inlets and outlets don't participate in the Schur solver since one of NF or FN for them is always 0
	// As a result we only have to work with items that have both an inlet and an outlet
#ifdef CADET_PARALLELIZE
	tbb::parallel_for(size_t(0), _inOutModels.size(), [=](size_t i)
#else
	for (unsigned int i = 0; i < _inOutModels.size(); ++i)
#endif
	{
		const unsigned int idxModel = _inOutModels[i];
		IUnitOperation* const m = _models[idxModel];
		const unsigned int offset = _dofOffset[idxModel];

		_jacNF[idxModel].multiplyVector(x, _tempState + offset);

		// Apply N_i^{-1} to tempState_i
		const int linSolve = m->linearSolve(t, alpha, outerTol, _tempState + offset, weight + offset, applyOffset(simState, offset));
		_errorIndicator[idxModel] = updateErrorIndicator(_errorIndicator[idxModel], linSolve);

		// Apply J_{f,i} and subtract results from z
		{
#ifdef CADET_PARALLELIZE
			SchurComplementMutex::scoped_lock l(_schurMutex);
#endif
			_jacFN[idxModel].multiplySubtract(_tempState + offset, z);
		}
	} CADET_PARFOR_END;

	return totalErrorIndicatorFromLocal(_errorIndicator);
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

void ModelSystem::setupParallelization(unsigned int numThreads)
{
	unsigned int tlsSize = 0;
	for (IUnitOperation const* m : _models)
		tlsSize = std::max(tlsSize, m->threadLocalMemorySize());

	_threadLocalStorage.resize(numThreads, tlsSize);
}

}  // namespace model

}  // namespace cadet
