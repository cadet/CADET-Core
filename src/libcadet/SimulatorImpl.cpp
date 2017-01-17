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

#include "cadet/Exceptions.hpp"
#include "cadet/SolutionRecorder.hpp"
#include "cadet/ParameterProvider.hpp"
#include "SimulatorImpl.hpp"
#include "SimulatableModel.hpp"
#include "UnitOperation.hpp"
#include "model/ModelSystemImpl.hpp"
#include "ParamIdUtil.hpp"

#include <idas/idas.h>
#include <idas/idas_impl.h>
#include "SundialsVector.hpp"

#include <vector>
#include <sstream>

#include "AutoDiff.hpp"
#include "LoggingUtils.hpp"
#include "Logging.hpp"

#ifdef _OPENMP
	#include <omp.h>
#endif

namespace
{
	template <class T>
	const std::vector<T> convertNVectorToStdVectorPtrs(N_Vector* vec, unsigned int numVec)
	{
		std::vector<T> sensState(numVec, nullptr);
		for (unsigned int i = 0; i < numVec; ++i)
			sensState[i] = NVEC_DATA(vec[i]);

		return sensState;
	}

	const std::vector<double*> convertNVectorToStdVectorPtrs(unsigned int& len, N_Vector* vec, unsigned int numVec)
	{
		len = NVEC_LENGTH(vec[0]);
		return convertNVectorToStdVectorPtrs<double*>(vec, numVec);
	}

	const std::vector<const double*> convertNVectorToStdVectorConstPtrs(unsigned int& len, N_Vector* vec, unsigned int numVec)
	{
		len = NVEC_LENGTH(vec[0]);
		return convertNVectorToStdVectorPtrs<const double*>(vec, numVec);
	}

	/**
	 * @brief Checks whether a given parameter @p id corresponds to a SECTION_TIMES parameter
	 * @param [in] id Parameter id to be checked
	 * @param [in] nSectionTimes Number of sections
	 * @return @c true if the given id corresponds to a SECTION_TIMES parameter, otherwise @c false
	 */
	inline bool isSectionTimeParameter(const cadet::ParameterId& id, unsigned int nSectionTimes)
	{
		return ((id.name == cadet::hashString("SECTION_TIMES")) && (id.section < nSectionTimes) && (id.component == cadet::CompIndep) &&
			(id.boundPhase == cadet::BoundPhaseIndep) && (id.reaction == cadet::ReactionIndep) && (id.unitOperation == cadet::UnitOpIndep));
	}
}

namespace cadet
{
	namespace log
	{
		inline std::ostream& operator<<(std::ostream& os, const N_Vector& nv)
		{
			double const* const ptrY = NVEC_DATA(nv);
			os << "[";
			for (unsigned int i = 0; i < NVEC_LENGTH(nv)-1; ++i)
				os << ptrY[i] << ",";
			os << ptrY[NVEC_LENGTH(nv)-1] << "]";
			return os;
		}
	}

	/**
	 * @brief IDAS error handler function
	 * @details Handles errors reported by the IDAS solver. See section 4.6.2 of the IDAS manual for details.
	 */
	void idasErrorHandler(int error_code, const char* module, const char* function, char* msg, void* eh_data)
	{
		cadet::Simulator* const sim = static_cast<cadet::Simulator*>(eh_data);

		std::ostringstream oss;
		oss << "In function '" << function << "' of module '" << module << "', error code '" << IDAGetReturnFlagName(error_code) << "':\n" << msg;

		// @todo Find an error handling system and put it here
		if (error_code < 0) 
		{
			// Error
			LOG(Error) << oss.str();
		}
		else
		{
			// Warning
			LOG(Warning) << oss.str();
		}
	}

	int residualDaeWrapper(double t, N_Vector y, N_Vector yDot, N_Vector res, void* userData)
	{
		cadet::Simulator* const sim = static_cast<cadet::Simulator*>(userData);
		const unsigned int secIdx = sim->getCurrentSection(t);
		const active timeFactor = sim->timeFactor();
		const int retVal = sim->_model->residualWithJacobian(sim->toRealTime(t), secIdx, timeFactor, NVEC_DATA(y), NVEC_DATA(yDot), NVEC_DATA(res), 
			sim->_vecADres, sim->_vecADy, sim->numSensitivityAdDirections());

		return retVal;
	}

	int linearSolveWrapper(IDAMem IDA_mem, N_Vector rhs, N_Vector weight, N_Vector y, N_Vector yDot, N_Vector res)
	{
		cadet::Simulator* const sim = static_cast<cadet::Simulator*>(IDA_mem->ida_lmem);
		const double t = static_cast<double>(sim->toRealTime(IDA_mem->ida_tn));
		const double alpha = IDA_mem->ida_cj;
		const double tol = IDA_mem->ida_epsNewt;
		const active timeFactor = sim->timeFactor();

		const int retVal =  sim->_model->linearSolve(t, static_cast<double>(timeFactor), alpha, tol, NVEC_DATA(rhs), NVEC_DATA(weight), NVEC_DATA(y), NVEC_DATA(yDot), NVEC_DATA(res));
		return retVal;
	}

	int residualSensWrapper(int ns, double t, N_Vector y, N_Vector yDot, N_Vector res, 
			N_Vector* yS, N_Vector* ySDot, N_Vector* resS,
			void *userData, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
	{
		cadet::Simulator* const sim = static_cast<cadet::Simulator*>(userData);
		const std::vector<const double*> sensY = convertNVectorToStdVectorPtrs<const double*>(yS, ns);
		const std::vector<const double*> sensYdot = convertNVectorToStdVectorPtrs<const double*>(ySDot, ns);
		std::vector<double*> sensRes = convertNVectorToStdVectorPtrs<double*>(resS, ns);
		const unsigned int secIdx = sim->getCurrentSection(t);
		const active timeFactor = sim->timeFactor();

		return sim->_model->residualSensFwd(ns, sim->toRealTime(t), secIdx, timeFactor, NVEC_DATA(y), NVEC_DATA(yDot), NVEC_DATA(res), 
			sensY, sensYdot, sensRes, sim->_vecADres, NVEC_DATA(tmp1), NVEC_DATA(tmp2), NVEC_DATA(tmp3));
	}

	Simulator::Simulator() : _model(nullptr), _solRecorder(nullptr), _idaMemBlock(nullptr), _vecStateY(nullptr), 
		_vecStateYdot(nullptr), _vecFwdYs(nullptr), _vecFwdYsDot(nullptr),
		_relTolS(1.0e-9), _absTol(1, 1.0e-12), _relTol(1.0e-9), _initStepSize(1, 1.0e-6), _maxSteps(10000), _curSec(0),
		_skipConsistencyStateY(false), _skipConsistencySensitivity(false), _consistentInitMode(ConsistentInitialization::Full), 
		_consistentInitModeSens(ConsistentInitialization::Full), _vecADres(nullptr), _vecADy(nullptr), _lastIntTime(0.0)
	{
#if defined(ACTIVE_ADOLC) || defined(ACTIVE_SFAD) || defined(ACTIVE_SETFAD)
		LOG(Debug) << "Resetting AD directions from " << ad::getDirections() << " to default " << SFAD_DEFAULT_DIR;
		ad::setDirections(SFAD_DEFAULT_DIR);
#endif		
	}

	Simulator::~Simulator() CADET_NOEXCEPT
	{
		clearModel();
	}

	void Simulator::clearModel() CADET_NOEXCEPT
	{
		delete[] _vecADy;
		delete[] _vecADres;

		if ((_sensitiveParams.slices() > 0) && _vecFwdYs)
		{
			NVec_DestroyArray(_vecFwdYs, _sensitiveParams.slices());
			NVec_DestroyArray(_vecFwdYsDot, _sensitiveParams.slices());
		}
		_sensitiveParams.clear();
		
		if (_vecStateYdot)
			NVec_Destroy(_vecStateYdot);
		if (_vecStateY)
			NVec_Destroy(_vecStateY);

		if (_idaMemBlock)
			IDAFree(&_idaMemBlock);		
	}

	void Simulator::initializeModel(IModelSystem& model)
	{
		// Clean up
		clearModel();

		// We only provide model::ModelSystem as implementation
		_model = static_cast<model::ModelSystem*>(&model);

		// Allocate and initialize state vectors
		const unsigned int nDOFs = _model->numDofs();
		_vecStateY = NVec_New(nDOFs);
		_vecStateYdot = NVec_New(nDOFs);

		// Propagate section times if available
		if (_sectionTimes.size() > 0)
		{
			bool* const secCont = new bool[_sectionContinuity.size()];
			std::copy(_sectionContinuity.begin(), _sectionContinuity.end(), secCont);

			// Convert from active to double
			double* const secTimes = new double[_sectionTimes.size()];
			for (unsigned int i = 0; i < _sectionTimes.size(); ++i)
				secTimes[i] = static_cast<double>(_sectionTimes[i]);

			_model->setSectionTimes(secTimes, secCont, _sectionTimes.size() - 1);

			delete[] secTimes;
			delete[] secCont;
		}

		// Initialize with all zeros, correct initial conditions will be set later
		NVec_Const(0.0, _vecStateY);
		NVec_Const(0.0, _vecStateYdot);

		// Create IDAS internal memory
		_idaMemBlock = IDACreate();

		// IDAS Step 4.1: Specify error handler function
		IDASetErrHandlerFn(_idaMemBlock, &idasErrorHandler, this);

		// IDAS Step 5: Initialize the solver
		_model->applyInitialCondition(NVEC_DATA(_vecStateY), NVEC_DATA(_vecStateYdot));

		// Use 0.0 as beginning of simulation time if we haven't set section times yet
		if (_transformedTimes.size() > 0)
			IDAInit(_idaMemBlock, &residualDaeWrapper, _transformedTimes[0], _vecStateY, _vecStateYdot);
		else
			IDAInit(_idaMemBlock, &residualDaeWrapper, 0.0, _vecStateY, _vecStateYdot);			

		// IDAS Step 6: Specify integration tolerances (S: scalar; V: array)
		updateMainErrorTolerances();

		// IDAS Step 7.1: Set optional inputs

		// Attach user data structure
		IDASetUserData(_idaMemBlock, this);

		// Set maximum number of steps
		IDASetMaxNumSteps(_idaMemBlock, _maxSteps);

//		// Set maximum step size ///todo make that user choosable
//		IDASetMaxStep(_idaMemBlock, 50.0);

		// Specify the linear solver.
		IDAMem IDA_mem = static_cast<IDAMem>(_idaMemBlock);

		IDA_mem->ida_linit          = nullptr;
		IDA_mem->ida_lsetup         = nullptr;
		IDA_mem->ida_lsolve         = &linearSolveWrapper;
		IDA_mem->ida_lperf          = nullptr;
		IDA_mem->ida_lfree          = nullptr;
		IDA_mem->ida_setupNonNull   = false;
		IDA_mem->ida_lmem           = this;

		// Allocate memory for AD if required
		if (_model->usesAD())
		{
			_vecADres = new active[nDOFs];
			_vecADy = new active[nDOFs];
		}
	}

	void Simulator::updateMainErrorTolerances()
	{
		if (!_idaMemBlock)
			return;

		if (_absTol.size() > 1)
		{
			if (!_model)
				return;

			N_Vector absTolTemp = NVec_New(_model->numDofs());
			const unsigned int pureDofs = _model->numPureDofs();

			// Check whether user has given us full absolute error for all (pure) DOFs
			if (_absTol.size() >= pureDofs)
			{
				// Copy error tolerances for pure data
				std::copy(_absTol.data(), _absTol.data() + pureDofs, NVEC_DATA(absTolTemp));

				// Calculate error tolerances for coupling DOFs and append them
				const std::vector<double> addAbsErrTol = _model->calculateErrorTolsForAdditionalDofs(_absTol.data(), _absTol.size());
				std::copy(addAbsErrTol.data(), addAbsErrTol.data() + addAbsErrTol.size(), NVEC_DATA(absTolTemp) + pureDofs);
			}
			else
			{
				// We've received an expandable error specification
				_model->expandErrorTol(_absTol.data(), _absTol.size(), NVEC_DATA(absTolTemp));
			}

			IDASVtolerances(_idaMemBlock, _relTol, absTolTemp);
			NVec_Destroy(absTolTemp);
		}
		else
			IDASStolerances(_idaMemBlock, _relTol, _absTol[0]);		
	}

	void Simulator::preFwdSensInit(unsigned int nSens)
	{
		// Turn off solution of sensitivity systems (this will be overridden by a call to IDASensInit below)
		// In fact, this has only an effect, if at first a computation with sensitivities is performed and then
		// (without clearing and reallocating internal memory by cs_free/cs_malloc) another computation without
		// sensitivities is started.
		IDASensToggleOff(_idaMemBlock);

		if (_vecFwdYs)
		{
			NVec_DestroyArray(_vecFwdYs, nSens);
			NVec_DestroyArray(_vecFwdYsDot, nSens);
			_vecFwdYs = nullptr;
			_vecFwdYsDot = nullptr;
		}

		// Allocate sensitivity state vectors
		if (nSens > 0)
		{
			_vecFwdYs     = NVec_CloneArray(nSens, _vecStateY);
			_vecFwdYsDot  = NVec_CloneArray(nSens, _vecStateYdot);

			// Allocate memory for AD if not already done
			if (!_vecADres)
				_vecADres = new active[_model->numDofs()];
		}
	}

	void Simulator::postFwdSensInit(unsigned int nSens)
	{
		// Initialize IDA sensitivity computation
		IDASensInit(_idaMemBlock, nSens, IDA_STAGGERED, &cadet::residualSensWrapper, _vecFwdYs, _vecFwdYsDot);

		// Set sensitivity integration tolerances
		IDASensSStolerances(_idaMemBlock, _relTolS, _absTolS.data());

		// Activate sensitivity error control
		IDASetSensErrCon(_idaMemBlock, true);
	}

	void Simulator::initializeFwdSensitivities()
	{
		const unsigned int nSens = _sensitiveParams.slices();
		preFwdSensInit(nSens);

		if (nSens == 0)
			return;

		for (unsigned int dir = 0; dir < nSens; ++dir)
		{
			// Initialize sensitivity vectors with 0.0
			NVec_Const(0.0, _vecFwdYs[dir]);
			NVec_Const(0.0, _vecFwdYsDot[dir]);
		}

		// Compute consistent initial conditions for sensitivity systems later
		_skipConsistencySensitivity = false;

		postFwdSensInit(nSens);
	}

	void Simulator::initializeFwdSensitivities(double const * const* const initSens, double const * const* const initSensDot)
	{
		const unsigned int nSens = _sensitiveParams.slices();
		preFwdSensInit(nSens);

		if (nSens == 0)
			return;

		if (!initSens)
			throw std::invalid_argument("Pointer to initial sensitivities is NULL");

		if (!initSensDot)
			throw std::invalid_argument("Pointer to initial sensitivity time derivatives is NULL");

		for (unsigned int dir = 0; dir < nSens; ++dir)
		{
			double* const yS = NVEC_DATA(_vecFwdYs[dir]);
			double const* const src = initSens[dir];
			double* const ySdot = NVEC_DATA(_vecFwdYsDot[dir]);
			double const* const srcDot = initSensDot[dir];

			// Initialize sensitivity vectors with given data
			std::copy(src, src + NVEC_LENGTH(_vecFwdYs[dir]), yS);
			std::copy(srcDot, srcDot + NVEC_LENGTH(_vecFwdYsDot[dir]), ySdot);
		}

		// Assume consistent values were given
		_skipConsistencySensitivity = true;

		postFwdSensInit(nSens);
	}

	void Simulator::setInitialConditionFwdSensitivities(double const * const* const initSens, double const * const* const initSensDot)
	{
		const unsigned int nSens = _sensitiveParams.slices();
		if (nSens == 0)
			return;

		if (!initSens && !initSensDot)
		{
			for (unsigned int dir = 0; dir < nSens; ++dir)
			{
				// Initialize sensitivity vectors with 0.0
				NVec_Const(0.0, _vecFwdYs[dir]);
				NVec_Const(0.0, _vecFwdYsDot[dir]);
			}

			// Compute consistent initial conditions for sensitivity systems later
			_skipConsistencySensitivity = false;
		}
		else
		{
			for (unsigned int dir = 0; dir < nSens; ++dir)
			{
				double* const yS = NVEC_DATA(_vecFwdYs[dir]);
				double const* const src = initSens[dir];
				double* const ySdot = NVEC_DATA(_vecFwdYsDot[dir]);
				double const* const srcDot = initSensDot[dir];

				// Initialize sensitivity vectors with given data
				std::copy(src, src + NVEC_LENGTH(_vecFwdYs[dir]), yS);
				std::copy(srcDot, srcDot + NVEC_LENGTH(_vecFwdYsDot[dir]), ySdot);
			}

			// Assume consistent values were given
			_skipConsistencySensitivity = true;
		}
	}

	void Simulator::setInitialCondition(IParameterProvider& paramProvider)
	{
		_model->applyInitialCondition(paramProvider, NVEC_DATA(_vecStateY), NVEC_DATA(_vecStateYdot));
		IDAReInit(_idaMemBlock, _transformedTimes[0], _vecStateY, _vecStateYdot);

		// Better check for consistency
		_skipConsistencyStateY = false;
	}

	void Simulator::setInitialCondition(double const* const initState)
	{
		// Copy initial state
		double* const y = NVEC_DATA(_vecStateY);
		std::copy(initState, initState + NVEC_LENGTH(_vecStateY), y);

		IDAReInit(_idaMemBlock, _transformedTimes[0], _vecStateY, _vecStateYdot);

		// We need to compute matching yDot for consistency
		_skipConsistencyStateY = false;
	}

	void Simulator::setInitialCondition(double const* const initState, double const* const initStateDot)
	{
		// Copy initial state
		double* const y = NVEC_DATA(_vecStateY);
		std::copy(initState, initState + NVEC_LENGTH(_vecStateY), y);

		// Copy initial time derivative state
		double* const yDot = NVEC_DATA(_vecStateYdot);
		std::copy(initStateDot, initStateDot + NVEC_LENGTH(_vecStateY), yDot);

		IDAReInit(_idaMemBlock, _transformedTimes[0], _vecStateY, _vecStateYdot);

		// Assume the initial state is consistent
		_skipConsistencyStateY = true;
	}

	void Simulator::skipConsistentInitialization()
	{
		_skipConsistencyStateY = true;
		_skipConsistencySensitivity = true;
	}

	void Simulator::setConsistentInitialization(ConsistentInitialization ci)
	{
		_consistentInitMode = ci;
	}

	void Simulator::setConsistentInitializationSens(ConsistentInitialization ci)
	{
		_consistentInitModeSens = ci;
	}

	std::unordered_map<ParameterId, double> Simulator::getAllParameterValues() const
	{
		std::unordered_map<ParameterId, double> data;
		if (_model)
			data = _model->getAllParameterValues();

		// Add section times
		const StringHash secTimesName = hashString("SECTION_TIMES");
		for (unsigned int i = 0; i < _sectionTimes.size(); ++i)
			data[makeParamId(secTimesName, UnitOpIndep, CompIndep, BoundPhaseIndep, ReactionIndep, i)] = static_cast<double>(_sectionTimes[i]);

		return data;
	}

	bool Simulator::hasParameter(const ParameterId& pId) const
	{
		if (isSectionTimeParameter(pId, _sectionTimes.size()))
			return true;

		if (_model)
			return _model->hasParameter(pId);

		return false;
	}

	void Simulator::setSensitiveParameter(const ParameterId& id)
	{
		setSensitiveParameter(id, 1.0e-5);
	}

	void Simulator::setSensitiveParameter(const ParameterId& id, double absTolS)
	{
		setSensitiveParameter(&id, nullptr, 1, absTolS);
	}

	void Simulator::setSensitiveParameter(ParameterId const* ids, double const* diffFactors, unsigned int numParams, double absTolS)
	{
		// Set AD directions
		const unsigned int adDir = numSensitivityAdDirections();
		for (unsigned int i = 0; i < numParams; ++i)
		{
			double localDiffFactor = 1.0;
			if (diffFactors)
				localDiffFactor = diffFactors[i];

			bool paramFound = setSectionTimesSensitive(ids[i], adDir, localDiffFactor);
			paramFound = _model->setSensitiveParameter(ids[i], adDir, localDiffFactor) || paramFound;

			if (!paramFound)
			{
				LOG(Warning) << "Warning: Unkown parameter " << ids[i] << " in parameter join was ignored";
			}
		}

		_sensitiveParams.pushBackSlice(ids, numParams);
		_absTolS.push_back(absTolS);
		
		if (diffFactors)
			_sensitiveParamsFactor.insert(_sensitiveParamsFactor.end(), diffFactors, diffFactors + numParams);
		else
			_sensitiveParamsFactor.insert(_sensitiveParamsFactor.end(), numParams, 1.0);
	}

	void Simulator::setSensitiveParameter(ParameterId const* ids, unsigned int numParams, double absTolS)
	{
		setSensitiveParameter(ids, nullptr, numParams, absTolS);
	}

	void Simulator::resetSensParams()
	{
		using std::to_string;

		unsigned int globalIdx = 0;
		for (unsigned int i = 0; i < _sensitiveParams.slices(); ++i)
		{
			ParameterId const* const ids = _sensitiveParams[i];

			for (unsigned int j = 0; j < _sensitiveParams.sliceSize(i); ++j, ++globalIdx)
			{
				bool paramFound = setSectionTimesSensitive(ids[j], i, _sensitiveParamsFactor[globalIdx]);
				paramFound = _model->setSensitiveParameter(ids[j], i, _sensitiveParamsFactor[globalIdx]) || paramFound;

				if (!paramFound)
					throw InvalidParameterException("Sensitive parameter " + to_string(ids[j]) + " disappeared");
			}
		}
	}

	bool Simulator::setSectionTimesSensitive(const ParameterId& id, unsigned int adDirection, double adValue)
	{
		if (isSectionTimeParameter(id, _sectionTimes.size()))
		{
			// Correct adValue
			LOG(Debug) << "Found parameter " << id << " in SECTION_TIMES: Dir " << adDirection << " is set to " << adValue << " [" << static_cast<double>(_sectionTimes[id.section]) << " < " << _sectionTimes.size() << "]";
			_sectionTimes[id.section].setADValue(adDirection, adValue);

			return true;
		}
		return false;
	}

	void Simulator::clearSensParams()
	{
		_sensitiveParams.clear();
		_sensitiveParamsFactor.clear();
		_absTolS.clear();

		_model->clearSensParams();
		for (unsigned int i = 0; i < _sectionTimes.size(); ++i)
			_sectionTimes[i].setADValue(0.0);

		initializeFwdSensitivities();
	}	

	unsigned int Simulator::numSensParams() const CADET_NOEXCEPT
	{
		return _sensitiveParams.slices();
	}

	void Simulator::setSensitiveParameterValue(const ParameterId& id, double value)
	{
		if (isSectionTimeParameter(id, _sectionTimes.size()))
		{
			_sectionTimes[id.section].setValue(value);

			// Update the time transformation
			calculateTimeTransformation(false);

			// Do not exit here, since model can also contain instances of SECTION_TIMES
		}

		if (!_model)
			return;

		util::SlicedVector<ParameterId>::size_type temp;
		util::SlicedVector<ParameterId>::size_type linearIndex;
		util::SlicedVector<ParameterId>::size_type idxSlice = _sensitiveParams.findElementAndSlice(id, temp, linearIndex);
		if (idxSlice < _sensitiveParams.slices())
		{
			ParameterId const* const fusedIds = _sensitiveParams[idxSlice];
			const util::SlicedVector<ParameterId>::size_type sliceOffset = _sensitiveParams.sliceOffset(idxSlice);

			// Noramlize parameter value
			value /= _sensitiveParamsFactor[linearIndex];

			// Take care of linear factors
			// @todo add constant offset to linear combination
			for (unsigned int i = 0; i < _sensitiveParams.sliceSize(idxSlice); ++i)
				_model->setSensitiveParameterValue(fusedIds[i], _sensitiveParamsFactor[sliceOffset + i] * value);
		}
	}

	void Simulator::setSensitiveParameterValue(unsigned int idx, double value)
	{
		if (idx >= _sensitiveParams.slices())
			return;

		ParameterId const* const paramIds = _sensitiveParams[idx];
		const util::SlicedVector<ParameterId>::size_type sliceOffset = _sensitiveParams.sliceOffset(idx);

		for (unsigned int i = 0; i < _sensitiveParams.sliceSize(idx); ++i)
		{
			const ParameterId& id = paramIds[i];
			if (isSectionTimeParameter(id, _sectionTimes.size()))
			{
				_sectionTimes[id.section].setValue(value);

				// Update the time transformation
				calculateTimeTransformation(false);

				// Do not exit here, since model can also contain instances of SECTION_TIMES
			}

			if (!_model)
				continue;

			// Take care of linear factors
			// @todo add constant offset to linear combination
			_model->setSensitiveParameterValue(id, _sensitiveParamsFactor[sliceOffset + i] * value);
		}
	}

	void Simulator::setSensitiveParameterFactors(unsigned int idx, double const* factors)
	{
		if (idx >= _sensitiveParams.slices())
			return;

		ParameterId const* const paramIds = _sensitiveParams[idx];
		const util::SlicedVector<ParameterId>::size_type sliceOffset = _sensitiveParams.sliceOffset(idx);

		for (unsigned int i = 0; i < _sensitiveParams.sliceSize(idx); ++i)
		{
			// Update the linear factor
			_sensitiveParamsFactor[sliceOffset + i] = factors[i];

			const ParameterId& id = paramIds[i];
			setSectionTimesSensitive(id, i, factors[i]);

			if (_model)
				_model->setSensitiveParameter(id, i, factors[i]);
		}
	}

	void Simulator::setParameterValue(const ParameterId& id, double value)
	{
		if (isSectionTimeParameter(id, _sectionTimes.size()))
		{
			_sectionTimes[id.section].setValue(value);

			// Update the time transformation
			calculateTimeTransformation(false);
			
			// Do not exit here, since model can also contain instances of SECTION_TIMES
		}

		if (_model)
			_model->setParameter(id, value);
	}

	void Simulator::setSolutionTimes(const std::vector<double>& solutionTimes)
	{
		_solutionTimes = solutionTimes;
		_solutionTimesOriginal = solutionTimes;

		// Transform user solution times
		transformSolutionTimes();
	}

	const std::vector<double>& Simulator::getSolutionTimes() const
	{
		return _solutionTimes;
	}

	void Simulator::setSectionTimes(const std::vector<double>& sectionTimes)
	{
		setSectionTimes(sectionTimes, std::vector<bool>(sectionTimes.size() - 1, false));
	}

	void Simulator::setSectionTimes(const std::vector<double>& sectionTimes, const std::vector<bool>& sectionContinuity)
	{
		// Ensure that at least one section is defined
		if (sectionTimes.size() < 2)
			throw std::invalid_argument("At least one section has to be specified!");

		// Ensure that all section start times are smaller than their end times
		for (unsigned int i = 0; i < sectionTimes.size() - 1; ++i)
			if (sectionTimes[i] > sectionTimes[i + 1])
				throw InvalidParameterException("The end time of each section must be greater than its start time (failed for section " + std::to_string(i) + ")!");

		// Send original section times to models
		if (_model)
		{
			bool* const secCont = new bool[sectionContinuity.size()];
			std::copy(sectionContinuity.begin(), sectionContinuity.end(), secCont);
			_model->setSectionTimes(sectionTimes.data(), secCont, sectionTimes.size() - 1);
			delete[] secCont;
		}

		// Copy section times into AD (active) data type
		_sectionTimes.clear();
		_sectionTimes.reserve(sectionTimes.size());
		for (unsigned int i = 0; i < sectionTimes.size(); ++i)
			_sectionTimes.push_back(sectionTimes[i]);

		_sectionContinuity = sectionContinuity;

		// Calculate the transformed section times
		calculateTimeTransformation(true);

		// Set AD sensitivities
		unsigned int globalIdx = 0;
		for (unsigned int i = 0; i < _sensitiveParams.slices(); ++i)
		{
			ParameterId const* const ids = _sensitiveParams[i];

			for (unsigned int j = 0; j < _sensitiveParams.sliceSize(i); ++j, ++globalIdx)
			{
				setSectionTimesSensitive(ids[j], i, _sensitiveParamsFactor[globalIdx]);
			}
		}
	}

	void Simulator::transformSolutionTimes()
	{
		for (unsigned int i = 0; i < _solutionTimes.size(); ++i)
			_solutionTimes[i] = toTransformedTime(_solutionTimesOriginal[i], _sectionTimes, _transformedTimes);
	}

	void Simulator::calculateTimeTransformation(bool resetSens)
	{
		// Transform section times
		_transformedTimes.clear();
		_transformedTimes.reserve(_sectionTimes.size());
		for (unsigned int i = 0; i < _sectionTimes.size(); ++i)
			_transformedTimes.push_back(static_cast<double>(_sectionTimes[i]));
//			_transformedTimes.push_back(i);

		// Transform user solution times
		transformSolutionTimes();

		// Set AD directions for section times again
		if (resetSens)
			resetSensParams();

		// Update IDAS
		IDAReInit(_idaMemBlock, _transformedTimes[0], _vecStateY, _vecStateYdot);		
	}

	void Simulator::setSolutionRecorder(ISolutionRecorder* recorder)
	{
		_solRecorder = recorder;
		if (_solRecorder)
		{
			_solRecorder->prepare(NVEC_LENGTH(_vecStateY), _sensitiveParams.slices(), _solutionTimes.size());
			_model->reportSolutionStructure(*_solRecorder);
		}
	}

	const active Simulator::timeFactor(unsigned int curSec) const
	{
//		return (_transformedTimes[curSec + 1] - _transformedTimes[curSec]) / static_cast<double>(_sectionTimes[curSec + 1] - _sectionTimes[curSec]);
		return (_transformedTimes[curSec + 1] - _transformedTimes[curSec]) / (_sectionTimes[curSec + 1] - _sectionTimes[curSec]);
	}

	template <typename time_t>
	double Simulator::toTransformedTime(double t, const std::vector<time_t>& oldTimes, const std::vector<double>& newTimes) const
	{
		if (t == static_cast<double>(oldTimes[0]))
			return newTimes[0];
		
		for (unsigned int i = 1; i < newTimes.size(); ++i)
		{
			if (t <= static_cast<double>(oldTimes[i]))
				return newTimes[i-1] + (t - static_cast<double>(oldTimes[i-1])) / (static_cast<double>(oldTimes[i]) - static_cast<double>(oldTimes[i-1])) * (newTimes[i] - newTimes[i-1]);
		}

		return -1.0;
	}

	active Simulator::toRealTime(double t, unsigned int curSec) const
	{
		for (unsigned int i = curSec; i < _transformedTimes.size()-1; ++i)
		{
			if ((t >= _transformedTimes[i]) && (t <= _transformedTimes[i + 1]))
			{
				return _sectionTimes[i] + (t - _transformedTimes[i]) / (_transformedTimes[i + 1] - _transformedTimes[i]) * (_sectionTimes[i + 1] - _sectionTimes[i]);
			}
		}
		return active(-1.0);
	}

	void Simulator::integrate()
	{
		// In this function the model is integrated by IDAS from the SUNDIALS package.
		// The authors of IDAS recommend to restart the time integrator when a discontinuity
		// is encountered (see https://computation.llnl.gov/casc/sundials/support/notes.html#disc).
		// The sectionTime (together with the sectionContinuity) array indicates such
		// discontinuitites and the solver is restarted accordingly. This also requires
		// the computation of consistent initial values for each restart.

		_timerIntegration.start();

		// Set number of AD directions
		// @todo This is problematic if multiple Simulators are run concurrently!
#if defined(ACTIVE_ADOLC) || defined(ACTIVE_SFAD) || defined(ACTIVE_SETFAD)
		LOG(Debug) << "Setting AD directions from " << ad::getDirections() << " to " << numSensitivityAdDirections() + _model->requiredADdirs();
		ad::setDirections(numSensitivityAdDirections() + _model->requiredADdirs());
#endif
		// Setup AD vectors by model
		// @todo Check if this is necessary (dirty flag)
		_model->prepareADvectors(_vecADres, _vecADy, numSensitivityAdDirections());

		std::vector<double>::const_iterator it;
		double tOut = 0.0;

		const bool writeAtUserTimes = _solutionTimes.size() > 0;
		const bool wantSensitivities = _sensitiveParams.slices() > 0;

		if (_solRecorder)
		{
			_solRecorder->notifyIntegrationStart(NVEC_LENGTH(_vecStateY), _sensitiveParams.slices(), _solutionTimes.size());
			_model->reportSolutionStructure(*_solRecorder);			
		}

		// Decide whether to use user specified solution output times (IDA_NORMAL)
		// or internal integrator steps (IDA_ONE_STEP)
		int idaTask = IDA_ONE_STEP;
		if (writeAtUserTimes)
		{
			idaTask = IDA_NORMAL;
		}

		LOG(Debug) << "Integration span: [" << _transformedTimes[0] << ", " << _transformedTimes.back() 
			<< "] transformed, [" << static_cast<double>(_sectionTimes[0]) << ", " << static_cast<double>(_sectionTimes.back()) << "] sections";
		
		if (writeAtUserTimes)
		{
			LOG(Debug) << "Solution time span: [" << _solutionTimes[0] << ", " << _solutionTimes.back() << "]";
		}

		double transformedT = _transformedTimes[0];
		_curSec = 0;
		const double tEnd = writeAtUserTimes ? _solutionTimes.back() : _transformedTimes.back();
		while (transformedT < tEnd)
		{
			// Get smallest index with t_i >= transformedT (t_i being a _transformedTimes element)
			// This will return i if transformedT == _transformedTimes[i], which effectively advances
			// the index if required
			_curSec = getNextSection(transformedT, _curSec);
			const double startTime = _transformedTimes[_curSec];

			// Determine continuous time slice
			unsigned int skip = 1; // Always finish the current section
			for (size_t i = _curSec; i < _transformedTimes.size() - 2; ++i)
			{
				if (!_sectionContinuity[i])
					break;

				// This is a continuous section transition, so we don't need to
				// restart the integrator and just integrate for a longer time
				++skip;
			}

			const double endTime = writeAtUserTimes ? std::min(_transformedTimes[_curSec + skip], tEnd) : _transformedTimes[_curSec + skip];
			transformedT = startTime;

			LOG(Debug) << " ###### SECTION " << _curSec << " from " << startTime << " to " << endTime;

			// IDAS Step 7.3: Set the initial step size
			const double stepSize = _initStepSize.size() > 1 ? _initStepSize[_curSec] : _initStepSize[0];
			IDASetInitStep(_idaMemBlock, stepSize);

			// IDAS Step 7.4: Set the stop time
			IDASetStopTime(_idaMemBlock, endTime);

			// Update Jacobian
			active realT = toRealTime(transformedT, _curSec);
			_model->notifyDiscontinuousSectionTransition(static_cast<double>(realT), _curSec);

			// Get time factor
			const active curTimeFactor = timeFactor(_curSec);

			// Compute consistent initial values
			LOG(Debug) << "---====--- CONSISTENCY ---====--- ";
			const double consPrev = _model->residualNorm(static_cast<double>(realT), _curSec, static_cast<double>(curTimeFactor), NVEC_DATA(_vecStateY), NVEC_DATA(_vecStateYdot));
			LOG(Debug) << " ==========> Consistency error prev: " << consPrev;

			if (!_skipConsistencyStateY && (_consistentInitMode != ConsistentInitialization::None))
			{
				if ((_consistentInitMode == ConsistentInitialization::Full) || ((_curSec == 0) && (_consistentInitMode == ConsistentInitialization::FullFirstOnly)))
				{
					_model->consistentInitialConditions(static_cast<double>(realT), _curSec, static_cast<double>(curTimeFactor), NVEC_DATA(_vecStateY), NVEC_DATA(_vecStateYdot), 
						_vecADres, _vecADy, numSensitivityAdDirections(), _algTol);

					const double consPost = _model->residualNorm(static_cast<double>(realT), _curSec, static_cast<double>(curTimeFactor), NVEC_DATA(_vecStateY), NVEC_DATA(_vecStateYdot));
					LOG(Debug) << " ==========> Consistency error post Full: " << consPost;
				}
				else if ((_consistentInitMode == ConsistentInitialization::Lean) || ((_curSec == 0) && (_consistentInitMode == ConsistentInitialization::LeanFirstOnly)))
				{
					_model->leanConsistentInitialConditions(static_cast<double>(realT), _curSec, static_cast<double>(curTimeFactor), NVEC_DATA(_vecStateY), NVEC_DATA(_vecStateYdot), 
						_vecADres, _vecADy, numSensitivityAdDirections(), _algTol);

					const double consPost = _model->residualNorm(static_cast<double>(realT), _curSec, static_cast<double>(curTimeFactor), NVEC_DATA(_vecStateY), NVEC_DATA(_vecStateYdot));
					LOG(Debug) << " ==========> Consistency error post Lean: " << consPost;
				}
				else
				{
					LOG(Debug) << " ==========> Consistent initialization NOT performed (mode " << to_string(_consistentInitMode) << ")";
				}

				LOG(Debug) << "y = " << _vecStateY << ";";
				LOG(Debug) << "yDot = " << _vecStateYdot << ";";
			}
			_skipConsistencyStateY = false;

			if ((_sensitiveParams.slices() > 0) && !_skipConsistencySensitivity && (_consistentInitModeSens != ConsistentInitialization::None))
			{
				const std::vector<const double*> sensYdbg = convertNVectorToStdVectorPtrs<const double*>(_vecFwdYs, _sensitiveParams.slices());
				const std::vector<const double*> sensYdotDbg = convertNVectorToStdVectorPtrs<const double*>(_vecFwdYsDot, _sensitiveParams.slices());

				std::vector<double> norms(sensYdbg.size(), 0.0);
				std::vector<double> temp(_model->numDofs(), 0.0);
				_model->residualSensFwdNorm(sensYdbg.size(), realT, _curSec, curTimeFactor, NVEC_DATA(_vecStateY), NVEC_DATA(_vecStateYdot),
					sensYdbg, sensYdotDbg, norms.data(), _vecADres, temp.data());

				LOG(Debug) << " ==========> Sens consistency error prev: " << norms;

				if ((_consistentInitModeSens == ConsistentInitialization::Full) || ((_curSec == 0) && (_consistentInitModeSens == ConsistentInitialization::FullFirstOnly)))
				{
					// Compute consistent initial conditions for sensitivity subsystems
					std::vector<double*> sensY = convertNVectorToStdVectorPtrs<double*>(_vecFwdYs, _sensitiveParams.slices());
					std::vector<double*> sensYdot = convertNVectorToStdVectorPtrs<double*>(_vecFwdYsDot, _sensitiveParams.slices());
					_model->consistentIntialSensitivity(realT, _curSec, curTimeFactor, NVEC_DATA(_vecStateY), NVEC_DATA(_vecStateYdot), sensY, sensYdot, _vecADres, _vecADy);

					_model->residualSensFwdNorm(sensYdbg.size(), realT, _curSec, curTimeFactor, NVEC_DATA(_vecStateY), NVEC_DATA(_vecStateYdot),
						sensYdbg, sensYdotDbg, norms.data(), _vecADres, temp.data());

					LOG(Debug) << " ==========> Sens consistency error post Full: " << norms;
				}
				else if ((_consistentInitModeSens == ConsistentInitialization::Lean) || ((_curSec == 0) && (_consistentInitModeSens == ConsistentInitialization::LeanFirstOnly)))
				{
					// Compute consistent initial conditions for sensitivity subsystems
					std::vector<double*> sensY = convertNVectorToStdVectorPtrs<double*>(_vecFwdYs, _sensitiveParams.slices());
					std::vector<double*> sensYdot = convertNVectorToStdVectorPtrs<double*>(_vecFwdYsDot, _sensitiveParams.slices());
					_model->leanConsistentIntialSensitivity(realT, _curSec, curTimeFactor, NVEC_DATA(_vecStateY), NVEC_DATA(_vecStateYdot), sensY, sensYdot, _vecADres, _vecADy);

					_model->residualSensFwdNorm(sensYdbg.size(), realT, _curSec, curTimeFactor, NVEC_DATA(_vecStateY), NVEC_DATA(_vecStateYdot),
						sensYdbg, sensYdotDbg, norms.data(), _vecADres, temp.data());

					LOG(Debug) << " ==========> Sens consistency error post Lean: " << norms;
				}
				else
				{
					LOG(Debug) << " ==========> Sens consistent initialization NOT performed (mode " << to_string(_consistentInitModeSens) << ")";
				}

				for (unsigned int j = 0; j < norms.size(); ++j)
				{
					LOG(Debug) << "sensY[" << j << "] = " << _vecFwdYs[j] << ";";
					LOG(Debug) << "sensYdot[" << j << "] = " << _vecFwdYsDot[j] << ";";
				}
			}
			_skipConsistencySensitivity = false;

			// IDAS Step 5.2: Re-initialization of the solver
			IDAReInit(_idaMemBlock, startTime, _vecStateY, _vecStateYdot);
			if (numSensParams() > 0)
				IDASensReInit(_idaMemBlock, IDA_STAGGERED, _vecFwdYs, _vecFwdYsDot);

			// Inititalize the IDA solver flag
			int solverFlag = IDA_SUCCESS;

			if (writeAtUserTimes)
			{
				// Write initial conditions only if desired by user
				if (_curSec == 0 && _solutionTimes.front() == transformedT)
					writeSolution(static_cast<double>(realT));

				// Initialize iterator and forward it to the first solution time that lies inside the current section
				it = _solutionTimes.begin();
				while ((*it) <= startTime) ++it;
			}
			else
			{
				// Always write initial conditions if solutions are written at integration times
				if (_curSec == 0) writeSolution(static_cast<double>(realT));

				// Here tOut - only during the first call to IDASolve - specifies the direction
				// and rough scale of the independent variable, see IDAS Guide p.33
				tOut = endTime;
			}

			// Main loop which integrates the system until reaching the end time of the current section
			// or until an error occures
			while ((solverFlag == IDA_SUCCESS) || (solverFlag == IDA_ROOT_RETURN))
			{
				// Update tOut if we write solutions at user specified times
				if (writeAtUserTimes)
				{
					// Check if user specified times are sufficiently long.
					// otherwise integrate till IDA_TSTOP_RETURN
					if (it == _solutionTimes.end())
						break;
					else
						tOut = *it;
				}

				// IDA Step 11: Advance solution in time
				solverFlag = IDASolve(_idaMemBlock, tOut, &transformedT, _vecStateY, _vecStateYdot, idaTask);
				LOG(Debug) << "Solve from " << transformedT << " to " << tOut << " => " << (solverFlag == IDA_SUCCESS ? "IDA_SUCCESS" : "") << (solverFlag == IDA_TSTOP_RETURN ? "IDA_TSTOP_RETURN" : "");
				realT = toRealTime(transformedT, _curSec);

				switch (solverFlag)
				{
				case IDA_SUCCESS:
					// tOut was reached

					// Extract sensitivity information from IDA (required for consistent initialization
					// and output of sensitivities)
					if (wantSensitivities)
					{
						IDAGetSens(_idaMemBlock, &transformedT, _vecFwdYs);
						IDAGetSensDky(_idaMemBlock, transformedT, 1, _vecFwdYsDot);
					}
					writeSolution(static_cast<double>(realT));
					++it;
					break;
				case IDA_ROOT_RETURN:
					// A root was found
					// Eventually call some routine
					break;
				case IDA_TSTOP_RETURN:
					// Extract sensitivity information from IDA (required for consistent initialization
					// and output of sensitivities)
					if (wantSensitivities)
					{
						IDAGetSens(_idaMemBlock, &transformedT, _vecFwdYs);
						IDAGetSensDky(_idaMemBlock, transformedT, 1, _vecFwdYsDot);
					}

					// Section end time was reached (in previous step)
					if (!writeAtUserTimes && (endTime == _transformedTimes.back()))
					{
						// Write a solution for the ultimate endTime in the last section,
						// when we write at integration times.
						writeSolution(static_cast<double>(realT));
					}
					break;
				default:
					_lastIntTime = _timerIntegration.stop();

					// An error occured
					LOG(Error) << "IDASolve returned " << IDAGetReturnFlagName(solverFlag) << " at t = " << static_cast<double>(realT);
					throw IntegrationException("Error in IDASolve!"); //todo might not be necessary
					break;
				} // switch

			} // while

		} // for (_sec ...)

		_lastIntTime = _timerIntegration.stop();
	}

	double const* Simulator::getLastSolution(unsigned int& len) const
	{
		len = NVEC_LENGTH(_vecStateY);
		return NVEC_DATA(_vecStateY);
	}

	double const* Simulator::getLastSolutionDerivative(unsigned int& len) const
	{
		len = NVEC_LENGTH(_vecStateYdot);
		return NVEC_DATA(_vecStateYdot);
	}

	const std::vector<double const*> Simulator::getLastSensitivities(unsigned int& len) const
	{
		return convertNVectorToStdVectorConstPtrs(len, _vecFwdYs, _sensitiveParams.slices());
	}

	const std::vector<double const*> Simulator::getLastSensitivityDerivatives(unsigned int& len) const
	{
		return convertNVectorToStdVectorConstPtrs(len, _vecFwdYsDot, _sensitiveParams.slices());
	}

	void Simulator::configureTimeIntegrator(double relTol, double absTol, double initStepSize, unsigned int maxSteps)
	{
		_absTol.clear();
		_absTol.push_back(absTol);

		_relTol = relTol;
		_maxSteps = maxSteps;
		
		_initStepSize.clear();
		_initStepSize.push_back(initStepSize);
	}

	void Simulator::configureTimeIntegrator(double relTol, double absTol, const std::vector<double>& initStepSizes, unsigned int maxSteps)
	{
		_absTol.clear();
		_absTol.push_back(absTol);

		_relTol = relTol;
		_maxSteps = maxSteps;
		_initStepSize = initStepSizes;
	}

	void Simulator::configure(IParameterProvider& paramProvider)
	{
		paramProvider.pushScope("time_integrator");
		
		_absTol.clear();
		if (paramProvider.isArray("ABSTOL"))
			_absTol = paramProvider.getDoubleArray("ABSTOL");
		else
			_absTol.push_back(paramProvider.getDouble("ABSTOL"));

		_relTol = paramProvider.getDouble("RELTOL");
		_algTol = paramProvider.getDouble("ALGTOL");
		_maxSteps = paramProvider.getInt("MAX_STEPS");

		_initStepSize.clear();
		if (paramProvider.isArray("INIT_STEP_SIZE"))
			_initStepSize = paramProvider.getDoubleArray("INIT_STEP_SIZE");
		else
			_initStepSize.push_back(paramProvider.getDouble("INIT_STEP_SIZE"));
		
		if (paramProvider.exists("RELTOL_SENS"))
			_relTolS = paramProvider.getDouble("RELTOL_SENS");
		else
			_relTolS = _relTol;

		paramProvider.popScope();

#ifdef _OPENMP
		if (paramProvider.exists("NTHREADS"))
			setNumThreads(paramProvider.getInt("NTHREADS"));
#endif

		_solutionTimes.clear();
		_solutionTimesOriginal.clear();
		if (paramProvider.exists("USER_SOLUTION_TIMES"))
		{
			_solutionTimes = paramProvider.getDoubleArray("USER_SOLUTION_TIMES");
			_solutionTimesOriginal = _solutionTimes;
		}

		if (paramProvider.exists("CONSISTENT_INIT_MODE"))
			_consistentInitMode = toConsistentInitialization(paramProvider.getInt("CONSISTENT_INIT_MODE"));

		if (paramProvider.exists("CONSISTENT_INIT_MODE_SENS"))
			_consistentInitModeSens = toConsistentInitialization(paramProvider.getInt("CONSISTENT_INIT_MODE_SENS"));

		// @todo: Read more configuration values
	}

	void Simulator::reconfigure(IParameterProvider& paramProvider)
	{
		configure(paramProvider);
	}

	void Simulator::setSensitivityErrorTolerance(double relTol, double const* absTol)
	{
		_relTolS = relTol;
		_absTolS.clear();

		const unsigned int nSens = _sensitiveParams.slices();
		if (_idaMemBlock && (nSens > 0) && _vecFwdYs && absTol)
		{
			_absTolS.insert(_absTolS.end(), absTol, absTol + nSens);

			// Set sensitivity integration tolerances
			IDASensSStolerances(_idaMemBlock, _relTolS, _absTolS.data());
		}
	}

	void Simulator::setRelativeErrorToleranceSens(double relTol)
	{
		_relTolS = relTol;
		if (_idaMemBlock && (_sensitiveParams.slices() > 0) && _vecFwdYs)
		{
			// Set sensitivity integration tolerances
			IDASensSStolerances(_idaMemBlock, _relTolS, _absTolS.data());
		}
	}

	void Simulator::setRelativeErrorTolerance(double relTol)
	{
		_relTol = relTol;
		updateMainErrorTolerances();
	}

	void Simulator::setAbsoluteErrorTolerance(double absTol)
	{
		_absTol.clear();
		_absTol.push_back(absTol);
		updateMainErrorTolerances();
	}

	void Simulator::setAbsoluteErrorTolerance(const std::vector<double>& absTol)
	{
		_absTol = absTol;
		updateMainErrorTolerances();
	}

	void Simulator::setInitialStepSize(double stepSize)
	{
		_initStepSize.clear();
		_initStepSize.push_back(stepSize);
	}

	void Simulator::setInitialStepSize(const std::vector<double>& stepSize)
	{
		_initStepSize = stepSize;
	}

	void Simulator::setMaximumSteps(unsigned int maxSteps)
	{
		_maxSteps = maxSteps;
		if (_idaMemBlock)
			IDASetMaxNumSteps(_idaMemBlock, _maxSteps);
	}


	bool Simulator::reconfigureModel(IParameterProvider& paramProvider)
	{
		if (!_model)
			return false;

		// Reconfigure the model
		const bool success = _model->reconfigure(paramProvider);

		// Set all AD directions for parameter sensitivities again
		resetSensParams();

		return success;
	}

	bool Simulator::reconfigureModel(IParameterProvider& paramProvider, unsigned int unitOpIdx)
	{
		if (!_model)
			return false;

		// Reconfigure the model
		const bool success = _model->reconfigureModel(paramProvider, unitOpIdx);

		// Set all AD directions for parameter sensitivities again
		resetSensParams();

		return success;
	}

	void Simulator::writeSolution(double t)
	{
		if (!_solRecorder)
			return;

		_solRecorder->beginTimestep(t);
		
		_solRecorder->beginSolution();
		_model->reportSolution(*_solRecorder, NVEC_DATA(_vecStateY));
		_solRecorder->endSolution();

		_solRecorder->beginSolutionDerivative();
		_model->reportSolution(*_solRecorder, NVEC_DATA(_vecStateYdot));
		_solRecorder->endSolutionDerivative();

		for (unsigned int i = 0; i < _sensitiveParams.slices(); ++i)
		{
			_solRecorder->beginSensitivity(*_sensitiveParams[i], i);
			_model->reportSolution(*_solRecorder, NVEC_DATA(_vecFwdYs[i]));
			_solRecorder->endSensitivity(*_sensitiveParams[i], i);

			_solRecorder->beginSensitivityDerivative(*_sensitiveParams[i], i);
			_model->reportSolution(*_solRecorder, NVEC_DATA(_vecFwdYsDot[i]));
			_solRecorder->endSensitivityDerivative(*_sensitiveParams[i], i);
		}

		_solRecorder->endTimestep();
	}

	unsigned int Simulator::getNextSection(double t, unsigned int startIdx) const
	{
		if (t < _transformedTimes[startIdx])
			return -1;

		for (unsigned int i = startIdx; i < _transformedTimes.size() - 1; ++i)
		{
			if (_transformedTimes[i] >= t)
				return i;
		}

		return -1;
	}

	unsigned int Simulator::getCurrentSection(double t) const
	{
		//TODO: Use binary search
		
		for (unsigned int i = _curSec; i < _transformedTimes.size() - 1; ++i)
		{
			if ((t >= _transformedTimes[i]) && (t <= _transformedTimes[i+1]))
				return i;
		}

		return -1;
	}

	IModelSystem* const Simulator::model() CADET_NOEXCEPT
	{
		return _model;
	}

	IModelSystem const* const Simulator::model() const CADET_NOEXCEPT
	{
		return _model;
	}

	unsigned int Simulator::numDofs() const CADET_NOEXCEPT
	{
		if (_model)
			return _model->numDofs();
		return 0;
	}

	void Simulator::setNumThreads(unsigned int nThreads) const
	{
#ifdef _OPENMP
		if (nThreads == 0)
			nThreads = omp_get_max_threads();
		omp_set_num_threads(nThreads);

		#ifdef CADET_SUNDIALS_OPENMP
			if (_vecStateY)
				NVec_SetThreads(_vecStateY, nThreads);
			if (_vecStateYdot)
				NVec_SetThreads(_vecStateYdot, nThreads);

			for (unsigned int i = 0; i < _sensitiveParams.slices(); ++i)
			{
				NVec_SetThreads(_vecFwdYs[i], nThreads);
				NVec_SetThreads(_vecFwdYsDot[i], nThreads);
			}
		#endif
#endif
	}

} // namespace cadet
