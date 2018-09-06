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

/**
 * @file 
 * Implements the LinearBinding class. This is the simplest, self-contained, and
 * complete example of an IBindingModel implementation. A slightly more complicated
 * example is given by LangmuirBinding.
 * 
 * The implementation in this file intentionally does not use common code as provided
 * by BindingModelBase or PureBindingModelBase. Thus, it can be used as an example
 * for learning how to develop custom binding models and getting an impression of
 * what happens behind the curtains. For serious binding models, have a look at the
 * LangmuirBinding model, which makes use of common code.
 * 
 * For implementing externally dependent binding models, the storage of the model 
 * parameters is encapsulated from the actual model implementation. The default
 * storage implementation just stores, configures, and registers the standard model
 * parameters and does not provide external dependence. A second implementation
 * additionally takes care of all parameters required for external dependence (i.e.,
 * additional parameters, external data source). The actual parameters (including
 * dependence on time, axial position, and external source) are cached in variables
 * which have the same name as the ones of the default storage. Using a templated
 * base implementation of the binding model, the different storage implementations
 * are plugged into the binding model, which always uses the variable names of the
 * default storage (default storage => original variables, dependent storage => cache).
 * In this way, support for external dependence can be added with minimal code duplication.
 * 
 * Each binding model implementation should provide a function that registers all
 * binding model variants (default, external dependence) in a map. This function
 * is expected to be in the cadet::model::binding namespace and usually found at
 * the very bottom of this file.
 */

#include "model/BindingModel.hpp"
#include "model/ExternalFunctionSupport.hpp"
#include "ParamIdUtil.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "model/Parameters.hpp"
#include "LocalVector.hpp"

#include <vector>
#include <unordered_map>
#include <functional>
#include <algorithm>
#include <iterator>

namespace cadet
{

namespace model
{

/**
 * @brief Handles linear binding model parameters that do not depend on external functions
 */
class LinearParamHandler : public ConstParamHandlerBase
{
public:

	/**
	 * @brief Holds actual parameter data
	 */
	struct ConstParams
	{
		std::vector<active> kA; //!< Adsorption rate
		std::vector<active> kD; //!< Desorption rate
	};

	typedef ConstParams params_t;

	/**
	 * @brief Returns name of the binding model
	 * @return Name of the binding model
	 */
	static const char* identifier() { return "LINEAR"; }

	LinearParamHandler() CADET_NOEXCEPT : _kA(&_localParams.kA), _kD(&_localParams.kD) { }

	/**
	 * @brief Reads parameters and verifies them
	 * @details See IBindingModel::configure() for details.
	 * @param [in] paramProvider IParameterProvider used for reading parameters
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 * @return @c true if the parameters were read and validated successfully, otherwise @c false
	 */
	inline bool configure(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBoundStates)
	{
		_kA.configure("LIN_KA", paramProvider, nComp, nBoundStates);
		_kD.configure("LIN_KD", paramProvider, nComp, nBoundStates);
		return validateConfig(nComp, nBoundStates);
	}

	/**
	 * @brief Registers all local parameters in a map for further use
	 * @param [in,out] parameters Map in which the parameters are stored
	 * @param [in] unitOpIdx Index of the unit operation used for registering the parameters
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 */
	inline void registerParameters(std::unordered_map<ParameterId, active*>& parameters, unsigned int unitOpIdx, unsigned int nComp, unsigned int const* nBoundStates)
	{
		_kA.registerParam("LIN_KA", parameters, unitOpIdx, nComp, nBoundStates);
		_kD.registerParam("LIN_KD", parameters, unitOpIdx, nComp, nBoundStates);
	}

	/**
	 * @brief Reserves space in the storage of the parameters
	 * @param [in] numElem Number of elements (total)
	 * @param [in] numSlices Number of slices / binding site types
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 */
	inline void reserve(unsigned int numElem, unsigned int numSlices, unsigned int nComp, unsigned int const* nBoundStates) \
	{
		_kA.reserve(numElem, numSlices, nComp, nBoundStates);
		_kD.reserve(numElem, numSlices, nComp, nBoundStates);
	}

	/**
	 * @brief Updates the parameter cache in order to take the external profile into account
	 * @param [in] t Current time
	 * @param [in] z Axial coordinate in the column
	 * @param [in] r Radial coordinate in the bead
	 * @param [in] secIdx Index of the current section
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 * @param [in,out] workSpace Memory buffer for updated data
	 * @return Externally dependent parameter values
	 */
	inline const params_t& update(double t, double z, double r, unsigned int secIdx, unsigned int nComp, unsigned int const* nBoundStates, void* workSpace) const
	{
		return _localParams;
	}

	/**
	 * @brief Calculates time derivative in case of external dependence
	 * @param [in] t Current time
	 * @param [in] z Axial coordinate in the column
	 * @param [in] r Radial coordinate in the bead
	 * @param [in] secIdx Index of the current section
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 * @param [in,out] workSpace Memory buffer for updated data
	 * @return Time derivatives of externally dependent parameters
	 */
	inline const params_t updateTimeDerivative(double t, double z, double r, unsigned int secIdx, unsigned int nComp, unsigned int const* nBoundStates, void* workSpace) const
	{
		return _localParams;
	}

protected:

	/**
	 * @brief Validates recently read parameters
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 * @return @c true if the parameters were validated successfully, otherwise @c false
	 */
	inline bool validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
	{
		if ((_kA.size() != _kD.size()) || (_kA.size() < nComp))
			throw InvalidParameterException("LIN_KA and LIN_KD have to have the same size");

		return true;
	}

	ConstParams _localParams; //!< Actual parameter data

	// Handlers provide configure(), reserve(), and registerParam() for parameters
	ScalarComponentDependentParameter _kA; //!< Handler for adsorption rate
	ScalarComponentDependentParameter _kD; //!< Handler for desorption rate
};

/**
 * @brief Handles linear binding model parameters that depend on an external function
 */
class ExtLinearParamHandler : public ExternalParamHandlerBase
{
public:

	/**
	 * @brief Holds actual parameter data
	 * @details The parameter data will be stored in a memory buffer. This requires
	 *          control over the location of the data, which is not provided by
	 *          std::vector and util::SlicedVector. There are "local" equivalents,
	 *          util::LocalVector and util::LocalSlicedVector, for these types that
	 *          allow to control the placement of the payload data. A single value (i.e.,
	 *          a single active or double) can simply be put in the struct.
	 */
	struct VariableParams
	{
		util::LocalVector<active> kA; //!< Adsorption rate
		util::LocalVector<active> kD; //!< Desorption rate
	};

	typedef VariableParams params_t;

	ExtLinearParamHandler() CADET_NOEXCEPT { }

	/**
	 * @brief Returns name of the binding model
	 * @return Name of the binding model
	 */
	static const char* identifier() { return "EXT_LINEAR"; }

	/**
	 * @brief Reads parameters and verifies them
	 * @details See IBindingModel::configure() for details.
	 * @param [in] paramProvider IParameterProvider used for reading parameters
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 * @return @c true if the parameters were read and validated successfully, otherwise @c false
	 */
	inline bool configure(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBoundStates)
	{
		_kA.configure("LIN_KA", paramProvider, nComp, nBoundStates);
		_kD.configure("LIN_KD", paramProvider, nComp, nBoundStates);
		
		// Number of externally dependent parameters (2) needs to be given to ExternalParamHandlerBase::configure()
		ExternalParamHandlerBase::configure(paramProvider, 2);
		return validateConfig(nComp, nBoundStates);
	}

	/**
	 * @brief Registers all local parameters in a map for further use
	 * @param [in,out] parameters Map in which the parameters are stored
	 * @param [in] unitOpIdx Index of the unit operation used for registering the parameters
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 */
	inline void registerParameters(std::unordered_map<ParameterId, active*>& parameters, unsigned int unitOpIdx, unsigned int nComp, unsigned int const* nBoundStates)
	{
		_kA.registerParam("LIN_KA", parameters, unitOpIdx, nComp, nBoundStates);
		_kD.registerParam("LIN_KD", parameters, unitOpIdx, nComp, nBoundStates);
	}

	/**
	 * @brief Reserves space in the storage of the parameters
	 * @param [in] numElem Number of elements (total)
	 * @param [in] numSlices Number of slices / binding site types
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 */
	inline void reserve(unsigned int numElem, unsigned int numSlices, unsigned int nComp, unsigned int const* nBoundStates) \
	{
		_kA.reserve(numElem, numSlices, nComp, nBoundStates);
		_kD.reserve(numElem, numSlices, nComp, nBoundStates);
	}

	/**
	 * @brief Updates the parameter cache in order to take the external profile into account
	 * @details The data and the returned value are constructed in the given @p workSpace memory buffer.
	 * @param [in] t Current time
	 * @param [in] z Axial coordinate in the column
	 * @param [in] r Radial coordinate in the bead
	 * @param [in] secIdx Index of the current section
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 * @param [in,out] workSpace Memory buffer for updated data
	 * @return Externally dependent parameter values
	 */
	inline const params_t& update(double t, double z, double r, unsigned int secIdx, unsigned int nComp, unsigned int const* nBoundStates, void* workSpace) const
	{
		// Construct params_t in workSpace
		params_t* const localParams = reinterpret_cast<params_t*>(workSpace);
		new (localParams) params_t;

		// Pointer to buffer for function evaluation
		double* const extFunBuffer = cadet::util::advancePointer<double>(workSpace, sizeof(params_t));

		// Pointer to buffer for actual parameter data (skip buffer for external time derivatives, see cacheSize() below)
		void* buffer = cadet::util::advancePointer<void>(workSpace, sizeof(params_t) + 2 * 2 * sizeof(double));

		// Evaluate external functions in buffer
		evaluateExternalFunctions(t, z, r, secIdx, 2, extFunBuffer);

		// Prepare the buffer for the data, update the data, and advance buffer pointer to next item
		_kA.prepareCache(localParams->kA, buffer);
		_kA.update(util::dataOfLocalVersion(localParams->kA), extFunBuffer[0], nComp, nBoundStates);
		buffer = util::advancePointer(buffer, util::memoryForDataOf(localParams->kA));

		_kD.prepareCache(localParams->kD, buffer);
		_kD.update(util::dataOfLocalVersion(localParams->kD), extFunBuffer[1], nComp, nBoundStates);

		return *localParams;
	}

	/**
	 * @brief Calculates time derivative in case of external dependence
	 * @details It is assumed that update() has just been called with the same @p workSpace before this function is called.
	 *          The time derivatives are constructed in the given @p workSpace memory buffer.
	 * @param [in] t Current time
	 * @param [in] z Axial coordinate in the column
	 * @param [in] r Radial coordinate in the bead
	 * @param [in] secIdx Index of the current section
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 * @param [in,out] workSpace Memory buffer for updated data
	 * @return Time derivatives of externally dependent parameters
	 */
	inline params_t updateTimeDerivative(double t, double z, double r, unsigned int secIdx, unsigned int nComp, unsigned int const* nBoundStates, void* workSpace) const
	{
		VariableParams p;

		// Pointers to buffers and already calculated external function values
		params_t* const localParams = reinterpret_cast<params_t*>(workSpace);
		double* const extFunBuffer = cadet::util::advancePointer<double>(workSpace, sizeof(params_t));
		double* const extDerivBuffer = extFunBuffer + 2;

		// Pointer to buffer for time derivatives
		void* buffer = util::ptrToEndOfData(localParams->kD);

		// Evaluate time derivatives of external functions
		evaluateTimeDerivativeExternalFunctions(t, z, r, secIdx, 2, extDerivBuffer);

		// Prepare the buffer for the data, update the data, and advance buffer pointer to next item
		_kA.prepareCache(p.kA, buffer);
		_kA.updateTimeDerivative(util::dataOfLocalVersion(p.kA), extFunBuffer[0], extDerivBuffer[0], nComp, nBoundStates);
		buffer = util::advancePointer(buffer, util::memoryForDataOf(p.kA));

		_kD.prepareCache(p.kD, buffer);
		_kD.updateTimeDerivative(util::dataOfLocalVersion(p.kD), extFunBuffer[1], extDerivBuffer[1], nComp, nBoundStates);

		return p;
	}

	/**
	 * @brief Returns how much memory is required for caching in bytes
	 * @details Memory size in bytes.
	 * @return Memory size in bytes
	 */
	inline std::size_t cacheSize() const CADET_NOEXCEPT
	{
		// Required buffer memory:
		//  + params_t object
		//  + buffer for external function evaluations (2 parameters)
		//  + buffer for external function time derivative evaluations (2 parameters)
		//  + buffer for actual parameter data (memory for _kA data + memory for _kD data)
		//  + buffer for parameter time derivatives (memory for _kA data + memory for _kD data)
		return sizeof(params_t) + 2 * 2 * sizeof(double) + 2 * (util::memoryForDataOf(_kA.base()) + util::memoryForDataOf(_kD.base()));
	}

protected:

	/**
	 * @brief Validates recently read parameters
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 * @return @c true if the parameters were validated successfully, otherwise @c false
	 */
	inline bool validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
	{
		if ((_kA.size() != _kD.size()) || (_kA.size() < nComp))
			throw InvalidParameterException("LIN_KA and LIN_KD have to have the same size");

		return true;
	}

	// Handlers provide configure(), reserve(), and registerParam() for parameters
	ExternalScalarComponentDependentParameter _kA; //!< Handler for adsorption rate
	ExternalScalarComponentDependentParameter _kD; //!< Handler for desorption rate
};


/**
 * @brief Defines the linear binding model
 * @details Implements the linear adsorption model \f$ \frac{\mathrm{d}q}{\mathrm{d}t} = k_a c_p - k_d q \f$.
 *          Multiple bound states are not supported. Components without bound state (i.e., non-binding components)
 *          are supported.
 *          
 *          See @cite Guiochon2006
 * @tparam ParamHandler_t Type that can add support for external function dependence
 */
template <class ParamHandler_t>
class LinearBindingBase : public IBindingModel
{
public:

	LinearBindingBase() : _nComp(0), _nBoundStates(nullptr) { }
	virtual ~LinearBindingBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }
	virtual const char* name() const CADET_NOEXCEPT { return ParamHandler_t::identifier(); }
	virtual bool requiresConfiguration() const CADET_NOEXCEPT { return true; }

	virtual bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
	{
		_nComp = nComp;
		_nBoundStates = nBound;
		if (hasMultipleBoundStates(nBound, nComp))
			throw InvalidParameterException("Linear binding model does not support multiple bound states");

		return true;
	}

	virtual bool configure(IParameterProvider& paramProvider, unsigned int unitOpIdx)
	{
		_parameters.clear();

		// Read binding dynamics (quasi-stationary, kinetic)
		_kineticBinding = paramProvider.getInt("IS_KINETIC");

		// Read parameters (k_a and k_d)
		_paramHandler.configure(paramProvider, _nComp, _nBoundStates);

		// Register parameters
		_paramHandler.registerParameters(_parameters, unitOpIdx, _nComp, _nBoundStates);

		return true;
	}

	virtual void fillBoundPhaseInitialParameters(ParameterId* params, unsigned int unitOpIdx) const CADET_NOEXCEPT
	{
		unsigned int ctr = 0;
		for (unsigned int c = 0; c < _nComp; ++c)
		{
			for (unsigned int bp = 0; bp < _nBoundStates[c]; ++bp, ++ctr)
				params[ctr] = makeParamId(hashString("INIT_Q"), unitOpIdx, c, bp, ReactionIndep, SectionIndep);
		}
	}

	virtual std::unordered_map<ParameterId, double> getAllParameterValues() const
	{
		std::unordered_map<ParameterId, double> data;
		std::transform(_parameters.begin(), _parameters.end(), std::inserter(data, data.end()),
		               [](const std::pair<const ParameterId, active*>& p) { return std::make_pair(p.first, static_cast<double>(*p.second)); });

		return data;
	}

	virtual bool hasParameter(const ParameterId& pId) const
	{
		return _parameters.find(pId) != _parameters.end();
	}

	virtual bool setParameter(const ParameterId& pId, int value)
	{
		return false;
	}

	virtual bool setParameter(const ParameterId& pId, double value)
	{
		auto paramHandle = _parameters.find(pId);
		if (paramHandle != _parameters.end())
		{
			paramHandle->second->setValue(value);
			return true;
		}

		return false;
	}

	virtual bool setParameter(const ParameterId& pId, bool value)
	{
		return false;
	}

	virtual void getAlgebraicBlock(unsigned int& idxStart, unsigned int& len) const
	{
		idxStart = 0;

		// If we have kinetic binding, then we don't have algebraic equations
		// Otherwise, all equations are algebraic
		if (_kineticBinding)
			len = 0;
		else
			len = numBoundStates(_nBoundStates, _nComp);
	}

	virtual active* getParameter(const ParameterId& pId)
	{
		auto paramHandle = _parameters.find(pId);
		if (paramHandle != _parameters.end())
		{
			return paramHandle->second;
		}

		return nullptr;
	}

	virtual unsigned int workspaceSize() const { return _paramHandler.cacheSize(); }

	virtual void consistentInitialState(double t, double z, double r, unsigned int secIdx, double* const vecStateY, double errorTol, 
		active* const adRes, active* const adY, unsigned int adEqOffset, unsigned int adDirOffset, const ad::IJacobianExtractor& jacExtractor, 
		double* const workingMemory, linalg::detail::DenseMatrixBase& workingMat) const
	{
		// If we have kinetic binding, there are no algebraic equations and we are done
		if (_kineticBinding)
			return;

		const typename ParamHandler_t::params_t& p = _paramHandler.update(t, z, r, secIdx, _nComp, _nBoundStates, workingMemory);

		// Compute the q_i from their corresponding c_{p,i}

		// Pointer to first component in liquid phase
		double const* yCp = vecStateY - _nComp;

		unsigned int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			// Solve  k_a * c_p - k_d * q == 0  for q to obtain  q = k_a / k_d * c_p
			vecStateY[bndIdx] = static_cast<double>(p.kA[i]) / static_cast<double>(p.kD[i]) * yCp[i];

			// Next bound component
			++bndIdx;
		}
	}

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { _paramHandler.setExternalFunctions(extFuns, size); }

	// The next four residual() function implementations, two analyticJacobian() function implementations, and
	// two jacobianAddDiscretized() function implementations are usually hidden behind
	// CADET_BINDINGMODEL_RESIDUAL_BOILERPLATE
	// which just expands to the eight implementations below.

	virtual int residual(const active& t, double z, double r, unsigned int secIdx, const active& timeFactor, 
		active const* y, double const* yDot, active* res, void* workSpace) const
	{
		return residualImpl<active, active, active>(t, z, r, secIdx, timeFactor, y, yDot, res, workSpace);
	}

	virtual int residual(double t, double z, double r, unsigned int secIdx, double timeFactor, 
		active const* y, double const* yDot, active* res, void* workSpace) const
	{
		return residualImpl<active, active, double>(t, z, r, secIdx, timeFactor, y, yDot, res, workSpace);
	}

	virtual int residual(const active& t, double z, double r, unsigned int secIdx, const active& timeFactor, 
		double const* y, double const* yDot, active* res, void* workSpace) const
	{
		return residualImpl<double, active, active>(t, z, r, secIdx, timeFactor, y, yDot, res, workSpace);
	}

	virtual int residual(double t, double z, double r, unsigned int secIdx, double timeFactor, 
		double const* y, double const* yDot, double* res, void* workSpace) const
	{
		return residualImpl<double, double, double>(t, z, r, secIdx, timeFactor, y, yDot, res, workSpace);
	}

	virtual void analyticJacobian(double t, double z, double r, unsigned int secIdx, double const* y, linalg::BandMatrix::RowIterator jac, void* workSpace) const
	{
		jacobianImpl(t, z, r, secIdx, y, jac, workSpace);
	}

	virtual void analyticJacobian(double t, double z, double r, unsigned int secIdx, double const* y, linalg::DenseBandedRowIterator jac, void* workSpace) const
	{
		jacobianImpl(t, z, r, secIdx, y, jac, workSpace);
	}

	virtual void jacobianAddDiscretized(double alpha, linalg::FactorizableBandMatrix::RowIterator jac) const
	{
		jacobianAddDiscretizedImpl(alpha, jac);
	}

	virtual void jacobianAddDiscretized(double alpha, linalg::DenseBandedRowIterator jac) const
	{
		jacobianAddDiscretizedImpl(alpha, jac);
	}

	virtual void multiplyWithDerivativeJacobian(double const* yDotS, double* const res, double timeFactor) const
	{
		// Multiplier is 0 if quasi-stationary and 1 if kinetic binding
		const double multiplier = _kineticBinding ? timeFactor : 0.0;

		const unsigned int eqSize = numBoundStates(_nBoundStates, _nComp);
		for (unsigned int i = 0; i < eqSize; ++i)
		{
			res[i] = multiplier * yDotS[i];
		}
	}

	virtual void timeDerivativeAlgebraicResidual(double t, double z, double r, unsigned int secIdx, double const* y, double* dResDt, void* workSpace) const
	{
		if (!hasAlgebraicEquations())
			return;

		if (!ParamHandler_t::dependsOnTime())
			return;

		// Update external function and compute time derivative of parameters
		const typename ParamHandler_t::params_t& p = _paramHandler.update(t, z, r, secIdx, _nComp, _nBoundStates, workSpace);
		const typename ParamHandler_t::params_t dpDt = _paramHandler.updateTimeDerivative(t, z, r, secIdx, _nComp, _nBoundStates, workSpace);

		// Pointer to first component in liquid phase
		double const* yCp = y - _nComp;

		unsigned int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			dResDt[bndIdx] = -(static_cast<double>(dpDt.kA[i]) * yCp[i] - static_cast<double>(dpDt.kD[i]) * y[bndIdx]);

			// Next bound component
			++bndIdx;
		}
	}

	virtual bool hasSalt() const CADET_NOEXCEPT { return false; }
	virtual bool supportsMultistate() const CADET_NOEXCEPT { return false; }
	virtual bool supportsNonBinding() const CADET_NOEXCEPT { return true; }
	virtual bool hasAlgebraicEquations() const CADET_NOEXCEPT { return !_kineticBinding; }
	virtual bool dependsOnTime() const CADET_NOEXCEPT { return ParamHandler_t::dependsOnTime(); }
	virtual bool requiresWorkspace() const CADET_NOEXCEPT { return ParamHandler_t::requiresWorkspace(); }

protected:
	int _nComp; //!< Number of components
	unsigned int const* _nBoundStates; //!< Array with number of bound states for each component
	bool _kineticBinding; //!< Determines whether binding is kinetic (@c true) or quasi-stationary (@c false)

	ParamHandler_t _paramHandler; //!< Parameters

	std::unordered_map<ParameterId, active*> _parameters; //!< Map used to translate ParameterIds to actual variables

	template <typename StateType, typename ResidualType, typename ParamType>
	int residualImpl(const ParamType& t, double z, double r, unsigned int secIdx, const ParamType& timeFactor,
		StateType const* y, double const* yDot, ResidualType* res, void* workSpace) const
	{
		const typename ParamHandler_t::params_t& p = _paramHandler.update(static_cast<double>(t), z, r, secIdx, _nComp, _nBoundStates, workSpace);

		// Implement k_a * c_{p,i} - k_d * q_i
		// Note that we actually need dq / dt - [k_a * c_{p,i} - k_d * q_i] = 0

		// Pointer to first component in liquid phase
		StateType const* yCp = y - _nComp;

		unsigned int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			res[bndIdx] = -(static_cast<ParamType>(p.kA[i]) * yCp[i] - static_cast<ParamType>(p.kD[i]) * y[bndIdx]);

			// Add time derivative if necessary
			if (_kineticBinding && yDot)
			{
				res[bndIdx] += timeFactor * yDot[bndIdx];
			}

			// Next bound component
			++bndIdx;
		}

		return 0;
	}

	template <typename RowIterator>
	inline void jacobianImpl(double t, double z, double r, unsigned int secIdx, double const* y, RowIterator jac, void* workSpace) const
	{
		const typename ParamHandler_t::params_t& p = _paramHandler.update(t, z, r, secIdx, _nComp, _nBoundStates, workSpace);

		int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			jac[0] = static_cast<double>(p.kD[i]); // dres / dq_i
			jac[i - bndIdx - _nComp] = -static_cast<double>(p.kA[i]); // dres / dc_{p,i}
			// The distance from liquid phase to solid phase is reduced for each non-binding component
			// since a bound state is neglected. The number of neglected bound states so far is i - bndIdx.
			// Thus, by going back nComp - (i - bndIdx) = -[ i - bndIdx - nComp ] we get to the corresponding
			// liquid phase component.

			++bndIdx;
			++jac;
		}
	}

	template <typename RowIterator>
	inline void jacobianAddDiscretizedImpl(double alpha, RowIterator jac) const
	{
		// We only add time derivatives for kinetic binding
		if (!_kineticBinding)
			return;

		// All equations are kinetic
		const unsigned int eqSize = numBoundStates(_nBoundStates, _nComp);
		for (unsigned int i = 0; i < eqSize; ++i, ++jac)
		{
			jac[0] += alpha;
		}
	}
};

typedef LinearBindingBase<LinearParamHandler> LinearBinding;
typedef LinearBindingBase<ExtLinearParamHandler> ExternalLinearBinding;

namespace binding
{
	void registerLinearModel(std::unordered_map<std::string, std::function<model::IBindingModel*()>>& bindings)
	{
		bindings[LinearBinding::identifier()] = []() { return new LinearBinding(); };
		bindings[ExternalLinearBinding::identifier()] = []() { return new ExternalLinearBinding(); };
	}
}  // namespace binding

}  // namespace model

}  // namespace cadet
