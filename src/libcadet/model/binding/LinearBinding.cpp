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
#include "model/binding/ExternalFunctionSupport.hpp"
#include "ParamIdUtil.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "ParamReaderHelper.hpp"

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
struct LinearParamHandler : public BindingParamHandlerBase
{
	static const char* identifier() { return "LINEAR"; }

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
		// Read k_a and k_d
		readParameterMatrix(kA, paramProvider, "LIN_KA", nComp, 1);
		readParameterMatrix(kD, paramProvider, "LIN_KD", nComp, 1);

		// Check parameters
		if ((kA.size() != kD.size()) || (kA.size() < nComp))
			throw InvalidParameterException("LIN_KA and LIN_KD have to have the same size");

		return true;
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
		registerComponentBoundStateDependentParam(hashString("LIN_KA"), parameters, kA, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("LIN_KD"), parameters, kD, unitOpIdx);
	}

	std::vector<active> kA; //!< Adsorption rate
	std::vector<active> kD; //!< Desorption rate
};

/**
 * @brief Handles linear binding model parameters that depend on an external function
 */
struct ExtLinearParamHandler : public ExternalBindingParamHandlerBase
{
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
		CADET_READPAR_MATRIX(kA, paramProvider, "LIN_KA", nComp, 1);
		CADET_READPAR_MATRIX(kD, paramProvider, "LIN_KD", nComp, 1);

		return ExternalBindingParamHandlerBase::configure(paramProvider, 2);
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
		CADET_REGPAR_COMPBND_VEC("LIN_KA", parameters, kA, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("LIN_KD", parameters, kD, unitOpIdx);
	}

	/**
	 * @brief Updates local parameter cache in order to take the external profile into account
	 * @details This function is declared const since the actual parameters are left unchanged by the method.
	 *         The cache is marked as mutable in order to make it writable.
	 * @param [in] t Current time
	 * @param [in] z Axial coordinate in the column
	 * @param [in] r Radial coordinate in the bead
	 * @param [in] secIdx Index of the current section
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 */
	inline void update(double t, double z, double r, unsigned int secIdx, unsigned int nComp, unsigned int const* nBoundStates) const
	{
		evaluateExternalFunctions(t, z, r, secIdx);
		for (unsigned int i = 0; i < nComp; ++i)
		{
			CADET_UPDATE_EXTDEP_VARIABLE_BRACES(kA, i, _extFunBuffer[0]);
			CADET_UPDATE_EXTDEP_VARIABLE_BRACES(kD, i, _extFunBuffer[1]);
		}
	}

	/**
	 * @brief Updates local parameter cache and calculates time derivative in case of external dependence
	 * @details This function is declared const since the actual parameters are left unchanged by the method.
	 *         The cache is marked as mutable in order to make it writable.
	 * @param [in] t Current time
	 * @param [in] z Axial coordinate in the column
	 * @param [in] r Radial coordinate in the bead
	 * @param [in] secIdx Index of the current section
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 */
	inline void updateTimeDerivative(double t, double z, double r, unsigned int secIdx, unsigned int nComp, unsigned int const* nBoundStates) const
	{
		const std::vector<double> extTimeDeriv = evaluateTimeDerivativeExternalFunctions(t, z, r, secIdx);
		for (unsigned int i = 0; i < nComp; ++i)
		{
			CADET_UPDATE_EXTDEP_VARIABLE_TDIFF_BRACES(kA, i, _extFunBuffer[0], extTimeDeriv[0]);
			CADET_UPDATE_EXTDEP_VARIABLE_TDIFF_BRACES(kD, i, _extFunBuffer[1], extTimeDeriv[1]);
		}
	}

	// Create all variables: The CADET_DEFINE_EXTDEP_VARIABLE(TYPE, BASENAME) macro creates the variables
	// <BASENAME> (mutable), <BASENAME>T0, <BASENAME>T1, <BASENAME>T2, and <BASENAME>T3
	// which are all of the same type TYPE. The first variable (<BASENAME>) is used as cache for the actual
	// value of the parameter which is computed from the coefficients <BASENAME>T0 - <BASENAME>T3 and the
	// current value of the external source (depending on time and axial position). It is marked as mutable
	// since it acts as a cache and has to be changed in update().

	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, kA)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, kD)
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

	virtual void configureModelDiscretization(unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
	{
		_nComp = nComp;
		_nBoundStates = nBound;
		if (hasMultipleBoundStates(nBound, nComp))
			throw InvalidParameterException("Linear binding model does not support multiple bound states");
	}

	virtual bool configure(IParameterProvider& paramProvider, unsigned int unitOpIdx)
	{
		// Read binding dynamics (quasi-stationary, kinetic)
		_kineticBinding = paramProvider.getInt("IS_KINETIC");

		// Read parameters (k_a and k_d)
		_p.configure(paramProvider, _nComp, _nBoundStates);

		// Register parameters
		_p.registerParameters(_parameters, unitOpIdx, _nComp, _nBoundStates);

		return true;
	}

	virtual bool reconfigure(IParameterProvider& paramProvider, unsigned int unitOpIdx)
	{
		// Configure again
		_parameters.clear();
		return configure(paramProvider, unitOpIdx);
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

	virtual unsigned int consistentInitializationWorkspaceSize() const { return 0; }

	virtual void consistentInitialState(double t, double z, double r, unsigned int secIdx, double* const vecStateY, double errorTol, 
		active* const adRes, active* const adY, unsigned int adEqOffset, unsigned int adDirOffset, const ad::IJacobianExtractor& jacExtractor, 
		double* const workingMemory, linalg::detail::DenseMatrixBase& workingMat) const
	{
		// If we have kinetic binding, there are no algebraic equations and we are done
		if (_kineticBinding)
			return;

		_p.update(t, z, r, secIdx, _nComp, _nBoundStates);

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
			vecStateY[bndIdx] = static_cast<double>(_p.kA[i]) / static_cast<double>(_p.kD[i]) * yCp[i];

			// Next bound component
			++bndIdx;
		}
	}

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { _p.setExternalFunctions(extFuns, size); }

	// The next four residual() function implementations, two analyticJacobian() function implementations, and
	// two jacobianAddDiscretized() function implementations are usually hidden behind
	// CADET_BINDINGMODEL_RESIDUAL_BOILERPLATE
	// which just expands to the eight implementations below.

	virtual int residual(const active& t, double z, double r, unsigned int secIdx, const active& timeFactor, 
		active const* y, double const* yDot, active* res) const
	{
		return residualImpl<active, active, active>(t, z, r, secIdx, timeFactor, y, yDot, res);
	}

	virtual int residual(double t, double z, double r, unsigned int secIdx, double timeFactor, 
		active const* y, double const* yDot, active* res) const
	{
		return residualImpl<active, active, double>(t, z, r, secIdx, timeFactor, y, yDot, res);
	}

	virtual int residual(const active& t, double z, double r, unsigned int secIdx, const active& timeFactor, 
		double const* y, double const* yDot, active* res) const
	{
		return residualImpl<double, active, active>(t, z, r, secIdx, timeFactor, y, yDot, res);
	}

	virtual int residual(double t, double z, double r, unsigned int secIdx, double timeFactor, 
		double const* y, double const* yDot, double* res) const
	{
		return residualImpl<double, double, double>(t, z, r, secIdx, timeFactor, y, yDot, res);
	}

	virtual void analyticJacobian(double t, double z, double r, unsigned int secIdx, double const* y, linalg::BandMatrix::RowIterator jac) const
	{
		jacobianImpl(t, z, r, secIdx, y, jac);
	}

	virtual void analyticJacobian(double t, double z, double r, unsigned int secIdx, double const* y, linalg::DenseBandedRowIterator jac) const
	{
		jacobianImpl(t, z, r, secIdx, y, jac);
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

	virtual void timeDerivativeAlgebraicResidual(double t, double z, double r, unsigned int secIdx, double const* y, double* dResDt) const
	{
		if (!hasAlgebraicEquations())
			return;

		if (!ParamHandler_t::dependsOnTime())
			return;

		// Update external function and compute time derivative of parameters
		_p.update(t, z, r, secIdx, _nComp, _nBoundStates);
		ParamHandler_t dpDt = _p;
		dpDt.updateTimeDerivative(t, z, r, secIdx, _nComp, _nBoundStates);

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

protected:
	int _nComp; //!< Number of components
	unsigned int const* _nBoundStates; //!< Array with number of bound states for each component
	bool _kineticBinding; //!< Determines whether binding is kinetic (@c true) or quasi-stationary (@c false)

	ParamHandler_t _p; //!< Parameters

	std::unordered_map<ParameterId, active*> _parameters; //!< Map used to translate ParameterIds to actual variables

	template <typename StateType, typename ResidualType, typename ParamType>
	int residualImpl(const ParamType& t, double z, double r, unsigned int secIdx, const ParamType& timeFactor,
		StateType const* y, double const* yDot, ResidualType* res) const
	{
		_p.update(static_cast<double>(t), z, r, secIdx, _nComp, _nBoundStates);

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

			res[bndIdx] = -(static_cast<ParamType>(_p.kA[i]) * yCp[i] - static_cast<ParamType>(_p.kD[i]) * y[bndIdx]);

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
	inline void jacobianImpl(double t, double z, double r, unsigned int secIdx, double const* y, RowIterator jac) const
	{
		_p.update(t, z, r, secIdx, _nComp, _nBoundStates);

		int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			jac[0] = static_cast<double>(_p.kD[i]); // dres / dq_i
			jac[i - bndIdx - _nComp] = -static_cast<double>(_p.kA[i]); // dres / dc_{p,i}
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
