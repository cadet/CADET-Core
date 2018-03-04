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

#include "model/binding/BindingModelBase.hpp"
#include "model/ExternalFunctionSupport.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "model/Parameters.hpp"
#include "LocalVector.hpp"

#include <vector>
#include <unordered_map>
#include <functional>

namespace cadet
{

namespace model
{

CADET_BINDINGPARAMS(AntiLangmuirParamHandler, ExtAntiLangmuirParamHandler, 
	(ScalarComponentDependentParameter, kA, "MCAL_KA") //!< Adsorption rate
	(ScalarComponentDependentParameter, kD, "MCAL_KD") //!< Desorption rate
	(ScalarComponentDependentParameter, qMax, "MCAL_QMAX") //!< Capacity
	(ScalarComponentDependentParameter, antiLangmuir, "MCAL_ANTILANGMUIR"), //!< Anti-Langmuir factor
);

inline const char* AntiLangmuirParamHandler::identifier() CADET_NOEXCEPT { return "MULTI_COMPONENT_ANTILANGMUIR"; }

inline bool AntiLangmuirParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_kA.size() != _kD.size()) || (_kA.size() != _qMax.size()) || (_kA.size() != _antiLangmuir.size()) || (_kA.size() < nComp))
		throw InvalidParameterException("MCAL_KA, MCAL_KD, MCAL_QMAX, and MCAL_ANTILANGMUIR have to have the same size");

	return true;
}

inline const char* ExtAntiLangmuirParamHandler::identifier() CADET_NOEXCEPT { return "EXT_MULTI_COMPONENT_ANTILANGMUIR"; }

inline bool ExtAntiLangmuirParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_kA.size() != _kD.size()) || (_kA.size() != _qMax.size()) || (_kA.size() != _antiLangmuir.size()) || (_kA.size() < nComp))
		throw InvalidParameterException("MCAL_KA, MCAL_KD, MCAL_QMAX, and MCAL_ANTILANGMUIR have to have the same size");

	return true;
}


/**
 * @brief Defines the multi component Anti-Langmuir binding model
 * @details Implements the Anti-Langmuir adsorption model: \f[ \begin{align} 
 *              \frac{\mathrm{d}q_i}{\mathrm{d}t} &= k_{a,i} c_{p,i} q_{\text{max},i} \left( 1 - \sum_j p_j \frac{q_j}{q_{\text{max},j}} \right) - k_{d,i} q_i,
 *          \end{align} \f]
 *          where the factor @f$ p_j \in \{ -1, 1\} @f$ determines the behaviour (@c -1 results in Anti-Langmuir, @c +1 in standard Langmuir).
 *          Multiple bound states are not supported. 
 *          Components without bound state (i.e., non-binding components) are supported.
 * @tparam ParamHandler_t Type that can add support for external function dependence
 */
template <class ParamHandler_t>
class AntiLangmuirBindingBase : public PureBindingModelBase
{
public:

	AntiLangmuirBindingBase() { }
	virtual ~AntiLangmuirBindingBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }
	virtual const char* name() const CADET_NOEXCEPT { return ParamHandler_t::identifier(); }

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { _paramHandler.setExternalFunctions(extFuns, size); }
	virtual bool dependsOnTime() const CADET_NOEXCEPT { return ParamHandler_t::dependsOnTime(); }

	CADET_PUREBINDINGMODELBASE_BOILERPLATE
	
protected:
	ParamHandler_t _paramHandler; //!< Handles parameters and their dependence on external functions

	virtual unsigned int paramCacheSize() const CADET_NOEXCEPT { return _paramHandler.cacheSize(); }

	virtual bool configureImpl(bool reconfigure, IParameterProvider& paramProvider, unsigned int unitOpIdx)
	{
		// Read parameters
		_paramHandler.configure(paramProvider, _nComp, _nBoundStates);

		// Register parameters
		_paramHandler.registerParameters(_parameters, unitOpIdx, _nComp, _nBoundStates);

		return true;
	}

	template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
	int residualImpl(const ParamType& t, double z, double r, unsigned int secIdx, const ParamType& timeFactor,
		StateType const* y, CpStateType const* yCp, double const* yDot, ResidualType* res, void* workSpace) const
	{
		const typename ParamHandler_t::params_t& p = _paramHandler.update(static_cast<double>(t), z, r, secIdx, _nComp, _nBoundStates, workSpace);

		// Protein equations: dq_i / dt - ( k_{a,i} * c_{p,i} * (1 - \sum q_i / q_{max,i}) - k_{d,i} * q_i) == 0
		//               <=>  dq_i / dt == k_{a,i} * c_{p,i} * (1 - \sum q_i / q_{max,i}) - k_{d,i} * q_i
		ResidualType qSum = 1.0;
		unsigned int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			qSum -= static_cast<ParamType>(p.antiLangmuir[i]) * y[bndIdx] / static_cast<ParamType>(p.qMax[i]);

			// Next bound component
			++bndIdx;
		}

		bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			// Residual
			res[bndIdx] = static_cast<ParamType>(p.kD[i]) * y[bndIdx] - static_cast<ParamType>(p.kA[i]) * yCp[i] * static_cast<ParamType>(p.qMax[i]) * qSum;

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
	void jacobianImpl(double t, double z, double r, unsigned int secIdx, double const* y, double const* yCp, RowIterator jac, void* workSpace) const
	{
		const typename ParamHandler_t::params_t& p = _paramHandler.update(t, z, r, secIdx, _nComp, _nBoundStates, workSpace);

		// Protein equations: dq_i / dt - ( k_{a,i} * c_{p,i} * (1 - \sum q_i / q_{max,i}) - k_{d,i} * q_i) == 0
		double qSum = 1.0;
		int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			qSum -= static_cast<double>(p.antiLangmuir[i]) * y[bndIdx] / static_cast<double>(p.qMax[i]);

			// Next bound component
			++bndIdx;
		}

		bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			const double ka = static_cast<double>(p.kA[i]);
			const double kd = static_cast<double>(p.kD[i]);

			// dres_i / dc_{p,i}
			jac[i - bndIdx - _nComp] = -ka * static_cast<double>(p.qMax[i]) * qSum;
			// Getting to c_{p,i}: -bndIdx takes us to q_0, another -nComp to c_{p,0} and a +i to c_{p,i}.
			//                     This means jac[i - bndIdx - nComp] corresponds to c_{p,i}.

			// Fill dres_i / dq_j
			int bndIdx2 = 0;
			for (int j = 0; j < _nComp; ++j)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[j] == 0)
					continue;

				// dres_i / dq_j
				jac[bndIdx2 - bndIdx] = ka * yCp[i] * static_cast<double>(p.antiLangmuir[j]) * static_cast<double>(p.qMax[i]) / static_cast<double>(p.qMax[j]);
				// Getting to q_j: -bndIdx takes us to q_0, another +bndIdx2 to q_j. This means jac[bndIdx2 - bndIdx] corresponds to q_j.

				++bndIdx2;
			}

			// Add to dres_i / dq_i
			jac[0] += kd;

			// Advance to next equation and Jacobian row
			++bndIdx;
			++jac;
		}
	}

};

typedef AntiLangmuirBindingBase<AntiLangmuirParamHandler> AntiLangmuirBinding;
typedef AntiLangmuirBindingBase<ExtAntiLangmuirParamHandler> ExternalAntiLangmuirBinding;

namespace binding
{
	void registerAntiLangmuirModel(std::unordered_map<std::string, std::function<model::IBindingModel*()>>& bindings)
	{
		bindings[AntiLangmuirBinding::identifier()] = []() { return new AntiLangmuirBinding(); };
		bindings[ExternalAntiLangmuirBinding::identifier()] = []() { return new ExternalAntiLangmuirBinding(); };
	}
}  // namespace binding

}  // namespace model

}  // namespace cadet
