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
#include "model/binding/ExternalFunctionSupport.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "ParamReaderHelper.hpp"

#include <vector>
#include <unordered_map>
#include <functional>

namespace cadet
{

namespace model
{

/**
 * @brief Handles Anti-Langmuir binding model parameters that do not depend on external functions
 */
struct AntiLangmuirParamHandler : public BindingParamHandlerBase
{
	static const char* identifier() { return "MULTI_COMPONENT_ANTILANGMUIR"; }

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
		readParameterMatrix(kA, paramProvider, "MCAL_KA", nComp, 1);
		readParameterMatrix(kD, paramProvider, "MCAL_KD", nComp, 1);
		readParameterMatrix(qMax, paramProvider, "MCAL_QMAX", nComp, 1);
		readParameterMatrix(antiLangmuir, paramProvider, "MCAL_ANTILANGMUIR", nComp, 1);

		// Check parameters
		if ((kA.size() != kD.size()) || (kA.size() != qMax.size()) || (kA.size() != antiLangmuir.size()) || (kA.size() < nComp))
			throw InvalidParameterException("MCAL_KA, MCAL_KD, MCAL_QMAX, and MCAL_ANTILANGMUIR have to have the same size");

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
		registerComponentBoundStateDependentParam(hashString("MCAL_KA"), parameters, kA, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("MCAL_KD"), parameters, kD, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("MCAL_QMAX"), parameters, qMax, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("MCAL_ANTILANGMUIR"), parameters, antiLangmuir, unitOpIdx);
	}

	std::vector<active> kA; //!< Adsorption rate
	std::vector<active> kD; //!< Desorption rate
	std::vector<active> qMax; //!< Capacity
	std::vector<active> antiLangmuir; //!< Anti-Langmuir factor
};

/**
 * @brief Handles Anti-Langmuir binding model parameters that depend on an external function
 */
struct ExtAntiLangmuirParamHandler : public ExternalBindingParamHandlerBase
{
	static const char* identifier() { return "EXT_MULTI_COMPONENT_ANTILANGMUIR"; }

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
		CADET_READPAR_MATRIX(kA, paramProvider, "MCAL_KA", nComp, 1);
		CADET_READPAR_MATRIX(kD, paramProvider, "MCAL_KD", nComp, 1);
		CADET_READPAR_MATRIX(qMax, paramProvider, "MCAL_QMAX", nComp, 1);
		CADET_READPAR_MATRIX(antiLangmuir, paramProvider, "MCAL_ANTILANGMUIR", nComp, 1);

		return ExternalBindingParamHandlerBase::configure(paramProvider, 4);
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
		CADET_REGPAR_COMPBND_VEC("MCAL_KA", parameters, kA, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("MCAL_KD", parameters, kD, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("MCAL_QMAX", parameters, qMax, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("MCAL_ANTILANGMUIR", parameters, antiLangmuir, unitOpIdx);
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
			CADET_UPDATE_EXTDEP_VARIABLE_BRACES(qMax, i, _extFunBuffer[2]);
			CADET_UPDATE_EXTDEP_VARIABLE_BRACES(antiLangmuir, i, _extFunBuffer[3]);
		}
	}

	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, kA)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, kD)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, qMax)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, antiLangmuir)
};

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

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { _p.setExternalFunctions(extFuns, size); }
	virtual bool dependsOnTime() const CADET_NOEXCEPT { return ParamHandler_t::dependsOnTime(); }

	CADET_PUREBINDINGMODELBASE_BOILERPLATE
	
protected:
	ParamHandler_t _p; //!< Handles parameters and their dependence on external functions

	virtual bool configureImpl(bool reconfigure, IParameterProvider& paramProvider, unsigned int unitOpIdx)
	{
		// Read parameters
		_p.configure(paramProvider, _nComp, _nBoundStates);

		// Register parameters
		_p.registerParameters(_parameters, unitOpIdx, _nComp, _nBoundStates);

		return true;
	}

	template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
	int residualImpl(const ParamType& t, double z, double r, unsigned int secIdx, const ParamType& timeFactor,
		StateType const* y, CpStateType const* yCp, double const* yDot, ResidualType* res) const
	{
		_p.update(static_cast<double>(t), z, r, secIdx, _nComp, _nBoundStates);

		// Protein equations: dq_i / dt - ( k_{a,i} * c_{p,i} * (1 - \sum q_i / q_{max,i}) - k_{d,i} * q_i) == 0
		//               <=>  dq_i / dt == k_{a,i} * c_{p,i} * (1 - \sum q_i / q_{max,i}) - k_{d,i} * q_i
		ResidualType qSum = 1.0;
		unsigned int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			qSum -= static_cast<ParamType>(_p.antiLangmuir[i]) * y[bndIdx] / static_cast<ParamType>(_p.qMax[i]);

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
			res[bndIdx] = static_cast<ParamType>(_p.kD[i]) * y[bndIdx] - static_cast<ParamType>(_p.kA[i]) * yCp[i] * static_cast<ParamType>(_p.qMax[i]) * qSum;

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
	void jacobianImpl(double t, double z, double r, unsigned int secIdx, double const* y, double const* yCp, RowIterator jac) const
	{
		_p.update(t, z, r, secIdx, _nComp, _nBoundStates);

		// Protein equations: dq_i / dt - ( k_{a,i} * c_{p,i} * (1 - \sum q_i / q_{max,i}) - k_{d,i} * q_i) == 0
		double qSum = 1.0;
		int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			qSum -= static_cast<double>(_p.antiLangmuir[i]) * y[bndIdx] / static_cast<double>(_p.qMax[i]);

			// Next bound component
			++bndIdx;
		}

		bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			const double ka = static_cast<double>(_p.kA[i]);
			const double kd = static_cast<double>(_p.kD[i]);

			// dres_i / dc_{p,i}
			jac[i - bndIdx - _nComp] = -ka * static_cast<double>(_p.qMax[i]) * qSum;
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
				jac[bndIdx2 - bndIdx] = ka * yCp[i] * static_cast<double>(_p.antiLangmuir[j]) * static_cast<double>(_p.qMax[i]) / static_cast<double>(_p.qMax[j]);
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
