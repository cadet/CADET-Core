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

#include <functional>
#include <unordered_map>
#include <string>
#include <vector>

namespace cadet
{

namespace model
{

/**
 * @brief Handles mobile phase modulator Langmuir binding model parameters that do not depend on external functions
 */
struct MPMLangmuirParamHandler : public BindingParamHandlerBase
{
	static const char* identifier() { return "MOBILE_PHASE_MODULATOR"; }

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
		readParameterMatrix(kA, paramProvider, "MPM_KA", nComp, 1);
		readParameterMatrix(kD, paramProvider, "MPM_KD", nComp, 1);
		readParameterMatrix(qMax, paramProvider, "MPM_QMAX", nComp, 1);
		readParameterMatrix(gamma, paramProvider, "MPM_GAMMA", nComp, 1);
		readParameterMatrix(beta, paramProvider, "MPM_BETA", nComp, 1);

		// Check parameters
		if ((kA.size() != kD.size()) || (kA.size() != qMax.size()) || (kA.size() != gamma.size())
			|| (kA.size() != beta.size()) || (kA.size() < nComp))
			throw InvalidParameterException("MPM_KA, MPM_KD, MPM_QMAX, MPM_GAMMA, and MPM_BETA have to have the same size");

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
		registerComponentBoundStateDependentParam(hashString("MPM_KA"), parameters, kA, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("MPM_KD"), parameters, kD, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("MPM_QMAX"), parameters, qMax, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("MPM_GAMMA"), parameters, gamma, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("MPM_BETA"), parameters, beta, unitOpIdx);
	}

	std::vector<active> kA; //!< Adsorption rate
	std::vector<active> kD; //!< Desorption rate
	std::vector<active> qMax; //!< Capacity
	std::vector<active> gamma; //!< Hydrophobicity
	std::vector<active> beta; //!< Ion exchange characteristics
};

/**
 * @brief Handles mobile phase modulator Langmuir binding model parameters that depend on an external function
 */
struct ExtMPMLangmuirParamHandler : public ExternalBindingParamHandlerBase
{
	static const char* identifier() { return "EXT_MOBILE_PHASE_MODULATOR"; }

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
		CADET_READPAR_MATRIX(kA, paramProvider, "MPM_KA", nComp, 1);
		CADET_READPAR_MATRIX(kD, paramProvider, "MPM_KD", nComp, 1);
		CADET_READPAR_MATRIX(qMax, paramProvider, "MPM_QMAX", nComp, 1);
		CADET_READPAR_MATRIX(gamma, paramProvider, "MPM_GAMMA", nComp, 1);
		CADET_READPAR_MATRIX(beta, paramProvider, "MPM_BETA", nComp, 1);

		return ExternalBindingParamHandlerBase::configure(paramProvider, 5);
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
		CADET_REGPAR_COMPBND_VEC("MPM_KA", parameters, kA, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("MPM_KD", parameters, kD, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("MPM_QMAX", parameters, qMax, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("MPM_GAMMA", parameters, gamma, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("MPM_BETA", parameters, beta, unitOpIdx);
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
			CADET_UPDATE_EXTDEP_VARIABLE_BRACES(gamma, i, _extFunBuffer[3]);
			CADET_UPDATE_EXTDEP_VARIABLE_BRACES(beta, i, _extFunBuffer[4]);
		}
	}

	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, kA)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, kD)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, qMax)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, gamma)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, beta)
};


/**
 * @brief Defines the mobile phase modulator Langmuir binding model
 * @details Implements the mobile phase modulator Langmuir adsorption model: \f[ \begin{align} 
 *              \frac{\mathrm{d}q_0}{\mathrm{d}t} &= 0 \\
 *              \frac{\mathrm{d}q_i}{\mathrm{d}t} &= k_{a,i} c_{p,i} e^{\gamma_i c_{p,0}} q_{\text{max},i} \left( 1 - \sum_j \frac{q_j}{q_{\text{max},j}} \right) - k_{d,i} c_{p,0}^{\beta_i} q_i
 *          \end{align} \f]
 *          While @f$ \gamma @f$ describes hydrophobicity, @f$ \beta @f$ accounts for ion-exchange characteristics.
 *          Multiple bound states are not supported. Component @c 0 is assumed to be salt, which is also assumed to be inert.
 *          Components without bound state (i.e., non-binding components) are supported.
 *          
 *          Note that the first equation is only used if salt (component @c 0) has a bound state.
 *          It is reasonable to set the number of bound states for salt to @c 0 in order to save time and memory.
 *          
 *          See @cite Melander1989 and @cite Karlsson2004.
 * @tparam ParamHandler_t Type that can add support for external function dependence
 */
template <class ParamHandler_t>
class MobilePhaseModulatorLangmuirBindingBase : public PureBindingModelBase
{
public:

	MobilePhaseModulatorLangmuirBindingBase() { }
	virtual ~MobilePhaseModulatorLangmuirBindingBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }
	virtual const char* name() const CADET_NOEXCEPT { return ParamHandler_t::identifier(); }

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { _p.setExternalFunctions(extFuns, size); }

	virtual bool hasSalt() const CADET_NOEXCEPT { return true; }
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
		StateType const* y, CpStateType const* yCp, double const* yDot, ResidualType* res, void* workSpace) const
	{
		_p.update(static_cast<double>(t), z, r, secIdx, _nComp, _nBoundStates);

		// Salt equation: dq_0 / dt == 0
		// Protein equations: dq_i / dt - ( k_{a,i} * exp(\gamma_i * c_{p,0}) * c_{p,i} * q_{max,i} * (1 - \sum q_i / q_{max,i}) - k_{d,i} * c_{p,0}^\beta_i * q_i) == 0
		//               <=>  dq_i / dt == k_{a,i} * exp(\gamma_i * c_{p,0}) * c_{p,i} * q_{max,i} * (1 - \sum q_i / q_{max,i}) - k_{d,i} * c_{p,0}^\beta_i * q_i
		ResidualType qSum = 1.0;
		unsigned int bndIdx = 0;
		if (_nBoundStates[0] == 1)
			bndIdx = 1;

		for (int i = 1; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			qSum -= y[bndIdx] / static_cast<ParamType>(_p.qMax[i]);

			// Next bound component
			++bndIdx;
		}

		// Handle salt equation
		bndIdx = 0;
		if (_nBoundStates[0] == 1)
		{
			if (_kineticBinding)
				res[0] = y[0];
			else
				res[0] = 0.0;

			// Add time derivative if necessary
			if (_kineticBinding && yDot)
				res[0] += timeFactor * yDot[bndIdx];

			bndIdx = 1;
		}

		// Handle protein equations
		for (int i = 1; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			// Residual
			res[bndIdx] = static_cast<ParamType>(_p.kD[i]) * pow(yCp[0], static_cast<ParamType>(_p.beta[i])) * y[bndIdx] - static_cast<ParamType>(_p.kA[i]) * exp(yCp[0] * static_cast<ParamType>(_p.gamma[i])) * yCp[i] * static_cast<ParamType>(_p.qMax[i]) * qSum;

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
		_p.update(t, z, r, secIdx, _nComp, _nBoundStates);
		
		// Salt equation
		int bndIdx = 0;
		if (_nBoundStates[0] == 1)
		{
			if (!_kineticBinding)
				jac[0] = 1.0;
			++jac;
			bndIdx = 1;
		}

		double qSum = 1.0;
		for (int i = 1; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			qSum -= y[bndIdx] / static_cast<double>(_p.qMax[i]);

			// Next bound component
			++bndIdx;
		}

		// Protein equations: dq_i / dt - ( k_{a,i} * exp(\gamma_i * c_{p,0}) * c_{p,i} * q_{max,i} * (1 - \sum q_i / q_{max,i}) - k_{d,i} * c_{p,0}^\beta_i * q_i) == 0
		bndIdx = 0;
		if (_nBoundStates[0] == 1)
			bndIdx = 1;

		for (int i = 1; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			const double gamma = static_cast<double>(_p.gamma[i]);
			const double beta = static_cast<double>(_p.beta[i]);
			const double qMax = static_cast<double>(_p.qMax[i]);
			const double ka = static_cast<double>(_p.kA[i]) * exp(gamma * yCp[0]);
			const double kdRaw = static_cast<double>(_p.kD[i]);

			// dres_i / dc_{p,0}
			jac[-bndIdx - _nComp] = -ka * yCp[i] * qMax * qSum * gamma + kdRaw * beta * y[bndIdx] * pow(yCp[0], beta - 1.0);
			// Getting to c_{p,0}: -bndIdx takes us to q_0, another -nComp to c_{p,0}.
			//                     This means jac[bndIdx - nComp] corresponds to c_{p,0}.

			// dres_i / dc_{p,i}
			jac[i - bndIdx - _nComp] = -ka * qMax * qSum;
			// Getting to c_{p,i}: -bndIdx takes us to q_0, another -nComp to c_{p,0} and a +i to c_{p,i}.
			//                     This means jac[i - bndIdx - nComp] corresponds to c_{p,i}.

			// Fill dres_i / dq_j
			int bndIdx2 = 0;
			if (_nBoundStates[0] == 1)
				bndIdx2 = 1;

			for (int j = 1; j < _nComp; ++j)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[j] == 0)
					continue;

				// dres_i / dq_j
				jac[bndIdx2 - bndIdx] = ka * yCp[i] * qMax / static_cast<double>(_p.qMax[j]);
				// Getting to q_j: -bndIdx takes us to q_0, another +bndIdx2 to q_j. This means jac[bndIdx2 - bndIdx] corresponds to q_j.

				++bndIdx2;
			}

			// Add to dres_i / dq_i
			jac[0] += kdRaw * pow(yCp[0], beta);

			// Advance to next equation and Jacobian row
			++bndIdx;
			++jac;
		}
	}
};


typedef MobilePhaseModulatorLangmuirBindingBase<MPMLangmuirParamHandler> MobilePhaseModulatorLangmuirBinding;
typedef MobilePhaseModulatorLangmuirBindingBase<ExtMPMLangmuirParamHandler> ExternalMobilePhaseModulatorLangmuirBinding;

namespace binding
{
	void registerMobilePhaseModulatorLangmuirModel(std::unordered_map<std::string, std::function<model::IBindingModel*()>>& bindings)
	{
		bindings[MobilePhaseModulatorLangmuirBinding::identifier()] = []() { return new MobilePhaseModulatorLangmuirBinding(); };
		bindings[ExternalMobilePhaseModulatorLangmuirBinding::identifier()] = []() { return new ExternalMobilePhaseModulatorLangmuirBinding(); };
	}
}  // namespace binding

}  // namespace model

}  // namespace cadet
