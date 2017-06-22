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

#include "model/binding/BindingModelBase.hpp"
#include "model/binding/ExternalFunctionSupport.hpp"
#include "model/binding/BindingModelMacros.hpp"
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
 * @brief Handles Kumar-Langmuir binding model parameters that do not depend on external functions
 */
struct KumarLangmuirParamHandler : public BindingParamHandlerBase
{
	static const char* identifier() { return "KUMAR_MULTI_COMPONENT_LANGMUIR"; }

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
		temperature = paramProvider.getDouble("KMCL_TEMP");
		readParameterMatrix(kA, paramProvider, "KMCL_KA", nComp, 1);
		readParameterMatrix(kD, paramProvider, "KMCL_KD", nComp, 1);
		readParameterMatrix(kAct, paramProvider, "KMCL_KACT", nComp, 1);
		readParameterMatrix(qMax, paramProvider, "KMCL_QMAX", nComp, 1);
		readParameterMatrix(nu, paramProvider, "KMCL_NU", nComp, 1);

		// Check parameters
		if ((kA.size() != kAct.size()) || (kA.size() != kD.size()) || (kA.size() != qMax.size()) || (kA.size() != nu.size()) || (kA.size() < nComp))
			throw InvalidParameterException("KMCL_KA, KMCL_KD, KMCL_KACT, KMCL_NU, and KMCL_QMAX have to have the same size");

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
		parameters[makeParamId(hashString("KMCL_TEMP"), unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &temperature;
		registerComponentBoundStateDependentParam(hashString("KMCL_KA"), parameters, kA, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("KMCL_KD"), parameters, kD, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("KMCL_KACT"), parameters, kAct, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("KMCL_QMAX"), parameters, qMax, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("KMCL_NU"), parameters, nu, unitOpIdx);
	}

	active temperature; //!< Temperature @f$ T @f$
	std::vector<active> nu; //!< Salt exponents / characteristic charges @f$ \nu @f$
	std::vector<active> kA; //!< Adsorption rate pre-exponential factor / frequency
	std::vector<active> kD; //!< Desorption rate
	std::vector<active> kAct; //!< Activation temperature @f$ k_{\text{act}} = \frac{E_a}{R} @f$
	std::vector<active> qMax; //!< Capacity
};

/**
 * @brief Handles Kumar-Langmuir binding model parameters that depend on an external function
 */
struct ExtKumarLangmuirParamHandler : public ExternalBindingParamHandlerBase
{
	static const char* identifier() { return "EXT_KUMAR_MULTI_COMPONENT_LANGMUIR"; }

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
		CADET_READPAR_SCALAR(temperature, paramProvider, "KMCL_TEMP");
		CADET_READPAR_MATRIX(kA, paramProvider, "KMCL_KA", nComp, 1);
		CADET_READPAR_MATRIX(kD, paramProvider, "KMCL_KD", nComp, 1);
		CADET_READPAR_MATRIX(kAct, paramProvider, "KMCL_KACT", nComp, 1);
		CADET_READPAR_MATRIX(qMax, paramProvider, "KMCL_QMAX", nComp, 1);
		CADET_READPAR_MATRIX(nu, paramProvider, "KMCL_NU", nComp, 1);

		return ExternalBindingParamHandlerBase::configure(paramProvider, 6);
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
		CADET_REGPAR_SCALAR("KMCL_TEMP", parameters, temperature, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("KMCL_KA", parameters, kA, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("KMCL_KD", parameters, kD, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("KMCL_KACT", parameters, kAct, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("KMCL_QMAX", parameters, qMax, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("KMCL_NU", parameters, nu, unitOpIdx);
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
			CADET_UPDATE_EXTDEP_VARIABLE_BRACES(kAct, i, _extFunBuffer[2]);
			CADET_UPDATE_EXTDEP_VARIABLE_BRACES(qMax, i, _extFunBuffer[3]);
			CADET_UPDATE_EXTDEP_VARIABLE_BRACES(nu, i, _extFunBuffer[4]);
		}
		CADET_UPDATE_EXTDEP_VARIABLE(temperature, _extFunBuffer[5]);
	}

	CADET_DEFINE_EXTDEP_VARIABLE(active, temperature)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, kA)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, kD)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, kAct)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, qMax)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, nu)
};


/**
 * @brief Defines the extended multi component Langmuir binding model used by Kumar et al.
 * @details Implements the extended Langmuir adsorption model: \f[ \begin{align} 
 *              \frac{\mathrm{d}q_i}{\mathrm{d}t} &= k_{a,i} \exp\left( \frac{k_{\text{act},i}}{T} \right) c_{p,i} q_{\text{max},i} \left( 1 - \sum_j \frac{q_j}{q_{\text{max},j}} \right) - \left( c_{p,0} \right)^{\nu_i} k_{d,i} q_i
 *          \end{align} \f]
 *          The first component @f$ c_{p,0} @f$ is assumed to be salt and should be set as non-binding (@c 0 bound states).
 *          Multiple bound states are not supported. Components without bound state (i.e., non-binding components) are supported.
 *          
 *          In this model, the true adsorption rate @f$ k_{a,\text{true}} @f$ is governed by the Arrhenius law in order to
 *          take temperature into account @f[ k_{a,\text{true}} = k_{a,i} \exp\left( \frac{k_{\text{act},i}}{T} \right). @f]
 *          Here, @f$ k_{a,i} @f$ is the frequency or pre-exponential factor and @f[ k_{\text{act},i} = \frac{E}{R} @f] is
 *          the activation temperature (@f$ E @f$ denotes the activation energy and @f$ R @f$ the Boltzmann gas constant).
 *          Desorption is modified by salt (component @c 0) which does not bind. The characteristic charge @f$ \nu @f$
 *          of the protein is taken into account by the power law.
 *          See @cite Kumar2015 for details.
 * @tparam ParamHandler_t Type that can add support for external function dependence
 */
template <class ParamHandler_t>
class KumarLangmuirBindingBase : public PureBindingModelBase
{
public:

	KumarLangmuirBindingBase() { }
	virtual ~KumarLangmuirBindingBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }
	virtual const char* name() const CADET_NOEXCEPT { return ParamHandler_t::identifier(); }

	virtual void configureModelDiscretization(unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
	{
		BindingModelBase::configureModelDiscretization(nComp, nBound, boundOffset);

		// Guarantee that salt has no bound state
		if (nBound[0] != 0)
			throw InvalidParameterException("Kumar-Langmuir binding model requires non-binding salt component");
	}

	virtual int residual(const active& t, double z, double r, unsigned int secIdx, const active& timeFactor,
		active const* y, double const* yDot, active* res) const;

	virtual int residual(double t, double z, double r, unsigned int secIdx, double timeFactor,
		active const* y, double const* yDot, active* res) const;

	virtual int residual(const active& t, double z, double r, unsigned int secIdx, const active& timeFactor,
		double const* y, double const* yDot, active* res) const;

	virtual int residual(double t, double z, double r, unsigned int secIdx, double timeFactor,
		double const* y, double const* yDot, double* res) const;

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { _p.setExternalFunctions(extFuns, size); }

	virtual bool hasSalt() const CADET_NOEXCEPT { return true; }	
	virtual bool dependsOnTime() const CADET_NOEXCEPT { return ParamHandler_t::dependsOnTime(); }

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

	virtual int residualCore(double t, double z, double r, unsigned int secIdx, double timeFactor,
		double const* y, double const* yCp, double const* yDot, double* res) const;
	virtual int residualCore(double t, double z, double r, unsigned int secIdx, double timeFactor,
		active const* y, double const* yCp, double const* yDot, active* res) const;

	virtual void analyticJacobianCore(double t, double z, double r, unsigned int secIdx, double const* y, 
		double const* yCp, linalg::BandMatrix::RowIterator jac) const;
	virtual void analyticJacobianCore(double t, double z, double r, unsigned int secIdx, double const* y, 
		double const* yCp, linalg::detail::DenseMatrixBase::RowIterator jac) const;

	template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
	int residualImpl(const ParamType& t, double z, double r, unsigned int secIdx, const ParamType& timeFactor,
		StateType const* y, CpStateType const* yCp, double const* yDot, ResidualType* res) const
	{
		_p.update(static_cast<double>(t), z, r, secIdx, _nComp, _nBoundStates);

		// Protein equations: dq_i / dt - ( k_{a,i} * exp( k_{act,i} / T ) * c_{p,i} * q_{max,i} * (1 - \sum_j q_j / q_{max,j}) - (c_{p,0})^{\nu_i} * k_{d,i} * q_i) == 0
		//               <=>  dq_i / dt == k_{a,i} * exp( k_{act,i} / T ) * c_{p,i} * q_{max,i} * (1 - \sum_j q_j / q_{max,j}) - (c_{p,0})^{\nu_i} * k_{d,i} * q_i
		ResidualType qSum = 1.0;
		unsigned int bndIdx = 0;
		for (int i = 1; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			qSum -= y[bndIdx] / static_cast<ParamType>(_p.qMax[i]);

			// Next bound component
			++bndIdx;
		}

		bndIdx = 0;
		for (int i = 1; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			// Residual
			const ResidualType ka = static_cast<ParamType>(_p.kA[i]) * exp(static_cast<ParamType>(_p.kAct[i]) / static_cast<ParamType>(_p.temperature));
			const ResidualType kd = pow(yCp[0], static_cast<ParamType>(_p.nu[i])) * static_cast<ParamType>(_p.kD[i]);
			res[bndIdx] = kd * y[bndIdx] - ka * yCp[i] * static_cast<ParamType>(_p.qMax[i]) * qSum;

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

		// Protein equations: dq_i / dt - ( k_{a,i} * exp( k_{act,i} / T ) * c_{p,i} * q_{max,i} * (1 - \sum_j q_j / q_{max,j}) - (c_{p,0})^{\nu_i} * k_{d,i} * q_i) == 0
		double qSum = 1.0;
		int bndIdx = 0;
		for (int i = 1; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			qSum -= y[bndIdx] / static_cast<double>(_p.qMax[i]);

			// Next bound component
			++bndIdx;
		}

		bndIdx = 0;
		for (int i = 1; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			const double ka = static_cast<double>(_p.kA[i]) * exp(static_cast<double>(_p.kAct[i]) / static_cast<double>(_p.temperature));
			const double kd = pow(yCp[0], static_cast<double>(_p.nu[i])) * static_cast<double>(_p.kD[i]);

			// dres_i / dc_{p,i}
			jac[i - bndIdx - _nComp] = -ka * static_cast<double>(_p.qMax[i]) * qSum;
			// Getting to c_{p,i}: -bndIdx takes us to q_0, another -nComp to c_{p,0} and a +i to c_{p,i}.
			//                     This means jac[i - bndIdx - nComp] corresponds to c_{p,i}.

			// dres_i / dc_{p,0}
			jac[i - bndIdx - _nComp - 1] = static_cast<double>(_p.nu[i]) * pow(yCp[0], static_cast<double>(_p.nu[i]) - 1.0) * y[bndIdx];

			// Fill dres_i / dq_j
			int bndIdx2 = 0;
			for (int j = 1; j < _nComp; ++j)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[j] == 0)
					continue;

				// dres_i / dq_j
				jac[bndIdx2 - bndIdx] = ka * yCp[i] * static_cast<double>(_p.qMax[i]) / static_cast<double>(_p.qMax[j]);
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

CADET_PUREBINDINGMODELBASE_TEMPLATED_BOILERPLATE_IMPL(KumarLangmuirBindingBase, ParamHandler_t)

typedef KumarLangmuirBindingBase<KumarLangmuirParamHandler> KumarLangmuirBinding;
typedef KumarLangmuirBindingBase<ExtKumarLangmuirParamHandler> ExternalKumarLangmuirBinding;

namespace binding
{
	void registerKumarLangmuirModel(std::unordered_map<std::string, std::function<model::IBindingModel*()>>& bindings)
	{
		bindings[KumarLangmuirBinding::identifier()] = []() { return new KumarLangmuirBinding(); };
		bindings[ExternalKumarLangmuirBinding::identifier()] = []() { return new ExternalKumarLangmuirBinding(); };
	}
}  // namespace binding

}  // namespace model

}  // namespace cadet
