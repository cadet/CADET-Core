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
 * @brief Handles Langmuir binding model parameters that do not depend on external functions
 */
struct LangmuirParamHandler : public BindingParamHandlerBase
{
	static const char* identifier() { return "MULTI_COMPONENT_LANGMUIR"; }

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
		readParameterMatrix(kA, paramProvider, "MCL_KA", nComp, 1);
		readParameterMatrix(kD, paramProvider, "MCL_KD", nComp, 1);
		readParameterMatrix(qMax, paramProvider, "MCL_QMAX", nComp, 1);

		// Check parameters
		if ((kA.size() != kD.size()) || (kA.size() != qMax.size()) || (kA.size() < nComp))
			throw InvalidParameterException("MCL_KA, MCL_KD, and MCL_QMAX have to have the same size");

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
		registerComponentBoundStateDependentParam(hashString("MCL_KA"), parameters, kA, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("MCL_KD"), parameters, kD, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("MCL_QMAX"), parameters, qMax, unitOpIdx);
	}

	std::vector<active> kA; //!< Adsorption rate
	std::vector<active> kD; //!< Desorption rate
	std::vector<active> qMax; //!< Capacity
};

/**
 * @brief Handles Langmuir binding model parameters that depend on an external function
 */
struct ExtLangmuirParamHandler : public ExternalBindingParamHandlerBase
{
	static const char* identifier() { return "EXT_MULTI_COMPONENT_LANGMUIR"; }

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
		CADET_READPAR_MATRIX(kA, paramProvider, "MCL_KA", nComp, 1);
		CADET_READPAR_MATRIX(kD, paramProvider, "MCL_KD", nComp, 1);
		CADET_READPAR_MATRIX(qMax, paramProvider, "MCL_QMAX", nComp, 1);

		return ExternalBindingParamHandlerBase::configure(paramProvider, 3);
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
		CADET_REGPAR_COMPBND_VEC("MCL_KA", parameters, kA, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("MCL_KD", parameters, kD, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("MCL_QMAX", parameters, qMax, unitOpIdx);
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
			CADET_UPDATE_EXTDEP_VARIABLE_TDIFF_BRACES(qMax, i, _extFunBuffer[2], extTimeDeriv[2]);
		}
	}

	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, kA)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, kD)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, qMax)
};


/**
 * @brief Defines the multi component Langmuir binding model
 * @details Implements the Langmuir adsorption model: \f[ \begin{align} 
 *              \frac{\mathrm{d}q_i}{\mathrm{d}t} &= k_{a,i} c_{p,i} q_{\text{max},i} \left( 1 - \sum_j \frac{q_j}{q_{\text{max},j}} \right) - k_{d,i} q_i
 *          \end{align} \f]
 *          Multiple bound states are not supported. 
 *          Components without bound state (i.e., non-binding components) are supported.
 *          
 *          See @cite Langmuir1916.
 * @tparam ParamHandler_t Type that can add support for external function dependence
 */
template <class ParamHandler_t>
class LangmuirBindingBase : public PureBindingModelBase
{
public:

	LangmuirBindingBase() { }
	virtual ~LangmuirBindingBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }
	virtual const char* name() const CADET_NOEXCEPT { return ParamHandler_t::identifier(); }

	virtual int residual(const active& t, double z, double r, unsigned int secIdx, const active& timeFactor,
		active const* y, double const* yDot, active* res) const;

	virtual int residual(double t, double z, double r, unsigned int secIdx, double timeFactor,
		active const* y, double const* yDot, active* res) const;

	virtual int residual(const active& t, double z, double r, unsigned int secIdx, const active& timeFactor,
		double const* y, double const* yDot, active* res) const;

	virtual int residual(double t, double z, double r, unsigned int secIdx, double timeFactor,
		double const* y, double const* yDot, double* res) const;

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { _p.setExternalFunctions(extFuns, size); }
	virtual bool dependsOnTime() const CADET_NOEXCEPT { return ParamHandler_t::dependsOnTime(); }

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

		// Protein equations: dq_i / dt - ( k_{a,i} * c_{p,i} * q_{max,i} * (1 - \sum_j q_j / q_{max,j}) - k_{d,i} * q_i) == 0
		//               <=>  dq_i / dt == k_{a,i} * c_{p,i} * q_{max,i} * (1 - \sum_j q_j / q_{max,j}) - k_{d,i} * q_i
		double qSum = 1.0;
		double qSumT = 0.0;
		unsigned int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			const double summand = y[bndIdx] / static_cast<double>(_p.qMax[i]);
			qSum -= summand;
			qSumT += summand / static_cast<double>(_p.qMax[i]) * static_cast<double>(dpDt.qMax[i]);

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
			dResDt[bndIdx] = static_cast<double>(dpDt.kD[i]) * y[bndIdx] 
				- yCp[i] * (static_cast<double>(dpDt.kA[i]) * static_cast<double>(_p.qMax[i]) * qSum
				           + static_cast<double>(_p.kA[i]) * static_cast<double>(dpDt.qMax[i]) * qSum
				           + static_cast<double>(_p.kA[i]) * static_cast<double>(_p.qMax[i]) * qSumT);

			// Next bound component
			++bndIdx;
		}
	}

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

		// Protein equations: dq_i / dt - ( k_{a,i} * c_{p,i} * q_{max,i} * (1 - \sum_j q_j / q_{max,j}) - k_{d,i} * q_i) == 0
		//               <=>  dq_i / dt == k_{a,i} * c_{p,i} * q_{max,i} * (1 - \sum_j q_j / q_{max,j}) - k_{d,i} * q_i
		ResidualType qSum = 1.0;
		unsigned int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			qSum -= y[bndIdx] / static_cast<ParamType>(_p.qMax[i]);

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

		// Protein equations: dq_i / dt - ( k_{a,i} * c_{p,i} * q_{max,i} * (1 - \sum_j q_j / q_{max,j}) - k_{d,i} * q_i) == 0
		double qSum = 1.0;
		int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			qSum -= y[bndIdx] / static_cast<double>(_p.qMax[i]);

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

CADET_PUREBINDINGMODELBASE_TEMPLATED_BOILERPLATE_IMPL(LangmuirBindingBase, ParamHandler_t)

typedef LangmuirBindingBase<LangmuirParamHandler> LangmuirBinding;
typedef LangmuirBindingBase<ExtLangmuirParamHandler> ExternalLangmuirBinding;

namespace binding
{
	void registerLangmuirModel(std::unordered_map<std::string, std::function<model::IBindingModel*()>>& bindings)
	{
		bindings[LangmuirBinding::identifier()] = []() { return new LangmuirBinding(); };
		bindings[ExternalLangmuirBinding::identifier()] = []() { return new ExternalLangmuirBinding(); };
	}
}  // namespace binding

}  // namespace model

}  // namespace cadet
