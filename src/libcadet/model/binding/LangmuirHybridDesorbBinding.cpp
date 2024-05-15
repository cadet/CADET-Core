// =============================================================================
//  CADET
//  
//  Copyright © 2008-2024: The CADET Authors
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
#include "SimulationTypes.hpp"

#include <functional>
#include <unordered_map>
#include <string>
#include <vector>

/*<codegen>
{
	"name": "LangmuirHDParamHandler",
	"externalName": "ExtLangmuirHDParamHandler",
	"parameters":
		[
			{ "type": "ScalarComponentDependentParameter", "varName": "kA", "confName": "MCL_KA"},
			{ "type": "ScalarComponentDependentParameter", "varName": "kD", "confName": "MCL_KD"},
			{ "type": "ScalarComponentDependentParameter", "varName": "qMax", "confName": "MCL_QMAX"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p1", "confName": "ANN_P1"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p2", "confName": "ANN_P2"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p3", "confName": "ANN_P3"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p4", "confName": "ANN_P4"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p5", "confName": "ANN_P5"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p6", "confName": "ANN_P6"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p7", "confName": "ANN_P7"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p8", "confName": "ANN_P8"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p9", "confName": "ANN_P9"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p10", "confName": "ANN_P10"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p11", "confName": "ANN_P11"},
			{ "type": "ScalarComponentDependentParameter", "varName": "P_norm", "confName": "ANN_P_NORM"},
			{ "type": "ScalarComponentDependentParameter", "varName": "q_norm", "confName": "ANN_Q_NORM"},
			{ "type": "ScalarComponentDependentParameter", "varName": "S_norm", "confName": "ANN_S_NORM"}
		]
}
</codegen>*/

/* Parameter description
 ------------------------
 kA = Adsorption rate
 kD = Desorption rate
 qMax = Capacity
 p1 = ANN param 1
 p2 = ANN param 2
 p3 = ANN param 3
 p4 = ANN param 4
 p5 = ANN param 5
 p6 = ANN param 6
 p7 = ANN param 7
 p8 = ANN param 8
 p9 = ANN param 9
 p10 = ANN param 10
 p11 = ANN param 11
 P_norm = Normative Vaule for Protien pore conc
 q_norm = Normative Vaule for Protien adsorbed conc
 S_norm = Normative Vaule for Salt pore conc
*/

namespace cadet
{

namespace model
{

inline const char* LangmuirHDParamHandler::identifier() CADET_NOEXCEPT { return "MULTI_COMPONENT_LANGMUIR_HD"; }

inline bool LangmuirHDParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_kA.size() != _kD.size()) || (_kA.size() != _qMax.size()) || (_kA.size() != _p1.size()) || (_kA.size() != _p2.size()) 
		|| (_kA.size() != _p3.size()) || (_kA.size() != _p4.size()) || (_kA.size() != _p5.size()) || (_kA.size() != _p6.size()) 
		|| (_kA.size() != _p7.size()) || (_kA.size() != _p8.size()) || (_kA.size() != _p9.size()) || (_kA.size() != _p10.size()) 
		|| (_kA.size() != _p11.size()) || (_kA.size() != _P_norm.size()) || (_kA.size() != _q_norm.size()) || (_kA.size() != _S_norm.size()) || (_kA.size() < nComp))
		throw InvalidParameterException("MCL_KA, MCL_KD, MCL_QMAX, and a ANN parameter have to have the same size");

	return true;
}

inline const char* ExtLangmuirHDParamHandler::identifier() CADET_NOEXCEPT { return "EXT_MULTI_COMPONENT_LANGMUIR_HD"; }

inline bool ExtLangmuirHDParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_kA.size() != _kD.size()) || (_kA.size() != _qMax.size()) || (_kA.size() != _p1.size()) || (_kA.size() != _p2.size())
		|| (_kA.size() != _p3.size()) || (_kA.size() != _p4.size()) || (_kA.size() != _p5.size()) || (_kA.size() != _p6.size())
		|| (_kA.size() != _p7.size()) || (_kA.size() != _p8.size()) || (_kA.size() != _p9.size()) || (_kA.size() != _p10.size())
		|| (_kA.size() != _p11.size()) || (_kA.size() != _P_norm.size()) || (_kA.size() != _q_norm.size()) || (_kA.size() != _S_norm.size()) || (_kA.size() < nComp))
		throw InvalidParameterException("MCL_KA, MCL_KD, MCL_QMAX, and a ANN parameter have to have the same size");

	return true;
}


/**
 * @brief Defines the multi component LangmuirHD binding model
 * @details Implements the LangmuirHD adsorption model: \f[ \begin{align} 
 *              \frac{\mathrm{d}q_i}{\mathrm{d}t} &= k_{a,i} c_{p,i} q_{\text{max},i} \left( 1 - \sum_j \frac{q_j}{q_{\text{max},j}} \right) - k_{d,i} q_i
 *          \end{align} \f]
 *          Multiple bound states are not supported. 
 *          Components without bound state (i.e., non-binding components) are supported.
 *          
 *          See @cite LangmuirHD1916.
 * @tparam ParamHandler_t Type that can add support for external function dependence
 */
template <class ParamHandler_t>
class LangmuirHybridDesorbBindingBase : public ParamHandlerBindingModelBase<ParamHandler_t>
{
public:

	LangmuirHybridDesorbBindingBase() { }
	virtual ~LangmuirHybridDesorbBindingBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }

	virtual bool implementsAnalyticJacobian() const CADET_NOEXCEPT { return true; }

	virtual void timeDerivativeQuasiStationaryFluxes(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yCp, double const* y, double* dResDt, LinearBufferAllocator workSpace) const
	{
		if (!this->hasQuasiStationaryReactions())
			return;

		if (!ParamHandler_t::dependsOnTime())
			return;

		// Update external function and compute time derivative of parameters
		typename ParamHandler_t::ParamsHandle p;
		typename ParamHandler_t::ParamsHandle dpDt;
		std::tie(p, dpDt) = _paramHandler.updateTimeDerivative(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Protein fluxes: -k_{a,i} * c_{p,i} * q_{max,i} * (1 - \sum_j q_j / q_{max,j}) + k_{d,i} * q_i
		double qSum = 1.0;
		double qSumT = 0.0;
		unsigned int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			const double summand = y[bndIdx] / static_cast<double>(p->qMax[i]);
			qSum -= summand;
			qSumT += summand / static_cast<double>(p->qMax[i]) * static_cast<double>(dpDt->qMax[i]);

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
			dResDt[bndIdx] = static_cast<double>(dpDt->kD[i]) * y[bndIdx] 
				- yCp[i] * (static_cast<double>(dpDt->kA[i]) * static_cast<double>(p->qMax[i]) * qSum
				           + static_cast<double>(p->kA[i]) * static_cast<double>(dpDt->qMax[i]) * qSum
				           + static_cast<double>(p->kA[i]) * static_cast<double>(p->qMax[i]) * qSumT);

			// Next bound component
			++bndIdx;
		}
	}

	CADET_BINDINGMODELBASE_BOILERPLATE

protected:
	using ParamHandlerBindingModelBase<ParamHandler_t>::_paramHandler;
	using ParamHandlerBindingModelBase<ParamHandler_t>::_reactionQuasistationarity;
	using ParamHandlerBindingModelBase<ParamHandler_t>::_nComp;
	using ParamHandlerBindingModelBase<ParamHandler_t>::_nBoundStates;

	template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
	int fluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, StateType const* y,
		CpStateType const* yCp, ResidualType* res, LinearBufferAllocator workSpace) const
	{
		typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Protein fluxes: -k_{a,i} * c_{p,i} * q_{max,i} * (1 - \sum_j q_j / q_{max,j}) + k_{d,i} * q_i
		ResidualType qSum = 1.0;
		unsigned int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			qSum -= y[bndIdx] / static_cast<ParamType>(p->qMax[i]);

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
			// ANN
			const double p1 = static_cast<double>(p->p1[i]);
			const double p2 = static_cast<double>(p->p2[i]);
			const double p3 = static_cast<double>(p->p3[i]);
			const double p4 = static_cast<double>(p->p4[i]);
			const double p5 = static_cast<double>(p->p5[i]);
			const double p6 = static_cast<double>(p->p6[i]);
			const double p7 = static_cast<double>(p->p7[i]);
			const double p8 = static_cast<double>(p->p8[i]);
			const double p9 = static_cast<double>(p->p9[i]);
			const double p10 = static_cast<double>(p->p10[i]);
			const double p11 = static_cast<double>(p->p11[i]);
			const double P_norm = static_cast<double>(p->P_norm[i]);
			const double q_norm = static_cast<double>(p->q_norm[i]);
			const double S_norm = static_cast<double>(p->S_norm[i]);

			// const double CF1 = exp(-((p1 * yCp[i] / P_norm) + (p3 * y[bndIdx] / q_norm) + (p5 * yCp[0] / S_norm) + p7));
			// const double CF2 = exp(-((p2 * yCp[i] / P_norm) + (p4 * y[bndIdx] / q_norm) + (p6 * yCp[0] / S_norm) + p8));

			// double ANN = (p11 + (p9 / (1 + CF1)) + (p10 / (1 + CF2)));
			
			res[bndIdx] = static_cast<ParamType>(p->kD[i]) * y[bndIdx] - static_cast<ParamType>(p->kA[i]) * yCp[i] * static_cast<ParamType>(p->qMax[i]) * qSum;

			// Next bound component
			++bndIdx;
		}

		return 0;
	}

	template <typename RowIterator>
	void jacobianImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* y, double const* yCp, int offsetCp, RowIterator jac, LinearBufferAllocator workSpace) const
	{
		typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Protein fluxes: -k_{a,i} * c_{p,i} * q_{max,i} * (1 - \sum_j q_j / q_{max,j}) + k_{d,i} * q_i
		double qSum = 1.0;
		int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			qSum -= y[bndIdx] / static_cast<double>(p->qMax[i]);

			// Next bound component
			++bndIdx;
		}

		bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;
			// Langmuir Params
			const double ka = static_cast<double>(p->kA[i]);
			const double kd = static_cast<double>(p->kD[i]);
			// ANN Params
			const double p1 = static_cast<double>(p->p1[i]);
			const double p2 = static_cast<double>(p->p2[i]);
			const double p3 = static_cast<double>(p->p3[i]);
			const double p4 = static_cast<double>(p->p4[i]);
			const double p5 = static_cast<double>(p->p5[i]);
			const double p6 = static_cast<double>(p->p6[i]);
			const double p7 = static_cast<double>(p->p7[i]);
			const double p8 = static_cast<double>(p->p8[i]);
			const double p9 = static_cast<double>(p->p9[i]);
			const double p10 = static_cast<double>(p->p10[i]);
			const double p11 = static_cast<double>(p->p11[i]);
			const double P_norm = static_cast<double>(p->P_norm[i]);
			const double q_norm = static_cast<double>(p->q_norm[i]);
			const double S_norm = static_cast<double>(p->S_norm[i]);
			// CF = Common Factor
			const double CF1 = exp(-((p1 * yCp[i] / P_norm) + (p3 * y[bndIdx] / q_norm) + (p5 * yCp[0] / S_norm) + p7));
			const double CF2 = exp(-((p2 * yCp[i] / P_norm) + (p4 * y[bndIdx] / q_norm) + (p6 * yCp[0] / S_norm) + p8));
			const double ANN = (p11 + (p9 / (1 + CF1)) + (p10 / (1 + CF2)));
			const double ANN_CF = (ANN / abs(ANN));


			// dres_i / dc_{p,i}
			const double dc_pi_ANN = (((p10 * p2 * CF2) / (P_norm * pow(CF2 + 1, 2))) + ((p1 * p9 * CF1) / (P_norm * pow(CF1 + 1, 2)))) * ANN_CF;
			jac[i - bndIdx - offsetCp] = kd * y[bndIdx] * dc_pi_ANN  - ka * static_cast<double>(p->qMax[i]) * qSum;
			//jac[i - bndIdx - offsetCp] = ka * static_cast<double>(p->qMax[i]) * qSum;
			// Getting to c_{p,i}: -bndIdx takes us to q_0, another -offsetCp to c_{p,0} and a +i to c_{p,i}.
			//                     This means jac[i - bndIdx - offsetCp] corresponds to c_{p,i}.


			// dres_i / dc_{p,0}
			const double dc_p0_ANN = (((p10 * p6 * CF2) / (S_norm * pow(CF2 + 1, 2))) + ((p5 * p9 * CF1) / (S_norm * pow(CF1 + 1, 2)))) * ANN_CF;
			jac[-bndIdx - offsetCp] = kd * y[bndIdx] * dc_p0_ANN;
			// Getting to c_{p,0}: -bndIdx takes us to q_0, another -offsetCp to c_{p,0}.
			//                     This means jac[bndIdx - offsetCp] corresponds to c_{p,0}.




			// Fill dres_i / dq_j
			int bndIdx2 = 0;
			for (int j = 0; j < _nComp; ++j)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[j] == 0)
					continue;

				// dres_i / dq_j
				jac[bndIdx2 - bndIdx] = ka * yCp[i] * static_cast<double>(p->qMax[i]) / static_cast<double>(p->qMax[j]);
				// Getting to q_j: -bndIdx takes us to q_0, another +bndIdx2 to q_j. This means jac[bndIdx2 - bndIdx] corresponds to q_j.

				++bndIdx2;
			}

			// Add to dres_i / dq_i
			const double dc_qi_ANN = (((p10 * p4 * CF2) / (q_norm * pow(CF2 + 1, 2))) + ((p3 * p9 * CF1) / (q_norm * pow(CF1 + 1, 2)))) * ANN_CF;
			jac[0] += kd*abs(ANN) + kd* y[bndIdx] * dc_qi_ANN;
			// jac[0] += kd;
			// Advance to next flux and Jacobian row
			++bndIdx;
			++jac;
		}
	}	
};

typedef LangmuirHybridDesorbBindingBase<LangmuirHDParamHandler> LangmuirHybridDesorbBinding;
typedef LangmuirHybridDesorbBindingBase<ExtLangmuirHDParamHandler> ExternalLangmuirHybridDesorbBinding;

namespace binding
{
	void registerLangmuirHybridDesorbModel(std::unordered_map<std::string, std::function<model::IBindingModel*()>>& bindings)
	{
		bindings[LangmuirHybridDesorbBinding::identifier()] = []() { return new LangmuirHybridDesorbBinding(); };
		bindings[ExternalLangmuirHybridDesorbBinding::identifier()] = []() { return new ExternalLangmuirHybridDesorbBinding(); };
	}
}  // namespace binding

}  // namespace model

}  // namespace cadet
