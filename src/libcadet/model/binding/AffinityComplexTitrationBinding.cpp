// =============================================================================
//  CADET
//  
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
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
#include "MathUtil.hpp"

#include <cmath>
#include <functional>
#include <unordered_map>
#include <string>
#include <vector>

namespace
{
	constexpr double ln10 = 2.30258509299404568402;
	constexpr double minSaltConcentration = 1e-30;
	constexpr double softplusTransitionFactor = 20.0;
	constexpr double minApparentCapacity = 1e-12;

	template <typename T>
	inline T stableActGate(T const& eta, T const& pKa, T const& saltAxis)
	{
		using std::pow;
		return 1.0 / (1.0 + pow(10.0, eta * (pKa - saltAxis)));
	}

	template <>
	inline double stableActGate<double>(double const& eta, double const& pKa, double const& saltAxis)
	{
		const double x = ln10 * eta * (pKa - saltAxis);
		if (x > 0.0)
		{
			const double e = std::exp(-x);
			return e / (1.0 + e);
		}

		const double e = std::exp(x);
		return 1.0 / (1.0 + e);
	}

	inline double stableActGateDerivSalt(double const eta, double const gate)
	{
		return ln10 * eta * gate * (1.0 - gate);
	}

	template <typename T>
	inline T softplusScaled(T const& x, T const& beta)
	{
		using std::exp;
		using std::log;
		return log(1.0 + exp(beta * x)) / beta;
	}

	template <>
	inline double softplusScaled<double>(double const& x, double const& beta)
	{
		const double bx = beta * x;
		if (bx > 0.0)
			return (bx + std::log1p(std::exp(-bx))) / beta;

		return std::log1p(std::exp(bx)) / beta;
	}

	template <typename T>
	inline T sigmoid(T const& x)
	{
		using std::exp;
		return 1.0 / (1.0 + exp(-x));
	}

	template <>
	inline double sigmoid<double>(double const& x)
	{
		if (x >= 0.0)
		{
			const double e = std::exp(-x);
			return 1.0 / (1.0 + e);
		}

		const double e = std::exp(x);
		return e / (1.0 + e);
	}
}

/*<codegen>
{
	"name": "AffinityComplexTitrationParamHandler",
	"externalName": "ExtAffinityComplexTitrationParamHandler",
	"parameters":
		[
			{ "type": "ScalarComponentDependentParameter", "varName": "kA", "confName": "ACT_KA"},
			{ "type": "ScalarComponentDependentParameter", "varName": "kD", "confName": "ACT_KD"},
			{ "type": "ScalarComponentDependentParameter", "varName": "qMax", "confName": "ACT_QMAX"},
			{ "type": "ScalarComponentDependentParameter", "varName": "etaA", "confName": "ACT_ETAA"},
			{ "type": "ScalarComponentDependentParameter", "varName": "etaG", "confName": "ACT_ETAG"},
			{ "type": "ScalarComponentDependentParameter", "varName": "pKaA", "confName": "ACT_PKAA"},
			{ "type": "ScalarComponentDependentParameter", "varName": "pKaG", "confName": "ACT_PKAG"}
		],
	"constantParameters":
		[
			{ "type": "ScalarBoolParameter", "varName": "usePSalt", "confName": "ACT_USE_P_SALT"}
		]
}
</codegen>*/

/* Parameter description
 ------------------------
 All these parameters also need to be defined for pH
 kA = Adsorption rate
 kD = Desorption rate
 qMax = Binding apacity
 etaA = Slope for the binding capacity changes against pH
 pKaA = Center for the binding capacity changes against pH
 etaG = Slope for the equilibrium constant changes against pH
 pKaG = Center for the equilibrium constant changes against pH
*/

namespace cadet
{

	namespace model
	{

		inline const char* AffinityComplexTitrationParamHandler::identifier() CADET_NOEXCEPT { return "AFFINITY_COMPLEX_TITRATION"; }

		inline bool AffinityComplexTitrationParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
		{
			if ((_kA.size() != _kD.size()) || (_kA.size() != _qMax.size()) || (_kA.size() < nComp) || (_kA.size() != _etaA.size()) || (_kA.size() != _pKaA.size()) || (_kA.size() != _etaG.size()) || (_kA.size() != _pKaG.size()))
				throw InvalidParameterException("ACT_KA, ACT_KD, ACT_QMAX, ACT_ETAA, ACT_PKAA, ACT_ETAG and ACT_PKAG have to have the same size");

			return true;
		}

		inline const char* ExtAffinityComplexTitrationParamHandler::identifier() CADET_NOEXCEPT { return "EXT_AFFINITY_COMPLEX_TITRATION"; }

		inline bool ExtAffinityComplexTitrationParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
		{
			if ((_kA.size() != _kD.size()) || (_kA.size() != _qMax.size()) || (_kA.size() < nComp) || (_kA.size() != _etaA.size()) || (_kA.size() != _pKaA.size()) || (_kA.size() != _etaG.size()) || (_kA.size() != _pKaG.size()))
				throw InvalidParameterException("ACT_KA, ACT_KD, ACT_QMAX, ACT_ETAA, ACT_PKAA, ACT_ETAG and ACT_PKAG have to have the same size");

			return true;
		}


		/**
		 * @brief Defines the affinity complex titration binding model
		 */
		template <class ParamHandler_t>
		class AffinityComplexTitrationBindingBase : public ParamHandlerBindingModelBase<ParamHandler_t>
		{
		public:

			AffinityComplexTitrationBindingBase() { }
			virtual ~AffinityComplexTitrationBindingBase() CADET_NOEXCEPT { }

			static const char* identifier() { return ParamHandler_t::identifier(); }

			virtual bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
			{
				const bool res = BindingModelBase::configureModelDiscretization(paramProvider, nComp, nBound, boundOffset);

				// Guarantee that pH is not a bound speies
				if (nBound[0] != 0)
					throw InvalidParameterException("Affinity complex titration binding model requires the first component (pSalt or salt concentration) to be non-binding");

				for (int i = 0; i < nComp; ++i)
				{
					if (nBound[i] > 1)
						throw InvalidParameterException("Currently the ACT isotherm model supports at most one bound state per component");
				}

				return res;
			}

			virtual bool implementsAnalyticJacobian() const CADET_NOEXCEPT { return true; }

			CADET_BINDINGMODELBASE_BOILERPLATE

		protected:
			using ParamHandlerBindingModelBase<ParamHandler_t>::_paramHandler;
			using ParamHandlerBindingModelBase<ParamHandler_t>::_reactionQuasistationarity;
			using ParamHandlerBindingModelBase<ParamHandler_t>::_nComp;
			using ParamHandlerBindingModelBase<ParamHandler_t>::_nBoundStates;

			// In the follwing the class method the binding model mass transfer behavior is implemented. 
			template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
			int fluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, StateType const* y,
				CpStateType const* yCp, ResidualType* res, LinearBufferAllocator workSpace) const
			{
				typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

				using CpStateParamType = typename DoubleActivePromoter<CpStateType, ParamType>::type;
				using ResParamType = typename DoubleActivePromoter<ResidualType, ParamType>::type;
				using std::log;

				// Protein fluxes
				ResidualType qSum = 0.0;
				unsigned int bndIdx = 0;
				const CpStateParamType saltAxis = _paramHandler.usePSalt().get() ? yCp[0] : -log(yCp[0] + static_cast<ParamType>(minSaltConcentration)) / static_cast<ParamType>(ln10);
				const ParamType saltAxisParam = static_cast<ParamType>(saltAxis);
				for (int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states (bound state index bndIdx is not advanced)
					if (_nBoundStates[i] == 0)
						continue;

					// y is only defined for components that have a bound state
					qSum += y[bndIdx];

					// Next bound component
					++bndIdx;
				}

				// the first component is pH
				bndIdx = 0;
				for (int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states (bound state index bndIdx is not advanced)
					if (_nBoundStates[i] == 0)
						continue;

					// Residual
					const ParamType pKaAaxis = static_cast<ParamType>(p->pKaA[i]);
					const ParamType pKaGaxis = static_cast<ParamType>(p->pKaG[i]);
					const ResParamType f_A = stableActGate(static_cast<ParamType>(p->etaA[i]), pKaAaxis, saltAxisParam);
					const ResParamType f_G = stableActGate(static_cast<ParamType>(p->etaG[i]), pKaGaxis, saltAxisParam);

					const ResParamType qApp = static_cast<ParamType>(p->qMax[i]) * f_A;
					const ResParamType qFree_local = qApp - qSum;
					const ResParamType qScale = (qApp > static_cast<ParamType>(minApparentCapacity)) ? qApp : static_cast<ParamType>(minApparentCapacity);
					const ResParamType beta = static_cast<ParamType>(softplusTransitionFactor) / qScale;
					const ResParamType qFree_eff = softplusScaled(qFree_local, beta);
					const ResParamType kA_local = static_cast<ParamType>(p->kA[i]) * f_G;
					res[bndIdx] = static_cast<ResidualType>(static_cast<ParamType>(p->kD[i]) * y[bndIdx] - kA_local * yCp[i] * qFree_eff);

					// Next bound component
					++bndIdx;
				}

				return 0;
			}

			template <typename RowIterator>
			void jacobianImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* y, double const* yCp, int offsetCp, RowIterator jac, LinearBufferAllocator workSpace) const
			{
				typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

				// Protein flux
				double qsum = 0.0;
				int bndIdx = 0;
				const double saltAxis = _paramHandler.usePSalt().get() ? yCp[0] : -std::log(yCp[0] + minSaltConcentration) / ln10;
				const double dSaltAxis_dC0 = _paramHandler.usePSalt().get() ? 1.0 : -1.0 / ((yCp[0] + minSaltConcentration) * ln10);
				for (int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states
					if (_nBoundStates[i] == 0)
						continue;

					qsum += y[bndIdx];

					// Next bound component
					++bndIdx;
				}

				bndIdx = 0;
				for (int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states (bound state index bndIdx is not advanced)
					if (_nBoundStates[i] == 0)
						continue;

					// local variables to aid the calculation of the jacobian
					const double pKaAaxis = static_cast<double>(p->pKaA[i]);
					const double pKaGaxis = static_cast<double>(p->pKaG[i]);
					const double f_A = stableActGate(static_cast<double>(p->etaA[i]), pKaAaxis, saltAxis);
					const double f_G = stableActGate(static_cast<double>(p->etaG[i]), pKaGaxis, saltAxis);

					const double f_A_deriv = stableActGateDerivSalt(static_cast<double>(p->etaA[i]), f_A);
					const double f_G_deriv = stableActGateDerivSalt(static_cast<double>(p->etaG[i]), f_G);

					const double qApp = static_cast<double>(p->qMax[i]) * f_A;
					const double qFree_local = qApp - qsum;
					const double qScale = (qApp > minApparentCapacity) ? qApp : minApparentCapacity;
					const double beta = softplusTransitionFactor / qScale;
					const double qFree_eff = softplusScaled(qFree_local, beta);
					const double dqFreeEff_dqFree = sigmoid(beta * qFree_local);

					const double qmax_times_ka = static_cast<double>(p->kA[i]) * static_cast<double>(p->qMax[i]);
					const double kA_times_fG = static_cast<double>(p->kA[i]) * f_G;

					// dres_i / dc_{p,i}
					jac[i - bndIdx - offsetCp] = -kA_times_fG * qFree_eff;
					// Getting to c_{p,i}: -bndIdx takes us to q_0, another -offsetCp to c_{p,0} and a +i to c_{p,i}.
					//                     This means jac[i - bndIdx - offsetCp] corresponds to c_{p,i}.

					// dres_i / d(salt axis), the first component is the selected salt representation
					jac[-bndIdx - offsetCp] = -(qmax_times_ka * f_G * f_A_deriv * dqFreeEff_dqFree + static_cast<double>(p->kA[i]) * f_G_deriv * qFree_eff) * yCp[i] * dSaltAxis_dC0;

					// Fill dres_i / dq_j
					int bndIdx2 = 0;
					for (int j = 0; j < _nComp; ++j)
					{
						// Skip components without bound states (bound state index bndIdx2 is not advanced)
						if (_nBoundStates[j] == 0)
							continue;

						// dres_i / dq_j, part 1
						jac[bndIdx2 - bndIdx] = kA_times_fG * yCp[i] * dqFreeEff_dqFree;
						// Getting to q_j: -bndIdx takes us to q_0, another +bndIdx2 to q_j. This means jac[bndIdx2 - bndIdx] corresponds to q_j.

						++bndIdx2;
					}

					// Add to dres_i / dq_i, part 2
					jac[0] += static_cast<double>(p->kD[i]);

					// Advance to next flux and Jacobian row
					++bndIdx;
					++jac;
				}
			}
		};

		// Do not forget to make changes in the following lines of code according to the guidelines given in the tutorial.
		typedef AffinityComplexTitrationBindingBase<AffinityComplexTitrationParamHandler> AffinityComplexTitrationBinding;
		typedef AffinityComplexTitrationBindingBase<ExtAffinityComplexTitrationParamHandler> ExternalAffinityComplexTitrationBinding;

		namespace binding
		{
			void registerAffinityComplexTitrationModel(std::unordered_map<std::string, std::function<model::IBindingModel* ()>>& bindings)
			{
				bindings[AffinityComplexTitrationBinding::identifier()] = []() { return new AffinityComplexTitrationBinding(); };
				bindings[ExternalAffinityComplexTitrationBinding::identifier()] = []() { return new ExternalAffinityComplexTitrationBinding(); };
			}
		}  // namespace binding

	}  // namespace model

}  // namespace cadet