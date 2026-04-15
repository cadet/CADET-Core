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

// Helper functions and constants for the AffinityComplexTitrationBinding model
namespace
{
	constexpr double ln10 = 2.30258509299404568402;
	constexpr double minIonConcentration = 1e-30;
	constexpr double softplusTransitionFactor = 100.0;
	constexpr double minApparentCapacity = 1e-20;

	// Here we use a numerically stable implementation of the logistic function f to avoid potential over and underflow.
	template <typename T>
	inline T stableActGate(T const& eta, T const& pKa, T const& saltAxis)
	{
		using std::exp;
		const T x = static_cast<T>(ln10) * eta * (pKa - saltAxis);
		if (x > static_cast<T>(0.0))
		{
			const T e = exp(-x);
			return e / (static_cast<T>(1.0) + e);
		}

		const T e = exp(x);
		return static_cast<T>(1.0) / (static_cast<T>(1.0) + e);
	}

	inline double stableActGateDerivIon(double const eta, double const gate)
	{
		return ln10 * eta * gate * (1.0 - gate);
	}

	// The softplus function is used to ensure that the apparent binding capacity does not become negative, which will cause IDAS to linger. 
	// It happens for certain param sets that decrease the capacity too abruptly in a multicomponent system. 
	// The transition factor beta determines how close to zero the apparent capacity can get before the function starts to deviate from the identity function.
	// An emipirical value of 100 is chosen here, which means that the apparent capacity can get as low as 1/100 = 1% of x before the softplus function starts to significantly increase it.
	// Approximation error near the clamp is roughly ln2 / beta = 0.7% x.
	template <typename T>
	inline T softplusScaled(T const& x, T const& beta)
	{
		using std::exp;
		using std::log;
		const T bx = beta * x;
		if (bx > static_cast<T>(0.0))
			return (bx + log(static_cast<T>(1.0) + exp(-bx))) / beta;

		return log(static_cast<T>(1.0) + exp(bx)) / beta;
	}

	template <typename T>
	inline T sigmoid(T const& x)
	{
		using std::exp;
		if (x >= static_cast<T>(0.0))
		{
			const T e = exp(-x);
			return static_cast<T>(1.0) / (static_cast<T>(1.0) + e);
		}

		const T e = exp(x);
		return e / (static_cast<T>(1.0) + e);
	}

	inline std::vector<double> getDoubleArrayWithFallback(cadet::IParameterProvider& paramProvider, const char* primaryName, const char* fallbackName)
	{
		if (paramProvider.exists(primaryName))
			return paramProvider.getDoubleArray(primaryName);
		if (paramProvider.exists(fallbackName))
			return paramProvider.getDoubleArray(fallbackName);
		return std::vector<double>();
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
			{ "type": "ScalarComponentDependentParameter", "varName": "etaG", "confName": "ACT_ETAG"}
		]
}
</codegen>*/

/* Parameter description
 ------------------------
 All these parameters also need to be defined for pH/pIon, but not used in the computation.
 kA = Adsorption rate
 kD = Desorption rate
 qMax = Binding apacity
 etaA = Slope for the binding capacity changes against pH/pIon
 pKaA = Center for the binding capacity changes against pH/pIon
 etaG = Slope for the equilibrium constant changes against pH/pIon
 pKaG = Center for the equilibrium constant changes against pH/pIon
*/

namespace cadet
{

	namespace model
	{

		inline const char* AffinityComplexTitrationParamHandler::identifier() CADET_NOEXCEPT { return "AFFINITY_COMPLEX_TITRATION"; }

		inline bool AffinityComplexTitrationParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
		{
			if ((_kA.size() != _kD.size()) || (_kA.size() != _qMax.size()) || (_kA.size() < nComp) || (_kA.size() != _etaA.size()) || (_kA.size() != _etaG.size()))
				throw InvalidParameterException("ACT_KA, ACT_KD, ACT_QMAX, ACT_ETAA and ACT_ETAG have to have the same size");

			return true;
		}

		inline const char* ExtAffinityComplexTitrationParamHandler::identifier() CADET_NOEXCEPT { return "EXT_AFFINITY_COMPLEX_TITRATION"; }

		inline bool ExtAffinityComplexTitrationParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
		{
			if ((_kA.size() != _kD.size()) || (_kA.size() != _qMax.size()) || (_kA.size() < nComp) || (_kA.size() != _etaA.size()) || (_kA.size() != _etaG.size()))
				throw InvalidParameterException("ACT_KA, ACT_KD, ACT_QMAX, ACT_ETAA and ACT_ETAG have to have the same size");

			return true;
		}


		/**
		 * @brief Defines the affinity complex titration binding model
		 */
		template <class ParamHandler_t>
		class AffinityComplexTitrationBindingBase : public ParamHandlerBindingModelBase<ParamHandler_t>
		{
		public:
			AffinityComplexTitrationBindingBase() : _useIonConc(false) {}

			virtual ~AffinityComplexTitrationBindingBase() CADET_NOEXCEPT {}

			static const char* identifier() { return ParamHandler_t::identifier(); }

			virtual bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
			{
				const bool res = BindingModelBase::configureModelDiscretization(paramProvider, nComp, nBound, boundOffset);

				// Guarantee that the first component is not a bound species
				if (nBound[0] != 0)
					throw InvalidParameterException("Affinity complex titration binding model requires the first component to be non-binding");

				for (int i = 0; i < nComp; ++i)
				{
					if (nBound[i] > 1)
						throw InvalidParameterException("Currently the ACT isotherm model supports at most one bound state per component");
				}

				// Read parameters related to ion concentration handling
				if (paramProvider.exists("ACT_USE_ION_CONC"))
					_useIonConc = paramProvider.getBool("ACT_USE_ION_CONC");
				else if (paramProvider.exists("EXT_ACT_USE_ION_CONC"))
					_useIonConc = paramProvider.getBool("EXT_ACT_USE_ION_CONC");
				else
					throw InvalidParameterException("ACT_USE_ION_CONC does not exist in model specification");

				// Read either pKa or cMid parameters depending on the value of _useIonConc
				if (!_useIonConc)
				{
					_cMidA.clear();
					_cMidG.clear();

					const std::vector<double> pKaA = getDoubleArrayWithFallback(paramProvider, "ACT_PKAA", "EXT_ACT_PKAA");
					const std::vector<double> pKaG = getDoubleArrayWithFallback(paramProvider, "ACT_PKAG", "EXT_ACT_PKAG");

					if (pKaA.empty() || pKaG.empty())
						throw InvalidParameterException("ACT_USE_ION_CONC=false requires ACT_PKAA and ACT_PKAG");
					if ((pKaA.size() < nComp) || (pKaG.size() < nComp))
						throw InvalidParameterException("ACT_PKAA and ACT_PKAG must have at least NCOMP entries");

					_pKaA.assign(pKaA.begin(), pKaA.begin() + nComp);
					_pKaG.assign(pKaG.begin(), pKaG.begin() + nComp);

				}
				else
				{
					_pKaA.clear();
					_pKaG.clear();

					const std::vector<double> cMidA = getDoubleArrayWithFallback(paramProvider, "ACT_CMID_A", "EXT_ACT_CMID_A");
					const std::vector<double> cMidG = getDoubleArrayWithFallback(paramProvider, "ACT_CMID_G", "EXT_ACT_CMID_G");
					if (cMidA.empty() || cMidG.empty())
						throw InvalidParameterException("ACT_USE_ION_CONC=true requires ACT_CMID_A and ACT_CMID_G");
					if ((cMidA.size() < nComp) || (cMidG.size() < nComp))
						throw InvalidParameterException("ACT_CMID_A and ACT_CMID_G must have at least NCOMP entries");

					_cMidA.assign(cMidA.begin(), cMidA.begin() + nComp);
					_cMidG.assign(cMidG.begin(), cMidG.begin() + nComp);

					for (unsigned int i = 0; i < nComp; ++i)
					{
						if ((cMidA[i] < 0.0) || (cMidG[i] < 0.0))
							throw InvalidParameterException("ACT_CMID_A and ACT_CMID_G entries must be >= 0 for ACT_USE_ION_CONC=true");
					}
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

			bool _useIonConc;
			std::vector<double> _pKaA;
			std::vector<double> _pKaG;
			std::vector<double> _cMidA;
			std::vector<double> _cMidG;

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
				const CpStateParamType IonAxis = !_useIonConc ? yCp[0] : -log(yCp[0] + static_cast<ParamType>(minIonConcentration)) / static_cast<ParamType>(ln10);

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

				// the first component is ion conc or negative log ion conc
				bndIdx = 0;
				for (unsigned int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states (bound state index bndIdx is not advanced)
					if (_nBoundStates[i] == 0)
						continue;

					// Residual
					const CpStateParamType pKaAaxis = !_useIonConc
						? static_cast<CpStateParamType>(_pKaA[i])
						: -log(static_cast<CpStateParamType>(_cMidA[i]) + static_cast<CpStateParamType>(minIonConcentration)) / static_cast<CpStateParamType>(ln10);

					const CpStateParamType pKaGaxis = !_useIonConc
						? static_cast<CpStateParamType>(_pKaG[i])
						: -log(static_cast<CpStateParamType>(_cMidG[i]) + static_cast<CpStateParamType>(minIonConcentration)) / static_cast<CpStateParamType>(ln10);

					const CpStateParamType etaAaxis = static_cast<CpStateParamType>(p->etaA[i]);
					const CpStateParamType etaGaxis = static_cast<CpStateParamType>(p->etaG[i]);

					const ResParamType f_A = stableActGate(etaAaxis, pKaAaxis, IonAxis);
					const ResParamType f_G = stableActGate(etaGaxis, pKaGaxis, IonAxis);

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
				const double IonAxis = !_useIonConc ? yCp[0] : -std::log(yCp[0] + minIonConcentration) / ln10;
				const double dIonAxis_dC0 = !_useIonConc ? 1.0 : -1.0 / ((yCp[0] + minIonConcentration) * ln10);

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

					// Local variables to aid the calculation of the Jacobian
					const double pKaAaxis = !_useIonConc ? _pKaA[i] : -std::log(_cMidA[i] + minIonConcentration) / ln10;
					const double pKaGaxis = !_useIonConc ? _pKaG[i] : -std::log(_cMidG[i] + minIonConcentration) / ln10;

					const double f_A = stableActGate(static_cast<double>(p->etaA[i]), pKaAaxis, IonAxis);
					const double f_G = stableActGate(static_cast<double>(p->etaG[i]), pKaGaxis, IonAxis);

					const double f_A_deriv = stableActGateDerivIon(static_cast<double>(p->etaA[i]), f_A);
					const double f_G_deriv = stableActGateDerivIon(static_cast<double>(p->etaG[i]), f_G);

					const double qApp = static_cast<double>(p->qMax[i]) * f_A;
					const double qFree_local = qApp - qsum;
					const double qScale = (qApp > minApparentCapacity) ? qApp : minApparentCapacity;
					const double beta = softplusTransitionFactor / qScale;
					const double bx = beta * qFree_local;
					const double qFree_eff = softplusScaled(qFree_local, beta);
					const double dqFreeEff_dqFree = sigmoid(bx);

					const double kA_i = static_cast<double>(p->kA[i]);
					const double kD_i = static_cast<double>(p->kD[i]);
					const double qMax_i = static_cast<double>(p->qMax[i]);

					const double kA_times_fG = kA_i * f_G;

					// dres_i / dc_{p,i}
					jac[i - bndIdx - offsetCp] = -kA_times_fG * qFree_eff;
					// Getting to c_{p,i}: -bndIdx takes us to q_0, another -offsetCp to c_{p,0} and a +i to c_{p,i}.
					//                     This means jac[i - bndIdx - offsetCp] corresponds to c_{p,i}.

					// dres_i / d(Ion axis), including the beta(qApp(IonAxis)) dependence
					const double dqApp_dIonAxis = qMax_i * f_A_deriv;

					double dsoftplus_dbeta = 0.0;
					double dbeta_dIonAxis = 0.0;

					if (qApp > minApparentCapacity)
					{
						const double logTerm = (bx > 0.0)
							? (bx + std::log1p(std::exp(-bx)))
							: std::log1p(std::exp(bx));

						dsoftplus_dbeta = (beta * qFree_local * dqFreeEff_dqFree - logTerm) / (beta * beta);
						dbeta_dIonAxis = -(beta / qApp) * dqApp_dIonAxis;
					}

					const double dqFreeEff_dIonAxis =
						dqFreeEff_dqFree * dqApp_dIonAxis + dsoftplus_dbeta * dbeta_dIonAxis;

					jac[-bndIdx - offsetCp] =
						-kA_i * yCp[i] * (f_G_deriv * qFree_eff + f_G * dqFreeEff_dIonAxis) * dIonAxis_dC0;

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
					jac[0] += kD_i;

					// Advance to next flux and Jacobian row
					++bndIdx;
					++jac;
				}
			}
		};

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