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

#include <functional>
#include <unordered_map>
#include <string>
#include <vector>

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
					throw InvalidParameterException("Affinity complex titration binding model requires the first component pH to be non-binding");

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

				// Protein fluxes
				ResidualType qSum = 0.0;
				unsigned int bndIdx = 0;
				for (int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states (bound state index bndIdx is not advanced)
					if (_nBoundStates[i] == 0)
						continue;

					// y is only defined for components that have a bound state
					qSum -= y[bndIdx] / static_cast<ParamType>(p->qMax[i]);

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
					const ResidualType qSum_local = qSum + 1.0 / (1.0 + pow(10.0, static_cast<ParamType>(p->etaA[i]) * (static_cast<ParamType>(p->pKaA[i]) - yCp[0])));
					const ResidualType kA_local = static_cast<ParamType>(p->kA[i]) / (1.0 + pow(10.0, static_cast<ParamType>(p->etaG[i]) * (static_cast<ParamType>(p->pKaG[i]) - yCp[0])));
					res[bndIdx] = static_cast<ParamType>(p->kD[i]) * y[bndIdx] - kA_local * yCp[i] * static_cast<ParamType>(p->qMax[i]) * qSum_local;

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
				for (int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states
					if (_nBoundStates[i] == 0)
						continue;

					qsum -= y[bndIdx] / static_cast<double>(p->qMax[i]);

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
					const double qMax_denominator = 1.0 + pow(10.0, static_cast<double>(p->etaA[i]) * (static_cast<double>(p->pKaA[i]) - yCp[0]));
					const double f_A = 1.0 / qMax_denominator;

					const double ka_denominator = 1.0 + pow(10.0, static_cast<double>(p->etaG[i]) * (static_cast<double>(p->pKaG[i]) - yCp[0]));
					const double f_G = 1.0 / ka_denominator;

					const double qSum_local = qsum + f_A;
					const double qmax_times_ka = static_cast<double>(p->kA[i]) * static_cast<double>(p->qMax[i]);

					const double f_A_deriv = std::log(10.0) * static_cast<double>(p->etaA[i]) * f_A * f_A * (qMax_denominator - 1.0);
					const double f_G_deriv = std::log(10.0) * static_cast<double>(p->etaG[i]) * f_G * f_G * (ka_denominator - 1.0);

					// dres_i / dc_{p,i}
					jac[i - bndIdx - offsetCp] = -qmax_times_ka * qSum_local * f_G;
					// Getting to c_{p,i}: -bndIdx takes us to q_0, another -offsetCp to c_{p,0} and a +i to c_{p,i}.
					//                     This means jac[i - bndIdx - offsetCp] corresponds to c_{p,i}.

					// dres_i / dpH, pH is the first c_{p,i}
					jac[-bndIdx - offsetCp] = -(f_A_deriv * f_G + f_A * f_G_deriv - f_G_deriv * qsum) * qmax_times_ka * yCp[i];

					// Fill dres_i / dq_j
					int bndIdx2 = 0;
					for (int j = 0; j < _nComp; ++j)
					{
						// Skip components without bound states (bound state index bndIdx2 is not advanced)
						if (_nBoundStates[j] == 0)
							continue;

						// dres_i / dq_j, part 1
						jac[bndIdx2 - bndIdx] = qmax_times_ka * f_G * yCp[i] / static_cast<double>(p->qMax[j]);
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