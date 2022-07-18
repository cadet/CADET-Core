// =============================================================================
//  CADET
//
//  Copyright © 2008-2021: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "model/binding/BindingModelBase.hpp"
#include "model/ExternalFunctionSupport.hpp"
#include "model/binding/BindingModelMacros.hpp"
#include "model/binding/RefConcentrationSupport.hpp"
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
	"name": "HICWANGParamHandler",
	"externalName": "ExtHICWANGParamHandler",
	"parameters":
		[
			{ "type": "ScalarParameter", "varName": "beta0", "confName": "BETA0"},
			{ "type": "ScalarParameter", "varName": "beta1", "confName": "BETA1"},
			{ "type": "ScalarComponentDependentParameter", "varName": "kEq", "confName": "KEQ"},
			{ "type": "ScalarComponentDependentParameter", "varName": "kKin", "confName": "KKIN"},
			{ "type": "ScalarComponentDependentParameter", "varName": "nu", "confName": "NU"},
			{ "type": "ScalarComponentDependentParameter", "varName": "qMax", "confName": "QMAX"}
		]
}
</codegen>*/

/* Parameter description
 ------------------------
 beta0 = bulk-like ordered water molecules at infinitely diluted salt concentration
 beta1 = influence of salt concentration on bulk-like ordered water molecules
 kEq = Equilibrium constant
 kKin = Kinetic constant
 nu = Number of binding sites
 qMax = Maximum binding capacity
*/

namespace cadet
{

	namespace model
	{

		inline const char* HICWANGParamHandler::identifier() CADET_NOEXCEPT { return "HICWANG"; }

		inline bool HICWANGParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
		{
			if ((_kEq.size() != _kKin.size()) || (_kEq.size() != _nu.size()) || (_kEq.size() != _qMax.size()) || (_kEq.size() < nComp))
				throw InvalidParameterException("KEQ, KKIN, NU, and QMAX have to have the same size");

			return true;
		}

		inline const char* ExtHICWANGParamHandler::identifier() CADET_NOEXCEPT { return "EXT_HICWANG"; }

		inline bool ExtHICWANGParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
		{
			if ((_kEq.size() != _kKin.size()) || (_kEq.size() != _nu.size()) || (_kEq.size() != _qMax.size()) || (_kEq.size() < nComp))
				throw InvalidParameterException("KEQ, KKIN, NU, and QMAX have to have the same size");

			return true;
		}


		/**
		 * @brief Defines the HIC Isotherm as described by Wang et al 2016
		 * @details Implements the "water on hydrophobic surfaces" model: \f[ \begin{align}
		 *				\beta&=\beta_0 e^{c_{p,0}\beta_1}\\
 		 *				\frac{\mathrm{d}q_i}{\mathrm{d}t} &= k_{eq,i} c_{p,i} \left( 1 - \sum_j \frac{q_j}{q_{max,j}} \right)^{\nu_i} - k_{kin,i} q_i^{1+\nu_i \beta}
 		 *			\end{align}  \f]
		 *          Component @c 0 is assumed to be salt without a bound state. Multiple bound states are not supported.
		 *          Components without bound state (i.e., salt and non-binding components) are supported.
		 *
		 *          See @cite Wang2016.
		 * @tparam ParamHandler_t Type that can add support for external function dependence
		 */
		template <class ParamHandler_t>
		class HICWangBase : public ParamHandlerBindingModelBase<ParamHandler_t>
		{
		public:

			HICWangBase() { }
			virtual ~HICWangBase() CADET_NOEXCEPT { }

			static const char* identifier() { return ParamHandler_t::identifier(); }

			virtual bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
			{
				const bool res = BindingModelBase::configureModelDiscretization(paramProvider, nComp, nBound, boundOffset);

				// Guarantee that salt has no bound state
				if (nBound[0] != 0)
					throw InvalidParameterException("HICWANG binding model requires exactly zero bound states for salt component");

				return res;
			}

			virtual bool hasSalt() const CADET_NOEXCEPT { return true; }
			virtual bool supportsMultistate() const CADET_NOEXCEPT { return false; }
			virtual bool supportsNonBinding() const CADET_NOEXCEPT { return true; }
			virtual bool hasQuasiStationaryReactions() const CADET_NOEXCEPT { return false; }
			virtual bool implementsAnalyticJacobian() const CADET_NOEXCEPT { return false; }

			virtual bool preConsistentInitialState(double t, unsigned int secIdx, const ColumnPosition& colPos, double* y, double const* yCp, LinearBufferAllocator workSpace) const
			{
				typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);
				return true;
			}

			virtual void postConsistentInitialState(double t, unsigned int secIdx, const ColumnPosition& colPos, double* y, double const* yCp, LinearBufferAllocator workSpace) const
			{
				preConsistentInitialState(t, secIdx, colPos, y, yCp, workSpace);
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
				using CpStateParamTypeParamType = typename DoubleActivePromoter<CpStateType, ParamType>::type;
				using StateParamType = typename DoubleActivePromoter<StateType, ParamType>::type;

				typename ParamHandler_t::ParamsHandle const _p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);


				auto beta0 = static_cast<ParamType>(_p->beta0);
				auto beta1 = static_cast<ParamType>(_p->beta1);

				auto beta = static_cast<ParamType>(beta0 * exp(beta1 * yCp[0]));

				ResidualType qSum = 1.0;
				ResidualType qProd = 0.0;
				unsigned int bndIdx = 0;
				for (int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states (bound state index bndIdx is not advanced)
					if (_nBoundStates[i] == 0)
						continue;

					auto q = static_cast<ParamType>(y[bndIdx]);
					qSum -= y[bndIdx] / static_cast<ParamType>(_p->qMax[i]);

					if (q < 0)
					{
						q = 0.0;
					}
					else
					{
						qProd += pow(q, static_cast<ParamType>(_p->nu[i]) * beta);
					}
					// Next bound component
					++bndIdx;
				}

				bndIdx = 0;

				for (int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states (bound state index bndIdx is not advanced)
					if (_nBoundStates[i] == 0)
						continue;
					auto qMax = static_cast<ParamType>(_p->qMax[i]);
					auto kKin = static_cast<ParamType>(_p->kKin[i]);
					auto kEq = static_cast<ParamType>(_p->kEq[i]);
					auto nu = static_cast<ParamType>(_p->nu[i]);

					if (y[bndIdx] < 0 || yCp[i] < 0)
					{
						res[bndIdx] = kKin * kEq * yCp[i] * (nu * y[bndIdx] / qMax - 1);
					}
					else
					{
						res[bndIdx] = kKin * (y[bndIdx] * static_cast<ParamType>(qProd) - kEq * static_cast<ParamType>(pow(qSum, nu)) * yCp[i]);
					}

					// Next bound component
					++bndIdx;
				}

				return 0;
			}

			template <typename RowIterator>
			void jacobianImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* y, double const* yCp, int offsetCp, RowIterator jac, LinearBufferAllocator workSpace) const
			{

			}
		};


		typedef HICWangBase<HICWANGParamHandler> HICWang;
		typedef HICWangBase<ExtHICWANGParamHandler> ExternalHICWang;

		namespace binding
		{
			void registerHICWangModel(std::unordered_map<std::string, std::function<model::IBindingModel* ()>>& bindings)
			{
				bindings[HICWang::identifier()] = []() { return new HICWang(); };
				bindings[ExternalHICWang::identifier()] = []() { return new ExternalHICWang(); };
			}
		}  // namespace binding

	}  // namespace model

}  // namespace cadet
