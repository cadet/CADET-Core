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
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "model/Parameters.hpp"
#include "LocalVector.hpp"
#include "SimulationTypes.hpp"
#include<cmath>
#include <vector>
#include <unordered_map>
#include <functional>
#include<iostream>
/*<codegen>
{
	"name": "EFunctionParamHandler",
	"externalName": "ExtEFunctionParamHandler",
	"parameters":
		[
			{ "type": "ScalarComponentDependentParameter", "varName": "eK", "confName": "EFUNC_K"},
			{ "type": "ScalarComponentDependentParameter", "varName": "eL", "confName": "EFUNC_L"},
			{ "type": "ScalarComponentDependentParameter", "varName": "eKkin", "confName": "EFUNC_KKIN"}
		]
}
</codegen>*/

/* Parameter description
 ------------------------
 eK	  : Some coefficient
 eL   : Saturation Limit
 eKkin : Kinetic coefficient
*/

namespace cadet
{

	namespace model
	{

		inline const char* EFunctionParamHandler::identifier() CADET_NOEXCEPT { return "EFUNCTION"; }

		inline bool EFunctionParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
		{
			if ((_eK.size() != _eL.size())  || (_eK.size() != _eKkin.size())  || (_eK.size() < nComp))
				throw InvalidParameterException("EFUNC_K, EFUNC_KKIN and EFUNC_L have to have the same size");

			return true;
		}

		inline const char* ExtEFunctionParamHandler::identifier() CADET_NOEXCEPT { return "EXT_EFUNCTION"; }

		inline bool ExtEFunctionParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
		{
			if ((_eK.size() != _eL.size()) || (_eK.size() != _eKkin.size()) || (_eK.size() < nComp))
				throw InvalidParameterException("EFUNC_K, EFUNC_KKIN and EFUNC_L have to have the same size");

			return true;
		}



		template <class ParamHandler_t>
		class EFunctionBindingBase : public ParamHandlerBindingModelBase<ParamHandler_t>
		{
		public:

			EFunctionBindingBase() { }
			virtual ~EFunctionBindingBase() CADET_NOEXCEPT { }

			static const char* identifier() { return ParamHandler_t::identifier(); }

			virtual bool implementsAnalyticJacobian() const CADET_NOEXCEPT { return true; }

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

				/* Protein Flux: dq/dt = k_kin * (q* – q) with
				q* = k_F * c^(1/n)
				*/
				unsigned int bndIdx = 0;
				const double threshold = 1e-14;

				for (int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states (bound state index bndIdx is not advanced)
					if (_nBoundStates[i] == 0)
						continue;
					const double eK_par = (double)static_cast<ParamType>(p->eK[i]);
					const double eL_par = (double)static_cast<ParamType>(p->eL[i]);
					const double ekkin_par = (double)static_cast<ParamType>(p->eKkin[i]);
					
					// Residual
					
					res[bndIdx] = ekkin_par * y[bndIdx] - ekkin_par * eL_par * (1 - std::exp(-(double)yCp[i]* eK_par) );

					// Next bound component
					++bndIdx;
				}

				return 0;
			}

			template <typename RowIterator>
			void jacobianImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* y, double const* yCp, int offsetCp, RowIterator jac, LinearBufferAllocator workSpace) const
			{
				typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

				const double threshold = 1e-14;

				int bndIdx = 0;
				for (int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states (bound state index bndIdx is not advanced)
					if (_nBoundStates[i] == 0)
						continue;
					// dres / dq_i
					const double eK_par = static_cast<double>(p->eK[i]);
					const double eL_par = static_cast<double>(p->eL[i]);
					const double ekkin_par = static_cast<double>(p->eKkin[i]);

					jac[0] = static_cast<double>(p->eKkin[i]);
					// dres / dc_{p,i}
										
					jac[i - bndIdx - offsetCp] = -ekkin_par * eL_par * eK_par *std::exp(-(double)yCp[i] * eK_par);


					// The distance from liquid phase to solid phase is reduced for each non-binding component
					// since a bound state is neglected. The number of neglected bound states so far is i - bndIdx.
					// Thus, by going back offsetCp - (i - bndIdx) = -[ i - bndIdx - offsetCp ] we get to the corresponding
					// liquid phase component.

					++bndIdx;
					++jac;
				}
			}

		};

		typedef EFunctionBindingBase<EFunctionParamHandler> EFunctionBinding;
		typedef EFunctionBindingBase<ExtEFunctionParamHandler> ExternalEFunctionBinding;

		namespace binding
		{
			void registerEFunctionModel(std::unordered_map<std::string, std::function<model::IBindingModel* ()>>& bindings)
			{
				bindings[EFunctionBinding::identifier()] = []() { return new EFunctionBinding(); };
				bindings[ExternalEFunctionBinding::identifier()] = []() { return new ExternalEFunctionBinding(); };
			}
		}  // namespace binding

	}  // namespace model

}  // namespace cadet
