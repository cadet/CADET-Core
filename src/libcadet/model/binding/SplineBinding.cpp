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
#include "model/binding/BindingModelMacros.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "model/Parameters.hpp"
#include "LocalVector.hpp"
#include "SimulationTypes.hpp"
#include "AdUtils.hpp"
#include "model/binding/spline.h"
#include <functional>
#include <unordered_map>
#include <string>
#include <vector>
#include <cmath>

/*<codegen>
{
	"name": "SplineParamHandler",
	"externalName": "ExtSplineParamHandler",
	"parameters":
		[
			{ "type": "ScalarComponentDependentParameter", "varName": "kKin", "confName": "SPLINE_KKIN"}
		]
}
</codegen>*/

namespace cadet
{

	namespace model
	{

		inline const char* SplineParamHandler::identifier() CADET_NOEXCEPT { return "SPLINE_INTERPOLATION"; }

		inline bool SplineParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
		{
			int nTotBnd = 0;
			for (int comp = 0; comp < nComp; comp++)
				nTotBnd += nBoundStates[comp];

			if (_kKin.size() < nTotBnd)
				throw InvalidParameterException("SPLINE_KKIN has to have NTOTALNBND entries");

			return true;
		}

		inline const char* ExtSplineParamHandler::identifier() CADET_NOEXCEPT { return "EXT_SPLINE_INTERPOLATION"; }

		inline bool ExtSplineParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
		{
			int nTotBnd = 0;
			for (int comp = 0; comp < nComp; comp++)
				nTotBnd += nBoundStates[comp];

			if (_kKin.size() < nTotBnd)
				throw InvalidParameterException("SPLINE_KKIN has to have NTOTALNBND entries");

			return true;
		}


		template <class ParamHandler_t>
		class SplineBindingBase : public ParamHandlerBindingModelBase<ParamHandler_t>
		{
		public:

			SplineBindingBase() :_porePhaseConc(), _splineParams() {}

			virtual ~SplineBindingBase() CADET_NOEXCEPT { }

			static const char* identifier() { return ParamHandler_t::identifier(); }

			virtual bool implementsAnalyticJacobian() const CADET_NOEXCEPT { return true; }

			virtual unsigned int workspaceSize(unsigned int nComp, unsigned int totalNumBoundStates, unsigned int const* nBoundStates) const CADET_NOEXCEPT
			{
				return ParamHandlerBindingModelBase<ParamHandler_t>::workspaceSize(nComp, totalNumBoundStates, nBoundStates)
					+ sizeof(active) * nComp + alignof(active);  // Add the memory for the qML buffer
			}

			CADET_BINDINGMODELBASE_BOILERPLATE

		protected:
			using ParamHandlerBindingModelBase<ParamHandler_t>::_paramHandler;
			using ParamHandlerBindingModelBase<ParamHandler_t>::_reactionQuasistationarity;
			using ParamHandlerBindingModelBase<ParamHandler_t>::_nComp;
			using ParamHandlerBindingModelBase<ParamHandler_t>::_nBoundStates;
			int _totBoundStates;
			std::vector<int> _bndStateOffset;

			// Storage for trained ANN curves for spline fitting
			std::vector< std::vector<double> > _porePhaseConc; // [_nComp][coeffs]
			std::vector< std::vector<double> > _splineParams; // [_totBoundStates][coeffs]

			/***************************************************************************************************/
			size_t find_interval(double x, const std::vector<double>& m_x) const
			{
				const std::vector<double>::const_iterator it = std::upper_bound(m_x.begin(), m_x.end(), x);
				return std::max(int(it - m_x.begin()) - 1, 0);
			}

			virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
			{
				const bool result = ParamHandlerBindingModelBase<ParamHandler_t>::configureImpl(paramProvider, unitOpIdx, parTypeIdx);

				_totBoundStates = 0;
				_bndStateOffset.resize(_nComp + 1);
				_bndStateOffset[0] = 0;

				for (int comp = 0; comp < _nComp + 1; ++comp)
				{
					if (comp > 0)
						_bndStateOffset[comp] = _bndStateOffset[comp - 1] + _nBoundStates[comp - 1];
				}
				_totBoundStates = _bndStateOffset[_nComp];

				// Input parameters

				_splineParams.resize(_totBoundStates);
				_porePhaseConc.resize(_nComp);

				std::vector<double> solidPhaseConc; // Temporary storage for solid phase concentration values for current component and bound state

				for (int comp = 0; comp < _nComp; ++comp)
				{
					if (_nBoundStates[comp] == 0)
						continue;

					_porePhaseConc[comp] = paramProvider.getDoubleArray(std::format("CP_VALS_COMP_{:03}", comp));

					for (int bnd = 0; bnd < _nBoundStates[comp]; ++bnd)
					{
						if (_nBoundStates[comp] == 1 && paramProvider.exists(std::format("CS_VALS_COMP_{:03}", comp)))
							solidPhaseConc = paramProvider.getDoubleArray(std::format("CS_VALS_COMP_{:03}", comp));
						else
						{
							std::string inputName = std::format("CS_VALS_COMP_{:03}", comp);
							solidPhaseConc = paramProvider.getDoubleArray(inputName + std::format("_BND_{:03}", bnd));
						}

						if (solidPhaseConc.size() != _porePhaseConc[comp].size())
							throw InvalidParameterException("CS_VALS and CP_VALS must have the same number of entries");

						tk::spline s;
						s.set_boundary(tk::spline::second_deriv, 0.0,
							tk::spline::first_deriv, 0.0);
						s.set_points(_porePhaseConc[comp], solidPhaseConc, tk::spline::cspline);
						s.make_monotonic();
						_splineParams[_bndStateOffset[comp] + bnd] = s.coeff();
					}
				}

				return result;
			}

			template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
			int fluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, StateType const* y,
				CpStateType const* yCp, ResidualType* res, LinearBufferAllocator workSpace) const
			{
				// Implements -kKin * (f(c_p) - q) = kKin * (q - f(c_p)) where f is the spline model

				typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

				BufferedArray<CpStateType> qMLarray = workSpace.array<CpStateType>(_totBoundStates);
				CpStateType* const qML = static_cast<CpStateType*>(qMLarray);

				splineModel<CpStateType>(qML, yCp, workSpace);

				int bndIdx = 0;

				for (int comp = 0; comp < _nComp; ++comp)
				{
					for (int bnd = 0; bnd < _nBoundStates[comp]; ++bnd, ++bndIdx)
					{
						res[bndIdx] = static_cast<ParamType>(p->kKin[bndIdx]) * (y[bndIdx] - qML[bndIdx]);
					}
				}

				return 0;
			}

			template <typename StateType>
			void splineModel(StateType* q, StateType const* cp, LinearBufferAllocator workSpace) const
			{
				int bndIdx = 0;

				for (int comp = 0; comp < _nComp; ++comp)
				{
					for (int bnd = 0; bnd < _nBoundStates[comp]; ++bnd, ++bndIdx)
					{
						const auto& coeffs = _splineParams[bndIdx];
						const int n_param = coeffs.size();
						const int n_pts = _porePhaseConc[comp].size();

						const int idx = find_interval(static_cast<double>(cp[comp]), _porePhaseConc[comp]);
						const StateType h = cp[comp] - _porePhaseConc[comp][idx];

						if (cp[comp] < _porePhaseConc[comp][0])
						{
							q[bndIdx] = (coeffs[1] * h + coeffs[2]) * h + coeffs[3];
						}
						else if (cp[comp] >= _porePhaseConc[comp][n_pts - 1])
						{
							q[bndIdx] = (coeffs[n_param - 3] * h + coeffs[n_param - 2]) * h + coeffs[n_param - 1];
						}
						else
						{
							q[bndIdx] = ((coeffs[4 * idx] * h + coeffs[4 * idx + 1]) * h + coeffs[4 * idx + 2]) * h + coeffs[4 * idx + 3];
						}
					}
				}
			}

			template <typename RowIterator>
			void jacobianImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
				double const* y, double const* yCp, int offsetCp,
				RowIterator jac, LinearBufferAllocator workSpace) const
			{
				// Implements the Jacobian of -kKin * (f(c_p) - q) wrt. c_p and q where f is the spline model

				typename ParamHandler_t::ParamsHandle const p =
					_paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

				int bndIdx = 0;

				for (int comp = 0; comp < _nComp; ++comp)
				{
					const double cp_val = yCp[comp];
					const int n_pts = _porePhaseConc[comp].size();

					for (int bnd = 0; bnd < _nBoundStates[comp]; ++bnd, ++bndIdx)
					{
						const auto& coeffs = _splineParams[bndIdx];
						const int n_param = coeffs.size();

						const int idx = find_interval(cp_val, _porePhaseConc[comp]);
						const double h = cp_val - _porePhaseConc[comp][idx];

						double dq_dcp = 0.0;

						if (cp_val < _porePhaseConc[comp][0])
						{
							dq_dcp = 2.0 * coeffs[1] * h + coeffs[2];
						}
						else if (cp_val >= _porePhaseConc[comp][n_pts - 1])
						{
							dq_dcp = 2.0 * coeffs[n_param - 3] * h + coeffs[n_param - 2];
						}
						else
						{
							dq_dcp = (3.0 * coeffs[4 * idx] * h + 2.0 * coeffs[4 * idx + 1]) * h + coeffs[4 * idx + 2];
						}

						jac[0] = static_cast<double>(p->kKin[bndIdx]);
						jac[comp - bnd - _bndStateOffset[comp] - offsetCp] =
							-static_cast<double>(p->kKin[bndIdx]) * dq_dcp;

						++jac;
					}
				}
			}

		};

		typedef SplineBindingBase<SplineParamHandler> SplineBinding;
		typedef SplineBindingBase<ExtSplineParamHandler> ExternalSplineBinding;

		namespace binding
		{
			void registerSplineModel(std::unordered_map<std::string, std::function<model::IBindingModel* ()>>& bindings)
			{
				bindings[SplineBinding::identifier()] = []() { return new SplineBinding(); };
				bindings[ExternalSplineBinding::identifier()] = []() { return new ExternalSplineBinding(); };
			}
		}  // namespace binding

	}  // namespace model

}  // namespace cadet
