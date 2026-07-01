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
#include <algorithm>

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

			if (_kKin.size() != nTotBnd)
				throw InvalidParameterException("SPLINE_KKIN has to have NTOTALBND entries");

			return true;
		}

		inline const char* ExtSplineParamHandler::identifier() CADET_NOEXCEPT { return "EXT_SPLINE_INTERPOLATION"; }

		inline bool ExtSplineParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
		{
			int nTotBnd = 0;
			for (int comp = 0; comp < nComp; comp++)
				nTotBnd += nBoundStates[comp];

			if (_kKin.size() != nTotBnd)
				throw InvalidParameterException("SPLINE_KKIN has to have NTOTALBND entries");

			return true;
		}


		template <class ParamHandler_t>
		class SplineBindingBase : public ParamHandlerBindingModelBase<ParamHandler_t>
		{
		public:

			SplineBindingBase() : _competitiveMode(false), _porePhaseConc(), _splineParams() {}

			virtual ~SplineBindingBase() CADET_NOEXCEPT { }

			static const char* identifier() { return ParamHandler_t::identifier(); }

			virtual bool implementsAnalyticJacobian() const CADET_NOEXCEPT { return true; }

			virtual bool requiresWorkspace() const CADET_NOEXCEPT { return true; }

			virtual unsigned int workspaceSize(unsigned int nComp, unsigned int totalNumBoundStates, unsigned int const* nBoundStates) const CADET_NOEXCEPT
			{
				return ParamHandlerBindingModelBase<ParamHandler_t>::workspaceSize(nComp, totalNumBoundStates, nBoundStates)
					+ sizeof(active) * (totalNumBoundStates + 4u * nComp)
					+ sizeof(size_t) * nComp + alignof(active); // qSpline plus competitive interpolation buffers
			}

			CADET_BINDINGMODELBASE_BOILERPLATE

		protected:
			using ParamHandlerBindingModelBase<ParamHandler_t>::_paramHandler;
			using ParamHandlerBindingModelBase<ParamHandler_t>::_reactionQuasistationarity;
			using ParamHandlerBindingModelBase<ParamHandler_t>::_nComp;
			using ParamHandlerBindingModelBase<ParamHandler_t>::_nBoundStates;
			int _totBoundStates;
			std::vector<int> _bndStateOffset;
			bool _competitiveMode;

			// Legacy one-dimensional spline storage q_i = f_i(c_i)
			std::vector< std::vector<double> > _porePhaseConc; // [_nComp][knots]
			std::vector< std::vector<double> > _splineParams; // [_totBoundStates][coeffs]
			// Structured-grid competitive spline storage q = f(c)
			tk::regular_grid_interpolator _competitiveSplineGrid;
			std::vector< std::vector<double> > _competitiveSplineValues; // [_totBoundStates][gridPoint]

			/***************************************************************************************************/
			size_t find_interval(double x, const std::vector<double>& m_x) const
			{
				const std::vector<double>::const_iterator it = std::upper_bound(m_x.begin(), m_x.end(), x);
				return std::max(int(it - m_x.begin()) - 1, 0);
			}

			bool almostEqual(double left, double right) const
			{
				const double scale = std::max(1.0, std::max(std::abs(left), std::abs(right)));
				return std::abs(left - right) <= 1e-12 * scale;
			}

			std::vector<double> uniqueSortedAxis(std::vector<double> values) const
			{
				std::sort(values.begin(), values.end());

				std::vector<double> axis;
				axis.reserve(values.size());

				for (std::vector<double>::const_iterator it = values.begin(); it != values.end(); ++it)
				{
					if (axis.empty() || !almostEqual(axis.back(), *it))
						axis.push_back(*it);
				}

				return axis;
			}

			size_t findAxisIndex(double value, const std::vector<double>& axis) const
			{
				const std::vector<double>::const_iterator it = std::lower_bound(axis.begin(), axis.end(), value);

				if (it != axis.end() && almostEqual(*it, value))
					return static_cast<size_t>(it - axis.begin());

				if (it != axis.begin() && almostEqual(*(it - 1), value))
					return static_cast<size_t>(it - axis.begin() - 1);

				throw InvalidParameterException("CP_VALS does not form a unique structured grid");
			}

			void configureCompetitiveSpline(IParameterProvider& paramProvider)
			{
				const std::vector<double> porePhaseConc = paramProvider.getDoubleArray("CP_VALS");
				const std::vector<double> solidPhaseConc = paramProvider.getDoubleArray("CS_VALS");

				if (porePhaseConc.empty() || solidPhaseConc.empty())
					throw InvalidParameterException("CP_VALS and CS_VALS must not be empty");

				if (porePhaseConc.size() % _nComp != 0)
					throw InvalidParameterException("CP_VALS must contain an integer number of pore-phase samples");

				if (solidPhaseConc.size() % _totBoundStates != 0)
					throw InvalidParameterException("CS_VALS must contain an integer number of bound-state samples");

				const size_t nSamples = porePhaseConc.size() / _nComp;
				if (solidPhaseConc.size() / _totBoundStates != nSamples)
					throw InvalidParameterException("CP_VALS and CS_VALS must contain the same number of samples");

				std::vector< std::vector<double> > axes(_nComp);
				for (int comp = 0; comp < _nComp; ++comp)
				{
					std::vector<double> compVals;
					compVals.reserve(nSamples);
					for (size_t sample = 0; sample < nSamples; ++sample)
						compVals.push_back(porePhaseConc[sample * _nComp + comp]);

					axes[comp] = uniqueSortedAxis(compVals);
					if (axes[comp].size() < 2)
						throw InvalidParameterException("Each spline input dimension requires at least two distinct CP_VALS entries");
				}

				size_t expectedSamples = 1;
				for (int comp = 0; comp < _nComp; ++comp)
					expectedSamples *= axes[comp].size();

				if (expectedSamples != nSamples)
					throw InvalidParameterException("CP_VALS does not cover a full Cartesian product grid");

				_competitiveSplineGrid.set_points(axes);
				_competitiveSplineValues.assign(_totBoundStates, std::vector<double>(_competitiveSplineGrid.num_points(), 0.0));

				std::vector<bool> seen(_competitiveSplineGrid.num_points(), false);
				for (size_t sample = 0; sample < nSamples; ++sample)
				{
					size_t flatIndex = 0;
					for (int comp = 0; comp < _nComp; ++comp)
					{
						const size_t axisIdx = findAxisIndex(porePhaseConc[sample * _nComp + comp], axes[comp]);
						flatIndex += axisIdx * _competitiveSplineGrid.stride(comp);
					}

					if (seen[flatIndex])
						throw InvalidParameterException("CP_VALS contains duplicate structured-grid points");
					seen[flatIndex] = true;

					for (int bndIdx = 0; bndIdx < _totBoundStates; ++bndIdx)
						_competitiveSplineValues[bndIdx][flatIndex] = solidPhaseConc[sample * _totBoundStates + bndIdx];
				}

				for (std::vector<bool>::const_iterator it = seen.begin(); it != seen.end(); ++it)
				{
					if (!(*it))
						throw InvalidParameterException("CP_VALS does not cover the full Cartesian grid");
				}
			}

			template <typename StateType>
			void competitiveSplineModel(StateType* q, StateType const* cp, LinearBufferAllocator workSpace) const
			{
				BufferedArray<size_t> lowerIndexArray = workSpace.array<size_t>(_nComp);
				size_t* const lowerIndices = static_cast<size_t*>(lowerIndexArray);
				BufferedArray<StateType> lowerWeightArray = workSpace.array<StateType>(_nComp);
				StateType* const lowerWeights = static_cast<StateType*>(lowerWeightArray);
				BufferedArray<StateType> upperWeightArray = workSpace.array<StateType>(_nComp);
				StateType* const upperWeights = static_cast<StateType*>(upperWeightArray);
				BufferedArray<StateType> invStepArray = workSpace.array<StateType>(_nComp);
				StateType* const invSteps = static_cast<StateType*>(invStepArray);

				_competitiveSplineGrid.compute_factors(cp, lowerIndices, lowerWeights, upperWeights, invSteps);

				for (int bndIdx = 0; bndIdx < _totBoundStates; ++bndIdx)
					_competitiveSplineGrid.evaluate(_competitiveSplineValues[bndIdx].data(), lowerIndices, lowerWeights, upperWeights, invSteps, q[bndIdx], static_cast<StateType*>(nullptr));
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

				_competitiveMode = false;

				// Input parameters
				_splineParams.resize(_totBoundStates);
				_porePhaseConc.resize(_nComp);
				_competitiveSplineValues.clear();

				if (paramProvider.exists("CP_VALS") || paramProvider.exists("CS_VALS"))
				{
					if (!(paramProvider.exists("CP_VALS") && paramProvider.exists("CS_VALS")))
						throw InvalidParameterException("Competitive spline mode requires both CP_VALS and CS_VALS");

					_competitiveMode = true;
					configureCompetitiveSpline(paramProvider);
					return result;
				}

				std::vector<double> solidPhaseConc;

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

				BufferedArray<CpStateType> qSplineArray = workSpace.array<CpStateType>(_totBoundStates);
				CpStateType* const qSpline = static_cast<CpStateType*>(qSplineArray);

				splineModel<CpStateType>(qSpline, yCp, workSpace);

				int bndIdx = 0;

				for (int comp = 0; comp < _nComp; ++comp)
				{
					for (int bnd = 0; bnd < _nBoundStates[comp]; ++bnd, ++bndIdx)
					{
						res[bndIdx] = static_cast<ParamType>(p->kKin[bndIdx]) * (y[bndIdx] - qSpline[bndIdx]);
					}
				}

				return 0;
			}

			template <typename StateType>
			void splineModel(StateType* q, StateType const* cp, LinearBufferAllocator workSpace) const
			{
				if (_competitiveMode)
				{
					competitiveSplineModel(q, cp, workSpace);
					return;
				}

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

				typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

				size_t* lowerIndices = nullptr;
				double* lowerWeights = nullptr;
				double* upperWeights = nullptr;
				double* invSteps = nullptr;
				double* dqDcpAll = nullptr;

				if (_competitiveMode)
				{
					BufferedArray<size_t> lowerIndexArray = workSpace.array<size_t>(_nComp);
					lowerIndices = static_cast<size_t*>(lowerIndexArray);
					BufferedArray<double> lowerWeightArray = workSpace.array<double>(_nComp);
					lowerWeights = static_cast<double*>(lowerWeightArray);
					BufferedArray<double> upperWeightArray = workSpace.array<double>(_nComp);
					upperWeights = static_cast<double*>(upperWeightArray);
					BufferedArray<double> invStepArray = workSpace.array<double>(_nComp);
					invSteps = static_cast<double*>(invStepArray);
					BufferedArray<double> derivativeArray = workSpace.array<double>(_nComp);
					dqDcpAll = static_cast<double*>(derivativeArray);

					_competitiveSplineGrid.compute_factors(yCp, lowerIndices, lowerWeights, upperWeights, invSteps);
				}

				int bndIdx = 0;

				for (int comp = 0; comp < _nComp; ++comp)
				{
					for (int bnd = 0; bnd < _nBoundStates[comp]; ++bnd, ++bndIdx)
					{
						jac[0] = static_cast<double>(p->kKin[bndIdx]);

						if (_competitiveMode)
						{
							double qValue = 0.0;
							_competitiveSplineGrid.evaluate(_competitiveSplineValues[bndIdx].data(), lowerIndices, lowerWeights, upperWeights, invSteps, qValue, dqDcpAll);

							for (int cpComp = 0; cpComp < _nComp; ++cpComp)
								jac[cpComp - bndIdx - offsetCp] = -static_cast<double>(p->kKin[bndIdx]) * dqDcpAll[cpComp];
						}
						else
						{
							const auto& coeffs = _splineParams[bndIdx];
							const int n_param = coeffs.size();
							const double cp_val = yCp[comp];
							const int n_pts = _porePhaseConc[comp].size();
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

							jac[comp - bndIdx - offsetCp] = -static_cast<double>(p->kKin[bndIdx]) * dq_dcp;
						}

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
