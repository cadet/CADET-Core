// SPDX-License-Identifier: AGPL-3.0-or-later
// =================================================================================
//  CADET
//
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Affero General Public
//  License v3.0 (or, at your option, any later version).
// =================================================================================

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
				// Buffers needed:
				//   fluxImpl:      qSpline[totBnd] + {lowerWeights, upperWeights, invSteps, cpBinding}[nComp] + lowerIndices[nComp] (size_t)
				//   jacobianImpl:  {lowerWeights, upperWeights, invSteps, dqDcp, cpBinding}[nComp] + lowerIndices[nComp] (size_t)
				// Use nComp as upper bound for nBindingComponents.
				return ParamHandlerBindingModelBase<ParamHandler_t>::workspaceSize(nComp, totalNumBoundStates, nBoundStates)
					+ sizeof(active) * (totalNumBoundStates + 5u * nComp)
					+ sizeof(size_t) * nComp + alignof(active);
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

			// Independent mode: per-component 1D spline storage, q_i = f_i(c_{p,i}).
			std::vector< std::vector<double> > _porePhaseConc; // [_nComp][knots]
			std::vector< std::vector<double> > _splineParams;  // [_totBoundStates][4*(nKnots-1) coeffs]
			// Competitive mode: multilinear interpolation on a Cartesian-product grid
			std::vector<int> _bindingComponents; // [_nComp] indices of components with nonzero bound states
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

			// Return a sorted copy of values with near-duplicate entries removed.
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

			// Return the index of value in axis, tolerating small floating-point errors.
			// Throws InvalidParameterException if value is not found.
			size_t findAxisIndex(double value, const std::vector<double>& axis) const
			{
				const std::vector<double>::const_iterator it = std::lower_bound(axis.begin(), axis.end(), value);

				if (it != axis.end() && almostEqual(*it, value))
					return static_cast<size_t>(it - axis.begin());

				if (it != axis.begin() && almostEqual(*(it - 1), value))
					return static_cast<size_t>(it - axis.begin() - 1);

				throw InvalidParameterException("CP_VALS_COMP does not form a unique structured grid");
			}

			void configureCompetitiveSpline(IParameterProvider& paramProvider)
			{
				const int nBindingComp = static_cast<int>(_bindingComponents.size());

				// Per-binding-component pore-phase concentration samples and grid axes.
				// Indexed 0..nBindingComp-1, not by original component index.
				std::vector< std::vector<double> > porePhaseConc(nBindingComp);
				std::vector< std::vector<double> > axes(nBindingComp);

				size_t nSamples = 0;
				for (int bindingIdx = 0; bindingIdx < nBindingComp; ++bindingIdx)
				{
					const int comp = _bindingComponents[bindingIdx];
					porePhaseConc[bindingIdx] = paramProvider.getDoubleArray(std::format("CP_VALS_COMP_{:03}", comp));

					if (porePhaseConc[bindingIdx].empty())
						throw InvalidParameterException("CP_VALS_COMP entries must not be empty");

					if (bindingIdx == 0)
						nSamples = porePhaseConc[bindingIdx].size();
					else if (porePhaseConc[bindingIdx].size() != nSamples)
						throw InvalidParameterException("All CP_VALS_COMP entries must contain the same number of samples");

					axes[bindingIdx] = uniqueSortedAxis(porePhaseConc[bindingIdx]);
					if (axes[bindingIdx].size() < 2)
						throw InvalidParameterException("Each spline input dimension requires at least two distinct CP_VALS_COMP entries");
				}

				size_t expectedSamples = 1;
				for (int bindingIdx = 0; bindingIdx < nBindingComp; ++bindingIdx)
					expectedSamples *= axes[bindingIdx].size();

				if (expectedSamples != nSamples)
					throw InvalidParameterException("CP_VALS_COMP does not cover a full Cartesian product grid");

				// Grid dimensions correspond to _bindingComponents[0..nBindingComp-1].
				_competitiveSplineGrid.set_points(axes);
				_competitiveSplineValues.assign(_totBoundStates, std::vector<double>(_competitiveSplineGrid.num_points(), 0.0));

				// Map each sample to its flat grid index.
				std::vector<size_t> sampleToFlatIndex(nSamples);
				std::vector<bool> seen(_competitiveSplineGrid.num_points(), false);
				for (size_t sample = 0; sample < nSamples; ++sample)
				{
					size_t flatIndex = 0;
					for (int bindingIdx = 0; bindingIdx < nBindingComp; ++bindingIdx)
					{
						const size_t axisIdx = findAxisIndex(porePhaseConc[bindingIdx][sample], axes[bindingIdx]);
						flatIndex += axisIdx * _competitiveSplineGrid.stride(bindingIdx);
					}

					if (seen[flatIndex])
						throw InvalidParameterException("CP_VALS_COMP contains duplicate structured-grid points");
					seen[flatIndex] = true;
					sampleToFlatIndex[sample] = flatIndex;
				}

				// Read solid-phase values; skip components with zero bound states.
				for (int comp = 0; comp < _nComp; ++comp)
				{
					for (int bnd = 0; bnd < _nBoundStates[comp]; ++bnd)
					{
						std::vector<double> solidPhaseConc;

						if (_nBoundStates[comp] == 1 && paramProvider.exists(std::format("CS_VALS_COMP_{:03}", comp)))
							solidPhaseConc = paramProvider.getDoubleArray(std::format("CS_VALS_COMP_{:03}", comp));
						else
						{
							const std::string inputName = std::format("CS_VALS_COMP_{:03}", comp);
							solidPhaseConc = paramProvider.getDoubleArray(inputName + std::format("_BND_{:03}", bnd));
						}

						if (solidPhaseConc.size() != nSamples)
							throw InvalidParameterException("CS_VALS_COMP and CP_VALS_COMP must have the same number of samples");

						const int bndIdx = _bndStateOffset[comp] + bnd;
						for (size_t sample = 0; sample < nSamples; ++sample)
							_competitiveSplineValues[bndIdx][sampleToFlatIndex[sample]] = solidPhaseConc[sample];
					}
				}

				for (std::vector<bool>::const_iterator it = seen.begin(); it != seen.end(); ++it)
				{
					if (!(*it))
						throw InvalidParameterException("CP_VALS_COMP does not cover the full Cartesian grid");
				}
			}

			void configureIndependentSpline(IParameterProvider& paramProvider)
			{
				// Configure independent mode: build one 1D cubic spline per bound state.
				// Boundary conditions: zero second derivative (left) → linear extrapolation;
				// zero first derivative (right) → flat extrapolation.
				std::vector<double> solidPhaseConc;

				for (int comp = 0; comp < _nComp; ++comp)
				{
					if (_nBoundStates[comp] == 0)
						continue;

					_porePhaseConc[comp] = paramProvider.getDoubleArray(std::format("CP_VALS_COMP_{:03}", comp));

					for (int bnd = 0; bnd < _nBoundStates[comp]; ++bnd)
					{
						std::string inputName = std::format("CS_VALS_COMP_{:03}", comp);
						solidPhaseConc = paramProvider.getDoubleArray(inputName + std::format("_BND_{:03}", bnd));

						if (solidPhaseConc.size() != _porePhaseConc[comp].size())
							throw InvalidParameterException("CS_VALS_COMP and CP_VALS_COMP must have the same number of entries");

						tk::spline s;
						s.set_boundary(tk::spline::second_deriv, 0.0,
							tk::spline::first_deriv, 0.0);
						s.set_points(_porePhaseConc[comp], solidPhaseConc, tk::spline::cspline);
						s.make_monotonic();
						_splineParams[_bndStateOffset[comp] + bnd] = s.coeff();
					}
				}
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

				// Collect components that participate in binding (used by competitive mode).
				_bindingComponents.clear();
				for (int comp = 0; comp < _nComp; ++comp)
				{
					if (_nBoundStates[comp] > 0)
						_bindingComponents.push_back(comp);
				}

				_splineParams.resize(_totBoundStates);
				_porePhaseConc.resize(_nComp);
				_competitiveSplineValues.clear();

				const std::string mode = paramProvider.getString("INTERPOLATION_MODE");

				if (mode == "COMPETITIVE_REGULAR_GRID")
				{
					_competitiveMode = true;
					configureCompetitiveSpline(paramProvider);
				}
				else if (mode == "INDEPENDENT")
				{
					_competitiveMode = false;
					configureIndependentSpline(paramProvider);
				}
				else
				{
					throw InvalidParameterException("INTERPOLATION_MODE must be INDEPENDENT or COMPETITIVE_REGULAR_GRID, got: " + mode);
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

			// Dispatch equilibrium evaluation to competitive (multilinear N-D) or
			// independent (cubic 1D) mode. Components with zero bound states are ignored
			// in both cases. In competitive mode, only binding components form the grid
			// dimensions; their cp values are extracted into a temporary buffer.
			template <typename StateType>
			void splineModel(StateType* q, StateType const* cp, LinearBufferAllocator workSpace) const
			{
				if (_competitiveMode)
				{
					const int nBindingComp = static_cast<int>(_bindingComponents.size());

					BufferedArray<size_t> lowerIndexArray = workSpace.array<size_t>(nBindingComp);
					size_t* const lowerIndices = static_cast<size_t*>(lowerIndexArray);
					BufferedArray<StateType> lowerWeightArray = workSpace.array<StateType>(nBindingComp);
					StateType* const lowerWeights = static_cast<StateType*>(lowerWeightArray);
					BufferedArray<StateType> upperWeightArray = workSpace.array<StateType>(nBindingComp);
					StateType* const upperWeights = static_cast<StateType*>(upperWeightArray);
					BufferedArray<StateType> invStepArray = workSpace.array<StateType>(nBindingComp);
					StateType* const invSteps = static_cast<StateType*>(invStepArray);

					// Extract pore-phase concentrations of binding components only.
					BufferedArray<StateType> cpBindingArray = workSpace.array<StateType>(nBindingComp);
					StateType* const cpBinding = static_cast<StateType*>(cpBindingArray);
					for (int bindingIdx = 0; bindingIdx < nBindingComp; ++bindingIdx)
						cpBinding[bindingIdx] = cp[_bindingComponents[bindingIdx]];

					// Locate the enclosing cell and compute interpolation weights once for all bound states.
					_competitiveSplineGrid.compute_factors(cpBinding, lowerIndices, lowerWeights, upperWeights, invSteps);

					for (int bndIdx = 0; bndIdx < _totBoundStates; ++bndIdx)
						_competitiveSplineGrid.evaluate(_competitiveSplineValues[bndIdx].data(), lowerIndices, lowerWeights, upperWeights, invSteps, q[bndIdx], static_cast<StateType*>(nullptr));

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
							// Left extrapolation: linear (zero-curvature boundary)
							q[bndIdx] = (coeffs[1] * h + coeffs[2]) * h + coeffs[3];
						}
						else if (cp[comp] >= _porePhaseConc[comp][n_pts - 1])
						{
							// Right extrapolation: flat (zero-slope boundary)
							q[bndIdx] = (coeffs[n_param - 3] * h + coeffs[n_param - 2]) * h + coeffs[n_param - 1];
						}
						else
						{
							// Interior: Horner evaluation of the cubic polynomial for interval idx
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
					const int nBindingComp = static_cast<int>(_bindingComponents.size());

					BufferedArray<size_t> lowerIndexArray = workSpace.array<size_t>(nBindingComp);
					lowerIndices = static_cast<size_t*>(lowerIndexArray);
					BufferedArray<double> lowerWeightArray = workSpace.array<double>(nBindingComp);
					lowerWeights = static_cast<double*>(lowerWeightArray);
					BufferedArray<double> upperWeightArray = workSpace.array<double>(nBindingComp);
					upperWeights = static_cast<double*>(upperWeightArray);
					BufferedArray<double> invStepArray = workSpace.array<double>(nBindingComp);
					invSteps = static_cast<double*>(invStepArray);
					BufferedArray<double> derivativeArray = workSpace.array<double>(nBindingComp);
					dqDcpAll = static_cast<double*>(derivativeArray);

					// Extract pore-phase concentrations of binding components only.
					BufferedArray<double> cpBindingArray = workSpace.array<double>(nBindingComp);
					double* const cpBinding = static_cast<double*>(cpBindingArray);
					for (int bindingIdx = 0; bindingIdx < nBindingComp; ++bindingIdx)
						cpBinding[bindingIdx] = yCp[_bindingComponents[bindingIdx]];

					_competitiveSplineGrid.compute_factors(cpBinding, lowerIndices, lowerWeights, upperWeights, invSteps);
				}

				int bndIdx = 0;

				for (int comp = 0; comp < _nComp; ++comp)
				{
					for (int bnd = 0; bnd < _nBoundStates[comp]; ++bnd, ++bndIdx)
					{
						// d(res_bndIdx)/d(q_bndIdx) = kKin
						jac[0] = static_cast<double>(p->kKin[bndIdx]);

						if (_competitiveMode)
						{
							const int nBindingComp = static_cast<int>(_bindingComponents.size());

							double qValue = 0.0;
							_competitiveSplineGrid.evaluate(_competitiveSplineValues[bndIdx].data(), lowerIndices, lowerWeights, upperWeights, invSteps, qValue, dqDcpAll);

							// Zero all cp Jacobian entries; only binding components have non-zero derivatives.
							for (int cpComp = 0; cpComp < _nComp; ++cpComp)
								jac[cpComp - bndIdx - offsetCp] = 0.0;

							// dqDcpAll[bindingIdx] = dq/d(cp[_bindingComponents[bindingIdx]])
							for (int bindingIdx = 0; bindingIdx < nBindingComp; ++bindingIdx)
							{
								const int cpComp = _bindingComponents[bindingIdx];
								jac[cpComp - bndIdx - offsetCp] = -static_cast<double>(p->kKin[bndIdx]) * dqDcpAll[bindingIdx];
							}
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
