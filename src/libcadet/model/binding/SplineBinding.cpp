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
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "model/Parameters.hpp"
#include "LocalVector.hpp"
#include "SimulationTypes.hpp"
#include "AdUtils.hpp"
#include "model/binding/spline.h"
//#include <iostream>
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
			{ "type": "ScalarComponentDependentParameter", "varName": "kKin", "confName": "ML_KKIN"}
		]
}
</codegen>*/

namespace cadet
{

	namespace model
	{

		inline const char* SplineParamHandler::identifier() CADET_NOEXCEPT { return "SPLINE"; }

		inline bool SplineParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
		{
			if (_kKin.size() < nComp)
				throw InvalidParameterException("ML_KKIN has to have NCOMP entries");

			return true;
		}

		inline const char* ExtSplineParamHandler::identifier() CADET_NOEXCEPT { return "EXT_SPLINE"; }

		inline bool ExtSplineParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
		{
			if (_kKin.size() < nComp)
				throw InvalidParameterException("ML_KKIN has to have NCOMP entries");

			return true;
		}


		template <class ParamHandler_t>
		class SplineBindingBase : public ParamHandlerBindingModelBase<ParamHandler_t>
		{
		public:

			SplineBindingBase() :pore_phase_concentration(),
				Spline_parameters(),
				solid_phase_concentration() {}

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

			std::vector<double> pore_phase_concentration; //Defining these vectors to store trained ANN curve for spline fitting
			std::vector<double> Spline_parameters;  //Defining these vectors to store trained ANN curve for spline fitting
			std::vector<double> solid_phase_concentration;

			/***************************************************************************************************/
			size_t find_closest(double x, const std::vector<double>& m_x) const
			{
				const std::vector<double>::const_iterator it = std::upper_bound(m_x.begin(), m_x.end(), x);
				return std::max(int(it - m_x.begin()) - 1, 0);
			}

			virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
			{
				const bool result = ParamHandlerBindingModelBase<ParamHandler_t>::configureImpl(paramProvider, unitOpIdx, parTypeIdx);

				//Input parameters

				// Read some ML parameters
				paramProvider.pushScope("model_weights");
				paramProvider.pushScope("spline_input_parameters");

				solid_phase_concentration = paramProvider.getDoubleArray("Q_VALS");
				pore_phase_concentration = paramProvider.getDoubleArray("C_VALS");

				tk::spline s;
				s.set_boundary(tk::spline::second_deriv, 0.0,
					tk::spline::second_deriv, 0.0);
				s.set_points(pore_phase_concentration, solid_phase_concentration, tk::spline::cspline);
				s.make_monotonic();
				Spline_parameters = s.coeff();

				paramProvider.popScope(); // model_weights
				paramProvider.popScope(); // adsorption

				return result;
			}

			template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
			int fluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, StateType const* y,
				CpStateType const* yCp, ResidualType* res, LinearBufferAllocator workSpace) const
			{
				typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

				// We have to implement -kKin * (f(c_p) - q) = kKin * (q - f(c_p))
				// where f is the ML model

				// y points to q
				// yCp points to c_p
				// Use workSpace to obtain scratch memory

				// Get a vector of doubles for the result of the ML model
				BufferedArray<CpStateType> qMLarray = workSpace.array<CpStateType>(_nComp);
				CpStateType* const qML = static_cast<CpStateType*>(qMLarray); //double const* qML

				// Run the ML model on the c_p
				mlModel<CpStateType>(qML, yCp, workSpace);

				// Compute res = kKin * (q - f(c_p))
				unsigned int bndIdx = 0;
				for (int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states (bound state index bndIdx is not advanced)
					if (_nBoundStates[i] == 0)
						continue;

					// Residual: kKin * (q - f(c_p))
					res[bndIdx] = static_cast<ParamType>(p->kKin[i]) * (y[bndIdx] - qML[bndIdx]);

					// Next bound component
					++bndIdx;
				}

				return 0;
			}

			template <typename StateType>
			void mlModel(StateType* q, StateType const* cp, LinearBufferAllocator workSpace) const
			{
				// JAZIB: Here goes the ML model
				// Fill q array (output) using cp array (input)
				/* ***************************************************
				****Forward propagation of the neural network********
				******************************************************/
				//				const double factor = 1 / (1 - 0.69);
				const size_t n = pore_phase_concentration.size();
				
				const size_t n_param = Spline_parameters.size();
				
				const size_t idx = find_closest(static_cast<double>(cp[0]), pore_phase_concentration);

				const StateType h = cp[0] - pore_phase_concentration[idx];

				if (cp[0] < pore_phase_concentration[0]) {
					// extrapolation to the left
					q[0] = (Spline_parameters[1] * h + Spline_parameters[2]) * h + Spline_parameters[3];
				}
				else if (cp[0] >= pore_phase_concentration[n - 1]) {
					// extrapolatwdwyvd ion to the right
					q[0] = (Spline_parameters[n_param - 3] * h + Spline_parameters[n_param - 2]) * h + Spline_parameters[n_param - 1];
				}
				else {
					// interpolation
					q[0] = ((Spline_parameters[4 * idx] * h + Spline_parameters[4 * idx + 1]) * h + Spline_parameters[4 * idx + 2]) * h + Spline_parameters[4 * idx + 3];
					//First_derivative = (3.0 * layer_0_kernel0[4 * idx] * h + 2.0 * layer_0_kernel0[4 * idx + 1]) * h + layer_0_kernel0[4 * idx + 2];
				}
			}

			template <typename RowIterator>
			void jacobianImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* y, double const* yCp, int offsetCp, RowIterator jac, LinearBufferAllocator workSpace) const
			{
				typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);
				//std::cout << "Time is : " << t << "\n";
				// We have to implement Jacobian of -kKin * (f(c_p) - q) wrt. c_p and q
				// where f is the ML model

				// y points to q
				// yCp points to c_p
				// Use workSpace to obtain scratch memory

				unsigned int bndIdx = 0;
				for (int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states (bound state index bndIdx is not advanced)
					if (_nBoundStates[i] == 0)
						continue;

					const double kkin = static_cast<double>(p->kKin[i]);

					for (int j = 0; j < _nComp; ++j)
					{
						//						const double factor = 1 / (1 - 0.69);
						const double pore_val = yCp[0];

						const size_t n = pore_phase_concentration.size();

						const size_t idx = find_closest(pore_val, pore_phase_concentration);

						const double h = yCp[0] - pore_phase_concentration[idx];

						//						std::vector< std::vector<double> >dout(1, std::vector<double>(1));
						//						dout[0][0] = 0.0;

						const double First_derivative = (3.0 * Spline_parameters[4 * idx] * h + 2.0 * Spline_parameters[4 * idx + 1]) * h + Spline_parameters[4 * idx + 2];

						//						dout[0][0] =  factor * First_derivative; //*keq
												// dres_i / dc_{p,j} = -kkin[i] * df_i / dc_{p,j}
												//dout[0][0] = keq * 148.422 / (1 + keq * yCp[0] * yCp[0]);
						jac[j - bndIdx - offsetCp] = -kkin * First_derivative;

					}

					// dres_i / dq_i
					jac[0] = kkin;

					// Advance to next flux and Jacobian row
					++bndIdx;
					++jac;
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
