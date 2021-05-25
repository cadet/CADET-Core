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
#include "C:\Users\hassan\Desktop\Spline Implementation\spline-master\src\spline.h"
#include<iostream>
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

			SplineBindingBase() : pore_phase_values(number_of_cp_points), Spline_parameters(trained_parameters){}

			virtual ~SplineBindingBase() CADET_NOEXCEPT { }

			static const char* identifier() { return ParamHandler_t::identifier(); }

			virtual bool implementsAnalyticJacobian() const CADET_NOEXCEPT { return true; }

			virtual unsigned int workspaceSize(unsigned int nComp, unsigned int totalNumBoundStates, unsigned int const* nBoundStates) const CADET_NOEXCEPT
			{
				return _paramHandler.cacheSize(nComp, totalNumBoundStates, nBoundStates)
					+ sizeof(double) * nComp + alignof(double) // Add the memory for the qML buffer
					+ sizeof(double) * nComp + alignof(double); // Add the memory for the c_p buffer	
					//+ sizeof(double) * number_of_input * number_of_nodes_W1 + alignof(double) 		// Add the memory for the Affine_layer_1 buffer	
					//+ sizeof(double) * number_of_input * number_of_nodes_W1 + alignof(double)    // Add the memory for the Affine_layer_1_Jacobian buffer
					//+ sizeof(double) * number_of_input * number_of_nodes_W2 + alignof(double)    // Add the memory for the Affine_layer_2 buffer
					//+ sizeof(double) * number_of_input * number_of_nodes_W2 + alignof(double)    // Add the memory for the Affine_layer_2_Jacobian buffer
					//+ sizeof(double) * number_of_input * number_of_outputs + alignof(double);    // Add the memory for the Affine_layer_3 buffer
			}

			CADET_BINDINGMODELBASE_BOILERPLATE

		protected:
			using ParamHandlerBindingModelBase<ParamHandler_t>::_paramHandler;
			using ParamHandlerBindingModelBase<ParamHandler_t>::_reactionQuasistationarity;
			using ParamHandlerBindingModelBase<ParamHandler_t>::_nComp;
			using ParamHandlerBindingModelBase<ParamHandler_t>::_nBoundStates;

			const int number_of_cp_points = 100;
			const int trained_parameters = (number_of_cp_points - 1) * 4;

			std::vector<double> pore_phase_values; //Defining these vectors to store trained ANN curve for spline fitting
			std::vector<double> Spline_parameters;  //Defining these vectors to store trained ANN curve for spline fitting


			/***************************************************************************************************/
			size_t find_closest(double x, std::vector<double>m_x) const
			{
				size_t idx;
				std::vector<double>::const_iterator it;
				it = std::upper_bound(m_x.begin(), m_x.end(), x);       // *it > x
				if (int(it - m_x.begin()) - 1 > 0)
				{
					idx = int(it - m_x.begin()) - 1;
				}
				else
				{
					idx = 0;
				}
				//size_t idx = int(it - m_x.begin()) - 1;
				//size_t idx = std::max(int(it - m_x.begin()) - 1, 0);   // m_x[idx] <= x
				return idx;
			}

			/*Question: One of the */
			virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
			{
				// Read parameters
				_paramHandler.configure(paramProvider, _nComp, _nBoundStates);

				// Register parameters
				_paramHandler.registerParameters(_parameters, unitOpIdx, parTypeIdx, _nComp, _nBoundStates);

				//Input parameters



				// Read some ML parameters
				paramProvider.pushScope("model_weights");
				paramProvider.pushScope("Spline_Input_Parameters");

				pore_phase_values = paramProvider.getDoubleArray("Input");
				Spline_parameters = paramProvider.getDoubleArray("Parameters");

				paramProvider.popScope(); // model_weights
				paramProvider.popScope(); // adsorption

				/**
				 * JAZIB: Read ML parameters and model data here (Done)
				 */

				return true;
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
				BufferedArray<double> qMLarray = workSpace.array<double>(_nComp);
				double* qML = static_cast<double*>(qMLarray); //double const* qML

				// Run the ML model on the c_p
				mlModel(qML, yCp, workSpace);

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

			void mlModel(double* q, double const* cp, LinearBufferAllocator workSpace) const
			{
				// JAZIB: Here goes the ML model
				// Fill q array (output) using cp array (input)
				/* ***************************************************
				****Forward propagation of the neural network********
				******************************************************/
				const double factor = 1 / (1 - 0.69);
				double pore_val = cp[0];
				size_t n = pore_phase_values.size();
				size_t n_param = Spline_parameters.size();

				size_t idx = find_closest(pore_val, pore_phase_values);

				double h = cp[0] - pore_phase_values[idx];
				double interpol;

				if (cp[0] < pore_phase_values[0]) {
					// extrapolation to the left
					interpol = (Spline_parameters[1] * h + Spline_parameters[2]) * h + Spline_parameters[3];
				}
				else if (cp[0] >= pore_phase_values[n - 1]) {
					// extrapolation to the right
					interpol = (Spline_parameters[n_param - 3] * h + Spline_parameters[n_param - 2]) * h + Spline_parameters[n_param - 1];
				}
				else {
					// interpolation

					interpol = ((Spline_parameters[4 * idx] * h + Spline_parameters[4 * idx + 1]) * h + Spline_parameters[4 * idx + 2]) * h + Spline_parameters[4 * idx + 3];
					//First_derivative = (3.0 * layer_0_kernel0[4 * idx] * h + 2.0 * layer_0_kernel0[4 * idx + 1]) * h + layer_0_kernel0[4 * idx + 2];
				}

				q[0] = interpol *factor;
				

			}

			void mlModel(double* q, active const* cp, LinearBufferAllocator workSpace) const
			{
				// Get a vector of doubles for c_p
				BufferedArray<double> cParray = workSpace.array<double>(_nComp);
				double* const cpDbl = static_cast<double*>(cParray);

				// Copy values of active to double (ignore / remove AD gradients)
				cadet::ad::copyFromAd(cp, cpDbl, _nComp);

				mlModel(q, cpDbl, workSpace);
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
						const double factor = 1 / (1 - 0.69);
						double First_derivative;
						double pore_val = yCp[0];

						size_t n = pore_phase_values.size();


						size_t idx = find_closest(pore_val, pore_phase_values);

						double h = yCp[0] - pore_phase_values[idx];

						std::vector< std::vector<double> >dout(1, std::vector<double>(1));
						dout[0][0] = 0.0;


						First_derivative = (3.0 * Spline_parameters[4 * idx] * h + 2.0 * Spline_parameters[4 * idx + 1]) * h + Spline_parameters[4 * idx + 2];
						/*tk::spline s;
						s.set_boundary(tk::spline::second_deriv, 0.0,
							tk::spline::second_deriv, 0.0);
						s.set_points(ANN_cp, ANN_q, tk::spline::cspline);
						s.make_monotonic();

						dout[0][0] = s.deriv(1, yCp[0]);*/
				
						dout[0][0] =  factor * First_derivative; //*keq
						// dres_i / dc_{p,j} = -kkin[i] * df_i / dc_{p,j}
						//dout[0][0] = keq * 148.422 / (1 + keq * yCp[0] * yCp[0]);
						jac[j - bndIdx - offsetCp] = -1 * (dout[0][0]);
						//std::cout << "Time " << t << " Jac is : " << jac[j - bndIdx - offsetCp] << "\n";
					}

					// dres_i / dq_i
					jac[0] = 1;

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
