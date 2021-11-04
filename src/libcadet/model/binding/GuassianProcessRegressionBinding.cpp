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
#include "model/binding/GPR_Class.h"
#include <functional>
#include <unordered_map>
#include <string>
#include <vector>
#include <cmath>

/*<codegen>
{
	"name": "GPRParamHandler",
	"externalName": "ExtGPRParamHandler",
	"parameters":
		[
			{ "type": "ScalarComponentDependentParameter", "varName": "kKin", "confName": "GPR_KKIN"}
		]
}
</codegen>*/

namespace cadet
{

	namespace model
	{

		inline const char* GPRParamHandler::identifier() CADET_NOEXCEPT { return "GPR"; }

		inline bool GPRParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
		{
			if (_kKin.size() < nComp)
				throw InvalidParameterException("ML_KKIN has to have NCOMP entries");

			return true;
		}

		inline const char* ExtGPRParamHandler::identifier() CADET_NOEXCEPT { return "EXT_GPR"; }

		inline bool ExtGPRParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
		{
			if (_kKin.size() < nComp)
				throw InvalidParameterException("ML_KKIN has to have NCOMP entries");

			return true;
		}


		template <class ParamHandler_t>
		class GPRBindingBase : public ParamHandlerBindingModelBase<ParamHandler_t>
		{
		public:

			GPRBindingBase() :pore_phase_concentration(),
				GPR_parameters(),
				solid_phase_concentration(){}

			virtual ~GPRBindingBase() CADET_NOEXCEPT { }

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

			std::vector<double> pore_phase_concentration; //Defining these vectors to store trained ANN curve for GPR fitting
			std::vector<double> GPR_parameters;  //Defining these vectors to store trained ANN curve for GPR fitting
			std::vector<double> solid_phase_concentration;
			std::vector<double> Kernel_Mat;
			std::vector<double> K_inv_y;
			double* offset;
			std::string kernel_name ;
			unsigned int ndim;
			/***************************************************************************************************/

			virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
			{
				const bool result = ParamHandlerBindingModelBase<ParamHandler_t>::configureImpl(paramProvider, unitOpIdx, parTypeIdx);

				//Input parameters
				
				// Read some ML parameters
				paramProvider.pushScope("training_data");
				

				solid_phase_concentration = paramProvider.getDoubleArray("Q_VALS");
				pore_phase_concentration = paramProvider.getDoubleArray("C_VALS");
				GPR_parameters = paramProvider.getDoubleArray("TRAINED_PARAMS");
				int kernel_name1 = paramProvider.getInt("KERNAL");
				ndim = paramProvider.getInt("NDIM");

				paramProvider.popScope(); // training_data

				std::vector<double> kernel_m(solid_phase_concentration.size() * solid_phase_concentration.size(), 0.0);
				std::vector<double> K_inv_y_temp(solid_phase_concentration.size(), 0.0);
				if (kernel_name1 == 1)
					kernel_name = "MLP";

				std::vector<double> test_pt = { 0.0 };

				GP::GPR_Class gp(solid_phase_concentration.size(), solid_phase_concentration.size(), ndim,
					GPR_parameters[0], GPR_parameters[1], GPR_parameters[2], GPR_parameters[3], GPR_parameters[4],
					GPR_parameters[5], GPR_parameters[6],
					kernel_name);
				if (kernel_name == "MLP")
				{
					gp.MLP_kernel(pore_phase_concentration.size(), pore_phase_concentration.size(), ndim,
						pore_phase_concentration.data(), pore_phase_concentration.data(), kernel_m.data());
				}

				gp.kernel_inv_y(pore_phase_concentration.data(), solid_phase_concentration.data(), kernel_m.data(), K_inv_y_temp.data());

				offset = gp.prediction(pore_phase_concentration.data(), solid_phase_concentration.data(),
					test_pt.data(), kernel_m.data(), K_inv_y_temp.data(), kernel_name);
				
				K_inv_y = K_inv_y_temp;
				Kernel_Mat = kernel_m;

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

			template <typename StateType>
			void mlModel(double* q, StateType const* cp, LinearBufferAllocator workSpace) const
			{
				// JAZIB: Here goes the ML model
				// Fill q array (output) using cp array (input)
				/* ***************************************************
				****Forward propagation of the neural network********
				******************************************************/
				const double* prediction;
				
				std::vector<double> matCp(_nComp, 0.0);

				for (unsigned int i = 0; i < _nComp; i++)
				{
					matCp[i]  = static_cast<double>(cp[0]);
				}
				GP::GPR_Class gp(solid_phase_concentration.size(), solid_phase_concentration.size(), ndim,
					GPR_parameters[0], GPR_parameters[1], GPR_parameters[2], GPR_parameters[3], GPR_parameters[4],
					GPR_parameters[5], GPR_parameters[6],
					kernel_name);

								
				prediction = gp.prediction(pore_phase_concentration.data(), solid_phase_concentration.data(),
					matCp.data(), Kernel_Mat.data(), K_inv_y.data(), kernel_name);


				for (unsigned int i = 0; i < _nComp; i++)
				{
					q[i] = prediction[i] - offset[0];
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

				
				std::vector<double> jacobian(_nComp, 0.0);
				std::vector<double> matCp(_nComp, 0.0);

				for (unsigned int i = 0; i < _nComp; i++)
				{
					matCp[i] = yCp[i];
				}

				GP::GPR_Class gp(solid_phase_concentration.size(), solid_phase_concentration.size(), ndim,
					GPR_parameters[0], GPR_parameters[1], GPR_parameters[2], GPR_parameters[3], GPR_parameters[4],
					GPR_parameters[5], GPR_parameters[6],
					kernel_name);

				gp.MLP_derivative(pore_phase_concentration.data(), matCp.data(), K_inv_y.data(), jacobian.data());

				unsigned int bndIdx = 0;
				for (int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states (bound state index bndIdx is not advanced)
					if (_nBoundStates[i] == 0)
						continue;

					const double kkin = static_cast<double>(p->kKin[i]);
					

					for (int j = 0; j < _nComp; ++j)
					{
						
						jac[j - bndIdx - offsetCp] = -kkin * jacobian[j];

					}

					// dres_i / dq_i
					jac[0] = kkin;

					// Advance to next flux and Jacobian row
					++bndIdx;
					++jac;
				}
			}
		};

		typedef GPRBindingBase<GPRParamHandler> GPRBinding;
		typedef GPRBindingBase<ExtGPRParamHandler> ExternalGPRBinding;

		namespace binding
		{
			void registerGPRModel(std::unordered_map<std::string, std::function<model::IBindingModel* ()>>& bindings)
			{
				bindings[GPRBinding::identifier()] = []() { return new GPRBinding(); };
				bindings[ExternalGPRBinding::identifier()] = []() { return new ExternalGPRBinding(); };
			}
		}  // namespace binding

	}  // namespace model

}  // namespace cadet
