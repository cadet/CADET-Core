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
#include "model/binding/Class_ANN.h"
#include<iostream>
#include <functional>
#include <unordered_map>
#include <string>
#include <vector>
#include <cmath>

/*<codegen>
{
	"name": "MachineLearningParamHandler",
	"externalName": "ExtMachineLearningParamHandler",
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

		inline const char* MachineLearningParamHandler::identifier() CADET_NOEXCEPT { return "MACHINE_LEARNING"; }

		inline bool MachineLearningParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
		{
			if (_kKin.size() < nComp)
				throw InvalidParameterException("ML_KKIN has to have NCOMP entries");

			return true;
		}

		inline const char* ExtMachineLearningParamHandler::identifier() CADET_NOEXCEPT { return "EXT_MACHINE_LEARNING"; }

		inline bool ExtMachineLearningParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
		{
			if (_kKin.size() < nComp)
				throw InvalidParameterException("ML_KKIN has to have NCOMP entries");

			return true;
		}


		template <class ParamHandler_t>
		class MachineLearningBindingBase : public ParamHandlerBindingModelBase<ParamHandler_t>
		{
		public:

			MachineLearningBindingBase() : bias0(), kernel0(),
				bias1(), kernel1(), bias2(), kernel2(), offset(), normalization_factor(), porosity_factor(), first_Affine_layer_pred(),
				second_Affine_layer_pred(), first_Affine_layer(), first_Affine_layer_temp(), second_Affine_layer(),
				Jacobian_1(), Jacobian_2() {}

			virtual ~MachineLearningBindingBase() CADET_NOEXCEPT { }

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


			unsigned int number_of_layers;
			unsigned int number_of_nodes; //It is asssumed that number of nodes remains same in each hidden layer
			unsigned int dimension;
			const unsigned int number_of_input{ 1 };
			const unsigned int number_of_output{ 1 };

			std::vector<double> bias0;
			std::vector<double> kernel0;

			std::vector<double> bias1;
			std::vector<double> kernel1;

			std::vector<double> bias2;
			std::vector<double> kernel2;

			std::vector<double> offset;
			std::vector<double> normalization_factor;
			std::vector<double> porosity_factor;
			std::string Implementation;
			mutable std::vector<double> first_Affine_layer_pred;
			mutable std::vector<double> second_Affine_layer_pred;

			mutable std::vector<double> first_Affine_layer;
			mutable std::vector<double> first_Affine_layer_temp;
			mutable std::vector<double> second_Affine_layer;
			mutable std::vector<double> Jacobian_1;
			mutable std::vector<double> Jacobian_2;

			/*Question: One of the */
			virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
			{
				// Read parameters
				_paramHandler.configure(paramProvider, _nComp, _nBoundStates);

				// Register parameters
				_paramHandler.registerParameters(_parameters, unitOpIdx, parTypeIdx, _nComp, _nBoundStates);

				number_of_layers = paramProvider.getInt("LAYERS");
				number_of_nodes = paramProvider.getInt("NODES");
				dimension = paramProvider.getInt("NDIM");
				normalization_factor = paramProvider.getDoubleArray("NORM_FACTOR");
				porosity_factor = paramProvider.getDoubleArray("POROS_FACTOR");
				Implementation = paramProvider.getString("IMPL");

				double val = 0.0;
				std::vector<double> input_temp = { val * normalization_factor[0] };
				std::vector<double> offset_temp(number_of_output, 0.0);
				// Read some ML parameters

				if (number_of_layers == 2)
				{
					paramProvider.pushScope("model_weights");
					paramProvider.pushScope("layer_0");

					bias0 = paramProvider.getDoubleArray("BIAS"); //Reading the biad vector from input file
					kernel0 = paramProvider.getDoubleArray("KERNEL");


					paramProvider.popScope(); // model_weights


					paramProvider.pushScope("layer_1");

					bias1 = paramProvider.getDoubleArray("BIAS");
					kernel1 = paramProvider.getDoubleArray("KERNEL");


					paramProvider.popScope(); // model_weights

					paramProvider.pushScope("layer_2");

					bias2 = paramProvider.getDoubleArray("BIAS");
					kernel2 = paramProvider.getDoubleArray("KERNEL");

					paramProvider.popScope(); // model_weights
					paramProvider.popScope(); // adsorption


					Class_ANN ANN(number_of_layers, number_of_nodes, dimension, number_of_input, number_of_output);

					ANN.prediction_two_layers_native(input_temp.data(), kernel0.data(), kernel1.data(), kernel2.data(),
						bias0.data(), bias1.data(), bias2.data(), offset_temp.data());

					offset = offset_temp;
					std::cout << "Test: " << offset[0] << "\n";
				}
				else if (number_of_layers == 1)
				{
					paramProvider.pushScope("model_weights");
					paramProvider.pushScope("layer_0");

					bias0 = paramProvider.getDoubleArray("BIAS"); //Reading the biad vector from input file
					kernel0 = paramProvider.getDoubleArray("KERNEL");


					paramProvider.popScope(); // model_weights


					paramProvider.pushScope("layer_1");

					bias1 = paramProvider.getDoubleArray("BIAS");
					kernel1 = paramProvider.getDoubleArray("KERNEL");

					paramProvider.popScope(); // model_weights
					paramProvider.popScope(); // adsorption

					for (unsigned int i = 0; i < dimension * number_of_output; i++)
						offset_temp[i] = bias1[i];

					Class_ANN ANN(number_of_layers, number_of_nodes, dimension, number_of_input, number_of_output);

					ANN.prediction_single_layer(input_temp.data(), kernel0.data(), kernel1.data(), bias0.data(),
						bias1.data(), offset_temp.data());

					offset = offset_temp;
					std::cout << "Test: " << offset[0] << "\n";

				}
				else
					throw InvalidParameterException("Number of layers: " + std::to_string(number_of_layers) + " is not supported. Select 1 or 2.");

				/**
				 * JAZIB: Read ML parameters and model data here (Done)
				 */

				return true;
			}

			inline void Matrix_vector_bias(const double* A, int rowA, int colA, const double* B, const double* C, double* result, bool flag) const
			{
				double sum;

				for (unsigned int i = 0; i < colA; ++i)
				{
					sum = 0;
					for (unsigned int j = 0; j < rowA; ++j)
						sum += A[i + j * colA] * B[j];
					result[i] = sum + C[i];
					// Applying Elu activation function
					if (flag)
					{
						if (result[i] < 0)
							result[i] = exp(result[i]) - 1;
					}
				}

			}

			inline void Matrix_vector_bias_backward(const double* A, int rowA, int colA, const double* B, double* result) const
			{
				double sum;

				for (unsigned int i = 0; i < rowA; ++i)
				{
					sum = 0;
					for (unsigned int j = 0; j < colA; ++j)
						sum += A[j + i * rowA] * B[j];
					result[i] = sum;
				}

			}

			inline void ReLu_vec(double* Affine, unsigned int size) const
			{
				for (unsigned int i = 0; i < size; i++)
				{
					if (Affine[i] < 0)
						Affine[i] = exp(Affine[i]) - 1;
				}
			}

			inline void Elu_backward(double* output, const double* matrix_1, unsigned int size) const
			{
				for (int i = 0; i < size; i++)
				{
					if (matrix_1[i] < 0)
						output[i] = exp(matrix_1[i]) * output[i];
				}
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

				double* prediction = new double[dimension * number_of_output];

				double* matCp = new double[_nComp];

				//std::vector<double> prediction(_nComp, 0.0);
				for (unsigned int i = 0; i < _nComp; i++)
				{
					matCp[i] = cp[i] * normalization_factor[i];
					prediction[i] = 0.0;
				}


				if (number_of_layers == 2)
				{
					first_Affine_layer_pred = bias0;
					second_Affine_layer_pred = bias1;

					if (Implementation == "Native")
					{
						Matrix_vector_bias(kernel0.data(), _nComp, number_of_nodes, matCp, bias0.data(), first_Affine_layer_pred.data(), true);

						Matrix_vector_bias(kernel1.data(), number_of_nodes, number_of_nodes, first_Affine_layer_pred.data(), bias1.data(), second_Affine_layer_pred.data(), true);

						Matrix_vector_bias(kernel2.data(), number_of_nodes, _nComp, second_Affine_layer_pred.data(), bias2.data(), prediction, false);

					}
					else if (Implementation == "Lapack")
					{

						for (unsigned int i = 0; i < _nComp; i++)
						{
							prediction[i] = bias2[i];
						}
						cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							dimension, number_of_nodes, number_of_input, 1.0, matCp, number_of_input,
							kernel0.data(), number_of_nodes, 1.0, first_Affine_layer_pred.data(), number_of_nodes);

						ReLu_vec(first_Affine_layer_pred.data(), number_of_nodes);

						cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							dimension, number_of_nodes, number_of_nodes, 1.0,
							first_Affine_layer_pred.data(), number_of_nodes,
							kernel1.data(), number_of_nodes, 1.0, second_Affine_layer_pred.data(), number_of_nodes);

						ReLu_vec(second_Affine_layer_pred.data(), number_of_nodes);

						cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							dimension, number_of_output, number_of_nodes, 1.0,
							second_Affine_layer_pred.data(), number_of_nodes, kernel2.data(),
							number_of_output, 1.0, prediction, number_of_output);

						/*Class_ANN ANN(number_of_layers, number_of_nodes, dimension, number_of_input, number_of_output);
						ANN.prediction_two_layers(matCp, kernel0.data(), kernel1.data(), kernel2.data(),
							first_Affine_layer_pred.data(), second_Affine_layer_pred.data(), prediction);*/
					}
					else
						throw InvalidParameterException("Select a correct implementation choice: 1. Native or 2.Lapack ");
				}
				else if (number_of_layers == 1)
				{
					first_Affine_layer_pred = bias0;
					if (Implementation == "Native")
					{
						Matrix_vector_bias(kernel0.data(), _nComp, number_of_nodes, matCp, bias0.data(), first_Affine_layer_pred.data(), true);

						Matrix_vector_bias(kernel1.data(), number_of_nodes, _nComp, first_Affine_layer_pred.data(), bias1.data(), prediction, false);
					}
					else if (Implementation == "Lapack")
					{
						for (unsigned int i = 0; i < dimension * number_of_output; i++)
							prediction[i] = bias1[i];

						cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							dimension, number_of_nodes, number_of_input, 1.0, matCp, number_of_input,
							kernel0.data(), number_of_nodes, 1.0, first_Affine_layer_pred.data(), number_of_nodes);

						ReLu_vec(first_Affine_layer_pred.data(), number_of_nodes);

						cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							dimension, number_of_output, number_of_nodes, 1.0,
							first_Affine_layer_pred.data(), number_of_nodes, kernel1.data(),
							number_of_output, 1.0, prediction, number_of_output);
						//ANN.prediction_single_layer(matCp.data(), kernel0.data(), kernel1.data(),
							//bias0.data(), bias1.data(), prediction);
					}
					
				}

				for (unsigned int i = 0; i < _nComp; i++)
				{
					q[i] = (prediction[i] - offset[0]) * porosity_factor[0]; //- offset[0]
					//std::cout << cp[i] << "\t" <<q[i] << "\n";
				}

				delete[] prediction;
				delete[] matCp;
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
				//std::cout << "Time : " << t << "\n";

				double* jacobian = new double[dimension * number_of_output];

				double* matCp = new double[_nComp];
				for (int i = 0; i < _nComp; ++i)
				{
					matCp[i] = yCp[i] * normalization_factor[0];
				}


				if (number_of_layers == 2)
				{
					first_Affine_layer = bias0;
					first_Affine_layer_temp = first_Affine_layer;
					second_Affine_layer = bias1;
					Jacobian_1 = kernel2;
					Jacobian_2 = Jacobian_1;

					if (Implementation == "Native")
					{

						//X = Ax + a
						Matrix_vector_bias(kernel0.data(), number_of_input, number_of_nodes, matCp, bias0.data(), first_Affine_layer.data(), false);

						first_Affine_layer_temp = first_Affine_layer;

						// X' = Elu(X)
						ReLu_vec(first_Affine_layer.data(), number_of_nodes);

						// Y = X'*B + b
						Matrix_vector_bias(kernel1.data(), number_of_nodes, number_of_nodes, first_Affine_layer.data(), bias1.data(), second_Affine_layer.data(), false);

						// Y' = C*Y
						Elu_backward(Jacobian_1.data(), second_Affine_layer.data(), dimension * number_of_nodes);


						//Z = Y'*B^T
						Matrix_vector_bias_backward(kernel1.data(), number_of_nodes, number_of_nodes, Jacobian_1.data(), Jacobian_2.data());


						// Z' = Z*Elu'(X)
						Elu_backward(Jacobian_2.data(), first_Affine_layer_temp.data(), number_of_nodes);

						// dNN/dx = XX*A^T
						Matrix_vector_bias_backward(kernel0.data(), number_of_input, number_of_nodes, Jacobian_2.data(), jacobian);

					}
					else if (Implementation == "Lapack")
					{
						
						cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							dimension, number_of_nodes, number_of_input, 1.0, matCp,
							number_of_input, kernel0.data(), number_of_nodes, 1.0, first_Affine_layer.data(), number_of_nodes);

						first_Affine_layer_temp = first_Affine_layer;

						ReLu_vec(first_Affine_layer.data(), number_of_nodes);

						cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							dimension, number_of_nodes, number_of_nodes, 1.0, first_Affine_layer.data(), number_of_nodes,
							kernel1.data(), number_of_nodes, 1.0, second_Affine_layer.data(), number_of_nodes);

						// ZZ = C * Elu'(Z)
						Elu_backward(Jacobian_1.data(), second_Affine_layer.data(), dimension * number_of_nodes);

						// YY = ZZ*B^T
						cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
							dimension, number_of_nodes, number_of_nodes, 1.0, Jacobian_1.data(),
							number_of_nodes, kernel1.data(), number_of_nodes, 0.0, Jacobian_2.data(), number_of_nodes);

						// XX = YY*Elu'(X)
						Elu_backward(Jacobian_2.data(), first_Affine_layer_temp.data(), number_of_nodes);

						// dNN/dx = XX*A^T
						cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							dimension, number_of_input, number_of_nodes, 1.0, Jacobian_2.data(),
							number_of_nodes, kernel0.data(), number_of_input, 0.0, jacobian, number_of_input);

						//Class_ANN ANN(number_of_layers, number_of_nodes, dimension, number_of_input, number_of_output);
						//ANN.jacobian_two_layers(matCp, kernel0.data(), kernel1.data(), kernel2.data(),
							//bias0.data(), bias1.data(), first_Affine_layer.data(), first_Affine_layer_temp.data(),
							//second_Affine_layer.data(), Jacobian_1.data(), Jacobian_2.data(), jacobian);
					}

					//ANN.jacobian_two_layers_native(matCp.data(), kernel0.data(), kernel1.data(), kernel2.data(),
						//bias0.data(), bias1.data(), bias2.data(), jacobian);
				}
				else if (number_of_layers == 1)
				{
					first_Affine_layer = bias0;
					Jacobian_1 = kernel1;
					
					if (Implementation == "Lapack")
					{
						cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							dimension, number_of_nodes, number_of_input, 1.0, matCp, number_of_input, kernel0.data(),
							number_of_nodes, 1.0, first_Affine_layer.data(), number_of_nodes);

						// XX = B*Elu'(X)
						Elu_backward(Jacobian_1.data(), first_Affine_layer.data(), dimension * number_of_nodes);

						// dNN/dx = XX*A^T
						cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
							dimension, number_of_input, number_of_nodes, 1.0, Jacobian_1.data(), number_of_nodes,
							kernel0.data(), number_of_input, 0.0, jacobian, number_of_input);
					}
					else if (Implementation == "Native")
					{
						//X = Ax + a
						Matrix_vector_bias(kernel0.data(), number_of_input, number_of_nodes, matCp, bias0.data(), first_Affine_layer.data(), false);

						// XX = B*Elu'(X)
						Elu_backward(Jacobian_1.data(), first_Affine_layer.data(), dimension * number_of_nodes);

						// dNN/dx = XX*A^T
						Matrix_vector_bias_backward(kernel0.data(), number_of_input, number_of_nodes, Jacobian_1.data(), jacobian);
					}
					//for (unsigned int i = 0; i < dimension * number_of_output; i++)
						//jacobian[i] = 0.0;
					//ANN.jacobian_single_layer(matCp.data(), kernel0.data(), kernel1.data(), bias0.data(), bias1.data(), jacobian);

				}


				unsigned int bndIdx = 0;
				for (int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states (bound state index bndIdx is not advanced)
					if (_nBoundStates[i] == 0)
						continue;

					const double kkin = static_cast<double>(p->kKin[i]);

					for (int j = 0; j < _nComp; ++j)
					{
						jac[j - bndIdx - offsetCp] = -kkin * (jacobian[j] * porosity_factor[0]);
						//std::cout << "jac is :" << jac[j - bndIdx - offsetCp] << "\n";
						//if (t< 10.0)
							//std::cout << "Time is :" << t  << "\n";
					}

					// dres_i / dq_i
					jac[0] = kkin;

					// Advance to next flux and Jacobian row
					++bndIdx;
					++jac;
				}
				delete[] jacobian;
				delete[] matCp;
			}
		};

		typedef MachineLearningBindingBase<MachineLearningParamHandler> MachineLearningBinding;
		typedef MachineLearningBindingBase<ExtMachineLearningParamHandler> ExternalMachineLearningBinding;

		namespace binding
		{
			void registerMachineLearningModel(std::unordered_map<std::string, std::function<model::IBindingModel* ()>>& bindings)
			{
				bindings[MachineLearningBinding::identifier()] = []() { return new MachineLearningBinding(); };
				bindings[ExternalMachineLearningBinding::identifier()] = []() { return new ExternalMachineLearningBinding(); };
			}
		}  // namespace binding

	}  // namespace model

}  // namespace cadet
