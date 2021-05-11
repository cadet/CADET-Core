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

			MachineLearningBindingBase() : bias_vector_1(number_of_input, std::vector<double>(number_of_nodes_W1)),
				first_hidden_layer_W1(number_of_input, std::vector<double>(number_of_nodes_W1)),
				bias_vector_2(number_of_input, std::vector<double>(number_of_nodes_W2)),
				second_hidden_layer_W2(number_of_nodes_W1, std::vector<double>(number_of_nodes_W2)),
				bias_vector_3(number_of_outputs, std::vector<double>(number_of_input)),
				last_layer(number_of_nodes_W2, std::vector<double>(number_of_outputs)),
				Affine_layer_1(number_of_input, std::vector<double>(number_of_nodes_W1, 0.0)),
				Affine_layer1_Jacobian(number_of_input, std::vector<double>(number_of_nodes_W1, 0.0)),
				Affine_layer_2(number_of_input, std::vector<double>(number_of_nodes_W2, 0.0)),
				Affine_layer_2_Jacobian(number_of_input, std::vector<double>(number_of_nodes_W2, 0.0)),
				Affine_layer_3(number_of_input, std::vector<double>(number_of_outputs, 0.0)) {}

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


			const int number_of_input{ 1 };
			const int number_of_hidden_layers{ 2 };
			const int number_of_outputs{ 1 };

			const int number_of_nodes_W1{ 75 };
			const int number_of_nodes_W2{ 75 };
			/*******************************************************/
			/* Defining 2D vectors to copy the 1D format input data*/
			/*******************************************************/
			std::vector< std::vector<double> > bias_vector_1;
			std::vector< std::vector<double> > first_hidden_layer_W1;
			std::vector< std::vector<double> > bias_vector_2;
			std::vector< std::vector<double> > second_hidden_layer_W2;
			std::vector< std::vector<double> > bias_vector_3;
			std::vector< std::vector<double> > last_layer;

			/***************************************************************************************************/
			/***Defining the vectors that will be used in the calculation of predicted value of isotherm*******/
			/***************************************************************************************************/
			mutable std::vector< std::vector<double> > Affine_layer_1;
			mutable std::vector< std::vector<double> > Affine_layer1_Jacobian;
			mutable std::vector< std::vector<double> > Affine_layer_2;
			mutable std::vector< std::vector<double> > Affine_layer_2_Jacobian;
			mutable std::vector< std::vector<double> > Affine_layer_3;

			/***************************************************************************************************/

		/* @breif: Converting 1D vectors read from HDF5 file to 2D vectors
			to facilitate in matrix multiplication for NN process */
			template<class T>
			void oneD_to_twoD(const std::vector<T>(&Matrix1D), std::vector<std::vector<T>>(&Matrix2D))
			{
				int rows = Matrix2D.size();
				int cols = Matrix2D[0].size();

				for (int i = 0; i < rows; i++)
					for (int j = 0; j < cols; j++)
						Matrix2D[i][j] = Matrix1D[j + i * cols];
			}

			// JAZIB: Add ML model parameters here
			// @brief: Initialising matrices fucntion with all zero
			template<class T>
			void Intialisation(T(&matrix)) const
			{
				int rows = matrix.size();
				int cols = matrix[0].size();

				static int count = 1;
				//cout << "Initialising Layer to zeros " << endl;
				for (int i = 0; i < rows; i++)
				{
					for (int j = 0; j < cols; j++)
					{
						matrix[i][j] = 0.0;
						//cout << matrix[i][j] << "\t";
					}
					//cout << endl;
				}
			}

			// @brief: creating a diagonal matrices of all 1's
			template<class t>
			void diagonal_one(t(&matrix)) const
			{
				static int count = 1;
				int rows = matrix.size();
				int cols = matrix[0].size();

				//cout << "diagonal matrix: " << endl;
				for (int i = 0; i < rows; i++)
				{
					for (int j = 0; j < cols; j++)
					{
						if (i == j)
							matrix[i][j] = 1.0;
						else
							matrix[i][j] = 0.0;
					}
				}

			}

			/* @Brief: Propagating the input layer to the first hidden layer and also using activation function right after
			Out: void function
			In params :
			Affine_layer1 contains the results of a* x + b
			Input layer : input values 1x1 is used for this case
			first_hidden_layer_W1: Trained weight matrix of first hidden layer
			bias_vector_1 : trained bias vector for 1st hidden layer*/
			template <class T>
			void Matrix_mul(T(&output), const T(&matrix_1), const T(&matrix_2), const T(&bias)) const
			{
				int rows_1 = matrix_1.size();
				int cols_1 = matrix_1[0].size();
				int rows_2 = matrix_2.size();
				int cols_2 = matrix_2[0].size();


				for (int i = 0; i < rows_1; ++i)
					for (int j = 0; j < cols_2; ++j)
						for (int k = 0; k < cols_1; ++k)
						{
							output[i][j] += matrix_1[i][k] * matrix_2[k][j];
						}
				//Adding bias vector
				for (int i = 0; i < rows_1; i++)
					for (int j = 0; j < cols_2; j++)
						output[i][j] = output[i][j] + bias[i][j];
			}


			// @brief: Activation fucntion: ReLu

			void ReLu(std::vector< std::vector<double>>(&Affine)) const
			{
				int rows = Affine.size();
				int cols = Affine[0].size();
				const int one = 1;

				for (int i = 0; i < rows; i++)
					for (int j = 0; j < cols; j++)
					{
						if (Affine[i][j] > 0)
							Affine[i][j] = one * Affine[i][j];
						else
							Affine[i][j] = (exp(Affine[i][j])) - 1;
					}
			}

			void ReLu_New(std::vector< std::vector<double>>(&Affine)) const
			{
				int rows = Affine.size();
				int cols = Affine[0].size();
				const int one = 1;

				for (int i = 0; i < rows; i++)
					for (int j = 0; j < cols; j++)
					{
						if (Affine[i][j] > 0)
							Affine[i][j] = one * Affine[i][j];
						else
							Affine[i][j] = 0.0;
					}
			}

			//template<class T>
			void Affine_Backward(std::vector< std::vector<double> >(&output), const std::vector< std::vector<double> >(&matrix_1), const std::vector< std::vector<double> >(&matrix_2)) const
			{
				int rows_1 = matrix_1.size();
				int cols_1 = matrix_1[0].size();
				int rows_2 = matrix_2.size();
				int cols_2 = matrix_2[0].size();


				std::vector< std::vector<double> > transpose(cols_2, std::vector<double>(rows_2, 0.0));
				/*Add try catch exception over here*/


				/*Taking transpose of the Matrix_1 to allow multiplication of Argument 1 and Argument 2*/
				for (int i = 0; i < rows_2; ++i)
					for (int j = 0; j < cols_2; ++j) {
						transpose[j][i] = matrix_2[i][j];
					}

				/*Copying dout into another matrix*/
				std::vector< std::vector<double> > dout_temp(rows_1, std::vector<double>(cols_1, 0.0));


				for (int i = 0; i < rows_1; i++)
					for (int j = 0; j < cols_1; j++)
						dout_temp[i][j] = output[i][j];

				/*Carrying out the dot product*/
				for (int i = 0; i < rows_1; ++i)
				{
					for (int j = 0; j < rows_2; ++j)
						for (int k = 0; k < cols_1; ++k)
						{
							output[i][j] += matrix_1[i][k] * transpose[k][j];
						}
				}

			}

			/*@breif: Backward pass of ReLU
		   : param output : Upstream Jacobian
		   : param matrix_1 : the cached input for this layer
		   : return : the jacobian matrix containing derivatives of the O neural network outputs with respect to
		   this layer's inputs, evaluated at x.
		   */
		   //template<class T>
			void relu_backward(std::vector< std::vector<double> >(&output), const std::vector< std::vector<double> >(&matrix_1)) const
			{
				int rows_1 = matrix_1.size();
				int cols_1 = matrix_1[0].size();
				const int one = 1;

				for (int i = 0; i < rows_1; i++)
					for (int j = 0; j < cols_1; j++)
					{
						if (matrix_1[i][j] > 0)
							output[i][j] = one * output[i][j];
						else
							output[i][j] = (exp(matrix_1[i][j])) * output[i][j];
					}
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
				paramProvider.pushScope("layer_0");


				const std::vector<double> bias0 = paramProvider.getDoubleArray("bias"); //Reading the biad vector from input file
				oneD_to_twoD(bias0, bias_vector_1); // Copying 1D bias vector into 2D vector format

				const std::vector<double> kernel0 = paramProvider.getDoubleArray("kernel");
				oneD_to_twoD(kernel0, first_hidden_layer_W1); // Copying 1D kernel vector into 2D vector format

				paramProvider.popScope(); // model_weights


				paramProvider.pushScope("layer_1");

				const std::vector<double> bias1 = paramProvider.getDoubleArray("bias");
				oneD_to_twoD(bias1, bias_vector_2); // Copying 1D bias vector into 2D vector format



				const std::vector<double> kernel1 = paramProvider.getDoubleArray("kernel");
				oneD_to_twoD(kernel1, second_hidden_layer_W2);

				paramProvider.popScope(); // model_weights

				paramProvider.pushScope("layer_2");

				const std::vector<double> bias2 = paramProvider.getDoubleArray("bias");
				oneD_to_twoD(bias2, bias_vector_3); // Copying 1D bias vector into 2D vector format

				const std::vector<double> kernel2 = paramProvider.getDoubleArray("kernel");
				oneD_to_twoD(kernel2, last_layer);

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
				const double epsilon_e = 0.399;
				double h = 0.0;
				double keq = 486.453;
				/*double qq = 155.01962744;*/
				double a = 1.3724461e-13;
				double b = 5.5858426e4;
				double epsil = 8.48484848e-4;



				double input_temp = cp[0] * keq;


				std::vector<std::vector<double>> matCp = { {input_temp} };


				Matrix_mul(Affine_layer_1, matCp, first_hidden_layer_W1, bias_vector_1);
				ReLu(Affine_layer_1);


				/* @Brief: Propagating the first hidden layer to the second hidden layer and also using activation function right after*/

				Matrix_mul(Affine_layer_2, Affine_layer_1, second_hidden_layer_W2, bias_vector_2);
				ReLu(Affine_layer_2);



				/* @Brief: Propagating the input layer to the first hidden layer */
				Matrix_mul(Affine_layer_3, Affine_layer_2, last_layer, bias_vector_3);

				if (cp[0] == 0.0)
				{
					h = Affine_layer_3[0][0]; //+ h * Affine_layer_3[0][0];//(a * pow(input_temp, 2) + b * input_temp);//+h * Affine_layer_3[0][0];
				}

				q[0] = (Affine_layer_3[0][0] + 8.0516016483306885e-01) * factor + 4.6422148705738839e-07;// 

				//q[0] = q[0] / (1.0 - 0.69);
				//q[0] = Affine_layer_3[0][0];
				//q[0] = keq * 148.422 * cp[0] / (1 + keq * cp[0]);
				//std::cout << "For Input: " << cp[0] << " Prediction is: " << q[0] << "\n";


				Intialisation(Affine_layer_1);
				Intialisation(Affine_layer_2);
				Intialisation(Affine_layer_3);

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
						const double epsilon_e = 0.399;
						double h = 0.0;
						double keq = 486.453;
						double a = 1.3724461e-13;
						double b = 5.5858426e4;
						double epsil = 1e-4;

						double input_temp = yCp[0] * keq;

						std::vector<std::vector<double>> matCp = { {input_temp} };
						//std::cout << "Input is: " << matCp[0][0] << "\n";

						h = input_temp / epsil;

						// JAZIB: Compute Jacobian of ML model
						/* ***************************************************
						**************Jacobian of the neural network**********
						******************************************************/

						std::vector<std::vector<double>>Ax_a(number_of_input, std::vector<double>(number_of_nodes_W1, 0.0));
						std::vector<std::vector<double>>HB_b(number_of_input, std::vector<double>(number_of_nodes_W2, 0.0));
						//*******************************************************
						// Computing H*B + b 

						Matrix_mul(Ax_a, matCp, first_hidden_layer_W1, bias_vector_1);
						ReLu(Ax_a);


						Matrix_mul(HB_b, Ax_a, second_hidden_layer_W2, bias_vector_2);
						ReLu(HB_b);

						//*****************************************************

						std::vector< std::vector<double> >dout(number_of_outputs, std::vector<double>(number_of_input));
						diagonal_one(dout);

						Affine_Backward(Affine_layer_2_Jacobian, dout, last_layer);

						relu_backward(Affine_layer_2_Jacobian, HB_b);

						Affine_Backward(Affine_layer1_Jacobian, Affine_layer_2_Jacobian, second_hidden_layer_W2);

						//********Again computing Ax+a for taking the elu derivative**********
						std::vector<std::vector<double>>Ax_a_second(number_of_input, std::vector<double>(number_of_nodes_W1, 0.0));
						Matrix_mul(Ax_a_second, matCp, first_hidden_layer_W1, bias_vector_1);
						//********************************************************************

						relu_backward(Affine_layer1_Jacobian, Ax_a_second);

						Intialisation(dout);
						Affine_Backward(dout, Affine_layer1_Jacobian, first_hidden_layer_W1);


						/*if (input_temp <= epsil)
						{
							dout[0][0] = (1-h)*(keq * 148.422 / (1 + keq * yCp[0] * yCp[0])) + h * keq * dout[0][0];//(2 * a * input_temp + b);//+h * keq * dout[0][0];
							//std::cout << "NN Jacobian is: " << dout[0][0] << "\n";
						}
						else
						{
							dout[0][0] = keq * dout[0][0];
						}*/

						/*if (input_temp <= 0.02212813)
						{
							h = (keq * qq) / pow(1 + keq * yCp[0], 2);
							//h = -164422620000000 * pow(input_temp, 5) + 10321455000000 * pow(input_temp, 4) - 245540280000 * pow(input_temp, 3) + 2779448100 * pow(input_temp, 2) - 15609898 * input_temp + 4.0575e4;
						}
						else
						{
							h = 5*a * pow(yCp[0], 4) + 4*b * pow(yCp[0], 3) + 3*c * pow(yCp[0], 2) + 2*d * pow(yCp[0], 1) + e ;
						}*/
						dout[0][0] = keq * factor * dout[0][0];
						// dres_i / dc_{p,j} = -kkin[i] * df_i / dc_{p,j}
						//dout[0][0] = keq * 148.422 / (1 + keq * yCp[0] * yCp[0]);
						jac[j - bndIdx - offsetCp] = -1 * (dout[0][0]);

						/*if (t > 3449.1)
						{
							std::cout << "time: " << t <<  "\n";
						}*/

						// Getting to c_{p,j}: -bndIdx takes us to q_0, another -offsetCp to c_{p,0} and a +j to c_{p,j}.
						//                     This means jac[j - bndIdx - offsetCp] corresponds to c_{p,j}.
						//std::cout << "Time " << t << " Jac is : " << jac[j - bndIdx - offsetCp]<<"\n";
						Intialisation(Affine_layer_2_Jacobian);
						Intialisation(Affine_layer1_Jacobian);
					}

					// dres_i / dq_i
					jac[0] = 1;

					// Advance to next flux and Jacobian row
					++bndIdx;
					++jac;
				}
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
