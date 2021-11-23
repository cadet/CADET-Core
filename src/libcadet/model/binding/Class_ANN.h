#include <cstdio>
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include "../../../mkl_files/mkl.h"
//#include "../../../../ThirdParty/sundials/include/sundials/sundials_lapack.h"
#include <chrono>
#include <stdexcept>

#pragma once
class Class_ANN
{
protected:
	unsigned int number_of_layers;
	unsigned int number_of_nodes; //It is asssumed that number of nodes remains same in each hidden layer
	unsigned int dimension;
	unsigned int number_of_input;
	unsigned int number_of_output;
	
public:
	Class_ANN(unsigned int num_layer, unsigned num_nodes, unsigned int dims, unsigned int num_inputs, unsigned int num_outputs):
		number_of_layers(num_layer), number_of_nodes(num_nodes), dimension(dims), number_of_input(num_inputs),
		number_of_output(num_outputs) {}
	
	void Elu(double* Affine, unsigned int size)
	{

		for (unsigned int i = 0; i < size; i++)
		{
			if (Affine[i] < 0)
				Affine[i] = exp(Affine[i]) - 1;
		}
	}

	void Elu_backward(double* output,const double* matrix_1, unsigned int size)
	{
		for (int i = 0; i < size; i++)
		{
			if (matrix_1[i] > 0)
				output[i] = output[i];
			else
				output[i] = exp(matrix_1[i]) * output[i];
		}
	}

	void Matrix_vector_bias(const double* A, int rowA, int colA, const double* B, const double* C, double* result, int alpha)
	{
		double sum;

		for (unsigned int i = 0; i < colA; ++i)
		{
			sum = 0;
			for (unsigned int j = 0; j < rowA; ++j)
				sum += A[i + j * colA] * B[j];
			result[i] = sum + alpha * C[i];
		}

	}	
	
	void Matrix_vector_bias_backward(const double* A, int rowA, int colA, const double* B, double* result)
	{
		double sum;

		for (unsigned int i = 0; i < rowA; ++i)
		{
			sum = 0.0;
			for (unsigned int j = 0; j < colA; ++j)
				sum += A[j + i * rowA] * B[j];
			result[i] = sum;
		}

	}

	void Initialization(std::vector<double>Vector)
	{
		for (unsigned int i = 0; i < Vector.size(); ++i)
			Vector[i] = 0.0;
	}

	void prediction_two_layers_native(const double* input, const double* weights_hidden_layer_1, const double* weights_hidden_layer_2,
		const double* weights_output_layer, const double* bias_hidden_layer_1, const double* bias_hidden_layer_2,
		const double* bias_output_layer, double* prediction)
	{
		//std::vector<double> func_vec(number_of_nodes * number_of_input, 0.0);
		//std::vector<double> func_vec_second(number_of_nodes * number_of_input, 0.0);
		double* func_vec = new double[number_of_nodes * number_of_input];
		double* func_vec_second = new double[number_of_nodes * number_of_input];

		double sum;

		// From input to first hidden layer
		for (unsigned int i = 0; i < number_of_nodes; ++i)
		{
			sum = 0.0;
			for (unsigned int j = 0; j < number_of_input; ++j)
				sum += weights_hidden_layer_1[i + j * number_of_nodes] * input[j];
			func_vec[i] = sum + bias_hidden_layer_1[i];
			if (func_vec[i] < 0) // Applying Elu activation function
				func_vec[i] = exp(func_vec[i]) - 1;
		}
		
		// From first hidden layer to second hidden layer
		for (unsigned int i = 0; i < number_of_nodes; ++i)
		{
			sum = 0.0;
			for (unsigned int j = 0; j < number_of_nodes; ++j)
				sum += weights_hidden_layer_2[i + j * number_of_nodes] * func_vec[j];
			func_vec_second[i] = sum + bias_hidden_layer_2[i];
			if (func_vec_second[i] < 0) // Applying Elu activation function
				func_vec_second[i] = exp(func_vec_second[i]) - 1;
		}

		// From second hidden layer to third hidden layer
		for (unsigned int i = 0; i < number_of_output; ++i)
		{
			sum = 0.0;
			for (unsigned int j = 0; j < number_of_nodes; ++j)
				sum += weights_output_layer[i + j * number_of_output] * func_vec_second[j];
			prediction[i] = sum + bias_output_layer[i];
		}
		
		delete [] func_vec;
		delete [] func_vec_second;
	}

	void prediction_two_layers(const double* input, const double* weights_hidden_layer_1, const double* weights_hidden_layer_2,
		const double* weights_output_layer, double* first_Affine_layer_pred,  double* second_Affine_layer_pred,
		double* prediction)
	{
		
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_nodes, number_of_input, 1.0, input, number_of_input, weights_hidden_layer_1, number_of_nodes, 1.0, first_Affine_layer_pred, number_of_nodes);

		Elu(first_Affine_layer_pred, number_of_nodes);

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_nodes, number_of_nodes, 1.0, first_Affine_layer_pred, number_of_nodes, weights_hidden_layer_2, number_of_nodes, 1.0, second_Affine_layer_pred, number_of_nodes);


		Elu(second_Affine_layer_pred, number_of_nodes);

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_output, number_of_nodes, 1.0, second_Affine_layer_pred, number_of_nodes, weights_output_layer, number_of_output, 1.0, prediction, number_of_output);
				
	}

	void prediction_single_layer(const double* input, const double* weights_hidden_layer_1,
		const double* weights_output_layer, const double* bias_hidden_layer_1, const double* bias_output_layer, double* prediction)
	{
		/*double* first_Affine_layer = new double[number_of_nodes];

		for (unsigned int i = 0; i < number_of_nodes; ++i)
		{
			first_Affine_layer[i] = bias_hidden_layer_1[i];
		}


		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_nodes, number_of_input, 1.0, input, number_of_input, weights_hidden_layer_1, number_of_nodes, 1.0, first_Affine_layer, number_of_nodes);

		Elu(first_Affine_layer, number_of_nodes);

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_output, number_of_nodes, 1.0, first_Affine_layer, number_of_nodes, weights_output_layer, number_of_output, 1.0, prediction, number_of_output);
		
		delete[] first_Affine_layer;*/
	}

	void jacobian_two_layers(const double* input, const double* weights_hidden_layer_1, const double* weights_hidden_layer_2,
		const double* weights_output_layer, const double* bias_hidden_layer_1,
		const double* bias_hidden_layer_2, double* first_Affine_layer, double * first_Affine_layer_temp,
		double* second_Affine_layer, double* Jacobian_1, double* Jacobian_2,  double* dout)
	{
		/*double* first_Affine_layer = new double[number_of_nodes];
		double* first_Affine_layer_temp = new double[number_of_nodes];
		double* second_Affine_layer = new double[number_of_nodes];
		double* Jacobian_1 = new double[number_of_nodes];
		double* Jacobian_2 = new double[number_of_nodes];*/
				
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_nodes, number_of_input, 1.0, input, number_of_input, weights_hidden_layer_1, number_of_nodes, 1.0, first_Affine_layer, number_of_nodes);
		
		for(unsigned int i =0 ; i< number_of_nodes; ++i)
			first_Affine_layer_temp[i] = first_Affine_layer[i];

		Elu(first_Affine_layer, number_of_nodes);


		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_nodes, number_of_nodes, 1.0, first_Affine_layer, number_of_nodes, weights_hidden_layer_2, number_of_nodes, 1.0, second_Affine_layer, number_of_nodes);


		// ZZ = C * Elu'(Z)
		Elu_backward(Jacobian_1, second_Affine_layer, dimension * number_of_nodes);



		// YY = ZZ*B^T
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
			dimension, number_of_nodes, number_of_nodes, 1.0, Jacobian_1, number_of_nodes, weights_hidden_layer_2, number_of_nodes, 0.0, Jacobian_2, number_of_nodes);

		// XX = YY*Elu'(X)
		Elu_backward(Jacobian_2, first_Affine_layer_temp, number_of_nodes);


		// dNN/dx = XX*A^T
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_input, number_of_nodes, 1.0, Jacobian_2, number_of_nodes, weights_hidden_layer_1, number_of_input, 0.0, dout, number_of_input);
				
	}

	void jacobian_two_layers_native(const double* input, const double* weights_hidden_layer_1, const double* weights_hidden_layer_2,
		const double* weights_output_layer, const double* bias_hidden_layer_1, const double* bias_hidden_layer_2,
		const double* bias_output_layer, double* jacobian)
	{
		double* func_vec = new double[number_of_nodes * number_of_input];
		double* func_vec_temp = new double[number_of_nodes * number_of_input];

		double* func_vec_second = new double[number_of_nodes * number_of_input];

		//std::vector<double> func_vec(number_of_nodes * number_of_input, 0.0);
		//std::vector<double> func_vec_temp(number_of_nodes * number_of_input, 0.0);

		//std::vector<double> func_vec_second(number_of_nodes * number_of_input, 0.0);
		double* Jacobian_1 = new double[number_of_nodes * number_of_input];
		double* Jacobian_2 = new double[number_of_nodes];
		
		//std::vector<double> Jacobian_1(number_of_input * number_of_nodes, 0.0);
		//std::vector<double> Jacobian_2(number_of_nodes, 0.0);


		for (unsigned int i = 0; i < number_of_input * number_of_nodes; ++i)
			Jacobian_1[i] = weights_output_layer[i];
		
		double sum = 0.0;
		//X = Ax + a

		// From input to first hidden layer
		for (unsigned int i = 0; i < number_of_nodes; ++i)
		{
			sum = 0.0;
			for (unsigned int j = 0; j < number_of_input; ++j)
				sum += weights_hidden_layer_1[i + j * number_of_nodes] * input[j];
			func_vec[i] = sum + bias_hidden_layer_1[i];
			func_vec_temp[i] = func_vec[i];
			if (func_vec[i] < 0) // Applying Elu activation function
				func_vec[i] = exp(func_vec[i]) - 1;
		}

		
		// Y = X'*B + b
		// From first hidden layer to second hidden layer
		for (unsigned int i = 0; i < number_of_nodes; ++i)
		{
			sum = 0;
			for (unsigned int j = 0; j < number_of_nodes; ++j)
				sum += weights_hidden_layer_2[i + j * number_of_nodes] * func_vec[j];
			func_vec_second[i] = sum + bias_hidden_layer_2[i];
			// Y' = C*Y
			if (func_vec_second[i]<0)
				Jacobian_1[i] = exp(func_vec_second[i]) * Jacobian_1[i];
		}

				
		//Z = Y'*B^T
		for (unsigned int i = 0; i < number_of_nodes; ++i)
		{
			sum = 0;
			for (unsigned int j = 0; j < number_of_nodes; ++j)
				sum += weights_hidden_layer_2[j + i * number_of_nodes] * Jacobian_1[j];
			Jacobian_2[i] = sum;
			// Z' = Z*Elu'(X)
			if (func_vec_temp[i] < 0) 
				Jacobian_2[i] = exp(func_vec_temp[i]) * Jacobian_2[i];
		}
		
		
		// dNN/dx = XX*A^T
		for (unsigned int i = 0; i < number_of_output; ++i)
		{
			sum = 0.0;
			for (unsigned int j = 0; j < number_of_nodes; ++j)
				sum += weights_hidden_layer_1[j + i * number_of_output] * Jacobian_2[j];
			jacobian[i] = sum;
		}
		
		delete [] func_vec;
		delete [] func_vec_temp;
		delete [] func_vec_second;
		delete [] Jacobian_1;
		delete [] Jacobian_2;
	}

	void jacobian_single_layer(const double* input, const double* weights_hidden_layer_1,
		const double* weights_output_layer, const double* bias_hidden_layer_1, const double* bias_output_layer, double* dout)
	{

		/*double* first_Affine_layer = new double[number_of_nodes];
		double* Jacobian_1 = new double[number_of_nodes];

		for (unsigned int i = 0; i < number_of_nodes; ++i)
		{
			first_Affine_layer[i] = bias_hidden_layer_1[i];
			Jacobian_1[i] = weights_output_layer[i];
		}


		// X= Ax+a
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_nodes, number_of_input, 1.0, input, number_of_input, weights_hidden_layer_1, number_of_nodes, 1.0, first_Affine_layer, number_of_nodes);

		// XX = B*Elu'(X)
		Elu_backward(Jacobian_1, first_Affine_layer, dimension * number_of_nodes);

		// dNN/dx = XX*A^T
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_input, number_of_nodes, 1.0, Jacobian_1, number_of_nodes, weights_hidden_layer_1, number_of_input, 0.0, dout, number_of_input);

		delete[] first_Affine_layer;
		delete[] Jacobian_1;*/
	}

};

//Class_ANN ANN(number_of_hidden_layers, number_of_nodes_W1, num_dims, number_of_input, number_of_outputs);
/*pred = ANN.prediction_two_layers(input.data(),layer_0_kernel0.data(),layer_1_kernel0.data(),layer_2_kernel0.data(),
		layer_0_bias0.data(),layer_1_bias0.data(),layer_2_bias0.data());*/