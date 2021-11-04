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
	Class_ANN(unsigned int num_layer, unsigned num_nodes, unsigned int dims, unsigned int num_inputs, unsigned int num_outputs) :
		number_of_layers(num_layer), number_of_nodes(num_nodes),
		dimension(dims), number_of_input(num_inputs), number_of_output(num_outputs) {}

	void Elu(double* Affine, unsigned int size)
	{

		for (unsigned int i = 0; i < size; i++)
		{
			if (Affine[i] > 0)
				Affine[i] = Affine[i];
			else
				Affine[i] = exp(Affine[i]) - 1;
		}
	}

	void Elu_backward(double* output, double* matrix_1, unsigned int size)
	{
		for (int i = 0; i < size; i++)
		{
			if (matrix_1[i] > 0)
				output[i] = output[i];
			else
				output[i] = exp(matrix_1[i]) * output[i];
		}
	}

	void Initialization(std::vector<double>Vector)
	{
		for (unsigned int i = 0; i < Vector.size(); ++i)
			Vector[i] = 0.0;
	}

	void prediction_two_layers(const double* input, std::vector<double> weights_hidden_layer_1, std::vector<double> weights_hidden_layer_2,
		std::vector<double> weights_output_layer, std::vector<double> bias_hidden_layer_1, std::vector<double> bias_hidden_layer_2, std::vector<double> bias_output_layer,
		double* prediction)
	{
		std::vector<double> first_Affine_layer = bias_hidden_layer_1;
		std::vector<double>second_Affine_layer = bias_hidden_layer_2;


		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_nodes, number_of_input, 1.0, input, number_of_input, weights_hidden_layer_1.data(), number_of_nodes, 1.0, first_Affine_layer.data(), number_of_nodes);

		Elu(first_Affine_layer.data(), number_of_nodes);

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_nodes, number_of_nodes, 1.0, first_Affine_layer.data(), number_of_nodes, weights_hidden_layer_2.data(), number_of_nodes, 1.0, second_Affine_layer.data(), number_of_nodes);


		Elu(second_Affine_layer.data(), number_of_nodes);

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_output, number_of_nodes, 1.0, second_Affine_layer.data(), number_of_nodes, weights_output_layer.data(), number_of_output, 1.0, prediction, number_of_output);
		
		Initialization(first_Affine_layer);
		Initialization(second_Affine_layer);
		
	}

	void prediction_single_layer(const double* input, std::vector<double> weights_hidden_layer_1,
		std::vector<double> weights_output_layer, std::vector<double> bias_hidden_layer_1, std::vector<double> bias_output_layer, double* prediction)
	{
		double* first_Affine_layer = bias_hidden_layer_1.data();


		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_nodes, number_of_input, 1.0, input, number_of_input, weights_hidden_layer_1.data(), number_of_nodes, 1.0, first_Affine_layer, number_of_nodes);

		Elu(first_Affine_layer, number_of_nodes);

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_output, number_of_nodes, 1.0, first_Affine_layer, number_of_nodes, weights_output_layer.data(), number_of_output, 1.0, prediction, number_of_output);
		
	}

	void jacobian_two_layers(const double* input, std::vector<double> weights_hidden_layer_1, std::vector<double> weights_hidden_layer_2,
		std::vector<double> weights_output_layer, std::vector<double> bias_hidden_layer_1, std::vector<double> bias_hidden_layer_2, std::vector<double> bias_output_layer, double* dout)
	{
		std::vector<double> first_Affine_layer = bias_hidden_layer_1;
		std::vector<double> second_Affine_layer = bias_hidden_layer_2;
		std::vector<double> Jacobian_1 = weights_output_layer;
		std::vector<double> Jacobian_2(number_of_nodes, 0.0);
		

		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_nodes, number_of_input, 1.0, input, number_of_input, weights_hidden_layer_1.data(), number_of_nodes, 1.0, first_Affine_layer.data(), number_of_nodes);
		std::vector<double> first_Affine_layer_temp = first_Affine_layer;

		Elu(first_Affine_layer.data(), number_of_nodes);


		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_nodes, number_of_nodes, 1.0, first_Affine_layer.data(), number_of_nodes, weights_hidden_layer_2.data(), number_of_nodes, 1.0, second_Affine_layer.data(), number_of_nodes);


		// ZZ = C * Elu'(Z)
		Elu_backward(Jacobian_1.data(), second_Affine_layer.data(), dimension * number_of_nodes);



		// YY = ZZ*B^T
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
			dimension, number_of_nodes, number_of_nodes, 1.0, Jacobian_1.data(), number_of_nodes, weights_hidden_layer_2.data(), number_of_nodes, 0.0, Jacobian_2.data(), number_of_nodes);

		// XX = YY*Elu'(X)
		Elu_backward(Jacobian_2.data(), first_Affine_layer_temp.data(), Jacobian_2.size());


		// dNN/dx = XX*A^T
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_input, number_of_nodes, 1.0, Jacobian_2.data(), number_of_nodes, weights_hidden_layer_1.data(), number_of_input, 0.0, dout, number_of_input);

		Initialization(first_Affine_layer);
		Initialization(second_Affine_layer);
		Initialization(Jacobian_1);
		Initialization(Jacobian_2);

	}

	void jacobian_single_layer(const double* input, std::vector<double> weights_hidden_layer_1,
		std::vector<double> weights_output_layer, std::vector<double> bias_hidden_layer_1, std::vector<double> bias_output_layer, double* dout)
	{
		std::vector<double> first_Affine_layer = bias_hidden_layer_1;
		std::vector<double> Jacobian_1 = weights_output_layer;

		// X= Ax+a
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_nodes, number_of_input, 1.0, input, number_of_input, weights_hidden_layer_1.data(), number_of_nodes, 1.0, first_Affine_layer.data(), number_of_nodes);

		// XX = B*Elu'(X)
		Elu_backward(Jacobian_1.data(), first_Affine_layer.data(), dimension * number_of_nodes);

		// dNN/dx = XX*A^T
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
			dimension, number_of_input, number_of_nodes, 1.0, Jacobian_1.data(), number_of_nodes, weights_hidden_layer_1.data(), number_of_input, 0.0, dout, number_of_input);

	}

};

//Class_ANN ANN(number_of_hidden_layers, number_of_nodes_W1, num_dims, number_of_input, number_of_outputs);
/*pred = ANN.prediction_two_layers(input.data(),layer_0_kernel0.data(),layer_1_kernel0.data(),layer_2_kernel0.data(),
		layer_0_bias0.data(),layer_1_bias0.data(),layer_2_bias0.data());*/