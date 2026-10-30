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
#pragma once

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

// =============================================================================
// ClassANN — Feedforward neural network with ELU activations.
//
// Architecture (N hidden layers, scalar output):
//   y = W_{N} * ELU(... ELU(W_0 * x + b_0) ...) + b_{N}
//
// Weight storage convention (matches original training export):
//   All weight matrices are stored COLUMN-MAJOR (Fortran order).
//   A weight matrix W of shape (n_out x n_in) is stored such that
//   W(row=i, col=j) = kernel[i + j * n_out].
//   This means Eigen::Map<MatrixXd>(kernel.data(), n_out, n_in) gives W directly.
//
// All hidden layers share the same width (num_nodes).
// All workspace buffers are pre-allocated at construction — no heap allocation
// inside forward/backward methods.
// =============================================================================

class ClassANN
{
public:
	// -------------------------------------------------------------------------
	// Constructor — pre-allocates all workspace buffers.
	// num_layer:   number of hidden layers (>= 1); all share width num_nodes
	// num_nodes:   nodes per hidden layer
	// num_inputs:  input dimensionality
	// num_outputs: output dimensionality (1 for single-component isotherm)
	// -------------------------------------------------------------------------
	ClassANN(unsigned int num_layer, unsigned int num_nodes,
	          unsigned int num_inputs, unsigned int num_outputs)
		: number_of_layers(num_layer)
		, number_of_nodes(num_nodes)
		, number_of_input(num_inputs)
		, number_of_output(num_outputs)
		, ws_z(num_layer, Eigen::VectorXd(num_nodes))
		, ws_h(num_layer, Eigen::VectorXd(num_nodes))
		, ws_delta_cur(num_nodes)
		, ws_delta_nxt(num_nodes)
	{
		if (num_layer < 1)
			throw std::invalid_argument(
				"ClassANN: at least 1 hidden layer is required.");
	}

	// -------------------------------------------------------------------------
	// General forward pass for N hidden layers.
	//
	// kernels[0..N-1]: hidden layer weight matrices, col-major (num_nodes x prev)
	// kernels[N]:      output layer weight matrix,   col-major (num_output x num_nodes)
	// biases[0..N-1]:  hidden layer biases  (num_nodes,)
	// biases[N]:       output layer bias    (num_output,)
	// input:           x  (num_input,)
	// output:          y  (num_output,)  — must be pre-allocated by caller
	// -------------------------------------------------------------------------
	void forward(
		const std::vector<const double*>& kernels,
		const std::vector<const double*>& biases,
		const double* input,
		double* output) const
	{
		using Map    = Eigen::Map<const Eigen::MatrixXd>;
		using MapVec = Eigen::Map<const Eigen::VectorXd>;

		const unsigned int N = number_of_layers;

		// First hidden layer: h[0] = ELU(W[0] * x + b[0])
		ws_h[0] = Map(kernels[0], number_of_nodes, number_of_input)
		          * MapVec(input, number_of_input)
		          + MapVec(biases[0], number_of_nodes);
		elu_inplace(ws_h[0]);

		// Remaining hidden layers: h[l] = ELU(W[l] * h[l-1] + b[l])
		for (unsigned int l = 1; l < N; ++l)
		{
			ws_h[l] = Map(kernels[l], number_of_nodes, number_of_nodes)
			          * ws_h[l - 1]
			          + MapVec(biases[l], number_of_nodes);
			elu_inplace(ws_h[l]);
		}

		// Output layer (no activation): y = W[N] * h[N-1] + b[N]
		Eigen::Map<Eigen::VectorXd>(output, number_of_output) =
		    Map(kernels[N], number_of_output, number_of_nodes) * ws_h[N - 1]
		    + MapVec(biases[N], number_of_output);
	}

	// -------------------------------------------------------------------------
	// General Jacobian (gradient dy/dx) for N hidden layers.
	//
	// kernels / biases: same layout as forward().
	//   biases[N] (output layer) is accepted but not used — keeps the
	//   interface symmetric with forward() so callers can reuse the same arrays.
	//
	// For scalar output (num_output=1) this gives a (num_input,) gradient vector.
	// -------------------------------------------------------------------------
	void jacobian(
		const std::vector<const double*>& kernels,
		const std::vector<const double*>& biases,
		const double* input,
		double* grad_out) const
	{
		using Map    = Eigen::Map<const Eigen::MatrixXd>;
		using MapVec = Eigen::Map<const Eigen::VectorXd>;

		const unsigned int N = number_of_layers;

		// Forward pass: collect pre-activations (ws_z) and post-activations (ws_h)
		ws_z[0] = Map(kernels[0], number_of_nodes, number_of_input)
		          * MapVec(input, number_of_input)
		          + MapVec(biases[0], number_of_nodes);
		ws_h[0] = ws_z[0];
		elu_inplace(ws_h[0]);

		for (unsigned int l = 1; l < N; ++l)
		{
			ws_z[l] = Map(kernels[l], number_of_nodes, number_of_nodes)
			          * ws_h[l - 1]
			          + MapVec(biases[l], number_of_nodes);
			ws_h[l] = ws_z[l];
			elu_inplace(ws_h[l]);
		}

		// Backward pass
		// Initial delta: W_out^T * 1  elementwise  ELU'(z[N-1])
		ws_delta_cur = Map(kernels[N], number_of_output, number_of_nodes).transpose()
		               * Eigen::VectorXd::Ones(number_of_output);
		elu_backward_inplace(ws_delta_cur, ws_z[N - 1]);

		// Propagate delta back through hidden layers N-2 down to 0
		for (int l = static_cast<int>(N) - 2; l >= 0; --l)
		{
			ws_delta_nxt = Map(kernels[l + 1], number_of_nodes, number_of_nodes).transpose()
			               * ws_delta_cur;
			elu_backward_inplace(ws_delta_nxt, ws_z[l]);
			ws_delta_cur.swap(ws_delta_nxt);
		}

		// Input-space gradient: grad = W_0^T * delta
		Eigen::Map<Eigen::VectorXd>(grad_out, number_of_input) =
		    Map(kernels[0], number_of_nodes, number_of_input).transpose() * ws_delta_cur;
	}

private:
	unsigned int number_of_layers;
	unsigned int number_of_nodes;
	unsigned int number_of_input;
	unsigned int number_of_output;

	// Pre-allocated workspace (mutable so forward/backward can be called on const objects)
	mutable std::vector<Eigen::VectorXd> ws_z;        // pre-activations,  ws_z[l] for hidden layer l
	mutable std::vector<Eigen::VectorXd> ws_h;        // post-activations, ws_h[l] for hidden layer l
	mutable Eigen::VectorXd             ws_delta_cur; // current backprop delta
	mutable Eigen::VectorXd             ws_delta_nxt; // scratch buffer for backprop swap

	// -------------------------------------------------------------------------
	// ELU activation applied in-place.
	// ELU(x) = x          if x >= 0
	//        = exp(x) - 1  if x <  0
	// -------------------------------------------------------------------------
	static void elu_inplace(Eigen::VectorXd& v)
	{
		for (int i = 0; i < v.size(); ++i)
			if (v[i] < 0.0)
				v[i] = std::exp(v[i]) - 1.0;
	}

	// -------------------------------------------------------------------------
	// Multiply gradient vector element-wise by ELU'(z_pre).
	// ELU'(z) = 1       if z >= 0
	//         = exp(z)   if z <  0
	//
	// grad:   current backprop gradient (modified in-place)
	// z_pre:  pre-activation values at the same layer
	// -------------------------------------------------------------------------
	static void elu_backward_inplace(Eigen::VectorXd& grad, const Eigen::VectorXd& z_pre)
	{
		for (int i = 0; i < grad.size(); ++i)
			if (z_pre[i] < 0.0)
				grad[i] *= std::exp(z_pre[i]);
		// if z_pre[i] >= 0: ELU'(z) = 1, no change needed
	}
};
