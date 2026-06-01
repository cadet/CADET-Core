#pragma once

#include <Eigen/Dense>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

// =============================================================================
// Class_ANN — Feedforward neural network with ELU activations.
//
// Architecture (single-component, scalar output):
//   1-hidden-layer:  y = W2 * ELU(W1 * x + b1) + b2
//   2-hidden-layer:  y = W3 * ELU(W2 * ELU(W1 * x + b1) + b2) + b3
//
// Weight storage convention (matches original training export):
//   All weight matrices are stored COLUMN-MAJOR (Fortran order).
//   A weight matrix W of shape (n_out x n_in) is stored such that
//   W(row=i, col=j) = kernel[i + j * n_out].
//   This means Eigen::Map<MatrixXd>(kernel.data(), n_out, n_in) gives W directly.
//
// MKL has been removed. All linear algebra uses Eigen.
// All workspace buffers are pre-allocated at construction — no heap allocation
// inside forward/backward methods.
// =============================================================================

class Class_ANN
{
public:
	// -------------------------------------------------------------------------
	// Constructor — pre-allocates all workspace buffers.
	// num_layer: number of hidden layers (1 or 2)
	// num_nodes: nodes per hidden layer (same for all hidden layers)
	// num_inputs: input dimensionality
	// num_outputs: output dimensionality (1 for single-component isotherm)
	// -------------------------------------------------------------------------
	Class_ANN(unsigned int num_layer, unsigned int num_nodes,
	          unsigned int num_inputs, unsigned int num_outputs)
		: number_of_layers(num_layer)
		, number_of_nodes(num_nodes)
		, number_of_input(num_inputs)
		, number_of_output(num_outputs)
		, ws_z1(num_nodes)
		, ws_z1_pre(num_nodes)
		, ws_z2(num_nodes)
		, ws_delta(num_nodes)
		, ws_delta2(num_nodes)
	{
		if (num_layer != 1 && num_layer != 2)
			throw std::invalid_argument(
				"Class_ANN: only 1 or 2 hidden layers are supported.");
	}

	// -------------------------------------------------------------------------
	// Forward pass — 1-hidden-layer network.
	//
	// y = W2 * ELU(W1 * x + b1) + b2
	//
	// kernel0: W1 col-major (num_nodes x num_input)
	// bias0:   b1 (num_nodes,)
	// kernel1: W2 col-major (num_output x num_nodes)
	// bias1:   b2 (num_output,)
	// input:   x  (num_input,)
	// output:  y  (num_output,)   -- must be pre-allocated by caller
	// -------------------------------------------------------------------------
	void forward_single_layer(
		const double* input,
		const double* kernel0, const double* bias0,
		const double* kernel1, const double* bias1,
		double* output) const
	{
		using Map    = Eigen::Map<const Eigen::MatrixXd>;
		using MapVec = Eigen::Map<const Eigen::VectorXd>;

		// z1 = W1 * x + b1
		ws_z1 = Map(kernel0, number_of_nodes, number_of_input)
		        * MapVec(input, number_of_input)
		        + MapVec(bias0, number_of_nodes);

		elu_inplace(ws_z1);   // h1 = ELU(z1), stored in ws_z1

		// y = W2 * h1 + b2
		Eigen::Map<Eigen::VectorXd>(output, number_of_output) =
		    Map(kernel1, number_of_output, number_of_nodes) * ws_z1
		    + MapVec(bias1, number_of_output);
	}

	// -------------------------------------------------------------------------
	// Forward pass — 2-hidden-layer network.
	//
	// y = W3 * ELU(W2 * ELU(W1 * x + b1) + b2) + b3
	//
	// kernel0: W1 col-major (num_nodes x num_input)
	// bias0:   b1 (num_nodes,)
	// kernel1: W2 col-major (num_nodes x num_nodes)
	// bias1:   b2 (num_nodes,)
	// kernel2: W3 col-major (num_output x num_nodes)
	// bias2:   b3 (num_output,)
	// input:   x  (num_input,)
	// output:  y  (num_output,)   -- must be pre-allocated by caller
	// -------------------------------------------------------------------------
	void forward_two_layers(
		const double* input,
		const double* kernel0, const double* bias0,
		const double* kernel1, const double* bias1,
		const double* kernel2, const double* bias2,
		double* output) const
	{
		using Map    = Eigen::Map<const Eigen::MatrixXd>;
		using MapVec = Eigen::Map<const Eigen::VectorXd>;

		// z1 = W1 * x + b1,  h1 = ELU(z1)
		ws_z1 = Map(kernel0, number_of_nodes, number_of_input)
		        * MapVec(input, number_of_input)
		        + MapVec(bias0, number_of_nodes);
		elu_inplace(ws_z1);

		// z2 = W2 * h1 + b2,  h2 = ELU(z2)
		ws_z2 = Map(kernel1, number_of_nodes, number_of_nodes) * ws_z1
		        + MapVec(bias1, number_of_nodes);
		elu_inplace(ws_z2);

		// y = W3 * h2 + b3  (no activation on output layer)
		Eigen::Map<Eigen::VectorXd>(output, number_of_output) =
		    Map(kernel2, number_of_output, number_of_nodes) * ws_z2
		    + MapVec(bias2, number_of_output);
	}

	// -------------------------------------------------------------------------
	// Jacobian — 1-hidden-layer network.
	//
	// dy/dx = W1^T * diag(ELU'(z1)) * W2^T
	// For scalar output (num_output=1) this gives a (num_input,) gradient vector.
	//
	// The pre-activation z1 is computed internally and stored in ws_z1_pre.
	// -------------------------------------------------------------------------
	void jacobian_single_layer(
		const double* input,
		const double* kernel0, const double* bias0,
		const double* kernel1,
		double* grad_out) const
	{
		using Map    = Eigen::Map<const Eigen::MatrixXd>;
		using MapVec = Eigen::Map<const Eigen::VectorXd>;

		// Forward: z1 = W1 * x + b1  (keep pre-activation for ELU')
		ws_z1_pre = Map(kernel0, number_of_nodes, number_of_input)
		            * MapVec(input, number_of_input)
		            + MapVec(bias0, number_of_nodes);

		// delta = W2^T * 1  (output layer weight transposed, scalar output)
		// For num_output=1: W2 has shape (1 x num_nodes), W2^T is (num_nodes x 1)
		ws_delta = Map(kernel1, number_of_output, number_of_nodes).transpose()
		           * Eigen::VectorXd::Ones(number_of_output);

		// Multiply element-wise by ELU'(z1)
		elu_backward_inplace(ws_delta, ws_z1_pre);

		// grad = W1^T * delta
		Eigen::Map<Eigen::VectorXd>(grad_out, number_of_input) =
		    Map(kernel0, number_of_nodes, number_of_input).transpose() * ws_delta;
	}

	// -------------------------------------------------------------------------
	// Jacobian — 2-hidden-layer network.
	//
	// dy/dx = W1^T * diag(ELU'(z1)) * W2^T * diag(ELU'(z2)) * W3^T
	//
	// Backprop steps:
	//   delta  = W3^T * 1              [elementwise * ELU'(z2)]
	//   delta2 = W2^T * delta          [elementwise * ELU'(z1)]
	//   grad   = W1^T * delta2
	// -------------------------------------------------------------------------
	void jacobian_two_layers(
		const double* input,
		const double* kernel0, const double* bias0,
		const double* kernel1, const double* bias1,
		const double* kernel2,
		double* grad_out) const
	{
		using Map    = Eigen::Map<const Eigen::MatrixXd>;
		using MapVec = Eigen::Map<const Eigen::VectorXd>;

		// Forward pass to get pre-activations z1 and z2
		// z1 = W1 * x + b1
		ws_z1_pre = Map(kernel0, number_of_nodes, number_of_input)
		            * MapVec(input, number_of_input)
		            + MapVec(bias0, number_of_nodes);

		// h1 = ELU(z1)
		ws_z1 = ws_z1_pre;
		elu_inplace(ws_z1);

		// z2 = W2 * h1 + b2  (keep pre-activation for ELU')
		ws_z2 = Map(kernel1, number_of_nodes, number_of_nodes) * ws_z1
		        + MapVec(bias1, number_of_nodes);

		// Backward pass
		// delta = W3^T (scalar output => W3 shape (1 x n_h), W3^T shape (n_h,))
		ws_delta = Map(kernel2, number_of_output, number_of_nodes).transpose()
		           * Eigen::VectorXd::Ones(number_of_output);

		// delta .*= ELU'(z2)
		elu_backward_inplace(ws_delta, ws_z2);

		// delta2 = W2^T * delta  .*  ELU'(z1)
		ws_delta2 = Map(kernel1, number_of_nodes, number_of_nodes).transpose() * ws_delta;
		elu_backward_inplace(ws_delta2, ws_z1_pre);

		// grad = W1^T * delta2
		Eigen::Map<Eigen::VectorXd>(grad_out, number_of_input) =
		    Map(kernel0, number_of_nodes, number_of_input).transpose() * ws_delta2;
	}

private:
	unsigned int number_of_layers;
	unsigned int number_of_nodes;
	unsigned int number_of_input;
	unsigned int number_of_output;

	// Pre-allocated workspace (mutable so forward/backward can be called on const objects)
	mutable Eigen::VectorXd ws_z1;       // pre- or post-activation layer 1
	mutable Eigen::VectorXd ws_z1_pre;   // pre-activation layer 1 (kept for backprop)
	mutable Eigen::VectorXd ws_z2;       // pre-activation layer 2 (kept for backprop)
	mutable Eigen::VectorXd ws_delta;    // backprop gradient buffer layer 2
	mutable Eigen::VectorXd ws_delta2;   // backprop gradient buffer layer 1

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
	// ELU'(z) = 1          if z >= 0
	//         = exp(z)      if z <  0
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
