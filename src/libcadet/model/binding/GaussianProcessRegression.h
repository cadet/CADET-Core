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

namespace GP
{
	class ClassGPR
	{
	public:
		ClassGPR(unsigned int m, unsigned int n, unsigned int k,
			double mlp_weight_var, double mlp_bias_var, double mlp_var, double lin_var,
			double rbf_var, double rbf_ls, double gaussian_var, const std::string& kernel_name)
			: kernelType(KernelType::Unknown)
			, mlp_weight_variance(mlp_weight_var)
			, mlp_bias_variance(mlp_bias_var)
			, mlp_variance(mlp_var)
			, gaussian_variance(gaussian_var)
			, rbf_variance(rbf_var)
			, rbf_lengthscale(rbf_ls)   // Expected to be ls^2 as exported from Python
			, linear_variance(lin_var)
			, M(m)
			, N(n)
			, K(k)
			// Precompute kernel constants
			, rbf_coeff(-1.0 / (2.0 * rbf_ls))
			, mlp_coeff(mlp_var * 2.0 / 3.14159265358979323846)
			, mlp_denom_const(mlp_bias_var + 1.0)
			// Allocate workspace buffers (sized for M x M — training kernel matrix)
			, workspace_sq(static_cast<std::size_t>(m) * m)
			, workspace_rbf(static_cast<std::size_t>(m) * m)
			, workspace_lin(static_cast<std::size_t>(m) * m)
			, workspace_mlp(static_cast<std::size_t>(m) * m)
			, workspace_k_star(m)
			, workspace_norms_x(m)
			, workspace_norms_y(m)
		{
			if (kernel_name == "MLP")
				kernelType = KernelType::MLP;
			else if (kernel_name == "RBF")
				kernelType = KernelType::RBF;
			else if (kernel_name == "RBF_Linear")
				kernelType = KernelType::RBF_Linear;
			else if (kernel_name == "MLP_Linear")
				kernelType = KernelType::MLP_Linear;
			else
				throw std::invalid_argument("GPR: unsupported kernel: " + kernel_name);
		}

		// -----------------------------------------------------------------------
		// Solve (K_train + gaussian_variance * I) * alpha = Y_train via Cholesky.
		// Uses column-major Eigen::MatrixXd — required for Eigen::LLT correctness.
		// Does NOT modify kernel_train (unlike original MKL dposv which destroyed it).
		// -----------------------------------------------------------------------
		void kernel_inv_y(const double* Y_train, const double* kernel_train,
			double* chol_solution) const
		{
			Eigen::MatrixXd Km(static_cast<int>(M), static_cast<int>(M));
			for (unsigned int i = 0; i < M; ++i)
				for (unsigned int j = 0; j < M; ++j)
					Km(i, j) = kernel_train[i * M + j];

			Km.diagonal().array() += gaussian_variance;

			Eigen::Map<const Eigen::VectorXd> y(Y_train, static_cast<int>(M));

			Eigen::LLT<Eigen::MatrixXd> llt(Km);
			if (llt.info() != Eigen::Success)
				throw std::runtime_error(
					"GPR: kernel matrix is not positive definite. "
					"Check gaussian_variance and training data.");

			Eigen::Map<Eigen::VectorXd>(chol_solution, static_cast<int>(M)) = llt.solve(y);
		}

		// -----------------------------------------------------------------------
		// GPR mean prediction: mu = k(X_train, x*)^T * alpha
		// X_train: (M x K), X_test: (1 x K), chol_solution: (M,)
		// -----------------------------------------------------------------------
		double prediction(const double* X_train, const double* X_test,
			const double* chol_solution) const
		{
			double* k_star = const_cast<double*>(workspace_k_star.data());

			switch (kernelType)
			{
				case KernelType::MLP:
					MLP_kernel(M, 1u, K, X_train, X_test, k_star);
					break;
				case KernelType::RBF:
					RBF_kernel(X_train, X_test, k_star, 1u);
					break;
				case KernelType::RBF_Linear:
					RBF_Linear_Kernel(X_train, X_test, k_star, 1u);
					break;
				case KernelType::MLP_Linear:
					MLP_Linear_Kernel(X_train, X_test, k_star, 1u);
					break;
				default:
					break;
			}

			return Eigen::Map<const Eigen::VectorXd>(k_star, static_cast<int>(M))
				.dot(Eigen::Map<const Eigen::VectorXd>(chol_solution, static_cast<int>(M)));
		}

		// -----------------------------------------------------------------------
		// Compute the full training kernel matrix K(X, X): (M x M)
		// Used once at configuration time.
		// -----------------------------------------------------------------------
		void GPR_kernel(const double* X, const double* Y, double* kernel) const
		{
			switch (kernelType)
			{
				case KernelType::MLP:
					MLP_kernel(M, M, K, X, Y, kernel);
					break;
				case KernelType::RBF:
					RBF_kernel(X, Y, kernel, M);
					break;
				case KernelType::RBF_Linear:
					RBF_Linear_Kernel(X, Y, kernel, M);
					break;
				case KernelType::MLP_Linear:
					MLP_Linear_Kernel(X, Y, kernel, M);
					break;
				default:
					break;
			}
		}

		// -----------------------------------------------------------------------
		// Gradient of GPR mean prediction w.r.t. x* (ND generalisation).
		// Returns a vector of length K: grad[d] = d(mu)/d(x*_d)
		//
		// dmu/dx*_d = sum_i alpha_i * dk(x_i, x*)/dx*_d
		//
		// X_train: (M x K), X_test: (1 x K), chol_solution: (M,)
		// -----------------------------------------------------------------------
		std::vector<double> GPR_derivative(const double* X_train, const double* X_test,
			const double* chol_solution) const
		{
			switch (kernelType)
			{
				case KernelType::MLP:
					return MLP_derivative(X_train, X_test, chol_solution);
				case KernelType::RBF:
					return RBF_derivative(X_train, X_test, chol_solution);
				case KernelType::RBF_Linear:
					return RBF_Linear_derivative(X_train, X_test, chol_solution);
				case KernelType::MLP_Linear:
					return MLP_Linear_derivative(X_train, X_test, chol_solution);
				default:
					return std::vector<double>(K, 0.0);
			}
		}

	private:

		enum class KernelType : unsigned int
		{
			MLP       = 0,
			RBF       = 1,
			RBF_Linear = 2,
			MLP_Linear = 3,
			Unknown   = 4
		};

		// -----------------------------------------------------------------------
		// Squared distance: out[i,j] = ||X[i] - Y[j]||^2
		// X: (train_data_size x K), Y: (test_data_size x K), row-major
		// -----------------------------------------------------------------------
		void sqdist(const double* X, const double* Y, double* out,
			unsigned int train_data_size, unsigned int test_data_size) const
		{
			for (unsigned int i = 0; i < train_data_size; ++i)
			{
				for (unsigned int j = 0; j < test_data_size; ++j)
				{
					double s = 0.0;
					for (unsigned int d = 0; d < K; ++d)
					{
						const double diff = X[i * K + d] - Y[j * K + d];
						s += diff * diff;
					}
					out[i * test_data_size + j] = s;
				}
			}
		}

		// -----------------------------------------------------------------------
		// RBF kernel: k(x,y) = rbf_variance * exp(-||x-y||^2 / (2 * rbf_lengthscale))
		// rbf_lengthscale is ls^2 as exported from scikit-learn.
		// X: (M x K), Y: (test_data_size x K), kernel: (M x test_data_size)
		// -----------------------------------------------------------------------
		void RBF_kernel(const double* X, const double* Y, double* kernel,
			unsigned int test_data_size) const
		{
			double* sq = const_cast<double*>(workspace_sq.data());
			sqdist(X, Y, sq, M, test_data_size);

			for (unsigned int i = 0; i < M; ++i)
				for (unsigned int j = 0; j < test_data_size; ++j)
					kernel[i * test_data_size + j] = rbf_variance
						* std::exp(rbf_coeff * sq[i * test_data_size + j]);
		}

		// -----------------------------------------------------------------------
		// Linear kernel: k(x,y) = linear_variance * (x . y)
		// X: (M x K), Y: (test_data_size x K), kernel: (M x test_data_size)
		// -----------------------------------------------------------------------
		void Linear_kernel(const double* X, const double* Y, double* kernel,
			unsigned int test_data_size) const
		{
			for (unsigned int i = 0; i < M; ++i)
			{
				for (unsigned int j = 0; j < test_data_size; ++j)
				{
					double dot = 0.0;
					for (unsigned int d = 0; d < K; ++d)
						dot += X[i * K + d] * Y[j * K + d];
					kernel[i * test_data_size + j] = linear_variance * dot;
				}
			}
		}

		// -----------------------------------------------------------------------
		// MLP (arc-cosine) kernel
		// k(x,y) = mlp_coeff * asin(num / (denomX * denomY))
		// num    = w * dot(x,y) + b
		// denomX = sqrt(w * ||x||^2 + b + 1),  denomY = sqrt(w * ||y||^2 + b + 1)
		// mlp_coeff = mlp_variance * 2/pi  (precomputed)
		// X: (x_row x x_col), Y: (y_row x x_col), kernel: (x_row x y_row)
		// -----------------------------------------------------------------------
		void MLP_kernel(unsigned int x_row, unsigned int y_row, unsigned int x_col,
			const double* X, const double* Y, double* kernel) const
		{
			double* norms_x = const_cast<double*>(workspace_norms_x.data());
			for (unsigned int i = 0; i < x_row; ++i)
			{
				double xnorm = 0.0;
				for (unsigned int d = 0; d < x_col; ++d)
					xnorm += X[i * x_col + d] * X[i * x_col + d];
				norms_x[i] = std::sqrt(mlp_weight_variance * xnorm + mlp_denom_const);
			}

			double* norms_y = const_cast<double*>(workspace_norms_y.data());
			for (unsigned int j = 0; j < y_row; ++j)
			{
				double ynorm = 0.0;
				for (unsigned int d = 0; d < x_col; ++d)
					ynorm += Y[j * x_col + d] * Y[j * x_col + d];
				norms_y[j] = std::sqrt(mlp_weight_variance * ynorm + mlp_denom_const);
			}

			for (unsigned int i = 0; i < x_row; ++i)
			{
				const double denomX = norms_x[i];
				for (unsigned int j = 0; j < y_row; ++j)
				{
					double dot = 0.0;
					for (unsigned int d = 0; d < x_col; ++d)
						dot += X[i * x_col + d] * Y[j * x_col + d];

					const double num = mlp_weight_variance * dot + mlp_bias_variance;
					const double arg = std::clamp(num / (denomX * norms_y[j]), -1.0, 1.0);
					kernel[i * y_row + j] = mlp_coeff * std::asin(arg);
				}
			}
		}

		// -----------------------------------------------------------------------
		// Composite kernels
		// -----------------------------------------------------------------------
		void RBF_Linear_Kernel(const double* X, const double* Y, double* out,
			unsigned int test_data_size) const
		{
			const std::size_t sz = static_cast<std::size_t>(M) * test_data_size;
			double* rbf = const_cast<double*>(workspace_rbf.data());
			double* lin = const_cast<double*>(workspace_lin.data());
			RBF_kernel(X, Y, rbf, test_data_size);
			Linear_kernel(X, Y, lin, test_data_size);
			for (std::size_t i = 0; i < sz; ++i)
				out[i] = rbf[i] + lin[i];
		}

		void MLP_Linear_Kernel(const double* X, const double* Y, double* out,
			unsigned int test_data_size) const
		{
			const std::size_t sz = static_cast<std::size_t>(M) * test_data_size;
			double* mlp = const_cast<double*>(workspace_mlp.data());
			double* lin = const_cast<double*>(workspace_lin.data());
			MLP_kernel(M, test_data_size, K, X, Y, mlp);
			Linear_kernel(X, Y, lin, test_data_size);
			for (std::size_t i = 0; i < sz; ++i)
				out[i] = mlp[i] + lin[i];
		}

		// -----------------------------------------------------------------------
		// Gradient of RBF kernel w.r.t. x* — ND generalisation.
		//
		// dk(x_i, x*)/dx*_d = -(x*_d - x_i_d) / rbf_lengthscale * k(x_i, x*)
		//
		// Returns grad of length K.
		// -----------------------------------------------------------------------
		std::vector<double> RBF_derivative(const double* X_train, const double* X_test,
			const double* chol_solution) const
		{
			double* k_star = const_cast<double*>(workspace_k_star.data());
			RBF_kernel(X_train, X_test, k_star, 1u);

			const double inv_ls = 1.0 / rbf_lengthscale;
			std::vector<double> grad(K, 0.0);

			for (unsigned int i = 0; i < M; ++i)
			{
				const double a = chol_solution[i] * k_star[i];
				for (unsigned int d = 0; d < K; ++d)
				{
					const double diff = X_test[d] - X_train[i * K + d];
					grad[d] += -(diff * inv_ls) * a;
				}
			}
			return grad;
		}

		// -----------------------------------------------------------------------
		// Gradient of Linear kernel w.r.t. x* — ND generalisation.
		//
		// dk(x_i, x*)/dx*_d = linear_variance * x_i_d
		//
		// Returns grad of length K.
		// -----------------------------------------------------------------------
		std::vector<double> Linear_derivative(const double* X_train,
			const double* chol_solution) const
		{
			std::vector<double> grad(K, 0.0);
			for (unsigned int i = 0; i < M; ++i)
			{
				const double a = linear_variance * chol_solution[i];
				for (unsigned int d = 0; d < K; ++d)
					grad[d] += a * X_train[i * K + d];
			}
			return grad;
		}

		// -----------------------------------------------------------------------
		// Gradient of MLP kernel w.r.t. x* — ND generalisation.
		//
		// u   = num / (g * h),   num = w * dot(x_i, x*) + b
		// g   = sqrt(w * ||x_i||^2 + b + 1)   (independent of x*)
		// h   = sqrt(w * ||x*||^2  + b + 1)
		//
		// d(num)/dx*_d = w * x_i_d
		// d(h)/dx*_d   = w * x*_d / h
		// d(u)/dx*_d   = d(num)/dx*_d / (g*h) - num * d(h)/dx*_d / (g * h^2)
		// dk/dx*_d     = mlp_coeff * d(u)/dx*_d / sqrt(1 - u^2)
		//
		// Returns grad of length K.
		// -----------------------------------------------------------------------
		std::vector<double> MLP_derivative(const double* X_train, const double* X_test,
			const double* chol_solution) const
		{
			// Precompute h and dh/dx*_d (both depend only on x*)
			double tnorm = 0.0;
			for (unsigned int d = 0; d < K; ++d)
				tnorm += X_test[d] * X_test[d];
			const double h = std::sqrt(mlp_weight_variance * tnorm + mlp_denom_const);

			// dh/dx*_d = w * x*_d / h  — stored as a vector
			std::vector<double> dh(K, 0.0);
			if (h > 0.0)
				for (unsigned int d = 0; d < K; ++d)
					dh[d] = mlp_weight_variance * X_test[d] / h;

			std::vector<double> grad(K, 0.0);

			for (unsigned int i = 0; i < M; ++i)
			{
				double xnorm = 0.0;
				double dot   = 0.0;
				for (unsigned int d = 0; d < K; ++d)
				{
					const double xv = X_train[i * K + d];
					xnorm += xv * xv;
					dot   += xv * X_test[d];
				}

				const double g   = std::sqrt(mlp_weight_variance * xnorm + mlp_denom_const);
				const double num = mlp_weight_variance * dot + mlp_bias_variance;
				const double gh  = g * h;

				if (gh <= 0.0)
					continue;

				const double u          = std::clamp(num / gh, -1.0, 1.0);
				const double oneMinusU2 = std::max(1.0e-14, 1.0 - u * u);
				const double scale      = mlp_coeff * chol_solution[i] / std::sqrt(oneMinusU2);

				for (unsigned int d = 0; d < K; ++d)
				{
					// dnum/dx*_d = w * x_i_d
					const double dnum_d = mlp_weight_variance * X_train[i * K + d];
					const double du_d   = dnum_d / gh - num * dh[d] / (g * h * h);
					grad[d] += scale * du_d;
				}
			}

			return grad;
		}

		// -----------------------------------------------------------------------
		// Composite derivatives
		// -----------------------------------------------------------------------
		std::vector<double> RBF_Linear_derivative(const double* X_train, const double* X_test,
			const double* chol_solution) const
		{
			std::vector<double> grad = RBF_derivative(X_train, X_test, chol_solution);
			const std::vector<double> lin = Linear_derivative(X_train, chol_solution);
			for (unsigned int d = 0; d < K; ++d)
				grad[d] += lin[d];
			return grad;
		}

		std::vector<double> MLP_Linear_derivative(const double* X_train, const double* X_test,
			const double* chol_solution) const
		{
			std::vector<double> grad = MLP_derivative(X_train, X_test, chol_solution);
			const std::vector<double> lin = Linear_derivative(X_train, chol_solution);
			for (unsigned int d = 0; d < K; ++d)
				grad[d] += lin[d];
			return grad;
		}

		// Kernel type
		KernelType kernelType;

		// Kernel hyperparameters
		double       mlp_weight_variance;
		double       mlp_bias_variance;
		double       mlp_variance;
		double       gaussian_variance;
		double       rbf_variance;
		double       rbf_lengthscale;   // = ls^2 as exported from Python
		double       linear_variance;

		unsigned int M;   // number of training points
		unsigned int N;   // unused — kept for API compatibility
		unsigned int K;   // input dimensionality

		// Precomputed constants
		double rbf_coeff;        // -1 / (2 * rbf_lengthscale)
		double mlp_coeff;        // mlp_variance * 2/pi
		double mlp_denom_const;  // mlp_bias_variance + 1.0

		// Pre-allocated mutable workspace buffers
		mutable std::vector<double> workspace_sq;       // M x M
		mutable std::vector<double> workspace_rbf;      // M x M
		mutable std::vector<double> workspace_lin;      // M x M
		mutable std::vector<double> workspace_mlp;      // M x M
		mutable std::vector<double> workspace_k_star;   // M
		mutable std::vector<double> workspace_norms_x;  // M
		mutable std::vector<double> workspace_norms_y;  // M
	};

}  // namespace GP
