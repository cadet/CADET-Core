#pragma once

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace GP
{
	class GPR_Class
	{
	public:
		GPR_Class(unsigned int m, unsigned int n, unsigned int k,
			double mlp_weight_var, double mlp_bias_var, double mlp_var, double lin_var,
			double rbf_var, double rbf_ls, double gaussian_var, const std::string& kernel_name)
			: mlp_weight_variance(mlp_weight_var)
			, mlp_bias_variance(mlp_bias_var)
			, mlp_variance(mlp_var)
			, gaussian_variance(gaussian_var)
			, rbf_variance(rbf_var)
			, rbf_lengthscale(rbf_ls)
			, linear_variance(lin_var)
			, kernel_choice(kernel_name)
			, M(m)
			, N(n)
			, K(k)
		{
		}

		void sqdist(const double* X, const double* Y, double* sqdist_t1_t2,
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
					sqdist_t1_t2[i * test_data_size + j] = s;
				}
			}
		}

		void RBF_kernel(const double* X, const double* Y, double* kernel, unsigned int test_data_size) const
		{
			std::vector<double> sqdist_values(static_cast<std::size_t>(M) * test_data_size, 0.0);
			sqdist(X, Y, sqdist_values.data(), M, test_data_size);

			for (unsigned int i = 0; i < M; ++i)
			{
				for (unsigned int j = 0; j < test_data_size; ++j)
				{
					const double exponent = -sqdist_values[i * test_data_size + j] / (2.0 * rbf_lengthscale);
					kernel[i * test_data_size + j] = rbf_variance * std::exp(exponent);
				}
			}
		}

		void Linear_kernel(const double* X, const double* Y, double* kernel, unsigned int test_data_size) const
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

		void MLP_kernel(unsigned int x_row, unsigned int y_col, unsigned int x_col,
			const double* X, const double* Y, double* kernel) const
		{
			constexpr double pi = 3.14159265358979323846;
			for (unsigned int i = 0; i < x_row; ++i)
			{
				for (unsigned int j = 0; j < y_col; ++j)
				{
					double dot = 0.0;
					double xnorm = 0.0;
					double ynorm = 0.0;
					for (unsigned int d = 0; d < x_col; ++d)
					{
						const double xv = X[i * x_col + d];
						const double yv = Y[j * x_col + d];
						dot += xv * yv;
						xnorm += xv * xv;
						ynorm += yv * yv;
					}

					const double numerator = mlp_weight_variance * dot + mlp_bias_variance;
					const double denomX = std::sqrt(mlp_weight_variance * xnorm + mlp_bias_variance + 1.0);
					const double denomY = std::sqrt(mlp_weight_variance * ynorm + mlp_bias_variance + 1.0);
					const double denom = denomX * denomY;
					const double arg = (denom > 0.0) ? std::clamp(numerator / denom, -1.0, 1.0) : 0.0;
					kernel[i * y_col + j] = mlp_variance * (2.0 / pi) * std::asin(arg);
				}
			}
		}

		void kernel_inv_y(const double* Y_train, const double* kernel_train, double* chol_solution) const
		{
			using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
			Eigen::Map<const RowMajorMatrix> kernel_map(kernel_train, M, M);
			Eigen::Map<const Eigen::VectorXd> y_map(Y_train, M);

			RowMajorMatrix kernel_mat = kernel_map;
			kernel_mat.diagonal().array() += gaussian_variance;

			Eigen::LLT<RowMajorMatrix> llt(kernel_mat);
			if (llt.info() != Eigen::Success)
				throw std::runtime_error("Cholesky decomposition failed for GPR kernel matrix");

			const Eigen::VectorXd alpha = llt.solve(y_map);
			if (llt.info() != Eigen::Success)
				throw std::runtime_error("Failed to solve GPR linear system");

			Eigen::Map<Eigen::VectorXd>(chol_solution, M) = alpha;
		}

		void RBF_Linear_Kernel(const double* X, const double* Y, double* out, unsigned int test_data_size) const
		{
			std::vector<double> rbf_kernel(static_cast<std::size_t>(M) * test_data_size, 0.0);
			std::vector<double> linear_kernel(static_cast<std::size_t>(M) * test_data_size, 0.0);
			RBF_kernel(X, Y, rbf_kernel.data(), test_data_size);
			Linear_kernel(X, Y, linear_kernel.data(), test_data_size);
			for (std::size_t i = 0; i < rbf_kernel.size(); ++i)
				out[i] = rbf_kernel[i] + linear_kernel[i];
		}

		void MLP_Linear_Kernel(const double* X, const double* Y, double* out, unsigned int test_data_size) const
		{
			std::vector<double> mlp_kernel(static_cast<std::size_t>(M) * test_data_size, 0.0);
			std::vector<double> linear_kernel(static_cast<std::size_t>(M) * test_data_size, 0.0);
			MLP_kernel(M, test_data_size, K, X, Y, mlp_kernel.data());
			Linear_kernel(X, Y, linear_kernel.data(), test_data_size);
			for (std::size_t i = 0; i < mlp_kernel.size(); ++i)
				out[i] = mlp_kernel[i] + linear_kernel[i];
		}

		double prediction(const double* X_train, const double* X_test, const double* chol_solution) const
		{
			std::vector<double> kernel_train_test(M, 0.0);

			if (kernel_choice == "MLP")
				MLP_kernel(M, 1u, K, X_train, X_test, kernel_train_test.data());
			else if (kernel_choice == "RBF")
				RBF_kernel(X_train, X_test, kernel_train_test.data(), 1u);
			else if (kernel_choice == "RBF_Linear")
				RBF_Linear_Kernel(X_train, X_test, kernel_train_test.data(), 1u);
			else if (kernel_choice == "MLP_Linear")
				MLP_Linear_Kernel(X_train, X_test, kernel_train_test.data(), 1u);
			else
				throw std::invalid_argument("Unsupported kernel selected in GPR prediction");

			const Eigen::Map<const Eigen::VectorXd> k_vec(kernel_train_test.data(), M);
			const Eigen::Map<const Eigen::VectorXd> alpha(chol_solution, M);
			return k_vec.dot(alpha);
		}

		double MLP_derivative(const double* X_train, const double* X_test, const double* chol_solution) const
		{
			constexpr double pi = 3.14159265358979323846;
			double pred_derivative = 0.0;

			for (unsigned int i = 0; i < M; ++i)
			{
				double dot = 0.0;
				double xnorm = 0.0;
				double tnorm = 0.0;
				for (unsigned int d = 0; d < K; ++d)
				{
					const double xv = X_train[i * K + d];
					const double tv = X_test[d];
					dot += xv * tv;
					xnorm += xv * xv;
					tnorm += tv * tv;
				}

				const double a = mlp_weight_variance;
				const double b = mlp_bias_variance;
				const double numerator = a * dot + b;
				const double g = std::sqrt(a * xnorm + b + 1.0);
				const double h = std::sqrt(a * tnorm + b + 1.0);
				const double denom = g * h;
				if (denom <= 0.0)
					continue;

				const double u = std::clamp(numerator / denom, -1.0, 1.0);

				// Derivative with respect to the first input dimension (single-component usage)
				const double x0 = X_train[i * K];
				const double t0 = X_test[0];
				const double hprime = (h > 0.0) ? (a * t0 / h) : 0.0;
				const double du = (a * x0 * denom - numerator * g * hprime) / (denom * denom);
				const double one_minus_u2 = std::max(1e-14, 1.0 - u * u);
				const double dk = mlp_variance * (2.0 / pi) * du / std::sqrt(one_minus_u2);
				pred_derivative += dk * chol_solution[i];
			}

			return pred_derivative;
		}

		double Linear_derivative(const double* X_train, const double* chol_solution) const
		{
			double linear_derivative = 0.0;
			for (unsigned int i = 0; i < M; ++i)
				linear_derivative += linear_variance * X_train[i * K] * chol_solution[i];
			return linear_derivative;
		}

		double RBF_derivative(const double* X_train, const double* X_test, const double* chol_solution) const
		{
			std::vector<double> kernel_train_test(M, 0.0);
			RBF_kernel(X_train, X_test, kernel_train_test.data(), 1u);

			double pred_derivative = 0.0;
			for (unsigned int i = 0; i < M; ++i)
			{
				const double x_tilda_star = X_test[0] - X_train[i * K];
				pred_derivative += -(x_tilda_star / rbf_lengthscale) * kernel_train_test[i] * chol_solution[i];
			}
			return pred_derivative;
		}

		double RBF_Linear_derivative(const double* X_train, const double* X_test, const double* chol_solution) const
		{
			return RBF_derivative(X_train, X_test, chol_solution) + Linear_derivative(X_train, chol_solution);
		}

	private:
		double mlp_weight_variance;
		double mlp_bias_variance;
		double mlp_variance;
		double gaussian_variance;
		double rbf_variance;
		double rbf_lengthscale;
		double linear_variance;
		std::string kernel_choice;
		unsigned int M;
		unsigned int N;
		unsigned int K;
	};
}
