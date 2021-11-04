#pragma once
#define _USE_MATH_DEFINES 
#define M_PI       3.14159265358979323846   // pi

#include <cstdio>
#include <iostream>
#include <cassert>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include "../../../mkl_files/mkl.h"
#include <chrono>
#include <stdexcept>

namespace
{
	namespace GP {
		class GPR_Class
		{
		protected:
			//std::vector<double> X_train, Y_train; //Training features and outputs or training dataset
			double mlp_weigth_variance;
			double mlp_bias_variance;
			double mlp_variance;
			double guassian_variance;
			double rbf_variance;
			double rbf_lengthscale;
			double linear_variance;
			unsigned int test_data_size;
			std::string kernel_choice;
			MKL_INT M, N, K; //M: Length of training dataset, N: length of test data, K: Dimension
		public:

			GPR_Class(unsigned int m, unsigned int n, unsigned int k,
				double mlp_weigth_var, double mlp_weigth_bs,double mlp_var, double lin_var,
				double rbf_var, double rbf_ls, double guass_var, std::string kernel_name) :
				M(m), N(n), K(k), 
				test_data_size(k), mlp_weigth_variance(mlp_weigth_var), mlp_bias_variance(mlp_weigth_bs), mlp_variance(mlp_var),
				rbf_variance(rbf_var), rbf_lengthscale(rbf_ls), linear_variance(lin_var),
				guassian_variance(guass_var), kernel_choice(kernel_name) {}

			//Method to calculate the sum of squared differences
			void sqdist(const double* X,const  double* Y, double* sqdist_t1_t2, unsigned int train_data_size, unsigned int test_data_size)
			{
				std::vector<double> sqrdTrain(train_data_size * K, 0.0);
				std::vector<double> sqrdTrainSum(train_data_size, 0.0);
				std::vector<double> sqrdTest(test_data_size * K, 0.0);
				std::vector<double> sqrdTestSum(test_data_size, 0.0);
				std::vector<double> sqdist1(train_data_size * test_data_size, 0.0);
				double temp = 0.0;
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					train_data_size, test_data_size, K, 1.0, X, K, Y, test_data_size, 1.0, sqdist1.data(), test_data_size);

				for (unsigned int i = 0; i < train_data_size; ++i)
					sqrdTrain[i] = std::pow(X[i], 2);

				for (unsigned int i = 0; i < test_data_size; ++i)
					sqrdTest[i] = std::pow(Y[i], 2);

				// Adding multi-dimensional support
				if (K > 1)
				{
					for (unsigned int i = 0; i < sqrdTrainSum.size(); ++i)
					{
						temp = sqrdTrain[i];
						for (unsigned int p = 1; p < K; ++p)
							temp = temp + sqrdTrain[i + p * train_data_size];
						sqrdTrainSum[i] = temp;
					}

					temp = 0;
					for (unsigned int i = 0; i < sqrdTestSum.size(); ++i)
					{
						temp = sqrdTest[i];
						for (unsigned int p = 1; p < K; ++p)
							temp = temp + sqrdTrain[i + p * test_data_size];
						sqrdTestSum[i] = temp;
					}
				}
				else
				{
					sqrdTestSum = sqrdTest;
					sqrdTrainSum = sqrdTrain;
				}

				for (unsigned int i = 0; i < test_data_size; ++i)
				{
					for (unsigned int j = 0; j < train_data_size; ++j)
					{
						sqdist_t1_t2[j + train_data_size * i] = sqrdTrainSum[j] + sqrdTestSum[i] - 2 * sqdist1[j + train_data_size * i];
					}
				}
			}

			// Method to compute radial basis function kernel
			void RBF_kernel(const double* X,const double* Y, double* Kernel, unsigned int test_data_size)
			{

				std::vector<double> sqdist1(M * test_data_size, 0.0);

				sqdist(X, Y, sqdist1.data(), M, test_data_size);

				for (unsigned int i = 0; i < test_data_size; ++i)
				{
					for (unsigned int j = 0; j < M; ++j)
					{
						Kernel[j + M * i] = rbf_variance * std::exp(-sqdist1[j + M * i] / (2 * rbf_lengthscale));
					}
				}
			}

			// Method to compute linear kernel
			void Linear_kernel(const double* X,const  double* Y, double* Kernel,
				unsigned int test_data_size)
			{
				double temp = 0;

				//Linear Kernel : K(X,Y) = \sum_i=0^dimensions (\sigma_i * X_i*Y_i)
				for (unsigned int i = 0; i < M; ++i)
				{
					for (unsigned int j = 0; j < test_data_size; ++j)
					{
						if (K > 1)
						{
							temp = linear_variance * X[i] * Y[j];
							for (unsigned int p = 1; p < K; ++p)
								temp = temp + linear_variance * X[i + p * M] * Y[j + p * test_data_size];
							Kernel[i + M * j] = temp;
						}
						else
							Kernel[i + M * j] = linear_variance * X[i] * Y[j];
					}
				}
			}


			// Method to crete Multi-layer perceptron kernel
			void MLP_kernel(unsigned int x_row, unsigned int y_col, unsigned int x_col,
				const double* X,const double* Y, double* Kernel)
			{
				/*Parameters:


				*/

				double* numerator = new double[x_row * y_col]; // \sigma_w^2*X*XtT + \sigma_b^2
				std::vector<double> denom1(x_row * x_col, 0.0); // \sqrt(\sigma_w^2*X*XT + \sigma_b^2 + 1)
				std::vector<double> denom1_MD(x_row * 1, 0.0); // For multi-dimensional support
				std::vector<double> denom2(x_col * y_col, 0.0);//  \sqrt(\sigma_w ^ 2 * Xt * XtT + \sigma_b ^ 2 + 1)
				double* fraction = new double[x_row * y_col]; // numerator/denom1*denom2

				for (unsigned int i = 0; i < (x_row * y_col); i++)
				{
					if (mlp_bias_variance != 0) {
						numerator[i] = 1.0;
						fraction[i] = 1.0;
					}
					else { numerator[i] = 0.0;  fraction[i] = 0.0; }
				}

				//Calculating the numerator
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					x_row, y_col, x_col, mlp_weigth_variance, X, x_col, Y, y_col, mlp_bias_variance, numerator, y_col);

				for (unsigned int i = 0; i < x_row * x_col; ++i)
				{
					denom1[i] = std::sqrt((mlp_weigth_variance * X[i] * X[i] + mlp_bias_variance + 1));
				}

				for (unsigned int i = 0; i < x_col * y_col; ++i)
				{
					denom2[i] = std::sqrt((mlp_weigth_variance * Y[i] * Y[i] + mlp_bias_variance + 1));
				}
				if (x_col > 1) //TODO later on
				{
					for (unsigned int i = 0; i < denom1_MD.size(); ++i)
					{
						denom1_MD[i] = denom1[i] + denom1[i + x_row];
					}
				}
				else
					denom1_MD = denom1;
				//Calculating the fraction inside the asin(fraction)
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					x_row, y_col, x_col, 1.0, denom1_MD.data(), x_col, denom2.data(), y_col, 0.0, fraction, y_col);

				for (unsigned i = 0; i < x_row * y_col; ++i)
				{
					Kernel[i] = mlp_variance * (2.0 / M_PI) * std::asin(numerator[i] / fraction[i]);
				}
				delete[] numerator;
				delete[] fraction;
			}

			// Method to calculate the \alpha = inv(K + \alpna_n*I)*y
			void kernel_inv_y(const double* X_train,const double* Y_train,const double* kernel_train, double* chol_solution)
			{
				MKL_INT info, n = M, nrhs = K, lda = M, ldb = K;
				std::vector<double> Kernel_upper_tri(M * M, 0.0);
				std::vector<double> Kernel_train_test(M * K, 0.0);

				for (unsigned i = 0; i < M * K; ++i)
					chol_solution[i] = Y_train[i];

				// Storing the Kernel matrix in upper traingular matrix form
				for (unsigned int i = 0; i < M; i++)
				{
					for (unsigned int j = 0; j < M; j++)
					{
						if (i > j)
							Kernel_upper_tri[j + i * M] = 0.0;
						if (i == j)
							Kernel_upper_tri[j + i * M] = kernel_train[j + i * M] + guassian_variance;
						else
							Kernel_upper_tri[j + i * M] = kernel_train[j + i * M];
					}
				}
				info = LAPACKE_dposv(LAPACK_ROW_MAJOR, 'U', n, nrhs, Kernel_upper_tri.data(), lda, chol_solution, ldb);
				if (info > 0) {
					std::cout << "The leading minor of order " << info << " is not positive ";
					std::cout << "definite in the kernel matrix;\nthe GPR based prediction cannot cannot not be computed.\n";
					exit(1);
				}
			}

			// Method to compute RBF+Linear Kernel
			void RBF_Linear_Kernel(const double* X,const double* Y, double* RBF_Linear_kernel,
				unsigned int test_data_size)
			{
				std::vector<double> RBF_Kernel(M * test_data_size, 0.0);
				std::vector<double> Linear_Kernel(M * test_data_size, 0.0);

				RBF_kernel(X, Y, RBF_Kernel.data(), test_data_size);
				Linear_kernel(X, Y, Linear_Kernel.data(), test_data_size);

				for (unsigned int i = 0; i < test_data_size; ++i)
				{
					for (unsigned int j = 0; j < M; ++j)
						RBF_Linear_kernel[j + test_data_size * i] = RBF_Kernel[j + test_data_size * i] + Linear_Kernel[j + test_data_size * i];
				}
			}

			// Method to compute MLP+Linear Kernel
			void MLP_Linear_Kernel(const double* X,const double* Y, double* MLP_Linear_kernel,
				unsigned int test_data_size)
			{
				std::vector<double> MLP_Kernel(M * test_data_size, 0.0);
				std::vector<double> Linear_Kernel(M * test_data_size, 0.0);

				MLP_kernel(M, test_data_size, K, X, Y, MLP_Kernel.data());
				Linear_kernel(X, Y, Linear_Kernel.data(), test_data_size);

				for (unsigned int i = 0; i < test_data_size; ++i)
				{
					for (unsigned int j = 0; j < M; ++j)
						MLP_Linear_kernel[j + test_data_size * i] = MLP_Kernel[j + test_data_size * i] + Linear_Kernel[j + test_data_size * i];
				}
			}

			// Method to make the prediction on unseen dataset using GPR
			double* prediction(const double* X_train,const double* Y_train,double* X_test,const double* kernel_train,const double* chol_solution,
								std::string kernal_name)
			{
				/*Parameters:
				X_train: Inputs of training dataset (pore phase concentration)
				Y_train: Outputs of training dataset (solid-phase concentration)
				X_test: Test input (evaulation of any pore phase concentration value provided by CADET solver)
				kernel_train: Kernel matrix for training dataset [dim]: number of training data points x number of training data points
				*/
				
				unsigned int data_size = M;
				unsigned int dimension = K;
				std::vector<double> Kernel_train_test(M * K, 0.0);
				double* predicted = new double[K];

				if (kernal_name == "MLP")
					MLP_kernel(M, K, K, X_train, X_test, Kernel_train_test.data());
				else if (kernal_name == "RBF")
					RBF_kernel(X_train, X_test, Kernel_train_test.data(), 1);
				else if (kernal_name == "RBF_Linear")
					RBF_Linear_Kernel(X_train, X_test, Kernel_train_test.data(), 1);
				else if (kernal_name == "MLP_Linear")
					MLP_Linear_Kernel(X_train, X_test, Kernel_train_test.data(), 1);
				else
				{
					std::cout << "Select a valid kernel. 1) MLP 2) RBF or 3) RBF_Linear \n";
					exit(0);
				}

				
				//Calculation of prediction
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					dimension, dimension, data_size, 1.0, Kernel_train_test.data(), data_size, chol_solution, dimension, 0.0, predicted, dimension);

				//std::cout << "Prediction is : \n" << predicted[0];
				return predicted;
			}

			// Function to determine the derivative of GPR based on MLP kernel
			void MLP_derivative(const double* X_train,const double* X_test,const double* chol_solution, double* pred_derivative)
			{
				unsigned int test_set_dim = K;
				//chol_solution-> inv(K+\sigmaI)y

				std::vector<double> derivative_numerator(M * K, 0.0);
				std::vector<double> derivative_denominator_1(M * K, 0.0);
				std::vector<double> derivative_denominator_2(test_set_dim * K, 0.0);
				std::vector<double> derivative_denominator(M * K, 0.0);
				std::vector<double> fraction(M * K, 0.0);
				std::vector<double> Kernel_Xt_Xs(M * K, 0.0);
				std::vector<double> Kernel_alpha(M * K, 0.0);


				MLP_kernel(M, K, K, X_train, X_test, Kernel_Xt_Xs.data());
				// Computing derivative_numerator, derivative_denominator_1
				for (unsigned int i = 0; i < M; ++i)
				{
					for (unsigned int j = 0; j < K; ++j)
					{
						derivative_numerator[i + M * j] = mlp_weigth_variance * mlp_bias_variance * X_train[i + M * j] + mlp_weigth_variance * X_train[i + M * j] - mlp_bias_variance * mlp_weigth_variance * X_test[j];
						derivative_denominator_1[i + M * j] = std::sqrt(mlp_weigth_variance * X_train[i + M * j] * X_train[i + M * j] + 1 + mlp_bias_variance);
						Kernel_Xt_Xs[i + M * j] = Kernel_Xt_Xs[i + M * j] / (mlp_variance * 2 / M_PI);
						Kernel_alpha[i + M * j] = (1 / std::cos(Kernel_Xt_Xs[i + M * j])) * chol_solution[i + M * j];
					}
				}

				// Computing derivative_denominator_2
				for (unsigned int i = 0; i < test_set_dim; ++i)
				{
					for (unsigned int j = 0; j < K; ++j)
					{
						derivative_denominator_2[i + test_set_dim * j] = std::pow(mlp_weigth_variance * X_test[i + test_set_dim * j] * X_test[i + test_set_dim * j] + mlp_bias_variance + 1, 1.5);
					}
				}

				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					M, K, K, 1.0, derivative_denominator_1.data(), K, derivative_denominator_2.data(), K, 0.0, derivative_denominator.data(), K);

				for (unsigned int i = 0; i < M; ++i)
				{
					for (unsigned int j = 0; j < K; ++j)
					{
						fraction[i + M * j] = mlp_variance * (2 / M_PI) * derivative_numerator[i + M * j] / derivative_denominator[i + M * j];
					}
				}
				//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					//x_row, y_col, x_col, 1.0, denom1_MD.data(), x_col, denom2.data(), y_col, 0.0, fraction, y_col);
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					K, K, M, 1.0, fraction.data(), M, Kernel_alpha.data(), K, 0.0, pred_derivative, K);
				//std::cout << "Deriv: \n";
				//std::cout << pred_derivative[0] << "\n";
			}

			// Method to compute the derivative of GPR based on Linear kernel
			void Linear_derivative(const double* X_train,const double* X_test,const double* chol_solution,
				double* linear_derivative)
			{
				// For the moment 1-D case only
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					K, K, M, linear_variance, X_train, M, chol_solution, K, 0.0, linear_derivative, K);
			}

			// Method to compute the derivative of GPR based on RBF kernel
			void RBF_derivative(const double* X_train,const double* X_test,const double* chol_solution, double* pred_derivative,
				unsigned int test_data_size)
			{
				std::vector<double> kernel_train_test(M * test_data_size, 0.0);
				std::vector<double> kernel_alpha(M * K, 0.0);
				std::vector<double>X_tilda_star(M * K, 0.0);

				RBF_kernel(X_train, X_test, kernel_train_test.data(), test_data_size);
				for (unsigned int i = 0; i < K; ++i)
				{
					for (unsigned int j = 0; j < M; ++j)
					{
						X_tilda_star[j + M * i] = X_test[i] - X_train[j + M * i];
						kernel_alpha[j + M * i] = kernel_train_test[j] * chol_solution[j + M * i];
					}
				}
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
					K, K, M, -rbf_lengthscale, X_tilda_star.data(), M, kernel_alpha.data(), K, 0.0, pred_derivative, K);
				std::cout << "\nDerivative is: \n";
				std::cout << pred_derivative[0] << "\n";
			}

			//Method to compute the derivative of GPR based on RBF+Linear kernel
			void RBF_Linear_derivative(const double* X_train,const double* X_test,const double* chol_solution, double* pred_derivative,
				unsigned int test_data_size)
			{
				std::vector<double> RBF_derivative(K * K, 0.0);
				std::vector<double> Linear_derivative(K * K, 0.0);

				this->RBF_derivative(X_train, X_test, chol_solution, RBF_derivative.data(), test_data_size);
				this->Linear_derivative(X_train, X_test, chol_solution, Linear_derivative.data());

				for (unsigned int i = 0; i < K * K; ++i)
					pred_derivative[i] = RBF_derivative[i] + Linear_derivative[i];
			}
		};
	}
}
//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	//x_row, y_col, x_col, 1.0, denom1_MD.data(), x_col, denom2.data(), y_col, 0.0, fraction, y_col);

/*double* alpha_kernal = new double[data_size];
				MKL_INT info, n = data_size, nrhs = dimension, lda = data_size, ldb = dimension;
				for (unsigned int i = 0; i < data_size; ++i)
					alpha_kernal[i] = 0.0;
				std::vector<double> Kernel_upper_tri(data_size * data_size, 0.0);
				std::vector<double> Kernel_train_test(data_size * dimension, 0.0);
				std::vector<double> chol_solution(data_size * dimension, 0.0);
				//double* chol_solution = new double[data_size * dimension];
				double* predicted = new double[dimension];

				for (unsigned i = 0; i < dimension; ++i)
					predicted[i] = 0.0;

				for (unsigned i = 0; i < data_size * dimension; ++i)
					chol_solution[i] = Y_train[i];

				// Storing the Kernel matrix in upper traingular matrix form
				for (unsigned int i = 0; i < data_size; i++)
				{
					for (unsigned int j = 0; j < data_size; j++)
					{
						if (i > j)
							Kernel_upper_tri[j + i * data_size] = 0.0;
						if (i == j)
							Kernel_upper_tri[j + i * data_size] = kernel_train[j + i * data_size] + guassian_variance;
						else
							Kernel_upper_tri[j + i * data_size] = kernel_train[j + i * data_size];
					}
				}
				info = LAPACKE_dposv(LAPACK_ROW_MAJOR, 'U', n, nrhs, Kernel_upper_tri.data(), lda, chol_solution.data(), ldb);
				if (info > 0) {
					std::cout<<"The leading minor of order "<<info<<" is not positive ";
					std::cout<<"definite in the kernel matrix;\nthe GPR based prediction cannot cannot not be computed.\n";
					exit(1);
				}
				std::cout << "\n";*/#pragma once
