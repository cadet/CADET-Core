// =============================================================================
//  CADET
//  
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "Weno.hpp"

#include <cmath>

namespace
{

	/**
	 * @brief Integral of the weighted monomial @f$ (a + b \xi)^m \xi^n @f$ over @f$ [\xi_L, \xi_R] @f$ for @f$ m \le 2 @f$
	 */
	double weightMomentIntegral(double xiL, double xiR, double a, double b, int m, int n)
	{
		// Binomial expansion of (a + b*xi)^m
		double coef[3] = { 1.0, 0.0, 0.0 };
		if (m == 1)
		{
			coef[0] = a;
			coef[1] = b;
		}
		else if (m == 2)
		{
			coef[0] = a * a;
			coef[1] = 2.0 * a * b;
			coef[2] = b * b;
		}

		double res = 0.0;
		for (int k = 0; k <= m; ++k)
		{
			const int p = n + k + 1;
			res += coef[k] * (std::pow(xiR, p) - std::pow(xiL, p)) / static_cast<double>(p);
		}
		return res;
	}

	/**
	 * @brief Solves the dense linear system A x = b in place via Gauss elimination with partial pivoting
	 * @param [in] n System size
	 * @param [in,out] A Row-major matrix, destroyed in the process
	 * @param [in,out] b Right hand side, contains the solution on exit
	 * @return @c true on success, @c false if the matrix is (numerically) singular
	 */
	bool gaussSolve(int n, double* A, double* b)
	{
		for (int col = 0; col < n; ++col)
		{
			int piv = col;
			for (int r = col + 1; r < n; ++r)
			{
				if (std::abs(A[r * n + col]) > std::abs(A[piv * n + col]))
					piv = r;
			}
			if (std::abs(A[piv * n + col]) < 1e-200)
				return false;

			if (piv != col)
			{
				for (int c = col; c < n; ++c)
					std::swap(A[piv * n + c], A[col * n + c]);
				std::swap(b[piv], b[col]);
			}

			const double d = A[col * n + col];
			for (int r = col + 1; r < n; ++r)
			{
				const double f = A[r * n + col] / d;
				if (f == 0.0)
					continue;
				for (int c = col; c < n; ++c)
					A[r * n + c] -= f * A[col * n + c];
				b[r] -= f * b[col];
			}
		}

		for (int r = n - 1; r >= 0; --r)
		{
			double s = b[r];
			for (int c = r + 1; c < n; ++c)
				s -= A[r * n + c] * b[c];
			b[r] = s / A[r * n + r];
		}
		return true;
	}

	/**
	 * @brief Inverts the dense matrix A (row-major, n x n) into inv via Gauss-Jordan with partial pivoting
	 * @return @c true on success, @c false if the matrix is (numerically) singular
	 */
	bool invertSmall(int n, const double* A, double* inv)
	{
		double work[9];
		for (int i = 0; i < n * n; ++i)
		{
			work[i] = A[i];
			inv[i] = 0.0;
		}
		for (int i = 0; i < n; ++i)
			inv[i * n + i] = 1.0;

		for (int col = 0; col < n; ++col)
		{
			int piv = col;
			for (int r = col + 1; r < n; ++r)
			{
				if (std::abs(work[r * n + col]) > std::abs(work[piv * n + col]))
					piv = r;
			}
			if (std::abs(work[piv * n + col]) < 1e-200)
				return false;

			if (piv != col)
			{
				for (int c = 0; c < n; ++c)
				{
					std::swap(work[piv * n + c], work[col * n + c]);
					std::swap(inv[piv * n + c], inv[col * n + c]);
				}
			}

			const double d = work[col * n + col];
			for (int c = 0; c < n; ++c)
			{
				work[col * n + c] /= d;
				inv[col * n + c] /= d;
			}

			for (int r = 0; r < n; ++r)
			{
				if (r == col)
					continue;
				const double f = work[r * n + col];
				if (f == 0.0)
					continue;
				for (int c = 0; c < n; ++c)
				{
					work[r * n + c] -= f * work[col * n + c];
					inv[r * n + c] -= f * inv[col * n + c];
				}
			}
		}
		return true;
	}

} // namespace

namespace cadet
{

// Initialization of static WENO coefficients
const double Weno::_wenoD2[2] = { 2.0/3.0, 1.0/3.0 };

const double Weno::_wenoC2[2*2] = { 0.5, -0.5,
                                    0.5,  1.5 };
const double Weno::_wenoJbvv2[2*3*3] =
{0, 2,   0,-2,   0, 0,
 0,-2,   2, 2,  -2, 0,
 0, 0,  -2, 0,   2, 0};

const double Weno::_wenoD3[3] = { 0.3, 0.6, 0.1 };

const double Weno::_wenoC3[3*3] = { 1.0/3.0, -1.0/6.0,  1.0/3.0,
                                    5.0/6.0,  5.0/6.0, -7.0/6.0,
                                   -1.0/6.0,  1.0/3.0, 11.0/6.0  };
const double Weno::_wenoJbvv3[3*5*5] = // Used to generate Jbv: vec(Jbv) = A*v
{0,0,8.0/3,   0,0,-19.0/3,        0,0,11.0/3,            0,0,0,             0,0,0,
 0,0,-19.0/3, 0,8.0/3,50.0/3,     0,-13.0/3,-31.0/3,     0,5.0/3,0,         0,0,0,
 0,0,11.0/3,  0,-13.0/3,-31.0/3,  20.0/3,26.0/3,20.0/3, -31.0/3,-13.0/3,0,  11.0/3,0,0,
 0,0,0,       0,5.0/3,0,         -31.0/3,-13.0/3,0,      50.0/3,8.0/3,0,   -19.0/3,0,0,
 0,0,0,       0,0,0,              11.0/3,0,0,           -19.0/3,0,0,         8.0/3,0,0};

bool Weno::computeGeomExactCellCoefficients(double a, double b, int weightPower, const std::vector<double>& faces,
	int cell, bool forwardFlow, int order, GeomExactCellCoefficients& out)
{
	const int r = order;
	const int dir = forwardFlow ? 1 : -1;

	// Local coordinates centered at the reconstructed face and scaled by the current cell width:
	// xi = (x - x_face) / h. The weight becomes (aP + bP * xi)^m; conditioning of the moment
	// systems is O(1) in these units. The reconstructed face sits at xi = 0, the current cell
	// spans [-1, 0] (forward flow) or [0, 1] (backward flow).
	const double faceX = forwardFlow ? faces[cell + 1] : faces[cell];
	const double h = faces[cell + 1] - faces[cell];
	const double aP = a + b * faceX;
	const double bP = b * h;

	// Geometry-weighted moments mu[q + r - 1][n] = A_{cell + dir*q}[xi^n] of all stencil cells
	// (flow offsets q = -(r-1), ..., r-1)
	double mu[5][5];
	for (int q = -(r - 1); q <= r - 1; ++q)
	{
		const int pc = cell + dir * q;
		const double xiL = (faces[pc] - faceX) / h;
		const double xiR = (faces[pc + 1] - faceX) / h;
		const double w0 = weightMomentIntegral(xiL, xiR, aP, bP, weightPower, 0);
		for (int n = 0; n <= 2 * r - 2; ++n)
			mu[q + r - 1][n] = weightMomentIntegral(xiL, xiR, aP, bP, weightPower, n) / w0;
	}

	// Candidate stencil weights: substencil rIdx covers flow offsets -rIdx, ..., -rIdx + r - 1;
	// requiring exact reproduction of the weighted averages of xi^n and evaluation at xi = 0
	// yields the Vandermonde-like system  sum_j mu_{cell(j)}[n] * C[rIdx][j] = delta_{n,0}.
	for (int rIdx = 0; rIdx < r; ++rIdx)
	{
		double A[9];
		double rhs[3];
		for (int n = 0; n < r; ++n)
		{
			rhs[n] = (n == 0) ? 1.0 : 0.0;
			for (int j = 0; j < r; ++j)
				A[n * r + j] = mu[(j - rIdx) + r - 1][n];
		}
		if (!gaussSolve(r, A, rhs))
			return false;
		for (int j = 0; j < r; ++j)
			out.C[rIdx][j] = rhs[j];
	}

	// Full-stencil weights (exactness on polynomials up to degree 2r - 2)
	const int fs = 2 * r - 1;
	double Af[25];
	double wf[5];
	for (int n = 0; n < fs; ++n)
	{
		wf[n] = (n == 0) ? 1.0 : 0.0;
		for (int j = 0; j < fs; ++j)
			Af[n * fs + j] = mu[j][n];
	}
	if (!gaussSolve(fs, Af, wf))
		return false;
	// wf[j] multiplies the stencil value at flow offset q = j - (r - 1)

	// Optimal weights: the most downwind stencil value (q = r - 1) appears only in substencil 0,
	// the most upwind one (q = -(r-1)) only in substencil r - 1; the remaining weight follows from
	// the central value. All other matching equations of the overdetermined system are verified below.
	double D[3] = { 0.0, 0.0, 0.0 };
	D[0] = wf[2 * (r - 1)] / out.C[0][r - 1];
	D[r - 1] = wf[0] / out.C[r - 1][0];
	if (r == 3)
		D[1] = (wf[r - 1] - D[0] * out.C[0][0] - D[2] * out.C[2][2]) / out.C[1][1];

	double scale = 1.0;
	for (int j = 0; j < fs; ++j)
		scale = std::max(scale, std::abs(wf[j]));
	for (int j = 0; j < fs; ++j)
	{
		const int q = j - (r - 1);
		double lhs = 0.0;
		for (int rIdx = 0; rIdx < r; ++rIdx)
		{
			const int jj = q + rIdx;
			if ((jj >= 0) && (jj < r))
				lhs += D[rIdx] * out.C[rIdx][jj];
		}
		if (std::abs(lhs - wf[j]) > 1e-8 * scale)
			return false;
	}

	double dSum = 0.0;
	for (int k = 0; k < r; ++k)
	{
		if (!(D[k] > 0.0))
			return false;
		dSum += D[k];
	}
	if (std::abs(dSum - 1.0) > 1e-8)
		return false;
	for (int k = 0; k < r; ++k)
		out.D[k] = D[k];

	// Smoothness indicator quadratic forms: with the weighted-Lagrange basis ell_j of substencil
	// rIdx (i.e., p = sum_j cbar_j ell_j), the Jiang-Shu indicator on the current cell I_i is
	//   beta = sum_{l=1}^{r-1} dx^{2l-1} int_{I_i} (d^l p / dx^l)^2 dx
	//        = sum_{m,n} B[rIdx][m][n] cbar_m cbar_n,
	// which is dimensionless in the local coordinates (the dx factors cancel exactly).
	const double xa = forwardFlow ? -1.0 : 0.0;
	const double xb = forwardFlow ? 0.0 : 1.0;
	const double i0 = xb - xa;
	const double i1 = 0.5 * (xb * xb - xa * xa);
	const double i2 = (xb * xb * xb - xa * xa * xa) / 3.0;

	for (int rIdx = 0; rIdx < r; ++rIdx)
	{
		double M[9];
		double Minv[9];
		for (int j = 0; j < r; ++j)
			for (int n = 0; n < r; ++n)
				M[j * r + n] = mu[(j - rIdx) + r - 1][n];
		if (!invertSmall(r, M, Minv))
			return false;

		// ell_j(xi) = sum_n Minv[n][j] * xi^n
		for (int j = 0; j < r; ++j)
		{
			const double c1j = Minv[1 * r + j];
			const double c2j = (r == 3) ? Minv[2 * r + j] : 0.0;
			for (int l = 0; l < r; ++l)
			{
				const double c1l = Minv[1 * r + l];
				const double c2l = (r == 3) ? Minv[2 * r + l] : 0.0;

				// int (c1j + 2 c2j xi)(c1l + 2 c2l xi) dxi + int (2 c2j)(2 c2l) dxi
				double bjl = c1j * c1l * i0 + 2.0 * (c1j * c2l + c2j * c1l) * i1 + 4.0 * c2j * c2l * i2;
				if (r == 3)
					bjl += 4.0 * c2j * c2l * i0;
				out.B[rIdx][j][l] = bjl;
			}
		}
	}

	return true;
}

void Weno::prepareGeometryExactCoefficients(double a, double b, int weightPower, const std::vector<double>& faces)
{
	if (faces.size() < 2)
		throw InvalidParameterException("prepareGeometryExactCoefficients requires at least two cell faces");
	if ((weightPower < 0) || (weightPower > 2))
		throw InvalidParameterException("prepareGeometryExactCoefficients supports weight exponents 0, 1, and 2 only");

	const int nCells = static_cast<int>(faces.size()) - 1;
	_geomForward.resize(nCells);
	_geomBackward.resize(nCells);

	for (int cell = 0; cell < nCells; ++cell)
	{
		for (int dirIdx = 0; dirIdx < 2; ++dirIdx)
		{
			const bool fwd = (dirIdx == 0);
			GeomExactCellCoefficients& out = fwd ? _geomForward[cell] : _geomBackward[cell];

			// Reduce order near boundaries such that the stencil always fits inside the domain
			const int flowIdx = fwd ? cell : (nCells - 1 - cell);
			int order = std::max(1, std::min(std::min(flowIdx + 1, _order), std::min(nCells - flowIdx, _order)));

			// Fall back to lower order if the optimal weights do not exist or are not positive
			// on this cell (possible on strongly distorted grids)
			while ((order > 1) && !computeGeomExactCellCoefficients(a, b, weightPower, faces, cell, fwd, order, out))
				--order;

			out.order = order;
		}
	}
}

} // namespace cadet
