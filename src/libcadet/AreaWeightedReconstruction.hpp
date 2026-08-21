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

#ifndef LIBCADET_AREAWEIGHTEDRECONSTRUCTION_HPP_
#define LIBCADET_AREAWEIGHTEDRECONSTRUCTION_HPP_

#include "AutoDiff.hpp"

#include <cmath>
#include <vector>

namespace cadet
{

template <typename Container>
struct ReverseMomentAccessor
{
	const Container& data;
	inline auto operator[](unsigned int idx) const -> decltype(data[0]) { return data[data.size() - 1 - idx]; }
};

template <typename ContainerType>
struct AreaMomentsBundle
{
	ContainerType centroids;
	ContainerType secondMoments;
	ContainerType thirdMoments;
	ContainerType fourthMoments;
};

inline void computeRadialAreaMoments(
	unsigned int nCol,
	const active* cellBounds,
	std::vector<active>& centroids,
	std::vector<active>& secondMoments,
	std::vector<active>& thirdMoments,
	std::vector<active>& fourthMoments)
{
	centroids.resize(nCol);
	secondMoments.resize(nCol);
	thirdMoments.resize(nCol);
	fourthMoments.resize(nCol);

	for (unsigned int j = 0; j < nCol; ++j)
	{
		const active rL = cellBounds[j];
		const active rR = cellBounds[j + 1];

		const active rL2 = rL * rL;
		const active rR2 = rR * rR;
		const active rL3 = rL2 * rL;
		const active rR3 = rR2 * rR;
		const active rL4 = rL2 * rL2;
		const active rR4 = rR2 * rR2;
		const active rL5 = rL4 * rL;
		const active rR5 = rR4 * rR;
		const active rL6 = rL3 * rL3;
		const active rR6 = rR3 * rR3;

		// V_j^{(k)} = integral rho^{k+1} drho = (rR^{k+2} - rL^{k+2}) / (k+2)
		const active V0 = (rR2 - rL2) * 0.5;
		const active V1 = (rR3 - rL3) / 3.0;
		const active V2 = (rR4 - rL4) * 0.25;
		const active V3 = (rR5 - rL5) * 0.2;
		const active V4 = (rR6 - rL6) / 6.0;

		centroids[j] = V1 / V0;
		secondMoments[j] = V2 / V0;
		thirdMoments[j] = V3 / V0;
		fourthMoments[j] = V4 / V0;
	}
}

inline void computeFrustumAreaMoments(
	unsigned int nCol,
	const active* cellFaces,
	active innerRadius,
	active outerRadius,
	active colLength,
	std::vector<active>& centroids,
	std::vector<active>& secondMoments,
	std::vector<active>& thirdMoments,
	std::vector<active>& fourthMoments)
{
	centroids.resize(nCol);
	secondMoments.resize(nCol);
	thirdMoments.resize(nCol);
	fourthMoments.resize(nCol);

	const active a = innerRadius;
	const active b = (outerRadius - innerRadius) / colLength;
	const active a2 = a * a;
	const active ab2 = 2.0 * a * b;
	const active b2 = b * b;

	for (unsigned int j = 0; j < nCol; ++j)
	{
		const active xL = cellFaces[j];
		const active xR = cellFaces[j + 1];

		// A(x) = (a + bx)^2 = a^2 + 2abx + b^2 x^2
		// V_j^{(k)} = a^2 * (xR^{k+1}-xL^{k+1})/(k+1) + 2ab * (xR^{k+2}-xL^{k+2})/(k+2) + b^2 * (xR^{k+3}-xL^{k+3})/(k+3)
		const active dx1 = xR - xL;
		const active dx2 = xR * xR - xL * xL;
		const active dx3 = xR * xR * xR - xL * xL * xL;
		const active dx4 = xR * xR * xR * xR - xL * xL * xL * xL;
		const active dx5 = xR * xR * xR * xR * xR - xL * xL * xL * xL * xL;
		const active dx6 = xR * xR * xR * xR * xR * xR - xL * xL * xL * xL * xL * xL;
		const active dx7 = xR * xR * xR * xR * xR * xR * xR - xL * xL * xL * xL * xL * xL * xL;

		const active V0 = a2 * dx1       + ab2 * dx2 * 0.5  + b2 * dx3 / 3.0;
		const active V1 = a2 * dx2 * 0.5 + ab2 * dx3 / 3.0  + b2 * dx4 * 0.25;
		const active V2 = a2 * dx3 / 3.0 + ab2 * dx4 * 0.25 + b2 * dx5 * 0.2;
		const active V3 = a2 * dx4 * 0.25 + ab2 * dx5 * 0.2  + b2 * dx6 / 6.0;
		const active V4 = a2 * dx5 * 0.2  + ab2 * dx6 / 6.0  + b2 * dx7 / 7.0;

		centroids[j] = V1 / V0;
		secondMoments[j] = V2 / V0;
		thirdMoments[j] = V3 / V0;
		fourthMoments[j] = V4 / V0;
	}
}

namespace detail
{

inline void solve3x3FirstRow(const double mu[3], const double nu[3], double coeff[3])
{
	// M = | 1  mu[0]  nu[0] |
	//     | 1  mu[1]  nu[1] |
	//     | 1  mu[2]  nu[2] |
	// Compute first row of M^{-1} via cofactors of first column
	const double cf0 = mu[1] * nu[2] - nu[1] * mu[2];
	const double cf1 = -(mu[0] * nu[2] - nu[0] * mu[2]);
	const double cf2 = mu[0] * nu[1] - nu[0] * mu[1];
	const double det = cf0 + cf1 + cf2;

	coeff[0] = cf0 / det;
	coeff[1] = cf1 / det;
	coeff[2] = cf2 / det;
}

inline void solve5x5FirstRow(double M[5][5], double coeff[5])
{
	double rhs[5] = {1.0, 0.0, 0.0, 0.0, 0.0};

	// Transpose M (we solve M^T y = e_0 to get first row of M^{-1})
	for (int i = 0; i < 5; ++i)
		for (int j = i + 1; j < 5; ++j)
			std::swap(M[i][j], M[j][i]);

	for (int col = 0; col < 5; ++col)
	{
		int maxRow = col;
		double maxVal = std::abs(M[col][col]);
		for (int row = col + 1; row < 5; ++row)
		{
			const double v = std::abs(M[row][col]);
			if (v > maxVal) { maxVal = v; maxRow = row; }
		}
		if (maxRow != col)
		{
			std::swap(rhs[col], rhs[maxRow]);
			for (int j = 0; j < 5; ++j)
				std::swap(M[col][j], M[maxRow][j]);
		}
		for (int row = col + 1; row < 5; ++row)
		{
			const double factor = M[row][col] / M[col][col];
			rhs[row] -= factor * rhs[col];
			for (int j = col + 1; j < 5; ++j)
				M[row][j] -= factor * M[col][j];
		}
	}

	for (int row = 4; row >= 0; --row)
	{
		coeff[row] = rhs[row];
		for (int j = row + 1; j < 5; ++j)
			coeff[row] -= M[row][j] * coeff[j];
		coeff[row] /= M[row][row];
	}
}

template <typename MomentsType>
inline void computeMoments(unsigned int cellPhysIdx, double xf, const MomentsType& moments,
	double& mu1, double& mu2, double& mu3, double& mu4)
{
	const double c1 = static_cast<double>(moments.centroids[cellPhysIdx]);
	const double c2 = static_cast<double>(moments.secondMoments[cellPhysIdx]);
	const double c3 = static_cast<double>(moments.thirdMoments[cellPhysIdx]);
	const double c4 = static_cast<double>(moments.fourthMoments[cellPhysIdx]);

	mu1 = c1 - xf;
	mu2 = c2 - 2.0 * xf * c1 + xf * xf;
	mu3 = c3 - 3.0 * xf * c2 + 3.0 * xf * xf * c1 - xf * xf * xf;
	mu4 = c4 - 4.0 * xf * c3 + 6.0 * xf * xf * c2 - 4.0 * xf * xf * xf * c1 + xf * xf * xf * xf;
}

} // namespace detail

} // namespace cadet

#endif // LIBCADET_AREAWEIGHTEDRECONSTRACTION_HPP_
