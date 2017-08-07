// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2017: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file 
 * Defines helper functions for checking Jacobians.
 */

#ifndef CADETTEST_JACOBIANHELPER_HPP_
#define CADETTEST_JACOBIANHELPER_HPP_

#include "UnitOperation.hpp"

namespace cadet
{

namespace test
{

/**
 * @brief Compares two given Jacobians column by column
 * @details A column is extracted by calling multiplyWithJacobian() on the model.
 * @param [in] modelA Model A
 * @param [in] modelB Model B
 * @param [in] y State vector
 * @param [in] yDot Time derivative of state vector
 * @param [in] dir Memory for extracting a column
 * @param [in] colA Memory for Jacobian column of @p modelA
 * @param [in] colB Memory for Jacobian column of @p modelB
 */
inline void compareJacobian(cadet::IUnitOperation* modelA, cadet::IUnitOperation* modelB, double const* y, double const* yDot, double* dir, double* colA, double* colB)
{
	const unsigned int n = modelA->numDofs();
	std::fill(dir, dir + n, 0.0);
	for (unsigned int col = 0; col < n; ++col)
	{
		dir[col] = 1.0;

		modelA->multiplyWithJacobian(0.0, 0u, 1.0, y, yDot, dir, 1.0, 0.0, colA);
		modelB->multiplyWithJacobian(0.0, 0u, 1.0, y, yDot, dir, 1.0, 0.0, colB);

		for (unsigned int row = 0; row < n; ++row)
		{
			CAPTURE(row);
			CAPTURE(col);
			CHECK(colA[row] == Approx(colB[row]));
		}

		dir[col] = 0.0;
	}
}

/**
 * @brief Checks a Jacobian against finite differences
 * @details Uses finite differences on @p modelA to determine the Jacobian. The two
 *          Jacobians are compared column by column. A column is extracted by calling
 *          multiplyWithJacobian() on the model.
 * @param [in] modelA Model A
 * @param [in] modelB Model B
 * @param [in] y State vector
 * @param [in] yDot Time derivative of state vector
 * @param [in] dir Memory for computing finite differences
 * @param [in] colA Memory for Jacobian column of @p modelA
 * @param [in] colB Memory for Jacobian column of @p modelB
 */
inline void compareJacobianFD(cadet::IUnitOperation* modelA, cadet::IUnitOperation* modelB, double const* y, double const* yDot, double* dir, double* colA, double* colB)
{
	const double h = 1e-6;
	const unsigned int n = modelA->numDofs();
	for (unsigned int col = 0; col < n; ++col)
	{
		std::copy(y, y + n, dir);
		if (y[col] != 0.0)
			dir[col] = y[col] * (1.0 + h);
		else
			dir[col] = h;

		modelA->residual(0.0, 0u, 1.0, dir, yDot, colA);

		if (y[col] != 0.0)
			dir[col] = y[col] * (1.0 - h);
		else
			dir[col] = -h;

		modelA->residual(0.0, 0u, 1.0, dir, yDot, colB);

		for (unsigned int j = 0; j < n; ++j)
		{
			if (y[col] != 0.0)
				colA[j] = (colA[j] - colB[j]) / (y[col] * 2.0 * h);
			else
				colA[j] = (colA[j] - colB[j]) / (2.0 * h);
		}

		std::fill(dir, dir + n, 0.0);
		dir[col] = 1.0;
		modelB->multiplyWithJacobian(0.0, 0u, 1.0, y, yDot, dir, 1.0, 0.0, colB);

		for (unsigned int row = 0; row < n; ++row)
		{
			CAPTURE(row);
			CAPTURE(col);
			CHECK(colA[row] == Approx(colB[row]));
		}
	}
}

/**
 * @brief Checks the pattern of the two given Jacobians including sign of the entries
 * @details Uses finite differences on @p modelA to determine the Jacobian structure. The
 *          Jacobians are compared column by column. A column is extracted by calling
 *          multiplyWithJacobian() on the model.
 * @param [in] modelA Model A
 * @param [in] modelB Model B
 * @param [in] y State vector
 * @param [in] yDot Time derivative of state vector
 * @param [in] dir Memory for computing finite differences
 * @param [in] colA Memory for Jacobian column of @p modelA
 * @param [in] colB Memory for Jacobian column of @p modelB
 */
inline void checkJacobianPatternFD(cadet::IUnitOperation* modelA, cadet::IUnitOperation* modelB, double const* y, double const* yDot, double* dir, double* colA, double* colB)
{
	const double h = 1e-5;
	const unsigned int n = modelA->numDofs();
	for (unsigned int col = 0; col < n; ++col)
	{
		std::copy(y, y + n, dir);
		if (y[col] != 0.0)
			dir[col] = y[col] * (1.0 + h);
		else
			dir[col] = h;

		modelA->residual(0.0, 0u, 1.0, dir, yDot, colA);

		if (y[col] != 0.0)
			dir[col] = y[col] * (1.0 - h);
		else
			dir[col] = -h;

		modelA->residual(0.0, 0u, 1.0, dir, yDot, colB);

		for (unsigned int j = 0; j < n; ++j)
		{
			if (y[col] != 0.0)
				colA[j] = (colA[j] - colB[j]) / (y[col] * 2.0 * h);
			else
				colA[j] = (colA[j] - colB[j]) / (2.0 * h);
		}

		std::fill(dir, dir + n, 0.0);
		dir[col] = 1.0;
		modelB->multiplyWithJacobian(0.0, 0u, 1.0, y, yDot, dir, 1.0, 0.0, colB);

		// Check for pattern including sign
		for (unsigned int row = 0; row < n; ++row)
		{
			CAPTURE(row);
			CAPTURE(col);
			if (colA[row] == 0.0)
				CHECK(colB[row] == 0.0);
			else if (colA[row] > 0.0)
				CHECK(colB[row] > 0.0);
			else if (colA[row] < 0.0)
				CHECK(colB[row] < 0.0);
			else if (std::isnan(colA[row]))
				CHECK(std::isnan(colB[row]));
		}
	}
}

/**
 * @brief Compares two given time derivative Jacobians column by column
 * @details A column is extracted by calling multiplyWithDerivativeJacobian() on the model.
 * @param [in] modelA Model A
 * @param [in] modelB Model B
 * @param [in] y State vector
 * @param [in] yDot Time derivative of state vector
 * @param [in] dir Memory for extracting a column
 * @param [in] colA Memory for Jacobian column of @p modelA
 * @param [in] colB Memory for Jacobian column of @p modelB
 */
inline void compareTimeDerivativeJacobian(cadet::IUnitOperation* modelA, cadet::IUnitOperation* modelB, double const* y, double const* yDot, double* dir, double* colA, double* colB)
{
	const unsigned int n = modelA->numDofs();
	std::fill(dir, dir + n, 0.0);
	for (unsigned int col = 0; col < n; ++col)
	{
		dir[col] = 1.0;

		modelA->multiplyWithDerivativeJacobian(0.0, 0u, 1.0, y, yDot, dir, colA);
		modelB->multiplyWithDerivativeJacobian(0.0, 0u, 1.0, y, yDot, dir, colB);

		for (unsigned int row = 0; row < n; ++row)
		{
			CAPTURE(row);
			CAPTURE(col);
			CHECK(colA[row] == Approx(colB[row]));
		}

		dir[col] = 0.0;
	}	
}

/**
 * @brief Checks a time derivative Jacobian against finite differences
 * @details Uses finite differences on @p modelA to determine the time derivative
 *          Jacobian. The two Jacobians are compared column by column. A column is
 *          extracted by calling multiplyWithDerivativeJacobian() on the model.
 * @param [in] modelA Model A
 * @param [in] modelB Model B
 * @param [in] y State vector
 * @param [in] yDot Time derivative of state vector
 * @param [in] dir Memory for computing finite differences
 * @param [in] colA Memory for Jacobian column of @p modelA
 * @param [in] colB Memory for Jacobian column of @p modelB
 */
inline void compareTimeDerivativeJacobianFD(cadet::IUnitOperation* modelA, cadet::IUnitOperation* modelB, double const* y, double const* yDot, double* dir, double* colA, double* colB)
{
	const double h = 1e-6;
	const unsigned int n = modelA->numDofs();
	for (unsigned int col = 0; col < n; ++col)
	{
		std::copy(yDot, yDot + n, dir);

		if (yDot[col] != 0.0)
			dir[col] = yDot[col] * (1.0 + h);
		else
			dir[col] = h;

		modelA->residual(0.0, 0u, 1.0, y, dir, colA);

		if (yDot[col] != 0.0)
			dir[col] = yDot[col] * (1.0 - h);
		else
			dir[col] = -h;

		modelA->residual(0.0, 0u, 1.0, y, dir, colB);

		for (unsigned int j = 0; j < n; ++j)
		{
			if (y[col] != 0.0)
				colA[j] = (colA[j] - colB[j]) / (yDot[col] * 2.0 * h);
			else
				colA[j] = (colA[j] - colB[j]) / (2.0 * h);
		}

		std::fill(dir, dir + n, 0.0);
		dir[col] = 1.0;
		modelB->multiplyWithDerivativeJacobian(0.0, 0u, 1.0, y, yDot, dir, colB);

		for (unsigned int row = 0; row < n; ++row)
		{
			CAPTURE(row);
			CAPTURE(col);
			CHECK(colA[row] == Approx(colB[row]));
		}
	}
}

} // namespace test
} // namespace cadet

#endif  // CADETTEST_JACOBIANHELPER_HPP_
