// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//
//  Copyright © 2008-2020: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file 
 * Defines helper functions used by all ModelSystemImpl-xyz.cpp files.
 */

namespace
{
	/**
	 * @brief Computes a total return code from a list of separate return codes
	 * @details A negative return code indicates a non-recoverable error. Positive
	 *          values indicate recoverable errors and a value of @c 0 indicates
	 *          no error.
	 * @param [in] err List of error codes to be fused into one
	 * @return Total error code summarizing all codes in the list
	 */
	inline int totalErrorIndicatorFromLocal(const std::vector<int>& err)
	{
		int totalError = 0;
		for (unsigned int i = 0; i < err.size(); ++i)
		{
			// Negative values are non-recoverable errors
			if (err[i] < 0)
				return err[i];

			// 0 = okay, positive values = recoverable error
			totalError = std::max(totalError, err[i]);
		}
		return totalError;
	}

	/**
	 * @brief Fuses two error codes into one
	 * @details A negative return code indicates a non-recoverable error. Positive
	 *          values indicate recoverable errors and a value of @c 0 indicates
	 *          no error.
	 * @param [in] curCode Current error code
	 * @param [in] nextCode Next error code that is fused into the current one
	 * @return Fused error code summarizing both inputs
	 */
	inline int updateErrorIndicator(int curCode, int nextCode)
	{
		if ((curCode < 0) || (nextCode < 0))
			return std::min(curCode, nextCode);
		return std::max(curCode, nextCode);
	}

	template <typename state_t>
	inline state_t applyOffset(const state_t& state, unsigned int offset)
	{
		return state_t{
			state.vecStateY + offset,
			(state.vecStateYdot) ? (state.vecStateYdot + offset) : nullptr
		};
	}

	template <>
	inline cadet::AdJacobianParams applyOffset(const cadet::AdJacobianParams& adJac, unsigned int offset)
	{
		return cadet::AdJacobianParams{
			(adJac.adRes) ? (adJac.adRes + offset) : nullptr,
			(adJac.adY) ? (adJac.adY + offset) : nullptr,
			adJac.adDirOffset
		};
	}
}
