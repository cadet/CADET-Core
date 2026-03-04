// =============================================================================
//  CADET
//  
//  Copyright © 2008-present: 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file 
 * Implements the High Resolution Koren method
 */

#ifndef LIBCADET_HIGHRESKOREN_HPP_
#define LIBCADET_HIGHRESKOREN_HPP_

#include "AutoDiff.hpp"
#include "MathUtil.hpp"
#include "Memory.hpp"
#include "common/CompilerSpecific.hpp"
#include "cadet/Exceptions.hpp"

#include <algorithm>
#include <vector>

namespace cadet
{

/**
 * @brief Implements the High Resolution Koren scheme for convection
 * @details This scheme uses a van Leer flux limiter to help the scheme
 * to stay in the TVD region. The limiter uses a smoothness monitor to
 * monitor the smoothness of the solution and adjust the scheme between
 * first and second order. \varepsilon in the van Leer flux limiter is
 * set to 1e-10. The BOUNDARY_MODEL is set to 0. 
 *
 * On non-equidistant grids, calling setGrid() enables a corrected
 * gradient-ratio limiter and proper face-distance weighting, maintaining
 * true second-order accuracy regardless of cell-size variation.
 */
class HighResolutionKoren
{
public:

	/**
	 * @brief Creates the HighResolutionKoren scheme
	 * @details The max order is 2. 
	 */
	HighResolutionKoren() : _epsilon(), _order(2), _forwardFlow(true) { }

	/**
    * @brief Returns the maximum order \f$ r \f$ of the implemented schemes
    * @return Maximum WENO order \f$ r \f$
    */
	CADET_CONSTEXPR static inline unsigned int maxOrder() CADET_NOEXCEPT { return 2; }

	/**
	 * @brief Returns the maximum stencil size for the implemented schemes
	 * @return Maximum stencil size
	 */
	CADET_CONSTEXPR static inline unsigned int maxStencilSize() CADET_NOEXCEPT { return 2 * maxOrder() - 1; }

	/**
	 * @brief Reconstructs a cell face value from volume averages
	 * @param [in] cellIdx Index of the current cell
	 * @param [in] numCells Number of cells
	 * @param [in] w Stencil that contains 3 volume averages from which the cell face values are reconstructed centered at the 
	 *               current cell (i.e., index 0 is the current cell, -2 the next to previous cell, 2 the next but one cell)
	 * @param [out] result Reconstructed cell face value
	 * @param [out] Dvm Gradient of the reconstructed cell face value. The array needs to be of size 3, which are the derivatives of right face flux wrt n_{i-1}, n_{i}, n_{i+1}. 
	 * @tparam StateType Type of the state variables
	 * @tparam StencilType Type of the stencil (can be a dedicated class with overloaded operator[] or a simple pointer)
	 * @return Order of the scheme that was used in the computation: 2 or 1
	 */
	template <typename StateType, typename StencilType>
	int reconstruct(unsigned int cellIdx, unsigned int numCells, const StencilType& w, StateType& result, double* const Dvm)
	{
		return reconstruct<StateType, StencilType, true>(cellIdx, numCells, w, result, Dvm);
	}

	/**
	 * @brief Reconstructs a cell face value from volume averages
	 * @param [in] cellIdx Index of the current cell
	 * @param [in] numCells Number of cells
	 * @param [in] w Stencil that contains 3 volume averages from which the cell face values are reconstructed centered at the 
	 *               current cell (i.e., index 0 is the current cell, -2 the next to previous cell, 2 the next but one cell)
	 * @param [out] result Reconstructed cell face value
	 * @tparam StateType Type of the state variables
	 * @tparam StencilType Type of the stencil (can be a dedicated class with overloaded operator[] or a simple pointer)
	 * @return Order of the scheme that was used in the computation: 2 or 1
	 */
	template <typename StateType, typename StencilType>
	int reconstruct(unsigned int cellIdx, unsigned int numCells, const StencilType& w, StateType& result)
	{
		return reconstruct<StateType, StencilType, false>(cellIdx, numCells, w, result, nullptr);
	}

	/**
    * @brief Sets the \f$ \varepsilon \f$ of the WENO emthod (prevents division by zero in the weights)
    * @param [in] HR Koren \f$ \varepsilon \f$
    */
	inline void epsilon(double eps) { _epsilon = eps; }

	/**
	 * @brief Returns the \f$ \varepsilon \f$ of the HR Koren emthod (prevents division by zero in the weights)
	 * @return HR Koren \f$ \varepsilon \f$
	 */
	inline double epsilon() const CADET_NOEXCEPT { return _epsilon; }

	/**
	 * @brief Sets the order
	 * @param [in] order Order of the HR Koren method
	 */
	inline void order(int order)
	{
		cadet_assert(order <= 2);
		cadet_assert(order > 0);
		_order = order;
	}

	/**
	 * @brief Returns the HR Koren order
	 * @return Order of the HR Koren method
	 */
	inline int order() const CADET_NOEXCEPT { return _order; }

	/**
	 * @brief Returns the number of upper diagonals required in the Jacobian
	 * @return Number of required Jacobian upper diagonals
	 */
	inline unsigned int upperBandwidth() const CADET_NOEXCEPT { return _order - 1; }

	/**
	 * @brief Returns the number of lower diagonals required in the Jacobian
	 * @return Number of required Jacobian lower diagonals
	 */
	inline unsigned int lowerBandwidth() const CADET_NOEXCEPT { return _order - 1; }

	/**
	 * @brief Returns the size of the stencil (i.e., the number of required elements)
	 * @return Size of the stencil. Either 3 or 1. 
	 */
	inline unsigned int stencilSize() const CADET_NOEXCEPT { return 2 * _order - 1; }

	/**
	 * @brief Configures the scheme for a non-equidistant grid
	 * @details Stores the grid geometry and pre-computes per-cell spacings.
	 *          After calling this, reconstruct() applies a gradient-ratio limiter
	 *          using actual cell-centre distances, and weights the face value by
	 *          the true distance from cell centre to face rather than h/2.
	 *          On an equidistant grid this is mathematically identical to the
	 *          standard equidistant path.
	 *
	 *          Call setForwardFlow() to select which face (right or left) is
	 *          reconstructed; the default is forward (right face).
	 *
	 * @param [in] cellCenters  Cell-centre coordinates, length @p nCells
	 * @param [in] cellFaces    Cell-face coordinates, length @p nCells + 1
	 * @param [in] nCells       Number of cells
	 */
	void setGrid(const double* cellCenters, const double* cellFaces, unsigned int nCells)
	{
		// Store copies so the caller's arrays do not need to outlive this object.
		_centersStore.assign(cellCenters, cellCenters + nCells);
		_facesStore.assign(cellFaces, cellFaces + nCells + 1);

		// Pre-compute h_L[i] = x_i - x_{i-1} and h_R[i] = x_{i+1} - x_i.
		// Boundary entries remain 0.0 (never accessed at order-1 cells).
		_hLeft.assign(nCells, 0.0);
		_hRight.assign(nCells, 0.0);
		for (unsigned int i = 1; i < nCells; ++i)
			_hLeft[i] = cellCenters[i] - cellCenters[i - 1];
		for (unsigned int i = 0; i + 1 < nCells; ++i)
			_hRight[i] = cellCenters[i + 1] - cellCenters[i];
	}

	/**
	 * @brief Selects which face the subsequent reconstruct() calls target
	 * @param [in] forward  @c true  = right face of each cell (positive / outward flow),
	 *                      @c false = left  face of each cell (negative / inward  flow)
	 */
	inline void setForwardFlow(bool forward) CADET_NOEXCEPT { _forwardFlow = forward; }

private:

	 /**
      * @brief Reconstructs a cell face value from volume averages
      * @details When setGrid() has been called this function uses the corrected
      *          non-equidistant van Leer gradient-ratio formulation; otherwise it
      *          falls back to the standard equidistant formula.
      *
      * @param [in] cellIdx Index of the current cell
      * @param [in] numCells Number of cells
      * @param [in] w Stencil that contains 3 volume averages from which the cell face values are reconstructed centered at the
      *               current cell (i.e., index 0 is the current cell, -2 the next to previous cell, 2 the next but one cell)
      * @param [out] result Reconstructed cell face value
      * @param [out] Dvm Gradient of the reconstructed cell face value. The array needs to be of size 3, which are the derivatives of right face flux wrt c_{i-1}, c_{i}, c_{i+1}.
      * @tparam StateType Type of the state variables
      * @tparam StencilType Type of the stencil (can be a dedicated class with overloaded operator[] or a simple pointer)
      * @tparam wantJac Determines if the gradient is computed (@c true) or not (@c false)
      * @return Order of the scheme that was used in the computation: 2 or 1
      */

	template <typename StateType, typename StencilType, bool wantJac>
	int reconstruct(unsigned int cellIdx, unsigned int numCells, const StencilType& w, StateType& result, double* const Dvm)
	{
#if defined(ACTIVE_SETFAD) || defined(ACTIVE_SFAD)
		using cadet::sqr;
		using sfad::sqr;
#endif

		// Select the effective order: reduce near boundaries so the stencil fits.
		const int order = std::min(std::min(static_cast<int>(cellIdx) + 1, _order),
		                           std::min(static_cast<int>(numCells - cellIdx), _order));

		// First-order upwind – no grid-geometry correction needed.
		if (order == 1)
		{
			result = w[0];
			if (wantJac)
				*Dvm = 1.0;
			return order;
		}

		// ----------------------------------------------------------------
		// Non-equidistant path (active when setGrid() has been called)
		// ----------------------------------------------------------------
		if (!_hLeft.empty())
		{
			// For forward flow (right face of cell i):
			//   hR = x_{i+1} - x_i   (stencil "forward" step in physical space)
			//   hL = x_i   - x_{i-1} (stencil "backward" step)
			//   fR = x_{i+1/2} - x_i (centre-to-face distance)
			//
			// For backward flow the stencil is mirrored: w[1] = c_{i-1}, w[-1] = c_{i+1},
			// and we reconstruct the LEFT face.  Physical h_L / h_R roles are swapped,
			// and fR becomes the distance from centre to the LEFT face.
			const double hR = _forwardFlow ? _hRight[cellIdx] : _hLeft[cellIdx];
			const double hL = _forwardFlow ? _hLeft[cellIdx]  : _hRight[cellIdx];
			const double fR = _forwardFlow
			    ? (_facesStore[cellIdx + 1] - _centersStore[cellIdx])
			    : (_centersStore[cellIdx]   - _facesStore[cellIdx]);

			// Scale factor k = 2*fR/hR (equals 1 on a uniform grid where fR = hR/2).
			const double k = 2.0 * fR / hR;

			// Gradient-ratio numerator and denominator (with regularisation ε).
			//   A ~ (backward gradient) * hR   → r numerator
			//   B ~ (forward  gradient) * hL   → r denominator
			// Dividing A/B gives the gradient ratio r without explicit division by h.
			StateType const A = (w[0] - w[-1]) * hR + _epsilon;
			StateType const B = (w[1] - w[ 0]) * hL + _epsilon;

			// van Leer flux limiter φ(r) = 2r/(1+r); vLf = φ/2 = r/(1+r) = A/(A+B).
			// TVD condition: limiter active only when r > 0.
			const bool rPositive = static_cast<double>(A) * static_cast<double>(B) > 0.0;
			StateType const vLf  = rPositive ? A / (A + B) : StateType(0.0);

			result = w[0] + k * vLf * (w[1] - w[0]);

			if (wantJac)
			{
				if (rPositive)
				{
					const double Av    = static_cast<double>(A);
					const double Bv    = static_cast<double>(B);
					const double ApB   = Av + Bv;
					const double vv    = static_cast<double>(vLf);
					const double Delta = static_cast<double>(w[1] - w[0]);

					// d(result)/d(w[-1])
					Dvm[0] = -k * hR * Delta * Bv / (ApB * ApB);
					// d(result)/d(w[0])
					Dvm[1] =  1.0 + k * Delta * (hR * Bv + Av * hL) / (ApB * ApB) - k * vv;
					// d(result)/d(w[1])
					Dvm[2] =  k * Av / ApB * (1.0 - hL * Delta / ApB);
				}
				else
				{
					// Upwind limit: result = w[0]
					Dvm[0] = 0.0;
					Dvm[1] = 1.0;
					Dvm[2] = 0.0;
				}
			}
			return order;
		}

		// ----------------------------------------------------------------
		// Equidistant fallback (original formulation)
		// ----------------------------------------------------------------

		// calculate the smoothness indicators
		StateType const r = (w[0] - w[-1] + _epsilon) / (w[1] - w[0] + _epsilon);

		// apply the van Leer flux limiter. Other flux limiters may be applied here.
		StateType const vLf = (cadet_likely(r > 0)) ? r / (1 + r) : 0.0;

		// calculate the flux (0.5 * 2.0 are merged)
		result = w[0] + vLf * (w[1] - w[0]);

		// define some jacobian related local variables
		StateType const r_1_1 = 1.0 / (1.0 + r);
		StateType const ci_complex = 1.0 - _epsilon / (w[1] - w[0] + _epsilon);

		if (wantJac)
		{ 
			// Dvm start with index 0 (drivative wrt c_{i-1})
			if (cadet_likely(r > 0))
			{
				Dvm[0] = static_cast<double>(r_1_1) * (static_cast<double>(vLf) - 1.0) * static_cast<double>(ci_complex);
				Dvm[1] = 1.0 - static_cast<double>(ci_complex) * static_cast<double>(vLf) * static_cast<double>(vLf) - static_cast<double>(ci_complex) * static_cast<double>(r_1_1) * static_cast<double>(vLf) + static_cast<double>(ci_complex) * static_cast<double>(vLf) - static_cast<double>(vLf) + static_cast<double>(r_1_1) * static_cast<double>(ci_complex);
				Dvm[2] = static_cast<double>(vLf) * static_cast<double>(vLf) * static_cast<double>(ci_complex) - static_cast<double>(vLf) * static_cast<double>(ci_complex) + static_cast<double>(vLf);
			}
			else
			{
				Dvm[0] = 0.0;
				Dvm[1] = 1.0; // becomes upwind when r<=0
				Dvm[2] = 0.0;
			}
		}

		return order;
	}

	int    _order;        //!< Selected order: 2 or 1 for HR Koren
	double _epsilon;      //!< Small number preventing divsion by zero
	bool   _forwardFlow;  //!< true = reconstruct right face, false = left face

	// Non-equidistant grid data (empty when equidistant path is used)
	std::vector<double> _centersStore; //!< Owned copy of cell-centre positions
	std::vector<double> _facesStore;   //!< Owned copy of cell-face positions
	std::vector<double> _hLeft;        //!< h_L[i] = x_i - x_{i-1}
	std::vector<double> _hRight;       //!< h_R[i] = x_{i+1} - x_i

};

} // namespace cadet

#endif  // LIBCADET_HIGHRESKOREN_HPP_
