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

/**
 * @file 
 * Implements the WENO method
 */

#ifndef LIBCADET_WENO_HPP_
#define LIBCADET_WENO_HPP_

#include "AutoDiff.hpp"
#include "MathUtil.hpp"
#include "Memory.hpp"
#include "common/CompilerSpecific.hpp"
#include "cadet/Exceptions.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace cadet
{

/**
 * @brief Implements the WENO scheme for convection
 * @details This scheme is based on upwind stencils and provides WENO methods 1-1, 2-3, and 3-5.
 *          In general, WENO achieves order \f$ r \f$ using a stencil with \f$ (2r-1) \f$ points
 *          that is subdivided into \f$ r \f$ smaller stencils having \f$ r \f$ points each.
 *          WENO combines all substencils with an estimate of the smoothness of the solution (also obtained from the
 *          substencils) in order to achieve a non-oscillatory high order reconstruction of the face values given
 *          volume averages (cell boundary fluxes in finite volume schemes).
 *          For details see \cite Liu1994 and \cite Jiang1996.
 */
class Weno
{
public:

	/**
	 * @brief Boundary treatment method determines how the reconstruction handles insufficient available elements (i.e., less elements available than stencil size)
	 */
	enum class BoundaryTreatment : int
	{
		ReduceOrder = 0, //!< Reduce the order of the WENO method such that the stencil is small enough
		ZeroWeights = 1,  //!< Set weights of WENO method to 0 for unavailable elements
		ZeroWeightsForPnotZero = 2, //!< Set weights of WENO method to 0 for unavailable elements, except for the first cell (order is reduced to 1)
		LargeGhostNodes = 3,
	};

	/**
	 * @brief Creates the WENO scheme
	 */
	Weno() : _epsilon(1e-10), _order(maxOrder()), _boundaryTreatment(BoundaryTreatment::ReduceOrder), _intermediateValues(3 * maxOrder() * sizeof(active)) { }

	/**
	 * @brief Returns the maximum order \f$ r \f$ of the implemented schemes
	 * @return Maximum WENO order \f$ r \f$
	 */
	CADET_CONSTEXPR static inline unsigned int maxOrder() CADET_NOEXCEPT { return 3; }

	/**
	 * @brief Returns the maximum stencil size for the implemented schemes
	 * @return Maximum stencil size
	 */
	CADET_CONSTEXPR static inline unsigned int maxStencilSize() CADET_NOEXCEPT { return 2 * maxOrder() - 1; }

	/**
	 * @brief Reconstructs a cell face value from volume averages
	 * @param [in] cellIdx Index of the current cell
	 * @param [in] numCells Number of cells
	 * @param [in] w Stencil that contains the \f$ 2r-1 \f$ volume averages from which the cell face values are reconstructed centered at the 
	 *               current cell (i.e., index 0 is the current cell, -2 the next to previous cell, 2 the next but one cell)
	 * @param [out] result Reconstructed cell face value
	 * @param [out] Dvm Gradient of the reconstructed cell face value (array has to be of size \f$ 2r-1\f$ where \f$ r \f$ is the WENO order)
	 * @tparam StateType Type of the state variables
	 * @tparam StencilType Type of the stencil (can be a dedicated class with overloaded operator[] or a simple pointer)
	 * @return Order of the WENO scheme that was used in the computation
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
	 * @param [in] w Stencil that contains the \f$ 2r-1 \f$ volume averages from which the cell face values are reconstructed centered at the 
	 *               current cell (i.e., index 0 is the current cell, -2 the next to previous cell, 2 the next but one cell)
	 * @param [out] result Reconstructed cell face value
	 * @tparam StateType Type of the state variables
	 * @tparam StencilType Type of the stencil (can be a dedicated class with overloaded operator[] or a simple pointer)
	 * @return Order of the WENO scheme that was used in the computation
	 */
	template <typename StateType, typename StencilType>
	int reconstruct(unsigned int cellIdx, unsigned int numCells, const StencilType& w, StateType& result)
	{
		return reconstruct<StateType, StencilType, false>(cellIdx, numCells, w, result, nullptr);
	}

	template <typename StateType, typename StencilType, typename FaceContainerType>
	int reconstruct(unsigned int cellIdx, unsigned int numCells, const StencilType& w, StateType& result, double* const Dvm, const FaceContainerType& cellFaces)
	{
		return reconstructNonEq<StateType, StencilType, FaceContainerType, true>(cellIdx, numCells, w, result, Dvm, cellFaces);
	}

	template <typename StateType, typename StencilType, typename FaceContainerType>
	int reconstruct(unsigned int cellIdx, unsigned int numCells, const StencilType& w, StateType& result, const FaceContainerType& cellFaces)
	{
		return reconstructNonEq<StateType, StencilType, FaceContainerType, false>(cellIdx, numCells, w, result, nullptr, cellFaces);
	}

	/**
	 * @brief Reconstructs a cell face value from geometry-weighted (radial) volume averages
	 * @details Geometry-exact WENO for radial (cylinder shell) flow on equidistant grids.
	 *          The FV degrees of freedom are the rho-weighted averages
	 *          \f$ \bar c_j = \frac{2}{\rho_{j+1/2}^2 - \rho_{j-1/2}^2} \int_{\rho_{j-1/2}}^{\rho_{j+1/2}} \rho c \,\mathrm{d}\rho. \f$
	 *          All stencil weights, optimal weights, and smoothness indicators are closed-form rational
	 *          functions of the single dimensionless parameter \f$ \zeta = \rho_{i+1/2} / \Delta\rho \f$
	 *          (outer face radius of the current cell divided by the cell width). They reproduce
	 *          rho-weighted averages of polynomials exactly (degree 2 per substencil, degree 4 on the
	 *          full stencil) and reduce to the classical Jiang-Shu coefficients as \f$ \zeta \to \infty \f$.
	 *          See Shadab et al., Computers & Fluids 190 (2019) 398-424, for the general framework.
	 *
	 *          Boundary cells always use order reduction (the radial optimal weights do not exist for
	 *          stencils that protrude the domain), irrespective of the configured boundary treatment.
	 *
	 * @param [in] cellIdx Index of the current cell (in flow direction)
	 * @param [in] numCells Number of cells
	 * @param [in] zeta Outer face radius of the current cell divided by the cell width, \f$ \rho_{i+1/2} / \Delta\rho \f$
	 * @param [in] forwardFlow @c true for outward flow (reconstruction at the outer face \f$ \rho_{i+1/2} \f$),
	 *             @c false for inward flow (reconstruction at the inner face \f$ \rho_{i-1/2} \f$)
	 * @param [in] w Stencil in flow orientation (index 0 is the current cell, positive indices are upstream
	 *             of the reconstructed face as in the axial scheme)
	 * @param [out] result Reconstructed cell face value
	 * @param [out] Dvm Gradient of the reconstructed cell face value (array has to be of size \f$ 2r-1 \f$)
	 * @tparam StateType Type of the state variables
	 * @tparam StencilType Type of the stencil (can be a dedicated class with overloaded operator[] or a simple pointer)
	 * @return Order of the WENO scheme that was used in the computation
	 */
	template <typename StateType, typename StencilType>
	int reconstructRadial(unsigned int cellIdx, unsigned int numCells, double zeta, bool forwardFlow, const StencilType& w, StateType& result, double* const Dvm)
	{
		return reconstructRadial<StateType, StencilType, true>(cellIdx, numCells, zeta, forwardFlow, w, result, Dvm);
	}

	/**
	 * @brief Reconstructs a cell face value from geometry-weighted (radial) volume averages
	 * @details Same as the Jacobian variant of reconstructRadial(), but does not compute the gradient.
	 */
	template <typename StateType, typename StencilType>
	int reconstructRadial(unsigned int cellIdx, unsigned int numCells, double zeta, bool forwardFlow, const StencilType& w, StateType& result)
	{
		return reconstructRadial<StateType, StencilType, false>(cellIdx, numCells, zeta, forwardFlow, w, result, nullptr);
	}

	/**
	 * @brief Precomputes geometry-exact WENO coefficients for a general weighted geometry and grid
	 * @details The FV degrees of freedom are the geometry-weighted cell averages
	 *          \f$ \bar c_j = \int_{x_{j-1/2}}^{x_{j+1/2}} A(x) c \,\mathrm{d}x / \int_{x_{j-1/2}}^{x_{j+1/2}} A(x) \,\mathrm{d}x \f$
	 *          with the cross-section weight \f$ A(x) = (a + b x)^m \f$ (radial cylinder shell:
	 *          \f$ a = 0, b = 1, m = 1 \f$ with \f$ x = \rho \f$; frustum: \f$ a = r_\text{in},
	 *          b = (r_\text{out} - r_\text{in})/L, m = 2 \f$ with \f$ x = z \f$; \f$ m = 0 \f$
	 *          reduces to the axial scheme). For every cell and both flow directions, the candidate
	 *          stencil weights, optimal weights, and smoothness-indicator quadratic forms are
	 *          computed from the moment systems of the weighted averaging operator (evaluated in
	 *          well-conditioned local coordinates) and stored. Works on arbitrary (also
	 *          non-equidistant) grids and is exact for the given grid.
	 *
	 *          The moment systems reproduce weighted averages of polynomials exactly (degree
	 *          \f$ r-1 \f$ per substencil, degree \f$ 2r-2 \f$ on the full stencil). The optimal
	 *          weights are obtained from the overdetermined matching system; all residual equations
	 *          are checked. If the optimal weights do not exist or are not positive on a cell
	 *          (possible on strongly distorted grids), the order is reduced for that cell.
	 *          Boundary cells always use order reduction.
	 *
	 *          Must be called after order() has been set and whenever the grid changes.
	 *
	 * @param [in] a Constant coefficient of the weight function
	 * @param [in] b Linear coefficient of the weight function
	 * @param [in] weightPower Exponent \f$ m \in \{0, 1, 2\} \f$ of the weight function
	 * @param [in] faces Cell face coordinates (ascending, size number of cells + 1)
	 */
	void prepareGeometryExactCoefficients(double a, double b, int weightPower, const std::vector<double>& faces);

	/**
	 * @brief Returns whether geometry-exact coefficients have been precomputed
	 */
	inline bool hasGeometryExactCoefficients() const CADET_NOEXCEPT { return !_geomForward.empty(); }

	/**
	 * @brief Reconstructs a cell face value from geometry-weighted volume averages using precomputed coefficients
	 * @details Geometry-exact WENO reconstruction based on the per-cell coefficients precomputed by
	 *          prepareGeometryExactCoefficients(). The reconstructed face is the downwind face of the
	 *          given cell (outer/right face for forward flow, inner/left face for backward flow).
	 * @param [in] cellIdx Physical index of the current cell (i.e., cell counted from the first face,
	 *             independent of the flow direction)
	 * @param [in] numCells Number of cells
	 * @param [in] forwardFlow @c true for flow in the direction of increasing coordinate
	 * @param [in] w Stencil in flow orientation (index 0 is the current cell, positive indices are
	 *             upstream of the reconstructed face as in the axial scheme)
	 * @param [out] result Reconstructed cell face value
	 * @param [out] Dvm Gradient of the reconstructed cell face value (array has to be of size \f$ 2r-1 \f$)
	 * @return Order of the WENO scheme that was used in the computation
	 */
	template <typename StateType, typename StencilType>
	int reconstructGeomExact(unsigned int cellIdx, unsigned int numCells, bool forwardFlow, const StencilType& w, StateType& result, double* const Dvm)
	{
		return reconstructGeomExact<StateType, StencilType, true>(cellIdx, numCells, forwardFlow, w, result, Dvm);
	}

	/**
	 * @brief Reconstructs a cell face value from geometry-weighted volume averages using precomputed coefficients
	 * @details Same as the Jacobian variant of reconstructGeomExact(), but does not compute the gradient.
	 */
	template <typename StateType, typename StencilType>
	int reconstructGeomExact(unsigned int cellIdx, unsigned int numCells, bool forwardFlow, const StencilType& w, StateType& result)
	{
		return reconstructGeomExact<StateType, StencilType, false>(cellIdx, numCells, forwardFlow, w, result, nullptr);
	}

	/**
	 * @brief Sets the \f$ \varepsilon \f$ of the WENO emthod (prevents division by zero in the weights)
	 * @param [in] eps Boundary treatment method
	 */
	inline void epsilon(double eps) { _epsilon = eps; }

	/**
	 * @brief Returns the \f$ \varepsilon \f$ of the WENO emthod (prevents division by zero in the weights)
	 * @return WENO \f$ \varepsilon \f$
	 */
	inline double epsilon() const CADET_NOEXCEPT { return _epsilon; }

	/**
	 * @brief Sets the WENO order
	 * @param [in] order Order of the WENO method
	 */
	inline void order(int order)
	{
		if (order <= 0 || order > static_cast<int>(maxOrder()))
			throw InvalidParameterException("WENO Order specified as " + std::to_string(order) + " but must be > 0 and <= " + std::to_string(static_cast<int>(maxOrder())));
		_order = order;
	}
	
	/**
	 * @brief Returns the WENO order
	 * @return Order of the WENO method
	 */
	inline int order() const CADET_NOEXCEPT { return _order; }

	/**
	 * @brief Sets the boundary treatment method
	 * @param [in] bndTreatment Boundary treatment method
	 */
	inline void boundaryTreatment(BoundaryTreatment bndTreatment) { _boundaryTreatment = bndTreatment; }

	/**
	 * @brief Returns the boundary treatment method
	 * @return Boundary treatment method
	 */
	inline BoundaryTreatment boundaryTreatment() const CADET_NOEXCEPT { return _boundaryTreatment; }

	/**
	 * @brief Sets the boundary treatment method
	 * @param [in] bndTreatment Boundary treatment method, is converted to the BoundaryTreatment enum
	 */
	inline void boundaryTreatment(int bndTreatment)
	{
		switch(bndTreatment)
		{
			case static_cast<typename std::underlying_type<BoundaryTreatment>::type>(BoundaryTreatment::ReduceOrder):
				_boundaryTreatment = BoundaryTreatment::ReduceOrder;
				return;
			case static_cast<typename std::underlying_type<BoundaryTreatment>::type>(BoundaryTreatment::ZeroWeights):
				_boundaryTreatment = BoundaryTreatment::ZeroWeights;
				return;
			case static_cast<typename std::underlying_type<BoundaryTreatment>::type>(BoundaryTreatment::ZeroWeightsForPnotZero):
				_boundaryTreatment = BoundaryTreatment::ZeroWeightsForPnotZero;
				return;
			case static_cast<typename std::underlying_type<BoundaryTreatment>::type>(BoundaryTreatment::LargeGhostNodes):
				_boundaryTreatment = BoundaryTreatment::LargeGhostNodes;
				return;
		}
		throw InvalidParameterException("Unknown boundary treatment type");
	}

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
	 * @return Size of the stencil
	 */
	inline unsigned int stencilSize() const CADET_NOEXCEPT { return 2 * _order - 1; }

private:

	/**
	 * @brief Reconstructs a cell face value from volume averages
	 * @param [in] cellIdx Index of the current cell
	 * @param [in] numCells Number of cells
	 * @param [in] w Stencil that contains the \f$ 2r-1 \f$ volume averages from which the cell face values are reconstructed centered at the 
	 *               current cell (i.e., index 0 is the current cell, -2 the next to previous cell, 2 the next but one cell)
	 * @param [out] result Reconstructed cell face value
	 * @param [out] Dvm Gradient of the reconstructed cell face value (array has to be of size \f$ 2r-1\f$ where \f$ r \f$ is the WENO order)
	 * @tparam StateType Type of the state variables
	 * @tparam StencilType Type of the stencil (can be a dedicated class with overloaded operator[] or a simple pointer)
	 * @tparam wantJac Determines if the gradient is computed (@c true) or not (@c false)
	 * @return Order of the WENO scheme that was used in the computation
	 */

	// Reconstruct for equidistant grids
	template <typename StateType, typename StencilType, bool wantJac>
	int reconstruct(unsigned int cellIdx, unsigned int numCells, const StencilType& w, StateType& result, double* const Dvm)
	{
#if defined(ACTIVE_SETFAD) || defined(ACTIVE_SFAD)
		using cadet::sqr;
		using sfad::sqr;
#endif

		// Local order of the scheme that is actually used (may be changed by treatment of boundaries)
		int order = _order;

		// Boundaries
		int bnd = 0;
		switch (_boundaryTreatment)
		{
		default:
		case BoundaryTreatment::ReduceOrder:
			// Lower WENO order such that maximum order is used at all points
			// This very statement selects the max. weno order for the current column cell
			// order = min(maxOrderleft, maxOrderright)
			order = std::min(std::min(static_cast<int>(cellIdx) + 1, _order), std::min(static_cast<int>(numCells - cellIdx), _order));
			break;

		case BoundaryTreatment::ZeroWeights:
			// Zero weights for cell averages outside the domain
			if (cellIdx < static_cast<unsigned int>(order - 1))
				bnd = -(order - 1 - cellIdx);
			else if (cellIdx > numCells - order)
				bnd = numCells - cellIdx;
			break;

		case BoundaryTreatment::ZeroWeightsForPnotZero:
			// Zero weights for p != 0
			if (cellIdx == 0)
				order = 1;
			else
			{
				if (cellIdx < static_cast<unsigned int>(order - 1))
					bnd = -(order - 1 - cellIdx);
				else if (cellIdx > numCells - order)
					bnd = numCells - cellIdx;
			}
			break;
/*
		case BoundaryTreatment::LargeGhostNodes:
			// Large ghost points
			if (cellIdx == 0)
			{
				w[-1] = 1e20;
				w[-2] = 1e50;
			}
			else if (cellIdx == numCells - 2)
				w[2] = 1e20;
			else if (cellIdx == numCells - 1)
				w[2] = 1e50;
			break;
*/
		}		

		// Total stencil size
		const int sl = 2 * order - 1;

		// Simple upwind scheme
		if (order == 1)
		{
			result = w[0];
			if (wantJac)
				*Dvm = 1.0;
			return order;
		}

		// Allocate memory for intermediate values: beta, alpha (= omega), and vr
		StateType* const work = _intermediateValues.create<StateType>(3 * order);
		StateType* const beta  = work;
		StateType* const alpha = work + order;
		StateType* const omega = work + order;
		StateType* const vr    = work + 2*order; // Reconstructed values

		const double* d = nullptr;
		const double* c = nullptr;
		const double* Jbvv = nullptr;

		// Calculate smoothness measures
		switch (order)
		{
			case 2:
				beta[0] = sqr(w[1] - w[0]);
				beta[1] = sqr(w[0] - w[-1]);
				d = _wenoD2;
				c = _wenoC2;
				Jbvv = _wenoJbvv2;
				break;
			case 3:
				beta[0] = 13.0/12.0 * sqr(w[ 0] - 2.0 * w[ 1] + w[2]) + 0.25 * sqr(3.0 * w[ 0] - 4.0 * w[ 1] +       w[2]);
				beta[1] = 13.0/12.0 * sqr(w[-1] - 2.0 * w[ 0] + w[1]) + 0.25 * sqr(      w[-1] -       w[ 1]             );
				beta[2] = 13.0/12.0 * sqr(w[-2] - 2.0 * w[-1] + w[0]) + 0.25 * sqr(      w[-2] - 4.0 * w[-1] + 3.0 * w[0]);
				d = _wenoD3;
				c = _wenoC3;
				Jbvv = _wenoJbvv3;
				break;
		}

		// Add eps to avoid divide-by-zeros
		for (int r = 0; r < order; ++r)
			beta[r] += _epsilon;

		// Calculate weights
		for (int r = 0; r < order; ++r)
			alpha[r] = d[r] / sqr(beta[r]);

		// Avoid boundaries
		if (cadet_unlikely(bnd != 0))
		{
			if (bnd < 0)
				// Beginning of interval
				for (int r = 0; r < -bnd; ++r)
					alpha[order - 1 - r] = 0.0;
			else
				// End of interval
				for (int r = 0; r < bnd; ++r)
					alpha[r] = 0.0;
		}

		// Normalize weights
		StateType alpha_sum = alpha[0];
		for (int r = 1; r < order; ++r)
			alpha_sum += alpha[r];
		for (int r = 0; r < order; ++r)
			omega[r] /= alpha_sum;

		// Calculate reconstructed values
		for (int r = 0; r < order; ++r)
		{
			vr[r] = 0.0;
			for (int j = 0; j < order; ++j)
				vr[r] += c[r + order * j] * w[-r+j];
		}

		// Weighted sum
		result = 0;
		for (int r = 0; r < order; ++r)
			result += vr[r] * omega[r];

		// Jacobian
		if (wantJac)
		{
			// Dependencies
			// 1. Constant vr in (*)

			// Start with "d(result)/d(omega)" = vr and
			// multiply with "d(omega)/d(alpha)" to get "d(result)/d(alpha)"
			double dot = 0.0;
			for (int r = 0; r < order; ++r)
				dot += static_cast<double>(vr[r]) * static_cast<double>(omega[r]); //StateType(vr[r] * omega[r]);
			for (int r = 0; r < order; ++r)
				vr[r] = (vr[r] - dot) / alpha_sum;

			// Multiply with "d(alpha)/d(beta)" to get "d(result)/d(beta)"
			for (int r = 0; r < order; ++r)
				vr[r] *= -2.0 * d[r] / pow(beta[r], 3.0);

			// Multiply with "d(beta)/d(v)" to get Dvm = "d(result)/d(v)"
			for (int j = 0; j < sl; ++j)
			{
				Dvm[j] = 0.0;
				for (int r = 0; r < order; ++r)
				{
					dot = 0.0;
					for (int i = 0; i < sl; ++i)
						dot += static_cast<double>(Jbvv[r + order * j + order * sl * i]) * static_cast<double>(w[i - order + 1]);
					// To do: re-arange Jbvv to reduce cache misses !
					Dvm[j] += static_cast<double>(vr[r]) * dot; // StateType(vr[r] * dot);
				}
			}

			// 2. Constant omega[r] in (*)
			for (int r = 0; r < order; ++r)
				for (int j = 0; j < order; ++j)
					Dvm[order - 1 + j - r] += static_cast<double>(omega[r]) * c[r + order * j];
		}

		_intermediateValues.destroy<StateType>();
		return order;
	}

	template <typename FaceContainerType>
	static inline auto cellWidths3FromFaces(unsigned int cellIdx, unsigned int numCells, const FaceContainerType& cellFaces)
	{
		const auto hCol = cellFaces[cellIdx + 1] - cellFaces[cellIdx];
		auto hLeft = hCol;
		auto hRight = hCol;
		if (cellIdx > 0)
			hLeft = cellFaces[cellIdx] - cellFaces[cellIdx - 1];
		if (cellIdx + 1 < numCells)
			hRight = cellFaces[cellIdx + 2] - cellFaces[cellIdx + 1];
		return std::array<decltype(hCol), 3>{hLeft, hCol, hRight};
	}

	template <typename FaceContainerType>
	static inline auto cellWidths5FromFaces(unsigned int cellIdx, unsigned int numCells, const FaceContainerType& cellFaces)
	{
		const auto hCol = cellFaces[cellIdx + 1] - cellFaces[cellIdx];
		auto hm2 = hCol;
		auto hm1 = hCol;
		auto hp1 = hCol;
		auto hp2 = hCol;

		if (cellIdx > 1)
			hm2 = cellFaces[cellIdx - 1] - cellFaces[cellIdx - 2];
		if (cellIdx > 0)
			hm1 = cellFaces[cellIdx] - cellFaces[cellIdx - 1];
		if (cellIdx + 1 < numCells)
			hp1 = cellFaces[cellIdx + 2] - cellFaces[cellIdx + 1];
		if (cellIdx + 2 < numCells)
			hp2 = cellFaces[cellIdx + 3] - cellFaces[cellIdx + 2];

		return std::array<decltype(hCol), 5>{hm2, hm1, hCol, hp1, hp2};
	}


	template <typename StateType, typename StencilType, typename FaceContainerType, bool wantJac>
	int reconstructNonEq(unsigned int cellIdx, unsigned int numCells, const StencilType& w, StateType& result, double* const Dvm, const FaceContainerType& cellFaces)
	{
		int order = std::min(std::min(static_cast<int>(cellIdx) + 1, _order), std::min(static_cast<int>(numCells - cellIdx), _order));
		order = std::max(order, 1);
		if (order <= 1)
		{
			result = StateType(w[0]);
			if (wantJac)
				Dvm[0] = 1.0;
			return order;
		}

		if (order == 2)
		{
			const auto hRaw = cellWidths3FromFaces(cellIdx, numCells, cellFaces);
			const StateType hm1 = static_cast<StateType>(hRaw[0]);
			const StateType h0 = static_cast<StateType>(hRaw[1]);
			const StateType hp1 = static_cast<StateType>(hRaw[2]);
			const StateType wm1 = static_cast<StateType>(w[-1]);
			const StateType w0 = static_cast<StateType>(w[0]);
			const StateType wp1 = static_cast<StateType>(w[1]);
			const auto hs = hm1 + h0 + hp1;
			const auto hl = hm1 + h0;
			const auto hr = h0 + hp1;
			const auto q0c = hp1 / hr;
			const auto q1c = 1.0 + h0 / hl;
			const auto C0 = hl / hs;
			const auto C1 = 1.0 - C0;

			const auto d0 = wp1 - w0;
			const auto d1 = w0 - wm1;
			const auto IS0c = (2.0 * h0 / hr) * (2.0 * h0 / hr);
			const auto IS1c = (2.0 * h0 / hl) * (2.0 * h0 / hl);
			const auto IS0 = IS0c * d0 * d0;
			const auto IS1 = IS1c * d1 * d1;
			const auto a0 = C0 / ((_epsilon + IS0) * (_epsilon + IS0));
			const auto a1 = C1 / ((_epsilon + IS1) * (_epsilon + IS1));
			const auto W0 = a0 / (a0 + a1);
			const auto W1 = 1.0 - W0;

			const auto q0 = q0c * w0 + (1.0 - q0c) * wp1;
			const auto q1 = q1c * w0 + (1.0 - q1c) * wm1;
			result = StateType(W0 * q0 + W1 * q1);

			if (wantJac)
			{
				const auto dq0_d0 = q0c;
				const auto dq0_dp1 = 1.0 - q0c;
				const auto dq1_d0 = q1c;
				const auto dq1_dm1 = 1.0 - q1c;

				const auto dIS0_d0 = -2.0 * IS0c * d0;
				const auto dIS0_dp1 = 2.0 * IS0c * d0;
				const auto dIS1_d0 = 2.0 * IS1c * d1;
				const auto dIS1_dm1 = -2.0 * IS1c * d1;

				const auto da0_d0 = -2.0 * C0 * (1.0 / ((_epsilon + IS0) * (_epsilon + IS0) * (_epsilon + IS0))) * dIS0_d0;
				const auto da0_dp1 = -2.0 * C0 * (1.0 / ((_epsilon + IS0) * (_epsilon + IS0) * (_epsilon + IS0))) * dIS0_dp1;
				const auto da1_d0 = -2.0 * C1 * (1.0 / ((_epsilon + IS1) * (_epsilon + IS1) * (_epsilon + IS1))) * dIS1_d0;
				const auto da1_dm1 = -2.0 * C1 * (1.0 / ((_epsilon + IS1) * (_epsilon + IS1) * (_epsilon + IS1))) * dIS1_dm1;

				const auto aSum = a0 + a1;
				const auto dW0_d0 = (da0_d0 * aSum - a0 * (da0_d0 + da1_d0)) / (aSum * aSum);
				const auto dW0_dp1 = (da0_dp1 * aSum - a0 * da0_dp1) / (aSum * aSum);
				const auto dW0_dm1 = (-a0 * da1_dm1) / (aSum * aSum);
				const auto dW1_d0 = -dW0_d0;
				const auto dW1_dp1 = -dW0_dp1;
				const auto dW1_dm1 = -dW0_dm1;

				Dvm[0] = static_cast<double>(dW0_dm1 * q0 + W0 * 0.0 + dW1_dm1 * q1 + W1 * dq1_dm1);
				Dvm[1] = static_cast<double>(dW0_d0 * q0 + W0 * dq0_d0 + dW1_d0 * q1 + W1 * dq1_d0);
				Dvm[2] = static_cast<double>(dW0_dp1 * q0 + W0 * dq0_dp1 + dW1_dp1 * q1 + W1 * 0.0);
			}
			return order;
		}

		const auto hRaw = cellWidths5FromFaces(cellIdx, numCells, cellFaces);

		const StateType hm2 = static_cast<StateType>(hRaw[0]);
		const StateType hm1 = static_cast<StateType>(hRaw[1]);
		const StateType h0 = static_cast<StateType>(hRaw[2]);
		const StateType hp1 = static_cast<StateType>(hRaw[3]);
		const StateType hp2 = static_cast<StateType>(hRaw[4]);

		const StateType wm2 = static_cast<StateType>(w[-2]);
		const StateType wm1 = static_cast<StateType>(w[-1]);
		const StateType w0 = static_cast<StateType>(w[0]);
		const StateType wp1 = static_cast<StateType>(w[1]);
		const StateType wp2 = static_cast<StateType>(w[2]);

		const auto dsum = hm2 + hm1 + h0 + hp1 + hp2;
		const auto dI0 = hm2 + hm1 + h0;
		const auto dI1 = hm1 + h0 + hp1;
		const auto dI2 = h0 + hp1 + hp2;

		const auto q0c1 = hp1 * (dI2 - h0) / (dI2 - hp2) / dI2;
		const auto q0c2 = hp1 * h0 / dI2 / (dI2 - h0);
		const auto q1c1 = h0 * (dI1 - hp1) / dI1 / (dI1 - hm1);
		const auto q1c2 = h0 * hp1 / dI1 / (dI1 - hp1);
		const auto q2c1 = h0 * (dI0 - hm2) / dI0 / (dI0 - h0);
		const auto q2c2 = 1.0 + h0 / (hm1 + h0) + h0 / dI0;

		const auto C0 = dI0 * (dI0 - hm2) / dsum / (dsum - hm2);
		const auto C1 = dI0 / dsum * (dI2 - h0) / (hm1 + dI2) * (1.0 + (dsum - hm2) / (dsum - hp2));
		const auto C2 = hp1 * (hp1 + hp2) / dsum / (dsum - hp2);

		const auto IS0pre = 4.0 * (h0 / dI2) * (h0 / dI2);
		const auto IS0c1 = IS0pre * (10.0 * h0 * h0 + hp1 * (h0 + hp1)) / ((hp1 + hp2) * (hp1 + hp2));
		const auto IS0c2 = IS0pre * (20.0 * h0 * h0 + 2.0 * hp1 * (h0 + hp1) + (2.0 * hp1 + h0) * dI2) / ((hp1 + hp2) * (h0 + hp1));
		const auto IS0c3 = IS0pre * (10.0 * h0 * h0 + (2.0 * dI2 - hp2) * (dI2 + hp1)) / ((h0 + hp1) * (h0 + hp1));

		const auto IS1pre = 4.0 * (h0 / dI1) * (h0 / dI1);
		const auto IS1c1 = IS1pre * (10.0 * h0 * h0 + hp1 * (h0 + hp1)) / ((hm1 + h0) * (hm1 + h0));
		const auto IS1c2 = IS1pre * (20.0 * h0 * h0 - hp1 * hm1 - (h0 + hp1) * (h0 + hm1)) / ((h0 + hp1) * (h0 + hm1));
		const auto IS1c3 = IS1pre * (10.0 * h0 * h0 + hm1 * (hm1 + h0)) / ((h0 + hp1) * (h0 + hp1));

		const auto IS2pre = 4.0 * (h0 / dI0) * (h0 / dI0);
		const auto IS2c1 = IS2pre * (10.0 * h0 * h0 + hm1 * (hm1 + h0)) / ((hm2 + hm1) * (hm2 + hm1));
		const auto IS2c2 = IS2pre * (20.0 * h0 * h0 + 2.0 * hm1 * (hm1 + h0) + dI0 * (2.0 * hm1 + h0)) / ((hm1 + h0) * (hm2 + hm1));
		const auto IS2c3 = IS2pre * (10.0 * h0 * h0 + (2.0 * dI0 - hm2) * (dI0 + hm1)) / ((hm1 + h0) * (hm1 + h0));

		const auto dm2 = wm2;
		const auto dm1 = wm1;
		const auto d0 = w0;
		const auto dp1 = wp1;
		const auto dp2 = wp2;

		const auto IS0 = IS0c1 * (dp2 - dp1) * (dp2 - dp1) + IS0c2 * (dp2 - dp1) * (d0 - dp1) + IS0c3 * (d0 - dp1) * (d0 - dp1);
		const auto IS1 = IS1c1 * (dm1 - d0) * (dm1 - d0) + IS1c2 * (dp1 - d0) * (dm1 - d0) + IS1c3 * (dp1 - d0) * (dp1 - d0);
		const auto IS2 = IS2c1 * (dm2 - dm1) * (dm2 - dm1) + IS2c2 * (d0 - dm1) * (dm2 - dm1) + IS2c3 * (d0 - dm1) * (d0 - dm1);

		const auto a0 = C0 / ((_epsilon + IS0) * (_epsilon + IS0));
		const auto a1 = C1 / ((_epsilon + IS1) * (_epsilon + IS1));
		const auto a2 = C2 / ((_epsilon + IS2) * (_epsilon + IS2));
		const auto aSum = a0 + a1 + a2;
		const auto W0 = a0 / aSum;
		const auto W1 = a1 / aSum;
		const auto W2 = 1.0 - W0 - W1;

		const auto q0 = dp1 + q0c1 * (d0 - dp1) + q0c2 * (dp1 - dp2);
		const auto q1 = d0 + q1c1 * (dp1 - d0) + q1c2 * (d0 - dm1);
		const auto q2 = dm1 + q2c1 * (dm2 - dm1) + q2c2 * (d0 - dm1);
		result = StateType(W0 * q0 + W1 * q1 + W2 * q2);

		if (wantJac)
		{
			const auto dIS0_dp2 = 2.0 * IS0c1 * (dp2 - dp1) + IS0c2 * (d0 - dp1);
			const auto dIS0_dp1 = -2.0 * IS0c1 * (dp2 - dp1) + IS0c2 * (2.0 * dp1 - dp2 - d0) - 2.0 * IS0c3 * (d0 - dp1);
			const auto dIS0_d0 = IS0c2 * (dp2 - dp1) + 2.0 * IS0c3 * (d0 - dp1);

			const auto dIS1_dp1 = IS1c2 * (dm1 - d0) + 2.0 * IS1c3 * (dp1 - d0);
			const auto dIS1_d0 = -2.0 * IS1c1 * (dm1 - d0) + IS1c2 * (2.0 * d0 - dp1 - dm1) - 2.0 * IS1c3 * (dp1 - d0);
			const auto dIS1_dm1 = 2.0 * IS1c1 * (dm1 - d0) + IS1c2 * (dp1 - d0);

			const auto dIS2_d0 = IS2c2 * (dm2 - dm1) + 2.0 * IS2c3 * (d0 - dm1);
			const auto dIS2_dm1 = -2.0 * IS2c1 * (dm2 - dm1) + IS2c2 * (2.0 * dm1 - d0 - dm2) - 2.0 * IS2c3 * (d0 - dm1);
			const auto dIS2_dm2 = 2.0 * IS2c1 * (dm2 - dm1) + IS2c2 * (d0 - dm1);

			const auto da0_dp2 = -2.0 * C0 * (1.0 / ((_epsilon + IS0) * (_epsilon + IS0) * (_epsilon + IS0))) * dIS0_dp2;
			const auto da0_dp1 = -2.0 * C0 * (1.0 / ((_epsilon + IS0) * (_epsilon + IS0) * (_epsilon + IS0))) * dIS0_dp1;
			const auto da0_d0 = -2.0 * C0 * (1.0 / ((_epsilon + IS0) * (_epsilon + IS0) * (_epsilon + IS0))) * dIS0_d0;

			const auto da1_dp1 = -2.0 * C1 * (1.0 / ((_epsilon + IS1) * (_epsilon + IS1) * (_epsilon + IS1))) * dIS1_dp1;
			const auto da1_d0 = -2.0 * C1 * (1.0 / ((_epsilon + IS1) * (_epsilon + IS1) * (_epsilon + IS1))) * dIS1_d0;
			const auto da1_dm1 = -2.0 * C1 * (1.0 / ((_epsilon + IS1) * (_epsilon + IS1) * (_epsilon + IS1))) * dIS1_dm1;

			const auto da2_d0 = -2.0 * C2 * (1.0 / ((_epsilon + IS2) * (_epsilon + IS2) * (_epsilon + IS2))) * dIS2_d0;
			const auto da2_dm1 = -2.0 * C2 * (1.0 / ((_epsilon + IS2) * (_epsilon + IS2) * (_epsilon + IS2))) * dIS2_dm1;
			const auto da2_dm2 = -2.0 * C2 * (1.0 / ((_epsilon + IS2) * (_epsilon + IS2) * (_epsilon + IS2))) * dIS2_dm2;

			const auto zero = aSum * 0.0;
			const auto invASum2 = 1.0 / (aSum * aSum);

			const auto sumDa_dp2 = da0_dp2 + zero + zero;
			const auto sumDa_dp1 = da0_dp1 + da1_dp1 + zero;
			const auto sumDa_d0 = da0_d0 + da1_d0 + da2_d0;
			const auto sumDa_dm1 = zero + da1_dm1 + da2_dm1;
			const auto sumDa_dm2 = zero + zero + da2_dm2;

			const auto dW0_dp2 = (da0_dp2 * aSum - a0 * sumDa_dp2) * invASum2;
			const auto dW0_dp1 = (da0_dp1 * aSum - a0 * sumDa_dp1) * invASum2;
			const auto dW0_d0 = (da0_d0 * aSum - a0 * sumDa_d0) * invASum2;
			const auto dW0_dm1 = (zero * aSum - a0 * sumDa_dm1) * invASum2;
			const auto dW0_dm2 = (zero * aSum - a0 * sumDa_dm2) * invASum2;

			const auto dW1_dp2 = (zero * aSum - a1 * sumDa_dp2) * invASum2;
			const auto dW1_dp1 = (da1_dp1 * aSum - a1 * sumDa_dp1) * invASum2;
			const auto dW1_d0 = (da1_d0 * aSum - a1 * sumDa_d0) * invASum2;
			const auto dW1_dm1 = (da1_dm1 * aSum - a1 * sumDa_dm1) * invASum2;
			const auto dW1_dm2 = (zero * aSum - a1 * sumDa_dm2) * invASum2;

			const auto dW2_dp2 = -dW0_dp2 - dW1_dp2;
			const auto dW2_dp1 = -dW0_dp1 - dW1_dp1;
			const auto dW2_d0 = -dW0_d0 - dW1_d0;
			const auto dW2_dm1 = -dW0_dm1 - dW1_dm1;
			const auto dW2_dm2 = -dW0_dm2 - dW1_dm2;

			const auto dq0_dp2 = -q0c2;
			const auto dq0_dp1 = 1.0 - q0c1 + q0c2;
			const auto dq0_d0 = q0c1;

			const auto dq1_dp1 = q1c1;
			const auto dq1_d0 = 1.0 - q1c1 + q1c2;
			const auto dq1_dm1 = -q1c2;

			const auto dq2_d0 = q2c2;
			const auto dq2_dm1 = 1.0 - q2c1 - q2c2;
			const auto dq2_dm2 = q2c1;

			Dvm[0] = static_cast<double>(dW0_dm2 * q0 + dW1_dm2 * q1 + dW2_dm2 * q2 + W2 * dq2_dm2);
			Dvm[1] = static_cast<double>(dW0_dm1 * q0 + dW1_dm1 * q1 + dW2_dm1 * q2 + W1 * dq1_dm1 + W2 * dq2_dm1);
			Dvm[2] = static_cast<double>(dW0_d0 * q0 + dW1_d0 * q1 + dW2_d0 * q2 + W0 * dq0_d0 + W1 * dq1_d0 + W2 * dq2_d0);
			Dvm[3] = static_cast<double>(dW0_dp1 * q0 + dW1_dp1 * q1 + dW2_dp1 * q2 + W0 * dq0_dp1 + W1 * dq1_dp1);
			Dvm[4] = static_cast<double>(dW0_dp2 * q0 + dW1_dp2 * q1 + dW2_dp2 * q2 + W0 * dq0_dp2);
		}
		return order;
	}


	/**
	 * @brief Evaluates the radial WENO coefficients of order 3 (fifth order scheme) at a given \f$ \zeta \f$
	 * @details Fills the value weights @p C, the optimal weights @p D, and the smoothness indicator
	 *          quadratic forms @p B in flow-oriented substencil indexing: substencil \f$ r \f$ uses the
	 *          stencil values \f$ u_j = w[j - r] \f$, \f$ j = 0, \dots, 2 \f$, and
	 *          \f$ q_r = \sum_j C[r][j] u_j \f$, \f$ \beta_r = \sum_{m,n} B[r][m][n] u_m u_n \f$.
	 *          All formulas are closed-form rational functions of \f$ \zeta = \rho_{i+1/2}/\Delta\rho \f$,
	 *          derived from exact reproduction of rho-weighted cell averages (symbolically verified;
	 *          they agree with Shadab et al. 2019, Appendix A.2, and reduce to the Jiang-Shu
	 *          coefficients for \f$ \zeta \to \infty \f$).
	 */
	static void radialCoefficientsOrder3(double zeta, bool forwardFlow, double C[3][3], double D[3], double B[3][3][3])
	{
		const double z = zeta;
		const double z2 = z * z;
		const double z3 = z2 * z;
		const double z4 = z2 * z2;

		// Recurring factors
		const double f2m5 = 2.0 * z - 5.0;
		const double f2m3 = 2.0 * z - 3.0;
		const double f2m1 = 2.0 * z - 1.0;
		const double f2p1 = 2.0 * z + 1.0;
		const double f2p3 = 2.0 * z + 3.0;
		const double pm3 = z2 - 3.0 * z + 1.0;  // vanishes only for zeta < 3 (root at 2.618...)
		const double pm1 = z2 - z - 1.0;
		const double pp1 = z2 + z - 1.0;
		const double p5 = 3.0 * z4 - 6.0 * z3 - 13.0 * z2 + 16.0 * z + 12.0;

		// Substencil value weights at the reconstructed face (physical formulas):
		// S0 = {i-2, i-1, i}, S1 = {i-1, i, i+1}, S2 = {i, i+1, i+2}
		// Outer face rho_{i+1/2} ("+" family, used for forward/outward flow):
		// wP[k][m] multiplies cell i - 2 + k + m (physically ascending order within stencil k)
		const double wP0[3] = {
			 f2m5 * (4.0 * z2 - 9.0 * z + 4.0) / (12.0 * f2m3 * pm3),
			-(14.0 * z2 - 45.0 * z + 23.0) / (12.0 * pm3),
			 f2m1 * (22.0 * z2 - 90.0 * z + 85.0) / (12.0 * f2m3 * pm3)
		};
		const double wP1[3] = {
			-f2m3 * (2.0 * z2 - 1.0) / (12.0 * f2m1 * pm1),
			 (10.0 * z2 - 9.0 * z - 11.0) / (12.0 * pm1),
			 f2p1 * (4.0 * z2 - 9.0 * z + 4.0) / (12.0 * f2m1 * pm1)
		};
		const double wP2[3] = {
			 f2m1 * (4.0 * z2 + 9.0 * z + 4.0) / (12.0 * f2p1 * pp1),
			 (10.0 * z2 + 9.0 * z - 11.0) / (12.0 * pp1),
			-f2p3 * (2.0 * z2 - 1.0) / (12.0 * f2p1 * pp1)
		};
		// Inner face rho_{i-1/2} ("-" family, used for backward/inward flow):
		const double wM0[3] = {
			-f2m5 * (2.0 * z2 - 4.0 * z + 1.0) / (12.0 * f2m3 * pm3),
			 (10.0 * z2 - 29.0 * z + 8.0) / (12.0 * pm3),
			 f2m1 * (4.0 * z2 - 17.0 * z + 17.0) / (12.0 * f2m3 * pm3)
		};
		const double wM1[3] = {
			 f2m3 * (4.0 * z2 + z - 1.0) / (12.0 * f2m1 * pm1),
			 (10.0 * z2 - 11.0 * z - 10.0) / (12.0 * pm1),
			-f2p1 * (2.0 * z2 - 4.0 * z + 1.0) / (12.0 * f2m1 * pm1)
		};
		const double wM2[3] = {
			 f2m1 * (22.0 * z2 + 46.0 * z + 17.0) / (12.0 * f2p1 * pp1),
			-(14.0 * z2 + 17.0 * z - 8.0) / (12.0 * pp1),
			 f2p3 * (4.0 * z2 + z - 1.0) / (12.0 * f2p1 * pp1)
		};

		// Optimal weights (physical stencil numbering k = 0, 1, 2)
		const double dP[3] = {
			2.0 * f2m3 * pm3 * (3.0 * z4 - 10.0 * z2 + 4.0) / (5.0 * f2m1 * (4.0 * z2 - 9.0 * z + 4.0) * p5),
			3.0 * pm1 * (48.0 * z4 * z2 - 154.0 * z4 * z - 83.0 * z4 + 500.0 * z3 - 191.0 * z2 - 192.0 * z + 96.0) / (10.0 * (2.0 * z2 - 1.0) * (4.0 * z2 - 9.0 * z + 4.0) * p5),
			3.0 * f2p1 * pp1 * (6.0 * z4 - 25.0 * z3 + 20.0 * z2 + 15.0 * z - 12.0) / (10.0 * f2m1 * (2.0 * z2 - 1.0) * p5)
		};
		const double dM[3] = {
			3.0 * f2m3 * pm3 * (6.0 * z4 + z3 - 19.0 * z2 - 4.0 * z + 4.0) / (10.0 * f2m1 * (2.0 * z2 - 4.0 * z + 1.0) * p5),
			3.0 * pm1 * (48.0 * z4 * z2 - 134.0 * z4 * z - 133.0 * z4 + 412.0 * z3 - 9.0 * z2 - 112.0 * z + 24.0) / (10.0 * (2.0 * z2 - 4.0 * z + 1.0) * (4.0 * z2 + z - 1.0) * p5),
			2.0 * f2p1 * pp1 * (3.0 * z4 - 12.0 * z3 + 8.0 * z2 + 8.0 * z - 3.0) / (5.0 * f2m1 * (4.0 * z2 + z - 1.0) * p5)
		};

		// Smoothness indicator quadratic forms (physical, symmetric; Bk[m][n] multiplies the
		// physically ascending substencil values of S_k). Shared by both flow directions.
		const double q0 = 432.0 * f2m3 * f2m3 * pm3 * pm3;
		const double q1 = 432.0 * f2m1 * f2m1 * pm1 * pm1;
		const double q2 = 432.0 * f2p1 * f2p1 * pp1 * pp1;

		double B0[3][3];
		double B1[3][3];
		double B2[3][3];

		B0[0][0] = f2m5 * f2m5 * (576.0 * z4 - 2340.0 * z3 + 3183.0 * z2 - 1638.0 * z + 283.0) / q0;
		B0[0][1] = -f2m5 * f2m3 * (1368.0 * z4 - 6948.0 * z3 + 11136.0 * z2 - 6375.0 * z + 1193.0) / q0;
		B0[0][2] = f2m5 * f2m1 * (792.0 * z4 - 4824.0 * z3 + 10113.0 * z2 - 8427.0 * z + 2164.0) / q0;
		B0[1][1] = f2m3 * f2m3 * (3600.0 * z4 - 21888.0 * z3 + 42108.0 * z2 - 26868.0 * z + 5431.0) / q0;
		B0[1][2] = -f2m3 * f2m1 * (2232.0 * z4 - 15804.0 * z3 + 38532.0 * z2 - 36549.0 * z + 10328.0) / q0;
		B0[2][2] = f2m1 * f2m1 * (1440.0 * z4 - 11628.0 * z3 + 34251.0 * z2 - 43512.0 * z + 20164.0) / q0;

		B1[0][0] = f2m3 * f2m3 * (576.0 * z4 + 36.0 * z3 - 381.0 * z2 - 12.0 * z + 64.0) / q1;
		B1[0][1] = -f2m3 * f2m1 * (936.0 * z4 - 900.0 * z3 - 1104.0 * z2 + 297.0 * z + 266.0) / q1;
		B1[0][2] = f2m3 * f2p1 * (360.0 * z4 - 720.0 * z3 + 141.0 * z2 + 219.0 * z - 74.0) / q1;
		B1[1][1] = f2m1 * f2m1 * (1872.0 * z4 - 3744.0 * z3 - 1236.0 * z2 + 3108.0 * z + 1303.0) / q1;
		B1[1][2] = -f2m1 * f2p1 * (936.0 * z4 - 2844.0 * z3 + 1812.0 * z2 + 867.0 * z - 505.0) / q1;
		B1[2][2] = f2p1 * f2p1 * (576.0 * z4 - 2340.0 * z3 + 3183.0 * z2 - 1638.0 * z + 283.0) / q1;

		B2[0][0] = f2m1 * f2m1 * (1440.0 * z4 + 5868.0 * z3 + 8007.0 * z2 + 4134.0 * z + 715.0) / q2;
		B2[0][1] = -f2m1 * f2p1 * (2232.0 * z4 + 6876.0 * z3 + 4512.0 * z2 - 2031.0 * z - 1261.0) / q2;
		B2[0][2] = f2m1 * f2p3 * (792.0 * z4 + 1656.0 * z3 + 393.0 * z2 - 495.0 * z - 182.0) / q2;
		B2[1][1] = f2p1 * f2p1 * (3600.0 * z4 + 7488.0 * z3 - 1956.0 * z2 - 6084.0 * z + 2383.0) / q2;
		B2[1][2] = -f2p1 * f2p3 * (1368.0 * z4 + 1476.0 * z3 - 1500.0 * z2 - 525.0 * z + 374.0) / q2;
		B2[2][2] = f2p3 * f2p3 * (576.0 * z4 + 36.0 * z3 - 381.0 * z2 - 12.0 * z + 64.0) / q2;

		// Symmetrize
		B0[1][0] = B0[0][1]; B0[2][0] = B0[0][2]; B0[2][1] = B0[1][2];
		B1[1][0] = B1[0][1]; B1[2][0] = B1[0][2]; B1[2][1] = B1[1][2];
		B2[1][0] = B2[0][1]; B2[2][0] = B2[0][2]; B2[2][1] = B2[1][2];

		// Map to flow-oriented substencil indexing: substencil r uses u_j = w[j - r].
		if (forwardFlow)
		{
			// Flow offset j - r corresponds to physical cell i + (j - r); substencil r covers
			// cells {i - r, ..., i - r + 2} = S_{2-r} in physically ascending order.
			for (int j = 0; j < 3; ++j)
			{
				C[0][j] = wP2[j];
				C[1][j] = wP1[j];
				C[2][j] = wP0[j];
			}
			D[0] = dP[2]; D[1] = dP[1]; D[2] = dP[0];
			for (int m = 0; m < 3; ++m)
			{
				for (int n = 0; n < 3; ++n)
				{
					B[0][m][n] = B2[m][n];
					B[1][m][n] = B1[m][n];
					B[2][m][n] = B0[m][n];
				}
			}
		}
		else
		{
			// Flow offset j - r corresponds to physical cell i - (j - r); substencil r covers
			// cells {i + r - 2, ..., i + r} = S_r in physically descending order.
			for (int j = 0; j < 3; ++j)
			{
				C[0][j] = wM0[2 - j];
				C[1][j] = wM1[2 - j];
				C[2][j] = wM2[2 - j];
			}
			D[0] = dM[0]; D[1] = dM[1]; D[2] = dM[2];
			for (int m = 0; m < 3; ++m)
			{
				for (int n = 0; n < 3; ++n)
				{
					B[0][m][n] = B0[2 - m][2 - n];
					B[1][m][n] = B1[2 - m][2 - n];
					B[2][m][n] = B2[2 - m][2 - n];
				}
			}
		}
	}

	/**
	 * @brief Evaluates the radial WENO coefficients of order 2 (third order scheme) at a given \f$ \zeta \f$
	 * @details Same conventions as radialCoefficientsOrder3(), with substencils of two cells:
	 *          S0 = {i-1, i}, S1 = {i, i+1}. The smoothness indicators are perfect squares,
	 *          \f$ \beta_r = \gamma_r (u_1 - u_0)^2 \f$ with geometry factors \f$ \gamma_r(\zeta) \f$.
	 */
	static void radialCoefficientsOrder2(double zeta, bool forwardFlow, double C[3][3], double D[3], double G[2])
	{
		const double z = zeta;
		const double z2 = z * z;

		const double f2m3 = 2.0 * z - 3.0;
		const double f2m1 = 2.0 * z - 1.0;
		const double f2p1 = 2.0 * z + 1.0;
		const double s0den = 3.0 * z2 - 6.0 * z + 2.0;  // roots below zeta = 2
		const double s1den = 3.0 * z2 - 1.0;
		const double pm1 = z2 - z - 1.0;

		// Physical substencil value weights (ascending order within each stencil)
		// Outer face rho_{i+1/2} ("+"):
		const double vP0[2] = { -f2m3 * (3.0 * z - 2.0) / (4.0 * s0den), f2m1 * (9.0 * z - 14.0) / (4.0 * s0den) };
		const double vP1[2] = { f2m1 * (3.0 * z + 2.0) / (4.0 * s1den), f2p1 * (3.0 * z - 2.0) / (4.0 * s1den) };
		// Inner face rho_{i-1/2} ("-"):
		const double vM0[2] = { f2m3 * (3.0 * z - 1.0) / (4.0 * s0den), f2m1 * (3.0 * z - 5.0) / (4.0 * s0den) };
		const double vM1[2] = { f2m1 * (9.0 * z + 5.0) / (4.0 * s1den), -f2p1 * (3.0 * z - 1.0) / (4.0 * s1den) };

		// Optimal weights (physical numbering)
		const double eP[2] = {
			(2.0 * z2 - 1.0) * s0den / (3.0 * f2m1 * (3.0 * z - 2.0) * pm1),
			s1den * (4.0 * z2 - 9.0 * z + 4.0) / (3.0 * f2m1 * (3.0 * z - 2.0) * pm1)
		};
		const double eM[2] = {
			s0den * (4.0 * z2 + z - 1.0) / (3.0 * f2m1 * (3.0 * z - 1.0) * pm1),
			s1den * (2.0 * z2 - 4.0 * z + 1.0) / (3.0 * f2m1 * (3.0 * z - 1.0) * pm1)
		};

		// Smoothness indicator geometry factors: beta_{S0} = g0 * (cbar_i - cbar_{i-1})^2, etc.
		const double g0 = 9.0 * f2m3 * f2m3 * f2m1 * f2m1 / (16.0 * s0den * s0den);
		const double g1 = 9.0 * f2m1 * f2m1 * f2p1 * f2p1 / (16.0 * s1den * s1den);

		if (forwardFlow)
		{
			// Substencil r covers cells {i - r, i - r + 1} = S_{1-r} in ascending order
			C[0][0] = vP1[0]; C[0][1] = vP1[1];
			C[1][0] = vP0[0]; C[1][1] = vP0[1];
			D[0] = eP[1]; D[1] = eP[0];
			G[0] = g1; G[1] = g0;
		}
		else
		{
			// Substencil r covers cells {i + r - 1, i + r} = S_r in descending order
			C[0][0] = vM0[1]; C[0][1] = vM0[0];
			C[1][0] = vM1[1]; C[1][1] = vM1[0];
			D[0] = eM[0]; D[1] = eM[1];
			G[0] = g0; G[1] = g1;
		}
	}

	/**
	 * @brief Reconstructs a cell face value from geometry-weighted (radial) volume averages
	 * @details See the public reconstructRadial() methods. Boundary cells always use order
	 *          reduction because the radial optimal weights do not exist for stencils that
	 *          protrude the domain.
	 * @tparam wantJac Determines if the gradient is computed (@c true) or not (@c false)
	 */
	template <typename StateType, typename StencilType, bool wantJac>
	int reconstructRadial(unsigned int cellIdx, unsigned int numCells, double zeta, bool forwardFlow, const StencilType& w, StateType& result, double* const Dvm)
	{
#if defined(ACTIVE_SETFAD) || defined(ACTIVE_SFAD)
		using cadet::sqr;
		using sfad::sqr;
#endif

		// Reduce order near boundaries such that the stencil always fits inside the domain
		const int order = std::max(1, std::min(std::min(static_cast<int>(cellIdx) + 1, _order), std::min(static_cast<int>(numCells - cellIdx), _order)));

		// Simple upwind scheme
		if (order == 1)
		{
			result = w[0];
			if (wantJac)
				*Dvm = 1.0;
			return order;
		}

		const int sl = 2 * order - 1;

		// Evaluate the zeta-dependent coefficients
		double C[3][3];
		double D[3];
		double B[3][3][3];  // beta quadratic forms (order 3)
		double G[2];        // beta geometry factors (order 2)
		if (order == 3)
			radialCoefficientsOrder3(zeta, forwardFlow, C, D, B);
		else
			radialCoefficientsOrder2(zeta, forwardFlow, C, D, G);

		// Allocate memory for intermediate values: beta, alpha (= omega), and vr
		StateType* const work = _intermediateValues.create<StateType>(3 * order);
		StateType* const beta = work;
		StateType* const alpha = work + order;
		StateType* const omega = work + order;
		StateType* const vr = work + 2 * order; // Reconstructed values

		// Substencil r uses the stencil values u_j = w[j - r], j = 0, ..., order-1
		if (order == 3)
		{
			for (int r = 0; r < 3; ++r)
			{
				beta[r] = 0.0;
				for (int m = 0; m < 3; ++m)
					for (int n = 0; n < 3; ++n)
						beta[r] += B[r][m][n] * w[m - r] * w[n - r];
			}
		}
		else
		{
			beta[0] = G[0] * sqr(w[1] - w[0]);
			beta[1] = G[1] * sqr(w[0] - w[-1]);
		}

		// Add eps to avoid divide-by-zeros
		for (int r = 0; r < order; ++r)
			beta[r] += _epsilon;

		// Calculate weights
		for (int r = 0; r < order; ++r)
			alpha[r] = D[r] / sqr(beta[r]);

		// Normalize weights
		StateType alpha_sum = alpha[0];
		for (int r = 1; r < order; ++r)
			alpha_sum += alpha[r];
		for (int r = 0; r < order; ++r)
			omega[r] /= alpha_sum;

		// Calculate reconstructed values
		for (int r = 0; r < order; ++r)
		{
			vr[r] = 0.0;
			for (int j = 0; j < order; ++j)
				vr[r] += C[r][j] * w[j - r];
		}

		// Weighted sum
		result = 0.0;
		for (int r = 0; r < order; ++r)
			result += vr[r] * omega[r];

		// Jacobian
		if (wantJac)
		{
			// d(result)/d(w_o) = sum_r [ q_r * d(omega_r)/d(w_o) + omega_r * d(q_r)/d(w_o) ]
			// with d(omega_r)/d(w_o) = (d(alpha_r)/d(w_o) - omega_r * sum_k d(alpha_k)/d(w_o)) / alpha_sum,
			//      d(alpha_r)/d(w_o) = -2 D_r / (eps + beta_r)^3 * d(beta_r)/d(w_o).
			const double aSum = static_cast<double>(alpha_sum);

			// dBeta[r][i]: derivative of beta_r w.r.t. stencil position i (i = 0 corresponds to w[-order+1])
			double dBeta[3][5] = { { 0.0 } };
			if (order == 3)
			{
				for (int r = 0; r < 3; ++r)
					for (int m = 0; m < 3; ++m)
					{
						double dot = 0.0;
						for (int n = 0; n < 3; ++n)
							dot += B[r][m][n] * static_cast<double>(w[n - r]);
						dBeta[r][m - r + order - 1] = 2.0 * dot;
					}
			}
			else
			{
				dBeta[0][2] = 2.0 * G[0] * static_cast<double>(w[1] - w[0]);
				dBeta[0][1] = -dBeta[0][2];
				dBeta[1][1] = 2.0 * G[1] * static_cast<double>(w[0] - w[-1]);
				dBeta[1][0] = -dBeta[1][1];
			}

			// dAlpha[r][i]
			double dAlpha[3][5] = { { 0.0 } };
			for (int r = 0; r < order; ++r)
			{
				const double betaR = static_cast<double>(beta[r]);
				const double fac = -2.0 * D[r] / (betaR * betaR * betaR);
				for (int i = 0; i < sl; ++i)
					dAlpha[r][i] = fac * dBeta[r][i];
			}

			for (int i = 0; i < sl; ++i)
			{
				double dASum = 0.0;
				for (int r = 0; r < order; ++r)
					dASum += dAlpha[r][i];

				double val = 0.0;
				for (int r = 0; r < order; ++r)
				{
					// d(omega_r)/d(w_i)
					const double dOmega = (dAlpha[r][i] - static_cast<double>(omega[r]) * dASum) / aSum;
					val += static_cast<double>(vr[r]) * dOmega;

					// d(q_r)/d(w_i): stencil position i corresponds to offset o = i - order + 1 and
					// contributes to substencil r via j = o + r
					const int j = i - order + 1 + r;
					if ((j >= 0) && (j < order))
						val += static_cast<double>(omega[r]) * C[r][j];
				}
				Dvm[i] = val;
			}
		}

		_intermediateValues.destroy<StateType>();
		return order;
	}

	/**
	 * @brief Per-cell precomputed geometry-exact WENO coefficients (one flow direction)
	 * @details Flow-oriented indexing as in the equidistant schemes: substencil \f$ r \f$ uses the
	 *          stencil values \f$ u_j = w[j - r] \f$, \f$ j = 0, \ldots, \text{order}-1 \f$, with
	 *          \f$ q_r = \sum_j C[r][j] u_j \f$ and \f$ \beta_r = \sum_{m,n} B[r][m][n] u_m u_n \f$.
	 */
	struct GeomExactCellCoefficients
	{
		int order; //!< WENO order used at this cell (after boundary and validity reduction)
		double C[3][3]; //!< Substencil value weights
		double D[3]; //!< Optimal weights
		double B[3][3][3]; //!< Smoothness indicator quadratic forms
	};

	/**
	 * @brief Computes the geometry-exact coefficients of one cell for one flow direction (defined in Weno.cpp)
	 * @return @c true if the coefficients are valid (consistent and positive optimal weights), @c false otherwise
	 */
	static bool computeGeomExactCellCoefficients(double a, double b, int weightPower, const std::vector<double>& faces,
		int cell, bool forwardFlow, int order, GeomExactCellCoefficients& out);

	std::vector<GeomExactCellCoefficients> _geomForward; //!< Geometry-exact coefficients for forward flow
	std::vector<GeomExactCellCoefficients> _geomBackward; //!< Geometry-exact coefficients for backward flow

	/**
	 * @brief Reconstructs a cell face value from geometry-weighted volume averages using precomputed coefficients
	 * @details See the public reconstructGeomExact() methods.
	 * @tparam wantJac Determines if the gradient is computed (@c true) or not (@c false)
	 */
	template <typename StateType, typename StencilType, bool wantJac>
	int reconstructGeomExact(unsigned int cellIdx, unsigned int numCells, bool forwardFlow, const StencilType& w, StateType& result, double* const Dvm)
	{
#if defined(ACTIVE_SETFAD) || defined(ACTIVE_SFAD)
		using cadet::sqr;
		using sfad::sqr;
#endif

		const GeomExactCellCoefficients& cc = forwardFlow ? _geomForward[cellIdx] : _geomBackward[cellIdx];
		const int order = cc.order;

		// Simple upwind scheme
		if (order == 1)
		{
			result = w[0];
			if (wantJac)
				*Dvm = 1.0;
			return order;
		}

		const int sl = 2 * order - 1;

		// Allocate memory for intermediate values: beta, alpha (= omega), and vr
		StateType* const work = _intermediateValues.create<StateType>(3 * order);
		StateType* const beta = work;
		StateType* const alpha = work + order;
		StateType* const omega = work + order;
		StateType* const vr = work + 2 * order; // Reconstructed values

		// Smoothness indicators: beta_r = u_r^T B_r u_r with u_j = w[j - r]
		for (int r = 0; r < order; ++r)
		{
			beta[r] = 0.0;
			for (int m = 0; m < order; ++m)
				for (int n = 0; n < order; ++n)
					beta[r] += cc.B[r][m][n] * w[m - r] * w[n - r];

			// Add eps to avoid divide-by-zeros
			beta[r] += _epsilon;
		}

		// Calculate weights
		for (int r = 0; r < order; ++r)
			alpha[r] = cc.D[r] / sqr(beta[r]);

		// Normalize weights
		StateType alpha_sum = alpha[0];
		for (int r = 1; r < order; ++r)
			alpha_sum += alpha[r];
		for (int r = 0; r < order; ++r)
			omega[r] /= alpha_sum;

		// Calculate reconstructed values
		for (int r = 0; r < order; ++r)
		{
			vr[r] = 0.0;
			for (int j = 0; j < order; ++j)
				vr[r] += cc.C[r][j] * w[j - r];
		}

		// Weighted sum
		result = 0.0;
		for (int r = 0; r < order; ++r)
			result += vr[r] * omega[r];

		// Jacobian
		if (wantJac)
		{
			// d(result)/d(w_o) = sum_r [ q_r * d(omega_r)/d(w_o) + omega_r * d(q_r)/d(w_o) ]
			const double aSum = static_cast<double>(alpha_sum);

			// dBeta[r][i]: derivative of beta_r w.r.t. stencil position i (i = 0 corresponds to w[-order+1])
			double dBeta[3][5] = { { 0.0 } };
			for (int r = 0; r < order; ++r)
				for (int m = 0; m < order; ++m)
				{
					double dot = 0.0;
					for (int n = 0; n < order; ++n)
						dot += cc.B[r][m][n] * static_cast<double>(w[n - r]);
					dBeta[r][m - r + order - 1] = 2.0 * dot;
				}

			double dAlpha[3][5] = { { 0.0 } };
			for (int r = 0; r < order; ++r)
			{
				const double betaR = static_cast<double>(beta[r]);
				const double fac = -2.0 * cc.D[r] / (betaR * betaR * betaR);
				for (int i = 0; i < sl; ++i)
					dAlpha[r][i] = fac * dBeta[r][i];
			}

			for (int i = 0; i < sl; ++i)
			{
				double dASum = 0.0;
				for (int r = 0; r < order; ++r)
					dASum += dAlpha[r][i];

				double val = 0.0;
				for (int r = 0; r < order; ++r)
				{
					const double dOmega = (dAlpha[r][i] - static_cast<double>(omega[r]) * dASum) / aSum;
					val += static_cast<double>(vr[r]) * dOmega;

					const int j = i - order + 1 + r;
					if ((j >= 0) && (j < order))
						val += static_cast<double>(omega[r]) * cc.C[r][j];
				}
				Dvm[i] = val;
			}
		}

		_intermediateValues.destroy<StateType>();
		return order;
	}

	double _epsilon; //!< Small number preventing divsion by zero
	int _order; //!< Selected WENO order
	BoundaryTreatment _boundaryTreatment; //!< Controls how to treat boundary cells
	ArrayPool _intermediateValues; //!< Buffer for intermediate and temporary values

	static const double _wenoD2[2];
	static const double _wenoC2[2*2];
	static const double _wenoJbvv2[2*3*3];

	static const double _wenoD3[3];
	static const double _wenoC3[3*3];
	static const double _wenoJbvv3[3*5*5];

};

} // namespace cadet

#endif  // LIBCADET_WENO_HPP_
