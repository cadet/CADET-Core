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
 *
 * ## Non-equidistant grid support
 *
 * When setGrid() is called the scheme precomputes per-cell FV-consistent sub-stencil
 * reconstruction coefficients and optimal weights using the **primitive-function / Lagrange
 * derivative** approach:
 *
 *   1. Define the cumulative sum \f$S_{j+1/2} = \sum_{m} h_m \bar{u}_m\f$ at each face.
 *   2. For each sub-stencil \f$r\f$ fit a degree-\f$(p)\f$ polynomial through the \f$p+1\f$
 *      consecutive \f$S\f$ values that span the sub-stencil cells, then differentiate at the
 *      target face \f$x_{i+1/2}\f$ (or \f$x_{i-1/2}\f$ for backward flow).
 *   3. Because \f$S\f$ is linear in the cell averages \f$\bar{u}\f$, the result is
 *      \f$\hat{u}^{(r)} = \sum_j c_r^{(i)}[j] \, \bar{u}_{\text{stencil}}\f$ with
 *      per-face coefficients \f$c_r^{(i)}[j]\f$.
 *   4. Optimal weights \f$d_r^{(i)}\f$ are obtained from the outer-cell equations of the
 *      full \f$(2r-1)\f$-cell reconstruction (a small closed-form linear system).
 *
 * The smoothness indicators \f$\beta_r\f$ retain the standard equidistant-grid formulae;
 * this introduces only an \f$O(h)\f$ error in the ENO stencil-selection step but does not
 * reduce the overall convergence order of the scheme.
 *
 * Both a forward-flow set (reconstructing the right face of each cell) and a backward-flow
 * set (reconstructing the left face) are precomputed.  Call setForwardFlow() before each
 * residual evaluation so that reconstruct() uses the correct set.
 *
 * The equidistant static coefficients (_wenoC, _wenoD) serve as the fall-back whenever
 * setGrid() has not been called or near-boundary order-reduction is active.
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
	Weno() : _epsilon(1e-10), _order(maxOrder()), _boundaryTreatment(BoundaryTreatment::ReduceOrder),
	         _intermediateValues(3 * maxOrder() * sizeof(active)), _forwardFlow(true) { }

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
		cadet_assert(order <= static_cast<int>(maxOrder()));
		cadet_assert(order > 0);
		if (order != _order)
		{
			// Invalidate precomputed data since it depends on order.
			_cellReconC_fwd.clear();
			_cellReconD_fwd.clear();
			_cellReconC_bwd.clear();
			_cellReconD_bwd.clear();
		}
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

	/**
	 * @brief Precomputes per-cell FV-consistent reconstruction data for a non-equidistant grid
	 * @details For each interior cell and each sub-stencil \f$r\f$, the coefficients
	 *          \f$c_r[j]\f$ and optimal weights \f$d_r\f$ are derived from Lagrange
	 *          differentiation of the primitive-function polynomial (see class documentation).
	 *          Both forward-flow (right-face) and backward-flow (left-face) sets are stored.
	 *          Near-boundary cells that require order reduction fall back to the equidistant
	 *          static tables in reconstruct().
	 * @param [in] cellCenters  Cell-centre coordinates, length @p nCells (unused, kept for interface compatibility)
	 * @param [in] cellFaces    Cell-face coordinates, length @p nCells + 1
	 * @param [in] nCells       Number of cells
	 */
	void setGrid(const double* /*cellCenters*/, const double* cellFaces, unsigned int nCells)
	{
		_facesStore.assign(cellFaces, cellFaces + nCells + 1);
		_nCellsStored = nCells;

		const int ord = _order;
		const int ord2 = ord * ord;

		_cellReconC_fwd.resize(nCells * ord2);
		_cellReconD_fwd.resize(nCells * ord);
		_cellReconC_bwd.resize(nCells * ord2);
		_cellReconD_bwd.resize(nCells * ord);

		// Temporary storage for one sub-stencil row and one full-stencil row.
		double tmp[2 * maxOrder() - 1]; // max length = 2*maxOrder-1

		for (unsigned int i = 0; i < nCells; ++i)
		{
			// Interior condition: all sub-stencils and the full (2*ord-1)-cell stencil
			// must fit within [0, nCells).  For both forward and backward flow this
			// requires  ord-1 <= i <= nCells-ord.
			const bool interior = ((int)i >= ord - 1) && ((int)i + ord <= (int)nCells);

			// ---- Forward flow: reconstruct right face cellFaces[i+1] ----
			//
			// Sub-stencil r uses cells i-r ... i-r+ord-1 (physical left to right).
			// Face array: cellFaces[i-r .. i-r+ord].  Target local index t = r+1.
			// Stencil access w[-r+j] maps j → physical cell i-r+j (left-to-right order).
			// lagrangeDerivCoeffs output tmp[j] is the coefficient of physical cell i-r+j
			// which equals the coefficient of w[-r+j], i.e. c[r + ord*j].
			if (interior)
			{
				double cFull[2 * maxOrder() - 1];

				for (int r = 0; r < ord; ++r)
				{
					const double* f = &_facesStore[i - r]; // f[0..ord]
					lagrangeDerivCoeffs(f, (unsigned int)ord, (unsigned int)(r + 1), tmp);
					// Store as c[r + ord*j] layout (stride ord over j).
					for (int j = 0; j < ord; ++j)
						_cellReconC_fwd[i * ord2 + r + ord * j] = tmp[j];
				}

				// Full (2*ord-1)-cell stencil; target = cellFaces[i+1], local index = ord.
				lagrangeDerivCoeffs(&_facesStore[i - (ord - 1)],
				                    (unsigned int)(2 * ord - 1), (unsigned int)ord, cFull);
				computeOptimalWeights(ord, _cellReconC_fwd.data() + i * ord2,
				                      cFull, _cellReconD_fwd.data() + i * ord);
			}
			else
			{
				// Near boundary: copy static equidistant tables as fallback.
				const double* ds = (ord == 2) ? _wenoD2 : _wenoD3;
				const double* cs = (ord == 2) ? _wenoC2 : _wenoC3;
				for (int r = 0; r < ord; ++r)
				{
					_cellReconD_fwd[i * ord + r] = ds[r];
					for (int j = 0; j < ord; ++j)
						_cellReconC_fwd[i * ord2 + r + ord * j] = cs[r + ord * j];
				}
			}

			// ---- Backward flow: reconstruct left face cellFaces[i] ----
			//
			// In backward flow the stencil is mirrored: w[k] = ū_{col-k}.
			// Sub-stencil r accesses w[-r+j] = ū_{col+r-j} for j=0..ord-1
			//   → physical cells col+r, col+r-1, ..., col+r-ord+1 (right to left in stencil).
			// The face array for lagrangeDerivCoeffs is sorted left-to-right:
			//   cellFaces[i+r-ord+1 .. i+r+1], target = cellFaces[i] (left face),
			//   local index t = i - (i+r-ord+1) = ord-1-r.
			// lagrangeDerivCoeffs output tmp[l] gives the coefficient of physical cell
			//   i+r-ord+1+l.  In stencil ordering j=ord-1-l, so:
			//   c[r + ord*j] = c[r + ord*(ord-1-l)] = tmp[l],
			//   equivalently c[r + ord*k] = tmp[ord-1-k].
			//
			// The full-stencil backward uses the same physical face array as forward but
			// targets the left face (local index = ord-1 in the [i-(ord-1)..i+ord] array).
			// The lagrange output tmp_full[l] is for physical cell i-(ord-1)+l; in stencil
			// ordering j = 2*(ord-1) - l, so cFull_bwd[j] = tmp_full[2*(ord-1)-j].
			if (interior)
			{
				for (int r = 0; r < ord; ++r)
				{
					const double* f = &_facesStore[i + r - ord + 1]; // f[0..ord]
					lagrangeDerivCoeffs(f, (unsigned int)ord, (unsigned int)(ord - 1 - r), tmp);
					// Reverse: stencil index k corresponds to lagrange index ord-1-k.
					for (int k = 0; k < ord; ++k)
						_cellReconC_bwd[i * ord2 + r + ord * k] = tmp[ord - 1 - k];
				}

				// Full stencil: faces cellFaces[i-(ord-1) .. i+ord], target local index = ord-1.
				double cFullLagrange[2 * maxOrder() - 1];
				lagrangeDerivCoeffs(&_facesStore[i - (ord - 1)],
				                    (unsigned int)(2 * ord - 1), (unsigned int)(ord - 1), cFullLagrange);
				// Reverse for backward-stencil ordering: cFull_bwd[j] = cFullLagrange[2*(ord-1)-j].
				double cFull_bwd[2 * maxOrder() - 1];
				for (int j = 0; j < 2 * ord - 1; ++j)
					cFull_bwd[j] = cFullLagrange[2 * (ord - 1) - j];

				computeOptimalWeights(ord, _cellReconC_bwd.data() + i * ord2,
				                      cFull_bwd, _cellReconD_bwd.data() + i * ord);
			}
			else
			{
				const double* ds = (ord == 2) ? _wenoD2 : _wenoD3;
				const double* cs = (ord == 2) ? _wenoC2 : _wenoC3;
				for (int r = 0; r < ord; ++r)
				{
					_cellReconD_bwd[i * ord + r] = ds[r];
					for (int j = 0; j < ord; ++j)
						_cellReconC_bwd[i * ord2 + r + ord * j] = cs[r + ord * j];
				}
			}
		}
	}

	/**
	 * @brief Selects which face subsequent reconstruct() calls target
	 * @param [in] forward  @c true  = right face of each cell (positive / outward flow),
	 *                      @c false = left  face of each cell (negative / inward  flow)
	 */
	inline void setForwardFlow(bool forward) CADET_NOEXCEPT { _forwardFlow = forward; }

private:

	/**
	 * @brief Computes Lagrange-derivative sub-stencil reconstruction coefficients
	 * @details Fits a degree-@p n polynomial through the primitive function S at the
	 *          @p n+1 face positions @p f[0..n], then differentiates at the face with
	 *          local index @p t.  Because S is linear in the cell averages the result
	 *          is a linear combination of those averages:
	 *          \f[ c[j] = h_j \sum_{l=j+1}^{n} \lambda_l, \quad
	 *              \lambda_l = L_l'(f_t) \f]
	 *          where \f$L_l\f$ is the \f$l\f$-th Lagrange basis and \f$h_j = f_{j+1}-f_j\f$.
	 *          The output array @p c is laid out so that it can be directly used in place of
	 *          a row of the static _wenoC tables: @p c[r + ord*j] convention is achieved by
	 *          the caller placing the output at the correct offset.
	 *
	 * @param [in]  f  Face positions, length @p n+1
	 * @param [in]  n  Number of cells in the sub-stencil
	 * @param [in]  t  Local index of the target face within @p f (0 ≤ t ≤ n)
	 * @param [out] c  Reconstruction coefficients for cells j=0..n-1, stored as c[j*ord + r]
	 *                 but written here as c[r], c[r + ord], c[r + 2*ord], ...
	 *                 — actually stored contiguously c[0..n-1] and the caller uses the
	 *                 correct stride.  Here we write c[r + ord*j] via pointer arithmetic
	 *                 so the output stride must be 1 (see call sites).
	 *
	 * @note The output is written as c[j] for j=0..n-1 (j indexes the cell within the
	 *       sub-stencil).  The calling code reads it back as c[r + ord*j], so the caller
	 *       must write this sub-stencil's row into the correct sub-array.
	 */
	static void lagrangeDerivCoeffs(const double* f, unsigned int n, unsigned int t, double* c)
	{
		// Step 1: compute λ[l] = L_l'(f[t]) for l = 0 .. n
		//
		// At a node f[t] the Lagrange derivative simplifies to:
		//   l == t : λ[t] = Σ_{m≠t} 1/(f[t] − f[m])
		//   l != t : λ[l] = Π_{m≠l,m≠t}(f[t]−f[m]) / Π_{m≠l}(f[l]−f[m])
		double lambda[2 * maxOrder()]; // n ≤ 2*maxOrder-1, so n+1 ≤ 2*maxOrder
		for (unsigned int l = 0; l <= n; ++l)
		{
			if (l == t)
			{
				double s = 0.0;
				for (unsigned int m = 0; m <= n; ++m)
					if (m != t) s += 1.0 / (f[t] - f[m]);
				lambda[l] = s;
			}
			else
			{
				double num = 1.0, den = 1.0;
				for (unsigned int m = 0; m <= n; ++m)
				{
					if (m != l && m != t) num *= f[t] - f[m];
					if (m != l)           den *= f[l] - f[m];
				}
				lambda[l] = num / den;
			}
		}

		// Step 2: c[j] = h[j] * Σ_{l=j+1}^n λ[l]   (suffix sum, right-to-left)
		// The output array c stores c[0], c[1], ..., c[n-1] contiguously.
		// Caller's indexing: _cellReconC[cellIdx*ord2 + r + ord*j] = c[j] for sub-stencil r.
		// However the caller passes a pointer to the base of sub-stencil r's column in the
		// _cellReconC array (stride ord), but here we just write indices 0..n-1.
		// The actual stride is set by the caller; see the call sites in setGrid().
		//
		// NOTE: The static tables use c[r + ord*j] indexing (r = row, j = column).
		// Here we compute one row at a time, writing c[j] for j=0..n-1.
		// The caller stores this row at the correct memory offset.
		double suffix = 0.0;
		for (int j = (int)n - 1; j >= 0; --j)
		{
			suffix += lambda[(unsigned int)(j + 1)];
			c[j] = (f[j + 1] - f[j]) * suffix;
		}
	}

	/**
	 * @brief Computes per-face optimal WENO weights from the outer-cell equations
	 * @details For a sub-stencil layout where sub-stencil r=ord-1 is the only one
	 *          contributing to the leftmost full-stencil cell, and r=0 is the only one
	 *          contributing to the rightmost:
	 *          \f[ d_{ord-1} = c_{\mathrm{full}}[0]   / c_{ord-1}[0], \quad
	 *              d_0       = c_{\mathrm{full}}[2(ord-1)] / c_0[ord-1], \quad
	 *              d_1 = 1 - \sum_{r \ne 1} d_r \f]
	 *
	 * @param [in]  ord      WENO order
	 * @param [in]  cellC    Per-cell sub-stencil coefficients for this cell, size ord*ord,
	 *                       indexed as cellC[r + ord*j] (same layout as _wenoC static tables)
	 * @param [in]  cFull    Full-(2*ord-1)-cell reconstruction coefficients, size 2*ord-1
	 * @param [out] d        Optimal weights d[0..ord-1]
	 */
	static void computeOptimalWeights(int ord, const double* cellC, const double* cFull, double* d)
	{
		// Outer equations are each solved by one sub-stencil exclusively.
		// r=ord-1: leftmost sub-stencil contributes to full-stencil cell j=0
		//   → d[ord-1] * cellC[(ord-1) + ord*0] = cFull[0]
		// r=0:     rightmost sub-stencil contributes to full-stencil cell j=2*(ord-1)
		//   → d[0] * cellC[0 + ord*(ord-1)] = cFull[2*(ord-1)]
		d[ord - 1] = cFull[0]            / cellC[(ord - 1)];            // cellC[(ord-1) + ord*0]
		d[0]       = cFull[2 * (ord - 1)] / cellC[0 + ord * (ord - 1)]; // cellC[0 + ord*(ord-1)]

		// All middle weights: constrained by Σ d_r = 1
		// For ord == 2: there are no middle weights, check is automatic.
		// For ord == 3: d[1] = 1 - d[0] - d[2].
		double sum = d[0] + d[ord - 1];
		for (int r = 1; r < ord - 1; ++r)
			d[r] = 0.0; // initialise; only relevant when ord > 3 (not currently implemented)
		if (ord >= 3)
			d[1] = 1.0 - sum;
		// For ord == 2: d[0] + d[1] should already equal 1 by construction; no middle weights.
	}

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

		// Calculate smoothness measures (equidistant formula retained — see class docs)
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

		// Override d and c with per-cell non-equidistant values when available
		// and when operating at full order (boundary cells fall back to static tables).
		const bool useNonEqui = !_cellReconC_fwd.empty() && (order == _order);
		if (useNonEqui)
		{
			const int ord2 = order * order;
			if (_forwardFlow)
			{
				d = _cellReconD_fwd.data() + cellIdx * order;
				c = _cellReconC_fwd.data() + cellIdx * ord2;
			}
			else
			{
				d = _cellReconD_bwd.data() + cellIdx * order;
				c = _cellReconC_bwd.data() + cellIdx * ord2;
			}
			// Jbvv stays as the static-table value: β formula is unchanged.
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
				dot += static_cast<double>(vr[r]) * static_cast<double>(omega[r]);
			for (int r = 0; r < order; ++r)
				vr[r] = (vr[r] - dot) / alpha_sum;

			// Multiply with "d(alpha)/d(beta)" to get "d(result)/d(beta)"
			// Uses per-cell d[r] (already pointed to by d above).
			for (int r = 0; r < order; ++r)
				vr[r] *= -2.0 * d[r] / pow(beta[r], 3.0);

			// Multiply with "d(beta)/d(v)" to get Dvm = "d(result)/d(v)"
			// Uses static Jbvv (β formula is equidistant and unchanged).
			for (int j = 0; j < sl; ++j)
			{
				Dvm[j] = 0.0;
				for (int r = 0; r < order; ++r)
				{
					dot = 0.0;
					for (int i = 0; i < sl; ++i)
						dot += static_cast<double>(Jbvv[r + order * j + order * sl * i]) * static_cast<double>(w[i - order + 1]);
					Dvm[j] += static_cast<double>(vr[r]) * dot;
				}
			}

			// 2. Constant omega[r] in (*)
			// Uses per-cell c (already pointed to by c above).
			for (int r = 0; r < order; ++r)
				for (int j = 0; j < order; ++j)
					Dvm[order - 1 + j - r] += static_cast<double>(omega[r]) * c[r + order * j];
		}

		_intermediateValues.destroy<StateType>();
		return order;
	}

	double _epsilon; //!< Small number preventing division by zero
	int _order;      //!< Selected WENO order
	BoundaryTreatment _boundaryTreatment; //!< Controls how to treat boundary cells
	ArrayPool _intermediateValues; //!< Buffer for intermediate and temporary values
	bool _forwardFlow; //!< true = reconstruct right face; false = reconstruct left face

	// Non-equidistant per-cell data (empty when equidistant static tables are used)
	std::vector<double> _cellReconC_fwd; //!< Sub-stencil coefficients, forward flow,  shape [nCells][order][order], c[r + order*j]
	std::vector<double> _cellReconD_fwd; //!< Optimal weights,             forward flow,  shape [nCells][order]
	std::vector<double> _cellReconC_bwd; //!< Sub-stencil coefficients, backward flow, shape [nCells][order][order], c[r + order*j]
	std::vector<double> _cellReconD_bwd; //!< Optimal weights,             backward flow, shape [nCells][order]
	std::vector<double> _facesStore;     //!< Owned copy of face positions, length nCells+1
	unsigned int _nCellsStored = 0;      //!< Number of cells for which data was precomputed

	static const double _wenoD2[2];
	static const double _wenoC2[2*2];
	static const double _wenoJbvv2[2*3*3];

	static const double _wenoD3[3];
	static const double _wenoC3[3*3];
	static const double _wenoJbvv3[3*5*5];

};

} // namespace cadet

#endif  // LIBCADET_WENO_HPP_
