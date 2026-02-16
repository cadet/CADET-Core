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
 * Implements the kernel of the axial convection dispersion transport operator.
 */

#ifndef LIBCADET_AXIALCONVECTIONDISPERSIONKERNEL_HPP_
#define LIBCADET_AXIALCONVECTIONDISPERSIONKERNEL_HPP_

#include "AutoDiff.hpp"
#include "Memory.hpp"
#include "Stencil.hpp"
#include "linalg/CompressedSparseMatrix.hpp"
#include "SimulationTypes.hpp"
#include "model/ParameterDependence.hpp"
#include "model/UnitOperation.hpp"

#include <algorithm>
#include <cmath>

namespace cadet
{

class Weno;
class HighResolutionKoren;

namespace model
{

namespace parts
{

namespace convdisp
{

template <typename T, typename Reconstruction_t>
struct AxialFlowParameters
{
	T u;
	active const* d_ax;
	T h;
	double* reconstructionDerivatives; //!< Holds derivatives of the reconstruction scheme
	Reconstruction_t* reconstruction; //!< The reconstruction scheme implementation
	ArrayPool* stencilMemory; //!< Provides memory for the stencil
	int strideCell;
	unsigned int nComp;
	unsigned int nCol;
	unsigned int offsetToInlet; //!< Offset to the first component of the inlet DOFs in the local state vector
	unsigned int offsetToBulk; //!< Offset to the first component of the first bulk cell in the local state vector
	IParameterParameterDependence* parDep;
	const IModel& model;
	bool gridEquidistant;
	std::vector<active> const* cellFaces; //!< Positions of the cell faces for non-equidistant grids
};


namespace impl
{
	template <class FaceContainerType>
	struct ReverseFaceAccessor
	{
		const FaceContainerType& faces;
		inline auto operator[](unsigned int idx) const -> decltype(faces[0]) { return faces[faces.size() - 1 - idx]; }
	};

	template <typename StateType, typename ResidualType, typename ParamType, typename ReconstrType, typename RowIteratorType, bool wantJac, bool wantRes = true>
	int residualForwardsAxialFlow(const SimulationTime& simTime, StateType const* y, double const* yDot, ResidualType* res, RowIteratorType jacBegin, const AxialFlowParameters<ParamType, ReconstrType>& p)
	{
		// True if the grid cell faces are provided and gridEquidistant is set to be false. 
		const bool nonEqGrid = p.cellFaces && !p.gridEquidistant;

		// h is the size of the cell
		const ParamType h2 = p.h * p.h;

		// The stencil caches parts of the state vector for better spatial coherence
		typedef CachingStencil<StateType, ArrayPool> StencilType;
		StencilType stencil(std::max(p.reconstruction->stencilSize(), 3u), *p.stencilMemory, std::max(p.reconstruction->order() - 1, 1));

		// The RowIterator is always centered on the main diagonal.
		// This means that jac[0] is the main diagonal, jac[-1] is the first lower diagonal,
		// and jac[1] is the first upper diagonal. We can also access the rows from left to
		// right beginning with the last lower diagonal moving towards the main diagonal and
		// continuing to the last upper diagonal by using the native() method.
		RowIteratorType jac;

		ResidualType* const resBulk = wantRes ? res + p.offsetToBulk : nullptr;
		StateType const* const yBulk = y + p.offsetToBulk;
		const ParamType colLength = nonEqGrid ? static_cast<ParamType>(p.cellFaces->back()) : static_cast<ParamType>(0.0);                // length of the column, needed for parameter dependence in non-equidistant grids

		for (unsigned int comp = 0; comp < p.nComp; ++comp)
		{
			if (wantJac)
				jac = jacBegin + comp;

			ResidualType* const resBulkComp = wantRes ? resBulk + comp : nullptr;
			StateType const* const yBulkComp = yBulk + comp;

			// Add time derivative to each cell
			if (yDot && wantRes)
			{
				double const* const yDotBulkComp = yDot + p.offsetToBulk + comp;
				for (unsigned int col = 0; col < p.nCol; ++col)
					resBulkComp[col * p.strideCell] = yDotBulkComp[col * p.strideCell];
			}
			else if (wantRes)
			{
				for (unsigned int col = 0; col < p.nCol; ++col)
					resBulkComp[col * p.strideCell] = 0.0;
			}

			// Fill stencil (left side with zeros, right side with states)
			for (int i = -std::max(p.reconstruction->order(), 2) + 1; i < 0; ++i)
				stencil[i] = 0.0;
			for (int i = 0; i < std::max(p.reconstruction->order(), 2); ++i)
				stencil[i] = yBulkComp[i * p.strideCell];

			// Reset reconstruction output
			StateType vm(0.0); // reconstructed value
			if (wantJac)
				std::fill(p.reconstructionDerivatives, p.reconstructionDerivatives + p.reconstruction->stencilSize(), 0.0);

			int wenoOrder = 0;
			const ParamType d_ax = static_cast<ParamType>(p.d_ax[comp]);

			// Iterate over all cells
			for (unsigned int col = 0; col < p.nCol; ++col)
			{
				const ParamType hCol = nonEqGrid ? static_cast<ParamType>((*p.cellFaces)[col + 1] - (*p.cellFaces)[col]) : p.h;              // size of the current cell
				const ParamType invHCol = static_cast<ParamType>(p.u) / hCol;

				// ------------------- Dispersion -------------------

				// Right side, leave out if we're in the last cell (boundary condition)
				if (cadet_likely(col < p.nCol - 1))
				{
					const ParamType hRight = nonEqGrid ? static_cast<ParamType>((*p.cellFaces)[col + 2] - (*p.cellFaces)[col + 1]) : p.h;     // size of the right cell
					const ParamType deltaZ = nonEqGrid ? static_cast<ParamType>(0.5) * (hCol + hRight) : p.h;                                 // size of the distance between the centers of the two cells, needed for non-equidistant grids
					const double relCoord = nonEqGrid ?
						static_cast<double>((*p.cellFaces)[col + 1]) / static_cast<double>(colLength) :
						static_cast<double>(col + 1) / p.nCol;                                                                                // relative coordinate of the cell face for parameter dependence

					const ParamType d_ax_right = d_ax * p.parDep->getValue(p.model, ColumnPosition{ relCoord, 0.0, 0.0 }, comp, ParTypeIndep, BoundStateIndep, static_cast<ParamType>(p.u));
					const ParamType coeff = nonEqGrid ? d_ax_right / (hCol * deltaZ) : d_ax_right / h2;

					if (wantRes)
						resBulkComp[col * p.strideCell] -= coeff * (stencil[1] - stencil[0]);
					// Jacobian entries
					if (wantJac)
					{
						const double coeffJac = static_cast<double>(coeff);
						jac[0] += coeffJac;
						jac[p.strideCell] -= coeffJac;
					}
				}

				// Left side, leave out if we're in the first cell (boundary condition)
				if (cadet_likely(col > 0))
				{
					const ParamType hLeft = nonEqGrid ? static_cast<ParamType>((*p.cellFaces)[col] - (*p.cellFaces)[col - 1]) : p.h;
					const ParamType deltaZ = nonEqGrid ? static_cast<ParamType>(0.5) * (hLeft + hCol) : p.h;
					const double relCoord = nonEqGrid ?
						static_cast<double>((*p.cellFaces)[col]) / static_cast<double>(colLength) :
						static_cast<double>(col) / p.nCol;

					const ParamType d_ax_left = d_ax * p.parDep->getValue(p.model, ColumnPosition{ relCoord, 0.0, 0.0 }, comp, ParTypeIndep, BoundStateIndep, static_cast<ParamType>(p.u));
					const ParamType coeff = nonEqGrid ? d_ax_left / (hCol * deltaZ) : d_ax_left / h2;

					if (wantRes)
						resBulkComp[col * p.strideCell] -= coeff * (stencil[-1] - stencil[0]);
					// Jacobian entries
					if (wantJac)
					{
						const double coeffJac = static_cast<double>(coeff);
						jac[0] += coeffJac;
						jac[-p.strideCell] -= coeffJac;
					}
				}

				// ------------------- Convection -------------------

				// Add convection through this cell's left face
				if (cadet_likely(col > 0))
				{
					// Remember that vm still contains the reconstructed value of the previous 
					// cell's *right* face, which is identical to this cell's *left* face!
					if (wantRes)
						resBulkComp[col * p.strideCell] -= invHCol * vm;

					// Jacobian entries
					if (wantJac)
					{
						for (int i = 0; i < 2 * wenoOrder - 1; ++i)
							// Note that we have an offset of -1 here (compared to the right cell face below), since
							// the reconstructed value depends on the previous stencil (which has now been moved by one cell)
							jac[(i - wenoOrder) * p.strideCell] -= static_cast<double>(invHCol) * p.reconstructionDerivatives[i];
					}
				}
				else if (wantRes)
				{
					// In the first cell we need to apply the boundary condition: inflow concentration
					resBulkComp[col * p.strideCell] -= invHCol * y[p.offsetToInlet + comp];
				}

				// Reconstruct concentration on this cell's right face
				// Non equidistant grid
				if (nonEqGrid)
				{
					if (wantJac)
						wenoOrder = p.reconstruction->template reconstruct<StateType, StencilType>(col, p.nCol, stencil, vm, p.reconstructionDerivatives, *p.cellFaces);
					else
						wenoOrder = p.reconstruction->template reconstruct<StateType, StencilType>(col, p.nCol, stencil, vm, *p.cellFaces);
				}
				// Equadistant grid
				else if (wantJac)
					wenoOrder = p.reconstruction->template reconstruct<StateType, StencilType>(col, p.nCol, stencil, vm, p.reconstructionDerivatives);
				else
					wenoOrder = p.reconstruction->template reconstruct<StateType, StencilType>(col, p.nCol, stencil, vm);

				// Right side
				if (wantRes)
					resBulkComp[col * p.strideCell] += invHCol * vm;
				// Jacobian entries
				if (wantJac)
				{
					for (int i = 0; i < 2 * wenoOrder - 1; ++i)
						jac[(i - wenoOrder + 1) * p.strideCell] += static_cast<double>(invHCol) * p.reconstructionDerivatives[i];
				}

				// Update stencil
				const unsigned int shift = std::max(p.reconstruction->order(), 2);
				if (cadet_likely(col + shift < p.nCol))
					stencil.advance(yBulkComp[(col + shift) * p.strideCell]);
				else
					stencil.advance(0.0);

				if (wantJac)
				{
					if (cadet_likely(col < p.nCol - 1))
						jac += p.strideCell;
				}
			}
		}

		// Film diffusion with flux into beads is added in residualFlux() function

		return 0;
	}

	template <typename StateType, typename ResidualType, typename ParamType, typename ReconstrType, typename RowIteratorType, bool wantJac, bool wantRes = true>
	int residualBackwardsAxialFlow(const SimulationTime& simTime, StateType const* y, double const* yDot, ResidualType* res, RowIteratorType jacBegin, const AxialFlowParameters<ParamType, ReconstrType>& p)
	{
		const bool nonEqGrid = p.cellFaces && !p.gridEquidistant;
		const ParamType h2 = p.h * p.h;

		// The stencil caches parts of the state vector for better spatial coherence
		typedef CachingStencil<StateType, ArrayPool> StencilType;
		StencilType stencil(std::max(p.reconstruction->stencilSize(), 3u), *p.stencilMemory, std::max(p.reconstruction->order() - 1, 1));

		// The RowIterator is always centered on the main diagonal.
		// This means that jac[0] is the main diagonal, jac[-1] is the first lower diagonal,
		// and jac[1] is the first upper diagonal. We can also access the rows from left to
		// right beginning with the last lower diagonal moving towards the main diagonal and
		// continuing to the last upper diagonal by using the native() method.
		RowIteratorType jac;

		ResidualType* const resBulk = wantRes ? res + p.offsetToBulk : nullptr;
		StateType const* const yBulk = y + p.offsetToBulk;
		const ParamType colLength = nonEqGrid ? static_cast<ParamType>(p.cellFaces->back()) : static_cast<ParamType>(0.0);

		for (unsigned int comp = 0; comp < p.nComp; ++comp)
		{
			if (wantJac)
				jac = jacBegin + p.strideCell * (p.nCol - 1) + comp;

			ResidualType* const resBulkComp = wantRes ? resBulk + comp : nullptr;
			StateType const* const yBulkComp = yBulk + comp;

			// Add time derivative to each cell
			if (yDot && wantRes)
			{
				double const* const yDotBulkComp = yDot + p.offsetToBulk + comp;
				for (unsigned int col = 0; col < p.nCol; ++col)
					resBulkComp[col * p.strideCell] = yDotBulkComp[col * p.strideCell];
			}
			else if (wantRes)
			{
				for (unsigned int col = 0; col < p.nCol; ++col)
					resBulkComp[col * p.strideCell] = 0.0;
			}

			// Fill stencil (left side with zeros, right side with states)
			for (int i = -std::max(p.reconstruction->order(), 2) + 1; i < 0; ++i)
				stencil[i] = 0.0;
			for (int i = 0; i < std::max(p.reconstruction->order(), 2); ++i)
				stencil[i] = yBulkComp[(p.nCol - static_cast<unsigned int>(i) - 1) * p.strideCell];

			// Reset reconstruction output
			StateType vm(0.0); // reconstructed value
			if (wantJac)
				std::fill(p.reconstructionDerivatives, p.reconstructionDerivatives + p.reconstruction->stencilSize(), 0.0);

			int wenoOrder = 0;
			const ParamType d_ax = static_cast<ParamType>(p.d_ax[comp]);

			// Iterate over all cells (backwards)
			// Note that col wraps around to unsigned int's maximum value after 0
			for (unsigned int col = p.nCol - 1; col < p.nCol; --col)
			{
				const ParamType hCol = nonEqGrid ? static_cast<ParamType>((*p.cellFaces)[col + 1] - (*p.cellFaces)[col]) : p.h;
				const ParamType invHCol = static_cast<ParamType>(p.u) / hCol;

				// ------------------- Dispersion -------------------

				// Right side, leave out if we're in the first cell (boundary condition)
				if (cadet_likely(col < p.nCol - 1))
				{
					const ParamType hRight = nonEqGrid ? static_cast<ParamType>((*p.cellFaces)[col + 2] - (*p.cellFaces)[col + 1]) : p.h;
					const ParamType deltaZ = nonEqGrid ? static_cast<ParamType>(0.5) * (hCol + hRight) : p.h;
					const double relCoord = nonEqGrid ?
						static_cast<double>((*p.cellFaces)[col + 1]) / static_cast<double>(colLength) :
						static_cast<double>(col + 1) / p.nCol;
					const ParamType d_ax_right = d_ax * p.parDep->getValue(p.model, ColumnPosition{ relCoord, 0.0, 0.0 }, comp, ParTypeIndep, BoundStateIndep, static_cast<ParamType>(p.u));
					const ParamType coeff = nonEqGrid ? d_ax_right / (hCol * deltaZ) : d_ax_right / h2;

					if (wantRes)
						resBulkComp[col * p.strideCell] -= coeff * (stencil[1] - stencil[0]);
					// Jacobian entries
					if (wantJac)
					{
						const double coeffJac = static_cast<double>(coeff);
						jac[0] += coeffJac;
						jac[p.strideCell] -= coeffJac;
					}
				}

				// Left side, leave out if we're in the last cell (boundary condition)
				if (cadet_likely(col > 0))
				{
					const ParamType hLeft = nonEqGrid ? static_cast<ParamType>((*p.cellFaces)[col] - (*p.cellFaces)[col - 1]) : p.h;
					const ParamType deltaZ = nonEqGrid ? static_cast<ParamType>(0.5) * (hLeft + hCol) : p.h;
					const double relCoord = nonEqGrid ?
						static_cast<double>((*p.cellFaces)[col]) / static_cast<double>(colLength) :
						static_cast<double>(col) / p.nCol;
					const ParamType d_ax_left = d_ax * p.parDep->getValue(p.model, ColumnPosition{ relCoord, 0.0, 0.0 }, comp, ParTypeIndep, BoundStateIndep, static_cast<ParamType>(p.u));
					const ParamType coeff = nonEqGrid ? d_ax_left / (hCol * deltaZ) : d_ax_left / h2;

					if (wantRes)
						resBulkComp[col * p.strideCell] -= coeff * (stencil[-1] - stencil[0]);
					// Jacobian entries
					if (wantJac)
					{
						const double coeffJac = static_cast<double>(coeff);
						jac[0] += coeffJac;
						jac[-p.strideCell] -= coeffJac;
					}
				}

				// ------------------- Convection -------------------

				// Add convection through this cell's right face
				if (cadet_likely(col < p.nCol - 1))
				{
					// Remember that vm still contains the reconstructed value of the previous 
					// cell's *left* face, which is identical to this cell's *right* face!
					if (wantRes)
						resBulkComp[col * p.strideCell] += invHCol * vm;

					// Jacobian entries
					if (wantJac)
					{
						for (int i = 0; i < 2 * wenoOrder - 1; ++i)
							// Note that we have an offset of +1 here (compared to the left cell face below), since
							// the reconstructed value depends on the previous stencil (which has now been moved by one cell)
							jac[(wenoOrder - i) * p.strideCell] += static_cast<double>(invHCol) * p.reconstructionDerivatives[i];
					}
				}
				else if (wantRes)
				{
					// In the last cell (z = L) we need to apply the boundary condition: inflow concentration
					resBulkComp[col * p.strideCell] += invHCol * y[p.offsetToInlet + comp];
				}

				// Reconstruct concentration on this cell's left face
				if (nonEqGrid)
				{
					const ReverseFaceAccessor<std::vector<active>> reverseFaces{ *p.cellFaces };
					const unsigned int flowCellIdx = p.nCol - 1 - col;
					if (wantJac)
						wenoOrder = p.reconstruction->template reconstruct<StateType, StencilType>(flowCellIdx, p.nCol, stencil, vm, p.reconstructionDerivatives, reverseFaces);
					else
						wenoOrder = p.reconstruction->template reconstruct<StateType, StencilType>(flowCellIdx, p.nCol, stencil, vm, reverseFaces);
				}
				else if (wantJac)
					wenoOrder = p.reconstruction->template reconstruct<StateType, StencilType>(col, p.nCol, stencil, vm, p.reconstructionDerivatives);
				else
					wenoOrder = p.reconstruction->template reconstruct<StateType, StencilType>(col, p.nCol, stencil, vm);

				// Left face
				if (wantRes)
					resBulkComp[col * p.strideCell] -= invHCol * vm;
				// Jacobian entries
				if (wantJac)
				{
					for (int i = 0; i < 2 * wenoOrder - 1; ++i)
						jac[(wenoOrder - i - 1) * p.strideCell] -= static_cast<double>(invHCol) * p.reconstructionDerivatives[i];
				}

				// Update stencil (be careful because of wrap-around, might cause reading memory very far away [although never used])
				const unsigned int shift = std::max(p.reconstruction->order(), 2);
				if (cadet_likely(col - shift < p.nCol))
					stencil.advance(yBulkComp[(col - shift) * p.strideCell]);
				else
					stencil.advance(0.0);

				if (wantJac)
				{
					if (cadet_likely(col > 0))
						jac -= p.strideCell;
				}
			}
		}

		// Film diffusion with flux into beads is added in residualFlux() function

		return 0;
	}

} // namespace impl


template <typename StateType, typename ResidualType, typename ParamType, typename ReconstrType, typename RowIteratorType, bool wantJac, bool wantRes = true>
int residualKernelAxial(const SimulationTime& simTime, StateType const* y, double const* yDot, ResidualType* res, RowIteratorType jacBegin, const AxialFlowParameters<ParamType, ReconstrType>& p)
{
	if (p.u >= 0.0)
		return impl::residualForwardsAxialFlow<StateType, ResidualType, ParamType, ReconstrType, RowIteratorType, wantJac, wantRes>(simTime, y, yDot, res, jacBegin, p);
	else
		return impl::residualBackwardsAxialFlow<StateType, ResidualType, ParamType, ReconstrType, RowIteratorType, wantJac, wantRes>(simTime, y, yDot, res, jacBegin, p);
}

void sparsityPatternAxial(linalg::SparsityPatternRowIterator itBegin, unsigned int nComp, unsigned int nCol, int strideCell, double u, Weno& weno);
void sparsityPatternAxial(linalg::SparsityPatternRowIterator itBegin, unsigned int nComp, unsigned int nCol, int strideCell, double u, HighResolutionKoren& koren);

} // namespace convdisp
} // namespace parts
} // namespace model
} // namespace cadet

#endif  // LIBCADET_AXIALCONVECTIONDISPERSIONKERNEL_HPP_
