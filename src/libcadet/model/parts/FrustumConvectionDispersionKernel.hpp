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
 * Implements the kernel of the frustum convection dispersion transport operator.
 */

#ifndef LIBCADET_FRUSTUMCONVECTIONDISPERSIONKERNEL_HPP_
#define LIBCADET_FRUSTUMCONVECTIONDISPERSIONKERNEL_HPP_

#include "AutoDiff.hpp"
#include "Memory.hpp"
#include "Weno.hpp"
#include "HighResKoren.hpp"
#include "Stencil.hpp"
#include "linalg/CompressedSparseMatrix.hpp"
#include "SimulationTypes.hpp"
#include "model/ParameterDependence.hpp"
#include "model/UnitOperation.hpp"

#include <algorithm>
#include <cmath>

namespace cadet
{

class IModel;

namespace model
{

namespace parts
{

namespace convdisp
{

template <typename T, typename Reconstruction_t>
struct FrustumFlowParameters
{
	T u;
	T colLength;
	active const* disp;
	active const* cellCenters; //!< Coordinates of cell midpoints
	active const* cellFaces; //!< Coordinates of cell faces
	active const* cellFaceRadiusSq; //!< Squared radii at cell faces
	active const* cellVolume; //!< Volumes of cells
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
	bool gridEquidistant; //!< Determines whether the grid is equidistant
	std::vector<active> const* gridFaces; //!< Positions of the cell faces for non-equidistant grids
};


namespace impl
{
	template <class FaceContainerType>
	struct ReverseFaceAccessorFrustum
	{
		const FaceContainerType& faces;
		inline auto operator[](unsigned int idx) const -> decltype(faces[0]) { return faces[faces.size() - 1 - idx]; }
	};

	template <typename StateType, typename ResidualType, typename ParamType, typename ReconstrType, typename RowIteratorType, bool wantJac, bool wantRes = true>
	int residualForwardsFrustumFlow(const SimulationTime& simTime, StateType const* y, double const* yDot, ResidualType* res, RowIteratorType jacBegin, const FrustumFlowParameters<ParamType, ReconstrType>& p)
	{
		const bool nonEqGrid = p.gridFaces && !p.gridEquidistant;

		// The stencil caches parts of the state vector for better spatial coherence
		typedef CachingStencil<StateType, ArrayPool> StencilType;
		StencilType stencil(std::max(p.reconstruction->stencilSize(), 3u), *p.stencilMemory, std::max(p.reconstruction->order() - 1, 1));

		RowIteratorType jac;

		ResidualType* const resBulk = wantRes ? res + p.offsetToBulk : nullptr;
		StateType const* const yBulk = y + p.offsetToBulk;

		for (unsigned int comp = 0; comp < p.nComp; ++comp)
		{
			if (wantJac)
				jac = jacBegin + comp;

			ResidualType* const resBulkComp = wantRes ? resBulk + comp : nullptr;
			StateType const* const yBulkComp = yBulk + comp;

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

			for (int i = -std::max(p.reconstruction->order(), 2) + 1; i < 0; ++i)
				stencil[i] = 0.0;
			for (int i = 0; i < std::max(p.reconstruction->order(), 2); ++i)
				stencil[i] = yBulkComp[i * p.strideCell];

			StateType vm(0.0);
			if (wantJac)
				std::fill(p.reconstructionDerivatives, p.reconstructionDerivatives + p.reconstruction->stencilSize(), 0.0);

			int recOrder = 0;
			const ParamType disp = static_cast<ParamType>(p.disp[comp]);
			const ParamType colLength = static_cast<ParamType>(p.colLength);
			const ParamType pi = static_cast<ParamType>(3.1415926535897932384626434);

			for (unsigned int col = 0; col < p.nCol; ++col)
			{
				const ParamType preFac = pi / static_cast<ParamType>(p.cellVolume[col]);
				const ParamType convCoeff = preFac * static_cast<ParamType>(p.u);

				if (cadet_likely(col < p.nCol - 1))
				{
					const double relCoord = static_cast<double>(p.cellFaces[col + 1]) / static_cast<double>(colLength);
					const ParamType dispRight = disp * p.parDep->getValue(p.model, ColumnPosition{relCoord, 0.0, 0.0}, comp, ParTypeIndep, BoundStateIndep,
						static_cast<ParamType>(p.u) / static_cast<ParamType>(p.cellFaceRadiusSq[col + 1]));
					const ParamType coeff = dispRight * preFac * static_cast<ParamType>(p.cellFaceRadiusSq[col + 1]) / (static_cast<ParamType>(p.cellCenters[col + 1]) - static_cast<ParamType>(p.cellCenters[col]));
					if (wantRes)
						resBulkComp[col * p.strideCell] -= coeff * (stencil[1] - stencil[0]);
					if (wantJac)
					{
						const double coeffJac = static_cast<double>(coeff);
						jac[0] += coeffJac;
						jac[p.strideCell] -= coeffJac;
					}
				}

				if (cadet_likely(col > 0))
				{
					const double relCoord = static_cast<double>(p.cellFaces[col]) / static_cast<double>(colLength);
					const ParamType dispLeft = disp * p.parDep->getValue(p.model, ColumnPosition{relCoord, 0.0, 0.0}, comp, ParTypeIndep, BoundStateIndep,
						static_cast<ParamType>(p.u) / static_cast<ParamType>(p.cellFaceRadiusSq[col]));
					const ParamType coeff = dispLeft * preFac * static_cast<ParamType>(p.cellFaceRadiusSq[col]) / (static_cast<ParamType>(p.cellCenters[col]) - static_cast<ParamType>(p.cellCenters[col - 1]));
					if (wantRes)
						resBulkComp[col * p.strideCell] -= coeff * (stencil[-1] - stencil[0]);
					if (wantJac)
					{
						const double coeffJac = static_cast<double>(coeff);
						jac[0] += coeffJac;
						jac[-p.strideCell] -= coeffJac;
					}
				}

				if (cadet_likely(col > 0))
				{
					if (wantRes)
						resBulkComp[col * p.strideCell] -= convCoeff * vm;
					if (wantJac)
					{
						for (int i = 0; i < 2 * recOrder - 1; ++i)
							jac[(i - recOrder) * p.strideCell] -= static_cast<double>(convCoeff) * p.reconstructionDerivatives[i];
					}
				}
				else if (wantRes)
				{
					resBulkComp[col * p.strideCell] -= convCoeff * y[p.offsetToInlet + comp];
				}

				if (nonEqGrid)
				{
					if (wantJac)
						recOrder = p.reconstruction->template reconstruct<StateType, StencilType>(col, p.nCol, stencil, vm, p.reconstructionDerivatives, *p.gridFaces);
					else
						recOrder = p.reconstruction->template reconstruct<StateType, StencilType>(col, p.nCol, stencil, vm, *p.gridFaces);
				}
				else
				{
					if (wantJac)
						recOrder = p.reconstruction->template reconstruct<StateType, StencilType>(col, p.nCol, stencil, vm, p.reconstructionDerivatives);
					else
						recOrder = p.reconstruction->template reconstruct<StateType, StencilType>(col, p.nCol, stencil, vm);
				}

				if (wantRes)
					resBulkComp[col * p.strideCell] += convCoeff * vm;
				if (wantJac)
				{
					for (int i = 0; i < 2 * recOrder - 1; ++i)
						jac[(i - recOrder + 1) * p.strideCell] += static_cast<double>(convCoeff) * p.reconstructionDerivatives[i];
				}

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

		return 0;
	}

	template <typename StateType, typename ResidualType, typename ParamType, typename ReconstrType, typename RowIteratorType, bool wantJac, bool wantRes = true>
	int residualBackwardsFrustumFlow(const SimulationTime& simTime, StateType const* y, double const* yDot, ResidualType* res, RowIteratorType jacBegin, const FrustumFlowParameters<ParamType, ReconstrType>& p)
	{
		const bool nonEqGrid = p.gridFaces && !p.gridEquidistant;

		typedef CachingStencil<StateType, ArrayPool> StencilType;
		StencilType stencil(std::max(p.reconstruction->stencilSize(), 3u), *p.stencilMemory, std::max(p.reconstruction->order() - 1, 1));

		RowIteratorType jac;

		ResidualType* const resBulk = wantRes ? res + p.offsetToBulk : nullptr;
		StateType const* const yBulk = y + p.offsetToBulk;

		for (unsigned int comp = 0; comp < p.nComp; ++comp)
		{
			if (wantJac)
				jac = jacBegin + p.strideCell * (p.nCol - 1) + comp;

			ResidualType* const resBulkComp = wantRes ? resBulk + comp : nullptr;
			StateType const* const yBulkComp = yBulk + comp;

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

			for (int i = -std::max(p.reconstruction->order(), 2) + 1; i < 0; ++i)
				stencil[i] = 0.0;
			for (int i = 0; i < std::max(p.reconstruction->order(), 2); ++i)
				stencil[i] = yBulkComp[(p.nCol - static_cast<unsigned int>(i) - 1) * p.strideCell];

			StateType vm(0.0);
			if (wantJac)
				std::fill(p.reconstructionDerivatives, p.reconstructionDerivatives + p.reconstruction->stencilSize(), 0.0);

			int recOrder = 0;
			const ParamType disp = static_cast<ParamType>(p.disp[comp]);
			const ParamType colLength = static_cast<ParamType>(p.colLength);
			const ParamType pi = static_cast<ParamType>(3.1415926535897932384626434);

			for (unsigned int col = p.nCol - 1; col < p.nCol; --col)
			{
				const ParamType preFac = pi / static_cast<ParamType>(p.cellVolume[col]);
				const ParamType convCoeff = preFac * static_cast<ParamType>(p.u);

				if (cadet_likely(col < p.nCol - 1))
				{
					const double relCoord = static_cast<double>(p.cellFaces[col + 1]) / static_cast<double>(colLength);
					const ParamType dispRight = disp * p.parDep->getValue(p.model, ColumnPosition{relCoord, 0.0, 0.0}, comp, ParTypeIndep, BoundStateIndep,
						static_cast<ParamType>(p.u) / static_cast<ParamType>(p.cellFaceRadiusSq[col + 1]));
					const ParamType coeff = dispRight * preFac * static_cast<ParamType>(p.cellFaceRadiusSq[col + 1]) / (static_cast<ParamType>(p.cellCenters[col + 1]) - static_cast<ParamType>(p.cellCenters[col]));
					if (wantRes)
						resBulkComp[col * p.strideCell] -= coeff * (stencil[-1] - stencil[0]);
					if (wantJac)
					{
						const double coeffJac = static_cast<double>(coeff);
						jac[0] += coeffJac;
						jac[p.strideCell] -= coeffJac;
					}
				}

				if (cadet_likely(col > 0))
				{
					const double relCoord = static_cast<double>(p.cellFaces[col]) / static_cast<double>(colLength);
					const ParamType dispLeft = disp * p.parDep->getValue(p.model, ColumnPosition{relCoord, 0.0, 0.0}, comp, ParTypeIndep, BoundStateIndep,
						static_cast<ParamType>(p.u) / static_cast<ParamType>(p.cellFaceRadiusSq[col]));
                    const ParamType coeff = dispLeft * preFac * static_cast<ParamType>(p.cellFaceRadiusSq[col]) / (static_cast<ParamType>(p.cellCenters[col]) - static_cast<ParamType>(p.cellCenters[col - 1]));
					if (wantRes)
						resBulkComp[col * p.strideCell] -= coeff * (stencil[1] - stencil[0]);
					if (wantJac)
					{
						const double coeffJac = static_cast<double>(coeff);
						jac[0] += coeffJac;
						jac[-p.strideCell] -= coeffJac;
					}
				}

				if (cadet_likely(col < p.nCol - 1))
				{
					if (wantRes)
						resBulkComp[col * p.strideCell] += convCoeff * vm;
					if (wantJac)
					{
						for (int i = 0; i < 2 * recOrder - 1; ++i)
							jac[(recOrder - i) * p.strideCell] += static_cast<double>(convCoeff) * p.reconstructionDerivatives[i];
					}
				}
				else if (wantRes)
				{
					resBulkComp[col * p.strideCell] += convCoeff * y[p.offsetToInlet + comp];
				}

				if (nonEqGrid)
				{
					const ReverseFaceAccessorFrustum<std::vector<active>> reverseFaces{*p.gridFaces};
					const unsigned int flowCellIdx = p.nCol - 1 - col;
					if (wantJac)
						recOrder = p.reconstruction->template reconstruct<StateType, StencilType>(flowCellIdx, p.nCol, stencil, vm, p.reconstructionDerivatives, reverseFaces);
					else
						recOrder = p.reconstruction->template reconstruct<StateType, StencilType>(flowCellIdx, p.nCol, stencil, vm, reverseFaces);
				}
				else
				{
					if (wantJac)
						recOrder = p.reconstruction->template reconstruct<StateType, StencilType>(col, p.nCol, stencil, vm, p.reconstructionDerivatives);
					else
						recOrder = p.reconstruction->template reconstruct<StateType, StencilType>(col, p.nCol, stencil, vm);
				}

				if (wantRes)
					resBulkComp[col * p.strideCell] -= convCoeff * vm;
				if (wantJac)
				{
					for (int i = 0; i < 2 * recOrder - 1; ++i)
						jac[(recOrder - i - 1) * p.strideCell] -= static_cast<double>(convCoeff) * p.reconstructionDerivatives[i];
				}

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

		return 0;
	}

} // namespace impl


template <typename StateType, typename ResidualType, typename ParamType, typename ReconstrType, typename RowIteratorType, bool wantJac, bool wantRes = true>
int residualKernelFrustum(const SimulationTime& simTime, StateType const* y, double const* yDot, ResidualType* res, RowIteratorType jacBegin, const FrustumFlowParameters<ParamType, ReconstrType>& p)
{
	if (p.u >= 0.0)
		return impl::residualForwardsFrustumFlow<StateType, ResidualType, ParamType, ReconstrType, RowIteratorType, wantJac, wantRes>(simTime, y, yDot, res, jacBegin, p);
	else
		return impl::residualBackwardsFrustumFlow<StateType, ResidualType, ParamType, ReconstrType, RowIteratorType, wantJac, wantRes>(simTime, y, yDot, res, jacBegin, p);
}

void sparsityPatternFrustum(linalg::SparsityPatternRowIterator itBegin, unsigned int nComp, unsigned int nCol, int strideCell, double u, Weno& weno);
void sparsityPatternFrustum(linalg::SparsityPatternRowIterator itBegin, unsigned int nComp, unsigned int nCol, int strideCell, double u, HighResolutionKoren& koren);

} // namespace convdisp
} // namespace parts
} // namespace model
} // namespace cadet

#endif  // LIBCADET_FRUSTUMCONVECTIONDISPERSIONKERNEL_HPP_
