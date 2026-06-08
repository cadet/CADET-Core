// SPDX-License-Identifier: AGPL-3.0-or-later
// =================================================================================
//  CADET
//
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Affero General Public
//  License v3.0 (or, at your option, any later version).
// =================================================================================

/**
 * @file 
 * Implements the kernel of the frustum convection dispersion transport operator.
 */

#include "model/parts/FrustumConvectionDispersionKernel.hpp"
#include "Weno.hpp"
#include "HighResKoren.hpp"
#include "linalg/CompressedSparseMatrix.hpp"

namespace cadet
{

namespace model
{

namespace parts
{

namespace convdisp
{

namespace impl
{

	class DummyStencil
	{
	public:
		DummyStencil() { }
		inline double operator[](const int idx) const { return 0.0; }
	};

	template <typename ReconstrType>
	void sparsityPatternFrustum(linalg::SparsityPatternRowIterator itBegin, unsigned int nComp, unsigned int nCol, int strideCell, double u, ReconstrType& reconstruction)
	{
		impl::DummyStencil stencil;

		if (u >= 0.0)
		{
			for (unsigned int comp = 0; comp < nComp; ++comp)
			{
				linalg::SparsityPatternRowIterator jac = itBegin + comp;
				double vm = 0.0;
				int recOrder = 0;

				for (unsigned int col = 0; col < nCol; ++col)
				{
					jac.centered(0);

					if (cadet_likely(col < nCol - 1))
						jac.centered(strideCell);
					if (cadet_likely(col > 0))
						jac.centered(-strideCell);

					if (cadet_likely(col > 0))
					{
						for (int i = 0; i < 2 * recOrder - 1; ++i)
							jac.centered((i - recOrder) * strideCell);
					}

					recOrder = reconstruction.template reconstruct<double, impl::DummyStencil>(col, nCol, stencil, vm);

					for (int i = 0; i < 2 * recOrder - 1; ++i)
						jac.centered((i - recOrder + 1) * strideCell);

					if (cadet_likely(col < nCol - 1))
						jac += strideCell;
				}
			}
		}
		else
		{
			for (unsigned int comp = 0; comp < nComp; ++comp)
			{
				linalg::SparsityPatternRowIterator jac = itBegin + strideCell * (nCol - 1) + comp;
				double vm = 0.0;
				int recOrder = 0;

				for (unsigned int col = nCol - 1; col < nCol; --col)
				{
					jac.centered(0);

					if (cadet_likely(col < nCol - 1))
						jac.centered(strideCell);
					if (cadet_likely(col > 0))
						jac.centered(-strideCell);

					if (cadet_likely(col < nCol - 1))
					{
						for (int i = 0; i < 2 * recOrder - 1; ++i)
							jac.centered((recOrder - i) * strideCell);
					}

					recOrder = reconstruction.template reconstruct<double, impl::DummyStencil>(col, nCol, stencil, vm);

					for (int i = 0; i < 2 * recOrder - 1; ++i)
						jac.centered((recOrder - i - 1) * strideCell);

					if (cadet_likely(col > 0))
						jac -= strideCell;
				}
			}
		}
	}

} // namespace impl

void sparsityPatternFrustum(linalg::SparsityPatternRowIterator itBegin, unsigned int nComp, unsigned int nCol, int strideCell, double u, Weno& weno)
{
	impl::sparsityPatternFrustum(itBegin, nComp, nCol, strideCell, u, weno);
}

void sparsityPatternFrustum(linalg::SparsityPatternRowIterator itBegin, unsigned int nComp, unsigned int nCol, int strideCell, double u, HighResolutionKoren& koren)
{
	impl::sparsityPatternFrustum(itBegin, nComp, nCol, strideCell, u, koren);
}

} // namespace convdisp
} // namespace parts
} // namespace model
} // namespace cadet
