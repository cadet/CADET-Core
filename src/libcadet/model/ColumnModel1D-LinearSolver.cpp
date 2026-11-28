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

#include "model/ColumnModel1D.hpp"
#include "model/BindingModel.hpp"
#include "model/parts/BindingCellKernel.hpp"
#include "AdUtils.hpp"

#include <algorithm>
#include <cmath>
#include <functional>

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include "ParallelSupport.hpp"

#ifdef CADET_PARALLELIZE
	#include <tbb/parallel_for.h>
	#include <tbb/flow_graph.h>

	typedef tbb::flow::continue_node< tbb::flow::continue_msg > node_t;
	typedef const tbb::flow::continue_msg & msg_t;
#endif

namespace cadet
{

namespace model
{

/**
 * @brief Computes the solution of the linear system involving the system Jacobian
 * @details The system \f[ \left( \frac{\partial F}{\partial y} + \alpha \frac{\partial F}{\partial \dot{y}} \right) x = b \f]
 *          has to be solved. The right hand side \f$ b \f$ is given by @p rhs, the Jacobians are evaluated at the
 *          point \f$(y, \dot{y})\f$ given by @p y and @p yDot. The residual @p res at this point, \f$ F(t, y, \dot{y}) \f$,
 *          may help with this. Error weights (see IDAS guide) are given in @p weight. The solution is returned in @p rhs.
 *
 *          The full Jacobian @f$ J = \left( \frac{\partial F}{\partial y} + \alpha \frac{\partial F}{\partial \dot{y}} \right) @f$ is given by
 *          @f[ \begin{align}
				J =
				\left[\begin{array}{c|ccc|c}
					 J_0 & & & & J_{0,f} \\
					\hline
					 & J_1 & & & J_{1,f} \\
					 & & \ddots & & \vdots \\
					 & & & J_{N_z} & J_{N_z,f} \\
					\hline
					 J_{f,0} & J_{f,1} & \dots & J_{f,N_z} & J_f
				\end{array}\right].
			\end{align} @f]
 *          By decomposing the Jacobian @f$ J @f$ into @f$ J = LU @f$, we get
 *          @f[ \begin{align}
				L &= \left[\begin{array}{c|ccc|c}
					  J_0     &         &        &           &   \\
					\hline
							  & J_1     &        &           &   \\
							  &         & \ddots &           &   \\
							  &         &        & J_{N_z}   &   \\
					\hline
					  J_{f,0} & J_{f,1} & \dots & J_{f,N_z} & I
				\end{array}\right], \\
				U &= \left[\begin{array}{c|ccc|c}
					  I &   &        &   & J_0^{-1} \, J_{0,f}       \\
					\hline
						& I &        &   & J_1^{-1} \, J_{1,f}       \\
						&   & \ddots &   & \vdots                    \\
						&   &        & I & J_{N_z}^{-1} \, J_{N_z,f} \\
					\hline
						&   &        &   & S
				\end{array}\right].
			\end{align} @f]
 *          Here, the Schur-complement @f$ S @f$ is given by
 *          @f[ \begin{align}
				S = J_f - J_{f,0} \, J_0^{-1} \, J_{0,f} - \sum_{p=1}^{N_z}{J_{f,p} \, J_p^{-1} \, J_{p,f}}.
			\end{align} @f]
 *          Note that @f$ J_f = I @f$ is the identity matrix and that the off-diagonal blocks @f$ J_{i,f} @f$
 *          and @f$ J_{f,i} @f$ for @f$ i = 0, \dots, N_{z} @f$ are sparse.
 *
 *          Exploiting the decomposition, the solution procedure @f$ x = J^{-1}b = \left( LU \right)^{-1}b = U^{-1} L^{-1} b @f$
 *          works as follows:
 *              -# Factorize the diagonal blocks @f$ J_0, \dots, J_{N_z} @f$
 *              -# Solve @f$ y = L^{-1} b @f$ by forward substitution. This is accomplished by first solving the diagonal
 *                 blocks independently, that is,
 *                 @f[ y_i = J_{i}^{-1} b_i. @f]
 *                 Then, calculate the flux-part @f$ y_f @f$ by substituting in the already calculated solutions @f$ y_i @f$:
 *                 @f[ y_f = b_f - \sum_{i=0}^{N_z} J_{f,i} y_i. @f]
 *              -# Solve the Schur-complement @f$ S x_f = y_f @f$ using an iterative method that only requires
 *                 matrix-vector products. The already inverted diagonal blocks @f$ J_i^{-1} @f$ come in handy here.
 *              -# Solve the rest of the @f$ U x = y @f$ system by backward substitution. To be more precise, compute
 *                 @f[ x_i = y_i - J_i^{-1} J_{i,f} y_f. @f]
 *
 *
 * @param [in] t Current time point
 * @param [in] alpha Value of \f$ \alpha \f$ (arises from BDF time discretization)
 * @param [in] outerTol Error tolerance for the solution of the linear system from outer Newton iteration
 * @param [in,out] rhs On entry the right hand side of the linear equation system, on exit the solution
 * @param [in] weight Vector with error weights
 * @param [in] simState State of the simulation (state vector and its time derivatives) at which the Jacobian is evaluated
 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
 */
template <typename ConvDispOperator>
int ColumnModel1D<ConvDispOperator>::linearSolve(double t, double alpha, double outerTol, double* const rhs, double const* const weight,
	const ConstSimulationState& simState)
{
	BENCH_SCOPE(_timerLinearSolve);

	Indexer idxr(_disc);

	Eigen::Map<VectorXd> r(rhs, numDofs()); // map rhs to Eigen object

	// ==== Step 1: Factorize diagonal Jacobian blocks

	// Factorize partial Jacobians only if required

	if (_factorizeJacobian)
	{

		// Assemble and factorize discretized bulk Jacobian
		assembleDiscretizedGlobalJacobian(alpha, idxr);

		if (_blockSolver)
		{
			if (cadet_unlikely(!blockFactorize(idxr)))
			{
				LOG(Error) << "Factorize() failed";
				return 1;
			}
		}
		else
		{
			_linearSolver->factorize(_globalJacDisc, idxr.offsetC(), numPureDofs());

			if (cadet_unlikely(_linearSolver->info() != Eigen::Success))
			{
				LOG(Error) << "Factorize() failed";
				return 1;
			}
		}

		// Do not factorize again at next call without changed Jacobians
		_factorizeJacobian = false;
	}

	// ==== Step 1.5: Solve J c_uo = b_uo - A * c_in = b_uo - A*b_in

	// Handle inlet DOFs:
	// Inlet at z = 0 for forward flow, at z = L for backward flow. The inlet block ends at the last
	// discrete point under backward flow; how many points it covers depends on the transport
	// operator, see ColumnModel1D::multiplyWithJacobian().
	const unsigned int offInlet = _convDispOp.forwardFlow() ? 0u : (_disc.nPoints - static_cast<unsigned int>(_jacInlet.rows())) * idxr.strideColNode();

	const auto& cm = _reaction.conservedMoieties("liquid");
	if (cm.hasRapidEquilibrium())
	{
		const auto& L = cm.conservedMoietyMatrix();
        for (unsigned int moiety = 0; moiety < cm.numMoieties(); ++moiety)
        {
            double inletValue = 0.0;
            for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
                inletValue += L(moiety, comp) * r[comp];

            for (int node = 0; node < static_cast<int>(_jacInlet.rows()); ++node)
            {
                r[idxr.offsetC() + offInlet + moiety * idxr.strideColComp() + node * idxr.strideColNode()] -= _jacInlet(node, 0) * inletValue;
            }
        }

	}
	else
	{
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			for (int node = 0; node < static_cast<int>(_jacInlet.rows()); ++node)
			{
				r[idxr.offsetC() + offInlet + comp * idxr.strideColComp() + node * idxr.strideColNode()] -= _jacInlet(node, 0) * r[comp];
			}
		}
	}

	// ==== Step 2: Solve system of pure DOFs
	// The result is stored in rhs (in-place solution)

	if (_blockSolver)
	{
		if (cadet_unlikely(!blockSolve(r.data() + idxr.offsetC(), idxr)))
		{
			LOG(Error) << "Solve() failed";
			return 1;
		}
	}
	else
	{
		_linearSolver->solve(r.data() + idxr.offsetC(), numPureDofs());

		if (cadet_unlikely(_linearSolver->info() != Eigen::Success))
		{
			LOG(Error) << "Solve() failed";
			return 1;
		}
	}

	// The full solution is now stored in rhs
	return 0;
}

/**
 * @brief Prepares the block elimination of the particle states
 * @details Determines which pure DOF belongs to which particle block and how wide the bulk and
 *          particle blocks are, and allocates the matrices the elimination works on. Called once
 *          the discretization is known, since none of this changes during a simulation.
 * @param [in] idxr Indexer
 */
template <typename ConvDispOperator>
void ColumnModel1D<ConvDispOperator>::setupBlockSolver(const Indexer& idxr)
{
	const int nPureDofs = static_cast<int>(numPureDofs());
	const int nBulk = static_cast<int>(_disc.nPoints) * idxr.strideColNode();

	_blkBegin.clear();
	_blkSize.clear();
	_blkNode.clear();

	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		for (unsigned int node = 0; node < _disc.nPoints; ++node)
		{
			_blkBegin.push_back(idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{node}) - idxr.offsetC());
			_blkSize.push_back(idxr.strideParBlock(type));
			_blkNode.push_back(static_cast<int>(node));
		}
	}

	const int nBlocks = static_cast<int>(_blkBegin.size());

	_blkOwner.assign(nPureDofs, -1);
	for (int blk = 0; blk < nBlocks; ++blk)
	{
		for (int i = 0; i < _blkSize[blk]; ++i)
			_blkOwner[_blkBegin[blk] + i] = blk;
	}

	// The bandwidths follow from the sparsity pattern, which is fixed once the discretization is
	_bulkBandwidth = 0;
	_parBandwidth = 0;
	int maxBlockSize = 0;

	for (int row = 0; row < nPureDofs; ++row)
	{
		const int rowBlk = _blkOwner[row];

		for (Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(_globalJac, row + idxr.offsetC()); it; ++it)
		{
			const int col = static_cast<int>(it.col()) - idxr.offsetC();
			if ((col < 0) || (col >= nPureDofs))
				continue;

			if ((rowBlk < 0) && (_blkOwner[col] < 0))
				_bulkBandwidth = std::max(_bulkBandwidth, std::abs(col - row));
			else if ((rowBlk >= 0) && (_blkOwner[col] == rowBlk))
				_parBandwidth = std::max(_parBandwidth, std::abs(col - row));
		}
	}

	for (int blk = 0; blk < nBlocks; ++blk)
		maxBlockSize = std::max(maxBlockSize, _blkSize[blk]);

	_schurBulk.resize(nBulk, _bulkBandwidth, _bulkBandwidth);

	_blkJac.resize(nBlocks);
	_blkToBulk.resize(nBlocks);
	_bulkToBlk.resize(nBlocks);
	_blkSchurFactor.resize(nBlocks);

	for (int blk = 0; blk < nBlocks; ++blk)
	{
		_blkJac[blk].resize(_blkSize[blk], _parBandwidth, _parBandwidth);
		_blkToBulk[blk].setZero(idxr.strideColNode(), _blkSize[blk]);
		_bulkToBlk[blk].setZero(_blkSize[blk], idxr.strideColNode());
		_blkSchurFactor[blk].setZero(_blkSize[blk], idxr.strideColNode());
	}

	_blkScratch.resize(maxBlockSize);
}

/**
 * @brief Eliminates the particle blocks and factorizes the remaining bulk system
 * @details Splits the assembled Jacobian into the bulk block, one block per particle, and the
 *          film transfer couplings between them. Each particle block is factorized on its own and
 *          applied to its coupling, and the resulting Schur complement is subtracted from the bulk
 *          block before that is factorized. Since every coupling touches a single axial node, the
 *          Schur complement only changes the diagonal blocks of the bulk Jacobian and the bulk
 *          stays banded.
 * @param [in] idxr Indexer
 * @return @c true if all factorizations succeeded, otherwise @c false
 */
template <typename ConvDispOperator>
bool ColumnModel1D<ConvDispOperator>::blockFactorize(const Indexer& idxr)
{
	const int nPureDofs = static_cast<int>(numPureDofs());
	const int nBlocks = static_cast<int>(_blkBegin.size());
	const int strideNode = idxr.strideColNode();

	_schurBulk.setAll(0.0);
	for (int blk = 0; blk < nBlocks; ++blk)
	{
		_blkJac[blk].setAll(0.0);
		_blkToBulk[blk].setZero();
		_bulkToBlk[blk].setZero();
	}

	for (int row = 0; row < nPureDofs; ++row)
	{
		const int rowBlk = _blkOwner[row];

		for (Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(_globalJacDisc, row + idxr.offsetC()); it; ++it)
		{
			const int col = static_cast<int>(it.col()) - idxr.offsetC();
			if ((col < 0) || (col >= nPureDofs))
				continue;

			const int colBlk = _blkOwner[col];
			const double val = it.value();

			if (rowBlk < 0)
			{
				if (colBlk < 0)
					_schurBulk.row(row)[col - row] = val;
				else
					_blkToBulk[colBlk](row - _blkNode[colBlk] * strideNode, col - _blkBegin[colBlk]) = val;
			}
			else if (colBlk < 0)
				_bulkToBlk[rowBlk](row - _blkBegin[rowBlk], col - _blkNode[rowBlk] * strideNode) = val;
			else
				_blkJac[rowBlk].row(row - _blkBegin[rowBlk])[col - row] = val;
		}
	}

	for (int blk = 0; blk < nBlocks; ++blk)
	{
		if (cadet_unlikely(!_blkJac[blk].factorize()))
			return false;

		// _blkSchurFactor = J_p^{-1} J_{p,bulk}, reused by the back substitution in blockSolve()
		for (int comp = 0; comp < strideNode; ++comp)
		{
			for (int i = 0; i < _blkSize[blk]; ++i)
				_blkScratch[i] = _bulkToBlk[blk](i, comp);

			if (cadet_unlikely(!_blkJac[blk].solve(_blkScratch.data())))
				return false;

			for (int i = 0; i < _blkSize[blk]; ++i)
				_blkSchurFactor[blk](i, comp) = _blkScratch[i];
		}

		const Eigen::MatrixXd update = _blkToBulk[blk] * _blkSchurFactor[blk];
		const int base = _blkNode[blk] * strideNode;

		for (int row = 0; row < strideNode; ++row)
		{
			for (int col = 0; col < strideNode; ++col)
				_schurBulk.row(base + row)[col - row] -= update(row, col);
		}
	}

	return _schurBulk.factorize();
}

/**
 * @brief Solves the system using the factors prepared by blockFactorize()
 * @param [in,out] rhs On entry the right hand side of the pure DOFs, on exit their solution
 * @param [in] idxr Indexer
 * @return @c true if all solves succeeded, otherwise @c false
 */
template <typename ConvDispOperator>
bool ColumnModel1D<ConvDispOperator>::blockSolve(double* const rhs, const Indexer& idxr)
{
	const int nBlocks = static_cast<int>(_blkBegin.size());
	const int strideNode = idxr.strideColNode();

	// Forward substitution: take the particle contributions over to the bulk right hand side
	for (int blk = 0; blk < nBlocks; ++blk)
	{
		double* const blkRhs = rhs + _blkBegin[blk];

		if (cadet_unlikely(!_blkJac[blk].solve(blkRhs)))
			return false;

		double* const bulkRhs = rhs + _blkNode[blk] * strideNode;
		for (int row = 0; row < strideNode; ++row)
		{
			double sum = 0.0;
			for (int i = 0; i < _blkSize[blk]; ++i)
				sum += _blkToBulk[blk](row, i) * blkRhs[i];

			bulkRhs[row] -= sum;
		}
	}

	if (cadet_unlikely(!_schurBulk.solve(rhs)))
		return false;

	// Back substitution: the bulk solution corrects the particle states
	for (int blk = 0; blk < nBlocks; ++blk)
	{
		double* const blkRhs = rhs + _blkBegin[blk];
		double const* const bulkSol = rhs + _blkNode[blk] * strideNode;

		for (int i = 0; i < _blkSize[blk]; ++i)
		{
			double sum = 0.0;
			for (int comp = 0; comp < strideNode; ++comp)
				sum += _blkSchurFactor[blk](i, comp) * bulkSol[comp];

			blkRhs[i] -= sum;
		}
	}

	return true;
}

/**
 * @brief Assembles bulk Jacobian @f$ J_i @f$ (@f$ i > 0 @f$) of the time-discretized equations
 * @details The system \f[ \left( \frac{\partial F}{\partial y} + \alpha \frac{\partial F}{\partial \dot{y}} \right) x = b \f]
 *          has to be solved. The system Jacobian of the original equations,
 *          \f[ \frac{\partial F}{\partial y}, \f]
 *          is already computed (by AD or manually in residualImpl() with @c wantJac = true). This function is responsible
 *          for adding
 *          \f[ \alpha \frac{\partial F}{\partial \dot{y}} \f]
 *          to the system Jacobian, which yields the Jacobian of the time-discretized equations
 *          \f[ F\left(t, y_0, \sum_{k=0}^N \alpha_k y_k \right) = 0 \f]
 *          when a BDF method is used. The time integrator needs to solve this equation for @f$ y_0 @f$, which requires
 *          the solution of the linear system mentioned above (@f$ \alpha_0 = \alpha @f$ given in @p alpha).
 *
 * @param [in] alpha Value of \f$ \alpha \f$ (arises from BDF time discretization)
 */
template <typename ConvDispOperator>
void ColumnModel1D<ConvDispOperator>::assembleDiscretizedGlobalJacobian(double alpha, Indexer idxr) {

	/* add static (per section) jacobian without inlet */
	_globalJacDisc = _globalJac;

	// Add time derivatives to particle shells
	for (unsigned int parType = 0; parType < _disc.nParType; parType++) {
		for (unsigned int colNode = 0; colNode < _disc.nPoints; colNode++) {

			linalg::BandedEigenSparseRowIterator jac(_globalJacDisc, idxr.offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode }));

			addTimeDerivativeToJacobianParticleShell(jac, idxr, alpha, parType);

			// compute time derivative of remaining points
			// Iterator jac has already been advanced to next shell
			for (unsigned int j = 1; j < _disc.nParPoints[parType]; ++j)
			{
				addTimeDerivativeToJacobianParticleShell(jac, idxr, alpha, parType);
				// Iterator jac has already been advanced to next shell
			}
		}
	}

	const auto& cm = _reaction.conservedMoieties("liquid");
	if (cm.hasRapidEquilibrium())
	{
		const auto& L = cm.conservedMoietyMatrix();
		const unsigned int nMoieties = cm.numMoieties();
		// todo refactor this
		for (unsigned int point = 0; point < _disc.nPoints; ++point)
		{
			const int pointOffset = idxr.offsetC() + point * idxr.strideColNode();

			for (unsigned int m = 0; m < nMoieties; ++m)
			{
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
					_globalJacDisc.coeffRef(pointOffset + m, pointOffset + comp) += alpha * L(m, comp);
			}
			// No mass-matrix entries for equilibrium rows
		}
	}
	else
	{
		// Add the remaining bulk time derivatives to global Jacobian
		linalg::BandedEigenSparseRowIterator jac(_globalJacDisc, idxr.offsetC());
		for (int i = jac.row(); i < idxr.offsetCp(); ++i, ++jac)
			jac[0] += alpha;
	}
}

/**
 * @brief Adds Jacobian @f$ \frac{\partial F}{\partial \dot{y}} @f$ to bead rows of system Jacobian
 * @details Actually adds @f$ \alpha \frac{\partial F}{\partial \dot{y}} @f$, which is useful
 *          for constructing the linear system in BDF time discretization.
 * @param [in,out] jac On entry, RowIterator of the particle block pointing to the beginning of a bead shell;
 *                     on exit, the iterator points to the end of the bead shell
 * @param [in] idxr Indexer
 * @param [in] alpha Value of \f$ \alpha \f$ (arises from BDF time discretization)
 * @param [in] parType Index of the particle type
 */
template <typename ConvDispOperator>
void ColumnModel1D<ConvDispOperator>::addTimeDerivativeToJacobianParticleShell(linalg::BandedEigenSparseRowIterator& jac, const Indexer& idxr, double alpha, unsigned int parType)
{
	parts::cell::addTimeDerivativeToJacobianParticleShell<linalg::BandedEigenSparseRowIterator, true>(jac, alpha, static_cast<double>(_particles[parType]->getPorosity()), _disc.nComp, _disc.nBound + _disc.nComp * parType,
		_particles[parType]->getPoreAccessFactor(), _disc.strideBound[parType], _disc.boundOffset + _disc.nComp * parType, _binding[parType]->reactionQuasiStationarity());
}

}  // namespace model

}  // namespace cadet
