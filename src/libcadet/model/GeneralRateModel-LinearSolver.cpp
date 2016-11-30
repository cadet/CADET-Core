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

#include "model/GeneralRateModel.hpp"
#include "linalg/DenseMatrix.hpp"
#include "linalg/BandMatrix.hpp"
#include "AdUtils.hpp"

#include <algorithm>
#include <functional>

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include "ParallelSupport.hpp"

#ifdef CADET_PARALLELIZE
	#include <tbb/tbb.h>

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
 * @param [in] timeFactor Factor which is premultiplied to the time derivatives originating from time transformation
 * @param [in] alpha Value of \f$ \alpha \f$ (arises from BDF time discretization)
 * @param [in] outerTol Error tolerance for the solution of the linear system from outer Newton iteration
 * @param [in,out] rhs On entry the right hand side of the linear equation system, on exit the solution
 * @param [in] weight Vector with error weights
 * @param [in] y Pointer to global state vector at which the Jacobian is evaluated
 * @param [in] yDot Pointer to global time derivative state vector at which the Jacobian is evaluated
 * @param [in] res Pointer to global residual vector at the point @p y, @p yDot
 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
 */
int GeneralRateModel::linearSolve(double t, double timeFactor, double alpha, double outerTol, double* const rhs, double const* const weight,
	double const* const y, double const* const yDot, double const* const res)
{
	BENCH_SCOPE(_timerLinearSolve);

	Indexer idxr(_disc);

	// ==== Step 1: Factorize diagonal Jacobian blocks

	// Factorize partial Jacobians only if required

#ifdef CADET_PARALLELIZE
	tbb::flow::graph g;
#endif

#ifdef CADET_PARALLELIZE
	node_t A(g, [&](msg_t)
	{
#endif
		if (_factorizeJacobian)
		{
			// Assemble and factorize discretized system Jacobians
			// Threads that are done with the bulk column blocks can proceed to the particle blocks

#ifdef CADET_PARALLELIZE
			tbb::parallel_for(size_t(0), size_t(_disc.nComp), [&](size_t comp)
#else
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
#endif
			{
				// Assemble
				assembleDiscretizedJacobianColumnBlock(comp, alpha, idxr, timeFactor);

				// Factorize
				const bool result = _jacCdisc[comp].factorize();
				if (cadet_unlikely(!result))
				{
					LOG(Error) << "Factorize() failed for comp " << comp;
				}
			} CADET_PARFOR_END;
		}
	CADET_PARNODE_END;


	// Process the particle blocks
#ifdef CADET_PARALLELIZE
	node_t B(g, [&](msg_t)
	{
#endif
		if (_factorizeJacobian)
		{
#ifdef CADET_PARALLELIZE
			tbb::parallel_for(size_t(0), size_t(_disc.nCol), [&](size_t pblk)
#else
			for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
#endif
			{
				// Assemble
				assembleDiscretizedJacobianParticleBlock(pblk, alpha, idxr, timeFactor);

				// Factorize
				const bool result = _jacPdisc[pblk].factorize();
				if (cadet_unlikely(!result))
				{
					{
						LOG(Error) << "Factorize() failed for par block " << pblk;
					}
				}
			} CADET_PARFOR_END;
		}
	CADET_PARNODE_END;

	// ====== Step 1.5: Solve J c_uo = b_uo - A * c_in = b_uo - A*b_in

	// rhs is passed twice but due to the values in jacA the writes happen to a different area of the rhs than the reads.
#ifdef CADET_PARALLELIZE
	node_t C(g, [&](msg_t) 
	{
#endif
		// Do not factorize again at next call without changed Jacobians
		_factorizeJacobian = false;

		_jacInlet.multiplySubtract(rhs, rhs + idxr.offsetC());
	CADET_PARNODE_END;

	// ==== Step 2: Solve diagonal Jacobian blocks J_i to get y_i = J_i^{-1} b_i
	// The result is stored in rhs (in-place solution)


	// Threads that are done with solving the bulk column blocks can proceed
	// to solving the particle blocks
#ifdef CADET_PARALLELIZE
	node_t D(g, [&](msg_t)
	{
#endif
#ifdef CADET_PARALLELIZE
		tbb::parallel_for(size_t(0), size_t(_disc.nComp), [&](size_t comp)
#else
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
#endif
		{
			const bool result = _jacCdisc[comp].solve(rhs + comp * idxr.strideColComp() + idxr.offsetC());
			if (cadet_unlikely(!result))
			{
				LOG(Error) << "Solve() failed for comp " << comp;
			}
		} CADET_PARFOR_END;
	CADET_PARNODE_END;

#ifdef CADET_PARALLELIZE
	node_t E(g, [&](msg_t)
	{
#endif
#ifdef CADET_PARALLELIZE
		tbb::parallel_for(size_t(0), size_t(_disc.nCol), [&](size_t pblk)
#else
		for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
#endif
		{
			const bool result = _jacPdisc[pblk].solve(rhs + idxr.offsetCp(pblk));
			if (cadet_unlikely(!result))
			{
				LOG(Error) << "Solve() failed for par block " << pblk;
			}
		} CADET_PARFOR_END;
	CADET_PARNODE_END;

	// Solve last row of L with backwards substitution: y_f = b_f - \sum_{i=0}^{N_z} J_{f,i} y_i
	// Note that we cannot easily parallelize this loop since the results of the sparse
	// matrix-vector multiplications are added in-place to rhs. We would need one copy of rhs
	// for each thread and later fuse them together (reduction statement).
#ifdef CADET_PARALLELIZE
	node_t F(g, [&](msg_t) 
	{
#endif
		_jacFC.multiplySubtract(rhs + idxr.offsetC(), rhs + idxr.offsetJf());

		for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
		{
			_jacFP[pblk].multiplySubtract(rhs + idxr.offsetCp(pblk), rhs + idxr.offsetJf());
		}

		// Now, rhs contains the full intermediate solution y = L^{-1} b

		// Initialize temporary storage by copying over the fluxes
		// Note that the rest of _tempState is zeroed out in schurComplementMatrixVector()
		std::copy(rhs + idxr.offsetJf(), rhs + numDofs(), _tempState + idxr.offsetJf());

		// ==== Step 3: Solve Schur-complement to get x_f = S^{-1} y_f
		// Column and particle parts remain unchanged.
		// The only thing to be done is the iterative (and approximate)
		// solution of the Schur complement system:
		//     S * x_f = y_f

		// Note that rhs is updated in-place with the solution of the Schur-complement
		// The temporary storage is only needed to hold the right hand side of the Schur-complement
		const double tolerance = std::sqrt(static_cast<double>(numDofs())) * outerTol * _schurSafety;

		BENCH_START(_timerGmres);
		const int gmresResult = _gmres.solve(tolerance, weight + idxr.offsetJf(), _tempState + idxr.offsetJf(), rhs + idxr.offsetJf());
		BENCH_STOP(_timerGmres);

			// Remove temporary results that are leftovers from schurComplementMatrixVector()
		std::fill(_tempState + idxr.offsetC(), _tempState + idxr.offsetJf(), 0.0);

		// At this point, rhs contains the intermediate solution [y_0, ..., y_{N_z}, x_f]

		// ==== Step 4: Solve U * x = y by backward substitution
		// The fluxes are already solved and remain unchanged

		// Compute tempState_0 = J_{0,f} * y_f
		_jacCF.multiplyAdd(rhs + idxr.offsetJf(), _tempState + idxr.offsetC());
	CADET_PARNODE_END;

	// Threads that are done with solving the bulk column blocks can proceed
	// to solving the particle blocks
#ifdef CADET_PARALLELIZE
	node_t G(g, [&](msg_t) 
	{
#endif
#ifdef CADET_PARALLELIZE
		tbb::parallel_for(size_t(0), size_t(_disc.nComp), [&](size_t comp)
#else
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
#endif
		{
			double* const localCol = _tempState + comp * idxr.strideColComp() + idxr.offsetC();
			double* const rhsCol = rhs + comp * idxr.strideColComp() + idxr.offsetC();

			// Apply J_0^{-1} to tempState_0
			const bool result = _jacCdisc[comp].solve(localCol);
			if (cadet_unlikely(!result))
			{
				LOG(Error) << "Solve() failed for comp " << comp;
			}

			// Compute rhs_0 = y_0 - J_0^{-1} * J_{0,f} * y_f = y_0 - tempState_0
			for (unsigned int i = 0; i < _disc.nCol; ++i)
				rhsCol[i] -= localCol[i];
		} CADET_PARFOR_END;
	CADET_PARNODE_END;

#ifdef CADET_PARALLELIZE
	node_t H(g, [&](msg_t) 
	{
#endif
#ifdef CADET_PARALLELIZE
		tbb::parallel_for(size_t(0), size_t(_disc.nCol), [&](size_t pblk)
#else
		for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
#endif
		{
			double* const localPar = _tempState + idxr.offsetCp(pblk);
			double* const rhsPar = rhs + idxr.offsetCp(pblk);

			// Compute tempState_i = J_{i,f} * y_f
			_jacPF[pblk].multiplyAdd(rhs + idxr.offsetJf(), localPar);
			// Apply J_i^{-1} to tempState_i
			const bool result = _jacPdisc[pblk].solve(localPar);
			if (cadet_unlikely(!result))
			{
				LOG(Error) << "Solve() failed for par block " << pblk;
			}

			// Compute rhs_i = y_i - J_i^{-1} * J_{i,f} * y_f = y_i - tempState_i
			for (int i = 0; i < idxr.strideParBlock(); ++i)
				rhsPar[i] -= localPar[i];
		} CADET_PARFOR_END;
	CADET_PARNODE_END;

#ifdef CADET_PARALLELIZE
	// Create TBB dependency graph
	make_edge(A, C);
	make_edge(B, C);

	make_edge(C, D);
	make_edge(C, E);
	make_edge(D, F);
	make_edge(E, F);
	make_edge(F, G);
	make_edge(F, H);

	// Start the graph running
	A.try_put(tbb::flow::continue_msg());
	B.try_put(tbb::flow::continue_msg());

	// Wait for results
	g.wait_for_all();
#endif

	// The full solution is now stored in rhs
	return 0;
}

/**
 * @brief Performs the matrix-vector product @f$ z = Sx @f$ with the Schur-complement @f$ S @f$ from the Jacobian
 * @details The Schur-complement @f$ S @f$ is given by
 *          @f[ \begin{align}
				S &= J_f - J_{f,0} \, J_0^{-1} \, J_{0,f} - \sum_{p=1}^{N_z}{J_{f,p} \, J_p^{-1} \, J_{p,f}} \\
				  &= I - \sum_{p=0}^{N_z}{J_{f,p} \, J_p^{-1} \, J_{p,f}},
			\end{align} @f]
 *          where @f$ J_f = I @f$ is the identity matrix and the off-diagonal blocks @f$ J_{i,f} @f$
 *          and @f$ J_{f,i} @f$ for @f$ i = 0, \dots, N_{z} @f$ are sparse.
 *
 *          The matrix-vector multiplication is executed in parallel as follows:
 *              -# Compute @f$ J_{f,i} \, J_i^{-1} \, J_{i,f} @f$ independently (in parallel with respect to index @f$ i @f$)
 *              -# Subtract the result from @f$ z @f$
 *
 * @param [in] x Vector @f$ x @f$ the matrix @f$ S @f$ is multiplied with
 * @param [out] z Result of the matrix-vector multiplication
 * @return @c 0 if successful, any other value in case of failure
 */
int GeneralRateModel::schurComplementMatrixVector(double const* x, double* z) const
{
	BENCH_SCOPE(_timerMatVec);

	// Copy x over to result z, which corresponds to the application of the identity matrix
	std::copy(x, x + _disc.nCol * _disc.nComp, z);

	Indexer idxr(_disc);
	std::fill(_tempState + idxr.offsetC(), _tempState + idxr.offsetJf(), 0.0);

#ifdef CADET_PARALLELIZE
	tbb::flow::graph g;
#endif

	// Solve bulk column blocks first

	// Apply J_{0,f}
	_jacCF.multiplyAdd(x, _tempState + idxr.offsetC());

#ifdef CADET_PARALLELIZE
	node_t A(g, [&](msg_t) 
	{
#endif

#ifdef CADET_PARALLELIZE
		tbb::parallel_for(size_t(0), size_t(_disc.nComp), [&](size_t comp)
#else
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
#endif
		{
			// Apply J_0^{-1} of each component
			const bool result = _jacCdisc[comp].solve(_tempState + comp * idxr.strideColComp() + idxr.offsetC());
			if (cadet_unlikely(!result))
			{
				LOG(Error) << "Solve() failed for comp " << comp;
			}
		} CADET_PARFOR_END;
	CADET_PARNODE_END;

#ifdef CADET_PARALLELIZE
	node_t B(g, [&](msg_t)
	{
#endif
		// Handle particle blocks
#ifdef CADET_PARALLELIZE
		tbb::parallel_for(size_t(0), size_t(_disc.nCol), [&](size_t pblk)
#else
		for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
#endif
		{
			// Get this thread's temporary memory block
			double* const tmp = _tempState + idxr.offsetCp(pblk);

			// Apply J_{i,f}
			_jacPF[pblk].multiplyAdd(x, tmp);
			// Apply J_{i}^{-1}
			const bool result = _jacPdisc[pblk].solve(tmp);
			if (cadet_unlikely(!result))
			{
				LOG(Error) << "Solve() failed for par block " << pblk;
			}
		} CADET_PARFOR_END;
	CADET_PARNODE_END;

#ifdef CADET_PARALLELIZE
	node_t C(g, [&](msg_t) 
	{
#endif
		// Apply J_{f,0} and subtract results from z
		_jacFC.multiplySubtract(_tempState + idxr.offsetC(), z);

		for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
		{
			// Apply J_{f,i} and subtract results from z
			_jacFP[pblk].multiplySubtract(_tempState + idxr.offsetCp() + idxr.strideParBlock() * pblk, z);
		}
	CADET_PARNODE_END;

#ifdef CADET_PARALLELIZE
	make_edge(A, C);
	make_edge(B, C);

	// Start the graph running
	A.try_put(tbb::flow::continue_msg());
	B.try_put(tbb::flow::continue_msg());

	// Wait for results
	g.wait_for_all();
#endif

	return 0;
}

/**
 * @brief Assembles the column void Jacobian block @f$ J_0 @f$ of the time-discretized equations
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
 * @param [in] comp Index of the current component
 * @param [in] alpha Value of \f$ \alpha \f$ (arises from BDF time discretization)
 * @param [in] idxr Indexer
 * @param [in] timeFactor Factor which is premultiplied to the time derivatives originating from time transformation
 */
void GeneralRateModel::assembleDiscretizedJacobianColumnBlock(unsigned int comp, double alpha, const Indexer& idxr, double timeFactor)
{
	linalg::FactorizableBandMatrix& fbm = _jacCdisc[comp];
	const linalg::BandMatrix& bm = _jacC[comp];

	// Copy normal matrix over to factorizable matrix
	fbm.copyOver(bm);

	// Add time derivatives
	addTimeDerivativeToJacobianColumnBlock(fbm, idxr, alpha, timeFactor);
}

/**
 * @brief Adds the derivatives with respect to @f$ \dot{y} @f$ of @f$ F(t, y, \dot{y}) @f$ to the Jacobian blocks
 * @details Given a FactorizableBandMtrix @p fbm, this functions computes 
 *          @f[ \begin{align*} \text{fbm} = \text{fbm} + \alpha \frac{\partial F}{\partial \dot{y}}. \end{align*} @f]
 *          The factor @f$ \alpha @f$ is useful when constructing the linear system in the time integration process.
 * @param [in,out] fbm FactorizableBandMatrix to which the derivatives with respect to @f$ \dot{y} @f$ are added
 * @param [in] idxr Indexer
 * @param [in] alpha Factor in front of @f$ \frac{\partial F}{\partial \dot{y}} @f$
 * @param [in] timeFactor Factor which is premultiplied to the time derivatives originating from time transformation
 */
void GeneralRateModel::addTimeDerivativeToJacobianColumnBlock(linalg::FactorizableBandMatrix& fbm, const Indexer& idxr, double alpha, double timeFactor)
{
	alpha *= timeFactor;

	linalg::FactorizableBandMatrix::RowIterator jac = fbm.row(0);
	for (unsigned int i = 0; i < _disc.nCol; ++i, ++jac)
	{
		// Add time derivative to main diagonal
		jac[0] += alpha;
	}
}

/**
 * @brief Assembles a particle Jacobian block @f$ J_i @f$ (@f$ i > 0 @f$) of the time-discretized equations
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
 * @param [in] pblk Index of the particle block
 * @param [in] alpha Value of \f$ \alpha \f$ (arises from BDF time discretization)
 * @param [in] idxr Indexer
 * @param [in] timeFactor Factor which is premultiplied to the time derivatives originating from time transformation
 */
void GeneralRateModel::assembleDiscretizedJacobianParticleBlock(unsigned int pblk, double alpha, const Indexer& idxr, double timeFactor)
{
	linalg::FactorizableBandMatrix& fbm = _jacPdisc[pblk];
	const linalg::BandMatrix& bm = _jacP[pblk];

	// Copy normal matrix over to factorizable matrix
	fbm.copyOver(bm);

	// Add time derivatives to particle shells
	const double invBetaP = 1.0 / static_cast<double>(_parPorosity) - 1.0;
	linalg::FactorizableBandMatrix::RowIterator jac = fbm.row(0);
	for (unsigned int j = 0; j < _disc.nPar; ++j)
	{
		// Mobile phase
		addMobilePhaseTimeDerivativeToJacobianParticleBlock(jac, idxr, alpha, invBetaP, timeFactor);

		// Stationary phase
		_binding->jacobianAddDiscretized(alpha * timeFactor, jac);

		// Advance pointers over all bound states
		jac += idxr.strideParBound();
	}
}

/**
 * @brief Adds Jacobian @f$ \frac{\partial F}{\partial \dot{y}} @f$ to bead mobile phase rows of system Jacobian
 * @details Actually adds @f$ \alpha \frac{\partial F}{\partial \dot{y}} @f$, which is useful
 *          for constructing the linear system in BDF time discretization.
 * @param [in,out] jac On entry, RowIterator of the particle block pointing to the beginning of a bead shell;
 *                     on exit, the iterator points to the end of the mobile phase
 * @param [in] idxr Indexer
 * @param [in] alpha Value of \f$ \alpha \f$ (arises from BDF time discretization)
 * @param [in] invBetaP Inverse porosity term @f$\frac{1}{\beta_p}@f$
 * @param [in] timeFactor Factor which is premultiplied to the time derivatives originating from time transformation
 */
void GeneralRateModel::addMobilePhaseTimeDerivativeToJacobianParticleBlock(linalg::FactorizableBandMatrix::RowIterator& jac, const Indexer& idxr, double alpha, double invBetaP, double timeFactor)
{
	// Compute total factor
	alpha *= timeFactor;

	// Mobile phase
	for (int comp = 0; comp < static_cast<int>(_disc.nComp); ++comp, ++jac)
	{
		// Add derviative with respect to dc_p / dt to Jacobian
		jac[0] += alpha;

		// Add derivative with respect to dq / dt to Jacobian
		for (int i = 0; i < static_cast<int>(_disc.nBound[comp]); ++i)
		{
			// Index explanation:
			//   -comp -> go back to beginning of liquid phase
			//   + strideParLiquid() skip to solid phase
			//   + offsetBoundComp() jump to component (skips all bound states of previous components)
			//   + i go to current bound state
			jac[idxr.strideParLiquid() - comp + idxr.offsetBoundComp(comp) + i] += alpha * invBetaP;
		}
	}
}


}  // namespace model

}  // namespace cadet