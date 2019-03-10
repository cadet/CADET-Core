// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//
//  Copyright © 2008-2019: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "model/GeneralRateModel.hpp"
#include "model/BindingModel.hpp"
#include "linalg/DenseMatrix.hpp"
#include "linalg/BandMatrix.hpp"
#include "ParamReaderHelper.hpp"
#include "AdUtils.hpp"
#include "model/parts/BindingConsistentInit.hpp"
#include "SimulationTypes.hpp"

#include <algorithm>
#include <functional>

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include "ParallelSupport.hpp"
#ifdef CADET_PARALLELIZE
	#include <tbb/tbb.h>
#endif

namespace cadet
{

namespace model
{

void GeneralRateModel::applyInitialCondition(const SimulationState& simState) const
{
	Indexer idxr(_disc);

	// Check whether full state vector is available as initial condition
	if (!_initState.empty())
	{
		std::fill(simState.vecStateY, simState.vecStateY + idxr.offsetC(), 0.0);
		std::copy(_initState.data(), _initState.data() + numPureDofs(), simState.vecStateY + idxr.offsetC());

		if (!_initStateDot.empty())
		{
			std::fill(simState.vecStateYdot, simState.vecStateYdot + idxr.offsetC(), 0.0);
			std::copy(_initStateDot.data(), _initStateDot.data() + numPureDofs(), simState.vecStateYdot + idxr.offsetC());
		}
		else
			std::fill(simState.vecStateYdot, simState.vecStateYdot + numDofs(), 0.0);

		return;
	}

	double* const stateYbulk = simState.vecStateY + idxr.offsetC();

	// Loop over column cells
	for (unsigned int col = 0; col < _disc.nCol; ++col)
	{
		// Loop over components in cell
		for (unsigned comp = 0; comp < _disc.nComp; ++comp)
			stateYbulk[col * idxr.strideColCell() + comp * idxr.strideColComp()] = static_cast<double>(_initC[comp]);
	}

	// Loop over particles
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		for (unsigned int col = 0; col < _disc.nCol; ++col)
		{
			const unsigned int offset = idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{col});

			// Loop over particle cells
			for (unsigned int shell = 0; shell < _disc.nParCell[type]; ++shell)
			{
				const unsigned int shellOffset = offset + shell * idxr.strideParShell(type);

				// Initialize c_p
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
					simState.vecStateY[shellOffset + comp] = static_cast<double>(_initCp[comp + _disc.nComp * type]);

				// Initialize q
				active const* const iq = _initQ.data() + _disc.nBoundBeforeType[_singleBinding ? 0 : type];
				for (unsigned int bnd = 0; bnd < _disc.strideBound[_singleBinding ? 0 : type]; ++bnd)
					simState.vecStateY[shellOffset + idxr.strideParLiquid() + bnd] = static_cast<double>(iq[bnd]);
			}
		}
	}
}

void GeneralRateModel::readInitialCondition(IParameterProvider& paramProvider)
{
	_initState.clear();
	_initStateDot.clear();

	// Check if INIT_STATE is present
	if (paramProvider.exists("INIT_STATE"))
	{
		const std::vector<double> initState = paramProvider.getDoubleArray("INIT_STATE");
		_initState = std::vector<double>(initState.begin(), initState.begin() + numPureDofs());

		// Check if INIT_STATE contains the full state and its time derivative
		if (initState.size() >= 2 * numPureDofs())
			_initStateDot = std::vector<double>(initState.begin() + numPureDofs(), initState.begin() + 2 * numPureDofs());
		return;
	}

	const std::vector<double> initC = paramProvider.getDoubleArray("INIT_C");
	std::vector<double> initQ;

	if (paramProvider.exists("INIT_Q"))
		initQ = paramProvider.getDoubleArray("INIT_Q");

	if (initC.size() < _disc.nComp)
		throw InvalidParameterException("INIT_C does not contain enough values for all components");

	if ((_disc.strideBound[_disc.nParType] > 0) && (((initQ.size() < _disc.strideBound[_disc.nParType]) && !_singleBinding) || (initQ.size() < _disc.strideBound[0] && _singleBinding)))
		throw InvalidParameterException("INIT_Q does not contain enough values for all bound states");

	// Check if INIT_CP is present
	if (paramProvider.exists("INIT_CP"))
	{
		const std::vector<double> initCp = paramProvider.getDoubleArray("INIT_CP");

		if (initCp.size() < _disc.nComp * _disc.nParType)
			throw InvalidParameterException("INIT_CP does not contain enough values for all components");

		ad::copyToAd(initCp.data(), _initCp.data(), _disc.nComp * _disc.nParType);
	}
	else
	{
		for (unsigned int type = 0; type < _disc.nParType; ++type)
			ad::copyToAd(initC.data(), _initCp.data() + _disc.nComp * type, _disc.nComp);
	}

	ad::copyToAd(initC.data(), _initC.data(), _disc.nComp);
	if (!initQ.empty())
		ad::copyToAd(initQ.data(), _initQ.data(), _disc.strideBound[_singleBinding ? 0 : _disc.nParType]);
}

/**
 * @brief Computes consistent initial values (state variables without their time derivatives)
 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] the initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ have
 *          to be consistent. This functions updates the initial state \f$ y_0 \f$ and overwrites the time
 *          derivative \f$ \dot{y}_0 \f$ such that they are consistent.
 *          
 *          The process works in two steps:
 *          <ol>
 *              <li>Solve all algebraic equations in the model (e.g., quasi-stationary isotherms, reaction equilibria).
 *                 Once all @f$ c_i @f$, @f$ c_{p,i} @f$, and @f$ q_i^{(j)} @f$ have been computed, solve for the
 *                 fluxes @f$ j_{f,i} @f$ (only linear equations).</li>
 *              <li>Compute the time derivatives of the state @f$ \dot{y} @f$ such that the residual is 0.
 *                 However, because of the algebraic equations, we need additional conditions to fully determine
 *                 @f$ \dot{y}@f$. By differentiating the algebraic equations with respect to time, we get the
 *                 missing linear equations (recall that the state vector @f$ y @f$ is fixed). The resulting system
 *                 has a similar structure as the system Jacobian.
 *                 @f[ \begin{align}
 *                  \left[\begin{array}{c|ccc|c}
 *                     \dot{J}_0  &         &        &           &   \\
 *                     \hline
 *                                & \dot{J}_1     &        &           &   \\
 *                                &         & \ddots &           &   \\
 *                                &         &        & \dot{J}_{N_z}   &   \\
 *                     \hline
 *                        J_{f,0} & J_{f,1} & \dots & J_{f,N_z} & I
 *                  \end{array}\right],
 *                 \end{align} @f]
 *                 where @f$ \dot{J}_i @f$ denotes the Jacobian with respect to @f$ \dot{y}@f$. Note that the
 *                 @f$ J_{i,f} @f$ matrices in the right column are missing.
 *                 
 *     The right hand side of the linear system is given by the negative residual without contribution 
 *     of @f$ \dot{y} @f$ for differential equations and 0 for algebraic equations 
 *     (@f$ -\frac{\partial F}{\partial t}@f$, to be more precise).
 *     
 *     The linear system is solved by backsubstitution. First, the diagonal blocks are solved in parallel.
 *     Then, the equations for the fluxes @f$ j_f @f$ are solved by substituting in the solution of the
 *     diagonal blocks.</li>
 *          </ol>
 *     This function performs step 1. See consistentInitialTimeDerivative() for step 2.
 *     
 * 	   This function is to be used with consistentInitialTimeDerivative(). Do not mix normal and lean
 *     consistent initialization!
 *     
 * @param [in] simTime Simulation time information (time point, section index, pre-factor of time derivatives)
 * @param [in,out] vecStateY State vector with initial values that are to be updated for consistency
 * @param [in,out] adJac Jacobian information for AD (AD vectors for residual and state, direction offset)
 * @param [in] errorTol Error tolerance for algebraic equations
 * @todo Decrease amount of allocated memory by partially using temporary vectors (state and Schur complement)
 */
void GeneralRateModel::consistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol)
{
	BENCH_SCOPE(_timerConsistentInit);

	Indexer idxr(_disc);

	// Step 1: Solve algebraic equations

	// Step 1a: Compute quasi-stationary binding model state
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		if (_binding[type]->hasAlgebraicEquations())
		{
			// TODO: Check memory consumption and offsets
			// Round up
			const unsigned int requiredMem = (_binding[type]->workspaceSize(_disc.nComp, _disc.strideBound[type], _disc.nBound + type * _disc.nComp) + sizeof(double) - 1) / sizeof(double);

			ad::BandedJacobianExtractor jacExtractor(_jacP[type * _disc.nCol].lowerBandwidth(), _jacP[type * _disc.nCol].lowerBandwidth(), _jacP[type * _disc.nCol].upperBandwidth());

			//Problem capturing variables here
#ifdef CADET_PARALLELIZE
			BENCH_SCOPE(_timerConsistentInitPar);
			tbb::parallel_for(size_t(0), size_t(_disc.nCol), [&](size_t pblk)
#else
			for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
#endif
			{
				// Reuse memory of band matrix for dense matrix
				linalg::DenseMatrixView jacobianMatrix(_jacPdisc[type * _disc.nCol + pblk].data(), _jacPdisc[type * _disc.nCol + pblk].pivot(), _disc.strideBound[type], _disc.strideBound[type]);

				// Midpoint of current column cell (z coordinate) - needed in externally dependent adsorption kinetic
				const double z = (0.5 + static_cast<double>(pblk)) / static_cast<double>(_disc.nCol);

				// This loop cannot be run in parallel without creating a Jacobian matrix for each thread which would increase memory usage
				const int localOffsetToParticle = idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{static_cast<unsigned int>(pblk)});
				for(size_t shell = 0; shell < size_t(_disc.nParCell[type]); ++shell)
				{
					const int localOffsetInParticle = static_cast<int>(shell) * idxr.strideParShell(type) + idxr.strideParLiquid();

					// Get pointer to q variables in a shell of particle pblk
					double* const qShell = vecStateY + localOffsetToParticle + localOffsetInParticle;
					active* const localAdRes = adJac.adRes ? adJac.adRes + localOffsetToParticle : nullptr;
					active* const localAdY = adJac.adY ? adJac.adY + localOffsetToParticle : nullptr;

					// We are essentially creating a 2d vector of blocks out of a linear strip of memory
					const unsigned int offset = _bindingWorkspaceOffset[type] + requiredMem * pblk;

					// Solve algebraic variables
					_binding[type]->consistentInitialState(simTime.t, simTime.secIdx, ColumnPosition{z, 0.0, static_cast<double>(_parCenterRadius[_disc.nParCellsBeforeType[type] + shell]) / static_cast<double>(_parRadius[type])}, qShell, qShell - idxr.strideParLiquid(), errorTol,
						localAdRes, localAdY, localOffsetInParticle, adJac.adDirOffset, jacExtractor, _tempState + offset, jacobianMatrix);
				}
			} CADET_PARFOR_END;
		}
	}

	// Step 1b: Compute fluxes j_f

	// Reset j_f to 0.0
	double* const jf = vecStateY + idxr.offsetJf();
	std::fill(jf, jf + _disc.nComp * _disc.nCol * _disc.nParType, 0.0);

	solveForFluxes(vecStateY, idxr);
}

/**
 * @brief Computes consistent initial time derivatives
 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] the initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ have
 *          to be consistent. This functions updates the initial state \f$ y_0 \f$ and overwrites the time
 *          derivative \f$ \dot{y}_0 \f$ such that they are consistent.
 *          
 *          The process works in two steps:
 *          <ol>
 *              <li>Solve all algebraic equations in the model (e.g., quasi-stationary isotherms, reaction equilibria).
 *                 Once all @f$ c_i @f$, @f$ c_{p,i} @f$, and @f$ q_i^{(j)} @f$ have been computed, solve for the
 *                 fluxes @f$ j_{f,i} @f$ (only linear equations).</li>
 *              <li>Compute the time derivatives of the state @f$ \dot{y} @f$ such that the residual is 0.
 *                 However, because of the algebraic equations, we need additional conditions to fully determine
 *                 @f$ \dot{y}@f$. By differentiating the algebraic equations with respect to time, we get the
 *                 missing linear equations (recall that the state vector @f$ y @f$ is fixed). The resulting system
 *                 has a similar structure as the system Jacobian.
 *                 @f[ \begin{align}
 *                  \left[\begin{array}{c|ccc|c}
 *                     \dot{J}_0  &         &        &           &   \\
 *                     \hline
 *                                & \dot{J}_1     &        &           &   \\
 *                                &         & \ddots &           &   \\
 *                                &         &        & \dot{J}_{N_z}   &   \\
 *                     \hline
 *                        J_{f,0} & J_{f,1} & \dots & J_{f,N_z} & I
 *                  \end{array}\right],
 *                 \end{align} @f]
 *                 where @f$ \dot{J}_i @f$ denotes the Jacobian with respect to @f$ \dot{y}@f$. Note that the
 *                 @f$ J_{i,f} @f$ matrices in the right column are missing.
 *                 
 *     The right hand side of the linear system is given by the negative residual without contribution 
 *     of @f$ \dot{y} @f$ for differential equations and 0 for algebraic equations 
 *     (@f$ -\frac{\partial F}{\partial t}@f$, to be more precise).
 *     
 *     The linear system is solved by backsubstitution. First, the diagonal blocks are solved in parallel.
 *     Then, the equations for the fluxes @f$ j_f @f$ are solved by substituting in the solution of the
 *     diagonal blocks.</li>
 *          </ol>
 *     This function performs step 2. See consistentInitialState() for step 1.
 *     
 * 	   This function is to be used with consistentInitialState(). Do not mix normal and lean
 *     consistent initialization!
 *     
 * @param [in] simTime Simulation time information (time point, section index, pre-factor of time derivatives)
 * @param [in] vecStateY Consistently initialized state vector
 * @param [in,out] vecStateYdot On entry, residual without taking time derivatives into account. On exit, consistent state time derivatives.
 */
void GeneralRateModel::consistentInitialTimeDerivative(const SimulationTime& simTime, double const* vecStateY, double* const vecStateYdot)
{
	BENCH_SCOPE(_timerConsistentInit);

	Indexer idxr(_disc);

	// Step 2: Compute the correct time derivative of the state vector

	// Step 2a: Assemble, factorize, and solve diagonal blocks of linear system

	// Note that the residual has not been negated, yet. We will do that now.
	for (unsigned int i = 0; i < numDofs(); ++i)
		vecStateYdot[i] = -vecStateYdot[i];

	// Handle bulk column block
	_convDispOp.solveTimeDerivativeSystem(simTime, vecStateYdot + idxr.offsetC());

	// Process the particle blocks
#ifdef CADET_PARALLELIZE
	BENCH_START(_timerConsistentInitPar);
	tbb::parallel_for(size_t(0), size_t(_disc.nCol * _disc.nParType), [&](size_t pblk)
#else
	for (unsigned int pblk = 0; pblk < _disc.nCol * _disc.nParType; ++pblk)
#endif
	{
		const unsigned int type = pblk / _disc.nCol;
		const unsigned int par = pblk % _disc.nCol;

		// Midpoint of current column cell (z coordinate) - needed in externally dependent adsorption kinetic
		const double z = (0.5 + static_cast<double>(par)) / static_cast<double>(_disc.nCol);

		// Assemble
		linalg::FactorizableBandMatrix& fbm = _jacPdisc[pblk];
		fbm.setAll(0.0);

		linalg::FactorizableBandMatrix::RowIterator jac = fbm.row(0);
		for (unsigned int j = 0; j < _disc.nParCell[type]; ++j)
		{
			// Mobile phase (advances jac accordingly)
			addMobilePhaseTimeDerivativeToJacobianParticleBlock(jac, idxr, 1.0, simTime.timeFactor, type);

			// Stationary phase
			// Populate matrix with time derivative Jacobian first
			_binding[type]->jacobianAddDiscretized(simTime.timeFactor, jac);

			// Overwrite rows corresponding to algebraic equations with the Jacobian and set right hand side to 0
			if (_binding[type]->hasAlgebraicEquations())
			{
				const unsigned int requiredMem = (_binding[type]->workspaceSize(_disc.nComp, _disc.strideBound[type], _disc.nBound + type * _disc.nComp) + sizeof(double) - 1) / sizeof(double);

				parts::BindingConsistentInitializer::consistentInitialTimeDerivative(_binding[type], simTime, jac,
					_jacP[pblk].row(j * static_cast<unsigned int>(idxr.strideParShell(type)) + static_cast<unsigned int>(idxr.strideParLiquid())),
					vecStateYdot + idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{par}) + static_cast<int>(j) * idxr.strideParShell(type) + idxr.strideParLiquid(),
					ColumnPosition{z, 0.0, static_cast<double>(_parCenterRadius[_disc.nParCellsBeforeType[type] + j]) / static_cast<double>(_parRadius[type])},
					_tempState + _bindingWorkspaceOffset[type] + requiredMem * par);
			}

			// Advance pointers over all bound states
			jac += idxr.strideParBound(type);
		}

		// Precondition
		double* const scaleFactors = _tempState + idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{par});
		fbm.rowScaleFactors(scaleFactors);
		fbm.scaleRows(scaleFactors);

		// Factorize
		const bool result = fbm.factorize();
		if (!result)
		{
			LOG(Error) << "Factorize() failed for par block " << pblk << " (type " << type << " col " << par << ")\n" << fbm;
		}

		// Solve
		const bool result2 = fbm.solve(scaleFactors, vecStateYdot + idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{par}));
		if (!result2)
		{
			LOG(Error) << "Solve() failed for par block " << pblk << " (type " << type << " col " << par << ")";
		}
	} CADET_PARFOR_END;

#ifdef CADET_PARALLELIZE
	BENCH_STOP(_timerConsistentInitPar);
#endif

	// Step 2b: Solve for fluxes j_f by backward substitution

	// Reset \dot{j}_f to 0.0
	double* const jfDot = vecStateYdot + idxr.offsetJf();
	std::fill(jfDot, jfDot + _disc.nComp * _disc.nCol * _disc.nParType, 0.0);

	solveForFluxes(vecStateYdot, idxr);
}



/**
 * @brief Computes approximately / partially consistent initial values (state variables without their time derivatives)
 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] the initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ have
 *          to be consistent. This functions updates the initial state \f$ y_0 \f$ and overwrites the time
 *          derivative \f$ \dot{y}_0 \f$ such that they are consistent.
 *          
 *          This function performs a relaxed consistent initialization: Only parts of the vectors are updated
 *          and, hence, consistency is not guaranteed. Since there is less work to do, it is probably faster than
 *          the standard process represented by consistentInitialState().
 *          
 *          The process works in two steps:
 *          <ol>
 *              <li>Keep state and time derivative vectors as they are (i.e., do not solve algebraic equations).
 *                 Only solve for the fluxes @f$ j_{f,i} @f$ (only linear equations).</li>
 *              <li>Compute the time derivatives of the state @f$ \dot{y} @f$ such that the residual is 0 in the column
 *                 bulk and flux blocks. The resulting equations are stated below:
 *                 @f[ \begin{align}
 *                  \left[\begin{array}{c|ccc|c}
 *                     \dot{J}_0  &         &        &           &   \\
 *                     \hline
 *                        J_{f,0} & J_{f,1} & \dots & J_{f,N_z} & I
 *                  \end{array}\right],
 *                 \end{align} @f]
 *                 where @f$ \dot{J}_0 @f$ denotes the bulk block Jacobian with respect to @f$ \dot{y}@f$.
 *                 
 *     The right hand side of the linear system is given by the negative residual without contribution 
 *     of @f$ \dot{y} @f$ for the bulk block and 0 for the flux block.
 *     
 *     The linear system is solved by backsubstitution. First, the bulk block is solved.
 *     Then, the equations for the fluxes @f$ j_f @f$ are solved by substituting in the solution of the
 *     bulk block and the unchanged particle block time derivative vectors.</li>
 *          </ol>
 *     This function performs step 1. See leanConsistentInitialTimeDerivative() for step 2.
 *     
 * 	   This function is to be used with leanConsistentInitialTimeDerivative(). Do not mix normal and lean
 *     consistent initialization!
 *     
 * @param [in] simTime Simulation time information (time point, section index, pre-factor of time derivatives)
 * @param [in,out] vecStateY State vector with initial values that are to be updated for consistency
 * @param [in,out] adJac Jacobian information for AD (AD vectors for residual and state, direction offset)
 * @param [in] errorTol Error tolerance for algebraic equations
 */
void GeneralRateModel::leanConsistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol)
{
	if ((_parDiffusion.size() > _disc.nComp * _disc.nParType) || (_parSurfDiffusion.size() > _disc.strideBound[_disc.nParType]))
		LOG(Warning) << "Lean consistent initialization is not appropriate for section-dependent pore and surface diffusion";

	BENCH_SCOPE(_timerConsistentInit);

	Indexer idxr(_disc);

	// Step 1: Compute fluxes j_f

	// Reset j_f to 0.0
	double* const jf = vecStateY + idxr.offsetJf();
	std::fill(jf, jf + _disc.nComp * _disc.nCol * _disc.nParType, 0.0);

	solveForFluxes(vecStateY, idxr);
}

/**
 * @brief Computes approximately / partially consistent initial time derivatives
 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] the initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ have
 *          to be consistent. This functions updates the initial state \f$ y_0 \f$ and overwrites the time
 *          derivative \f$ \dot{y}_0 \f$ such that they are consistent.
 *          
 *          This function performs a relaxed consistent initialization: Only parts of the vectors are updated
 *          and, hence, consistency is not guaranteed. Since there is less work to do, it is probably faster than
 *          the standard process represented by consistentInitialTimeDerivative().
 *          
 *          The process works in two steps:
 *          <ol>
 *              <li>Keep state and time derivative vectors as they are (i.e., do not solve algebraic equations).
 *                 Only solve for the fluxes @f$ j_{f,i} @f$ (only linear equations).</li>
 *              <li>Compute the time derivatives of the state @f$ \dot{y} @f$ such that the residual is 0 in the column
 *                 bulk and flux blocks. The resulting equations are stated below:
 *                 @f[ \begin{align}
 *                  \left[\begin{array}{c|ccc|c}
 *                     \dot{J}_0  &         &        &           &   \\
 *                     \hline
 *                        J_{f,0} & J_{f,1} & \dots & J_{f,N_z} & I
 *                  \end{array}\right],
 *                 \end{align} @f]
 *                 where @f$ \dot{J}_0 @f$ denotes the bulk block Jacobian with respect to @f$ \dot{y}@f$.
 *                 
 *     The right hand side of the linear system is given by the negative residual without contribution 
 *     of @f$ \dot{y} @f$ for the bulk block and 0 for the flux block.
 *     
 *     The linear system is solved by backsubstitution. First, the bulk block is solved.
 *     Then, the equations for the fluxes @f$ j_f @f$ are solved by substituting in the solution of the
 *     bulk block and the unchanged particle block time derivative vectors.</li>
 *          </ol>
 *     This function performs step 2. See leanConsistentInitialState() for step 1.
 *     
 * 	   This function is to be used with leanConsistentInitialState(). Do not mix normal and lean
 *     consistent initialization!
 *     
 * @param [in] t Current time point
 * @param [in] timeFactor Used for time transformation (pre factor of time derivatives) and to compute parameter derivatives with respect to section length
 * @param [in] vecStateY (Lean) consistently initialized state vector
 * @param [in,out] vecStateYdot On entry, inconsistent state time derivatives. On exit, partially consistent state time derivatives.
 * @param [in] res On entry, residual without taking time derivatives into account. The data is overwritten during execution of the function.
 */
void GeneralRateModel::leanConsistentInitialTimeDerivative(double t, double timeFactor, double const* const vecStateY, double* const vecStateYdot, double* const res)
{
	if ((_parDiffusion.size() > _disc.nComp * _disc.nParType) || (_parSurfDiffusion.size() > _disc.strideBound[_disc.nParType]))
		LOG(Warning) << "Lean consistent initialization is not appropriate for section-dependent pore and surface diffusion";

	BENCH_SCOPE(_timerConsistentInit);

	Indexer idxr(_disc);

	// Step 2: Compute the correct time derivative of the state vector

	// Step 2a: Assemble, factorize, and solve column bulk block of linear system

	// Note that the residual is not negated as required at this point. We will fix that later.

	double* const resSlice = res + idxr.offsetC();

	// Handle bulk block
	_convDispOp.solveTimeDerivativeSystem(SimulationTime{t, 0u, timeFactor}, resSlice);

	// Note that we have solved with the *positive* residual as right hand side
	// instead of the *negative* one. Fortunately, we are dealing with linear systems,
	// which means that we can just negate the solution.
	double* const yDotSlice = vecStateYdot + idxr.offsetC();
	for (unsigned int i = 0; i < _disc.nCol * _disc.nComp; ++i)
		yDotSlice[i] = -resSlice[i];

	// Step 2b: Solve for fluxes j_f by backward substitution

	// Reset \dot{j}_f to 0.0
	double* const jfDot = vecStateYdot + idxr.offsetJf();
	std::fill(jfDot, jfDot + _disc.nComp * _disc.nCol * _disc.nParType, 0.0);

	solveForFluxes(vecStateYdot, idxr);
}

void GeneralRateModel::initializeSensitivityStates(const std::vector<double*>& vecSensY) const
{
	Indexer idxr(_disc);
	for (unsigned int param = 0; param < vecSensY.size(); ++param)
	{
		double* const stateYbulk = vecSensY[param] + idxr.offsetC();

		// Loop over column cells
		for (unsigned int col = 0; col < _disc.nCol; ++col)
		{
			// Loop over components in cell
			for (unsigned comp = 0; comp < _disc.nComp; ++comp)
				stateYbulk[col * idxr.strideColCell() + comp * idxr.strideColComp()] = _initC[comp].getADValue(param);
		}

		// Loop over particles
		for (unsigned int type = 0; type < _disc.nParType; ++type)
		{
			for (unsigned int col = 0; col < _disc.nCol; ++col)
			{
				const unsigned int offset = idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{col});

				// Loop over particle cells
				for (unsigned int shell = 0; shell < _disc.nParCell[type]; ++shell)
				{
					const unsigned int shellOffset = offset + shell * idxr.strideParShell(type);
					double* const stateYparticle = vecSensY[param] + shellOffset;
					double* const stateYparticleSolid = stateYparticle + idxr.strideParLiquid();
					
					// Initialize c_p
					for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
						stateYparticle[comp] = _initCp[comp + type * _disc.nComp].getADValue(param);

					// Initialize q
					for (unsigned int bnd = 0; bnd < _disc.strideBound[type]; ++bnd)
						stateYparticleSolid[bnd] = _initQ[bnd + _disc.nBoundBeforeType[type]].getADValue(param);
				}
			}
		}
	}
}

/**
 * @brief Computes consistent initial values and time derivatives of sensitivity subsystems
 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] and initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$,
 *          the sensitivity system for a parameter @f$ p @f$ reads
 *          \f[ \frac{\partial F}{\partial y}(t, y, \dot{y}) s + \frac{\partial F}{\partial \dot{y}}(t, y, \dot{y}) \dot{s} + \frac{\partial F}{\partial p}(t, y, \dot{y}) = 0. \f]
 *          The initial values of this linear DAE, @f$ s_0 = \frac{\partial y_0}{\partial p} @f$ and @f$ \dot{s}_0 = \frac{\partial \dot{y}_0}{\partial p} @f$
 *          have to be consistent with the sensitivity DAE. This functions updates the initial sensitivity\f$ s_0 \f$ and overwrites the time
 *          derivative \f$ \dot{s}_0 \f$ such that they are consistent.
 *          
 *          The process follows closely the one of consistentInitialConditions() and, in fact, is a linearized version of it.
 *          This is necessary because the initial conditions of the sensitivity system \f$ s_0 \f$ and \f$ \dot{s}_0 \f$ are
 *          related to the initial conditions \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ of the original DAE by differentiating them
 *          with respect to @f$ p @f$: @f$ s_0 = \frac{\partial y_0}{\partial p} @f$ and @f$ \dot{s}_0 = \frac{\partial \dot{y}_0}{\partial p}. @f$
 *          <ol>
 *              <li>Solve all algebraic equations in the model (e.g., quasi-stationary isotherms, reaction equilibria).
 *                 Once all @f$ c_i @f$, @f$ c_{p,i} @f$, and @f$ q_i^{(j)} @f$ have been computed, solve for the
 *                 fluxes @f$ j_{f,i} @f$. Let @f$ \mathcal{I}_a @f$ be the index set of algebraic equations, then, at this point, we have
 *                 \f[ \left( \frac{\partial F}{\partial y}(t, y_0, \dot{y}_0) s + \frac{\partial F}{\partial p}(t, y_0, \dot{y}_0) \right)_{\mathcal{I}_a} = 0. \f]</li>
 *              <li>Compute the time derivatives of the sensitivity @f$ \dot{s} @f$ such that the differential equations hold.
 *                 However, because of the algebraic equations, we need additional conditions to fully determine
 *                 @f$ \dot{s}@f$. By differentiating the algebraic equations with respect to time, we get the
 *                 missing linear equations (recall that the sensitivity vector @f$ s @f$ is fixed). The resulting system
 *                 has a similar structure as the system Jacobian.
 *                 @f[ \begin{align}
 *                  \left[\begin{array}{c|ccc|c}
 *                     \dot{J}_0  &         &        &           &   \\
 *                     \hline
 *                                & \dot{J}_1     &        &           &   \\
 *                                &         & \ddots &           &   \\
 *                                &         &        & \dot{J}_{N_z}   &   \\
 *                     \hline
 *                        J_{f,0} & J_{f,1} & \dots & J_{f,N_z} & I
 *                  \end{array}\right],
 *                 \end{align} @f]
 *                 where @f$ \dot{J}_i @f$ denotes the Jacobian with respect to @f$ \dot{y}@f$. Note that the
 *                 @f$ J_{i,f} @f$ matrices in the right column are missing.
 *                
 *     Let @f$ \mathcal{I}_d @f$ denote the index set of differential equations.
 *     The right hand side of the linear system is given by @f[ -\frac{\partial F}{\partial y}(t, y, \dot{y}) s - \frac{\partial F}{\partial p}(t, y, \dot{y}), @f]
 *     which is 0 for algebraic equations (@f$ -\frac{\partial^2 F}{\partial t \partial p}@f$, to be more precise).
 *     
 *     The linear system is solved by backsubstitution. First, the diagonal blocks are solved in parallel.
 *     Then, the equations for the fluxes @f$ j_f @f$ are solved by substituting in the solution of the
 *     diagonal blocks.</li>
 *          </ol>
 *     This function requires the parameter sensitivities to be computed beforehand and up-to-date Jacobians.
 * @param [in] simTime Simulation time information (time point, section index, pre-factor of time derivatives)
 * @param [in] simState Consistent state of the simulation (state vector and its time derivative)
 * @param [in,out] vecSensY Sensitivity subsystem state vectors
 * @param [in,out] vecSensYdot Time derivative state vectors of the sensitivity subsystems to be initialized
 * @param [in] adRes Pointer to residual vector of AD datatypes with parameter sensitivities
 * @todo Decrease amount of allocated memory by partially using temporary vectors (state and Schur complement)
 */
void GeneralRateModel::consistentInitialSensitivity(const ActiveSimulationTime& simTime, const ConstSimulationState& simState,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes)
{
	BENCH_SCOPE(_timerConsistentInit);

	Indexer idxr(_disc);

	for (unsigned int param = 0; param < vecSensY.size(); ++param)
	{
		double* const sensY = vecSensY[param];
		double* const sensYdot = vecSensYdot[param];

		// Copy parameter derivative dF / dp from AD and negate it
		for (unsigned int i = _disc.nComp; i < numDofs(); ++i)
			sensYdot[i] = -adRes[i].getADValue(param);

		// Step 1: Solve algebraic equations

		// Step 1a: Compute quasi-stationary binding model state
		for (unsigned int type = 0; type < _disc.nParType; ++type)
		{
			if (_binding[type]->hasAlgebraicEquations())
			{
#ifdef CADET_PARALLELIZE
				BENCH_SCOPE(_timerConsistentInitPar);
				tbb::parallel_for(size_t(0), size_t(_disc.nCol), [&](size_t pblk)
#else
				for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
#endif
				{
					// Get algebraic block
					unsigned int algStart = 0;
					unsigned int algLen = 0;
					_binding[type]->getAlgebraicBlock(algStart, algLen);

					// Reuse memory of band matrix for dense matrix
					linalg::DenseMatrixView jacobianMatrix(_jacPdisc[type * _disc.nCol + pblk].data(), _jacPdisc[type * _disc.nCol + pblk].pivot(), algLen, algLen);

					for (unsigned int shell = 0; shell < _disc.nParCell[type]; ++shell)
					{
						const unsigned int jacRowOffset = shell * static_cast<unsigned int>(idxr.strideParShell(type)) + static_cast<unsigned int>(idxr.strideParLiquid());
						const int localCpOffset = idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{static_cast<unsigned int>(pblk)}) + static_cast<int>(shell) * idxr.strideParShell(type);

						parts::BindingConsistentInitializer::consistentInitialSensitivityState(algStart, algLen, jacobianMatrix,
							_jacP[type * _disc.nCol + pblk], localCpOffset, jacRowOffset, idxr.strideParLiquid(), _disc.strideBound[type],
							sensY, sensYdot, _tempState + idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{static_cast<unsigned int>(pblk)}));
					}
				} CADET_PARFOR_END;
			}
		}

		// Step 1b: Compute fluxes j_f, right hand side is -dF / dp
		std::copy(sensYdot + idxr.offsetJf(), sensYdot + numDofs(), sensY + idxr.offsetJf());

		solveForFluxes(sensY, idxr);

		// Step 2: Compute the correct time derivative of the state vector

		// Step 2a: Assemble, factorize, and solve diagonal blocks of linear system

		// Compute right hand side by adding -dF / dy * s = -J * s to -dF / dp which is already stored in sensYdot
		multiplyWithJacobian(toSimple(simTime), simState, sensY, -1.0, 1.0, sensYdot);

		// Note that we have correctly negated the right hand side

		// Handle bulk block
		_convDispOp.solveTimeDerivativeSystem(toSimple(simTime), sensYdot + idxr.offsetC());

		// Process the particle blocks
#ifdef CADET_PARALLELIZE
		BENCH_START(_timerConsistentInitPar);
		tbb::parallel_for(size_t(0), size_t(_disc.nCol * _disc.nParType), [&](size_t pblk)
#else
		for (unsigned int pblk = 0; pblk < _disc.nCol * _disc.nParType; ++pblk)
#endif
		{
			const unsigned int type = pblk / _disc.nCol;
			const unsigned int par = pblk % _disc.nCol;

			// Assemble
			linalg::FactorizableBandMatrix& fbm = _jacPdisc[pblk];
			fbm.setAll(0.0);

			linalg::FactorizableBandMatrix::RowIterator jac = fbm.row(0);
			for (unsigned int j = 0; j < _disc.nParCell[type]; ++j)
			{
				// Mobile phase (advances jac accordingly)
				addMobilePhaseTimeDerivativeToJacobianParticleBlock(jac, idxr, 1.0, static_cast<double>(simTime.timeFactor), type);

				// Stationary phase
				// Populate matrix with time derivative Jacobian first
				_binding[type]->jacobianAddDiscretized(static_cast<double>(simTime.timeFactor), jac);

				// Overwrite rows corresponding to algebraic equations with the Jacobian and set right hand side to 0
				if (_binding[type]->hasAlgebraicEquations())
				{
					parts::BindingConsistentInitializer::consistentInitialSensitivityTimeDerivative(_binding[type], jac,
						_jacP[pblk].row(j * static_cast<unsigned int>(idxr.strideParShell(type)) + static_cast<unsigned int>(idxr.strideParLiquid())),
						sensYdot + idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{par}) + static_cast<int>(j) * idxr.strideParShell(type) + idxr.strideParLiquid());
				}

				// Advance pointers over all bound states
				jac += idxr.strideParBound(type);
			}   

			// Precondition
			double* const scaleFactors = _tempState + idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{par});
			fbm.rowScaleFactors(scaleFactors);
			fbm.scaleRows(scaleFactors);

			// Factorize
			const bool result = fbm.factorize();
			if (!result)
			{
				LOG(Error) << "Factorize() failed for par block " << pblk << " (type " << type << " col " << par << ")";
			}

			// Solve
			const bool result2 = fbm.solve(scaleFactors, sensYdot + idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{par}));
			if (!result2)
			{
				LOG(Error) << "Solve() failed for par block " << pblk << " (type " << type << " col " << par << ")";
			}
		} CADET_PARFOR_END;

#ifdef CADET_PARALLELIZE
		BENCH_STOP(_timerConsistentInitPar);
#endif

		// Step 2b: Solve for fluxes j_f by backward substitution
		solveForFluxes(sensYdot, idxr);
	}
}

/**
 * @brief Computes approximately / partially consistent initial values and time derivatives of sensitivity subsystems
 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] and initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$,
 *          the sensitivity system for a parameter @f$ p @f$ reads
 *          \f[ \frac{\partial F}{\partial y}(t, y, \dot{y}) s + \frac{\partial F}{\partial \dot{y}}(t, y, \dot{y}) \dot{s} + \frac{\partial F}{\partial p}(t, y, \dot{y}) = 0. \f]
 *          The initial values of this linear DAE, @f$ s_0 = \frac{\partial y_0}{\partial p} @f$ and @f$ \dot{s}_0 = \frac{\partial \dot{y}_0}{\partial p} @f$
 *          have to be consistent with the sensitivity DAE. This functions updates the initial sensitivity\f$ s_0 \f$ and overwrites the time
 *          derivative \f$ \dot{s}_0 \f$ such that they are consistent.
 *          
 *          The process follows closely the one of leanConsistentInitialConditions() and, in fact, is a linearized version of it.
 *          This is necessary because the initial conditions of the sensitivity system \f$ s_0 \f$ and \f$ \dot{s}_0 \f$ are
 *          related to the initial conditions \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ of the original DAE by differentiating them
 *          with respect to @f$ p @f$: @f$ s_0 = \frac{\partial y_0}{\partial p} @f$ and @f$ \dot{s}_0 = \frac{\partial \dot{y}_0}{\partial p}. @f$
 *          <ol>
 *              <li>Keep state and time derivative vectors as they are (i.e., do not solve algebraic equations).
 *                 Only solve for the fluxes @f$ j_{f,i} @f$ (only linear equations).</li>
 *              <li>Compute the time derivatives of the sensitivity @f$ \dot{s} @f$ such that the differential equations hold.
 *                 However, because of the algebraic equations, we need additional conditions to fully determine
 *                 @f$ \dot{s}@f$. By differentiating the algebraic equations with respect to time, we get the
 *                 missing linear equations (recall that the sensitivity vector @f$ s @f$ is fixed). The resulting
 *                 equations are stated below:
 *                 @f[ \begin{align}
 *                  \left[\begin{array}{c|ccc|c}
 *                     \dot{J}_0  &         &        &           &   \\
 *                     \hline
 *                        J_{f,0} & J_{f,1} & \dots & J_{f,N_z} & I
 *                  \end{array}\right],
 *                 \end{align} @f]
 *                 where @f$ \dot{J}_0 @f$ denotes the bulk block Jacobian with respect to @f$ \dot{y}@f$.
 *                
 *     Let @f$ \mathcal{I}_d @f$ denote the index set of differential equations.
 *     The right hand side of the linear system is given by @f[ -\frac{\partial F}{\partial y}(t, y, \dot{y}) s - \frac{\partial F}{\partial p}(t, y, \dot{y}), @f]
 *     which is 0 for algebraic equations (@f$ -\frac{\partial^2 F}{\partial t \partial p}@f$, to be more precise).
 *     
 *     The linear system is solved by backsubstitution. First, the bulk block is solved.
 *     Then, the equations for the fluxes @f$ j_f @f$ are solved by substituting in the solution of the
 *     bulk block and the unchanged particle block time derivative vectors.</li>
 *          </ol>
 *     This function requires the parameter sensitivities to be computed beforehand and up-to-date Jacobians.
 * @param [in] simTime Simulation time information (time point, section index, pre-factor of time derivatives)
 * @param [in] simState Consistent state of the simulation (state vector and its time derivative)
 * @param [in,out] vecSensY Sensitivity subsystem state vectors
 * @param [in,out] vecSensYdot Time derivative state vectors of the sensitivity subsystems to be initialized
 * @param [in] adRes Pointer to residual vector of AD datatypes with parameter sensitivities
 * @todo Decrease amount of allocated memory by partially using temporary vectors (state and Schur complement)
 */
void GeneralRateModel::leanConsistentInitialSensitivity(const ActiveSimulationTime& simTime, const ConstSimulationState& simState,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes)
{
	if ((_parDiffusion.size() > _disc.nComp * _disc.nParType) || (_parSurfDiffusion.size() > _disc.strideBound[_disc.nParType]))
		LOG(Warning) << "Lean consistent initialization is not appropriate for section-dependent pore and surface diffusion";

	BENCH_SCOPE(_timerConsistentInit);

	Indexer idxr(_disc);

	for (unsigned int param = 0; param < vecSensY.size(); ++param)
	{
		double* const sensY = vecSensY[param];
		double* const sensYdot = vecSensYdot[param];

		// Copy parameter derivative from AD to tempState and negate it
		// We need to use _tempState in order to keep sensYdot unchanged at this point
		for (unsigned int i = 0; i < idxr.offsetCp(); ++i)
			_tempState[i] = -adRes[i].getADValue(param);

		std::fill(_tempState + idxr.offsetCp(), _tempState + idxr.offsetJf(), 0.0);

		for (unsigned int i = idxr.offsetJf(); i < numDofs(); ++i)
			_tempState[i] = -adRes[i].getADValue(param);

		// Step 1: Compute fluxes j_f, right hand side is -dF / dp
		std::copy(_tempState + idxr.offsetJf(), _tempState + numDofs(), sensY + idxr.offsetJf());

		solveForFluxes(sensY, idxr);

		// Step 2: Compute the correct time derivative of the state vector

		// Step 2a: Assemble, factorize, and solve diagonal blocks of linear system

		// Compute right hand side by adding -dF / dy * s = -J * s to -dF / dp which is already stored in _tempState
		multiplyWithJacobian(toSimple(simTime), simState, sensY, -1.0, 1.0, _tempState);

		// Copy relevant parts to sensYdot for use as right hand sides
		std::copy(_tempState + idxr.offsetC(), _tempState + idxr.offsetCp(), sensYdot + idxr.offsetC());
		std::copy(_tempState + idxr.offsetJf(), _tempState + numDofs(), sensYdot);

		// Handle bulk block
		_convDispOp.solveTimeDerivativeSystem(toSimple(simTime), sensYdot + idxr.offsetC());

		// Step 2b: Solve for fluxes j_f by backward substitution
		solveForFluxes(sensYdot, idxr);
	}
}

/**
 * @brief Solves the algebraic flux equations for the fluxes @f$ j_f @f$
 * @details The equation to be solved is @f$ j_f - k_f * (c - c_p) == y @f$, where @f$ y @f$
 *          is a given vector.
 * @param [in,out] vecState On entry the state vector with @f$ y @f$ in its flux variables @f$ j_f @f$,
 *                 on exit the solution @f$ j_f. @f$
 * @param [in] idxr Indexer
 */
void GeneralRateModel::solveForFluxes(double* const vecState, const Indexer& idxr) const
{
	// We have j_f - k_f * (c - c_p) == 0 
	// Thus, jacFC contains -k_f and jacFP +k_f.
	// We just need to subtract both -k_f * c and k_f * c_p to get j_f == k_f * (c - c_p)

	double* const jf = vecState + idxr.offsetJf();

	// Note that we cannot parallelize this loop since we are updating the fluxes in-place
	_jacFC.multiplySubtract(vecState + idxr.offsetC(), jf);
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		linalg::DoubleSparseMatrix const* const jacFPtype = _jacFP + type * _disc.nCol;
		for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
			jacFPtype[pblk].multiplySubtract(vecState + idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{pblk}), jf);
	}
}

}  // namespace model

}  // namespace cadet
