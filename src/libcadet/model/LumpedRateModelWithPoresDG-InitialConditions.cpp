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

#include "model/LumpedRateModelWithPoresDG.hpp"
#include "model/BindingModel.hpp"
#include "ParamReaderHelper.hpp"
#include "AdUtils.hpp"
#include "model/parts/BindingCellKernel.hpp"
#include "SimulationTypes.hpp"
#include "SensParamUtil.hpp"
#include "linalg/Subset.hpp"

#include <algorithm>
#include <functional>

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include "ParallelSupport.hpp"
#ifdef CADET_PARALLELIZE
#include <tbb/parallel_for.h>
#endif

namespace cadet
{

	namespace model
	{

		int LumpedRateModelWithPoresDG::multiplexInitialConditions(const cadet::ParameterId& pId, unsigned int adDirection, double adValue)
		{
			if (_singleBinding)
			{
				if ((pId.name == hashString("INIT_CP")) && (pId.section == SectionIndep) && (pId.boundState == BoundStateIndep) && (pId.particleType == ParTypeIndep) && (pId.component != CompIndep) && (pId.reaction == ReactionIndep))
				{
					_sensParams.insert(&_initCp[pId.component]);
					for (unsigned int t = 0; t < _disc.nParType; ++t)
						_initCp[t * _disc.nComp + pId.component].setADValue(adDirection, adValue);
				}
				else if (pId.name == hashString("INIT_CP"))
					return -1;

				if ((pId.name == hashString("INIT_Q")) && (pId.section == SectionIndep) && (pId.boundState != BoundStateIndep) && (pId.particleType == ParTypeIndep) && (pId.component != CompIndep) && (pId.reaction == ReactionIndep))
				{
					_sensParams.insert(&_initQ[_disc.nBoundBeforeType[0] + _disc.boundOffset[pId.component] + pId.boundState]);
					for (unsigned int t = 0; t < _disc.nParType; ++t)
						_initQ[_disc.nBoundBeforeType[t] + _disc.boundOffset[t * _disc.nComp + pId.component] + pId.boundState].setADValue(adDirection, adValue);
				}
				else if (pId.name == hashString("INIT_Q"))
					return -1;
			}
			else
			{
				if ((pId.name == hashString("INIT_CP")) && (pId.section == SectionIndep) && (pId.boundState == BoundStateIndep) && (pId.particleType != ParTypeIndep) && (pId.component != CompIndep) && (pId.reaction == ReactionIndep))
				{
					_sensParams.insert(&_initCp[pId.particleType * _disc.nComp + pId.component]);
					_initCp[pId.particleType * _disc.nComp + pId.component].setADValue(adDirection, adValue);
				}
				else if (pId.name == hashString("INIT_CP"))
					return -1;

				if ((pId.name == hashString("INIT_Q")) && (pId.section == SectionIndep) && (pId.boundState != BoundStateIndep) && (pId.particleType != ParTypeIndep) && (pId.component != CompIndep) && (pId.reaction == ReactionIndep))
				{
					_sensParams.insert(&_initQ[_disc.nBoundBeforeType[pId.particleType] + _disc.boundOffset[pId.particleType * _disc.nComp + pId.component] + pId.boundState]);
					_initQ[_disc.nBoundBeforeType[pId.particleType] + _disc.boundOffset[pId.particleType * _disc.nComp + pId.component] + pId.boundState].setADValue(adDirection, adValue);
				}
				else if (pId.name == hashString("INIT_Q"))
					return -1;
			}
			return 0;
		}

		int LumpedRateModelWithPoresDG::multiplexInitialConditions(const cadet::ParameterId& pId, double val, bool checkSens)
		{
			if (_singleBinding)
			{
				if ((pId.name == hashString("INIT_CP")) && (pId.section == SectionIndep) && (pId.boundState == BoundStateIndep) && (pId.particleType == ParTypeIndep) && (pId.component != CompIndep) && (pId.reaction == ReactionIndep))
				{
					if (checkSens && !contains(_sensParams, &_initCp[pId.component]))
						return -1;

					for (unsigned int t = 0; t < _disc.nParType; ++t)
						_initCp[t * _disc.nComp + pId.component].setValue(val);
				}
				else if (pId.name == hashString("INIT_CP"))
					return -1;

				if ((pId.name == hashString("INIT_Q")) && (pId.section == SectionIndep) && (pId.boundState != BoundStateIndep) && (pId.particleType == ParTypeIndep) && (pId.component != CompIndep) && (pId.reaction == ReactionIndep))
				{
					if (checkSens && !contains(_sensParams, &_initQ[_disc.nBoundBeforeType[0] + _disc.boundOffset[pId.component] + pId.boundState]))
						return -1;

					for (unsigned int t = 0; t < _disc.nParType; ++t)
						_initQ[_disc.nBoundBeforeType[t] + _disc.boundOffset[t * _disc.nComp + pId.component] + pId.boundState].setValue(val);
				}
				else if (pId.name == hashString("INIT_Q"))
					return -1;
			}
			else
			{
				if ((pId.name == hashString("INIT_CP")) && (pId.section == SectionIndep) && (pId.boundState == BoundStateIndep) && (pId.particleType != ParTypeIndep) && (pId.component != CompIndep) && (pId.reaction == ReactionIndep))
				{
					if (checkSens && !contains(_sensParams, &_initCp[pId.particleType * _disc.nComp + pId.component]))
						return -1;

					_initCp[pId.particleType * _disc.nComp + pId.component].setValue(val);
				}
				else if (pId.name == hashString("INIT_CP"))
					return -1;

				if ((pId.name == hashString("INIT_Q")) && (pId.section == SectionIndep) && (pId.boundState != BoundStateIndep) && (pId.particleType != ParTypeIndep) && (pId.component != CompIndep) && (pId.reaction == ReactionIndep))
				{
					if (checkSens && !contains(_sensParams, &_initQ[_disc.nBoundBeforeType[pId.particleType] + _disc.boundOffset[pId.particleType * _disc.nComp + pId.component] + pId.boundState]))
						return -1;

					_initQ[_disc.nBoundBeforeType[pId.particleType] + _disc.boundOffset[pId.particleType * _disc.nComp + pId.component] + pId.boundState].setValue(val);
				}
				else if (pId.name == hashString("INIT_Q"))
					return -1;
			}
			return 0;
		}

		void LumpedRateModelWithPoresDG::applyInitialCondition(const SimulationState& simState) const
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
			for (unsigned int point = 0; point < _disc.nPoints; ++point)
			{
				// Loop over components in cell
				for (unsigned comp = 0; comp < _disc.nComp; ++comp)
					stateYbulk[point * idxr.strideColNode() + comp * idxr.strideColComp()] = static_cast<double>(_initC[comp]);
			}

			// Loop over particles
			for (unsigned int type = 0; type < _disc.nParType; ++type)
			{
				for (unsigned int point = 0; point < _disc.nPoints; ++point)
				{
					const unsigned int offset = idxr.offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ point });

					// Initialize c_p
					for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
						simState.vecStateY[offset + comp] = static_cast<double>(_initCp[comp + _disc.nComp * type]);

					// Initialize q
					active const* const iq = _initQ.data() + _disc.nBoundBeforeType[type];
					for (unsigned int bnd = 0; bnd < _disc.strideBound[type]; ++bnd)
						simState.vecStateY[offset + idxr.strideParLiquid() + bnd] = static_cast<double>(iq[bnd]);
				}
			}
		}

		void LumpedRateModelWithPoresDG::readInitialCondition(IParameterProvider& paramProvider)
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

			if (initC.size() < _disc.nComp)
				throw InvalidParameterException("INIT_C does not contain enough values for all components");

			ad::copyToAd(initC.data(), _initC.data(), _disc.nComp);

			// Check if INIT_CP is present
			if (paramProvider.exists("INIT_CP"))
			{
				const std::vector<double> initCp = paramProvider.getDoubleArray("INIT_CP");

				if (((initCp.size() < _disc.nComp) && _singleBinding) || ((initCp.size() < _disc.nComp * _disc.nParType) && !_singleBinding))
					throw InvalidParameterException("INIT_CP does not contain enough values for all components");

				if (!_singleBinding)
					ad::copyToAd(initCp.data(), _initCp.data(), _disc.nComp * _disc.nParType);
				else
				{
					for (unsigned int t = 0; t < _disc.nParType; ++t)
						ad::copyToAd(initCp.data(), _initCp.data() + t * _disc.nComp, _disc.nComp);
				}
			}
			else
			{
				for (unsigned int t = 0; t < _disc.nParType; ++t)
					ad::copyToAd(initC.data(), _initCp.data() + t * _disc.nComp, _disc.nComp);
			}

			std::vector<double> initQ;
			if (paramProvider.exists("INIT_Q"))
				initQ = paramProvider.getDoubleArray("INIT_Q");

			if (initQ.empty() || (_disc.strideBound[_disc.nParType] == 0))
				return;

			if ((_disc.strideBound[_disc.nParType] > 0) && (((initQ.size() < _disc.strideBound[_disc.nParType]) && !_singleBinding) || ((initQ.size() < _disc.strideBound[0]) && _singleBinding)))
				throw InvalidParameterException("INIT_Q does not contain enough values for all bound states");

			if (!_singleBinding)
				ad::copyToAd(initQ.data(), _initQ.data(), _disc.strideBound[_disc.nParType]);
			else
			{
				for (unsigned int t = 0; t < _disc.nParType; ++t)
					ad::copyToAd(initQ.data(), _initQ.data() + _disc.nBoundBeforeType[t], _disc.strideBound[t]);
			}
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
		 *                 No need to solve for fluxes @f$ j_{f,i} @f$ (since they are not part of the state for DG).</li>
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
		 *     No need to solve for fluxes @f$ j_{f,i} @f$ (since they are not part of the state for DG)</li>
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
		void LumpedRateModelWithPoresDG::consistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerConsistentInit);

			Indexer idxr(_disc);

			// Step 1: Solve algebraic equations

			// Step 1a: Compute quasi-stationary binding model state
			for (unsigned int type = 0; type < _disc.nParType; ++type)
			{
				if (!_binding[type]->hasQuasiStationaryReactions())
					continue;

				// Copy quasi-stationary binding mask to a local array that also includes the mobile phase
				std::vector<int> qsMask(_disc.nComp + _disc.strideBound[type], false);
				int const* const qsMaskSrc = _binding[type]->reactionQuasiStationarity();
				std::copy_n(qsMaskSrc, _disc.strideBound[type], qsMask.data() + _disc.nComp);

				// Activate mobile phase components that have at least one active bound state
				unsigned int bndStartIdx = 0;
				unsigned int numActiveComp = 0;
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				{
					for (unsigned int bnd = 0; bnd < _disc.nBound[_disc.nComp * type + comp]; ++bnd)
					{
						if (qsMaskSrc[bndStartIdx + bnd])
						{
							++numActiveComp;
							qsMask[comp] = true;
							break;
						}
					}

					bndStartIdx += _disc.nBound[_disc.nComp * type + comp];
				}

				const linalg::ConstMaskArray mask{ qsMask.data(), static_cast<int>(_disc.nComp + _disc.strideBound[type]) };
				const int probSize = linalg::numMaskActive(mask);

				//Problem capturing variables here
#ifdef CADET_PARALLELIZE
				BENCH_SCOPE(_timerConsistentInitPar);
				tbb::parallel_for(std::size_t(0), static_cast<std::size_t>(_disc.nPoints), [&](std::size_t pblk)
#else
				for (unsigned int pblk = 0; pblk < _disc.nPoints; ++pblk)
#endif
				{
					LinearBufferAllocator tlmAlloc = threadLocalMem.get();

					// Reuse memory of band matrix for dense matrix
					linalg::DenseMatrixView fullJacobianMatrix(_globalJacDisc.valuePtr() + _globalJacDisc.outerIndexPtr()[idxr.offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ static_cast<unsigned int>(pblk) }) - idxr.offsetC()], nullptr, mask.len, mask.len);

					// z coordinate (column length normed to 1) of current node - needed in externally dependent adsorption kinetic
					const double z = _convDispOp.relativeCoordinate(pblk);

					// Get workspace memory
					BufferedArray<double> nonlinMemBuffer = tlmAlloc.array<double>(_nonlinearSolver->workspaceSize(probSize));
					double* const nonlinMem = static_cast<double*>(nonlinMemBuffer);

					BufferedArray<double> solutionBuffer = tlmAlloc.array<double>(probSize);
					double* const solution = static_cast<double*>(solutionBuffer);

					BufferedArray<double> fullResidualBuffer = tlmAlloc.array<double>(mask.len);
					double* const fullResidual = static_cast<double*>(fullResidualBuffer);

					BufferedArray<double> fullXBuffer = tlmAlloc.array<double>(mask.len);
					double* const fullX = static_cast<double*>(fullXBuffer);

					BufferedArray<double> jacobianMemBuffer = tlmAlloc.array<double>(probSize * probSize);
					linalg::DenseMatrixView jacobianMatrix(static_cast<double*>(jacobianMemBuffer), _globalJacDisc.outerIndexPtr() + pblk * probSize, probSize, probSize);

					BufferedArray<double> conservedQuantsBuffer = tlmAlloc.array<double>(numActiveComp);
					double* const conservedQuants = static_cast<double*>(conservedQuantsBuffer);

					const parts::cell::CellParameters cellResParams
					{
						_disc.nComp,
						_disc.nBound + _disc.nComp * type,
						_disc.boundOffset + _disc.nComp * type,
						_disc.strideBound[type],
						_binding[type]->reactionQuasiStationarity(),
						_parPorosity[type],
						_poreAccessFactor.data() + _disc.nComp * type,
						_binding[type],
						(_dynReaction[type] && (_dynReaction[type]->numReactionsCombined() > 0)) ? _dynReaction[type] : nullptr
					};

					const int localOffsetToParticle = idxr.offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ static_cast<unsigned int>(pblk) });
					const int localOffsetInParticle = idxr.strideParLiquid();

					// Get pointer to q variables in a shell of particle pblk
					double* const qShell = vecStateY + localOffsetToParticle + localOffsetInParticle;
					active* const localAdRes = adJac.adRes ? adJac.adRes + localOffsetToParticle : nullptr;
					active* const localAdY = adJac.adY ? adJac.adY + localOffsetToParticle : nullptr;

					const ColumnPosition colPos{ z, 0.0, static_cast<double>(_parRadius[type]) * 0.5 };

					// Determine whether nonlinear solver is required
					if (!_binding[type]->preConsistentInitialState(simTime.t, simTime.secIdx, colPos, qShell, qShell - idxr.strideParLiquid(), tlmAlloc))
						CADET_PAR_CONTINUE;

					// Extract initial values from current state
					linalg::selectVectorSubset(qShell - _disc.nComp, mask, solution);

					// Save values of conserved moieties
					const double epsQ = 1.0 - static_cast<double>(_parPorosity[type]);
					linalg::conservedMoietiesFromPartitionedMask(mask, _disc.nBound + type * _disc.nComp, _disc.nComp, qShell - _disc.nComp, conservedQuants, static_cast<double>(_parPorosity[type]), epsQ);

					std::function<bool(double const* const, linalg::detail::DenseMatrixBase&)> jacFunc;
					if (localAdY && localAdRes) // todo fix AD consistent initialization with req. binding
					{
						jacFunc = [&](double const* const x, linalg::detail::DenseMatrixBase& mat)
						{
							// Copy over state vector to AD state vector (without changing directional values to keep seed vectors)
							// and initialize residuals with zero (also resetting directional values)
							ad::copyToAd(qShell - _disc.nComp, localAdY, mask.len);
							// @todo Check if this is necessary
							ad::resetAd(localAdRes, mask.len);

							// Prepare input vector by overwriting masked items
							linalg::applyVectorSubset(x, mask, localAdY);

							// Call residual function
							parts::cell::residualKernel<active, active, double, parts::cell::CellParameters, linalg::DenseBandedRowIterator, false, true>(
								simTime.t, simTime.secIdx, colPos, localAdY, nullptr, localAdRes, fullJacobianMatrix.row(0), cellResParams, tlmAlloc
								);

//							// todo check analytical jacobian
//#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
//							std::copy_n(qShell - _disc.nComp, mask.len, fullX);
//							linalg::applyVectorSubset(x, mask, fullX);
//
//							// Compute analytic Jacobian
//							parts::cell::residualKernel<double, double, double, parts::cell::CellParameters, linalg::DenseBandedRowIterator, true, true>(
//								simTime.t, simTime.secIdx, colPos, fullX, nullptr, fullResidual, fullJacobianMatrix.row(0), cellResParams, tlmAlloc
//								);
//
//							// Compare
//							const double diff = ad::compareDenseJacobianWithBandedAd(
//								adJac.adRes + idxr.offsetCp(ParticleTypeIndex{ type }), pblk * idxr.strideParBlock(type), adJac.adDirOffset, _jacP[type].lowerBandwidth(),
//								_jacP[type].lowerBandwidth(), _jacP[type].upperBandwidth(), fullJacobianMatrix
//							);
//							LOG(Debug) << "MaxDiff: " << diff;
//#endif
							// Extract Jacobian from AD
							// Read particle Jacobian entries from dedicated AD directions
							int offsetParticleTypeDirs = adJac.adDirOffset + _convDispOp.requiredADdirs();
							const active* const adRes = adJac.adRes;

							for (unsigned int type = 0; type < _disc.nParType; type++)
							{
								for (unsigned int par = 0; par < _disc.nPoints; par++)
								{
									const int eqOffset_res = idxr.offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ par });
									const int eqOffset_mat = idxr.offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ par }) - idxr.offsetC();
									for (unsigned int phase = 0; phase < idxr.strideParBlock(type); phase++)
									{
										for (unsigned int phaseTo = 0; phaseTo < idxr.strideParBlock(type); phaseTo++)
										{
											_globalJac.coeffRef(eqOffset_mat + phase, eqOffset_mat + phaseTo) = adRes[eqOffset_res + phase].getADValue(offsetParticleTypeDirs + phaseTo);
										}
									}
								}
								offsetParticleTypeDirs += idxr.strideParBlock(type);
							}

							// Extract Jacobian from full Jacobian
							mat.setAll(0.0);
							linalg::copyMatrixSubset(fullJacobianMatrix, mask, mask, mat);

							// Replace upper part with conservation relations
							mat.submatrixSetAll(0.0, 0, 0, numActiveComp, probSize);

							unsigned int bndIdx = 0;
							unsigned int rIdx = 0;
							unsigned int bIdx = 0;
							for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
							{
								if (!mask.mask[comp])
								{
									bndIdx += _disc.nBound[_disc.nComp * type + comp];
									continue;
								}

								mat.native(rIdx, rIdx) = static_cast<double>(_parPorosity[type]);

								for (unsigned int bnd = 0; bnd < _disc.nBound[_disc.nComp * type + comp]; ++bnd, ++bndIdx)
								{
									if (mask.mask[bndIdx])
									{
										mat.native(rIdx, bIdx + numActiveComp) = epsQ;
										++bIdx;
									}
								}

								++rIdx;
							}

							return true;
						};
					}
					else
					{
						jacFunc = [&](double const* const x, linalg::detail::DenseMatrixBase& mat)
						{
							// Prepare input vector by overwriting masked items
							std::copy_n(qShell - _disc.nComp, mask.len, fullX);
							linalg::applyVectorSubset(x, mask, fullX);

							// Call residual function
							parts::cell::residualKernel<double, double, double, parts::cell::CellParameters, linalg::DenseBandedRowIterator, true, true>(
								simTime.t, simTime.secIdx, colPos, fullX, nullptr, fullResidual, fullJacobianMatrix.row(0), cellResParams, tlmAlloc
								);

							// Extract Jacobian from full Jacobian
							mat.setAll(0.0);
							linalg::copyMatrixSubset(fullJacobianMatrix, mask, mask, mat);

							// Replace upper part with conservation relations
							mat.submatrixSetAll(0.0, 0, 0, numActiveComp, probSize);

							unsigned int bndIdx = 0;
							unsigned int rIdx = 0;
							unsigned int bIdx = 0;
							for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
							{
								if (!mask.mask[comp])
								{
									bndIdx += _disc.nBound[_disc.nComp * type + comp];
									continue;
								}

								mat.native(rIdx, rIdx) = static_cast<double>(_parPorosity[type]);

								for (unsigned int bnd = 0; bnd < _disc.nBound[_disc.nComp * type + comp]; ++bnd, ++bndIdx)
								{
									if (mask.mask[bndIdx])
									{
										mat.native(rIdx, bIdx + numActiveComp) = epsQ;
										++bIdx;
									}
								}

								++rIdx;
							}

							return true;
						};
					}

					// Apply nonlinear solver
					_nonlinearSolver->solve(
						[&](double const* const x, double* const r)
						{
							// Prepare input vector by overwriting masked items
							std::copy_n(qShell - _disc.nComp, mask.len, fullX);
							linalg::applyVectorSubset(x, mask, fullX);

							// Call residual function
							parts::cell::residualKernel<double, double, double, parts::cell::CellParameters, linalg::DenseBandedRowIterator, false, true>(
								simTime.t, simTime.secIdx, colPos, fullX, nullptr, fullResidual, fullJacobianMatrix.row(0), cellResParams, tlmAlloc
								);

							// Extract values from residual
							linalg::selectVectorSubset(fullResidual, mask, r);

							// Calculate residual of conserved moieties
							std::fill_n(r, numActiveComp, 0.0);
							unsigned int bndIdx = _disc.nComp;
							unsigned int rIdx = 0;
							unsigned int bIdx = 0;
							for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
							{
								if (!mask.mask[comp])
								{
									bndIdx += _disc.nBound[_disc.nComp * type + comp];
									continue;
								}

								r[rIdx] = static_cast<double>(_parPorosity[type]) * x[rIdx] - conservedQuants[rIdx];

								for (unsigned int bnd = 0; bnd < _disc.nBound[_disc.nComp * type + comp]; ++bnd, ++bndIdx)
								{
									if (mask.mask[bndIdx])
									{
										r[rIdx] += epsQ * x[bIdx + numActiveComp];
										++bIdx;
									}
								}

								++rIdx;
							}

							return true;
						},
						jacFunc, errorTol, solution, nonlinMem, jacobianMatrix, probSize);

					// Apply solution
					linalg::applyVectorSubset(solution, mask, qShell - idxr.strideParLiquid());

					// Refine / correct solution
					_binding[type]->postConsistentInitialState(simTime.t, simTime.secIdx, colPos, qShell, qShell - idxr.strideParLiquid(), tlmAlloc);

				} CADET_PARFOR_END;
			}

			// reset jacobian pattern //@todo can this be avoided?
			setGlobalJacPattern(_globalJacDisc, _dynReactionBulk);
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
		 *                 No need to solve for fluxes @f$ j_{f,i} @f$ (since they are not part of the state for DG).</li>
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
		 *     No need to solve for fluxes @f$ j_{f,i} @f$ (since they are not part of the state for DG).</li>
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
		void LumpedRateModelWithPoresDG::consistentInitialTimeDerivative(const SimulationTime& simTime, double const* vecStateY, double* const vecStateYdot, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerConsistentInit);

			Eigen::Map<const VectorXd> y(vecStateY, numDofs());
			Eigen::Map<VectorXd> yDot(vecStateYdot, numDofs());

			Indexer idxr(_disc);

			// Step 2: Compute the correct time derivative of the state vector

			// Step 2a: Assemble, factorize, and solve diagonal blocks of linear system

			// Note that the residual has not been negated, yet. We will do that now.
			for (unsigned int i = 0; i < numDofs(); ++i)
				vecStateYdot[i] = -vecStateYdot[i];

			// Handle bulk column block
			// Assemble
			double* vPtr = _globalJacDisc.valuePtr();
			for (int k = 0; k < _globalJacDisc.nonZeros(); k++) {
				vPtr[k] = 0.0;
			}
			_convDispOp.addTimeDerivativeToJacobian(1.0, _globalJacDisc);

			// Process the particle blocks
#ifdef CADET_PARALLELIZE
			BENCH_START(_timerConsistentInitPar);
			tbb::parallel_for(std::size_t(0), static_cast<std::size_t>(_disc.nParType), [&](std::size_t type)
#else
			for (unsigned int type = 0; type < _disc.nParType; ++type)
#endif
			{
				LinearBufferAllocator tlmAlloc = threadLocalMem.get();
				double* const dFluxDt = _tempState + idxr.offsetCp(ParticleTypeIndex{ static_cast<unsigned int>(type) });

				for (unsigned int pblk = 0; pblk < _disc.nPoints; ++pblk)
				{
					// z coordinate (column length normed to 1) of current node - needed in externally dependent adsorption kinetic
					const double z = _convDispOp.relativeCoordinate(pblk);

					// Assemble
					linalg::BandedEigenSparseRowIterator jacPar(_globalJacDisc, idxr.offsetCp(ParticleTypeIndex{ static_cast<unsigned int>(type) }, ParticleIndex{ pblk }) - idxr.offsetC());

					// Mobile and stationary phase (advances jac accordingly)
					addTimeDerivativeToJacobianParticleBlock(jacPar, idxr, 1.0, type);

					if (!_binding[type]->hasQuasiStationaryReactions())
						continue;

					// Get iterators to beginning of solid phase
					linalg::BandedEigenSparseRowIterator jacSolidOrig(_globalJac, idxr.offsetCp(ParticleTypeIndex{ static_cast<unsigned int>(type) }, ParticleIndex{ pblk }) - idxr.offsetC() + idxr.strideParLiquid());
					linalg::BandedEigenSparseRowIterator jacSolid(_globalJacDisc, idxr.offsetCp(ParticleTypeIndex{ static_cast<unsigned int>(type) }, ParticleIndex{ pblk }) - idxr.offsetC() + idxr.strideParLiquid());

					int const* const mask = _binding[type]->reactionQuasiStationarity();
					double* const qShellDot = vecStateYdot + idxr.offsetCp(ParticleTypeIndex{ static_cast<unsigned int>(type) }, ParticleIndex{ pblk }) + idxr.strideParLiquid();

					// Obtain derivative of fluxes wrt. time
					std::fill_n(dFluxDt, _disc.strideBound[type], 0.0);
					if (_binding[type]->dependsOnTime())
					{
						_binding[type]->timeDerivativeQuasiStationaryFluxes(simTime.t, simTime.secIdx,
							ColumnPosition{ z, 0.0, static_cast<double>(_parRadius[type]) * 0.5 },
							qShellDot - _disc.nComp, qShellDot, dFluxDt, tlmAlloc);
					}

					// Copy row from original Jacobian and set right hand side
					for (int i = 0; i < idxr.strideParBound(type); ++i, ++jacSolid, ++jacSolidOrig)
					{
						if (!mask[i])
							continue;

						jacSolid.copyRowFrom(jacSolidOrig);
						qShellDot[i] = -dFluxDt[i];
					}
				}

			} CADET_PARFOR_END;

#ifdef CADET_PARALLELIZE
			BENCH_STOP(_timerConsistentInitPar);
#endif

			// Factorize
			_linearSolver->factorize(_globalJacDisc);

			if (cadet_unlikely(_linearSolver->info() != Eigen::Success))
			{
				LOG(Error) << "Factorize() failed";
			}
			// Solve
			yDot.segment(idxr.offsetC(), numPureDofs()) = _linearSolver->solve(yDot.segment(idxr.offsetC(), numPureDofs()));

			if (cadet_unlikely(_linearSolver->info() != Eigen::Success))
			{
				LOG(Error) << "Solve() failed";
			}
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
		 *                 No need to solve for fluxes @f$ j_{f,i} @f$ (since they are not part of the state for DG).</li>
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
		 *     No need to solve for fluxes @f$ j_{f,i} @f$ (since they are not part of the state for DG).</li>
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
		void LumpedRateModelWithPoresDG::leanConsistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerConsistentInit);
			// Nothing to do here since there are no dedicated flux state variables for the DG LRMP discretization.
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
		 *                 No need to solve for fluxes @f$ j_{f,i} @f$ (since they are not part of the state for DG).</li>
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
		 *     No need to solve for fluxes @f$ j_{f,i} @f$ (since they are not part of the state for DG).</li>
		 *          </ol>
		 *     This function performs step 2. See leanConsistentInitialState() for step 1.
		 *
		 * 	   This function is to be used with leanConsistentInitialState(). Do not mix normal and lean
		 *     consistent initialization!
		 *
		 * @param [in] t Current time point
		 * @param [in] vecStateY (Lean) consistently initialized state vector
		 * @param [in,out] vecStateYdot On entry, inconsistent state time derivatives. On exit, partially consistent state time derivatives.
		 * @param [in] res On entry, residual without taking time derivatives into account. The data is overwritten during execution of the function.
		 */
		void LumpedRateModelWithPoresDG::leanConsistentInitialTimeDerivative(double t, double const* const vecStateY, double* const vecStateYdot, double* const res, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerConsistentInit);

			Indexer idxr(_disc);

			double* const resSlice = res + idxr.offsetC();

			// Step 2: Compute the correct time derivative of the state vector, i.e.
			// assemble, factorize, and solve column bulk block of linear system

			// Note that the residual is not negated as required at this point. We will fix that later.

			solveBulkTimeDerivativeSystem(SimulationTime{ t, 0u }, resSlice);

			// Note that we have solved with the *positive* residual as right hand side
			// instead of the *negative* one. Fortunately, we are dealing with linear systems,
			// which means that we can just negate the solution.
			double* const yDotSlice = vecStateYdot + idxr.offsetC();
			for (unsigned int i = 0; i < _disc.nElem * _disc.nComp; ++i)
				yDotSlice[i] = -resSlice[i];
		}

		void LumpedRateModelWithPoresDG::solveBulkTimeDerivativeSystem(const SimulationTime& simTime, double* const rhs) {
			
			Indexer idxr(_disc);

			// Assemble
			double* vals = _globalJacDisc.valuePtr();
			for (int entry = 0; entry < _globalJacDisc.nonZeros(); entry++)
				vals[entry] = 0.0;

			double alpha = 1.0;
			const int gapCell = idxr.strideColNode() - static_cast<int>(_disc.nComp) * idxr.strideColComp();
			linalg::BandedEigenSparseRowIterator jac(_globalJacDisc, 0);

			for (unsigned int point = 0; point < _disc.nPoints; ++point, jac += gapCell) {
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++jac) {
					// dc_b / dt in transport equation
					jac[0] += alpha;
				}
			}

			const int bulkRows = idxr.offsetCp() - idxr.offsetC();
			_linearSolver->analyzePattern(_globalJacDisc.block(0, 0, bulkRows, bulkRows));
			_linearSolver->factorize(_globalJacDisc.block(0, 0, bulkRows, bulkRows));

			if (_linearSolver->info() != Success) {
				LOG(Error) << "factorization failed in sensitivity initialization";
			}

			Eigen::Map<Eigen::VectorXd> ret_vec(rhs, bulkRows);
			ret_vec = _linearSolver->solve(ret_vec);

			// Use the factors to solve the linear system 
			if (_linearSolver->info() != Success) {
				LOG(Error) << "solve failed in sensitivity initialization";
			}

			// reset linear solver to global Jacobian
			_linearSolver->analyzePattern(_globalJacDisc);
		}

		void LumpedRateModelWithPoresDG::initializeSensitivityStates(const std::vector<double*>& vecSensY) const
		{
			Indexer idxr(_disc);
			for (std::size_t param = 0; param < vecSensY.size(); ++param)
			{
				double* const stateYbulk = vecSensY[param] + idxr.offsetC();

				// Loop over column cells
				for (unsigned int point = 0; point < _disc.nPoints; ++point)
				{
					// Loop over components in cell
					for (unsigned comp = 0; comp < _disc.nComp; ++comp)
						stateYbulk[point * idxr.strideColNode() + comp * idxr.strideColComp()] = _initC[comp].getADValue(param);
				}

				// Loop over particles
				for (unsigned int type = 0; type < _disc.nParType; ++type)
				{
					for (unsigned int point = 0; point < _disc.nPoints; ++point)
					{
						const unsigned int offset = idxr.offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ point });
						double* const stateYparticle = vecSensY[param] + offset;
						double* const stateYparticleSolid = stateYparticle + idxr.strideParLiquid();

						// Initialize c_p
						for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
							stateYparticle[comp] = _initCp[comp + _disc.nComp * type].getADValue(param);

						// Initialize q
						for (unsigned int bnd = 0; bnd < _disc.strideBound[type]; ++bnd)
							stateYparticleSolid[bnd] = _initQ[bnd + _disc.nBoundBeforeType[type]].getADValue(param);
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
		 *                 No need to solve for fluxes @f$ j_{f,i} @f$ (since they are not part of the state for DG). Let @f$ \mathcal{I}_a @f$ be the index set of algebraic equations, then, at this point, we have
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
		 *     No need to solve for fluxes @f$ j_{f,i} @f$ (since they are not part of the state for DG).</li>
		 *          </ol>
		 *     This function requires the parameter sensitivities to be computed beforehand and up-to-date Jacobians.
		 * @param [in] simTime Simulation time information (time point, section index, pre-factor of time derivatives)
		 * @param [in] simState Consistent state of the simulation (state vector and its time derivative)
		 * @param [in,out] vecSensY Sensitivity subsystem state vectors
		 * @param [in,out] vecSensYdot Time derivative state vectors of the sensitivity subsystems to be initialized
		 * @param [in] adRes Pointer to residual vector of AD datatypes with parameter sensitivities
		 * @todo Decrease amount of allocated memory by partially using temporary vectors (state and Schur complement)
		 */
		void LumpedRateModelWithPoresDG::consistentInitialSensitivity(const SimulationTime& simTime, const ConstSimulationState& simState,
			std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerConsistentInit);
		
			Indexer idxr(_disc);
		
			for (std::size_t param = 0; param < vecSensY.size(); ++param)
			{
				double* const sensY = vecSensY[param];
				double* const sensYdot = vecSensYdot[param];
		
				// Copy parameter derivative dF / dp from AD and negate it
				for (unsigned int i = _disc.nComp; i < numDofs(); ++i)
					sensYdot[i] = -adRes[i].getADValue(param);
		
				// Step 1: Compute quasi-stationary binding model state
				for (unsigned int type = 0; type < _disc.nParType; ++type)
				{
					if (!_binding[type]->hasQuasiStationaryReactions())
						continue;
		
					int const* const qsMask = _binding[type]->reactionQuasiStationarity();
					const linalg::ConstMaskArray mask{ qsMask, static_cast<int>(_disc.strideBound[type]) };
					const int probSize = linalg::numMaskActive(mask);
		
		#ifdef CADET_PARALLELIZE
					BENCH_SCOPE(_timerConsistentInitPar);
					tbb::parallel_for(std::size_t(0), static_cast<std::size_t>(_disc.nPoints), [&](std::size_t pblk)
		#else
					for (unsigned int pblk = 0; pblk < _disc.nPoints; ++pblk)
		#endif
					{
						// Reuse memory of band matrix for dense matrix
						linalg::DenseMatrixView jacobianMatrix(_globalJacDisc.valuePtr() + _globalJacDisc.outerIndexPtr()[idxr.offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ static_cast<unsigned int>(pblk) }) - idxr.offsetC()] + pblk * probSize * probSize, nullptr, probSize, probSize);

						//linalg::DenseMatrixView jacobianMatrix(_jacPdisc[type].data() + pblk * probSize * probSize, _jacPdisc[type].pivot() + pblk * probSize, probSize, probSize);
		
						// Get workspace memory
						LinearBufferAllocator tlmAlloc = threadLocalMem.get();
		
						BufferedArray<double> rhsBuffer = tlmAlloc.array<double>(probSize);
						double* const rhs = static_cast<double*>(rhsBuffer);
		
						BufferedArray<double> rhsUnmaskedBuffer = tlmAlloc.array<double>(idxr.strideParBound(type));
						double* const rhsUnmasked = static_cast<double*>(rhsUnmaskedBuffer);
		
						double* const maskedMultiplier = _tempState + idxr.offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ static_cast<unsigned int>(pblk) });
						double* const scaleFactors = _tempState + idxr.offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ static_cast<unsigned int>(pblk) });
		
						const int jacRowOffset = idxr.offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ static_cast<unsigned int>(pblk) }) - idxr.offsetC();
						const int jacColOffset = jacRowOffset;
						const int localQOffset = idxr.offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ static_cast<unsigned int>(pblk) }) + idxr.strideParLiquid();
		
						// Extract subproblem Jacobian from full Jacobian
						jacobianMatrix.setAll(0.0);
						linalg::copyMatrixSubset(_globalJac, mask, mask, jacRowOffset, jacColOffset, jacobianMatrix);
		
						// Construct right hand side
						linalg::selectVectorSubset(sensYdot + localQOffset, mask, rhs);
		
						// Zero out masked elements
						std::copy_n(sensY + localQOffset - idxr.strideParLiquid(), _disc.nComp + _disc.strideBound[type], maskedMultiplier);
						linalg::fillVectorSubset(maskedMultiplier + _disc.nComp, mask, 0.0);
		
						// Assemble right hand side
						Eigen::Map<VectorXd> maskedMultiplier_eigen(maskedMultiplier, idxr.strideParBlock(type));
						Eigen::Map<VectorXd> rhsUnmasked_eigen(rhsUnmasked, idxr.strideParBlock(type));
						rhsUnmasked_eigen = _globalJac.block(jacRowOffset, jacColOffset, idxr.strideParBlock(type), idxr.strideParBlock(type)) * maskedMultiplier_eigen;
						linalg::vectorSubsetAdd(rhsUnmasked, mask, -1.0, 1.0, rhs);

						// Precondition
						jacobianMatrix.rowScaleFactors(scaleFactors);
						jacobianMatrix.scaleRows(scaleFactors);
		
						// Solve
						jacobianMatrix.factorize();
						jacobianMatrix.solve(scaleFactors, rhs);
		
						// Write back
						linalg::applyVectorSubset(rhs, mask, sensY + localQOffset);
					} CADET_PARFOR_END;
				}
		
				// Step 2: Compute the correct time derivative of the state vector: Assemble, factorize, and solve diagonal blocks of linear system
		
				// Compute right hand side by adding -dF / dy * s = -J * s to -dF / dp which is already stored in sensYdot
				multiplyWithJacobian(simTime, simState, sensY, -1.0, 1.0, sensYdot);
		
				// Note that we have correctly negated the right hand side
		
				// Assemble bulk block
				double* vPtr = _globalJacDisc.valuePtr();
				for (int k = 0; k < _globalJacDisc.nonZeros(); k++) {
					vPtr[k] = 0.0;
				}
				_convDispOp.addTimeDerivativeToJacobian(1.0, _globalJacDisc);

				// Process the particle blocks
		#ifdef CADET_PARALLELIZE
				BENCH_START(_timerConsistentInitPar);
				tbb::parallel_for(std::size_t(0), static_cast<std::size_t>(_disc.nParType), [&](std::size_t type)
		#else
				for (unsigned int type = 0; type < _disc.nParType; ++type)
		#endif
				{
					for (unsigned int pblk = 0; pblk < _disc.nPoints; ++pblk)
					{
						// Assemble
						linalg::BandedEigenSparseRowIterator jacPar(_globalJacDisc, idxr.offsetCp(ParticleTypeIndex{ static_cast<unsigned int>(type) }, ParticleIndex{ pblk }) - idxr.offsetC());

						// Mobile and solid phase
						addTimeDerivativeToJacobianParticleBlock(jacPar, idxr, 1.0, type);
						// Iterator jac has already been advanced to next shell
		
						// Overwrite rows corresponding to algebraic equations with the Jacobian and set right hand side to 0
						if (_binding[type]->hasQuasiStationaryReactions())
						{
							// Get iterators to beginning of solid phase
							linalg::BandedEigenSparseRowIterator jacSolidOrig(_globalJac, idxr.offsetCp(ParticleTypeIndex{ static_cast<unsigned int>(type) }, ParticleIndex{ pblk }) - idxr.offsetC() + idxr.strideParLiquid());
							linalg::BandedEigenSparseRowIterator jacSolid(_globalJacDisc, idxr.offsetCp(ParticleTypeIndex{ static_cast<unsigned int>(type) }, ParticleIndex{ pblk }) - idxr.offsetC() + idxr.strideParLiquid());

							int const* const mask = _binding[type]->reactionQuasiStationarity();
							double* const qShellDot = sensYdot + idxr.offsetCp(ParticleTypeIndex{ static_cast<unsigned int>(type) }, ParticleIndex{ pblk }) + idxr.strideParLiquid();
		
							// Copy row from original Jacobian and set right hand side
							for (int i = 0; i < idxr.strideParBound(type); ++i, ++jacSolid, ++jacSolidOrig)
							{
								if (!mask[i])
									continue;
		
								jacSolid.copyRowFrom(jacSolidOrig);
		
								// Right hand side is -\frac{\partial^2 res(t, y, \dot{y})}{\partial p \partial t}
								// If the residual is not explicitly depending on time, this expression is 0
								// @todo This is wrong if external functions are used. Take that into account!
								qShellDot[i] = 0.0;
							}
						}
					}
		
				} CADET_PARFOR_END;
		
				Eigen::Map<VectorXd> yDot(sensYdot, numPureDofs());

				// Factorize
				_linearSolver->factorize(_globalJacDisc);

				if (cadet_unlikely(_linearSolver->info() != Eigen::Success))
				{
					LOG(Error) << "Factorize() failed";
				}
				// Solve
				yDot.segment(0, numPureDofs()) = _linearSolver->solve(yDot.segment(0, numPureDofs()));

				if (cadet_unlikely(_linearSolver->info() != Eigen::Success))
				{
					LOG(Error) << "Solve() failed";
				}

		#ifdef CADET_PARALLELIZE
				BENCH_STOP(_timerConsistentInitPar);
		#endif
		
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
		 *                 No need to solve for fluxes @f$ j_{f,i} @f$ (since they are not part of the state for DG).</li>
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
		 *     No need to solve for fluxes @f$ j_{f,i} @f$ (since they are not part of the state for DG).</li>
		 *          </ol>
		 *     This function requires the parameter sensitivities to be computed beforehand and up-to-date Jacobians.
		 * @param [in] simTime Simulation time information (time point, section index, pre-factor of time derivatives)
		 * @param [in] simState Consistent state of the simulation (state vector and its time derivative)
		 * @param [in,out] vecSensY Sensitivity subsystem state vectors
		 * @param [in,out] vecSensYdot Time derivative state vectors of the sensitivity subsystems to be initialized
		 * @param [in] adRes Pointer to residual vector of AD datatypes with parameter sensitivities
		 * @todo Decrease amount of allocated memory by partially using temporary vectors (state and Schur complement)
		 */
		void LumpedRateModelWithPoresDG::leanConsistentInitialSensitivity(const SimulationTime& simTime, const ConstSimulationState& simState,
			std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes, util::ThreadLocalStorage& threadLocalMem)
		{
			BENCH_SCOPE(_timerConsistentInit);

			Indexer idxr(_disc);

			for (std::size_t param = 0; param < vecSensY.size(); ++param)
			{
				double* const sensY = vecSensY[param];
				double* const sensYdot = vecSensYdot[param];

				// Copy parameter derivative from AD to tempState and negate it
				// We need to use _tempState in order to keep sensYdot unchanged at this point
				for (int i = 0; i < idxr.offsetCp(); ++i)
					_tempState[i] = -adRes[i].getADValue(param);

				// Compute the correct time derivative of the state vector, i.e.
				// assemble, factorize, and solve diagonal blocks of linear system

				// Compute right hand side by adding -dF / dy * s = -J * s to -dF / dp which is already stored in _tempState
				multiplyWithJacobian(simTime, simState, sensY, -1.0, 1.0, _tempState);

				// Copy relevant parts to sensYdot for use as right hand sides
				std::copy(_tempState + idxr.offsetC(), _tempState + idxr.offsetCp(), sensYdot + idxr.offsetC());

				// Handle bulk block
				solveBulkTimeDerivativeSystem(simTime, sensYdot + idxr.offsetC());
			}
		}

	}  // namespace model

}  // namespace cadet
