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

#include "model/FiniteBath.hpp"
#include "model/BindingModel.hpp"
#include "linalg/Subset.hpp"
#include "linalg/Norms.hpp"
#include "ParamReaderHelper.hpp"
#include "AdUtils.hpp"
#include "model/parts/BindingCellKernel.hpp"
#include "SimulationTypes.hpp"
#include "SensParamUtil.hpp"

#include <algorithm>
#include <functional>

#include "LoggingUtils.hpp"
#include "Logging.hpp"

namespace cadet
{

namespace model
{

int FiniteBath::multiplexInitialConditions(const cadet::ParameterId& pId, unsigned int adDirection, double adValue)
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

		if ((pId.name == hashString("INIT_CS")) && (pId.section == SectionIndep) && (pId.boundState != BoundStateIndep) && (pId.particleType == ParTypeIndep) && (pId.component != CompIndep) && (pId.reaction == ReactionIndep))
		{
			_sensParams.insert(&_initCs[_disc.nBoundBeforeType[0] + _disc.boundOffset[pId.component] + pId.boundState]);
			for (unsigned int t = 0; t < _disc.nParType; ++t)
				_initCs[_disc.nBoundBeforeType[t] + _disc.boundOffset[t * _disc.nComp + pId.component] + pId.boundState].setADValue(adDirection, adValue);
		}
		else if (pId.name == hashString("INIT_CS"))
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

		if ((pId.name == hashString("INIT_CS")) && (pId.section == SectionIndep) && (pId.boundState != BoundStateIndep) && (pId.particleType != ParTypeIndep) && (pId.component != CompIndep) && (pId.reaction == ReactionIndep))
		{
			_sensParams.insert(&_initCs[_disc.nBoundBeforeType[pId.particleType] + _disc.boundOffset[pId.particleType * _disc.nComp + pId.component] + pId.boundState]);
			_initCs[_disc.nBoundBeforeType[pId.particleType] + _disc.boundOffset[pId.particleType * _disc.nComp + pId.component] + pId.boundState].setADValue(adDirection, adValue);
		}
		else if (pId.name == hashString("INIT_CS"))
			return -1;
	}
	return 0;
}

int FiniteBath::multiplexInitialConditions(const cadet::ParameterId& pId, double val, bool checkSens)
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

		if ((pId.name == hashString("INIT_CS")) && (pId.section == SectionIndep) && (pId.boundState != BoundStateIndep) && (pId.particleType == ParTypeIndep) && (pId.component != CompIndep) && (pId.reaction == ReactionIndep))
		{
			if (checkSens && !contains(_sensParams, &_initCs[_disc.nBoundBeforeType[0] + _disc.boundOffset[pId.component] + pId.boundState]))
				return -1;

			for (unsigned int t = 0; t < _disc.nParType; ++t)
				_initCs[_disc.nBoundBeforeType[t] + _disc.boundOffset[t * _disc.nComp + pId.component] + pId.boundState].setValue(val);
		}
		else if (pId.name == hashString("INIT_CS"))
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

		if ((pId.name == hashString("INIT_CS")) && (pId.section == SectionIndep) && (pId.boundState != BoundStateIndep) && (pId.particleType != ParTypeIndep) && (pId.component != CompIndep) && (pId.reaction == ReactionIndep))
		{
			if (checkSens && !contains(_sensParams, &_initCs[_disc.nBoundBeforeType[pId.particleType] + _disc.boundOffset[pId.particleType * _disc.nComp + pId.component] + pId.boundState]))
				return -1;

			_initCs[_disc.nBoundBeforeType[pId.particleType] + _disc.boundOffset[pId.particleType * _disc.nComp + pId.component] + pId.boundState].setValue(val);
		}
		else if (pId.name == hashString("INIT_CS"))
			return -1;
	}
	return 0;
}

void FiniteBath::applyInitialCondition(const SimulationState& simState) const
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

	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		stateYbulk[comp] = static_cast<double>(_initC[comp]);

	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		const unsigned int offset = idxr.offsetCp(ParticleTypeIndex{ type });

		for (unsigned int shell = 0; shell < _disc.nParPoints[type]; ++shell)
		{
			const unsigned int shellOffset = offset + shell * idxr.strideParNode(type);

			// Initialize c^p
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				simState.vecStateY[shellOffset + comp] = static_cast<double>(_initCp[comp + _disc.nComp * type]);

			// Initialize c^s
			active const* const iq = _initCs.data() + _disc.nBoundBeforeType[type];
			for (unsigned int bnd = 0; bnd < _disc.strideBound[type]; ++bnd)
				simState.vecStateY[shellOffset + idxr.strideParLiquid() + bnd] = static_cast<double>(iq[bnd]);
		}
	}
}

void FiniteBath::readInitialCondition(IParameterProvider& paramProvider)
{
	_initState.clear();
	_initStateDot.clear();

	// Check if INIT_STATE is present
	if (paramProvider.exists("INIT_STATE"))
	{
		const std::vector<double> initState = paramProvider.getDoubleArray("INIT_STATE");
		_initState = std::vector<double>(initState.begin(), initState.begin() + numPureDofs());

		// Check if INIT_STATE contains the full state vector and its time derivative
		if (initState.size() >= 2 * numPureDofs())
			_initStateDot = std::vector<double>(initState.begin() + numPureDofs(), initState.begin() + 2 * numPureDofs());
		return;
	}

	const std::vector<double> initC = paramProvider.getDoubleArray("INIT_C");

	if (initC.size() < _disc.nComp)
		throw InvalidParameterException("INIT_C does not contain enough values for all components");

	ad::copyToAd(initC.data(), _initC.data(), _disc.nComp);

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		std::ostringstream oss;
		oss << "particle_type_" << std::setfill('0') << std::setw(3) << parType;

		if (!paramProvider.exists(oss.str()))
			continue;

		paramProvider.pushScope(oss.str());

		// Check if INIT_CP is present, otherwise copy from INIT_C
		if (paramProvider.exists("INIT_CP"))
		{
			const std::vector<double> initCp = paramProvider.getDoubleArray("INIT_CP");

			if (initCp.size() < _disc.nComp)
				throw InvalidParameterException("INIT_CP does not contain enough values for all components");

			ad::copyToAd(initCp.data(), _initCp.data() + parType * _disc.nComp, _disc.nComp);

			if (_singleBinding && (parType > 0))
			{
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				{
					if (initCp[comp] != static_cast<double>(_initCp[(parType - 1) * _disc.nComp + comp]))
						throw InvalidParameterException("Binding models were specified as particle type independent (see field BINDING_PARTYPE_DEPENDENT), but INIT_CP is different for particle type " + std::to_string(parType - 1) + " and " + std::to_string(parType));
				}
			}
		}
		else
			ad::copyToAd(initC.data(), _initCp.data() + parType * _disc.nComp, _disc.nComp);

		std::vector<double> initCs;
		if (paramProvider.exists("INIT_CS"))
			initCs = paramProvider.getDoubleArray("INIT_CS");

		if (initCs.empty() || (_disc.strideBound[parType] == 0))
		{
			paramProvider.popScope();
			continue;
		}

		if (initCs.size() < _disc.strideBound[parType])
			throw InvalidParameterException("INIT_CS does not contain enough values for all bound states");

		ad::copyToAd(initCs.data(), _initCs.data() + _disc.nBoundBeforeType[parType], _disc.strideBound[parType]);

		if (_singleBinding && (parType > 0))
		{
			for (unsigned int bnd = 0; bnd < _disc.strideBound[parType]; ++bnd)
			{
				if (initCs[bnd] != static_cast<double>(_initCs[_disc.nBoundBeforeType[parType - 1] + bnd]))
					throw InvalidParameterException("Binding models were specified as particle type independent (see field BINDING_PARTYPE_DEPENDENT), but INIT_CS is different for particle type " + std::to_string(parType - 1) + " and " + std::to_string(parType));
			}
		}

		paramProvider.popScope();
	}
}

void FiniteBath::consistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerConsistentInit);

	consistentInitialBulkLiquidEquilibrium(simTime, vecStateY, errorTol, threadLocalMem);

	// Initialize quasi-stationary binding states
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		if (!_binding[type]->hasQuasiStationaryReactions())
			continue;

		consistentInitialBindingEquilibrium(simTime, vecStateY, adJac, errorTol, threadLocalMem, type);
	}

	// The binding equilibrium above writes plain values into the AD state vector, which resets the
	// seed vectors (see sfad::Fwd::operator=), so they have to be restored for the next AD evaluation
	if (adJac.adY)
		prepareADvectors(adJac);
}

void FiniteBath::consistentInitialBulkLiquidEquilibrium(const SimulationTime& simTime, double* const vecStateY, double errorTol, util::ThreadLocalStorage& threadLocalMem)
{
	Indexer idxr(_disc);

	const auto& cm = _reaction.conservedMoieties("liquid");
	if (!cm.isEnabled() || (cm.numEquilibriumReactions() == 0))
		return;

	const auto& L = cm.conservedMoietyMatrix();
	const unsigned int nMoieties = cm.numMoieties();
	const unsigned int nEq = cm.numEquilibriumReactions();
	const unsigned int probSize = _disc.nComp;

	if (nMoieties + nEq != probSize)
		throw InvalidParameterException("FiniteBath consistent initialization: Invalid equilibrium initialization size");

	linalg::DenseMatrix jacobianMatrix;
	jacobianMatrix.resize(probSize, probSize);

	LinearBufferAllocator tlmAlloc = threadLocalMem.get();
	double* const cLocal = vecStateY + idxr.offsetC();
	const ColumnPosition colPos{ 0.0, 0.0, 0.0 };

	BufferedArray<double> solutionBuffer = tlmAlloc.array<double>(probSize);
	double* const solution = static_cast<double*>(solutionBuffer);

	std::copy_n(cLocal, probSize, solution);

	BufferedArray<double> conservedBuffer = tlmAlloc.array<double>(nMoieties);
	double* const conserved = static_cast<double*>(conservedBuffer);

	cm.applyToVector(conserved, cLocal, probSize);

	BufferedArray<double> baNonlinMem = tlmAlloc.array<double>(_nonlinearSolver->workspaceSize(probSize));
	double* const nonlinMem = static_cast<double*>(baNonlinMem);

	std::function<bool(double const* const, linalg::detail::DenseMatrixBase&)> jacFunc;
	jacFunc = [&](double const* const x, linalg::detail::DenseMatrixBase& mat)
	{
		mat.setAll(0.0);
		// Upper part: Conserved moieties
		for (unsigned int m = 0; m < nMoieties; ++m)
		{
			for (unsigned int i = 0; i < probSize; ++i)
				mat.native(m, i) = L(m, i);
		}

		unsigned int eqIdx = 0;
		for (auto const* reaction : _reaction.getDynReactionVector("liquid"))
		{
			if (!reaction)
				continue;

			reaction->analyticEquilibriumJacobian(simTime.t, simTime.secIdx, colPos, probSize, x, eqIdx, nMoieties,
				mat.row(nMoieties), tlmAlloc.manageRemainingMemory());
		}

		return eqIdx == nEq;
	};

	std::function<bool(double const* const x, double* const r)> resFunc;
	resFunc = [&](double const* const x, double* const r)
	{
		// Upper part: Conserved moieties
		cm.applyToVector(r, x, probSize);
		for (unsigned int m = 0; m < nMoieties; ++m)
			r[m] -= conserved[m];

		unsigned int eqIdx = 0;
		for (auto const* reaction : _reaction.getDynReactionVector("liquid"))
		{
			if (!reaction)
				continue;
			reaction->residualEquilibriumFlux(simTime.t, simTime.secIdx, colPos, probSize, x, r + nMoieties, eqIdx, tlmAlloc.manageRemainingMemory());
		}
		return eqIdx == nEq;
	};

	const bool success = _nonlinearSolver->solve(resFunc, jacFunc, errorTol, solution, nonlinMem, jacobianMatrix, probSize);

	// Composite solvers can report a failed intermediate method even if a
	// subsequent method has produced a valid final iterate.
	const bool validFinalIterate = success
		|| (resFunc(solution, nonlinMem) && (linalg::linfNorm(nonlinMem, probSize) <= errorTol));

	if (validFinalIterate)
		std::copy_n(solution, probSize, cLocal);
	else
		LOG(Error) << "Consistent liquid equilibrium initialization failed";
}

void FiniteBath::consistentInitialBindingEquilibrium(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, util::ThreadLocalStorage& threadLocalMem, unsigned int type)
{
	Indexer idxr(_disc);

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

	// AD directions of this particle type, see prepareADvectors()
	unsigned int parDirOffset = adJac.adDirOffset + _disc.nComp;
	for (unsigned int t = 0; t < type; ++t)
		parDirOffset += idxr.strideParBlock(t);

	LinearBufferAllocator tlmAlloc = threadLocalMem.get();

	BufferedArray<double> nonlinMemBuffer = tlmAlloc.array<double>(_nonlinearSolver->workspaceSize(probSize));
	double* const nonlinMem = static_cast<double*>(nonlinMemBuffer);

	BufferedArray<double> solutionBuffer = tlmAlloc.array<double>(probSize);
	double* const solution = static_cast<double*>(solutionBuffer);

	BufferedArray<double> fullResidualBuffer = tlmAlloc.array<double>(mask.len);
	double* const fullResidual = static_cast<double*>(fullResidualBuffer);

	BufferedArray<double> fullXBuffer = tlmAlloc.array<double>(mask.len);
	double* const fullX = static_cast<double*>(fullXBuffer);

	BufferedArray<double> fullJacobianBuffer = tlmAlloc.array<double>(mask.len * mask.len);
	linalg::DenseMatrixView fullJacobianMatrix(static_cast<double*>(fullJacobianBuffer), nullptr, mask.len, mask.len);

	BufferedArray<double> conservedQuantsBuffer = tlmAlloc.array<double>(numActiveComp);
	double* const conservedQuants = static_cast<double*>(conservedQuantsBuffer);

	// The nonlinear solver factorizes this matrix, so it has to own its pivot array
	linalg::DenseMatrix jacobianMatrix;
	jacobianMatrix.resize(probSize, probSize);

	const parts::cell::CellParameters cellResParams = makeCellResidualParams(type, mask.mask + _disc.nComp);

	const double epsQ = 1.0 - static_cast<double>(_particles[type]->getPorosity());
	const int localOffsetToParticle = idxr.offsetCp(ParticleTypeIndex{ type });

	for (unsigned int node = 0; node < _disc.nParPoints[type]; ++node)
	{
		const int localOffsetInParticle = static_cast<int>(node) * idxr.strideParNode(type);

		// Get pointer to q variables in a shell of the particle
		double* const qShell = vecStateY + localOffsetToParticle + localOffsetInParticle + idxr.strideParLiquid();
		active* const localAdRes = adJac.adRes ? adJac.adRes + localOffsetToParticle + localOffsetInParticle : nullptr;
		active* const localAdY = adJac.adY ? adJac.adY + localOffsetToParticle + localOffsetInParticle : nullptr;

		// r (particle) coordinate of current node
		const double r = _particles[type]->relativeCoordinate(node);
		const ColumnPosition colPos{ 0.0, 0.0, r };

		// Determine whether nonlinear solver is required
		if (!_binding[type]->preConsistentInitialState(simTime.t, simTime.secIdx, colPos, qShell, qShell - idxr.strideParLiquid(), tlmAlloc))
			continue;

		// Extract initial values from current state
		linalg::selectVectorSubset(qShell - _disc.nComp, mask, solution);

		// Save values of conserved moieties
		linalg::conservedMoietiesFromPartitionedMask(mask, _disc.nBound + type * _disc.nComp, _disc.nComp, qShell - _disc.nComp, conservedQuants, static_cast<double>(_particles[type]->getPorosity()), epsQ);

		// Replaces the upper part of the masked Jacobian by the conservation relations
		auto applyConservationRows = [&](linalg::detail::DenseMatrixBase& mat)
		{
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

				mat.native(rIdx, rIdx) = static_cast<double>(_particles[type]->getPorosity());

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
		};

		std::function<bool(double const* const, linalg::detail::DenseMatrixBase&)> jacFunc;

		if (localAdY && localAdRes)
		{
			jacFunc = [&](double const* const x, linalg::detail::DenseMatrixBase& mat)
			{
				// Copy over state vector to AD state vector (without changing directional values to keep seed vectors)
				// and initialize residuals with zero (also resetting directional values)
				ad::copyToAd(qShell - _disc.nComp, localAdY, mask.len);
				ad::resetAd(localAdRes, mask.len);

				// Prepare input vector by overwriting masked items
				linalg::applyVectorSubset(x, mask, localAdY);

				parts::cell::residualKernel<active, active, double, parts::cell::CellParameters, linalg::DenseBandedRowIterator, false, true>(
					simTime.t, simTime.secIdx, colPos, localAdY, nullptr, localAdRes, fullJacobianMatrix.row(0), cellResParams, tlmAlloc
					);

				// Extract the dense cell Jacobian from the dedicated AD directions of this particle type
				for (int row = 0; row < mask.len; ++row)
				{
					for (int col = 0; col < mask.len; ++col)
						fullJacobianMatrix.native(row, col) = localAdRes[row].getADValue(parDirOffset + localOffsetInParticle + col);
				}

				// Extract Jacobian from full Jacobian
				mat.setAll(0.0);
				linalg::copyMatrixSubset(fullJacobianMatrix, mask, mask, mat);

				applyConservationRows(mat);

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

				// The residual kernel only adds to the Jacobian, so the scratch matrix has to be cleared
				fullJacobianMatrix.setAll(0.0);

				parts::cell::residualKernel<double, double, double, parts::cell::CellParameters, linalg::DenseBandedRowIterator, true, true>(
					simTime.t, simTime.secIdx, colPos, fullX, nullptr, fullResidual, fullJacobianMatrix.row(0), cellResParams, tlmAlloc
					);

				// Extract Jacobian from full Jacobian
				mat.setAll(0.0);
				linalg::copyMatrixSubset(fullJacobianMatrix, mask, mask, mat);

				applyConservationRows(mat);

				return true;
			};
		}

		// Apply nonlinear solver
		_nonlinearSolver->solve(
			[&](double const* const x, double* const res)
			{
				// Prepare input vector by overwriting masked items
				std::copy_n(qShell - _disc.nComp, mask.len, fullX);
				linalg::applyVectorSubset(x, mask, fullX);

				parts::cell::residualKernel<double, double, double, parts::cell::CellParameters, linalg::DenseBandedRowIterator, false, true>(
					simTime.t, simTime.secIdx, colPos, fullX, nullptr, fullResidual, fullJacobianMatrix.row(0), cellResParams, tlmAlloc
					);

				// Extract values from residual
				linalg::selectVectorSubset(fullResidual, mask, res);

				// Calculate residual of conserved moieties
				std::fill_n(res, numActiveComp, 0.0);
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

					res[rIdx] = static_cast<double>(_particles[type]->getPorosity()) * x[rIdx] - conservedQuants[rIdx];

					for (unsigned int bnd = 0; bnd < _disc.nBound[_disc.nComp * type + comp]; ++bnd, ++bndIdx)
					{
						if (mask.mask[bndIdx])
						{
							res[rIdx] += epsQ * x[bIdx + numActiveComp];
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
	}
}

void FiniteBath::consistentInitialTimeDerivative(const SimulationTime& simTime, double const* vecStateY, double* const vecStateYdot, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerConsistentInit);

	Indexer idxr(_disc);

	// Rows are copied from the system Jacobian below, which requires both matrices to have the same
	// sparsity pattern. An AD Jacobian may have added entries to it, so take the pattern over here.
	_globalJacDisc = _globalJac;

	double* entries = _globalJacDisc.valuePtr();
	for (int entry = 0; entry < _globalJacDisc.nonZeros(); ++entry)
		entries[entry] = 0.0;

	Eigen::Map<Eigen::VectorXd> yDot(vecStateYdot, numDofs());

	// Note that the residual has not been negated, yet. We will do that now.
	for (unsigned int i = 0; i < numDofs(); ++i)
		vecStateYdot[i] = -vecStateYdot[i];

	consistentInitialBulkTimeDerivative(simTime, vecStateY, vecStateYdot, threadLocalMem);

	for (unsigned int type = 0; type < _disc.nParType; ++type)
		consistentInitialBindingTimeDerivative(simTime, vecStateYdot, threadLocalMem, type);

	_linearSolver->factorize(_globalJacDisc.block(idxr.offsetC(), idxr.offsetC(), numPureDofs(), numPureDofs()));
	if (cadet_unlikely(_linearSolver->info() != Eigen::Success))
	{
		LOG(Error) << "Factorize() failed";
	}

	yDot.segment(idxr.offsetC(), numPureDofs()) = _linearSolver->solve(yDot.segment(idxr.offsetC(), numPureDofs()));
	if (cadet_unlikely(_linearSolver->info() != Eigen::Success))
	{
		LOG(Error) << "Solve() failed";
	}
}

void FiniteBath::consistentInitialBulkTimeDerivative(const SimulationTime& simTime, double const* vecStateY, double* const vecStateYdot, util::ThreadLocalStorage& threadLocalMem)
{
	Indexer idxr(_disc);

	const auto& cm = _reaction.conservedMoieties("liquid");
	if (cm.isEnabled() && cm.numEquilibriumReactions() > 0)
	{
		const auto& L = cm.conservedMoietyMatrix();
		const unsigned int nMoieties = cm.numMoieties();
		const unsigned int nEq = cm.numEquilibriumReactions();

		LinearBufferAllocator tlmAlloc = threadLocalMem.get();

		double const* const cLocal = vecStateY + idxr.offsetC();
		double* const cDotLocal = vecStateYdot + idxr.offsetC();

		const ColumnPosition colPos{ 0.0, 0.0, 0.0 };

		// Differential rows: dF / dyDot = L
		for (unsigned int m = 0; m < nMoieties; ++m)
		{
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				_globalJacDisc.coeffRef(idxr.offsetC() + m, idxr.offsetC() + comp) = L(m, comp);
		}

		// Algebraic rows: dg_eq / dc * cDot = -dg_eq / dt
		linalg::BandedEigenSparseRowIterator eqJac(_globalJacDisc, idxr.offsetC() + nMoieties);
		unsigned int eqIdx = 0;
		for (auto const* reaction : _reaction.getDynReactionVector("liquid"))
		{
			if (!reaction)
				continue;

			reaction->analyticEquilibriumJacobian(simTime.t, simTime.secIdx, colPos, _disc.nComp, cLocal, eqIdx, nMoieties,
				eqJac, tlmAlloc.manageRemainingMemory());
		}

		cadet_assert(eqIdx == nEq);
		// Explicitly time-dependent equilibrium equations are not supported yet
		std::fill_n(cDotLocal + nMoieties, nEq, 0.0);
	}
	else
	{
		linalg::BandedEigenSparseRowIterator jacBlk(_globalJacDisc, idxr.offsetC());
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++jacBlk)
			jacBlk[0] = 1.0;
	}
}

void FiniteBath::consistentInitialBindingTimeDerivative(const SimulationTime& simTime, double* const vecStateYdot, util::ThreadLocalStorage& threadLocalMem, unsigned int type)
{
	Indexer idxr(_disc);

	linalg::BandedEigenSparseRowIterator jacPar(_globalJacDisc, idxr.offsetCp(ParticleTypeIndex{ type }));

	LinearBufferAllocator tlmAlloc = threadLocalMem.get();
	double* const dFluxDt = _tempState + idxr.offsetCp(ParticleTypeIndex{ type });

	for (unsigned int j = 0; j < _disc.nParPoints[type]; ++j)
	{
		addTimeDerivativeToJacobianParticleShell(jacPar, idxr, 1.0, type); // Iterator jacPar advances to next node

		if (!_binding[type]->hasQuasiStationaryReactions())
			continue;

		const int shellOffset = idxr.offsetCp(ParticleTypeIndex{ type }) + static_cast<int>(j) * idxr.strideParNode(type) + idxr.strideParLiquid();

		linalg::BandedEigenSparseRowIterator jacSolidOrig(_globalJac, shellOffset);
		linalg::BandedEigenSparseRowIterator jacSolid(_globalJacDisc, shellOffset);

		int const* const mask = _binding[type]->reactionQuasiStationarity();
		double* const qShellDot = vecStateYdot + shellOffset;

		// Obtain derivative of fluxes wrt. time
		std::fill_n(dFluxDt, _disc.strideBound[type], 0.0);
		if (_binding[type]->dependsOnTime())
		{
			const double r = _particles[type]->relativeCoordinate(j);

			_binding[type]->timeDerivativeQuasiStationaryFluxes(simTime.t, simTime.secIdx,
				ColumnPosition{ 0.0, 0.0, r },
				qShellDot - _disc.nComp, qShellDot, dFluxDt, tlmAlloc);
		}

		// Copy row from original Jacobian (without time derivatives) and set right hand side
		for (int i = 0; i < idxr.strideParBound(type); ++i, ++jacSolid, ++jacSolidOrig)
		{
			if (!mask[i])
				continue;

			jacSolid.copyRowFrom(jacSolidOrig);
			qShellDot[i] = -dFluxDt[i];
		}
	}
}

void FiniteBath::leanConsistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerConsistentInit);

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		_particles[parType]->leanConsistentInitialStateValidity();

	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		if (!_binding[type]->hasQuasiStationaryReactions())
			continue;

		leanConsistentInitialBindingEquilibrium(simTime, vecStateY, threadLocalMem, type);
	}
}

void FiniteBath::leanConsistentInitialBindingEquilibrium(const SimulationTime& simTime, double* const vecStateY, util::ThreadLocalStorage& threadLocalMem, unsigned int type)
{
	Indexer idxr(_disc);

	LinearBufferAllocator tlmAlloc = threadLocalMem.get();

	const int localOffsetToParticle = idxr.offsetCp(ParticleTypeIndex{ type });
	for (unsigned int shell = 0; shell < _disc.nParPoints[type]; ++shell)
	{
		const int localOffsetInParticle = static_cast<int>(shell) * idxr.strideParNode(type) + idxr.strideParLiquid();
		double* const qShell = vecStateY + localOffsetToParticle + localOffsetInParticle;

		const double r = _particles[type]->relativeCoordinate(shell);
		const ColumnPosition colPos{ 0.0, 0.0, r };

		// Perform consistent initialization that does not require a full fledged nonlinear solver
		if (!_binding[type]->preConsistentInitialState(simTime.t, simTime.secIdx, colPos, qShell, qShell - idxr.strideParLiquid(), tlmAlloc))
			continue;
	}
}

void FiniteBath::leanConsistentInitialTimeDerivative(double t, double const* const vecStateY, double* const vecStateYdot, double* const res, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerConsistentInit);

	Indexer idxr(_disc);

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		_particles[parType]->leanConsistentInitialTimeDerivativeValidity();

	// Note that the residual is not negated as required at this point. We will fix that later.
	double* const resSlice = res + idxr.offsetC();

	solveBulkTimeDerivativeSystem(resSlice);

	// Note that we have solved with the *positive* residual as right hand side
	// instead of the *negative* one. Fortunately, we are dealing with linear systems,
	// which means that we can just negate the solution.
	double* const yDotSlice = vecStateYdot + idxr.offsetC();
	for (unsigned int i = 0; i < _disc.nComp; ++i)
		yDotSlice[i] = -resSlice[i];
}

void FiniteBath::initializeSensitivityStates(const std::vector<double*>& vecSensY) const
{
	for (std::size_t param = 0; param < vecSensY.size(); ++param)
	{
		initializeSensitivityBulkStates(vecSensY[param], param);

		for (unsigned int type = 0; type < _disc.nParType; ++type)
			initializeSensitivityParticleStates(vecSensY[param], param, type);
	}
}

void FiniteBath::initializeSensitivityBulkStates(double* const sensY, std::size_t param) const
{
	Indexer idxr(_disc);
	double* const stateYbulk = sensY + idxr.offsetC();

	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		stateYbulk[comp] = _initC[comp].getADValue(param);
}

void FiniteBath::initializeSensitivityParticleStates(double* const sensY, std::size_t param, unsigned int type) const
{
	Indexer idxr(_disc);

	const unsigned int offset = idxr.offsetCp(ParticleTypeIndex{ type });

	for (unsigned int shell = 0; shell < _disc.nParPoints[type]; ++shell)
	{
		const unsigned int shellOffset = offset + shell * idxr.strideParNode(type);
		double* const stateYparticle = sensY + shellOffset;
		double* const stateYparticleSolid = stateYparticle + idxr.strideParLiquid();

		// Initialize c^p
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			stateYparticle[comp] = _initCp[comp + type * _disc.nComp].getADValue(param);

		// Initialize c^s
		for (unsigned int bnd = 0; bnd < _disc.strideBound[type]; ++bnd)
			stateYparticleSolid[bnd] = _initCs[bnd + _disc.nBoundBeforeType[type]].getADValue(param);
	}
}

void FiniteBath::consistentInitialSensitivity(const SimulationTime& simTime, const ConstSimulationState& simState,
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

		// Step 1: Solve algebraic equations
		for (unsigned int type = 0; type < _disc.nParType; ++type)
		{
			if (!_binding[type]->hasQuasiStationaryReactions())
				continue;

			consistentInitialSensitivityBindingEquilibrium(sensY, sensYdot, threadLocalMem, type);
		}

		// Step 2: Compute the correct time derivative of the sensitivity vector

		// Compute right hand side by adding -dF / dy * s = -J * s to -dF / dp which is already stored in sensYdot
		multiplyWithJacobian(simTime, simState, sensY, -1.0, 1.0, sensYdot);

		consistentInitialSensitivityBulkTimeDerivative();

		for (unsigned int type = 0; type < _disc.nParType; ++type)
			consistentInitialSensitivityBindingTimeDerivative(sensYdot, type);

		Eigen::Map<Eigen::VectorXd> yDot(sensYdot + idxr.offsetC(), numPureDofs());

		_linearSolver->factorize(_globalJacDisc.block(idxr.offsetC(), idxr.offsetC(), numPureDofs(), numPureDofs()));

		if (cadet_unlikely(_linearSolver->info() != Eigen::Success))
		{
			LOG(Error) << "Factorize() failed";
		}

		yDot.segment(0, numPureDofs()) = _linearSolver->solve(yDot.segment(0, numPureDofs()));

		if (cadet_unlikely(_linearSolver->info() != Eigen::Success))
		{
			LOG(Error) << "Solve() failed";
		}
	}
}

void FiniteBath::consistentInitialSensitivityBindingEquilibrium(double* const sensY, double const* const sensYdot, util::ThreadLocalStorage& threadLocalMem, unsigned int type)
{
	Indexer idxr(_disc);

	int const* const qsMask = _binding[type]->reactionQuasiStationarity();
	const linalg::ConstMaskArray mask{ qsMask, static_cast<int>(_disc.strideBound[type]) };
	const int probSize = linalg::numMaskActive(mask);

	if (probSize == 0)
		return;

	LinearBufferAllocator tlmAlloc = threadLocalMem.get();

	BufferedArray<double> rhsBuffer = tlmAlloc.array<double>(probSize);
	double* const rhs = static_cast<double*>(rhsBuffer);

	BufferedArray<double> rhsUnmaskedBuffer = tlmAlloc.array<double>(idxr.strideParBound(type));
	double* const rhsUnmasked = static_cast<double*>(rhsUnmaskedBuffer);

	BufferedArray<double> maskedMultiplierBuffer = tlmAlloc.array<double>(idxr.strideParNode(type));
	double* const maskedMultiplier = static_cast<double*>(maskedMultiplierBuffer);

	BufferedArray<double> scaleFactorBuffer = tlmAlloc.array<double>(probSize);
	double* const scaleFactors = static_cast<double*>(scaleFactorBuffer);

	// The matrix is factorized, so it has to own its pivot array
	linalg::DenseMatrix jacobianMatrix;
	jacobianMatrix.resize(probSize, probSize);

	for (unsigned int shell = 0; shell < _disc.nParPoints[type]; ++shell)
	{
		const int shellOffset = idxr.offsetCp(ParticleTypeIndex{ type }) + static_cast<int>(shell) * idxr.strideParNode(type);
		const int qOffset = shellOffset + idxr.strideParLiquid();

		// Extract the quasi-stationary subproblem from the full Jacobian
		jacobianMatrix.setAll(0.0);
		linalg::copyMatrixSubset(_globalJac, mask, mask, qOffset, qOffset, jacobianMatrix);

		// Construct right hand side, which holds -dF / dp - dF / dy * s at this point
		linalg::selectVectorSubset(sensYdot + qOffset, mask, rhs);

		// Zero out the unknown (masked) elements, the remainder of the shell is already known
		std::copy_n(sensY + shellOffset, idxr.strideParNode(type), maskedMultiplier);
		linalg::fillVectorSubset(maskedMultiplier + _disc.nComp, mask, 0.0);

		// Subtract the contribution of the known part of the shell from the right hand side
		Eigen::Map<Eigen::VectorXd> maskedMultiplier_eigen(maskedMultiplier, idxr.strideParNode(type));
		Eigen::Map<Eigen::VectorXd> rhsUnmasked_eigen(rhsUnmasked, idxr.strideParBound(type));
		rhsUnmasked_eigen = _globalJac.block(qOffset, shellOffset, idxr.strideParBound(type), idxr.strideParNode(type)) * maskedMultiplier_eigen;
		linalg::vectorSubsetAdd(rhsUnmasked, mask, -1.0, 1.0, rhs);

		// Precondition
		jacobianMatrix.rowScaleFactors(scaleFactors);
		jacobianMatrix.scaleRows(scaleFactors);

		jacobianMatrix.factorize();
		jacobianMatrix.solve(scaleFactors, rhs);

		// Write back
		linalg::applyVectorSubset(rhs, mask, sensY + qOffset);
	}
}

void FiniteBath::consistentInitialSensitivityBulkTimeDerivative()
{
	Indexer idxr(_disc);

	// Rows are copied from the system Jacobian later on, which requires both matrices to have the
	// same sparsity pattern, see consistentInitialTimeDerivative()
	_globalJacDisc = _globalJac;

	double* vPtr = _globalJacDisc.valuePtr();
	for (int k = 0; k < _globalJacDisc.nonZeros(); ++k)
		vPtr[k] = 0.0;

	linalg::BandedEigenSparseRowIterator jac(_globalJacDisc, idxr.offsetC());
	for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++jac)
		jac[0] += 1.0;
}

void FiniteBath::consistentInitialSensitivityBindingTimeDerivative(double* const sensYdot, unsigned int type)
{
	Indexer idxr(_disc);

	linalg::BandedEigenSparseRowIterator jacPar(_globalJacDisc, idxr.offsetCp(ParticleTypeIndex{ type }));

	for (unsigned int j = 0; j < _disc.nParPoints[type]; ++j)
	{
		// Populate matrix with time derivative Jacobian first
		addTimeDerivativeToJacobianParticleShell(jacPar, idxr, 1.0, type);
		// Iterator jacPar has already been advanced to next shell

		// Overwrite rows corresponding to algebraic equations with the Jacobian and set right hand side to 0
		if (!_binding[type]->hasQuasiStationaryReactions())
			continue;

		const int shellOffset = idxr.offsetCp(ParticleTypeIndex{ type }) + static_cast<int>(j) * idxr.strideParNode(type) + idxr.strideParLiquid();

		linalg::BandedEigenSparseRowIterator jacSolidOrig(_globalJac, shellOffset);
		linalg::BandedEigenSparseRowIterator jacSolid(_globalJacDisc, shellOffset);

		int const* const mask = _binding[type]->reactionQuasiStationarity();
		double* const qShellDot = sensYdot + shellOffset;

		for (int i = 0; i < idxr.strideParBound(type); ++i, ++jacSolid, ++jacSolidOrig)
		{
			if (!mask[i])
				continue;

			jacSolid.copyRowFrom(jacSolidOrig);

			// Right hand side is -\frac{\partial^2 res(t, y, \dot{y})}{\partial p \partial t}
			// If the residual is not explicitly depending on time, this expression is 0
			qShellDot[i] = 0.0;
		}
	}
}

void FiniteBath::solveBulkTimeDerivativeSystem(double* const rhs)
{
	// The bulk time derivative Jacobian of a finite bath is the identity, so the system is solved
	// by simply leaving the right hand side as it is.
}

void FiniteBath::leanConsistentInitialSensitivity(const SimulationTime& simTime, const ConstSimulationState& simState,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerConsistentInit);

	Indexer idxr(_disc);

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		_particles[parType]->leanConsistentInitialStateValidity();

	for (std::size_t param = 0; param < vecSensY.size(); ++param)
	{
		double const* const sensY = vecSensY[param];
		double* const sensYdot = vecSensYdot[param];

		// Copy parameter derivative from AD to tempState and negate it
		for (int i = 0; i < idxr.offsetCp(); ++i)
			_tempState[i] = -adRes[i].getADValue(param);

		std::fill(_tempState + idxr.offsetCp(), _tempState + numDofs(), 0.0);

		// Compute right hand side by adding -dF / dy * s = -J * s to -dF / dp which is already stored in _tempState
		multiplyWithJacobian(simTime, simState, sensY, -1.0, 1.0, _tempState);

		std::copy(_tempState + idxr.offsetC(), _tempState + idxr.offsetCp(), sensYdot + idxr.offsetC());

		solveBulkTimeDerivativeSystem(sensYdot + idxr.offsetC());
	}
}

}  // namespace model

}  // namespace cadet
