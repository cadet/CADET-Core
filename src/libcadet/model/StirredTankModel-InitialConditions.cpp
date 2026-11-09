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

#include "model/StirredTankModel.hpp"
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

int StirredTankModel::multiplexInitialConditions(const cadet::ParameterId& pId, unsigned int adDirection, double adValue)
{
	// INIT_CP is not a parameter of a tank whose particles own no liquid phase, see configure()
	if ((_disc.strideParLiquid == 0) && (pId.name == hashString("INIT_CP")))
		return -1;

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

int StirredTankModel::multiplexInitialConditions(const cadet::ParameterId& pId, double val, bool checkSens)
{
	// INIT_CP is not a parameter of a tank whose particles own no liquid phase, see configure()
	if ((_disc.strideParLiquid == 0) && (pId.name == hashString("INIT_CP")))
		return -1;

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

void StirredTankModel::applyInitialCondition(const SimulationState& simState) const
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

			// Initialize c^p, which a particle in rapid equilibrium with the bulk liquid does not have
			if (idxr.strideParLiquid() > 0)
			{
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
					simState.vecStateY[shellOffset + comp] = static_cast<double>(_initCp[comp + _disc.nComp * type]);
			}

			// Initialize c^s
			active const* const iq = _initCs.data() + _disc.nBoundBeforeType[type];
			for (unsigned int bnd = 0; bnd < _disc.strideBound[type]; ++bnd)
				simState.vecStateY[shellOffset + idxr.strideParLiquid() + bnd] = static_cast<double>(iq[bnd]);
		}
	}

	// Initialize the bulk liquid volume
	simState.vecStateY[idxr.offsetV()] = static_cast<double>(_initLiquidVolume);
}

void StirredTankModel::readInitialCondition(IParameterProvider& paramProvider)
{
	_initState.clear();
	_initStateDot.clear();

	// Check if INIT_STATE is present
	if (paramProvider.exists("INIT_STATE"))
	{
		const std::vector<double> initState = paramProvider.getDoubleArray("INIT_STATE");

		if (initState.size() < numPureDofs())
			throw InvalidParameterException("INIT_STATE does not contain enough values for all degrees of freedom (expected at least " + std::to_string(numPureDofs()) + ", got " + std::to_string(initState.size()) + ")");

		_initState = std::vector<double>(initState.begin(), initState.begin() + numPureDofs());

		// Check if INIT_STATE contains the full state vector and its time derivative
		if (initState.size() >= 2 * numPureDofs())
			_initStateDot = std::vector<double>(initState.begin() + numPureDofs(), initState.begin() + 2 * numPureDofs());
		return;
	}

	const std::vector<double> initC = paramProvider.getDoubleArray("INIT_C");

	if (initC.size() < _disc.nComp)
		throw InvalidParameterException("INIT_C does not contain enough values for all components");

	checkNonNegativeField(initC, "INIT_C");

	ad::copyToAd(initC.data(), _initC.data(), _disc.nComp);

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		std::ostringstream oss;
		oss << "particle_type_" << std::setfill('0') << std::setw(3) << parType;

		if (!paramProvider.exists(oss.str()))
			continue;

		paramProvider.pushScope(oss.str());

		// A particle in rapid equilibrium with the bulk liquid owns no liquid phase. Silently ignoring an
		// initial condition for a state that does not exist would hide the mistake, so it is rejected.
		if (_disc.strideParLiquid == 0)
		{
			if (paramProvider.exists("INIT_CP"))
				throw InvalidParameterException("INIT_CP was specified for particle type " + std::to_string(parType) + ", which is in rapid equilibrium with the bulk liquid and therefore has no liquid phase of its own");
		}
		// Check if INIT_CP is present, otherwise copy from INIT_C
		else if (paramProvider.exists("INIT_CP"))
		{
			const std::vector<double> initCp = paramProvider.getDoubleArray("INIT_CP");

			if (initCp.size() < _disc.nComp)
				throw InvalidParameterException("INIT_CP does not contain enough values for all components");

			checkNonNegativeField(initCp, "INIT_CP");

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

		checkNonNegativeField(initCs, "INIT_CS");

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

	// The liquid volume is read under the name it was specified with, see configure()
	if (paramProvider.exists("INIT_LIQUID_VOLUME"))
		_initLiquidVolume.setValue(paramProvider.getDouble("INIT_LIQUID_VOLUME"));
	else if (paramProvider.exists("LIQUID_VOLUME"))
		_initLiquidVolume.setValue(paramProvider.getDouble("LIQUID_VOLUME"));
}

void StirredTankModel::consistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerConsistentInit);

	// An empty tank has an algebraic rather than a differential bulk balance, which determines the
	// liquid phase concentrations outright
	if (consistentInitialEmptyTank(simTime, vecStateY, errorTol, threadLocalMem))
		return;

	consistentInitialBulkLiquidEquilibrium(simTime, vecStateY, errorTol, threadLocalMem);

	// A particle in rapid equilibrium binds to the bulk liquid, which every such particle type
	// shares, so their bound states cannot be initialized one type at a time against a liquid phase
	// of their own. They are solved jointly with the bulk liquid instead.
	if (_disc.strideParLiquid == 0)
		consistentInitialLiquidEquilibriumBinding(simTime, vecStateY, errorTol, threadLocalMem);
	else
	{
		// Initialize quasi-stationary binding states
		for (unsigned int type = 0; type < _disc.nParType; ++type)
		{
			if (!_binding[type]->hasQuasiStationaryReactions())
				continue;

			consistentInitialBindingEquilibrium(simTime, vecStateY, adJac, errorTol, threadLocalMem, type);
		}
	}

	// The binding equilibrium above writes plain values into the AD state vector, which resets the
	// seed vectors (see sfad::Fwd::operator=), so they have to be restored for the next AD evaluation
	if (adJac.adY)
		prepareADvectors(adJac);
}

bool StirredTankModel::consistentInitialEmptyTank(const SimulationTime& simTime, double* const vecStateY, double errorTol, util::ThreadLocalStorage& threadLocalMem)
{
	Indexer idxr(_disc);

	if (vecStateY[idxr.offsetV()] != 0.0)
		return false;

	// With a vanishing liquid volume the bulk balance loses its time derivative,
	//   \dot{V}^l c_i + V^p [sum_j d_j sum_m \dot{q}_{j,i,m}] = c_{in,i} F_in - c_i F_out,
	// which is algebraic in c. A particle that has pores would contribute a film transfer term here
	// instead, which is a linear system rather than the formula below.
	if (_disc.strideParLiquid != 0)
		throw InvalidParameterException("Consistent initialization of a stirred tank without liquid is only supported for particles in rapid equilibrium with the bulk liquid");

	if (_reaction.conservedMoieties("liquid").hasRapidEquilibrium())
		throw InvalidParameterException("Consistent initialization: reactions in rapid equilibrium require a non-zero liquid volume");

	const double flowIn = static_cast<double>(_flowRateIn);
	const double flowOut = static_cast<double>(_flowRateOut);

	// Volume: \dot{V} = F_in - F_out - F_filter
	const double vDot = flowIn - flowOut - static_cast<double>(_curFlowRateFilter);

	// A vanishing denominator means \dot{V} = F_out = 0, and then F_in = F_filter. The liquid would
	// be drawn off as fast as it enters, leaving the dry components behind, which is assumed not to
	// be intended. The initial conditions are left in place.
	const double denom = vDot + flowOut;
	if (denom == 0.0)
		return true;

	// Without bound states the balance is the explicit formula
	//   c_i = c_{in,i} F_in / (\dot{V}^l + F_out),
	// which is also the starting point of the solve below
	double* const c = vecStateY + idxr.offsetC();

	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		c[comp] = vecStateY[comp] * flowIn / denom;

	const double vSolid = static_cast<double>(_constSolidVolume);
	const unsigned int nTotalBound = (_disc.nParType > 0) ? _disc.strideBound[_disc.nParType] : 0u;

	if ((nTotalBound == 0) || (vSolid == 0.0))
		return true;

	// The bound state rates are what the balance couples the liquid phase to, and a quasi-stationary
	// state does not have one: its rate follows from differentiating its equilibrium condition, which
	// brings in \dot{c} and raises the index of the system. CSTRModel assumed every bound state to be
	// dynamic without saying so.
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		if (_binding[type]->hasQuasiStationaryReactions())
			throw InvalidParameterException("Consistent initialization: quasi-stationary binding in particle type " + std::to_string(type) + " requires a non-zero liquid volume");
	}

	// A dynamic bound state has \dot{q}_{j,i,m} = -f_{j,i,m}(c, q), see IBindingModel::flux(), so
	//   (\dot{V}^l + F_out) c_i - V^p [sum_j d_j sum_m f_{j,i,m}(c, q)] = c_{in,i} F_in
	// determines the liquid phase. The bound states are given, the binding fluxes are not generally
	// linear in c, so this is solved rather than evaluated.
	LinearBufferAllocator tlmAlloc = threadLocalMem.get();
	const ColumnPosition colPos{ 0.0, 0.0, 0.0 };
	const unsigned int probSize = _disc.nComp;

	const unsigned int maxStrideBound = *std::max_element(_disc.strideBound, _disc.strideBound + _disc.nParType);
	const unsigned int blockSize = _disc.nComp + maxStrideBound;

	BufferedArray<double> nonlinMemBuffer = tlmAlloc.array<double>(_nonlinearSolver->workspaceSize(probSize));
	double* const nonlinMem = static_cast<double*>(nonlinMemBuffer);

	// The binding model reaches the liquid phase by counting back from the bound states, so the two
	// are handed to it as one block
	BufferedArray<double> cellBuffer = tlmAlloc.array<double>(blockSize);
	double* const cell = static_cast<double*>(cellBuffer);
	double* const cellQ = cell + _disc.nComp;

	BufferedArray<double> fluxBuffer = tlmAlloc.array<double>(maxStrideBound);
	double* const flux = static_cast<double*>(fluxBuffer);

	BufferedArray<double> fluxJacBuffer = tlmAlloc.array<double>(static_cast<std::size_t>(blockSize) * blockSize);
	linalg::DenseMatrixView fluxJacobian(static_cast<double*>(fluxJacBuffer), nullptr, blockSize, blockSize);

	// The nonlinear solver factorizes this matrix, so it has to own its pivot array
	linalg::DenseMatrix jacobianMatrix;
	jacobianMatrix.resize(probSize, probSize);

	_nonlinearSolver->solve(
		[&](double const* const x, double* const res)
		{
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				res[comp] = denom * x[comp] - flowIn * vecStateY[comp];

			std::copy_n(x, _disc.nComp, cell);

			for (unsigned int type = 0; type < _disc.nParType; ++type)
			{
				const unsigned int strideBound = _disc.strideBound[type];
				if (strideBound == 0)
					continue;

				std::copy_n(vecStateY + idxr.offsetCp(ParticleTypeIndex{ type }), strideBound, cellQ);
				std::fill_n(flux, strideBound, 0.0);

				_binding[type]->flux(simTime.t, simTime.secIdx, colPos, cellQ, cell, flux, tlmAlloc.manageRemainingMemory());

				const double factor = vSolid * static_cast<double>(_parTypeVolFrac[type]);
				unsigned int const* const nBound = _disc.nBound + type * _disc.nComp;

				unsigned int idx = 0;
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				{
					for (unsigned int bnd = 0; bnd < nBound[comp]; ++bnd, ++idx)
						res[comp] -= factor * flux[idx];
				}
			}

			return true;
		},
		[&](double const* const x, linalg::detail::DenseMatrixBase& mat)
		{
			mat.setAll(0.0);
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				mat.native(comp, comp) = denom;

			std::copy_n(x, _disc.nComp, cell);

			for (unsigned int type = 0; type < _disc.nParType; ++type)
			{
				const unsigned int strideBound = _disc.strideBound[type];
				if (strideBound == 0)
					continue;

				std::copy_n(vecStateY + idxr.offsetCp(ParticleTypeIndex{ type }), strideBound, cellQ);
				fluxJacobian.setAll(0.0);

				// offsetCp counts backwards from a bound state row to the liquid phase, see
				// IBindingModel::analyticJacobian
				_binding[type]->analyticJacobian(simTime.t, simTime.secIdx, colPos, cellQ, _disc.nComp, fluxJacobian.row(_disc.nComp), tlmAlloc.manageRemainingMemory());

				const double factor = vSolid * static_cast<double>(_parTypeVolFrac[type]);
				unsigned int const* const nBound = _disc.nBound + type * _disc.nComp;

				unsigned int idx = 0;
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				{
					for (unsigned int bnd = 0; bnd < nBound[comp]; ++bnd, ++idx)
					{
						for (unsigned int toComp = 0; toComp < _disc.nComp; ++toComp)
							mat.native(comp, toComp) -= factor * fluxJacobian.native(_disc.nComp + idx, toComp);
					}
				}
			}

			return true;
		},
		errorTol, c, nonlinMem, jacobianMatrix, probSize);

	return true;
}

void StirredTankModel::consistentInitialLiquidEquilibriumBinding(const SimulationTime& simTime, double* const vecStateY, double errorTol, util::ThreadLocalStorage& threadLocalMem)
{
	Indexer idxr(_disc);

	if (_disc.nParType == 0)
		return;

	// Without particle liquid phases the bulk liquid and the bound states of every type are
	// contiguous, so they form a single masked system
	double* const c = vecStateY + idxr.offsetC();
	const unsigned int nTotalBound = _disc.strideBound[_disc.nParType];
	const int maskLen = static_cast<int>(_disc.nComp + nTotalBound);

	LinearBufferAllocator tlmAlloc = threadLocalMem.get();
	const ColumnPosition colPos{ 0.0, 0.0, 0.0 };

	std::vector<int> qsMask(maskLen, false);
	bool hasQuasiStationary = false;

	int* localMask = qsMask.data() + _disc.nComp;
	for (unsigned int type = 0; type < _disc.nParType; localMask += _disc.strideBound[type], ++type)
	{
		if (!_binding[type]->hasQuasiStationaryReactions())
			continue;

		double* const typeQ = vecStateY + idxr.offsetCp(ParticleTypeIndex{ type });

		if (!_binding[type]->preConsistentInitialState(simTime.t, simTime.secIdx, colPos, typeQ, c, tlmAlloc))
			continue;

		hasQuasiStationary = true;

		int const* const qsMaskSrc = _binding[type]->reactionQuasiStationarity();
		std::copy_n(qsMaskSrc, _disc.strideBound[type], localMask);

		// A component is part of the system as soon as one of its bound states is quasi-stationary
		unsigned int idx = 0;
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			const unsigned int nBoundComp = _disc.nBound[type * _disc.nComp + comp];

			for (unsigned int state = 0; state < nBoundComp; ++state)
			{
				if (qsMaskSrc[idx + state])
				{
					qsMask[comp] = true;
					break;
				}
			}

			idx += nBoundComp;
		}
	}

	if (!hasQuasiStationary)
		return;

	const linalg::ConstMaskArray mask{ qsMask.data(), maskLen };
	const int probSize = linalg::numMaskActive(mask);
	const unsigned int numActiveComp = linalg::numMaskActive(mask, _disc.nComp);

	// The conserved quantity of a component is the amount of it held by the tank, which the
	// equilibrium redistributes between the bulk liquid and the bound states but does not change
	const double vLiquid = vecStateY[idxr.offsetV()];
	const double vSolid = static_cast<double>(_constSolidVolume);
	const double porosity = vLiquid / (vLiquid + vSolid);
	const double epsQ = 1.0 - porosity;

	BufferedArray<double> nonlinMemBuffer = tlmAlloc.array<double>(_nonlinearSolver->workspaceSize(probSize));
	double* const nonlinMem = static_cast<double*>(nonlinMemBuffer);

	BufferedArray<double> solutionBuffer = tlmAlloc.array<double>(probSize);
	double* const solution = static_cast<double*>(solutionBuffer);

	BufferedArray<double> fullXBuffer = tlmAlloc.array<double>(maskLen);
	double* const fullX = static_cast<double*>(fullXBuffer);

	BufferedArray<double> fullResidualBuffer = tlmAlloc.array<double>(maskLen);
	double* const fullResidual = static_cast<double*>(fullResidualBuffer);

	BufferedArray<double> fullJacobianBuffer = tlmAlloc.array<double>(static_cast<std::size_t>(maskLen) * maskLen);
	linalg::DenseMatrixView fullJacobianMatrix(static_cast<double*>(fullJacobianBuffer), nullptr, maskLen, maskLen);

	BufferedArray<double> conservedQuantsBuffer = tlmAlloc.array<double>(std::max(numActiveComp, 1u));
	double* const conservedQuants = static_cast<double*>(conservedQuantsBuffer);

	// Extract initial values from current state
	linalg::selectVectorSubset(c, mask, solution);

	{
		unsigned int rIdx = 0;
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			if (!mask.mask[comp])
				continue;

			conservedQuants[rIdx] = porosity * c[comp];

			for (unsigned int type = 0; type < _disc.nParType; ++type)
			{
				const double factor = epsQ * static_cast<double>(_parTypeVolFrac[type]);
				const unsigned int offset = _disc.nComp + _disc.parTypeOffset[type] + _disc.boundOffset[type * _disc.nComp + comp];

				for (unsigned int state = 0; state < _disc.nBound[type * _disc.nComp + comp]; ++state)
				{
					if (!mask.mask[offset + state])
						continue;

					conservedQuants[rIdx] += factor * c[offset + state];
				}
			}

			++rIdx;
		}
	}

	// Evaluates the binding flux of every particle type against the shared bulk liquid, and
	// optionally its Jacobian, into the unmasked scratch arrays
	auto evaluateBindingFluxes = [&](bool wantJac)
	{
		std::fill_n(fullResidual, maskLen, 0.0);
		if (wantJac)
			fullJacobianMatrix.setAll(0.0);

		for (unsigned int type = 0; type < _disc.nParType; ++type)
		{
			const int rowOffset = static_cast<int>(_disc.nComp + _disc.parTypeOffset[type]);
			double* const typeQ = fullX + rowOffset;

			_binding[type]->flux(simTime.t, simTime.secIdx, colPos, typeQ, fullX, fullResidual + rowOffset, tlmAlloc.manageRemainingMemory());

			if (wantJac)
			{
				// offsetCp counts backwards from a bound state row to the liquid phase, see
				// IBindingModel::analyticJacobian
				_binding[type]->analyticJacobian(simTime.t, simTime.secIdx, colPos, typeQ, rowOffset, fullJacobianMatrix.row(rowOffset), tlmAlloc.manageRemainingMemory());
			}
		}
	};

	// Replaces the component rows of the masked Jacobian by the conservation relations
	auto applyConservationRows = [&](linalg::detail::DenseMatrixBase& mat)
	{
		mat.submatrixSetAll(0.0, 0, 0, numActiveComp, probSize);

		unsigned int rIdx = 0;
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			if (!mask.mask[comp])
				continue;

			mat.native(rIdx, rIdx) = porosity;

			for (unsigned int type = 0; type < _disc.nParType; ++type)
			{
				const double factor = epsQ * static_cast<double>(_parTypeVolFrac[type]);
				const unsigned int offset = _disc.nComp + _disc.parTypeOffset[type] + _disc.boundOffset[type * _disc.nComp + comp];
				unsigned int bIdx = numActiveComp + linalg::numMaskActive(mask, _disc.nComp, _disc.parTypeOffset[type] + _disc.boundOffset[type * _disc.nComp + comp]);

				for (unsigned int state = 0; state < _disc.nBound[type * _disc.nComp + comp]; ++state)
				{
					if (!mask.mask[offset + state])
						continue;

					mat.native(rIdx, bIdx) = factor;
					++bIdx;
				}
			}

			++rIdx;
		}
	};

	// The nonlinear solver factorizes this matrix, so it has to own its pivot array
	linalg::DenseMatrix jacobianMatrix;
	jacobianMatrix.resize(probSize, probSize);

	_nonlinearSolver->solve(
		[&](double const* const x, double* const res)
		{
			std::copy_n(c, maskLen, fullX);
			linalg::applyVectorSubset(x, mask, fullX);

			evaluateBindingFluxes(false);

			linalg::selectVectorSubset(fullResidual, mask, res);

			// Replace the component rows by the conservation relations
			unsigned int rIdx = 0;
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				if (!mask.mask[comp])
					continue;

				res[rIdx] = porosity * x[rIdx] - conservedQuants[rIdx];

				for (unsigned int type = 0; type < _disc.nParType; ++type)
				{
					const double factor = epsQ * static_cast<double>(_parTypeVolFrac[type]);
					const unsigned int offset = _disc.nComp + _disc.parTypeOffset[type] + _disc.boundOffset[type * _disc.nComp + comp];
					unsigned int bIdx = numActiveComp + linalg::numMaskActive(mask, _disc.nComp, _disc.parTypeOffset[type] + _disc.boundOffset[type * _disc.nComp + comp]);

					for (unsigned int state = 0; state < _disc.nBound[type * _disc.nComp + comp]; ++state)
					{
						if (!mask.mask[offset + state])
							continue;

						res[rIdx] += factor * x[bIdx];
						++bIdx;
					}
				}

				++rIdx;
			}

			return true;
		},
		[&](double const* const x, linalg::detail::DenseMatrixBase& mat)
		{
			std::copy_n(c, maskLen, fullX);
			linalg::applyVectorSubset(x, mask, fullX);

			evaluateBindingFluxes(true);

			mat.setAll(0.0);
			linalg::copyMatrixSubset(fullJacobianMatrix, mask, mask, mat);

			applyConservationRows(mat);

			return true;
		},
		errorTol, solution, nonlinMem, jacobianMatrix, probSize);

	// Apply solution
	linalg::applyVectorSubset(solution, mask, c);

	// Refine / correct solution
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		if (!_binding[type]->hasQuasiStationaryReactions())
			continue;

		double* const typeQ = vecStateY + idxr.offsetCp(ParticleTypeIndex{ type });
		_binding[type]->postConsistentInitialState(simTime.t, simTime.secIdx, colPos, typeQ, c, tlmAlloc);
	}
}

void StirredTankModel::consistentInitialBulkLiquidEquilibrium(const SimulationTime& simTime, double* const vecStateY, double errorTol, util::ThreadLocalStorage& threadLocalMem)
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
		throw InvalidParameterException("CSTR consistent initialization: Invalid equilibrium initialization size");

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

void StirredTankModel::consistentInitialBindingEquilibrium(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, util::ThreadLocalStorage& threadLocalMem, unsigned int type)
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

void StirredTankModel::consistentInitialTimeDerivative(const SimulationTime& simTime, double const* vecStateY, double* const vecStateYdot, util::ThreadLocalStorage& threadLocalMem)
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

	// An empty tank whose inflow is drawn off by the filter as fast as it enters has no rate to
	// determine, see consistentInitialEmptyTank(). The negated residual vanishes along with the
	// flow rates, so it is already the time derivative that is left in place here.
	if ((vecStateY[idxr.offsetV()] == 0.0) && (static_cast<double>(_flowRateIn) == static_cast<double>(_curFlowRateFilter)))
		return;

	consistentInitialBulkTimeDerivative(simTime, vecStateY, vecStateYdot, threadLocalMem);

	for (unsigned int type = 0; type < _disc.nParType; ++type)
		consistentInitialBindingTimeDerivative(simTime, vecStateYdot, threadLocalMem, type);

	_linearSolver->factorize(_globalJacDisc, idxr.offsetC(), numPureDofs());
	if (cadet_unlikely(_linearSolver->info() != Eigen::Success))
	{
		LOG(Error) << "Factorize() failed";
	}

	_linearSolver->solve(yDot.data() + idxr.offsetC(), numPureDofs());
	if (cadet_unlikely(_linearSolver->info() != Eigen::Success))
	{
		LOG(Error) << "Solve() failed";
	}
}

void StirredTankModel::consistentInitialBulkTimeDerivative(const SimulationTime& simTime, double const* vecStateY, double* const vecStateYdot, util::ThreadLocalStorage& threadLocalMem)
{
	Indexer idxr(_disc);

	// An empty tank has an algebraic rather than a differential bulk balance, which is
	// differentiated once below to obtain an equation for \dot{c}
	const double vLiquid = vecStateY[idxr.offsetV()];

	const auto& cm = _reaction.conservedMoieties("liquid");
	if (cm.isEnabled() && cm.numEquilibriumReactions() > 0)
	{
		if (vLiquid == 0.0)
			throw InvalidParameterException("Consistent initialization: reactions in rapid equilibrium require a non-zero liquid volume");

		const auto& L = cm.conservedMoietyMatrix();
		const unsigned int nMoieties = cm.numMoieties();
		const unsigned int nEq = cm.numEquilibriumReactions();

		LinearBufferAllocator tlmAlloc = threadLocalMem.get();

		double const* const cLocal = vecStateY + idxr.offsetC();
		double* const cDotLocal = vecStateYdot + idxr.offsetC();

		const ColumnPosition colPos{ 0.0, 0.0, 0.0 };

		// Differential rows: dF / dyDot = V^l * L, plus the liquid volume column L * c
		const double v = vecStateY[idxr.offsetV()];

		for (unsigned int m = 0; m < nMoieties; ++m)
		{
			double volumeEntry = 0.0;
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				_globalJacDisc.coeffRef(idxr.offsetC() + m, idxr.offsetC() + comp) = v * L(m, comp);
				volumeEntry += L(m, comp) * cLocal[comp];
			}

			_globalJacDisc.coeffRef(idxr.offsetC() + m, idxr.offsetV()) = volumeEntry;
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
	else if (vLiquid == 0.0)
	{
		// Without liquid the bulk balance reads
		//   \dot{V}^l c_i + V^p [sum_j d_j sum_m \dot{q}_{j,i,m}] = c_{in,i} F_in - c_i F_out
		// and holds no \dot{c}_i at all, see consistentInitialEmptyTank(). The residual therefore
		// determines no liquid phase rate and
		//   (\dot{V}^l + F_out) \dot{c}_i + V^p [sum_j d_j sum_m \dot{q}_{j,i,m}] = 0
		// is a choice rather than a consequence: it is what differentiating the balance gives if the
		// bound states are taken as frozen over the step. A rate that follows from the balance needs
		// \ddot{q}, which brings in \dot{c} through the binding model and raises the index again.
		// The bound state terms of the equation are contributed by the rapid equilibrium block below.
		if (_disc.strideParLiquid != 0)
			throw InvalidParameterException("Consistent initialization of a stirred tank without liquid is only supported for particles in rapid equilibrium with the bulk liquid");

		// \dot{V}^l + F_out = F_in - F_filter
		const double denom = static_cast<double>(_flowRateIn) - static_cast<double>(_curFlowRateFilter);
		double* const cDotLocal = vecStateYdot + idxr.offsetC();

		linalg::BandedEigenSparseRowIterator jacBlk(_globalJacDisc, idxr.offsetC());
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++jacBlk)
		{
			jacBlk[0] = denom;

			// The inlet concentration rate is not available here, which violates the assumption
			// that every outlet DOF is dynamic. Repairing this requires a different consistent
			// initialization concept, so \dot{c}_{in} is taken as zero.
			cDotLocal[comp] = 0.0;
		}
	}
	else
	{
		// The bulk balance is extensive, so dF / dyDot holds V^l on the diagonal and c_i in the
		// liquid volume column
		double const* const cLocal = vecStateY + idxr.offsetC();
		const int volumeOffset = idxr.offsetV() - idxr.offsetC();

		linalg::BandedEigenSparseRowIterator jacBlk(_globalJacDisc, idxr.offsetC());
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++jacBlk)
		{
			jacBlk[0] = vLiquid;
			jacBlk[volumeOffset - static_cast<int>(comp)] = cLocal[comp];
		}
	}

	// What the bound states of a particle in rapid equilibrium take up leaves the bulk liquid, so
	// their time derivatives appear in the bulk rows weighted with the particle volume
	if (_disc.strideParLiquid == 0)
	{
		const double vSolid = static_cast<double>(_constSolidVolume);

		for (unsigned int type = 0; type < _disc.nParType; ++type)
		{
			const double factor = vSolid * static_cast<double>(_parTypeVolFrac[type]);
			const int offsetQ = idxr.offsetCp(ParticleTypeIndex{ type });
			unsigned int const* const nBound = _disc.nBound + type * _disc.nComp;
			unsigned int const* const boundOffset = _disc.boundOffset + type * _disc.nComp;

			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				for (unsigned int bnd = 0; bnd < nBound[comp]; ++bnd)
					_globalJacDisc.coeffRef(idxr.offsetC() + comp, offsetQ + boundOffset[comp] + bnd) = factor;
			}
		}
	}

	// Volume: dF / dyDot = 1
	_globalJacDisc.coeffRef(idxr.offsetV(), idxr.offsetV()) = 1.0;
}

void StirredTankModel::consistentInitialBindingTimeDerivative(const SimulationTime& simTime, double* const vecStateYdot, util::ThreadLocalStorage& threadLocalMem, unsigned int type)
{
	Indexer idxr(_disc);

	linalg::BandedEigenSparseRowIterator jacPar(_globalJacDisc, idxr.offsetCp(ParticleTypeIndex{ type }));

	LinearBufferAllocator tlmAlloc = threadLocalMem.get();
	double* const dFluxDt = _tempState + idxr.offsetCp(ParticleTypeIndex{ type });

	// A particle in rapid equilibrium holds bound states only, and the liquid phase its fluxes depend
	// on is the bulk liquid rather than one of its own
	if (_disc.strideParLiquid == 0)
	{
		const int offsetQ = idxr.offsetCp(ParticleTypeIndex{ type });
		int const* const mask = _binding[type]->reactionQuasiStationarity();
		double* const qDot = vecStateYdot + offsetQ;

		// dF / dyDot of the bound states, which a quasi-stationary state does not have
		linalg::BandedEigenSparseRowIterator jacQ(_globalJacDisc, offsetQ);
		for (int i = 0; i < idxr.strideParBound(type); ++i, ++jacQ)
		{
			if (!mask[i])
				jacQ[0] += 1.0;
		}

		if (!_binding[type]->hasQuasiStationaryReactions())
			return;

		// Obtain derivative of fluxes wrt. time
		std::fill_n(dFluxDt, _disc.strideBound[type], 0.0);
		if (_binding[type]->dependsOnTime())
		{
			_binding[type]->timeDerivativeQuasiStationaryFluxes(simTime.t, simTime.secIdx,
				ColumnPosition{ 0.0, 0.0, _particles[type]->relativeCoordinate(0) },
				vecStateYdot + idxr.offsetC(), qDot, dFluxDt, tlmAlloc);
		}

		// Copy row from original Jacobian (without time derivatives) and set right hand side
		linalg::BandedEigenSparseRowIterator jacSolidOrig(_globalJac, offsetQ);
		linalg::BandedEigenSparseRowIterator jacSolid(_globalJacDisc, offsetQ);

		for (int i = 0; i < idxr.strideParBound(type); ++i, ++jacSolid, ++jacSolidOrig)
		{
			if (!mask[i])
				continue;

			jacSolid.copyRowFrom(jacSolidOrig);
			qDot[i] = -dFluxDt[i];
		}

		return;
	}

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

void StirredTankModel::leanConsistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerConsistentInit);

	// The bulk liquid equilibrium equations are algebraic, so they hold no time derivative that
	// leanConsistentInitialTimeDerivative() could use to absorb an inconsistent liquid phase.
	// Leaving it in place, as done below, therefore keeps a residual that the relaxed
	// initialization cannot repair later on.
	if (_reaction.conservedMoieties("liquid").hasRapidEquilibrium())
		throw InvalidParameterException("CSTR lean consistent initialization: reactions in rapid equilibrium are not supported yet");

	// An empty tank has an algebraic rather than a differential bulk balance, which determines the
	// liquid phase concentrations outright. The bound states are assumed to be approximately
	// correct already, which is what the relaxed initialization builds on.
	if (consistentInitialEmptyTank(simTime, vecStateY, errorTol, threadLocalMem))
		return;

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		_particles[parType]->leanConsistentInitialStateValidity();

	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		if (!_binding[type]->hasQuasiStationaryReactions())
			continue;

		leanConsistentInitialBindingEquilibrium(simTime, vecStateY, threadLocalMem, type);
	}
}

void StirredTankModel::leanConsistentInitialBindingEquilibrium(const SimulationTime& simTime, double* const vecStateY, util::ThreadLocalStorage& threadLocalMem, unsigned int type)
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

void StirredTankModel::leanConsistentInitialTimeDerivative(double t, double const* const vecStateY, double* const vecStateYdot, double* const res, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerConsistentInit);

	Indexer idxr(_disc);

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		_particles[parType]->leanConsistentInitialTimeDerivativeValidity();

	double* const resSlice = res + idxr.offsetC();
	const int volumeIdx = idxr.offsetV() - idxr.offsetC();

	// Note that the residual is not negated as required at this point, so negate the entries that
	// the bulk solve reads. The bound state rates are left where the relaxed initialization found
	// them, which is what makes it relaxed.
	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		resSlice[comp] = -resSlice[comp];

	resSlice[volumeIdx] = -resSlice[volumeIdx];

	if (vecStateY[idxr.offsetV()] == 0.0)
	{
		leanConsistentInitialEmptyTankTimeDerivative(vecStateYdot);
		vecStateYdot[idxr.offsetV()] = resSlice[volumeIdx];
		return;
	}

	solveBulkTimeDerivativeSystem(resSlice, vecStateY, vecStateYdot);

	std::copy_n(resSlice, _disc.nComp, vecStateYdot + idxr.offsetC());
	vecStateYdot[idxr.offsetV()] = resSlice[volumeIdx];
}

void StirredTankModel::leanConsistentInitialEmptyTankTimeDerivative(double* const vecStateYdot)
{
	Indexer idxr(_disc);

	// Without liquid the bulk balance holds no \dot{c}_i, and the rate
	//   \dot{c}_i = -V^p [sum_j d_j sum_m \dot{q}_{j,i,m}] / (\dot{V}^l + F_out)
	// is the choice consistentInitialBulkTimeDerivative() makes, repeated here so that the relaxed
	// initialization agrees with the full one. A particle that has pores would contribute a film
	// transfer term here instead, which is a linear system rather than the explicit formula below.
	if (_disc.strideParLiquid != 0)
		throw InvalidParameterException("Consistent initialization of a stirred tank without liquid is only supported for particles in rapid equilibrium with the bulk liquid");

	// A vanishing denominator leaves the time derivatives in place, see consistentInitialEmptyTank()
	const double denom = static_cast<double>(_flowRateIn) - static_cast<double>(_curFlowRateFilter);
	if (denom == 0.0)
		return;

	const double qFactor = static_cast<double>(_constSolidVolume) / denom;
	double* const cDot = vecStateYdot + idxr.offsetC();

	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
	{
		// The inlet concentration rate is not available here, see consistentInitialBulkTimeDerivative()
		cDot[comp] = 0.0;
	}

	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		unsigned int const* const boundOffset = _disc.boundOffset + type * _disc.nComp;
		unsigned int const* const nBound = _disc.nBound + type * _disc.nComp;
		double const* const typeQdot = vecStateYdot + idxr.offsetCp(ParticleTypeIndex{ type });

		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			double qDotSum = 0.0;
			double const* const localQdot = typeQdot + boundOffset[comp];
			for (unsigned int bnd = 0; bnd < nBound[comp]; ++bnd)
				qDotSum += localQdot[bnd];

			cDot[comp] -= qFactor * static_cast<double>(_parTypeVolFrac[type]) * qDotSum;
		}
	}
}

void StirredTankModel::initializeSensitivityStates(const std::vector<double*>& vecSensY) const
{
	for (std::size_t param = 0; param < vecSensY.size(); ++param)
	{
		initializeSensitivityBulkStates(vecSensY[param], param);

		for (unsigned int type = 0; type < _disc.nParType; ++type)
			initializeSensitivityParticleStates(vecSensY[param], param, type);
	}
}

void StirredTankModel::initializeSensitivityBulkStates(double* const sensY, std::size_t param) const
{
	Indexer idxr(_disc);
	double* const stateYbulk = sensY + idxr.offsetC();

	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		stateYbulk[comp] = _initC[comp].getADValue(param);

	// The initial liquid volume is a parameter, so a sensitivity with respect to it starts at one
	sensY[idxr.offsetV()] = _initLiquidVolume.getADValue(param);
}

void StirredTankModel::initializeSensitivityParticleStates(double* const sensY, std::size_t param, unsigned int type) const
{
	Indexer idxr(_disc);

	const unsigned int offset = idxr.offsetCp(ParticleTypeIndex{ type });

	for (unsigned int shell = 0; shell < _disc.nParPoints[type]; ++shell)
	{
		const unsigned int shellOffset = offset + shell * idxr.strideParNode(type);
		double* const stateYparticle = sensY + shellOffset;
		double* const stateYparticleSolid = stateYparticle + idxr.strideParLiquid();

		// Initialize c^p, which a particle in rapid equilibrium with the bulk liquid does not have
		if (idxr.strideParLiquid() > 0)
		{
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				stateYparticle[comp] = _initCp[comp + type * _disc.nComp].getADValue(param);
		}

		// Initialize c^s
		for (unsigned int bnd = 0; bnd < _disc.strideBound[type]; ++bnd)
			stateYparticleSolid[bnd] = _initCs[bnd + _disc.nBoundBeforeType[type]].getADValue(param);
	}
}

void StirredTankModel::consistentInitialSensitivity(const SimulationTime& simTime, const ConstSimulationState& simState,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerConsistentInit);

	// The bulk liquid equilibrium equations are algebraic, so dF / dyDot has nEq zero rows and the
	// state sensitivity has to satisfy the linearized equilibrium constraint. Neither is handled
	// here: step 1 below only solves the quasi-stationary binding equilibria, and
	// consistentInitialSensitivityBulkTimeDerivative() in step 2 puts V^l on every bulk liquid row
	// instead of the conserved moiety matrix.
	if (_reaction.conservedMoieties("liquid").hasRapidEquilibrium())
		throw InvalidParameterException("CSTR consistent sensitivity initialization: reactions in rapid equilibrium are not supported yet");

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

		consistentInitialSensitivityBulkTimeDerivative(simState.vecStateY);

		for (unsigned int type = 0; type < _disc.nParType; ++type)
			consistentInitialSensitivityBindingTimeDerivative(sensYdot, type);

		Eigen::Map<Eigen::VectorXd> yDot(sensYdot + idxr.offsetC(), numPureDofs());

		_linearSolver->factorize(_globalJacDisc, idxr.offsetC(), numPureDofs());

		if (cadet_unlikely(_linearSolver->info() != Eigen::Success))
		{
			LOG(Error) << "Factorize() failed";
		}

		_linearSolver->solve(yDot.data(), numPureDofs());

		if (cadet_unlikely(_linearSolver->info() != Eigen::Success))
		{
			LOG(Error) << "Solve() failed";
		}
	}
}

void StirredTankModel::consistentInitialSensitivityBindingEquilibrium(double* const sensY, double const* const sensYdot, util::ThreadLocalStorage& threadLocalMem, unsigned int type)
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

	// A particle in rapid equilibrium binds to the bulk liquid, so the known part of its subproblem
	// reaches from the bulk liquid to the liquid volume rather than covering a particle node
	const int knownStart = (_disc.strideParLiquid == 0) ? idxr.offsetC() : 0;
	const int knownLen = (_disc.strideParLiquid == 0) ? (idxr.offsetV() + 1 - idxr.offsetC()) : idxr.strideParNode(type);

	BufferedArray<double> maskedMultiplierBuffer = tlmAlloc.array<double>(knownLen);
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

		// Zero out the unknown (masked) elements, the remainder is already known
		const int blockStart = (_disc.strideParLiquid == 0) ? knownStart : shellOffset;

		std::copy_n(sensY + blockStart, knownLen, maskedMultiplier);
		linalg::fillVectorSubset(maskedMultiplier + (qOffset - blockStart), mask, 0.0);

		// Subtract the contribution of the known part from the right hand side
		Eigen::Map<Eigen::VectorXd> maskedMultiplier_eigen(maskedMultiplier, knownLen);
		Eigen::Map<Eigen::VectorXd> rhsUnmasked_eigen(rhsUnmasked, idxr.strideParBound(type));
		rhsUnmasked_eigen = _globalJac.block(qOffset, blockStart, idxr.strideParBound(type), knownLen) * maskedMultiplier_eigen;
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

void StirredTankModel::consistentInitialSensitivityBulkTimeDerivative(double const* const y)
{
	Indexer idxr(_disc);

	// Rows are copied from the system Jacobian later on, which requires both matrices to have the
	// same sparsity pattern, see consistentInitialTimeDerivative()
	_globalJacDisc = _globalJac;

	double* vPtr = _globalJacDisc.valuePtr();
	for (int k = 0; k < _globalJacDisc.nonZeros(); ++k)
		vPtr[k] = 0.0;

	// The bulk balance is extensive, so dF / dyDot holds V^l on the diagonal and c_i in the liquid
	// volume column
	double const* const c = y + idxr.offsetC();
	const double v = y[idxr.offsetV()];
	const int volumeOffset = idxr.offsetV() - idxr.offsetC();

	linalg::BandedEigenSparseRowIterator jac(_globalJacDisc, idxr.offsetC());
	for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++jac)
	{
		jac[0] += v;
		jac[volumeOffset - static_cast<int>(comp)] += c[comp];
	}

	// What the bound states of a particle in rapid equilibrium take up leaves the bulk liquid, so
	// their time derivatives appear in the bulk rows weighted with the particle volume, see
	// consistentInitialBulkTimeDerivative()
	if (_disc.strideParLiquid == 0)
	{
		const double vSolid = static_cast<double>(_constSolidVolume);

		for (unsigned int type = 0; type < _disc.nParType; ++type)
		{
			const double factor = vSolid * static_cast<double>(_parTypeVolFrac[type]);
			const int offsetQ = idxr.offsetCp(ParticleTypeIndex{ type });
			unsigned int const* const nBound = _disc.nBound + type * _disc.nComp;
			unsigned int const* const boundOffset = _disc.boundOffset + type * _disc.nComp;

			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				for (unsigned int bnd = 0; bnd < nBound[comp]; ++bnd)
					_globalJacDisc.coeffRef(idxr.offsetC() + comp, offsetQ + boundOffset[comp] + bnd) += factor;
			}
		}
	}

	// Volume: dF / dyDot = 1
	_globalJacDisc.coeffRef(idxr.offsetV(), idxr.offsetV()) += 1.0;
}

void StirredTankModel::consistentInitialSensitivityBindingTimeDerivative(double* const sensYdot, unsigned int type)
{
	Indexer idxr(_disc);

	linalg::BandedEigenSparseRowIterator jacPar(_globalJacDisc, idxr.offsetCp(ParticleTypeIndex{ type }));

	// A particle in rapid equilibrium holds bound states only, see consistentInitialBindingTimeDerivative()
	if (_disc.strideParLiquid == 0)
	{
		const int offsetQ = idxr.offsetCp(ParticleTypeIndex{ type });
		int const* const mask = _binding[type]->reactionQuasiStationarity();

		linalg::BandedEigenSparseRowIterator jacQ(_globalJacDisc, offsetQ);
		for (int i = 0; i < idxr.strideParBound(type); ++i, ++jacQ)
		{
			if (!mask[i])
				jacQ[0] += 1.0;
		}

		if (!_binding[type]->hasQuasiStationaryReactions())
			return;

		linalg::BandedEigenSparseRowIterator jacSolidOrig(_globalJac, offsetQ);
		linalg::BandedEigenSparseRowIterator jacSolid(_globalJacDisc, offsetQ);
		double* const qDot = sensYdot + offsetQ;

		for (int i = 0; i < idxr.strideParBound(type); ++i, ++jacSolid, ++jacSolidOrig)
		{
			if (!mask[i])
				continue;

			jacSolid.copyRowFrom(jacSolidOrig);
			qDot[i] = 0.0;
		}

		return;
	}

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

void StirredTankModel::solveBulkTimeDerivativeSystem(double* const rhs, double const* const y, double const* const yDot)
{
	Indexer idxr(_disc);

	// The bulk balance is extensive, so its time derivative Jacobian on the bulk concentrations and
	// the liquid volume is
	//   [ V^l I   c ]
	//   [   0    1  ]
	// which is solved by taking the volume rate from the last row and substituting it into the others.
	double const* const c = y + idxr.offsetC();
	const double v = y[idxr.offsetV()];
	const int volumeIdx = idxr.offsetV() - idxr.offsetC();

	const double vDot = rhs[volumeIdx];

	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		rhs[comp] = (rhs[comp] - c[comp] * vDot) / v;

	// Bound states in rapid equilibrium with the bulk liquid add V^p d_j \dot{q} to every bulk row.
	// Their rates are known here, as the relaxed initialization leaves them where it found them.
	if (_disc.strideParLiquid != 0)
		return;

	const double vSolid = static_cast<double>(_constSolidVolume);

	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		const double factor = vSolid * static_cast<double>(_parTypeVolFrac[type]) / v;
		unsigned int const* const boundOffset = _disc.boundOffset + type * _disc.nComp;
		unsigned int const* const nBound = _disc.nBound + type * _disc.nComp;
		double const* const typeQdot = yDot + idxr.offsetCp(ParticleTypeIndex{ type });

		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			double qDotSum = 0.0;
			double const* const localQdot = typeQdot + boundOffset[comp];
			for (unsigned int bnd = 0; bnd < nBound[comp]; ++bnd)
				qDotSum += localQdot[bnd];

			rhs[comp] -= factor * qDotSum;
		}
	}
}

void StirredTankModel::leanConsistentInitialSensitivity(const SimulationTime& simTime, const ConstSimulationState& simState,
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
		_tempState[idxr.offsetV()] = -adRes[idxr.offsetV()].getADValue(param);

		// Compute right hand side by adding -dF / dy * s = -J * s to -dF / dp which is already stored in _tempState
		multiplyWithJacobian(simTime, simState, sensY, -1.0, 1.0, _tempState);

		std::copy(_tempState + idxr.offsetC(), _tempState + idxr.offsetCp(), sensYdot + idxr.offsetC());
		sensYdot[idxr.offsetV()] = _tempState[idxr.offsetV()];

		solveBulkTimeDerivativeSystem(sensYdot + idxr.offsetC(), simState.vecStateY, sensYdot);
	}
}

}  // namespace model

}  // namespace cadet
