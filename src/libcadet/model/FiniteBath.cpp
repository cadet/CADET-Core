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
#include "BindingModelFactory.hpp"
#include "ReactionModelFactory.hpp"
#include "ParticleModelFactory.hpp"
#include "ParamReaderHelper.hpp"
#include "ParamReaderScopes.hpp"
#include "cadet/Exceptions.hpp"
#include "cadet/ExternalFunction.hpp"
#include "cadet/SolutionRecorder.hpp"
#include "ConfigurationHelper.hpp"
#include "model/BindingModel.hpp"
#include "model/ReactionModel.hpp"
#include "model/ParameterDependence.hpp"
#include "model/parts/BindingCellKernel.hpp"
#include "SimulationTypes.hpp"
#include "linalg/Norms.hpp"
#include "linalg/Subset.hpp"

#include "AdUtils.hpp"
#include "SensParamUtil.hpp"
#include "ParallelSupport.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <algorithm>
#include <functional>
#include <numeric>
#include <iterator>

using namespace Eigen;

namespace cadet
{

namespace model
{

namespace
{
	/**
	 * @brief Returns the name of the parameter group that holds the configuration of a particle type
	 */
	inline std::string particleGroupName(unsigned int parType)
	{
		std::ostringstream oss;
		oss << "particle_type_" << std::setfill('0') << std::setw(3) << parType;
		return oss.str();
	}
}

FiniteBath::FiniteBath(UnitOpIdx unitOpIdx) : UnitOperationBase(unitOpIdx),
	_linearSolver(nullptr), _liquidVolume(0.0), _bulkPorosity(1.0), _flowRateIn(0.0), _flowRateOut(0.0),
	_analyticJac(true), _jacobianAdDirs(0), _factorizeJacobian(false), _tempState(nullptr)
{
}

FiniteBath::~FiniteBath() CADET_NOEXCEPT
{
	delete[] _tempState;

	for (IParticleModel* pm : _particles)
		delete pm;

	_particles.clear();

	_binding.clear(); // binding models are deleted in the respective particle model

	_reaction.clearDynamicReactionModels();
	delete _linearSolver;
}

unsigned int FiniteBath::numDofs() const CADET_NOEXCEPT
{
	// Inlet DOFs, bulk liquid phase and particle phases
	return 2 * _disc.nComp + _disc.parTypeOffset[_disc.nParType];
}

unsigned int FiniteBath::numPureDofs() const CADET_NOEXCEPT
{
	return numDofs() - _disc.nComp;
}

bool FiniteBath::usesAD() const CADET_NOEXCEPT
{
#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	// We always need AD if we want to check the analytical Jacobian
	return true;
#else
	// We only need AD if we are not computing the Jacobian analytically
	return !_analyticJac;
#endif
}

bool FiniteBath::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper)
{
	_disc.nComp = paramProvider.getInt("NCOMP");

	if (paramProvider.exists("NPARTYPE"))
	{
		if (paramProvider.getInt("NPARTYPE") < 0)
			throw InvalidParameterException("Number of particle types must be >= 0!");

		_disc.nParType = paramProvider.getInt("NPARTYPE");
	}
	else
		_disc.nParType = 0;

	if ((_disc.nParType == 0) && paramProvider.exists("particle_type_000"))
		throw InvalidParameterException("NPARTYPE is set to 0, but group particle_type_000 exists.");

	const bool firstConfigCall = _tempState == nullptr; // used to not multiply allocate memory

	// ==== Read discretization
	bool analyticJac = true;
	{
		if (!paramProvider.exists("discretization"))
			throw InvalidParameterException("Group discretization is missing");

		paramProvider.pushScope("discretization");

		if (firstConfigCall)
			_linearSolver = cadet::linalg::setLinearSolver(paramProvider.exists("LINEAR_SOLVER") ? paramProvider.getString("LINEAR_SOLVER") : "SparseLU");

#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
		analyticJac = paramProvider.getBool("USE_ANALYTIC_JACOBIAN");
#else
		// Default to AD Jacobian when analytic Jacobian is to be checked
		analyticJac = false;
#endif

		// Create nonlinear solver for consistent initialization
		configureNonlinearSolver(paramProvider);

		paramProvider.popScope(); // discretization
	}

	// ==== Create and configure the particle models
	Indexer idxr(_disc);
	_particles = std::vector<IParticleModel*>(_disc.nParType, nullptr);

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		const std::string parGroup = particleGroupName(parType);
		if (!paramProvider.exists(parGroup))
			throw InvalidParameterException("Group " + parGroup + " is missing");

		paramProvider.pushScope(parGroup);

		const bool filmDiffusion = paramProvider.exists("HAS_FILM_DIFFUSION") ? paramProvider.getBool("HAS_FILM_DIFFUSION") : true;
		const bool poreDiffusion = paramProvider.exists("HAS_PORE_DIFFUSION") ? paramProvider.getBool("HAS_PORE_DIFFUSION") : false;
		const bool surfaceDiffusion = paramProvider.exists("HAS_SURFACE_DIFFUSION") ? paramProvider.getBool("HAS_SURFACE_DIFFUSION") : false;

		const std::string particleType = ParticleModel(filmDiffusion, poreDiffusion, surfaceDiffusion).getParticleTransportType();

		if (particleType == "EQUILIBRIUM_PARTICLE")
			throw InvalidParameterException("Unit operation FINITE_BATH does not support particles in rapid equilibrium with the bulk liquid (see field HAS_FILM_DIFFUSION) in particle type " + std::to_string(parType) + ", use unit operation CSTR instead");

		// The liquid in a finite bath is not flowing, so a film diffusion that depends on the interstitial velocity is meaningless
		if (paramProvider.exists("FILM_DIFFUSION_DEP"))
		{
			const std::string filmDiffDep = paramProvider.getString("FILM_DIFFUSION_DEP");
			if ((filmDiffDep != "NONE") && (filmDiffDep != "IDENTITY"))
				throw InvalidParameterException("Unit operation FINITE_BATH does not support a velocity dependent film diffusion, but FILM_DIFFUSION_DEP is set to " + filmDiffDep + " in particle type " + std::to_string(parType));
		}

		paramProvider.popScope(); // particle_type_xxx

		_particles[parType] = helper.createParticleModel(particleType);

		if (!_particles[parType])
			throw InvalidParameterException("Unknown particle model " + particleType);
	}

	bool particleConfSuccess = true;
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		particleConfSuccess = _particles[parType]->configureModelDiscretization(paramProvider, helper, _disc.nComp, parType, _disc.nParType, idxr.strideColComp()) && particleConfSuccess;

	// ==== Precompute the particle discretization dependent sizes and offsets
	if (firstConfigCall)
	{
		_disc.nParPoints = new unsigned int[std::max(_disc.nParType, 1u)];
		_disc.boundOffset = new unsigned int[std::max(_disc.nComp * _disc.nParType, 1u)];
		_disc.nBound = new unsigned int[std::max(_disc.nComp * _disc.nParType, 1u)];
		_disc.strideBound = new unsigned int[_disc.nParType + 1];
		_disc.nBoundBeforeType = new unsigned int[std::max(_disc.nParType, 1u)];
		_disc.parTypeOffset = new unsigned int[_disc.nParType + 1];
	}

	for (unsigned int type = 0; type < _disc.nParType; ++type)
		_disc.nParPoints[type] = _particles[type]->nDiscPoints();

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			_disc.nBound[parType * _disc.nComp + comp] = _particles[parType]->nBound()[comp];
	}

	const unsigned int nTotalBound = std::accumulate(_disc.nBound, _disc.nBound + _disc.nComp * _disc.nParType, 0u);

	_disc.strideBound[_disc.nParType] = nTotalBound;
	if (_disc.nParType > 0)
		_disc.nBoundBeforeType[0] = 0;

	for (unsigned int j = 0; j < _disc.nParType; ++j)
	{
		unsigned int* const ptrOffset = _disc.boundOffset + j * _disc.nComp;
		unsigned int* const ptrBound = _disc.nBound + j * _disc.nComp;

		ptrOffset[0] = 0;
		for (unsigned int i = 1; i < _disc.nComp; ++i)
			ptrOffset[i] = ptrOffset[i - 1] + ptrBound[i - 1];

		_disc.strideBound[j] = ptrOffset[_disc.nComp - 1] + ptrBound[_disc.nComp - 1];

		if (j != _disc.nParType - 1)
			_disc.nBoundBeforeType[j + 1] = _disc.nBoundBeforeType[j] + _disc.strideBound[j];
	}

	_disc.parTypeOffset[0] = 0;
	for (unsigned int j = 1; j < _disc.nParType + 1; ++j)
		_disc.parTypeOffset[j] = _disc.parTypeOffset[j - 1] + (_disc.nComp + _disc.strideBound[j - 1]) * _disc.nParPoints[j - 1];

	// Allocate space for initial conditions
	_initC.resize(_disc.nComp);
	_initCp.resize(_disc.nComp * _disc.nParType);
	_initCs.resize(nTotalBound);

	// ==== Construct and configure the bulk liquid reaction model
	bool reactionConfSuccess = true;
	_reaction.clearDynamicReactionModels();

	if (paramProvider.exists("NREAC_LIQUID"))
	{
		const int nReactions = paramProvider.getInt("NREAC_LIQUID");
		reactionConfSuccess = _reaction.configureDiscretization("liquid",
			nReactions,
			_disc.nComp,
			_disc.nParType > 0 ? _disc.nBound : nullptr,
			_disc.nParType > 0 ? _disc.boundOffset : nullptr,
			paramProvider,
			helper) && reactionConfSuccess;
	}
	else
		_reaction.empty();

	// ==== Binding models are owned by the particle models, only pointers are copied here
	_binding = std::vector<IBindingModel*>(_disc.nParType, nullptr);
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		_binding[parType] = _particles[parType]->getBinding();

		if (parType > 0)
		{
			if (_binding[parType] && (_singleBinding != !_particles[parType]->bindingParDep()))
				throw InvalidParameterException("Binding particle type dependence must be the same for all particle types, check field BINDING_PARTYPE_DEPENDENT");
		}
		else
			_singleBinding = _binding[parType] ? !_particles[parType]->bindingParDep() : true;
	}

	// Allocate memory
	if (firstConfigCall)
		_tempState = new double[numDofs()];

	_globalJac.resize(numDofs(), numDofs());
	_globalJacDisc.resize(numDofs(), numDofs());
	// pattern is set in configure(), after the particle models are fully configured

	useAnalyticJacobian(analyticJac);

	return particleConfSuccess && reactionConfSuccess;
}

bool FiniteBath::configure(IParameterProvider& paramProvider)
{
	_parameters.clear();

	// Read the geometry of the vessel
	_liquidVolume = paramProvider.getDouble("LIQUID_VOLUME");
	if (static_cast<double>(_liquidVolume) <= 0.0)
		throw InvalidParameterException("Field LIQUID_VOLUME has to be positive");

	if (_disc.nParType > 0)
	{
		if (!paramProvider.exists("BULK_POROSITY"))
			throw InvalidParameterException("The required parameter \"BULK_POROSITY\" was not found");

		_bulkPorosity = paramProvider.getDouble("BULK_POROSITY");

		if ((static_cast<double>(_bulkPorosity) <= 0.0) || (static_cast<double>(_bulkPorosity) > 1.0))
			throw InvalidParameterException("Field BULK_POROSITY has to be in (0, 1]");
	}
	else
		_bulkPorosity = paramProvider.exists("BULK_POROSITY") ? paramProvider.getDouble("BULK_POROSITY") : 1.0;

	_parameters[makeParamId(hashString("LIQUID_VOLUME"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_liquidVolume;
	_parameters[makeParamId(hashString("BULK_POROSITY"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_bulkPorosity;

	// Check whether PAR_TYPE_VOLFRAC is required or not
	if ((_disc.nParType > 1) && !paramProvider.exists("PAR_TYPE_VOLFRAC"))
		throw InvalidParameterException("The required parameter \"PAR_TYPE_VOLFRAC\" was not found");

	if (_disc.nParType > 1)
	{
		readScalarParameterOrArray(_parTypeVolFrac, paramProvider, "PAR_TYPE_VOLFRAC", 1);

		if (_parTypeVolFrac.size() != _disc.nParType)
			throw InvalidParameterException("Number of elements in field PAR_TYPE_VOLFRAC does not match number of particle types");

		const double volFracSum = std::accumulate(_parTypeVolFrac.begin(), _parTypeVolFrac.end(), 0.0,
			[](double a, const active& b) -> double { return a + static_cast<double>(b); });
		if (std::abs(1.0 - volFracSum) > 1e-10)
			throw InvalidParameterException("Sum of field PAR_TYPE_VOLFRAC differs from 1.0 (is " + std::to_string(volFracSum) + ")");

		registerParam1DArray(_parameters, _parTypeVolFrac, [=, this](bool multi, unsigned int type) { return makeParamId(hashString("PAR_TYPE_VOLFRAC"), _unitOpIdx, CompIndep, type, BoundStateIndep, ReactionIndep, SectionIndep); });
	}
	else if (_disc.nParType == 1)
		_parTypeVolFrac = std::vector<active>(1, 1.0);

	// Register initial conditions parameters
	registerParam1DArray(_parameters, _initC, [=, this](bool multi, unsigned int comp) { return makeParamId(hashString("INIT_C"), _unitOpIdx, comp, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep); });

	if (_singleBinding)
	{
		for (unsigned int c = 0; c < _disc.nComp; ++c)
			_parameters[makeParamId(hashString("INIT_CP"), _unitOpIdx, c, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_initCp[c];
	}
	else
		registerParam2DArray(_parameters, _initCp, [=, this](bool multi, unsigned int type, unsigned int comp) { return makeParamId(hashString("INIT_CP"), _unitOpIdx, comp, type, BoundStateIndep, ReactionIndep, SectionIndep); }, _disc.nComp);

	if (!_binding.empty())
	{
		const unsigned int maxBoundStates = *std::max_element(_disc.strideBound, _disc.strideBound + _disc.nParType);
		std::vector<ParameterId> initParams(maxBoundStates);

		if (_singleBinding)
		{
			_binding[0]->fillBoundPhaseInitialParameters(initParams.data(), _unitOpIdx, ParTypeIndep);

			active* const iq = _initCs.data() + _disc.nBoundBeforeType[0];
			for (unsigned int i = 0; i < _disc.strideBound[0]; ++i)
				_parameters[initParams[i]] = iq + i;
		}
		else
		{
			for (unsigned int type = 0; type < _disc.nParType; ++type)
			{
				_binding[type]->fillBoundPhaseInitialParameters(initParams.data(), _unitOpIdx, type);

				active* const iq = _initCs.data() + _disc.nBoundBeforeType[type];
				for (unsigned int i = 0; i < _disc.strideBound[type]; ++i)
					_parameters[initParams[i]] = iq + i;
			}
		}
	}

	// Reconfigure particle models
	bool particleConfSuccess = true;
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		particleConfSuccess = _particles[parType]->configure(_unitOpIdx, paramProvider, _parameters, _disc.nParType, _disc.nBoundBeforeType, _disc.strideBound[_disc.nParType]) && particleConfSuccess;

	// Reconfigure bulk liquid reaction model
	bool dynReactionConfSuccess = true;
	if (paramProvider.exists("NREAC_LIQUID"))
	{
		dynReactionConfSuccess = _reaction.configure("liquid", 0, _unitOpIdx, paramProvider) && dynReactionConfSuccess;
		dynReactionConfSuccess = _reaction.configureConservedMoities("liquid", _disc.nComp, 1e-14) && dynReactionConfSuccess;
	}

	const bool hasBulkReaction = _reaction.getDynReactionVector("liquid")[0] != nullptr;

	setJacobianPattern(_globalJac, 0, hasBulkReaction);
	resizeConservedMoietyJacobianBuffer();
	_globalJacDisc = _globalJac;
	_linearSolver->analyzePattern(_globalJacDisc.block(_disc.nComp, _disc.nComp, numPureDofs(), numPureDofs()));

	return particleConfSuccess && dynReactionConfSuccess;
}

unsigned int FiniteBath::threadLocalMemorySize() const CADET_NOEXCEPT
{
	LinearMemorySizer lms;

	// Memory for residualImpl()
	for (unsigned int i = 0; i < _disc.nParType; ++i)
	{
		if (_binding[i] && _binding[i]->requiresWorkspace())
			lms.fitBlock(_binding[i]->workspaceSize(_disc.nComp, _disc.strideBound[i], _disc.nBound + i * _disc.nComp));
	}
	_reaction.setWorkspaceRequirements("liquid", _disc.nComp, 0, lms);

	const unsigned int maxStrideBound = (_disc.nParType > 0)
		? *std::max_element(_disc.strideBound, _disc.strideBound + _disc.nParType)
		: 0u;
	lms.add<active>(_disc.nComp + maxStrideBound);
	lms.add<double>((maxStrideBound + _disc.nComp) * (maxStrideBound + _disc.nComp));

	lms.commit();
	const std::size_t resImplSize = lms.bufferSize();

	// Memory for consistentInitialState()
	lms.add<double>(_nonlinearSolver->workspaceSize(_disc.nComp + maxStrideBound) * sizeof(double));
	lms.add<double>(_disc.nComp + maxStrideBound);
	lms.add<double>(_disc.nComp + maxStrideBound);
	lms.add<double>(_disc.nComp + maxStrideBound);
	lms.add<double>((_disc.nComp + maxStrideBound) * (_disc.nComp + maxStrideBound));
	lms.add<double>(_disc.nComp);

	lms.addBlock(resImplSize);
	lms.commit();

	// Memory for consistentInitialSensitivity()
	lms.add<double>(_disc.nComp + maxStrideBound);
	lms.add<double>(maxStrideBound);
	lms.commit();

	return lms.bufferSize();
}

unsigned int FiniteBath::numAdDirsForJacobian() const CADET_NOEXCEPT
{
	// The bulk block is seeded densely (it is only nComp wide), every particle block gets its own
	// dedicated directions, just as in ColumnModel1D.
	Indexer idxr(_disc);

	int sumParBandwidth = 0;
	for (unsigned int type = 0; type < _disc.nParType; ++type)
		sumParBandwidth += idxr.strideParBlock(type);

	return _disc.nComp + sumParBandwidth;
}

void FiniteBath::useAnalyticJacobian(const bool analyticJac)
{
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	_analyticJac = analyticJac;
	_jacobianAdDirs = _analyticJac ? 0 : numAdDirsForJacobian();
#else
	// If CADET_CHECK_ANALYTIC_JACOBIAN is active, we always enable AD for comparison and use it in simulation
	_analyticJac = false;
	_jacobianAdDirs = numAdDirsForJacobian();
#endif
}

unsigned int FiniteBath::requiredADdirs() const CADET_NOEXCEPT
{
	const unsigned int numDirsBinding = maxBindingAdDirs();
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	return numDirsBinding + _jacobianAdDirs;
#else
	return numDirsBinding + numAdDirsForJacobian();
#endif
}

void FiniteBath::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, const ConstSimulationState& simState, const AdJacobianParams& adJac)
{
	// The liquid volume is constant, which requires the flow rates to cancel
	if (std::abs(static_cast<double>(_flowRateIn) - static_cast<double>(_flowRateOut)) > 1e-12 * std::max(1.0, std::abs(static_cast<double>(_flowRateIn))))
		throw InvalidParameterException("Inlet and outlet flow rate of unit " + std::to_string(_unitOpIdx) + " differ in section " + std::to_string(secIdx) + ", which is incompatible with the constant liquid volume of a finite bath");

	const bool hasReaction = _reaction.getDynReactionVector("liquid")[0];
	setJacobianPattern(_globalJac, secIdx, hasReaction);
	resizeConservedMoietyJacobianBuffer();
	_globalJacDisc = _globalJac;

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		_particles[parType]->notifyDiscontinuousSectionTransition(t, secIdx);

	_disc.curSection = secIdx;
}

void FiniteBath::setFlowRates(active const* in, active const* out) CADET_NOEXCEPT
{
	_flowRateIn = in[0];
	_flowRateOut = out[0];
}

void FiniteBath::reportSolution(ISolutionRecorder& recorder, double const* const solution) const
{
	Exporter expr(_disc, *this, solution);
	recorder.beginUnitOperation(_unitOpIdx, *this, expr);
	recorder.endUnitOperation();
}

void FiniteBath::reportSolutionStructure(ISolutionRecorder& recorder) const
{
	Exporter expr(_disc, *this, nullptr);
	recorder.unitOperationStructure(_unitOpIdx, *this, expr);
}

void FiniteBath::prepareADvectors(const AdJacobianParams& adJac) const
{
	// Early out if AD is disabled
	if (!adJac.adY)
		return;

	Indexer idxr(_disc);

	// The bulk block only has nComp entries and may be dense because of liquid phase reactions,
	// so we simply give every bulk DOF its own direction.
	active* const adVecBulk = adJac.adY + idxr.offsetC();
	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
	{
		adVecBulk[comp].fillADValue(adJac.adDirOffset, 0.0);
		adVecBulk[comp].setADValue(adJac.adDirOffset + comp, 1.0);
	}

	// Every particle type gets its own dedicated set of directions, the block is treated as dense
	unsigned int adDirOffset = adJac.adDirOffset + _disc.nComp;

	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		active* const adVec = adJac.adY + idxr.offsetCp(ParticleTypeIndex{ type });

		for (int eq = 0; eq < idxr.strideParBlock(type); ++eq)
		{
			adVec[eq].fillADValue(adJac.adDirOffset, 0.0);
			adVec[eq].setADValue(adDirOffset + eq, 1.0);
		}

		adDirOffset += idxr.strideParBlock(type);
	}
}

void FiniteBath::extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset)
{
	Indexer idxr(_disc);

	// The entries are written through row iterators, which only touch the declared sparsity pattern.
	// Inserting entries would make the pattern of the Jacobian and of its time discretized counterpart
	// diverge, and an entry outside the pattern is zero by construction anyway.

	// Extract the bulk block, which holds the flow, the liquid reaction and the diagonal of the film diffusion
	{
		linalg::BandedEigenSparseRowIterator jacBulk(_globalJac, idxr.offsetC());

		for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++jacBulk)
		{
			for (unsigned int toComp = 0; toComp < _disc.nComp; ++toComp)
				jacBulk[static_cast<int>(toComp) - static_cast<int>(comp)] = adRes[idxr.offsetC() + comp].getADValue(adDirOffset + toComp);
		}
	}

	// Read the particle Jacobian entries from the dedicated AD directions
	unsigned int offsetParticleTypeDirs = adDirOffset + _disc.nComp;

	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		const int offsetPar = idxr.offsetCp(ParticleTypeIndex{ type });
		linalg::BandedEigenSparseRowIterator jacPar(_globalJac, offsetPar);

		for (int phase = 0; phase < idxr.strideParBlock(type); ++phase, ++jacPar)
		{
			for (int phaseTo = 0; phaseTo < idxr.strideParBlock(type); ++phaseTo)
				jacPar[phaseTo - phase] = adRes[offsetPar + phase].getADValue(offsetParticleTypeDirs + phaseTo);
		}

		offsetParticleTypeDirs += idxr.strideParBlock(type);
	}

	// The film diffusion couples the bulk and the particle blocks, which are seeded with different
	// directions, so these entries are added analytically
	const active velocity = 0.0;

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		_particles[parType]->calcFilmDiffJacobian(_disc.curSection, idxr.offsetCp(ParticleTypeIndex{ parType }), idxr.offsetC(), 1, _disc.nParType, static_cast<double>(_bulkPorosity), &_parTypeVolFrac[0], &velocity, _globalJac, true);

	const auto& cm = _reaction.conservedMoieties("liquid");
	if (cm.isEnabled() && (cm.numEquilibriumReactions() > 0) && (_disc.nParType > 0))
	{
		cm.applyToMatrix(_globalJac, _disc.nComp, idxr.offsetC(), idxr.offsetCp(), _globalJac.cols(),
			_cMJacobianEntries.data(), _cMJacobianEntries.size());
	}

	if (!_globalJac.isCompressed())
		_globalJac.makeCompressed();
}

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN

void FiniteBath::checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const
{
	Indexer idxr(_disc);

	double maxDiff = 0.0;

	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
	{
		for (unsigned int toComp = 0; toComp < _disc.nComp; ++toComp)
		{
			const double adVal = adRes[idxr.offsetC() + comp].getADValue(adDirOffset + toComp);
			const double anaVal = _globalJac.coeff(idxr.offsetC() + comp, idxr.offsetC() + toComp);
			maxDiff = std::max(maxDiff, std::abs(adVal - anaVal));
		}
	}

	unsigned int offsetParticleTypeDirs = adDirOffset + _disc.nComp;
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		const int offsetPar = idxr.offsetCp(ParticleTypeIndex{ type });

		for (int phase = 0; phase < idxr.strideParBlock(type); ++phase)
		{
			for (int phaseTo = 0; phaseTo < idxr.strideParBlock(type); ++phaseTo)
			{
				const double adVal = adRes[offsetPar + phase].getADValue(offsetParticleTypeDirs + phaseTo);
				const double anaVal = _globalJac.coeff(offsetPar + phase, offsetPar + phaseTo);
				maxDiff = std::max(maxDiff, std::abs(adVal - anaVal));
			}
		}

		offsetParticleTypeDirs += idxr.strideParBlock(type);
	}

	LOG(Debug) << "-> Max diff: " << maxDiff;
}

#endif

int FiniteBath::jacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	_factorizeJacobian = true;

	if (_analyticJac)
		return residualImpl<double, double, double, true, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
	else
		return residualWithJacobian(simTime, ConstSimulationState{ simState.vecStateY, nullptr }, nullptr, adJac, threadLocalMem);
}

int FiniteBath::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	return residualImpl<double, double, double, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
}

int FiniteBath::residualWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	return residual(simTime, simState, res, adJac, threadLocalMem, true, false);
}

int FiniteBath::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res,
	const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem, bool updateJacobian, bool paramSensitivity)
{
	if (updateJacobian)
	{
		_factorizeJacobian = true;

#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
		if (_analyticJac)
		{
			if (paramSensitivity)
			{
				const int retCode = residualImpl<double, active, active, true>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, adJac.adRes, threadLocalMem);

				if (res)
					ad::copyFromAd(adJac.adRes, res, numDofs());

				return retCode;
			}
			else
				return residualImpl<double, double, double, true>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
		}
		else
		{
			ad::copyToAd(simState.vecStateY, adJac.adY, numDofs());
			ad::resetAd(adJac.adRes, numDofs());

			int retCode = 0;
			if (paramSensitivity)
				retCode = residualImpl<active, active, active, false>(simTime.t, simTime.secIdx, adJac.adY, simState.vecStateYdot, adJac.adRes, threadLocalMem);
			else
				retCode = residualImpl<active, active, double, false>(simTime.t, simTime.secIdx, adJac.adY, simState.vecStateYdot, adJac.adRes, threadLocalMem);

			if (res)
				ad::copyFromAd(adJac.adRes, res, numDofs());

			extractJacobianFromAD(adJac.adRes, adJac.adDirOffset);

			return retCode;
		}
#else
		ad::copyToAd(simState.vecStateY, adJac.adY, numDofs());
		ad::resetAd(adJac.adRes, numDofs());

		int retCode = 0;
		if (paramSensitivity)
			retCode = residualImpl<active, active, active, false>(simTime.t, simTime.secIdx, adJac.adY, simState.vecStateYdot, adJac.adRes, threadLocalMem);
		else
			retCode = residualImpl<active, active, double, false>(simTime.t, simTime.secIdx, adJac.adY, simState.vecStateYdot, adJac.adRes, threadLocalMem);

		if (res)
		{
			retCode = residualImpl<double, double, double, true>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);

			checkAnalyticJacobianAgainstAd(adJac.adRes, adJac.adDirOffset);
		}

		extractJacobianFromAD(adJac.adRes, adJac.adDirOffset);

		return retCode;
#endif
	}
	else
	{
		if (paramSensitivity)
		{
			ad::resetAd(adJac.adRes, numDofs());

			const int retCode = residualImpl<double, active, active, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, adJac.adRes, threadLocalMem);

			if (res)
				ad::copyFromAd(adJac.adRes, res, numDofs());

			return retCode;
		}
		else
			return residualImpl<double, double, double, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
	}
}

template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
int FiniteBath::residualImpl(double t, unsigned int secIdx, StateType const* const y, double const* const yDot, ResidualType* const res, util::ThreadLocalStorage& threadLocalMem)
{
	if (wantRes)
	{
		Eigen::Map<Eigen::Vector<ResidualType, Dynamic>> resi(res, numDofs());
		resi.setZero();
	}

	LinearBufferAllocator tlmAlloc = threadLocalMem.get();
	Indexer idxr(_disc);

	if (wantJac)
	{
		_globalJac.coeffs().setZero();

		// Static (per section) part of the Jacobian: outflow, particle diffusion and film diffusion
		linalg::BandedEigenSparseRowIterator jacBulk(_globalJac, idxr.offsetC());
		const double flowOut = static_cast<double>(_flowRateOut) / static_cast<double>(_liquidVolume);
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++jacBulk)
			jacBulk[0] = flowOut;

		for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
			_particles[parType]->calcParticleDiffJacobian(secIdx, 0, idxr.offsetCp(ParticleTypeIndex{ parType }), _globalJac);

		const active velocity = 0.0;
		for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
			_particles[parType]->calcFilmDiffJacobian(secIdx, idxr.offsetCp(ParticleTypeIndex{ parType }), idxr.offsetC(), 1, _disc.nParType, static_cast<double>(_bulkPorosity), &_parTypeVolFrac[0], &velocity, _globalJac);
	}

	residualBulk<StateType, ResidualType, ParamType, wantJac, wantRes>(t, secIdx, y, yDot, res, threadLocalMem);

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		linalg::BandedEigenSparseRowIterator jacIt;

		if (wantJac)
			jacIt = linalg::BandedEigenSparseRowIterator(_globalJac, idxr.offsetCp(ParticleTypeIndex{ parType }));

		const active velocity = 0.0;
		model::columnPackingParameters packing
		{
			_parTypeVolFrac[parType],
			_bulkPorosity,
			ColumnPosition{ 0.0, 0.0, 0.0 },
			velocity
		};

		_particles[parType]->residual(t, secIdx,
			y + idxr.offsetCp(ParticleTypeIndex{ parType }),
			y + idxr.offsetC(),
			yDot ? yDot + idxr.offsetCp(ParticleTypeIndex{ parType }) : nullptr,
			res ? res + idxr.offsetCp(ParticleTypeIndex{ parType }) : nullptr,
			res ? res + idxr.offsetC() : nullptr,
			packing, jacIt, tlmAlloc,
			typename cadet::ParamSens<ParamType>::enabled()
		);
	}

	const auto& cm = _reaction.conservedMoieties("liquid");
	if (cm.isEnabled() && cm.numEquilibriumReactions() > 0)
	{
		const unsigned int nMoieties = cm.numMoieties();
		const unsigned int nEq = cm.numEquilibriumReactions();

		if (wantRes)
		{
			BufferedArray<ResidualType> oldRes = tlmAlloc.array<ResidualType>(_disc.nComp);
			ResidualType* const resC = res + idxr.offsetC();

			std::copy_n(resC, _disc.nComp, static_cast<ResidualType*>(oldRes));
			cm.applyToVector(resC, static_cast<ResidualType*>(oldRes), _disc.nComp);

			ResidualType* const eqRes = resC + nMoieties;
			std::fill_n(eqRes, nEq, 0.0);

			unsigned int eqIdx = 0;
			for (auto const* reaction : _reaction.getDynReactionVector("liquid"))
			{
				if (!reaction)
					continue;

				LinearBufferAllocator subAlloc = tlmAlloc.manageRemainingMemory();
				reaction->residualEquilibriumFlux(
					t, secIdx, ColumnPosition{ 0.0, 0.0, 0.0 },
					_disc.nComp, y + idxr.offsetC(), eqRes, eqIdx, subAlloc);
			}

			cadet_assert(eqIdx == nEq);
		}

		if (wantJac)
		{
			LinearBufferAllocator subAlloc = tlmAlloc.manageRemainingMemory();

			cm.applyToMatrix(_globalJac, _disc.nComp, idxr.offsetC(), 0, _globalJac.cols(),
				_cMJacobianEntries.data(), _cMJacobianEntries.size());

			// Replace the remaining component rows by equilibrium equations
			for (unsigned int r = 0; r < nEq; ++r)
				_globalJac.innerVector(idxr.offsetC() + nMoieties + r) *= 0.0;

			linalg::BandedEigenSparseRowIterator eqJac(_globalJac, idxr.offsetC() + nMoieties);
			unsigned int eqIdxJac = 0;
			for (auto const* reaction : _reaction.getDynReactionVector("liquid"))
			{
				if (!reaction)
					continue;

				reaction->analyticEquilibriumJacobian(
					t, secIdx, ColumnPosition{ 0.0, 0.0, 0.0 },
					_disc.nComp, reinterpret_cast<double const*>(y + idxr.offsetC()),
					eqIdxJac, nMoieties, eqJac, subAlloc);
			}

			cadet_assert(eqIdxJac == nEq);
		}
	}

	if (!wantRes)
		return 0;

	// Handle inlet DOFs, which are simply copied to the residual
	for (unsigned int i = 0; i < _disc.nComp; ++i)
		res[i] = y[i];

	return 0;
}

template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
int FiniteBath::residualBulk(double t, unsigned int secIdx, StateType const* yBase, double const* yDotBase, ResidualType* resBase, util::ThreadLocalStorage& threadLocalMem)
{
	Indexer idxr(_disc);
	LinearBufferAllocator tlmAlloc = threadLocalMem.get();

	StateType const* const y = yBase + idxr.offsetC();

	if (wantRes)
	{
		const ParamType invVolume = 1.0 / static_cast<ParamType>(_liquidVolume);
		const ParamType flowIn = static_cast<ParamType>(_flowRateIn) * invVolume;
		const ParamType flowOut = static_cast<ParamType>(_flowRateOut) * invVolume;

		ResidualType* const res = resBase + idxr.offsetC();
		double const* const yDot = yDotBase ? yDotBase + idxr.offsetC() : nullptr;

		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			res[comp] += -flowIn * yBase[comp] + flowOut * y[comp];

			if (cadet_likely(yDot))
				res[comp] += yDot[comp];
		}
	}

	if (_reaction.getDynReactionVector("liquid").size() == 0)
		return 0;

	const ColumnPosition colPos{ 0.0, 0.0, 0.0 };
	linalg::BandedEigenSparseRowIterator jac(_globalJac, idxr.offsetC());

	for (auto i = 0; i < _reaction.getDynReactionVector("liquid").size(); ++i)
	{
		if (!_reaction.getDynReactionVector("liquid")[i])
			continue;

		if (wantRes)
			_reaction.getDynReactionVector("liquid")[i]->residualFluxAdd(t, secIdx, colPos, _disc.nComp, y, resBase + idxr.offsetC(), -1.0, tlmAlloc);

		if (wantJac)
			_reaction.getDynReactionVector("liquid")[i]->analyticJacobianAdd(t, secIdx, colPos, _disc.nComp, reinterpret_cast<double const*>(y), -1.0, jac, tlmAlloc);
	}

	return 0;
}

parts::cell::CellParameters FiniteBath::makeCellResidualParams(unsigned int parType, int const* qsReaction) const
{
	return parts::cell::CellParameters
		{
			_disc.nComp,
			_disc.nBound + _disc.nComp * parType,
			_disc.boundOffset + _disc.nComp * parType,
			_disc.strideBound[parType],
			qsReaction,
			_particles[parType]->getPorosity(),
			_particles[parType]->getPoreAccessFactor(),
			_binding[parType],
			nullptr
		};
}

int FiniteBath::residualSensFwdWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidualSens);

	return residual(simTime, simState, nullptr, adJac, threadLocalMem, true, true);
}

int FiniteBath::residualSensFwdAdOnly(const SimulationTime& simTime, const ConstSimulationState& simState, active* const adRes, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidualSens);

	return residualImpl<double, active, active, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, adRes, threadLocalMem);
}

int FiniteBath::residualSensFwdCombine(const SimulationTime& simTime, const ConstSimulationState& simState,
	const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS, active const* adRes,
	double* const tmp1, double* const tmp2, double* const tmp3)
{
	BENCH_SCOPE(_timerResidualSens);

	for (std::size_t param = 0; param < yS.size(); ++param)
	{
		// Directional derivative (dF / dy) * s
		multiplyWithJacobian(SimulationTime{ 0.0, 0u }, ConstSimulationState{ nullptr, nullptr }, yS[param], 1.0, 0.0, tmp1);

		// Directional derivative (dF / dyDot) * sDot
		multiplyWithDerivativeJacobian(SimulationTime{ 0.0, 0u }, ConstSimulationState{ nullptr, nullptr }, ySdot[param], tmp2);

		double* const ptrResS = resS[param];

		for (unsigned int i = 0; i < numDofs(); ++i)
			ptrResS[i] = tmp1[i] + tmp2[i] + adRes[i].getADValue(param);
	}

	return 0;
}

void FiniteBath::multiplyWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* yS, double alpha, double beta, double* ret)
{
	Indexer idxr(_disc);

	// Handle identity matrix of inlet DOFs
	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		ret[comp] = alpha * yS[comp] + beta * ret[comp];

	// Main Jacobian
	Eigen::Map<Eigen::VectorXd> ret_vec(ret + idxr.offsetC(), numPureDofs());
	Eigen::Map<const Eigen::VectorXd> yS_vec(yS + idxr.offsetC(), numPureDofs());
	ret_vec = alpha * _globalJac.block(idxr.offsetC(), idxr.offsetC(), numPureDofs(), numPureDofs()) * yS_vec + beta * ret_vec;

	// Inlet coupling
	const double jacInlet = inletJacobian();

	const auto& cm = _reaction.conservedMoieties("liquid");
	if (cm.isEnabled() && cm.numEquilibriumReactions() > 0)
	{
		const auto& L = cm.conservedMoietyMatrix();
		for (unsigned int moiety = 0; moiety < cm.numMoieties(); ++moiety)
		{
			double inletDirection = 0.0;
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				inletDirection += L(moiety, comp) * yS[comp];

			ret[idxr.offsetC() + moiety] += alpha * jacInlet * inletDirection;
		}
	}
	else
	{
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			ret[idxr.offsetC() + comp] += alpha * jacInlet * yS[comp];
	}
}

void FiniteBath::multiplyWithDerivativeJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* sDot, double* ret)
{
	Indexer idxr(_disc);

	// Bulk liquid phase
	std::copy_n(sDot + idxr.offsetC(), _disc.nComp, ret + idxr.offsetC());

	// Particles
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		const double invBetaP = (1.0 / static_cast<double>(_particles[type]->getPorosity()) - 1.0);
		unsigned int const* const nBound = _disc.nBound + type * _disc.nComp;
		unsigned int const* const boundOffset = _disc.boundOffset + type * _disc.nComp;
		int const* const qsReaction = _binding[type]->reactionQuasiStationarity();

		const int offsetCpType = idxr.offsetCp(ParticleTypeIndex{ type });
		for (unsigned int shell = 0; shell < _disc.nParPoints[type]; ++shell)
		{
			const int offsetCpShell = offsetCpType + shell * idxr.strideParNode(type);
			double const* const mobileSdot = sDot + offsetCpShell;
			double* const mobileRet = ret + offsetCpShell;

			parts::cell::multiplyWithDerivativeJacobianKernel<true>(mobileSdot, mobileRet, _disc.nComp, nBound, boundOffset, _disc.strideBound[type], qsReaction, 1.0, invBetaP);
		}
	}

	const auto& cm = _reaction.conservedMoieties("liquid");
	if (cm.isEnabled() && cm.numEquilibriumReactions() > 0)
		cm.applyToDerivativeVector(ret + idxr.offsetC(), sDot + idxr.offsetC(), _disc.nComp);

	// Handle inlet DOFs (all algebraic)
	std::fill_n(ret, _disc.nComp, 0.0);
}

void FiniteBath::setExternalFunctions(IExternalFunction** extFuns, unsigned int size)
{
	for (IBindingModel* bm : _binding)
	{
		if (bm)
			bm->setExternalFunctions(extFuns, size);
	}
}

bool FiniteBath::setParameter(const ParameterId& pId, double value)
{
	if (pId.unitOperation == _unitOpIdx)
	{
		const int mpIc = multiplexInitialConditions(pId, value, false);
		if (mpIc > 0)
			return true;
		else if (mpIc < 0)
			return false;

		bool paramExists = false;
		for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		{
			const bool paramExistsNow = _particles[parType]->setParameter(pId, value);
			paramExists = paramExists || paramExistsNow;

			if (paramExists)
			{
				if ((pId.name == hashString("PAR_RADIUS")) || (pId.name == hashString("PAR_CORERADIUS")))
					_particles[parType]->updateRadialDisc();

				if ((pId.particleType != ParTypeIndep && parType == pId.particleType) || (pId.particleType == ParTypeIndep && parType == _disc.nParType - 1))
					return true;
			}
		}

		if (model::setParameter(pId, value, _reaction.getDynReactionVector("liquid"), false))
			return true;
	}

	return UnitOperationBase::setParameter(pId, value);
}

bool FiniteBath::setParameter(const ParameterId& pId, int value)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return false;

	bool paramExists = false;
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		const bool paramExistsNow = _particles[parType]->setParameter(pId, value);
		paramExists = paramExists || paramExistsNow;

		if (paramExists)
		{
			if ((pId.particleType != ParTypeIndep && parType == pId.particleType) || (pId.particleType == ParTypeIndep && parType == _disc.nParType - 1))
				return true;
		}
	}

	if (pId.unitOperation == _unitOpIdx)
	{
		if (model::setParameter(pId, value, _reaction.getDynReactionVector("liquid"), false))
			return true;
	}

	return UnitOperationBase::setParameter(pId, value);
}

bool FiniteBath::setParameter(const ParameterId& pId, bool value)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return false;

	bool paramExists = false;
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		const bool paramExistsNow = _particles[parType]->setParameter(pId, value);
		paramExists = paramExists || paramExistsNow;

		if (paramExists)
		{
			if ((pId.particleType != ParTypeIndep && parType == pId.particleType) || (pId.particleType == ParTypeIndep && parType == _disc.nParType - 1))
				return true;
		}
	}

	if (pId.unitOperation == _unitOpIdx)
	{
		if (model::setParameter(pId, value, _reaction.getDynReactionVector("liquid"), false))
			return true;
	}

	return UnitOperationBase::setParameter(pId, value);
}

void FiniteBath::setSensitiveParameterValue(const ParameterId& pId, double value)
{
	if (pId.unitOperation == _unitOpIdx)
	{
		if (multiplexInitialConditions(pId, value, true) != 0)
			return;

		bool paramExists = false;
		for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		{
			const bool paramExistsNow = _particles[parType]->setSensitiveParameterValue(_sensParams, pId, value);
			paramExists = paramExists || paramExistsNow;

			if (paramExists)
			{
				if ((pId.name == hashString("PAR_RADIUS")) || (pId.name == hashString("PAR_CORERADIUS")))
					_particles[parType]->updateRadialDisc();

				if ((pId.particleType != ParTypeIndep && parType == pId.particleType) || (pId.particleType == ParTypeIndep && parType == _disc.nParType - 1))
					return;
			}
		}

		if (model::setSensitiveParameterValue(pId, value, _sensParams, _reaction.getDynReactionVector("liquid"), false))
			return;
	}

	UnitOperationBase::setSensitiveParameterValue(pId, value);
}

bool FiniteBath::setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue)
{
	if (pId.unitOperation == _unitOpIdx)
	{
		const int mpIc = multiplexInitialConditions(pId, adDirection, adValue);
		if (mpIc > 0)
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}
		else if (mpIc < 0)
			return false;

		bool paramExists = false;
		for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		{
			const bool paramExistsNow = _particles[parType]->setSensitiveParameter(_sensParams, pId, adDirection, adValue);
			paramExists = paramExists || paramExistsNow;

			if (paramExists)
			{
				if ((pId.name == hashString("PAR_RADIUS")) || (pId.name == hashString("PAR_CORERADIUS")))
					_particles[parType]->updateRadialDisc();

				if ((pId.particleType != ParTypeIndep && parType == pId.particleType) || (pId.particleType == ParTypeIndep && parType == _disc.nParType - 1))
				{
					LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
					return true;
				}
			}
		}

		if (model::setSensitiveParameter(pId, adDirection, adValue, _sensParams, _reaction.getDynReactionVector("liquid"), false))
		{
			LOG(Debug) << "Found parameter " << pId << " in DynamicBulkReactionModel: Dir " << adDirection << " is set to " << adValue;
			return true;
		}
	}

	return UnitOperationBase::setSensitiveParameter(pId, adDirection, adValue);
}

std::unordered_map<ParameterId, double> FiniteBath::getAllParameterValues() const
{
	std::unordered_map<ParameterId, double> data = UnitOperationBase::getAllParameterValues();

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		_particles[parType]->getAllParameterValues(data);

	return data;
}

double FiniteBath::getParameterDouble(const ParameterId& pId) const
{
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		const double val = _particles[parType]->getParameterDouble(pId);
		if (val)
			return val;
	}

	return UnitOperationBase::getParameterDouble(pId);
}

bool FiniteBath::hasParameter(const ParameterId& pId) const
{
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		if (_particles[parType]->hasParameter(pId))
			return true;
	}

	return UnitOperationBase::hasParameter(pId);
}

void FiniteBath::resizeConservedMoietyJacobianBuffer()
{
	const auto& cm = _reaction.conservedMoieties("liquid");
	if (!cm.isEnabled() || (cm.numEquilibriumReactions() == 0))
		return;

	Indexer idxr(_disc);
	std::size_t numEntries = 0;

	for (unsigned int state = 0; state < _disc.nComp; ++state)
	{
		for (SparseMatrix<double, RowMajor>::InnerIterator it(_globalJac, idxr.offsetC() + state); it; ++it)
			++numEntries;
	}

	_cMJacobianEntries.resize(numEntries);
}

int FiniteBath::Exporter::writeMobilePhase(double* buffer) const
{
	std::copy_n(_idx.c(_data), _disc.nComp, buffer);
	return _disc.nComp;
}

int FiniteBath::Exporter::writeSolidPhase(double* buffer) const
{
	int numWritten = 0;
	for (unsigned int i = 0; i < _disc.nParType; ++i)
	{
		const int n = writeSolidPhase(i, buffer);
		buffer += n;
		numWritten += n;
	}
	return numWritten;
}

int FiniteBath::Exporter::writeParticleMobilePhase(double* buffer) const
{
	int numWritten = 0;
	for (unsigned int i = 0; i < _disc.nParType; ++i)
	{
		const int n = writeParticleMobilePhase(i, buffer);
		buffer += n;
		numWritten += n;
	}
	return numWritten;
}

int FiniteBath::Exporter::writeSolidPhase(unsigned int parType, double* buffer) const
{
	cadet_assert(parType < _disc.nParType);

	const unsigned int stride = _disc.nComp + _disc.strideBound[parType];
	double const* ptr = _data + _idx.offsetCp(ParticleTypeIndex{ parType }) + _disc.nComp;

	for (unsigned int j = 0; j < _disc.nParPoints[parType]; ++j)
	{
		std::copy_n(ptr, _disc.strideBound[parType], buffer);
		buffer += _disc.strideBound[parType];
		ptr += stride;
	}

	return _disc.nParPoints[parType] * _disc.strideBound[parType];
}

int FiniteBath::Exporter::writeParticleMobilePhase(unsigned int parType, double* buffer) const
{
	cadet_assert(parType < _disc.nParType);

	const unsigned int stride = _disc.nComp + _disc.strideBound[parType];
	double const* ptr = _data + _idx.offsetCp(ParticleTypeIndex{ parType });

	for (unsigned int j = 0; j < _disc.nParPoints[parType]; ++j)
	{
		std::copy_n(ptr, _disc.nComp, buffer);
		buffer += _disc.nComp;
		ptr += stride;
	}

	return _disc.nParPoints[parType] * _disc.nComp;
}

int FiniteBath::Exporter::writeInlet(unsigned int port, double* buffer) const
{
	cadet_assert(port == 0);
	std::copy_n(_data, _disc.nComp, buffer);
	return _disc.nComp;
}

int FiniteBath::Exporter::writeInlet(double* buffer) const
{
	std::copy_n(_data, _disc.nComp, buffer);
	return _disc.nComp;
}

int FiniteBath::Exporter::writeOutlet(unsigned int port, double* buffer) const
{
	cadet_assert(port == 0);
	std::copy_n(_idx.c(_data), _disc.nComp, buffer);
	return _disc.nComp;
}

int FiniteBath::Exporter::writeOutlet(double* buffer) const
{
	std::copy_n(_idx.c(_data), _disc.nComp, buffer);
	return _disc.nComp;
}

int FiniteBath::linearSolve(double t, double alpha, double outerTol, double* const rhs, double const* const weight,
	const ConstSimulationState& simState)
{
	BENCH_SCOPE(_timerLinearSolve);

	Indexer idxr(_disc);

	Eigen::Map<VectorXd> r(rhs, numDofs());

	if (_factorizeJacobian)
	{
		assembleDiscretizedGlobalJacobian(alpha, idxr);

		BENCH_START(_timerFactorize);
		_linearSolver->factorize(_globalJacDisc.block(idxr.offsetC(), idxr.offsetC(), numPureDofs(), numPureDofs()));
		BENCH_STOP(_timerFactorize);

		if (cadet_unlikely(_linearSolver->info() != Eigen::Success))
		{
			LOG(Error) << "Factorize() failed";
		}

		_factorizeJacobian = false;
	}

	// Handle inlet DOFs: solve J c_uo = b_uo - A * b_in
	const double jacInlet = inletJacobian();

	const auto& cm = _reaction.conservedMoieties("liquid");
	if (cm.isEnabled() && cm.numEquilibriumReactions() > 0)
	{
		const auto& L = cm.conservedMoietyMatrix();
		for (unsigned int moiety = 0; moiety < cm.numMoieties(); ++moiety)
		{
			double inletValue = 0.0;
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				inletValue += L(moiety, comp) * r[comp];

			r[idxr.offsetC() + moiety] -= jacInlet * inletValue;
		}
	}
	else
	{
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			r[idxr.offsetC() + comp] -= jacInlet * r[comp];
	}

	r.segment(idxr.offsetC(), numPureDofs()) = _linearSolver->solve(r.segment(idxr.offsetC(), numPureDofs()));

	if (cadet_unlikely(_linearSolver->info() != Eigen::Success))
	{
		LOG(Error) << "Solve() failed";
	}

	return 0;
}

void FiniteBath::assembleDiscretizedGlobalJacobian(double alpha, Indexer idxr)
{
	// Add the static (per section) Jacobian without inlet
	_globalJacDisc = _globalJac;

	// Add time derivatives to the particle shells
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		linalg::BandedEigenSparseRowIterator jac(_globalJacDisc, idxr.offsetCp(ParticleTypeIndex{ parType }));

		for (unsigned int j = 0; j < _disc.nParPoints[parType]; ++j)
		{
			addTimeDerivativeToJacobianParticleShell(jac, idxr, alpha, parType);
			// Iterator jac has already been advanced to next shell
		}
	}

	const auto& cm = _reaction.conservedMoieties("liquid");
	if (cm.isEnabled() && cm.numEquilibriumReactions() > 0)
	{
		const auto& L = cm.conservedMoietyMatrix();
		const unsigned int nMoieties = cm.numMoieties();

		for (unsigned int m = 0; m < nMoieties; ++m)
		{
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				_globalJacDisc.coeffRef(idxr.offsetC() + m, idxr.offsetC() + comp) += alpha * L(m, comp);
		}
		// No mass-matrix entries for equilibrium rows
	}
	else
	{
		linalg::BandedEigenSparseRowIterator jac(_globalJacDisc, idxr.offsetC());
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++jac)
			jac[0] += alpha;
	}
}

void FiniteBath::addTimeDerivativeToJacobianParticleShell(linalg::BandedEigenSparseRowIterator& jac, const Indexer& idxr, double alpha, unsigned int parType)
{
	parts::cell::addTimeDerivativeToJacobianParticleShell<linalg::BandedEigenSparseRowIterator, true>(jac, alpha, static_cast<double>(_particles[parType]->getPorosity()), _disc.nComp, _disc.nBound + _disc.nComp * parType,
		_particles[parType]->getPoreAccessFactor(), _disc.strideBound[parType], _disc.boundOffset + _disc.nComp * parType, _binding[parType]->reactionQuasiStationarity());
}

}  // namespace model

}  // namespace cadet

#include "model/FiniteBath-InitialConditions.cpp"
