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

StirredTankModel::StirredTankModel(UnitOpIdx unitOpIdx) : UnitOperationBase(unitOpIdx),
	_linearSolver(nullptr), _initLiquidVolume(0.0), _constSolidVolume(0.0), _constantVolume(false),
	_flowRateIn(0.0), _flowRateOut(0.0), _curFlowRateFilter(0.0),
	_analyticJac(true), _jacobianAdDirs(0), _factorizeJacobian(false), _tempState(nullptr)
{
}

StirredTankModel::~StirredTankModel() CADET_NOEXCEPT
{
	delete[] _tempState;

	for (IParticleModel* pm : _particles)
		delete pm;

	_particles.clear();

	_binding.clear(); // binding models are deleted in the respective particle model

	_reaction.clearDynamicReactionModels();
	delete _linearSolver;
}

unsigned int StirredTankModel::numDofs() const CADET_NOEXCEPT
{
	// Inlet DOFs, bulk liquid phase, particle phases and the bulk liquid volume
	return 2 * _disc.nComp + _disc.parTypeOffset[_disc.nParType] + 1;
}

unsigned int StirredTankModel::numPureDofs() const CADET_NOEXCEPT
{
	return numDofs() - _disc.nComp;
}

bool StirredTankModel::usesAD() const CADET_NOEXCEPT
{
#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	// We always need AD if we want to check the analytical Jacobian
	return true;
#else
	// We only need AD if we are not computing the Jacobian analytically
	return !_analyticJac;
#endif
}

bool StirredTankModel::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper)
{
	const int nComp = paramProvider.getInt("NCOMP");
	if (nComp < 1)
		throw InvalidParameterException("Field NCOMP has to be at least 1, but is " + std::to_string(nComp));

	_disc.nComp = nComp;

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
	// A well mixed vessel whose particles are not spatially resolved has nothing to discretize, so
	// the group is optional and everything in it falls back to its default
	bool analyticJac = true;
	{
		const bool hasDiscretization = paramProvider.exists("discretization");

		if (hasDiscretization)
			paramProvider.pushScope("discretization");

		if (firstConfigCall)
			_linearSolver = cadet::linalg::setLinearSolver(paramProvider.exists("LINEAR_SOLVER") ? paramProvider.getString("LINEAR_SOLVER") : "SparseLU");

#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
		analyticJac = paramProvider.exists("USE_ANALYTIC_JACOBIAN") ? paramProvider.getBool("USE_ANALYTIC_JACOBIAN") : true;
#else
		// Default to AD Jacobian when analytic Jacobian is to be checked
		analyticJac = false;
#endif

		// Create nonlinear solver for consistent initialization
		configureNonlinearSolver(paramProvider);

		if (hasDiscretization)
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

		// The particle model of a particle type follows from its transport fields alone. The field that
		// decides between rapid equilibrium and a particle liquid phase changes the state layout of the
		// whole unit operation, so it is stated rather than defaulted.
		if (!paramProvider.exists("HAS_FILM_DIFFUSION"))
			throw InvalidParameterException("The required parameter \"HAS_FILM_DIFFUSION\" was not found in particle type " + std::to_string(parType));

		const bool filmDiffusion = paramProvider.getBool("HAS_FILM_DIFFUSION");
		const bool poreDiffusion = paramProvider.exists("HAS_PORE_DIFFUSION") ? paramProvider.getBool("HAS_PORE_DIFFUSION") : false;
		const bool surfaceDiffusion = paramProvider.exists("HAS_SURFACE_DIFFUSION") ? paramProvider.getBool("HAS_SURFACE_DIFFUSION") : false;

		const std::string particleType = ParticleModel(filmDiffusion, poreDiffusion, surfaceDiffusion).getParticleTransportType();

		// A particle in rapid equilibrium with the bulk liquid owns no liquid phase, so its state block
		// holds bound states only. Mixing the two kinds would make the block layout depend on the
		// particle type, which the indexer of this unit operation does not express.
		const bool isEquilibrium = particleType == "LIQUID_EQUILIBRIUM_PARTICLE";

		if (parType == 0)
			_disc.strideParLiquid = isEquilibrium ? 0u : _disc.nComp;
		else if (isEquilibrium != (_disc.strideParLiquid == 0u))
			throw InvalidParameterException("Particle types in rapid equilibrium with the bulk liquid (see field HAS_FILM_DIFFUSION) must not be mixed with other particle types, check particle type " + std::to_string(parType));

		// The liquid in a stirred tank is not flowing, so a film diffusion that depends on the interstitial velocity is meaningless
		if (paramProvider.exists("FILM_DIFFUSION_DEP"))
		{
			const std::string filmDiffDep = paramProvider.getString("FILM_DIFFUSION_DEP");
			if ((filmDiffDep != "NONE") && (filmDiffDep != "IDENTITY"))
				throw InvalidParameterException("Unit operation CSTR does not support a velocity dependent film diffusion, but FILM_DIFFUSION_DEP is set to " + filmDiffDep + " in particle type " + std::to_string(parType));
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

	// The particle models state their own block size, which is what makes room for a particle that
	// owns no liquid phase and therefore contributes bound states only
	_disc.parTypeOffset[0] = 0;
	for (unsigned int j = 1; j < _disc.nParType + 1; ++j)
		_disc.parTypeOffset[j] = _disc.parTypeOffset[j - 1] + _particles[j - 1]->strideParBlock();

	// A particle in rapid equilibrium binds to the bulk liquid, which is not adjacent to its bound
	// states, so it has to be told how far back the liquid phase lies
	for (unsigned int j = 0; j < _disc.nParType; ++j)
		_particles[j]->setBulkRowOffset(idxr.offsetCp(ParticleTypeIndex{ j }) - idxr.offsetC());

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

bool StirredTankModel::configure(IParameterProvider& paramProvider)
{
	_parameters.clear();

	// Read the flow rate of liquid drawn without components
	_flowRateFilter.clear();
	if (paramProvider.exists("FLOWRATE_FILTER"))
	{
		readScalarParameterOrArray(_flowRateFilter, paramProvider, "FLOWRATE_FILTER", 1);

		for (std::size_t i = 0; i < _flowRateFilter.size(); ++i)
		{
			if (static_cast<double>(_flowRateFilter[i]) < 0.0)
				throw InvalidParameterException("Field FLOWRATE_FILTER has to be non-negative, but entry " + std::to_string(i) + " is " + std::to_string(static_cast<double>(_flowRateFilter[i])));
		}

		registerParam1DArray(_parameters, _flowRateFilter, [=, this](bool multi, unsigned int sec) { return makeParamId(hashString("FLOWRATE_FILTER"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, multi ? sec : SectionIndep); });
	}

	// Read the geometry of the vessel. The volumes are given by one of two mutually exclusive
	// parameterizations:
	//  - variable volume: INIT_LIQUID_VOLUME is the initial bulk liquid volume, which evolves with the
	//    net flow rate, and CONST_SOLID_VOLUME is the particle volume including the pore liquid
	//  - constant volume: LIQUID_VOLUME is the prescribed bulk liquid volume and the porosity the
	//    ratio of bulk liquid volume to total volume, which together fix the particle volume
	//
	// The porosity is named after what it describes, consistently with the 1D column units: particles
	// that have pores hold liquid of their own, so the ratio is the one of the bulk liquid alone and
	// the field is COL_POROSITY. Particles in rapid equilibrium hold no liquid, so all of it is bulk
	// liquid and the field is TOTAL_POROSITY. Without particles there is nothing to lump, so
	// COL_POROSITY applies as well.
	const bool usesColPorosity = (_disc.nParType == 0) || (_disc.strideParLiquid != 0);
	const char* const porosityField = usesColPorosity ? "COL_POROSITY" : "TOTAL_POROSITY";

	const bool hasVariableVolume = paramProvider.exists("INIT_LIQUID_VOLUME") || paramProvider.exists("CONST_SOLID_VOLUME");
	const bool hasConstantVolume = paramProvider.exists("LIQUID_VOLUME") || paramProvider.exists(porosityField);

	if (hasVariableVolume && hasConstantVolume)
		throw InvalidParameterException("Fields LIQUID_VOLUME and " + std::string(porosityField) + " must not be combined with INIT_LIQUID_VOLUME or CONST_SOLID_VOLUME");

	if (!hasVariableVolume && !hasConstantVolume)
		throw InvalidParameterException("Field INIT_LIQUID_VOLUME or LIQUID_VOLUME required");

	_constantVolume = hasConstantVolume;

	if (_constantVolume)
	{
		if (!paramProvider.exists("LIQUID_VOLUME"))
			throw InvalidParameterException("Field LIQUID_VOLUME required when " + std::string(porosityField) + " is given");

		_initLiquidVolume = paramProvider.getDouble("LIQUID_VOLUME");

		if (static_cast<double>(_initLiquidVolume) <= 0.0)
			throw InvalidParameterException("Field LIQUID_VOLUME has to be positive");

		if ((_disc.nParType > 0) && !paramProvider.exists(porosityField))
			throw InvalidParameterException("The required parameter \"" + std::string(porosityField) + "\" was not found");

		const double porosity = paramProvider.exists(porosityField) ? paramProvider.getDouble(porosityField) : 1.0;

		if ((porosity <= 0.0) || (porosity > 1.0))
			throw InvalidParameterException("Field " + std::string(porosityField) + " has to be in (0, 1]");

		// V^p = (V^l + V^p) * (1 - eps) -> V^p = V^l * (1 - eps) / eps
		_constSolidVolume = static_cast<double>(_initLiquidVolume) * (1.0 - porosity) / porosity;
	}
	else
	{
		if (!paramProvider.exists("INIT_LIQUID_VOLUME"))
			throw InvalidParameterException("Field INIT_LIQUID_VOLUME required when CONST_SOLID_VOLUME is given");

		_initLiquidVolume = paramProvider.getDouble("INIT_LIQUID_VOLUME");

		if (static_cast<double>(_initLiquidVolume) < 0.0)
			throw InvalidParameterException("Field INIT_LIQUID_VOLUME has to be non-negative");

		if (_disc.nParType > 0)
		{
			if (!paramProvider.exists("CONST_SOLID_VOLUME"))
				throw InvalidParameterException("The required parameter \"CONST_SOLID_VOLUME\" was not found");
		}

		_constSolidVolume = paramProvider.exists("CONST_SOLID_VOLUME") ? paramProvider.getDouble("CONST_SOLID_VOLUME") : 0.0;

		if (static_cast<double>(_constSolidVolume) < 0.0)
			throw InvalidParameterException("Field CONST_SOLID_VOLUME has to be non-negative");
	}

	// The liquid volume is registered under the name that was used to specify it, since changing it
	// means something different in either parameterization
	if (_constantVolume)
		_parameters[makeParamId(hashString("LIQUID_VOLUME"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_initLiquidVolume;
	else
		_parameters[makeParamId(hashString("INIT_LIQUID_VOLUME"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_initLiquidVolume;

	_parameters[makeParamId(hashString("CONST_SOLID_VOLUME"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_constSolidVolume;

	// Check whether PAR_TYPE_VOLFRAC is required or not
	if ((_disc.nParType > 1) && !paramProvider.exists("PAR_TYPE_VOLFRAC"))
		throw InvalidParameterException("The required parameter \"PAR_TYPE_VOLFRAC\" was not found");

	if (_disc.nParType > 1)
	{
		readScalarParameterOrArray(_parTypeVolFrac, paramProvider, "PAR_TYPE_VOLFRAC", 1);

		if (_parTypeVolFrac.size() != _disc.nParType)
			throw InvalidParameterException("Number of elements in field PAR_TYPE_VOLFRAC does not match number of particle types");

		// A single volume fraction is a fraction of the solid phase, so it is in [0, 1] on its own.
		// The sum below does not imply this, since it is also met by entries that cancel each other.
		for (std::size_t i = 0; i < _parTypeVolFrac.size(); ++i)
		{
			const double volFrac = static_cast<double>(_parTypeVolFrac[i]);
			if ((volFrac < 0.0) || (volFrac > 1.0))
				throw InvalidParameterException("Field PAR_TYPE_VOLFRAC has to be in [0.0, 1.0], but entry " + std::to_string(i) + " is " + std::to_string(volFrac));
		}

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

	// A particle in rapid equilibrium with the bulk liquid owns no liquid phase, so it has no initial
	// particle liquid concentration either
	if (_disc.strideParLiquid > 0)
	{
		if (_singleBinding)
		{
			for (unsigned int c = 0; c < _disc.nComp; ++c)
				_parameters[makeParamId(hashString("INIT_CP"), _unitOpIdx, c, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_initCp[c];
		}
		else
			registerParam2DArray(_parameters, _initCp, [=, this](bool multi, unsigned int type, unsigned int comp) { return makeParamId(hashString("INIT_CP"), _unitOpIdx, comp, type, BoundStateIndep, ReactionIndep, SectionIndep); }, _disc.nComp);
	}

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
		particleConfSuccess = _particles[parType]->configure(_unitOpIdx, paramProvider, _parameters, _disc.nParType, _disc.nBoundBeforeType, _disc.strideBound[_disc.nParType], false) && particleConfSuccess;

	// Reconfigure bulk liquid reaction model
	bool dynReactionConfSuccess = true;
	if (paramProvider.exists("NREAC_LIQUID"))
	{
		dynReactionConfSuccess = _reaction.configure("liquid", 0, _unitOpIdx, paramProvider, true) && dynReactionConfSuccess;
	}

	// A cross phase reaction of a particle type in rapid equilibrium acts on the bulk liquid, so it
	// couples the bulk components just like a reaction of the bulk itself
	const bool hasBulkReaction = (_reaction.getDynReactionVector("liquid")[0] != nullptr) || hasEquilibriumParticleCrossPhaseReaction();

	setJacobianPattern(_globalJac, 0, hasBulkReaction);
	resizeConservedMoietyJacobianBuffer();
	_globalJacDisc = _globalJac;
	_linearSolver->analyzePattern(_globalJacDisc, _disc.nComp, numPureDofs());

	return particleConfSuccess && dynReactionConfSuccess;
}

unsigned int StirredTankModel::threadLocalMemorySize() const CADET_NOEXCEPT
{
	LinearMemorySizer lms;

	// Memory for residualImpl()
	for (unsigned int i = 0; i < _disc.nParType; ++i)
	{
		if (_binding[i] && _binding[i]->requiresWorkspace())
			lms.fitBlock(_binding[i]->workspaceSize(_disc.nComp, _disc.strideBound[i], _disc.nBound + i * _disc.nComp));

		// The solid phase and cross phase reactions of a particle type in rapid equilibrium are
		// evaluated here, see residualEquilibriumParticleReactions()
		const ReactionSystem& parReaction = *_particles[i]->getReaction();
		parReaction.setWorkspaceRequirements("cross_phase", _disc.nComp, _disc.strideBound[i], lms);
		parReaction.setWorkspaceRequirements("solid", _disc.nComp, _disc.strideBound[i], lms);
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

unsigned int StirredTankModel::numAdDirsForJacobian() const CADET_NOEXCEPT
{
	// The bulk block is seeded densely (it is only nComp wide), every particle block gets its own
	// dedicated directions, just as in ColumnModel1D. The liquid volume gets one of its own, since
	// the extensive bulk balance depends on it and nothing else supplies that column.
	Indexer idxr(_disc);

	int sumParBandwidth = 0;
	for (unsigned int type = 0; type < _disc.nParType; ++type)
		sumParBandwidth += idxr.strideParBlock(type);

	return _disc.nComp + sumParBandwidth + 1;
}

void StirredTankModel::useAnalyticJacobian(const bool analyticJac)
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

unsigned int StirredTankModel::requiredADdirs() const CADET_NOEXCEPT
{
	const unsigned int numDirsBinding = maxBindingAdDirs();
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	return numDirsBinding + _jacobianAdDirs;
#else
	return numDirsBinding + numAdDirsForJacobian();
#endif
}

void StirredTankModel::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, const ConstSimulationState& simState, const AdJacobianParams& adJac)
{
	if (_flowRateFilter.size() > 1)
		_curFlowRateFilter = _flowRateFilter[secIdx];
	else if (_flowRateFilter.size() == 1)
		_curFlowRateFilter = _flowRateFilter[0];

	// A prescribed constant liquid volume requires a vanishing net flow rate
	if (_constantVolume)
	{
		const double flowIn = static_cast<double>(_flowRateIn);
		const double flowOut = static_cast<double>(_flowRateOut);
		const double flowFilter = static_cast<double>(_curFlowRateFilter);

		if (std::abs(flowIn - flowOut - flowFilter) > 1e-12 * std::max(1.0, std::abs(flowIn)))
			throw InvalidParameterException("Inlet, outlet and filter flow rate of unit " + std::to_string(_unitOpIdx) + " do not cancel in section " + std::to_string(secIdx) + ", which is required for the constant liquid volume given by LIQUID_VOLUME");
	}

	const bool hasReaction = (_reaction.getDynReactionVector("liquid")[0] != nullptr) || hasEquilibriumParticleCrossPhaseReaction();
	setJacobianPattern(_globalJac, secIdx, hasReaction);
	resizeConservedMoietyJacobianBuffer();
	_globalJacDisc = _globalJac;

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		_particles[parType]->notifyDiscontinuousSectionTransition(t, secIdx);

	_disc.curSection = secIdx;
}

void StirredTankModel::setFlowRates(active const* in, active const* out) CADET_NOEXCEPT
{
	_flowRateIn = in[0];
	_flowRateOut = out[0];
}

void StirredTankModel::reportSolution(ISolutionRecorder& recorder, double const* const solution) const
{
	Exporter expr(_disc, *this, solution);
	recorder.beginUnitOperation(_unitOpIdx, *this, expr);
	recorder.endUnitOperation();
}

void StirredTankModel::reportSolutionStructure(ISolutionRecorder& recorder) const
{
	Exporter expr(_disc, *this, nullptr);
	recorder.unitOperationStructure(_unitOpIdx, *this, expr);
}

void StirredTankModel::prepareADvectors(const AdJacobianParams& adJac) const
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

	// The liquid volume gets the direction behind the particle blocks
	active& adVolume = adJac.adY[idxr.offsetV()];
	adVolume.fillADValue(adJac.adDirOffset, 0.0);
	adVolume.setADValue(adDirOffset, 1.0);
}

void StirredTankModel::extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset)
{
	Indexer idxr(_disc);

	// The entries are written through row iterators, which only touch the declared sparsity pattern.
	// Inserting entries would make the pattern of the Jacobian and of its time discretized counterpart
	// diverge, and an entry outside the pattern is zero by construction anyway.

	// Extract the bulk block, which holds the flow, the liquid reaction and the diagonal of the film diffusion
	{
		linalg::BandedEigenSparseRowIterator jacBulk(_globalJac, idxr.offsetC());

		// The direction of the liquid volume sits behind the bulk and the particle blocks
		unsigned int volumeDir = adDirOffset + _disc.nComp;
		for (unsigned int type = 0; type < _disc.nParType; ++type)
			volumeDir += idxr.strideParBlock(type);

		const int volumeOffset = idxr.offsetV() - idxr.offsetC();

		for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++jacBulk)
		{
			for (unsigned int toComp = 0; toComp < _disc.nComp; ++toComp)
				jacBulk[static_cast<int>(toComp) - static_cast<int>(comp)] = adRes[idxr.offsetC() + comp].getADValue(adDirOffset + toComp);

			jacBulk[volumeOffset - static_cast<int>(comp)] = adRes[idxr.offsetC() + comp].getADValue(volumeDir);
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

		// A particle in rapid equilibrium binds to the bulk liquid, so its bound states depend on the
		// bulk concentrations. That coupling sits outside the particle block and is not covered by the
		// film diffusion below, which this particle does not have, so it is read from the bulk
		// directions here.
		if (_disc.strideParLiquid == 0)
		{
			const int bulkOffset = idxr.offsetC() - offsetPar;
			linalg::BandedEigenSparseRowIterator jacQ(_globalJac, offsetPar);

			for (int phase = 0; phase < idxr.strideParBlock(type); ++phase, ++jacQ)
			{
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
					jacQ[bulkOffset + static_cast<int>(comp) - phase] = adRes[offsetPar + phase].getADValue(adDirOffset + comp);
			}

			// The reverse coupling comes from the solid phase and cross phase reactions of such a
			// particle type, which the hosting unit operation evaluates against the bulk liquid, see
			// residualEquilibriumParticleReactions()
			linalg::BandedEigenSparseRowIterator jacC(_globalJac, idxr.offsetC());

			for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++jacC)
			{
				const int parOffset = offsetPar - static_cast<int>(idxr.offsetC() + comp);

				for (int phaseTo = 0; phaseTo < idxr.strideParBlock(type); ++phaseTo)
					jacC[parOffset + phaseTo] = adRes[idxr.offsetC() + comp].getADValue(offsetParticleTypeDirs + phaseTo);
			}
		}

		offsetParticleTypeDirs += idxr.strideParBlock(type);
	}

	// The film diffusion couples the bulk and the particle blocks, which are seeded with different
	// directions, so these entries are added analytically
	const active velocity = 0.0;
	const bulkFilmTransfer& bulkFilm = filmTransferOperator();

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		_particles[parType]->calcFilmDiffJacobian(_disc.curSection, idxr.offsetCp(ParticleTypeIndex{ parType }), idxr.offsetC(), 1, _disc.nParType, bulkFilm, &_parTypeVolFrac[0], &velocity, _globalJac, true);

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

void StirredTankModel::checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const
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

int StirredTankModel::jacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	_factorizeJacobian = true;

	if (_analyticJac)
		return residualImpl<double, double, double, true, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
	else
		return residualWithJacobian(simTime, ConstSimulationState{ simState.vecStateY, nullptr }, nullptr, adJac, threadLocalMem);
}

int StirredTankModel::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	return residualImpl<double, double, double, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
}

int StirredTankModel::residualWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	return residual(simTime, simState, res, adJac, threadLocalMem, true, false);
}

int StirredTankModel::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res,
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
int StirredTankModel::residualImpl(double t, unsigned int secIdx, StateType const* const y, double const* const yDot, ResidualType* const res, util::ThreadLocalStorage& threadLocalMem)
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

		// Static (per section) part of the Jacobian: outflow, particle diffusion and film diffusion.
		// The bulk balance is extensive, so the outflow enters undivided and the volume rate is added
		// by residualBulk(), which knows the current state.
		linalg::BandedEigenSparseRowIterator jacBulk(_globalJac, idxr.offsetC());
		const double flowOut = static_cast<double>(_flowRateOut);
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++jacBulk)
			jacBulk[0] = flowOut;

		for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
			_particles[parType]->calcParticleDiffJacobian(secIdx, 0, idxr.offsetCp(ParticleTypeIndex{ parType }), _globalJac);

		const active velocity = 0.0;
		const bulkFilmTransfer& bulkFilm = filmTransferOperator();

		for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
			_particles[parType]->calcFilmDiffJacobian(secIdx, idxr.offsetCp(ParticleTypeIndex{ parType }), idxr.offsetC(), 1, _disc.nParType, bulkFilm, &_parTypeVolFrac[0], &velocity, _globalJac);
	}

	residualBulk<StateType, ResidualType, ParamType, wantJac, wantRes>(t, secIdx, y, yDot, res, threadLocalMem);

	const bulkFilmTransfer& bulkFilm = filmTransferOperator();

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		linalg::BandedEigenSparseRowIterator jacIt;

		if (wantJac)
			jacIt = linalg::BandedEigenSparseRowIterator(_globalJac, idxr.offsetCp(ParticleTypeIndex{ parType }));

		const active velocity = 0.0;
		model::columnPackingParameters packing
		{
			_parTypeVolFrac[parType],
			bulkFilm,
			0,
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

	// A particle type with pores evaluates its own solid phase and cross phase reactions against its
	// own liquid phase. One in rapid equilibrium has none, so its reactions reach the bulk liquid and
	// are evaluated here, where the volumes that weigh them are known.
	if (_disc.strideParLiquid == 0)
	{
		for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
			residualEquilibriumParticleReactions<StateType, ResidualType, ParamType, wantJac, wantRes>(t, secIdx, y, res, parType, threadLocalMem);
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
int StirredTankModel::residualBulk(double t, unsigned int secIdx, StateType const* yBase, double const* yDotBase, ResidualType* resBase, util::ThreadLocalStorage& threadLocalMem)
{
	Indexer idxr(_disc);
	LinearBufferAllocator tlmAlloc = threadLocalMem.get();

	StateType const* const y = yBase + idxr.offsetC();

	// The bulk liquid volume and its rate of change are states, so the rate is taken from yDot. It
	// equals F_in - F_out - F_filter only once the volume equation below is satisfied, which is not
	// the case during a Newton iteration.
	const StateType& v = yBase[idxr.offsetV()];
	const double vDot = yDotBase ? yDotBase[idxr.offsetV()] : 0.0;

	const ParamType flowIn = static_cast<ParamType>(_flowRateIn);
	const ParamType flowOut = static_cast<ParamType>(_flowRateOut);

	if (wantRes)
	{
		ResidualType* const res = resBase + idxr.offsetC();
		double const* const yDot = yDotBase ? yDotBase + idxr.offsetC() : nullptr;

		// Concentrations, in extensive form:
		//   V^l \dot{c}_i + \dot{V}^l c_i - F_in c_{in,i} + F_out c_i + (film transfer) - V^l f_react,i = 0
		// The film transfer term is added by the particle models, weighted with V^p, and the reaction
		// term is added below.
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			res[comp] += -flowIn * yBase[comp] + flowOut * y[comp];

			if (cadet_likely(yDot))
				res[comp] += v * yDot[comp] + vDot * y[comp];
		}

		// Volume: \dot{V} - F_in + F_out + F_filter = 0
		resBase[idxr.offsetV()] = vDot - flowIn + flowOut + static_cast<ParamType>(_curFlowRateFilter);
	}

	if (wantJac)
	{
		// dRes_c / dc_i picks up the volume rate, the remaining entries are set by the static part of
		// the Jacobian and by the particle models
		linalg::BandedEigenSparseRowIterator jacBulk(_globalJac, idxr.offsetC());
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++jacBulk)
			jacBulk[0] += vDot;

		// dRes_c / dV^l = \dot{c}_i, the reaction part is added below
		if (yDotBase)
		{
			double const* const yDot = yDotBase + idxr.offsetC();
			const int volumeOffset = idxr.offsetV() - idxr.offsetC();
			linalg::BandedEigenSparseRowIterator jacVol(_globalJac, idxr.offsetC());

			for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++jacVol)
				jacVol[volumeOffset - static_cast<int>(comp)] += yDot[comp];
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

		if (wantRes || wantJac)
		{
			// The reaction term enters the extensive balance scaled with the bulk liquid volume, and
			// its unscaled flux is the derivative of that term with respect to the volume
			LinearBufferAllocator subAlloc = tlmAlloc.manageRemainingMemory();
			BufferedArray<ResidualType> flux = subAlloc.array<ResidualType>(_disc.nComp);
			std::fill_n(static_cast<ResidualType*>(flux), _disc.nComp, 0.0);

			_reaction.getDynReactionVector("liquid")[i]->residualFluxAdd(t, secIdx, colPos, _disc.nComp, y, static_cast<ResidualType*>(flux), -1.0, subAlloc);

			if (wantRes)
			{
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
					resBase[idxr.offsetC() + comp] += v * flux[comp];
			}

			if (wantJac)
			{
				const int volumeOffset = idxr.offsetV() - idxr.offsetC();
				linalg::BandedEigenSparseRowIterator jacVol(_globalJac, idxr.offsetC());

				// dRes_c / dV^l += f_react,i
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++jacVol)
					jacVol[volumeOffset - static_cast<int>(comp)] += static_cast<double>(flux[comp]);

				_reaction.getDynReactionVector("liquid")[i]->analyticJacobianAdd(t, secIdx, colPos, _disc.nComp, reinterpret_cast<double const*>(y), -static_cast<double>(v), jac, subAlloc);
			}
		}
	}

	return 0;
}

bool StirredTankModel::hasEquilibriumParticleCrossPhaseReaction() const CADET_NOEXCEPT
{
	if (_disc.strideParLiquid != 0)
		return false;

	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		for (auto const* reaction : _particles[type]->getReaction()->getDynReactionVector("cross_phase"))
		{
			if (reaction)
				return true;
		}
	}

	return false;
}

/**
 * @brief Adds the solid phase and cross phase reactions of a particle type in rapid equilibrium
 * @details A particle in rapid equilibrium has no liquid phase of its own, so its cross phase
 *          reactions act on the bulk liquid. The terms are extensive, which is why they are
 *          evaluated here rather than in @c LiquidEquilibriumParticle: the particle model knows
 *          neither @f$ V^\ell @f$ nor that the balance is written in holdups.
 *
 *          The bulk balance already carries @f$ V^p \sum_m \dot{q}_{i,m} @f$ from the binding
 *          term, which assumes that everything the bound states take up comes out of the liquid.
 *          A reaction between bound states does not, so its flux is subtracted there again. That
 *          is what sends the solid flux into the bulk rows as well.
 * @param [in] t Current time point
 * @param [in] secIdx Index of the current section
 * @param [in] y Pointer to the first element of the state vector
 * @param [in,out] res Pointer to the first element of the residual
 * @param [in] type Index of the particle type
 * @param [in,out] threadLocalMem Thread local memory
 */
template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
int StirredTankModel::residualEquilibriumParticleReactions(double t, unsigned int secIdx, StateType const* const y, ResidualType* const res, unsigned int type, util::ThreadLocalStorage& threadLocalMem)
{
	Indexer idxr(_disc);

	const unsigned int strideBound = _disc.strideBound[type];
	if (strideBound == 0)
		return 0;

	const ReactionSystem& reaction = *_particles[type]->getReaction();
	const std::vector<IDynamicReactionModel*>& crossPhase = reaction.getDynReactionVector("cross_phase");
	const std::vector<IDynamicReactionModel*>& solid = reaction.getDynReactionVector("solid");

	if (crossPhase.empty() && solid.empty())
		return 0;

	LinearBufferAllocator tlmAlloc = threadLocalMem.get();

	int const* const qsReaction = _binding[type]->reactionQuasiStationarity();
	unsigned int const* const nBound = _disc.nBound + type * _disc.nComp;

	const int qRowOffset = idxr.offsetCp(ParticleTypeIndex{ type }) - idxr.offsetC();
	const int volumeOffset = idxr.offsetV() - idxr.offsetC();

	StateType const* const c = y + idxr.offsetC();
	StateType const* const q = y + idxr.offsetCp(ParticleTypeIndex{ type });
	ResidualType* const resC = res ? res + idxr.offsetC() : nullptr;
	ResidualType* const resQ = res ? res + idxr.offsetCp(ParticleTypeIndex{ type }) : nullptr;

	// The liquid phase the cross phase reaction acts on is the bulk liquid, so it carries V^l. What
	// the reaction produces in the solid is held by the particles, which carry V^p d_j.
	const StateType& v = y[idxr.offsetV()];
	const ParamType vSolidVolFrac = static_cast<ParamType>(_constSolidVolume) * static_cast<ParamType>(_parTypeVolFrac[type]);

	const ColumnPosition colPos{ 0.0, 0.0, 0.0 };
	const int jacCols = static_cast<int>(_disc.nComp + strideBound);

	for (int phase = 0; phase < 2; ++phase)
	{
		const bool isCrossPhase = (phase == 0);
		const std::vector<IDynamicReactionModel*>& reactions = isCrossPhase ? crossPhase : solid;

		for (auto const* dynReaction : reactions)
		{
			if (!dynReaction)
				continue;

			LinearBufferAllocator subAlloc = tlmAlloc.manageRemainingMemory();

			BufferedArray<ResidualType> fluxBuffer = subAlloc.array<ResidualType>(_disc.nComp + strideBound);
			ResidualType* const fluxLiquid = static_cast<ResidualType*>(fluxBuffer);
			ResidualType* const fluxSolid = fluxLiquid + _disc.nComp;

			std::fill_n(fluxLiquid, _disc.nComp + strideBound, 0.0);

			// The unscaled flux is also the derivative of the term it enters with respect to the volume
			// that scales it, so it is needed for the Jacobian as well
			if (isCrossPhase)
				dynReaction->residualCombinedAdd(t, secIdx, colPos, c, q, fluxLiquid, fluxSolid, -1.0, subAlloc);
			else
				dynReaction->residualFluxAdd(t, secIdx, colPos, strideBound, q, fluxSolid, -1.0, subAlloc);

			if (wantRes && resC)
			{
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
					resC[comp] += v * fluxLiquid[comp];

				unsigned int idx = 0;
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				{
					for (unsigned int bnd = 0; bnd < nBound[comp]; ++bnd, ++idx)
					{
						resC[comp] += vSolidVolFrac * fluxSolid[idx];

						// A quasi-stationary bound state has no rate equation for the reaction to enter
						if (!qsReaction[idx])
							resQ[idx] += fluxSolid[idx];
					}
				}
			}

			if (!wantJac)
				continue;

			// Rows 0 to NCOMP - 1 hold the liquid phase fluxes and the remaining ones the solid phase
			// fluxes, both differentiated with respect to the bulk liquid followed by the bound states
			BufferedArray<double> jacBuffer = subAlloc.array<double>(static_cast<std::size_t>(jacCols) * jacCols);
			linalg::DenseMatrixView jacFlux(static_cast<double*>(jacBuffer), nullptr, jacCols, jacCols);
			jacFlux.setAll(0.0);

			if (isCrossPhase)
			{
				dynReaction->analyticJacobianCombinedAdd(t, secIdx, colPos, reinterpret_cast<double const*>(c), reinterpret_cast<double const*>(q),
					-1.0, jacFlux.row(0), jacFlux.row(_disc.nComp), subAlloc);
			}
			else
			{
				dynReaction->analyticJacobianAdd(t, secIdx, colPos, strideBound, reinterpret_cast<double const*>(q),
					-1.0, jacFlux.row(_disc.nComp), subAlloc);
			}

			linalg::BandedEigenSparseRowIterator jacBulk(_globalJac, idxr.offsetC());
			linalg::BandedEigenSparseRowIterator jacPar(_globalJac, idxr.offsetCp(ParticleTypeIndex{ type }));

			const double vLiquid = static_cast<double>(v);
			const double solidFactor = static_cast<double>(vSolidVolFrac);

			unsigned int idx = 0;
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				const int compShift = -static_cast<int>(comp);

				// dRes_c / dc and dRes_c / dq of the liquid phase flux, which carries V^l. The unscaled
				// flux is what that term contributes to dRes_c / dV^l.
				(jacBulk + comp)[volumeOffset + compShift] += static_cast<double>(fluxLiquid[comp]);
				(jacBulk + comp).addArray(jacFlux.rowPtr(comp), compShift, _disc.nComp, vLiquid);
				(jacBulk + comp).addArray(jacFlux.rowPtr(comp) + _disc.nComp, qRowOffset + compShift, strideBound, vLiquid);

				for (unsigned int bnd = 0; bnd < nBound[comp]; ++bnd, ++idx)
				{
					double const* const solidRow = jacFlux.rowPtr(_disc.nComp + idx);

					(jacBulk + comp).addArray(solidRow, compShift, _disc.nComp, solidFactor);
					(jacBulk + comp).addArray(solidRow + _disc.nComp, qRowOffset + compShift, strideBound, solidFactor);

					if (qsReaction[idx])
						continue;

					const int boundShift = -static_cast<int>(idx);

					(jacPar + idx).addArray(solidRow, -qRowOffset + boundShift, _disc.nComp, 1.0);
					(jacPar + idx).addArray(solidRow + _disc.nComp, boundShift, strideBound, 1.0);
				}
			}
		}
	}

	return 0;
}

parts::cell::CellParameters StirredTankModel::makeCellResidualParams(unsigned int parType, int const* qsReaction) const
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

int StirredTankModel::residualSensFwdWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidualSens);

	return residual(simTime, simState, nullptr, adJac, threadLocalMem, true, true);
}

int StirredTankModel::residualSensFwdAdOnly(const SimulationTime& simTime, const ConstSimulationState& simState, active* const adRes, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidualSens);

	return residualImpl<double, active, active, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, adRes, threadLocalMem);
}

int StirredTankModel::residualSensFwdCombine(const SimulationTime& simTime, const ConstSimulationState& simState,
	const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS, active const* adRes,
	double* const tmp1, double* const tmp2, double* const tmp3)
{
	BENCH_SCOPE(_timerResidualSens);

	for (std::size_t param = 0; param < yS.size(); ++param)
	{
		// Directional derivative (dF / dy) * s
		multiplyWithJacobian(simTime, simState, yS[param], 1.0, 0.0, tmp1);

		// Directional derivative (dF / dyDot) * sDot. The mass matrix of the bulk balance holds the
		// liquid volume and the bulk concentrations, so the state is required here.
		multiplyWithDerivativeJacobian(simTime, simState, ySdot[param], tmp2);

		double* const ptrResS = resS[param];

		for (unsigned int i = 0; i < numDofs(); ++i)
			ptrResS[i] = tmp1[i] + tmp2[i] + adRes[i].getADValue(param);
	}

	return 0;
}

void StirredTankModel::multiplyWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* yS, double alpha, double beta, double* ret)
{
	Indexer idxr(_disc);

	// Handle identity matrix of inlet DOFs
	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		ret[comp] = (beta == 0.0) ? alpha * yS[comp] : alpha * yS[comp] + beta * ret[comp];

	// Main Jacobian
	Eigen::Map<Eigen::VectorXd> ret_vec(ret + idxr.offsetC(), numPureDofs());
	Eigen::Map<const Eigen::VectorXd> yS_vec(yS + idxr.offsetC(), numPureDofs());
	if (beta == 0.0)
		ret_vec = alpha * _globalJac.block(idxr.offsetC(), idxr.offsetC(), numPureDofs(), numPureDofs()) * yS_vec;
	else
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

void StirredTankModel::multiplyWithDerivativeJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* sDot, double* ret)
{
	Indexer idxr(_disc);

	// Bulk liquid phase. The balance is extensive, so the mass matrix row of component i is
	//   V^l \dot{c}_i + c_i \dot{V}^l
	double const* const c = simState.vecStateY + idxr.offsetC();
	const double v = simState.vecStateY[idxr.offsetV()];

	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		ret[idxr.offsetC() + comp] = v * sDot[idxr.offsetC() + comp] + c[comp] * sDot[idxr.offsetV()];

	// Volume: \dot{V} - F_in + F_out + F_filter = 0
	ret[idxr.offsetV()] = sDot[idxr.offsetV()];

	// Particles
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		unsigned int const* const nBound = _disc.nBound + type * _disc.nComp;
		unsigned int const* const boundOffset = _disc.boundOffset + type * _disc.nComp;
		int const* const qsReaction = _binding[type]->reactionQuasiStationarity();

		const int offsetCpType = idxr.offsetCp(ParticleTypeIndex{ type });

		if (_disc.strideParLiquid == 0)
		{
			// A particle in rapid equilibrium holds bound states only, and whatever they take up
			// leaves the bulk liquid. Quasi-stationary states have no time derivative of their own,
			// but they do change the amount held by the solid, so they contribute to the bulk row.
			const double vSolidVolFrac = static_cast<double>(_constSolidVolume) * static_cast<double>(_parTypeVolFrac[type]);
			double* const retQ = ret + offsetCpType;
			double const* const sDotQ = sDot + offsetCpType;

			for (unsigned int idx = 0; idx < _disc.strideBound[type]; ++idx)
				retQ[idx] = qsReaction[idx] ? 0.0 : sDotQ[idx];

			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				for (unsigned int bnd = 0; bnd < nBound[comp]; ++bnd)
					ret[idxr.offsetC() + comp] += vSolidVolFrac * sDotQ[boundOffset[comp] + bnd];
			}

			continue;
		}

		const double invBetaP = (1.0 / static_cast<double>(_particles[type]->getPorosity()) - 1.0);
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

void StirredTankModel::setExternalFunctions(IExternalFunction** extFuns, unsigned int size)
{
	for (IBindingModel* bm : _binding)
	{
		if (bm)
			bm->setExternalFunctions(extFuns, size);
	}
}

bool StirredTankModel::setParameter(const ParameterId& pId, double value)
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

bool StirredTankModel::setParameter(const ParameterId& pId, int value)
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

bool StirredTankModel::setParameter(const ParameterId& pId, bool value)
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

void StirredTankModel::setSensitiveParameterValue(const ParameterId& pId, double value)
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

bool StirredTankModel::setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue)
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

std::unordered_map<ParameterId, double> StirredTankModel::getAllParameterValues() const
{
	std::unordered_map<ParameterId, double> data = UnitOperationBase::getAllParameterValues();

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		_particles[parType]->getAllParameterValues(data);

	return data;
}

double StirredTankModel::getParameterDouble(const ParameterId& pId) const
{
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		const double val = _particles[parType]->getParameterDouble(pId);
		if (val)
			return val;
	}

	return UnitOperationBase::getParameterDouble(pId);
}

bool StirredTankModel::hasParameter(const ParameterId& pId) const
{
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		if (_particles[parType]->hasParameter(pId))
			return true;
	}

	return UnitOperationBase::hasParameter(pId);
}

void StirredTankModel::resizeConservedMoietyJacobianBuffer()
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

int StirredTankModel::Exporter::writeMobilePhase(double* buffer) const
{
	std::copy_n(_idx.c(_data), _disc.nComp, buffer);
	return _disc.nComp;
}

int StirredTankModel::Exporter::writeSolidPhase(double* buffer) const
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

int StirredTankModel::Exporter::writeParticleMobilePhase(double* buffer) const
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

int StirredTankModel::Exporter::writeSolidPhase(unsigned int parType, double* buffer) const
{
	cadet_assert(parType < _disc.nParType);

	const unsigned int stride = _idx.strideParNode(parType);
	double const* ptr = _data + _idx.offsetCp(ParticleTypeIndex{ parType }) + _idx.strideParLiquid();

	for (unsigned int j = 0; j < _disc.nParPoints[parType]; ++j)
	{
		std::copy_n(ptr, _disc.strideBound[parType], buffer);
		buffer += _disc.strideBound[parType];
		ptr += stride;
	}

	return _disc.nParPoints[parType] * _disc.strideBound[parType];
}

int StirredTankModel::Exporter::writeParticleMobilePhase(unsigned int parType, double* buffer) const
{
	cadet_assert(parType < _disc.nParType);

	const unsigned int stride = _idx.strideParNode(parType);
	const unsigned int strideLiquid = _idx.strideParLiquid();
	double const* ptr = _data + _idx.offsetCp(ParticleTypeIndex{ parType });

	for (unsigned int j = 0; j < _disc.nParPoints[parType]; ++j)
	{
		std::copy_n(ptr, strideLiquid, buffer);
		buffer += strideLiquid;
		ptr += stride;
	}

	return _disc.nParPoints[parType] * strideLiquid;
}

int StirredTankModel::Exporter::writeInlet(unsigned int port, double* buffer) const
{
	cadet_assert(port == 0);
	std::copy_n(_data, _disc.nComp, buffer);
	return _disc.nComp;
}

int StirredTankModel::Exporter::writeInlet(double* buffer) const
{
	std::copy_n(_data, _disc.nComp, buffer);
	return _disc.nComp;
}

int StirredTankModel::Exporter::writeOutlet(unsigned int port, double* buffer) const
{
	cadet_assert(port == 0);
	std::copy_n(_idx.c(_data), _disc.nComp, buffer);
	return _disc.nComp;
}

int StirredTankModel::Exporter::writeOutlet(double* buffer) const
{
	std::copy_n(_idx.c(_data), _disc.nComp, buffer);
	return _disc.nComp;
}

int StirredTankModel::linearSolve(double t, double alpha, double outerTol, double* const rhs, double const* const weight,
	const ConstSimulationState& simState)
{
	BENCH_SCOPE(_timerLinearSolve);

	Indexer idxr(_disc);

	Eigen::Map<VectorXd> r(rhs, numDofs());

	if (_factorizeJacobian)
	{
		assembleDiscretizedGlobalJacobian(alpha, idxr, simState.vecStateY);

		BENCH_START(_timerFactorize);
		_linearSolver->factorize(_globalJacDisc, idxr.offsetC(), numPureDofs());
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

	_linearSolver->solve(r.data() + idxr.offsetC(), numPureDofs());

	if (cadet_unlikely(_linearSolver->info() != Eigen::Success))
	{
		LOG(Error) << "Solve() failed";
	}

	return 0;
}

void StirredTankModel::assembleDiscretizedGlobalJacobian(double alpha, Indexer idxr, double const* const y)
{
	// Add the static (per section) Jacobian without inlet
	_globalJacDisc = _globalJac;

	// Add time derivatives to the particle shells
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		linalg::BandedEigenSparseRowIterator jac(_globalJacDisc, idxr.offsetCp(ParticleTypeIndex{ parType }));

		if (_disc.strideParLiquid == 0)
		{
			// A particle in rapid equilibrium holds bound states only. Its own rows carry a time
			// derivative unless the state is quasi-stationary, and the amount it takes up leaves the
			// bulk liquid, which the extensive balance receives weighted with the particle volume.
			int const* const qsReaction = _binding[parType]->reactionQuasiStationarity();
			unsigned int const* const nBound = _disc.nBound + parType * _disc.nComp;
			unsigned int const* const boundOffset = _disc.boundOffset + parType * _disc.nComp;
			const double vSolidVolFrac = static_cast<double>(_constSolidVolume) * static_cast<double>(_parTypeVolFrac[parType]);
			const int offsetCpType = idxr.offsetCp(ParticleTypeIndex{ parType });

			for (unsigned int idx = 0; idx < _disc.strideBound[parType]; ++idx, ++jac)
			{
				if (!qsReaction[idx])
					jac[0] += alpha;
			}

			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				for (unsigned int bnd = 0; bnd < nBound[comp]; ++bnd)
					_globalJacDisc.coeffRef(idxr.offsetC() + comp, offsetCpType + boundOffset[comp] + bnd) += alpha * vSolidVolFrac;
			}

			continue;
		}

		for (unsigned int j = 0; j < _disc.nParPoints[parType]; ++j)
		{
			addTimeDerivativeToJacobianParticleShell(jac, idxr, alpha, parType);
			// Iterator jac has already been advanced to next shell
		}
	}

	// The bulk balance is extensive, so its mass matrix row of component i is
	//   V^l \dot{c}_i + c_i \dot{V}^l
	double const* const c = y + idxr.offsetC();
	const double v = y[idxr.offsetV()];

	const auto& cm = _reaction.conservedMoieties("liquid");
	if (cm.isEnabled() && cm.numEquilibriumReactions() > 0)
	{
		const auto& L = cm.conservedMoietyMatrix();
		const unsigned int nMoieties = cm.numMoieties();

		for (unsigned int m = 0; m < nMoieties; ++m)
		{
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				_globalJacDisc.coeffRef(idxr.offsetC() + m, idxr.offsetC() + comp) += alpha * v * L(m, comp);

			double volumeEntry = 0.0;
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				volumeEntry += L(m, comp) * c[comp];

			_globalJacDisc.coeffRef(idxr.offsetC() + m, idxr.offsetV()) += alpha * volumeEntry;
		}
		// No mass-matrix entries for equilibrium rows
	}
	else
	{
		const int volumeOffset = idxr.offsetV() - idxr.offsetC();
		linalg::BandedEigenSparseRowIterator jac(_globalJacDisc, idxr.offsetC());

		for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++jac)
		{
			jac[0] += alpha * v;
			jac[volumeOffset - static_cast<int>(comp)] += alpha * c[comp];
		}
	}

	// Volume: \dot{V} - F_in + F_out + F_filter = 0
	_globalJacDisc.coeffRef(idxr.offsetV(), idxr.offsetV()) += alpha;
}

void StirredTankModel::addTimeDerivativeToJacobianParticleShell(linalg::BandedEigenSparseRowIterator& jac, const Indexer& idxr, double alpha, unsigned int parType)
{
	parts::cell::addTimeDerivativeToJacobianParticleShell<linalg::BandedEigenSparseRowIterator, true>(jac, alpha, static_cast<double>(_particles[parType]->getPorosity()), _disc.nComp, _disc.nBound + _disc.nComp * parType,
		_particles[parType]->getPoreAccessFactor(), _disc.strideBound[parType], _disc.boundOffset + _disc.nComp * parType, _binding[parType]->reactionQuasiStationarity());
}

}  // namespace model

}  // namespace cadet

#include "model/StirredTankModel-InitialConditions.cpp"
