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

#include "model/ColumnModel1D.hpp"
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

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <algorithm>
#include <functional>
#include <numeric>
#include <iterator>


#include "ParallelSupport.hpp"
#ifdef CADET_PARALLELIZE
	#include <tbb/parallel_for.h>
#endif

namespace cadet
{

namespace model
{

constexpr double SurfVolRatioSphere = 3.0;
constexpr double SurfVolRatioCylinder = 2.0;
constexpr double SurfVolRatioSlab = 1.0;


ColumnModel1D::ColumnModel1D(UnitOpIdx unitOpIdx) : UnitOperationBase(unitOpIdx),
	_globalJac(), _globalJacDisc(), _jacInlet(),
	_analyticJac(true), _jacobianAdDirs(0), _factorizeJacobian(false), _tempState(nullptr),
	_initC(0), _initCp(0), _initCs(0), _initState(0), _initStateDot(0)
{
}

ColumnModel1D::~ColumnModel1D() CADET_NOEXCEPT
{
	delete[] _tempState;

	for (IParticleModel* pm : _particles)
		delete pm;

	_particles.clear();

	_binding.clear(); // binding models are deleted in the respective particle model
	//_dynReaction.clear(); // particle reaction models are deleted in the respective particle model

	_reaction.clearDynamicReactionModels();
	delete _linearSolver;
}

unsigned int ColumnModel1D::numDofs() const CADET_NOEXCEPT
{
	// Column bulk DOFs: nPoints * nComp
	// Particle DOFs: nPoints * nParType particles each having nComp (liquid phase) + sum boundStates (solid phase) DOFs
	//                in each shell; there are nParCell shells for each particle type
	// Inlet DOFs: nComp
	return _disc.nPoints * _disc.nComp + _disc.parTypeOffset[_disc.nParType] + _disc.nComp;
}

unsigned int ColumnModel1D::numPureDofs() const CADET_NOEXCEPT
{
	// Column bulk DOFs: nPoints * nComp
	// Particle DOFs: nPoints particles each having nComp (liquid phase) + sum boundStates (solid phase) DOFs
	//                in each shell; there are nPar shells
	return _disc.nPoints * _disc.nComp  + _disc.parTypeOffset[_disc.nParType];
}


bool ColumnModel1D::usesAD() const CADET_NOEXCEPT
{
#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	// We always need AD if we want to check the analytical Jacobian
	return true;
#else
	// We only need AD if we are not computing the Jacobian analytically
	return !_analyticJac;
#endif
}

bool ColumnModel1D::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper)
{
	// Read unit type as name to allow model configuration via GRM, LRMP. Here, the particle types are set correspondingly
	std::string unitName = paramProvider.getString("UNIT_TYPE");

	// ==== Read discretization
	_disc.nComp = paramProvider.getInt("NCOMP");

	_disc.nParType = paramProvider.exists("NPARTYPE") ? paramProvider.getInt("NPARTYPE") : 0;
	if (_disc.nParType < 0)
		throw InvalidParameterException("Number of particle types must be >= 0!");

	if (_disc.nParType == 0 && paramProvider.exists("particle_type_000"))
		throw InvalidParameterException("NPARTYPE is set to 0, but group particle_type_000 exists.");

	paramProvider.pushScope("discretization");

	const bool firstConfigCall = _tempState == nullptr; // used to not multiply allocate memory

	if (firstConfigCall)
		_linearSolver = cadet::linalg::setLinearSolver(paramProvider.exists("LINEAR_SOLVER") ? paramProvider.getString("LINEAR_SOLVER") : "SparseLU");

	if (paramProvider.getString("SPATIAL_METHOD") != "DG" && _disc.nParType == 0)
		throw InvalidParameterException("Column with NPARTYPE = 0 is only available for DG discretization");

	if (paramProvider.exists("POLYDEG"))
		_disc.polyDeg = paramProvider.getInt("POLYDEG");
	else
		_disc.polyDeg = 4u; // default value
	if (_disc.polyDeg < 1)
		throw InvalidParameterException("Polynomial degree must be at least 1!");
	else if (_disc.polyDeg < 3)
		LOG(Warning) << "Polynomial degree > 2 in bulk discretization (cf. POLYDEG) is always recommended for performance reasons.";

	_disc.nNodes = _disc.polyDeg + 1;

	if (paramProvider.exists("NELEM"))
		_disc.nElem = paramProvider.getInt("NELEM");
	else if (paramProvider.exists("NCOL"))
		_disc.nElem = std::max(1u, paramProvider.getInt("NCOL") / _disc.nNodes); // number of elements is rounded down
	else
		throw InvalidParameterException("Specify field NELEM (or NCOL)");

	if (_disc.nElem < 1)
		throw InvalidParameterException("Number of column elements must be at least 1!");

	_disc.nPoints = _disc.nNodes * _disc.nElem;

	int polynomial_integration_mode = 0;
	if (paramProvider.exists("EXACT_INTEGRATION"))
		polynomial_integration_mode = paramProvider.getInt("EXACT_INTEGRATION");
	_disc.exactInt = static_cast<bool>(polynomial_integration_mode); // only integration mode 0 applies the inexact collocated diagonal LGL mass matrix

	paramProvider.popScope();

	// Create and configure particle model
	Indexer idxr(_disc);
	_particles = std::vector<IParticleModel*>(_disc.nParType, nullptr);

	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		const std::string parGroup = "particle_type_" + std::string(3 - std::to_string(parType).length(), '0') + std::to_string(parType);
		if (!paramProvider.exists(parGroup))
			throw InvalidParameterException("Unit type was specified as " + unitName + ", but group " + parGroup + " is missing");
		paramProvider.pushScope(parGroup);

		if (unitName == "COLUMN_MODEL_1D")
		{
			const bool filmDiffusion = paramProvider.getBool("HAS_FILM_DIFFUSION");
			const bool poreDiffusion = paramProvider.exists("HAS_PORE_DIFFUSION") ? paramProvider.getBool("HAS_PORE_DIFFUSION") : false;
			const bool surfaceDiffusion = paramProvider.exists("HAS_SURFACE_DIFFUSION") ? paramProvider.getBool("HAS_SURFACE_DIFFUSION") : false;
			const std::string particleType = ParticleModel(filmDiffusion, poreDiffusion, surfaceDiffusion).getParticleTransportType();

			_particles[parType] = helper.createParticleModel(particleType);

			if (!_particles[parType])
				throw InvalidParameterException("Unknown particle model " + particleType);
		}
		else
		{
			std::string particleType = "NONE";
			if (paramProvider.exists("HAS_FILM_DIFFUSION"))
			{
				const bool filmDiffusion = paramProvider.getBool("HAS_FILM_DIFFUSION");
				const bool poreDiffusion = paramProvider.exists("HAS_PORE_DIFFUSION") ? paramProvider.getBool("HAS_PORE_DIFFUSION") : false;
				const bool surfaceDiffusion = paramProvider.exists("HAS_SURFACE_DIFFUSION") ? paramProvider.getBool("HAS_SURFACE_DIFFUSION") : false;
				particleType = ParticleModel(filmDiffusion, poreDiffusion, surfaceDiffusion).getParticleTransportType();
			}

			if (unitName == "GENERAL_RATE_MODEL")
			{
				particleType = particleType == "NONE" ? "GENERAL_RATE_PARTICLE" : particleType;

				if (particleType == "GENERAL_RATE_PARTICLE")
					_particles[parType] = helper.createParticleModel("GENERAL_RATE_PARTICLE");
				else
					throw InvalidParameterException("Unit type was specified as " + unitName + ", which is inconsistent with specified particle model " + particleType);
			}
			else if (unitName == "LUMPED_RATE_MODEL_WITH_PORES")
			{
				particleType = particleType == "NONE" ? "HOMOGENEOUS_PARTICLE" : particleType;

				if (particleType == "HOMOGENEOUS_PARTICLE")
					_particles[parType] = helper.createParticleModel("HOMOGENEOUS_PARTICLE");
				else
					throw InvalidParameterException("Unit type was specified as " + unitName + ", which is inconsistent with specified particle model " + particleType);
			}
			else
				throw InvalidParameterException("Failed to configure unit type " + unitName);
		}

		paramProvider.popScope(); // particle_type_xxx
	}

	if (std::any_of(_particles.begin(), _particles.end(), [](IParticleModel* pm) { return !static_cast<bool>(pm); }))
		throw InvalidParameterException("Particle model configuration failed");

	bool particleConfSuccess = true;
	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		particleConfSuccess = particleConfSuccess && _particles[parType]->configureModelDiscretization(paramProvider, helper, _disc.nComp, parType, _disc.nParType, idxr.strideColComp());
	}
	paramProvider.pushScope("discretization");

	if (firstConfigCall)
	{
		_disc.nParPoints = new unsigned int[_disc.nParType];
		for (int type = 0; type < _disc.nParType; type++)
		{
			_disc.nParPoints[type] = _particles[type]->nDiscPoints();
		}
	}

	_disc.nBound = new unsigned int[_disc.nParType * _disc.nComp];
	for (int parType = 0; parType < _disc.nParType; parType++)
		for (int comp = 0; comp < _disc.nComp; comp++)
			_disc.nBound[parType * _disc.nComp + comp] = _particles[parType]->nBound()[comp];

	const unsigned int nTotalBound = std::accumulate(_disc.nBound, _disc.nBound + _disc.nComp * _disc.nParType, 0u);

	// Precompute offsets and total number of bound states (DOFs in solid phase)
	if (firstConfigCall)
	{
		_disc.boundOffset = new unsigned int[_disc.nComp * _disc.nParType];
		_disc.strideBound = new unsigned int[_disc.nParType + 1];
		_disc.nBoundBeforeType = new unsigned int[std::max(_disc.nParType, 1u)];
	}
	_disc.strideBound[_disc.nParType] = nTotalBound;
	_disc.nBoundBeforeType[0] = 0;
	for (unsigned int j = 0; j < _disc.nParType; ++j)
	{
		unsigned int* const ptrOffset = _disc.boundOffset + j * _disc.nComp;
		unsigned int* const ptrBound = _disc.nBound + j * _disc.nComp;

		ptrOffset[0] = 0;
		for (unsigned int i = 1; i < _disc.nComp; ++i)
		{
			ptrOffset[i] = ptrOffset[i - 1] + ptrBound[i - 1];
		}
		_disc.strideBound[j] = ptrOffset[_disc.nComp - 1] + ptrBound[_disc.nComp - 1];

		if (j != _disc.nParType - 1)
			_disc.nBoundBeforeType[j + 1] = _disc.nBoundBeforeType[j] + _disc.strideBound[j];
	}

	// Precompute offsets of particle type DOFs
	if (firstConfigCall)
	{
		_disc.parTypeOffset = new unsigned int[_disc.nParType + 1];
	}
	_disc.parTypeOffset[0] = 0;
	unsigned int nTotalParPoints = 0;
	for (unsigned int j = 1; j < _disc.nParType + 1; ++j)
	{
		_disc.parTypeOffset[j] = _disc.parTypeOffset[j-1] + (_disc.nComp + _disc.strideBound[j-1]) * _disc.nParPoints[j-1] * _disc.nPoints;
		nTotalParPoints += _disc.nParPoints[j-1];
	}

	// Determine whether analytic Jacobian should be used but don't set it right now.
	// We need to setup Jacobian matrices first.
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	const bool analyticJac = paramProvider.getBool("USE_ANALYTIC_JACOBIAN");
#else
	const bool analyticJac = false;
#endif

	// Allocate space for initial conditions
	_initC.resize(_disc.nComp);
	_initCp.resize(_disc.nComp * _disc.nParType);
	_initCs.resize(nTotalBound);

	// Create nonlinear solver for consistent initialization
	configureNonlinearSolver(paramProvider);

	paramProvider.popScope();

	// ==== Construct and configure convection dispersion operator

	unsigned int strideColNode = _disc.nComp;
	const bool transportSuccess = _convDispOp.configureModelDiscretization(paramProvider, helper, _disc.nComp, polynomial_integration_mode, _disc.nElem, _disc.polyDeg, strideColNode);

	_disc.curSection = -1;

	// ==== Construct and configure dynamic reaction model
	bool reactionConfSuccess = true;
	_reaction.clearDynamicReactionModels();

	// Bulk liquid phase reactions
	if (paramProvider.exists("NREAC_LIQUID"))
	{
		int nReactions = paramProvider.getInt("NREAC_LIQUID");
		reactionConfSuccess = _reaction.configureDiscretization("liquid",
			nReactions,
			_disc.nComp,
			_disc.nBound,
			_disc.boundOffset,
			paramProvider,
			helper) && reactionConfSuccess;
	}
	else
	{
		_reaction.empty();
	}

	// ==== Construct and configure binding and particle reaction -> done in particle model, only pointers are copied here.
	_binding = std::vector<IBindingModel*>(_disc.nParType, nullptr);
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		_binding[parType] = _particles[parType]->getBinding();

		// Check if binding and reaction particle type dependence is the same for all particle types
		if (parType > 0)
		{
			if (_binding[parType])
			{
				if (_singleBinding != !_particles[parType]->bindingParDep())
					throw InvalidParameterException("Binding particle type dependence must be the same for all particle types, check field BINDING_PARTYPE_DEPENDENT");
			}
		}
		else // if no particle reaction or binding exists in first particle type, default to single mode
		{
			_singleBinding = _binding[parType] ? !_particles[parType]->bindingParDep() : true;
		}
	}

	// Allocate memory
	if (firstConfigCall)
		_tempState = new double[numDofs()];

	if (_disc.exactInt)
		_jacInlet.resize(_disc.nNodes, 1); // first cell depends on inlet concentration (same for every component)
	else
		_jacInlet.resize(1, 1); // first node depends on inlet concentration (same for every component)

	// set jacobian pattern
	_globalJacDisc.resize(numDofs(), numDofs());
	_globalJac.resize(numDofs(), numDofs());
	// pattern is set in configure(), after surface diffusion is read

	// Set whether analytic Jacobian is used
	useAnalyticJacobian(analyticJac);

	return transportSuccess && particleConfSuccess && reactionConfSuccess;
}
 
bool ColumnModel1D::configure(IParameterProvider& paramProvider)
{
	_parameters.clear();

	const bool transportSuccess = _convDispOp.configure(_unitOpIdx, paramProvider, _parameters);

	// Read geometry parameters handled by unit operation
	_colPorosity = paramProvider.getDouble("COL_POROSITY");

	// Check whether PAR_TYPE_VOLFRAC is required or not
	if ((_disc.nParType > 1) && !paramProvider.exists("PAR_TYPE_VOLFRAC"))
		throw InvalidParameterException("The required parameter \"PAR_TYPE_VOLFRAC\" was not found");

	_axiallyConstantParTypeVolFrac = true;

	if (_disc.nParType > 1)
	{
		readScalarParameterOrArray(_parTypeVolFrac, paramProvider, "PAR_TYPE_VOLFRAC", 1);
		if (_parTypeVolFrac.size() == _disc.nParType)
		{

			// Expand to all axial cells
			_parTypeVolFrac.resize(_disc.nPoints * _disc.nParType, 1.0);
			for (unsigned int i = 1; i < _disc.nPoints; ++i)
				std::copy(_parTypeVolFrac.begin(), _parTypeVolFrac.begin() + _disc.nParType, _parTypeVolFrac.begin() + _disc.nParType * i);
		}
		else
			_axiallyConstantParTypeVolFrac = false;

		if (_disc.nParType * _disc.nPoints != _parTypeVolFrac.size())
			throw InvalidParameterException("Number of elements in field PAR_TYPE_VOLFRAC does not match number of particle types times number of axial cells");

		// Check that particle volume fractions sum to 1.0
		for (unsigned int i = 0; i < _disc.nPoints; ++i)
		{
			const double volFracSum = std::accumulate(_parTypeVolFrac.begin() + i * _disc.nParType, _parTypeVolFrac.begin() + (i + 1) * _disc.nParType, 0.0,
				[](double a, const active& b) -> double { return a + static_cast<double>(b); });
			if (std::abs(1.0 - volFracSum) > 1e-10)
				throw InvalidParameterException("Sum of field PAR_TYPE_VOLFRAC differs from 1.0 (is " + std::to_string(volFracSum) + ") in axial cell " + std::to_string(i));
		}

		if (_axiallyConstantParTypeVolFrac)
		{
			// Register only the first nParType items
			for (unsigned int i = 0; i < _disc.nParType; ++i)
				_parameters[makeParamId(hashString("PAR_TYPE_VOLFRAC"), _unitOpIdx, CompIndep, i, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parTypeVolFrac[i];
		}
		else
			registerParam2DArray(_parameters, _parTypeVolFrac, [=](bool multi, unsigned cell, unsigned int type) { return makeParamId(hashString("PAR_TYPE_VOLFRAC"), _unitOpIdx, CompIndep, type, BoundStateIndep, ReactionIndep, cell); }, _disc.nParType);
	}
	else if (_disc.nParType == 1)
		_parTypeVolFrac = std::vector<active>(_disc.nPoints, 1.0);

	// Read vectorial parameters (which may also be section dependent; transport)
	// Add parameters to map
	_parameters[makeParamId(hashString("COL_POROSITY"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_colPorosity;

	// Register initial conditions parameters
	registerParam1DArray(_parameters, _initC, [=](bool multi, unsigned int comp) { return makeParamId(hashString("INIT_C"), _unitOpIdx, comp, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep); });

	if (_singleBinding)
	{
		for (unsigned int c = 0; c < _disc.nComp; ++c)
			_parameters[makeParamId(hashString("INIT_CP"), _unitOpIdx, c, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_initCp[c];
	}
	else
		registerParam2DArray(_parameters, _initCp, [=](bool multi, unsigned int type, unsigned int comp) { return makeParamId(hashString("INIT_CP"), _unitOpIdx, comp, type, BoundStateIndep, ReactionIndep, SectionIndep); }, _disc.nComp);


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

	// Reconfigure particle model
	bool particleConfSuccess = true;
	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		particleConfSuccess = particleConfSuccess && _particles[parType]->configure(_unitOpIdx, paramProvider, _parameters, _disc.nParType, _disc.nBoundBeforeType, _disc.strideBound[_disc.nParType]);
	}

	// Reconfigure bulk liquid reaction model
	bool dynReactionConfSuccess = true;
	if (paramProvider.exists("NREAC_LIQUID"))
		dynReactionConfSuccess = _reaction.configure("liquid", 0, _unitOpIdx, paramProvider) && dynReactionConfSuccess;

	// jaobian pattern set after binding and particle surface diffusion are configured
	bool hasBulkReaction = false;
	if (_reaction.getDynReactionVector("liquid")[0] != nullptr)
		hasBulkReaction = true;

	setJacobianPattern(_globalJac, 0, hasBulkReaction);
	_globalJacDisc = _globalJac;
	// the solver repetitively solves the linear system with a static pattern of the jacobian (set above). 
	// The goal of analyzePattern() is to reorder the nonzero elements of the matrix, such that the factorization step creates less fill-in
	_linearSolver->analyzePattern(_globalJacDisc.block(_disc.nComp, _disc.nComp, numPureDofs(), numPureDofs()));

	return transportSuccess && particleConfSuccess && dynReactionConfSuccess;
}

unsigned int ColumnModel1D::threadLocalMemorySize() const CADET_NOEXCEPT
{
	LinearMemorySizer lms;

	// Memory for residualImpl()
	for (unsigned int i = 0; i < _disc.nParType; ++i)
	{
		if (_binding[i] && _binding[i]->requiresWorkspace())
			lms.fitBlock(_binding[i]->workspaceSize(_disc.nComp, _disc.strideBound[i], _disc.nBound + i * _disc.nComp));
	}
	// Bulk reaction
	_reaction.setWorkspaceRequirements("liquid", _disc.nComp, 0, lms);

	const unsigned int maxStrideBound = *std::max_element(_disc.strideBound, _disc.strideBound + _disc.nParType);
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

	// Memory for consistentInitialSensitivity
	lms.add<double>(_disc.nComp + maxStrideBound);
	lms.add<double>(maxStrideBound);
	lms.commit();

	return lms.bufferSize();
}

unsigned int ColumnModel1D::numAdDirsForJacobian() const CADET_NOEXCEPT
{
	// The global DG Jacobian is banded around the main diagonal and has additional (also banded) entries for film diffusion.
	// To feasibly seed and reconstruct the Jacobian, we create dedicated active directions for the bulk and each particle type (see @ todo)
	Indexer idxr(_disc);
	
	int sumParBandwidth = 0;
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		sumParBandwidth += idxr.strideParBlock(type);
	}

	return _convDispOp.requiredADdirs() + sumParBandwidth;
}

void ColumnModel1D::useAnalyticJacobian(const bool analyticJac)
{
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	_analyticJac = analyticJac;
	if (!_analyticJac)
		_jacobianAdDirs = numAdDirsForJacobian();
	else
		_jacobianAdDirs = 0;
#else
	// If CADET_CHECK_ANALYTIC_JACOBIAN is active, we always enable AD for comparison and use it in simulation
	_analyticJac = false;
	_jacobianAdDirs = numAdDirsForJacobian();
#endif
}

void ColumnModel1D::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, const ConstSimulationState& simState, const AdJacobianParams& adJac)
{
	Indexer idxr(_disc);

	// todo: only reset jacobian pattern if it changes, i.e. once in configuration and then only for changes in SurfDiff+kinetic binding.
	bool hasReaction = _reaction.getDynReactionVector("liquid")[0];
	setJacobianPattern(_globalJac, 0, hasReaction);
	_globalJacDisc = _globalJac;

	_convDispOp.notifyDiscontinuousSectionTransition(t, secIdx, _jacInlet);

	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		_particles[parType]->notifyDiscontinuousSectionTransition(t, secIdx);
	}

	_disc.curSection = secIdx;
}

void ColumnModel1D::setFlowRates(active const* in, active const* out) CADET_NOEXCEPT
{
	_convDispOp.setFlowRates(in[0], out[0], _colPorosity);
}

void ColumnModel1D::reportSolution(ISolutionRecorder& recorder, double const* const solution) const
{
	Exporter expr(_disc, *this, solution);
	recorder.beginUnitOperation(_unitOpIdx, *this, expr);
	recorder.endUnitOperation();
}

void ColumnModel1D::reportSolutionStructure(ISolutionRecorder& recorder) const
{
	Exporter expr(_disc, *this, nullptr);
	recorder.unitOperationStructure(_unitOpIdx, *this, expr);
}

unsigned int ColumnModel1D::requiredADdirs() const CADET_NOEXCEPT
{
	const unsigned int numDirsBinding = maxBindingAdDirs();
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	return numDirsBinding + _jacobianAdDirs;
#else
	// If CADET_CHECK_ANALYTIC_JACOBIAN is active, we always need the AD directions for the Jacobian
	return numDirsBinding + numAdDirsForJacobian();
#endif
}

void ColumnModel1D::prepareADvectors(const AdJacobianParams& adJac) const
{
	// Early out if AD is disabled
	if (!adJac.adY)
		return;

	Indexer idxr(_disc);

	// The global DG Jacobian is banded around the main diagonal and has additional (also banded, but offset) entries for film diffusion,
	// i.e. banded AD vector seeding is not sufficient (as it is for the FV Jacobians, see @puttmann2016 and the DG LRM Jacobian).
	// The compressed vectorial AD seeding and Jacobian construction is described in the following.
	// The global DG Jacobian is banded around the main diagonal and has additional (also banded) entries for film diffusion.
	// To feasibly seed and reconstruct the Jacobian (we need information for decompression), we create dedicated active directions for
	// the bulk and each particle type. Particle AD directions are treated as dense (per particle block) since only nCells > 6 would
	// justify band compression, which rarely ever happens.

	// We begin by seeding the (banded around main diagonal) bulk Jacobian block
	// We have differing Jacobian structures for exact integration and collocation DG scheme, i.e. we need different seed vectors
	// collocation DG: 2 * N_n * (N_c + N_q) + 1 = total bandwidth (main diagonal entries maximally depend on the next and last N_n liquid phase entries of same component)
	//    ex. int. DG: 4 * N_n * (N_c + N_q) + 1 = total bandwidth (main diagonal entries maximally depend on the next and last 2*N_n liquid phase entries of same component)
	const int lowerBandwidth = (_disc.exactInt) ? 2 * _disc.nNodes * idxr.strideColNode() : _disc.nNodes * idxr.strideColNode();
	const int upperBandwidth = lowerBandwidth;
	const int bulkRows = idxr.offsetCp() - idxr.offsetC();
	ad::prepareAdVectorSeedsForBandMatrix(adJac.adY + _disc.nComp, adJac.adDirOffset, bulkRows, lowerBandwidth, upperBandwidth, lowerBandwidth);

	// We now seed the particle Jacobian blocks using the individual AD directions for each particle type.
	unsigned int adDirOffset = adJac.adDirOffset + _convDispOp.requiredADdirs();

	for (unsigned int type = 0; type < _disc.nParType; type++)
	{
		for (unsigned int parBlock = 0; parBlock < _disc.nPoints; parBlock++)
		{
			// move adVec pointer to start of current particle block 
			active* _adVec = adJac.adY + idxr.offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ parBlock });

			for (int eq = 0; eq < idxr.strideParBlock(type); ++eq)
			{
				// Clear previously set directions
				_adVec[eq].fillADValue(adJac.adDirOffset, 0.0);
				// Set direction
				_adVec[eq].setADValue(adDirOffset + eq, 1.0);
			}
		}
		if (type < _disc.nParType - 1u) // move to dedicated DoFs of next particle type
			adDirOffset += idxr.strideParBlock(type);
	}
}
/**
 * @brief Extracts the system Jacobian from band compressed AD seed vectors
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
void ColumnModel1D::extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset)
{
	Indexer idxr(_disc);

	const active* const adVec = adRes + idxr.offsetC();

	/* Extract bulk phase equations entries */
	const int lowerBandwidth = (_disc.exactInt) ? 2 * _disc.nNodes * idxr.strideColNode() : _disc.nNodes * idxr.strideColNode();
	const int upperBandwidth = lowerBandwidth;
	const int stride = lowerBandwidth + 1 + upperBandwidth;
	int diagDir = lowerBandwidth;
	const int bulkDoFs = idxr.offsetCp() - idxr.offsetC();
	const int eqOffset = 0;
	const int matOffset = idxr.offsetC();
	ad::extractBandedBlockEigenJacobianFromAd(adVec, adDirOffset, diagDir, lowerBandwidth, upperBandwidth, eqOffset, bulkDoFs, _globalJac, matOffset);

	/* Handle particle liquid and solid phase equations entries */
	// Read particle Jacobian entries from dedicated AD directions
	int offsetParticleTypeDirs = adDirOffset + _convDispOp.requiredADdirs();

	for (unsigned int type = 0; type < _disc.nParType; type++)
	{
		for (unsigned int par = 0; par < _disc.nPoints; par++)
		{
			const int eqOffset_res = idxr.offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ par });
			const int eqOffset_mat = idxr.offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ par });
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

	/* Add analytically derived flux entries (only those that are part of the outlier bands) */
	// todo extract these entries instead of analytical calculation?
	for (unsigned int parType = 0; parType < _disc.nParType; parType++)
	{
		_particles[parType]->calcFilmDiffJacobian(_disc.curSection, idxr.offsetCp(ParticleTypeIndex{static_cast<unsigned int>(parType)}), idxr.offsetC(), _disc.nPoints, _disc.nParType, static_cast<double>(_colPorosity), &_parTypeVolFrac[0], _globalJac, true);
	}

	if (!_globalJac.isCompressed())
		_globalJac.makeCompressed();
}

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN

/**
 * @brief Compares the analytical Jacobian with a Jacobian derived by AD
 * @details The analytical Jacobian is assumed to be stored in the corresponding band matrices.
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
void ColumnModel1D::checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const
{
	Indexer idxr(_disc);

	const active* const adVec = adRes + idxr.offsetC();

	/* Extract bulk phase equations entries */
	const int lowerBandwidth = (_disc.exactInt) ? 2 * _disc.nNodes * idxr.strideColNode() : _disc.nNodes * idxr.strideColNode();
	const int upperBandwidth = lowerBandwidth;
	const int stride = lowerBandwidth + 1 + upperBandwidth;
	int diagDir = lowerBandwidth;
	const int bulkDoFs = idxr.offsetCp() - idxr.offsetC();
	const int eqOffset = 0;
	const int matOffset = idxr.offsetC();

	double JacMaxError = 0.0;
	int row = 0;
	int col = 0;

	for (int eq = eqOffset; eq < eqOffset + bulkDoFs; ++eq)
	{
		// Start with lowest subdiagonal and stay in the range of the columns:
		// diagDir might be too big for the matrix and, hence, dir ranges between
		// diagDir - lowerBandwidth <= dir <= diagDir + upperBandwidth
		int dir = diagDir - lowerBandwidth + eq % stride;

		// Loop over diagonals
		for (int diag = 0; diag < stride; ++diag)
		{
			if (eq - lowerBandwidth + diag >= eqOffset && // left block boundary
				eq - lowerBandwidth + diag < eqOffset + bulkDoFs && // right block boundary
				adVec[eq].getADValue(adDirOffset + dir) != 0.0 // keep pattern
				)
			{
				row = matOffset + eq;
				col = matOffset + eq - lowerBandwidth + diag;

				const double localError = std::abs(_globalJac.coeff(row, col) - adVec[eq].getADValue(adDirOffset + dir));

				if (localError > 1e-14)
				{
					JacMaxError = std::max(JacMaxError, localError);
					if (row < idxr.offsetCp() && col < idxr.offsetCp())
						LOG(Debug) << "Error in bulk Jacobian: " << localError << " at row, col: " << row << ", " << col;
					else if ((row >= idxr.offsetCp() && col < idxr.offsetCp()) || (row < idxr.offsetCp() && col >= idxr.offsetCp()))
						LOG(Debug) << "Error in film diffusion Jacobian: " << localError << " at row, col: " << row << ", " << col;
					else
						LOG(Debug) << "Error in particle Jacobian: " << localError << " at row, col: " << row << ", " << col;

					LOG(Debug) << "AD Jacobian value: " << adVec[eq].getADValue(adDirOffset + dir);
					LOG(Debug) << "Analytical Jacobian value: " << _globalJac.coeff(row, col);
				}
			}
			// Wrap around at end of row and jump to lowest subdiagonal
			if (dir == diagDir + upperBandwidth)
				dir = diagDir - lowerBandwidth;
			else
				++dir;
		}
	}

	/* Handle particle liquid and solid phase equations entries */
	// Read particle Jacobian entries from dedicated AD directions
	int offsetParticleTypeDirs = adDirOffset + _convDispOp.requiredADdirs();

	for (unsigned int type = 0; type < _disc.nParType; type++)
	{
		for (unsigned int par = 0; par < _disc.nPoints; par++)
		{
			const int eqOffset_res = idxr.offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ par });
			const int eqOffset_mat = idxr.offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ par });
			for (unsigned int phase = 0; phase < idxr.strideParBlock(type); phase++)
			{
				for (unsigned int phaseTo = 0; phaseTo < idxr.strideParBlock(type); phaseTo++)
				{
					row = eqOffset_mat + phase;
					col = eqOffset_mat + phaseTo;

					const double localError = std::abs(_globalJac.coeff(row, col) - adRes[eqOffset_res + phase].getADValue(offsetParticleTypeDirs + phaseTo));

					LOG(Debug) << "Error in particle type " + std::to_string(type) + " Jacobian: " << localError << " at row, col: " << row << ", " << col;
					LOG(Debug) << "AD Jacobian value: " << adRes[eqOffset_res + phase].getADValue(offsetParticleTypeDirs + phaseTo);
					LOG(Debug) << "Analytical Jacobian value: " << _globalJac.coeff(row, col);

				}
			}
		}
		offsetParticleTypeDirs += idxr.strideParBlock(type);
	}

	// Analytical Jacobain is always used for film diffusion, nothing to do here
}

#endif

int ColumnModel1D::jacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	_factorizeJacobian = true;

	if (_analyticJac)
		return residualImpl<double, double, double, true, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
	else
		return residualWithJacobian(simTime, ConstSimulationState{ simState.vecStateY, nullptr }, nullptr, adJac, threadLocalMem);
}

int ColumnModel1D::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	// Evaluate residual do not compute Jacobian or parameter sensitivities
	return residualImpl<double, double, double, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
}

int ColumnModel1D::residualWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	//FDJac = calcFDJacobian(static_cast<const double*>(simState.vecStateY), static_cast<const double*>(simState.vecStateYdot), simTime, threadLocalMem, 2.0); // debug code

	// Evaluate residual, use AD for Jacobian if required but do not evaluate parameter derivatives
	return residual(simTime, simState, res, adJac, threadLocalMem, true, false);
}

int ColumnModel1D::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res,
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

				// Copy AD residuals to original residuals vector
				if (res)
					ad::copyFromAd(adJac.adRes, res, numDofs());

				return retCode;
			}
			else
				return residualImpl<double, double, double, true>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
		}
		else
		{
			// Copy over state vector to AD state vector (without changing directional values to keep seed vectors)
			// and initialize residuals with zero (also resetting directional values)
			ad::copyToAd(simState.vecStateY, adJac.adY, numDofs());
			// @todo Check if this is necessary
			ad::resetAd(adJac.adRes, numDofs());

			// Evaluate with AD enabled
			int retCode = 0;
			if (paramSensitivity)
				retCode = residualImpl<active, active, active, false>(simTime.t, simTime.secIdx, adJac.adY, simState.vecStateYdot, adJac.adRes, threadLocalMem);
			else
				retCode = residualImpl<active, active, double, false>(simTime.t, simTime.secIdx, adJac.adY, simState.vecStateYdot, adJac.adRes, threadLocalMem);

			// Copy AD residuals to original residuals vector
			if (res)
				ad::copyFromAd(adJac.adRes, res, numDofs());

			// Extract Jacobian
			extractJacobianFromAD(adJac.adRes, adJac.adDirOffset);

			return retCode;
		}
#else
		// Copy over state vector to AD state vector (without changing directional values to keep seed vectors)
		// and initialize residuals with zero (also resetting directional values)
		ad::copyToAd(simState.vecStateY, adJac.adY, numDofs());
		// @todo Check if this is necessary
		ad::resetAd(adJac.adRes, numDofs());

		// Evaluate with AD enabled
		int retCode = 0;
		if (paramSensitivity)
			retCode = residualImpl<active, active, active, false>(simTime.t, simTime.secIdx, adJac.adY, simState.vecStateYdot, adJac.adRes, threadLocalMem);
		else
			retCode = residualImpl<active, active, double, false>(simTime.t, simTime.secIdx, adJac.adY, simState.vecStateYdot, adJac.adRes, threadLocalMem);

		// Only do comparison if we have a residuals vector (which is not always the case)
		if (res)
		{
			// Evaluate with analytical Jacobian which is stored in the band matrices
			retCode = residualImpl<double, double, double, true>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);

			// Compare AD with anaytic Jacobian
			checkAnalyticJacobianAgainstAd(adJac.adRes, adJac.adDirOffset);
		}

		// Extract Jacobian
		extractJacobianFromAD(adJac.adRes, adJac.adDirOffset);

		return retCode;
#endif
	}
	else
	{
		if (paramSensitivity)
		{
			// initialize residuals with zero
			// @todo Check if this is necessary
			ad::resetAd(adJac.adRes, numDofs());

			const int retCode = residualImpl<double, active, active, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, adJac.adRes, threadLocalMem);

			// Copy AD residuals to original residuals vector
			if (res)
				ad::copyFromAd(adJac.adRes, res, numDofs());

			return retCode;
		}
		else
			return residualImpl<double, double, double, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
	}
}

template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
int ColumnModel1D::residualImpl(double t, unsigned int secIdx, StateType const* const y, double const* const yDot, ResidualType* const res, util::ThreadLocalStorage& threadLocalMem)
{
	if (wantRes)
	{
		Eigen::Map<Eigen::Vector<ResidualType, Dynamic>> resi(res, numDofs());
		resi.setZero();
	}
	if (wantJac)
		_globalJac.coeffs().setZero();

	LinearBufferAllocator tlmAlloc = threadLocalMem.get();
	Indexer idxr(_disc);

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(std::size_t(0), static_cast<std::size_t>(_disc.nPoints * _disc.nParType + 1), [&](std::size_t pblk)
#else
	for (unsigned int pblk = 0; pblk < _disc.nPoints * _disc.nParType + 1; ++pblk)
#endif
	{
		if (cadet_unlikely(pblk == 0))
		{
			if (wantJac)
			{
				// estimate new static (per section) jacobian
				bool success = calcTransportJacobian(secIdx);

				if (cadet_unlikely(!success))
					LOG(Error) << "Jacobian pattern did not fit the analytical transport Jacobian assembly";
			}

			residualBulk<StateType, ResidualType, ParamType, wantJac, wantRes>(t, secIdx, y, yDot, res, threadLocalMem);

			continue;
		}

		const unsigned int parType = (pblk - 1) / _disc.nPoints;
		const unsigned int colNode = (pblk - 1) % _disc.nPoints;

		linalg::BandedEigenSparseRowIterator jacIt;

		if (wantJac)
		{
			jacIt = linalg::BandedEigenSparseRowIterator(
				_globalJac,
				idxr.offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode })
			);
		}

		model::columnPackingParameters packing
		{
			_parTypeVolFrac[parType + _disc.nParType * colNode],
			_colPorosity,
			ColumnPosition{ _convDispOp.relativeCoordinate(colNode), 0.0, 0.0 }
		};


		_particles[parType]->residual(t, secIdx,
			y + idxr.offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode }),
			y + idxr.offsetC() + colNode * idxr.strideColNode(),
			yDot ? yDot + idxr.offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode }) : nullptr,
			res ? res + idxr.offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode }) : nullptr,
			res ? res + idxr.offsetC() + colNode * idxr.strideColNode() : nullptr,
			packing, jacIt, tlmAlloc,
			typename cadet::ParamSens<ParamType>::enabled()
		);
	}

	if (!wantRes)
		return 0;

	BENCH_STOP(_timerResidualPar);

	// Handle inlet DOFs, which are simply copied to the residual
	for (unsigned int i = 0; i < _disc.nComp; ++i)
	{
		res[i] = y[i];
	}

	return 0;
}

template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
int ColumnModel1D::residualBulk(double t, unsigned int secIdx, StateType const* yBase, double const* yDotBase, ResidualType* resBase, util::ThreadLocalStorage& threadLocalMem)
{
	if (wantRes)
		_convDispOp.residual(*this, t, secIdx, yBase, yDotBase, resBase, typename cadet::ParamSens<ParamType>::enabled());

	if (_reaction.getDynReactionVector("liquid").size() == 0)
		return 0;

	Indexer idxr(_disc);
	LinearBufferAllocator tlmAlloc = threadLocalMem.get();
	StateType const* y = yBase + idxr.offsetC();

	if (wantJac && !wantRes) // only compute Jacobian
	{
		for (unsigned int col = 0; col < _disc.nPoints; ++col, y += idxr.strideColNode())
		{
			const ColumnPosition colPos{ _convDispOp.relativeCoordinate(col), 0.0, 0.0 };
			linalg::BandedEigenSparseRowIterator jac(_globalJac, idxr.offsetC() + col * idxr.strideColNode());
			
			for (auto i = 0; i < _reaction.getDynReactionVector("liquid").size(); i++)
			{
				if (!_reaction.getDynReactionVector("liquid")[i])
					continue;

				_reaction.getDynReactionVector("liquid")[i]->analyticJacobianAdd(t, secIdx, colPos, _disc.nComp, reinterpret_cast<double const*>(y), -1.0, jac, tlmAlloc);
			}
		}

		return 0;
	}

	ResidualType* res = resBase + idxr.offsetC();

	for (unsigned int col = 0; col < _disc.nPoints; ++col, y += idxr.strideColNode(), res += idxr.strideColNode())
	{
		const ColumnPosition colPos{ _convDispOp.relativeCoordinate(col), 0.0, 0.0};
		linalg::BandedEigenSparseRowIterator jac(_globalJac, idxr.offsetC() + col * idxr.strideColNode());

		for (auto i = 0; i < _reaction.getDynReactionVector("liquid").size(); i++)
		{
			if (!_reaction.getDynReactionVector("liquid")[i])
				continue;

			_reaction.getDynReactionVector("liquid")[i]->residualFluxAdd(t, secIdx, colPos, _disc.nComp, y, res, -1.0, tlmAlloc);

			if (wantJac)
			{
				// static_cast should be sufficient here, but this statement is also analyzed when wantJac = false
				_reaction.getDynReactionVector("liquid")[i]->analyticJacobianAdd(t, secIdx, colPos, _disc.nComp, reinterpret_cast<double const*>(y), -1.0, jac, tlmAlloc);
			}
		}
	}


	return 0;
}

parts::cell::CellParameters ColumnModel1D::makeCellResidualParams(unsigned int parType, int const* qsReaction) const
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

int ColumnModel1D::residualSensFwdWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidualSens);

	// Evaluate residual for all parameters using AD in vector mode and at the same time update the
	// Jacobian (in one AD run, if analytic Jacobians are disabled)
	return residual(simTime, simState, nullptr, adJac, threadLocalMem, true, true);
}

int ColumnModel1D::residualSensFwdAdOnly(const SimulationTime& simTime, const ConstSimulationState& simState, active* const adRes, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidualSens);

	// Evaluate residual for all parameters using AD in vector mode
	return residualImpl<double, active, active, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, adRes, threadLocalMem);
}

int ColumnModel1D::residualSensFwdCombine(const SimulationTime& simTime, const ConstSimulationState& simState,
	const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS, active const* adRes,
	double* const tmp1, double* const tmp2, double* const tmp3)
{
	BENCH_SCOPE(_timerResidualSens);

	// tmp1 stores result of (dF / dy) * s
	// tmp2 stores result of (dF / dyDot) * sDot

	for (std::size_t param = 0; param < yS.size(); ++param)
	{
		// Directional derivative (dF / dy) * s
		multiplyWithJacobian(SimulationTime{0.0, 0u}, ConstSimulationState{nullptr, nullptr}, yS[param], 1.0, 0.0, tmp1);

		// Directional derivative (dF / dyDot) * sDot
		multiplyWithDerivativeJacobian(SimulationTime{0.0, 0u}, ConstSimulationState{nullptr, nullptr}, ySdot[param], tmp2);

		double* const ptrResS = resS[param];

		BENCH_START(_timerResidualSensPar);

		// Complete sens residual is the sum:
		// TODO: Chunk TBB loop
#ifdef CADET_PARALLELIZE
		tbb::parallel_for(std::size_t(0), static_cast<std::size_t>(numDofs()), [&](std::size_t i)
#else
		for (unsigned int i = 0; i < numDofs(); ++i)
#endif
		{
			ptrResS[i] = tmp1[i] + tmp2[i] + adRes[i].getADValue(param);
		} CADET_PARFOR_END;

		BENCH_STOP(_timerResidualSensPar);
	}

	return 0;
}
/**
 * @brief Multiplies the given vector with the system Jacobian (i.e., @f$ \frac{\partial F}{\partial y}\left(t, y, \dot{y}\right) @f$)
 * @details Actually, the operation @f$ z = \alpha \frac{\partial F}{\partial y} x + \beta z @f$ is performed.
 *
 *          Note that residual() or one of its cousins has to be called with the requested point @f$ (t, y, \dot{y}) @f$ once
 *          before calling multiplyWithJacobian() as this implementation ignores the given @f$ (t, y, \dot{y}) @f$.
 * @param [in] simTime Current simulation time point
 * @param [in] simState Simulation state vectors
 * @param [in] yS Vector @f$ x @f$ that is transformed by the Jacobian @f$ \frac{\partial F}{\partial y} @f$
 * @param [in] alpha Factor @f$ \alpha @f$ in front of @f$ \frac{\partial F}{\partial y} @f$
 * @param [in] beta Factor @f$ \beta @f$ in front of @f$ z @f$
 * @param [in,out] ret Vector @f$ z @f$ which stores the result of the operation
 */
void ColumnModel1D::multiplyWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* yS, double alpha, double beta, double* ret)
{
	Indexer idxr(_disc);

	// Handle identity matrix of inlet DOFs
	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
	{
		ret[comp] = alpha * yS[comp] + beta * ret[comp];
	}

	// Main Jacobian
	Eigen::Map<Eigen::VectorXd> ret_vec(ret + idxr.offsetC(), numPureDofs());
	Eigen::Map<const Eigen::VectorXd> yS_vec(yS + idxr.offsetC(), numPureDofs());
	ret_vec = alpha * _globalJac.block(idxr.offsetC(), idxr.offsetC(), numPureDofs(), numPureDofs()) * yS_vec + beta * ret_vec;

	// Map inlet DOFs to the column inlet (first bulk cells)
	// Inlet at z = 0 for forward flow, at z = L for backward flow.
	unsigned int offInlet = _convDispOp.forwardFlow() ? 0 : (_disc.nElem - 1u) * idxr.strideColCell();

	for (unsigned int comp = 0; comp < _disc.nComp; comp++) {
		for (unsigned int node = 0; node < (_disc.exactInt ? _disc.nNodes : 1); node++) {
			ret[idxr.offsetC() + offInlet + comp * idxr.strideColComp() + node * idxr.strideColNode()] += alpha * _jacInlet(node, 0) * yS[comp];
		}
	}
}

/**
 * @brief Multiplies the time derivative Jacobian @f$ \frac{\partial F}{\partial \dot{y}}\left(t, y, \dot{y}\right) @f$ with a given vector
 * @details The operation @f$ z = \frac{\partial F}{\partial \dot{y}} x @f$ is performed.
 *          The matrix-vector multiplication is performed matrix-free (i.e., no matrix is explicitly formed).
 * @param [in] simTime Current simulation time point
 * @param [in] simState Simulation state vectors
 * @param [in] sDot Vector @f$ x @f$ that is transformed by the Jacobian @f$ \frac{\partial F}{\partial \dot{y}} @f$
 * @param [out] ret Vector @f$ z @f$ which stores the result of the operation
 */
void ColumnModel1D::multiplyWithDerivativeJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* sDot, double* ret)
{
	Indexer idxr(_disc);

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(std::size_t(0), static_cast<std::size_t>(_disc.nPoints * _disc.nParType + 1), [&](std::size_t idx)
#else
	for (unsigned int idx = 0; idx < _disc.nPoints * _disc.nParType + 1; ++idx)
#endif
	{
		if (cadet_unlikely(idx == 0))
		{
			_convDispOp.multiplyWithDerivativeJacobian(simTime, sDot, ret);
		}
		else
		{
			const unsigned int idxParLoop = idx - 1;
			const unsigned int pblk = idxParLoop % _disc.nPoints;
			const unsigned int type = idxParLoop / _disc.nPoints;

			const double invBetaP = (1.0 / static_cast<double>(_particles[type]->getPorosity()) - 1.0);
			unsigned int const* const nBound = _disc.nBound + type * _disc.nComp;
			unsigned int const* const boundOffset = _disc.boundOffset + type * _disc.nComp;
			int const* const qsReaction = _binding[type]->reactionQuasiStationarity();

			// Particle shells
			const int offsetCpType = idxr.offsetCp(ParticleTypeIndex{ type }, ParticleIndex{ pblk });
			for (unsigned int shell = 0; shell < _disc.nParPoints[type]; ++shell)
			{
				const int offsetCpShell = offsetCpType + shell * idxr.strideParNode(type);
				double const* const mobileSdot = sDot + offsetCpShell;
				double* const mobileRet = ret + offsetCpShell;

				parts::cell::multiplyWithDerivativeJacobianKernel<true>(mobileSdot, mobileRet, _disc.nComp, nBound, boundOffset, _disc.strideBound[type], qsReaction, 1.0, invBetaP);
			}
		}
	} CADET_PARFOR_END;

	// Handle inlet DOFs (all algebraic)
	std::fill_n(ret, _disc.nComp, 0.0);
}

void ColumnModel1D::setExternalFunctions(IExternalFunction** extFuns, unsigned int size)
{
	for (IBindingModel* bm : _binding)
	{
		if (bm)
			bm->setExternalFunctions(extFuns, size);
	}
}

unsigned int ColumnModel1D::localOutletComponentIndex(unsigned int port) const CADET_NOEXCEPT
{
	// Inlets are duplicated so need to be accounted for
	if (_convDispOp.forwardFlow())
		// Forward Flow: outlet is last cell
		return _disc.nComp + (_disc.nPoints - 1) * _disc.nComp;
	else
		// Backward flow: Outlet is first cell
		return _disc.nComp;
}

unsigned int ColumnModel1D::localInletComponentIndex(unsigned int port) const CADET_NOEXCEPT
{
	// Always 0 due to dedicated inlet DOFs
	return 0;
}

unsigned int ColumnModel1D::localOutletComponentStride(unsigned int port) const CADET_NOEXCEPT
{
	return 1;
}

unsigned int ColumnModel1D::localInletComponentStride(unsigned int port) const CADET_NOEXCEPT
{
	return 1;
}

void ColumnModel1D::expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut)
{
	// @todo Write this function
}

bool ColumnModel1D::setParameter(const ParameterId& pId, double value)
{
	if (pId.unitOperation == _unitOpIdx)
	{
		const int mpIc = multiplexInitialConditions(pId, value, false);
		if (mpIc > 0)
			return true;
		else if (mpIc < 0)
			return false;

		// Intercept changes to PAR_TYPE_VOLFRAC when not specified per axial cell (but once globally)
		if (_axiallyConstantParTypeVolFrac && (pId.name == hashString("PAR_TYPE_VOLFRAC")))
		{
			if ((pId.section != SectionIndep) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep))
				return false;
			if (pId.particleType >= _disc.nParType)
				return false;

			for (unsigned int i = 0; i < _disc.nPoints; ++i)
				_parTypeVolFrac[i * _disc.nParType + pId.particleType].setValue(value);

			return true;
		}

		// Parameters with particle type independence will be set for all particle types that have this parameter.
		// E.g. for a parameter type independent surface diffusion, the value is set for all particle types that have this parameter.
		bool paramExists = false;
		for (int parType = 0; parType < _disc.nParType; parType++)
		{
			const bool paramExistsNow = _particles[parType]->setParameter(pId, value);
			paramExists = paramExists || paramExistsNow;

			if (paramExists)
			{
				// Check whether particle radius or core radius has changed and update radial discretization if necessary
				if ((pId.name == hashString("PAR_RADIUS")) || (pId.name == hashString("PAR_CORERADIUS")))
					_particles[parType]->updateRadialDisc();

				if ((pId.particleType != ParTypeIndep && parType == pId.particleType) || (pId.particleType == ParTypeIndep && parType == _disc.nParType - 1))
				{
					return true;
				}
			}
		}

		if (_convDispOp.setParameter(pId, value))
			return true;

		if (model::setParameter(pId, value, _reaction.getDynReactionVector("liquid"), false))
			return true;
	}

	return UnitOperationBase::setParameter(pId, value);
}

bool ColumnModel1D::setParameter(const ParameterId& pId, int value)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return false;

	// Parameters with particle type independence will be set for all particle types that have this parameter.
	// E.g. for a parameter type independent surface diffusion, the value is set for all particle types that have this parameter.
	bool paramExists = false;
	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		const bool paramExistsNow = _particles[parType]->setParameter(pId, value);
		paramExists = paramExists || paramExistsNow;

		if (paramExists)
		{
			if ((pId.particleType != ParTypeIndep && parType == pId.particleType) || (pId.particleType == ParTypeIndep && parType == _disc.nParType - 1))
			{
				return true;
			}
		}
	}

	if (pId.unitOperation == _unitOpIdx)
	{
		if (model::setParameter(pId, value, _reaction.getDynReactionVector("liquid"), false))
			return true;
	}

	return UnitOperationBase::setParameter(pId, value);
}

bool ColumnModel1D::setParameter(const ParameterId& pId, bool value)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return false;

	// Parameters with particle type independence will be set for all particle types that have this parameter.
	// E.g. for a parameter type independent surface diffusion, the value is set for all particle types that have this parameter.
	bool paramExists = false;
	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		const bool paramExistsNow = _particles[parType]->setParameter(pId, value);
		paramExists = paramExists || paramExistsNow;

		if (paramExists)
		{
			if ((pId.particleType != ParTypeIndep && parType == pId.particleType) || (pId.particleType == ParTypeIndep && parType == _disc.nParType - 1))
			{
				return true;
			}
		}
	}

	if (pId.unitOperation == _unitOpIdx)
	{
		if (model::setParameter(pId, value, _reaction.getDynReactionVector("liquid"), false))
			return true;
	}

	return UnitOperationBase::setParameter(pId, value);
}

void ColumnModel1D::setSensitiveParameterValue(const ParameterId& pId, double value)
{
	if (pId.unitOperation == _unitOpIdx)
	{
		if (multiplexInitialConditions(pId, value, true) != 0)
			return;

		// Intercept changes to PAR_TYPE_VOLFRAC when not specified per axial cell (but once globally)
		if (_axiallyConstantParTypeVolFrac && (pId.name == hashString("PAR_TYPE_VOLFRAC")))
		{
			if ((pId.section != SectionIndep) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep))
				return;
			if (pId.particleType >= _disc.nParType)
				return;

			if (!contains(_sensParams, &_parTypeVolFrac[pId.particleType]))
				return;

			for (unsigned int i = 0; i < _disc.nPoints; ++i)
				_parTypeVolFrac[i * _disc.nParType + pId.particleType].setValue(value);

			return;
		}

		// Parameters with particle type independence will be set for all particle types that have this parameter.
		// E.g. for a parameter type independent surface diffusion, the AD value is set for all particle types that have this parameter.
		bool paramExists = false;
		for (int parType = 0; parType < _disc.nParType; parType++)
		{
			const bool paramExistsNow = _particles[parType]->setSensitiveParameterValue(_sensParams, pId, value);
			paramExists = paramExists || paramExistsNow;

			if (paramExists)
			{
				// Check whether particle radius or core radius has changed and update radial discretization if necessary
				if ((pId.name == hashString("PAR_RADIUS")) || (pId.name == hashString("PAR_CORERADIUS")))
					_particles[parType]->updateRadialDisc();

				if ((pId.particleType != ParTypeIndep && parType == pId.particleType) || (pId.particleType == ParTypeIndep && parType == _disc.nParType - 1))
				{
					return;
				}
			}
		}

		if (_convDispOp.setSensitiveParameterValue(_sensParams, pId, value))
			return;

		if (model::setSensitiveParameterValue(pId, value, _sensParams, _reaction.getDynReactionVector("liquid"), false))
			return;
	}

	return UnitOperationBase::setSensitiveParameterValue(pId, value);
}

bool ColumnModel1D::setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue)
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

		// Intercept changes to PAR_TYPE_VOLFRAC when not specified per axial cell (but once globally)
		if (_axiallyConstantParTypeVolFrac && (pId.name == hashString("PAR_TYPE_VOLFRAC")))
		{
			if ((pId.section != SectionIndep) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep))
				return false;
			if (pId.particleType >= _disc.nParType)
				return false;

			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;

			// Register parameter and set AD seed / direction
			_sensParams.insert(&_parTypeVolFrac[pId.particleType]);
			for (unsigned int i = 0; i < _disc.nPoints; ++i)
				_parTypeVolFrac[i * _disc.nParType + pId.particleType].setADValue(adDirection, adValue);

			return true;
		}

		// Parameter sensitivities with particle type independence will be set for all particle types that have this parameter.
		// E.g. for a parameter type independent surface diffusion sensitivity, the sensitivity is set for all particle types that have this parameter.
		bool paramExists = false;
		for (int parType = 0; parType < _disc.nParType; parType++)
		{
			const bool paramExistsNow = _particles[parType]->setSensitiveParameter(_sensParams, pId, adDirection, adValue);
			paramExists = paramExists || paramExistsNow;

			if (paramExists)
			{
				// Check whether particle radius or core radius has been set active and update radial discretization if necessary
				// Note that we need to recompute the radial discretization variables (_parCellSize, _parCenterRadius, _parOuterSurfAreaPerVolume, _parInnerSurfAreaPerVolume)
				// because their gradient has changed (although their nominal value has not changed).
				if ((pId.name == hashString("PAR_RADIUS")) || (pId.name == hashString("PAR_CORERADIUS")))
				{
					_particles[parType]->updateRadialDisc();
				}

				// continue loop for particle type independent parameters to set the respective parameter sensitivity in all particle types
				if ((pId.particleType != ParTypeIndep && parType == pId.particleType) || (pId.particleType == ParTypeIndep && parType == _disc.nParType - 1))
				{
					LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
					return true;
				}
			}
		}

		if (_convDispOp.setSensitiveParameter(_sensParams, pId, adDirection, adValue))
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		if (model::setSensitiveParameter(pId, adDirection, adValue, _sensParams, _reaction.getDynReactionVector("liquid"), false))
		{
			LOG(Debug) << "Found parameter " << pId << " in DynamicBulkReactionModel: Dir " << adDirection << " is set to " << adValue;
			return true;
		}
	}

	return UnitOperationBase::setSensitiveParameter(pId, adDirection, adValue);
}

std::unordered_map<ParameterId, double> ColumnModel1D::getAllParameterValues() const
{
	std::unordered_map<ParameterId, double> data = UnitOperationBase::getAllParameterValues();

	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		_particles[parType]->getAllParameterValues(data);
	}

	return data;
}

double ColumnModel1D::getParameterDouble(const ParameterId& pId) const
{
	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		double val = _particles[parType]->getParameterDouble(pId);
		if (val)
			return val;
	}

	return UnitOperationBase::getParameterDouble(pId);
}

bool ColumnModel1D::hasParameter(const ParameterId& pId) const
{
	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		if (_particles[parType]->hasParameter(pId))
			return true;
	}

	return UnitOperationBase::hasParameter(pId);
}

int ColumnModel1D::Exporter::writeMobilePhase(double* buffer) const
{
	const int blockSize = numMobilePhaseDofs();
	std::copy_n(_idx.c(_data), blockSize, buffer);
	return blockSize;
}

int ColumnModel1D::Exporter::writeSolidPhase(double* buffer) const
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

int ColumnModel1D::Exporter::writeParticleMobilePhase(double* buffer) const
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

int ColumnModel1D::Exporter::writeSolidPhase(unsigned int parType, double* buffer) const
{
	cadet_assert(parType < _disc.nParType);

	const unsigned int stride = _disc.nComp + _disc.strideBound[parType];
	double const* ptr = _data + _idx.offsetCp(ParticleTypeIndex{ parType }) + _disc.nComp;
	for (unsigned int i = 0; i < _disc.nPoints; ++i)
	{
		for (unsigned int j = 0; j < _disc.nParPoints[parType]; ++j)
		{
			std::copy_n(ptr, _disc.strideBound[parType], buffer);
			buffer += _disc.strideBound[parType];
			ptr += stride;
		}
	}
	return _disc.nPoints * _disc.nParPoints[parType] * _disc.strideBound[parType];
}

int ColumnModel1D::Exporter::writeParticleMobilePhase(unsigned int parType, double* buffer) const
{
	cadet_assert(parType < _disc.nParType);

	const unsigned int stride = _disc.nComp + _disc.strideBound[parType];
	double const* ptr = _data + _idx.offsetCp(ParticleTypeIndex{ parType });
	for (unsigned int i = 0; i < _disc.nPoints; ++i)
	{
		for (unsigned int j = 0; j < _disc.nParPoints[parType]; ++j)
		{
			std::copy_n(ptr, _disc.nComp, buffer);
			buffer += _disc.nComp;
			ptr += stride;
		}
	}
	return _disc.nPoints * _disc.nParPoints[parType] * _disc.nComp;
}

int ColumnModel1D::Exporter::writeParticleFlux(double* buffer) const
{
	return 0;
}

int ColumnModel1D::Exporter::writeParticleFlux(unsigned int parType, double* buffer) const
{
	return 0;
}

int ColumnModel1D::Exporter::writeInlet(unsigned int port, double* buffer) const
{
	cadet_assert(port == 0);
	std::copy_n(_data, _disc.nComp, buffer);
	return _disc.nComp;
}

int ColumnModel1D::Exporter::writeInlet(double* buffer) const
{
	std::copy_n(_data, _disc.nComp, buffer);
	return _disc.nComp;
}

int ColumnModel1D::Exporter::writeOutlet(unsigned int port, double* buffer) const
{
	cadet_assert(port == 0);

	if (_model._convDispOp.forwardFlow())
		std::copy_n(&_idx.c(_data, _disc.nPoints - 1, 0), _disc.nComp, buffer);
	else
		std::copy_n(&_idx.c(_data, 0, 0), _disc.nComp, buffer);

	return _disc.nComp;
}

int ColumnModel1D::Exporter::writeOutlet(double* buffer) const
{
	if (_model._convDispOp.forwardFlow())
		std::copy_n(&_idx.c(_data, _disc.nPoints - 1, 0), _disc.nComp, buffer);
	else
		std::copy_n(&_idx.c(_data, 0, 0), _disc.nComp, buffer);

	return _disc.nComp;
}

}  // namespace model

}  // namespace cadet

#include "model/ColumnModel1D-InitialConditions.cpp"
#include "model/ColumnModel1D-LinearSolver.cpp"
