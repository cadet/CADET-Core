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

#include "model/GeneralRateModelDG.hpp"
#include "BindingModelFactory.hpp"
#include "ReactionModelFactory.hpp"
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


GeneralRateModelDG::GeneralRateModelDG(UnitOpIdx unitOpIdx) : UnitOperationBase(unitOpIdx),
	_globalJac(), _globalJacDisc(), _jacInlet(), _dynReactionBulk(nullptr),
	_analyticJac(true), _jacobianAdDirs(0), _factorizeJacobian(false), _tempState(nullptr),
	_initC(0), _initCp(0), _initQ(0), _initState(0), _initStateDot(0)
{
}

GeneralRateModelDG::~GeneralRateModelDG() CADET_NOEXCEPT
{
	delete[] _tempState;

	_binding.clear(); // binding models are deleted in the respective particle model
	_dynReaction.clear(); // particle reaction models are deleted in the respective particle model
	delete[] _particle;

	delete _dynReactionBulk;

	delete _linearSolver;
}

unsigned int GeneralRateModelDG::numDofs() const CADET_NOEXCEPT
{
	// Column bulk DOFs: nPoints * nComp
	// Particle DOFs: nPoints * nParType particles each having nComp (liquid phase) + sum boundStates (solid phase) DOFs
	//                in each shell; there are nParCell shells for each particle type
	// Inlet DOFs: nComp
	return _disc.nPoints * _disc.nComp + _disc.parTypeOffset[_disc.nParType] + _disc.nComp;
}

unsigned int GeneralRateModelDG::numPureDofs() const CADET_NOEXCEPT
{
	// Column bulk DOFs: nPoints * nComp
	// Particle DOFs: nPoints particles each having nComp (liquid phase) + sum boundStates (solid phase) DOFs
	//                in each shell; there are nPar shells
	return _disc.nPoints * _disc.nComp  + _disc.parTypeOffset[_disc.nParType];
}


bool GeneralRateModelDG::usesAD() const CADET_NOEXCEPT
{
#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	// We always need AD if we want to check the analytical Jacobian
	return true;
#else
	// We only need AD if we are not computing the Jacobian analytically
	return !_analyticJac;
#endif
}

bool GeneralRateModelDG::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper)
{
	// ==== Read discretization
	_disc.nComp = paramProvider.getInt("NCOMP");

	std::vector<int> nBound;
	const bool newNBoundInterface = paramProvider.exists("NBOUND");
	const bool newNPartypeInterface = paramProvider.exists("NPARTYPE");

	paramProvider.pushScope("discretization");

	const bool firstConfigCall = _tempState == nullptr; // used to not multiply allocate memory

	if (firstConfigCall)
		_linearSolver = cadet::linalg::setLinearSolver(paramProvider.exists("LINEAR_SOLVER") ? paramProvider.getString("LINEAR_SOLVER") : "SparseLU");

	if (!newNBoundInterface && paramProvider.exists("NBOUND")) // done here and in this order for backwards compatibility
		nBound = paramProvider.getIntArray("NBOUND");
	else
	{
		paramProvider.popScope();
		nBound = paramProvider.getIntArray("NBOUND");
		paramProvider.pushScope("discretization");
	}
	if (nBound.size() < _disc.nComp)
		throw InvalidParameterException("Field NBOUND contains too few elements (NCOMP = " + std::to_string(_disc.nComp) + " required)");

	if (paramProvider.exists("POLYDEG"))
		_disc.polyDeg = paramProvider.getInt("POLYDEG");
	else
		_disc.polyDeg = 4u; // default value
	if (paramProvider.getInt("POLYDEG") < 1)
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
	
	if (!newNPartypeInterface && paramProvider.exists("NPARTYPE")) // done here and in this order for backwards compatibility
		_disc.nParType = paramProvider.getInt("NPARTYPE");
	else if (newNPartypeInterface)
	{
		paramProvider.popScope();
		_disc.nParType = paramProvider.getInt("NPARTYPE");
		paramProvider.pushScope("discretization");
	}
	else // Infer number of particle types
		_disc.nParType = nBound.size() / _disc.nComp;

	paramProvider.popScope();
	Indexer idxr(_disc);
	_particle = new parts::GeneralRateParticle[_disc.nParType];
	bool particleConfSuccess = true;
	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		particleConfSuccess = particleConfSuccess && _particle[parType].configureModelDiscretization(paramProvider, helper, _disc.nComp, parType, _disc.nParType, idxr.strideColComp());
	}
	paramProvider.pushScope("discretization");

	if (firstConfigCall)
	{
		_disc.nParPoints = new unsigned int[_disc.nParType];
		for (int type = 0; type < _disc.nParType; type++)
		{
			_disc.nParPoints[type] = _particle[type].nDiscPoints();
		}
	}

	_disc.newStaticJac = true;

	if (firstConfigCall)
		_disc.nBound = new unsigned int[_disc.nComp * _disc.nParType];
	if (nBound.size() < _disc.nComp * _disc.nParType)
	{
		// Multiplex number of bound states to all particle types
		for (unsigned int i = 0; i < _disc.nParType; ++i)
			std::copy_n(nBound.begin(), _disc.nComp, _disc.nBound + i * _disc.nComp);
	}
	else
		std::copy_n(nBound.begin(), _disc.nComp * _disc.nParType, _disc.nBound);

	const unsigned int nTotalBound = std::accumulate(_disc.nBound, _disc.nBound + _disc.nComp * _disc.nParType, 0u);

	// Precompute offsets and total number of bound states (DOFs in solid phase)
	if (firstConfigCall)
	{
		_disc.boundOffset = new unsigned int[_disc.nComp * _disc.nParType];
		_disc.strideBound = new unsigned int[_disc.nParType + 1];
		_disc.nBoundBeforeType = new unsigned int[_disc.nParType];
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
	_initQ.resize(nTotalBound);

	// Create nonlinear solver for consistent initialization
	configureNonlinearSolver(paramProvider);

	paramProvider.popScope();

	// ==== Construct and configure parameter dependencies

	unsigned int strideColNode = _disc.nComp;
	const bool transportSuccess = _convDispOp.configureModelDiscretization(paramProvider, helper, _disc.nComp, polynomial_integration_mode, _disc.nElem, _disc.polyDeg, strideColNode);

	_disc.curSection = -1;

	// ==== Construct and configure dynamic reaction model
	bool reactionConfSuccess = true;

	_dynReactionBulk = nullptr;
	if (paramProvider.exists("REACTION_MODEL"))
	{
		const std::string dynReactName = paramProvider.getString("REACTION_MODEL");
		_dynReactionBulk = helper.createDynamicReactionModel(dynReactName);
		if (!_dynReactionBulk)
			throw InvalidParameterException("Unknown dynamic reaction model " + dynReactName);

		if (_dynReactionBulk->usesParamProviderInDiscretizationConfig())
			paramProvider.pushScope("reaction_bulk");

		reactionConfSuccess = _dynReactionBulk->configureModelDiscretization(paramProvider, _disc.nComp, nullptr, nullptr);

		if (_dynReactionBulk->usesParamProviderInDiscretizationConfig())
			paramProvider.popScope();
	}

	// ==== Construct and configure binding and particle reaction -> done in particle model, only pointers are copied here.
	_binding = std::vector<IBindingModel*>(_disc.nParType, nullptr);
	_dynReaction = std::vector<IDynamicReactionModel*>(_disc.nParType, nullptr);

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		_binding[parType] = _particle[parType].getBinding();
		_singleBinding = _particle[parType].singleBinding();
		if (parType > 0 && _singleBinding != _particle[parType].singleBinding())
			throw InvalidParameterException("Configuration of binding went wrong");

		_dynReaction[parType] = _particle[parType].getReaction();
		_singleDynReaction = _particle[parType].singleReaction();
		if (parType > 0 && _singleDynReaction != _particle[parType].singleReaction())
			throw InvalidParameterException("Configuration of particle reaction went wrong");
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
 
bool GeneralRateModelDG::configure(IParameterProvider& paramProvider)
{
	_parameters.clear();

	const bool transportSuccess = _convDispOp.configure(_unitOpIdx, paramProvider, _parameters);

	// Read geometry parameters handled by unit operation
	_colPorosity = paramProvider.getDouble("COL_POROSITY");

	// Check whether PAR_TYPE_VOLFRAC is required or not
	if ((_disc.nParType > 1) && !paramProvider.exists("PAR_TYPE_VOLFRAC"))
		throw InvalidParameterException("The required parameter \"PAR_TYPE_VOLFRAC\" was not found");

	// Let PAR_TYPE_VOLFRAC default to 1.0 for backwards compatibility
	if (paramProvider.exists("PAR_TYPE_VOLFRAC"))
	{
		readScalarParameterOrArray(_parTypeVolFrac, paramProvider, "PAR_TYPE_VOLFRAC", 1);
		if (_parTypeVolFrac.size() == _disc.nParType)
		{
			_axiallyConstantParTypeVolFrac = true;

			// Expand to all axial cells
			_parTypeVolFrac.resize(_disc.nPoints * _disc.nParType, 1.0);
			for (unsigned int i = 1; i < _disc.nPoints; ++i)
				std::copy(_parTypeVolFrac.begin(), _parTypeVolFrac.begin() + _disc.nParType, _parTypeVolFrac.begin() + _disc.nParType * i);
		}
		else
			_axiallyConstantParTypeVolFrac = false;
	}
	else
	{
		_parTypeVolFrac.resize(_disc.nPoints, 1.0);
		_axiallyConstantParTypeVolFrac = false;
	}

	// Check whether all sizes are matched
	if (_disc.nParType * _disc.nPoints != _parTypeVolFrac.size())
		throw InvalidParameterException("Number of elements in field PAR_TYPE_VOLFRAC does not match number of particle types times number of axial cells");

	// Check that particle volume fractions sum to 1.0
	for (unsigned int i = 0; i < _disc.nPoints; ++i)
	{
		const double volFracSum = std::accumulate(_parTypeVolFrac.begin() + i * _disc.nParType, _parTypeVolFrac.begin() + (i+1) * _disc.nParType, 0.0,
			[](double a, const active& b) -> double { return a + static_cast<double>(b); });
		if (std::abs(1.0 - volFracSum) > 1e-10)
			throw InvalidParameterException("Sum of field PAR_TYPE_VOLFRAC differs from 1.0 (is " + std::to_string(volFracSum) + ") in axial cell " + std::to_string(i));
	}

	// Read vectorial parameters (which may also be section dependent; transport)
	_filmDiffusionMode = readAndRegisterMultiplexCompTypeSecParam(paramProvider, _parameters, _filmDiffusion, "FILM_DIFFUSION", _disc.nParType, _disc.nComp, _unitOpIdx);

	if ((_filmDiffusion.size() < _disc.nComp * _disc.nParType) || (_filmDiffusion.size() % (_disc.nComp * _disc.nParType) != 0))
		throw InvalidParameterException("Number of elements in field FILM_DIFFUSION is not a positive multiple of NCOMP * NPARTYPE (" + std::to_string(_disc.nComp * _disc.nParType) + ")");

	if (paramProvider.exists("PORE_ACCESSIBILITY"))
		_poreAccessFactorMode = readAndRegisterMultiplexCompTypeSecParam(paramProvider, _parameters, _poreAccessFactor, "PORE_ACCESSIBILITY", _disc.nParType, _disc.nComp, _unitOpIdx);
	else
	{
		_poreAccessFactorMode = MultiplexMode::ComponentType;
		_poreAccessFactor = std::vector<cadet::active>(_disc.nComp * _disc.nParType, 1.0);
	}

	if (_disc.nComp * _disc.nParType != _poreAccessFactor.size())
		throw InvalidParameterException("Number of elements in field PORE_ACCESSIBILITY differs from NCOMP * NPARTYPE (" + std::to_string(_disc.nComp * _disc.nParType) + ")");

	// Add parameters to map
	_parameters[makeParamId(hashString("COL_POROSITY"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_colPorosity;

	if (_axiallyConstantParTypeVolFrac)
	{
		// Register only the first nParType items
		for (unsigned int i = 0; i < _disc.nParType; ++i)
			_parameters[makeParamId(hashString("PAR_TYPE_VOLFRAC"), _unitOpIdx, CompIndep, i, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parTypeVolFrac[i];
	}
	else
		registerParam2DArray(_parameters, _parTypeVolFrac, [=](bool multi, unsigned cell, unsigned int type) { return makeParamId(hashString("PAR_TYPE_VOLFRAC"), _unitOpIdx, CompIndep, type, BoundStateIndep, ReactionIndep, cell); }, _disc.nParType);

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

			active* const iq = _initQ.data() + _disc.nBoundBeforeType[0];
			for (unsigned int i = 0; i < _disc.strideBound[0]; ++i)
				_parameters[initParams[i]] = iq + i;
		}
		else
		{
			for (unsigned int type = 0; type < _disc.nParType; ++type)
			{
				_binding[type]->fillBoundPhaseInitialParameters(initParams.data(), _unitOpIdx, type);

				active* const iq = _initQ.data() + _disc.nBoundBeforeType[type];
				for (unsigned int i = 0; i < _disc.strideBound[type]; ++i)
					_parameters[initParams[i]] = iq + i;
			}
		}
	}

	// Reconfigure particle model
	bool particleConfSuccess = true;
	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		particleConfSuccess = particleConfSuccess && _particle[parType].configure(_unitOpIdx, paramProvider, _parameters, _disc.nParType, _disc.nBoundBeforeType, _disc.strideBound[_disc.nParType]);
	}

	// Reconfigure reaction model
	bool dynReactionConfSuccess = true;
	if (_dynReactionBulk && _dynReactionBulk->requiresConfiguration())
	{
		paramProvider.pushScope("reaction_bulk");
		dynReactionConfSuccess = _dynReactionBulk->configure(paramProvider, _unitOpIdx, ParTypeIndep);
		paramProvider.popScope();
	}

	// jaobian pattern set after binding and particle surface diffusion are configured
	setJacobianPattern_GRM(_globalJac, 0, _dynReactionBulk);
	_globalJacDisc = _globalJac;
	// the solver repetitively solves the linear system with a static pattern of the jacobian (set above). 
	// The goal of analyzePattern() is to reorder the nonzero elements of the matrix, such that the factorization step creates less fill-in
	_linearSolver->analyzePattern(_globalJacDisc.block(_disc.nComp, _disc.nComp, numPureDofs(), numPureDofs()));

	return transportSuccess && particleConfSuccess && dynReactionConfSuccess;
}

unsigned int GeneralRateModelDG::threadLocalMemorySize() const CADET_NOEXCEPT
{
	LinearMemorySizer lms;

	// Memory for residualImpl()
	for (unsigned int i = 0; i < _disc.nParType; ++i)
	{
		if (_binding[i] && _binding[i]->requiresWorkspace())
			lms.fitBlock(_binding[i]->workspaceSize(_disc.nComp, _disc.strideBound[i], _disc.nBound + i * _disc.nComp));

		if (_dynReaction[i] && _dynReaction[i]->requiresWorkspace())
			lms.fitBlock(_dynReaction[i]->workspaceSize(_disc.nComp, _disc.strideBound[i], _disc.nBound + i * _disc.nComp));
	}

	if (_dynReactionBulk && _dynReactionBulk->requiresWorkspace())
		lms.fitBlock(_dynReactionBulk->workspaceSize(_disc.nComp, 0, nullptr));

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

unsigned int GeneralRateModelDG::numAdDirsForJacobian() const CADET_NOEXCEPT
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

void GeneralRateModelDG::useAnalyticJacobian(const bool analyticJac)
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

void GeneralRateModelDG::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, const ConstSimulationState& simState, const AdJacobianParams& adJac)
{
	Indexer idxr(_disc);

	// todo: only reset jacobian pattern if it changes, i.e. once in configuration and then only for changes in SurfDiff+kinetic binding.
 	setJacobianPattern_GRM(_globalJac, 0, _dynReactionBulk);
	_globalJacDisc = _globalJac;

	_convDispOp.notifyDiscontinuousSectionTransition(t, secIdx, _jacInlet);

	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		_particle[parType].notifyDiscontinuousSectionTransition(t, secIdx, getSectionDependentSlice(_filmDiffusion, _disc.nComp * _disc.nParType, secIdx), &_poreAccessFactor[0]);
	}

	_disc.curSection = secIdx;
	_disc.newStaticJac = true;
}

void GeneralRateModelDG::setFlowRates(active const* in, active const* out) CADET_NOEXCEPT
{
	_convDispOp.setFlowRates(in[0], out[0], _colPorosity);
}

void GeneralRateModelDG::reportSolution(ISolutionRecorder& recorder, double const* const solution) const
{
	Exporter expr(_disc, *this, solution);
	recorder.beginUnitOperation(_unitOpIdx, *this, expr);
	recorder.endUnitOperation();
}

void GeneralRateModelDG::reportSolutionStructure(ISolutionRecorder& recorder) const
{
	Exporter expr(_disc, *this, nullptr);
	recorder.unitOperationStructure(_unitOpIdx, *this, expr);
}

unsigned int GeneralRateModelDG::requiredADdirs() const CADET_NOEXCEPT
{
	const unsigned int numDirsBinding = maxBindingAdDirs();
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	return numDirsBinding + _jacobianAdDirs;
#else
	// If CADET_CHECK_ANALYTIC_JACOBIAN is active, we always need the AD directions for the Jacobian
	return numDirsBinding + numAdDirsForJacobian();
#endif
}

void GeneralRateModelDG::prepareADvectors(const AdJacobianParams& adJac) const
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
void GeneralRateModelDG::extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset)
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
		_particle[parType].calcFilmDiffJacobian(_disc.curSection, idxr.offsetCp(ParticleTypeIndex{static_cast<unsigned int>(parType)}), idxr.offsetC(), _disc.nPoints, _disc.nParType, static_cast<double>(_colPorosity), &_parTypeVolFrac[0], _globalJac, true);
	}
}

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN

/**
 * @brief Compares the analytical Jacobian with a Jacobian derived by AD
 * @details The analytical Jacobian is assumed to be stored in the corresponding band matrices.
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
void GeneralRateModelDG::checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const
{
	// todo write this function
	//Indexer idxr(_disc);

	//LOG(Debug) << "AD dir offset: " << adDirOffset << " DiagDirCol: " << _convDispOp.jacobian().lowerBandwidth() << " DiagDirPar: " << _jacP[0].lowerBandwidth();

	//// Column
	//const double maxDiffCol = _convDispOp.checkAnalyticJacobianAgainstAd(adRes, adDirOffset);

	//// Particles
	//double maxDiffPar = 0.0;
	//for (unsigned int type = 0; type < _disc.nParType; ++type)
	//{
	//	for (unsigned int pblk = 0; pblk < _disc.nPoints; ++pblk)
	//	{
	//		linalg::BandMatrix& jacMat = _jacP[_disc.nPoints * type + pblk];
	//		const double localDiff = ad::compareBandedJacobianWithAd(adRes + idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{pblk}), adDirOffset, jacMat.lowerBandwidth(), jacMat);
	//		LOG(Debug) << "-> Par type " << type << " block " << pblk << " diff: " << localDiff;
	//		maxDiffPar = std::max(maxDiffPar, localDiff);
	//	}
	//}
}

#endif

int GeneralRateModelDG::jacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	_factorizeJacobian = true;

	if (_analyticJac)
		return residualImpl<double, double, double, true, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
	else
		return residualWithJacobian(simTime, ConstSimulationState{ simState.vecStateY, nullptr }, nullptr, adJac, threadLocalMem);
}

int GeneralRateModelDG::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	// Evaluate residual do not compute Jacobian or parameter sensitivities
	return residualImpl<double, double, double, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
}

int GeneralRateModelDG::residualWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	//FDJac = calcFDJacobian(static_cast<const double*>(simState.vecStateY), static_cast<const double*>(simState.vecStateYdot), simTime, threadLocalMem, 2.0); // debug code

	// Evaluate residual, use AD for Jacobian if required but do not evaluate parameter derivatives
	return residual(simTime, simState, res, adJac, threadLocalMem, true, false);
}

int GeneralRateModelDG::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res,
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
int GeneralRateModelDG::residualImpl(double t, unsigned int secIdx, StateType const* const y, double const* const yDot, ResidualType* const res, util::ThreadLocalStorage& threadLocalMem)
{
	if (wantRes)
	{
		double* const resPtr = reinterpret_cast<double* const>(res);
		Eigen::Map<Eigen::VectorXd> resi(resPtr, numDofs());
		resi.setZero();
	}

	if (wantJac)
	{
		if (!wantRes || _disc.newStaticJac)
		{
			// estimate new static (per section) jacobian
			bool success = calcStaticAnaJacobian_GRM(secIdx);

			_disc.newStaticJac = false;

			if (cadet_unlikely(!success)) {
				LOG(Error) << "Jacobian pattern did not fit the Jacobian estimation";
			}
		}
	}

	residualBulk<StateType, ResidualType, ParamType, wantJac, wantRes>(t, secIdx, y, yDot, res, threadLocalMem);

	BENCH_START(_timerResidualPar);

	LinearBufferAllocator tlmAlloc = threadLocalMem.get();
	Indexer idxr(_disc);

	for (unsigned int pblk = 0; pblk < _disc.nPoints * _disc.nParType; ++pblk)
	{
		const unsigned int parType = pblk / _disc.nPoints;
		const unsigned int colNode = pblk % _disc.nPoints;

		linalg::BandedEigenSparseRowIterator jacIt(_globalJac, idxr.offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode }));
		ColumnPosition colPos{ _convDispOp.relativeCoordinate(colNode), 0.0, 0.0 }; // Relative position of current node - needed in externally dependent adsorption kinetic

		_particle[parType].residual<wantJac, wantRes>(t, secIdx,
			y + idxr.offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode }),
			y + idxr.offsetC() + colNode * idxr.strideColNode(),
			yDot ? yDot + idxr.offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode }) : nullptr,
			res ? res + idxr.offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colNode }) : nullptr,
			colPos, jacIt, tlmAlloc,
			typename cadet::ParamSens<ParamType>::enabled()
		);
	}

	if (!wantRes)
		return 0;

	BENCH_STOP(_timerResidualPar);

	residualFlux<StateType, ResidualType, ParamType>(t, secIdx, y, yDot, res);

	// Handle inlet DOFs, which are simply copied to the residual
	for (unsigned int i = 0; i < _disc.nComp; ++i)
	{
		res[i] = y[i];
	}

	return 0;
}

template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
int GeneralRateModelDG::residualBulk(double t, unsigned int secIdx, StateType const* yBase, double const* yDotBase, ResidualType* resBase, util::ThreadLocalStorage& threadLocalMem)
{
	if (wantRes)
		_convDispOp.residual(*this, t, secIdx, yBase, yDotBase, resBase, typename cadet::ParamSens<ParamType>::enabled());

	if (!_dynReactionBulk || (_dynReactionBulk->numReactionsLiquid() == 0))
		return 0;

	Indexer idxr(_disc);
	LinearBufferAllocator tlmAlloc = threadLocalMem.get();

	// Dynamic reactions
	if (_dynReactionBulk) {

		StateType const* y = yBase + idxr.offsetC();

		if (wantJac && !wantRes) // only compute Jacobian
		{
			for (unsigned int col = 0; col < _disc.nPoints; ++col, y += idxr.strideColNode())
			{
				const ColumnPosition colPos{ (0.5 + static_cast<double>(col)) / static_cast<double>(_disc.nPoints), 0.0, 0.0 };

				linalg::BandedEigenSparseRowIterator jac(_globalJac, col * idxr.strideColNode());
				// static_cast should be sufficient here, but this statement is also analyzed when wantJac = false
				_dynReactionBulk->analyticJacobianLiquidAdd(t, secIdx, colPos, reinterpret_cast<double const*>(y), -1.0, jac, tlmAlloc);
			}

			return 0;
		}

		ResidualType* res = resBase + idxr.offsetC();

		for (unsigned int col = 0; col < _disc.nPoints; ++col, y += idxr.strideColNode(), res += idxr.strideColNode())
		{
			const ColumnPosition colPos{ _convDispOp.relativeCoordinate(col), 0.0, 0.0};
			_dynReactionBulk->residualLiquidAdd(t, secIdx, colPos, y, res, -1.0, tlmAlloc);

			if (wantJac)
			{
				linalg::BandedEigenSparseRowIterator jac(_globalJac, col * idxr.strideColNode());
				// static_cast should be sufficient here, but this statement is also analyzed when wantJac = false
				_dynReactionBulk->analyticJacobianLiquidAdd(t, secIdx, colPos, reinterpret_cast<double const*>(y), -1.0, jac, tlmAlloc);
			}
		}
	}

	return 0;
}

template <typename StateType, typename ResidualType, typename ParamType>
int GeneralRateModelDG::residualFlux(double t, unsigned int secIdx, StateType const* yBase, double const* yDotBase, ResidualType* resBase)
{
	Indexer idxr(_disc);

	const ParamType invBetaC = 1.0 / static_cast<ParamType>(_colPorosity) - 1.0;

	// Get offsets
	ResidualType* const resCol = resBase + idxr.offsetC();
	StateType const* const yCol = yBase + idxr.offsetC();

	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		ResidualType* const resParType = resBase + idxr.offsetCp(ParticleTypeIndex{type});
		StateType const* const yParType = yBase + idxr.offsetCp(ParticleTypeIndex{type});

		const ParamType epsP = static_cast<ParamType>(_particle[type].getPorosity());

		// Ordering of diffusion:
		// sec0type0comp0, sec0type0comp1, sec0type0comp2, sec0type1comp0, sec0type1comp1, sec0type1comp2,
		// sec1type0comp0, sec1type0comp1, sec1type0comp2, sec1type1comp0, sec1type1comp1, sec1type1comp2, ...
		active const* const filmDiff = getSectionDependentSlice(_filmDiffusion, _disc.nComp * _disc.nParType, secIdx) + type * _disc.nComp;

		const ParamType surfaceToVolumeRatio = _particle[type].surfaceToVolumeRatio<ParamType>();

		const ParamType jacCF_val = invBetaC * surfaceToVolumeRatio;
		const ParamType jacPF_val = -1.0 / epsP;

		// Add flux to column void / bulk volume
		for (unsigned int i = 0; i < _disc.nPoints * _disc.nComp; ++i)
		{
			const unsigned int colNode = i / _disc.nComp;
			const unsigned int comp = i - colNode * _disc.nComp;
			// + 1/Beta_c * (surfaceToVolumeRatio_{p,j}) * d_j * (k_f * [c_l - c_p])
			resCol[i] += static_cast<ParamType>(filmDiff[comp]) * jacCF_val * static_cast<ParamType>(_parTypeVolFrac[type + colNode * _disc.nParType])
				        * (yCol[i] - yParType[colNode * idxr.strideParBlock(type) + (_disc.nParPoints[type] - 1) * idxr.strideParNode(type) + comp]);
		}

		//  Bead boundary condition is computed in particle residual.

	}

	return 0;
}

parts::cell::CellParameters GeneralRateModelDG::makeCellResidualParams(unsigned int parType, int const* qsReaction) const
{
	return parts::cell::CellParameters
		{
			_disc.nComp,
			_disc.nBound + _disc.nComp * parType,
			_disc.boundOffset + _disc.nComp * parType,
			_disc.strideBound[parType],
			qsReaction,
			_particle[parType].getPorosity(),
			_poreAccessFactor.data() + _disc.nComp * parType,
			_binding[parType],
			(_dynReaction[parType] && (_dynReaction[parType]->numReactionsCombined() > 0)) ? _dynReaction[parType] : nullptr
		};
}

int GeneralRateModelDG::residualSensFwdWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidualSens);

	// Evaluate residual for all parameters using AD in vector mode and at the same time update the
	// Jacobian (in one AD run, if analytic Jacobians are disabled)
	return residual(simTime, simState, nullptr, adJac, threadLocalMem, true, true);
}

int GeneralRateModelDG::residualSensFwdAdOnly(const SimulationTime& simTime, const ConstSimulationState& simState, active* const adRes, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidualSens);

	// Evaluate residual for all parameters using AD in vector mode
	return residualImpl<double, active, active, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, adRes, threadLocalMem);
}

int GeneralRateModelDG::residualSensFwdCombine(const SimulationTime& simTime, const ConstSimulationState& simState,
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
void GeneralRateModelDG::multiplyWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* yS, double alpha, double beta, double* ret)
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
void GeneralRateModelDG::multiplyWithDerivativeJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* sDot, double* ret)
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

			const double invBetaP = (1.0 / static_cast<double>(_particle[type].getPorosity()) - 1.0);
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

void GeneralRateModelDG::setExternalFunctions(IExternalFunction** extFuns, unsigned int size)
{
	for (IBindingModel* bm : _binding)
	{
		if (bm)
			bm->setExternalFunctions(extFuns, size);
	}
}

unsigned int GeneralRateModelDG::localOutletComponentIndex(unsigned int port) const CADET_NOEXCEPT
{
	// Inlets are duplicated so need to be accounted for
	if (_convDispOp.forwardFlow())
		// Forward Flow: outlet is last cell
		return _disc.nComp + (_disc.nPoints - 1) * _disc.nComp;
	else
		// Backward flow: Outlet is first cell
		return _disc.nComp;
}

unsigned int GeneralRateModelDG::localInletComponentIndex(unsigned int port) const CADET_NOEXCEPT
{
	// Always 0 due to dedicated inlet DOFs
	return 0;
}

unsigned int GeneralRateModelDG::localOutletComponentStride(unsigned int port) const CADET_NOEXCEPT
{
	return 1;
}

unsigned int GeneralRateModelDG::localInletComponentStride(unsigned int port) const CADET_NOEXCEPT
{
	return 1;
}

void GeneralRateModelDG::expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut)
{
	// @todo Write this function
}

bool GeneralRateModelDG::setParameter(const ParameterId& pId, double value)
{
	if (pId.unitOperation == _unitOpIdx)
	{
		if (multiplexCompTypeSecParameterValue(pId, hashString("PORE_ACCESSIBILITY"), _poreAccessFactorMode, _poreAccessFactor, _disc.nParType, _disc.nComp, value, nullptr))
			return true;
		if (multiplexCompTypeSecParameterValue(pId, hashString("FILM_DIFFUSION"), _filmDiffusionMode, _filmDiffusion, _disc.nParType, _disc.nComp, value, nullptr))
			return true;
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

		for (int parType = 0; parType < _disc.nParType; parType++)
		{
			if (_particle[parType].setParameter(pId, value))
				return true;
		}

		if (_convDispOp.setParameter(pId, value))
			return true;

		if (model::setParameter(pId, value, std::vector<IDynamicReactionModel*>{ _dynReactionBulk }, true))
			return true;
	}

	const bool result = UnitOperationBase::setParameter(pId, value);

	// Check whether particle radius or core radius has changed and update radial discretization if necessary
	if (result && ((pId.name == hashString("PAR_RADIUS")) || (pId.name == hashString("PAR_CORERADIUS"))))
	{
		for (int parType = 0; parType < _disc.nParType; parType++)
		{
			_particle[parType].updateRadialDisc();
		}
	}

	return result;
}

bool GeneralRateModelDG::setParameter(const ParameterId& pId, int value)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return false;

	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		if (_particle[parType].setParameter(pId, value))
			return true;
	}

	if (pId.unitOperation == _unitOpIdx)
	{
		if (model::setParameter(pId, value, std::vector<IDynamicReactionModel*>{ _dynReactionBulk }, true))
			return true;
	}

	return UnitOperationBase::setParameter(pId, value);
}

bool GeneralRateModelDG::setParameter(const ParameterId& pId, bool value)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return false;

	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		if (_particle[parType].setParameter(pId, value))
			return true;
	}

	if (pId.unitOperation == _unitOpIdx)
	{
		if (model::setParameter(pId, value, std::vector<IDynamicReactionModel*>{ _dynReactionBulk }, true))
			return true;
	}

	return UnitOperationBase::setParameter(pId, value);
}

void GeneralRateModelDG::setSensitiveParameterValue(const ParameterId& pId, double value)
{
	if (pId.unitOperation == _unitOpIdx)
	{
		if (multiplexCompTypeSecParameterValue(pId, hashString("PORE_ACCESSIBILITY"), _poreAccessFactorMode, _poreAccessFactor, _disc.nParType, _disc.nComp, value, &_sensParams))
			return;
		if (multiplexCompTypeSecParameterValue(pId, hashString("FILM_DIFFUSION"), _filmDiffusionMode, _filmDiffusion, _disc.nParType, _disc.nComp, value, &_sensParams))
			return;
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

		for (int parType = 0; parType < _disc.nParType; parType++)
		{
			if (_particle[parType].setSensitiveParameterValue(_sensParams, pId, value))
				return;
		}

		if (_convDispOp.setSensitiveParameterValue(_sensParams, pId, value))
			return;

		if (model::setSensitiveParameterValue(pId, value, _sensParams, std::vector<IDynamicReactionModel*>{ _dynReactionBulk }, true))
			return;
	}

	UnitOperationBase::setSensitiveParameterValue(pId, value);

	// Check whether particle radius or core radius has changed and update radial discretization if necessary
	if ((pId.name == hashString("PAR_RADIUS")) || (pId.name == hashString("PAR_CORERADIUS")))
	{
		for (int parType = 0; parType < _disc.nParType; parType++)
		{
			_particle[parType].updateRadialDisc();
		}
	}
}

bool GeneralRateModelDG::setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue)
{
	if (pId.unitOperation == _unitOpIdx)
	{
		if (multiplexCompTypeSecParameterAD(pId, hashString("PORE_ACCESSIBILITY"), _poreAccessFactorMode, _poreAccessFactor, _disc.nParType, _disc.nComp, adDirection, adValue, _sensParams))
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		if (multiplexCompTypeSecParameterAD(pId, hashString("FILM_DIFFUSION"), _filmDiffusionMode, _filmDiffusion, _disc.nParType, _disc.nComp, adDirection, adValue, _sensParams))
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}

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

		for (int parType = 0; parType < _disc.nParType; parType++)
		{
			if (_particle[parType].setSensitiveParameter(_sensParams, pId, adDirection, adValue))
			{
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

		if (model::setSensitiveParameter(pId, adDirection, adValue, _sensParams, std::vector<IDynamicReactionModel*> { _dynReactionBulk }, true))
		{
			LOG(Debug) << "Found parameter " << pId << " in DynamicBulkReactionModel: Dir " << adDirection << " is set to " << adValue;
			return true;
		}
	}

	const bool result = UnitOperationBase::setSensitiveParameter(pId, adDirection, adValue);

	// Check whether particle radius or core radius has been set active and update radial discretization if necessary
	// Note that we need to recompute the radial discretization variables (_parCellSize, _parCenterRadius, _parOuterSurfAreaPerVolume, _parInnerSurfAreaPerVolume)
	// because their gradient has changed (although their nominal value has not changed).
	if ((pId.name == hashString("PAR_RADIUS")) || (pId.name == hashString("PAR_CORERADIUS")))
	{
		for (int parType = 0; parType < _disc.nParType; parType++)
		{
			_particle[parType].updateRadialDisc();
		}
	}

	return result;
}

std::unordered_map<ParameterId, double> GeneralRateModelDG::getAllParameterValues() const
{
	std::unordered_map<ParameterId, double> data = UnitOperationBase::getAllParameterValues();

	std::vector<IParameterStateDependence*> parDepSurfDiffusion(_disc.nParType, nullptr);
	bool singleParDepSurfDiffusion = false;

	for (int type = 0; type < _disc.nParType; type++)
	{
		parDepSurfDiffusion[type] = _particle[type].getParDepSurfDiffusion();
		singleParDepSurfDiffusion = _particle[0].singleParDepSurfDiffusion();
		if (_particle[type].singleParDepSurfDiffusion() != singleParDepSurfDiffusion)
			throw InvalidParameterException("Something went wrong configuring the surface diffusion parameter dependence");
	}

	model::getAllParameterValues(data, parDepSurfDiffusion, singleParDepSurfDiffusion);

	return data;
}

double GeneralRateModelDG::getParameterDouble(const ParameterId& pId) const
{
	double val = 0.0;

	std::vector<IParameterStateDependence*> parDepSurfDiffusion(_disc.nParType, nullptr);
	bool singleParDepSurfDiffusion = false;

	for (int type = 0; type < _disc.nParType; type++)
	{
		parDepSurfDiffusion[type] = _particle[type].getParDepSurfDiffusion();
		singleParDepSurfDiffusion = _particle[0].singleParDepSurfDiffusion();
		if (_particle[type].singleParDepSurfDiffusion() != singleParDepSurfDiffusion)
			throw InvalidParameterException("Something went wrong configuring the surface diffusion parameter dependence");
	}

	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		if (model::getParameterDouble(pId, parDepSurfDiffusion, singleParDepSurfDiffusion, val))
			return val;
	}

	// Not found
	return UnitOperationBase::getParameterDouble(pId);
}

bool GeneralRateModelDG::hasParameter(const ParameterId& pId) const
{
	std::vector<IParameterStateDependence*> parDepSurfDiffusion(_disc.nParType, nullptr);
	bool singleParDepSurfDiffusion = false;

	for (int type = 0; type < _disc.nParType; type++)
	{
		parDepSurfDiffusion[type] = _particle[type].getParDepSurfDiffusion();
		singleParDepSurfDiffusion = _particle[0].singleParDepSurfDiffusion();
		if (_particle[type].singleParDepSurfDiffusion() != singleParDepSurfDiffusion)
			throw InvalidParameterException("Something went wrong configuring the surface diffusion parameter dependence");
	}

	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		if (model::hasParameter(pId, parDepSurfDiffusion, singleParDepSurfDiffusion))
			return true;
	}

	return UnitOperationBase::hasParameter(pId);
}

int GeneralRateModelDG::Exporter::writeMobilePhase(double* buffer) const
{
	const int blockSize = numMobilePhaseDofs();
	std::copy_n(_idx.c(_data), blockSize, buffer);
	return blockSize;
}

int GeneralRateModelDG::Exporter::writeSolidPhase(double* buffer) const
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

int GeneralRateModelDG::Exporter::writeParticleMobilePhase(double* buffer) const
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

int GeneralRateModelDG::Exporter::writeSolidPhase(unsigned int parType, double* buffer) const
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

int GeneralRateModelDG::Exporter::writeParticleMobilePhase(unsigned int parType, double* buffer) const
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

int GeneralRateModelDG::Exporter::writeParticleFlux(double* buffer) const
{
	return 0;
}

int GeneralRateModelDG::Exporter::writeParticleFlux(unsigned int parType, double* buffer) const
{
	return 0;
}

int GeneralRateModelDG::Exporter::writeInlet(unsigned int port, double* buffer) const
{
	cadet_assert(port == 0);
	std::copy_n(_data, _disc.nComp, buffer);
	return _disc.nComp;
}

int GeneralRateModelDG::Exporter::writeInlet(double* buffer) const
{
	std::copy_n(_data, _disc.nComp, buffer);
	return _disc.nComp;
}

int GeneralRateModelDG::Exporter::writeOutlet(unsigned int port, double* buffer) const
{
	cadet_assert(port == 0);

	if (_model._convDispOp.forwardFlow())
		std::copy_n(&_idx.c(_data, _disc.nPoints - 1, 0), _disc.nComp, buffer);
	else
		std::copy_n(&_idx.c(_data, 0, 0), _disc.nComp, buffer);

	return _disc.nComp;
}

int GeneralRateModelDG::Exporter::writeOutlet(double* buffer) const
{
	if (_model._convDispOp.forwardFlow())
		std::copy_n(&_idx.c(_data, _disc.nPoints - 1, 0), _disc.nComp, buffer);
	else
		std::copy_n(&_idx.c(_data, 0, 0), _disc.nComp, buffer);

	return _disc.nComp;
}

}  // namespace model

}  // namespace cadet

#include "model/GeneralRateModelDG-InitialConditions.cpp"
#include "model/GeneralRateModelDG-LinearSolver.cpp"

