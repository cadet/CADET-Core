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

#include "model/RadialGeneralRateModelDG.hpp"
#include "model/particle/GeneralRateParticle.hpp"
#include "BindingModelFactory.hpp"
#include "ParamReaderHelper.hpp"
#include "cadet/Exceptions.hpp"
#include "cadet/ExternalFunction.hpp"
#include "cadet/SolutionRecorder.hpp"
#include "ConfigurationHelper.hpp"
#include "model/BindingModel.hpp"
#include "model/ParameterDependence.hpp"
#include "SimulationTypes.hpp"
#include "linalg/DenseMatrix.hpp"
#include "linalg/Norms.hpp"
#include "linalg/Subset.hpp"
#include "model/parts/BindingCellKernel.hpp"
#include "model/parts/DGToolbox.hpp"

#include "AdUtils.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <algorithm>
#include <functional>
#include <numeric>

#include "ParallelSupport.hpp"
#ifdef CADET_PARALLELIZE
#include <tbb/parallel_for.h>
#endif

#define EIGEN_USE_MKL_ALL

using namespace Eigen;

namespace cadet
{

namespace model
{

RadialGeneralRateModelDG::RadialGeneralRateModelDG(UnitOpIdx unitOpIdx) : UnitOperationBase(unitOpIdx),
	_dynReactionBulk(nullptr),
	_analyticJac(true), _jacobianAdDirs(0), _factorizeJacobian(false), _tempState(nullptr),
	_initC(0), _initState(0), _initStateDot(0)
{
	// Multiple particle types are supported
	_singleBinding = false;
	_singleDynReaction = false;
}

RadialGeneralRateModelDG::~RadialGeneralRateModelDG() CADET_NOEXCEPT
{
	delete[] _tempState;
	delete _linearSolver;

	for (IParticleModel* pm : _particles)
		delete pm;
	_particles.clear();

	delete[] _disc.nBound;
	delete[] _disc.boundOffset;
	delete[] _disc.nParPoints;
	delete[] _disc.parTypeOffset;
}

unsigned int RadialGeneralRateModelDG::numDofs() const CADET_NOEXCEPT
{
	// Inlet DOFs: nComp
	// Bulk DOFs: nPoints * nComp (mobile phase c)
	// Particle DOFs: sum over particle types and bulk points
	return _disc.nComp + _disc.nPoints * _disc.nComp + _disc.parTypeOffset[_disc.nParType];
}

unsigned int RadialGeneralRateModelDG::numPureDofs() const CADET_NOEXCEPT
{
	// Bulk DOFs: nPoints * nComp
	// Particle DOFs: sum over particle types
	return _disc.nPoints * _disc.nComp + _disc.parTypeOffset[_disc.nParType];
}

bool RadialGeneralRateModelDG::usesAD() const CADET_NOEXCEPT
{
#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	// We always need AD if we want to check the analytical Jacobian
	return true;
#else
	// We only need AD if we are not computing the Jacobian analytically
	return !_analyticJac;
#endif
}

bool RadialGeneralRateModelDG::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper)
{
	const bool firstConfigCall = _tempState == nullptr;

	// Read discretization
	_disc.nComp = paramProvider.getInt("NCOMP");

	// Read number of particle types
	_disc.nParType = paramProvider.exists("NPARTYPE") ? paramProvider.getInt("NPARTYPE") : 1;
	if (_disc.nParType < 1)
		throw InvalidParameterException("Number of particle types must be at least 1");

	// Read bound states from first particle type (assume all have same)
	paramProvider.pushScope("particle_type_000");

	std::vector<int> nBound = paramProvider.getIntArray("NBOUND");
	if (nBound.size() < _disc.nComp)
		throw InvalidParameterException("Field NBOUND contains too few elements (NCOMP = " + std::to_string(_disc.nComp) + " required)");

	if (firstConfigCall)
		_disc.nBound = new unsigned int[_disc.nComp];
	std::copy_n(nBound.begin(), _disc.nComp, _disc.nBound);

	paramProvider.popScope();

	paramProvider.pushScope("discretization");

	if (firstConfigCall)
		_linearSolver = cadet::linalg::setLinearSolver(paramProvider.exists("LINEAR_SOLVER") ? paramProvider.getString("LINEAR_SOLVER") : "SparseLU");

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
		_disc.nElem = std::max(1u, paramProvider.getInt("NCOL") / _disc.nNodes);
	else
		throw InvalidParameterException("Specify field NELEM (or NCOL)");

	if (_disc.nElem < 1)
		throw InvalidParameterException("Number of radial elements must be at least 1!");

	_disc.nPoints = _disc.nNodes * _disc.nElem;

	// Radial DG always uses exact integration
	int polynomial_integration_mode = 1;

	// Precompute offsets and total number of bound states
	if (firstConfigCall)
		_disc.boundOffset = new unsigned int[_disc.nComp];
	_disc.boundOffset[0] = 0;
	for (unsigned int i = 1; i < _disc.nComp; ++i)
	{
		_disc.boundOffset[i] = _disc.boundOffset[i - 1] + _disc.nBound[i - 1];
	}
	_disc.strideBound = _disc.boundOffset[_disc.nComp - 1] + _disc.nBound[_disc.nComp - 1];

	// Determine whether analytic Jacobian should be used
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	const bool analyticJac = paramProvider.getBool("USE_ANALYTIC_JACOBIAN");
#else
	const bool analyticJac = false;
#endif

	// Allocate space for initial conditions
	_initC.resize(_disc.nElem * _disc.nNodes * _disc.nComp);

	// Create nonlinear solver for consistent initialization
	configureNonlinearSolver(paramProvider);

	paramProvider.popScope();

	// Configure convection-dispersion operator (bulk phase stride is just nComp)
	const unsigned int strideNode = _disc.nComp;
	const bool transportSuccess = _convDispOp.configureModelDiscretization(paramProvider, helper, _disc.nComp, polynomial_integration_mode, _disc.nElem, _disc.polyDeg, strideNode);

	_disc.curSection = -1;

	// Allocate and configure particle models
	if (firstConfigCall)
	{
		_disc.nParPoints = new unsigned int[_disc.nParType];
		_disc.parTypeOffset = new unsigned int[_disc.nParType + 1];
	}

	_particles.resize(_disc.nParType, nullptr);
	_initCp.resize(_disc.nParType);
	_initCs.resize(_disc.nParType);

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		// Create particle model (GENERAL_RATE_PARTICLE handles its own particle_type_XXX scope push/pop)
		_particles[parType] = helper.createParticleModel("GENERAL_RATE_PARTICLE");
		if (!_particles[parType])
			throw InvalidParameterException("Failed to create particle model for type " + std::to_string(parType));

		// Configure particle model discretization (expects paramProvider at unit_001 scope)
		if (!_particles[parType]->configureModelDiscretization(paramProvider, helper, _disc.nComp, parType, _disc.nParType, _disc.nComp))
			throw InvalidParameterException("Failed to configure particle model discretization for type " + std::to_string(parType));

		_disc.nParPoints[parType] = _particles[parType]->nDiscPoints();

		// Allocate initial conditions for this particle type
		_initCp[parType].resize(_disc.nPoints * _disc.nParPoints[parType] * _disc.nComp);
		_initCs[parType].resize(_disc.nPoints * _disc.nParPoints[parType] * _disc.strideBound);
	}

	// Compute particle type offsets
	_disc.parTypeOffset[0] = 0;
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		// Each bulk point has nParPoints[parType] particle shells with (nComp + strideBound) DOFs each
		const unsigned int parDofs = _disc.nPoints * _disc.nParPoints[parType] * (_disc.nComp + _disc.strideBound);
		_disc.parTypeOffset[parType + 1] = _disc.parTypeOffset[parType] + parDofs;
	}
	_disc.nTotalParPoints = _disc.parTypeOffset[_disc.nParType];

	// Jacobian allocation
	_jacInlet.resize(_disc.nNodes, 1);

	// Total Jacobian size: bulk + all particle DOFs
	const unsigned int jacSize = _disc.nPoints * _disc.nComp + _disc.parTypeOffset[_disc.nParType];
	_jac.resize(jacSize, jacSize);
	_jacDisc.resize(jacSize, jacSize);

	// Set whether analytic Jacobian is used
	useAnalyticJacobian(analyticJac);

	// Configure binding models (from base class)
	clearBindingModels();
	_binding.resize(_disc.nParType, nullptr);

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		const std::string parGroup = "particle_type_" + std::string(3 - std::to_string(parType).length(), '0') + std::to_string(parType);
		paramProvider.pushScope(parGroup);

		if (paramProvider.exists("adsorption_model"))
		{
			paramProvider.pushScope("adsorption");

			const std::string bindModelName = paramProvider.getString("ADSORPTION_MODEL");
			_binding[parType] = helper.createBindingModel(bindModelName);
			if (!_binding[parType])
				throw InvalidParameterException("Unknown binding model " + bindModelName);

			if (!_binding[parType]->configureModelDiscretization(paramProvider, _disc.nComp, _disc.nBound, _disc.boundOffset))
				throw InvalidParameterException("Failed to configure binding model for particle type " + std::to_string(parType));

			paramProvider.popScope();
		}
		else
		{
			_binding[parType] = helper.createBindingModel("NONE");
			_binding[parType]->configureModelDiscretization(paramProvider, _disc.nComp, _disc.nBound, _disc.boundOffset);
		}

		paramProvider.popScope();
	}

	// Configure dynamic reaction models
	clearDynamicReactionModels();
	_dynReaction.resize(_disc.nParType, nullptr);

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		const std::string parGroup = "particle_type_" + std::string(3 - std::to_string(parType).length(), '0') + std::to_string(parType);
		paramProvider.pushScope(parGroup);

		if (paramProvider.exists("reaction_model"))
		{
			paramProvider.pushScope("reaction");

			const std::string dynReactionName = paramProvider.getString("REACTION_MODEL");
			_dynReaction[parType] = helper.createDynamicReactionModel(dynReactionName);
			if (!_dynReaction[parType])
				throw InvalidParameterException("Unknown dynamic reaction model " + dynReactionName);

			if (!_dynReaction[parType]->configureModelDiscretization(paramProvider, _disc.nComp, _disc.nBound, _disc.boundOffset))
				throw InvalidParameterException("Failed to configure dynamic reaction model");

			paramProvider.popScope();
		}

		paramProvider.popScope();
	}

	// Allocate temporaries
	unsigned int maxWorkspace = 0;
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		const unsigned int ws = std::max({
			_disc.nComp + _disc.strideBound,
			_binding[parType] ? _binding[parType]->workspaceSize(_disc.nComp, _disc.strideBound, _disc.nBound) : 0u,
			_dynReaction[parType] ? _dynReaction[parType]->workspaceSize(_disc.nComp, _disc.strideBound, _disc.nBound) : 0u
		});
		maxWorkspace = std::max(maxWorkspace, ws);
	}

	const unsigned int nDof = numDofs();
	if (firstConfigCall || (_tempState == nullptr))
	{
		delete[] _tempState;
		_tempState = new double[nDof + _disc.nPoints * maxWorkspace];
	}

	return transportSuccess;
}

bool RadialGeneralRateModelDG::configure(IParameterProvider& paramProvider)
{
	_parameters.clear();

	// Read geometry parameters
	_colPorosity = paramProvider.getDouble("COL_POROSITY");

	// Register column porosity
	_parameters[makeParamId(hashString("COL_POROSITY"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_colPorosity;

	// Read particle type volume fractions
	if (_disc.nParType > 1)
	{
		readScalarParameterOrArray(_parTypeVolFrac, paramProvider, "PAR_TYPE_VOLFRAC", 1);
		if (_parTypeVolFrac.size() == _disc.nParType)
		{
			// Expand to all radial points
			_parTypeVolFrac.resize(_disc.nPoints * _disc.nParType, 1.0);
			for (unsigned int i = 1; i < _disc.nPoints; ++i)
				std::copy(_parTypeVolFrac.begin(), _parTypeVolFrac.begin() + _disc.nParType, _parTypeVolFrac.begin() + _disc.nParType * i);
		}

		if (_disc.nParType * _disc.nPoints != _parTypeVolFrac.size())
			throw InvalidParameterException("Number of elements in field PAR_TYPE_VOLFRAC does not match number of particle types times number of radial points");

		// Check that particle volume fractions sum to 1.0
		for (unsigned int i = 0; i < _disc.nPoints; ++i)
		{
			const double volFracSum = std::accumulate(_parTypeVolFrac.begin() + i * _disc.nParType, _parTypeVolFrac.begin() + (i + 1) * _disc.nParType, 0.0,
				[](double a, const active& b) -> double { return a + static_cast<double>(b); });
			if (std::abs(1.0 - volFracSum) > 1e-10)
				throw InvalidParameterException("Sum of field PAR_TYPE_VOLFRAC differs from 1.0 (is " + std::to_string(volFracSum) + ") in radial point " + std::to_string(i));
		}

		// Register parameters
		for (unsigned int i = 0; i < _disc.nParType; ++i)
			_parameters[makeParamId(hashString("PAR_TYPE_VOLFRAC"), _unitOpIdx, CompIndep, i, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parTypeVolFrac[i];
	}
	else if (_disc.nParType == 1)
		_parTypeVolFrac = std::vector<active>(_disc.nPoints, 1.0);

	// Configure particle models (GeneralRateParticle handles its own particle_type_XXX scope push/pop)
	unsigned int nBoundBeforeType = 0;
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		if (!_particles[parType]->configure(_unitOpIdx, paramProvider, _parameters, _disc.nParType, &nBoundBeforeType, _disc.strideBound))
			throw InvalidParameterException("Failed to configure particle model for type " + std::to_string(parType));

		nBoundBeforeType += _disc.strideBound;
	}

	// Configure binding models
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		if (_binding[parType] && _binding[parType]->requiresConfiguration())
		{
			const std::string parGroup = "particle_type_" + std::string(3 - std::to_string(parType).length(), '0') + std::to_string(parType);
			paramProvider.pushScope(parGroup);
			paramProvider.pushScope("adsorption");

			if (!_binding[parType]->configure(paramProvider, _unitOpIdx, parType))
				throw InvalidParameterException("Failed to configure binding model for particle type " + std::to_string(parType));

			paramProvider.popScope();
			paramProvider.popScope();
		}
	}

	// Configure dynamic reaction models
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		if (_dynReaction[parType])
		{
			const std::string parGroup = "particle_type_" + std::string(3 - std::to_string(parType).length(), '0') + std::to_string(parType);
			paramProvider.pushScope(parGroup);
			paramProvider.pushScope("reaction");

			if (!_dynReaction[parType]->configure(paramProvider, _unitOpIdx, parType))
				throw InvalidParameterException("Failed to configure dynamic reaction model for particle type " + std::to_string(parType));

			paramProvider.popScope();
			paramProvider.popScope();
		}
	}

	// Configure convection-dispersion operator
	const bool transportSuccess = _convDispOp.configure(_unitOpIdx, paramProvider, _parameters);

	// Set Jacobian pattern
	setPattern(_jac, false);
	setPattern(_jacDisc, true);

	return transportSuccess;
}

unsigned int RadialGeneralRateModelDG::threadLocalMemorySize() const CADET_NOEXCEPT
{
	return 0;
}

unsigned int RadialGeneralRateModelDG::requiredADdirs() const CADET_NOEXCEPT
{
	return _jacobianAdDirs;
}

void RadialGeneralRateModelDG::useAnalyticJacobian(const bool analyticJac)
{
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	_analyticJac = analyticJac;
	if (!_analyticJac)
	{
		// AD Jacobian required - compute from conv-disp bandwidth
		const unsigned int bulkBandwidth = _convDispOp.requiredADdirs();
		// Particle phase has local coupling within each particle
		unsigned int maxParStride = 0;
		for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		{
			const unsigned int parStride = _disc.nParPoints[parType] * (_disc.nComp + _disc.strideBound);
			maxParStride = std::max(maxParStride, parStride);
		}
		// Film diffusion couples bulk and particle surface
		_jacobianAdDirs = std::max(bulkBandwidth, 2 * _disc.nNodes * maxParStride + 1);
	}
	else
		_jacobianAdDirs = 0;
#else
	// Use AD Jacobian if analytic Jacobian is to be checked
	_analyticJac = false;
	const unsigned int bulkBandwidth = _convDispOp.requiredADdirs();
	unsigned int maxParStride = 0;
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		const unsigned int parStride = _disc.nParPoints[parType] * (_disc.nComp + _disc.strideBound);
		maxParStride = std::max(maxParStride, parStride);
	}
	_jacobianAdDirs = std::max(bulkBandwidth, 2 * _disc.nNodes * maxParStride + 1);
#endif
}

void RadialGeneralRateModelDG::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, const ConstSimulationState& simState, const AdJacobianParams& adJac)
{
	// Update section index
	updateSection(secIdx);

	// Update convection-dispersion operator for new section
	_convDispOp.notifyDiscontinuousSectionTransition(t, secIdx, _jacInlet);

	// Update particle models
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		_particles[parType]->notifyDiscontinuousSectionTransition(t, secIdx);
	}
}

void RadialGeneralRateModelDG::setFlowRates(active const* in, active const* out) CADET_NOEXCEPT
{
	_convDispOp.setFlowRates(in[0], out[0], _colPorosity);
}

void RadialGeneralRateModelDG::reportSolution(ISolutionRecorder& recorder, double const* const solution) const
{
	Exporter expr(_disc, *this, solution);
	recorder.beginUnitOperation(_unitOpIdx, *this, expr);
	recorder.endUnitOperation();
}

void RadialGeneralRateModelDG::reportSolutionStructure(ISolutionRecorder& recorder) const
{
	Exporter expr(_disc, *this, nullptr);
	recorder.unitOperationStructure(_unitOpIdx, *this, expr);
}

int RadialGeneralRateModelDG::jacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	if (_analyticJac)
		return residualImpl<double, double, double, true, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
	else
		return residualWithJacobian(simTime, ConstSimulationState{ simState.vecStateY, nullptr }, nullptr, adJac, threadLocalMem);
}

int RadialGeneralRateModelDG::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	return residualImpl<double, double, double, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
}

int RadialGeneralRateModelDG::residualWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	return residual(simTime, simState, res, adJac, threadLocalMem, true, false);
}

int RadialGeneralRateModelDG::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem, bool updateJacobian, bool paramSensitivity)
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
int RadialGeneralRateModelDG::residualImpl(double t, unsigned int secIdx, StateType const* const y_, double const* const yDot_, ResidualType* const res_, util::ThreadLocalStorage& threadLocalMem)
{
	Indexer idxr(_disc);

	// Compute Jacobian if requested
	if (wantJac) {
		if (!wantRes || _disc.newStaticJac) {
			_convDispOp.calcTransportJacobian(_jac, _jacInlet, 0);
			_disc.newStaticJac = false;
		}
	}

	// Initialize residual to zero
	if (wantRes)
	{
		std::fill(res_, res_ + numDofs(), ResidualType(0.0));
	}

	// Compute bulk convection dispersion residual
	if (wantRes)
		_convDispOp.residual(*this, t, secIdx, y_, yDot_, res_, typename cadet::ParamSens<ParamType>::enabled());

	BENCH_START(_timerResidualPar);

	// Compute particle residuals
	// Loop over particle types and bulk points
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
#ifdef CADET_PARALLELIZE
		tbb::parallel_for(std::size_t(0), static_cast<std::size_t>(_disc.nPoints), [&](std::size_t colNode)
#else
		for (unsigned int colNode = 0; colNode < _disc.nPoints; ++colNode)
#endif
		{
			LinearBufferAllocator tlmAlloc = threadLocalMem.get();

			linalg::BandedEigenSparseRowIterator jacIt(_jac, idxr.offsetCp(parType) - idxr.offsetC() + colNode * idxr.strideParBlock(parType));

			// Get column packing parameters (porosity, volume fraction, position)
			model::columnPackingParameters packing
			{
				_parTypeVolFrac[parType + _disc.nParType * colNode],
				_colPorosity,
				ColumnPosition{ _convDispOp.relativeCoordinate(colNode), 0.0, 0.0 }
			};

			// Call particle residual
			// Note: The particle model handles film diffusion BC at particle surface
			if constexpr (std::is_same_v<StateType, double> && std::is_same_v<ResidualType, double>)
			{
				_particles[parType]->residual(t, secIdx,
					y_ + idxr.offsetCp(parType) + colNode * idxr.strideParBlock(parType),
					y_ + idxr.offsetC() + colNode * idxr.strideColNode(),
					yDot_ ? yDot_ + idxr.offsetCp(parType) + colNode * idxr.strideParBlock(parType) : nullptr,
					wantRes ? res_ + idxr.offsetCp(parType) + colNode * idxr.strideParBlock(parType) : nullptr,
					wantRes ? res_ + idxr.offsetC() + colNode * idxr.strideColNode() : nullptr,
					packing, jacIt, tlmAlloc,
					WithoutParamSensitivity()
				);
			}
			else if constexpr (std::is_same_v<StateType, double> && std::is_same_v<ResidualType, active>)
			{
				_particles[parType]->residual(t, secIdx,
					y_ + idxr.offsetCp(parType) + colNode * idxr.strideParBlock(parType),
					y_ + idxr.offsetC() + colNode * idxr.strideColNode(),
					yDot_ ? yDot_ + idxr.offsetCp(parType) + colNode * idxr.strideParBlock(parType) : nullptr,
					wantRes ? res_ + idxr.offsetCp(parType) + colNode * idxr.strideParBlock(parType) : nullptr,
					wantRes ? res_ + idxr.offsetC() + colNode * idxr.strideColNode() : nullptr,
					packing, jacIt, tlmAlloc,
					WithParamSensitivity()
				);
			}
			else if constexpr (std::is_same_v<StateType, active> && std::is_same_v<ResidualType, active>)
			{
				_particles[parType]->residual(t, secIdx,
					y_ + idxr.offsetCp(parType) + colNode * idxr.strideParBlock(parType),
					y_ + idxr.offsetC() + colNode * idxr.strideColNode(),
					yDot_ ? yDot_ + idxr.offsetCp(parType) + colNode * idxr.strideParBlock(parType) : nullptr,
					wantRes ? res_ + idxr.offsetCp(parType) + colNode * idxr.strideParBlock(parType) : nullptr,
					wantRes ? res_ + idxr.offsetC() + colNode * idxr.strideColNode() : nullptr,
					packing, jacIt, tlmAlloc,
					WithoutParamSensitivity()
				);
			}

		} CADET_PARFOR_END;
	}

	BENCH_STOP(_timerResidualPar);

	// Add time derivatives for bulk phase
	if (wantRes && yDot_)
	{
		for (unsigned int node = 0; node < _disc.nPoints; ++node)
		{
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				const unsigned int bulkIdx = idxr.offsetC() + comp + node * idxr.strideColNode();
				res_[bulkIdx] += yDot_[bulkIdx];
			}
		}
	}

	return 0;
}

int RadialGeneralRateModelDG::residualSensFwdAdOnly(const SimulationTime& simTime, const ConstSimulationState& simState, active* const adRes, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidualSens);

	return residualImpl<double, active, active, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, adRes, threadLocalMem);
}

int RadialGeneralRateModelDG::residualSensFwdWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidualSens);

	return residual(simTime, simState, nullptr, adJac, threadLocalMem, true, true);
}

int RadialGeneralRateModelDG::residualSensFwdCombine(const SimulationTime& simTime, const ConstSimulationState& simState,
	const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS, active const* adRes,
	double* const tmp1, double* const tmp2, double* const tmp3)
{
	BENCH_SCOPE(_timerResidualSens);

	for (std::size_t param = 0; param < yS.size(); ++param)
	{
		multiplyWithJacobian(SimulationTime{ 0.0, 0u }, ConstSimulationState{ nullptr, nullptr }, yS[param], 1.0, 0.0, tmp1);
		multiplyWithDerivativeJacobian(SimulationTime{ 0.0, 0u }, ConstSimulationState{ nullptr, nullptr }, ySdot[param], tmp2);

		double* const ptrResS = resS[param];

		BENCH_START(_timerResidualSensPar);

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

void RadialGeneralRateModelDG::multiplyWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* yS, double alpha, double beta, double* ret)
{
	Indexer idxr(_disc);

	// Inlet DOFs
	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
	{
		ret[comp] = alpha * yS[comp] + beta * ret[comp];
	}

	// Pure DOFs (bulk + particles)
	Eigen::Map<Eigen::VectorXd> ret_vec(ret + idxr.offsetC(), numPureDofs());
	Eigen::Map<const Eigen::VectorXd> yS_vec(yS + idxr.offsetC(), numPureDofs());
	ret_vec = alpha * _jac * yS_vec + beta * ret_vec;

	// Inlet contribution
	unsigned int offInlet = _convDispOp.forwardFlow() ? 0 : (_disc.nElem - 1u) * idxr.strideColCell();

	for (unsigned int comp = 0; comp < _disc.nComp; comp++) {
		for (unsigned int node = 0; node < _disc.nNodes; node++) {
			ret[idxr.offsetC() + offInlet + node * idxr.strideColNode() + comp] += alpha * _jacInlet.coeff(node, 0) * yS[comp];
		}
	}
}

void RadialGeneralRateModelDG::multiplyWithDerivativeJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* sDot, double* ret)
{
	Indexer idxr(_disc);

	std::fill_n(ret, numDofs(), 0.0);

	// Bulk phase: dc/dt
	for (unsigned int node = 0; node < _disc.nPoints; ++node)
	{
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			const unsigned int idx = idxr.offsetC() + comp + node * idxr.strideColNode();
			ret[idx] = sDot[idx];
		}
	}

	// Particle phase: dcp/dt + Fp * sum_j dq_j/dt
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		const double invBetaP = 1.0 / static_cast<double>(_particles[parType]->getPorosity()) - 1.0;

		for (unsigned int colNode = 0; colNode < _disc.nPoints; ++colNode)
		{
			for (unsigned int parNode = 0; parNode < _disc.nParPoints[parType]; ++parNode)
			{
				const unsigned int parOffset = idxr.offsetCp(parType) + colNode * idxr.strideParBlock(parType) + parNode * idxr.strideParNode(parType);

				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				{
					const unsigned int cpIdx = parOffset + comp;
					ret[cpIdx] = sDot[cpIdx];

					// Add bound state contributions
					if (_binding[parType] && !_binding[parType]->reactionQuasiStationarity()[comp])
					{
						for (unsigned int bnd = 0; bnd < _disc.nBound[comp]; ++bnd)
						{
							const unsigned int qIdx = parOffset + _disc.nComp + _disc.boundOffset[comp] + bnd;
							ret[cpIdx] += invBetaP * sDot[qIdx];
						}
					}
				}

				// Bound states: dq/dt
				for (unsigned int bnd = 0; bnd < _disc.strideBound; ++bnd)
				{
					const unsigned int qIdx = parOffset + _disc.nComp + bnd;
					ret[qIdx] = sDot[qIdx];
				}
			}
		}
	}
}

void RadialGeneralRateModelDG::setExternalFunctions(IExternalFunction** extFuns, unsigned int size)
{
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		if (_binding[parType])
			_binding[parType]->setExternalFunctions(extFuns, size);
		if (_dynReaction[parType])
			_dynReaction[parType]->setExternalFunctions(extFuns, size);
	}
}

unsigned int RadialGeneralRateModelDG::localOutletComponentIndex(unsigned int port) const CADET_NOEXCEPT
{
	Indexer idxr(_disc);
	// Outlet is at end of column for forward flow, beginning for backward
	if (_convDispOp.forwardFlow())
		return idxr.offsetC() + (_disc.nPoints - 1) * idxr.strideColNode();
	else
		return idxr.offsetC();
}

unsigned int RadialGeneralRateModelDG::localOutletComponentStride(unsigned int port) const CADET_NOEXCEPT
{
	return 1;
}

unsigned int RadialGeneralRateModelDG::localInletComponentIndex(unsigned int port) const CADET_NOEXCEPT
{
	return 0;
}

unsigned int RadialGeneralRateModelDG::localInletComponentStride(unsigned int port) const CADET_NOEXCEPT
{
	return 1;
}

void RadialGeneralRateModelDG::expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut)
{
	// TODO: Implement proper error tolerance expansion
	std::fill_n(expandOut, numDofs(), errorSpec[0]);
}

int RadialGeneralRateModelDG::linearSolve(double t, double alpha, double tol, double* const rhs, double const* const weight, const ConstSimulationState& simState)
{
	BENCH_SCOPE(_timerLinearSolve);

	Indexer idxr(_disc);

	if (_factorizeJacobian)
	{
		// Assemble discretized Jacobian
		assembleDiscretizedJacobian(alpha, idxr);

		// Factorize
		_linearSolver->analyzePattern(_jacDisc);
		_linearSolver->factorize(_jacDisc);

		_factorizeJacobian = false;
	}

	// Solve
	Eigen::Map<Eigen::VectorXd> r(rhs + idxr.offsetC(), numPureDofs());
	r = _linearSolver->solve(r);

	return 0;
}

void RadialGeneralRateModelDG::assembleDiscretizedJacobian(double alpha, const Indexer& idxr)
{
	_jacDisc = _jac;

	// Add time derivatives to bulk phase
	for (unsigned int node = 0; node < _disc.nPoints; ++node)
	{
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			const unsigned int idx = comp + node * idxr.strideColNode();
			_jacDisc.coeffRef(idx, idx) += alpha;
		}
	}

	// Add time derivatives to particle phases
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		const double invBetaP = 1.0 / static_cast<double>(_particles[parType]->getPorosity()) - 1.0;
		const unsigned int poreOffset = _disc.nPoints * _disc.nComp + _disc.parTypeOffset[parType];

		for (unsigned int colNode = 0; colNode < _disc.nPoints; ++colNode)
		{
			for (unsigned int parNode = 0; parNode < _disc.nParPoints[parType]; ++parNode)
			{
				const unsigned int parOffset = poreOffset + colNode * idxr.strideParBlock(parType) - _disc.parTypeOffset[parType] + parNode * idxr.strideParNode(parType);

				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				{
					const unsigned int cpIdx = parOffset + comp;
					_jacDisc.coeffRef(cpIdx, cpIdx) += alpha;

					// Coupling with bound states
					for (unsigned int bnd = 0; bnd < _disc.nBound[comp]; ++bnd)
					{
						const unsigned int qIdx = parOffset + _disc.nComp + _disc.boundOffset[comp] + bnd;
						_jacDisc.coeffRef(cpIdx, qIdx) += alpha * invBetaP;
					}
				}

				// Bound state time derivatives
				for (unsigned int bnd = 0; bnd < _disc.strideBound; ++bnd)
				{
					const unsigned int qIdx = parOffset + _disc.nComp + bnd;
					_jacDisc.coeffRef(qIdx, qIdx) += alpha;
				}
			}
		}
	}
}

void RadialGeneralRateModelDG::prepareADvectors(const AdJacobianParams& adJac) const
{
	ad::prepareAdVectorSeedsForBandMatrix(adJac.adY + Indexer(_disc).offsetC(), adJac.adDirOffset, numPureDofs(), _jac.outerSize(), 0, numPureDofs());
}

void RadialGeneralRateModelDG::extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset)
{
	Indexer idxr(_disc);
	const int nDOFs = numPureDofs();
	const double* const adVec = reinterpret_cast<const double*>(adRes) + idxr.offsetC();

	for (int row = 0; row < _jac.rows(); row++)
	{
		for (Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(_jac, row); it; ++it)
		{
			const int col = it.col();
			it.valueRef() = adVec[row * (adDirOffset + nDOFs + 1) + adDirOffset + col];
		}
	}
}

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
void RadialGeneralRateModelDG::checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const
{
	// TODO: Implement Jacobian check
}
#endif

void RadialGeneralRateModelDG::applyInitialCondition(const SimulationState& simState) const
{
	Indexer idxr(_disc);

	// Inlet DOFs
	std::fill_n(simState.vecStateY, _disc.nComp, 0.0);

	// Bulk phase
	for (unsigned int i = 0; i < _disc.nPoints * _disc.nComp; ++i)
		simState.vecStateY[idxr.offsetC() + i] = static_cast<double>(_initC[i]);

	// Particle phases
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		for (unsigned int colNode = 0; colNode < _disc.nPoints; ++colNode)
		{
			for (unsigned int parNode = 0; parNode < _disc.nParPoints[parType]; ++parNode)
			{
				const unsigned int parOffset = idxr.offsetCp(parType) + colNode * idxr.strideParBlock(parType) + parNode * idxr.strideParNode(parType);
				const unsigned int initOffset = (colNode * _disc.nParPoints[parType] + parNode);

				// Pore phase
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
					simState.vecStateY[parOffset + comp] = static_cast<double>(_initCp[parType][initOffset * _disc.nComp + comp]);

				// Solid phase
				for (unsigned int bnd = 0; bnd < _disc.strideBound; ++bnd)
					simState.vecStateY[parOffset + _disc.nComp + bnd] = static_cast<double>(_initCs[parType][initOffset * _disc.strideBound + bnd]);
			}
		}
	}

	// Time derivatives
	if (simState.vecStateYdot)
		std::fill_n(simState.vecStateYdot, numDofs(), 0.0);
}

void RadialGeneralRateModelDG::readInitialCondition(IParameterProvider& paramProvider)
{
	Indexer idxr(_disc);

	// Bulk phase initial conditions
	std::vector<double> initC = paramProvider.getDoubleArray("INIT_C");
	if (initC.size() < _disc.nComp)
		throw InvalidParameterException("INIT_C has fewer than NCOMP (" + std::to_string(_disc.nComp) + ") entries");

	// Fill initial conditions for all nodes
	if (initC.size() >= _disc.nPoints * _disc.nComp)
	{
		// Per-node initialization
		for (unsigned int i = 0; i < _disc.nPoints * _disc.nComp; ++i)
			_initC[i] = initC[i];
	}
	else
	{
		// Per-component initialization
		for (unsigned int node = 0; node < _disc.nPoints; ++node)
		{
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				_initC[node * _disc.nComp + comp] = initC[comp];
		}
	}

	// Particle phase initial conditions per particle type
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		const std::string parGroup = "particle_type_" + std::string(3 - std::to_string(parType).length(), '0') + std::to_string(parType);
		paramProvider.pushScope(parGroup);

		// Pore phase
		if (paramProvider.exists("INIT_CP"))
		{
			std::vector<double> initCp = paramProvider.getDoubleArray("INIT_CP");
			const unsigned int nTotalPar = _disc.nPoints * _disc.nParPoints[parType];

			if (initCp.size() >= nTotalPar * _disc.nComp)
			{
				for (unsigned int i = 0; i < nTotalPar * _disc.nComp; ++i)
					_initCp[parType][i] = initCp[i];
			}
			else
			{
				for (unsigned int i = 0; i < nTotalPar; ++i)
				{
					for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
						_initCp[parType][i * _disc.nComp + comp] = initCp[comp % initCp.size()];
				}
			}
		}
		else
		{
			// Default: pore phase = bulk phase (at corresponding radial position)
			for (unsigned int colNode = 0; colNode < _disc.nPoints; ++colNode)
			{
				for (unsigned int parNode = 0; parNode < _disc.nParPoints[parType]; ++parNode)
				{
					const unsigned int initOffset = colNode * _disc.nParPoints[parType] + parNode;
					for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
						_initCp[parType][initOffset * _disc.nComp + comp] = _initC[colNode * _disc.nComp + comp];
				}
			}
		}

		// Solid phase
		if (_disc.strideBound > 0)
		{
			std::vector<double> initCs = paramProvider.getDoubleArray("INIT_Q");
			const unsigned int nTotalPar = _disc.nPoints * _disc.nParPoints[parType];

			if (initCs.size() >= nTotalPar * _disc.strideBound)
			{
				for (unsigned int i = 0; i < nTotalPar * _disc.strideBound; ++i)
					_initCs[parType][i] = initCs[i];
			}
			else
			{
				for (unsigned int i = 0; i < nTotalPar; ++i)
				{
					for (unsigned int bnd = 0; bnd < _disc.strideBound; ++bnd)
						_initCs[parType][i * _disc.strideBound + bnd] = initCs[bnd % initCs.size()];
				}
			}
		}

		paramProvider.popScope();
	}
}

void RadialGeneralRateModelDG::consistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerConsistentInit);

	// Apply initial conditions
	applyInitialCondition(SimulationState{ vecStateY, nullptr });

	// Solve for consistent initial state if binding is quasi-stationary
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		if (_binding[parType] && _binding[parType]->hasQuasiStationaryReactions())
		{
			// TODO: Implement quasi-stationary binding initialization
		}
	}
}

void RadialGeneralRateModelDG::consistentInitialTimeDerivative(const SimulationTime& simTime, double const* vecStateY, double* const vecStateYdot, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerConsistentInit);

	// Set time derivatives to zero initially
	std::fill_n(vecStateYdot, numDofs(), 0.0);
}

void RadialGeneralRateModelDG::initializeSensitivityStates(const std::vector<double*>& vecSensY) const
{
	for (double* sensY : vecSensY)
		std::fill_n(sensY, numDofs(), 0.0);
}

void RadialGeneralRateModelDG::consistentInitialSensitivity(const SimulationTime& simTime, const ConstSimulationState& simState,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes, util::ThreadLocalStorage& threadLocalMem)
{
	for (unsigned int i = 0; i < vecSensY.size(); ++i)
	{
		std::fill_n(vecSensY[i], numDofs(), 0.0);
		std::fill_n(vecSensYdot[i], numDofs(), 0.0);
	}
}

void RadialGeneralRateModelDG::leanConsistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, util::ThreadLocalStorage& threadLocalMem)
{
	consistentInitialState(simTime, vecStateY, adJac, errorTol, threadLocalMem);
}

void RadialGeneralRateModelDG::leanConsistentInitialTimeDerivative(double t, double const* const vecStateY, double* const vecStateYdot, double* const res, util::ThreadLocalStorage& threadLocalMem)
{
	std::fill_n(vecStateYdot, numDofs(), 0.0);
}

void RadialGeneralRateModelDG::leanConsistentInitialSensitivity(const SimulationTime& simTime, const ConstSimulationState& simState,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes, util::ThreadLocalStorage& threadLocalMem)
{
	consistentInitialSensitivity(simTime, simState, vecSensY, vecSensYdot, adRes, threadLocalMem);
}

bool RadialGeneralRateModelDG::setParameter(const ParameterId& pId, double value)
{
	if (pId.unitOperation == _unitOpIdx)
	{
		if (_convDispOp.setParameter(pId, value))
			return true;
	}

	const bool found = UnitOperationBase::setParameter(pId, value);
	if (!found && (pId.unitOperation == _unitOpIdx))
	{
		for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		{
			if (_particles[parType] && _particles[parType]->setParameter(pId, value))
				return true;
			if (_binding[parType] && _binding[parType]->setParameter(pId, value))
				return true;
			if (_dynReaction[parType] && _dynReaction[parType]->setParameter(pId, value))
				return true;
		}
	}

	return found;
}

bool RadialGeneralRateModelDG::setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue)
{
	if (pId.unitOperation == _unitOpIdx)
	{
		if (_convDispOp.setSensitiveParameter(_sensParams, pId, adDirection, adValue))
			return true;
	}

	return UnitOperationBase::setSensitiveParameter(pId, adDirection, adValue);
}

void RadialGeneralRateModelDG::setSensitiveParameterValue(const ParameterId& pId, double value)
{
	if (pId.unitOperation == _unitOpIdx)
	{
		if (_convDispOp.setSensitiveParameterValue(_sensParams, pId, value))
			return;
	}

	UnitOperationBase::setSensitiveParameterValue(pId, value);
}

// Jacobian pattern functions
void RadialGeneralRateModelDG::setPattern(Eigen::SparseMatrix<double, RowMajor>& mat, bool stateDer)
{
	std::vector<T> tripletList;

	Indexer idxr(_disc);

	// Estimate number of entries
	unsigned int conv_disp_entries = _convDispOp.nJacEntries(false);
	unsigned int particle_entries = 0;
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		// Film diffusion coupling + particle diffusion + binding
		particle_entries += _disc.nPoints * _particles[parType]->jacobianNNZperParticle();
	}

	tripletList.reserve(conv_disp_entries + particle_entries);

	// Add convection-dispersion pattern (bulk phase)
	convDispJacPattern(tripletList);

	// Add particle Jacobian pattern
	particleJacPattern(tripletList, _disc.curSection >= 0 ? _disc.curSection : 0);

	if (stateDer)
		stateDerPattern(tripletList);

	mat.setFromTriplets(tripletList.begin(), tripletList.end());
}

void RadialGeneralRateModelDG::convDispJacPattern(std::vector<T>& tripletList)
{
	_convDispOp.convDispJacPattern(tripletList);
}

void RadialGeneralRateModelDG::particleJacPattern(std::vector<T>& tripletList, unsigned int secIdx)
{
	Indexer idxr(_disc);
	const unsigned int bulkOffset = _disc.nPoints * _disc.nComp;

	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		const unsigned int parOffset = _disc.parTypeOffset[parType];

		for (unsigned int colNode = 0; colNode < _disc.nPoints; ++colNode)
		{
			const unsigned int localParOffset = bulkOffset + parOffset + colNode * idxr.strideParBlock(parType) - _disc.parTypeOffset[parType];
			const unsigned int bulkNodeOffset = colNode * idxr.strideColNode();

			_particles[parType]->setParJacPattern(tripletList, localParOffset, bulkNodeOffset, colNode, secIdx);
		}
	}
}

void RadialGeneralRateModelDG::stateDerPattern(std::vector<T>& tripletList)
{
	Indexer idxr(_disc);
	const unsigned int bulkOffset = _disc.nPoints * _disc.nComp;

	// Bulk phase time derivative
	for (unsigned int point = 0; point < _disc.nPoints; point++)
	{
		for (unsigned int comp = 0; comp < _disc.nComp; comp++)
		{
			const unsigned int idx = comp + point * idxr.strideColNode();
			tripletList.push_back(T(idx, idx, 0.0));
		}
	}

	// Particle phase time derivatives
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		const unsigned int parOffset = _disc.parTypeOffset[parType];

		for (unsigned int colNode = 0; colNode < _disc.nPoints; ++colNode)
		{
			for (unsigned int parNode = 0; parNode < _disc.nParPoints[parType]; ++parNode)
			{
				const unsigned int localParOffset = bulkOffset + parOffset + colNode * idxr.strideParBlock(parType) - _disc.parTypeOffset[parType] + parNode * idxr.strideParNode(parType);

				// Pore phase time derivative and coupling to bound states
				for (unsigned int comp = 0; comp < _disc.nComp; comp++)
				{
					const unsigned int cpIdx = localParOffset + comp;
					tripletList.push_back(T(cpIdx, cpIdx, 0.0));

					// Coupling with bound states
					for (unsigned int bnd = 0; bnd < _disc.nBound[comp]; ++bnd)
					{
						const unsigned int qIdx = localParOffset + _disc.nComp + _disc.boundOffset[comp] + bnd;
						tripletList.push_back(T(cpIdx, qIdx, 0.0));
					}
				}

				// Bound state time derivatives
				for (unsigned int bnd = 0; bnd < _disc.strideBound; ++bnd)
				{
					const unsigned int qIdx = localParOffset + _disc.nComp + bnd;
					tripletList.push_back(T(qIdx, qIdx, 0.0));
				}
			}
		}
	}
}

// Exporter implementations
unsigned int RadialGeneralRateModelDG::Exporter::numParticleMobilePhaseDofs() const CADET_NOEXCEPT
{
	unsigned int total = 0;
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		total += _disc.nComp * _disc.nPoints * _disc.nParPoints[parType];
	return total;
}

unsigned int RadialGeneralRateModelDG::Exporter::numSolidPhaseDofs() const CADET_NOEXCEPT
{
	unsigned int total = 0;
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
		total += _disc.strideBound * _disc.nPoints * _disc.nParPoints[parType];
	return total;
}

int RadialGeneralRateModelDG::Exporter::writeMobilePhase(double* buffer) const
{
	const double* data = _data + _idx.offsetC();
	for (unsigned int i = 0; i < _disc.nPoints; ++i)
	{
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			buffer[comp * _disc.nPoints + i] = data[i * _idx.strideColNode() + comp];
	}
	return _disc.nComp * _disc.nPoints;
}

int RadialGeneralRateModelDG::Exporter::writeSolidPhase(double* buffer) const
{
	int written = 0;
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		written += writeSolidPhase(parType, buffer + written);
	}
	return written;
}

int RadialGeneralRateModelDG::Exporter::writeParticleMobilePhase(double* buffer) const
{
	int written = 0;
	for (unsigned int parType = 0; parType < _disc.nParType; ++parType)
	{
		written += writeParticleMobilePhase(parType, buffer + written);
	}
	return written;
}

int RadialGeneralRateModelDG::Exporter::writeSolidPhase(unsigned int parType, double* buffer) const
{
	const double* data = _data + _idx.offsetCp(parType);
	const unsigned int nParPoints = _disc.nParPoints[parType];

	for (unsigned int colNode = 0; colNode < _disc.nPoints; ++colNode)
	{
		for (unsigned int parNode = 0; parNode < nParPoints; ++parNode)
		{
			for (unsigned int bnd = 0; bnd < _disc.strideBound; ++bnd)
			{
				const unsigned int srcIdx = colNode * _idx.strideParBlock(parType) + parNode * _idx.strideParNode(parType) + _disc.nComp + bnd;
				buffer[bnd * _disc.nPoints * nParPoints + colNode * nParPoints + parNode] = data[srcIdx];
			}
		}
	}
	return _disc.strideBound * _disc.nPoints * nParPoints;
}

int RadialGeneralRateModelDG::Exporter::writeParticleMobilePhase(unsigned int parType, double* buffer) const
{
	const double* data = _data + _idx.offsetCp(parType);
	const unsigned int nParPoints = _disc.nParPoints[parType];

	for (unsigned int colNode = 0; colNode < _disc.nPoints; ++colNode)
	{
		for (unsigned int parNode = 0; parNode < nParPoints; ++parNode)
		{
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				const unsigned int srcIdx = colNode * _idx.strideParBlock(parType) + parNode * _idx.strideParNode(parType) + comp;
				buffer[comp * _disc.nPoints * nParPoints + colNode * nParPoints + parNode] = data[srcIdx];
			}
		}
	}
	return _disc.nComp * _disc.nPoints * nParPoints;
}

int RadialGeneralRateModelDG::Exporter::writeInlet(unsigned int port, double* buffer) const
{
	cadet_assert(port == 0);
	std::copy_n(_data, _disc.nComp, buffer);
	return _disc.nComp;
}

int RadialGeneralRateModelDG::Exporter::writeInlet(double* buffer) const
{
	return writeInlet(0, buffer);
}

int RadialGeneralRateModelDG::Exporter::writeOutlet(unsigned int port, double* buffer) const
{
	cadet_assert(port == 0);

	if (_model._convDispOp.forwardFlow())
	{
		const double* data = _data + _idx.offsetC() + (_disc.nPoints - 1) * _idx.strideColNode();
		std::copy_n(data, _disc.nComp, buffer);
	}
	else
	{
		const double* data = _data + _idx.offsetC();
		std::copy_n(data, _disc.nComp, buffer);
	}
	return _disc.nComp;
}

int RadialGeneralRateModelDG::Exporter::writeOutlet(double* buffer) const
{
	return writeOutlet(0, buffer);
}

int RadialGeneralRateModelDG::Exporter::writePrimaryCoordinates(double* coords) const
{
	// Compute DG node coordinates for radial bulk discretization
	for (unsigned int i = 0; i < _disc.nElem; i++) {
		for (unsigned int j = 0; j < _disc.nNodes; j++) {
			// Mapping for radial coordinates: rho_left + 0.5 * deltaRho * (1 + xi_j)
			coords[i * _disc.nNodes + j] = _model._convDispOp.elemLeftBound(i) +
				0.5 * (static_cast<double>(_model._convDispOp.columnLength()) / static_cast<double>(_disc.nElem)) *
				(1.0 + _model._convDispOp.LGLnodes()[j]);
		}
	}
	return _disc.nPoints;
}

int RadialGeneralRateModelDG::Exporter::writeParticleCoordinates(unsigned int parType, double* coords) const
{
	return _model._particles[parType]->writeParticleCoordinates(coords);
}

} // namespace model
} // namespace cadet
