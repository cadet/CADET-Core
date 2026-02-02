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

#include "model/LumpedRateModelWithPores.hpp"
#include "BindingModelFactory.hpp"
#include "ParamReaderHelper.hpp"
#include "ParamReaderScopes.hpp"
#include "cadet/Exceptions.hpp"
#include "cadet/ExternalFunction.hpp"
#include "cadet/SolutionRecorder.hpp"
#include "ConfigurationHelper.hpp"
#include "model/BindingModel.hpp"
#include "model/ReactionModel.hpp"
#include "model/parts/BindingCellKernel.hpp"
#include "model/ParameterDependence.hpp"
#include "SimulationTypes.hpp"
#include "linalg/DenseMatrix.hpp"
#include "linalg/BandMatrix.hpp"
#include "linalg/Norms.hpp"

#include "Stencil.hpp"
#include "AdUtils.hpp"
#include "SensParamUtil.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <algorithm>
#include <functional>

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

template <typename Operator>
int schurComplementMultiplierLRMPores(void* userData, double const* x, double* z)
{
	typedef LumpedRateModelWithPores<Operator> LRMP;
	LRMP* const lrm = static_cast<LRMP*>(userData);
	return lrm->schurComplementMatrixVector(x, z);
}


template <typename ConvDispOperator>
LumpedRateModelWithPores<ConvDispOperator>::LumpedRateModelWithPores(UnitOpIdx unitOpIdx) : UnitOperationBase(unitOpIdx),
 _filmDiffDep(nullptr), _jacP(0), _jacPdisc(0), _jacPF(0), _jacFP(0), _jacInlet(), _analyticJac(true),
	_jacobianAdDirs(0), _factorizeJacobian(false), _tempState(nullptr), _initC(0), _initCp(0), _initCs(0),
	_initState(0), _initStateDot(0)
{
}

template <typename ConvDispOperator>
LumpedRateModelWithPores<ConvDispOperator>::~LumpedRateModelWithPores() CADET_NOEXCEPT
{
	delete[] _tempState;
	delete _filmDiffDep;
	clearDynamicReactionModels();
}

template <typename ConvDispOperator>
unsigned int LumpedRateModelWithPores<ConvDispOperator>::numDofs() const CADET_NOEXCEPT
{
	// Column bulk DOFs: nCol * nComp
	// Particle DOFs: nCol * nParType particles each having nComp (liquid phase) + sum boundStates (solid phase) DOFs
	// Flux DOFs: nCol * nComp * nParType
	// Inlet DOFs: nComp
	return _disc.nComp + _disc.nComp * _disc.nCol * (1 + _disc.nParType) + _disc.parTypeOffset[_disc.nParType];
}

template <typename ConvDispOperator>
unsigned int LumpedRateModelWithPores<ConvDispOperator>::numPureDofs() const CADET_NOEXCEPT
{
	// Column bulk DOFs: nCol * nComp
	// Particle DOFs: nCol * nParType particles each having nComp (liquid phase) + sum boundStates (solid phase) DOFs
	// Flux DOFs: nCol * nComp * nParType
	return _disc.nComp * _disc.nCol * (1 + _disc.nParType) + _disc.parTypeOffset[_disc.nParType];
}


template <typename ConvDispOperator>
bool LumpedRateModelWithPores<ConvDispOperator>::usesAD() const CADET_NOEXCEPT
{
#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	// We always need AD if we want to check the analytical Jacobian
	return true;
#else
	// We only need AD if we are not computing the Jacobian analytically
	return !_analyticJac;
#endif
}

template <typename ConvDispOperator>
bool LumpedRateModelWithPores<ConvDispOperator>::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper)
{
	// ==== Read discretization
	_disc.nComp = paramProvider.getInt("NCOMP");
	_disc.nParType = paramProvider.exists("NPARTYPE") ? paramProvider.getInt("NPARTYPE") : 1;
	if (_disc.nParType < 1)
		throw InvalidParameterException("Number of particle types must be >= 1 for the GENERAL_RATE_MODEL with arrow-head Finite Volume discretization");

	paramProvider.pushScope("discretization");
	_disc.nCol = paramProvider.getInt("NCOL");
	paramProvider.popScope();

	// Precompute offsets and total number of bound states (DOFs in solid phase)
	_disc.nBound = new unsigned int[_disc.nComp * _disc.nParType];

	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		paramProvider.pushScope("particle_type_" + std::string(3 - std::to_string(parType).length(), '0') + std::to_string(parType));

		std::vector<int> nBound = paramProvider.getIntArray("NBOUND");
		if (nBound.size() != _disc.nComp)
			throw InvalidParameterException("Field NBOUND contains too few elements (NCOMP = " + std::to_string(_disc.nComp) + " required)");
		std::copy_n(nBound.begin(), _disc.nComp, _disc.nBound + parType * _disc.nComp);

		paramProvider.popScope();
	}

	const unsigned int nTotalBound = std::accumulate(_disc.nBound, _disc.nBound + _disc.nComp * _disc.nParType, 0u);

	_disc.boundOffset = new unsigned int[_disc.nComp * _disc.nParType];
	_disc.strideBound = new unsigned int[_disc.nParType + 1];
	_disc.nBoundBeforeType = new unsigned int[_disc.nParType];
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

	// Configure particles
	_parGeomSurfToVol = std::vector<double>(_disc.nParType, SurfVolRatioSphere);
	clearBindingModels();
	_binding = std::vector<IBindingModel*>(_disc.nParType, nullptr);
	bool bindingConfSuccess = true;
	
	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		paramProvider.pushScope("particle_type_" + std::string(3 - std::to_string(parType).length(), '0') + std::to_string(parType));

		// ==== Construct and configure binding model

		std::string bindModelName = "NONE";
		if (paramProvider.exists("ADSORPTION_MODEL"))
			bindModelName = paramProvider.getString("ADSORPTION_MODEL");

		if (paramProvider.exists("BINDING_PARTYPE_DEPENDENT"))
			_singleBinding = !paramProvider.getInt("BINDING_PARTYPE_DEPENDENT");
		else
			_singleBinding = _disc.nParType == 1;

		_binding[parType] = helper.createBindingModel(bindModelName);
		if (!_binding[parType])
			throw InvalidParameterException("Unknown binding model " + bindModelName);

		MultiplexedScopeSelector scopeGuard(paramProvider, "adsorption", _binding[parType]->usesParamProviderInDiscretizationConfig());
		bindingConfSuccess = _binding[parType]->configureModelDiscretization(paramProvider, _disc.nComp, _disc.nBound + parType * _disc.nComp, _disc.boundOffset + parType * _disc.nComp) && bindingConfSuccess;

		// ==== Construct and configure dynamic particle reaction model

		// Set particle geometry
		if (paramProvider.exists("PAR_GEOM"))
		{
			std::string geom = paramProvider.getString("PAR_GEOM");
			if (geom == "SPHERE")
				_parGeomSurfToVol[parType] = SurfVolRatioSphere;
			else if (geom == "CYLINDER")
				_parGeomSurfToVol[parType] = SurfVolRatioCylinder;
			else if (geom == "SLAB")
				_parGeomSurfToVol[parType] = SurfVolRatioSlab;
			else
				throw InvalidParameterException("Unknown particle geometry type " + geom);
		}

		paramProvider.popScope();
	}

	// Precompute offsets of particle type DOFs
	_disc.parTypeOffset = new unsigned int[_disc.nParType + 1];
	_disc.parTypeOffset[0] = 0;
	for (unsigned int j = 1; j < _disc.nParType + 1; ++j)
	{
		_disc.parTypeOffset[j] = _disc.parTypeOffset[j-1] + (_disc.nComp + _disc.strideBound[j-1]) * _disc.nCol;
	}

	paramProvider.pushScope("discretization");
	// Determine whether analytic Jacobian should be used but don't set it right now.
	// We need to setup Jacobian matrices first.
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	const bool analyticJac = paramProvider.getBool("USE_ANALYTIC_JACOBIAN");
#else
	const bool analyticJac = false;
#endif

	// Initialize and configure GMRES for solving the Schur-complement
	_gmres.initialize(_disc.nCol * _disc.nComp * _disc.nParType, paramProvider.getInt("MAX_KRYLOV"), linalg::toOrthogonalization(paramProvider.getInt("GS_TYPE")), paramProvider.getInt("MAX_RESTARTS"));
	_gmres.matrixVectorMultiplier(&schurComplementMultiplierLRMPores<ConvDispOperator>, this);
	_schurSafety = paramProvider.getDouble("SCHUR_SAFETY");

	// Allocate space for initial conditions
	_initC.resize(_disc.nComp);
	_initCp.resize(_disc.nComp * _disc.nParType);
	_initCs.resize(nTotalBound);

	// Create nonlinear solver for consistent initialization
	configureNonlinearSolver(paramProvider);

	paramProvider.popScope();

	const bool transportSuccess = _convDispOp.configureModelDiscretization(paramProvider, helper, _disc.nComp, _disc.nCol);

	paramProvider.pushScope("particle_type_000");
	if (paramProvider.exists("FILM_DIFFUSION_DEP"))
	{
		const std::string paramDepName = paramProvider.getString("FILM_DIFFUSION_DEP");
		_filmDiffDep = helper.createParameterParameterDependence(paramDepName);
		if (!_filmDiffDep)
			throw InvalidParameterException("Unknown parameter dependence " + paramDepName + " in FILM_DIFFUSION_DEP");

		_filmDiffDep->configureModelDiscretization(paramProvider);
	}
	else
		_filmDiffDep = helper.createParameterParameterDependence("CONSTANT_ONE");
	paramProvider.popScope();

	// Allocate memory
	Indexer idxr(_disc);

	_jacInlet.resize(_disc.nComp);

	_jacP.resize(_disc.nParType);
	_jacPdisc.resize(_disc.nParType);
	for (unsigned int i = 0; i < _disc.nParType; ++i)
	{
		_jacPdisc[i].resize(_disc.nCol * (_disc.nComp + _disc.strideBound[i]), _disc.nComp + _disc.strideBound[i] - 1, _disc.nComp + _disc.strideBound[i] - 1);
		_jacP[i].resize(_disc.nCol * (_disc.nComp + _disc.strideBound[i]), _disc.nComp + _disc.strideBound[i] - 1, _disc.nComp + _disc.strideBound[i] - 1);
	}

	_jacPF.resize(_disc.nParType);
	_jacFP.resize(_disc.nParType);
	for (unsigned int i = 0; i < _disc.nParType; ++i)
	{
		_jacPF[i].resize(_disc.nComp * _disc.nCol);
		_jacFP[i].resize(_disc.nComp * _disc.nCol);
	}

	_jacCF.resize(_disc.nComp * _disc.nCol * _disc.nParType);
	_jacFC.resize(_disc.nComp * _disc.nCol * _disc.nParType);

	// Set whether analytic Jacobian is used
	useAnalyticJacobian(analyticJac);


	// ==== Construct and configure dynamic reaction model
	clearDynamicReactionModels();
	_reacParticle = std::vector<ReactionSystem*>(_disc.nParType, nullptr);
	ReactionSystem::create(_reacParticle);
	
	bool reactionConfSuccess = true;
	if (_disc.nParType > 0)
	{
		for (unsigned int par = 0; par < _disc.nParType; par++)
		{
			char particleScope[32];
			snprintf(particleScope, sizeof(particleScope), "particle_type_%03d", par);
			paramProvider.pushScope(particleScope); // particle_type_xxx

			if (paramProvider.exists("NREAC_CROSS_PHASE"))
			{
				int nReactions = paramProvider.getInt("NREAC_CROSS_PHASE");
				reactionConfSuccess = _reacParticle[par]->configureDiscretization("cross_phase",
					nReactions,
					_disc.nComp,
					_disc.nBound + par * _disc.nComp,
					_disc.boundOffset + par * _disc.nComp,
					paramProvider,
					helper) && reactionConfSuccess;

			}
			if (paramProvider.exists("NREAC_LIQUID"))
			{
				int nReactions = paramProvider.getInt("NREAC_LIQUID");
				reactionConfSuccess = _reacParticle[par]->configureDiscretization("liquid",
					nReactions,
					_disc.nComp,
					_disc.nBound + par * _disc.nComp,
					_disc.boundOffset + par * _disc.nComp,
					paramProvider,
					helper) && reactionConfSuccess;
			}
			if (paramProvider.exists("NREAC_SOLID"))
			{
				int nReactions = paramProvider.getInt("NREAC_SOLID");
				reactionConfSuccess = _reacParticle[par]->configureDiscretization("solid",
					nReactions,
					_disc.nComp,
					_disc.nBound + par * _disc.nComp,
					_disc.boundOffset + par * _disc.nComp,
					paramProvider,
					helper) && reactionConfSuccess;

			}
			paramProvider.popScope(); // particle_type_xxx

		}

	}
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

	// Setup the memory for tempState based on state vector
	_tempState = new double[numDofs()];

	return transportSuccess && bindingConfSuccess && reactionConfSuccess;
}

template <typename ConvDispOperator>
bool LumpedRateModelWithPores<ConvDispOperator>::configure(IParameterProvider& paramProvider)
{
	_parameters.clear();

	const bool transportSuccess = _convDispOp.configure(_unitOpIdx, paramProvider, _parameters);

	if (_filmDiffDep)
	{
		if (!_filmDiffDep->configure(paramProvider, _unitOpIdx, ParTypeIndep, BoundStateIndep, "FILM_DIFFUSION_DEP"))
			throw InvalidParameterException("Failed to configure film diffusion parameter dependency (FILM_DIFFUSION_DEP)");
	}

	// Read geometry parameters
	_colPorosity = paramProvider.getDouble("COL_POROSITY");
	_parRadius.resize(_disc.nParType);
	_parPorosity.resize(_disc.nParType);

	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		paramProvider.pushScope("particle_type_" + std::string(3 - std::to_string(parType).length(), '0') + std::to_string(parType));

		_parRadius[parType] = paramProvider.getDouble("PAR_RADIUS");

		if (paramProvider.exists("PAR_RADIUS_PARTYPE_DEPENDENT"))
			_singleParRadius = !paramProvider.getBool("PAR_RADIUS_PARTYPE_DEPENDENT");

		if (_singleParRadius && parType == 0)
			_parameters[makeParamId(hashStringRuntime("PAR_RADIUS"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parRadius[parType];
		else if (!_singleParRadius)
			_parameters[makeParamId(hashStringRuntime("PAR_RADIUS"), _unitOpIdx, CompIndep, parType, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parRadius[parType];

		_parPorosity[parType] = paramProvider.getDouble("PAR_POROSITY");

		if (paramProvider.exists("PAR_POROSITY_PARTYPE_DEPENDENT"))
			_singleParPorosity = !paramProvider.getBool("PAR_RADIUS_PARTYPE_DEPENDENT");

		if (_singleParPorosity && parType == 0)
			_parameters[makeParamId(hashStringRuntime("PAR_POROSITY"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parPorosity[parType];
		else if (!_singleParPorosity)
			_parameters[makeParamId(hashStringRuntime("PAR_POROSITY"), _unitOpIdx, CompIndep, parType, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parPorosity[parType];

		paramProvider.popScope();
	}

	// Check whether PAR_TYPE_VOLFRAC is required or not
	if ((_disc.nParType > 1) && !paramProvider.exists("PAR_TYPE_VOLFRAC"))
		throw InvalidParameterException("The required parameter \"PAR_TYPE_VOLFRAC\" was not found");

	if (paramProvider.exists("PAR_TYPE_VOLFRAC"))
	{
		readScalarParameterOrArray(_parTypeVolFrac, paramProvider, "PAR_TYPE_VOLFRAC", 1);
		if (_parTypeVolFrac.size() == _disc.nParType)
		{
			_axiallyConstantParTypeVolFrac = true;

			// Expand to all axial cells
			_parTypeVolFrac.resize(_disc.nCol * _disc.nParType, 1.0);
			for (unsigned int i = 1; i < _disc.nCol; ++i)
				std::copy(_parTypeVolFrac.begin(), _parTypeVolFrac.begin() + _disc.nParType, _parTypeVolFrac.begin() + _disc.nParType * i);
		}
		else
			_axiallyConstantParTypeVolFrac = false;
	}
	else
	{
		_parTypeVolFrac.resize(_disc.nCol, 1.0);
		_axiallyConstantParTypeVolFrac = false;
	}

	// Check that particle volume fractions sum to 1.0
	for (unsigned int i = 0; i < _disc.nCol; ++i)
	{
		const double volFracSum = std::accumulate(_parTypeVolFrac.begin() + i * _disc.nParType, _parTypeVolFrac.begin() + (i+1) * _disc.nParType, 0.0,
			[](double a, const active& b) -> double { return a + static_cast<double>(b); });
		if (std::abs(1.0 - volFracSum) > 1e-10)
			throw InvalidParameterException("Sum of field PAR_TYPE_VOLFRAC differs from 1.0 (is " + std::to_string(volFracSum) + ") in axial cell " + std::to_string(i));
	}

	// Read vectorial parameters (which may also be section dependent; transport)
	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		paramProvider.pushScope("particle_type_" + std::string(3 - std::to_string(parType).length(), '0') + std::to_string(parType));

		bool filmDiffParTypeDep = paramProvider.exists("FILM_DIFFUSION_PARTYPE_DEPENDENT") ? paramProvider.getInt("FILM_DIFFUSION_PARTYPE_DEPENDENT") : true;
		MultiplexMode tmpMode = parType > 0 ? _filmDiffusionMode : MultiplexMode::Independent;
		_filmDiffusionMode = readAndRegisterMultiplexCompSecParam(paramProvider, _parameters, _filmDiffusion, "FILM_DIFFUSION", _disc.nComp, parType, filmDiffParTypeDep, _unitOpIdx, _filmDiffusion.size());
		if (parType > 0 && tmpMode != _filmDiffusionMode)
			throw InvalidParameterException("Inconsistent film diffusion multiplex mode accross particle types " + std::to_string(parType - 1) + " and " + std::to_string(parType));

		if (paramProvider.exists("PORE_ACCESSIBILITY"))
		{
			if (parType > 0)
				tmpMode = _poreAccessFactorMode;
			bool poreAccessParTypeDep = paramProvider.exists("PORE_ACCESSIBILITY_PARTYPE_DEPENDENT") ? paramProvider.getBool("PORE_ACCESSIBILITY_PARTYPE_DEPENDENT") : true;
			_poreAccessFactorMode = readAndRegisterMultiplexCompSecParam(paramProvider, _parameters, _poreAccessFactor, "PORE_ACCESSIBILITY", _disc.nComp, parType, poreAccessParTypeDep, _unitOpIdx);
			if (parType > 0 && tmpMode != _poreAccessFactorMode)
				throw InvalidParameterException("Inconsistent pore access factor mode accross particle types " + std::to_string(parType - 1) + " and " + std::to_string(parType));
		}
		else if (parType == 0)
		{
			_poreAccessFactorMode = MultiplexMode::ComponentType;
			_poreAccessFactor = std::vector<cadet::active>(_disc.nComp * _disc.nParType, 1.0);
		}

		paramProvider.popScope();
	}

	// Check whether all sizes are matched
	if (_disc.nParType != _parRadius.size())
		throw InvalidParameterException("Number of elements in field PAR_RADIUS does not match number of particle types");
	if (_disc.nParType * _disc.nCol != _parTypeVolFrac.size())
		throw InvalidParameterException("Number of elements in field PAR_TYPE_VOLFRAC does not match number of particle types");
	if (_disc.nParType != _parPorosity.size())
		throw InvalidParameterException("Number of elements in field PAR_POROSITY does not match number of particle types");
	if ((_filmDiffusion.size() < _disc.nComp * _disc.nParType) || (_filmDiffusion.size() % (_disc.nComp * _disc.nParType) != 0))
		throw InvalidParameterException("Number of elements in field FILM_DIFFUSION is not a positive multiple of NCOMP * NPARTYPE (" + std::to_string(_disc.nComp * _disc.nParType) + ")");
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
		registerParam2DArray(_parameters, _parTypeVolFrac, [=, this](bool multi, unsigned cell, unsigned int type) { return makeParamId(hashString("PAR_TYPE_VOLFRAC"), _unitOpIdx, CompIndep, type, BoundStateIndep, ReactionIndep, cell); }, _disc.nParType);

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

	// Reconfigure binding model
	bool bindingConfSuccess = true;
	bool dynReactionConfSuccess = true;
	if (!_binding.empty())
	{
		if (_singleBinding)
		{
			if (_binding[0] && _binding[0]->requiresConfiguration())
			{
				paramProvider.pushScope("particle_type_000");
				paramProvider.pushScope("adsorption");
				bindingConfSuccess = _binding[0]->configure(paramProvider, _unitOpIdx, ParTypeIndep);
				paramProvider.popScope();
				paramProvider.popScope();
			}
		}
		else
		{
			for (unsigned int type = 0; type < _disc.nParType; ++type)
			{
				if (!_binding[type] || !_binding[type]->requiresConfiguration())
					continue;

				paramProvider.pushScope("particle_type_" + std::string(3 - std::to_string(type).length(), '0') + std::to_string(type));
				paramProvider.pushScope("adsorption");
				bindingConfSuccess = bindingConfSuccess && _binding[type]->configure(paramProvider, _unitOpIdx, type);
				paramProvider.popScope();
				paramProvider.popScope();
			}
		}
	}
	// Reconfigure reaction model
	if (paramProvider.exists("NREAC_LIQUID"))
		dynReactionConfSuccess = _reaction.configure("liquid", 0, _unitOpIdx, paramProvider) && dynReactionConfSuccess;

	for (unsigned int par = 0; par < _disc.nParType; par++)
	{
		char particleScope[32];
		snprintf(particleScope, sizeof(particleScope), "particle_type_%03d", par);

		if (paramProvider.exists(particleScope))
		{
			paramProvider.pushScope(particleScope);
			if (paramProvider.exists("NREAC_CROSS_PHASE"))
				dynReactionConfSuccess = _reacParticle[par]->configure("cross_phase", par, _unitOpIdx, paramProvider) && dynReactionConfSuccess;
			if (paramProvider.exists("NREAC_LIQUID"))
				dynReactionConfSuccess = _reacParticle[par]->configure("liquid", par, _unitOpIdx, paramProvider) && dynReactionConfSuccess;
			if (paramProvider.exists("NREAC_SOLID"))
				dynReactionConfSuccess = _reacParticle[par]->configure("solid", par, _unitOpIdx, paramProvider) && dynReactionConfSuccess;

			paramProvider.popScope();
		}
	}

	return transportSuccess && bindingConfSuccess && dynReactionConfSuccess;
}

template <typename ConvDispOperator>
unsigned int LumpedRateModelWithPores<ConvDispOperator>::threadLocalMemorySize() const CADET_NOEXCEPT
{
	LinearMemorySizer lms;

	// Handle all reactions
	for (unsigned int i = 0; i < _disc.nParType; ++i)
	{
		_reacParticle[i]->setWorkspaceRequirements("cross_phase",  _disc.nComp, _disc.strideBound[i], lms);
		_reacParticle[i]->setWorkspaceRequirements("liquid", _disc.nComp, _disc.strideBound[i], lms);
		_reacParticle[i]->setWorkspaceRequirements("solid",  _disc.nComp, _disc.strideBound[i], lms);

	}
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

template <typename ConvDispOperator>
unsigned int LumpedRateModelWithPores<ConvDispOperator>::numAdDirsForJacobian() const CADET_NOEXCEPT
{
	// We need as many directions as the highest bandwidth of the diagonal blocks:
	// The bandwidth of the column block depends on the size of the WENO stencil, whereas
	// the bandwidth of the particle blocks are given by the number of components and bound states.

	// Get maximum stride of particle type blocks
	int maxStride = 0;
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		maxStride = std::max(maxStride, _jacP[type].stride());
	}

	return std::max(_convDispOp.requiredADdirs(), maxStride);
}

template <typename ConvDispOperator>
void LumpedRateModelWithPores<ConvDispOperator>::useAnalyticJacobian(const bool analyticJac)
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

template <typename ConvDispOperator>
void LumpedRateModelWithPores<ConvDispOperator>::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, const ConstSimulationState& simState, const AdJacobianParams& adJac)
{
	// Setup flux Jacobian blocks at the beginning of the simulation or in case of
	// section dependent film or particle diffusion coefficients
	if ((secIdx == 0) || isSectionDependent(_filmDiffusionMode))
		assembleOffdiagJac(t, secIdx);

	Indexer idxr(_disc);

	// AxialConvectionDispersionOperator tells us whether flow direction has changed
	if (!_convDispOp.notifyDiscontinuousSectionTransition(t, secIdx, adJac))
		return;

	// Setup the matrix connecting inlet DOFs to first column cells
	_jacInlet.clear();
	const double v = _convDispOp.inletJacobianFactor();

	if (_convDispOp.forwardFlow())
	{
		// Forwards flow

		// Place entries for inlet DOF to first column cell conversion
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			_jacInlet.addElement(comp * idxr.strideColComp(), comp, -v);
	}
	else
	{
		// Backwards flow

		// Place entries for inlet DOF to last column cell conversion
		const unsigned int offset = (_disc.nCol - 1) * idxr.strideColCell();
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			_jacInlet.addElement(offset + comp * idxr.strideColComp(), comp, v);
	}
}

template <typename ConvDispOperator>
void LumpedRateModelWithPores<ConvDispOperator>::setFlowRates(active const* in, active const* out) CADET_NOEXCEPT
{
	_convDispOp.setFlowRates(in[0], out[0], _colPorosity);
}

template <typename ConvDispOperator>
void LumpedRateModelWithPores<ConvDispOperator>::reportSolution(ISolutionRecorder& recorder, double const* const solution) const
{
	Exporter expr(_disc, *this, solution);
	recorder.beginUnitOperation(_unitOpIdx, *this, expr);
	recorder.endUnitOperation();
}

template <typename ConvDispOperator>
void LumpedRateModelWithPores<ConvDispOperator>::reportSolutionStructure(ISolutionRecorder& recorder) const
{
	Exporter expr(_disc, *this, nullptr);
	recorder.unitOperationStructure(_unitOpIdx, *this, expr);
}


template <typename ConvDispOperator>
unsigned int LumpedRateModelWithPores<ConvDispOperator>::requiredADdirs() const CADET_NOEXCEPT
{
	const unsigned int numDirsBinding = maxBindingAdDirs();
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	return numDirsBinding + _jacobianAdDirs;
#else
	// If CADET_CHECK_ANALYTIC_JACOBIAN is active, we always need the AD directions for the Jacobian
	return numDirsBinding + numAdDirsForJacobian();
#endif
}

template <typename ConvDispOperator>
void LumpedRateModelWithPores<ConvDispOperator>::prepareADvectors(const AdJacobianParams& adJac) const
{
	// Early out if AD is disabled
	if (!adJac.adY)
		return;

	Indexer idxr(_disc);

	// Column block
	_convDispOp.prepareADvectors(adJac);

	// Particle block
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		const unsigned int lowerParBandwidth = _jacP[type].lowerBandwidth();
		const unsigned int upperParBandwidth = _jacP[type].upperBandwidth();

		ad::prepareAdVectorSeedsForBandMatrix(adJac.adY + idxr.offsetCp(ParticleTypeIndex{type}), adJac.adDirOffset, idxr.strideParBlock(type) * _disc.nCol, lowerParBandwidth, upperParBandwidth, lowerParBandwidth);
	}
}

/**
 * @brief Extracts the system Jacobian from band compressed AD seed vectors
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
template <typename ConvDispOperator>
void LumpedRateModelWithPores<ConvDispOperator>::extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset)
{
	Indexer idxr(_disc);

	// Column
	_convDispOp.extractJacobianFromAD(adRes, adDirOffset);

	// Particles
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		linalg::BandMatrix& jacMat = _jacP[type];
		ad::extractBandedJacobianFromAd(adRes + idxr.offsetCp(ParticleTypeIndex{type}), adDirOffset, jacMat.lowerBandwidth(), jacMat);
	}
}

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN

/**
 * @brief Compares the analytical Jacobian with a Jacobian derived by AD
 * @details The analytical Jacobian is assumed to be stored in the corresponding band matrices.
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
template <typename ConvDispOperator>
void LumpedRateModelWithPores<ConvDispOperator>::checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const
{
	Indexer idxr(_disc);

	LOG(Debug) << "AD dir offset: " << adDirOffset << " DiagDirCol: " << _convDispOp.jacobian().lowerBandwidth();

	// Column
	const double maxDiffCol = _convDispOp.checkAnalyticJacobianAgainstAd(adRes, adDirOffset);

	// Particles
	double maxDiffPar = 0.0;
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		const linalg::BandMatrix& jacMat = _jacP[type];
		const double localDiff = ad::compareBandedJacobianWithAd(adRes + idxr.offsetCp(ParticleTypeIndex{type}), adDirOffset, jacMat.lowerBandwidth(), jacMat);
		LOG(Debug) << "-> Par type " << type << " diff: " << localDiff;
		maxDiffPar = std::max(maxDiffPar, localDiff);
	}
}

#endif

template <typename ConvDispOperator>
int LumpedRateModelWithPores<ConvDispOperator>::jacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	_factorizeJacobian = true;

	if (_analyticJac)
		return residualImpl<double, double, double, true, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
	else
		return residualWithJacobian(simTime, ConstSimulationState{ simState.vecStateY, nullptr }, nullptr, adJac, threadLocalMem);
}

template <typename ConvDispOperator>
int LumpedRateModelWithPores<ConvDispOperator>::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	// Evaluate residual do not compute Jacobian or parameter sensitivities
	return residualImpl<double, double, double, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
}

template <typename ConvDispOperator>
int LumpedRateModelWithPores<ConvDispOperator>::residualWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	// Evaluate residual, use AD for Jacobian if required but do not evaluate parameter derivatives
	return residual(simTime, simState, res, adJac, threadLocalMem, true, false);
}

template <typename ConvDispOperator>
int LumpedRateModelWithPores<ConvDispOperator>::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res,
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
			// Compute Jacobian via AD

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
		// Compute Jacobian via AD

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

template <typename ConvDispOperator>
template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
int LumpedRateModelWithPores<ConvDispOperator>::residualImpl(double t, unsigned int secIdx, StateType const* const y, double const* const yDot, ResidualType* const res, util::ThreadLocalStorage& threadLocalMem)
{
	// Reset Jacobian
	if (wantJac)
	{
		for (unsigned int type = 0; type < _disc.nParType; ++type)
			_jacP[type].setAll(0.0);
	}

	BENCH_START(_timerResidualPar);

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(std::size_t(0), static_cast<std::size_t>(_disc.nCol * _disc.nParType + 1), [&](std::size_t pblk)
#else
	for (unsigned int pblk = 0; pblk < _disc.nCol * _disc.nParType + 1; ++pblk)
#endif
	{
		if (cadet_unlikely(pblk == 0))
			residualBulk<StateType, ResidualType, ParamType, wantJac, wantRes>(t, secIdx, y, yDot, res, threadLocalMem);
		else
		{
			const unsigned int type = (pblk - 1) / _disc.nCol;
			const unsigned int par = (pblk - 1) % _disc.nCol;
			residualParticle<StateType, ResidualType, ParamType, wantJac, wantRes>(t, type, par, secIdx, y, yDot, res, threadLocalMem);
		}
	} CADET_PARFOR_END;

	BENCH_STOP(_timerResidualPar);

	if (!wantRes)
		return 0;

	residualFlux<StateType, ResidualType, ParamType>(t, secIdx, y, yDot, res);

	// Handle inlet DOFs, which are simply copied to res
	for (unsigned int i = 0; i < _disc.nComp; ++i)
	{
		res[i] = y[i];
	}

	return 0;
}

template <typename ConvDispOperator>
template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
int LumpedRateModelWithPores<ConvDispOperator>::residualBulk(double t, unsigned int secIdx, StateType const* yBase, double const* yDotBase, ResidualType* resBase, util::ThreadLocalStorage& threadLocalMem)
{
	if (wantRes)
		_convDispOp.residual(*this, t, secIdx, yBase, yDotBase, resBase, wantJac, typename ParamSens<ParamType>::enabled());
	else
		_convDispOp.jacobian(*this, t, secIdx, yBase, nullptr, nullptr);

	if (_reaction.getDynReactionVector("liquid").empty() || !_reaction.getDynReactionVector("liquid")[0])
		return 0;

	// Get offsets
	Indexer idxr(_disc);
	StateType const* y = yBase + idxr.offsetC();
	ResidualType* res = wantRes ? resBase + idxr.offsetC() : nullptr;
	LinearBufferAllocator tlmAlloc = threadLocalMem.get();

	for (unsigned int col = 0; col < _disc.nCol; ++col, y += idxr.strideColCell(), res += idxr.strideColCell())
	{
		const ColumnPosition colPos{ (0.5 + static_cast<double>(col)) / static_cast<double>(_disc.nCol), 0.0, 0.0 };
		
		for (auto i = 0; i < _reaction.getDynReactionVector("liquid").size(); i++)
		{

			if (!_reaction.getDynReactionVector("liquid")[i])
				continue;

		if (wantRes)
				_reaction.getDynReactionVector("liquid")[i]->residualFluxAdd(t, secIdx, colPos, _disc.nComp, y, res, -1.0, tlmAlloc);

		if (wantJac)
		{
			// static_cast should be sufficient here, but this statement is also analyzed when wantJac = false
				_reaction.getDynReactionVector("liquid")[i]->analyticJacobianAdd(t, secIdx, colPos, _disc.nComp, reinterpret_cast<double const*>(y), -1.0, _convDispOp.jacobian().row(col * idxr.strideColCell()), tlmAlloc);
			}
		}
	}

	return 0;
}

template <typename ConvDispOperator>
template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
int LumpedRateModelWithPores<ConvDispOperator>::residualParticle(double t, unsigned int parType, unsigned int colCell, unsigned int secIdx, StateType const* yBase, double const* yDotBase, ResidualType* resBase, util::ThreadLocalStorage& threadLocalMem)
{
	Indexer idxr(_disc);

	// Go to the particle block of the given type and column cell
	StateType const* y = yBase + idxr.offsetCp(ParticleTypeIndex{parType}, ParticleIndex{colCell});
	double const* yDot = yDotBase ? yDotBase + idxr.offsetCp(ParticleTypeIndex{parType}, ParticleIndex{colCell}) : nullptr;
	ResidualType* res = wantRes ? resBase + idxr.offsetCp(ParticleTypeIndex{parType}, ParticleIndex{colCell}) : nullptr;

	// Prepare parameters
	const ParamType radius = static_cast<ParamType>(_parRadius[parType]);

	// Midpoint of current column cell (z coordinate) - needed in externally dependent adsorption kinetic
	const double z = _convDispOp.relativeCoordinate(colCell);

	const parts::cell::CellParameters cellResParams
		{
			_disc.nComp,
			_disc.nBound + _disc.nComp * parType,
			_disc.boundOffset + _disc.nComp * parType,
			_disc.strideBound[parType],
			_binding[parType]->reactionQuasiStationarity(),
			_parPorosity[parType],
			_poreAccessFactor.data() + _disc.nComp * parType,
			_binding[parType],
			_reacParticle[parType]
		};

	// Handle time derivatives, binding, dynamic reactions
	if (wantRes)
		parts::cell::residualKernel<StateType, ResidualType, ParamType, parts::cell::CellParameters, linalg::BandMatrix::RowIterator, wantJac, true>(
			t, secIdx, ColumnPosition{ z, 0.0, static_cast<double>(radius) * 0.5 }, y, yDot, res,
			_jacP[parType].row(colCell * idxr.strideParBlock(parType)), cellResParams, threadLocalMem.get()
		);
	else
		parts::cell::residualKernel<StateType, ResidualType, ParamType, parts::cell::CellParameters, linalg::BandMatrix::RowIterator, wantJac, false, false>(
			t, secIdx, ColumnPosition{ z, 0.0, static_cast<double>(radius) * 0.5 }, y, yDot, res,
			_jacP[parType].row(colCell * idxr.strideParBlock(parType)), cellResParams, threadLocalMem.get()
		);

	return 0;
}

template <typename ConvDispOperator>
template <typename StateType, typename ResidualType, typename ParamType>
int LumpedRateModelWithPores<ConvDispOperator>::residualFlux(double t, unsigned int secIdx, StateType const* yBase, double const* yDotBase, ResidualType* resBase)
{
	Indexer idxr(_disc);

	const ParamType invBetaC = 1.0 / static_cast<ParamType>(_colPorosity) - 1.0;

	// Get offsets
	ResidualType* const resCol = resBase + idxr.offsetC();
	ResidualType* const resFlux = resBase + idxr.offsetJf();

	StateType const* const yCol = yBase + idxr.offsetC();
	StateType const* const yFlux = yBase + idxr.offsetJf();

	// J_f block (identity matrix), adds flux state to flux equation
	for (unsigned int i = 0; i < _disc.nComp * _disc.nCol * _disc.nParType; ++i)
		resFlux[i] = yFlux[i];

	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		ResidualType* const resParType = resBase + idxr.offsetCp(ParticleTypeIndex{type});
		ResidualType* const resFluxType = resBase + idxr.offsetJf(ParticleTypeIndex{type});

		StateType const* const yParType = yBase + idxr.offsetCp(ParticleTypeIndex{type});
		StateType const* const yFluxType = yBase + idxr.offsetJf(ParticleTypeIndex{type});

		const ParamType epsP = static_cast<ParamType>(_parPorosity[type]);
		const ParamType radius = static_cast<ParamType>(_parRadius[type]);
		active const* const filmDiff = getSectionDependentSlice(_filmDiffusion, _disc.nComp * _disc.nParType, secIdx) + type * _disc.nComp;
		active const* const poreAccFactor = _poreAccessFactor.data() + type * _disc.nComp;

		const ParamType jacCF_val = invBetaC * _parGeomSurfToVol[type] / radius;
		const ParamType jacPF_val = -_parGeomSurfToVol[type] / (epsP * radius);

		// J_{0,f} block, adds flux to column void / bulk volume equations
		for (unsigned int i = 0; i < _disc.nCol * _disc.nComp; ++i)
		{
			const unsigned int colCell = i / _disc.nComp;
			const unsigned int comp = i % _disc.nComp;

			const double relPos = _convDispOp.relativeCoordinate(colCell);
			const active curVelocity = _convDispOp.currentVelocity(relPos);
			const active modifier = _filmDiffDep->getValue(*this, ColumnPosition{relPos, 0.0, 0.0}, comp, ParTypeIndep, BoundStateIndep, curVelocity);

			resCol[i] += jacCF_val * static_cast<ParamType>(filmDiff[comp]) * static_cast<ParamType>(modifier) * static_cast<ParamType>(_parTypeVolFrac[type + _disc.nParType * colCell]) * yFluxType[i];
		}

		// J_{f,0} block, adds bulk volume state c_i to flux equation
		for (unsigned int bnd = 0; bnd < _disc.nCol; ++bnd)
		{
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				const unsigned int eq = bnd * idxr.strideColCell() + comp * idxr.strideColComp();
				resFluxType[eq] -= yCol[eq];
			}
		}

		// J_{p,f} block, adds flux to particle / bead volume equations
		for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
		{
			const double relPos = _convDispOp.relativeCoordinate(pblk);
			const active curVelocity = _convDispOp.currentVelocity(relPos);

			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				const active modifier = _filmDiffDep->getValue(*this, ColumnPosition{relPos, 0.0, 0.0}, comp, ParTypeIndep, BoundStateIndep, curVelocity);

				const unsigned int eq = pblk * idxr.strideColCell() + comp * idxr.strideColComp();
				resParType[pblk * idxr.strideParBlock(type) + comp] += jacPF_val / static_cast<ParamType>(poreAccFactor[comp]) * static_cast<ParamType>(filmDiff[comp]) * static_cast<ParamType>(modifier) * yFluxType[eq];
			}
		}

		// J_{f,p} block, adds outer bead shell state c_{p,i} to flux equation
		for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
		{
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				const unsigned int eq = pblk * idxr.strideColCell() + comp * idxr.strideColComp();
				resFluxType[eq] += yParType[comp + pblk * idxr.strideParBlock(type)];
			}
		}
	}

	return 0;
}

/**
 * @brief Assembles off diagonal Jacobian blocks
 * @details Assembles the fixed blocks @f$ J_{0,f}, \dots, J_{N_p,f} @f$ and @f$ J_{f,0}, \dots, J_{f, N_p}. @f$
 *          The blocks are fixed for each section.
 * @param [in] t Current time
 * @param [in] secIdx Index of the current section
 */
template <typename ConvDispOperator>
void LumpedRateModelWithPores<ConvDispOperator>::assembleOffdiagJac(double t, unsigned int secIdx)
{
	// Clear matrices for new assembly
	_jacCF.clear();
	_jacFC.clear();
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		_jacPF[type].clear();
		_jacFP[type].clear();
	}

	Indexer idxr(_disc);

	const double invBetaC = 1.0 / static_cast<double>(_colPorosity) - 1.0;

	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		const unsigned int typeOffset = type * _disc.nCol * _disc.nComp;

		const double epsP = static_cast<double>(_parPorosity[type]);
		const double radius = static_cast<double>(_parRadius[type]);
		const double jacCF_val = invBetaC * _parGeomSurfToVol[type] / radius;
		const double jacPF_val = -_parGeomSurfToVol[type] / (radius * epsP);

		active const* const filmDiff = getSectionDependentSlice(_filmDiffusion, _disc.nComp * _disc.nParType, secIdx) + type * _disc.nComp;
		active const* const poreAccFactor = _poreAccessFactor.data() + type * _disc.nComp;

		// Note that the J_f block, which is the identity matrix, is treated in the linear solver

		// J_{0,f} block, adds flux to column void / bulk volume equations
		for (unsigned int eq = 0; eq < _disc.nCol * _disc.nComp; ++eq)
		{
			const unsigned int colCell = eq / _disc.nComp;
			const unsigned int comp = eq % _disc.nComp;

			const double relPos = _convDispOp.relativeCoordinate(colCell);
			const double curVelocity = static_cast<double>(_convDispOp.currentVelocity(relPos));
			const double modifier = _filmDiffDep->getValue(*this, ColumnPosition{relPos, 0.0, 0.0}, comp, ParTypeIndep, BoundStateIndep, curVelocity);

			// Main diagonal corresponds to j_{f,i} (flux) state variable
			_jacCF.addElement(eq, eq + typeOffset, jacCF_val * static_cast<double>(filmDiff[comp]) * modifier * static_cast<double>(_parTypeVolFrac[type + _disc.nParType * colCell]));
		}

		// J_{f,0} block, adds bulk volume state c_i to flux equation
		for (unsigned int eq = 0; eq < _disc.nCol * _disc.nComp; ++eq)
		{
			_jacFC.addElement(eq + typeOffset, eq, -1.0);
		}

		// J_{p,f} block, implements bead boundary condition in outer bead shell equation
		for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
		{
			const double relPos = _convDispOp.relativeCoordinate(pblk);
			const double curVelocity = static_cast<double>(_convDispOp.currentVelocity(relPos));

			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				const unsigned int eq = typeOffset + pblk * idxr.strideColCell() + comp * idxr.strideColComp();
				const unsigned int col = pblk * idxr.strideParBlock(type) + comp;

				const double modifier = _filmDiffDep->getValue(*this, ColumnPosition{relPos, 0.0, 0.0}, comp, ParTypeIndep, BoundStateIndep, curVelocity);
				_jacPF[type].addElement(col, eq, jacPF_val / static_cast<double>(poreAccFactor[comp]) * static_cast<double>(filmDiff[comp]) * modifier);
			}
		}

		// J_{f,p} block, adds outer bead shell state c_{p,i} to flux equation
		for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
		{
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				const unsigned int eq = typeOffset + pblk * idxr.strideColCell() + comp * idxr.strideColComp();
				const unsigned int col = pblk * idxr.strideParBlock(type) + comp;
				_jacFP[type].addElement(eq, col, 1.0);
			}
		}
	}
}

template <typename ConvDispOperator>
int LumpedRateModelWithPores<ConvDispOperator>::residualSensFwdWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState,
	const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidualSens);

	// Evaluate residual for all parameters using AD in vector mode and at the same time update the
	// Jacobian (in one AD run, if analytic Jacobians are disabled)
	return residual(simTime, simState, nullptr, adJac, threadLocalMem, true, true);
}

template <typename ConvDispOperator>
int LumpedRateModelWithPores<ConvDispOperator>::residualSensFwdAdOnly(const SimulationTime& simTime, const ConstSimulationState& simState, active* const adRes, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidualSens);

	// Evaluate residual for all parameters using AD in vector mode
	return residualImpl<double, active, active, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, adRes, threadLocalMem);
}

template <typename ConvDispOperator>
int LumpedRateModelWithPores<ConvDispOperator>::residualSensFwdCombine(const SimulationTime& simTime, const ConstSimulationState& simState,
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
template <typename ConvDispOperator>
void LumpedRateModelWithPores<ConvDispOperator>::multiplyWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* yS, double alpha, double beta, double* ret)
{
	Indexer idxr(_disc);

	// Handle identity matrix of inlet DOFs
	for (unsigned int i = 0; i < _disc.nComp; ++i)
	{
		ret[i] = alpha * yS[i] + beta * ret[i];
	}

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(std::size_t(0), static_cast<std::size_t>(_disc.nParType + 1), [&](std::size_t idx)
#else
	for (unsigned int idx = 0; idx < _disc.nParType + 1; ++idx)
#endif
	{
		if (cadet_unlikely(idx == 0))
		{
			// Interstitial block
			_convDispOp.jacobian().multiplyVector(yS + idxr.offsetC(), alpha, beta, ret + idxr.offsetC());
			_jacCF.multiplyVector(yS + idxr.offsetJf(), alpha, 1.0, ret + idxr.offsetC());
		}
		else
		{
			// Particle blocks
			const unsigned int type = idx - 1;
			const int localOffset = idxr.offsetCp(ParticleTypeIndex{type});
			_jacP[type].multiplyVector(yS + localOffset, alpha, beta, ret + localOffset);
			_jacPF[type].multiplyVector(yS + idxr.offsetJf(), alpha, 1.0, ret + localOffset);
		}
	} CADET_PARFOR_END;

	// Handle flux equation

	// Set fluxes(ret) = fluxes(yS)
	// This applies the identity matrix in the bottom right corner of the Jaocbian (flux equation)
	for (unsigned int i = idxr.offsetJf(); i < numDofs(); ++i)
		ret[i] = alpha * yS[i] + beta * ret[i];

	double* const retJf = ret + idxr.offsetJf();
	_jacFC.multiplyVector(yS + idxr.offsetC(), alpha, 1.0, retJf);
	for (unsigned int type = 0; type < _disc.nParType; ++type)
	{
		_jacFP[type].multiplyVector(yS + idxr.offsetCp(ParticleTypeIndex{type}), alpha, 1.0, retJf);
	}

	// Map inlet DOFs to the column inlet (first bulk cells)
	_jacInlet.multiplyAdd(yS, ret + idxr.offsetC(), alpha);
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
template <typename ConvDispOperator>
void LumpedRateModelWithPores<ConvDispOperator>::multiplyWithDerivativeJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* sDot, double* ret)
{
	Indexer idxr(_disc);

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(std::size_t(0), static_cast<std::size_t>(_disc.nCol * _disc.nParType + 1), [&](std::size_t idx)
#else
	for (unsigned int idx = 0; idx < _disc.nCol * _disc.nParType + 1; ++idx)
#endif
	{
		if (cadet_unlikely(idx == 0))
		{
			_convDispOp.multiplyWithDerivativeJacobian(simTime, sDot, ret);
		}
		else
		{
			const unsigned int idxParLoop = idx - 1;
			const unsigned int pblk = idxParLoop % _disc.nCol;
			const unsigned int type = idxParLoop / _disc.nCol;

			// Particle
			double const* const localSdot = sDot + idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{pblk});
			double* const localRet = ret + idxr.offsetCp(ParticleTypeIndex{type}, ParticleIndex{pblk});

			unsigned int const* const nBound = _disc.nBound + type * _disc.nComp;
			unsigned int const* const boundOffset = _disc.boundOffset + type * _disc.nComp;

			// Mobile phase
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				// Add derivative with respect to dc_p / dt to Jacobian
				localRet[comp] = localSdot[comp];

				const double invBetaP = (1.0 - static_cast<double>(_parPorosity[type])) / (static_cast<double>(_poreAccessFactor[type * _disc.nComp + comp]) * static_cast<double>(_parPorosity[type]));

				// Add derivative with respect to dq / dt to Jacobian (normal equations)
				for (unsigned int i = 0; i < nBound[comp]; ++i)
				{
					// Index explanation:
					//   nComp -> skip mobile phase
					//   + boundOffset[comp] skip bound states of all previous components
					//   + i go to current bound state
					localRet[comp] += invBetaP * localSdot[_disc.nComp + boundOffset[comp] + i];
				}
			}

			// Solid phase
			double const* const solidSdot = localSdot + _disc.nComp;
			double* const solidRet = localRet + _disc.nComp;
			int const* const qsReaction = _binding[type]->reactionQuasiStationarity();

			for (unsigned int bnd = 0; bnd < _disc.strideBound[type]; ++bnd)
			{
				// Add derivative with respect to dynamic states to Jacobian
				if (qsReaction[bnd])
					solidRet[bnd] = 0.0;
				else
					solidRet[bnd] = solidSdot[bnd];
			}
		}
	} CADET_PARFOR_END;

	// Handle fluxes (all algebraic)
	double* const dFdyDot = ret + idxr.offsetJf();
	std::fill(dFdyDot, dFdyDot + _disc.nCol * _disc.nComp * _disc.nParType, 0.0);

	// Handle inlet DOFs (all algebraic)
	std::fill_n(ret, _disc.nComp, 0.0);
}

template <typename ConvDispOperator>
void LumpedRateModelWithPores<ConvDispOperator>::setExternalFunctions(IExternalFunction** extFuns, unsigned int size)
{
	for (IBindingModel* bm : _binding)
	{
		if (bm)
			bm->setExternalFunctions(extFuns, size);
	}
}

template <typename ConvDispOperator>
unsigned int LumpedRateModelWithPores<ConvDispOperator>::localOutletComponentIndex(unsigned int port) const CADET_NOEXCEPT
{
	// Inlets are duplicated so need to be accounted for
	if (_convDispOp.forwardFlow())
		// Forward Flow: outlet is last cell
		return _disc.nComp + (_disc.nCol - 1) * _disc.nComp;
	else
		// Backward flow: Outlet is first cell
		return _disc.nComp;
}

template <typename ConvDispOperator>
unsigned int LumpedRateModelWithPores<ConvDispOperator>::localInletComponentIndex(unsigned int port) const CADET_NOEXCEPT
{
	return 0;
}

template <typename ConvDispOperator>
unsigned int LumpedRateModelWithPores<ConvDispOperator>::localOutletComponentStride(unsigned int port) const CADET_NOEXCEPT
{
	return 1;
}

template <typename ConvDispOperator>
unsigned int LumpedRateModelWithPores<ConvDispOperator>::localInletComponentStride(unsigned int port) const CADET_NOEXCEPT
{
	return 1;
}

template <typename ConvDispOperator>
void LumpedRateModelWithPores<ConvDispOperator>::expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut)
{
	// @todo Write this function
}

template <typename ConvDispOperator>
bool LumpedRateModelWithPores<ConvDispOperator>::setParameter(const ParameterId& pId, double value)
{
	if (pId.unitOperation == _unitOpIdx)
	{
		// Intercept changes to PAR_TYPE_VOLFRAC when not specified per axial cell (but once globally)
		if (_axiallyConstantParTypeVolFrac && (pId.name == hashString("PAR_TYPE_VOLFRAC")))
		{
			if ((pId.section != SectionIndep) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep))
				return false;
			if (pId.particleType >= _disc.nParType)
				return false;

			for (unsigned int i = 0; i < _disc.nCol; ++i)
				_parTypeVolFrac[i * _disc.nParType + pId.particleType].setValue(value);

			return true;
		}

		for (int parType = 0; parType < _disc.nParType; parType++)
		{
			if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_RADIUS"), _singleParRadius, _parRadius[parType], parType, value, nullptr))
				return true;
			if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_POROSITY"), _singleParPorosity, _parPorosity[parType], parType, value, nullptr))
				return true;

			if (singleTypeMultiplexCompTypeSecParameterValue(pId, hashString("FILM_DIFFUSION"), _filmDiffusionMode, _filmDiffusion, _disc.nComp, parType, value, nullptr, (_filmDiffusion.size() / _disc.nParType) * parType))
				return true;
			if (singleTypeMultiplexCompTypeSecParameterValue(pId, hashString("PORE_ACCESSIBILITY"), _poreAccessFactorMode, _poreAccessFactor, _disc.nComp, parType, value, nullptr, (_poreAccessFactor.size() / _disc.nParType) * parType))
				return true;
		}

		const int mpIc = multiplexInitialConditions(pId, value, false);
		if (mpIc > 0)
			return true;
		else if (mpIc < 0)
			return false;

		if (_convDispOp.setParameter(pId, value))
			return true;

		if (_filmDiffDep)
		{
			if (_filmDiffDep->hasParameter(pId))
			{
				_filmDiffDep->setParameter(pId, value);
				return true;
			}
		}

		if (model::setParameter(pId, value, _reaction.getDynReactionVector("liquid"), false))
			return true;
	}

	return UnitOperationBase::setParameter(pId, value);
}

template <typename ConvDispOperator>
bool LumpedRateModelWithPores<ConvDispOperator>::setParameter(const ParameterId& pId, int value)
{
	if (pId.unitOperation == _unitOpIdx)
	{
		if (model::setParameter(pId, value, _reaction.getDynReactionVector("liquid"), false))
			return true;
	}
	{
		if (model::setParameter(pId, value, _reaction.getDynReactionVector("liquid"), false))
			return true;
	}

	return UnitOperationBase::setParameter(pId, value);
}

template <typename ConvDispOperator>
bool LumpedRateModelWithPores<ConvDispOperator>::setParameter(const ParameterId& pId, bool value)
{
	if (pId.unitOperation == _unitOpIdx)
	{
		if (model::setParameter(pId, value, _reaction.getDynReactionVector("liquid"), false))
			return true;
	}

	return UnitOperationBase::setParameter(pId, value);
}

template <typename ConvDispOperator>
void LumpedRateModelWithPores<ConvDispOperator>::setSensitiveParameterValue(const ParameterId& pId, double value)
{
	if (pId.unitOperation == _unitOpIdx)
	{
		// Intercept changes to PAR_TYPE_VOLFRAC when not specified per axial cell (but once globally)
		if (_axiallyConstantParTypeVolFrac && (pId.name == hashString("PAR_TYPE_VOLFRAC")))
		{
			if ((pId.section != SectionIndep) || (pId.component != CompIndep) || (pId.boundState != BoundStateIndep) || (pId.reaction != ReactionIndep))
				return;
			if (pId.particleType >= _disc.nParType)
				return;

			if (!contains(_sensParams, &_parTypeVolFrac[pId.particleType]))
				return;

			for (unsigned int i = 0; i < _disc.nCol; ++i)
				_parTypeVolFrac[i * _disc.nParType + pId.particleType].setValue(value);

			return;
		}

		for (int parType = 0; parType < _disc.nParType; parType++)
		{
			if (singleTypeMultiplexCompTypeSecParameterValue(pId, hashString("PORE_ACCESSIBILITY"), _poreAccessFactorMode, _poreAccessFactor, _disc.nComp, parType, value, &_sensParams, (_poreAccessFactor.size() / _disc.nParType) * parType))
				return;
			if (singleTypeMultiplexCompTypeSecParameterValue(pId, hashString("FILM_DIFFUSION"), _filmDiffusionMode, _filmDiffusion, _disc.nComp, parType, value, &_sensParams, (_filmDiffusion.size() / _disc.nParType) * parType))
				return;
			if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_RADIUS"), _singleParRadius, _parRadius[parType], parType, value, &_sensParams))
				return;
			if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_POROSITY"), _singleParPorosity, _parPorosity[parType], parType, value, &_sensParams))
				return;
		}

		if (_convDispOp.setSensitiveParameterValue(_sensParams, pId, value))
			return;

		if (_filmDiffDep)
		{
			active* const param = _filmDiffDep->getParameter(pId);
			if (param)
			{
				param->setValue(value);
			}
		}

		if (model::setSensitiveParameterValue(pId, value, _sensParams, _reaction.getDynReactionVector("liquid"), false))
			return;
	}

	UnitOperationBase::setSensitiveParameterValue(pId, value);
}

template <typename ConvDispOperator>
bool LumpedRateModelWithPores<ConvDispOperator>::setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue)
{
	if (pId.unitOperation == _unitOpIdx)
	{
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
			for (unsigned int i = 0; i < _disc.nCol; ++i)
				_parTypeVolFrac[i * _disc.nParType + pId.particleType].setADValue(adDirection, adValue);

			return true;
		}

		for (int parType = 0; parType < _disc.nParType; parType++)
		{
			if (singleTypeMultiplexCompTypeSecParameterAD(pId, hashString("PORE_ACCESSIBILITY"), _poreAccessFactorMode, _poreAccessFactor, _disc.nComp, parType, adDirection, adValue, _sensParams, (_poreAccessFactor.size() / _disc.nParType) * parType))
			{
				LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
				return true;
			}

			if (singleTypeMultiplexCompTypeSecParameterAD(pId, hashString("FILM_DIFFUSION"), _filmDiffusionMode, _filmDiffusion, _disc.nComp, parType, adDirection, adValue, _sensParams, (_filmDiffusion.size() / _disc.nParType) * parType))
			{
				LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
				return true;
			}

			if (singleTypeMultiplexTypeParameterAD(pId, hashString("PAR_RADIUS"), _singleParRadius, _parRadius[parType], parType, adDirection, adValue, _sensParams))
			{
				LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
				return true;
			}

			if (singleTypeMultiplexTypeParameterAD(pId, hashString("PAR_POROSITY"), _singleParPorosity, _parPorosity[parType], parType, adDirection, adValue, _sensParams))
			{
				LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
				return true;
			}
		}

		const int mpIc = multiplexInitialConditions(pId, adDirection, adValue);
		if (mpIc > 0)
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}
		else if (mpIc < 0)
			return false;

		if (_convDispOp.setSensitiveParameter(_sensParams, pId, adDirection, adValue))
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		if (_filmDiffDep)
		{
			active* const param = _filmDiffDep->getParameter(pId);
			if (param)
			{
				LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
				param->setADValue(adDirection, adValue);
				return true;
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


template <typename ConvDispOperator>
int LumpedRateModelWithPores<ConvDispOperator>::Exporter::writeMobilePhase(double* buffer) const
{
	const int blockSize = _disc.nComp * _disc.nCol;
	std::copy_n(_idx.c(_data), blockSize, buffer);
	return blockSize;
}

template <typename ConvDispOperator>
int LumpedRateModelWithPores<ConvDispOperator>::Exporter::writeSolidPhase(double* buffer) const
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

template <typename ConvDispOperator>
int LumpedRateModelWithPores<ConvDispOperator>::Exporter::writeParticleMobilePhase(double* buffer) const
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

template <typename ConvDispOperator>
int LumpedRateModelWithPores<ConvDispOperator>::Exporter::writeSolidPhase(unsigned int parType, double* buffer) const
{
	cadet_assert(parType < _disc.nParType);

	const unsigned int stride = _disc.nComp + _disc.strideBound[parType];
	double const* ptr = _data + _idx.offsetCp(ParticleTypeIndex{parType}) + _idx.strideParLiquid();
	for (unsigned int i = 0; i < _disc.nCol; ++i)
	{
		std::copy_n(ptr, _disc.strideBound[parType], buffer);
		buffer += _disc.strideBound[parType];
		ptr += stride;
	}
	return _disc.nCol * _disc.strideBound[parType];
}

template <typename ConvDispOperator>
int LumpedRateModelWithPores<ConvDispOperator>::Exporter::writeParticleMobilePhase(unsigned int parType, double* buffer) const
{
	cadet_assert(parType < _disc.nParType);

	const unsigned int stride = _disc.nComp + _disc.strideBound[parType];
	double const* ptr = _data + _idx.offsetCp(ParticleTypeIndex{parType});
	for (unsigned int i = 0; i < _disc.nCol; ++i)
	{
		std::copy_n(ptr, _disc.nComp, buffer);
		buffer += _disc.nComp;
		ptr += stride;
	}
	return _disc.nCol * _disc.nComp;
}

template <typename ConvDispOperator>
int LumpedRateModelWithPores<ConvDispOperator>::Exporter::writeParticleFlux(double* buffer) const
{
	const int blockSize = numParticleFluxDofs();
	std::copy_n(_idx.jf(_data), blockSize, buffer);
	return blockSize;
}

template <typename ConvDispOperator>
int LumpedRateModelWithPores<ConvDispOperator>::Exporter::writeParticleFlux(unsigned int parType, double* buffer) const
{
	const unsigned int blockSize = _disc.nComp * _disc.nCol;
	std::copy_n(_idx.jf(_data) + blockSize * parType, blockSize, buffer);
	return blockSize;
}

template <typename ConvDispOperator>
int LumpedRateModelWithPores<ConvDispOperator>::Exporter::writeInlet(unsigned int port, double* buffer) const
{
	cadet_assert(port == 0);
	std::copy_n(_data, _disc.nComp, buffer);
	return _disc.nComp;
}

template <typename ConvDispOperator>
int LumpedRateModelWithPores<ConvDispOperator>::Exporter::writeInlet(double* buffer) const
{
	std::copy_n(_data, _disc.nComp, buffer);
	return _disc.nComp;
}

template <typename ConvDispOperator>
int LumpedRateModelWithPores<ConvDispOperator>::Exporter::writeOutlet(unsigned int port, double* buffer) const
{
	cadet_assert(port == 0);

	if (_model._convDispOp.forwardFlow())
		std::copy_n(&_idx.c(_data, _disc.nCol - 1, 0), _disc.nComp, buffer);
	else
		std::copy_n(&_idx.c(_data, 0, 0), _disc.nComp, buffer);

	return _disc.nComp;
}

template <typename ConvDispOperator>
int LumpedRateModelWithPores<ConvDispOperator>::Exporter::writeOutlet(double* buffer) const
{
	if (_model._convDispOp.forwardFlow())
		std::copy_n(&_idx.c(_data, _disc.nCol - 1, 0), _disc.nComp, buffer);
	else
		std::copy_n(&_idx.c(_data, 0, 0), _disc.nComp, buffer);

	return _disc.nComp;
}

}  // namespace model

}  // namespace cadet

#include "model/LumpedRateModelWithPores-InitialConditions.cpp"
#include "model/LumpedRateModelWithPores-LinearSolver.cpp"

namespace cadet
{

namespace model
{

// Template instantiations
template class LumpedRateModelWithPores<parts::AxialConvectionDispersionOperator>;
template class LumpedRateModelWithPores<parts::RadialConvectionDispersionOperator>;

IUnitOperation* createAxialFVLRMP(UnitOpIdx uoId)
{
	typedef LumpedRateModelWithPores<parts::AxialConvectionDispersionOperator> AxialLRMP;

	return new AxialLRMP(uoId);
}

IUnitOperation* createRadialFVLRMP(UnitOpIdx uoId)
{
	typedef LumpedRateModelWithPores<parts::RadialConvectionDispersionOperator> RadialLRMP;

	return new RadialLRMP(uoId);
}

}  // namespace model

}  // namespace cadet
