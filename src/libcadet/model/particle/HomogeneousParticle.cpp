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

#include "model/particle/HomogeneousParticle.hpp"
#include "model/parts/ParticleDiffusionOperatorDG.hpp"

#include "cadet/Exceptions.hpp"
#include "BindingModelFactory.hpp"
#include "ReactionModelFactory.hpp"
#include "ParamReaderHelper.hpp"
#include "ParamReaderScopes.hpp"
#include "AdUtils.hpp"
#include "SimulationTypes.hpp"
#include "model/ParameterDependence.hpp"
#include "model/parts/BindingCellKernel.hpp"
#include "SensParamUtil.hpp"
#include "ConfigurationHelper.hpp"
#include "model/BindingModel.hpp"
#include "model/ReactionModel.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

namespace cadet
{

namespace model
{

	/**
	 * @brief Creates a HomogeneousParticle
	 */
	HomogeneousParticle::HomogeneousParticle() : _boundOffset(nullptr)
	{
	}

	HomogeneousParticle::~HomogeneousParticle() CADET_NOEXCEPT
	{
		delete[] _boundOffset;
	}

	bool HomogeneousParticle::configureModelDiscretization_old(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp)
	{
		_strideBulkComp = strideBulkComp;

		_parTypeIdx = parTypeIdx;
		_nComp = nComp;

		// Read particle geometry and default to "SPHERICAL"
		_parGeomSurfToVol = _SurfVolRatioSphere;
		if (paramProvider.exists("PAR_GEOM"))
		{
			std::vector<std::string> pg = paramProvider.getStringArray("PAR_GEOM");
			if ((pg.size() == 1) && (nParType > 1))
			{
				// Multiplex using first value
				pg.resize(nParType, pg[0]);
			}
			else if (pg.size() < nParType)
				throw InvalidParameterException("Field PAR_GEOM contains too few elements (" + std::to_string(nParType) + " required)");

			if (pg[_parTypeIdx] == "SPHERE")
				_parGeomSurfToVol = _SurfVolRatioSphere;
			else if (pg[_parTypeIdx] == "CYLINDER")
				_parGeomSurfToVol = _SurfVolRatioCylinder;
			else if (pg[_parTypeIdx] == "SLAB")
				_parGeomSurfToVol = _SurfVolRatioSlab;
			else
				throw InvalidParameterException("Unknown particle geometry type \"" + pg[_parTypeIdx] + "\" at index " + std::to_string(_parTypeIdx) + " of field PAR_GEOM");
		}

		std::vector<int> nBound;
		const bool newNBoundInterface = paramProvider.exists("NBOUND");

		paramProvider.pushScope("discretization");

		if (!newNBoundInterface && paramProvider.exists("NBOUND")) // done here and in this order for backwards compatibility
		{
			nBound = paramProvider.getIntArray("NBOUND");
			paramProvider.popScope();
		}
		else
		{
			paramProvider.popScope();
			nBound = paramProvider.getIntArray("NBOUND");
		}
		if (nBound.size() < _nComp)
			throw InvalidParameterException("Field NBOUND contains too few elements (NCOMP = " + std::to_string(_nComp) + " required)");

		std::vector<int> stridesParTypeBound(nParType + 1);
		std::vector<int> nBoundBeforeType(nParType);
		if (!_nBound)
			_nBound = std::make_shared<unsigned int[]>(_nComp);

		if (nBound.size() < _nComp * nParType)
		{
			std::copy_n(nBound.begin(), _nComp, _nBound.get());

			stridesParTypeBound[0] = std::accumulate(nBound.begin(), nBound.begin() + _nComp, 0);
			nBoundBeforeType[0] = 0;

			for (int type = 1; type < nParType; type++)
			{
				stridesParTypeBound[type] = stridesParTypeBound[0];
				nBoundBeforeType[type] += nBoundBeforeType[type - 1] + _nBound[type - 1];
			}
		}
		else
		{
			std::copy_n(nBound.begin() + _parTypeIdx * _nComp, _nComp, _nBound.get());

			stridesParTypeBound[0] = std::accumulate(nBound.begin(), nBound.begin() + _nComp, 0);
			nBoundBeforeType[0] = 0;

			for (int type = 1; type < nParType; type++)
			{
				stridesParTypeBound[type] = std::accumulate(nBound.begin() + type * _nComp, nBound.begin() + (type + 1) * _nComp, 0);
				nBoundBeforeType[type] += nBoundBeforeType[type - 1] + nBound[type - 1];
			}
		}

		// Precompute offsets and total number of bound states (DOFs in solid phase)
		if (!_boundOffset)
			_boundOffset = new unsigned int[_nComp];

		_boundOffset[0] = 0.0;
		_strideBound = std::accumulate(_nBound.get(), _nBound.get() + _nComp, 0u);

		for (unsigned int i = 1; i < _nComp; ++i)
		{
			_boundOffset[i] = _boundOffset[i - 1] + _nBound[i - 1];
		}

		// ==== Construct and configure binding model
		_binding = nullptr;
		std::vector<std::string> bindModelNames = { "NONE" };
		if (paramProvider.exists("ADSORPTION_MODEL"))
			bindModelNames = paramProvider.getStringArray("ADSORPTION_MODEL");

		if (paramProvider.exists("ADSORPTION_MODEL_MULTIPLEX"))
			_bindingParDep = (paramProvider.getInt("ADSORPTION_MODEL_MULTIPLEX") == 1);
		else
			{
			// Infer multiplex mode
			_bindingParDep = (bindModelNames.size() == 1);
			}

		if (!_bindingParDep && (bindModelNames.size() < nParType))
			throw InvalidParameterException("Field ADSORPTION_MODEL contains too few elements (" + std::to_string(nParType) + " required)");
		else if (_bindingParDep && (bindModelNames.size() != 1))
			throw InvalidParameterException("Field ADSORPTION_MODEL requires (only) 1 element");

		bool bindingConfSuccess = true;

		_binding = helper.createBindingModel(bindModelNames[_bindingParDep ? 0 : _parTypeIdx]);
		if (!_binding)
			throw InvalidParameterException("Unknown binding model " + bindModelNames[_bindingParDep ? 0 : _parTypeIdx]);

		MultiplexedScopeSelector scopeGuard(paramProvider, "adsorption", _bindingParDep, _parTypeIdx, nParType == 1, _binding->usesParamProviderInDiscretizationConfig());
		bindingConfSuccess = _binding->configureModelDiscretization(paramProvider, _nComp, _nBound.get(), _boundOffset);

		// ==== Construct and configure dynamic reaction model
		bool reactionConfSuccess = true;

		_dynReaction = nullptr;

		if (paramProvider.exists("REACTION_MODEL_PARTICLES"))
		{
			const std::vector<std::string> dynReactModelNames = paramProvider.getStringArray("REACTION_MODEL_PARTICLES");

			if (paramProvider.exists("REACTION_MODEL_PARTICLES_MULTIPLEX"))
				_reactionParDep = (paramProvider.getInt("REACTION_MODEL_PARTICLES_MULTIPLEX") == 1);
			else
			{
				// Infer multiplex mode
				_reactionParDep = (dynReactModelNames.size() == 1);
			}

			if (!_reactionParDep && (dynReactModelNames.size() < nParType))
				throw InvalidParameterException("Field REACTION_MODEL_PARTICLES contains too few elements (" + std::to_string(nParType) + " required)");
			else if (_reactionParDep && (dynReactModelNames.size() != 1))
				throw InvalidParameterException("Field REACTION_MODEL_PARTICLES requires (only) 1 element");

			_dynReaction = helper.createDynamicReactionModel(dynReactModelNames[_reactionParDep ? 0 : _parTypeIdx]);

			if (!_dynReaction)
				throw InvalidParameterException("Unknown dynamic reaction model " + dynReactModelNames[_reactionParDep ? 0 : _parTypeIdx]);

			MultiplexedScopeSelector scopeGuard(paramProvider, "reaction_particle", _reactionParDep, _parTypeIdx, nParType == 1, _dynReaction->usesParamProviderInDiscretizationConfig());
			reactionConfSuccess = _dynReaction->configureModelDiscretization(paramProvider, _nComp, _nBound.get(), _boundOffset) && reactionConfSuccess;
		}

		return bindingConfSuccess && reactionConfSuccess;
	}

	bool HomogeneousParticle::configure_old(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound)
	{
		// Read geometry parameters
		std::vector<double> parRadii(nParType);
		_parRadiusParTypeDep = readScalarParameterOrArray(parRadii, paramProvider, "PAR_RADIUS", nParType);

		if (_parRadiusParTypeDep)
		{
			_parRadius = parRadii[0];
			if (_parTypeIdx == 0)
				parameters[makeParamId(hashStringRuntime("PAR_RADIUS"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parRadius;
		}
		else
		{
			_parRadius = parRadii[_parTypeIdx];
			parameters[makeParamId(hashStringRuntime("PAR_RADIUS"), unitOpIdx, CompIndep, _parTypeIdx, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parRadius;
		}

		std::vector<double> parPorosities(nParType);
		_parPorosityParTypeDep = readScalarParameterOrArray(parPorosities, paramProvider, "PAR_POROSITY", nParType);
		if (_parPorosityParTypeDep)
		{
			_parPorosity = parPorosities[0];
			if (_parTypeIdx == 0)
				parameters[makeParamId(hashStringRuntime("PAR_POROSITY"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parPorosity;
		}
		else
		{
			_parPorosity = parPorosities[_parTypeIdx];
			parameters[makeParamId(hashStringRuntime("PAR_POROSITY"), unitOpIdx, CompIndep, _parTypeIdx, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parPorosity;
		}

		if (nParType != parRadii.size())
			throw InvalidParameterException("Number of elements in field PAR_RADIUS does not match number of particle types");
		if (nParType != parPorosities.size())
			throw InvalidParameterException("Number of elements in field PAR_POROSITY does not match number of particle types");

		if(_parRadius <= 0.0)
			throw InvalidParameterException("Particle radius is not bigger than zero for particle type " + std::to_string(_parTypeIdx));
		if (_parPorosity <= 0.0 || _parPorosity > 1.0)
			throw InvalidParameterException("Particle porosity is not within (0, 1] for particle type " + std::to_string(_parTypeIdx));

		// Read and register film diffusion, poreAccesFactor
		_filmDiffusionMode = readAndRegisterSingleTypeMultiplexCompTypeSecParam(paramProvider, parameters, _filmDiffusion, "FILM_DIFFUSION", nParType, _nComp, _parTypeIdx, unitOpIdx);

		if (paramProvider.exists("PORE_ACCESSIBILITY"))
			_poreAccessFactorMode = readAndRegisterSingleTypeMultiplexCompTypeSecParam(paramProvider, parameters, _poreAccessFactor, "PORE_ACCESSIBILITY", nParType, _nComp, _parTypeIdx, unitOpIdx);
		else
		{
			_poreAccessFactorMode = MultiplexMode::ComponentType;
			_poreAccessFactor = std::vector<cadet::active>(_nComp, 1.0);
		}
		if (_poreAccessFactorMode == MultiplexMode::ComponentSectionType || _poreAccessFactorMode == MultiplexMode::ComponentSection)
		{
			throw InvalidParameterException("Section dependence not supported for PORE_ACCESSIBILITY");
		}

		// Reconfigure binding model
		bool bindingConfSuccess = true;
		if (_binding)
		{
			if (_binding->requiresConfiguration())
			{
				if (_bindingParDep)
				{
					MultiplexedScopeSelector scopeGuard(paramProvider, "adsorption", true);
					bindingConfSuccess = _binding->configure(paramProvider, unitOpIdx, ParTypeIndep);
				}
				else
				{
					MultiplexedScopeSelector scopeGuard(paramProvider, "adsorption", _parTypeIdx, nParType == 1, true);
					bindingConfSuccess = _binding->configure(paramProvider, unitOpIdx, _parTypeIdx);
				}
			}
		}

		// Reconfigure reaction model
		bool dynReactionConfSuccess = true;
		if (_dynReaction && _dynReaction->requiresConfiguration())
		{
			if (_reactionParDep)
			{
				MultiplexedScopeSelector scopeGuard(paramProvider, "reaction_particle", true);
				dynReactionConfSuccess = _dynReaction->configure(paramProvider, unitOpIdx, ParTypeIndep) && dynReactionConfSuccess;
			}
			else
			{
				MultiplexedScopeSelector scopeGuard(paramProvider, "reaction_particle", _parTypeIdx, nParType == 1, true);
				dynReactionConfSuccess = _dynReaction->configure(paramProvider, unitOpIdx, _parTypeIdx) && dynReactionConfSuccess;
			}
		}

		return bindingConfSuccess && dynReactionConfSuccess;
	}

	bool HomogeneousParticle::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp)
	{
		_strideBulkComp = strideBulkComp;

		_parTypeIdx = parTypeIdx;
		_nComp = nComp;

		std::ostringstream parTypeIdxString;
		parTypeIdxString << std::setfill('0') << std::setw(3) << std::setprecision(0) << _parTypeIdx;
		paramProvider.pushScope("particle_type_" + parTypeIdxString.str());

		// Read particle geometry and default to Sphere
		_parGeomSurfToVol = _SurfVolRatioSphere;
		if (paramProvider.exists("PAR_GEOM"))
		{
			std::vector<std::string> pg = paramProvider.getStringArray("PAR_GEOM");
			if (pg.size() > 1)
				throw InvalidParameterException("Only one geometry must be specified, multiple are given for particle type " + std::to_string(_parTypeIdx));
			else

				if (pg[0] == "SPHERE")
					_parGeomSurfToVol = _SurfVolRatioSphere;
				else if (pg[0] == "CYLINDER")
					_parGeomSurfToVol = _SurfVolRatioCylinder;
				else if (pg[0] == "SLAB")
					_parGeomSurfToVol = _SurfVolRatioSlab;
				else
					throw InvalidParameterException("Unknown particle geometry type \"" + pg[0] + "\" for particle type " + std::to_string(_parTypeIdx));
		}

		std::vector<int> nBound = paramProvider.getIntArray("NBOUND");
		if (nBound.size() != _nComp)
			throw InvalidParameterException("Field NBOUND does not contain NCOMP = " + std::to_string(_nComp) + " entries for particle type " + std::to_string(_parTypeIdx));

		if (!_nBound)
			_nBound = std::make_shared<unsigned int[]>(_nComp);
		std::copy_n(nBound.begin(), _nComp, _nBound.get());

		// Precompute offsets and total number of bound states (DOFs in solid phase)
		if (!_boundOffset)
			_boundOffset = new unsigned int[_nComp];

		_boundOffset[0] = 0.0;
		_strideBound = std::accumulate(_nBound.get(), _nBound.get() + _nComp, 0u);

		for (unsigned int i = 1; i < _nComp; ++i)
		{
			_boundOffset[i] = _boundOffset[i - 1] + _nBound[i - 1];
		}

		// ==== Construct and configure binding model

		_binding = nullptr;
		std::vector<std::string> bindModelNames = { "NONE" };
		bool bindingConfSuccess = true;

		if (paramProvider.exists("ADSORPTION_MODEL"))
		{
			bindModelNames = paramProvider.getStringArray("ADSORPTION_MODEL");

			if (_bindingParDep && (bindModelNames.size() != 1))
				throw InvalidParameterException("Field ADSORPTION_MODEL requires (only) 1 element");

			if (paramProvider.exists("adsorption"))
			{
				paramProvider.pushScope("adsorption");
				_bindingParDep = paramProvider.exists("BINDING_PARTYPE_DEPENDENT") ? paramProvider.getInt("BINDING_PARTYPE_DEPENDENT") : true;
				paramProvider.popScope();
			}
		}

		_binding = helper.createBindingModel(bindModelNames[0]);
		if (!_binding)
			throw InvalidParameterException("Unknown binding model " + bindModelNames[0]);

		_bindingParDep = true;

		if (_binding->usesParamProviderInDiscretizationConfig())
		{
			paramProvider.pushScope("adsorption");

			if (paramProvider.exists("BINDING_PARTYPE_DEPENDENT"))
				_bindingParDep = paramProvider.getBool("BINDING_PARTYPE_DEPENDENT");

			paramProvider.popScope();
		}
		else if (bindModelNames[0] == "NONE")
			_nBound = std::make_shared<unsigned int[]>(_nComp, 0);
		else
			throw InvalidParameterException("Binding model " + bindModelNames[0] + " was specified, but group \"adsorption\" is missing for particle type " + std::to_string(_parTypeIdx));

		if (_binding->usesParamProviderInDiscretizationConfig())
			paramProvider.pushScope("adsorption");

		bindingConfSuccess = _binding->configureModelDiscretization(paramProvider, _nComp, _nBound.get(), _boundOffset);

		if (_binding->usesParamProviderInDiscretizationConfig())
			paramProvider.popScope();

		// ==== Construct and configure dynamic reaction model

		_dynReaction = nullptr;
		_reactionParDep = true;
		bool reactionConfSuccess = true;

		if (paramProvider.exists("REACTION_MODEL"))
		{
			const std::vector<std::string> dynReactModelNames = paramProvider.getStringArray("REACTION_MODEL");

			if (dynReactModelNames.size() != 1)
				throw InvalidParameterException("Field REACTION_MODEL_PARTICLES requires (only) 1 element");

			_dynReaction = helper.createDynamicReactionModel(dynReactModelNames[0]);

			if (!_dynReaction)
				throw InvalidParameterException("Unknown dynamic reaction model " + dynReactModelNames[0]);

			if (paramProvider.exists("reaction"))
			{
				paramProvider.pushScope("reaction");
				_reactionParDep = paramProvider.exists("REACTION_PARTYPE_DEPENDENT") ? paramProvider.getInt("REACTION_PARTYPE_DEPENDENT") : true;
				paramProvider.popScope();
			}

			reactionConfSuccess = _dynReaction->configureModelDiscretization(paramProvider, _nComp, _nBound.get(), _boundOffset) && reactionConfSuccess;
		}

		paramProvider.popScope(); // particle_type_{:03}

		return bindingConfSuccess && reactionConfSuccess;
	}

	bool HomogeneousParticle::configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound)
	{
		std::ostringstream parTypeIdxString;
		parTypeIdxString << std::setfill('0') << std::setw(3) << std::setprecision(0) << _parTypeIdx;
		paramProvider.pushScope("particle_type_" + parTypeIdxString.str());

		// Read geometry parameters

		_parRadius = paramProvider.getDouble("PAR_RADIUS");
		
		_parRadiusParTypeDep =  true;
		if (paramProvider.exists("PAR_RADIUS_PARTYPE_DEPENDENT"))
			_parRadiusParTypeDep = paramProvider.getBool("PAR_RADIUS_PARTYPE_DEPENDENT");

		if (!_parRadiusParTypeDep && _parTypeIdx == 0)
			parameters[makeParamId(hashStringRuntime("PAR_RADIUS"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parRadius;
		else if (_parRadiusParTypeDep)
			parameters[makeParamId(hashStringRuntime("PAR_RADIUS"), unitOpIdx, CompIndep, _parTypeIdx, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parRadius;

		if (_parRadius <= 0.0)
			throw InvalidParameterException("Particle radius is not bigger than zero for particle type " + std::to_string(_parTypeIdx));

		_parPorosity = paramProvider.getDouble("PAR_POROSITY");

		_parPorosityParTypeDep = true;
		if (paramProvider.exists("PAR_POROSITY_PARTYPE_DEPENDENT"))
			_parPorosityParTypeDep = paramProvider.getBool("PAR_POROSITY_PARTYPE_DEPENDENT");

		if (!_parPorosityParTypeDep && _parTypeIdx == 0)
			parameters[makeParamId(hashStringRuntime("PAR_POROSITY"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parPorosity;
		else if (_parPorosityParTypeDep)
			parameters[makeParamId(hashStringRuntime("PAR_POROSITY"), unitOpIdx, CompIndep, _parTypeIdx, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parPorosity;

		if (_parPorosity <= 0.0 || _parPorosity > 1.0)
			throw InvalidParameterException("Particle porosity is not within (0, 1] for particle type " + std::to_string(_parTypeIdx));

		bool filmDiffParTypeDep = paramProvider.exists("FILM_DIFFUSION_PARTYPE_DEPENDENT") ? paramProvider.getBool("FILM_DIFFUSION_PARTYPE_DEPENDENT") : true;
		_filmDiffusionMode = newIF_readAndRegisterMultiplexCompSecParam(paramProvider, parameters, _filmDiffusion, "FILM_DIFFUSION", _nComp, _parTypeIdx, filmDiffParTypeDep, unitOpIdx);

		if (paramProvider.exists("PORE_ACCESSIBILITY"))
		{
			bool poreAccessParTypeDep = paramProvider.exists("PORE_ACCESSIBILITY_PARTYPE_DEPENDENT") ? paramProvider.getBool("PORE_ACCESSIBILITY_PARTYPE_DEPENDENT") : true;
			_poreAccessFactorMode = newIF_readAndRegisterMultiplexCompSecParam(paramProvider, parameters, _poreAccessFactor, "PORE_ACCESSIBILITY", _nComp, _parTypeIdx, poreAccessParTypeDep, unitOpIdx);
		}
		else
		{
			_poreAccessFactorMode = MultiplexMode::ComponentType;
			_poreAccessFactor = std::vector<cadet::active>(_nComp, 1.0);
		}

		//// Done in the unit operation: Register initial conditions parameters
		//registerParam1DArray(parameters, _initC, [=](bool multi, unsigned int comp) { return makeParamId(hashString("INIT_C"), unitOpIdx, comp, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep); });

		//if (__bindingParDep)
		//{
		//	for (unsigned int c = 0; c < nComp; ++c)
		//		parameters[makeParamId(hashString("INIT_CP"), unitOpIdx, c, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_initCp[c];
		//}
		//else
		//	registerParam2DArray(parameters, _initCp, [=](bool multi, unsigned int type, unsigned int comp) { return makeParamId(hashString("INIT_CP"), unitOpIdx, comp, type, BoundStateIndep, ReactionIndep, SectionIndep); }, nComp);


		//if (!_binding.empty())
		//{
		//	const unsigned int maxBoundStates = *std::max_element(_strideBound, _strideBound + _nParType);
		//	std::vector<ParameterId> initParams(maxBoundStates);

		//	if (__bindingParDep)
		//	{
		//		_binding[0]->fillBoundPhaseInitialParameters(initParams.data(), unitOpIdx, ParTypeIndep);

		//		active* const iq = _initQ.data() + _nBoundBeforeType[0];
		//		for (unsigned int i = 0; i < _strideBound[0]; ++i)
		//			parameters[initParams[i]] = iq + i;
		//	}
		//	else
		//	{
		//		for (unsigned int type = 0; type < _nParType; ++type)
		//		{
		//			_binding[type]->fillBoundPhaseInitialParameters(initParams.data(), unitOpIdx, type);

		//			active* const iq = _initQ.data() + _nBoundBeforeType[type];
		//			for (unsigned int i = 0; i < _strideBound[type]; ++i)
		//				parameters[initParams[i]] = iq + i;
		//		}
		//	}
		//}

		// Reconfigure binding model
		bool bindingConfSuccess = true;
		if (_binding)
		{
			if (_binding->requiresConfiguration())
			{
				if (_bindingParDep)
				{
					MultiplexedScopeSelector scopeGuard(paramProvider, "adsorption", true);
					bindingConfSuccess = _binding->configure(paramProvider, unitOpIdx, ParTypeIndep);
				}
				else
				{
					MultiplexedScopeSelector scopeGuard(paramProvider, "adsorption", _parTypeIdx, nParType == 1, true);
					bindingConfSuccess = _binding->configure(paramProvider, unitOpIdx, _parTypeIdx);
				}
			}
		}

		// Reconfigure reaction model
		bool dynReactionConfSuccess = true;
		if (_dynReaction && _dynReaction->requiresConfiguration())
		{
			if (_reactionParDep)
			{
				MultiplexedScopeSelector scopeGuard(paramProvider, "reaction_particle", true);
				dynReactionConfSuccess = _dynReaction->configure(paramProvider, unitOpIdx, ParTypeIndep) && dynReactionConfSuccess;
			}
			else
			{
				MultiplexedScopeSelector scopeGuard(paramProvider, "reaction_particle", _parTypeIdx, nParType == 1, true);
				dynReactionConfSuccess = _dynReaction->configure(paramProvider, unitOpIdx, _parTypeIdx) && dynReactionConfSuccess;
			}
		}

		paramProvider.popScope(); // particle_type_{:03}

		return bindingConfSuccess && dynReactionConfSuccess;
	}

	bool HomogeneousParticle::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx)
	{
		return true;
	}

	int HomogeneousParticle::writeParticleCoordinates(double* coords) const
	{
		coords[0] = static_cast<double>(_parRadius) / 0.5;
		return 1;
	}

	parts::cell::CellParameters HomogeneousParticle::makeCellResidualParams(int const* qsReaction, unsigned int const* nBound) const
	{
		return parts::cell::CellParameters
		{
			_nComp,
			nBound,
			_boundOffset,
			_strideBound,
			qsReaction,
			getPorosity(),
			getPoreAccessFactor(),
			_binding,
			(_dynReaction && (_dynReaction->numReactionsCombined() > 0)) ? _dynReaction : nullptr
		};
	}

	int HomogeneousParticle::residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, double* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity)
	{
		if (resPar)
		{
			if (jacIt.data())
				return residualImpl<double, double, double, true, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
			else
				return residualImpl<double, double, double, false, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
		}
		else if (jacIt.data())
			return residualImpl<double, double, double, true, false>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
		else
			return -1;
	}
	int HomogeneousParticle::residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity)
	{
		if (jacIt.data())
			return residualImpl<double, active, active, true, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
		else
			return residualImpl<double, active, active, false, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
	}
	int HomogeneousParticle::residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity)
	{
		if (jacIt.data())
			return residualImpl<active, active, double, true, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
		else
			return residualImpl<active, active, double, false, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
	}
	int HomogeneousParticle::residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity)
	{
		if (jacIt.data())
			return residualImpl<active, active, active, true, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
		else
			return residualImpl<active, active, active, false, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
	}

	template <typename StateType, typename ResidualType, typename ParamType, bool wantNonLinJac, bool wantRes>
	int HomogeneousParticle::residualImpl(double t, unsigned int secIdx, StateType const* yPar, StateType const* yBulk, double const* yDotPar, ResidualType* resPar, ResidualType* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc)
	{
		int const* const qsBinding = _binding->reactionQuasiStationarity();
		const parts::cell::CellParameters cellResParams = makeCellResidualParams(qsBinding, _nBound.get());

		linalg::BandedEigenSparseRowIterator jacBase = jacIt;

		// Handle time derivatives, binding, dynamic reactions

		// r (particle) coordinate of current node (particle radius normed to 1) - needed in externally dependent adsorption kinetic
		packing.colPos.particle = relativeCoordinate(0);

		if (wantRes)
			parts::cell::residualKernel<StateType, ResidualType, ParamType, parts::cell::CellParameters, linalg::BandedEigenSparseRowIterator, wantNonLinJac, true>(
				t, secIdx, packing.colPos, yPar, yDotPar ? yDotPar : nullptr, resPar ? resPar : nullptr, jacIt, cellResParams, tlmAlloc
			);
		else
			parts::cell::residualKernel<StateType, ResidualType, ParamType, parts::cell::CellParameters, linalg::BandedEigenSparseRowIterator, wantNonLinJac, false, false>(
				t, secIdx, packing.colPos, yPar, yDotPar ? yDotPar : nullptr, resPar ? resPar : nullptr, jacIt, cellResParams, tlmAlloc
			);

		const ParamType jacPF_val = -static_cast<ParamType>(surfaceToVolumeRatio()) / static_cast<ParamType>(_parPorosity);
		const ParamType invBetaC = 1.0 / static_cast<ParamType>(packing.colPorosity) - 1.0;
		const ParamType jacCF_val = invBetaC * static_cast<ParamType>(surfaceToVolumeRatio());

		// Film diffusion flux
		if (wantRes)
		{
			const active* const filmDiff = getSectionDependentSlice(_filmDiffusion, _nComp, secIdx);

			for (unsigned int comp = 0; comp < _nComp; ++comp)
			{
				// flux into particle
				resPar[comp] += jacPF_val / static_cast<ParamType>(_poreAccessFactor[comp]) * static_cast<ParamType>(filmDiff[comp]) * (yBulk[comp * _strideBulkComp] - yPar[comp]);
				// flux into bulk
				resBulk[comp] += jacCF_val * static_cast<ParamType>(filmDiff[comp]) * static_cast<ParamType>(packing.parTypeVolFrac) * (yBulk[comp] - yPar[comp]);
			}
		}

		return true;
	}

	void HomogeneousParticle::setParJacPattern(std::vector<Eigen::Triplet<double>>& tripletList, const unsigned int offsetPar, const unsigned int offsetBulk, unsigned int colNode, unsigned int secIdx) const
	{
		// binding Jacobian pattern
		// add dense nComp * nBound blocks, since all solid and liquid entries can be coupled through binding.
		for (unsigned int parState = 0; parState < _nComp + _strideBound; parState++) {
			for (unsigned int toParState = 0; toParState < _nComp + _strideBound; toParState++) {
				tripletList.push_back(Eigen::Triplet<double>(offsetPar + parState, offsetPar + toParState, 0.0));
			}
		}

		// flux Jacobian pattern

		for (unsigned int comp = 0; comp < _nComp; comp++)
		{
			tripletList.push_back(Eigen::Triplet<double>(offsetPar + comp, offsetBulk + comp, 0.0));
			tripletList.push_back(Eigen::Triplet<double>(offsetBulk + comp, offsetPar + comp, 0.0));
		}
	}


	unsigned int HomogeneousParticle::jacobianNNZperParticle() const
	{
		return (_nComp + _strideBound) * (_nComp + _strideBound) + _nComp * 4; // reaction, binding patter + film diffusion pattern for one particle
	}

	int HomogeneousParticle::calcParticleDiffJacobian(const int secIdx, const int colNode, const int offsetLocalCp, Eigen::SparseMatrix<double, RowMajor>& globalJac)
	{
		return 1;
	}
	
	int HomogeneousParticle::calcFilmDiffJacobian(unsigned int secIdx, const int offsetCp, const int offsetC, const int nBulkPoints, const int nParType, const double colPorosity, const active* const parTypeVolFrac, Eigen::SparseMatrix<double, RowMajor>& globalJac, bool crossDepsOnly)
	{
		const double invBetaC = 1.0 / static_cast<double>(colPorosity) - 1.0;

		const double jacCF_val = invBetaC * static_cast<double>(surfaceToVolumeRatio());
		const double jacPF_val = -static_cast<double>(surfaceToVolumeRatio()) / static_cast<double>(_parPorosity);

		// Ordering of diffusion:
		// sec0type0comp0, sec0type0comp1, sec0type0comp2, sec0type1comp0, sec0type1comp1, sec0type1comp2,
		// sec1type0comp0, sec1type0comp1, sec1type0comp2, sec1type1comp0, sec1type1comp1, sec1type1comp2, ...
		active const* const filmDiff = getSectionDependentSlice(_filmDiffusion, _nComp, secIdx);
		active const* const poreAccFactor = _poreAccessFactor.data();

		linalg::BandedEigenSparseRowIterator jacC(globalJac, offsetC);
		linalg::BandedEigenSparseRowIterator jacP(globalJac, offsetCp);

		for (unsigned int colNode = 0; colNode < nBulkPoints; colNode++, jacP += _strideBound)
		{
			for (unsigned int comp = 0; comp < _nComp; comp++, ++jacC, ++jacP) {

				// add Cl on Cl entries (added since already set in bulk jacobian)
				// row: already at bulk phase. already at current node and component.
				// col: already at bulk phase. already at current node and component.
				if (!crossDepsOnly)
					jacC[0] += jacCF_val * static_cast<double>(filmDiff[comp]) * static_cast<double>(parTypeVolFrac[_parTypeIdx + nParType * colNode]);
				// add Cl on Cp entries
				// row: already at bulk phase. already at current node and component.
				// col: jump to particle phase
				jacC[jacP.row() - jacC.row()] = -jacCF_val * static_cast<double>(filmDiff[comp]) * static_cast<double>(parTypeVolFrac[_parTypeIdx + nParType * colNode]);

				// add Cp on Cp entries
				// row: already at particle. already at current node and liquid state.
				// col: already at particle. already at current node and liquid state.
				if (!crossDepsOnly)
					jacP[0] = -jacPF_val / static_cast<double>(poreAccFactor[comp]) * static_cast<double>(filmDiff[comp]);
				// add Cp on Cl entries
				// row: already at particle. already at current node and liquid state.
				// col: go to flux of current parType and adjust for offsetC. jump over previous colNodes and add component offset
				jacP[jacC.row() - jacP.row()] = jacPF_val / static_cast<double>(poreAccFactor[comp]) * static_cast<double>(filmDiff[comp]);
			}
		}

		return 1;
	}

	bool HomogeneousParticle::setParameter(const ParameterId& pId, double value)
	{
		if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_RADIUS"), _parRadiusParTypeDep, _parRadius, _parTypeIdx, value, nullptr))
			return true;
		if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_POROSITY"), _parPorosityParTypeDep, _parPorosity, _parTypeIdx, value, nullptr))
			return true;
		if (singleTypeMultiplexCompTypeSecParameterValue(pId, hashString("PORE_ACCESSIBILITY"), _poreAccessFactorMode, _poreAccessFactor, _nComp, _parTypeIdx, value, nullptr))
			return true;
		if (singleTypeMultiplexCompTypeSecParameterValue(pId, hashString("FILM_DIFFUSION"), _filmDiffusionMode, _filmDiffusion, _nComp, _parTypeIdx, value, nullptr))
			return true;

		return false;
	}

	bool HomogeneousParticle::setParameter(const ParameterId& pId, int value)
	{
		return false;
	}

	bool HomogeneousParticle::setParameter(const ParameterId& pId, bool value)
	{
		return false;
	}

	bool HomogeneousParticle::setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& pId, double value)
	{
		if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_RADIUS"), _parRadiusParTypeDep, _parRadius, _parTypeIdx, value, &sensParams))
			return true;

		if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_POROSITY"), _parPorosityParTypeDep, _parPorosity, _parTypeIdx, value, &sensParams))
			return true;

		if (singleTypeMultiplexCompTypeSecParameterValue(pId, hashString("PORE_ACCESSIBILITY"), _poreAccessFactorMode, _poreAccessFactor, _nComp, _parTypeIdx, value, &sensParams))
			return true;

		if (singleTypeMultiplexCompTypeSecParameterValue(pId, hashString("FILM_DIFFUSION"), _filmDiffusionMode, _filmDiffusion, _nComp, _parTypeIdx, value, &sensParams))
			return true;

		return false;
	}

	bool HomogeneousParticle::setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue)
	{
		if (singleTypeMultiplexTypeParameterAD(pId, hashString("PAR_RADIUS"), _parRadiusParTypeDep, _parRadius, _parTypeIdx, adDirection, adValue, sensParams))
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		if (singleTypeMultiplexTypeParameterAD(pId, hashString("PAR_POROSITY"), _parPorosityParTypeDep, _parPorosity, _parTypeIdx, adDirection, adValue, sensParams))
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		if (singleTypeMultiplexCompTypeSecParameterAD(pId, hashString("PORE_ACCESSIBILITY"), _poreAccessFactorMode, _poreAccessFactor, _nComp, _parTypeIdx, adDirection, adValue, sensParams))
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		if (singleTypeMultiplexCompTypeSecParameterAD(pId, hashString("FILM_DIFFUSION"), _filmDiffusionMode, _filmDiffusion, _nComp, _parTypeIdx, adDirection, adValue, sensParams))
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		return false;
	}

	std::unordered_map<ParameterId, double> HomogeneousParticle::getAllParameterValues(std::unordered_map<ParameterId, double>& data) const
	{
		return data;
	}

	double HomogeneousParticle::getParameterDouble(const ParameterId& pId) const
	{
		return static_cast<double>(false);
	}

	bool HomogeneousParticle::hasParameter(const ParameterId& pId) const
	{
		return false;
	}

	void registerHomogeneousParticleModel(std::unordered_map<std::string, std::function<model::IParticleModel* ()>>& particles)
	{
		particles[HomogeneousParticle::identifier()] = []() { return new HomogeneousParticle(); };
	}

}  // namespace model

}  // namespace cadet
