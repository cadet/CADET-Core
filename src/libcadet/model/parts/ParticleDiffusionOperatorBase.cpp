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

#include "model/parts/ParticleDiffusionOperatorBase.hpp"
#include "model/ParameterDependence.hpp"
#include "ConfigurationHelper.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

namespace cadet
{

namespace model
{

namespace parts
{
	ParticleDiffusionOperatorBase::ParticleDiffusionOperatorBase() : _boundOffset(nullptr), _reqBinding(nullptr), _parDepSurfDiffusion(nullptr)
	{
	}

	ParticleDiffusionOperatorBase::~ParticleDiffusionOperatorBase() CADET_NOEXCEPT
	{
		delete[] _boundOffset;
		delete[] _parDepSurfDiffusion;
	}
	

	bool ParticleDiffusionOperatorBase::configureModelDiscretization_old(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp)
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
			nBound = paramProvider.getIntArray("NBOUND");
		else
		{
			paramProvider.popScope();
			nBound = paramProvider.getIntArray("NBOUND");
			paramProvider.pushScope("discretization");
		}
		if (nBound.size() < _nComp)
			throw InvalidParameterException("Field NBOUND contains too few elements (NCOMP = " + std::to_string(_nComp) + " required)");

		std::vector<int> stridesParTypeBound(nParType + 1);
		std::vector<int> nBoundBeforeType(nParType);
		if (!_nBound)
			_nBound = std::make_shared<unsigned int []>(_nComp);

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

		// ==== Construct and configure parameter dependencies
		paramProvider.popScope();

		delete _parDepSurfDiffusion;
		bool parSurfDiffDepConfSuccess = true;
		if (paramProvider.exists("PAR_SURFDIFFUSION_DEP"))
		{
			const std::vector<std::string> psdDepNames = paramProvider.getStringArray("PAR_SURFDIFFUSION_DEP");
			if ((psdDepNames.size() == 1) || (nParType == 1))
				_paramDepSurfDiffTypeDep =  false;

			if (!_paramDepSurfDiffTypeDep && (psdDepNames.size() < nParType))
				throw InvalidParameterException("Field PAR_SURFDIFFUSION_DEP contains too few elements (" + std::to_string(nParType) + " required)");
			else if (_paramDepSurfDiffTypeDep && (psdDepNames.size() != 1))
				throw InvalidParameterException("Field PAR_SURFDIFFUSION_DEP requires (only) 1 element");

			if (_paramDepSurfDiffTypeDep)
			{
				if ((psdDepNames[0] == "") || (psdDepNames[0] == "NONE") || (psdDepNames[0] == "DUMMY"))
				{
					_hasParDepSurfDiffusion = false;
					_paramDepSurfDiffTypeDep =  false;
					_parDepSurfDiffusion = nullptr;
				}
				else
				{
					_parDepSurfDiffusion = helper.createParameterStateDependence(psdDepNames[0]);
					if (!_parDepSurfDiffusion)
						throw InvalidParameterException("Unknown parameter dependence " + psdDepNames[0]);

					parSurfDiffDepConfSuccess = _parDepSurfDiffusion->configureModelDiscretization(paramProvider, _nComp, _nBound.get(), _boundOffset);
					_hasParDepSurfDiffusion = true;
				}
			}
			else
			{
				if (!(psdDepNames[0] == "") || (psdDepNames[0] == "NONE") || (psdDepNames[0] == "DUMMY"))
				{
					_parDepSurfDiffusion = helper.createParameterStateDependence(psdDepNames[_parTypeIdx]);
					if (!_parDepSurfDiffusion)
						throw InvalidParameterException("Unknown parameter dependence " + psdDepNames[parTypeIdx]);

					parSurfDiffDepConfSuccess = _parDepSurfDiffusion->configureModelDiscretization(paramProvider, _nComp, _nBound.get(), _boundOffset) && parSurfDiffDepConfSuccess;
				}

				_hasParDepSurfDiffusion = _parDepSurfDiffusion;
			}
		}
		else
		{
			_hasParDepSurfDiffusion = false;
			_paramDepSurfDiffTypeDep =  false;
			_parDepSurfDiffusion = nullptr;
		}

		// Check whether surface diffusion is present
		_hasSurfaceDiffusion = false;
		if (paramProvider.exists("PAR_SURFDIFFUSION"))
		{
			const std::vector<double> surfDiff = paramProvider.getDoubleArray("PAR_SURFDIFFUSION");
			// Assume particle surface diffusion if a parameter dependence is present
			if (_parDepSurfDiffusion)
			{
				_hasSurfaceDiffusion = true;
			}
			else
			{
				double const* const lsd = surfDiff.data() + nBoundBeforeType[_parTypeIdx];

				// Check surface diffusion coefficients
				for (unsigned int j = 0; j < _strideBound; ++j)
				{
					if (lsd[j] != 0.0)
					{
						_hasSurfaceDiffusion = true;
						break;
					}
				}
			}
		}

		paramProvider.pushScope("discretization");

		return parSurfDiffDepConfSuccess;
	}

	bool ParticleDiffusionOperatorBase::configure_old(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound, const int* reqBinding, const bool hasDynamicReactions)
	{
		_reqBinding = reqBinding;
		_hasDynamicReactions = hasDynamicReactions;

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
		// Let PAR_CORERADIUS default to 0.0 for backwards compatibility
		if (paramProvider.exists("PAR_CORERADIUS"))
		{
			std::vector<double> parCoreRadii(nParType);
			_parCoreRadiusParTypeDep = readScalarParameterOrArray(parCoreRadii, paramProvider, "PAR_CORERADIUS", nParType);
			if (_parCoreRadiusParTypeDep)
			{
				_parCoreRadius = parCoreRadii[0];
				if (_parTypeIdx == 0)
					parameters[makeParamId(hashStringRuntime("PAR_CORERADIUS"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parCoreRadius;
			}
			else
			{
				_parCoreRadius = parCoreRadii[_parTypeIdx];
				parameters[makeParamId(hashStringRuntime("PAR_CORERADIUS"), unitOpIdx, CompIndep, _parTypeIdx, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parCoreRadius;
			}
		}
		else
		{
			_parCoreRadiusParTypeDep =  false;
			_parCoreRadius = 0.0;
		}

		_filmDiffusionMode = readAndRegisterSingleTypeMultiplexCompTypeSecParam(paramProvider, parameters, _filmDiffusion, "FILM_DIFFUSION", nParType, _nComp, _parTypeIdx, unitOpIdx);

		if (paramProvider.exists("PORE_ACCESSIBILITY"))
			_poreAccessFactorMode = readAndRegisterSingleTypeMultiplexCompTypeSecParam(paramProvider, parameters, _poreAccessFactor, "PORE_ACCESSIBILITY", nParType, _nComp, _parTypeIdx, unitOpIdx);
		else
		{
			_poreAccessFactorMode = MultiplexMode::ComponentType;
			_poreAccessFactor = std::vector<cadet::active>(_nComp, 1.0);
		}

		_invBetaP.resize(_nComp);
		for (int comp = 0; comp < _nComp; comp++)
			_invBetaP[comp] = (1.0 - _parPorosity) / (_poreAccessFactor[comp] * _parPorosity);

		_parDiffusionMode = readAndRegisterSingleTypeMultiplexCompTypeSecParam(paramProvider, parameters, _parDiffusion, "PAR_DIFFUSION", nParType, _nComp, _parTypeIdx, unitOpIdx);

		if (paramProvider.exists("PAR_SURFDIFFUSION"))
			_parSurfDiffusionMode = readAndRegisterSingleTypeMultiplexBndCompTypeSecParam(paramProvider, parameters, _parSurfDiffusion, "PAR_SURFDIFFUSION", nTotalBound, _nComp, _strideBound, _nBound.get(), nBoundBeforeType, _parTypeIdx, unitOpIdx);
		else
		{
			_parSurfDiffusionMode = MultiplexMode::Component;
			_parSurfDiffusion.resize(_strideBound, 0.0);
		}

		bool parSurfDiffDepConfSuccess = true;
		if (_hasParDepSurfDiffusion)
		{
			if (_paramDepSurfDiffTypeDep && _parDepSurfDiffusion)
			{
				parSurfDiffDepConfSuccess = _parDepSurfDiffusion->configure(paramProvider, unitOpIdx, ParTypeIndep, "PAR_SURFDIFFUSION");
			}
			else if (!_paramDepSurfDiffTypeDep && _parDepSurfDiffusion)
			{
				parSurfDiffDepConfSuccess = _parDepSurfDiffusion->configure(paramProvider, unitOpIdx, _parTypeIdx, "PAR_SURFDIFFUSION") && parSurfDiffDepConfSuccess;
			}
		}

		return true;
	}

	bool ParticleDiffusionOperatorBase::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp)
	{
		_strideBulkComp = strideBulkComp;

		_parTypeIdx = parTypeIdx;
		_nComp = nComp;

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

		// Precompute offsets and total number of bound states (DOFs in solid phase)
		if (!_boundOffset)
			_boundOffset = new unsigned int[_nComp];

		_boundOffset[0] = 0.0;
		_strideBound = std::accumulate(_nBound.get(), _nBound.get() + _nComp, 0u);

		for (unsigned int i = 1; i < _nComp; ++i)
		{
			_boundOffset[i] = _boundOffset[i - 1] + _nBound[i - 1];
		}

		// ==== Construct and configure parameter dependencies

		delete _parDepSurfDiffusion;
		bool parSurfDiffDepConfSuccess = true;
		if (paramProvider.exists("PAR_SURFDIFFUSION_DEP"))
		{
			const std::string psdDepName = paramProvider.getString("PAR_SURFDIFFUSION_DEP");
			
			_paramDepSurfDiffTypeDep = paramProvider.exists("PAR_SURFDIFFUSIN_DEP_PARTYPE_DEPENDENT") ? paramProvider.getBool("PAR_SURFDIFFUSIN_DEP_PARTYPE_DEPENDENT") : true;

			if ((psdDepName == "") || (psdDepName == "NONE") || (psdDepName == "DUMMY"))
			{
				_hasParDepSurfDiffusion = false;
				_parDepSurfDiffusion = nullptr;
			}
			else
			{
				_parDepSurfDiffusion = helper.createParameterStateDependence(psdDepName);
				if (!_parDepSurfDiffusion)
					throw InvalidParameterException("Unknown parameter dependence " + psdDepName);

				parSurfDiffDepConfSuccess = _parDepSurfDiffusion->configureModelDiscretization(paramProvider, _nComp, _nBound.get(), _boundOffset);
				_hasParDepSurfDiffusion = true;
			}
		}
		else
		{
			_hasParDepSurfDiffusion = false;
			_paramDepSurfDiffTypeDep =  false;
			_parDepSurfDiffusion = nullptr;
		}

		// Check whether surface diffusion is present
		_hasSurfaceDiffusion = false;
		if (paramProvider.exists("PAR_SURFDIFFUSION"))
		{
			const std::vector<double> surfDiff = paramProvider.getDoubleArray("PAR_SURFDIFFUSION");
			// Assume particle surface diffusion if a parameter dependence is present
			if (_parDepSurfDiffusion)
			{
				_hasSurfaceDiffusion = true;
			}
			else
			{
				double const* const lsd = surfDiff.data();

				// Check surface diffusion coefficients
				for (unsigned int j = 0; j < _strideBound; ++j)
				{
					if (lsd[j] != 0.0)
					{
						_hasSurfaceDiffusion = true;
						break;
					}
				}
			}
		}

		return parSurfDiffDepConfSuccess;
	}

	bool ParticleDiffusionOperatorBase::configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound, const int* reqBinding, const bool hasDynamicReactions)
	{
		_reqBinding = reqBinding;
		_hasDynamicReactions = hasDynamicReactions;

		// Read geometry parameters
		_parRadius = paramProvider.getDouble("PAR_RADIUS");

		_parRadiusParTypeDep = true;
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
			_parRadiusParTypeDep = paramProvider.getBool("PAR_POROSITY_PARTYPE_DEPENDENT");

		if (!_parRadiusParTypeDep && _parTypeIdx == 0)
			parameters[makeParamId(hashStringRuntime("PAR_POROSITY"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parPorosity;
		else if (_parRadiusParTypeDep)
			parameters[makeParamId(hashStringRuntime("PAR_POROSITY"), unitOpIdx, CompIndep, _parTypeIdx, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parPorosity;

		if (_parPorosity <= 0.0 || _parPorosity > 1.0)
			throw InvalidParameterException("Particle porosity is not within (0, 1] for particle type " + std::to_string(_parTypeIdx));

		if (paramProvider.exists("PAR_CORERADIUS"))
		{
			_parCoreRadius = paramProvider.getDouble("PAR_CORERADIUS");

			_parCoreRadiusParTypeDep = true;
			if (paramProvider.exists("PAR_CORERADIUS_PARTYPE_DEPENDENT"))
				_parCoreRadiusParTypeDep = paramProvider.getBool("PAR_CORERADIUS_PARTYPE_DEPENDENT");

			if (!_parRadiusParTypeDep && _parTypeIdx == 0)
				parameters[makeParamId(hashStringRuntime("PAR_CORERADIUS"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parCoreRadius;
			else if (_parCoreRadiusParTypeDep)
				parameters[makeParamId(hashStringRuntime("PAR_CORERADIUS"), unitOpIdx, CompIndep, _parTypeIdx, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parCoreRadius;

			if (_parRadius <= 0.0)
				throw InvalidParameterException("Particle radius is not bigger than zero for particle type " + std::to_string(_parTypeIdx));

		}
		else
		{
			_parCoreRadiusParTypeDep = true;
			_parCoreRadius = 0.0;
		}

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

		_invBetaP.resize(_nComp);
		for (int comp = 0; comp < _nComp; comp++)
			_invBetaP[comp] = (1.0 - _parPorosity) / (_poreAccessFactor[comp] * _parPorosity);

		bool parDiffParTypeDep = paramProvider.exists("PAR_DIFFUSION_PARTYPE_DEPENDENT") ? paramProvider.getBool("PAR_DIFFUSION_PARTYPE_DEPENDENT") : true;
		_parDiffusionMode = newIF_readAndRegisterMultiplexCompSecParam(paramProvider, parameters, _parDiffusion, "PAR_DIFFUSION", _nComp, _parTypeIdx, parDiffParTypeDep, unitOpIdx);

		if (paramProvider.exists("PAR_SURFDIFFUSION"))
		{
			bool parSurfDiffParTypeDep = paramProvider.exists("PAR_SURFDIFFUSION_PARTYPE_DEPENDENT") ? paramProvider.getBool("PAR_SURFDIFFUSION_PARTYPE_DEPENDENT") : true;
			_parSurfDiffusionMode = newIF_readAndRegisterMultiplexBndCompSecParam(paramProvider, parameters, _parSurfDiffusion, "PAR_SURFDIFFUSION", nTotalBound, _nComp, _strideBound, _nBound.get(), _parTypeIdx, parSurfDiffParTypeDep, unitOpIdx);
		}
		else
		{
			_parSurfDiffusionMode = MultiplexMode::Component;
			_parSurfDiffusion.resize(_strideBound, 0.0);
		}

		bool parSurfDiffDepConfSuccess = true;
		if (_hasParDepSurfDiffusion)
		{
			if (!_paramDepSurfDiffTypeDep && _parDepSurfDiffusion)
			{
				parSurfDiffDepConfSuccess = _parDepSurfDiffusion->configure(paramProvider, unitOpIdx, ParTypeIndep, "PAR_SURFDIFFUSION");
			}
			else if (_paramDepSurfDiffTypeDep && _parDepSurfDiffusion)
			{
				parSurfDiffDepConfSuccess = _parDepSurfDiffusion->configure(paramProvider, unitOpIdx, _parTypeIdx, "PAR_SURFDIFFUSION") && parSurfDiffDepConfSuccess;
			}
		}

		return parSurfDiffDepConfSuccess;
	}

	bool ParticleDiffusionOperatorBase::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx)
	{
		return true;
	}

	bool ParticleDiffusionOperatorBase::setParameter(const ParameterId& pId, double value)
	{
		if (singleTypeMultiplexCompTypeSecParameterValue(pId, hashString("PAR_DIFFUSION"), _parDiffusionMode, _parDiffusion, _nComp, _parTypeIdx, value, nullptr))
			return true;
		if (singleTypeMultiplexBndCompTypeSecParameterValue(pId, hashString("PAR_SURFDIFFUSION"), _parSurfDiffusionMode, _parSurfDiffusion, _nComp, _strideBound, _boundOffset, _parTypeIdx, value, nullptr))
			return true;

		if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_RADIUS"), _parRadiusParTypeDep, _parRadius, _parTypeIdx, value, nullptr))
			return true;
		if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_CORERADIUS"), _parCoreRadiusParTypeDep, _parCoreRadius, _parTypeIdx, value, nullptr))
			return true;
		if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_POROSITY"), _parPorosityParTypeDep, _parPorosity, _parTypeIdx, value, nullptr))
			return true;
		if (singleTypeMultiplexCompTypeSecParameterValue(pId, hashString("PORE_ACCESSIBILITY"), _poreAccessFactorMode, _poreAccessFactor, _nComp, _parTypeIdx, value, nullptr))
			return true;
		if (singleTypeMultiplexCompTypeSecParameterValue(pId, hashString("FILM_DIFFUSION"), _filmDiffusionMode, _filmDiffusion, _nComp, _parTypeIdx, value, nullptr))
			return true;

		if (!_paramDepSurfDiffTypeDep || pId.particleType == _parTypeIdx)
			if (_parDepSurfDiffusion && _parDepSurfDiffusion->setParameter(pId, value))
				return true;

		return false;
	}

	bool ParticleDiffusionOperatorBase::setParameter(const ParameterId& pId, int value)
	{
		if (!_paramDepSurfDiffTypeDep || pId.particleType == _parTypeIdx)
			if (_parDepSurfDiffusion && _parDepSurfDiffusion->setParameter(pId, value))
				return true;

		return false;
	}

	bool ParticleDiffusionOperatorBase::setParameter(const ParameterId& pId, bool value)
	{
		if (!_paramDepSurfDiffTypeDep || pId.particleType == _parTypeIdx)
			if (_parDepSurfDiffusion && _parDepSurfDiffusion->setParameter(pId, value))
				return true;

		return false;
	}

	bool ParticleDiffusionOperatorBase::setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& pId, double value)
	{
		if (singleTypeMultiplexCompTypeSecParameterValue(pId, hashString("PAR_DIFFUSION"), _parDiffusionMode, _parDiffusion, _nComp, _parTypeIdx, value, &sensParams))
			return true;
		if (singleTypeMultiplexBndCompTypeSecParameterValue(pId, hashString("PAR_SURFDIFFUSION"), _parSurfDiffusionMode, _parSurfDiffusion, _nComp, _strideBound, _boundOffset, _parTypeIdx, value, &sensParams))
			return true;

		if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_RADIUS"), _parRadiusParTypeDep, _parRadius, _parTypeIdx, value, &sensParams))
			return true;
		if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_CORERADIUS"), _parCoreRadiusParTypeDep, _parCoreRadius, _parTypeIdx, value, &sensParams))
			return true;
		if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_POROSITY"), _parPorosityParTypeDep, _parPorosity, _parTypeIdx, value, &sensParams))
			return true;

		if (singleTypeMultiplexCompTypeSecParameterValue(pId, hashString("PORE_ACCESSIBILITY"), _poreAccessFactorMode, _poreAccessFactor, _nComp, _parTypeIdx, value, &sensParams))
			return true;
		if (singleTypeMultiplexCompTypeSecParameterValue(pId, hashString("FILM_DIFFUSION"), _filmDiffusionMode, _filmDiffusion, _nComp, _parTypeIdx, value, &sensParams))
			return true;

		if (model::setSensitiveParameterValue(pId, value, sensParams, std::vector< IParameterStateDependence*>(1, _parDepSurfDiffusion), _paramDepSurfDiffTypeDep))
			return true;

		return false;
	}

	bool ParticleDiffusionOperatorBase::setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue)
	{
		if (singleTypeMultiplexCompTypeSecParameterAD(pId, hashString("PAR_DIFFUSION"), _parDiffusionMode, _parDiffusion, _nComp, _parTypeIdx, adDirection, adValue, sensParams))
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		if (singleTypeMultiplexBndCompTypeSecParameterAD(pId, hashString("PAR_SURFDIFFUSION"), _parSurfDiffusionMode, _parSurfDiffusion, _nComp, _strideBound, _boundOffset, _parTypeIdx, adDirection, adValue, sensParams))
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		if (model::setSensitiveParameter(pId, adDirection, adValue, sensParams, std::vector< IParameterStateDependence*>(1, _parDepSurfDiffusion), _paramDepSurfDiffTypeDep))
		{
			LOG(Debug) << "Found parameter " << pId << " in surface diffusion parameter dependence: Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		if (singleTypeMultiplexTypeParameterAD(pId, hashString("PAR_RADIUS"), _parRadiusParTypeDep, _parRadius, _parTypeIdx, adDirection, adValue, sensParams))
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		if (singleTypeMultiplexTypeParameterAD(pId, hashString("PAR_CORERADIUS"), _parCoreRadiusParTypeDep, _parCoreRadius, _parTypeIdx, adDirection, adValue, sensParams))
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


} // namespace parts
} // namespace model
} // namespace cadet