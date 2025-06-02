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

namespace parts
{
	/**
	 * @brief Creates a ParticleDiffusionOperatorDG
	 */
	ParticleDiffusionOperatorDG::ParticleDiffusionOperatorDG() : _binding(nullptr), _boundOffset(nullptr), _nBound(nullptr), _localFlux(nullptr), _deltaR(nullptr), _Ir(nullptr), _DGjacParDispBlocks(nullptr), _minus_InvMM_ST(nullptr), _parInvMM(nullptr), _parDepSurfDiffusion(nullptr)
	{
	}

	ParticleDiffusionOperatorDG::~ParticleDiffusionOperatorDG() CADET_NOEXCEPT
	{
		delete _binding;
		delete _dynReaction;
		delete[] _boundOffset;
		delete[] _nBound;
		delete[] _localFlux;
		delete[] _deltaR;
		delete[] _Ir;
		delete[] _DGjacParDispBlocks;
		delete[] _minus_InvMM_ST;
		delete[] _parInvMM;
		delete[] _parDepSurfDiffusion;
	}

	bool ParticleDiffusionOperatorDG::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp)
	{
		_strideBulkComp = strideBulkComp;

		_parTypeIdx = parTypeIdx;
		_nComp = nComp;

		_filmDiffusion.resize(_nComp); // filled in notifyDiscontinuousSectionTransition
		_poreAccessFactor.resize(_nComp); // filled in notifyDiscontinuousSectionTransition
		_invBetaP.resize(_nComp); // filled in notifyDiscontinuousSectionTransition
		
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

		std::vector<int> parPolyDegs(nParType);
		std::vector<int> parNelements(nParType);

		int nTotalParPoints = 0;

		if (paramProvider.exists("PAR_POLYDEG"))
		{
			std::vector<int> nParPoints(nParType);

			parPolyDegs = paramProvider.getIntArray("PAR_POLYDEG");

			if ((std::any_of(parPolyDegs.begin(), parPolyDegs.end(), [](int value) { return value < 1; })))
				throw InvalidParameterException("Particle polynomial degree (PAR_POLYDEG) must be at least 1!");
			parNelements = paramProvider.getIntArray("PAR_NELEM");

			if ((std::any_of(parNelements.begin(), parNelements.end(), [](int value) { return value < 1; })))
				throw InvalidParameterException("Particle number of elements (PAR_NELEM) must be at least 1!");

			if (parPolyDegs.size() == 1)
				_parPolyDeg = parPolyDegs[0];
			else if (parPolyDegs.size() < nParType)
				throw InvalidParameterException("Field PAR_POLYDEG must have 1 or _nParType (" + std::to_string(nParType) + ") entries");
			else
				_parPolyDeg = parPolyDegs[_parTypeIdx];
			if (parNelements.size() == 1)
				_nParElem = parNelements[0];
			else if (parNelements.size() < nParType)
				throw InvalidParameterException("Field PAR_NELEM must have 1 or _nParType (" + std::to_string(nParType) + ") entries");
			else
				_nParElem = parNelements[_parTypeIdx];

			std::transform(parPolyDegs.begin(), parPolyDegs.end(), parNelements.begin(), nParPoints.begin(), [](int x, int y) { return (x + 1) * y; });
			nTotalParPoints = std::accumulate(nParPoints.begin(), nParPoints.end(), 0);
		}
		else if (paramProvider.exists("NPAR"))
		{
			const std::vector<int> nParPoints = paramProvider.getIntArray("NPAR");
			if ((nParPoints.size() > 1) && (nParPoints.size() < nParType))
				throw InvalidParameterException("Field NPAR must have 1 or NPARTYPE (" + std::to_string(nParType) + ") entries");

			if (nParPoints.size() == 1)
				_parPolyDeg = std::max(1, std::min(nParPoints[0] - 1, 4));
			else
				_parPolyDeg = std::max(1, std::min(nParPoints[_parTypeIdx] - 1, 4));

			nTotalParPoints = std::accumulate(nParPoints.begin(), nParPoints.end(), 0);
		}
		else
			throw InvalidParameterException("Specify field PAR_POLYDEG (or NPAR)");

		if (paramProvider.exists("PAR_GSM"))
		{
			std::vector<bool> parGSMs = paramProvider.getBoolArray("PAR_GSM");
			if (parGSMs.size() == 1)
				_parGSM = parGSMs[0];
			else if (parGSMs.size() < nParType)
				throw InvalidParameterException("Field PAR_GSM must have 1 or NPARTYPE (" + std::to_string(nParType) + ") entries");
			else
				_parGSM = parGSMs[_parTypeIdx];
			if (_parGSM && _nParElem != 1)
				throw InvalidParameterException("Field PAR_NELEM must equal one to use a GSM discretization in the corresponding particle type");
		}
		else // Use GSM as default for particle discretization
		{
			_parGSM = (_nParElem == 1);
		}

		initializeDG();

		std::vector<int> stridesParTypeBound(nParType + 1);
		std::vector<int> nBoundBeforeType(nParType);
		if (!_nBound)
			_nBound = new unsigned int[_nComp];

		if (nBound.size() < _nComp * nParType)
		{
			std::copy_n(nBound.begin(), _nComp, _nBound);

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
			std::copy_n(nBound.begin() + _parTypeIdx * _nComp, _nComp, _nBound);

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
		_strideBound = std::accumulate(_nBound, _nBound + _nComp, 0u);

		for (unsigned int i = 1; i < _nComp; ++i)
		{
			_boundOffset[i] = _boundOffset[i - 1] + _nBound[i - 1];
		}

		// Configure particle discretization
		_parElementSize.resize(_nParElem);
		_parCenterRadius.resize(_nParElem);
		_parOuterSurfAreaPerVolume.resize(_nParElem);
		_parInnerSurfAreaPerVolume.resize(_nParElem);

		// Read particle discretization mode and default to "EQUIDISTANT_PAR"
		_parDiscMode = ParticleDiscretizationMode::Equidistant;
		std::vector<std::string> pdt = paramProvider.getStringArray("PAR_DISC_TYPE");
		if ((pdt.size() == 1) && (nParType > 1))
		{
			// Multiplex using first value
			pdt.resize(nParType, pdt[0]);
		}
		else if (pdt.size() < nParType)
			throw InvalidParameterException("Field PAR_DISC_TYPE contains too few elements (" + std::to_string(nParType) + " required)");

		if (pdt[_parTypeIdx] == "EQUIVOLUME_PAR")
			_parDiscMode = ParticleDiscretizationMode::Equivolume;
		else if (pdt[_parTypeIdx] == "USER_DEFINED_PAR")
			_parDiscMode = ParticleDiscretizationMode::UserDefined;

		// Read particle geometry and default to "SPHERICAL"
		paramProvider.popScope();
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
		paramProvider.pushScope("discretization");

		if (paramProvider.exists("PAR_DISC_VECTOR"))
		{
			std::vector<double> pdv = paramProvider.getDoubleArray("PAR_DISC_VECTOR");
			if (pdv.size() < nTotalParPoints + nParType)
				throw InvalidParameterException("Field PAR_DISC_VECTOR contains too few elements (Sum [NPAR + 1] = " + std::to_string(nTotalParPoints + nParType) + " required)");
		}

		paramProvider.popScope();

		// ==== Construct and configure parameter dependencies
		delete _parDepSurfDiffusion;
		bool parSurfDiffDepConfSuccess = true;
		if (paramProvider.exists("PAR_SURFDIFFUSION_DEP"))
		{
			const std::vector<std::string> psdDepNames = paramProvider.getStringArray("PAR_SURFDIFFUSION_DEP");
			if ((psdDepNames.size() == 1) || (nParType == 1))
				_singleParDepSurfDiffusion = true;

			if (!_singleParDepSurfDiffusion && (psdDepNames.size() < nParType))
				throw InvalidParameterException("Field PAR_SURFDIFFUSION_DEP contains too few elements (" + std::to_string(nParType) + " required)");
			else if (_singleParDepSurfDiffusion && (psdDepNames.size() != 1))
				throw InvalidParameterException("Field PAR_SURFDIFFUSION_DEP requires (only) 1 element");

			if (_singleParDepSurfDiffusion)
			{
				if ((psdDepNames[0] == "") || (psdDepNames[0] == "NONE") || (psdDepNames[0] == "DUMMY"))
				{
					_hasParDepSurfDiffusion = false;
					_singleParDepSurfDiffusion = true;
					_parDepSurfDiffusion = nullptr;
				}
				else
				{
					_parDepSurfDiffusion = helper.createParameterStateDependence(psdDepNames[0]);
					if (!_parDepSurfDiffusion)
						throw InvalidParameterException("Unknown parameter dependence " + psdDepNames[0]);

					parSurfDiffDepConfSuccess = _parDepSurfDiffusion->configureModelDiscretization(paramProvider, _nComp, _nBound, _boundOffset);
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

					parSurfDiffDepConfSuccess = _parDepSurfDiffusion->configureModelDiscretization(paramProvider, _nComp, _nBound, _boundOffset) && parSurfDiffDepConfSuccess;
				}

				_hasParDepSurfDiffusion = _parDepSurfDiffusion;
			}
		}
		else
		{
			_hasParDepSurfDiffusion = false;
			_singleParDepSurfDiffusion = true;
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
				for (unsigned int j = 0; j < stridesParTypeBound[_parTypeIdx]; ++j)
				{
					if (lsd[j] != 0.0)
					{
						_hasSurfaceDiffusion = true;
						break;
					}
				}
			}
		}

		// ==== Construct and configure binding model
		_binding = nullptr;
		std::vector<std::string> bindModelNames = { "NONE" };
		if (paramProvider.exists("ADSORPTION_MODEL"))
			bindModelNames = paramProvider.getStringArray("ADSORPTION_MODEL");

		if (paramProvider.exists("ADSORPTION_MODEL_MULTIPLEX"))
			_singleBinding = (paramProvider.getInt("ADSORPTION_MODEL_MULTIPLEX") == 1);
		else
		{
			// Infer multiplex mode
			_singleBinding = (bindModelNames.size() == 1);
		}

		if (!_singleBinding && (bindModelNames.size() < nParType))
			throw InvalidParameterException("Field ADSORPTION_MODEL contains too few elements (" + std::to_string(nParType) + " required)");
		else if (_singleBinding && (bindModelNames.size() != 1))
			throw InvalidParameterException("Field ADSORPTION_MODEL requires (only) 1 element");

		bool bindingConfSuccess = true;

		_binding = helper.createBindingModel(bindModelNames[_singleBinding ? 0 : _parTypeIdx]);
		if (!_binding)
			throw InvalidParameterException("Unknown binding model " + bindModelNames[_singleBinding ? 0 : _parTypeIdx]);

		MultiplexedScopeSelector scopeGuard(paramProvider, "adsorption", _singleBinding, _parTypeIdx, nParType == 1, _binding->usesParamProviderInDiscretizationConfig());
		bindingConfSuccess = _binding->configureModelDiscretization(paramProvider, _nComp, _nBound, _boundOffset);

		// ==== Construct and configure dynamic reaction model
		bool reactionConfSuccess = true;

		_dynReaction = nullptr;

		if (paramProvider.exists("REACTION_MODEL_PARTICLES"))
		{
			const std::vector<std::string> dynReactModelNames = paramProvider.getStringArray("REACTION_MODEL_PARTICLES");

			if (paramProvider.exists("REACTION_MODEL_PARTICLES_MULTIPLEX"))
				_singleDynReaction = (paramProvider.getInt("REACTION_MODEL_PARTICLES_MULTIPLEX") == 1);
			else
			{
				// Infer multiplex mode
				_singleDynReaction = (dynReactModelNames.size() == 1);
			}

			if (!_singleDynReaction && (dynReactModelNames.size() < nParType))
				throw InvalidParameterException("Field REACTION_MODEL_PARTICLES contains too few elements (" + std::to_string(nParType) + " required)");
			else if (_singleDynReaction && (dynReactModelNames.size() != 1))
				throw InvalidParameterException("Field REACTION_MODEL_PARTICLES requires (only) 1 element");

			_dynReaction = helper.createDynamicReactionModel(dynReactModelNames[_singleDynReaction ? 0 : _parTypeIdx]);

			if (!_dynReaction)
				throw InvalidParameterException("Unknown dynamic reaction model " + dynReactModelNames[_singleDynReaction ? 0 : _parTypeIdx]);

			MultiplexedScopeSelector scopeGuard(paramProvider, "reaction_particle", _singleDynReaction, _parTypeIdx, nParType == 1, _dynReaction->usesParamProviderInDiscretizationConfig());
			reactionConfSuccess = _dynReaction->configureModelDiscretization(paramProvider, _nComp, _nBound, _boundOffset) && reactionConfSuccess;
		}

		return parSurfDiffDepConfSuccess && bindingConfSuccess && reactionConfSuccess;
	}

	bool ParticleDiffusionOperatorDG::configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound)
	{
		// Read geometry parameters
		std::vector<double> parRadii(nParType);
		_singleParRadius = readScalarParameterOrArray(parRadii, paramProvider, "PAR_RADIUS", nParType);
		
		if (_singleParRadius)
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
		_singleParPorosity = readScalarParameterOrArray(parPorosities, paramProvider, "PAR_POROSITY", nParType);
		if (_singleParPorosity)
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
			_singleParCoreRadius = readScalarParameterOrArray(parCoreRadii, paramProvider, "PAR_CORERADIUS", nParType);
			if (_singleParCoreRadius)
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
			_singleParCoreRadius = true;
			_parCoreRadius = 0.0;
		}

		// Check whether PAR_TYPE_VOLFRAC is required or not
		if ((nParType > 1) && !paramProvider.exists("PAR_TYPE_VOLFRAC"))
			throw InvalidParameterException("The required parameter \"PAR_TYPE_VOLFRAC\" was not found");

		// todo: PAR_TYPE_VOLFRAC remains parameter of the unit, not the particles?
		//// Let PAR_TYPE_VOLFRAC default to 1.0 for backwards compatibility
		//if (paramProvider.exists("PAR_TYPE_VOLFRAC"))
		//{
		//	readScalarParameterOrArray(_parTypeVolFrac, paramProvider, "PAR_TYPE_VOLFRAC", 1);
		//	if (_parTypeVolFrac.size() == _nParType)
		//	{
		//		_axiallyConstantParTypeVolFrac = true;

		//		// Expand to all axial elements
		//		_parTypeVolFrac.resize(nPoints * _nParType, 1.0);
		//		for (unsigned int i = 1; i < nPoints; ++i)
		//			std::copy(_parTypeVolFrac.begin(), _parTypeVolFrac.begin() + _nParType, _parTypeVolFrac.begin() + _nParType * i);
		//	}
		//	else
		//		_axiallyConstantParTypeVolFrac = false;
		//}
		//else
		//{
		//	_parTypeVolFrac.resize(nPoints, 1.0);
		//	_axiallyConstantParTypeVolFrac = false;
		//}

		//// Check that particle volume fractions sum to 1.0
		//for (unsigned int i = 0; i < nPoints; ++i)
		//{
		//	const double volFracSum = std::accumulate(_parTypeVolFrac.begin() + i * _nParType, _parTypeVolFrac.begin() + (i + 1) * _nParType, 0.0,
		//		[](double a, const active& b) -> double { return a + static_cast<double>(b); });
		//	if (std::abs(1.0 - volFracSum) > 1e-10)
		//		throw InvalidParameterException("Sum of field PAR_TYPE_VOLFRAC differs from 1.0 (is " + std::to_string(volFracSum) + ") in axial elem " + std::to_string(i));
		//}

		// Read vectorial parameters (which may also be section dependent; transport)
		//_filmDiffusionMode = readAndRegisterMultiplexCompTypeSecParam(paramProvider, parameters, _filmDiffusion, "FILM_DIFFUSION", _nParType, nComp, unitOpIdx); // todo film diffusion remains unit operation parameter?

		_parDiffusionMode = readAndRegisterSingleTypeMultiplexCompTypeSecParam(paramProvider, parameters, _parDiffusion, "PAR_DIFFUSION", nParType, _nComp, _parTypeIdx, unitOpIdx);

		if (paramProvider.exists("PAR_SURFDIFFUSION"))
			_parSurfDiffusionMode = readAndRegisterSingleTypeMultiplexBndCompTypeSecParam(paramProvider, parameters, _parSurfDiffusion, "PAR_SURFDIFFUSION", nTotalBound, _nComp, _strideBound, _nBound, nBoundBeforeType, _parTypeIdx, unitOpIdx);
		else
		{
			_parSurfDiffusionMode = MultiplexMode::Component;
			_parSurfDiffusion.resize(_strideBound, 0.0);
		}

		bool parSurfDiffDepConfSuccess = true;
		if (_hasParDepSurfDiffusion)
		{
			if (_singleParDepSurfDiffusion && _parDepSurfDiffusion)
			{
				parSurfDiffDepConfSuccess = _parDepSurfDiffusion->configure(paramProvider, unitOpIdx, ParTypeIndep, "PAR_SURFDIFFUSION");
			}
			else if (!_singleParDepSurfDiffusion && _parDepSurfDiffusion)
			{
				parSurfDiffDepConfSuccess = _parDepSurfDiffusion->configure(paramProvider, unitOpIdx, _parTypeIdx, "PAR_SURFDIFFUSION") && parSurfDiffDepConfSuccess;
			}
		}

		//if ((_filmDiffusion.size() < nComp * _nParType) || (_filmDiffusion.size() % (nComp * _nParType) != 0))
		//	throw InvalidParameterException("Number of elements in field FILM_DIFFUSION is not a positive multiple of NCOMP * _nParType (" + std::to_string(nComp * _nParType) + ")");
		if ((_parDiffusion.size() < _nComp) || (_parDiffusion.size() % (_nComp) != 0))
			throw InvalidParameterException("Number of elements in field PAR_DIFFUSION is not a positive multiple of NCOMP * NPARTYPE (" + std::to_string(_nComp * nParType) + ")");
		if ((_parSurfDiffusion.size() < _strideBound) || ((nTotalBound > 0) && (_parSurfDiffusion.size() % _strideBound != 0)))
			throw InvalidParameterException("Number of elements in field PAR_SURFDIFFUSION is not a positive multiple of NTOTALBND (" + std::to_string(nTotalBound) + ")");

		//if (paramProvider.exists("PORE_ACCESSIBILITY"))
		//	_poreAccessFactorMode = readAndRegisterMultiplexCompTypeSecParam(paramProvider, parameters, _poreAccessFactor, "PORE_ACCESSIBILITY", _nParType, nComp, unitOpIdx);
		//else
		//{
		//	_poreAccessFactorMode = MultiplexMode::ComponentType;
		//	_poreAccessFactor = std::vector<cadet::active>(nComp * _nParType, 1.0);
		//}

		//if (nComp * _nParType != _poreAccessFactor.size())
		//	throw InvalidParameterException("Number of elements in field PORE_ACCESSIBILITY differs from NCOMP * _nParType (" + std::to_string(nComp * _nParType) + ")");

		//// Add parameters to map
		//parameters[makeParamId(hashString("COL_POROSITY"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_colPorosity;

		//if (_axiallyConstantParTypeVolFrac)
		//{
		//	// Register only the first _nParType items
		//	for (unsigned int i = 0; i < _nParType; ++i)
		//		parameters[makeParamId(hashString("PAR_TYPE_VOLFRAC"), unitOpIdx, CompIndep, i, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parTypeVolFrac[i];
		//}
		//else
		//	registerParam2DArray(parameters, _parTypeVolFrac, [=](bool multi, unsigned elem, unsigned int type) { return makeParamId(hashString("PAR_TYPE_VOLFRAC"), unitOpIdx, CompIndep, type, BoundStateIndep, ReactionIndep, elem); }, _nParType);

		// Calculate the particle radial discretization variables (_parElementSize, _parCenterRadius, etc.)
		if (_deltaR == nullptr)
			_deltaR = new active[_nParElem];
		updateRadialDisc();

		//// Register initial conditions parameters
		//registerParam1DArray(parameters, _initC, [=](bool multi, unsigned int comp) { return makeParamId(hashString("INIT_C"), unitOpIdx, comp, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep); });

		//if (_singleBinding)
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

		//	if (_singleBinding)
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
				if (_singleBinding)
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
			if (_singleDynReaction)
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

		return parSurfDiffDepConfSuccess && bindingConfSuccess && dynReactionConfSuccess;
	}

	void ParticleDiffusionOperatorDG::setEquidistantRadialDisc()
	{
		active* const ptrCenterRadius = _parCenterRadius.data();
		active* const ptrOuterSurfAreaPerVolume = _parOuterSurfAreaPerVolume.data();
		active* const ptrInnerSurfAreaPerVolume = _parInnerSurfAreaPerVolume.data();

		const active radius = _parRadius - _parCoreRadius;
		const active _Dr = radius / static_cast<double>(_nParElem);
		std::fill(_parElementSize.data(), _parElementSize.data() + _nParElem, _Dr);

		if (_parGeomSurfToVol == _SurfVolRatioSphere)
		{
			for (unsigned int elem = 0; elem < _nParElem; ++elem)
			{
				const active r_out = _parRadius - static_cast<double>(elem) * _Dr;
				const active r_in = _parRadius - static_cast<double>(elem + 1) * _Dr;

				ptrCenterRadius[elem] = _parRadius - (0.5 + static_cast<double>(elem)) * _Dr;

				// Compute denominator -> corresponding to elem volume
				const active vol = pow(r_out, 3.0) - pow(r_in, 3.0);

				ptrOuterSurfAreaPerVolume[elem] = 3.0 * sqr(r_out) / vol;
				ptrInnerSurfAreaPerVolume[elem] = 3.0 * sqr(r_in) / vol;
			}
		}
		else if (_parGeomSurfToVol == _SurfVolRatioCylinder)
		{
			for (unsigned int elem = 0; elem < _nParElem; ++elem)
			{
				const active r_out = _parRadius - static_cast<double>(elem) * _Dr;
				const active r_in = _parRadius - static_cast<double>(elem + 1) * _Dr;

				ptrCenterRadius[elem] = _parRadius - (0.5 + static_cast<double>(elem)) * _Dr;

				// Compute denominator -> corresponding to elem volume
				const active vol = sqr(r_out) - sqr(r_in);

				ptrOuterSurfAreaPerVolume[elem] = 2.0 * r_out / vol;
				ptrInnerSurfAreaPerVolume[elem] = 2.0 * r_in / vol;
			}
		}
		else if (_parGeomSurfToVol == _SurfVolRatioSlab)
		{
			for (unsigned int elem = 0; elem < _nParElem; ++elem)
			{
				const active r_out = _parRadius - static_cast<double>(elem) * _Dr;
				const active r_in = _parRadius - static_cast<double>(elem + 1) * _Dr;

				ptrCenterRadius[elem] = _parRadius - (0.5 + static_cast<double>(elem)) * _Dr;

				// Compute denominator -> corresponding to elem volume
				const active vol = r_out - r_in;

				ptrOuterSurfAreaPerVolume[elem] = 1.0 / vol;
				ptrInnerSurfAreaPerVolume[elem] = 1.0 / vol;
			}
		}
	}
	/**
	 * @brief Computes the radial nodes in the beads in such a way that all elements have the same volume
	 */
	void ParticleDiffusionOperatorDG::setEquivolumeRadialDisc()
	{
		active* const ptrElemSize = _parElementSize.data();
		active* const ptrCenterRadius = _parCenterRadius.data();
		active* const ptrOuterSurfAreaPerVolume = _parOuterSurfAreaPerVolume.data();
		active* const ptrInnerSurfAreaPerVolume = _parInnerSurfAreaPerVolume.data();

		if (_parGeomSurfToVol == _SurfVolRatioSphere)
		{
			active r_out = _parRadius;
			active r_in = _parCoreRadius;
			const active volumePerelement = (pow(_parRadius, 3.0) - pow(_parCoreRadius, 3.0)) / static_cast<double>(_nParElem);

			for (unsigned int elem = 0; elem < _nParElem; ++elem)
			{
				if (elem != (_nParElem - 1))
					r_in = pow(pow(r_out, 3.0) - volumePerelement, (1.0 / 3.0));
				else
					r_in = _parCoreRadius;

				ptrElemSize[elem] = r_out - r_in;
				ptrCenterRadius[elem] = (r_out + r_in) * 0.5;

				ptrOuterSurfAreaPerVolume[elem] = 3.0 * sqr(r_out) / volumePerelement;
				ptrInnerSurfAreaPerVolume[elem] = 3.0 * sqr(r_in) / volumePerelement;
				// Note that the DG particle elements are oppositely ordered compared to the FV particle elements
				_deltaR[_nParElem - (elem + 1)] = r_out - r_in;

				// For the next elem: r_out == r_in of the current elem
				r_out = r_in;
			}
		}
		else if (_parGeomSurfToVol == _SurfVolRatioCylinder)
		{
			active r_out = _parRadius;
			active r_in = _parCoreRadius;
			const active volumePerelement = (sqr(_parRadius) - sqr(_parCoreRadius)) / static_cast<double>(_nParElem);

			for (unsigned int elem = 0; elem < _nParElem; ++elem)
			{
				if (elem != (_nParElem - 1))
					r_in = sqrt(sqr(r_out) - volumePerelement);
				else
					r_in = _parCoreRadius;

				ptrElemSize[elem] = r_out - r_in;
				ptrCenterRadius[elem] = (r_out + r_in) * 0.5;

				ptrOuterSurfAreaPerVolume[elem] = 2.0 * r_out / volumePerelement;
				ptrInnerSurfAreaPerVolume[elem] = 2.0 * r_in / volumePerelement;
				// Note that the DG particle elements are oppositely ordered compared to the FV particle elements
				_deltaR[_nParElem - (elem + 1)] = r_out - r_in;

				// For the next elem: r_out == r_in of the current elem
				r_out = r_in;
			}
		}
		else if (_parGeomSurfToVol == _SurfVolRatioSlab)
		{
			active r_out = _parRadius;
			active r_in = _parCoreRadius;
			const active volumePerelement = (_parRadius - _parCoreRadius) / static_cast<double>(_nParElem);

			for (unsigned int elem = 0; elem < _nParElem; ++elem)
			{
				if (elem != (_nParElem - 1))
					r_in = r_out - volumePerelement;
				else
					r_in = _parCoreRadius;

				ptrElemSize[elem] = r_out - r_in;
				ptrCenterRadius[elem] = (r_out + r_in) * 0.5;

				ptrOuterSurfAreaPerVolume[elem] = 1.0 / volumePerelement;
				ptrInnerSurfAreaPerVolume[elem] = 1.0 / volumePerelement;
				// Note that the DG particle elements are oppositely ordered compared to the FV particle elements
				_deltaR[_nParElem - (elem + 1)] = r_out - r_in;

				// For the next elem: r_out == r_in of the current elem
				r_out = r_in;
			}
		}
	}

	/**
	 * @brief Computes all helper quantities for radial bead discretization from given radial elem boundaries
	 * @details Calculates surface areas per volume for every element and the radial element centers.
	 */
	void ParticleDiffusionOperatorDG::setUserdefinedRadialDisc()
	{
		active* const ptrElemSize = _parElementSize.data();
		active* const ptrCenterRadius = _parCenterRadius.data();
		active* const ptrOuterSurfAreaPerVolume = _parOuterSurfAreaPerVolume.data();
		active* const ptrInnerSurfAreaPerVolume = _parInnerSurfAreaPerVolume.data();

		// Care for the right ordering and include 0.0 / 1.0 if not already in the vector.
		std::vector<active> orderedInterfaces = std::vector<active>(_parDiscVector.begin(), _parDiscVector.begin() + _nParElem + 1);

		// Sort in descending order
		std::sort(orderedInterfaces.begin(), orderedInterfaces.end(), std::greater<active>());

		// Force first and last element to be 1.0 and 0.0, respectively
		orderedInterfaces[0] = 1.0;
		orderedInterfaces.back() = 0.0;

		// Map [0, 1] -> [core radius, particle radius] via linear interpolation
		for (unsigned int elem = 0; elem < _nParElem; ++elem)
			orderedInterfaces[elem] = static_cast<double>(orderedInterfaces[elem]) * (_parRadius - _parCoreRadius) + _parCoreRadius;

		if (_parGeomSurfToVol == _SurfVolRatioSphere)
		{
			for (unsigned int elem = 0; elem < _nParElem; ++elem)
			{
				ptrElemSize[elem] = orderedInterfaces[elem] - orderedInterfaces[elem + 1];
				ptrCenterRadius[elem] = (orderedInterfaces[elem] + orderedInterfaces[elem + 1]) * 0.5;

				// Compute denominator -> corresponding to elem volume
				const active vol = pow(orderedInterfaces[elem], 3.0) - pow(orderedInterfaces[elem + 1], 3.0);

				ptrOuterSurfAreaPerVolume[elem] = 3.0 * sqr(orderedInterfaces[elem]) / vol;
				ptrInnerSurfAreaPerVolume[elem] = 3.0 * sqr(orderedInterfaces[elem + 1]) / vol;
				// Note that the DG particle elements are oppositely ordered compared to the FV particle elements
				_deltaR[_nParElem - (elem + 1)] = ptrOuterSurfAreaPerVolume[elem] - ptrInnerSurfAreaPerVolume[elem];
			}
		}
		else if (_parGeomSurfToVol == _SurfVolRatioCylinder)
		{
			for (unsigned int elem = 0; elem < _nParElem; ++elem)
			{
				ptrElemSize[elem] = orderedInterfaces[elem] - orderedInterfaces[elem + 1];
				ptrCenterRadius[elem] = (orderedInterfaces[elem] + orderedInterfaces[elem + 1]) * 0.5;

				// Compute denominator -> corresponding to elem volume
				const active vol = sqr(orderedInterfaces[elem]) - sqr(orderedInterfaces[elem + 1]);

				ptrOuterSurfAreaPerVolume[elem] = 2.0 * orderedInterfaces[elem] / vol;
				ptrInnerSurfAreaPerVolume[elem] = 2.0 * orderedInterfaces[elem + 1] / vol;
				// Note that the DG particle elements are oppositely ordered compared to the FV particle elements
				_deltaR[_nParElem - (elem + 1)] = ptrOuterSurfAreaPerVolume[elem] - ptrInnerSurfAreaPerVolume[elem];
			}
		}
		else if (_parGeomSurfToVol == _SurfVolRatioSlab)
		{
			for (unsigned int elem = 0; elem < _nParElem; ++elem)
			{
				ptrElemSize[elem] = orderedInterfaces[elem] - orderedInterfaces[elem + 1];
				ptrCenterRadius[elem] = (orderedInterfaces[elem] + orderedInterfaces[elem + 1]) * 0.5;

				// Compute denominator -> corresponding to elem volume
				const active vol = orderedInterfaces[elem] - orderedInterfaces[elem + 1];

				ptrOuterSurfAreaPerVolume[elem] = 1.0 / vol;
				ptrInnerSurfAreaPerVolume[elem] = 1.0 / vol;
				// Note that the DG particle elements are oppositely ordered compared to the FV particle elements
				_deltaR[_nParElem - (elem + 1)] = ptrOuterSurfAreaPerVolume[elem] - ptrInnerSurfAreaPerVolume[elem];
			}
		}
	}

	// todo: parameter sensitivities for particle radius. Here, we have the problem that every DG operator becomes an active type.
	// We have sums of three matrices each multiplied by its own metric term, which could be applied iteratively to the residual/solution.
	// Not needed for Slab, and only two matrices required for Cylinder.
	// This approach should only be used when necessary, i.e. solely when particle radius parameter sensitivity is required.
	void ParticleDiffusionOperatorDG::updateRadialDisc()
	{
		if (_parDiscMode == ParticleDiscretizationMode::Equidistant)
		{
			for (int elem = 0; elem < _nParElem; elem++)
			{
				_deltaR[elem] = (_parRadius - _parCoreRadius) / _nParElem;
			}
			setEquidistantRadialDisc();
		}
		else if (_parDiscMode == ParticleDiscretizationMode::Equivolume)
			setEquivolumeRadialDisc();
		else if (_parDiscMode == ParticleDiscretizationMode::UserDefined)
			setUserdefinedRadialDisc();

		/*		metrics		*/
		// estimate element dependent operators

		for (int elem = 0; elem < _nParElem; elem++)
		{
			for (int node = 0; node < _nParNode; node++)
				_Ir[elem][node] = _deltaR[elem] / 2.0 * (_parNodes[node] + 1.0);

			active r_L = _parCoreRadius + elem * _deltaR[elem]; // left boundary of current elem

			_Ir[elem] = _Ir[elem] + VectorXd::Ones(_nParNode) * r_L;

			if (_parGeomSurfToVol == _SurfVolRatioSphere)
				_Ir[elem] = _Ir[elem].array().square();
			else if (_parGeomSurfToVol == _SurfVolRatioSlab)
				_Ir[elem].setOnes(); // no metric terms for slab

			// compute mass matrices for exact integration based on particle geometry, via transformation to normalized Jacobi polynomials with weight function w
			if (_parGeomSurfToVol == _SurfVolRatioSphere) // r^2 =  r_i^2 + (1 + \xi) * r_i * DeltaR_i / 2.0 + (1 + \xi)^2 * (DeltaR_i / 2.0)^2
			{
				_parInvMM[elem] = parts::dgtoolbox::mMatrix(_parPolyDeg, _parNodes, 0.0, 2.0) * pow((static_cast<double>(_deltaR[elem]) / 2.0), 2.0);
				if (elem > 0 || _parCoreRadius != 0.0) // following contributions are zero for f_Irst elem when R_c = 0 (no particle core)
					_parInvMM[elem] += parts::dgtoolbox::mMatrix(_parPolyDeg, _parNodes, 0.0, 1.0) * (static_cast<double>(_deltaR[elem]) * static_cast<double>(r_L))
					+ parts::dgtoolbox::mMatrix(_parPolyDeg, _parNodes, 0.0, 0.0) * pow(static_cast<double>(r_L), 2.0);

				_parInvMM[elem] = _parInvMM[elem].inverse();
				_minus_InvMM_ST[elem] = -_parInvMM[elem] * _parPolyDerM.transpose() * _parInvMM[elem].inverse();

				// particle GSM specific second order stiffness matrix (single element, i.e. _nParElem = 1)
				_secondOrderStiffnessM = std::pow(static_cast<double>(_parCoreRadius), 2.0) * parts::dgtoolbox::secondOrderStiffnessMatrix(_parPolyDeg, 0.0, 0.0, _parNodes);
				_secondOrderStiffnessM += static_cast<double>(_deltaR[0]) * static_cast<double>(_parCoreRadius) * parts::dgtoolbox::secondOrderStiffnessMatrix(_parPolyDeg, 0.0, 1.0, _parNodes);
				_secondOrderStiffnessM += std::pow(static_cast<double>(_deltaR[0]) / 2.0, 2.0) * parts::dgtoolbox::secondOrderStiffnessMatrix(_parPolyDeg, 0.0, 2.0, _parNodes);
			}
			else if (_parGeomSurfToVol == _SurfVolRatioCylinder) // r = r_i + (1 + \xi) * DeltaR_i / 2.0
			{
				_parInvMM[elem] = parts::dgtoolbox::mMatrix(_parPolyDeg, _parNodes, 0.0, 1.0) * (static_cast<double>(_deltaR[elem]) / 2.0);
				if (elem > 0 || _parCoreRadius != 0.0) // following contribution is zero for f_Irst elem when R_c = 0 (no particle core)
					_parInvMM[elem] += parts::dgtoolbox::mMatrix(_parPolyDeg, _parNodes, 0.0, 0.0) * static_cast<double>(r_L);

				_parInvMM[elem] = _parInvMM[elem].inverse();
				_minus_InvMM_ST[elem] = -_parInvMM[elem] * _parPolyDerM.transpose() * _parInvMM[elem].inverse();

				// particle GSM specific second order stiffness matrix (single element, i.e. _nParElem = 1)
				_secondOrderStiffnessM = static_cast<double>(_parCoreRadius) * parts::dgtoolbox::secondOrderStiffnessMatrix(_parPolyDeg, 0.0, 0.0, _parNodes);
				_secondOrderStiffnessM += static_cast<double>(_deltaR[0]) / 2.0 * parts::dgtoolbox::secondOrderStiffnessMatrix(_parPolyDeg, 0.0, 1.0, _parNodes);
			}
			else if (_parGeomSurfToVol == _SurfVolRatioSlab) // r = 1
			{
				_minus_InvMM_ST[elem] = -_parInvMM[elem] * _parPolyDerM.transpose() * _parInvMM[elem].inverse();

				_secondOrderStiffnessM = parts::dgtoolbox::secondOrderStiffnessMatrix(_parPolyDeg, 0.0, 0.0, _parNodes);
				_parInvMM[elem] = parts::dgtoolbox::invMMatrix(_parPolyDeg, _parNodes, 0.0, 0.0);
			}
		}

		_minus_parInvMM_Ar = -_parInvMM[0] * _secondOrderStiffnessM;
	}

	/**
	 * @brief allocates memory for DG operators and computes those that are metric independent. Also allocates required containers needed for the DG discretization.
	 */
	void ParticleDiffusionOperatorDG::initializeDG()
	{
		/* Allocate space for DG operators and containers */

		const bool firstConfigCall = _parInvMM == nullptr; // used to not multiply allocate memory

		// particles
		if (firstConfigCall)
		{
			_localFlux = new active[_nComp];
		}

		_nParNode = _parPolyDeg + 1u;
		_nParPoints = _nParNode * _nParElem;
		_g_p.resize(_nParPoints);
		_g_p.setZero();
		_g_pSum.resize(_nParPoints);
		_g_pSum.setZero();
		_surfaceFluxParticle.resize(_nParElem + 1);
		_surfaceFluxParticle.setZero();
		_parNodes.resize(_nParNode);
		_parNodes.setZero();
		_parInvWeights.resize(_nParNode);
		_parInvWeights.setZero();
		_parPolyDerM.resize(_nParNode, _nParNode);
		_parPolyDerM.setZero();
		_parInvMM_Leg.resize(_nParNode, _nParNode);
		_parInvMM_Leg.setZero();

		if (firstConfigCall)
		{
			_Ir = new Vector<active, Dynamic>[_nParElem];
			for (int elem = 0; elem < _nParElem; elem++)
			{
				_Ir[elem].resize(_nParNode);
				_Ir[elem].setZero();
			}
			_minus_InvMM_ST = new MatrixXd[_nParElem];
			_parInvMM = new MatrixXd[_nParElem];
		}

		/* compute metric independent DG operators for bulk and particles. Note that metric dependent DG operators are computet in updateRadialDisc(). */

		parts::dgtoolbox::lglNodesWeights(_parPolyDeg, _parNodes, _parInvWeights, true);
		_parPolyDerM = parts::dgtoolbox::derivativeMatrix(_parPolyDeg, _parNodes);
		_parInvMM_Leg = parts::dgtoolbox::invMMatrix(_parPolyDeg, _parNodes, 0.0, 0.0);
	}


	/**
	 * @brief Notifies the operator that a discontinuous section transition is in progress
	 * @param [in] t Current time point
	 * @param [in] secIdx Index of the new section that is about to be integrated
	 * @return @c true if flow direction has changed, otherwise @c false
	 */
	bool ParticleDiffusionOperatorDG::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active const* const filmDiff, active const* const poreAccessFactor)
	{
			for (int comp = 0; comp < _nComp; comp++)
			{
				_filmDiffusion[comp] = filmDiff[_parTypeIdx * _nComp + comp];
				_poreAccessFactor[comp] = poreAccessFactor[_parTypeIdx * _nComp + comp];
				_invBetaP[comp] = (1.0 - _parPorosity) / (_poreAccessFactor[comp] * _parPorosity);
			}

		//_curSection = secIdx;
		//_newStaticJac = true;

		// todo update operators

		initializeDGjac(_parGeomSurfToVol);

		return true;
	}
	/**
	 * @brief calculates the physical radial/particle coordinates of the DG discretization with double! interface nodes
	 */
	int ParticleDiffusionOperatorDG::getParticleCoordinates(double* coords) const
	{
		active const* const pcr = _parCenterRadius.data();

		// Note that the DG particle shells are oppositely ordered compared to the FV particle shells
		for (unsigned int par = 0; par < _nParPoints; par++) {

			unsigned int cell = std::floor(par / _nParNode);

			double r_L = static_cast<double>(pcr[cell]) - 0.5 * static_cast<double>(_deltaR[cell]);
			coords[par] = r_L + 0.5 * static_cast<double>(_deltaR[cell]) * (1.0 + _parNodes[par % _nParNode]);
		}

		return _nParPoints;
	}
	/**
	 * @brief Computes the residual of the transport equations
	 * @param [in] model Model that owns the operator
	 * @param [in] t Current time point
	 * @param [in] secIdx Index of the current section
	 * @param [in] yPar Pointer to particle phase entry in unit state vector
	 * @param [in] yBulk Pointer to corresponding bulk phase entry in unit state vector
	 * @param [in] yDotPar Pointer to particle phase derivative entry in unit state vector
	 * @param [out] resPar Pointer Pointer to particle phase entry in unit residual vector
	 * @param [out] colPos column position of the particle (particle coordinate zero)
	 * @param [in] jacIt Matrix iterator pointing to the particle phase entry in the unit Jacobian
	 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
	 */
	template<bool wantJac, bool wantRes>
	int ParticleDiffusionOperatorDG::residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity)
	{
		return residualImpl<double, double, double, wantJac, wantRes>(t, secIdx, yPar, yBulk, yDotPar, resPar, colPos, jacIt, tlmAlloc);
	}

	template<bool wantJac, bool wantRes>
	int ParticleDiffusionOperatorDG::residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity)
	{
		return residualImpl<active, active, double, wantJac, wantRes>(t, secIdx, yPar, yBulk, yDotPar, resPar, colPos, jacIt, tlmAlloc);
	}

	template<bool wantJac, bool wantRes>
	int ParticleDiffusionOperatorDG::residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity)
	{
		return residualImpl<double, active, active, wantJac, wantRes>(t, secIdx, yPar, yBulk, yDotPar, resPar, colPos, jacIt, tlmAlloc);
	}

	template<bool wantJac, bool wantRes>
	int ParticleDiffusionOperatorDG::residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity)
	{
		return residualImpl<active, active, active, wantJac, wantRes>(t, secIdx, yPar, yBulk, yDotPar, resPar, colPos, jacIt, tlmAlloc);
	}

	template int ParticleDiffusionOperatorDG::residual<true, true>(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity);
	template int ParticleDiffusionOperatorDG::residual<false, true>(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity);
	template int ParticleDiffusionOperatorDG::residual<true, false>(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity);

	template int ParticleDiffusionOperatorDG::residual<true, true>(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity);
	template int ParticleDiffusionOperatorDG::residual<false, true>(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity);

	template int ParticleDiffusionOperatorDG::residual<true, true>(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity);
	template int ParticleDiffusionOperatorDG::residual<false, true>(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity);

	template int ParticleDiffusionOperatorDG::residual<true, true>(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity);
	template int ParticleDiffusionOperatorDG::residual<false, true>(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity);

	/**
	 * @brief promotes doubles to actives
	 * @detail promotes consecutive doubles to consecutive actives (with zero gradients) based on input double pointer
	 */
	void vectorPromoter(double* state, const unsigned int nVals)
	{
		const int nDirs = ad::getDirections();
		const int stride = (1 + nDirs);
		const int ADsize = stride * nVals;
		double buff = 0.0;

		for (int val = 1; val <= nVals; val++) // start with last entry to avoid overwriting
		{
			buff = state[nVals - val];
			std::fill(state + ADsize - val * stride, state + ADsize - (val - 1) * stride, 0.0);
			state[ADsize - val * stride] = buff;
		}
	}

	cell::CellParameters ParticleDiffusionOperatorDG::makeCellResidualParams(int const* qsReaction) const
	{
		return cell::CellParameters
		{
			_nComp,
			_nBound,
			_boundOffset,
			_strideBound,
			qsReaction,
			_parPorosity,
			_poreAccessFactor.data(),
			_binding,
			(_dynReaction && (_dynReaction->numReactionsCombined() > 0)) ? _dynReaction : nullptr
		};
	}

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
	int ParticleDiffusionOperatorDG::particleDiffusionImpl(double t, unsigned int secIdx, StateType const* yPar, StateType const* yBulk, double const* yDotPar, ResidualType* resPar, const int* const qsBinding, const bool hasDynamicReactions, linalg::BandedEigenSparseRowIterator& jacBase)
	{
		// Add the DG discretized solid entries of the jacobian that get overwritten by the binding kernel.
		// These entries only exist for the GRM with surface diffusion
		if (wantJac && hasDynamicReactions && _hasSurfaceDiffusion)
			addSolidDGentries(secIdx, jacBase, qsBinding);

		if (!wantRes)
			return 0;

		/* Mobile phase RHS	*/

		// Get film diffusion flux at current node to compute boundary condition
		for (unsigned int comp = 0; comp < _nComp; comp++) {
			_localFlux[comp] = _filmDiffusion[comp] * (yBulk[comp * _strideBulkComp] - yPar[(_nParPoints - 1) * strideParNode() + comp]);
		}

		active const* const parDiff = getSectionDependentSlice(_parDiffusion, _nComp, secIdx);

		// Ordering of particle surface diffusion:
		// bnd0comp0, bnd1comp0, bnd0comp1, bnd1comp1, bnd0comp2, bnd1comp2
		active const* const parSurfDiff = getSectionDependentSlice(_parSurfDiffusion, _strideBound, secIdx);

		const int nNodes = _nParNode;
		const int nElem = _nParElem;
		const int nPoints = _nParPoints;
		const int nComp = _nComp;

		if (_parGSM) // GSM implementation
		{
			for (unsigned int comp = 0; comp < nComp; comp++)
			{

				// ====================================================================================//
				// solve GSM-discretized particle mass balance   									   //
				// ====================================================================================//

				Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> resCp(resPar + comp, nPoints, InnerStride<Dynamic>(strideParNode()));
				Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> Cp(yPar + comp, nPoints, InnerStride<Dynamic>(strideParNode()));

				// Use auxiliary variable to get c^p + \sum 1 / \Beta_p c^s
				Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<>> sum_cp_cs(reinterpret_cast<ResidualType*>(&_g_pSum[0]), nPoints, InnerStride<>(1));
				sum_cp_cs = static_cast<ParamType>(parDiff[comp]) * Cp.template cast<ResidualType>();

				for (int bnd = 0; bnd < _nBound[comp]; bnd++)
				{
					if (parSurfDiff[_boundOffset[comp] + bnd] != 0.0) // some bound states might still not be effected by surface diffusion
					{
						Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> c_s(yPar + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd, nPoints, InnerStride<Dynamic>(strideParNode()));
						sum_cp_cs += static_cast<ParamType>(_invBetaP[comp]) * static_cast<ParamType>(parSurfDiff[_boundOffset[comp] + bnd]) * c_s;

						/* For kinetic bindings with surface diffusion: add the additional DG-discretized particle mass balance equations to residual */

						if (!qsBinding[bnd])
						{
							// Eigen access to current bound state residual
							Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> resCs(resPar + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd,
								nPoints, InnerStride<Dynamic>(strideParNode()));

							// Use auxiliary variable to get \Beta_p D_s c^s
							Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<>> c_s_modified(reinterpret_cast<ResidualType*>(&_g_p[0]), nPoints, InnerStride<>(1));

							// Apply squared inverse mapping and surface diffusion
							c_s_modified = 2.0 / static_cast<ParamType>(_deltaR[0]) * 2.0 / static_cast<ParamType>(_deltaR[0]) *
								static_cast<ParamType>(parSurfDiff[_boundOffset[comp] + bnd]) * c_s;

							Eigen::Map<const Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> c_s_modified_const(&c_s_modified[0], nPoints, InnerStride<Dynamic>(1));
							parGSMVolumeIntegral<ResidualType, ResidualType>(c_s_modified_const, resCs);

							// Leave out the surface integral as we only have one element, i.e. we apply BC with zeros
						}
					}
				}

				// Apply squared inverse mapping
				sum_cp_cs *= 2.0 / static_cast<ParamType>(_deltaR[0]) * 2.0 / static_cast<ParamType>(_deltaR[0]);

				Eigen::Map<const Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> sum_cp_cs_const(&sum_cp_cs[0], nPoints, InnerStride<Dynamic>(1));
				parGSMVolumeIntegral<ResidualType, ResidualType>(sum_cp_cs_const, resCp);

				// Pass sum_cp_cs_const to match the DGSEM interface; nullptr might also be feasible
				parSurfaceIntegral<ResidualType, ResidualType>(sum_cp_cs_const, resCp, nNodes, 1u, false, comp);

			}
		}
		else // DGSEM implementation
		{
			for (unsigned int comp = 0; comp < nComp; comp++)
			{
				// =====================================================================================================//
				// Solve auxiliary systems  d_p g_p + d_s beta_p sum g_s= d (d_p c_p + d_s beta_p sum c_s) / d xi		//
				// =====================================================================================================//

				// Component-wise! strides
				unsigned int strideelem = nNodes;
				unsigned int strideNode = 1u;

				// Reset cache for auxiliary variable
				Eigen::Map<Vector<StateType, Dynamic>, 0, InnerStride<>> g_p(reinterpret_cast<StateType*>(&_g_p[0]), nPoints, InnerStride<>(1));
				Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<>> g_pSum(reinterpret_cast<ResidualType*>(&_g_pSum[0]), nPoints, InnerStride<>(1));
				g_p.setZero();
				g_pSum.setZero();

				Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> cp(yPar + comp, _nParPoints, InnerStride<Dynamic>(strideParNode()));

				// Handle surface diffusion: Compute auxiliary variable; For kinetic bindings: add additional mass balance to residual of respective bound state
				if (_hasSurfaceDiffusion)
				{
					for (int bnd = 0; bnd < _nBound[comp]; bnd++)
					{
						if (parSurfDiff[_boundOffset[comp] + bnd] != 0.0) // some bound states might still not be effected by surface diffusion
						{
							// Get solid phase vector
							Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> c_s(yPar + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd,
								nPoints, InnerStride<Dynamic>(strideParNode()));
							// Compute g_s = d c_s / d xi
							solve_auxiliary_DG<StateType>(c_s, strideelem, strideNode, comp);
							// Apply invBeta_p, d_s and add to sum -> gSum += d_s * invBeta_p * (D c - M^-1 B [c - c^*])
							g_pSum += g_p.template cast<ResidualType>() * static_cast<ParamType>(_invBetaP[comp]) * static_cast<ParamType>(parSurfDiff[_boundOffset[comp] + bnd]);

							/* For kinetic bindings with surface diffusion: add the additional DG-discretized particle mass balance equations to residual */

							if (!qsBinding[bnd])
							{
								// Eigen access to current bound state residual
								Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> resCs(resPar + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd,
									nPoints, InnerStride<Dynamic>(strideParNode()));

								// Promote auxiliary variable storage from double to active if required
								// @todo is there a more efficient or elegant solution?
								if (std::is_same<ResidualType, active>::value && std::is_same<StateType, double>::value)
									vectorPromoter(reinterpret_cast<double*>(&g_p[0]), nPoints); // reinterpret_cast only required because statement is scanned also when StateType != double

								// Access auxiliary variable as ResidualType
								Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<>> g_p_ResType(reinterpret_cast<ResidualType*>(&_g_p[0]), nPoints, InnerStride<>(1));

								applyParInvMap<ResidualType, ParamType>(g_p_ResType);
								g_p_ResType *= static_cast<ParamType>(parSurfDiff[_boundOffset[comp] + bnd]);

								// Eigen access to auxiliary variable of current bound state
								Eigen::Map<const Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> g_p_ResType_const(&g_p_ResType[0], nPoints, InnerStride<Dynamic>(1));

								// Add - D_r * gs to the residual, including metric part.->res = invMap^2* [ -D_r * (d_s c^s) ]
								parVolumeIntegral<ResidualType, ResidualType>(false, g_p_ResType_const, resCs);

								// Add M^-1 B (gs - gs^*) to the residual -> res =  invMap^2 * [ - D_r * (d_s c^s) + M^-1 B (gs - gs^*) ]
								parSurfaceIntegral<ResidualType, ResidualType>(g_p_ResType_const, resCs, strideelem, strideNode, false, comp, true);
							}
						}
					}
				}

				// Compute g_p = d c_p / d xi
				solve_auxiliary_DG<StateType>(cp, strideelem, strideNode, comp);

				// Add particle diffusion part to auxiliary variable sum -> gSum += d_p * (D c - M^-1 B [c - c^*])
				g_pSum += g_p * static_cast<ParamType>(parDiff[comp]);

				// apply squared inverse mapping to sum of bound state auxiliary variables -> gSum = - invMap^2 * (d_p * c^p + sum_mi d_s invBeta_p c^s)
				applyParInvMap<ResidualType, ParamType>(g_pSum);

				// ====================================================================================//
				// solve DG-discretized particle mass balance   									     //
				// ====================================================================================//

				  /* Solve DG-discretized particle mass balance equation */

				Eigen::Map<const Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> g_pSum_const(&g_pSum[0], nPoints, InnerStride<Dynamic>(1));

				// Eigen access to particle liquid residual
				Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> resCp(resPar + comp, nPoints, InnerStride<Dynamic>(strideParNode()));

				// Add - D_r * (g_sum) to the residual, including metric part. -> res = - D_r * (d_p * c^p + invBeta_p sum_mi d_s c^s)
				parVolumeIntegral<ResidualType, ResidualType>(false, g_pSum_const, resCp);

				// Add M^-1 B (g_sum - g_sum^*) to the residual -> res = - D_r * (d_p * c^p + invBeta_p sum_mi d_s c^s) + M^-1 B (g_sum - g_sum^*)
				parSurfaceIntegral<ResidualType, ResidualType>(g_pSum_const, resCp, strideelem, strideNode, false, comp);

			}
		}

		return true;
	}


	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
	int ParticleDiffusionOperatorDG::residualImpl(double t, unsigned int secIdx, StateType const* yPar, StateType const* yBulk, double const* yDotPar, ResidualType* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc)
	{
		int const* const qsBinding = _binding->reactionQuasiStationarity();
		const cell::CellParameters cellResParams = makeCellResidualParams(qsBinding);

		linalg::BandedEigenSparseRowIterator jacBase = jacIt;

		// Handle time derivatives, binding, dynamic reactions: residualKernel computes discrete point wise,
		// so we loop over each discrete particle point
		for (unsigned int par = 0; par < _nParPoints; ++par)
		{
			// local pointers to current particle node, needed in residualKernel
			StateType const* local_y = yPar + par * strideParNode();
			double const* local_yDot = yDotPar ? yDotPar + par * strideParNode() : nullptr;
			ResidualType* local_res = resPar ? resPar + par * strideParNode() : nullptr;

			// r (particle) coordinate of current node (particle radius normed to 1) - needed in externally dependent adsorption kinetic
			colPos.particle = relativeCoordinate(par);

			if (wantRes)
				cell::residualKernel<StateType, ResidualType, ParamType, cell::CellParameters, linalg::BandedEigenSparseRowIterator, wantJac, true>(
					t, secIdx, colPos, local_y, local_yDot, local_res, jacIt, cellResParams, tlmAlloc
				);
			else
				cell::residualKernel<StateType, ResidualType, ParamType, cell::CellParameters, linalg::BandedEigenSparseRowIterator, wantJac, false, false>(
					t, secIdx, colPos, local_y, local_yDot, local_res, jacIt, cellResParams, tlmAlloc
				);

			// Move rowiterator to next particle node
			jacIt += strideParNode();
		}

		particleDiffusionImpl<StateType, ResidualType, ParamType, wantJac, wantRes>(t, secIdx, yPar, yBulk, yDotPar, resPar, qsBinding, _binding->hasDynamicReactions(), jacBase);

		return true;
	}

	template<typename ResidualType, typename ParamType>
	void ParticleDiffusionOperatorDG::applyParInvMap(Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<>>& state)
	{
		for (int elem = 0; elem < _nParElem; elem++) {
			state.segment(elem * _nParNode, _nParNode) *= 2.0 / static_cast<ParamType>(_deltaR[elem]) * 2.0 / static_cast<ParamType>(_deltaR[elem]);
		}
	}

	template<typename StateType, typename ResidualType>
	void ParticleDiffusionOperatorDG::parGSMVolumeIntegral(Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state, Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>>& stateDer)
	{
		int nNodes = _nParNode;

		stateDer.segment(0, nNodes)
			-= (_minus_parInvMM_Ar.template cast<StateType>() * state.segment(0, nNodes)).template cast<ResidualType>();
	}

	template<typename StateType, typename ResidualType>
	void ParticleDiffusionOperatorDG::parVolumeIntegral(const bool aux, Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state, Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>>& stateDer)
	{
		int nNodes = _nParNode;

		/* no additional metric term for auxiliary equation or particle equation with exact integration scheme
		   -> res = - D * (d_p * c^p + invBeta_p sum_mi d_s c^s) */
		if (aux || (_parGeomSurfToVol == _SurfVolRatioSlab)) {
			// comp-elem-node state vector: use of Eigen lib performance
			for (unsigned int elem = 0; elem < _nParElem; elem++) {
				stateDer.segment(elem * nNodes, nNodes)
					-= (_parPolyDerM.template cast<StateType>() * state.segment(elem * nNodes, nNodes)).template cast<ResidualType>();
			}
		}
		else if (_parGeomSurfToVol != _SurfVolRatioSlab) {
			// comp-elem-node state vector: use of Eigen lib performance
			for (unsigned int elem = 0; elem < _nParElem; elem++) {
				stateDer.segment(elem * nNodes, nNodes)
					-= (_minus_InvMM_ST[elem].template cast<StateType>() * state.segment(elem * nNodes, nNodes)).template cast<ResidualType>();
			}
		}
	}
	/*
	 * @brief calculates the interface fluxes g* of particle mass balance equation and implements the respective boundary conditions
	 * @param [in] aux bool if interface flux for auxiliary equation
	 * @param [in] addParDisc bool if interface flux for additional particle DG-discretized equation
	*/
	template<typename StateType>
	void ParticleDiffusionOperatorDG::InterfaceFluxParticle(Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state,
		const unsigned int strideelem, const unsigned int strideNode, const bool aux, const int comp, const bool addParDisc)
	{
		Eigen::Map<Vector<StateType, Dynamic>, 0, InnerStride<>> _surfFluxPar(reinterpret_cast<StateType*>(&_surfaceFluxParticle[0]), _nParElem + 1, InnerStride<>(1));

		// reset surface flux storage as it is used multiple times
		_surfFluxPar.setZero();

		// numerical flux: state* = 0.5 (state^+ + state^-)

		// calculate inner interface fluxes
		for (unsigned int elem = 1u; elem < _nParElem; elem++) {
			_surfFluxPar[elem] // left interfaces
				= 0.5 * (state[elem * strideelem - strideNode] + // outer/left node
					state[elem * strideelem]); // inner/right node
		}

		// calculate boundary interface fluxes.
		if (aux) { // ghost nodes given by state^- := state^+ for auxiliary equation
			_surfFluxPar[0] = state[0];

			_surfFluxPar[_nParElem] = state[_nParElem * strideelem - strideNode];
		}
		else if (addParDisc) {
			_surfFluxPar[0] = 0.0;

			_surfFluxPar[_nParElem] = 0.0;
		}
		else {

			// film diffusion BC
			_surfFluxPar[_nParElem] = static_cast<StateType>(_localFlux[comp])
				/ (static_cast<double>(_parPorosity) * static_cast<double>(_poreAccessFactor[comp]))
				* (2.0 / static_cast<double>(_deltaR[0])); // inverse squared mapping was also applied, so we apply Map * invMap^2 = invMap

			// inner particle BC
			_surfFluxPar[0] = 0.0;

		}
	}
	/**
	 * @brief calculates the particle surface Integral (type- and component-wise)
	 * @param [in] state relevant state vector
	 * @param [in] stateDer state derivative vector the solution is added to
	 * @param [in] aux true for auxiliary equation, false for main equation
	 * @param [in] strideelem component-wise elem stride
	 * @param [in] strideNodecomponent-wise node stride
	 * @param [in] comp current component
	*/
	template<typename StateType, typename ResidualType>
	void ParticleDiffusionOperatorDG::parSurfaceIntegral(Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state,
		Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>>& stateDer, unsigned const int strideelem, unsigned const int strideNode,
		const bool aux, const int comp, const bool addParDisc)
	{
		Eigen::Map<Vector<StateType, Dynamic>, 0, InnerStride<>> _surfFluxPar(reinterpret_cast<StateType*>(&_surfaceFluxParticle[0]), _nParElem + 1, InnerStride<>(1));

		// calc numerical flux values
		InterfaceFluxParticle<StateType>(state, strideelem, strideNode, aux, comp, addParDisc);

		// strong surface integral -> M^-1 B [state - state*]
		for (unsigned int elem = 0; elem < _nParElem; elem++) {

			for (unsigned int Node = 0; Node < _nParNode; Node++) {
				if (aux) { // strong surface integral -> M^-1 B [state - state*] and B has two non-zero entries, -1 and 1
					stateDer[elem * strideelem + Node * strideNode]
						-= _parInvMM_Leg(Node, 0) * (state[elem * strideelem] - _surfFluxPar[elem])
						- _parInvMM_Leg(Node, _parPolyDeg) * (state[elem * strideelem + _parPolyDeg * strideNode] - _surfFluxPar[elem + 1u]);
				}
				else {
					if (_parGeomSurfToVol == _SurfVolRatioSlab) { // strong surface integral -> M^-1 B [state - state*] and B has two non-zero entries, -1 and 1
						stateDer[elem * strideelem + Node * strideNode]
							-= static_cast<ResidualType>(
								_parInvMM[elem](Node, 0) * (state[elem * strideelem] - _surfFluxPar[elem])
								- _parInvMM[elem](Node, _parPolyDeg) * (state[elem * strideelem + _parPolyDeg * strideNode] - _surfFluxPar[elem + 1u])
								);
					}
					else if (_parGeomSurfToVol == _SurfVolRatioCylinder) { // weak surface integral -> M^-1 B [- state*] and B has two non-zero entries, which depend on metrics
						stateDer[elem * strideelem + Node * strideNode]
							-= static_cast<ResidualType>(
								_Ir[elem][0] * _parInvMM[elem](Node, 0) * (-_surfFluxPar[elem])
								+ _Ir[elem][_nParNode - 1] * _parInvMM[elem](Node, _parPolyDeg) * _surfFluxPar[elem + 1u]
								);
					}
					else if (_parGeomSurfToVol == _SurfVolRatioSphere) { // weak surface integral -> M^-1 B [- state*] and B has two non-zero entries, which depend on metrics
						stateDer[elem * strideelem + Node * strideNode]
							-= static_cast<ResidualType>(
								_Ir[elem][0] * _parInvMM[elem](Node, 0) * (-_surfFluxPar[elem])
								+ _Ir[elem][_nParNode - 1] * _parInvMM[elem](Node, _parPolyDeg) * _surfFluxPar[elem + 1u]
								);
					}
				}
			}
		}
	}
	/**
	 * @brief solves the auxiliary system g = d c / d xi
	 * @detail computes g = Dc - M^-1 B [c - c^*] and stores this in _g_p
	*/
	template<typename StateType>
	void ParticleDiffusionOperatorDG::solve_auxiliary_DG(Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<>>& conc, unsigned int strideelem, unsigned int strideNode, int comp)
	{
		Eigen::Map<Vector<StateType, Dynamic>, 0, InnerStride<>> g_p(reinterpret_cast<StateType*>(&_g_p[0]), _nParPoints, InnerStride<>(1));
		Eigen::Map<Vector<StateType, Dynamic>, 0, InnerStride<>> _surfFluxPar(reinterpret_cast<StateType*>(&_surfaceFluxParticle[0]), _nParElem + 1, InnerStride<>(1));
		_surfFluxPar.setZero(); // reset surface flux storage as it is used multiple times
		g_p.setZero(); // reset auxiliary variable g

		// ========================================================================================//
		// solve auxiliary systems g = d c / d xi	 =>		g_p = Dc - M^-1 B [c - c^*]			   //
		// ========================================================================================//

		parVolumeIntegral<StateType, StateType>(true, conc, g_p); // volumne integral in strong DG form: - D c

		parSurfaceIntegral<StateType>(conc, g_p, strideelem, strideNode, true, comp); // surface integral in strong DG form: M^-1 B [c - c^*]

		g_p *= -1.0; // auxiliary factor -1
	}

	// ==========================================================================================================================================================  //
	// ========================================						DG particle Jacobian							=============================================  //
	// ==========================================================================================================================================================  //

	/**
	 * @brief calculates the particle dispersion jacobian Pattern of the DG scheme for the given particle type and bead
	*/
	void ParticleDiffusionOperatorDG::calcParticleJacobianPattern(std::vector<ParticleDiffusionOperatorDG::T>& tripletList, unsigned int offset, unsigned int colNode, unsigned int secIdx)
	{
		// Ordering of particle surface diffusion:
		// bnd0comp0, bnd1comp0, bnd0comp1, bnd1comp1, bnd0comp2, bnd1comp2
		active const* const _parSurfDiff = getSectionDependentSlice(_parSurfDiffusion, _strideBound, secIdx);

		int const* const qsBinding = _binding->reactionQuasiStationarity();

		// (global) strides
		unsigned int selem = _nParNode * strideParNode();
		unsigned int sNode = strideParNode();
		unsigned int sComp = 1u;
		unsigned int nNodes = _nParNode;

		// case: one elem  -> diffBlock \in R^(nParNodes x nParNodes), GBlock = parPolyDerM
		if (_nParElem == 1) {

			// fill the jacobian: add dispersion block for each unbound and bound component, adjusted for the respective coefficients
			for (unsigned int comp = 0; comp < _nComp; comp++) {
				for (unsigned int i = 0; i < nNodes; i++) {
					for (unsigned int j = 0; j < nNodes; j++) {
						// handle liquid state
						// row: add component offset and go node strides from there for each dispersion block entry
						// col: add component offset and go node strides from there for each dispersion block entry
						tripletList.push_back(T(offset + comp * sComp + i * sNode,
							offset + comp * sComp + j * sNode, 0.0));

						// handle surface diffusion of bound states.
						if (_hasSurfaceDiffusion) {

							for (unsigned int bnd = 0; bnd < _nBound[comp]; bnd++) {
								if (_parSurfDiff[_boundOffset[comp] + bnd] != 0.0) {
									// row: add current component offset and go node strides from there for each dispersion block entry
									// col: jump oover liquid states, add current bound state offset and go node strides from there for each dispersion block entry
									tripletList.push_back(T(offset + comp * sComp + i * sNode,
										offset + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd + j * sNode, 0.0));

									/* add surface diffusion dispersion block to solid */
									if (!qsBinding[offsetBoundComp(ComponentIndex{ comp }) + bnd]) {
										// row: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry
										// col: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry
										tripletList.push_back(T(offset + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd + i * sNode,
											offset + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd + j * sNode, 0.0));

									}
								}
							}
						}
					}
				}
			}
		}
		else {
			/*			boundary elements			*/

			/*			 left boundary elem				*/

			unsigned int special = 0u; if (_nParElem < 3u) special = 1u; // limits the iterator for special case nelements = 3 (dependence on additional entry)
			// fill the jacobian: add dispersion block for each unbound and bound component, adjusted for the respective coefficients
			for (unsigned int comp = 0; comp < _nComp; comp++) {
				for (unsigned int i = 0; i < nNodes; i++) {
					for (unsigned int j = nNodes + 1; j < 3 * nNodes + 2 - special; j++) {
						// handle liquid state
						// row: add component offset and go node strides from there for each dispersion block entry
						// col: add component offset and go node strides from there for each dispersion block entry. adjust for j start
						tripletList.push_back(T(offset + comp * sComp + i * sNode,
							offset + comp * sComp + j * sNode - (nNodes + 1) * sNode,
							0.0));

						// handle surface diffusion of bound states. binding is handled in residualKernel().
						if (_hasSurfaceDiffusion) {

							for (unsigned int bnd = 0; bnd < _nBound[comp]; bnd++) {
								if (_parSurfDiff[_boundOffset[comp] + bnd] != 0.0) {
									// row: add current component offset and go node strides from there for each dispersion block entry
									// col: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry. adjust for j start
									tripletList.push_back(T(offset + comp * sComp + i * sNode,
										offset + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd + j * sNode - (nNodes + 1) * sNode,
										0.0));

									/* add surface diffusion dispersion block to solid */
									if (!qsBinding[offsetBoundComp(ComponentIndex{ comp }) + bnd]) {
										// row: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry
										// col: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry. adjust for j start
										tripletList.push_back(T(offset + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd + i * sNode,
											offset + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd + j * sNode - (nNodes + 1) * sNode,
											0.0));

									}
								}
							}
						}
					}
				}
			}

			/*			 right boundary elem				*/

			// fill the jacobian: add dispersion block for each unbound and bound component, adjusted for the respective coefficients
			for (unsigned int comp = 0; comp < _nComp; comp++) {
				for (unsigned int i = 0; i < nNodes; i++) {

					for (unsigned int j = special; j < 2 * nNodes + 1; j++) {
						// handle liquid state
						// row: add component offset and jump over previous elements. Go node strides from there for each dispersion block entry
						// col: add component offset and jump over previous elements. Go back one elem (and node or adjust for start) and go node strides from there for each dispersion block entry.
						tripletList.push_back(T(offset + comp * sComp + (_nParElem - 1) * selem + i * sNode,
							offset + comp * sComp + (_nParElem - 1) * selem - selem - sNode + j * sNode,
							0.0));

						// handle surface diffusion of bound states. binding is handled in residualKernel().
						if (_hasSurfaceDiffusion) {

							for (unsigned int bnd = 0; bnd < _nBound[comp]; bnd++) {
								if (_parSurfDiff[_boundOffset[comp] + bnd] != 0.0) {
									// row: add component offset and jump over previous elements. Go node strides from there for each dispersion block entry
									// col: jump over liquid states, add current bound state offset and jump over previous elements. Go back one elem (and node or adjust for start) and go node strides from there for each dispersion block entry.
									tripletList.push_back(T(offset + comp * sComp + (_nParElem - 1) * selem + i * sNode,
										offset + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd + (_nParElem - 2) * selem - sNode + j * sNode,
										0.0));

									/* add surface diffusion dispersion block to solid */
									if (!qsBinding[offsetBoundComp(ComponentIndex{ comp }) + bnd]) {
										// row: jump over previous elements and over liquid states, add current bound state offset. go node strides from there for each dispersion block entry
										// col: jump over previous elements and over liquid states, add current bound state offset. go back one elem (and node or adjust for start) and go node strides from there for each dispersion block entry.
										tripletList.push_back(T(offset + (_nParElem - 1) * selem + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd + i * sNode,
											offset + (_nParElem - 2) * selem + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) - sNode + bnd + j * sNode,
											0.0));
									}
								}
							}
						}
					}
				}
			}
			if (_nParElem == 3) {
				for (unsigned int comp = 0; comp < _nComp; comp++) {
					for (unsigned int i = 0; i < nNodes; i++) {
						for (unsigned int j = 1; j < 3 * nNodes + 2 - 1; j++) {
							// handle liquid state
							// row: add component offset and jump over previous elem. Go node strides from there for each dispersion block entry
							// col: add component offset. Go node strides from there for each dispersion block entry. adjust for j start
							tripletList.push_back(T(offset + comp * sComp + selem + i * sNode,
								offset + comp * sComp + j * sNode - sNode,
								0.0));

							// handle surface diffusion of bound states. binding is handled in residualKernel().
							if (_hasSurfaceDiffusion) {

								for (unsigned int bnd = 0; bnd < _nBound[comp]; bnd++) {
									if (_parSurfDiff[_boundOffset[comp] + bnd] != 0.0) {
										// row: add component offset and jump over previous elements. Go node strides from there for each dispersion block entry
										// col: jump over liquid states, add current bound state offset and jump over previous elem. go back one elem and go node strides from there for each dispersion block entry. adjust for j start
										tripletList.push_back(T(offset + comp * sComp + selem + i * sNode,
											offset + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd + j * sNode - sNode,
											0.0));

										/* add surface diffusion dispersion block to solid */
										if (!qsBinding[offsetBoundComp(ComponentIndex{ comp }) + bnd]) {
											// row: jump over previous elements and over liquid states, add current bound state offset. go node strides from there for each dispersion block entry
											// col: jump over liquid states, add current bound state offset and jump over previous elem. go node strides from there for each dispersion block entry. adjust for j start
											tripletList.push_back(T(offset + selem + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd + i * sNode,
												offset + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd + j * sNode - sNode,
												0.0));

										}
									}
								}
							}

						}
					}
				}
			}// special case nelements == 3
			/*	boundary elem neighbours (exist only if nelements >= 4)	*/
			if (_nParElem >= 4) {

				for (unsigned int comp = 0; comp < _nComp; comp++) {
					for (unsigned int i = 0; i < nNodes; i++) {
						for (unsigned int j = 1; j < 3 * nNodes + 2; j++) {
							// handle liquid state
							// row: add component offset and jump over previous elem. Go node strides from there for each dispersion block entry
							// col: add component offset. Go node strides from there for each dispersion block entry. adjust for j start
							tripletList.push_back(T(offset + comp * sComp + selem + i * sNode,
								offset + comp * sComp + j * sNode - sNode,
								0.0));

							// handle surface diffusion of bound states. binding is handled in residualKernel().
							if (_hasSurfaceDiffusion) {

								for (unsigned int bnd = 0; bnd < _nBound[comp]; bnd++) {
									if (_parSurfDiff[_boundOffset[comp] + bnd] != 0.0) {
										// row: add component offset and jump over previous elements. Go node strides from there for each dispersion block entry
										// col: jump over liquid states, add current bound state offset and jump over previous elements. Go back one elem and go node strides from there for each dispersion block entry. adjust for j start
										tripletList.push_back(T(offset + comp * sComp + selem + i * sNode,
											offset + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd + j * sNode - sNode,
											0.0));

										/* add surface diffusion dispersion block to solid */
										if (!qsBinding[offsetBoundComp(ComponentIndex{ comp }) + bnd]) {
											// row: jump over previous elements and over liquid states, add current bound state offset. go node strides from there for each dispersion block entry
											// col: jump over previous elements and over liquid states, add current bound state offset. go back one elem and go node strides from there for each dispersion block entry. adjust for j start
											tripletList.push_back(T(offset + selem + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd + i * sNode,
												offset + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd + j * sNode - sNode,
												0.0));

										}
									}
								}
							}

						}
					}
				}

				for (unsigned int comp = 0; comp < _nComp; comp++) {
					for (unsigned int i = 0; i < nNodes; i++) {
						for (unsigned int j = 0; j < 3 * nNodes + 2 - 1; j++) {
							// handle liquid state
							// row: add component offset and jump over previous elements. Go node strides from there for each dispersion block entry
							// col: add component offset and jump over previous elements. Go back one elem and node. Go node strides from there for each dispersion block entry.
							tripletList.push_back(T(offset + comp * sComp + (_nParElem - 2) * selem + i * sNode,
								offset + comp * sComp + (_nParElem - 2) * selem - selem - sNode + j * sNode,
								0.0));

							// handle surface diffusion of bound states. binding is handled in residualKernel().
							if (_hasSurfaceDiffusion) {

								for (unsigned int bnd = 0; bnd < _nBound[comp]; bnd++) {
									if (_parSurfDiff[_boundOffset[comp] + bnd] != 0.0) {
										// row: add component offset and jump over previous elements. Go node strides from there for each dispersion block entry
										// col: jump over liquid states, add current bound state offset and jump over previous elements. Go back one elem and node and go node strides from there for each dispersion block entry
										tripletList.push_back(T(offset + comp * sComp + (_nParElem - 2) * selem + i * sNode,
											offset + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd + (_nParElem - 2) * selem - selem - sNode + j * sNode,
											0.0));

										/* add surface diffusion dispersion block to solid */
										if (!qsBinding[offsetBoundComp(ComponentIndex{ comp }) + bnd]) {
											// row: jump over previous elements and over liquid states, add current bound state offset. go node strides from there for each dispersion block entry
											// col: jump over previous elements and over liquid states, add current bound state offset. go back one elem and node and go node strides from there for each dispersion block entry
											tripletList.push_back(T(offset + (_nParElem - 2) * selem + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd + i * sNode,
												offset + (_nParElem - 2) * selem + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) - selem - sNode + bnd + j * sNode,
												0.0));

										}
									}
								}
							}
						}
					}
				}
			}

			/* Inner elements (exist only if nelements >= 5) */

			if (_nParElem >= 5) {

				for (unsigned int elem = 2; elem < _nParElem - 2; elem++) {

					for (unsigned int comp = 0; comp < _nComp; comp++) {
						for (unsigned int i = 0; i < nNodes; i++) {
							for (unsigned int j = 0; j < 3 * nNodes + 2; j++) {
								// handle liquid state
								// row: add component offset and jump over previous elements. Go node strides from there for each dispersion block entry
								// col: add component offset and jump over previous elements. Go back one elem and node. Go node strides from there for each dispersion block entry.
								tripletList.push_back(T(offset + comp * sComp + elem * selem + i * sNode,
									offset + comp * sComp + elem * selem - selem - sNode + j * sNode,
									0.0));

								// handle surface diffusion of bound states. binding is handled in residualKernel().
								if (_hasSurfaceDiffusion) {

									for (unsigned int bnd = 0; bnd < _nBound[comp]; bnd++) {
										if (_parSurfDiff[_boundOffset[comp] + bnd] != 0.0) {
											// row: add component offset and jump over previous elements. Go node strides from there for each dispersion block entry
											// col: jump over liquid states, add current bound state offset and jump over previous elements. Go back one elem and node and go node strides from there for each dispersion block entry
											tripletList.push_back(T(offset + comp * sComp + elem * selem + i * sNode,
												offset + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd + elem * selem - selem - sNode + j * sNode,
												0.0));

											/* add surface diffusion dispersion block to solid */
											if (!qsBinding[offsetBoundComp(ComponentIndex{ comp }) + bnd]) {
												// row: jump over previous elements and over liquid states, add current bound state offset. go node strides from there for each dispersion block entry
												// col: jump over previous elements and over liquid states, add current bound state offset. go back one elem and node and go node strides from there for each dispersion block entry
												tripletList.push_back(T(offset + elem * selem + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd + i * sNode,
													offset + elem * selem + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd - selem - sNode + j * sNode,
													0.0));

											}
										}
									}
								}
							}
						}
					}
				}

			}

		} // if nelements > 1
	}
	
	unsigned int ParticleDiffusionOperatorDG::calcParDispNNZ()
	{
		return _nComp * ((3u * _nParElem - 2u) * _nParNode * _nParNode + (2u * _nParElem - 3u) * _nParNode);
	}

	/**
	 *@brief adds the time derivative entries from particle equations
	 *@detail since the main diagonal entries are already set, we actually only set the solid phase time derivative entries for the discretized particle mass balance equations
	 */
	void ParticleDiffusionOperatorDG::parTimeDerJacPattern_GRM(std::vector<ParticleDiffusionOperatorDG::T>& tripletList, unsigned int offset, unsigned int colNode, unsigned int secIdx)
	{
		active const* const _parSurfDiff = getSectionDependentSlice(_parSurfDiffusion, _strideBound, secIdx);

		for (unsigned int parNode = 0; parNode < _nParPoints; parNode++)
		{
			for (unsigned int comp = 0; comp < _nComp; comp++)
			{
				for (unsigned int bnd = 0; bnd < _nBound[comp]; bnd++) {
					// row: jump over previous nodes add current component offset
					// col: jump over previous nodes, liquid phase and previous bound states
					tripletList.push_back(T(offset + parNode * strideParNode() + comp,
						offset + parNode * strideParNode() + strideParLiquid() + offsetBoundComp(ComponentIndex{ comp }) + bnd,
						0.0));
				}
			}
		}
	}
	/**
	 * @brief sets the sparsity pattern of the binding Jacobian
	 */
	void ParticleDiffusionOperatorDG::parBindingPattern_GRM(std::vector<ParticleDiffusionOperatorDG::T>& tripletList, const int offset, const unsigned int colNode)
	{
		// every bound state might depend on every bound and liquid state
		for (int parNode = 0; parNode < _nParPoints; parNode++)
		{
			for (int bnd = 0; bnd < _strideBound; bnd++)
			{
				for (int conc = 0; conc < strideParNode(); conc++) {
					// row: jump over previous nodes and liquid states and add current bound state offset
					// col: jump over previous nodes and add current concentration offset (liquid and bound)
					tripletList.push_back(T(offset + parNode * strideParNode() + strideParLiquid() + bnd,
						offset + parNode * strideParNode() + conc, 0.0));
				}
			}
		}
	}
	/**
	 * @brief calculates the DG Jacobian auxiliary block
	 * @param [in] exInt true if exact integration DG scheme
	 * @param [in] elemIdx elem index
	 */
	MatrixXd ParticleDiffusionOperatorDG::getParGBlock(unsigned int elemIdx)
	{
		// Auxiliary Block [ d g(c) / d c ], additionally depends on boundary entries of neighbouring elements
		MatrixXd gBlock = MatrixXd::Zero(_nParNode, _nParNode + 2);
		gBlock.block(0, 1, _nParNode, _nParNode) = _parPolyDerM;

		if (elemIdx == 0 || elemIdx == _nParElem + 1)
		{ // elemIdx out of bounds
			return MatrixXd::Zero(_nParNode, _nParNode + 2);
		}
		if (elemIdx != 1 && elemIdx != _nParElem)
		{ // inner elem
			gBlock.block(0, 0, _nParNode, 1) -= 0.5 * _parInvMM_Leg.block(0, 0, _nParNode, 1);
			gBlock.block(0, 1, _nParNode, 1) += 0.5 * _parInvMM_Leg.block(0, 0, _nParNode, 1);
			gBlock.block(0, _nParNode, _nParNode, 1) -= 0.5 * _parInvMM_Leg.block(0, _nParNode - 1, _nParNode, 1);
			gBlock.block(0, _nParNode + 1, _nParNode, 1) += 0.5 * _parInvMM_Leg.block(0, _nParNode - 1, _nParNode, 1);
		}
		else if (elemIdx == 1u)
		{ // left boundary elem
			if (elemIdx == _nParElem) // special case one elem
				return gBlock * 2.0 / static_cast<double>(_deltaR[(elemIdx - 1)]);
			gBlock.block(0, _nParNode, _nParNode, 1) -= 0.5 * _parInvMM_Leg.block(0, _nParNode - 1, _nParNode, 1);
			gBlock.block(0, _nParNode + 1, _nParNode, 1) += 0.5 * _parInvMM_Leg.block(0, _nParNode - 1, _nParNode, 1);
		}
		else if (elemIdx == _nParElem)
		{ // right boundary elem
			gBlock.block(0, 0, _nParNode, 1) -= 0.5 * _parInvMM_Leg.block(0, 0, _nParNode, 1);
			gBlock.block(0, 1, _nParNode, 1) += 0.5 * _parInvMM_Leg.block(0, 0, _nParNode, 1);
		}
		gBlock *= 2.0 / static_cast<double>(_deltaR[(elemIdx - 1)]);

		return gBlock;
	}

	/**
	 * @brief calculates the num. flux part of a dispersion DG Jacobian block
	 * @param [in] elemIdx elem index
	 * @param [in] leftG left neighbour auxiliary block
	 * @param [in] middleG neighbour auxiliary block
	 * @param [in] rightG neighbour auxiliary block
	 */
	MatrixXd ParticleDiffusionOperatorDG::parAuxBlockGstar(unsigned int elemIdx, MatrixXd leftG, MatrixXd middleG, MatrixXd rightG) {

		// auxiliary block [ d g^* / d c ], depends on whole previous and subsequent elem plus f_Irst entries of subsubsequent elements
		MatrixXd gStarDC = MatrixXd::Zero(_nParNode, 3 * _nParNode + 2);
		// NOTE: N = polyDeg
		// indices  gStarDC    :     0   ,   1   , ..., nNodes; nNodes+1, ..., 2 * nNodes;	2*nNodes+1, ..., 3 * nNodes; 3*nNodes+1
		// derivative index j  : -(N+1)-1, -(N+1),... ,  -1   ;   0     , ...,		N	 ;	  N + 1	  , ..., 2N + 2    ; 2(N+1) +1
		// auxiliary block [d g^* / d c]
		if (elemIdx != 1)
		{
			gStarDC.block(0, _nParNode, 1, _nParNode + 2) += middleG.block(0, 0, 1, _nParNode + 2);
			gStarDC.block(0, 0, 1, _nParNode + 2) += leftG.block(_nParNode - 1, 0, 1, _nParNode + 2);
		}
		if (elemIdx != _nParElem)
		{
			gStarDC.block(_nParNode - 1, _nParNode, 1, _nParNode + 2) += middleG.block(_nParNode - 1, 0, 1, _nParNode + 2);
			gStarDC.block(_nParNode - 1, 2 * _nParNode, 1, _nParNode + 2) += rightG.block(0, 0, 1, _nParNode + 2);
		}
		gStarDC *= 0.5;

		return gStarDC;
	}

	/**
	 * @brief calculates the lifting matrix B
	 * @param [in] elem element index
	 * @param [in] parGeomSurfToVol particle geometry
	 */
	MatrixXd ParticleDiffusionOperatorDG::getParBMatrix(int elem, double parGeomSurfToVol) {
		// also known as "lifting" matrix and includes metric dependent terms for particle discretization
		MatrixXd B = MatrixXd::Zero(_nParNode, _nParNode);
		if (parGeomSurfToVol == _SurfVolRatioSlab)
		{
			B(0, 0) = -1.0;
			B(_nParNode - 1, _nParNode - 1) = 1.0;
		}
		else
		{
			B(0, 0) = -static_cast<double>(_Ir[(elem - 1)][0]);
			B(_nParNode - 1, _nParNode - 1) = static_cast<double>(_Ir[(elem - 1)][_nParNode - 1]);
		}

		return B;
	}

	/**
	 * @brief calculates the dispersion part of the DG jacobian
	 * @param [in] parGeomSurfToVol particle geometry
	 */
	MatrixXd ParticleDiffusionOperatorDG::GSMjacobianParDispBlock(double parGeomSurfToVol)
	{
		MatrixXd dispBlock;

		// We have to match the DGSEM interface, where the dispersion block [ d RHS_disp / d c ] depends on whole previous and subsequent elem plus f_Irst entries of subsubsequent elements
		dispBlock = MatrixXd::Zero(_nParNode, 3 * _nParNode + 2);

		dispBlock.block(0, _nParNode + 1, _nParNode, _nParNode) = _minus_parInvMM_Ar;
		dispBlock *= 2.0 / static_cast<double>(_deltaR[0]) * 2.0 / static_cast<double>(_deltaR[0]);

		return -dispBlock; // *-1 for residual
	}

	/**
	 * @brief calculates the dispersion part of the DG jacobian
	 * @param [in] elemIdx elem index
	 * @param [in] parGeomSurfToVol particle geometry
	 */
	MatrixXd ParticleDiffusionOperatorDG::DGjacobianParDispBlock(unsigned int elemIdx, double parGeomSurfToVol)
	{
		MatrixXd dispBlock;
		// Inner dispersion block [ d RHS_disp / d c ], depends on whole previous and subsequent elem plus f_Irst entries of subsubsequent elements
		dispBlock = MatrixXd::Zero(_nParNode, 3 * _nParNode + 2);

		MatrixXd B = getParBMatrix(elemIdx, parGeomSurfToVol); // "Lifting" matrix
		MatrixXd gBlock = getParGBlock(elemIdx); // current elem auxiliary block matrix
		MatrixXd gStarDC = parAuxBlockGstar(elemIdx, getParGBlock(elemIdx - 1), gBlock, getParGBlock(elemIdx + 1)); // Numerical flux block

		if (parGeomSurfToVol != _SurfVolRatioSlab) // weak form DGSEM required
			dispBlock.block(0, _nParNode, _nParNode, _nParNode + 2) = _minus_InvMM_ST[(elemIdx - 1)] * gBlock;
		else // strong form DGSEM
			dispBlock.block(0, _nParNode, _nParNode, _nParNode + 2) = (_parPolyDerM - _parInvMM[(elemIdx - 1)] * B) * gBlock;

		dispBlock += _parInvMM[(elemIdx - 1)] * B * gStarDC;
		dispBlock *= 2.0 / static_cast<double>(_deltaR[(elemIdx - 1)]);

		return -dispBlock; // *-1 for residual
	}

	void ParticleDiffusionOperatorDG::initializeDGjac(double parGeomSurfToVol)
	{
		// particle jacobian blocks (each is unique)
		_DGjacParDispBlocks = new MatrixXd[_nParElem];

		for (unsigned int block = 0; block < _nParElem; block++)
		{
			if (_parGSM)
				_DGjacParDispBlocks[block] = GSMjacobianParDispBlock(parGeomSurfToVol);
			else
				_DGjacParDispBlocks[block] = DGjacobianParDispBlock(block + 1u, parGeomSurfToVol);
		}
	}
	/**
	 * @brief adds jacobian entries which have been overwritten by the binding kernel (only use for surface diffusion combined with kinetic binding)
	 * @detail only adds the entries d RHS_i / d c^s_i, which lie on the diagonal
	 */
	 int ParticleDiffusionOperatorDG::addSolidDGentries(const int secIdx, linalg::BandedEigenSparseRowIterator& jacBase, const int* const reqBinding)
	 {
	 	active const* const parSurfDiff = getSectionDependentSlice(_parSurfDiffusion, _strideBound, secIdx);

	 	// Get jacobian iterator at first solid entry of first particle of current type
	 	linalg::BandedEigenSparseRowIterator jac = jacBase + strideParLiquid();

	 	for (unsigned int elem = 0; elem < _nParElem; elem++)
	 		addDiagonalSolidJacobianEntries(_DGjacParDispBlocks[elem].block(0, _nParNode + 1, _nParNode, _nParNode), jac, reqBinding, parSurfDiff);

	 	return 1;
	 }
	/**
	 * @brief adds a state block into the system jacobian.
	 * @param [in] block (sub)block whose diagonal entries are to be added
	 * @param [in] jac row iterator at first (i.e. upper) entry
	 * @param [in] reqBinding pointer to binding kinetics
	 * @param [in] type particle type
	 */
	void ParticleDiffusionOperatorDG::addDiagonalSolidJacobianEntries(Eigen::MatrixXd block, linalg::BandedEigenSparseRowIterator& jac, const int* const reqBinding, const active* const surfDiffPtr)
	{
		for (unsigned int i = 0; i < block.rows(); i++, jac += strideParLiquid())
		{
			for (unsigned int comp = 0; comp < _nComp; comp++)
			{
				for (unsigned int bnd = 0; bnd < _nBound[comp]; bnd++, ++jac)
				{

					if (static_cast<double>(surfDiffPtr[offsetBoundComp(ComponentIndex{ comp }) + bnd]) != 0.0
						&& !reqBinding[offsetBoundComp(ComponentIndex{ comp }) + bnd]) {
						// row, col: at current node and bound state
						jac[0] += block(i, i)
							* static_cast<double>(surfDiffPtr[offsetBoundComp(ComponentIndex{ comp }) + bnd]);
					}
				}
			}
		}
	}
	/**
	 * @brief adds a state block into the particle jacobian.
	 * @param [in] block (sub)block to be added
	 * @param [in] jac row iterator at first (i.e. upper) entry
	 * @param [in] parDiff pointer to particle diffusion parameters
	 * @param [in] surfDiff pointer to particle surface diffusion parameters
	 * @param [in] beta_p pointer to particle porosity parameters
	 * @param [in] nonKinetic pointer to binding kinetics parameters
	 * @param [in] type particle type
	 * @param [in] nBlocks number of blocks, i.e. elements/elements, to be inserted
	 * @param [in] offRowToCol column to row offset (i.e. start at upper left corner of block)
	 */
	void ParticleDiffusionOperatorDG::insertParJacBlock(Eigen::MatrixXd block, linalg::BandedEigenSparseRowIterator& jac, const active* const parDiff, const active* const surfDiff, const active* const beta_p, const int* nonKinetic, unsigned int nBlocks, int offRowToCol)
	{
		for (unsigned int elem = 0; elem < nBlocks; elem++) {
			for (unsigned int i = 0; i < block.rows(); i++) {
				for (unsigned int comp = 0; comp < _nComp; comp++, ++jac) {
					for (unsigned int j = 0; j < block.cols(); j++) {
						/* liquid on liquid blocks */
						// row: at current node and component; col: jump to node j
						jac[(j - i) * strideParNode() + offRowToCol] = block(i, j) * static_cast<double>(parDiff[comp]);
					}
					/* liquid on solid blocks */
					for (unsigned int bnd = 0; bnd < _nBound[comp]; bnd++) {
						if (static_cast<double>(surfDiff[offsetBoundComp(ComponentIndex{ comp }) + bnd]) != 0.0) {
							for (unsigned int j = 0; j < block.cols(); j++) {
								// row: at current node and component; col: jump to node j and to current bound state
								jac[(j - i) * strideParNode() + offRowToCol + strideParLiquid() - comp
									+ offsetBoundComp(ComponentIndex{ comp }) + bnd
								]
									= block(i, j) * static_cast<double>(beta_p[comp])
									* static_cast<double>(surfDiff[offsetBoundComp(ComponentIndex{ comp }) + bnd]);
							}
						}
					}
				}
				/* solid on solid blocks */
				for (unsigned int comp = 0; comp < _nComp; comp++) {
					for (unsigned int bnd = 0; bnd < _nBound[comp]; bnd++, ++jac) {
						if (static_cast<double>(surfDiff[offsetBoundComp(ComponentIndex{ comp }) + bnd]) != 0.0
							&& !nonKinetic[offsetBoundComp(ComponentIndex{ comp }) + bnd]) {
							for (unsigned int j = 0; j < block.cols(); j++) {
								// row: at current node and bound state; col: jump to node j
								jac[(j - i) * strideParNode() + offRowToCol + bnd]
									= block(i, j)
									* static_cast<double>(surfDiff[offsetBoundComp(ComponentIndex{ comp }) + bnd]);
							}
						}
					}
				}
			}
		}
	}
	/**
	 * @brief analytically calculates the static (per section) particle diffusion Jacobian
	 * @return 1 if jacobain calculation fits the predefined pattern of the jacobian, 0 if not.
	 */
	int ParticleDiffusionOperatorDG::calcStaticAnaParticleDiffJacobian(const int secIdx, const int colNode, const int offsetLocalCp, Eigen::SparseMatrix<double, RowMajor>& globalJac)
	{
		// Prepare parameters
		const active* const parDiff = getSectionDependentSlice(_parDiffusion, _nComp, secIdx);

		// Ordering of particle surface diffusion:
		// bnd0comp0, bnd1comp0, bnd0comp1, bnd1comp1, bnd0comp2, bnd1comp2
		const active* const  parSurfDiff = getSectionDependentSlice(_parSurfDiffusion, _strideBound, secIdx);

		const active* const invBetaP = &_invBetaP[0];

		const int* const reqBinding = _binding->reactionQuasiStationarity();

		// (global) strides
		unsigned int selem = _nParNode * strideParNode();
		unsigned int sNode = strideParNode();
		unsigned int sComp = 1u;
		unsigned int nNodes = _nParNode;

		/* Special case */
		if (_nParElem == 1) {

			linalg::BandedEigenSparseRowIterator jacIt(globalJac, offsetLocalCp); // row iterator starting at first elem, first component

			insertParJacBlock(_DGjacParDispBlocks[0].block(0, nNodes + 1, nNodes, nNodes), jacIt, parDiff, parSurfDiff, invBetaP, reqBinding, 1u, 0);
			return 1;
		}

		/* Special case */
		if (_nParElem == 2) {

			linalg::BandedEigenSparseRowIterator jacIt(globalJac, offsetLocalCp); // row iterator starting at first elem, first component

			insertParJacBlock(_DGjacParDispBlocks[0].block(0, nNodes + 1, nNodes, 2 * nNodes), jacIt, parDiff, parSurfDiff, invBetaP, reqBinding, 1u, 0);
			// right Bacobian block, iterator is already moved to second elem
			insertParJacBlock(_DGjacParDispBlocks[1].block(0, 1, nNodes, 2 * nNodes), jacIt, parDiff, parSurfDiff, invBetaP, reqBinding, 1u, -strideParElem());
			return 1;
		}

		/* Special case */
		if (_nParElem == 3) {

			linalg::BandedEigenSparseRowIterator jacIt(globalJac, offsetLocalCp + strideParElem()); // row iterator starting at first elem, first component

			insertParJacBlock(_DGjacParDispBlocks[1].block(0, 1, nNodes, 3 * nNodes), jacIt, parDiff, parSurfDiff, invBetaP, reqBinding, 1u, -strideParElem());
		}

		/* Inner elements (exist only if nelements >= 5) */
		if (_nParElem >= 5) {

			linalg::BandedEigenSparseRowIterator jacIt(globalJac, offsetLocalCp + strideParElem() * 2); // row iterator starting at third elem, first component

			// insert all (nElem - 4) inner elem blocks
			for (unsigned int elem = 2; elem < _nParElem - 2; elem++)
				insertParJacBlock(_DGjacParDispBlocks[elem], jacIt, parDiff, parSurfDiff, invBetaP, reqBinding, 1u, -(strideParElem() + strideParNode()));
		}

		/*	boundary elem neighbours (exist only if nelements >= 4)	*/
		if (_nParElem >= 4) {

			linalg::BandedEigenSparseRowIterator jacIt(globalJac, offsetLocalCp + strideParElem()); // row iterator starting at second elem, first component

			insertParJacBlock(_DGjacParDispBlocks[1].block(0, 1, nNodes, 3 * nNodes + 1), jacIt, parDiff, parSurfDiff, invBetaP, reqBinding, 1u, -strideParElem());

			jacIt += (_nParElem - 4) * strideParElem(); // move iterator to preultimate elem (already at third elem)
			insertParJacBlock(_DGjacParDispBlocks[_nParElem - 2u].block(0, 0, nNodes, 3 * nNodes + 1), jacIt, parDiff, parSurfDiff, invBetaP, reqBinding, 1u, -(strideParElem() + strideParNode()));
		}

		/*			boundary elements (exist only if nelements >= 3)			*/
		if (_nParElem >= 3) {

			linalg::BandedEigenSparseRowIterator jacIt(globalJac, offsetLocalCp); // row iterator starting at first elem, first component

			insertParJacBlock(_DGjacParDispBlocks[0].block(0, nNodes + 1, nNodes, 2 * nNodes + 1), jacIt, parDiff, parSurfDiff, invBetaP, reqBinding, 1u, 0);

			jacIt += (_nParElem - 2) * strideParElem(); // move iterator to last elem (already at second elem)
			insertParJacBlock(_DGjacParDispBlocks[_nParElem - 1u].block(0, 0, nNodes, 2 * nNodes + 1), jacIt, parDiff, parSurfDiff, invBetaP, reqBinding, 1u, -(strideParElem() + strideParNode()));
		}

		return true;
	}

	int ParticleDiffusionOperatorDG::calcFilmDiffJacobian(unsigned int secIdx, const int offsetCp, const int offsetC, const int nBulkPoints, const int nParType, const double colPorosity, const active* const parTypeVolFrac, Eigen::SparseMatrix<double, RowMajor>& globalJac, bool outliersOnly)
	{
		// lifting matrix entry for exact integration scheme depends on metrics for sphere and cylinder
		double exIntLiftContribution = static_cast<double>(_Ir[_nParElem - 1][_nParNode - 1]);
		if (_parGeomSurfToVol == _SurfVolRatioSlab)
			exIntLiftContribution = 1.0;

		// Ordering of diffusion:
		// sec0type0comp0, sec0type0comp1, sec0type0comp2, sec0type1comp0, sec0type1comp1, sec0type1comp2,
		// sec1type0comp0, sec1type0comp1, sec1type0comp2, sec1type1comp0, sec1type1comp1, sec1type1comp2, ...
		active const* const filmDiff = getSectionDependentSlice(_filmDiffusion, _nComp, secIdx);

		linalg::BandedEigenSparseRowIterator jacCl(globalJac, offsetC);
		linalg::BandedEigenSparseRowIterator jacCp(globalJac, offsetCp + (_nParPoints - 1) * strideParNode()); // iterator at the outer particle boundary

		for (unsigned int blk = 0; blk < nBulkPoints; blk++)
		{
			for (unsigned int comp = 0; comp < _nComp; comp++, ++jacCp, ++jacCl) {
				// add Cl on Cl entries (added since these entries are also touched by bulk jacobian)
				// row: already at bulk phase. already at current node and component.
				// col: already at bulk phase. already at current node and component.
				if (!outliersOnly)
					jacCl[0] += static_cast<double>(filmDiff[comp]) * (1.0 - colPorosity) / colPorosity
					* _parGeomSurfToVol / static_cast<double>(_parRadius)
					* static_cast<double>(parTypeVolFrac[_parTypeIdx + blk * nParType]);
				// add Cl on Cp entries (added since these entries are also touched by bulk jacobian)
				// row: already at bulk phase. already at current node and component.
				// col: go to current particle phase entry.
				jacCl[jacCp.row() - jacCl.row()] = -static_cast<double>(filmDiff[comp]) * (1.0 - colPorosity) / colPorosity
					* _parGeomSurfToVol / static_cast<double>(_parRadius)
					* static_cast<double>(parTypeVolFrac[_parTypeIdx + blk * nParType]);


				unsigned int entry = jacCp.row();
				for (int node = _parPolyDeg; node >= 0; node--, jacCp -= strideParNode()) {
					// row: already at particle. Already at current node and liquid state.
					// col: original entry at outer node.
					if (!outliersOnly) // Cp on Cb
						jacCp[entry - jacCp.row()]
						+= static_cast<double>(filmDiff[comp]) * 2.0 / static_cast<double>(_deltaR[0]) * _parInvMM[_nParElem - 1](node, _nParNode - 1) * exIntLiftContribution / static_cast<double>(_parPorosity) / static_cast<double>(_poreAccessFactor[comp]);
					// row: already at particle. Already at current node and liquid state.
					// col: go to current bulk phase.
					jacCp[jacCl.row() - jacCp.row()]
						= -static_cast<double>(filmDiff[comp]) * 2.0 / static_cast<double>(_deltaR[0]) * _parInvMM[_nParElem - 1](node, _nParNode - 1) * exIntLiftContribution / static_cast<double>(_parPorosity) / static_cast<double>(_poreAccessFactor[comp]);
				}
				// set back iterator to first node as required by component loop
				jacCp += _nParNode * strideParNode();
			}
			if (blk < nBulkPoints - 1) // execute iteration statement only when condition is true in next loop.
				jacCp += _strideBound + (_nParPoints - 1) * strideParNode();
		}

		return 1;
	}

	bool ParticleDiffusionOperatorDG::setParameter(const ParameterId& pId, double value)
	{
		if (singleTypeMultiplexCompTypeSecParameterValue(pId, hashString("PAR_DIFFUSION"), _parDiffusionMode, _parDiffusion, _nComp, _parTypeIdx, value, nullptr))
			return true;
		if (singleTypeMultiplexBndCompTypeSecParameterValue(pId, hashString("PAR_SURFDIFFUSION"), _parSurfDiffusionMode, _parSurfDiffusion, _nComp, _strideBound, _boundOffset, _parTypeIdx, value, nullptr))
			return true;

		if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_RADIUS"), _singleParRadius, _parRadius, _parTypeIdx, value, nullptr))
			return true;
		if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_CORERADIUS"), _singleParCoreRadius, _parCoreRadius, _parTypeIdx, value, nullptr))
			return true;
		if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_POROSITY"), _singleParPorosity, _parPorosity, _parTypeIdx, value, nullptr))
			return true;

		if (_singleParDepSurfDiffusion)
			if (_parDepSurfDiffusion && _parDepSurfDiffusion->setParameter(pId, value))
				return true;

		return false;
	}

	bool ParticleDiffusionOperatorDG::setParameter(const ParameterId& pId, int value)
	{
		if (_singleParDepSurfDiffusion)
			if (_parDepSurfDiffusion && _parDepSurfDiffusion->setParameter(pId, value))
				return true;

		return false;
	}

	bool ParticleDiffusionOperatorDG::setParameter(const ParameterId& pId, bool value)
	{
		if (_singleParDepSurfDiffusion)
			if (_parDepSurfDiffusion && _parDepSurfDiffusion->setParameter(pId, value))
				return true;

		return false;
	}

	bool ParticleDiffusionOperatorDG::setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& pId, double value)
	{
		if (singleTypeMultiplexCompTypeSecParameterValue(pId, hashString("PAR_DIFFUSION"), _parDiffusionMode, _parDiffusion, _nComp, _parTypeIdx, value, &sensParams))
			return true;
		if (singleTypeMultiplexBndCompTypeSecParameterValue(pId, hashString("PAR_SURFDIFFUSION"), _parSurfDiffusionMode, _parSurfDiffusion, _nComp, _strideBound, _boundOffset, _parTypeIdx, value, &sensParams))
			return true;

		if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_RADIUS"), _singleParRadius, _parRadius, _parTypeIdx, value, &sensParams))
			return true;
		if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_CORERADIUS"), _singleParCoreRadius, _parCoreRadius, _parTypeIdx, value, &sensParams))
			return true;
		if (singleTypeMultiplexTypeParameterValue(pId, hashString("PAR_POROSITY"), _singleParPorosity, _parPorosity, _parTypeIdx, value, &sensParams))
			return true;

		if (model::setSensitiveParameterValue(pId, value, sensParams, std::vector< IParameterStateDependence*>(1, _parDepSurfDiffusion), _singleParDepSurfDiffusion))
			return true;

		return false;
	}

	bool ParticleDiffusionOperatorDG::setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue)
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

		if (model::setSensitiveParameter(pId, adDirection, adValue, sensParams, std::vector< IParameterStateDependence*>(1, _parDepSurfDiffusion), _singleParDepSurfDiffusion))
		{
			LOG(Debug) << "Found parameter " << pId << " in surface diffusion parameter dependence: Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		if (singleTypeMultiplexTypeParameterAD(pId, hashString("PAR_RADIUS"), _singleParRadius, _parRadius, _parTypeIdx, adDirection, adValue, sensParams))
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		if (singleTypeMultiplexTypeParameterAD(pId, hashString("PAR_CORERADIUS"), _singleParCoreRadius, _parCoreRadius, _parTypeIdx, adDirection, adValue, sensParams))
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		if (singleTypeMultiplexTypeParameterAD(pId, hashString("PAR_POROSITY"), _singleParPorosity, _parPorosity, _parTypeIdx, adDirection, adValue, sensParams))
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		return false;
	}


}  // namespace parts

}  // namespace model

}  // namespace cadet
