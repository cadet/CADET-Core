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

#include "ParamReaderHelper.hpp"
#include "AdUtils.hpp"
#include "SimulationTypes.hpp"
#include "model/ParameterDependence.hpp"
#include "SensParamUtil.hpp"
#include "ConfigurationHelper.hpp"

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
	ParticleDiffusionOperatorDG::ParticleDiffusionOperatorDG() : _nParElem(nullptr), _nParPointsBeforeType(nullptr),
		_parPolyDeg(nullptr), _nParNode(nullptr), _nParPoints(nullptr), _parExactInt(nullptr), _parGSM(nullptr), //_parTypeOffset(nullptr),
		_nBound(nullptr), _boundOffset(nullptr), _strideBound(nullptr), _nBoundBeforeType(nullptr), _offsetSurfDiff(nullptr), _deltaR(nullptr), _parNodes(nullptr),
		_parPolyDerM(nullptr), _minus_InvMM_ST(nullptr), _minus_parInvMM_Ar(nullptr), _parInvWeights(nullptr), _parInvMM(nullptr), _parInvMM_Leg(nullptr),
		_Ir(nullptr), _Dr(nullptr), _secondOrderStiffnessM(nullptr), _DGjacParDispBlocks(nullptr), _g_p(nullptr), _g_pSum(nullptr),
		_surfaceFluxParticle(nullptr), _localFlux(nullptr)
	{
	}

	ParticleDiffusionOperatorDG::~ParticleDiffusionOperatorDG() CADET_NOEXCEPT
	{
		delete[] _nParElem;
		delete[] _nParPointsBeforeType;
		delete[] _parPolyDeg;
		delete[] _nParNode;
		delete[] _nParPoints;
		delete[] _parExactInt;
		delete[] _parGSM;
		//delete[] _parTypeOffset;
		delete[] _nBound;
		delete[] _boundOffset;
		delete[] _strideBound;
		delete[] _nBoundBeforeType;

		delete[] _offsetSurfDiff;

		delete[] _deltaR;
		delete[] _parNodes;
		delete[] _parPolyDerM;
		delete[] _minus_InvMM_ST;
		delete[] _minus_parInvMM_Ar;
		delete[] _parInvWeights;
		delete[] _parInvMM;
		delete[] _parInvMM_Leg;
		delete[] _Ir;
		delete[] _Dr;
		delete[] _secondOrderStiffnessM;

		delete[] _DGjacParDispBlocks;

		delete[] _g_p;
		delete[] _g_pSum;
		delete[] _surfaceFluxParticle;
		delete[] _localFlux;
	}

	void ParticleDiffusionOperatorDG::clearParDepSurfDiffusion()
	{
		if (_singleParDepSurfDiffusion)
		{
			if (!_parDepSurfDiffusion.empty())
				delete _parDepSurfDiffusion[0];
		}
		else
		{
			for (IParameterStateDependence* pd : _parDepSurfDiffusion)
				delete pd;
		}

		_parDepSurfDiffusion.clear();
	}

	bool ParticleDiffusionOperatorDG::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int nParType, const int strideBulkComp)
	{
		_strideBulkComp = strideBulkComp;

		_nParType = nParType;
		_nComp = nComp;

		_filmDiffusion.resize(_nComp * _nParType); // filled in notifyDiscontinuousSectionTransition
		_poreAccessFactor.resize(_nComp * _nParType); // filled in notifyDiscontinuousSectionTransition

		std::vector<int> nBound;
		const bool newNBoundInterface = paramProvider.exists("NBOUND");

		paramProvider.pushScope("discretization");

		const bool firstConfigCall = _parPolyDeg == nullptr;

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

		std::vector<int> parPolyDegs(_nParType);
		std::vector<int> parNelems(_nParType);
		std::vector<bool> parExactInts(_nParType, true);

		if (firstConfigCall) // avoid memory leaks from reinitialization
		{
			_parPolyDeg = new unsigned int[_nParType];
			_nParElem = new unsigned int[_nParType];
			_parExactInt = new bool[_nParType];
			_parGSM = new bool[_nParType];
		}

		if (paramProvider.exists("PAR_POLYDEG"))
		{
			parPolyDegs = paramProvider.getIntArray("PAR_POLYDEG");

			if ((std::any_of(parPolyDegs.begin(), parPolyDegs.end(), [](int value) { return value < 1; })))
				throw InvalidParameterException("Particle polynomial degrees must be at least 1!");
			parNelems = paramProvider.getIntArray("PAR_NELEM");

			if ((std::any_of(parNelems.begin(), parNelems.end(), [](int value) { return value < 1; })))
				throw InvalidParameterException("Particle number of elements must be at least 1!");

			if (paramProvider.exists("PAR_EXACT_INTEGRATION"))
				parExactInts = paramProvider.getBoolArray("PAR_EXACT_INTEGRATION");

			if ((std::any_of(parExactInts.begin(), parExactInts.end(), [](bool value) { return !value; })))
				LOG(Warning) << "Inexact integration method (cf. PAR_EXACT_INTEGRATION) in particles might add severe! stiffness to the system and disables consistent initialization!";

			if (parPolyDegs.size() == 1)
			{
				// Multiplex number of particle elements to all particle types
				for (unsigned int i = 0; i < _nParType; ++i)
					std::fill(_parPolyDeg, _parPolyDeg + _nParType, parPolyDegs[0]);
			}
			else if (parPolyDegs.size() < _nParType)
				throw InvalidParameterException("Field PAR_POLYDEG must have 1 or _nParType (" + std::to_string(_nParType) + ") entries");
			else
				std::copy_n(parPolyDegs.begin(), _nParType, _parPolyDeg);
			if (parNelems.size() == 1)
			{
				// Multiplex number of particle elements to all particle types
				for (unsigned int i = 0; i < _nParType; ++i)
					std::fill(_nParElem, _nParElem + _nParType, parNelems[0]);
			}
			else if (parNelems.size() < _nParType)
				throw InvalidParameterException("Field PAR_NELEM must have 1 or _nParType (" + std::to_string(_nParType) + ") entries");
			else
				std::copy_n(parNelems.begin(), _nParType, _nParElem);
			if (parExactInts.size() == 1)
			{
				// Multiplex number of particle elements to all particle types
				for (unsigned int i = 0; i < _nParType; ++i)
					std::fill(_parExactInt, _parExactInt + _nParType, parExactInts[0]);
			}
			else if (parExactInts.size() < _nParType)
				throw InvalidParameterException("Field PAR_EXACT_INTEGRATION must have 1 or _nParType (" + std::to_string(_nParType) + ") entries");
			else
				std::copy_n(parExactInts.begin(), _nParType, _parExactInt);
		}
		else if (paramProvider.exists("NPAR"))
		{
			const std::vector<int> _nParPoints = paramProvider.getIntArray("NPAR");
			if ((_nParPoints.size() > 1) && (_nParPoints.size() < _nParType))
				throw InvalidParameterException("Field NPAR must have 1 or _nParType (" + std::to_string(_nParType) + ") entries");

			for (unsigned int par = 0; par < _nParType; par++)
				_parPolyDeg[par] = std::max(1, std::min(_nParPoints[par] - 1, 4));
		}
		else
			throw InvalidParameterException("Specify field PAR_POLYDEG (or NPAR)");

		if (paramProvider.exists("PAR_GSM"))
		{
			std::vector<bool> parGSMs = paramProvider.getBoolArray("PAR_GSM");
			if (parGSMs.size() == 1)
				std::fill(_parGSM, _parGSM + _nParType, parGSMs[0]);
			else
				std::copy_n(parGSMs.begin(), _nParType, _parGSM);
			for (int type = 0; type < _nParType; type++)
				if (_parGSM[type] && _nParElem[type] != 1)
					throw InvalidParameterException("Field PAR_NELEM must equal one to use a GSM discretization in the corresponding particle type");
		}
		else // Use GSM as default for particle discretization
		{
			for (int type = 0; type < _nParType; type++)
				_parGSM[type] = (_nParElem[type] == 1);
		}

		initializeDG();

		if (firstConfigCall)
			_nBound = new unsigned int[_nComp * _nParType];
		if (nBound.size() < _nComp * _nParType)
		{
			// Multiplex number of bound states to all particle types
			for (unsigned int i = 0; i < _nParType; ++i)
				std::copy_n(nBound.begin(), _nComp, _nBound + i * _nComp);
		}
		else
			std::copy_n(nBound.begin(), _nComp * _nParType, _nBound);

		const unsigned int nTotalBound = std::accumulate(_nBound, _nBound + _nComp * _nParType, 0u);

		// Precompute offsets and total number of bound states (DOFs in solid phase)
		if (firstConfigCall)
		{
			_boundOffset = new unsigned int[_nComp * _nParType];
			_strideBound = new unsigned int[_nParType + 1];
			_nBoundBeforeType = new unsigned int[_nParType];
		}
		_strideBound[_nParType] = nTotalBound;
		_nBoundBeforeType[0] = 0;
		for (unsigned int j = 0; j < _nParType; ++j)
		{
			unsigned int* const ptrOffset = _boundOffset + j * _nComp;
			unsigned int* const ptrBound = _nBound + j * _nComp;

			ptrOffset[0] = 0;
			for (unsigned int i = 1; i < _nComp; ++i)
			{
				ptrOffset[i] = ptrOffset[i - 1] + ptrBound[i - 1];
			}
			_strideBound[j] = ptrOffset[_nComp - 1] + ptrBound[_nComp - 1];

			if (j != _nParType - 1)
				_nBoundBeforeType[j + 1] = _nBoundBeforeType[j] + _strideBound[j];
		}

		// Precompute offsets of particle type DOFs
		if (firstConfigCall)
		{
			//_parTypeOffset = new unsigned int[_nParType + 1];
			_nParPointsBeforeType = new unsigned int[_nParType + 1];
		}
		//_parTypeOffset[0] = 0;
		_nParPointsBeforeType[0] = 0;
		unsigned int nTotalParPoints = 0;
		for (unsigned int j = 1; j < _nParType + 1; ++j)
		{
			//_parTypeOffset[j] = _parTypeOffset[j - 1] + (_nComp + _strideBound[j - 1]) * _nParPoints[j - 1] * _nPoints;
			_nParPointsBeforeType[j] = _nParPointsBeforeType[j - 1] + _nParPoints[j - 1];
			nTotalParPoints += _nParPoints[j - 1];
		}
		_nParPointsBeforeType[_nParType] = nTotalParPoints;

		// Configure particle discretization
		_parElementSize.resize(_offsetMetric[_nParType]);
		_parCenterRadius.resize(_offsetMetric[_nParType]);
		_parOuterSurfAreaPerVolume.resize(_offsetMetric[_nParType]);
		_parInnerSurfAreaPerVolume.resize(_offsetMetric[_nParType]);

		// Read particle discretization mode and default to "EQUIDISTANT_PAR"
		_parDiscMode = std::vector<ParticleDiscretizationMode>(_nParType, ParticleDiscretizationMode::Equidistant);
		std::vector<std::string> pdt = paramProvider.getStringArray("PAR_DISC_TYPE");
		if ((pdt.size() == 1) && (_nParType > 1))
		{
			// Multiplex using first value
			pdt.resize(_nParType, pdt[0]);
		}
		else if (pdt.size() < _nParType)
			throw InvalidParameterException("Field PAR_DISC_TYPE contains too few elements (" + std::to_string(_nParType) + " required)");

		for (unsigned int i = 0; i < _nParType; ++i)
		{
			if (pdt[i] == "EQUIVOLUME_PAR")
				_parDiscMode[i] = ParticleDiscretizationMode::Equivolume;
			else if (pdt[i] == "USER_DEFINED_PAR")
				_parDiscMode[i] = ParticleDiscretizationMode::UserDefined;
		}

		// Read particle geometry and default to "SPHERICAL"
		paramProvider.popScope();
		_parGeomSurfToVol = std::vector<double>(_nParType, _SurfVolRatioSphere);
		if (paramProvider.exists("PAR_GEOM"))
		{
			std::vector<std::string> pg = paramProvider.getStringArray("PAR_GEOM");
			if ((pg.size() == 1) && (_nParType > 1))
			{
				// Multiplex using first value
				pg.resize(_nParType, pg[0]);
			}
			else if (pg.size() < _nParType)
				throw InvalidParameterException("Field PAR_GEOM contains too few elements (" + std::to_string(_nParType) + " required)");

			for (unsigned int i = 0; i < _nParType; ++i)
			{
				if (pg[i] == "SPHERE")
					_parGeomSurfToVol[i] = _SurfVolRatioSphere;
				else if (pg[i] == "CYLINDER")
					_parGeomSurfToVol[i] = _SurfVolRatioCylinder;
				else if (pg[i] == "SLAB")
					_parGeomSurfToVol[i] = _SurfVolRatioSlab;
				else
					throw InvalidParameterException("Unknown particle geometry type \"" + pg[i] + "\" at index " + std::to_string(i) + " of field PAR_GEOM");
			}
		}
		paramProvider.pushScope("discretization");

		if (paramProvider.exists("PAR_DISC_VECTOR"))
		{
			_parDiscVector = paramProvider.getDoubleArray("PAR_DISC_VECTOR");
			if (_parDiscVector.size() < nTotalParPoints + _nParType)
				throw InvalidParameterException("Field PAR_DISC_VECTOR contains too few elements (Sum [NPAR + 1] = " + std::to_string(nTotalParPoints + _nParType) + " required)");
		}

		// Determine whether surface diffusion optimization is applied (decreases Jacobian size) //@TODO?
		const bool optimizeParticleJacobianBandwidth = paramProvider.exists("OPTIMIZE_PAR_BANDWIDTH") ? paramProvider.getBool("OPTIMIZE_PAR_BANDWIDTH") : true;

		paramProvider.popScope();

		// ==== Construct and configure parameter dependencies
		clearParDepSurfDiffusion();
		bool parSurfDiffDepConfSuccess = true;
		if (paramProvider.exists("PAR_SURFDIFFUSION_DEP"))
		{
			const std::vector<std::string> psdDepNames = paramProvider.getStringArray("PAR_SURFDIFFUSION_DEP");
			if ((psdDepNames.size() == 1) || (_nParType == 1))
				_singleParDepSurfDiffusion = true;

			if (!_singleParDepSurfDiffusion && (psdDepNames.size() < _nParType))
				throw InvalidParameterException("Field PAR_SURFDIFFUSION_DEP contains too few elements (" + std::to_string(_nParType) + " required)");
			else if (_singleParDepSurfDiffusion && (psdDepNames.size() != 1))
				throw InvalidParameterException("Field PAR_SURFDIFFUSION_DEP requires (only) 1 element");

			if (_singleParDepSurfDiffusion)
			{
				if ((psdDepNames[0] == "") || (psdDepNames[0] == "NONE") || (psdDepNames[0] == "DUMMY"))
				{
					_hasParDepSurfDiffusion = false;
					_singleParDepSurfDiffusion = true;
					_parDepSurfDiffusion = std::vector<IParameterStateDependence*>(_nParType, nullptr);
				}
				else
				{
					IParameterStateDependence* const pd = helper.createParameterStateDependence(psdDepNames[0]);
					if (!pd)
						throw InvalidParameterException("Unknown parameter dependence " + psdDepNames[0]);

					_parDepSurfDiffusion = std::vector<IParameterStateDependence*>(_nParType, pd);
					parSurfDiffDepConfSuccess = pd->configureModelDiscretization(paramProvider, _nComp, _nBound, _boundOffset);
					_hasParDepSurfDiffusion = true;
				}
			}
			else
			{
				_parDepSurfDiffusion = std::vector<IParameterStateDependence*>(_nParType, nullptr);

				for (unsigned int i = 0; i < _nParType; ++i)
				{
					if ((psdDepNames[0] == "") || (psdDepNames[0] == "NONE") || (psdDepNames[0] == "DUMMY"))
						continue;

					_parDepSurfDiffusion[i] = helper.createParameterStateDependence(psdDepNames[i]);
					if (!_parDepSurfDiffusion[i])
						throw InvalidParameterException("Unknown parameter dependence " + psdDepNames[i]);

					parSurfDiffDepConfSuccess = _parDepSurfDiffusion[i]->configureModelDiscretization(paramProvider, _nComp, _nBound + i * _nComp, _boundOffset + i * _nComp) && parSurfDiffDepConfSuccess;
				}

				_hasParDepSurfDiffusion = std::any_of(_parDepSurfDiffusion.cbegin(), _parDepSurfDiffusion.cend(), [](IParameterStateDependence const* pd) -> bool { return pd; });
			}
		}
		else
		{
			_hasParDepSurfDiffusion = false;
			_singleParDepSurfDiffusion = true;
			_parDepSurfDiffusion = std::vector<IParameterStateDependence*>(_nParType, nullptr);
		}

		if (optimizeParticleJacobianBandwidth)
		{
			// Check whether surface diffusion is present
			_hasSurfaceDiffusion = std::vector<bool>(_nParType, false);
			if (paramProvider.exists("PAR_SURFDIFFUSION"))
			{
				const std::vector<double> surfDiff = paramProvider.getDoubleArray("PAR_SURFDIFFUSION");
				for (unsigned int i = 0; i < _nParType; ++i)
				{
					// Assume particle surface diffusion if a parameter dependence is present
					if (_parDepSurfDiffusion[i])
					{
						_hasSurfaceDiffusion[i] = true;
						continue;
					}

					double const* const lsd = surfDiff.data() + _nBoundBeforeType[i];

					// Check surface diffusion coefficients of each particle type
					for (unsigned int j = 0; j < _strideBound[i]; ++j)
					{
						if (lsd[j] != 0.0)
						{
							_hasSurfaceDiffusion[i] = true;
							break;
						}
					}
				}
			}
		}
		else
		{
			// Assume that surface diffusion is present
			_hasSurfaceDiffusion = std::vector<bool>(_nParType, true);
		}

		//_curSection = -1;

		//// ==== Construct and configure binding model
		//clearBindingModels();
		//_binding = std::vector<IBindingModel*>(_nParType, nullptr);

		//std::vector<std::string> bindModelNames = { "NONE" };
		//if (paramProvider.exists("ADSORPTION_MODEL"))
		//	bindModelNames = paramProvider.getStringArray("ADSORPTION_MODEL");

		//if (paramProvider.exists("ADSORPTION_MODEL_MULTIPLEX"))
		//	_singleBinding = (paramProvider.getInt("ADSORPTION_MODEL_MULTIPLEX") == 1);
		//else
		//{
		//	// Infer multiplex mode
		//	_singleBinding = (bindModelNames.size() == 1);
		//}

		//if (!_singleBinding && (bindModelNames.size() < _nParType))
		//	throw InvalidParameterException("Field ADSORPTION_MODEL contains too few elements (" + std::to_string(_nParType) + " required)");
		//else if (_singleBinding && (bindModelNames.size() != 1))
		//	throw InvalidParameterException("Field ADSORPTION_MODEL requires (only) 1 element");

		//bool bindingConfSuccess = true;
		//for (unsigned int i = 0; i < _nParType; ++i)
		//{
		//	if (_singleBinding && (i > 0))
		//	{
		//		// Reuse first binding model
		//		_binding[i] = _binding[0];
		//	}
		//	else
		//	{
		//		_binding[i] = helper.createBindingModel(bindModelNames[i]);
		//		if (!_binding[i])
		//			throw InvalidParameterException("Unknown binding model " + bindModelNames[i]);

		//		MultiplexedScopeSelector scopeGuard(paramProvider, "adsorption", _singleBinding, i, _nParType == 1, _binding[i]->usesParamProviderInDiscretizationConfig());
		//		bindingConfSuccess = _binding[i]->configureModelDiscretization(paramProvider, _nComp, _nBound + i * _nComp, _boundOffset + i * _nComp) && bindingConfSuccess;
		//	}
		//}

		//// ==== Construct and configure dynamic reaction model
		//bool reactionConfSuccess = true;

		//_dynReactionBulk = nullptr;
		//if (paramProvider.exists("REACTION_MODEL"))
		//{
		//	const std::string dynReactName = paramProvider.getString("REACTION_MODEL");
		//	_dynReactionBulk = helper.createDynamicReactionModel(dynReactName);
		//	if (!_dynReactionBulk)
		//		throw InvalidParameterException("Unknown dynamic reaction model " + dynReactName);

		//	if (_dynReactionBulk->usesParamProviderInDiscretizationConfig())
		//		paramProvider.pushScope("reaction_bulk");

		//	reactionConfSuccess = _dynReactionBulk->configureModelDiscretization(paramProvider, _nComp, nullptr, nullptr);

		//	if (_dynReactionBulk->usesParamProviderInDiscretizationConfig())
		//		paramProvider.popScope();
		//}

		//clearDynamicReactionModels();
		//_dynReaction = std::vector<IDynamicReactionModel*>(_nParType, nullptr);

		//if (paramProvider.exists("REACTION_MODEL_PARTICLES"))
		//{
		//	const std::vector<std::string> dynReactModelNames = paramProvider.getStringArray("REACTION_MODEL_PARTICLES");

		//	if (paramProvider.exists("REACTION_MODEL_PARTICLES_MULTIPLEX"))
		//		_singleDynReaction = (paramProvider.getInt("REACTION_MODEL_PARTICLES_MULTIPLEX") == 1);
		//	else
		//	{
		//		// Infer multiplex mode
		//		_singleDynReaction = (dynReactModelNames.size() == 1);
		//	}

		//	if (!_singleDynReaction && (dynReactModelNames.size() < _nParType))
		//		throw InvalidParameterException("Field REACTION_MODEL_PARTICLES contains too few elements (" + std::to_string(_nParType) + " required)");
		//	else if (_singleDynReaction && (dynReactModelNames.size() != 1))
		//		throw InvalidParameterException("Field REACTION_MODEL_PARTICLES requires (only) 1 element");

		//	for (unsigned int i = 0; i < _nParType; ++i)
		//	{
		//		if (_singleDynReaction && (i > 0))
		//		{
		//			// Reuse first binding model
		//			_dynReaction[i] = _dynReaction[0];
		//		}
		//		else
		//		{
		//			_dynReaction[i] = helper.createDynamicReactionModel(dynReactModelNames[i]);
		//			if (!_dynReaction[i])
		//				throw InvalidParameterException("Unknown dynamic reaction model " + dynReactModelNames[i]);

		//			MultiplexedScopeSelector scopeGuard(paramProvider, "reaction_particle", _singleDynReaction, i, _nParType == 1, _dynReaction[i]->usesParamProviderInDiscretizationConfig());
		//			reactionConfSuccess = _dynReaction[i]->configureModelDiscretization(paramProvider, _nComp, _nBound + i * _nComp, _boundOffset + i * _nComp) && reactionConfSuccess;
		//		}
		//	}
		//}

		//// Allocate memory
		//if (firstConfigCall)
		//	_tempState = new double[numDofs()];

		return parSurfDiffDepConfSuccess /*&& bindingConfSuccess*/ /*&& reactionConfSuccess*/;
	}

	bool ParticleDiffusionOperatorDG::configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters)
	{
		const bool firstConfigCall = _offsetSurfDiff == nullptr; // used to not multiply allocate memory

		// Read geometry parameters
		_singleParRadius = readAndRegisterMultiplexTypeParam(paramProvider, parameters, _parRadius, "PAR_RADIUS", _nParType, unitOpIdx);
		_singleParPorosity = readAndRegisterMultiplexTypeParam(paramProvider, parameters, _parPorosity, "PAR_POROSITY", _nParType, unitOpIdx);

		// Let PAR_CORERADIUS default to 0.0 for backwards compatibility
		if (paramProvider.exists("PAR_CORERADIUS"))
			_singleParCoreRadius = readAndRegisterMultiplexTypeParam(paramProvider, parameters, _parCoreRadius, "PAR_CORERADIUS", _nParType, unitOpIdx);
		else
		{
			_singleParCoreRadius = true;
			_parCoreRadius = std::vector<active>(_nParType, 0.0);
		}

		// Check whether PAR_TYPE_VOLFRAC is required or not
		if ((_nParType > 1) && !paramProvider.exists("PAR_TYPE_VOLFRAC"))
			throw InvalidParameterException("The required parameter \"PAR_TYPE_VOLFRAC\" was not found");

		// todo: PAR_TYPE_VOLFRAC remains parameter of the unit, not the particles?
		//// Let PAR_TYPE_VOLFRAC default to 1.0 for backwards compatibility
		//if (paramProvider.exists("PAR_TYPE_VOLFRAC"))
		//{
		//	readScalarParameterOrArray(_parTypeVolFrac, paramProvider, "PAR_TYPE_VOLFRAC", 1);
		//	if (_parTypeVolFrac.size() == _nParType)
		//	{
		//		_axiallyConstantParTypeVolFrac = true;

		//		// Expand to all axial cells
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

		// Check whether all sizes are matched
		if (_nParType != _parRadius.size())
			throw InvalidParameterException("Number of elements in field PAR_RADIUS does not match number of particle types");
		//if (_nParType * nPoints != _parTypeVolFrac.size())
		//	throw InvalidParameterException("Number of elements in field PAR_TYPE_VOLFRAC does not match number of particle types times number of axial cells");
		if (_nParType != _parPorosity.size())
			throw InvalidParameterException("Number of elements in field PAR_POROSITY does not match number of particle types");
		if (_nParType != _parCoreRadius.size())
			throw InvalidParameterException("Number of elements in field PAR_CORERADIUS does not match number of particle types");

		//// Check that particle volume fractions sum to 1.0
		//for (unsigned int i = 0; i < nPoints; ++i)
		//{
		//	const double volFracSum = std::accumulate(_parTypeVolFrac.begin() + i * _nParType, _parTypeVolFrac.begin() + (i + 1) * _nParType, 0.0,
		//		[](double a, const active& b) -> double { return a + static_cast<double>(b); });
		//	if (std::abs(1.0 - volFracSum) > 1e-10)
		//		throw InvalidParameterException("Sum of field PAR_TYPE_VOLFRAC differs from 1.0 (is " + std::to_string(volFracSum) + ") in axial cell " + std::to_string(i));
		//}

		// Read vectorial parameters (which may also be section dependent; transport)
		//_filmDiffusionMode = readAndRegisterMultiplexCompTypeSecParam(paramProvider, parameters, _filmDiffusion, "FILM_DIFFUSION", _nParType, nComp, unitOpIdx); // todo film diffusion remains unit operation parameter?
		_parDiffusionMode = readAndRegisterMultiplexCompTypeSecParam(paramProvider, parameters, _parDiffusion, "PAR_DIFFUSION", _nParType, _nComp, unitOpIdx);

		if (firstConfigCall)
			_offsetSurfDiff = new unsigned int[_strideBound[_nParType]];
		if (paramProvider.exists("PAR_SURFDIFFUSION"))
			_parSurfDiffusionMode = readAndRegisterMultiplexBndCompTypeSecParam(paramProvider, parameters, _parSurfDiffusion, "PAR_SURFDIFFUSION", _nParType, _nComp, _strideBound, _nBound, unitOpIdx);
		else
		{
			_parSurfDiffusionMode = MultiplexMode::Component;
			_parSurfDiffusion.resize(_strideBound[_nParType], 0.0);
		}

		bool parSurfDiffDepConfSuccess = true;
		if (_hasParDepSurfDiffusion)
		{
			if (_singleParDepSurfDiffusion && _parDepSurfDiffusion[0])
			{
				parSurfDiffDepConfSuccess = _parDepSurfDiffusion[0]->configure(paramProvider, unitOpIdx, ParTypeIndep, "PAR_SURFDIFFUSION");
			}
			else if (!_singleParDepSurfDiffusion)
			{
				for (unsigned int i = 0; i < _nParType; ++i)
				{
					if (!_parDepSurfDiffusion[i])
						continue;

					parSurfDiffDepConfSuccess = _parDepSurfDiffusion[i]->configure(paramProvider, unitOpIdx, i, "PAR_SURFDIFFUSION") && parSurfDiffDepConfSuccess;
				}
			}
		}

		//if ((_filmDiffusion.size() < nComp * _nParType) || (_filmDiffusion.size() % (nComp * _nParType) != 0))
		//	throw InvalidParameterException("Number of elements in field FILM_DIFFUSION is not a positive multiple of NCOMP * _nParType (" + std::to_string(nComp * _nParType) + ")");
		if ((_parDiffusion.size() < _nComp * _nParType) || (_parDiffusion.size() % (_nComp * _nParType) != 0))
			throw InvalidParameterException("Number of elements in field PAR_DIFFUSION is not a positive multiple of NCOMP * _nParType (" + std::to_string(_nComp * _nParType) + ")");
		if ((_parSurfDiffusion.size() < _strideBound[_nParType]) || ((_strideBound[_nParType] > 0) && (_parSurfDiffusion.size() % _strideBound[_nParType] != 0)))
			throw InvalidParameterException("Number of elements in field PAR_SURFDIFFUSION is not a positive multiple of NTOTALBND (" + std::to_string(_strideBound[_nParType]) + ")");

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
		//	// Register only the f_Irst _nParType items
		//	for (unsigned int i = 0; i < _nParType; ++i)
		//		parameters[makeParamId(hashString("PAR_TYPE_VOLFRAC"), unitOpIdx, CompIndep, i, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parTypeVolFrac[i];
		//}
		//else
		//	registerParam2DArray(parameters, _parTypeVolFrac, [=](bool multi, unsigned cell, unsigned int type) { return makeParamId(hashString("PAR_TYPE_VOLFRAC"), unitOpIdx, CompIndep, type, BoundStateIndep, ReactionIndep, cell); }, _nParType);

		// Calculate the particle radial discretization variables (_parElementSize, _parCenterRadius, etc.)
		if (firstConfigCall)
			_deltaR = new active[_offsetMetric[_nParType]];
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

		//// Reconfigure binding model
		//bool bindingConfSuccess = true;
		//if (!_binding.empty())
		//{
		//	if (_singleBinding)
		//	{
		//		if (_binding[0] && _binding[0]->requ_IresConfiguration())
		//		{
		//			MultiplexedScopeSelector scopeGuard(paramProvider, "adsorption", true);
		//			bindingConfSuccess = _binding[0]->configure(paramProvider, unitOpIdx, ParTypeIndep);
		//		}
		//	}
		//	else
		//	{
		//		for (unsigned int type = 0; type < _nParType; ++type)
		//		{
		//			if (!_binding[type] || !_binding[type]->requ_IresConfiguration())
		//				continue;

		//			MultiplexedScopeSelector scopeGuard(paramProvider, "adsorption", type, _nParType == 1, true);
		//			bindingConfSuccess = _binding[type]->configure(paramProvider, unitOpIdx, type) && bindingConfSuccess;
		//		}
		//	}
		//}

		//// Reconfigure reaction model
		//bool dynReactionConfSuccess = true;
		//if (_dynReactionBulk && _dynReactionBulk->requ_IresConfiguration())
		//{
		//	paramProvider.pushScope("reaction_bulk");
		//	dynReactionConfSuccess = _dynReactionBulk->configure(paramProvider, unitOpIdx, ParTypeIndep);
		//	paramProvider.popScope();
		//}

		//if (_singleDynReaction)
		//{
		//	if (_dynReaction[0] && _dynReaction[0]->requ_IresConfiguration())
		//	{
		//		MultiplexedScopeSelector scopeGuard(paramProvider, "reaction_particle", true);
		//		dynReactionConfSuccess = _dynReaction[0]->configure(paramProvider, unitOpIdx, ParTypeIndep) && dynReactionConfSuccess;
		//	}
		//}
		//else
		//{
		//	for (unsigned int type = 0; type < _nParType; ++type)
		//	{
		//		if (!_dynReaction[type] || !_dynReaction[type]->requ_IresConfiguration())
		//			continue;

		//		MultiplexedScopeSelector scopeGuard(paramProvider, "reaction_particle", type, _nParType == 1, true);
		//		dynReactionConfSuccess = _dynReaction[type]->configure(paramProvider, unitOpIdx, type) && dynReactionConfSuccess;
		//	}
		//}

		// jaobian pattern set after binding and particle surface diffusion are configured
		//setJacobianPattern_GRM(_globalJac, 0, _dynReactionBulk); // todo set pattern for particles in global unit operation Jacobain

		return parSurfDiffDepConfSuccess/* && bindingConfSuccess && dynReactionConfSuccess*/;
	}

	void ParticleDiffusionOperatorDG::setEquidistantRadialDisc(unsigned int parType)
	{
		active* const ptrCenterRadius = _parCenterRadius.data() + _offsetMetric[parType];
		active* const ptrOuterSurfAreaPerVolume = _parOuterSurfAreaPerVolume.data() + _offsetMetric[parType];
		active* const ptrInnerSurfAreaPerVolume = _parInnerSurfAreaPerVolume.data() + _offsetMetric[parType];

		const active radius = _parRadius[parType] - _parCoreRadius[parType];
		const active _Dr = radius / static_cast<double>(_nParElem[parType]);
		std::fill(_parElementSize.data() + _offsetMetric[parType], _parElementSize.data() + _offsetMetric[parType] + _nParElem[parType], _Dr);

		if (_parGeomSurfToVol[parType] == _SurfVolRatioSphere)
		{
			for (unsigned int cell = 0; cell < _nParElem[parType]; ++cell)
			{
				const active r_out = _parRadius[parType] - static_cast<double>(cell) * _Dr;
				const active r_in = _parRadius[parType] - static_cast<double>(cell + 1) * _Dr;

				ptrCenterRadius[cell] = _parRadius[parType] - (0.5 + static_cast<double>(cell)) * _Dr;

				// Compute denominator -> corresponding to cell volume
				const active vol = pow(r_out, 3.0) - pow(r_in, 3.0);

				ptrOuterSurfAreaPerVolume[cell] = 3.0 * sqr(r_out) / vol;
				ptrInnerSurfAreaPerVolume[cell] = 3.0 * sqr(r_in) / vol;
			}
		}
		else if (_parGeomSurfToVol[parType] == _SurfVolRatioCylinder)
		{
			for (unsigned int cell = 0; cell < _nParElem[parType]; ++cell)
			{
				const active r_out = _parRadius[parType] - static_cast<double>(cell) * _Dr;
				const active r_in = _parRadius[parType] - static_cast<double>(cell + 1) * _Dr;

				ptrCenterRadius[cell] = _parRadius[parType] - (0.5 + static_cast<double>(cell)) * _Dr;

				// Compute denominator -> corresponding to cell volume
				const active vol = sqr(r_out) - sqr(r_in);

				ptrOuterSurfAreaPerVolume[cell] = 2.0 * r_out / vol;
				ptrInnerSurfAreaPerVolume[cell] = 2.0 * r_in / vol;
			}
		}
		else if (_parGeomSurfToVol[parType] == _SurfVolRatioSlab)
		{
			for (unsigned int cell = 0; cell < _nParElem[parType]; ++cell)
			{
				const active r_out = _parRadius[parType] - static_cast<double>(cell) * _Dr;
				const active r_in = _parRadius[parType] - static_cast<double>(cell + 1) * _Dr;

				ptrCenterRadius[cell] = _parRadius[parType] - (0.5 + static_cast<double>(cell)) * _Dr;

				// Compute denominator -> corresponding to cell volume
				const active vol = r_out - r_in;

				ptrOuterSurfAreaPerVolume[cell] = 1.0 / vol;
				ptrInnerSurfAreaPerVolume[cell] = 1.0 / vol;
			}
		}
	}
	/**
	 * @brief Computes the radial nodes in the beads in such a way that all elements have the same volume
	 */
	void ParticleDiffusionOperatorDG::setEquivolumeRadialDisc(unsigned int parType)
	{
		active* const ptrCellSize = _parElementSize.data() + _offsetMetric[parType];
		active* const ptrCenterRadius = _parCenterRadius.data() + _offsetMetric[parType];
		active* const ptrOuterSurfAreaPerVolume = _parOuterSurfAreaPerVolume.data() + _offsetMetric[parType];
		active* const ptrInnerSurfAreaPerVolume = _parInnerSurfAreaPerVolume.data() + _offsetMetric[parType];

		if (_parGeomSurfToVol[parType] == _SurfVolRatioSphere)
		{
			active r_out = _parRadius[parType];
			active r_in = _parCoreRadius[parType];
			const active volumePerelement = (pow(_parRadius[parType], 3.0) - pow(_parCoreRadius[parType], 3.0)) / static_cast<double>(_nParElem[parType]);

			for (unsigned int cell = 0; cell < _nParElem[parType]; ++cell)
			{
				if (cell != (_nParElem[parType] - 1))
					r_in = pow(pow(r_out, 3.0) - volumePerelement, (1.0 / 3.0));
				else
					r_in = _parCoreRadius[parType];

				ptrCellSize[cell] = r_out - r_in;
				ptrCenterRadius[cell] = (r_out + r_in) * 0.5;

				ptrOuterSurfAreaPerVolume[cell] = 3.0 * sqr(r_out) / volumePerelement;
				ptrInnerSurfAreaPerVolume[cell] = 3.0 * sqr(r_in) / volumePerelement;
				// Note that the DG particle elements are oppositely ordered compared to the FV particle elements
				_deltaR[_offsetMetric[parType] + _nParElem[parType] - (cell + 1)] = r_out - r_in;

				// For the next cell: r_out == r_in of the current cell
				r_out = r_in;
			}
		}
		else if (_parGeomSurfToVol[parType] == _SurfVolRatioCylinder)
		{
			active r_out = _parRadius[parType];
			active r_in = _parCoreRadius[parType];
			const active volumePerelement = (sqr(_parRadius[parType]) - sqr(_parCoreRadius[parType])) / static_cast<double>(_nParElem[parType]);

			for (unsigned int cell = 0; cell < _nParElem[parType]; ++cell)
			{
				if (cell != (_nParElem[parType] - 1))
					r_in = sqrt(sqr(r_out) - volumePerelement);
				else
					r_in = _parCoreRadius[parType];

				ptrCellSize[cell] = r_out - r_in;
				ptrCenterRadius[cell] = (r_out + r_in) * 0.5;

				ptrOuterSurfAreaPerVolume[cell] = 2.0 * r_out / volumePerelement;
				ptrInnerSurfAreaPerVolume[cell] = 2.0 * r_in / volumePerelement;
				// Note that the DG particle elements are oppositely ordered compared to the FV particle elements
				_deltaR[_offsetMetric[parType] + _nParElem[parType] - (cell + 1)] = r_out - r_in;

				// For the next cell: r_out == r_in of the current cell
				r_out = r_in;
			}
		}
		else if (_parGeomSurfToVol[parType] == _SurfVolRatioSlab)
		{
			active r_out = _parRadius[parType];
			active r_in = _parCoreRadius[parType];
			const active volumePerelement = (_parRadius[parType] - _parCoreRadius[parType]) / static_cast<double>(_nParElem[parType]);

			for (unsigned int cell = 0; cell < _nParElem[parType]; ++cell)
			{
				if (cell != (_nParElem[parType] - 1))
					r_in = r_out - volumePerelement;
				else
					r_in = _parCoreRadius[parType];

				ptrCellSize[cell] = r_out - r_in;
				ptrCenterRadius[cell] = (r_out + r_in) * 0.5;

				ptrOuterSurfAreaPerVolume[cell] = 1.0 / volumePerelement;
				ptrInnerSurfAreaPerVolume[cell] = 1.0 / volumePerelement;
				// Note that the DG particle elements are oppositely ordered compared to the FV particle elements
				_deltaR[_offsetMetric[parType] + _nParElem[parType] - (cell + 1)] = r_out - r_in;

				// For the next cell: r_out == r_in of the current cell
				r_out = r_in;
			}
		}
	}

	/**
	 * @brief Computes all helper quantities for radial bead discretization from given radial cell boundaries
	 * @details Calculates surface areas per volume for every element and the radial element centers.
	 */
	void ParticleDiffusionOperatorDG::setUserdefinedRadialDisc(unsigned int parType)
	{
		active* const ptrCellSize = _parElementSize.data() + _offsetMetric[parType];
		active* const ptrCenterRadius = _parCenterRadius.data() + _offsetMetric[parType];
		active* const ptrOuterSurfAreaPerVolume = _parOuterSurfAreaPerVolume.data() + _offsetMetric[parType];
		active* const ptrInnerSurfAreaPerVolume = _parInnerSurfAreaPerVolume.data() + _offsetMetric[parType];

		// Care for the right ordering and include 0.0 / 1.0 if not already in the vector.
		std::vector<active> orderedInterfaces = std::vector<active>(_parDiscVector.begin() + _offsetMetric[parType] + parType,
			_parDiscVector.begin() + _offsetMetric[parType] + parType + _nParElem[parType] + 1);

		// Sort in descending order
		std::sort(orderedInterfaces.begin(), orderedInterfaces.end(), std::greater<active>());

		// Force f_Irst and last element to be 1.0 and 0.0, respectively
		orderedInterfaces[0] = 1.0;
		orderedInterfaces.back() = 0.0;

		// Map [0, 1] -> [core radius, particle radius] via linear interpolation
		for (unsigned int cell = 0; cell < _nParElem[parType]; ++cell)
			orderedInterfaces[cell] = static_cast<double>(orderedInterfaces[cell]) * (_parRadius[parType] - _parCoreRadius[parType]) + _parCoreRadius[parType];

		if (_parGeomSurfToVol[parType] == _SurfVolRatioSphere)
		{
			for (unsigned int cell = 0; cell < _nParElem[parType]; ++cell)
			{
				ptrCellSize[cell] = orderedInterfaces[cell] - orderedInterfaces[cell + 1];
				ptrCenterRadius[cell] = (orderedInterfaces[cell] + orderedInterfaces[cell + 1]) * 0.5;

				// Compute denominator -> corresponding to cell volume
				const active vol = pow(orderedInterfaces[cell], 3.0) - pow(orderedInterfaces[cell + 1], 3.0);

				ptrOuterSurfAreaPerVolume[cell] = 3.0 * sqr(orderedInterfaces[cell]) / vol;
				ptrInnerSurfAreaPerVolume[cell] = 3.0 * sqr(orderedInterfaces[cell + 1]) / vol;
				// Note that the DG particle elements are oppositely ordered compared to the FV particle elements
				_deltaR[_offsetMetric[parType] + _nParElem[parType] - (cell + 1)] = ptrOuterSurfAreaPerVolume[cell] - ptrInnerSurfAreaPerVolume[cell];
			}
		}
		else if (_parGeomSurfToVol[parType] == _SurfVolRatioCylinder)
		{
			for (unsigned int cell = 0; cell < _nParElem[parType]; ++cell)
			{
				ptrCellSize[cell] = orderedInterfaces[cell] - orderedInterfaces[cell + 1];
				ptrCenterRadius[cell] = (orderedInterfaces[cell] + orderedInterfaces[cell + 1]) * 0.5;

				// Compute denominator -> corresponding to cell volume
				const active vol = sqr(orderedInterfaces[cell]) - sqr(orderedInterfaces[cell + 1]);

				ptrOuterSurfAreaPerVolume[cell] = 2.0 * orderedInterfaces[cell] / vol;
				ptrInnerSurfAreaPerVolume[cell] = 2.0 * orderedInterfaces[cell + 1] / vol;
				// Note that the DG particle elements are oppositely ordered compared to the FV particle elements
				_deltaR[_offsetMetric[parType] + _nParElem[parType] - (cell + 1)] = ptrOuterSurfAreaPerVolume[cell] - ptrInnerSurfAreaPerVolume[cell];
			}
		}
		else if (_parGeomSurfToVol[parType] == _SurfVolRatioSlab)
		{
			for (unsigned int cell = 0; cell < _nParElem[parType]; ++cell)
			{
				ptrCellSize[cell] = orderedInterfaces[cell] - orderedInterfaces[cell + 1];
				ptrCenterRadius[cell] = (orderedInterfaces[cell] + orderedInterfaces[cell + 1]) * 0.5;

				// Compute denominator -> corresponding to cell volume
				const active vol = orderedInterfaces[cell] - orderedInterfaces[cell + 1];

				ptrOuterSurfAreaPerVolume[cell] = 1.0 / vol;
				ptrInnerSurfAreaPerVolume[cell] = 1.0 / vol;
				// Note that the DG particle elements are oppositely ordered compared to the FV particle elements
				_deltaR[_offsetMetric[parType] + _nParElem[parType] - (cell + 1)] = ptrOuterSurfAreaPerVolume[cell] - ptrInnerSurfAreaPerVolume[cell];
			}
		}
	}

	// todo: parameter sensitivities for particle radius. Here, we have the problem that every DG operator becomes an active type.
	// alternatively (only for exact integration), we could store more matrices and compute the metric dependend calculations in the residual.
	// inexact integration approach is deprecated anyways but would requ_Ire active type DG operators since every entry is multiplied by an individual metric term,
	// wherease for the exact integration approach we have sums of three matrices each multiplied by its own metric term, whcih could be applied iteratively to the residual/solution.
	// Not needed for Slab, and only two matrices (exact integration approach) requ_Ired for Cylinder.
	// This approach should only be used when necessary, i.e. solely when particle radius parameter sensitivity is requ_Ired.
	void ParticleDiffusionOperatorDG::updateRadialDisc()
	{

		for (unsigned int parType = 0; parType < _nParType; ++parType)
		{
			if (_parDiscMode[parType] == ParticleDiscretizationMode::Equidistant)
			{
				for (int cell = 0; cell < _nParElem[parType]; cell++)
				{
					_deltaR[_offsetMetric[parType] + cell] = (_parRadius[parType] - _parCoreRadius[parType]) / _nParElem[parType];
				}
				setEquidistantRadialDisc(parType);
			}
			else if (_parDiscMode[parType] == ParticleDiscretizationMode::Equivolume)
				setEquivolumeRadialDisc(parType);
			else if (_parDiscMode[parType] == ParticleDiscretizationMode::UserDefined)
				setUserdefinedRadialDisc(parType);
		}

		/*		metrics		*/
		// estimate cell dependent D_r

		for (int parType = 0; parType < _nParType; parType++)
		{
			for (int cell = 0; cell < _nParElem[parType]; cell++)
			{
				for (int node = 0; node < _nParNode[parType]; node++)
					_Ir[_offsetMetric[parType] + cell][node] = _deltaR[_offsetMetric[parType] + cell] / 2.0 * (_parNodes[parType][node] + 1.0);

				active r_L = _parCoreRadius[parType] + cell * _deltaR[_offsetMetric[parType] + cell]; // left boundary of current cell

				_Ir[_offsetMetric[parType] + cell] = _Ir[_offsetMetric[parType] + cell] + VectorXd::Ones(_nParNode[parType]) * r_L;

				if (_parGeomSurfToVol[parType] == _SurfVolRatioSphere)
					_Ir[_offsetMetric[parType] + cell] = _Ir[_offsetMetric[parType] + cell].array().square();
				else if (_parGeomSurfToVol[parType] == _SurfVolRatioSlab)
					_Ir[_offsetMetric[parType] + cell].setOnes(); // no metric terms for slab

				// (D_r)_{i, j} = D_{i, j} * (r_j / r_i) [only needed for inexact integration]
				_Dr[_offsetMetric[parType] + cell] = _parPolyDerM[parType];
				_Dr[_offsetMetric[parType] + cell].array().rowwise() *= _Ir[_offsetMetric[parType] + cell].array().template cast<double>().transpose();
				_Dr[_offsetMetric[parType] + cell].array().colwise() *= _Ir[_offsetMetric[parType] + cell].array().template cast<double>().cwiseInverse();

				// compute mass matrices for exact integration based on particle geometry, via transformation to normalized Jacobi polynomials with weight function w
				if (_parGeomSurfToVol[parType] == _SurfVolRatioSphere) // r^2 =  r_i^2 + (1 + \xi) * r_i * DeltaR_i / 2.0 + (1 + \xi)^2 * (DeltaR_i / 2.0)^2
				{
					_parInvMM[_offsetMetric[parType] + cell] = parts::dgtoolbox::mMatrix(_parPolyDeg[parType], _parNodes[parType], 0.0, 2.0) * pow((static_cast<double>(_deltaR[_offsetMetric[parType] + cell]) / 2.0), 2.0);
					if (cell > 0 || _parCoreRadius[parType] != 0.0) // following contributions are zero for f_Irst cell when R_c = 0 (no particle core)
						_parInvMM[_offsetMetric[parType] + cell] += parts::dgtoolbox::mMatrix(_parPolyDeg[parType], _parNodes[parType], 0.0, 1.0) * (static_cast<double>(_deltaR[_offsetMetric[parType] + cell]) * static_cast<double>(r_L))
						+ parts::dgtoolbox::mMatrix(_parPolyDeg[parType], _parNodes[parType], 0.0, 0.0) * pow(static_cast<double>(r_L), 2.0);

					_parInvMM[_offsetMetric[parType] + cell] = _parInvMM[_offsetMetric[parType] + cell].inverse();
					_minus_InvMM_ST[_offsetMetric[parType] + cell] = -_parInvMM[_offsetMetric[parType] + cell] * _parPolyDerM[parType].transpose() * _parInvMM[_offsetMetric[parType] + cell].inverse();

					// particle GSM specific second order stiffness matrix (single element, i.e. _nParElem = 1)
					_secondOrderStiffnessM[parType] = std::pow(static_cast<double>(_parCoreRadius[parType]), 2.0) * parts::dgtoolbox::secondOrderStiffnessMatrix(_parPolyDeg[parType], 0.0, 0.0, _parNodes[parType]);
					_secondOrderStiffnessM[parType] += static_cast<double>(_deltaR[_offsetMetric[parType]]) * static_cast<double>(_parCoreRadius[parType]) * parts::dgtoolbox::secondOrderStiffnessMatrix(_parPolyDeg[parType], 0.0, 1.0, _parNodes[parType]);
					_secondOrderStiffnessM[parType] += std::pow(static_cast<double>(_deltaR[_offsetMetric[parType]]) / 2.0, 2.0) * parts::dgtoolbox::secondOrderStiffnessMatrix(_parPolyDeg[parType], 0.0, 2.0, _parNodes[parType]);
				}
				else if (_parGeomSurfToVol[parType] == _SurfVolRatioCylinder) // r = r_i + (1 + \xi) * DeltaR_i / 2.0
				{
					_parInvMM[_offsetMetric[parType] + cell] = parts::dgtoolbox::mMatrix(_parPolyDeg[parType], _parNodes[parType], 0.0, 1.0) * (static_cast<double>(_deltaR[_offsetMetric[parType] + cell]) / 2.0);
					if (cell > 0 || _parCoreRadius[parType] != 0.0) // following contribution is zero for f_Irst cell when R_c = 0 (no particle core)
						_parInvMM[_offsetMetric[parType] + cell] += parts::dgtoolbox::mMatrix(_parPolyDeg[parType], _parNodes[parType], 0.0, 0.0) * static_cast<double>(r_L);

					_parInvMM[_offsetMetric[parType] + cell] = _parInvMM[_offsetMetric[parType] + cell].inverse();
					_minus_InvMM_ST[_offsetMetric[parType] + cell] = -_parInvMM[_offsetMetric[parType] + cell] * _parPolyDerM[parType].transpose() * _parInvMM[_offsetMetric[parType] + cell].inverse();

					// particle GSM specific second order stiffness matrix (single element, i.e. _nParElem = 1)
					_secondOrderStiffnessM[parType] = static_cast<double>(_parCoreRadius[parType]) * parts::dgtoolbox::secondOrderStiffnessMatrix(_parPolyDeg[parType], 0.0, 0.0, _parNodes[parType]);
					_secondOrderStiffnessM[parType] += static_cast<double>(_deltaR[_offsetMetric[parType]]) / 2.0 * parts::dgtoolbox::secondOrderStiffnessMatrix(_parPolyDeg[parType], 0.0, 1.0, _parNodes[parType]);
				}
				else if (_parGeomSurfToVol[parType] == _SurfVolRatioSlab) // r = 1
				{
					_minus_InvMM_ST[_offsetMetric[parType] + cell] = -_parInvMM[_offsetMetric[parType] + cell] * _parPolyDerM[parType].transpose() * _parInvMM[_offsetMetric[parType] + cell].inverse();

					_secondOrderStiffnessM[parType] = parts::dgtoolbox::secondOrderStiffnessMatrix(_parPolyDeg[parType], 0.0, 0.0, _parNodes[parType]);
					_parInvMM[_offsetMetric[parType] + cell] = parts::dgtoolbox::invMMatrix(_parPolyDeg[parType], _parNodes[parType], 0.0, 0.0);
				}
			}

			_minus_parInvMM_Ar[parType] = -_parInvMM[_offsetMetric[parType]] * _secondOrderStiffnessM[parType];
		}
	}

	/**
	 * @brief allocates memory for DG operators and computes those that are metric independent. Also allocates requ_Ired containers needed for the DG discretization.
	 */
	void ParticleDiffusionOperatorDG::initializeDG()
	{
		/* Allocate space for DG operators and containers */

		const bool firstConfigCall = _nParNode == nullptr; // used to not multiply allocate memory

		// particles
		if (firstConfigCall)
		{
			_nParNode = new unsigned int[_nParType];
			_nParPoints = new unsigned int[_nParType];
			_g_p = new Vector<active, Dynamic>[_nParType];
			_g_pSum = new Vector<active, Dynamic>[_nParType];
			_surfaceFluxParticle = new Vector<active, Dynamic>[_nParType];
			_parNodes = new VectorXd[_nParType];
			_parInvWeights = new VectorXd[_nParType];
			_parInvMM_Leg = new MatrixXd[_nParType];
			_parPolyDerM = new MatrixXd[_nParType];
			_localFlux = new active[_nComp];
		}

		for (int parType = 0; parType < _nParType; parType++)
		{
			_nParNode[parType] = _parPolyDeg[parType] + 1u;
			_nParPoints[parType] = _nParNode[parType] * _nParElem[parType];
			_g_p[parType].resize(_nParPoints[parType]);
			_g_p[parType].setZero();
			_g_pSum[parType].resize(_nParPoints[parType]);
			_g_pSum[parType].setZero();
			_surfaceFluxParticle[parType].resize(_nParElem[parType] + 1);
			_surfaceFluxParticle[parType].setZero();
			_parNodes[parType].resize(_nParNode[parType]);
			_parNodes[parType].setZero();
			_parInvWeights[parType].resize(_nParNode[parType]);
			_parInvWeights[parType].setZero();
			_parPolyDerM[parType].resize(_nParNode[parType], _nParNode[parType]);
			_parPolyDerM[parType].setZero();
			_parInvMM_Leg[parType].resize(_nParNode[parType], _nParNode[parType]);
			_parInvMM_Leg[parType].setZero();
		}

		_offsetMetric = VectorXi::Zero(_nParType + 1);
		for (int parType = 1; parType <= _nParType; parType++) {
			_offsetMetric[parType] = _offsetMetric[parType - 1] + _nParElem[parType - 1];
		}

		if (firstConfigCall)
		{
			_Dr = new MatrixXd[_offsetMetric[_nParType]];
			_Ir = new Vector<active, Dynamic>[_offsetMetric[_nParType]];
			for (int parType = 0; parType < _nParType; parType++)
			{
				for (int elem = 0; elem < _nParElem[parType]; elem++)
				{
					_Ir[_offsetMetric[parType] + elem].resize(_nParNode[parType]);
					_Ir[_offsetMetric[parType] + elem].setZero();
					_Dr[_offsetMetric[parType] + elem].resize(_nParNode[parType], _nParNode[parType]);
					_Dr[_offsetMetric[parType] + elem].setZero();
				}
			}
			_minus_InvMM_ST = new MatrixXd[_offsetMetric[_nParType]];
			_parInvMM = new MatrixXd[_offsetMetric[_nParType]];
			_secondOrderStiffnessM = new MatrixXd[_nParType];
			_minus_parInvMM_Ar = new MatrixXd[_nParType];
		}

		/* compute metric independent DG operators for bulk and particles. Note that metric dependent DG operators are computet in updateRadialDisc(). */

		for (int parType = 0; parType < _nParType; parType++)
		{
			parts::dgtoolbox::lglNodesWeights(_parPolyDeg[parType], _parNodes[parType], _parInvWeights[parType], true);
			_parPolyDerM[parType] = parts::dgtoolbox::derivativeMatrix(_parPolyDeg[parType], _parNodes[parType]);
			_parInvMM_Leg[parType] = parts::dgtoolbox::invMMatrix(_parPolyDeg[parType], _parNodes[parType], 0.0, 0.0);
		}
	}


	/**
	 * @brief Notifies the operator that a discontinuous section transition is in progress
	 * @param [in] t Current time point
	 * @param [in] secIdx Index of the new section that is about to be integrated
	 * @return @c true if flow direction has changed, otherwise @c false
	 */
	bool ParticleDiffusionOperatorDG::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active const* const filmDiff, active const* const poreAccessFactor)
	{
		// calculate offsets between surface diffusion storage and state vector order
		orderSurfDiff();

		for (int fd = 0; fd < _nComp * _nParType; fd++)
		{
			_filmDiffusion[fd] = filmDiff[fd];
			_poreAccessFactor[fd] = poreAccessFactor[fd];
		}

		//_curSection = secIdx;
		//_newStaticJac = true;

		// todo update operators and Jacobian blocks

		return true;
	}

	/**
	 * @brief Computes the residual of the transport equations
	 * @param [in] model Model that owns the operator
	 * @param [in] t Current time point
	 * @param [in] secIdx Index of the current section
	 * @param [in] y Pointer to unit operation's state vector
	 * @param [in] yDot Pointer to unit operation's time derivative state vector
	 * @param [out] res Pointer to unit operation's residual vector
	 * @param [in] jac Matrix that holds the Jacobian
	 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
	 */
	int ParticleDiffusionOperatorDG::residual(double t, unsigned int parType, unsigned int colNode, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, int const* const qsBinding, WithoutParamSensitivity)
	{
		return residualImpl<double, double, double>(t, parType, colNode, secIdx, yPar, yBulk, yDotPar, resPar, qsBinding);
	}

	int ParticleDiffusionOperatorDG::residual(double t, unsigned int parType, unsigned int colNode, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, int const* const qsBinding, WithoutParamSensitivity)
	{
		return residualImpl<active, active, double>(t, parType, colNode, secIdx, yPar, yBulk, yDotPar, resPar, qsBinding);
	}

	int ParticleDiffusionOperatorDG::residual(double t, unsigned int parType, unsigned int colNode, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, int const* const qsBinding, WithParamSensitivity)
	{
		return residualImpl<double, active, active>(t, parType, colNode, secIdx, yPar, yBulk, yDotPar, resPar, qsBinding);
	}

	int ParticleDiffusionOperatorDG::residual(double t, unsigned int parType, unsigned int colNode, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, int const* const qsBinding, WithParamSensitivity)
	{
		return residualImpl<active, active, active>(t, parType, colNode, secIdx, yPar, yBulk, yDotPar, resPar, qsBinding);
	}

	template <typename StateType, typename ResidualType, typename ParamType>
	int ParticleDiffusionOperatorDG::residualImpl(double t, unsigned int parType, unsigned int colNode, unsigned int secIdx, StateType const* yPar, StateType const* yBulk, double const* yDotPar, ResidualType* resPar, int const* const qsBinding)
	{
		/* Mobile phase RHS	*/

		// Get film diffusion flux at current node to compute boundary condition
		for (unsigned int comp = 0; comp < _nComp; comp++) {
			_localFlux[comp] = _filmDiffusion[comp] * (yBulk[comp * _strideBulkComp] - yPar[comp]);
		}

		active const* const parDiff = getSectionDependentSlice(_parDiffusion, _nComp * _nParType, secIdx) + parType * _nComp;

		// Ordering of particle surface diffusion:
		// bnd0comp0, bnd0comp1, bnd0comp2, bnd1comp0, bnd1comp1, bnd1comp2
		active const* const parSurfDiff = getSectionDependentSlice(_parSurfDiffusion, _strideBound[_nParType], secIdx) + _nBoundBeforeType[parType];

		const int nNodes = _nParNode[parType];
		const int nElem = _nParElem[parType];
		const int nPoints = _nParPoints[parType];
		const int nComp = _nComp;

		if (_parGSM[parType]) // GSM implementation
		{
			for (unsigned int comp = 0; comp < nComp; comp++)
			{

				// ====================================================================================//
				// solve GSM-discretized particle mass balance   									   //
				// ====================================================================================//

				Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> resCp(resPar + comp, nPoints, InnerStride<Dynamic>(strideParNode(parType)));
				Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> Cp(yPar + comp, nPoints, InnerStride<Dynamic>(strideParNode(parType)));

				// Use auxiliary variable to get c^p + \sum 1 / \Beta_p c^s
				Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<>> sum_cp_cs(reinterpret_cast<ResidualType*>(&_g_pSum[parType][0]), nPoints, InnerStride<>(1));
				sum_cp_cs = static_cast<ParamType>(parDiff[comp]) * Cp.template cast<ResidualType>();

				for (int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++)
				{
					if (parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)] != 0.0) // some bound states might still not be effected by surface diffusion
					{
						Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> c_s(yPar + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd, nPoints, InnerStride<Dynamic>(strideParNode(parType)));
						ParamType invBetaP = (1.0 - static_cast<ParamType>(_parPorosity[parType])) / (static_cast<ParamType>(_poreAccessFactor[_nComp * parType + comp]) * static_cast<ParamType>(_parPorosity[parType]));
						sum_cp_cs += invBetaP * static_cast<ParamType>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)]) * c_s;

						/* For kinetic bindings with surface diffusion: add the additional DG-discretized particle mass balance equations to residual */

						if (!qsBinding[bnd])
						{
							// Eigen access to current bound state residual
							Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> resCs(resPar + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd,
								nPoints, InnerStride<Dynamic>(strideParNode(parType)));

							// Use auxiliary variable to get \Beta_p D_s c^s
							Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<>> c_s_modified(reinterpret_cast<ResidualType*>(&_g_p[parType][0]), nPoints, InnerStride<>(1));

							// Apply squared inverse mapping and surface diffusion
							c_s_modified = 2.0 / static_cast<ParamType>(_deltaR[_offsetMetric[parType]]) * 2.0 / static_cast<ParamType>(_deltaR[_offsetMetric[parType]]) *
								static_cast<ParamType>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)]) * c_s;

							Eigen::Map<const Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> c_s_modified_const(&c_s_modified[0], nPoints, InnerStride<Dynamic>(1));
							parGSMVolumeIntegral<ResidualType, ResidualType>(parType, c_s_modified_const, resCs);

							// Leave out the surface integral as we only have one element, i.e. we apply BC with zeros
						}
					}
				}

				// Apply squared inverse mapping
				sum_cp_cs *= 2.0 / static_cast<ParamType>(_deltaR[_offsetMetric[parType]]) * 2.0 / static_cast<ParamType>(_deltaR[_offsetMetric[parType]]);

				Eigen::Map<const Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> sum_cp_cs_const(&sum_cp_cs[0], nPoints, InnerStride<Dynamic>(1));
				parGSMVolumeIntegral<ResidualType, ResidualType>(parType, sum_cp_cs_const, resCp);

				// Pass sum_cp_cs_const to match the DGSEM interface; nullptr might also be feasible
				parSurfaceIntegral<ResidualType, ResidualType>(parType, sum_cp_cs_const, resCp, nNodes, 1u, false, comp);

			}
		}
		else // DGSEM implementation
		{
			for (unsigned int comp = 0; comp < nComp; comp++)
			{
				// Component dependent (through access factor) inverse Beta_P
				ParamType invBetaP = (1.0 - static_cast<ParamType>(_parPorosity[parType])) / (static_cast<ParamType>(_poreAccessFactor[_nComp * parType + comp]) * static_cast<ParamType>(_parPorosity[parType]));

				// =====================================================================================================//
				// Solve auxiliary systems  d_p g_p + d_s beta_p sum g_s= d (d_p c_p + d_s beta_p sum c_s) / d xi		//
				// =====================================================================================================//
				
				// Component-wise! strides
				unsigned int strideCell = nNodes;
				unsigned int strideNode = 1u;

				// Reset cache for auxiliary variable
				Eigen::Map<Vector<StateType, Dynamic>, 0, InnerStride<>> g_p(reinterpret_cast<StateType*>(&_g_p[parType][0]), nPoints, InnerStride<>(1));
				Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<>> g_pSum(reinterpret_cast<ResidualType*>(&_g_pSum[parType][0]), nPoints, InnerStride<>(1));
				g_p.setZero();
				g_pSum.setZero();

				Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> cp(yPar + comp, _nParPoints[parType], InnerStride<Dynamic>(strideParNode(parType)));

				// Handle surface diffusion: Compute auxiliary variable; For kinetic bindings: add additional mass balance to residual of respective bound state
				if (_hasSurfaceDiffusion[parType])
				{
					for (int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++)
					{
						if (parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)] != 0.0) // some bound states might still not be effected by surface diffusion
						{
							// Get solid phase vector
							Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> c_s(yPar + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd,
								nPoints, InnerStride<Dynamic>(strideParNode(parType)));
							// Compute g_s = d c_s / d xi
							solve_auxiliary_DG<StateType>(parType, c_s, strideCell, strideNode, comp);
							// Apply invBeta_p, d_s and add to sum -> gSum += d_s * invBeta_p * (D c - M^-1 B [c - c^*])
							g_pSum += g_p.template cast<ResidualType>() * invBetaP * static_cast<ParamType>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)]);

							/* For kinetic bindings with surface diffusion: add the additional DG-discretized particle mass balance equations to residual */

							if (!qsBinding[bnd])
							{
								// Eigen access to current bound state residual
								Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> resCs(resPar + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd,
									nPoints, InnerStride<Dynamic>(strideParNode(parType)));

								// Promote auxiliary variable storage from double to active if required
								// @todo is there a more efficient or elegant solution?
								if (std::is_same<ResidualType, active>::value && std::is_same<StateType, double>::value)
									vectorPromoter(reinterpret_cast<double*>(&g_p[0]), nPoints); // reinterpret_cast only required because statement is scanned also when StateType != double

								// Access auxiliary variable as ResidualType
								Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<>> g_p_ResType(reinterpret_cast<ResidualType*>(&_g_p[parType][0]), nPoints, InnerStride<>(1));

								applyParInvMap<ResidualType, ParamType>(g_p_ResType, parType);
								g_p_ResType *= static_cast<ParamType>(parSurfDiff[getOffsetSurfDiff(parType, comp, bnd)]);

								// Eigen access to auxiliary variable of current bound state
								Eigen::Map<const Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> g_p_ResType_const(&g_p_ResType[0], nPoints, InnerStride<Dynamic>(1));

								// Add - D_r * gs to the residual, including metric part.->res = invMap^2* [ -D_r * (d_s c^s) ]
								parVolumeIntegral<ResidualType, ResidualType>(parType, false, g_p_ResType_const, resCs);

								// Add M^-1 B (gs - gs^*) to the residual -> res =  invMap^2 * [ - D_r * (d_s c^s) + M^-1 B (gs - gs^*) ]
								parSurfaceIntegral<ResidualType, ResidualType>(parType, g_p_ResType_const, resCs, strideCell, strideNode, false, comp, true);
							}
						}
					}
				}

				// Compute g_p = d c_p / d xi
				solve_auxiliary_DG<StateType>(parType, cp, strideCell, strideNode, comp);

				// Add particle diffusion part to auxiliary variable sum -> gSum += d_p * (D c - M^-1 B [c - c^*])
				g_pSum += g_p * static_cast<ParamType>(parDiff[comp]);

				// apply squared inverse mapping to sum of bound state auxiliary variables -> gSum = - invMap^2 * (d_p * c^p + sum_mi d_s invBeta_p c^s)
				applyParInvMap<ResidualType, ParamType>(g_pSum, parType);

				// ====================================================================================//
				// solve DG-discretized particle mass balance   									     //
				// ====================================================================================//

				  /* Solve DG-discretized particle mass balance equation */

				Eigen::Map<const Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> g_pSum_const(&g_pSum[0], nPoints, InnerStride<Dynamic>(1));

				// Eigen access to particle liquid residual
				Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> resCp(resPar + comp, nPoints, InnerStride<Dynamic>(strideParNode(parType)));

				// Add - D_r * (g_sum) to the residual, including metric part. -> res = - D_r * (d_p * c^p + invBeta_p sum_mi d_s c^s)
				parVolumeIntegral<ResidualType, ResidualType>(parType, false, g_pSum_const, resCp);

				// Add M^-1 B (g_sum - g_sum^*) to the residual -> res = - D_r * (d_p * c^p + invBeta_p sum_mi d_s c^s) + M^-1 B (g_sum - g_sum^*)
				parSurfaceIntegral<ResidualType, ResidualType>(parType, g_pSum_const, resCp, strideCell, strideNode, false, comp);

			}
		}

		return true;
	}

	/**
	 * @brief calculates the DG Jacobian auxiliary block
	 * @param [in] exInt true if exact integration DG scheme
	 * @param [in] elemIdx elem index
	 */
	MatrixXd ParticleDiffusionOperatorDG::getParGBlock(unsigned int elemIdx, unsigned int parType)
	{
		// Auxiliary Block [ d g(c) / d c ], additionally depends on boundary entries of neighbouring elements
		MatrixXd gBlock = MatrixXd::Zero(_nParNode[parType], _nParNode[parType] + 2);
		gBlock.block(0, 1, _nParNode[parType], _nParNode[parType]) = _parPolyDerM[parType];
		if (_parExactInt[parType]) {
			if (elemIdx == 0 || elemIdx == _nParElem[parType] + 1) { // elemIdx out of bounds
				return MatrixXd::Zero(_nParNode[parType], _nParNode[parType] + 2);
			}
			if (elemIdx != 1 && elemIdx != _nParElem[parType]) { // inner elem
				gBlock.block(0, 0, _nParNode[parType], 1) -= 0.5 * _parInvMM_Leg[parType].block(0, 0, _nParNode[parType], 1);
				gBlock.block(0, 1, _nParNode[parType], 1) += 0.5 * _parInvMM_Leg[parType].block(0, 0, _nParNode[parType], 1);
				gBlock.block(0, _nParNode[parType], _nParNode[parType], 1) -= 0.5 * _parInvMM_Leg[parType].block(0, _nParNode[parType] - 1, _nParNode[parType], 1);
				gBlock.block(0, _nParNode[parType] + 1, _nParNode[parType], 1) += 0.5 * _parInvMM_Leg[parType].block(0, _nParNode[parType] - 1, _nParNode[parType], 1);
			}
			else if (elemIdx == 1u) { // left boundary elem
				if (elemIdx == _nParElem[parType]) // special case one elem
					return gBlock * 2.0 / static_cast<double>(_deltaR[_offsetMetric[parType] + (elemIdx - 1)]);
				gBlock.block(0, _nParNode[parType], _nParNode[parType], 1) -= 0.5 * _parInvMM_Leg[parType].block(0, _nParNode[parType] - 1, _nParNode[parType], 1);
				gBlock.block(0, _nParNode[parType] + 1, _nParNode[parType], 1) += 0.5 * _parInvMM_Leg[parType].block(0, _nParNode[parType] - 1, _nParNode[parType], 1);
			}
			else if (elemIdx == _nParElem[parType]) { // right boundary elem
				gBlock.block(0, 0, _nParNode[parType], 1) -= 0.5 * _parInvMM_Leg[parType].block(0, 0, _nParNode[parType], 1);
				gBlock.block(0, 1, _nParNode[parType], 1) += 0.5 * _parInvMM_Leg[parType].block(0, 0, _nParNode[parType], 1);
			}
			gBlock *= 2.0 / static_cast<double>(_deltaR[_offsetMetric[parType] + (elemIdx - 1)]);
		}
		else {
			// inexact integration not maintained due to inferior performance. Code is part of calcParticleCollocationDGSEMJacobian()
		}

		return gBlock;
	}

	/**
	 * @brief calculates the num. flux part of a dispersion DG Jacobian block
	 * @param [in] elemIdx elem index
	 * @param [in] leftG left neighbour auxiliary block
	 * @param [in] middleG neighbour auxiliary block
	 * @param [in] rightG neighbour auxiliary block
	 */
	MatrixXd ParticleDiffusionOperatorDG::parAuxBlockGstar(unsigned int elemIdx, unsigned int parType, MatrixXd leftG, MatrixXd middleG, MatrixXd rightG) {

		// auxiliary block [ d g^* / d c ], depends on whole previous and subsequent elem plus f_Irst entries of subsubsequent elements
		MatrixXd gStarDC = MatrixXd::Zero(_nParNode[parType], 3 * _nParNode[parType] + 2);
		// NOTE: N = polyDeg
		// indices  gStarDC    :     0   ,   1   , ..., nNodes; nNodes+1, ..., 2 * nNodes;	2*nNodes+1, ..., 3 * nNodes; 3*nNodes+1
		// derivative index j  : -(N+1)-1, -(N+1),... ,  -1   ;   0     , ...,		N	 ;	  N + 1	  , ..., 2N + 2    ; 2(N+1) +1
		// auxiliary block [d g^* / d c]
		if (elemIdx != 1) {
			gStarDC.block(0, _nParNode[parType], 1, _nParNode[parType] + 2) += middleG.block(0, 0, 1, _nParNode[parType] + 2);
			gStarDC.block(0, 0, 1, _nParNode[parType] + 2) += leftG.block(_nParNode[parType] - 1, 0, 1, _nParNode[parType] + 2);
		}
		if (elemIdx != _nParElem[parType]) {
			gStarDC.block(_nParNode[parType] - 1, _nParNode[parType], 1, _nParNode[parType] + 2) += middleG.block(_nParNode[parType] - 1, 0, 1, _nParNode[parType] + 2);
			gStarDC.block(_nParNode[parType] - 1, 2 * _nParNode[parType], 1, _nParNode[parType] + 2) += rightG.block(0, 0, 1, _nParNode[parType] + 2);
		}
		gStarDC *= 0.5;

		return gStarDC;
	}

	/**
	 * @brief calculates the lifting matrix B
	 * @param [in] parType particle type index
	 * @param [in] elem element index
	 * @param [in] parGeomSurfToVol particle geometry
	 */
	MatrixXd ParticleDiffusionOperatorDG::getParBMatrix(int parType, int elem, double parGeomSurfToVol) {
		// also known as "lifting" matrix and includes metric dependent terms for particle discretization
		MatrixXd B = MatrixXd::Zero(_nParNode[parType], _nParNode[parType]);
		if (parGeomSurfToVol == _SurfVolRatioSlab) {
			B(0, 0) = -1.0;
			B(_nParNode[parType] - 1, _nParNode[parType] - 1) = 1.0;
		}
		else {
			B(0, 0) = -static_cast<double>(_Ir[_offsetMetric[parType] + (elem - 1)][0]);
			B(_nParNode[parType] - 1, _nParNode[parType] - 1) = static_cast<double>(_Ir[_offsetMetric[parType] + (elem - 1)][_nParNode[parType] - 1]);
		}

		return B;
	}

	/**
	 * @brief calculates the dispersion part of the DG jacobian
	 * @param [in] parType particle type index
	 * @param [in] parGeomSurfToVol particle geometry
	 */
	MatrixXd ParticleDiffusionOperatorDG::GSMjacobianParDispBlock(unsigned int parType, double parGeomSurfToVol) {

		MatrixXd dispBlock;

		// We have to match the DGSEM interface, where the dispersion block [ d RHS_disp / d c ] depends on whole previous and subsequent elem plus f_Irst entries of subsubsequent elements
		dispBlock = MatrixXd::Zero(_nParNode[parType], 3 * _nParNode[parType] + 2);

		dispBlock.block(0, _nParNode[parType] + 1, _nParNode[parType], _nParNode[parType]) = _minus_parInvMM_Ar[parType];
		dispBlock *= 2.0 / static_cast<double>(_deltaR[_offsetMetric[parType]]) * 2.0 / static_cast<double>(_deltaR[_offsetMetric[parType]]);

		return -dispBlock; // *-1 for residual
	}

	/**
	 * @brief calculates the dispersion part of the DG jacobian
	 * @param [in] elemIdx elem index
	 * @param [in] parType particle type index
	 * @param [in] parGeomSurfToVol particle geometry
	 */
	MatrixXd ParticleDiffusionOperatorDG::DGjacobianParDispBlock(unsigned int elemIdx, unsigned int parType, double parGeomSurfToVol) {

		MatrixXd dispBlock;
		// Inner dispersion block [ d RHS_disp / d c ], depends on whole previous and subsequent elem plus f_Irst entries of subsubsequent elements
		dispBlock = MatrixXd::Zero(_nParNode[parType], 3 * _nParNode[parType] + 2);

		if (_parExactInt[parType])
		{
			MatrixXd B = getParBMatrix(parType, elemIdx, parGeomSurfToVol); // "Lifting" matrix
			MatrixXd gBlock = getParGBlock(elemIdx, parType); // current elem auxiliary block matrix
			MatrixXd gStarDC = parAuxBlockGstar(elemIdx, parType, getParGBlock(elemIdx - 1, parType), gBlock, getParGBlock(elemIdx + 1, parType)); // Numerical flux block

			if (parGeomSurfToVol != _SurfVolRatioSlab) // weak form DGSEM requ_Ired
				dispBlock.block(0, _nParNode[parType], _nParNode[parType], _nParNode[parType] + 2) = _minus_InvMM_ST[_offsetMetric[parType] + (elemIdx - 1)] * gBlock;
			else // strong form DGSEM
				dispBlock.block(0, _nParNode[parType], _nParNode[parType], _nParNode[parType] + 2) = (_parPolyDerM[parType] - _parInvMM[_offsetMetric[parType] + (elemIdx - 1)] * B) * gBlock;

			dispBlock += _parInvMM[_offsetMetric[parType] + (elemIdx - 1)] * B * gStarDC;
			dispBlock *= 2.0 / static_cast<double>(_deltaR[_offsetMetric[parType] + (elemIdx - 1)]);
		}
		else
		{
			// inexact integration is not maintained due to inferior performance. Code is in calcParticleCollocationDGSEMJacobian
		}

		return -dispBlock; // *-1 for residual
	}

	void ParticleDiffusionOperatorDG::initializeDGjac(std::vector<double> parGeomSurfToVol) {

		// particle jacobian blocks (each is unique)
		_DGjacParDispBlocks = new MatrixXd[std::accumulate(_nParElem, _nParElem + _nParType, 0)];

		for (unsigned int type = 0; type < _nParType; type++)
		{
			for (unsigned int block = 0; block < _nParElem[type]; block++)
			{
				if (_parGSM[type])
					_DGjacParDispBlocks[_offsetMetric[type] + block] = GSMjacobianParDispBlock(type, parGeomSurfToVol[type]);
				else
					_DGjacParDispBlocks[_offsetMetric[type] + block] = DGjacobianParDispBlock(block + 1u, type, parGeomSurfToVol[type]);
			}
		}
	}


}  // namespace parts

}  // namespace model

}  // namespace cadet
