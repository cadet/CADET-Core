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
		_parPolyDeg(nullptr), _nParNode(nullptr), _nParPoints(nullptr), _parGSM(nullptr), //_parTypeOffset(nullptr),
		_nBound(nullptr), _boundOffset(nullptr), _strideBound(nullptr), _nBoundBeforeType(nullptr), _deltaR(nullptr), _parNodes(nullptr),
		_parPolyDerM(nullptr), _minus_InvMM_ST(nullptr), _minus_parInvMM_Ar(nullptr), _parInvWeights(nullptr), _parInvMM(nullptr), _parInvMM_Leg(nullptr),
		_Ir(nullptr), _secondOrderStiffnessM(nullptr), _DGjacParDispBlocks(nullptr), _g_p(nullptr), _g_pSum(nullptr),
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
		delete[] _parGSM;
		//delete[] _parTypeOffset;
		delete[] _nBound;
		delete[] _boundOffset;
		delete[] _strideBound;
		delete[] _nBoundBeforeType;

		delete[] _deltaR;
		delete[] _parNodes;
		delete[] _parPolyDerM;
		delete[] _minus_InvMM_ST;
		delete[] _minus_parInvMM_Ar;
		delete[] _parInvWeights;
		delete[] _parInvMM;
		delete[] _parInvMM_Leg;
		delete[] _Ir;
		delete[] _secondOrderStiffnessM;

		delete[] _DGjacParDispBlocks;

		delete[] _g_p;
		delete[] _g_pSum;
		delete[] _surfaceFluxParticle;
		delete[] _localFlux;

		clearParDepSurfDiffusion();
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
		_invBetaP.resize(_nComp * _nParType); // filled in notifyDiscontinuousSectionTransition
		
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
		std::vector<int> parNelements(_nParType);

		if (firstConfigCall) // avoid memory leaks from reinitialization
		{
			_parPolyDeg = new unsigned int[_nParType];
			_nParElem = new unsigned int[_nParType];
			_parGSM = new bool[_nParType];
		}

		if (paramProvider.exists("PAR_POLYDEG"))
		{
			parPolyDegs = paramProvider.getIntArray("PAR_POLYDEG");

			if ((std::any_of(parPolyDegs.begin(), parPolyDegs.end(), [](int value) { return value < 1; })))
				throw InvalidParameterException("Particle polynomial degrees must be at least 1!");
			parNelements = paramProvider.getIntArray("PAR_NELEM");

			if ((std::any_of(parNelements.begin(), parNelements.end(), [](int value) { return value < 1; })))
				throw InvalidParameterException("Particle number of elements must be at least 1!");

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
			if (parNelements.size() == 1)
			{
				// Multiplex number of particle elements to all particle types
				for (unsigned int i = 0; i < _nParType; ++i)
					std::fill(_nParElem, _nParElem + _nParType, parNelements[0]);
			}
			else if (parNelements.size() < _nParType)
				throw InvalidParameterException("Field PAR_NELEM must have 1 or _nParType (" + std::to_string(_nParType) + ") entries");
			else
				std::copy_n(parNelements.begin(), _nParType, _nParElem);
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
		const bool firstConfigCall = _deltaR == nullptr; // used to not multiply allocate memory

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

		// Check whether all sizes are matched
		if (_nParType != _parRadius.size())
			throw InvalidParameterException("Number of elements in field PAR_RADIUS does not match number of particle types");
		//if (_nParType * nPoints != _parTypeVolFrac.size())
		//	throw InvalidParameterException("Number of elements in field PAR_TYPE_VOLFRAC does not match number of particle types times number of axial elements");
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
		//		throw InvalidParameterException("Sum of field PAR_TYPE_VOLFRAC differs from 1.0 (is " + std::to_string(volFracSum) + ") in axial elem " + std::to_string(i));
		//}

		// Read vectorial parameters (which may also be section dependent; transport)
		//_filmDiffusionMode = readAndRegisterMultiplexCompTypeSecParam(paramProvider, parameters, _filmDiffusion, "FILM_DIFFUSION", _nParType, nComp, unitOpIdx); // todo film diffusion remains unit operation parameter?
		_parDiffusionMode = readAndRegisterMultiplexCompTypeSecParam(paramProvider, parameters, _parDiffusion, "PAR_DIFFUSION", _nParType, _nComp, unitOpIdx);

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
		//	// Register only the first _nParType items
		//	for (unsigned int i = 0; i < _nParType; ++i)
		//		parameters[makeParamId(hashString("PAR_TYPE_VOLFRAC"), unitOpIdx, CompIndep, i, BoundStateIndep, ReactionIndep, SectionIndep)] = &_parTypeVolFrac[i];
		//}
		//else
		//	registerParam2DArray(parameters, _parTypeVolFrac, [=](bool multi, unsigned elem, unsigned int type) { return makeParamId(hashString("PAR_TYPE_VOLFRAC"), unitOpIdx, CompIndep, type, BoundStateIndep, ReactionIndep, elem); }, _nParType);

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
		//		if (_binding[0] && _binding[0]->requiresConfiguration())
		//		{
		//			MultiplexedScopeSelector scopeGuard(paramProvider, "adsorption", true);
		//			bindingConfSuccess = _binding[0]->configure(paramProvider, unitOpIdx, ParTypeIndep);
		//		}
		//	}
		//	else
		//	{
		//		for (unsigned int type = 0; type < _nParType; ++type)
		//		{
		//			if (!_binding[type] || !_binding[type]->requiresConfiguration())
		//				continue;

		//			MultiplexedScopeSelector scopeGuard(paramProvider, "adsorption", type, _nParType == 1, true);
		//			bindingConfSuccess = _binding[type]->configure(paramProvider, unitOpIdx, type) && bindingConfSuccess;
		//		}
		//	}
		//}

		//// Reconfigure reaction model
		//bool dynReactionConfSuccess = true;
		//if (_dynReactionBulk && _dynReactionBulk->requiresConfiguration())
		//{
		//	paramProvider.pushScope("reaction_bulk");
		//	dynReactionConfSuccess = _dynReactionBulk->configure(paramProvider, unitOpIdx, ParTypeIndep);
		//	paramProvider.popScope();
		//}

		//if (_singleDynReaction)
		//{
		//	if (_dynReaction[0] && _dynReaction[0]->requiresConfiguration())
		//	{
		//		MultiplexedScopeSelector scopeGuard(paramProvider, "reaction_particle", true);
		//		dynReactionConfSuccess = _dynReaction[0]->configure(paramProvider, unitOpIdx, ParTypeIndep) && dynReactionConfSuccess;
		//	}
		//}
		//else
		//{
		//	for (unsigned int type = 0; type < _nParType; ++type)
		//	{
		//		if (!_dynReaction[type] || !_dynReaction[type]->requiresConfiguration())
		//			continue;

		//		MultiplexedScopeSelector scopeGuard(paramProvider, "reaction_particle", type, _nParType == 1, true);
		//		dynReactionConfSuccess = _dynReaction[type]->configure(paramProvider, unitOpIdx, type) && dynReactionConfSuccess;
		//	}
		//}

		// jaobian pattern set after binding and particle surface diffusion are configured
		//setJacobianPattern_GRM(globalJac, 0, _dynReactionBulk); // todo set pattern for particles in global unit operation Jacobain

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
			for (unsigned int elem = 0; elem < _nParElem[parType]; ++elem)
			{
				const active r_out = _parRadius[parType] - static_cast<double>(elem) * _Dr;
				const active r_in = _parRadius[parType] - static_cast<double>(elem + 1) * _Dr;

				ptrCenterRadius[elem] = _parRadius[parType] - (0.5 + static_cast<double>(elem)) * _Dr;

				// Compute denominator -> corresponding to elem volume
				const active vol = pow(r_out, 3.0) - pow(r_in, 3.0);

				ptrOuterSurfAreaPerVolume[elem] = 3.0 * sqr(r_out) / vol;
				ptrInnerSurfAreaPerVolume[elem] = 3.0 * sqr(r_in) / vol;
			}
		}
		else if (_parGeomSurfToVol[parType] == _SurfVolRatioCylinder)
		{
			for (unsigned int elem = 0; elem < _nParElem[parType]; ++elem)
			{
				const active r_out = _parRadius[parType] - static_cast<double>(elem) * _Dr;
				const active r_in = _parRadius[parType] - static_cast<double>(elem + 1) * _Dr;

				ptrCenterRadius[elem] = _parRadius[parType] - (0.5 + static_cast<double>(elem)) * _Dr;

				// Compute denominator -> corresponding to elem volume
				const active vol = sqr(r_out) - sqr(r_in);

				ptrOuterSurfAreaPerVolume[elem] = 2.0 * r_out / vol;
				ptrInnerSurfAreaPerVolume[elem] = 2.0 * r_in / vol;
			}
		}
		else if (_parGeomSurfToVol[parType] == _SurfVolRatioSlab)
		{
			for (unsigned int elem = 0; elem < _nParElem[parType]; ++elem)
			{
				const active r_out = _parRadius[parType] - static_cast<double>(elem) * _Dr;
				const active r_in = _parRadius[parType] - static_cast<double>(elem + 1) * _Dr;

				ptrCenterRadius[elem] = _parRadius[parType] - (0.5 + static_cast<double>(elem)) * _Dr;

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
	void ParticleDiffusionOperatorDG::setEquivolumeRadialDisc(unsigned int parType)
	{
		active* const ptrElemSize = _parElementSize.data() + _offsetMetric[parType];
		active* const ptrCenterRadius = _parCenterRadius.data() + _offsetMetric[parType];
		active* const ptrOuterSurfAreaPerVolume = _parOuterSurfAreaPerVolume.data() + _offsetMetric[parType];
		active* const ptrInnerSurfAreaPerVolume = _parInnerSurfAreaPerVolume.data() + _offsetMetric[parType];

		if (_parGeomSurfToVol[parType] == _SurfVolRatioSphere)
		{
			active r_out = _parRadius[parType];
			active r_in = _parCoreRadius[parType];
			const active volumePerelement = (pow(_parRadius[parType], 3.0) - pow(_parCoreRadius[parType], 3.0)) / static_cast<double>(_nParElem[parType]);

			for (unsigned int elem = 0; elem < _nParElem[parType]; ++elem)
			{
				if (elem != (_nParElem[parType] - 1))
					r_in = pow(pow(r_out, 3.0) - volumePerelement, (1.0 / 3.0));
				else
					r_in = _parCoreRadius[parType];

				ptrElemSize[elem] = r_out - r_in;
				ptrCenterRadius[elem] = (r_out + r_in) * 0.5;

				ptrOuterSurfAreaPerVolume[elem] = 3.0 * sqr(r_out) / volumePerelement;
				ptrInnerSurfAreaPerVolume[elem] = 3.0 * sqr(r_in) / volumePerelement;
				// Note that the DG particle elements are oppositely ordered compared to the FV particle elements
				_deltaR[_offsetMetric[parType] + _nParElem[parType] - (elem + 1)] = r_out - r_in;

				// For the next elem: r_out == r_in of the current elem
				r_out = r_in;
			}
		}
		else if (_parGeomSurfToVol[parType] == _SurfVolRatioCylinder)
		{
			active r_out = _parRadius[parType];
			active r_in = _parCoreRadius[parType];
			const active volumePerelement = (sqr(_parRadius[parType]) - sqr(_parCoreRadius[parType])) / static_cast<double>(_nParElem[parType]);

			for (unsigned int elem = 0; elem < _nParElem[parType]; ++elem)
			{
				if (elem != (_nParElem[parType] - 1))
					r_in = sqrt(sqr(r_out) - volumePerelement);
				else
					r_in = _parCoreRadius[parType];

				ptrElemSize[elem] = r_out - r_in;
				ptrCenterRadius[elem] = (r_out + r_in) * 0.5;

				ptrOuterSurfAreaPerVolume[elem] = 2.0 * r_out / volumePerelement;
				ptrInnerSurfAreaPerVolume[elem] = 2.0 * r_in / volumePerelement;
				// Note that the DG particle elements are oppositely ordered compared to the FV particle elements
				_deltaR[_offsetMetric[parType] + _nParElem[parType] - (elem + 1)] = r_out - r_in;

				// For the next elem: r_out == r_in of the current elem
				r_out = r_in;
			}
		}
		else if (_parGeomSurfToVol[parType] == _SurfVolRatioSlab)
		{
			active r_out = _parRadius[parType];
			active r_in = _parCoreRadius[parType];
			const active volumePerelement = (_parRadius[parType] - _parCoreRadius[parType]) / static_cast<double>(_nParElem[parType]);

			for (unsigned int elem = 0; elem < _nParElem[parType]; ++elem)
			{
				if (elem != (_nParElem[parType] - 1))
					r_in = r_out - volumePerelement;
				else
					r_in = _parCoreRadius[parType];

				ptrElemSize[elem] = r_out - r_in;
				ptrCenterRadius[elem] = (r_out + r_in) * 0.5;

				ptrOuterSurfAreaPerVolume[elem] = 1.0 / volumePerelement;
				ptrInnerSurfAreaPerVolume[elem] = 1.0 / volumePerelement;
				// Note that the DG particle elements are oppositely ordered compared to the FV particle elements
				_deltaR[_offsetMetric[parType] + _nParElem[parType] - (elem + 1)] = r_out - r_in;

				// For the next elem: r_out == r_in of the current elem
				r_out = r_in;
			}
		}
	}

	/**
	 * @brief Computes all helper quantities for radial bead discretization from given radial elem boundaries
	 * @details Calculates surface areas per volume for every element and the radial element centers.
	 */
	void ParticleDiffusionOperatorDG::setUserdefinedRadialDisc(unsigned int parType)
	{
		active* const ptrElemSize = _parElementSize.data() + _offsetMetric[parType];
		active* const ptrCenterRadius = _parCenterRadius.data() + _offsetMetric[parType];
		active* const ptrOuterSurfAreaPerVolume = _parOuterSurfAreaPerVolume.data() + _offsetMetric[parType];
		active* const ptrInnerSurfAreaPerVolume = _parInnerSurfAreaPerVolume.data() + _offsetMetric[parType];

		// Care for the right ordering and include 0.0 / 1.0 if not already in the vector.
		std::vector<active> orderedInterfaces = std::vector<active>(_parDiscVector.begin() + _offsetMetric[parType] + parType,
			_parDiscVector.begin() + _offsetMetric[parType] + parType + _nParElem[parType] + 1);

		// Sort in descending order
		std::sort(orderedInterfaces.begin(), orderedInterfaces.end(), std::greater<active>());

		// Force first and last element to be 1.0 and 0.0, respectively
		orderedInterfaces[0] = 1.0;
		orderedInterfaces.back() = 0.0;

		// Map [0, 1] -> [core radius, particle radius] via linear interpolation
		for (unsigned int elem = 0; elem < _nParElem[parType]; ++elem)
			orderedInterfaces[elem] = static_cast<double>(orderedInterfaces[elem]) * (_parRadius[parType] - _parCoreRadius[parType]) + _parCoreRadius[parType];

		if (_parGeomSurfToVol[parType] == _SurfVolRatioSphere)
		{
			for (unsigned int elem = 0; elem < _nParElem[parType]; ++elem)
			{
				ptrElemSize[elem] = orderedInterfaces[elem] - orderedInterfaces[elem + 1];
				ptrCenterRadius[elem] = (orderedInterfaces[elem] + orderedInterfaces[elem + 1]) * 0.5;

				// Compute denominator -> corresponding to elem volume
				const active vol = pow(orderedInterfaces[elem], 3.0) - pow(orderedInterfaces[elem + 1], 3.0);

				ptrOuterSurfAreaPerVolume[elem] = 3.0 * sqr(orderedInterfaces[elem]) / vol;
				ptrInnerSurfAreaPerVolume[elem] = 3.0 * sqr(orderedInterfaces[elem + 1]) / vol;
				// Note that the DG particle elements are oppositely ordered compared to the FV particle elements
				_deltaR[_offsetMetric[parType] + _nParElem[parType] - (elem + 1)] = ptrOuterSurfAreaPerVolume[elem] - ptrInnerSurfAreaPerVolume[elem];
			}
		}
		else if (_parGeomSurfToVol[parType] == _SurfVolRatioCylinder)
		{
			for (unsigned int elem = 0; elem < _nParElem[parType]; ++elem)
			{
				ptrElemSize[elem] = orderedInterfaces[elem] - orderedInterfaces[elem + 1];
				ptrCenterRadius[elem] = (orderedInterfaces[elem] + orderedInterfaces[elem + 1]) * 0.5;

				// Compute denominator -> corresponding to elem volume
				const active vol = sqr(orderedInterfaces[elem]) - sqr(orderedInterfaces[elem + 1]);

				ptrOuterSurfAreaPerVolume[elem] = 2.0 * orderedInterfaces[elem] / vol;
				ptrInnerSurfAreaPerVolume[elem] = 2.0 * orderedInterfaces[elem + 1] / vol;
				// Note that the DG particle elements are oppositely ordered compared to the FV particle elements
				_deltaR[_offsetMetric[parType] + _nParElem[parType] - (elem + 1)] = ptrOuterSurfAreaPerVolume[elem] - ptrInnerSurfAreaPerVolume[elem];
			}
		}
		else if (_parGeomSurfToVol[parType] == _SurfVolRatioSlab)
		{
			for (unsigned int elem = 0; elem < _nParElem[parType]; ++elem)
			{
				ptrElemSize[elem] = orderedInterfaces[elem] - orderedInterfaces[elem + 1];
				ptrCenterRadius[elem] = (orderedInterfaces[elem] + orderedInterfaces[elem + 1]) * 0.5;

				// Compute denominator -> corresponding to elem volume
				const active vol = orderedInterfaces[elem] - orderedInterfaces[elem + 1];

				ptrOuterSurfAreaPerVolume[elem] = 1.0 / vol;
				ptrInnerSurfAreaPerVolume[elem] = 1.0 / vol;
				// Note that the DG particle elements are oppositely ordered compared to the FV particle elements
				_deltaR[_offsetMetric[parType] + _nParElem[parType] - (elem + 1)] = ptrOuterSurfAreaPerVolume[elem] - ptrInnerSurfAreaPerVolume[elem];
			}
		}
	}

	// todo: parameter sensitivities for particle radius. Here, we have the problem that every DG operator becomes an active type.
	// We have sums of three matrices each multiplied by its own metric term, which could be applied iteratively to the residual/solution.
	// Not needed for Slab, and only two matrices required for Cylinder.
	// This approach should only be used when necessary, i.e. solely when particle radius parameter sensitivity is required.
	void ParticleDiffusionOperatorDG::updateRadialDisc()
	{
		for (unsigned int parType = 0; parType < _nParType; ++parType)
		{
			if (_parDiscMode[parType] == ParticleDiscretizationMode::Equidistant)
			{
				for (int elem = 0; elem < _nParElem[parType]; elem++)
				{
					_deltaR[_offsetMetric[parType] + elem] = (_parRadius[parType] - _parCoreRadius[parType]) / _nParElem[parType];
				}
				setEquidistantRadialDisc(parType);
			}
			else if (_parDiscMode[parType] == ParticleDiscretizationMode::Equivolume)
				setEquivolumeRadialDisc(parType);
			else if (_parDiscMode[parType] == ParticleDiscretizationMode::UserDefined)
				setUserdefinedRadialDisc(parType);
		}

		/*		metrics		*/
		// estimate element dependent operators

		for (int parType = 0; parType < _nParType; parType++)
		{
			for (int elem = 0; elem < _nParElem[parType]; elem++)
			{
				for (int node = 0; node < _nParNode[parType]; node++)
					_Ir[_offsetMetric[parType] + elem][node] = _deltaR[_offsetMetric[parType] + elem] / 2.0 * (_parNodes[parType][node] + 1.0);

				active r_L = _parCoreRadius[parType] + elem * _deltaR[_offsetMetric[parType] + elem]; // left boundary of current elem

				_Ir[_offsetMetric[parType] + elem] = _Ir[_offsetMetric[parType] + elem] + VectorXd::Ones(_nParNode[parType]) * r_L;

				if (_parGeomSurfToVol[parType] == _SurfVolRatioSphere)
					_Ir[_offsetMetric[parType] + elem] = _Ir[_offsetMetric[parType] + elem].array().square();
				else if (_parGeomSurfToVol[parType] == _SurfVolRatioSlab)
					_Ir[_offsetMetric[parType] + elem].setOnes(); // no metric terms for slab

				// compute mass matrices for exact integration based on particle geometry, via transformation to normalized Jacobi polynomials with weight function w
				if (_parGeomSurfToVol[parType] == _SurfVolRatioSphere) // r^2 =  r_i^2 + (1 + \xi) * r_i * DeltaR_i / 2.0 + (1 + \xi)^2 * (DeltaR_i / 2.0)^2
				{
					_parInvMM[_offsetMetric[parType] + elem] = parts::dgtoolbox::mMatrix(_parPolyDeg[parType], _parNodes[parType], 0.0, 2.0) * pow((static_cast<double>(_deltaR[_offsetMetric[parType] + elem]) / 2.0), 2.0);
					if (elem > 0 || _parCoreRadius[parType] != 0.0) // following contributions are zero for f_Irst elem when R_c = 0 (no particle core)
						_parInvMM[_offsetMetric[parType] + elem] += parts::dgtoolbox::mMatrix(_parPolyDeg[parType], _parNodes[parType], 0.0, 1.0) * (static_cast<double>(_deltaR[_offsetMetric[parType] + elem]) * static_cast<double>(r_L))
						+ parts::dgtoolbox::mMatrix(_parPolyDeg[parType], _parNodes[parType], 0.0, 0.0) * pow(static_cast<double>(r_L), 2.0);

					_parInvMM[_offsetMetric[parType] + elem] = _parInvMM[_offsetMetric[parType] + elem].inverse();
					_minus_InvMM_ST[_offsetMetric[parType] + elem] = -_parInvMM[_offsetMetric[parType] + elem] * _parPolyDerM[parType].transpose() * _parInvMM[_offsetMetric[parType] + elem].inverse();

					// particle GSM specific second order stiffness matrix (single element, i.e. _nParElem = 1)
					_secondOrderStiffnessM[parType] = std::pow(static_cast<double>(_parCoreRadius[parType]), 2.0) * parts::dgtoolbox::secondOrderStiffnessMatrix(_parPolyDeg[parType], 0.0, 0.0, _parNodes[parType]);
					_secondOrderStiffnessM[parType] += static_cast<double>(_deltaR[_offsetMetric[parType]]) * static_cast<double>(_parCoreRadius[parType]) * parts::dgtoolbox::secondOrderStiffnessMatrix(_parPolyDeg[parType], 0.0, 1.0, _parNodes[parType]);
					_secondOrderStiffnessM[parType] += std::pow(static_cast<double>(_deltaR[_offsetMetric[parType]]) / 2.0, 2.0) * parts::dgtoolbox::secondOrderStiffnessMatrix(_parPolyDeg[parType], 0.0, 2.0, _parNodes[parType]);
				}
				else if (_parGeomSurfToVol[parType] == _SurfVolRatioCylinder) // r = r_i + (1 + \xi) * DeltaR_i / 2.0
				{
					_parInvMM[_offsetMetric[parType] + elem] = parts::dgtoolbox::mMatrix(_parPolyDeg[parType], _parNodes[parType], 0.0, 1.0) * (static_cast<double>(_deltaR[_offsetMetric[parType] + elem]) / 2.0);
					if (elem > 0 || _parCoreRadius[parType] != 0.0) // following contribution is zero for f_Irst elem when R_c = 0 (no particle core)
						_parInvMM[_offsetMetric[parType] + elem] += parts::dgtoolbox::mMatrix(_parPolyDeg[parType], _parNodes[parType], 0.0, 0.0) * static_cast<double>(r_L);

					_parInvMM[_offsetMetric[parType] + elem] = _parInvMM[_offsetMetric[parType] + elem].inverse();
					_minus_InvMM_ST[_offsetMetric[parType] + elem] = -_parInvMM[_offsetMetric[parType] + elem] * _parPolyDerM[parType].transpose() * _parInvMM[_offsetMetric[parType] + elem].inverse();

					// particle GSM specific second order stiffness matrix (single element, i.e. _nParElem = 1)
					_secondOrderStiffnessM[parType] = static_cast<double>(_parCoreRadius[parType]) * parts::dgtoolbox::secondOrderStiffnessMatrix(_parPolyDeg[parType], 0.0, 0.0, _parNodes[parType]);
					_secondOrderStiffnessM[parType] += static_cast<double>(_deltaR[_offsetMetric[parType]]) / 2.0 * parts::dgtoolbox::secondOrderStiffnessMatrix(_parPolyDeg[parType], 0.0, 1.0, _parNodes[parType]);
				}
				else if (_parGeomSurfToVol[parType] == _SurfVolRatioSlab) // r = 1
				{
					_minus_InvMM_ST[_offsetMetric[parType] + elem] = -_parInvMM[_offsetMetric[parType] + elem] * _parPolyDerM[parType].transpose() * _parInvMM[_offsetMetric[parType] + elem].inverse();

					_secondOrderStiffnessM[parType] = parts::dgtoolbox::secondOrderStiffnessMatrix(_parPolyDeg[parType], 0.0, 0.0, _parNodes[parType]);
					_parInvMM[_offsetMetric[parType] + elem] = parts::dgtoolbox::invMMatrix(_parPolyDeg[parType], _parNodes[parType], 0.0, 0.0);
				}
			}

			_minus_parInvMM_Ar[parType] = -_parInvMM[_offsetMetric[parType]] * _secondOrderStiffnessM[parType];
		}
	}

	/**
	 * @brief allocates memory for DG operators and computes those that are metric independent. Also allocates required containers needed for the DG discretization.
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
			_Ir = new Vector<active, Dynamic>[_offsetMetric[_nParType]];
			for (int parType = 0; parType < _nParType; parType++)
			{
				for (int elem = 0; elem < _nParElem[parType]; elem++)
				{
					_Ir[_offsetMetric[parType] + elem].resize(_nParNode[parType]);
					_Ir[_offsetMetric[parType] + elem].setZero();
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
		for (int dep = 0; dep < _nComp * _nParType; dep++)
		{
			_filmDiffusion[dep] = filmDiff[dep];
			_poreAccessFactor[dep] = poreAccessFactor[dep];
		}
		for (int parType = 0; parType < _nParType; parType++)
		{
			for (int comp = 0; comp < _nComp; comp++)
			{
				_invBetaP[parType * _nComp + comp] = (1.0 - _parPorosity[parType]) / (_poreAccessFactor[_nComp * parType + comp] * _parPorosity[parType]);
			}
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
	int ParticleDiffusionOperatorDG::getParticleCoordinates(unsigned int parType, double* coords) const
	{
		active const* const pcr = _parCenterRadius.data() + _offsetMetric[parType];

		// Note that the DG particle shells are oppositely ordered compared to the FV particle shells
		for (unsigned int par = 0; par < _nParPoints[parType]; par++) {

			unsigned int cell = std::floor(par / _nParNode[parType]);

			double r_L = static_cast<double>(pcr[cell]) - 0.5 * static_cast<double>(_deltaR[_offsetMetric[parType] + cell]);
			coords[par] = r_L + 0.5 * static_cast<double>(_deltaR[_offsetMetric[parType] + cell]) * (1.0 + _parNodes[parType][par % _nParNode[parType]]);
		}

		return _nParPoints[parType];
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
	int ParticleDiffusionOperatorDG::residual(double t, unsigned int secIdx, unsigned int parType, unsigned int colNode, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, int const* const qsBinding, WithoutParamSensitivity)
	{
		return residualImpl<double, double, double>(t, secIdx, parType, colNode, yPar, yBulk, yDotPar, resPar, qsBinding);
	}

	int ParticleDiffusionOperatorDG::residual(double t, unsigned int secIdx, unsigned int parType, unsigned int colNode, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, int const* const qsBinding, WithoutParamSensitivity)
	{
		return residualImpl<active, active, double>(t, secIdx, parType, colNode, yPar, yBulk, yDotPar, resPar, qsBinding);
	}

	int ParticleDiffusionOperatorDG::residual(double t, unsigned int secIdx, unsigned int parType, unsigned int colNode, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, int const* const qsBinding, WithParamSensitivity)
	{
		return residualImpl<double, active, active>(t, secIdx, parType, colNode, yPar, yBulk, yDotPar, resPar, qsBinding);
	}

	int ParticleDiffusionOperatorDG::residual(double t, unsigned int secIdx, unsigned int parType, unsigned int colNode, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, int const* const qsBinding, WithParamSensitivity)
	{
		return residualImpl<active, active, active>(t, secIdx, parType, colNode, yPar, yBulk, yDotPar, resPar, qsBinding);
	}

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

	template <typename StateType, typename ResidualType, typename ParamType>
	int ParticleDiffusionOperatorDG::residualImpl(double t, unsigned int secIdx, unsigned int parType, unsigned int colNode, StateType const* yPar, StateType const* yBulk, double const* yDotPar, ResidualType* resPar, int const* const qsBinding)
	{
		/* Mobile phase RHS	*/

		// Get film diffusion flux at current node to compute boundary condition
		for (unsigned int comp = 0; comp < _nComp; comp++) {
			_localFlux[comp] = _filmDiffusion[comp + parType * _nComp] * (yBulk[comp * _strideBulkComp] - yPar[(_nParPoints[parType] - 1) * strideParNode(parType) + comp]);
		}

		active const* const parDiff = getSectionDependentSlice(_parDiffusion, _nComp * _nParType, secIdx) + parType * _nComp;

		// Ordering of particle surface diffusion:
		// bnd0comp0, bnd1comp0, bnd0comp1, bnd1comp1, bnd0comp2, bnd1comp2
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
					if (parSurfDiff[_boundOffset[parType * _nComp + comp] + bnd] != 0.0) // some bound states might still not be effected by surface diffusion
					{
						Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> c_s(yPar + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd, nPoints, InnerStride<Dynamic>(strideParNode(parType)));
						sum_cp_cs += static_cast<ParamType>(_invBetaP[_nComp * parType + comp]) * static_cast<ParamType>(parSurfDiff[_boundOffset[parType * _nComp + comp] + bnd]) * c_s;

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
								static_cast<ParamType>(parSurfDiff[_boundOffset[parType * _nComp + comp] + bnd]) * c_s;

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
				// =====================================================================================================//
				// Solve auxiliary systems  d_p g_p + d_s beta_p sum g_s= d (d_p c_p + d_s beta_p sum c_s) / d xi		//
				// =====================================================================================================//
				
				// Component-wise! strides
				unsigned int strideelem = nNodes;
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
						if (parSurfDiff[_boundOffset[parType * _nComp + comp] + bnd] != 0.0) // some bound states might still not be effected by surface diffusion
						{
							// Get solid phase vector
							Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>> c_s(yPar + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd,
								nPoints, InnerStride<Dynamic>(strideParNode(parType)));
							// Compute g_s = d c_s / d xi
							solve_auxiliary_DG<StateType>(parType, c_s, strideelem, strideNode, comp);
							// Apply invBeta_p, d_s and add to sum -> gSum += d_s * invBeta_p * (D c - M^-1 B [c - c^*])
							g_pSum += g_p.template cast<ResidualType>() * static_cast<ParamType>(_invBetaP[_nComp * parType + comp]) * static_cast<ParamType>(parSurfDiff[_boundOffset[parType * _nComp + comp] + bnd]);

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
								g_p_ResType *= static_cast<ParamType>(parSurfDiff[_boundOffset[parType * _nComp + comp] + bnd]);

								// Eigen access to auxiliary variable of current bound state
								Eigen::Map<const Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>> g_p_ResType_const(&g_p_ResType[0], nPoints, InnerStride<Dynamic>(1));

								// Add - D_r * gs to the residual, including metric part.->res = invMap^2* [ -D_r * (d_s c^s) ]
								parVolumeIntegral<ResidualType, ResidualType>(parType, false, g_p_ResType_const, resCs);

								// Add M^-1 B (gs - gs^*) to the residual -> res =  invMap^2 * [ - D_r * (d_s c^s) + M^-1 B (gs - gs^*) ]
								parSurfaceIntegral<ResidualType, ResidualType>(parType, g_p_ResType_const, resCs, strideelem, strideNode, false, comp, true);
							}
						}
					}
				}

				// Compute g_p = d c_p / d xi
				solve_auxiliary_DG<StateType>(parType, cp, strideelem, strideNode, comp);

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
				parSurfaceIntegral<ResidualType, ResidualType>(parType, g_pSum_const, resCp, strideelem, strideNode, false, comp);

			}
		}

		return true;
	}

	template<typename ResidualType, typename ParamType>
	void ParticleDiffusionOperatorDG::applyParInvMap(Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<>>& state, unsigned int parType)
	{
		for (int elem = 0; elem < _nParElem[parType]; elem++) {
			state.segment(elem * _nParNode[parType], _nParNode[parType]) *= 2.0 / static_cast<ParamType>(_deltaR[_offsetMetric[parType] + elem]) * 2.0 / static_cast<ParamType>(_deltaR[_offsetMetric[parType] + elem]);
		}
	}

	template<typename StateType, typename ResidualType>
	void ParticleDiffusionOperatorDG::parGSMVolumeIntegral(const int parType, Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state, Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>>& stateDer)
	{
		int nNodes = _nParNode[parType];

		stateDer.segment(0, nNodes)
			-= (_minus_parInvMM_Ar[parType].template cast<StateType>() * state.segment(0, nNodes)).template cast<ResidualType>();
	}

	template<typename StateType, typename ResidualType>
	void ParticleDiffusionOperatorDG::parVolumeIntegral(const int parType, const bool aux, Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state, Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>>& stateDer)
	{
		int nNodes = _nParNode[parType];

		/* no additional metric term for auxiliary equation or particle equation with exact integration scheme
		   -> res = - D * (d_p * c^p + invBeta_p sum_mi d_s c^s) */
		if (aux || (_parGeomSurfToVol[parType] == _SurfVolRatioSlab)) {
			// comp-elem-node state vector: use of Eigen lib performance
			for (unsigned int elem = 0; elem < _nParElem[parType]; elem++) {
				stateDer.segment(elem * nNodes, nNodes)
					-= (_parPolyDerM[parType].template cast<StateType>() * state.segment(elem * nNodes, nNodes)).template cast<ResidualType>();
			}
		}
		else if (_parGeomSurfToVol[parType] != _SurfVolRatioSlab) {
			// comp-elem-node state vector: use of Eigen lib performance
			for (unsigned int elem = 0; elem < _nParElem[parType]; elem++) {
				stateDer.segment(elem * nNodes, nNodes)
					-= (_minus_InvMM_ST[_offsetMetric[parType] + elem].template cast<StateType>() * state.segment(elem * nNodes, nNodes)).template cast<ResidualType>();
			}
		}
	}
	/*
	 * @brief calculates the interface fluxes g* of particle mass balance equation and implements the respective boundary conditions
	 * @param [in] aux bool if interface flux for auxiliary equation
	 * @param [in] addParDisc bool if interface flux for additional particle DG-discretized equation
	*/
	template<typename StateType>
	void ParticleDiffusionOperatorDG::InterfaceFluxParticle(int parType, Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state,
		const unsigned int strideelem, const unsigned int strideNode, const bool aux, const int comp, const bool addParDisc)
	{
		Eigen::Map<Vector<StateType, Dynamic>, 0, InnerStride<>> _surfFluxPar(reinterpret_cast<StateType*>(&_surfaceFluxParticle[parType][0]), _nParElem[parType] + 1, InnerStride<>(1));

		// reset surface flux storage as it is used multiple times
		_surfFluxPar.setZero();

		// numerical flux: state* = 0.5 (state^+ + state^-)

		// calculate inner interface fluxes
		for (unsigned int elem = 1u; elem < _nParElem[parType]; elem++) {
			_surfFluxPar[elem] // left interfaces
				= 0.5 * (state[elem * strideelem - strideNode] + // outer/left node
					state[elem * strideelem]); // inner/right node
		}

		// calculate boundary interface fluxes.
		if (aux) { // ghost nodes given by state^- := state^+ for auxiliary equation
			_surfFluxPar[0] = state[0];

			_surfFluxPar[_nParElem[parType]] = state[_nParElem[parType] * strideelem - strideNode];
		}
		else if (addParDisc) {
			_surfFluxPar[0] = 0.0;

			_surfFluxPar[_nParElem[parType]] = 0.0;
		}
		else {

			// film diffusion BC
			_surfFluxPar[_nParElem[parType]] = static_cast<StateType>(_localFlux[comp])
				/ (static_cast<double>(_parPorosity[parType]) * static_cast<double>(_poreAccessFactor[parType * _nComp + comp]))
				* (2.0 / static_cast<double>(_deltaR[_offsetMetric[parType]])); // inverse squared mapping was also applied, so we apply Map * invMap^2 = invMap

			// inner particle BC
			_surfFluxPar[0] = 0.0;

		}
	}
	/**
	 * @brief calculates the particle surface Integral (type- and component-wise)
	 * @param [in] parType current particle type
	 * @param [in] state relevant state vector
	 * @param [in] stateDer state derivative vector the solution is added to
	 * @param [in] aux true for auxiliary equation, false for main equation
	 * @param [in] strideelem component-wise elem stride
	 * @param [in] strideNodecomponent-wise node stride
	 * @param [in] comp current component
	*/
	template<typename StateType, typename ResidualType>
	void ParticleDiffusionOperatorDG::parSurfaceIntegral(int parType, Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<Dynamic>>& state,
		Eigen::Map<Vector<ResidualType, Dynamic>, 0, InnerStride<Dynamic>>& stateDer, unsigned const int strideelem, unsigned const int strideNode,
		const bool aux, const int comp, const bool addParDisc)
	{
		Eigen::Map<Vector<StateType, Dynamic>, 0, InnerStride<>> _surfFluxPar(reinterpret_cast<StateType*>(&_surfaceFluxParticle[parType][0]), _nParElem[parType] + 1, InnerStride<>(1));

		// calc numerical flux values
		InterfaceFluxParticle<StateType>(parType, state, strideelem, strideNode, aux, comp, addParDisc);

		// strong surface integral -> M^-1 B [state - state*]
		for (unsigned int elem = 0; elem < _nParElem[parType]; elem++) {

			for (unsigned int Node = 0; Node < _nParNode[parType]; Node++) {
				if (aux) { // strong surface integral -> M^-1 B [state - state*] and B has two non-zero entries, -1 and 1
					stateDer[elem * strideelem + Node * strideNode]
						-= _parInvMM_Leg[parType](Node, 0) * (state[elem * strideelem] - _surfFluxPar[elem])
						- _parInvMM_Leg[parType](Node, _parPolyDeg[parType]) * (state[elem * strideelem + _parPolyDeg[parType] * strideNode] - _surfFluxPar[elem + 1u]);
				}
				else {
					if (_parGeomSurfToVol[parType] == _SurfVolRatioSlab) { // strong surface integral -> M^-1 B [state - state*] and B has two non-zero entries, -1 and 1
						stateDer[elem * strideelem + Node * strideNode]
							-= static_cast<ResidualType>(
								_parInvMM[parType](Node, 0) * (state[elem * strideelem] - _surfFluxPar[elem])
								- _parInvMM[parType](Node, _parPolyDeg[parType]) * (state[elem * strideelem + _parPolyDeg[parType] * strideNode] - _surfFluxPar[elem + 1u])
								);
					}
					else if (_parGeomSurfToVol[parType] == _SurfVolRatioCylinder) { // weak surface integral -> M^-1 B [- state*] and B has two non-zero entries, which depend on metrics
						stateDer[elem * strideelem + Node * strideNode]
							-= static_cast<ResidualType>(
								_Ir[_offsetMetric[parType] + elem][0] * _parInvMM[_offsetMetric[parType] + elem](Node, 0) * (-_surfFluxPar[elem])
								+ _Ir[_offsetMetric[parType] + elem][_nParNode[parType] - 1] * _parInvMM[_offsetMetric[parType] + elem](Node, _parPolyDeg[parType]) * _surfFluxPar[elem + 1u]
								);
					}
					else if (_parGeomSurfToVol[parType] == _SurfVolRatioSphere) { // weak surface integral -> M^-1 B [- state*] and B has two non-zero entries, which depend on metrics
						stateDer[elem * strideelem + Node * strideNode]
							-= static_cast<ResidualType>(
								_Ir[_offsetMetric[parType] + elem][0] * _parInvMM[_offsetMetric[parType] + elem](Node, 0) * (-_surfFluxPar[elem])
								+ _Ir[_offsetMetric[parType] + elem][_nParNode[parType] - 1] * _parInvMM[_offsetMetric[parType] + elem](Node, _parPolyDeg[parType]) * _surfFluxPar[elem + 1u]
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
	void ParticleDiffusionOperatorDG::solve_auxiliary_DG(int parType, Eigen::Map<const Vector<StateType, Dynamic>, 0, InnerStride<>>& conc, unsigned int strideelem, unsigned int strideNode, int comp)
	{
		Eigen::Map<Vector<StateType, Dynamic>, 0, InnerStride<>> g_p(reinterpret_cast<StateType*>(&_g_p[parType][0]), _nParPoints[parType], InnerStride<>(1));
		Eigen::Map<Vector<StateType, Dynamic>, 0, InnerStride<>> _surfFluxPar(reinterpret_cast<StateType*>(&_surfaceFluxParticle[parType][0]), _nParElem[parType] + 1, InnerStride<>(1));
		_surfFluxPar.setZero(); // reset surface flux storage as it is used multiple times
		g_p.setZero(); // reset auxiliary variable g

		// ========================================================================================//
		// solve auxiliary systems g = d c / d xi	 =>		g_p = Dc - M^-1 B [c - c^*]			   //
		// ========================================================================================//

		parVolumeIntegral<StateType, StateType>(parType, true, conc, g_p); // volumne integral in strong DG form: - D c

		parSurfaceIntegral<StateType>(parType, conc, g_p, strideelem, strideNode, true, comp); // surface integral in strong DG form: M^-1 B [c - c^*]

		g_p *= -1.0; // auxiliary factor -1
	}

	// ==========================================================================================================================================================  //
	// ========================================						DG particle Jacobian							=============================================  //
	// ==========================================================================================================================================================  //

	/**
	 * @brief calculates the particle dispersion jacobian Pattern of the DG scheme for the given particle type and bead
	*/
	void ParticleDiffusionOperatorDG::calcParticleJacobianPattern(std::vector<ParticleDiffusionOperatorDG::T>& tripletList, unsigned int offset, unsigned int parType, unsigned int colNode, unsigned int secIdx, int const* const qsBinding)
	{
		// Ordering of particle surface diffusion:
		// bnd0comp0, bnd1comp0, bnd0comp1, bnd1comp1, bnd0comp2, bnd1comp2
		active const* const _parSurfDiff = getSectionDependentSlice(_parSurfDiffusion, _strideBound[_nParType], secIdx) + _nBoundBeforeType[parType];

		// (global) strides
		unsigned int selem = _nParNode[parType] * strideParNode(parType);
		unsigned int sNode = strideParNode(parType);
		unsigned int sComp = 1u;
		unsigned int nNodes = _nParNode[parType];

		// case: one elem  -> diffBlock \in R^(nParNodes x nParNodes), GBlock = parPolyDerM
		if (_nParElem[parType] == 1) {

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
						if (_hasSurfaceDiffusion[parType]) {

							for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
								if (_parSurfDiff[_boundOffset[parType * _nComp + comp] + bnd] != 0.0) {
									// row: add current component offset and go node strides from there for each dispersion block entry
									// col: jump oover liquid states, add current bound state offset and go node strides from there for each dispersion block entry
									tripletList.push_back(T(offset + comp * sComp + i * sNode,
										offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode, 0.0));

									/* add surface diffusion dispersion block to solid */
									if (!qsBinding[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
										// row: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry
										// col: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry
										tripletList.push_back(T(offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
											offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode, 0.0));

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

			unsigned int special = 0u; if (_nParElem[parType] < 3u) special = 1u; // limits the iterator for special case nelements = 3 (dependence on additional entry)
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
						if (_hasSurfaceDiffusion[parType]) {

							for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
								if (_parSurfDiff[_boundOffset[parType * _nComp + comp] + bnd] != 0.0) {
									// row: add current component offset and go node strides from there for each dispersion block entry
									// col: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry. adjust for j start
									tripletList.push_back(T(offset + comp * sComp + i * sNode,
										offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode - (nNodes + 1) * sNode,
										0.0));

									/* add surface diffusion dispersion block to solid */
									if (!qsBinding[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
										// row: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry
										// col: jump over liquid states, add current bound state offset and go node strides from there for each dispersion block entry. adjust for j start
										tripletList.push_back(T(offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
											offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode - (nNodes + 1) * sNode,
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
						tripletList.push_back(T(offset + comp * sComp + (_nParElem[parType] - 1) * selem + i * sNode,
							offset + comp * sComp + (_nParElem[parType] - 1) * selem - selem - sNode + j * sNode,
							0.0));

						// handle surface diffusion of bound states. binding is handled in residualKernel().
						if (_hasSurfaceDiffusion[parType]) {

							for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
								if (_parSurfDiff[_boundOffset[parType * _nComp + comp] + bnd] != 0.0) {
									// row: add component offset and jump over previous elements. Go node strides from there for each dispersion block entry
									// col: jump over liquid states, add current bound state offset and jump over previous elements. Go back one elem (and node or adjust for start) and go node strides from there for each dispersion block entry.
									tripletList.push_back(T(offset + comp * sComp + (_nParElem[parType] - 1) * selem + i * sNode,
										offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + (_nParElem[parType] - 2) * selem - sNode + j * sNode,
										0.0));

									/* add surface diffusion dispersion block to solid */
									if (!qsBinding[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
										// row: jump over previous elements and over liquid states, add current bound state offset. go node strides from there for each dispersion block entry
										// col: jump over previous elements and over liquid states, add current bound state offset. go back one elem (and node or adjust for start) and go node strides from there for each dispersion block entry.
										tripletList.push_back(T(offset + (_nParElem[parType] - 1) * selem + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
											offset + (_nParElem[parType] - 2) * selem + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) - sNode + bnd + j * sNode,
											0.0));
									}
								}
							}
						}
					}
				}
			}
			if (_nParElem[parType] == 3) {
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
							if (_hasSurfaceDiffusion[parType]) {

								for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
									if (_parSurfDiff[_boundOffset[parType * _nComp + comp] + bnd] != 0.0) {
										// row: add component offset and jump over previous elements. Go node strides from there for each dispersion block entry
										// col: jump over liquid states, add current bound state offset and jump over previous elem. go back one elem and go node strides from there for each dispersion block entry. adjust for j start
										tripletList.push_back(T(offset + comp * sComp + selem + i * sNode,
											offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode - sNode,
											0.0));

										/* add surface diffusion dispersion block to solid */
										if (!qsBinding[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
											// row: jump over previous elements and over liquid states, add current bound state offset. go node strides from there for each dispersion block entry
											// col: jump over liquid states, add current bound state offset and jump over previous elem. go node strides from there for each dispersion block entry. adjust for j start
											tripletList.push_back(T(offset + selem + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
												offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode - sNode,
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
			if (_nParElem[parType] >= 4) {

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
							if (_hasSurfaceDiffusion[parType]) {

								for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
									if (_parSurfDiff[_boundOffset[parType * _nComp + comp] + bnd] != 0.0) {
										// row: add component offset and jump over previous elements. Go node strides from there for each dispersion block entry
										// col: jump over liquid states, add current bound state offset and jump over previous elements. Go back one elem and go node strides from there for each dispersion block entry. adjust for j start
										tripletList.push_back(T(offset + comp * sComp + selem + i * sNode,
											offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode - sNode,
											0.0));

										/* add surface diffusion dispersion block to solid */
										if (!qsBinding[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
											// row: jump over previous elements and over liquid states, add current bound state offset. go node strides from there for each dispersion block entry
											// col: jump over previous elements and over liquid states, add current bound state offset. go back one elem and go node strides from there for each dispersion block entry. adjust for j start
											tripletList.push_back(T(offset + selem + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
												offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + j * sNode - sNode,
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
							tripletList.push_back(T(offset + comp * sComp + (_nParElem[parType] - 2) * selem + i * sNode,
								offset + comp * sComp + (_nParElem[parType] - 2) * selem - selem - sNode + j * sNode,
								0.0));

							// handle surface diffusion of bound states. binding is handled in residualKernel().
							if (_hasSurfaceDiffusion[parType]) {

								for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
									if (_parSurfDiff[_boundOffset[parType * _nComp + comp] + bnd] != 0.0) {
										// row: add component offset and jump over previous elements. Go node strides from there for each dispersion block entry
										// col: jump over liquid states, add current bound state offset and jump over previous elements. Go back one elem and node and go node strides from there for each dispersion block entry
										tripletList.push_back(T(offset + comp * sComp + (_nParElem[parType] - 2) * selem + i * sNode,
											offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + (_nParElem[parType] - 2) * selem - selem - sNode + j * sNode,
											0.0));

										/* add surface diffusion dispersion block to solid */
										if (!qsBinding[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
											// row: jump over previous elements and over liquid states, add current bound state offset. go node strides from there for each dispersion block entry
											// col: jump over previous elements and over liquid states, add current bound state offset. go back one elem and node and go node strides from there for each dispersion block entry
											tripletList.push_back(T(offset + (_nParElem[parType] - 2) * selem + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
												offset + (_nParElem[parType] - 2) * selem + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) - selem - sNode + bnd + j * sNode,
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

			if (_nParElem[parType] >= 5) {

				for (unsigned int elem = 2; elem < _nParElem[parType] - 2; elem++) {

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
								if (_hasSurfaceDiffusion[parType]) {

									for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
										if (_parSurfDiff[_boundOffset[parType * _nComp + comp] + bnd] != 0.0) {
											// row: add component offset and jump over previous elements. Go node strides from there for each dispersion block entry
											// col: jump over liquid states, add current bound state offset and jump over previous elements. Go back one elem and node and go node strides from there for each dispersion block entry
											tripletList.push_back(T(offset + comp * sComp + elem * selem + i * sNode,
												offset + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + elem * selem - selem - sNode + j * sNode,
												0.0));

											/* add surface diffusion dispersion block to solid */
											if (!qsBinding[offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd]) {
												// row: jump over previous elements and over liquid states, add current bound state offset. go node strides from there for each dispersion block entry
												// col: jump over previous elements and over liquid states, add current bound state offset. go back one elem and node and go node strides from there for each dispersion block entry
												tripletList.push_back(T(offset + elem * selem + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd + i * sNode,
													offset + elem * selem + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd - selem - sNode + j * sNode,
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
	
	unsigned int ParticleDiffusionOperatorDG::calcParDispNNZ(int parType)
	{
		return _nComp * ((3u * _nParElem[parType] - 2u) * _nParNode[parType] * _nParNode[parType] + (2u * _nParElem[parType] - 3u) * _nParNode[parType]);
	}

	/**
	 *@brief adds the time derivative entries from particle equations
	 *@detail since the main diagonal entries are already set, we actually only set the solid phase time derivative entries for the discretized particle mass balance equations
	 */
	void ParticleDiffusionOperatorDG::parTimeDerJacPattern_GRM(std::vector<ParticleDiffusionOperatorDG::T>& tripletList, unsigned int offset, unsigned int parType, unsigned int colNode, unsigned int secIdx)
	{
		active const* const _parSurfDiff = getSectionDependentSlice(_parSurfDiffusion, _strideBound[_nParType], secIdx) + _nBoundBeforeType[parType];

		for (unsigned int parNode = 0; parNode < _nParPoints[parType]; parNode++)
		{
			for (unsigned int comp = 0; comp < _nComp; comp++)
			{
				for (unsigned int bnd = 0; bnd < _nBound[parType * _nComp + comp]; bnd++) {
					// row: jump over previous nodes add current component offset
					// col: jump over previous nodes, liquid phase and previous bound states
					tripletList.push_back(T(offset + parNode * strideParNode(parType) + comp,
						offset + parNode * strideParNode(parType) + strideParLiquid() + offsetBoundComp(ParticleTypeIndex{ parType }, ComponentIndex{ comp }) + bnd,
						0.0));
				}
			}
		}
	}
	/**
	 * @brief sets the sparsity pattern of the binding Jacobian
	 */
	void ParticleDiffusionOperatorDG::parBindingPattern_GRM(std::vector<ParticleDiffusionOperatorDG::T>& tripletList, const int offset, const unsigned int parType, const unsigned int colNode)
	{
		// every bound state might depend on every bound and liquid state
		for (int parNode = 0; parNode < _nParPoints[parType]; parNode++)
		{
			for (int bnd = 0; bnd < _strideBound[parType]; bnd++)
			{
				for (int conc = 0; conc < strideParNode(parType); conc++) {
					// row: jump over previous nodes and liquid states and add current bound state offset
					// col: jump over previous nodes and add current concentration offset (liquid and bound)
					tripletList.push_back(T(offset + parNode * strideParNode(parType) + strideParLiquid() + bnd,
						offset + parNode * strideParNode(parType) + conc, 0.0));
				}
			}
		}
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

		if (elemIdx == 0 || elemIdx == _nParElem[parType] + 1)
		{ // elemIdx out of bounds
			return MatrixXd::Zero(_nParNode[parType], _nParNode[parType] + 2);
		}
		if (elemIdx != 1 && elemIdx != _nParElem[parType])
		{ // inner elem
			gBlock.block(0, 0, _nParNode[parType], 1) -= 0.5 * _parInvMM_Leg[parType].block(0, 0, _nParNode[parType], 1);
			gBlock.block(0, 1, _nParNode[parType], 1) += 0.5 * _parInvMM_Leg[parType].block(0, 0, _nParNode[parType], 1);
			gBlock.block(0, _nParNode[parType], _nParNode[parType], 1) -= 0.5 * _parInvMM_Leg[parType].block(0, _nParNode[parType] - 1, _nParNode[parType], 1);
			gBlock.block(0, _nParNode[parType] + 1, _nParNode[parType], 1) += 0.5 * _parInvMM_Leg[parType].block(0, _nParNode[parType] - 1, _nParNode[parType], 1);
		}
		else if (elemIdx == 1u)
		{ // left boundary elem
			if (elemIdx == _nParElem[parType]) // special case one elem
				return gBlock * 2.0 / static_cast<double>(_deltaR[_offsetMetric[parType] + (elemIdx - 1)]);
			gBlock.block(0, _nParNode[parType], _nParNode[parType], 1) -= 0.5 * _parInvMM_Leg[parType].block(0, _nParNode[parType] - 1, _nParNode[parType], 1);
			gBlock.block(0, _nParNode[parType] + 1, _nParNode[parType], 1) += 0.5 * _parInvMM_Leg[parType].block(0, _nParNode[parType] - 1, _nParNode[parType], 1);
		}
		else if (elemIdx == _nParElem[parType])
		{ // right boundary elem
			gBlock.block(0, 0, _nParNode[parType], 1) -= 0.5 * _parInvMM_Leg[parType].block(0, 0, _nParNode[parType], 1);
			gBlock.block(0, 1, _nParNode[parType], 1) += 0.5 * _parInvMM_Leg[parType].block(0, 0, _nParNode[parType], 1);
		}
		gBlock *= 2.0 / static_cast<double>(_deltaR[_offsetMetric[parType] + (elemIdx - 1)]);

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
		if (elemIdx != 1)
		{
			gStarDC.block(0, _nParNode[parType], 1, _nParNode[parType] + 2) += middleG.block(0, 0, 1, _nParNode[parType] + 2);
			gStarDC.block(0, 0, 1, _nParNode[parType] + 2) += leftG.block(_nParNode[parType] - 1, 0, 1, _nParNode[parType] + 2);
		}
		if (elemIdx != _nParElem[parType])
		{
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
		if (parGeomSurfToVol == _SurfVolRatioSlab)
		{
			B(0, 0) = -1.0;
			B(_nParNode[parType] - 1, _nParNode[parType] - 1) = 1.0;
		}
		else
		{
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
	MatrixXd ParticleDiffusionOperatorDG::GSMjacobianParDispBlock(unsigned int parType, double parGeomSurfToVol)
	{
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
	MatrixXd ParticleDiffusionOperatorDG::DGjacobianParDispBlock(unsigned int elemIdx, unsigned int parType, double parGeomSurfToVol)
	{
		MatrixXd dispBlock;
		// Inner dispersion block [ d RHS_disp / d c ], depends on whole previous and subsequent elem plus f_Irst entries of subsubsequent elements
		dispBlock = MatrixXd::Zero(_nParNode[parType], 3 * _nParNode[parType] + 2);

		MatrixXd B = getParBMatrix(parType, elemIdx, parGeomSurfToVol); // "Lifting" matrix
		MatrixXd gBlock = getParGBlock(elemIdx, parType); // current elem auxiliary block matrix
		MatrixXd gStarDC = parAuxBlockGstar(elemIdx, parType, getParGBlock(elemIdx - 1, parType), gBlock, getParGBlock(elemIdx + 1, parType)); // Numerical flux block

		if (parGeomSurfToVol != _SurfVolRatioSlab) // weak form DGSEM required
			dispBlock.block(0, _nParNode[parType], _nParNode[parType], _nParNode[parType] + 2) = _minus_InvMM_ST[_offsetMetric[parType] + (elemIdx - 1)] * gBlock;
		else // strong form DGSEM
			dispBlock.block(0, _nParNode[parType], _nParNode[parType], _nParNode[parType] + 2) = (_parPolyDerM[parType] - _parInvMM[_offsetMetric[parType] + (elemIdx - 1)] * B) * gBlock;

		dispBlock += _parInvMM[_offsetMetric[parType] + (elemIdx - 1)] * B * gStarDC;
		dispBlock *= 2.0 / static_cast<double>(_deltaR[_offsetMetric[parType] + (elemIdx - 1)]);

		return -dispBlock; // *-1 for residual
	}

	void ParticleDiffusionOperatorDG::initializeDGjac(std::vector<double> parGeomSurfToVol)
	{
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
	/**
	 * @brief adds jacobian entries which have been overwritten by the binding kernel (only use for surface diffusion combined with kinetic binding)
	 * @detail only adds the entries d RHS_i / d c^s_i, which lie on the diagonal
	 */
	int ParticleDiffusionOperatorDG::addSolidDGentries(const int secIdx, const int parType, const int nBulk, const int* const reqBinding, const int offsetCp, Eigen::SparseMatrix<double, RowMajor>& globalJac)
	{
		active const* const parSurfDiff = getSectionDependentSlice(_parSurfDiffusion, _strideBound[_nParType], secIdx) + _nBoundBeforeType[parType];

		for (unsigned int blk = 0; blk < nBulk; blk++) {
			// Get jacobian iterator at first solid entry of first particle of current type
			linalg::BandedEigenSparseRowIterator jac(globalJac, offsetCp + blk * strideParBlock(parType) + strideParLiquid());

			for (unsigned int elem = 0; elem < _nParElem[parType]; elem++)
				addDiagonalSolidJacobianEntries(_DGjacParDispBlocks[_offsetMetric[parType] + elem].block(0, _nParNode[parType] + 1, _nParNode[parType], _nParNode[parType]), jac, reqBinding, parSurfDiff, parType);
		}

		return 1;
	}
	/**
	 * @brief adds a state block into the system jacobian.
	 * @param [in] block (sub)block whose diagonal entries are to be added
	 * @param [in] jac row iterator at first (i.e. upper) entry
	 * @param [in] reqBinding pointer to binding kinetics
	 * @param [in] type particle type
	 */
	void ParticleDiffusionOperatorDG::addDiagonalSolidJacobianEntries(Eigen::MatrixXd block, linalg::BandedEigenSparseRowIterator& jac, const int* const reqBinding, const active* const surfDiffPtr, unsigned int type)
	{
		for (unsigned int i = 0; i < block.rows(); i++, jac += strideParLiquid())
		{
			for (unsigned int comp = 0; comp < _nComp; comp++)
			{
				for (unsigned int bnd = 0; bnd < _nBound[type * _nComp + comp]; bnd++, ++jac)
				{

					if (static_cast<double>(surfDiffPtr[offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd]) != 0.0
						&& !reqBinding[offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd]) {
						// row, col: at current node and bound state
						jac[0] += block(i, i)
							* static_cast<double>(surfDiffPtr[offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd]);
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
	void ParticleDiffusionOperatorDG::insertParJacBlock(Eigen::MatrixXd block, linalg::BandedEigenSparseRowIterator& jac, const active* const parDiff, const active* const surfDiff, const active* const beta_p, const int* nonKinetic, unsigned int type, unsigned int nBlocks, int offRowToCol)
	{
		for (unsigned int elem = 0; elem < nBlocks; elem++) {
			for (unsigned int i = 0; i < block.rows(); i++) {
				for (unsigned int comp = 0; comp < _nComp; comp++, ++jac) {
					for (unsigned int j = 0; j < block.cols(); j++) {
						/* liquid on liquid blocks */
						// row: at current node and component; col: jump to node j
						jac[(j - i) * strideParNode(type) + offRowToCol] = block(i, j) * static_cast<double>(parDiff[comp]);
					}
					/* liquid on solid blocks */
					for (unsigned int bnd = 0; bnd < _nBound[type * _nComp + comp]; bnd++) {
						if (static_cast<double>(surfDiff[offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd]) != 0.0) {
							for (unsigned int j = 0; j < block.cols(); j++) {
								// row: at current node and component; col: jump to node j and to current bound state
								jac[(j - i) * strideParNode(type) + offRowToCol + strideParLiquid() - comp
									+ offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd
								]
									= block(i, j) * static_cast<double>(beta_p[comp])
									* static_cast<double>(surfDiff[offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd]);
							}
						}
					}
				}
				/* solid on solid blocks */
				for (unsigned int comp = 0; comp < _nComp; comp++) {
					for (unsigned int bnd = 0; bnd < _nBound[type * _nComp + comp]; bnd++, ++jac) {
						if (static_cast<double>(surfDiff[offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd]) != 0.0
							&& !nonKinetic[offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd]) {
							for (unsigned int j = 0; j < block.cols(); j++) {
								// row: at current node and bound state; col: jump to node j
								jac[(j - i) * strideParNode(type) + offRowToCol + bnd]
									= block(i, j)
									* static_cast<double>(surfDiff[offsetBoundComp(ParticleTypeIndex{ type }, ComponentIndex{ comp }) + bnd]);
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
	int ParticleDiffusionOperatorDG::calcStaticAnaParticleDiffJacobian(const int secIdx, const int parType, const int colNode, const int* const reqBinding, const int offsetLocalCp, Eigen::SparseMatrix<double, RowMajor>& globalJac)
	{
		// Prepare parameters
		const active* const parDiff = getSectionDependentSlice(_parDiffusion, _nComp * _nParType, secIdx) + parType * _nComp;

		// Ordering of particle surface diffusion:
		// bnd0comp0, bnd1comp0, bnd0comp1, bnd1comp1, bnd0comp2, bnd1comp2
		const active* const  parSurfDiff = getSectionDependentSlice(_parSurfDiffusion, _strideBound[_nParType], secIdx) + _nBoundBeforeType[parType];

		const active* const invBetaP = &_invBetaP[parType * _nComp];

		// (global) strides
		unsigned int selem = _nParNode[parType] * strideParNode(parType);
		unsigned int sNode = strideParNode(parType);
		unsigned int sComp = 1u;
		unsigned int nNodes = _nParNode[parType];

		/* Special case */
		if (_nParElem[parType] == 1) {

			linalg::BandedEigenSparseRowIterator jacIt(globalJac, offsetLocalCp); // row iterator starting at first elem, first component

			insertParJacBlock(_DGjacParDispBlocks[_offsetMetric[parType] + 0].block(0, nNodes + 1, nNodes, nNodes), jacIt, parDiff, parSurfDiff, invBetaP, reqBinding, parType, 1u, 0);
			return 1;
		}

		/* Special case */
		if (_nParElem[parType] == 2) {

			linalg::BandedEigenSparseRowIterator jacIt(globalJac, offsetLocalCp); // row iterator starting at first elem, first component

			insertParJacBlock(_DGjacParDispBlocks[_offsetMetric[parType] + 0].block(0, nNodes + 1, nNodes, 2 * nNodes), jacIt, parDiff, parSurfDiff, invBetaP, reqBinding, parType, 1u, 0);
			// right Bacobian block, iterator is already moved to second elem
			insertParJacBlock(_DGjacParDispBlocks[_offsetMetric[parType] + 1].block(0, 1, nNodes, 2 * nNodes), jacIt, parDiff, parSurfDiff, invBetaP, reqBinding, parType, 1u, -strideParElem(parType));
			return 1;
		}

		/* Special case */
		if (_nParElem[parType] == 3) {

			linalg::BandedEigenSparseRowIterator jacIt(globalJac, offsetLocalCp + strideParElem(parType)); // row iterator starting at first elem, first component

			insertParJacBlock(_DGjacParDispBlocks[_offsetMetric[parType] + 1].block(0, 1, nNodes, 3 * nNodes), jacIt, parDiff, parSurfDiff, invBetaP, reqBinding, parType, 1u, -strideParElem(parType));
		}

		/* Inner elements (exist only if nelements >= 5) */
		if (_nParElem[parType] >= 5) {

			linalg::BandedEigenSparseRowIterator jacIt(globalJac, offsetLocalCp + strideParElem(parType) * 2); // row iterator starting at third elem, first component

			// insert all (nElem - 4) inner elem blocks
			for (unsigned int elem = 2; elem < _nParElem[parType] - 2; elem++)
				insertParJacBlock(_DGjacParDispBlocks[_offsetMetric[parType] + elem], jacIt, parDiff, parSurfDiff, invBetaP, reqBinding, parType, 1u, -(strideParElem(parType) + strideParNode(parType)));
		}

		/*	boundary elem neighbours (exist only if nelements >= 4)	*/
		if (_nParElem[parType] >= 4) {

			linalg::BandedEigenSparseRowIterator jacIt(globalJac, offsetLocalCp + strideParElem(parType)); // row iterator starting at second elem, first component

			insertParJacBlock(_DGjacParDispBlocks[_offsetMetric[parType] + 1].block(0, 1, nNodes, 3 * nNodes + 1), jacIt, parDiff, parSurfDiff, invBetaP, reqBinding, parType, 1u, -strideParElem(parType));

			jacIt += (_nParElem[parType] - 4) * strideParElem(parType); // move iterator to preultimate elem (already at third elem)
			insertParJacBlock(_DGjacParDispBlocks[_offsetMetric[parType] + _nParElem[parType] - 2u].block(0, 0, nNodes, 3 * nNodes + 1), jacIt, parDiff, parSurfDiff, invBetaP, reqBinding, parType, 1u, -(strideParElem(parType) + strideParNode(parType)));
		}

		/*			boundary elements (exist only if nelements >= 3)			*/
		if (_nParElem[parType] >= 3) {

			linalg::BandedEigenSparseRowIterator jacIt(globalJac, offsetLocalCp); // row iterator starting at first elem, first component

			insertParJacBlock(_DGjacParDispBlocks[_offsetMetric[parType] + 0].block(0, nNodes + 1, nNodes, 2 * nNodes + 1), jacIt, parDiff, parSurfDiff, invBetaP, reqBinding, parType, 1u, 0);

			jacIt += (_nParElem[parType] - 2) * strideParElem(parType); // move iterator to last elem (already at second elem)
			insertParJacBlock(_DGjacParDispBlocks[_offsetMetric[parType] + _nParElem[parType] - 1u].block(0, 0, nNodes, 2 * nNodes + 1), jacIt, parDiff, parSurfDiff, invBetaP, reqBinding, parType, 1u, -(strideParElem(parType) + strideParNode(parType)));
		}

		return true;
	}

	int ParticleDiffusionOperatorDG::calcFilmDiffJacobian(unsigned int secIdx, const int parType, const int offsetCp, const int offsetC, const int nBulkPoints, const double colPorosity, const active* const parTypeVolFrac, Eigen::SparseMatrix<double, RowMajor>& globalJac, bool outliersOnly)
	{
		// lifting matrix entry for exact integration scheme depends on metrics for sphere and cylinder
		double exIntLiftContribution = static_cast<double>(_Ir[_offsetMetric[parType] + _nParElem[parType] - 1][_nParNode[parType] - 1]);
		if (_parGeomSurfToVol[parType] == _SurfVolRatioSlab)
			exIntLiftContribution = 1.0;

		// Ordering of diffusion:
		// sec0type0comp0, sec0type0comp1, sec0type0comp2, sec0type1comp0, sec0type1comp1, sec0type1comp2,
		// sec1type0comp0, sec1type0comp1, sec1type0comp2, sec1type1comp0, sec1type1comp1, sec1type1comp2, ...
		active const* const filmDiff = getSectionDependentSlice(_filmDiffusion, _nComp * _nParType, secIdx) + parType * _nComp;

		linalg::BandedEigenSparseRowIterator jacCl(globalJac, offsetC);
		linalg::BandedEigenSparseRowIterator jacCp(globalJac, offsetCp + (_nParPoints[parType] - 1) * strideParNode(parType)); // iterator at the outer particle boundary

		for (unsigned int blk = 0; blk < nBulkPoints; blk++)
		{
			for (unsigned int comp = 0; comp < _nComp; comp++, ++jacCp, ++jacCl) {
				// add Cl on Cl entries (added since these entries are also touched by bulk jacobian)
				// row: already at bulk phase. already at current node and component.
				// col: already at bulk phase. already at current node and component.
				if (!outliersOnly)
					jacCl[0] += static_cast<double>(filmDiff[comp]) * (1.0 - colPorosity) / colPorosity
					* _parGeomSurfToVol[parType] / static_cast<double>(_parRadius[parType])
					* static_cast<double>(parTypeVolFrac[parType + blk * _nParType]);
				// add Cl on Cp entries (added since these entries are also touched by bulk jacobian)
				// row: already at bulk phase. already at current node and component.
				// col: go to current particle phase entry.
				jacCl[jacCp.row() - jacCl.row()] = -static_cast<double>(filmDiff[comp]) * (1.0 - colPorosity) / colPorosity
					* _parGeomSurfToVol[parType] / static_cast<double>(_parRadius[parType])
					* static_cast<double>(parTypeVolFrac[parType + blk * _nParType]);


				unsigned int entry = jacCp.row();
				for (int node = _parPolyDeg[parType]; node >= 0; node--, jacCp -= strideParNode(parType)) {
					// row: already at particle. Already at current node and liquid state.
					// col: original entry at outer node.
					if (!outliersOnly) // Cp on Cb
						jacCp[entry - jacCp.row()]
						+= static_cast<double>(filmDiff[comp]) * 2.0 / static_cast<double>(_deltaR[_offsetMetric[parType]]) * _parInvMM[_offsetMetric[parType] + _nParElem[parType] - 1](node, _nParNode[parType] - 1) * exIntLiftContribution / static_cast<double>(_parPorosity[parType]) / static_cast<double>(_poreAccessFactor[parType * _nComp + comp]);
					// row: already at particle. Already at current node and liquid state.
					// col: go to current bulk phase.
					jacCp[jacCl.row() - jacCp.row()]
						= -static_cast<double>(filmDiff[comp]) * 2.0 / static_cast<double>(_deltaR[_offsetMetric[parType]]) * _parInvMM[_offsetMetric[parType] + _nParElem[parType] - 1](node, _nParNode[parType] - 1) * exIntLiftContribution / static_cast<double>(_parPorosity[parType]) / static_cast<double>(_poreAccessFactor[parType * _nComp + comp]);
				}
				// set back iterator to first node as required by component loop
				jacCp += _nParNode[parType] * strideParNode(parType);
			}
			if (blk < nBulkPoints - 1) // execute iteration statement only when condition is true in next loop.
				jacCp += _strideBound[parType] + (_nParPoints[parType] - 1) * strideParNode(parType);
		}

		return 1;
	}

	bool ParticleDiffusionOperatorDG::setParameter(const ParameterId& pId, double value)
	{
		if (multiplexCompTypeSecParameterValue(pId, hashString("PAR_DIFFUSION"), _parDiffusionMode, _parDiffusion, _nParType, _nComp, value, nullptr))
			return true;
		if (multiplexBndCompTypeSecParameterValue(pId, hashString("PAR_SURFDIFFUSION"), _parSurfDiffusionMode, _parSurfDiffusion, _nParType, _nComp, _strideBound, _nBound, _boundOffset, _nBoundBeforeType, value, nullptr))
			return true;

		if (multiplexTypeParameterValue(pId, hashString("PAR_RADIUS"), _singleParRadius, _parRadius, value, nullptr))
			return true;
		if (multiplexTypeParameterValue(pId, hashString("PAR_CORERADIUS"), _singleParCoreRadius, _parCoreRadius, value, nullptr))
			return true;
		if (multiplexTypeParameterValue(pId, hashString("PAR_POROSITY"), _singleParPorosity, _parPorosity, value, nullptr))
			return true;

		if (model::setParameter(pId, value, _parDepSurfDiffusion, _singleParDepSurfDiffusion))
			return true;

		return false;
	}

	bool ParticleDiffusionOperatorDG::setParameter(const ParameterId& pId, int value)
	{
		if (model::setParameter(pId, value, _parDepSurfDiffusion, _singleParDepSurfDiffusion))
			return true;
	}

	bool ParticleDiffusionOperatorDG::setParameter(const ParameterId& pId, bool value)
	{
		if (model::setParameter(pId, value, _parDepSurfDiffusion, _singleParDepSurfDiffusion))
			return true;
	}

	bool ParticleDiffusionOperatorDG::setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& pId, double value)
	{
		if (multiplexCompTypeSecParameterValue(pId, hashString("PAR_DIFFUSION"), _parDiffusionMode, _parDiffusion, _nParType, _nComp, value, &sensParams))
			return true;
		if (multiplexBndCompTypeSecParameterValue(pId, hashString("PAR_SURFDIFFUSION"), _parSurfDiffusionMode, _parSurfDiffusion, _nParType, _nComp, _strideBound, _nBound, _boundOffset, _nBoundBeforeType, value, &sensParams))
			return true;

		if (multiplexTypeParameterValue(pId, hashString("PAR_RADIUS"), _singleParRadius, _parRadius, value, &sensParams))
			return true;
		if (multiplexTypeParameterValue(pId, hashString("PAR_CORERADIUS"), _singleParCoreRadius, _parCoreRadius, value, &sensParams))
			return true;
		if (multiplexTypeParameterValue(pId, hashString("PAR_POROSITY"), _singleParPorosity, _parPorosity, value, &sensParams))
			return true;

		if (model::setSensitiveParameterValue(pId, value, sensParams, _parDepSurfDiffusion, _singleParDepSurfDiffusion))
			return true;

		return false;
	}

	bool ParticleDiffusionOperatorDG::setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue)
	{
		if (multiplexCompTypeSecParameterAD(pId, hashString("PAR_DIFFUSION"), _parDiffusionMode, _parDiffusion, _nParType, _nComp, adDirection, adValue, sensParams))
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		if (multiplexBndCompTypeSecParameterAD(pId, hashString("PAR_SURFDIFFUSION"), _parSurfDiffusionMode, _parSurfDiffusion, _nParType, _nComp, _strideBound, _nBound, _boundOffset, _nBoundBeforeType, adDirection, adValue, sensParams))
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		if (model::setSensitiveParameter(pId, adDirection, adValue, sensParams, _parDepSurfDiffusion, _singleParDepSurfDiffusion))
		{
			LOG(Debug) << "Found parameter " << pId << " in surface diffusion parameter dependence: Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		if (multiplexTypeParameterAD(pId, hashString("PAR_RADIUS"), _singleParRadius, _parRadius, adDirection, adValue, sensParams))
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		if (multiplexTypeParameterAD(pId, hashString("PAR_CORERADIUS"), _singleParCoreRadius, _parCoreRadius, adDirection, adValue, sensParams))
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		if (multiplexTypeParameterAD(pId, hashString("PAR_POROSITY"), _singleParPorosity, _parPorosity, adDirection, adValue, sensParams))
		{
			LOG(Debug) << "Found parameter " << pId << ": Dir " << adDirection << " is set to " << adValue;
			return true;
		}

		return false;
	}


}  // namespace parts

}  // namespace model

}  // namespace cadet
