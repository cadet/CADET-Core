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

namespace cadet
{

namespace model
{

namespace parts
{
	ParticleDiffusionOperatorBase::ParticleDiffusionOperatorBase() : _boundOffset(nullptr), _nBound(nullptr), _reqBinding(nullptr), _parDepSurfDiffusion(nullptr)
	{
	}

	ParticleDiffusionOperatorBase::~ParticleDiffusionOperatorBase() CADET_NOEXCEPT
	{
		delete[] _boundOffset;
		delete[] _nBound;
		delete[] _parDepSurfDiffusion;
	}
	

	bool ParticleDiffusionOperatorBase::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp)
	{
		_strideBulkComp = strideBulkComp;

		_parTypeIdx = parTypeIdx;
		_nComp = nComp;

		_filmDiffusion.resize(_nComp); // filled in notifyDiscontinuousSectionTransition
		_poreAccessFactor.resize(_nComp); // filled in notifyDiscontinuousSectionTransition
		_invBetaP.resize(_nComp); // filled in notifyDiscontinuousSectionTransition

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

		// ==== Construct and configure parameter dependencies
		paramProvider.popScope();

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

	bool ParticleDiffusionOperatorBase::configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound, const int* reqBinding, const bool hasDynamicReactions)
	{
		_reqBinding = reqBinding;
		_hasDynamicReactions = hasDynamicReactions;

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
	}

} // namespace parts
} // namespace model
} // namespace cadet