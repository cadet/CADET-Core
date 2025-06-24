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

#include "model/ColumnModel2D.hpp"
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
#include "SimulationTypes.hpp"
#include "linalg/Norms.hpp"

#include "AdUtils.hpp"
#include "SensParamUtil.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <algorithm>
#include <functional>
#include <numeric>

#include "ParallelSupport.hpp"
#ifdef CADET_PARALLELIZE
#include <tbb/parallel_for.h>
#endif

namespace
{

cadet::model::MultiplexMode readAndRegisterMultiplexParam(cadet::IParameterProvider& paramProvider, std::unordered_map<cadet::ParameterId, cadet::active*>& parameters, std::vector<cadet::active>& values, const std::string& name, unsigned int nAxial, unsigned int nRad, unsigned int nParType, cadet::UnitOpIdx uoi)
{
	cadet::model::MultiplexMode mode = cadet::model::MultiplexMode::Independent;
	readScalarParameterOrArray(values, paramProvider, name, 1);
	if (paramProvider.exists(name + "_MULTIPLEX"))
	{
		const int modeConfig = paramProvider.getInt(name + "_MULTIPLEX");
		if (modeConfig == 0)
		{
			mode = cadet::model::MultiplexMode::Independent;
			if (values.size() != nParType)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(nParType) + ")");
		}
		else if (modeConfig == 1)
		{
			mode = cadet::model::MultiplexMode::Radial;
			if (values.size() != nRad * nParType)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(nRad * nParType) + ")");
		}
		else if (modeConfig == 2)
		{
			mode = cadet::model::MultiplexMode::Axial;
			if (values.size() != nAxial * nParType)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(nAxial * nParType) + ")");
		}
		else if (modeConfig == 3)
		{
			mode = cadet::model::MultiplexMode::AxialRadial;
			if (values.size() != nAxial * nRad * nParType)
				throw cadet::InvalidParameterException("Number of elements in field " + name + " inconsistent with " + name + "_MULTIPLEX (should be " + std::to_string(nAxial * nRad * nParType) + ")");
		}
	}
	else
	{
		if (values.size() == nParType)
			mode = cadet::model::MultiplexMode::Independent;
		else if (values.size() == nRad * nParType)
			mode = cadet::model::MultiplexMode::Radial;
		else if (values.size() == nAxial * nParType)
			mode = cadet::model::MultiplexMode::Axial;
		else if (values.size() == nRad * nAxial * nParType)
			mode = cadet::model::MultiplexMode::AxialRadial;
		else
			throw cadet::InvalidParameterException("Could not infer multiplex mode of field " + name + ", set " + name + "_MULTIPLEX or change number of elements");
	}

	const cadet::StringHash nameHash = cadet::hashStringRuntime(name);
	switch (mode)
	{
		case cadet::model::MultiplexMode::Independent:
			{
				std::vector<cadet::active> p(nAxial * nRad * nParType);
				for (unsigned int s = 0; s < nAxial * nRad; ++s)
					std::copy(values.begin(), values.end(), p.begin() + s * nParType);

				values = std::move(p);

				for (unsigned int s = 0; s < nParType; ++s)
					parameters[cadet::makeParamId(nameHash, uoi, cadet::CompIndep, s, cadet::BoundStateIndep, cadet::ReactionIndep, cadet::SectionIndep)] = &values[s];
			}
			break;
		case cadet::model::MultiplexMode::Radial:
			{
				std::vector<cadet::active> p(nAxial * nRad * nParType);
				for (unsigned int s = 0; s < nAxial; ++s)
					std::copy(values.begin(), values.end(), p.begin() + s * nParType * nRad);

				values = std::move(p);

				for (unsigned int s = 0; s < nRad; ++s)
				{
					for (unsigned int i = 0; i < nParType; ++i)
						parameters[cadet::makeParamId(nameHash, uoi, cadet::CompIndep, i, cadet::BoundStateIndep, s, cadet::SectionIndep)] = &values[s * nParType + i];
				}
			}
			break;
		case cadet::model::MultiplexMode::Axial:
			{
				std::vector<cadet::active> p(nAxial * nRad * nParType);
				for (unsigned int i = 0; i < nAxial; ++i)
				{
					for (unsigned int j = 0; j < nRad; ++j)
						std::copy(values.begin() + i * nParType, values.begin() + (i+1) * nParType, p.begin() + i * nRad * nParType + j * nParType);
				}

				values = std::move(p);

				for (unsigned int s = 0; s < nAxial; ++s)
				{
					for (unsigned int i = 0; i < nParType; ++i)
						parameters[cadet::makeParamId(nameHash, uoi, cadet::CompIndep, i, cadet::BoundStateIndep, cadet::ReactionIndep, s)] = &values[s * nParType * nRad + i];
				}
			}
			break;
		case cadet::model::MultiplexMode::AxialRadial:
			cadet::registerParam3DArray(parameters, values, [=](bool multi, unsigned int ax, unsigned int rad, unsigned int pt) { return cadet::makeParamId(nameHash, uoi, cadet::CompIndep, pt, cadet::BoundStateIndep, rad, ax); }, nParType, nRad);
			break;
		case cadet::model::MultiplexMode::RadialSection:
		case cadet::model::MultiplexMode::Component:
		case cadet::model::MultiplexMode::ComponentRadial:
		case cadet::model::MultiplexMode::ComponentRadialSection:
		case cadet::model::MultiplexMode::ComponentSection:
		case cadet::model::MultiplexMode::Section:
		case cadet::model::MultiplexMode::Type:
		case cadet::model::MultiplexMode::ComponentType:
		case cadet::model::MultiplexMode::ComponentSectionType:
			cadet_assert(false);
			break;
	}

	return mode;
}

bool multiplexParameterValue(const cadet::ParameterId& pId, cadet::StringHash nameHash, cadet::model::MultiplexMode mode, std::vector<cadet::active>& data, unsigned int nAxial, unsigned int nRad, unsigned int nParType, double value, std::unordered_set<cadet::active*> const* sensParams)
{
	if (pId.name != nameHash)
		return false;

	switch (mode)
	{
		case cadet::model::MultiplexMode::Independent:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.particleType]))
					return false;

				for (unsigned int i = 0; i < nAxial * nRad; ++i)
					data[i * nParType + pId.particleType].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::Radial:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction == cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.reaction * nParType + pId.particleType]))
					return false;

				for (unsigned int i = 0; i < nAxial; ++i)
					data[i * nRad * nParType + pId.reaction * nParType + pId.particleType].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::Axial:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.section * nParType * nRad + pId.particleType]))
					return false;

				for (unsigned int i = 0; i < nRad; ++i)
					data[pId.section * nParType * nRad + i * nParType + pId.particleType].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::AxialRadial:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction == cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				if (sensParams && !cadet::contains(*sensParams, &data[pId.section * nParType * nRad + pId.reaction * nParType + pId.particleType]))
					return false;

				data[pId.section * nParType * nRad + pId.reaction * nParType + pId.particleType].setValue(value);

				return true;
			}
		case cadet::model::MultiplexMode::RadialSection:
		case cadet::model::MultiplexMode::Component:
		case cadet::model::MultiplexMode::ComponentRadial:
		case cadet::model::MultiplexMode::ComponentRadialSection:
		case cadet::model::MultiplexMode::ComponentSection:
		case cadet::model::MultiplexMode::Section:
		case cadet::model::MultiplexMode::Type:
		case cadet::model::MultiplexMode::ComponentType:
		case cadet::model::MultiplexMode::ComponentSectionType:
			cadet_assert(false);
			break;
	}

	return false;
}

bool multiplexParameterAD(const cadet::ParameterId& pId, cadet::StringHash nameHash, cadet::model::MultiplexMode mode, std::vector<cadet::active>& data, unsigned int nAxial, unsigned int nRad, unsigned int nParType, unsigned int adDirection, double adValue, std::unordered_set<cadet::active*>& sensParams)
{
	if (pId.name != nameHash)
		return false;

	switch (mode)
	{
		case cadet::model::MultiplexMode::Independent:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.particleType]);

				for (unsigned int i = 0; i < nAxial * nRad; ++i)
					data[i * nParType + pId.particleType].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::Radial:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction == cadet::ReactionIndep) || (pId.section != cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.reaction * nParType + pId.particleType]);

				for (unsigned int i = 0; i < nAxial; ++i)
					data[i * nRad * nParType + pId.reaction * nParType + pId.particleType].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::Axial:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction != cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.section * nParType * nRad + pId.particleType]);

				for (unsigned int i = 0; i < nRad; ++i)
					data[pId.section * nParType * nRad + i * nParType + pId.particleType].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::AxialRadial:
			{
				if ((pId.component != cadet::CompIndep) || (pId.particleType == cadet::ParTypeIndep) || (pId.boundState != cadet::BoundStateIndep)
					|| (pId.reaction == cadet::ReactionIndep) || (pId.section == cadet::SectionIndep))
					return false;

				sensParams.insert(&data[pId.section * nParType * nRad + pId.reaction * nParType + pId.particleType]);
				data[pId.section * nParType * nRad + pId.reaction * nParType + pId.particleType].setADValue(adDirection, adValue);

				return true;
			}
		case cadet::model::MultiplexMode::RadialSection:
		case cadet::model::MultiplexMode::Component:
		case cadet::model::MultiplexMode::ComponentRadial:
		case cadet::model::MultiplexMode::ComponentRadialSection:
		case cadet::model::MultiplexMode::ComponentSection:
		case cadet::model::MultiplexMode::Section:
		case cadet::model::MultiplexMode::Type:
		case cadet::model::MultiplexMode::ComponentType:
		case cadet::model::MultiplexMode::ComponentSectionType:
			cadet_assert(false);
			break;
	}

	return false;
}


}  // namespace


namespace cadet
{

namespace model
{

ColumnModel2D::ColumnModel2D(UnitOpIdx unitOpIdx) : UnitOperationBase(unitOpIdx),
	_dynReactionBulk(nullptr), _jacInlet(),	_analyticJac(true), _jacobianAdDirs(0), _factorizeJacobian(false), _tempState(nullptr),
	_initC(0), _singleRadiusInitC(true), _initCp(0), _singleRadiusInitCp(true), _initCs(0), _singleRadiusInitCs(true), _initState(0), _initStateDot(0)
{
}

ColumnModel2D::~ColumnModel2D() CADET_NOEXCEPT
{
	delete[] _tempState;

	for (IParticleModel* pm : _particles)
		delete pm;

	_particles.clear();

	_binding.clear(); // binding models are deleted in the respective particle model
	_dynReaction.clear(); // particle reaction models are deleted in the respective particle model

	delete _dynReactionBulk;

	delete[] _disc.parTypeOffset;
	delete[] _disc.nBound;
	delete[] _disc.boundOffset;
	delete[] _disc.strideBound;
	delete[] _disc.nBoundBeforeType;

	delete _linearSolver;
}

unsigned int ColumnModel2D::numDofs() const CADET_NOEXCEPT
{
	// Column bulk DOFs: axNPoints * radNPoints * nComp
	// Particle DOFs: axNPoints * radNPoints * nParType particles each having nComp (liquid phase) + sum boundStates (solid phase) DOFs
	// Inlet DOFs: nComp * nRad
	return _disc.axNPoints * _disc.radNPoints * _disc.nComp + _disc.parTypeOffset[_disc.nParType] + _disc.nComp * _disc.radNPoints;
}

unsigned int ColumnModel2D::numPureDofs() const CADET_NOEXCEPT
{
	// Column bulk DOFs: axNPoints * radNPoints * nComp
	// Particle DOFs: axNPoints * radNPoints * nParType particles each having nComp (liquid phase) + sum boundStates (solid phase) DOFs
	return _disc.axNPoints * _disc.radNPoints * _disc.nComp + _disc.parTypeOffset[_disc.nParType];
}


bool ColumnModel2D::usesAD() const CADET_NOEXCEPT
{
#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	// We always need AD if we want to check the analytical Jacobian
	return true;
#else
	// We only need AD if we are not computing the Jacobian analytically
	return !_analyticJac;
#endif
}

bool ColumnModel2D::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper)
{
	const bool firstConfigCall = _tempState == nullptr; // used to avoid multiply allocation

	_disc.nComp = paramProvider.getInt("NCOMP");

	_disc.nParType = paramProvider.exists("NPARTYPE") ? paramProvider.getInt("NPARTYPE") : 0;
	if (_disc.nParType < 0)
		throw InvalidParameterException("Number of particle types must be >= 0!");

	if (_disc.nParType == 0 && paramProvider.exists("particle_type_000"))
		throw InvalidParameterException("NPARTYPE is set to 0, but group particle_type_000 exists.");

	// Create and configure particle model
	_particles = std::vector<IParticleModel*>(_disc.nParType, nullptr);

	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		paramProvider.pushScope("particle_type_" + std::string(3 - std::to_string(parType).length(), '0') + std::to_string(parType));
		std::string particleModelName = paramProvider.getString("PARTICLE_TYPE");

		_particles[parType] = helper.createParticleModel(particleModelName);
		if (!_particles[parType])
			throw InvalidParameterException("Unknown particle model " + particleModelName);

		paramProvider.popScope();
	}

	bool particleConfSuccess = true;
	Indexer idxr(_disc);

	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		particleConfSuccess = particleConfSuccess && _particles[parType]->configureModelDiscretization(paramProvider, helper, _disc.nComp, parType, _disc.nParType, idxr.strideColComp());
	}

	if (firstConfigCall)
	{
		_disc.nParPoints = new unsigned int[_disc.nParType];
		for (int type = 0; type < _disc.nParType; type++)
		{
			_disc.nParPoints[type] = _particles[type]->nDiscPoints();
		}
	}

	_disc.newStaticJac = true;

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

	// ==== Construct and configure convection dispersion operator

	const unsigned int strideRadNode = _disc.nComp;
	const bool transportSuccess = _convDispOp.configureModelDiscretization(paramProvider, helper, _disc.nComp, strideRadNode);

	_disc.axNPoints = _convDispOp.axNPoints();
	_disc.radNPoints = _convDispOp.radNPoints();
	_disc.radNNodes = _convDispOp.radNNodes();
	_disc.radNElem = _convDispOp.radNElem();
	_disc.nBulkPoints = _disc.axNPoints * _disc.radNPoints;

	// Precompute offsets of particle type DOFs
	if (firstConfigCall)
		_disc.parTypeOffset = new unsigned int[_disc.nParType + 1];

	_disc.parTypeOffset[0] = 0;

	for (unsigned int j = 1; j < _disc.nParType + 1; ++j)
		_disc.parTypeOffset[j] = _disc.parTypeOffset[j - 1] + (_disc.nComp + _disc.strideBound[j - 1]) * _disc.nParPoints[j - 1] * _disc.nBulkPoints;

	paramProvider.pushScope("discretization");

	if (firstConfigCall)
		_linearSolver = cadet::linalg::setLinearSolver(paramProvider.exists("LINEAR_SOLVER") ? paramProvider.getString("LINEAR_SOLVER") : "SparseLU");

	// Determine whether analytic Jacobian should be used but don't set it right now.
	// We need to setup Jacobian matrices first.
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	const bool analyticJac = paramProvider.getBool("USE_ANALYTIC_JACOBIAN");
#else
	const bool analyticJac = false;
#endif

	// Create nonlinear solver for consistent initialization
	configureNonlinearSolver(paramProvider);

	paramProvider.popScope();

	// Allocate memory
	_initC.resize(_disc.nComp * _disc.radNElem);
	_initCp.resize(_disc.nComp * _disc.nParType * _disc.radNElem);
	_initCs.resize(nTotalBound * _disc.radNElem);

	// Set whether analytic Jacobian is used
	useAnalyticJacobian(analyticJac);

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
		_binding[parType] = _particles[parType]->getBinding();

		_dynReaction[parType] = _particles[parType]->getReaction();

		// Check if binding and reaction particle type dependence is the same for all particle types
		if (parType > 0)
		{
			if (_binding[parType])
			{
				if (_singleBinding != !_particles[parType]->bindingParDep())
					throw InvalidParameterException("Binding particle type dependence must be the same for all particle types, check field BINDING_PARTYPE_DEPENDENT");
			}

			if (_dynReaction[parType])
			{
				if (_singleDynReaction != !_particles[parType]->reactionParDep())
					throw InvalidParameterException("Reaction particle type dependence must be the same for all particle types, check field REACTION_PARTYPE_DEPENDENT");
			}
		}
		else // if no particle reaction or binding exists in first particle type, default to single mode
		{
			_singleBinding = _binding[parType] ? !_particles[parType]->bindingParDep() : true;
			_singleDynReaction = _dynReaction[parType] ? !_particles[parType]->reactionParDep() : true;
		}
	}

	// Setup the memory for tempState and Jacobians
	if (firstConfigCall)
	{
		_tempState = new double[numDofs()];
		_jacInlet.resize(_convDispOp.axNNodes() * _disc.radNPoints * _disc.nComp, _disc.radNPoints * _disc.nComp);
		// Allocate Jacobian memory; pattern will be set and analyzed in configure()
		_globalJac.resize(numPureDofs(), numPureDofs());
		_globalJacDisc.resize(numPureDofs(), numPureDofs());
	}

	_jacInlet.setZero();

	return transportSuccess && particleConfSuccess && reactionConfSuccess;
}

bool ColumnModel2D::configure(IParameterProvider& paramProvider)
{
	_parameters.clear();

	const bool transportSuccess = _convDispOp.configure(_unitOpIdx, paramProvider, _parameters);

	// Check whether PAR_TYPE_VOLFRAC is required or not
	if ((_disc.nParType > 1) && !paramProvider.exists("PAR_TYPE_VOLFRAC"))
		throw InvalidParameterException("The required parameter \"PAR_TYPE_VOLFRAC\" was not found");

	if (paramProvider.exists("PAR_TYPE_VOLFRAC"))
		_parTypeVolFracMode = readAndRegisterMultiplexParam(paramProvider, _parameters, _parTypeVolFrac, "PAR_TYPE_VOLFRAC", _disc.axNPoints, _disc.radNPoints, _disc.nParType, _unitOpIdx);
	else if (_disc.nParType == 1)
	{
		_parTypeVolFrac.resize(_disc.axNPoints * _disc.radNPoints, 1.0);
		_parTypeVolFracMode = MultiplexMode::Independent;
	}

	if (_disc.nParType * _disc.axNPoints * _disc.radNPoints != _parTypeVolFrac.size())
		throw InvalidParameterException("Number of elements in field PAR_TYPE_VOLFRAC does not match number of particle types times (RAD_POLYDEG+1)*RAD_NELEM * (AX_POLYDEG+1)*AX_NELEM");

	// Check that particle volume fractions sum up to 1.0
	if (_disc.nParType > 0)
	{
		for (unsigned int i = 0; i < _disc.nBulkPoints; ++i)
		{
			const double volFracSum = std::accumulate(_parTypeVolFrac.begin() + i * _disc.nParType, _parTypeVolFrac.begin() + (i + 1) * _disc.nParType, 0.0,
				[](double a, const active& b) -> double { return a + static_cast<double>(b); });
			if (std::abs(1.0 - volFracSum) > 1e-12)
				throw InvalidParameterException("Sum of field PAR_TYPE_VOLFRAC differs from 1.0 (is " + std::to_string(volFracSum) + ") in axial element " + std::to_string(i / _disc.radNElem) + " radial element " + std::to_string(i % _disc.radNElem));
		}
	}

	// Register initial conditions parameters
	registerParam1DArray(_parameters, _initC, [=](bool multi, unsigned int comp) { return makeParamId(hashString("INIT_C"), _unitOpIdx, comp, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep); });
	if (_disc.radNElem > 1)
		registerParam2DArray(_parameters, _initC, [=](bool multi, unsigned int rad, unsigned int comp) { return makeParamId(hashString("INIT_C"), _unitOpIdx, comp, ParTypeIndep, BoundStateIndep, rad, SectionIndep); }, _disc.nComp);

	if (_singleBinding)
	{
		for (unsigned int c = 0; c < _disc.nComp; ++c)
			_parameters[makeParamId(hashString("INIT_CP"), _unitOpIdx, c, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_initCp[c];

		if (_disc.radNElem > 1)
		{
			for (unsigned int r = 0; r < _disc.radNElem; ++r)
			{
				for (unsigned int c = 0; c < _disc.nComp; ++c)
					_parameters[makeParamId(hashString("INIT_CP"), _unitOpIdx, c, ParTypeIndep, BoundStateIndep, r, SectionIndep)] = &_initCp[r * _disc.nComp * _disc.nParType + c];
			}
		}
	}
	else
	{
		registerParam2DArray(_parameters, _initCp, [=](bool multi, unsigned int type, unsigned int comp) { return makeParamId(hashString("INIT_CP"), _unitOpIdx, comp, type, BoundStateIndep, ReactionIndep, SectionIndep); }, _disc.nComp);
		if (_disc.radNElem > 1)
			registerParam3DArray(_parameters, _initCp, [=](bool multi, unsigned int rad, unsigned int type, unsigned int comp) { return makeParamId(hashString("INIT_CP"), _unitOpIdx, comp, type, BoundStateIndep, rad, SectionIndep); }, _disc.nComp, _disc.nParType);
	}

	if (!_binding.empty())
	{
		const unsigned int maxBoundStates = *std::max_element(_disc.strideBound, _disc.strideBound + _disc.nParType);
		std::vector<ParameterId> initParams(maxBoundStates);

		// Register radially independent
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

		// Register radially dependent
		if (_disc.radNElem > 1)
		{
			if (_singleBinding)
			{
				_binding[0]->fillBoundPhaseInitialParameters(initParams.data(), _unitOpIdx, ParTypeIndep);

				for (unsigned int r = 0; r < _disc.radNElem; ++r)
				{
					for (ParameterId& pId : initParams)
						pId.reaction = r;

					active* const iq = _initCs.data() + _disc.nBoundBeforeType[0] + r * _disc.strideBound[_disc.nParType];
					for (unsigned int i = 0; i < _disc.strideBound[0]; ++i)
						_parameters[initParams[i]] = iq + i;
				}
			}
			else
			{
				for (unsigned int r = 0; r < _disc.radNElem; ++r)
				{
					for (unsigned int type = 0; type < _disc.nParType; ++type)
					{
						_binding[type]->fillBoundPhaseInitialParameters(initParams.data(), _unitOpIdx, type);
						for (ParameterId& pId : initParams)
							pId.reaction = r;

						active* const iq = _initCs.data() + _disc.nBoundBeforeType[type] + r * _disc.strideBound[_disc.nParType];
						for (unsigned int i = 0; i < _disc.strideBound[type]; ++i)
							_parameters[initParams[i]] = iq + i;
					}
				}
			}
		}
	}

	// Reconfigure particle model
	bool particleConfSuccess = true;
	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		particleConfSuccess = particleConfSuccess && _particles[parType]->configure(_unitOpIdx, paramProvider, _parameters, _disc.nParType, _disc.nBoundBeforeType, _disc.strideBound[_disc.nParType]);
	}

	// Reconfigure reaction model
	bool dynReactionConfSuccess = true;
	if (_dynReactionBulk && _dynReactionBulk->requiresConfiguration())
	{
		paramProvider.pushScope("reaction_bulk");
		dynReactionConfSuccess = _dynReactionBulk->configure(paramProvider, _unitOpIdx, ParTypeIndep);
		paramProvider.popScope();
	}

	setJacobianPattern(_globalJac, 0, _dynReactionBulk);
	_globalJacDisc = _globalJac;

	// the solver repetitively solves the linear system with a static pattern of the jacobian (set above). 
	// The goal of analyzePattern() is to reorder the nonzero elements of the matrix, such that the factorization step creates less fill-in
	_linearSolver->analyzePattern(_globalJacDisc);


	return transportSuccess && particleConfSuccess && dynReactionConfSuccess;
}


unsigned int ColumnModel2D::threadLocalMemorySize() const CADET_NOEXCEPT
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

unsigned int ColumnModel2D::numAdDirsForJacobian() const CADET_NOEXCEPT
{
	// todo compressed/seeded AD
	return numDofs();
}

void ColumnModel2D::useAnalyticJacobian(const bool analyticJac)
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

void ColumnModel2D::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, const ConstSimulationState& simState, const AdJacobianParams& adJac)
{
	Indexer idxr(_disc);

	_disc.newStaticJac = true;

	// ConvectionDispersionOperator tells us whether flow direction has changed
	if (!_convDispOp.notifyDiscontinuousSectionTransition(t, secIdx));

	for (int parType = 0; parType < _disc.nParType; parType++)
	{
		_particles[parType]->notifyDiscontinuousSectionTransition(t, secIdx);
	}
}

void ColumnModel2D::setFlowRates(active const* in, active const* out) CADET_NOEXCEPT
{
	_convDispOp.setFlowRates(in, out);
}

void ColumnModel2D::reportSolution(ISolutionRecorder& recorder, double const* const solution) const
{
	Exporter expr(_disc, *this, solution);
	recorder.beginUnitOperation(_unitOpIdx, *this, expr);
	recorder.endUnitOperation();
}

void ColumnModel2D::reportSolutionStructure(ISolutionRecorder& recorder) const
{
	Exporter expr(_disc, *this, nullptr);
	recorder.unitOperationStructure(_unitOpIdx, *this, expr);
}


unsigned int ColumnModel2D::requiredADdirs() const CADET_NOEXCEPT
{
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	return _jacobianAdDirs;
#else
	// If CADET_CHECK_ANALYTIC_JACOBIAN is active, we always need the AD directions for the Jacobian
	return numAdDirsForJacobian();
#endif
}

void ColumnModel2D::prepareADvectors(const AdJacobianParams& adJac) const
{
	// Early out if AD is disabled
	if (!adJac.adY)
		return;

	Indexer idxr(_disc);

	const int adDirOffset = adJac.adDirOffset;
	active * adVec = adJac.adY;

	// Dense seeding implemented here, could be improved by figuring out the sparse seeding / colouring

	// Start with diagonal Jacobian element
	for (int eq = 0; eq < numDofs(); ++eq)
	{
		adVec[eq].fillADValue(adDirOffset, 0.0); // Clear previously set directions
		adVec[eq].setADValue(adDirOffset + eq, 1.0);
	}
}

/**
 * @brief Extracts the system Jacobian from AD seed vectors
 * @param [in] adRes Residual vector of AD datatypes with AD seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
void ColumnModel2D::extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset)
{
	Indexer idxr(_disc);

	active const* const adResUnit = adRes + adDirOffset + idxr.offsetC();

	for (int i = 0; i < _convDispOp.axNNodes() * _disc.radNPoints * _disc.nComp; i++)
	{
		for (int j = 0; j < _disc.radNPoints * _disc.nComp; j++)
		{
			const double val = adResUnit[i].getADValue(j + adDirOffset);
			if(std::abs(val) > 1e-15)
				_jacInlet(i, j) = val;
		}
	}

	for (int i = 0; i < _globalJac.rows(); i++)
	{
		for (int j = 0; j < _globalJac.cols(); j++)
		{
			const double val = adResUnit[i].getADValue(j + idxr.offsetC() + adDirOffset);
			if (std::abs(val) > 1e-15)
			_globalJac.coeffRef(i, j) = val;
		}
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
void ColumnModel2D::checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const
{
	Indexer idxr(_disc);

	const int offCp = idxr.offsetCp() - idxr.offsetC();

	double inletJacMaxError = 0.0;
	int row = 0;
	int col = 0;

	active const* const adResUnit = adRes + adDirOffset + idxr.offsetC();

	for (int i = 0; i < _convDispOp.axNNodes() * _disc.radNPoints * _disc.nComp; i++)
	{
		for (int j = 0; j < _disc.radNPoints * _disc.nComp; j++)
		{
			const double localError = std::abs(_jacInlet(i, j) - adResUnit[i].getADValue(j + adDirOffset));
			if (localError > 1e-14)
			{
				inletJacMaxError = std::max(inletJacMaxError, localError);
				row = i;
				col = j;
				LOG(Debug) << "Error in Inlet Jacobian: " << localError << " at row, col: " << i << ", " << j;
			}
		}
	}

	double mainJacMaxError = 0.0;

	for (int i = 0; i < _globalJac.rows(); i++)
	{
		for (int j = 0; j < _globalJac.cols(); j++)
		{
			const double localError = std::abs(_globalJac.coeff(i, j) - adResUnit[i].getADValue(j + idxr.offsetC() + adDirOffset));
			if (localError > 1e-14)
			{
				mainJacMaxError = std::max(mainJacMaxError, localError);
				row = i;
				col = j;

				if(i < offCp && j < offCp)
					LOG(Debug) << "Error in bulk Jacobian: " << localError << " at row, col: " << i << ", " << j;
				else if((i >= offCp && j < offCp) || (i < offCp && j >= offCp))
					LOG(Debug) << "Error in film diffusion Jacobian: " << localError << " at row, col: " << i << ", " << j;
				else
					LOG(Debug) << "Error in particle Jacobian: " << localError << " at row, col: " << i << ", " << j;
			
				LOG(Debug) << "AD Jacobian value: " << adResUnit[i].getADValue(j + idxr.offsetC() + adDirOffset);
				LOG(Debug) << "Analytical Jacobian value: " << _globalJac.coeff(i, j);
			}
		}
	}
}

#endif

int ColumnModel2D::jacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	_factorizeJacobian = true;

	if (_analyticJac)
		return residualImpl<double, double, double, true, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
	else
		return residualWithJacobian(simTime, ConstSimulationState{ simState.vecStateY, nullptr }, nullptr, adJac, threadLocalMem);
}

int ColumnModel2D::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	// Evaluate residual do not compute Jacobian or parameter sensitivities
	return residualImpl<double, double, double, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
}

int ColumnModel2D::residualWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidual);

	// Evaluate residual, use AD for Jacobian if required but do not evaluate parameter derivatives
	return residual(simTime, simState, res, adJac, threadLocalMem, true, false);
}

int ColumnModel2D::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res,
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

template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
int ColumnModel2D::residualImpl(double t, unsigned int secIdx, StateType const* const y, double const* const yDot, ResidualType* const res, util::ThreadLocalStorage& threadLocalMem)
{
	if (wantRes)
	{
		double* const resPtr = reinterpret_cast<double* const>(res);
		Eigen::Map<Eigen::VectorXd> resi(resPtr, numDofs());
		resi.setZero();
	}
	if (wantJac)
	{
		std::fill_n(_globalJac.valuePtr(), _globalJac.nonZeros(), 0.0);
	}

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(std::size_t(0), static_cast<std::size_t>(_disc.axNPoints * _disc.radNPoints * _disc.nParType + 1), [&](std::size_t pblk)
#else
	for (unsigned int pblk = 0; pblk < _disc.axNPoints * _disc.radNPoints * _disc.nParType + 1; ++pblk)
#endif
	{
		if (cadet_unlikely(pblk == 0))
		{
			if (wantJac || _disc.newStaticJac)
			{
				bool success = calcTransportJacobian(secIdx);

				_disc.newStaticJac = false;

				if (cadet_unlikely(!success)) {
					LOG(Error) << "Jacobian pattern did not fit the Jacobian estimation";
				}
			}

			residualBulk<StateType, ResidualType, ParamType, wantJac, wantRes>(t, secIdx, y, yDot, res, threadLocalMem);
		}
		else
		{
			const unsigned int parType = (pblk - 1) / (_disc.axNPoints * _disc.radNPoints);
			const unsigned int colPoint = (pblk - 1) % (_disc.axNPoints * _disc.radNPoints);
			const unsigned int axPoint = std::floor(colPoint / _disc.radNPoints);
			const unsigned int radPoint = colPoint % _disc.radNPoints;
			const unsigned int radElem = std::floor((colPoint - axPoint * _disc.radNPoints) / _disc.radNNodes);

			LinearBufferAllocator tlmAlloc = threadLocalMem.get();
			Indexer idxr(_disc);

			linalg::BandedEigenSparseRowIterator jacIt(_globalJac, idxr.offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colPoint }) - idxr.offsetC());
			model::columnPackingParameters packing
			{
				_parTypeVolFrac[colPoint * _disc.nParType + parType],
				_convDispOp.columnPorosity(radElem),
				ColumnPosition{ _convDispOp.relativeAxialCoordinate(axPoint), _convDispOp.relativeRadialCoordinate(radPoint), 0.0 }
			};

			_particles[parType]->residual(t, secIdx,
				y + idxr.offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colPoint }),
				y + idxr.offsetC() + axPoint * idxr.strideColAxialNode() + radPoint * idxr.strideColRadialNode(),
				yDot ? yDot + idxr.offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colPoint }) : nullptr,
				res ? res + idxr.offsetCp(ParticleTypeIndex{ parType }, ParticleIndex{ colPoint }) : nullptr,
				res ? res + idxr.offsetC() + axPoint * idxr.strideColAxialNode() + radPoint * idxr.strideColRadialNode() : nullptr,
				packing, jacIt, tlmAlloc,
				typename cadet::ParamSens<ParamType>::enabled()
			);
		}
	} CADET_PARFOR_END;

	// Handle inlet DOFs, which are simply copied to res
	for (unsigned int i = 0; i < _disc.nComp * _disc.radNPoints; ++i)
	{
		res[i] = y[i];
	}

	return 0;
}

template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
int ColumnModel2D::residualBulk(double t, unsigned int secIdx, StateType const* yBase, double const* yDotBase, ResidualType* resBase, util::ThreadLocalStorage& threadLocalMem)
{
	if (wantRes)
		_convDispOp.residual(*this, t, secIdx, yBase, yDotBase, resBase, wantJac, typename ParamSens<ParamType>::enabled());

	if (!_dynReactionBulk || (_dynReactionBulk->numReactionsLiquid() == 0))
		return 0;

	// Get offsets
	Indexer idxr(_disc);
	StateType const* y = yBase + idxr.offsetC();
	ResidualType* res = resBase + idxr.offsetC();
	LinearBufferAllocator tlmAlloc = threadLocalMem.get();

	for (unsigned int axNode = 0; axNode < _disc.axNPoints; ++axNode)
	{
		for (unsigned int radNode = 0; radNode < _disc.radNPoints; ++radNode, y += idxr.strideColRadialNode(), res += idxr.strideColRadialNode())
		{
			const double r = _convDispOp.relativeRadialCoordinate(radNode);
			const double z = _convDispOp.relativeAxialCoordinate(axNode);;

			const ColumnPosition colPos{ z, r, 0.0 };
			_dynReactionBulk->residualLiquidAdd(t, secIdx, colPos, y, res, -1.0, tlmAlloc);

			if (wantJac)
			{
				const int rowIdx = axNode * idxr.strideColAxialNode() + radNode * idxr.strideColRadialNode();
				linalg::BandedEigenSparseRowIterator jac(_globalJac, rowIdx);

				// static_cast should be sufficient here, but this statement is also analyzed when wantJac = false
				_dynReactionBulk->analyticJacobianLiquidAdd(t, secIdx, colPos, reinterpret_cast<double const*>(y), -1.0, jac, tlmAlloc);
			}
		}
	}

	return 0;
}

int ColumnModel2D::residualSensFwdWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidualSens);

	// Evaluate residual for all parameters using AD in vector mode and at the same time update the
	// Jacobian (in one AD run, if analytic Jacobians are disabled)
	return residual(simTime, simState, nullptr, adJac, threadLocalMem, true, true);
}

int ColumnModel2D::residualSensFwdAdOnly(const SimulationTime& simTime, const ConstSimulationState& simState, active* const adRes, util::ThreadLocalStorage& threadLocalMem)
{
	BENCH_SCOPE(_timerResidualSens);

	// Evaluate residual for all parameters using AD in vector mode
	return residualImpl<double, active, active, false>(simTime.t, simTime.secIdx, simState.vecStateY, simState.vecStateYdot, adRes, threadLocalMem);
}

int ColumnModel2D::residualSensFwdCombine(const SimulationTime& simTime, const ConstSimulationState& simState,
	const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS, active const* adRes,
	double* const tmp1, double* const tmp2, double* const tmp3)
{
	BENCH_SCOPE(_timerResidualSens);

	// tmp1 stores result of (dF / dy) * s
	// tmp2 stores result of (dF / dyDot) * sDot

	const SimulationTime cst{simTime.t, simTime.secIdx};
	const ConstSimulationState css{nullptr, nullptr};

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(std::size_t(0), static_cast<std::size_t>(yS.size()), [&](std::size_t param)
#else
	for (std::size_t param = 0; param < yS.size(); ++param)
#endif
	{

		// Directional derivative (dF / dy) * s
		multiplyWithJacobian(cst, css, yS[param], 1.0, 0.0, tmp1);

		// Directional derivative (dF / dyDot) * sDot
		multiplyWithDerivativeJacobian(cst, css, ySdot[param], tmp2);

		double* const ptrResS = resS[param];

		BENCH_START(_timerResidualSensPar);

		// Complete sens residual is the sum:
		for (unsigned int i = 0; i < numDofs(); ++i)
		{
			ptrResS[i] = tmp1[i] + tmp2[i] + adRes[i].getADValue(param);
		}

		BENCH_STOP(_timerResidualSensPar);
	} CADET_PARFOR_END;

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
void ColumnModel2D::multiplyWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* yS, double alpha, double beta, double* ret)
{
	Indexer idxr(_disc);

	// Handle identity matrix of inlet DOFs
	for (unsigned int i = 0; i < _disc.nComp * _disc.radNPoints; ++i)
	{
		ret[i] = alpha * yS[i] + beta * ret[i];
	}

	// Main Jacobian
	Eigen::Map<Eigen::VectorXd> ret_vec(ret + idxr.offsetC(), numPureDofs());
	Eigen::Map<const Eigen::VectorXd> yS_vec(yS + idxr.offsetC(), numPureDofs());
	ret_vec = alpha * _globalJac * yS_vec + beta * ret_vec;

	// Map inlet DOFs to the column inlet (first bulk nodes)
	Eigen::Map<const Eigen::VectorXd> yInlet(yS + idxr.offsetC(), _disc.radNPoints * _disc.nComp);
	Eigen::Map<Eigen::VectorXd> retInlet(ret + idxr.offsetC(), _convDispOp.axNNodes() * _disc.radNPoints * _disc.nComp);

	retInlet += _jacInlet * yInlet;
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
void ColumnModel2D::multiplyWithDerivativeJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* sDot, double* ret)
{
	Indexer idxr(_disc);
	std::fill_n(ret, numDofs(), 0.0);

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(std::size_t(0), static_cast<std::size_t>(_disc.nBulkPoints * _disc.nParType + 1), [&](std::size_t idx)
#else
	for (unsigned int idx = 0; idx < _disc.nBulkPoints * _disc.nParType + 1; ++idx)
#endif
	{
		if (cadet_unlikely(idx == 0))
		{
			_convDispOp.multiplyWithDerivativeJacobian(simTime, sDot, ret);
		}
		else
		{
			const unsigned int idxParLoop = idx - 1;
			const unsigned int pblk = idxParLoop % _disc.nBulkPoints;
			const unsigned int type = idxParLoop / _disc.nBulkPoints;

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

void ColumnModel2D::setExternalFunctions(IExternalFunction** extFuns, unsigned int size)
{
	for (IBindingModel* bm : _binding)
	{
		if (bm)
			bm->setExternalFunctions(extFuns, size);
	}
}

unsigned int ColumnModel2D::localOutletComponentIndex(unsigned int port) const CADET_NOEXCEPT
{
	const int radZone = std::floor(port / _disc.radNNodes);

	// Inlets are duplicated so need to be accounted for
	if (static_cast<double>(_convDispOp.currentVelocity(radZone)) >= 0.0)
		// Forward Flow: outlet is last node
		return _disc.nComp * _disc.radNPoints + (_disc.axNPoints - 1) * _disc.nComp * _disc.radNPoints + port * _disc.nComp;
	else
		// Backward flow: Outlet is first node
		return _disc.nComp * _disc.radNPoints + _disc.nComp * port;
}

unsigned int ColumnModel2D::localInletComponentIndex(unsigned int port) const CADET_NOEXCEPT
{
	// Always 0 due to dedicated inlet DOFs
	return _disc.nComp * port;
}

unsigned int ColumnModel2D::localOutletComponentStride(unsigned int port) const CADET_NOEXCEPT
{
	return 1;
}

unsigned int ColumnModel2D::localInletComponentStride(unsigned int port) const CADET_NOEXCEPT
{
	return 1;
}

void ColumnModel2D::expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut)
{
	// @todo Write this function
}

bool ColumnModel2D::setParameter(const ParameterId& pId, double value)
{
	if (pId.unitOperation == _unitOpIdx)
	{
		if (multiplexParameterValue(pId, hashString("PAR_TYPE_VOLFRAC"), _parTypeVolFracMode, _parTypeVolFrac, _disc.axNPoints, _disc.radNPoints, _disc.nParType, value, nullptr))
			return true;
		const int mpIc = multiplexInitialConditions(pId, value, false);
		if (mpIc > 0)
			return true;
		else if (mpIc < 0)
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

		if (model::setParameter(pId, value, std::vector<IDynamicReactionModel*>{ _dynReactionBulk }, true))
			return true;
	}

	const bool result = UnitOperationBase::setParameter(pId, value);

	return result;
}

bool ColumnModel2D::setParameter(const ParameterId& pId, int value)
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
		if (model::setParameter(pId, value, std::vector<IDynamicReactionModel*>{ _dynReactionBulk }, true))
			return true;
	}

	return UnitOperationBase::setParameter(pId, value);
}

bool ColumnModel2D::setParameter(const ParameterId& pId, bool value)
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
		if (model::setParameter(pId, value, std::vector<IDynamicReactionModel*>{ _dynReactionBulk }, true))
			return true;
	}

	return UnitOperationBase::setParameter(pId, value);
}

void ColumnModel2D::setSensitiveParameterValue(const ParameterId& pId, double value)
{
	if (pId.unitOperation == _unitOpIdx)
	{
		if (multiplexParameterValue(pId, hashString("PAR_TYPE_VOLFRAC"), _parTypeVolFracMode, _parTypeVolFrac, _disc.axNPoints, _disc.radNPoints, _disc.nParType, value, &_sensParams))
			return;
		if (multiplexInitialConditions(pId, value, true) != 0)
			return;

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

		if (model::setSensitiveParameterValue(pId, value, _sensParams, std::vector<IDynamicReactionModel*>{ _dynReactionBulk }, true))
			return;
	}

	UnitOperationBase::setSensitiveParameterValue(pId, value);
}

bool ColumnModel2D::setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue)
{
	if (pId.unitOperation == _unitOpIdx)
	{
		if (multiplexParameterAD(pId, hashString("PAR_TYPE_VOLFRAC"), _parTypeVolFracMode, _parTypeVolFrac, _disc.axNPoints, _disc.radNPoints, _disc.nParType, adDirection, adValue, _sensParams))
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

		if (model::setSensitiveParameter(pId, adDirection, adValue, _sensParams, std::vector<IDynamicReactionModel*> { _dynReactionBulk }, true))
		{
			LOG(Debug) << "Found parameter " << pId << " in DynamicBulkReactionModel: Dir " << adDirection << " is set to " << adValue;
			return true;
		}
	}

	const bool result = UnitOperationBase::setSensitiveParameter(pId, adDirection, adValue);

	return result;
}


int ColumnModel2D::Exporter::writeMobilePhase(double* buffer) const
{
	const int blockSize = numMobilePhaseDofs();
	std::copy_n(_idx.c(_data), blockSize, buffer);
	return blockSize;
}

int ColumnModel2D::Exporter::writeSolidPhase(double* buffer) const
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

int ColumnModel2D::Exporter::writeParticleMobilePhase(double* buffer) const
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

int ColumnModel2D::Exporter::writeSolidPhase(unsigned int parType, double* buffer) const
{
	cadet_assert(parType < _disc.nParType);

	const unsigned int stride = _disc.nComp + _disc.strideBound[parType];
	double const* ptr = _data + _idx.offsetCp(ParticleTypeIndex{parType}) + _idx.strideParLiquid();
	for (unsigned int i = 0; i < _disc.axNPoints * _disc.radNPoints; ++i)
	{
		for (unsigned int j = 0; j < _disc.nParPoints[parType]; ++j)
		{
			std::copy_n(ptr, _disc.strideBound[parType], buffer);
			buffer += _disc.strideBound[parType];
			ptr += stride;
		}
	}
	return _disc.nBulkPoints * _disc.nParPoints[parType] * _disc.strideBound[parType];
}

int ColumnModel2D::Exporter::writeParticleMobilePhase(unsigned int parType, double* buffer) const
{
	cadet_assert(parType < _disc.nParType);

	const unsigned int stride = _disc.nComp + _disc.strideBound[parType];
	double const* ptr = _data + _idx.offsetCp(ParticleTypeIndex{parType});
	for (unsigned int i = 0; i < _disc.axNPoints * _disc.radNPoints; ++i)
	{
		for (unsigned int j = 0; j < _disc.nParPoints[parType]; ++j)
		{
			std::copy_n(ptr, _disc.nComp, buffer);
			buffer += _disc.nComp;
			ptr += stride;
		}
	}
	return _disc.nBulkPoints * _disc.nParPoints[parType] * _disc.nComp;
}

int ColumnModel2D::Exporter::writeInlet(unsigned int port, double* buffer) const
{
	cadet_assert(port < _disc.radNPoints);
	std::copy_n(_data + port * _disc.nComp, _disc.nComp, buffer);
	return _disc.nComp;
}

int ColumnModel2D::Exporter::writeInlet(double* buffer) const
{
	std::copy_n(_data, _disc.nComp * _disc.radNPoints, buffer);
	return _disc.nComp * _disc.radNPoints;
}

int ColumnModel2D::Exporter::writeOutlet(unsigned int port, double* buffer) const
{
	cadet_assert(port < _disc.radNPoints);

	const int radZone = std::floor(port / _disc.radNNodes);

	if (_model._convDispOp.currentVelocity(radZone) >= 0)
		std::copy_n(&_idx.c(_data, _disc.axNPoints - 1, port, 0), _disc.nComp, buffer);
	else
		std::copy_n(&_idx.c(_data, 0, port, 0), _disc.nComp, buffer);

	return _disc.nComp;
}

int ColumnModel2D::Exporter::writeOutlet(double* buffer) const
{
	for (int i = 0; i < _disc.radNPoints; ++i)
	{
		writeOutlet(i, buffer);
		buffer += _disc.nComp;
	}
	return _disc.nComp * _disc.radNPoints;
}


void registerColumnModel2D(std::unordered_map<std::string, std::function<IUnitOperation*(UnitOpIdx, IParameterProvider&)>>& models)
{
	models[ColumnModel2D::identifier()] = [](UnitOpIdx uoId, IParameterProvider&) { return new ColumnModel2D(uoId); };
	models["COLUMN_MODEL_2D"] = [](UnitOpIdx uoId, IParameterProvider&) { return new ColumnModel2D(uoId); };
}

}  // namespace model

}  // namespace cadet
