// SPDX-License-Identifier: AGPL-3.0-or-later
// =================================================================================
//  CADET
//
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Affero General Public
//  License v3.0 (or, at your option, any later version).
// =================================================================================

#include "model/particle/LiquidEquilibriumParticle.hpp"

#include "cadet/Exceptions.hpp"
#include "BindingModelFactory.hpp"
#include "ReactionModelFactory.hpp"
#include "ParamReaderHelper.hpp"
#include "ParamReaderScopes.hpp"
#include "AdUtils.hpp"
#include "SimulationTypes.hpp"
#include "model/parts/BindingCellKernel.hpp"
#include "SensParamUtil.hpp"
#include "ConfigurationHelper.hpp"
#include "model/BindingModel.hpp"
#include "model/ReactionModel.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <algorithm>
#include <numeric>
#include <sstream>
#include <iomanip>

namespace cadet
{

namespace model
{

	namespace
	{
		inline void bindingFlux(IBindingModel* binding, double t, unsigned int secIdx, const ColumnPosition& colPos, active const* yCp, active const* yQ, active* res, LinearBufferAllocator buffer, WithParamSensitivity)
		{
			binding->flux(t, secIdx, colPos, yQ, yCp, res, buffer, WithParamSensitivity());
		}

		inline void bindingFlux(IBindingModel* binding, double t, unsigned int secIdx, const ColumnPosition& colPos, active const* yCp, active const* yQ, active* res, LinearBufferAllocator buffer, WithoutParamSensitivity)
		{
			binding->flux(t, secIdx, colPos, yQ, yCp, res, buffer, WithoutParamSensitivity());
		}

		inline void bindingFlux(IBindingModel* binding, double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yCp, double const* yQ, active* res, LinearBufferAllocator buffer, WithParamSensitivity)
		{
			binding->flux(t, secIdx, colPos, yQ, yCp, res, buffer);
		}

		inline void bindingFlux(IBindingModel* binding, double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yCp, double const* yQ, double* res, LinearBufferAllocator buffer, WithoutParamSensitivity)
		{
			binding->flux(t, secIdx, colPos, yQ, yCp, res, buffer);
		}
	}

	LiquidEquilibriumParticle::LiquidEquilibriumParticle() : _zero(0.0), _strideBulkComp(1), _bulkRowOffset(0), _boundOffset(nullptr), _strideBound(0)
	{
	}

	LiquidEquilibriumParticle::~LiquidEquilibriumParticle() CADET_NOEXCEPT
	{
		delete[] _boundOffset;
	}

	bool LiquidEquilibriumParticle::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp)
	{
		_strideBulkComp = strideBulkComp;
		_parTypeIdx = parTypeIdx;
		_nComp = nComp;

		// Default to the layout of a column node, where the liquid phase directly precedes the bound states
		_bulkRowOffset = nComp;

		std::ostringstream parTypeIdxString;
		parTypeIdxString << std::setfill('0') << std::setw(3) << std::setprecision(0) << _parTypeIdx;
		paramProvider.pushScope("particle_type_" + parTypeIdxString.str());

		std::vector<int> nBound = paramProvider.getIntArray("NBOUND");
		if (nBound.size() != _nComp)
			throw InvalidParameterException("Field NBOUND does not contain NCOMP = " + std::to_string(_nComp) + " entries for particle type " + std::to_string(_parTypeIdx));

		if (!_nBound)
			_nBound = std::make_shared<unsigned int[]>(_nComp);
		std::copy_n(nBound.begin(), _nComp, _nBound.get());

		// Precompute offsets and total number of bound states
		if (!_boundOffset)
			_boundOffset = new unsigned int[_nComp];

		_boundOffset[0] = 0u;
		for (unsigned int i = 1; i < _nComp; ++i)
			_boundOffset[i] = _boundOffset[i - 1] + _nBound[i - 1];

		_strideBound = std::accumulate(_nBound.get(), _nBound.get() + _nComp, 0u);

		// ==== Construct and configure binding model

		_binding = nullptr;
		std::string bindModelName = "NONE";
		bool bindingConfSuccess = true;

		if (paramProvider.exists("ADSORPTION_MODEL"))
			bindModelName = paramProvider.getString("ADSORPTION_MODEL");

		_binding = helper.createBindingModel(bindModelName);
		if (!_binding)
			throw InvalidParameterException("Unknown binding model " + bindModelName);

		// Unlike the particle models that have pores, NBOUND is kept without a binding model, which
		// leaves bound states that solid phase and cross phase reactions can act on

		{
			MultiplexedScopeSelector scopeGuard(paramProvider, "adsorption", _binding->usesParamProviderInDiscretizationConfig());
			bindingConfSuccess = _binding->configureModelDiscretization(paramProvider, _nComp, _nBound.get(), _boundOffset);
			_bindingParDep = paramProvider.exists("BINDING_PARTYPE_DEPENDENT") ? paramProvider.getBool("BINDING_PARTYPE_DEPENDENT") : true;
		}

		// ==== Construct and configure reaction models
		// A reaction acting on the liquid phase belongs to the hosting unit operation, which owns that
		// phase, so only the phases this particle owns are configured here.

		_reactionParDep = true;
		bool reactionConfSuccess = true;

		_reaction.clearDynamicReactionModels();

		if (paramProvider.exists("NREAC_CROSS_PHASE"))
			reactionConfSuccess = _reaction.configureDiscretization("cross_phase", paramProvider.getInt("NREAC_CROSS_PHASE"), _nComp, _nBound.get(), _boundOffset, paramProvider, helper) && reactionConfSuccess;

		if (paramProvider.exists("NREAC_SOLID"))
			reactionConfSuccess = _reaction.configureDiscretization("solid", paramProvider.getInt("NREAC_SOLID"), _nComp, _nBound.get(), _boundOffset, paramProvider, helper) && reactionConfSuccess;

		if (paramProvider.exists("REACTION_PARTYPE_DEPENDENT"))
			_reactionParDep = paramProvider.getBool("REACTION_PARTYPE_DEPENDENT");

		paramProvider.popScope(); // particle_type_{:03}

		return bindingConfSuccess && reactionConfSuccess;
	}

	bool LiquidEquilibriumParticle::configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound, bool uoImplementsReqReactions)
	{
		std::ostringstream parTypeIdxString;
		parTypeIdxString << std::setfill('0') << std::setw(3) << std::setprecision(0) << _parTypeIdx;
		paramProvider.pushScope("particle_type_" + parTypeIdxString.str());

		// The particle has neither a geometry nor transport parameters, so there is nothing to read here
		// beyond the binding and reaction models.

		bool bindingConfSuccess = true;
		if (_binding && _binding->requiresConfiguration())
		{
			MultiplexedScopeSelector scopeGuard(paramProvider, "adsorption", _binding->requiresConfiguration());
			bindingConfSuccess = _binding->configure(paramProvider, unitOpIdx, _bindingParDep ? static_cast<ParticleTypeIdx>(_parTypeIdx) : ParTypeIndep);
		}

		bool reactionConfSuccess = true;
		if (paramProvider.exists("NREAC_CROSS_PHASE"))
			reactionConfSuccess = _reaction.configure("cross_phase", _parTypeIdx, unitOpIdx, paramProvider, uoImplementsReqReactions) && reactionConfSuccess;
		if (paramProvider.exists("NREAC_SOLID"))
			reactionConfSuccess = _reaction.configure("solid", _parTypeIdx, unitOpIdx, paramProvider, uoImplementsReqReactions) && reactionConfSuccess;

		paramProvider.popScope(); // particle_type_{:03}

		return bindingConfSuccess && reactionConfSuccess;
	}

	bool LiquidEquilibriumParticle::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx)
	{
		return true;
	}

	int LiquidEquilibriumParticle::writeParticleCoordinates(double* coords) const
	{
		// The particle is not spatially resolved
		coords[0] = 0.0;
		return 1;
	}

	parts::cell::CellParameters LiquidEquilibriumParticle::makeCellResidualParams(int const* qsReaction, unsigned int const* nBound) const
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
			&_reaction
		};
	}

	template <typename StateType, typename ResidualType, typename ParamType, bool wantNonLinJac, bool wantRes>
	int LiquidEquilibriumParticle::residualImpl(double t, unsigned int secIdx, StateType const* yPar, StateType const* yBulk, double const* yDotPar, ResidualType* resPar, ResidualType* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc)
	{
		int const* const qsBinding = _binding->reactionQuasiStationarity();

		packing.colPos.particle = relativeCoordinate(0);

		// The particle is in rapid equilibrium with the liquid phase it is suspended in, so the binding
		// model is evaluated against that liquid phase rather than a particle liquid phase of its own.
		if (wantRes)
		{
			bindingFlux(_binding, t, secIdx, packing.colPos, yBulk, yPar, resPar, tlmAlloc, typename ParamSens<ParamType>::enabled());

			if (_binding->hasDynamicReactions() && yDotPar)
			{
				for (unsigned int idx = 0; idx < _strideBound; ++idx)
				{
					// Quasi-stationary bound states have no time derivative
					if (qsBinding[idx])
						continue;

					resPar[idx] += yDotPar[idx];
				}
			}

			// Whatever the bound states take up leaves the liquid phase, which is a state of the hosting
			// unit operation. The weight is the one the unit operation hands to every particle type for
			// the film transfer term, so this term is scaled exactly like the flux of a particle that
			// does have a film resistance. Quasi-stationary states contribute too: the amount held by
			// the solid changes regardless of how fast the binding is.
			if (resBulk && yDotPar)
			{
				for (unsigned int comp = 0; comp < _nComp; ++comp)
				{
					double qDotSum = 0.0;
					for (unsigned int bnd = 0; bnd < _nBound[comp]; ++bnd)
						qDotSum += yDotPar[_boundOffset[comp] + bnd];

					const ParamType bulkWeight = packing.bulkFilm.nodeLocal() ? static_cast<ParamType>(*packing.bulkFilm.weight(packing.bulkPoint, comp)) : ParamType(1.0);

					resBulk[comp] += bulkWeight * static_cast<ParamType>(packing.parTypeVolFrac) * qDotSum;
				}
			}
		}

		if (wantNonLinJac)
			_binding->analyticJacobian(t, secIdx, packing.colPos, reinterpret_cast<double const*>(yPar), _bulkRowOffset, jacIt, tlmAlloc);

		return 1;
	}

	int LiquidEquilibriumParticle::residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, double* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity)
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

	int LiquidEquilibriumParticle::residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity)
	{
		if (jacIt.data())
			return residualImpl<double, active, active, true, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
		else
			return residualImpl<double, active, active, false, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
	}

	int LiquidEquilibriumParticle::residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity)
	{
		if (jacIt.data())
			return residualImpl<active, active, double, true, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
		else
			return residualImpl<active, active, double, false, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
	}

	int LiquidEquilibriumParticle::residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity)
	{
		if (jacIt.data())
			return residualImpl<active, active, active, true, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
		else
			return residualImpl<active, active, active, false, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
	}

	void LiquidEquilibriumParticle::setParJacPattern(std::vector<Eigen::Triplet<double>>& tripletList, const unsigned int offsetPar, const unsigned int offsetBulk, unsigned int colNode, unsigned int secIdx, const bulkFilmTransfer& bulkFilm) const
	{
		// Binding and reaction couple every bound state of this type to every other bound state and to
		// every liquid phase component, so both blocks are dense.
		for (unsigned int state = 0; state < _strideBound; ++state)
		{
			for (unsigned int toState = 0; toState < _strideBound; ++toState)
				tripletList.push_back(Eigen::Triplet<double>(offsetPar + state, offsetPar + toState, 0.0));

			for (unsigned int comp = 0; comp < _nComp; ++comp)
				tripletList.push_back(Eigen::Triplet<double>(offsetPar + state, offsetBulk + comp, 0.0));
		}

		// The bound states enter the liquid phase balance through their time derivative. The bulk side
		// operator may spread that over further discrete points, see bulkFilmTransfer.
		for (unsigned int comp = 0; comp < _nComp; ++comp)
		{
			for (int k = 0; k < bulkFilm.stencilSize; ++k)
			{
				const int rowShift = bulkFilm.offsets ? bulkFilm.offset(colNode)[k] * bulkFilm.strideBulkPoint : 0;

				for (unsigned int state = 0; state < _strideBound; ++state)
					tripletList.push_back(Eigen::Triplet<double>(offsetBulk + rowShift + comp, offsetPar + state, 0.0));
			}
		}
	}

	unsigned int LiquidEquilibriumParticle::jacobianNNZperParticle() const
	{
		return _strideBound * _strideBound + 2 * _strideBound * _nComp;
	}

	bool LiquidEquilibriumParticle::setParameter(const ParameterId& pId, double value)
	{
		return false;
	}

	bool LiquidEquilibriumParticle::setParameter(const ParameterId& pId, int value)
	{
		return false;
	}

	bool LiquidEquilibriumParticle::setParameter(const ParameterId& pId, bool value)
	{
		return false;
	}

	bool LiquidEquilibriumParticle::setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& pId, double value)
	{
		return false;
	}

	bool LiquidEquilibriumParticle::setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue)
	{
		return false;
	}

	std::unordered_map<ParameterId, double> LiquidEquilibriumParticle::getAllParameterValues(std::unordered_map<ParameterId, double>& data) const
	{
		return data;
	}

	double LiquidEquilibriumParticle::getParameterDouble(const ParameterId& pId) const
	{
		return static_cast<double>(false);
	}

	bool LiquidEquilibriumParticle::hasParameter(const ParameterId& pId) const
	{
		return false;
	}

	void registerLiquidEquilibriumParticleModel(std::unordered_map<std::string, std::function<model::IParticleModel* ()>>& particles)
	{
		particles[LiquidEquilibriumParticle::identifier()] = []() { return new LiquidEquilibriumParticle(); };
	}

}  // namespace model

}  // namespace cadet
