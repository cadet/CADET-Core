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

#include "model/particle/GeneralRateParticle.hpp"
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
	 * @brief Creates a GeneralRateParticle
	 */
	GeneralRateParticle::GeneralRateParticle() : _parDiffOp(nullptr)
	{
	}

	GeneralRateParticle::~GeneralRateParticle() CADET_NOEXCEPT
	{
		delete _parDiffOp;
	}

	bool GeneralRateParticle::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp)
	{
		_parTypeIdx = parTypeIdx;
		_nComp = nComp;

		// ==== Construct and configure particle transport and discretization

		paramProvider.pushScope("discretization");

		if(paramProvider.exists("PAR_SPATIAL_METHOD"))
		{
			const std::string parSpatialMethod = paramProvider.getString("PAR_SPATIAL_METHOD");
			if (parSpatialMethod != "DG")
				throw InvalidParameterException("Unsupported PAR_SPATIAL_METHOD '" + parSpatialMethod + "' for GeneralRateParticle. Only 'DG' is supported.");

			_parDiffOp = new parts::ParticleDiffusionOperatorDG();
		}
		else
			_parDiffOp = new parts::ParticleDiffusionOperatorDG();

		paramProvider.popScope();

		const bool particleTransportConfigSuccess = _parDiffOp->configureModelDiscretization(paramProvider, helper, nComp, parTypeIdx, nParType, strideBulkComp);

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
		bindingConfSuccess = _binding->configureModelDiscretization(paramProvider, _nComp, _parDiffOp->nBound(), _parDiffOp->offsetBoundComp());

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
			reactionConfSuccess = _dynReaction->configureModelDiscretization(paramProvider, _nComp, _parDiffOp->nBound(), _parDiffOp->offsetBoundComp()) && reactionConfSuccess;
		}

		return particleTransportConfigSuccess && bindingConfSuccess && reactionConfSuccess;
	}

	bool GeneralRateParticle::configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound)
	{
		//// Done in the unit operation: Register initial conditions parameters
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

		// Reconfigure particle transport and discretization
		const bool parTransportConfigSuccess = _parDiffOp->configure(unitOpIdx, paramProvider, parameters, nParType, nBoundBeforeType, nTotalBound, _binding->reactionQuasiStationarity(), _binding->hasDynamicReactions());

		return parTransportConfigSuccess && bindingConfSuccess && dynReactionConfSuccess;
	}

	bool GeneralRateParticle::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active const* const filmDiff, active const* const poreAccessFactor)
	{
		return _parDiffOp->notifyDiscontinuousSectionTransition(t, secIdx, filmDiff, poreAccessFactor);
	}

	int GeneralRateParticle::getParticleCoordinates(double* coords) const
	{
		return _parDiffOp->getParticleCoordinates(coords);
	}

	parts::cell::CellParameters GeneralRateParticle::makeCellResidualParams(int const* qsReaction, unsigned int const* nBound) const
	{
		return parts::cell::CellParameters
		{
			_nComp,
			nBound,
			_parDiffOp->offsetBoundComp(),
			_parDiffOp->strideBound(),
			qsReaction,
			getPorosity(),
			getPoreAccessfactor(),
			_binding,
			(_dynReaction && (_dynReaction->numReactionsCombined() > 0)) ? _dynReaction : nullptr
		};
	}

	int GeneralRateParticle::residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity)
	{
		if (resPar)
		{
			if (jacIt.data())
				return residualImpl<double, double, double, true, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, colPos, jacIt, tlmAlloc);
			else
				return residualImpl<double, double, double, false, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, colPos, jacIt, tlmAlloc);
		}
		else if (jacIt.data())
			return residualImpl<double, double, double, true, false>(t, secIdx, yPar, yBulk, yDotPar, resPar, colPos, jacIt, tlmAlloc);
		else
			return -1;
	}
	int GeneralRateParticle::residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity)
	{
		if (jacIt.data())
			return residualImpl<double, active, active, true, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, colPos, jacIt, tlmAlloc);
		else
			return residualImpl<double, active, active, false, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, colPos, jacIt, tlmAlloc);
	}
	int GeneralRateParticle::residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity)
	{
		if (jacIt.data())
			return residualImpl<active, active, double, true, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, colPos, jacIt, tlmAlloc);
		else
			return residualImpl<active, active, double, false, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, colPos, jacIt, tlmAlloc);
	}
	int GeneralRateParticle::residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity)
	{
		if (jacIt.data())
			return residualImpl<active, active, active, true, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, colPos, jacIt, tlmAlloc);
		else
			return residualImpl<active, active, active, false, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, colPos, jacIt, tlmAlloc);
	}

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
	int GeneralRateParticle::residualImpl(double t, unsigned int secIdx, StateType const* yPar, StateType const* yBulk, double const* yDotPar, ResidualType* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc)
	{
		int const* const qsBinding = _binding->reactionQuasiStationarity();
		const parts::cell::CellParameters cellResParams = makeCellResidualParams(qsBinding, _parDiffOp->nBound());

		linalg::BandedEigenSparseRowIterator jacBase = jacIt;

		// Handle time derivatives, binding, dynamic reactions: residualKernel computes discrete point wise,
		// so we loop over each discrete particle point
		for (unsigned int par = 0; par < nDiscPoints(); ++par)
		{
			// local pointers to current particle node, needed in residualKernel
			StateType const* local_y = yPar + par * stridePoint();
			double const* local_yDot = yDotPar ? yDotPar + par * stridePoint() : nullptr;
			ResidualType* local_res = resPar ? resPar + par * stridePoint() : nullptr;

			// r (particle) coordinate of current node (particle radius normed to 1) - needed in externally dependent adsorption kinetic
			colPos.particle = relativeCoordinate(par);

			if (wantRes)
				parts::cell::residualKernel<StateType, ResidualType, ParamType, parts::cell::CellParameters, linalg::BandedEigenSparseRowIterator, wantJac, true>(
					t, secIdx, colPos, local_y, local_yDot, local_res, jacIt, cellResParams, tlmAlloc
				);
			else
				parts::cell::residualKernel<StateType, ResidualType, ParamType, parts::cell::CellParameters, linalg::BandedEigenSparseRowIterator, wantJac, false, false>(
					t, secIdx, colPos, local_y, local_yDot, local_res, jacIt, cellResParams, tlmAlloc
				);

			// Move rowiterator to next particle node
			jacIt += stridePoint();
		}

		ResidualType* wantResPtr = wantRes ? resPar : nullptr;
		linalg::BandedEigenSparseRowIterator jacJojo = wantJac ? jacBase : linalg::BandedEigenSparseRowIterator{};
		return _parDiffOp->residual(t, secIdx, yPar, yBulk, yDotPar, wantResPtr, jacJojo, typename ParamSens<ParamType>::enabled());
	}

	unsigned int GeneralRateParticle::jacobianNNZperParticle() const
	{
		return _parDiffOp->jacobianNNZperParticle();
	}

	int GeneralRateParticle::calcStaticAnaParticleDiffJacobian(const int secIdx, const int colNode, const int offsetLocalCp, Eigen::SparseMatrix<double, RowMajor>& globalJac)
	{
		return _parDiffOp->calcStaticAnaParticleDiffJacobian(secIdx, colNode, offsetLocalCp, globalJac);
	}
	
	int GeneralRateParticle::calcFilmDiffJacobian(unsigned int secIdx, const int offsetCp, const int offsetC, const int nBulkPoints, const int nParType, const double colPorosity, const active* const parTypeVolFrac, Eigen::SparseMatrix<double, RowMajor>& globalJac, bool outliersOnly)
	{
		return _parDiffOp->calcFilmDiffJacobian(secIdx, offsetCp, offsetC, nBulkPoints, nParType, colPorosity, parTypeVolFrac, globalJac, outliersOnly);
	}

	bool GeneralRateParticle::setParameter(const ParameterId& pId, double value)
	{
		return _parDiffOp->setParameter(pId, value);
	}

	bool GeneralRateParticle::setParameter(const ParameterId& pId, int value)
	{
		return _parDiffOp->setParameter(pId, value);
	}

	bool GeneralRateParticle::setParameter(const ParameterId& pId, bool value)
	{
		return _parDiffOp->setParameter(pId, value);
	}

	bool GeneralRateParticle::setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& pId, double value)
	{
		return _parDiffOp->setSensitiveParameterValue(sensParams, pId, value);
	}

	bool GeneralRateParticle::setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue)
	{
		return _parDiffOp->setSensitiveParameter(sensParams, pId, adDirection, adValue);
	}

	std::unordered_map<ParameterId, double> GeneralRateParticle::getAllParameterValues(std::unordered_map<ParameterId, double>& data) const
	{
		model::getAllParameterValues(data, std::vector<IParameterStateDependence*>{_parDiffOp->getParDepSurfDiffusion()}, _parDiffOp->singleParDepSurfDiffusion());

		return data;
	}

	double GeneralRateParticle::getParameterDouble(const ParameterId& pId) const
	{
		double val = 0.0;

		if (model::getParameterDouble(pId, std::vector<IParameterStateDependence*>{_parDiffOp->getParDepSurfDiffusion()}, _parDiffOp->singleParDepSurfDiffusion(), val))
			return val;
		else
			return static_cast<double>(false);
	}

	bool GeneralRateParticle::hasParameter(const ParameterId& pId) const
	{
		if (model::hasParameter(pId, std::vector<IParameterStateDependence*>{_parDiffOp->getParDepSurfDiffusion()}, _parDiffOp->singleParDepSurfDiffusion()))
			return true;
	}

}  // namespace model

}  // namespace cadet
