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

		std::ostringstream parTypeIdxString;
		parTypeIdxString << std::setfill('0') << std::setw(3) << std::setprecision(0) << _parTypeIdx;
		paramProvider.pushScope("particle_type_" + parTypeIdxString.str());

		// ==== Construct and configure binding model

		_binding = nullptr;
		std::string bindModelName = "NONE";
		bool bindingConfSuccess = true;

		if (paramProvider.exists("ADSORPTION_MODEL"))
			bindModelName = paramProvider.getString("ADSORPTION_MODEL");

		_binding = helper.createBindingModel(bindModelName);
		if (!_binding)
			throw InvalidParameterException("Unknown binding model " + bindModelName);

		if (bindModelName == "NONE")
			_nBound = std::make_shared<unsigned int[]>(_nComp, 0);
		else
		{
			std::vector<int> nBound = paramProvider.getIntArray("NBOUND");
			if (nBound.size() != _nComp)
				throw InvalidParameterException("Field NBOUND does not contain NCOMP = " + std::to_string(_nComp) + " entries for particle type " + std::to_string(_parTypeIdx));

			if (!_nBound)
				_nBound = std::make_shared<unsigned int[]>(_nComp);
			std::copy_n(nBound.begin(), _nComp, _nBound.get());
		}

		// ==== Construct and configure particle transport and discretization

		paramProvider.pushScope("discretization");

		if (paramProvider.exists("SPATIAL_METHOD"))
		{
			const std::string parSpatialMethod = paramProvider.getString("SPATIAL_METHOD");
			if (parSpatialMethod != "DG")
				throw InvalidParameterException("Unsupported SPATIAL_METHOD '" + parSpatialMethod + "' for GeneralRateParticle. Only 'DG' is supported for now.");

			_parDiffOp = new parts::ParticleDiffusionOperatorDG();
		}
		else
			_parDiffOp = new parts::ParticleDiffusionOperatorDG();

		paramProvider.popScope();

		_parDiffOp->setNBound(_nBound);

		const bool particleTransportConfigSuccess = _parDiffOp->configureModelDiscretization(paramProvider, helper, nComp, parTypeIdx, nParType, strideBulkComp);

		// ==== Configure binding model
		{
			MultiplexedScopeSelector scopeGuard(paramProvider, "adsorption", _binding->usesParamProviderInDiscretizationConfig());
			bindingConfSuccess = _binding->configureModelDiscretization(paramProvider, _nComp, _nBound.get(), _parDiffOp->offsetBoundComp());
			_bindingParDep = paramProvider.exists("BINDING_PARTYPE_DEPENDENT") ? paramProvider.getBool("BINDING_PARTYPE_DEPENDENT") : true;
		}

		// ==== Construct and configure dynamic reaction model

		_reactionParDep = true;
		bool reactionConfSuccess = true;

		_reaction.clearDynamicReactionModels();
		bool hasReaction = false;

		unsigned int noOfSet = 0;
		if (paramProvider.exists("NREAC_CROSS_PHASE"))
		{
			hasReaction = true;
			int nReac = paramProvider.getInt("NREAC_CROSS_PHASE");

			reactionConfSuccess = _reaction.configureDiscretization("cross_phase",
				nReac,
				_nComp,
				_nBound.get(),
				&noOfSet,
				paramProvider,
				helper) && reactionConfSuccess;

		}
		if (paramProvider.exists("NREAC_LIQUID"))
		{
			hasReaction = true;
			int nReac = paramProvider.getInt("NREAC_LIQUID");

			reactionConfSuccess = _reaction.configureDiscretization("liquid",
				nReac,
				_nComp,
				_nBound.get(),
				&noOfSet,
				paramProvider,
				helper) && reactionConfSuccess;
		}
		if (paramProvider.exists("NREAC_SOLID"))
		{
			hasReaction = true;
			int nReac = paramProvider.getInt("NREAC_SOLID");

			reactionConfSuccess = _reaction.configureDiscretization("solid",
				nReac,
				_nComp,
				_nBound.get(),
				&noOfSet,
				paramProvider,
				helper) && reactionConfSuccess;

		}

		if (!hasReaction)
		{
			_reaction.empty();
		}

		paramProvider.popScope(); // particle_type_{:03}

		return particleTransportConfigSuccess && bindingConfSuccess && reactionConfSuccess;
	}

	bool GeneralRateParticle::configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound)
	{
		std::ostringstream parTypeIdxString;
		parTypeIdxString << std::setfill('0') << std::setw(3) << std::setprecision(0) << _parTypeIdx;
		paramProvider.pushScope("particle_type_" + parTypeIdxString.str());

		// Reconfigure binding model
		bool bindingConfSuccess = true;
		if (_binding)
		{
			MultiplexedScopeSelector scopeGuard(paramProvider, "adsorption", _binding->requiresConfiguration());
			bindingConfSuccess = _binding->configure(paramProvider, unitOpIdx, ParTypeIndep);
		}

		// Reconfigure reaction model
		bool dynReactionConfSuccess = true;

		if (paramProvider.exists("NREAC_CROSS_PHASE"))
			dynReactionConfSuccess = _reaction.configure("cross_phase", _parTypeIdx, unitOpIdx, paramProvider) && dynReactionConfSuccess;
		if (paramProvider.exists("NREAC_LIQUID"))
			dynReactionConfSuccess = _reaction.configure("liquid", _parTypeIdx, unitOpIdx, paramProvider) && dynReactionConfSuccess;
		if (paramProvider.exists("NREAC_SOLID"))
			dynReactionConfSuccess = _reaction.configure("solid", _parTypeIdx, unitOpIdx, paramProvider) && dynReactionConfSuccess;

		// Reconfigure particle transport and discretization
		const bool parTransportConfigSuccess = _parDiffOp->configure(unitOpIdx, paramProvider, parameters, nParType, nBoundBeforeType, nTotalBound, _binding->reactionQuasiStationarity());
		
		paramProvider.popScope(); // particle_type_{:03}

		return parTransportConfigSuccess && bindingConfSuccess && dynReactionConfSuccess;
	}

	bool GeneralRateParticle::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx)
	{
		return _parDiffOp->notifyDiscontinuousSectionTransition(t, secIdx);
	}

	int GeneralRateParticle::writeParticleCoordinates(double* coords) const
	{
		return _parDiffOp->writeParticleCoordinates(coords);
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
			getPoreAccessFactor(),
			_binding,
			&_reaction
		};
	}

	int GeneralRateParticle::residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, double* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity)
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
	int GeneralRateParticle::residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity)
	{
		if (jacIt.data())
			return residualImpl<double, active, active, true, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
		else
			return residualImpl<double, active, active, false, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
	}
	int GeneralRateParticle::residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity)
	{
		if (jacIt.data())
			return residualImpl<active, active, double, true, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
		else
			return residualImpl<active, active, double, false, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
	}
	int GeneralRateParticle::residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity)
	{
		if (jacIt.data())
			return residualImpl<active, active, active, true, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
		else
			return residualImpl<active, active, active, false, true>(t, secIdx, yPar, yBulk, yDotPar, resPar, resBulk, packing, jacIt, tlmAlloc);
	}

	template <typename StateType, typename ResidualType, typename ParamType, bool wantNonLinJac, bool wantRes>
	int GeneralRateParticle::residualImpl(double t, unsigned int secIdx, StateType const* yPar, StateType const* yBulk, double const* yDotPar, ResidualType* resPar, ResidualType* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc)
	{
		int const* const qsBinding = _binding->reactionQuasiStationarity();
		const parts::cell::CellParameters cellResParams = makeCellResidualParams(qsBinding, _nBound.get());

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
			packing.colPos.particle = relativeCoordinate(par);

			if (wantRes)
				parts::cell::residualKernel<StateType, ResidualType, ParamType, parts::cell::CellParameters, linalg::BandedEigenSparseRowIterator, wantNonLinJac, true>(
					t, secIdx, packing.colPos, local_y, local_yDot, local_res, jacIt, cellResParams, tlmAlloc
				);
			else
				parts::cell::residualKernel<StateType, ResidualType, ParamType, parts::cell::CellParameters, linalg::BandedEigenSparseRowIterator, wantNonLinJac, false, false>(
					t, secIdx, packing.colPos, local_y, local_yDot, local_res, jacIt, cellResParams, tlmAlloc
				);

			// Move rowiterator to next particle node
			jacIt += stridePoint();
		}

		// particle diffusion, including film diffusion boundary condition
		ResidualType* wantResPtr = wantRes ? resPar : nullptr;
		linalg::BandedEigenSparseRowIterator jacSafe = wantNonLinJac ? jacBase : linalg::BandedEigenSparseRowIterator{};
		_parDiffOp->residual(t, secIdx, yPar, yBulk, yDotPar, wantResPtr, jacSafe, typename ParamSens<ParamType>::enabled());

		if (wantRes)
		{
			// film diffusion bulk eq. term
			active const* const filmDiff = _parDiffOp->getFilmDiffusion(secIdx);
			const ParamType invBetaC = 1.0 / static_cast<ParamType>(packing.colPorosity) - 1.0;
			const ParamType jacCF_val = invBetaC * static_cast<ParamType>(surfaceToVolumeRatio());
			const ParamType jacPF_val = -1.0 / static_cast<ParamType>(getPorosity());

			// Add flux to column void / bulk volume
			for (unsigned int comp = 0; comp < _nComp; ++comp)
			{
				// + 1/Beta^c * (surfaceToVolumeRatio^p_j) * d_j * (k_f * [c^b - c^p])
				resBulk[comp] += static_cast<ParamType>(filmDiff[comp]) * jacCF_val * static_cast<ParamType>(packing.parTypeVolFrac) * (yBulk[comp] - yPar[(nDiscPoints() - 1) * stridePoint() + comp]);
			}
		}

		return 0;
	}

	unsigned int GeneralRateParticle::jacobianNNZperParticle() const
	{
		const int bindingNNZ = _parDiffOp->nDiscPoints() * (_parDiffOp->strideBound() + _nComp) * (_parDiffOp->strideBound() + _nComp);

		return _parDiffOp->jacobianNNZperParticle() + bindingNNZ;
	}

	int GeneralRateParticle::calcParticleDiffJacobian(const int secIdx, const int colNode, const int offsetLocalCp, Eigen::SparseMatrix<double, RowMajor>& globalJac)
	{
		return _parDiffOp->calcParticleDiffJacobian(secIdx, colNode, offsetLocalCp, globalJac);
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
		model::getAllParameterValues(data, std::vector<IParameterStateDependence*>{_parDiffOp->getParDepSurfDiffusion()}, _parDiffOp->paramDepSurfDiffParTypeIndep());
		
		return data;
	}

	double GeneralRateParticle::getParameterDouble(const ParameterId& pId) const
	{
		double val = 0.0;

		if (model::getParameterDouble(pId, std::vector<IParameterStateDependence*>{_parDiffOp->getParDepSurfDiffusion()}, _parDiffOp->paramDepSurfDiffParTypeIndep(), val))
			return val;
		else
			return static_cast<double>(false);
	}

	bool GeneralRateParticle::hasParameter(const ParameterId& pId) const
	{
		if (model::hasParameter(pId, std::vector<IParameterStateDependence*>{_parDiffOp->getParDepSurfDiffusion()}, _parDiffOp->paramDepSurfDiffParTypeIndep()))
			return true;

		return false;
	}

	void registerGeneralRateParticleModel(std::unordered_map<std::string, std::function<model::IParticleModel* ()>>& particles)
	{
		particles[GeneralRateParticle::identifier()] = []() { return new GeneralRateParticle(); };
	}

}  // namespace model

}  // namespace cadet
