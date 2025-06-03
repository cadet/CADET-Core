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
	 * @brief Creates a GeneralRateParticle
	 */
	GeneralRateParticle::GeneralRateParticle()
	{
	}

	GeneralRateParticle::~GeneralRateParticle() CADET_NOEXCEPT
	{
		delete _parDiffOp;
	}

	bool GeneralRateParticle::configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp)
	{
		_parDiffOp = new parts::ParticleDiffusionOperatorDG();

		return _parDiffOp->configureModelDiscretization(paramProvider, helper, nComp, parTypeIdx, nParType, strideBulkComp);
	}

	bool GeneralRateParticle::configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound)
	{
		return _parDiffOp->configure(unitOpIdx, paramProvider, parameters, nParType, nBoundBeforeType, nTotalBound);
	}

	/**
	 * @brief Notifies the operator that a discontinuous section transition is in progress
	 * @param [in] t Current time point
	 * @param [in] secIdx Index of the new section that is about to be integrated
	 * @return @c true if flow direction has changed, otherwise @c false
	 */
	bool GeneralRateParticle::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active const* const filmDiff, active const* const poreAccessFactor)
	{
		return _parDiffOp->notifyDiscontinuousSectionTransition(t, secIdx, filmDiff, poreAccessFactor);
	}
	/**
	 * @brief calculates the physical radial/particle coordinates of the DG discretization with double! interface nodes
	 */
	int GeneralRateParticle::getParticleCoordinates(double* coords) const
	{
		return _parDiffOp->getParticleCoordinates(coords);
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
	int GeneralRateParticle::residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity)
	{
		return residualImpl<double, double, double, wantJac, wantRes>(t, secIdx, yPar, yBulk, yDotPar, resPar, colPos, jacIt, tlmAlloc);
	}

	template<bool wantJac, bool wantRes>
	int GeneralRateParticle::residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity)
	{
		return residualImpl<active, active, double, wantJac, wantRes>(t, secIdx, yPar, yBulk, yDotPar, resPar, colPos, jacIt, tlmAlloc);
	}

	template<bool wantJac, bool wantRes>
	int GeneralRateParticle::residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity)
	{
		return residualImpl<double, active, active, wantJac, wantRes>(t, secIdx, yPar, yBulk, yDotPar, resPar, colPos, jacIt, tlmAlloc);
	}

	template<bool wantJac, bool wantRes>
	int GeneralRateParticle::residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity)
	{
		return residualImpl<active, active, active, wantJac, wantRes>(t, secIdx, yPar, yBulk, yDotPar, resPar, colPos, jacIt, tlmAlloc);
	}

	template int GeneralRateParticle::residual<true, true>(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity);
	template int GeneralRateParticle::residual<false, true>(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity);
	template int GeneralRateParticle::residual<true, false>(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity);

	template int GeneralRateParticle::residual<true, true>(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity);
	template int GeneralRateParticle::residual<false, true>(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity);

	template int GeneralRateParticle::residual<true, true>(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity);
	template int GeneralRateParticle::residual<false, true>(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity);

	template int GeneralRateParticle::residual<true, true>(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity);
	template int GeneralRateParticle::residual<false, true>(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
	int GeneralRateParticle::residualImpl(double t, unsigned int secIdx, StateType const* yPar, StateType const* yBulk, double const* yDotPar, ResidualType* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc)
	{
		return _parDiffOp->residualImpl<StateType, ResidualType, ParamType, wantJac, wantRes>(t, secIdx, yPar, yBulk, yDotPar, resPar, colPos, jacIt, tlmAlloc);
	}

	unsigned int GeneralRateParticle::calcParDiffNNZ()
	{
		return _parDiffOp->calcParDiffNNZ();
	}

	/**
	 * @brief analytically calculates the static (per section) particle diffusion Jacobian
	 * @return 1 if jacobain calculation fits the predefined pattern of the jacobian, 0 if not.
	 */
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

}  // namespace parts

}  // namespace model

}  // namespace cadet
