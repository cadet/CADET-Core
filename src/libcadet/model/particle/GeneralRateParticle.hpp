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

/**
 * @file
 * Defines the general rate particle model as e.g. used in the classical General Rate Model of chromatography
 */

#ifndef LIBCADET_GENERALRATEPARTICLE_HPP_
#define LIBCADET_GENERALRATEPARTICLE_HPP_

#include "model/particle/ParticleModel.hpp"
#include "model/parts/ParticleDiffusionOperatorBase.hpp"
#include "model/BindingModel.hpp"
#include "cadet/StrongTypes.hpp"
#include "ParamIdUtil.hpp"
#include "AutoDiff.hpp"
#include "Memory.hpp"
#include "model/ModelUtils.hpp"
#include "SimulationTypes.hpp"
#include "ParamReaderHelper.hpp"
#include "linalg/BandedEigenSparseRowIterator.hpp"
#include "model/parts/DGToolbox.hpp"
#include "model/ParameterMultiplexing.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <unordered_map>
#include <unordered_set>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>
#include <vector>

using namespace Eigen;

namespace cadet
{

	class IConfigHelper;

namespace model
{

class IDynamicReactionModel;
class IParameterStateDependence;

namespace parts
{
	namespace cell
	{
		struct CellParameters;
	}
}
	/**
	 * @brief General Rate Particle Model
	 * @details Implements the equation
	 *
	 @f[ \begin{align}
	 - \left. \left( \varepsilon^{\mathrm{p}} D^{\mathrm{p}}_{i} \frac{\partial c^{\mathrm{p}}_{i}}{\partial r} + (1 - \varepsilon^{\mathrm{p}}) D^{\mathrm{s}}_{i} \frac{\partial c^{\mathrm{s}}_{i}}{\partial r} \right) \right|_{r=0}
	 &= 0, \\
	 \varepsilon^{\mathrm{p}} \left. \left( \varepsilon^{\mathrm{p}}  D^{\mathrm{p}}_{i} \frac{\partial c^{\mathrm{p}}_{i}}{\partial r} + (1 - \varepsilon^{\mathrm{p}} ) D^{\mathrm{s}}_{i} \frac{\partial c^{\mathrm{s}}_{i}}{\partial r} \right)\right|_{r = R^{\mathrm{p}}_{}}
	 &= k^{\mathrm{f}}_{i} \left. \left( c^{\mathrm{b}}_i - c^{\mathrm{p}}_{i} \right|_{r = R^{\mathrm{p}}_{}} \right)
	 \end{align} @f]
	 * Additionally implements the variants for cylindrical and slab-shaped particles
	 * Methods are described in @cite Breuer2023
	 *
	 * This class does not store the Jacobian. It only fills existing matrices given to its residual() functions.
	 * It assumes that there is no offset to the particle entries in the local state vector
	 */
	class GeneralRateParticle : public IParticleModel
	{
	public:

		GeneralRateParticle();
		~GeneralRateParticle() CADET_NOEXCEPT;

		static inline const char* identifier() CADET_NOEXCEPT { return "GENERAL_RATE_PARTICLE"; }

		bool configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp) override;
		bool configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound) override;

		void updateRadialDisc() { _parDiffOp->updateRadialDisc(); }

		parts::cell::CellParameters makeCellResidualParams(int const* qsReaction, unsigned int const* nBound) const override;

		bool notifyDiscontinuousSectionTransition(double t, unsigned int secIdx) override;

		int residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, double* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity) override;
		int residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity) override;
		int residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity) override;
		int residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity) override;

		double relativeCoordinate(const unsigned int nodeIdx) const CADET_NOEXCEPT  override { return _parDiffOp->relativeCoordinate(nodeIdx); }

		active surfaceToVolumeRatio() const CADET_NOEXCEPT override
		{
			return _parDiffOp->surfaceToVolumeRatio();
		}

		parts::ParticleDiffusionOperatorBase* _parDiffOp;

		inline const active& getPorosity() const CADET_NOEXCEPT  override { return _parDiffOp->getPorosity(); }
		inline const active* getPoreAccessFactor() const CADET_NOEXCEPT  override { return _parDiffOp->getPoreAccessFactor(); }
		inline const active* getFilmDiffusion(const unsigned int secIdx) const CADET_NOEXCEPT { return _parDiffOp->getFilmDiffusion(secIdx);  }

		inline int nDiscPoints() const CADET_NOEXCEPT  override { return _parDiffOp->nDiscPoints(); }
		inline int strideParBlock() const CADET_NOEXCEPT  override { return nDiscPoints() * stridePoint(); }

		int writeParticleCoordinates(double* coords) const override;

		void setParJacPattern(std::vector<Eigen::Triplet<double>>& tripletList, const unsigned int offsetPar, const unsigned int offsetBulk, unsigned int colNode, unsigned int secIdx) const override
		{
			_parDiffOp->setParticleJacobianPattern(tripletList, offsetPar, offsetBulk, colNode, secIdx);
		}

		unsigned int jacobianNNZperParticle() const override;
		int calcParticleDiffJacobian(const int secIdx, const int colNode, const int offsetLocalCp, Eigen::SparseMatrix<double, RowMajor>& globalJac) override;
		int calcFilmDiffJacobian(unsigned int secIdx, const int offsetCp, const int offsetC, const int nBulkPoints, const int nParType, const double colPorosity, const active* const parTypeVolFrac, Eigen::SparseMatrix<double, RowMajor>& globalJac, bool crossDepsOnly = false) override;

		bool setParameter(const ParameterId& pId, double value) override;
		bool setParameter(const ParameterId& pId, int value) override;
		bool setParameter(const ParameterId& pId, bool value) override;
		bool setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue) override;
		bool setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& pId, double value) override;

		bool hasParameter(const ParameterId& pId) const override;
		double getParameterDouble(const ParameterId& pId) const override;
		std::unordered_map<ParameterId, double> getAllParameterValues(std::unordered_map<ParameterId, double>& data) const override;

		bool leanConsistentInitialStateValidity() const override
		{
			if (isSectionDependent(_parDiffOp->parDiffMode()) || isSectionDependent(_parDiffOp->parSurfDiffMode()))
			{
				LOG(Warning) << "Lean consistent initialization is not appropriate for section-dependent pore and surface diffusion in particle type " + std::to_string(_parTypeIdx);
				return false;
			}
			else
				return true;
		}

		bool leanConsistentInitialTimeDerivativeValidity() const override
		{
			if (isSectionDependent(_parDiffOp->parDiffMode()) || isSectionDependent(_parDiffOp->parSurfDiffMode()))
			{
				LOG(Warning) << "Lean consistent initialization is not appropriate for section-dependent pore and surface diffusion in particle type " + std::to_string(_parTypeIdx);
				return false;
			}
			else
				return true;
		}

		virtual bool isParticleLumped() const CADET_NOEXCEPT { return false; }

	protected:

		/**
		 * @brief stride over one discrete point
		 */
		inline int stridePoint() const CADET_NOEXCEPT { return static_cast<int>(_nComp) + _parDiffOp->strideBound(); }

		template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
		int residualImpl(double t, unsigned int secIdx, StateType const* yPar, StateType const* yBulk, double const* yDotPar, ResidualType* resPar, ResidualType* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc);
	};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_GENERALRATEPARTICLE_HPP_
