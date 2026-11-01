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

/**
 * @file
 * Defines the liquid equilibrium particle model as e.g. used in the stirred tank
 */

#ifndef LIBCADET_LIQUIDEQUILIBRIUMPARTICLE_HPP_
#define LIBCADET_LIQUIDEQUILIBRIUMPARTICLE_HPP_

#include "model/particle/ParticleModel.hpp"
#include "model/BindingModel.hpp"
#include "cadet/StrongTypes.hpp"
#include "ParamIdUtil.hpp"
#include "AutoDiff.hpp"
#include "Memory.hpp"
#include "model/ModelUtils.hpp"
#include "SimulationTypes.hpp"
#include "ParamReaderHelper.hpp"
#include "linalg/BandedEigenSparseRowIterator.hpp"
#include "model/ParameterMultiplexing.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <unordered_map>
#include <unordered_set>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <memory>
#include <vector>

namespace cadet
{

	class IConfigHelper;

namespace model
{

class IDynamicReactionModel;

namespace parts
{
	namespace cell
	{
		struct CellParameters;
	}
}

	/**
	 * @brief Liquid equilibrium particle model
	 * @details The particle is in rapid equilibrium with the liquid phase it is suspended in, that is,
	 *          it has neither pores nor a film diffusion resistance. It therefore contributes bound
	 *          states only, and its binding model is evaluated against that liquid phase directly:
	 @f[ \begin{align}
	 \frac{\partial q_i}{\partial t} &= f_{\text{iso}}(c, q)
	 \end{align} @f]
	 *          The whole particle volume is solid, so the porosity is zero and the bound state
	 *          concentrations refer to the particle volume.
	 *
	 *          The liquid phase it binds to is a state of the hosting unit operation rather than of
	 *          the particle, so the amount of substance that the bound states take up has to appear
	 *          in that unit operation's liquid balance. This particle adds it to @c resBulk itself,
	 *          weighted exactly like the film transfer term of the particles that do have a film
	 *          resistance, see @c bulkFilmTransfer.
	 */
	class LiquidEquilibriumParticle : public IParticleModel
	{
	public:

		LiquidEquilibriumParticle();
		~LiquidEquilibriumParticle() CADET_NOEXCEPT;

		static inline const char* identifier() CADET_NOEXCEPT { return "LIQUID_EQUILIBRIUM_PARTICLE"; }

		bool configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp) override;
		bool configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound, bool uoImplementsReqReactions) override;

		void updateRadialDisc() override { }

		parts::cell::CellParameters makeCellResidualParams(int const* qsReaction, unsigned int const* nBound) const override;

		bool notifyDiscontinuousSectionTransition(double t, unsigned int secIdx) override;

		int residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, double* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity) override;
		int residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity) override;
		int residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity) override;
		int residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity) override;

		double relativeCoordinate(const unsigned int nodeIdx) const CADET_NOEXCEPT override { return 0.5; }

		/**
		 * @brief Surface to volume ratio
		 * @details Transport into the particle is not resolved, so no surface area enters the equations.
		 */
		active surfaceToVolumeRatio() const CADET_NOEXCEPT override { return _zero; }

		/**
		 * @brief Particle porosity
		 * @details The particle has no pores, so its volume is entirely solid and the bound state
		 *          concentrations refer to the particle volume.
		 */
		inline const active& getPorosity() const CADET_NOEXCEPT override { return _zero; }
		inline const active* getPoreAccessFactor() const CADET_NOEXCEPT override { return nullptr; }
		inline const active* getFilmDiffusion(const unsigned int secIdx) const CADET_NOEXCEPT override { return nullptr; }

		inline int nDiscPoints() const CADET_NOEXCEPT override { return 1; }

		/**
		 * @brief Stride of the state block of this particle
		 * @details There is no particle liquid phase, so the block holds the bound states only. This is
		 *          what a hosting unit operation has to honour when it lays out its state vector, and it
		 *          is why a particle type of this kind cannot simply be slotted into a model that assumes
		 *          every type to own a liquid phase of NCOMP entries.
		 */
		inline int strideParBlock() const CADET_NOEXCEPT override { return static_cast<int>(_strideBound); }

		int writeParticleCoordinates(double* coords) const override;

		void bulkFilmTransferCoefficients(unsigned int secIdx, active const* const pointVelocity, const int nBulkPoints, active* filmDiffusion) const override { }

		void setParJacPattern(std::vector<Eigen::Triplet<double>>& tripletList, const unsigned int offsetPar, const unsigned int offsetBulk, unsigned int colNode, unsigned int secIdx, const bulkFilmTransfer& bulkFilm) const override;

		unsigned int jacobianNNZperParticle() const override;
		int calcParticleDiffJacobian(const int secIdx, const int colNode, const int offsetLocalCp, Eigen::SparseMatrix<double, Eigen::RowMajor>& globalJac) override { return 1; }
		int calcFilmDiffJacobian(unsigned int secIdx, const int offsetCp, const int offsetC, const int nBulkPoints, const int nParType, const bulkFilmTransfer& bulkFilm, const active* const parTypeVolFrac, const active* const pointVelocity, Eigen::SparseMatrix<double, Eigen::RowMajor>& globalJac, bool crossDepsOnly = false) override { return 1; }

		bool setParameter(const ParameterId& pId, double value) override;
		bool setParameter(const ParameterId& pId, int value) override;
		bool setParameter(const ParameterId& pId, bool value) override;
		bool setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue) override;
		bool setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& pId, double value) override;

		bool hasParameter(const ParameterId& pId) const override;
		double getParameterDouble(const ParameterId& pId) const override;
		std::unordered_map<ParameterId, double> getAllParameterValues(std::unordered_map<ParameterId, double>& data) const override;

		bool leanConsistentInitialStateValidity() const override { return true; }
		bool leanConsistentInitialTimeDerivativeValidity() const override { return true; }

		bool isParticleLumped() const CADET_NOEXCEPT override { return true; }

		/**
		 * @brief Sets the distance from the bound states of this particle type back to the liquid phase
		 * @details The bound states are in rapid equilibrium with a liquid phase that is not adjacent to
		 *          them in the state vector whenever the hosting unit operation holds more than one
		 *          particle type. The distance is the same for the state vector and for the rows of the
		 *          Jacobian, and it follows the convention of IBindingModel::analyticJacobian, whose
		 *          offsetCp counts backwards. Defaults to the layout of a column node, where the liquid
		 *          phase directly precedes the bound states.
		 */
		void setBulkRowOffset(int offset) CADET_NOEXCEPT override { _bulkRowOffset = offset; }
		inline int bulkRowOffset() const CADET_NOEXCEPT { return _bulkRowOffset; }

	protected:

		active _zero; //!< Holds 0.0, returned where the interface hands out a reference to a quantity this particle does not have

		/* strides and offsets */
		int _strideBulkComp; //!< Component stride in bulk state vector
		int _bulkRowOffset; //!< Distance from the bound states back to the liquid phase, see setBulkRowOffset()
		unsigned int* _boundOffset; //!< Array with offset to the first bound state of each component in the solid phase
		unsigned int _strideBound; //!< Total number of bound states

		/* residual implementation */
		template <typename StateType, typename ResidualType, typename ParamType, bool wantNonLinJac, bool wantRes>
		int residualImpl(double t, unsigned int secIdx, StateType const* yPar, StateType const* yBulk, double const* yDotPar, ResidualType* resPar, ResidualType* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc);
	};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_LIQUIDEQUILIBRIUMPARTICLE_HPP_
