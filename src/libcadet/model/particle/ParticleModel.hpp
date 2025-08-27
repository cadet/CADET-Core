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
 * Defines the particle model interface
 */

#ifndef LIBCADET_IPARTICLEMODEL_HPP_
#define LIBCADET_IPARTICLEMODEL_HPP_

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

#include <bitset>
#include <unordered_map>
#include <unordered_set>
#include <Eigen/Sparse>
#include <array>
#include <vector>

namespace cadet
{

	class IConfigHelper;

	namespace model
	{

		class IDynamicReactionModel;
		class IParameterStateDependence;

		class ParticleModel {
		private:
			// Bit positions for each particle transport aspect
			enum class TransportBit : size_t {
				FilmDiffusion = 0,
				PoreDiffusion = 1,
				SurfaceDiffusion = 2
			};

			std::bitset<3> _transportFlags;

			void validateConfiguration() const
			{
				if (!hasFilmDiffusion() && (hasPoreDiffusion() || hasSurfaceDiffusion()))
					throw std::invalid_argument("Invalid Particle configuration: instantaneous/no film diffusion is not supported in combination with pore diffusion and surface diffusion");
			}

		public:

			ParticleModel(bool filmDiffusion, bool poreDiffusion, bool surfaceDiffusion) {
				_transportFlags[static_cast<size_t>(TransportBit::FilmDiffusion)] = filmDiffusion;
				_transportFlags[static_cast<size_t>(TransportBit::PoreDiffusion)] = poreDiffusion;
				_transportFlags[static_cast<size_t>(TransportBit::SurfaceDiffusion)] = surfaceDiffusion;

				validateConfiguration();
			}

			bool hasFilmDiffusion() const {
				return _transportFlags[static_cast<size_t>(TransportBit::FilmDiffusion)];
			}

			bool hasPoreDiffusion() const {
				return _transportFlags[static_cast<size_t>(TransportBit::PoreDiffusion)];
			}

			bool hasSurfaceDiffusion() const {
				return _transportFlags[static_cast<size_t>(TransportBit::SurfaceDiffusion)];
			}
			
			std::string getParticleTransportType() const
			{
				validateConfiguration();

				if (!hasFilmDiffusion())
					return "EQUILIBRIUM_PARTICLE";
				else if (!hasPoreDiffusion() && !hasSurfaceDiffusion())
					return "HOMOGENEOUS_PARTICLE";
				else if (hasPoreDiffusion() || hasSurfaceDiffusion())
					return "GENERAL_RATE_PARTICLE";
				else
					return "UNKNOWN";
			}
		};

		namespace parts
		{
			namespace cell
			{
				struct CellParameters;
			}
		}

		struct columnPackingParameters
		{
			const active& parTypeVolFrac;
			const active& colPorosity;
			ColumnPosition colPos;
		};

		/**
		 * @brief Implements the particle model interface
		 */
		class IParticleModel
		{
		public:

			IParticleModel() : _binding(nullptr), _dynReaction(nullptr)
			{
			}
			virtual ~IParticleModel() CADET_NOEXCEPT
			{
				delete _binding;
				delete _dynReaction;
			}

			virtual bool configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp) = 0;
			virtual bool configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound) = 0;

			virtual void updateRadialDisc() = 0;

			virtual parts::cell::CellParameters makeCellResidualParams(int const* qsReaction, unsigned int const* nBound) const = 0;
			/**
			 * @brief Notifies the operator that a discontinuous section transition is in progress
			 * @param [in] t Current time point
			 * @param [in] secIdx Index of the new section that is about to be integrated
			 * @return @c true if flow direction has changed, otherwise @c false
			 */
			virtual bool notifyDiscontinuousSectionTransition(double t, unsigned int secIdx) = 0;

			/**
			 * @brief Computes the residual of the particle equations and updates nonlinear Jacobian entries
			 * @param [in] model Model that owns the operator
			 * @param [in] t Current time point
			 * @param [in] secIdx Index of the current section
			 * @param [in] yPar Pointer to particle phase entry in unit state vector
			 * @param [in] yBulk Pointer to corresponding bulk phase entry in unit state vector
			 * @param [in] yDotPar Pointer to particle phase derivative entry in unit state vector
			 * @param [out] resPar Pointer Pointer to particle phase entry in unit residual vector, nullptr if no residual shall be computed
			 * @param [out] packing column packing parameters
			 * @param [in] jacIt Row iterator pointing to the particle phase entry in the unit Jacobian, uninitialized if no Jacobian shall be computed, otherwise only non-linear Jacobian entries will be added
			 * @param [in] tlmAlloc memory allocator
			 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
			 */
			virtual int residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, double* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity) = 0;
			virtual int residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity) = 0;
			virtual int residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity) = 0;
			virtual int residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity) = 0;

			virtual double relativeCoordinate(const unsigned int nodeIdx) const CADET_NOEXCEPT = 0;

			virtual active surfaceToVolumeRatio() const CADET_NOEXCEPT = 0;

			inline IBindingModel* getBinding() const CADET_NOEXCEPT { return _binding; }
			inline bool bindingParDep() const CADET_NOEXCEPT { return _bindingParDep; }
			inline IDynamicReactionModel* getReaction() const CADET_NOEXCEPT { return _dynReaction; }
			inline bool reactionParDep() const CADET_NOEXCEPT { return _reactionParDep; }

			virtual inline const active& getPorosity() const CADET_NOEXCEPT = 0;
			virtual inline const active* getPoreAccessFactor() const CADET_NOEXCEPT = 0;
			virtual inline const active* getFilmDiffusion(const unsigned int secIdx) const CADET_NOEXCEPT = 0;

			/**
			 * @brief total number discrete points per particle
			 */
			virtual inline int nDiscPoints() const CADET_NOEXCEPT = 0;
			/**
			 * @brief stride over one particle
			 */
			virtual inline int strideParBlock() const CADET_NOEXCEPT = 0;
			/**
			 * @brief computes the physical particle coordinates for every spatial point
			 */
			virtual int writeParticleCoordinates(double* coords) const = 0;

			/**
			 * @brief sets the particle sparsity pattern wrt the global Jacobian
			 */
			virtual void setParJacPattern(std::vector<Eigen::Triplet<double>>& tripletList, const unsigned int offsetPar, const unsigned int offsetBulk, unsigned int colNode, unsigned int secIdx) const = 0;

			virtual unsigned int jacobianNNZperParticle() const = 0;
			/**
			 * @brief analytically calculates the static (per section) particle diffusion Jacobian
			 * @return 1 if jacobain calculation fits the predefined pattern of the jacobian, 0 if not.
			 */
			virtual int calcParticleDiffJacobian(const int secIdx, const int colNode, const int offsetLocalCp, Eigen::SparseMatrix<double, Eigen::RowMajor>& globalJac) = 0;
			/**
			 * @brief calculates the film diffusion Jacobian
			 * @return 1 if jacobain calculation fits the predefined pattern of the jacobian, 0 if not
			 * @param [in] parTypeVolFrac pointer to the particle type volume fractions in column-position-major order
			 */
			virtual int calcFilmDiffJacobian(unsigned int secIdx, const int offsetCp, const int offsetC, const int nBulkPoints, const int nParType, const double colPorosity, const active* const parTypeVolFrac, Eigen::SparseMatrix<double, Eigen::RowMajor>& globalJac, bool crossDepsOnly = false) = 0;

			virtual bool setParameter(const ParameterId& pId, double value) = 0;
			virtual bool setParameter(const ParameterId& pId, int value) = 0;
			virtual bool setParameter(const ParameterId& pId, bool value) = 0;
			virtual bool setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue) = 0;
			virtual bool setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& pId, double value) = 0;

			virtual bool hasParameter(const ParameterId& pId) const = 0;
			virtual double getParameterDouble(const ParameterId& pId) const = 0;
			virtual std::unordered_map<ParameterId, double> getAllParameterValues(std::unordered_map<ParameterId, double>& data) const = 0;

			virtual bool leanConsistentInitialStateValidity() const = 0;

			virtual bool leanConsistentInitialTimeDerivativeValidity() const = 0;

			unsigned int* nBound() CADET_NOEXCEPT { return _nBound.get(); }

			virtual bool isParticleLumped() const CADET_NOEXCEPT = 0;

			protected:

				unsigned int _parTypeIdx; //!< Particle type index (wrt the unit operation that owns this particle model)
				unsigned int _nComp; //!< Number of components

				IBindingModel* _binding; //!< Binding model
				std::shared_ptr<unsigned int[]> _nBound; //!< Array with number of bound states for each component
				bool _bindingParDep; //!< Whether the binding model parameters depend on the particle type
				IDynamicReactionModel* _dynReaction; //!< Dynamic reaction model
				bool _reactionParDep; //!< Whether the binding model parameters depend on the particle type
		};

	} // namespace model
} // namespace cadet

#endif  // LIBCADET_IPARTICLEMODEL_HPP_
