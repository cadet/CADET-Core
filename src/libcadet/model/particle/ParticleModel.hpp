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

#include <unordered_map>
#include <unordered_set>
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
			virtual bool notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active const* const filmDiff, active const* const poreAccessFactor) = 0;

			/**
			 * @brief Computes the residual of the particle equations and updates nonlinear Jacobian entries
			 * @param [in] model Model that owns the operator
			 * @param [in] t Current time point
			 * @param [in] secIdx Index of the current section
			 * @param [in] yPar Pointer to particle phase entry in unit state vector
			 * @param [in] yBulk Pointer to corresponding bulk phase entry in unit state vector
			 * @param [in] yDotPar Pointer to particle phase derivative entry in unit state vector
			 * @param [out] resPar Pointer Pointer to particle phase entry in unit residual vector, nullptr if no residual shall be computed
			 * @param [out] colPos column position of the particle (particle coordinate zero)
			 * @param [in] jacIt Row iterator pointing to the particle phase entry in the unit Jacobian, uninitialized if no Jacobian shall be computed, otherwise only non-linear Jacobian entries will be added
			 * @param [in] tlmAlloc memory allocator
			 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
			 */
			virtual int residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity) = 0;
			virtual int residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity) = 0;
			virtual int residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity) = 0;
			virtual int residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity) = 0;

			unsigned int _parTypeIdx; //!< Particle type index (wrt the unit operation that owns this particle model)

			unsigned int _nComp; //!< Number of components

			IBindingModel* _binding; //!< Binding model
			bool _singleBinding; //!< Determines whether only a single binding model is present in the whole unit
			IDynamicReactionModel* _dynReaction; //!< Dynamic reaction model
			bool _singleDynReaction; //!< Determines whether only a single particle reaction model is present in the whole unit

			virtual double relativeCoordinate(const unsigned int nodeIdx) const CADET_NOEXCEPT = 0;

			virtual active surfaceToVolumeRatio() const CADET_NOEXCEPT = 0;

			inline IBindingModel* getBinding() const CADET_NOEXCEPT { return _binding; }
			inline bool singleBinding() const CADET_NOEXCEPT { return _singleBinding; }
			inline IDynamicReactionModel* getReaction() const CADET_NOEXCEPT { return _dynReaction; }
			inline bool singleReaction() const CADET_NOEXCEPT { return _singleDynReaction; }

			virtual inline const active& getPorosity() const CADET_NOEXCEPT = 0;
			virtual inline const active* getPoreAccessfactor() const CADET_NOEXCEPT = 0;

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
			virtual int calcStaticAnaParticleDiffJacobian(const int secIdx, const int colNode, const int offsetLocalCp, Eigen::SparseMatrix<double, RowMajor>& globalJac) = 0;
			virtual int calcFilmDiffJacobian(unsigned int secIdx, const int offsetCp, const int offsetC, const int nBulkPoints, const int nParType, const double colPorosity, const active* const parTypeVolFrac, Eigen::SparseMatrix<double, RowMajor>& globalJac, bool crossDepsOnly = false) = 0;

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

		};

	} // namespace model
} // namespace cadet

#endif  // LIBCADET_IPARTICLEMODEL_HPP_
