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
 * Defines the particle dispersion transport operator according to the discontinuous Galerkin discretization.
 */

#ifndef LIBCADET_GENERALRATEPARTICLE_HPP_
#define LIBCADET_GENERALRATEPARTICLE_HPP_

#include "model/parts/ParticleDiffusionOperatorDG.hpp"
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
	//constexpr double _SurfVolRatioSphere = 3.0; //!< Surface to volume ratio for a spherical particle
	//constexpr double _SurfVolRatioCylinder = 2.0; //!< Surface to volume ratio for a cylindrical particle
	//constexpr double _SurfVolRatioSlab = 1.0; //!< Surface to volume ratio for a slab-shaped particle

	/**
	 * @brief Particle dispersion transport operator
	 * @details Implements the equation
	 *
	 * @f[\begin{align}
	 \frac{\partial c^{\mathrm{p}}_{i}}{\partial t} + \frac{1 - \varepsilon^{\mathrm{p}}}{\varepsilon^{\mathrm{p}}} \frac{\partial c^{\mathrm{s}}_{i}}{\partial t} &=
	 \frac{1}{r^2} \frac{\partial }{\partial r} \left( r^2 D_{i}^{\mathrm{p}} \frac{\partial c^{\mathrm{p}}_{i}}{\partial r} \right)
	 - \frac{1 - \varepsilon^{\mathrm{p}}}{\varepsilon^{\mathrm{p}}} \frac{1}{r^2} \frac{\partial }{\partial r} \left( r^2 D_{i}^{\mathrm{s}} \frac{\partial c^{\mathrm{s}}_{i}}{\partial r} \right) \\
	 \end{align} @f]
	 * with Danckwerts boundary conditions (see @cite Danckwerts1953)
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
	class GeneralRateParticle
	{
	public:

		GeneralRateParticle();
		~GeneralRateParticle() CADET_NOEXCEPT;

		bool configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp);
		bool configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound);

		void updateRadialDisc() { _parDiffOp->updateRadialDisc(); }

		parts::cell::CellParameters makeCellResidualParams(int const* qsReaction, unsigned int const* nBound) const;

		//void clearParDepSurfDiffusion();

		bool notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active const* const filmDiff, active const* const poreAccessFactor);

		template<bool wantJac, bool wantRes>
		int residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity);
		template<bool wantJac, bool wantRes>
		int residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity);
		template<bool wantJac, bool wantRes>
		int residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity);
		template<bool wantJac, bool wantRes>
		int residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity);

		unsigned int _parTypeIdx; //!< Particle type index (wrt the unit operation that owns this particle model)

		unsigned int _nComp; //!< Number of components

		IBindingModel* _binding; //!< Binding model
		bool _singleBinding; //!< Determines whether only a single binding model is present in the whole unit
		IDynamicReactionModel* _dynReaction; //!< Dynamic reaction model
		bool _singleDynReaction; //!< Determines whether only a single particle reaction model is present in the whole unit

		double relativeCoordinate(const unsigned int nodeIdx) const CADET_NOEXCEPT { return _parDiffOp->relativeCoordinate(nodeIdx); }

		template<typename ParamType>
		ParamType surfaceToVolumeRatio() const CADET_NOEXCEPT
		{
			return _parDiffOp->surfaceToVolumeRatio<ParamType>();
		}

		parts::ParticleDiffusionOperatorDG* _parDiffOp;

		unsigned int* nBound() CADET_NOEXCEPT { return _parDiffOp->_nBound; }; //!< Array with number of bound states for each component
		inline IBindingModel* getBinding() const CADET_NOEXCEPT { return _binding; }
		inline bool singleBinding() const CADET_NOEXCEPT { return _singleBinding; }
		inline IDynamicReactionModel* getReaction() const CADET_NOEXCEPT { return _dynReaction; }
		inline bool singleReaction() const CADET_NOEXCEPT { return _singleDynReaction; }

		inline active& getPorosity() const CADET_NOEXCEPT { return _parDiffOp->_parPorosity; }
		inline active* getPoreAccessfactor() const CADET_NOEXCEPT { return &_parDiffOp->_poreAccessFactor[0]; }
		inline IParameterStateDependence* getParDepSurfDiffusion() const CADET_NOEXCEPT { return _parDiffOp->_parDepSurfDiffusion; }
		inline bool singleParDepSurfDiffusion() const CADET_NOEXCEPT { return _parDiffOp->_singleParDepSurfDiffusion; }
		inline MultiplexMode parDiffMode() const CADET_NOEXCEPT { return _parDiffOp->_parDiffusionMode; }
		inline MultiplexMode parSurfDiffMode() const CADET_NOEXCEPT { return _parDiffOp->_parSurfDiffusionMode; }

		/**
		 * @brief array with offsets to the first bound state of each component in the solid phase
		 */
		inline unsigned int* offsetBoundComp() const CADET_NOEXCEPT { return _parDiffOp->_boundOffset; };
		/**
		 * @brief offset to the first bound state
		 */
		inline unsigned int offsetBoundComp(ComponentIndex comp) const CADET_NOEXCEPT { return offsetBoundComp()[comp.value]; }
		/**
		 * @brief total number of bound states
		 */
		inline unsigned int strideBound() const CADET_NOEXCEPT { return _parDiffOp->_strideBound; };
		/**
		 * @brief total number discrete points per particle
		 */
		inline int nDiscPoints() const CADET_NOEXCEPT { return _parDiffOp->_nParPoints; }
		/**
		 * @brief stride over all components in the liquid phase
		 */
		inline int strideLiquid() const CADET_NOEXCEPT { return static_cast<int>(_nComp); }
		/**
		 * @brief stride over one discrete point
		 */
		inline int stridePoint() const CADET_NOEXCEPT { return strideLiquid() + strideBound(); }
		/**
		 * @brief stride over one particle
		 */
		inline int strideParBlock() const CADET_NOEXCEPT { return nDiscPoints() * stridePoint(); }

		int getParticleCoordinates(double* coords) const;

		typedef Eigen::Triplet<double> T;

		/**
		 * @brief sets the particle sparsity pattern wrt the global Jacobian
		 */
		void setParJacPattern(std::vector<T>& tripletList, const unsigned int offsetPar, const unsigned int offsetBulk, unsigned int colNode, unsigned int secIdx)
		{
			_parDiffOp->calcParticleJacobianPattern(tripletList, offsetPar, offsetBulk, colNode, secIdx);

			_parDiffOp->parTimeDerJacPattern_GRM(tripletList, offsetPar, colNode, secIdx);

			_parDiffOp->parBindingPattern_GRM(tripletList, offsetPar, colNode);
		}

		unsigned int calcParDiffNNZ();
		int calcStaticAnaParticleDiffJacobian(const int secIdx, const int colNode, const int offsetLocalCp, Eigen::SparseMatrix<double, RowMajor>& globalJac);
		int calcFilmDiffJacobian(unsigned int secIdx, const int offsetCp, const int offsetC, const int nBulkPoints, const int nParType, const double colPorosity, const active* const parTypeVolFrac, Eigen::SparseMatrix<double, RowMajor>& globalJac, bool outliersOnly = false);

		bool setParameter(const ParameterId& pId, double value);
		bool setParameter(const ParameterId& pId, int value);
		bool setParameter(const ParameterId& pId, bool value);
		bool setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue);
		bool setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& pId, double value);

	protected:

		template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
		int residualImpl(double t, unsigned int secIdx, StateType const* yPar, StateType const* yBulk, double const* yDotPar, ResidualType* resPar, ColumnPosition colPos, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc);
	};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_GENERALRATEPARTICLE_HPP_
