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
 * Defines the homogeneous particle model as e.g. used in the Lumped Rate Model with Pores
 */

#ifndef LIBCADET_HOMOGENEOUSPARTICLE_HPP_
#define LIBCADET_HOMOGENEOUSPARTICLE_HPP_

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
	 * @brief Homogeneous Particle Model
	 @f[ \begin{align}
	 \frac{\partial c_i}{\partial t} &= - u \frac{\partial c_i}{\partial z} + D_{\text{ax},i} \frac{\partial^2 c_i}{\partial z^2} - \frac{1 - \varepsilon_c}{\varepsilon_c} \frac{3 k_{f,i}}{r_p} j_{f,i} \\
	 \frac{\partial c_{p,i}}{\partial t} + \frac{1 - \varepsilon_p}{\varepsilon_p} \frac{\partial q_{i}}{\partial t} &= \frac{3 k_{f,i}}{\varepsilon_p r_p} (c_i - c_{p,i}) \\
	 a \frac{\partial q_i}{\partial t} &= f_{\text{iso}}(c_p, q)
	 \end{align} @f]
	 * Additionally implements the variants for cylindrical and slab-shaped particles
	 */
	class HomogeneousParticle : public IParticleModel
	{
	public:

		HomogeneousParticle();
		~HomogeneousParticle() CADET_NOEXCEPT;

		static inline const char* identifier() CADET_NOEXCEPT { return "HOMOGENEOUS_PARTICLE"; }

		bool configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp) override;
		bool configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound) override;
		bool configureModelDiscretization_old(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp) override;
		bool configure_old(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound) override;

		void updateRadialDisc() { }

		parts::cell::CellParameters makeCellResidualParams(int const* qsReaction, unsigned int const* nBound) const override;

		bool notifyDiscontinuousSectionTransition(double t, unsigned int secIdx) override;

		int residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, double* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity) override;
		int residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity) override;
		int residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithoutParamSensitivity) override;
		int residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, active* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc, WithParamSensitivity) override;

		double relativeCoordinate(const unsigned int nodeIdx) const CADET_NOEXCEPT  override { return 0.5; }

		active surfaceToVolumeRatio() const CADET_NOEXCEPT override
		{
			return _parGeomSurfToVol / _parRadius;
		}

		inline const active& getPorosity() const CADET_NOEXCEPT  override { return _parPorosity; }
		inline const active* getPoreAccessfactor() const CADET_NOEXCEPT  override { return &_poreAccessFactor[0]; }
		inline const active* getFilmDiffusion(const unsigned int secIdx) const CADET_NOEXCEPT { return getSectionDependentSlice(_filmDiffusion, _nComp, secIdx); }

		inline int nDiscPoints() const CADET_NOEXCEPT  override { return 1; }
		inline int strideParBlock() const CADET_NOEXCEPT  override { return nDiscPoints() * stridePoint(); }

		int writeParticleCoordinates(double* coords) const override;

		void setParJacPattern(std::vector<Eigen::Triplet<double>>& tripletList, const unsigned int offsetPar, const unsigned int offsetBulk, unsigned int colNode, unsigned int secIdx) const override;

		unsigned int jacobianNNZperParticle() const override;
		int calcStaticAnaParticleDiffJacobian(const int secIdx, const int colNode, const int offsetLocalCp, Eigen::SparseMatrix<double, RowMajor>& globalJac) override;
		int calcFilmDiffJacobian(unsigned int secIdx, const int offsetCp, const int offsetC, const int nBulkPoints, const int nParType, const double colPorosity, const active* const parTypeVolFrac, Eigen::SparseMatrix<double, RowMajor>& globalJac, bool crossDepsOnly = false) override;

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

	protected:

		/* diffusion */
		std::vector<active> _filmDiffusion; //!< Particle diffusion coefficient \f$ D_p \f$
		MultiplexMode _filmDiffusionMode;
		std::vector<active> _poreAccessFactor; //!< Pore accessibility factor \f$ F_{\text{acc}} \f$
		MultiplexMode _poreAccessFactorMode;

		/* geometry */
		double _SurfVolRatioSphere = 3.0; //!< Surface to volume ratio for a spherical particle
		double _SurfVolRatioCylinder = 2.0; //!< Surface to volume ratio for a cylindrical particle
		double _SurfVolRatioSlab = 1.0; //!< Surface to volume ratio for a slab-shaped particle
		double _parGeomSurfToVol; //!< Particle surface to volume ratio factor (i.e., 3.0 for spherical, 2.0 for cylindrical, 1.0 for hexahedral)
		active _parRadius; //!< Particle radius \f$ r_p \f$
		bool _parRadiusParTypeDep; //!< Determines whether or not the radius is particle type dependent / specific. Used in parameter sensitivities
		active _parPorosity; //!< Particle porosity (internal porosity) \f$ \varepsilon_p \f$
		bool _parPorosityParTypeDep; //!< Determines whether or not the porosity is particle type dependent / specific. Used in parameter sensitivities

		/* strides and offsets */
		int _strideBulkComp; //!< Component stride in bulk state vector
		unsigned int* _boundOffset; //!< Array with offset to the first bound state of each component in the solid phase
		unsigned int _strideBound; //!< Total number of bound states
		/**
		 * @brief stride over one discrete point
		 */
		inline int stridePoint() const CADET_NOEXCEPT { return static_cast<int>(_nComp) + _strideBound; }

		/* residual implementation */
		template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
		int residualImpl(double t, unsigned int secIdx, StateType const* yPar, StateType const* yBulk, double const* yDotPar, ResidualType* resPar, ResidualType* resBulk, columnPackingParameters packing, linalg::BandedEigenSparseRowIterator& jacIt, LinearBufferAllocator tlmAlloc);
	};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_HOMOGENEOUSPARTICLE_HPP_
