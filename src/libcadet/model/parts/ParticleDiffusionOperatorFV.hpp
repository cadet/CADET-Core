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

#ifndef LIBCADET_PARTCICLEDIFFUSIONOPERATORFV_HPP_
#define LIBCADET_PARTCICLEDIFFUSIONOPERATORFV_HPP_

#include "model/parts/ParticleDiffusionOperatorBase.hpp"
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

namespace parts
{
	namespace cell
	{
		struct CellParameters;
	}

	/**
	 * @brief Particle dispersion transport operator based on a FVM
	 * @details Implements the equation
	 *
	 @f[ \begin{align}
	 - \left. \left( \varepsilon^{\mathrm{p}} D^{\mathrm{p}}_{i} \frac{\partial c^{\mathrm{p}}_{i}}{\partial r} + (1 - \varepsilon^{\mathrm{p}}) D^{\mathrm{s}}_{i} \frac{\partial c^{\mathrm{s}}_{i}}{\partial r} \right) \right|_{r=0}
	 &= 0, \\
	 \varepsilon^{\mathrm{p}} \left. \left( \varepsilon^{\mathrm{p}}  D^{\mathrm{p}}_{i} \frac{\partial c^{\mathrm{p}}_{i}}{\partial r} + (1 - \varepsilon^{\mathrm{p}} ) D^{\mathrm{s}}_{i} \frac{\partial c^{\mathrm{s}}_{i}}{\partial r} \right)\right|_{r = R^{\mathrm{p}}_{}}
	 &= k^{\mathrm{f}}_{i} \left. \left( c^{\mathrm{b}}_i - c^{\mathrm{p}}_{i} \right|_{r = R^{\mathrm{p}}_{}} \right)
	 \end{align} @f]
	 * Additionally implements the variants for cylindrical and slab-shaped particles
	 * Methods are described in @cite Lieres2010
	 *
	 * This class does not store the Jacobian. It only fills existing matrices given to its residual() functions.
	 * It assumes that there is no offset to the particle entries in the local state vector
	 */
	class ParticleDiffusionOperatorFV : public ParticleDiffusionOperatorBase
	{
	public:

		ParticleDiffusionOperatorFV();
		~ParticleDiffusionOperatorFV() CADET_NOEXCEPT;

		bool configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp);
		bool configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound, const int* reqBinding);

		void setEquidistantRadialDisc();
		void setEquivolumeRadialDisc();
		void setUserdefinedRadialDisc();
		void updateRadialDisc();

		void clearParDepSurfDiffusion();

		bool notifyDiscontinuousSectionTransition(double t, unsigned int secIdx);

		int residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, linalg::BandedEigenSparseRowIterator& jacIt, WithoutParamSensitivity);
		int residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, linalg::BandedEigenSparseRowIterator& jacIt, WithParamSensitivity);
		int residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, linalg::BandedEigenSparseRowIterator& jacIt, WithoutParamSensitivity);
		int residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, linalg::BandedEigenSparseRowIterator& jacIt, WithParamSensitivity);

		std::vector<active> _parOuterSurfAreaPerVolume; //!< Particle element outer sphere surface to volume ratio
		std::vector<active> _parInnerSurfAreaPerVolume; //!< Particle element inner sphere surface to volume ratio

		/* Model discretization */

		enum class ParticleDiscretizationMode : int
		{
			/**
			 * Equidistant distribution of element edges
			 */
			Equidistant,

			/**
			 * Volumes of elements are uniform
			 */
			Equivolume,

			/**
			 * element edges specified by user
			 */
			UserDefined
		};

		double relativeCoordinate(const unsigned int cellIdx) const CADET_NOEXCEPT
		{
			return static_cast<double>(std::accumulate(&_deltaR[0], &_deltaR[cellIdx], _deltaR[cellIdx] * 0.5) / (_parRadius - _parCoreRadius));
		}

		int calcFilmDiffJacobian(unsigned int secIdx, const int offsetCp, const int offsetC, const int nBulkPoints, const int nParType, const double colPorosity, const active* const parTypeVolFrac, Eigen::SparseMatrix<double, RowMajor>& globalJac, bool outliersOnly = false);

		int writeParticleCoordinates(double* coords) const;

		typedef Eigen::Triplet<double> T;

		int particleJacobianBandwidth(unsigned int& lowerBandwidth, unsigned int& upperBandWidth) const;

		void setParticleJacobianPattern(std::vector<T>& tripletList, unsigned int offsetPar, unsigned int offsetBulk, unsigned int colNode, unsigned int secIdx);

		unsigned int jacobianNNZperParticle() const;

		int calcParticleDiffJacobian(const int secIdx, const int colNode, const int offsetLocalCp, Eigen::SparseMatrix<double, RowMajor>& globalJac);

		bool setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue);

	protected:

		void parBindingPattern(std::vector<Eigen::Triplet<double>>& tripletList, const int offset, const unsigned int colNode);

		template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
		int residualImpl(double t, unsigned int secIdx, StateType const* yPar, StateType const* yBulk, double const* yDotPar, ResidualType* resPar, linalg::BandedEigenSparseRowIterator& jacBase);

		ParticleDiscretizationMode _parDiscMode; //!< Particle discretization mode

		std::vector<double> _parDiscVector; //!< Particle discretization element boundary coodinates

		std::vector<active> _parCenterRadius; //!< Particle node-centered position for each particle node

		/* FV specific operators */

		std::vector<active> _deltaR; //!< particle cell spacing
		ArrayPool _discParFlux; //!< Storage for discretized @f$ k_f @f$ value
		int _boundaryOrderFV; //!< Order of the bulk-particle boundary discretization
	};

} // namespace parts
} // namespace model
} // namespace cadet

#endif  // LIBCADET_PARTCICLEDIFFUSIONOPERATORFV_HPP_
