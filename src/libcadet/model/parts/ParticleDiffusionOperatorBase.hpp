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
 * Defines the particle diffusion operator interface.
 */

#ifndef LIBCADET_PARTCICLEDIFFUSIONOPERATORBASE_HPP_
#define LIBCADET_PARTCICLEDIFFUSIONOPERATORBASE_HPP_

#include "cadet/StrongTypes.hpp"
#include "ParamIdUtil.hpp"
#include "AutoDiff.hpp"
#include "Memory.hpp"
#include "model/ModelUtils.hpp"
#include "SimulationTypes.hpp"
#include "ParamReaderHelper.hpp"
#include "model/ParameterMultiplexing.hpp"
#include "linalg/BandedEigenSparseRowIterator.hpp"

#include <unordered_map>
#include <unordered_set>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>
#include <vector>

namespace cadet
{

class IConfigHelper;

namespace model
{

class IParameterStateDependence;

namespace parts
{
	constexpr double _SurfVolRatioSphere = 3.0; //!< Surface to volume ratio for a spherical particle
	constexpr double _SurfVolRatioCylinder = 2.0; //!< Surface to volume ratio for a cylindrical particle
	constexpr double _SurfVolRatioSlab = 1.0; //!< Surface to volume ratio for a slab-shaped particle

	class ParticleDiffusionOperatorBase
	{
	public:

		ParticleDiffusionOperatorBase();

		virtual ~ParticleDiffusionOperatorBase() CADET_NOEXCEPT;

		/**
		^* @brief Sets fixed parameters of the particle diffusion operator (e.g., the number of discretization points, components, bound states)
		 * @details This function is called prior to configure() by the underlying model.
		 *          It can only be called once. Whereas non-structural model parameters
		 *          (e.g., rate constants) are configured by configure(), this function
		 *          sets structural parameters (e.g., number of components and bound states).
		 *
		 * @param [in] paramProvider Parameter provider
		 * @param [in] nComp Number of components
		 * @param [in] nBound Array of size @p nComp which contains the number of bound states for each component
		 * @param [in] boundOffset Array of size @p nComp with offsets to the first bound state of each component beginning from the solid phase
		 */
		virtual bool configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper, const int nComp, const int parTypeIdx, const int nParType, const int strideBulkComp) = 0;

		/**
		 * @brief Configures the model by extracting all non-structural parameters (e.g., model parameters) from the given @p paramProvider
		 * @details The scope of the cadet::IParameterProvider is left unchanged on return.
		 *
		 *          The structure of the model is left unchanged, that is, the number of degrees of
		 *          freedom stays the same (e.g., number of bound states is left unchanged). Only
		 *          true (non-structural) model parameters are read and changed.
		 *
		 *          This function may only be called if configureModelDiscretization() has been called
		 *          in the past. Contrary to configureModelDiscretization(), it can be called multiple
		 *          times.
		 * @param [in] paramProvider Parameter provider
		 * @return @c true if the configuration was successful, otherwise @c false
		 */
		virtual bool configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, const int nParType, const unsigned int* nBoundBeforeType, const int nTotalBound, const int* reqBinding, const bool hasDynamicReactions) = 0;

		/**
		 * @brief updates radial discretization operators based on configured parameters
		 */
		virtual void updateRadialDisc() = 0;

		/**
		 * @brief Notifies the operator that a discontinuous section transition is in progress
		 * @param [in] t Current time point
		 * @param [in] secIdx Index of the new section that is about to be integrated
		 * @param [in] filmDiff pointer to film diffusion parameter of unit operation
		 * @param [in] poreAccessFactor pointer to pore access factor parameter of unit operation
		 * @return @c true if flow direction has changed, otherwise @c false
		 */
		virtual bool notifyDiscontinuousSectionTransition(double t, unsigned int secIdx) = 0;

		/**
		 * @brief Computes the residual of the transport equations
		 * @param [in] model Model that owns the operator
		 * @param [in] t Current time point
		 * @param [in] secIdx Index of the current section
		 * @param [in] yPar Pointer to particle phase entry in unit state vector
		 * @param [in] yBulk Pointer to corresponding bulk phase entry in unit state vector
		 * @param [in] yDotPar Pointer to particle phase derivative entry in unit state vector
		 * @param [out] resPar Pointer Pointer to particle phase entry in unit residual vector, nullptr if no residual shall be computed
		 * @param [in] jacIt Row iterator pointing to the particle phase entry in the unit Jacobian, uninitialized if no Jacobian shall be computed
		 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
		 */
		virtual int residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, double* resPar, linalg::BandedEigenSparseRowIterator& jacIt, WithoutParamSensitivity) = 0;
		virtual int residual(double t, unsigned int secIdx, double const* yPar, double const* yBulk, double const* yDotPar, active* resPar, linalg::BandedEigenSparseRowIterator& jacIt, WithParamSensitivity) = 0;
		virtual int residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, linalg::BandedEigenSparseRowIterator& jacIt, WithoutParamSensitivity) = 0;
		virtual int residual(double t, unsigned int secIdx, active const* yPar, active const* yBulk, double const* yDotPar, active* resPar, linalg::BandedEigenSparseRowIterator& jacIt, WithParamSensitivity) = 0;

		virtual int calcFilmDiffJacobian(unsigned int secIdx, const int offsetCp, const int offsetC, const int nBulkPoints, const int nParType, const double colPorosity, const active* const parTypeVolFrac, Eigen::SparseMatrix<double, Eigen::RowMajor>& globalJac, bool outliersOnly = false) = 0;
		virtual int calcParticleDiffJacobian(const int secIdx, const int colNode, const int offsetLocalCp, Eigen::SparseMatrix<double, Eigen::RowMajor>& globalJac) = 0;
		
		/**
		 *@brief adds the solid time derivative and binding pattern to the list
		 */
		virtual void setParticleJacobianPattern(std::vector<Eigen::Triplet<double>>& tripletList, unsigned int offsetPar, unsigned int offsetBulk, unsigned int colNode, unsigned int secIdx) = 0;

		/**
		 * @brief calculates and returns the physical particle coordinates according to the discretization
		 */
		virtual int writeParticleCoordinates(double* coords) const = 0;
		/**
		 * @brief calculates and returns the relative particle coordinate in [0, 1] for the given node index
		 */
		virtual double relativeCoordinate(const unsigned int nodeIdx) const = 0;

		active surfaceToVolumeRatio() const CADET_NOEXCEPT
		{
			return _parGeomSurfToVol / _parRadius;
		}

		/**
		 * @brief Number of non-zero Jacobian entries occuring from a single particle of this type
		 * @detail includes intra-particle entries and film diffusion entries
		 */
		virtual unsigned int jacobianNNZperParticle() const = 0;
		
		virtual bool setParameter(const ParameterId& pId, double value);
		virtual bool setParameter(const ParameterId& pId, int value);
		virtual bool setParameter(const ParameterId& pId, bool value);
		virtual bool setSensitiveParameter(std::unordered_set<active*>& sensParams, const ParameterId& pId, unsigned int adDirection, double adValue);
		virtual bool setSensitiveParameterValue(const std::unordered_set<active*>& sensParams, const ParameterId& pId, double value);
		
		/**
		 * @brief array with number of bound states for each component
		 */
		std::shared_ptr<unsigned int[]> nBound() CADET_NOEXCEPT { return _nBound; };
		void setNBound(std::shared_ptr<unsigned int[]> nBound) CADET_NOEXCEPT { _nBound = nBound; };
		/**
		 * @brief array with offsets to the first bound state of each component in the solid phase
		 */
		inline unsigned int* offsetBoundComp() const CADET_NOEXCEPT { return _boundOffset; };
		/**
		 * @brief offset to the first bound state
		 */
		inline unsigned int offsetBoundComp(ComponentIndex comp) const CADET_NOEXCEPT { return offsetBoundComp()[comp.value]; }
		/**
		 * @brief total number of bound states
		 */
		inline unsigned int strideBound() const CADET_NOEXCEPT { return _strideBound; };
		/**
		 * @brief total number discrete points per particle
		 */
		inline int nDiscPoints() const CADET_NOEXCEPT { return _nParPoints; }

		inline const active& getPorosity() const CADET_NOEXCEPT { return _parPorosity; }
		inline const active* getPoreAccessFactor() const CADET_NOEXCEPT { return &_poreAccessFactor[0]; }
		inline const active* getFilmDiffusion(const unsigned int secIdx) const CADET_NOEXCEPT { return getSectionDependentSlice(_filmDiffusion, _nComp, secIdx); }
		inline IParameterStateDependence* getParDepSurfDiffusion() const CADET_NOEXCEPT { return _parDepSurfDiffusion; }
		inline bool paramDepSurfDiffParTypeIndep() const CADET_NOEXCEPT { return !_paramDepSurfDiffTypeDep; }
		inline MultiplexMode parDiffMode() const CADET_NOEXCEPT { return _parDiffusionMode; }
		inline MultiplexMode parSurfDiffMode() const CADET_NOEXCEPT { return _parSurfDiffusionMode; }

	protected:

		/* component system */
		unsigned int _nComp; //!< Number of components
		std::shared_ptr<unsigned int[]> _nBound; //!< Array with number of bound states for each component

		/* geometry */
		double _parGeomSurfToVol; //!< Particle surface to volume ratio factor (i.e., 3.0 for spherical, 2.0 for cylindrical, 1.0 for hexahedral)
		active _parRadius; //!< Particle radius \f$ r_p \f$
		bool _parRadiusParTypeDep; //!< Determines whether radius is particle type dependent, needed for sensitivitites
		active _parCoreRadius; //!< Particle core radius \f$ r_c \f$
		bool _parCoreRadiusParTypeDep; //!< Determines whether core radius is particle type dependent, needed for sensitivitites
		active _parPorosity; //!< Particle porosity (internal porosity) \f$ \varepsilon_p \f$
		bool _parPorosityParTypeDep; //!< Determines whether porosity is particle type dependent, needed for sensitivitites
		std::vector<active> _poreAccessFactor; //!< Pore accessibility factor \f$ F_{\text{acc}} \f$
		MultiplexMode _poreAccessFactorMode; //!< Determines the multiplex of the pore access factor, needed for sensitivitites
		std::vector<active> _invBetaP; //!< Ratio of solid to liquid particle volume

		/* diffusion rates */
		std::vector<active> _filmDiffusion; //!< Film diffusion coefficient \f$ k_f \f$
		MultiplexMode _filmDiffusionMode; //!< Determines the multiplex of film diffusion, needed for sensitivitites
		std::vector<active> _parDiffusion; //!< Particle diffusion coefficient \f$ D_p \f$
		MultiplexMode _parDiffusionMode; //!< Determines the multiplex of particle diffusion, needed for sensitivitites
		std::vector<active> _parSurfDiffusion; //!< Particle surface diffusion coefficient \f$ D_s \f$
		MultiplexMode _parSurfDiffusionMode; //!< Determines the multiplex of surface diffusion, needed for sensitivitites
		IParameterStateDependence* _parDepSurfDiffusion; //!< Parameter dependence of surface diffusion
		bool _paramDepSurfDiffTypeDep; //!< Determines whether parameter dependence of surface diffusion is particle type dependent
		bool _hasParDepSurfDiffusion; //!< Determines whether particle surface diffusion parameter dependencies are present
		bool _hasSurfaceDiffusion; //!< Determines whether surface diffusion is present

		/* binding and reactions */
		const int* _reqBinding; //!< Array of size @p _strideBound with flags whether binding is in rapid equilibrium for each bound state
		bool _hasDynamicReactions; //! Determines whether or not the binding has any dynamic reactions

		/* discretization - state vector strides, offsets */
		unsigned int _parTypeIdx; //!< Particle type index (wrt the unit operation that owns this particle model)
		int _strideBulkComp; //!< Component stride in bulk state vector
		unsigned int* _boundOffset; //!< Array with offset to the first bound state of each component in the solid phase
		unsigned int _strideBound; //!< Total number of bound states
		unsigned int _nParPoints; //!< Total number of discrete points per particle
		inline int strideBulkComp() const CADET_NOEXCEPT { return _strideBulkComp; }
		inline int strideParComp() const CADET_NOEXCEPT { return 1; }
		inline int strideParLiquid() const CADET_NOEXCEPT { return static_cast<int>(_nComp); }
		inline int strideParBound() const CADET_NOEXCEPT { return static_cast<int>(_strideBound); }
		inline int strideParPoint() const CADET_NOEXCEPT { return strideParLiquid() + strideParBound(); }
		inline int strideParBlock() const CADET_NOEXCEPT { return static_cast<int>(_nParPoints) * strideParPoint(); }

	};

} // namespace parts
} // namespace model
} // namespace cadet

#endif  // LIBCADET_PARTCICLEDIFFUSIONOPERATORBASE_HPP_