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
			* @brief Sets fixed parameters of the particle diffusion operator (e.g., the number of discretization points, components, bound states)
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

		unsigned int _parTypeIdx; //!< Particle type index (wrt the unit operation that owns this particle model)
		unsigned int _nComp; //!< Number of components
		int _strideBulkComp; //!< Component stride in bulk state vector
		double _parGeomSurfToVol; //!< Particle surface to volume ratio factor (i.e., 3.0 for spherical, 2.0 for cylindrical, 1.0 for hexahedral)
		std::vector<active> _filmDiffusion; //!< Particle diffusion coefficient \f$ D_p \f$
		//MultiplexMode _filmDiffusionMode;
		std::vector<active> _poreAccessFactor; //!< Pore accessibility factor \f$ F_{\text{acc}} \f$
		//MultiplexMode _poreAccessFactorMode;
		std::vector<active> _invBetaP; //!< Ratio of solid to liquid particle volume
		unsigned int* _nBound; //!< Array with number of bound states for each component
		unsigned int* _boundOffset; //!< Array with offset to the first bound state of each component in the solid phase
		unsigned int _strideBound; //!< Total number of bound states
		const int* _reqBinding; //!< Array of size @p _strideBound with flags whether binding is in rapid equilibrium for each bound state
		bool _hasDynamicReactions; //! Determines whether or not the binding has any dynamic reactions
		active _parRadius; //!< Particle radius \f$ r_p \f$
		bool _singleParRadius;
		active _parCoreRadius; //!< Particle core radius \f$ r_c \f$
		bool _singleParCoreRadius;
		active _parPorosity; //!< Particle porosity (internal porosity) \f$ \varepsilon_p \f$
		bool _singleParPorosity;
		std::vector<active> _parDiffusion; //!< Particle diffusion coefficient \f$ D_p \f$
		MultiplexMode _parDiffusionMode;
		std::vector<active> _parSurfDiffusion; //!< Particle surface diffusion coefficient \f$ D_s \f$
		MultiplexMode _parSurfDiffusionMode;
		IParameterStateDependence* _parDepSurfDiffusion; //!< Parameter dependencies for particle surface diffusion
		bool _singleParDepSurfDiffusion; //!< Determines whether a single parameter dependence for particle surface diffusion is used
		bool _hasParDepSurfDiffusion; //!< Determines whether particle surface diffusion parameter dependencies are present
		bool _hasSurfaceDiffusion; //!< Determines whether surface diffusion is present
	};

} // namespace parts
} // namespace model
} // namespace cadet

#endif  // LIBCADET_PARTCICLEDIFFUSIONOPERATORBASE_HPP_