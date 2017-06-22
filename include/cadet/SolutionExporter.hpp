// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2017: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file 
 * Defines interfaces which export the solution to the user space.
 */

#ifndef LIBCADET_SOLUTIONEXPORTER_HPP_
#define LIBCADET_SOLUTIONEXPORTER_HPP_

#include <cstdint>

#include "cadet/LibExportImport.hpp"
#include "cadet/cadetCompilerInfo.hpp"

namespace cadet
{

/**
 * @brief Elements of state vector ordering
 */
enum class StateOrdering : uint8_t
{
	Component,
	AxialCell,
	RadialCell,
	Phase
};

/**
 * @brief Interface providing functionality for exporting the solution to the user space
 */
class CADET_API ISolutionExporter
{
public:

	/**
	 * @brief Returns whether the associated model has multiple bound states
	 * @details Components with multiple bound states have several entries in the solid phase.
	 * @return @c true if multiple bound states are present, otherwise @c false
	 */
	virtual bool hasMultipleBoundStates() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns whether the associated model has a flux into the particles
	 * @return @c true if flux into the particles is present, otherwise @c false
	 */
	virtual bool hasParticleFlux() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns whether the associated model has a particle mobile phase
	 * @return @c true if particle mobile phase is present, otherwise @c false
	 */
	virtual bool hasParticleMobilePhase() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns whether the associated model has non-binding components
	 * @details Non-binding components do not have a solid phase representation.
	 * @return @c true if non-binding components are present, otherwise @c false
	 */
	virtual bool hasNonBindingComponents() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the number of components
	 * @return Number of components
	 */
	virtual unsigned int numComponents() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the number of axial cells
	 * @return Number of axial cells
	 */
	virtual unsigned int numAxialCells() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the number of radial cells
	 * @return Number of radial cells
	 */
	virtual unsigned int numRadialCells() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the number of total bound states for all components
	 * @return Total number of bound states for all components
	 */
	virtual unsigned int numBoundStates() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns an array with the number of bound states for each component
	 * @return Pointer to array holding the number of bound states for each component
	 */
	virtual unsigned int const* numBoundStatesPerComponent() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the number of bound states of the given component
	 * @param [in] comp Index of the component
	 * @return Number of bound states of component @p comp
	 */
	virtual unsigned int numBoundStates(unsigned int comp) const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the number of main DOFs in the column / interstitial bulk volume
	 * @return Number of main DOFs
	 */
	virtual unsigned int numColumnDofs() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the number of all DOFs in the particles, if particles are supported
	 * @return Total number of particle DOFs
	 */
	virtual unsigned int numParticleDofs() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the number of bulk-bead flux DOFs, if particle fluxes are supported
	 * @return Total number of particle flux DOFs
	 */
	virtual unsigned int numFluxDofs() const CADET_NOEXCEPT = 0;


	/**
	 * @brief Returns a concentration value
	 * @details This refers to the main mobile phase (e.g., interstitial concentration).
	 *          Every unit operation model possesses a main mobile phase.
	 * 
	 * @param [in] component Component index
	 * @param [in] axialCell Index of the axial cell
	 * 
	 * @return Concentration at the given point
	 */
	virtual double concentration(unsigned int component, unsigned int axialCell) const = 0;

	/**
	 * @brief Returns a flux from the main mobile phase to the particle
	 * @details Not every model has fluxes.
	 * 
	 * @param [in] component Component index
	 * @param [in] axialCell Index of the axial cell
	 * 
	 * @return Flux at the given point
	 */
	virtual double flux(unsigned int component, unsigned int axialCell) const = 0;

	/**
	 * @brief Returns a concentration in the particle mobile phase
	 * @details Not every model has particle mobile phases.
	 * 
	 * @param [in] component Component index
	 * @param [in] axialCell Index of the axial cell
	 * @param [in] radialCell Index of the radial cell
	 * 
	 * @return Concentration at the given point
	 */
	virtual double mobilePhase(unsigned int component, unsigned int axialCell, unsigned int radialCell) const = 0;

	/**
	 * @brief Returns a concentration in the solid phase
	 * @details Not every model has solid phases.
	 * 
	 * @param [in] component Component index
	 * @param [in] axialCell Index of the axial cell
	 * @param [in] radialCell Index of the radial cell
	 * @param [in] boundState Index of the bound state
	 * 
	 * @return Concentration at the given point
	 */
	virtual double solidPhase(unsigned int component, unsigned int axialCell, unsigned int radialCell, unsigned int boundState) const = 0;


	/**
	 * @brief Provides direct access to the underlying main mobile phase state vector
	 * @details The ordering of the data inside the state vector is provided by concentrationOrdering()
	 * @return Pointer to the first element of the state vector
	 */
	virtual double const* concentration() const = 0;

	/**
	 * @brief Provides direct access to the underlying flux state vector
	 * @details The ordering of the data inside the state vector is provided by fluxOrdering()
	 * @return Pointer to the first element of the state vector or @c NULL if the model does not support it
	 */
	virtual double const* flux() const = 0;

	/**
	 * @brief Provides direct access to the underlying particle mobile phase state vector
	 * @details The ordering of the data inside the state vector is provided by mobilePhaseOrdering()
	 * @return Pointer to the first element of the state vector or @c NULL if the model does not support it
	 */
	virtual double const* mobilePhase() const = 0;

	/**
	 * @brief Provides direct access to the underlying solid phase state vector
	 * @details The ordering of the data inside the state vector is provided by solidPhaseOrdering()
	 * @return Pointer to the first element of the state vector or @c NULL if the model does not support it
	 */
	virtual double const* solidPhase() const = 0;

	/**
	 * @brief Provides direct access to the inlet state vector
	 * @details The inlet state vector only contains one value for each main mobile phase component (see numComponents()).
	 *          The stride required for the access is returned in @p stride.
	 * @param [out] stride Stride of the vector access
	 * @return Pointer to the first element of the inlet state vector
	 */
	virtual double const* inlet(unsigned int& stride) const = 0;

	/**
	 * @brief Provides direct access to the outlet state vector
	 * @details The outlet state vector only contains one value for each main mobile phase component (see numComponents()).
	 *          The stride required for the access is returned in @p stride.
	 * @param [out] stride Stride of the vector access
	 * @return Pointer to the first element of the outlet state vector
	 */
	virtual double const* outlet(unsigned int& stride) const = 0;


	/**
	 * @brief Returns an array describing the ordering of the main mobile phase state vector
	 * @details A pointer to the first element of the state vector ordering array is returned. The length
	 *          of the array is returned in the parameter @p len. The elements of the array indicate the
	 *          order of loops required to extract the data.
	 * 
	 * @param [out] len Length of the returned ordering vector
	 * @return Pointer to first element of the ordering vector
	 */
	virtual StateOrdering const* concentrationOrdering(unsigned int& len) const = 0;

	/**
	 * @brief Returns an array describing the ordering of the flux state vector
	 * @details A pointer to the first element of the state vector ordering array is returned. The length
	 *          of the array is returned in the parameter @p len. The elements of the array indicate the
	 *          order of loops required to extract the data.
	 * 
	 * @param [out] len Length of the returned ordering vector
	 * @return Pointer to first element of the ordering vector or @c NULL if fluxes are not supported
	 */
	virtual StateOrdering const* fluxOrdering(unsigned int& len) const = 0;

	/**
	 * @brief Returns an array describing the ordering of the mobile phase state vector
	 * @details A pointer to the first element of the state vector ordering array is returned. The length
	 *          of the array is returned in the parameter @p len. The elements of the array indicate the
	 *          order of loops required to extract the data.
	 * 
	 * @param [out] len Length of the returned ordering vector
	 * @return Pointer to first element of the ordering vector or @c NULL if mobile phases are not supported
	 */
	virtual StateOrdering const* mobilePhaseOrdering(unsigned int& len) const = 0;

	/**
	 * @brief Returns an array describing the ordering of the solid phase state vector
	 * @details A pointer to the first element of the state vector ordering array is returned. The length
	 *          of the array is returned in the parameter @p len. The elements of the array indicate the
	 *          order of loops required to extract the data.
	 * 
	 * @param [out] len Length of the returned ordering vector
	 * @return Pointer to first element of the ordering vector or @c NULL if solid phases are not supported
	 */
	virtual StateOrdering const* solidPhaseOrdering(unsigned int& len) const = 0;
};

} // namespace cadet

#endif  // LIBCADET_SOLUTIONEXPORTER_HPP_
