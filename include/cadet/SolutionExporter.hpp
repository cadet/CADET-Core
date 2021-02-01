// =============================================================================
//  CADET
//  
//  Copyright © 2008-2020: The CADET Authors
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
	ParticleType,
	ParticleShell,
	BoundState
};

/**
 * @brief Interface providing functionality for exporting the solution to the user space
 */
class CADET_API ISolutionExporter
{
public:

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
	 * @brief Returns whether the associated model has a solid phase
	 * @return @c true if solid phase is present, otherwise @c false
	 */
	virtual bool hasSolidPhase() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns whether the associated model has volume DOFs
	 * @details Models that can accumulate (e.g., CSTR) usually have varying volume.
	 * @return @c true if volume DOFs are present, otherwise @c false
	 */
	virtual bool hasVolume() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the number of components
	 * @return Number of components
	 */
	virtual unsigned int numComponents() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the number of axial column cells
	 * @return Number of axial column cells
	 */
	virtual unsigned int numAxialCells() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the number of radial column cells
	 * @return Number of radial column cells
	 */
	virtual unsigned int numRadialCells() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the number of inlet ports
	 * @return Number of inlet ports
	 */
	virtual unsigned int numInletPorts() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the number of outlet ports
	 * @return Number of outlet ports
	 */
	virtual unsigned int numOutletPorts() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the number of particle types
	 * @return Number of particle types
	 */
	virtual unsigned int numParticleTypes() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the number of particle shells
	 * @param [in] parType Particle type index
	 * @return Number of particle shells
	 */
	virtual unsigned int numParticleShells(unsigned int parType) const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the total number of bound states for all components of a given particle type
	 * @param [in] parType Particle type index
	 * @return Total number of bound states for all components of a particle type
	 */
	virtual unsigned int numBoundStates(unsigned int parType) const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the total number of DOFs in the interstitial bulk volume
	 * @return Number of main / bulk DOFs
	 */
	virtual unsigned int numBulkDofs() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the total number of mobile phase DOFs in the particles, if particles are supported
	 * @details The total number of DOFs is returned, i.e., the sum of all particle cells' mobile phase DOFs.
	 *          This includes all particle types.
	 * 
	 * @param [in] parType Particle type index
	 * @return Total number of particle mobile phase DOFs
	 */
	virtual unsigned int numParticleMobilePhaseDofs(unsigned int parType) const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the total number of solid phase DOFs
	 * @details The total number of DOFs is returned, i.e., the sum of all column and particle cells' solid phase DOFs.
	 *          This includes all particle types.
	 * 
	 * @param [in] parType Particle type index
	 * @return Total number of solid phase DOFs
	 */
	virtual unsigned int numSolidPhaseDofs(unsigned int parType) const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the number of bulk-bead flux DOFs, if particle fluxes are supported
	 * @return Total number of particle flux DOFs
	 */
	virtual unsigned int numFluxDofs() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the number of volume DOFs, if volume DOFs are supported
	 * @return Number of volume DOFs
	 */
	virtual unsigned int numVolumeDofs() const CADET_NOEXCEPT = 0;


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
	 * @brief Provides direct access to the underlying particle mobile phase state vector part of the given particle type
	 * @details The ordering of the (full) data inside the state vector is provided by mobilePhaseOrdering()
	 * @param [in] parType Particle type index
	 * @return Pointer to the first element of the state vector or @c NULL if the model does not support it
	 */
	virtual double const* particleMobilePhase(unsigned int parType) const = 0;

	/**
	 * @brief Provides direct access to the underlying solid phase state vector part of the given particle type
	 * @details The ordering of the (full) data inside the state vector is provided by solidPhaseOrdering()
	 * @param [in] parType Particle type index
	 * @return Pointer to the first element of the state vector or @c NULL if the model does not support it
	 */
	virtual double const* solidPhase(unsigned int parType) const = 0;

	/**
	 * @brief Provides direct access to the underlying volume slice of the state vector
	 * @return Pointer to the first element of the volume slice or @c NULL if the model does not support volume DOFs
	 */
	virtual double const* volume() const = 0;

	/**
	 * @brief Provides direct access to the inlet state vector
	 * @details The inlet state vector only contains one value for each main mobile phase component (see numComponents()).
	 *          The stride required for the access is returned in @p stride.
	 * @param [in] port Index of the port
	 * @param [out] stride Stride of the vector access
	 * @return Pointer to the first element of the inlet state vector
	 */
	virtual double const* inlet(unsigned int port, unsigned int& stride) const = 0;

	/**
	 * @brief Provides direct access to the outlet state vector
	 * @details The outlet state vector only contains one value for each main mobile phase component (see numComponents()).
	 *          The stride required for the access is returned in @p stride.
	 * @param [in] port Index of the port
	 * @param [out] stride Stride of the vector access
	 * @return Pointer to the first element of the outlet state vector
	 */
	virtual double const* outlet(unsigned int port, unsigned int& stride) const = 0;


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

	/**
	 * @brief Returns the number of elements between two bulk mobile phase DOF blocks
	 * @details Stride between two bulk mobile phase DOF blocks.
	 * @return Number of elements between two bulk mobile phase DOF blocks
	 */
	virtual unsigned int bulkMobilePhaseStride() const = 0;

	/**
	 * @brief Returns the number of elements between two particle mobile phase DOF blocks
	 * @details Stride between two particle mobile phase DOF blocks of the given particle type.
	 * @param [in] parType Particle type index
	 * @return Number of elements between two particle mobile phase DOF blocks
	 */
	virtual unsigned int particleMobilePhaseStride(unsigned int parType) const = 0;

	/**
	 * @brief Returns the number of elements between two solid phase DOF blocks
	 * @details Stride between two solid phase DOF blocks of the given particle type.
	 * @param [in] parType Particle type index
	 * @return Number of elements between two solid phase DOF blocks
	 */
	virtual unsigned int solidPhaseStride(unsigned int parType) const = 0;

	/**
	 * @brief Returns axial coordinates of the DOFs
	 * @param [out] Pointer to array that is filled with axial coordinates of DOFs
	 */
	virtual void axialCoordinates(double* coords) const = 0;

	/**
	 * @brief Returns radial coordinates of the DOFs
	 * @param [out] Pointer to array that is filled with radial coordinates of DOFs
	 */
	virtual void radialCoordinates(double* coords) const = 0;

	/**
	 * @brief Returns particle coordinates of the DOFs
	 * @param [in] parType Particle type index
	 * @param [out] Pointer to array that is filled with particle coordinates of DOFs
	 */
	virtual void particleCoordinates(unsigned int parType, double* coords) const = 0;
};

} // namespace cadet

#endif  // LIBCADET_SOLUTIONEXPORTER_HPP_
