// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2019: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file 
 * Defines useful types for the simulation.
 */

#ifndef LIBCADET_SIMULATIONTYPES_HPP_
#define LIBCADET_SIMULATIONTYPES_HPP_

#include "AutoDiff.hpp"

namespace cadet
{

/**
 * @brief Position inside a column unit operation
 */
struct ColumnPosition
{
	double axial; //!< Axial bulk coordinate z
	double radial; //!< Radial bulk coordinate rho
	double particle; //!< Radial particle coordinate r
};

/**
 * @brief Common parameters for Jacobians via AD
 */
struct AdJacobianParams
{
	active* adRes; //!< Residual vector
	active* adY; //!< State vector
	unsigned int adDirOffset; //!< Number of reserved AD directions (e.g., for sensitivities)
};

/**
 * @brief Simulation time point
 */
struct SimulationTime
{
	double t; //!< Time of the simulation
	unsigned int secIdx; //!< Index of the current section
	double timeFactor; //!< Used for time transformation (pre factor of time derivatives) and to compute parameter derivatives with respect to section length
};

/**
 * @brief Simulation time point
 */
struct ActiveSimulationTime
{
	const active& t; //!< Time of the simulation
	unsigned int secIdx; //!< Index of the current section
	const active& timeFactor; //!< Used for time transformation (pre factor of time derivatives) and to compute parameter derivatives with respect to section length
};

/**
 * @brief State vectors of the simulation that can be updated
 */
struct SimulationState
{
	double* vecStateY; //!< State vector
	double* vecStateYdot; //!< Time derivative of state vector
};

/**
 * @brief State vectors of the simulation that cannot be changed
 */
struct ConstSimulationState
{
	double const* vecStateY; //!< State vector
	double const* vecStateYdot; //!< Time derivative of state vector
};

/**
 * @brief Converts a mutable simulation state to a constant one
 * @param [in] simState Mutable simulation state
 * @return Constant simulation state
 */
inline ConstSimulationState toConst(const SimulationState& simState)
{
	return ConstSimulationState{simState.vecStateY, simState.vecStateYdot};
}

/**
 * @brief Converts a simple simulation time point to an AD-enabled (i.e., active) one
 * @param [in] simTime Simple simulation time
 * @return AD-enabled simulation time
 */
inline ActiveSimulationTime toActive(const SimulationTime& simTime)
{
	return ActiveSimulationTime{active(simTime.t), simTime.secIdx, active(simTime.timeFactor)};
}

/**
 * @brief Converts an AD-enabled (i.e., active) simulation time point to a simple one
 * @param [in] simTime AD-enabled simulation time
 * @return Simple simulation time
 */
inline SimulationTime toSimple(const ActiveSimulationTime& simTime)
{
	return SimulationTime{static_cast<double>(simTime.t), simTime.secIdx, static_cast<double>(simTime.timeFactor)};
}

} // namespace cadet

#endif  // LIBCADET_SIMULATIONTYPES_HPP_
