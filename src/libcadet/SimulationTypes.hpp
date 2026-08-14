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
	double primary; //!< primary bulk coordinate
	double secondary; //!< secondary bulk coordinate
	double particle; //!< Radial particle coordinate
};

/**
 * @brief Quadrature rule for averaging a spatially varying exchange coefficient onto the bulk points
 * @details Bulk-particle exchange terms (film diffusion) carry a coefficient that may vary in space.
 *          The bulk discretization provides the positions at which the coefficient has to be evaluated
 *          and the weights that average these samples onto a bulk point,
 *          @f[ \tilde{k}_i = \sum_q w_{iq} k\left( \rho_{iq} \right), @f]
 *          so that the phase owning the coefficient can evaluate it without knowing the bulk
 *          discretization, and the bulk discretization does not need to know the coefficient.
 *          The weights of a bulk point sum to one, hence a constant coefficient is reproduced exactly.
 *          A single quadrature point with weight one recovers evaluation at the bulk point itself.
 */
struct CouplingQuadrature
{
	int nQuadPoints; //!< Number of quadrature points per bulk point
	ColumnPosition const* colPos; //!< Positions in relative coordinates, bulk point major
	active const* velocity; //!< Interstitial velocity at the positions, bulk point major
	double const* weights; //!< Averaging weights, bulk point major

	/**
	 * @brief Returns the rule belonging to a single bulk point
	 */
	inline CouplingQuadrature slice(int point) const CADET_NOEXCEPT
	{
		return CouplingQuadrature{ nQuadPoints, colPos + point * nQuadPoints, velocity + point * nQuadPoints, weights + point * nQuadPoints };
	}
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
 * @brief Type tag for including parameter sensitivities
 */
struct WithParamSensitivity {};

/**
 * @brief Type tag for excluding parameter sensitivities
 */
struct WithoutParamSensitivity {};

/**
 * @brief Type used for choosing the parameter sensitivity tag type
 * @details Chooses WithoutParamSensitivity if @p T is @c double and
 *          WithParamSensitivity if @p T is @c active.
 * @tparam T Parameter data type
 */
template <typename T>
struct ParamSens { };

template <>
struct ParamSens<double> { typedef WithoutParamSensitivity enabled; };

template <>
struct ParamSens<active> { typedef WithParamSensitivity enabled; };

} // namespace cadet

#endif  // LIBCADET_SIMULATIONTYPES_HPP_
