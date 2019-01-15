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
 * Defines tests for unit operations with particles.
 */

#ifndef CADETTEST_PARTICLEHELPER_HPP_
#define CADETTEST_PARTICLEHELPER_HPP_

#include "cadet/ParameterId.hpp"

#include <vector>

namespace cadet
{

class JsonParameterProvider;

namespace test
{

namespace particle
{

	/**
	 * @brief Extends a model to multiple particle types by replicating the first type
	 * @param [in,out] jpp ParameterProvider to change the number of axial cells in
	 * @param [in] unit Index of unit operation
	 * @param [in] nTypes Total number of particle types
	 * @param [in] volFrac Array with volume fractions of particle types
	 */
	void extendModelToManyParticleTypes(cadet::JsonParameterProvider& jpp, UnitOpIdx unit, unsigned int nTypes, double const* const volFrac);

	/**
	 * @brief Extends a model to multiple particle types by replicating the first type
	 * @param [in,out] jpp ParameterProvider to change the number of axial cells in
	 * @param [in] nTypes Total number of particle types
	 * @param [in] volFrac Array with volume fractions of particle types
	 */
	void extendModelToManyParticleTypes(cadet::JsonParameterProvider& jpp, unsigned int nTypes, double const* const volFrac);

	/**
	 * @brief Sets the volume fractions of the particle types
	 * @param [in,out] jpp ParameterProvider to change the number of axial cells in
	 * @param [in] unit Index of unit operation
	 * @param [in] volFrac Array with volume fractions of particle types
	 */
	void setParticleTypeVolumeFractions(cadet::JsonParameterProvider& jpp, UnitOpIdx unit, const std::vector<double>& volFrac);

	/**
	 * @brief Sets the volume fractions of the particle types
	 * @param [in,out] jpp ParameterProvider to change the number of axial cells in
	 * @param [in] volFrac Array with volume fractions of particle types
	 */
	void setParticleTypeVolumeFractions(cadet::JsonParameterProvider& jpp, const std::vector<double>& volFrac);

	/**
	 * @brief Checks whether results of simulation with one particle type matches the one with two (identical) particle types
	 * @details The LWE example is taken and simulated. The existing particle type is replicated such that there are
	 *          two identical particle types. Another simulation is run with the two particle types.
	 *          The results of these simulations must match.
	 * 
	 * @param [in] uoType Unit operation type
	 * @param [in] absTols Array with absolute error tolerances
	 * @param [in] relTols Array with relative error tolerances
	 */
	void testOneVsTwoIdenticalParticleTypes(const char* uoType, double absTol, double relTol);

	/**
	 * @brief Checks whether, when using two separate identical particle types, results of one type match the other
	 * @details The LWE example is taken and the existing particle type is replicated such that there are
	 *          two identical particle types. Two simulations are run each one using only one particle type.
	 *          The results of these simulations must match.
	 * 
	 * @param [in] uoType Unit operation type
	 * @param [in] absTols Array with absolute error tolerances
	 * @param [in] relTols Array with relative error tolerances
	 */
	void testSeparateIdenticalParticleTypes(const char* uoType, double absTol, double relTol);

	/**
	 * @brief Checks whether a linear binding model with multiple identical particle types produces the same as result as a single type model
	 * @details The linear benchmark problem is run. Then, two additional identical particle types are added.
	 *          The results of the multi type simulation must match the ones of the single type simulation.
	 *          This check is conducted with both dynamic and quasi-stationary binding.
	 * 
	 * @param [in] uoType Unit operation type
	 * @param [in] absTols Array with absolute error tolerances
	 * @param [in] relTols Array with relative error tolerances
	 */
	void testLinearMixedParticleTypes(const char* uoType, double absTol, double relTol);

	/**
	 * @brief Checks the full analytic Jacobian against AD for a model with multiple particle types
	 * @param [in] uoType Unit operation type
	 */
	void testJacobianMixedParticleTypes(const std::string& uoType);

	/**
	 * @brief Checks the (analytic) time derivative Jacobian against FD for a model with multiple particle types
	 * @details Uses centered finite differences.
	 * @param [in] uoType Unit operation type
	 * @param [in] h Step size of centered finite differences
	 * @param [in] absTol Absolute error tolerance
	 * @param [in] relTol Relative error tolerance
	 */
	void testTimeDerivativeJacobianMixedParticleTypesFD(const std::string& uoType, double h, double absTol, double relTol);

} // namespace particle
} // namespace test
} // namespace cadet

#endif  // CADETTEST_PARTICLEHELPER_HPP_
