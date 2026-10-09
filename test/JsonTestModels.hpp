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
 * Defines a ParameterProvider that uses JSON.
 */

#ifndef CADETTEST_JSONTESTMODELS_HPP_
#define CADETTEST_JSONTESTMODELS_HPP_

#include "common/JsonParameterProvider.hpp"

cadet::JsonParameterProvider createColumnWithSMA(const std::string& uoType, const std::string& spatialScheme);
cadet::JsonParameterProvider createColumnWithTwoCompLinearBinding(const std::string& uoType, const std::string& spatialScheme);
cadet::JsonParameterProvider createColumnLinearBenchmark(bool dynamicBinding, bool nonBinding, const std::string& uoType, const std::string& spatialScheme);
nlohmann::json createLWEJson(const std::string& uoType, const std::string& spatialMethod);
cadet::JsonParameterProvider createLWE(const std::string& uoType, const std::string& spatialScheme);
cadet::JsonParameterProvider createPulseInjectionColumn(const std::string& uoType, const std::string& spatialScheme, bool dynamicBinding);
cadet::JsonParameterProvider createLinearBenchmark(bool dynamicBinding, bool nonBinding, const std::string& uoType, const std::string& spatialScheme);
cadet::JsonParameterProvider createCSTR(unsigned int nComp);
cadet::JsonParameterProvider createCSTRBenchmark(unsigned int nSec, double endTime, double interval);
nlohmann::json createColumn2ParType1GeneralRate1HomoParticleBothWithTwoCompLinearJson(const std::string& uoType, const std::string& spatialMethod);

/**
 * @brief Creates a finite bath with two components and one particle type
 * @param [in] particleType Either @c HOMOGENEOUS_PARTICLE or @c GENERAL_RATE_PARTICLE
 * @param [in] parMethod Spatial discretization of general rate particles, either @c FV or @c DG
 * @param [in] binding Determines whether a linear binding model with one bound state per component is used
 */
nlohmann::json createFiniteBathJson(const std::string& particleType, const std::string& parMethod, bool binding = true);
cadet::JsonParameterProvider createFiniteBath(const std::string& particleType, const std::string& parMethod, bool binding = true);
cadet::JsonParameterProvider createFiniteBathBenchmark(const std::string& particleType, const std::string& parMethod, bool binding, unsigned int nSec, double endTime, double interval);

#endif  // CADETTEST_JSONTESTMODELS_HPP_
