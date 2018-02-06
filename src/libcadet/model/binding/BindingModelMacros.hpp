// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2018: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file 
 * Defines some useful macros for IBindingModel implementations.
 */

#ifndef LIBCADET_BINDINGMODELMACROS_HPP_
#define LIBCADET_BINDINGMODELMACROS_HPP_

/**
 * @brief Inserts implementations of all residual() method variants
 * @details An IBindingModel implementation has to provide residual() methods for different variants of state and
 *          parameter type. This macro saves some time by providing those implementations. It assumes that the 
 *          implementation provides a templatized residualImpl() function which realizes all required variants.
 *          
 *          The implementation is inserted inline in the class declaration.
 */
#define CADET_BINDINGMODEL_RESIDUAL_BOILERPLATE                                                                     \
	virtual int residual(const active& t, double z, double r, unsigned int secIdx, const active& timeFactor,        \
		active const* y, double const* yDot, active* res) const                                                     \
	{                                                                                                               \
		return residualImpl<active, active, active, active>(t, z, r, secIdx, timeFactor, y, y - _nComp, yDot, res); \
	}                                                                                                               \
	                                                                                                                \
	virtual int residual(double t, double z, double r, unsigned int secIdx, double timeFactor,                      \
		active const* y, double const* yDot, active* res) const                                                     \
	{                                                                                                               \
		return residualImpl<active, active, active, double>(t, z, r, secIdx, timeFactor, y, y - _nComp, yDot, res); \
	}                                                                                                               \
	                                                                                                                \
	virtual int residual(const active& t, double z, double r, unsigned int secIdx, const active& timeFactor,        \
		double const* y, double const* yDot, active* res) const                                                     \
	{                                                                                                               \
		return residualImpl<double, double, active, active>(t, z, r, secIdx, timeFactor, y, y - _nComp, yDot, res); \
	}                                                                                                               \
	                                                                                                                \
	virtual int residual(double t, double z, double r, unsigned int secIdx, double timeFactor,                      \
		double const* y, double const* yDot, double* res) const                                                     \
	{                                                                                                               \
		return residualImpl<double, double, double, double>(t, z, r, secIdx, timeFactor, y, y - _nComp, yDot, res); \
	}


/**
 * @brief Inserts implementations of all residual(), analyticJacobian(), and jacobianAddDiscretized() method variants
 * @details An IBindingModel implementation has to provide residual(), analyticJacobian(), and jacobianAddDiscretized()
 *          methods for different variants of state and parameter type. This macro saves some time by providing those
 *          implementations. It assumes that the implementation provides a templatized residualImpl(), jacobianImpl(),
 *          and jacobianAddDiscretizedImpl() function which realize all required variants.
 *          
 *          The implementation is inserted inline in the class declaration.
 */
#define CADET_BINDINGMODEL_RESIDUAL_JACOBIAN_BOILERPLATE                                                            \
	CADET_BINDINGMODEL_RESIDUAL_BOILERPLATE                                                                         \
	virtual void analyticJacobian(double t, double z, double r, unsigned int secIdx, double const* y,               \
		linalg::BandMatrix::RowIterator jac) const                                                                  \
	{                                                                                                               \
		jacobianImpl(t, z, r, secIdx, y, y - _nComp, jac);                                                          \
	}                                                                                                               \
	                                                                                                                \
	virtual void analyticJacobian(double t, double z, double r, unsigned int secIdx, double const* y,               \
		linalg::DenseBandedRowIterator jac) const                                                                   \
	{                                                                                                               \
		jacobianImpl(t, z, r, secIdx, y, y - _nComp, jac);                                                          \
	}                                                                                                               \
	                                                                                                                \
	virtual void jacobianAddDiscretized(double alpha, linalg::FactorizableBandMatrix::RowIterator jac) const        \
	{                                                                                                               \
		jacobianAddDiscretizedImpl(alpha, jac);                                                                     \
	}                                                                                                               \
	                                                                                                                \
	virtual void jacobianAddDiscretized(double alpha, linalg::DenseBandedRowIterator jac) const                     \
	{                                                                                                               \
		jacobianAddDiscretizedImpl(alpha, jac);                                                                     \
	}


/**
 * @brief Inserts implementations of all residual(), analyticJacobian(), and jacobianAddDiscretized() method variants
 * @details An IBindingModel implementation has to provide residual(), analyticJacobian(), and jacobianAddDiscretized()
 *          methods for different variants of state and parameter type. This macro saves some time by providing those
 *          implementations. It assumes that the implementation provides a templatized residualImpl(), jacobianImpl(),
 *          and jacobianAddDiscretizedImpl() function which realize all required variants.
 *          
 *          The implementation is inserted as a standalone function definition outside of class declaration.
 * 
 * @param CLASSNAME Name of the IBindingModel implementation (including template)
 * @param TEMPLATELINE Line before each function that may contain a template<typename TEMPLATENAME> modifier
 */
#define CADET_BINDINGMODEL_RESIDUAL_JACOBIAN_BOILERPLATE_IMPL_BASE(CLASSNAME, TEMPLATELINE)                         \
	TEMPLATELINE                                                                                                    \
	int CLASSNAME::residual(const active& t, double z, double r, unsigned int secIdx, const active& timeFactor,     \
		active const* y, double const* yDot, active* res) const                                                     \
	{                                                                                                               \
		return residualImpl<active, active, active, active>(t, z, r, secIdx, timeFactor, y, y - _nComp, yDot, res); \
	}                                                                                                               \
	                                                                                                                \
	TEMPLATELINE                                                                                                    \
	int CLASSNAME::residual(double t, double z, double r, unsigned int secIdx, double timeFactor,                   \
		active const* y, double const* yDot, active* res) const                                                     \
	{                                                                                                               \
		return residualImpl<active, active, active, double>(t, z, r, secIdx, timeFactor, y, y - _nComp, yDot, res); \
	}                                                                                                               \
	                                                                                                                \
	TEMPLATELINE                                                                                                    \
	int CLASSNAME::residual(const active& t, double z, double r, unsigned int secIdx, const active& timeFactor,     \
		double const* y, double const* yDot, active* res) const                                                     \
	{                                                                                                               \
		return residualImpl<double, double, active, active>(t, z, r, secIdx, timeFactor, y, y - _nComp, yDot, res); \
	}                                                                                                               \
	                                                                                                                \
	TEMPLATELINE                                                                                                    \
	int CLASSNAME::residual(double t, double z, double r, unsigned int secIdx, double timeFactor,                   \
		double const* y, double const* yDot, double* res) const                                                     \
	{                                                                                                               \
		return residualImpl<double, double, double, double>(t, z, r, secIdx, timeFactor, y, y - _nComp, yDot, res); \
	}                                                                                                               \
	                                                                                                                \
	TEMPLATELINE                                                                                                    \
	void CLASSNAME::analyticJacobian(double t, double z, double r, unsigned int secIdx, double const* y,            \
		linalg::BandMatrix::RowIterator jac) const                                                                  \
	{                                                                                                               \
		jacobianImpl(t, z, r, secIdx, y, y - _nComp, jac);                                                          \
	}                                                                                                               \
	                                                                                                                \
	TEMPLATELINE                                                                                                    \
	void CLASSNAME::analyticJacobian(double t, double z, double r, unsigned int secIdx, double const* y,            \
		linalg::DenseBandedRowIterator jac) const                                                                   \
	{                                                                                                               \
		jacobianImpl(t, z, r, secIdx, y, y - _nComp, jac);                                                          \
	}                                                                                                               \
	                                                                                                                \
	TEMPLATELINE                                                                                                    \
	void CLASSNAME::jacobianAddDiscretized(double alpha, linalg::FactorizableBandMatrix::RowIterator jac) const     \
	{                                                                                                               \
		jacobianAddDiscretizedImpl(alpha, jac);                                                                     \
	}                                                                                                               \
	                                                                                                                \
	TEMPLATELINE                                                                                                    \
	void CLASSNAME::jacobianAddDiscretized(double alpha, linalg::DenseBandedRowIterator jac) const                  \
	{                                                                                                               \
		jacobianAddDiscretizedImpl(alpha, jac);                                                                     \
	}

/**
 * @brief Inserts implementations of all residual() method variants which forward to residualImpl() template function
 * @details An IBindingModel implementation has to provide residual() methods for different variants of state and
 *          parameter type. This macro saves some time by providing those implementations. It assumes that the
 *          implementation provides a templatized residualImpl() function which realizes all required variants.
 *          
 *          The implementation is inserted as a standalone function definition outside of class declaration.
 * 
 * @param CLASSNAME Name of the IBindingModel implementation
 */
#define CADET_BINDINGMODEL_RESIDUAL_BOILERPLATE_IMPL(CLASSNAME)                                                    \
	CADET_BINDINGMODEL_RESIDUAL_BOILERPLATE_IMPL_BASE(CLASSNAME,)

/**
 * @brief Inserts implementations of all residual() method variants which forward to residualImpl() template function
 * @details An IBindingModel implementation has to provide residual() methods for different variants of state and
 *          parameter type. This macro saves some time by providing those implementations. It assumes that the
 *          implementation provides a templatized residualImpl() function which realizes all required variants.
 *          
 *          The implementation is inserted as a standalone function definition outside of class declaration.
 * 
 * @param CLASSNAME Name of the IBindingModel implementation
 * @param TEMPLATENAME Name of the template parameter that handles externally dependent binding models (optional)
 */
#define CADET_BINDINGMODEL_RESIDUAL_TEMPLATED_BOILERPLATE_IMPL(CLASSNAME,TEMPLATENAME)                             \
	CADET_BINDINGMODEL_RESIDUAL_BOILERPLATE_IMPL_BASE(CLASSNAME<TEMPLATENAME>, template<typename TEMPLATENAME>)

#endif  // LIBCADET_BINDINGMODELMACROS_HPP_
