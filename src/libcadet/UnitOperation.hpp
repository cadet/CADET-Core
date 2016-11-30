// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2016: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file 
 * Defines an unit operation model interface.
 */

#ifndef LIBCADET_IUNITOPERATION_HPP_
#define LIBCADET_IUNITOPERATION_HPP_

#include "cadet/Model.hpp"
#include "SimulatableModel.hpp"

namespace cadet
{

class IExternalFunction;

/**
 * @brief Defines an unit operation model interface
 * @details This interface is mainly used to connect different unit operations with each other
 *          in the IModelSystem. The connection is either made by accessing inlet / outlet DOFs
 *          from the local slice of the global state vector or, in case of an inlet, by taking
 *          some data (local array) as source. The implementing class only has to provide one
 *          of those two means depending on the result of numDofs(). If <tt>numDofs() == 0</tt>,
 *          then getData() and getDataActive() are used instead of the index() and stride()
 *          functions, which are used for <tt>numDofs() > 0</tt>.
 */
class IUnitOperation : public IModel, public ISimulatableModel
{
public:
	virtual unsigned int numPureDofs() const CADET_NOEXCEPT { return numDofs(); }

	/**
	 * @brief Returns the number of components
	 * @details It is assumed that the number of components is also the number of inputs
	 *          and outputs of the unit operation.
	 * @return Number of components
	 */
	virtual unsigned int numComponents() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns whether this unit operation possesses an inlet
	 * @return @c true if the unit operation can take in a stream, otherwise @c false
	 */
	virtual bool hasInlet() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns whether this unit operation possesses an outlet
	 * @return @c true if the unit operation can output a stream, otherwise @c false
	 */
	virtual bool hasOutlet() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns a pointer to a consecutive array with outlet data
	 * @details Is used as data source if <tt>numDofs() == 0</tt>.
	 * @return Pointer to a consecutive array with outlet data
	 */
	virtual double const* const getData() const CADET_NOEXCEPT = 0;
	virtual active const* const getDataActive() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the factor that is used to connect an outlet with an inlet
	 * @param [in] compIdx Index of the component
	 * @param [in] secIdx Index of the section
	 * @return Factor
	 */
	virtual active inletConnectionFactorActive(unsigned int compIdx, unsigned int secIdx) const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the factor that is used to connect an outlet with an inlet
	 * @param [in] compIdx Index of the component
	 * @param [in] secIdx Index of the section
	 * @return Factor
	 */
	virtual double inletConnectionFactor(unsigned int compIdx, unsigned int secIdx) const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the local index of the first outlet component in the managed state vector slice
	 * @details As each unit operation manages a slice of the global state vector on its own, the
	 *          position of inlet and outlet components as DOF in the slice may differ. This function
	 *          returns the local index in the slice of the first outlet component.
	 * @return Local index of the first outlet component in the state vector slice
	 */
	virtual unsigned int localOutletComponentIndex() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the local index stride in the managed state vector slice
	 * @details This function returns the number of elements between two components in the slice of
	 *          the global state vector.
	 * @return Local stride in the state vector slice
	 */
	virtual unsigned int localOutletComponentStride() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the local index of the first inlet component in the managed state vector slice
	 * @details As each unit operation manages a slice of the global state vector on its own, the
	 *          position of inlet and outlet components as DOF in the slice may differ. This function
	 *          returns the local index in the slice of the first inlet component.
	 * @return Local index of the first inlet component in the state vector slice
	 */
	virtual unsigned int localInletComponentIndex() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the local index stride in the managed state vector slice
	 * @details This function returns the number of elements between two components in the slice of
	 *          the global state vector.
	 * @return Local stride in the state vector slice
	 */
	virtual unsigned int localInletComponentStride() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Evaluates the residuals with AD to compute the sensitivity derivatives
	 * @details Evaluates @f$ \frac{\partial F}{\partial p} @f$, where @f$ p @f$ are the sensitive parameters
	 *          and @f$ F @f$ is the residual function of this model.
	 * @param [in] t Current time point
	 * @param [in] secIdx Index of the current section
	 * @param [in] timeFactor Used for time transformation (pre factor of time derivatives) and to compute parameter derivatives with respect to section length
	 * @param [in] y Pointer to global state vector
	 * @param [in] yDot Pointer to global time derivative state vector
	 * @param [in,out] adRes Pointer to global residual vector of AD datatypes for computing the sensitivity derivatives
	 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
	 */
	virtual int residualSensFwdAdOnly(const active& t, unsigned int secIdx, const active& timeFactor,
		double const* const y, double const* const yDot, active* const adRes) = 0;

	/**
	 * @brief Computes the residual of the forward sensitivity systems using the result of residualSensFwdAdOnly()
	 * @details Assembles and evaluates the residuals of the sensitivity systems
	 *          @f[ \frac{F}{\partial y} s + \frac{F}{\partial \dot{y}} \dot{s} + \frac{\partial F}{\partial p_i} = 0. @f]
	 * @param [in] timeFactor Used for time transformation (pre factor of time derivatives) and to compute parameter derivatives with respect to section length
	 * @param [in] yS Pointers to global sensitivity state vectors
	 * @param [in] ySdot Pointers to global sensitivity time derivative state vectors
	 * @param [out] resS Pointers to global sensitivity residuals
	 * @param [in] adRes Pointer to global residual vector of AD datatypes with the sensitivity derivatives from residualSensFwdAdOnly()
	 * @param [in] tmp1 Temporary storage in the size of global state vector @p y
	 * @param [in] tmp2 Temporary storage in the size of global state vector of @p y
	 * @param [in] tmp3 Temporary storage in the size of global state vector of @p y
	 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
	 */
	virtual int residualSensFwdCombine(const active& timeFactor, const std::vector<const double*>& yS, const std::vector<const double*>& ySdot,
		const std::vector<double*>& resS, active const* adRes, double* const tmp1, double* const tmp2, double* const tmp3) = 0;

	/**
	 * @brief Evaluates the residuals with AD to compute the parameter sensitivities and at the same time updates the Jacobian
	 * @details Evaluates @f$ \frac{\partial F}{\partial p} @f$, where @f$ p @f$ are the sensitive parameters
	 *          and @f$ F @f$ is the residual function of this model. At the same time updates the Jacobian of the system.
	 * @param [in] t Current time point
	 * @param [in] secIdx Index of the current section
	 * @param [in] timeFactor Used for time transformation (pre factor of time derivatives) and to compute parameter derivatives with respect to section length
	 * @param [in] y Pointer to global state vector
	 * @param [in] yDot Pointer to global time derivative state vector
	 * @param [in,out] adRes Pointer to global residual vector of AD datatypes for computing the sensitivity derivatives
	 * @param [in,out] adY Pointer to global state vector of AD datatypes that can be used for computing the Jacobian (or @c nullptr if AD is disabled)
	 * @param [in] numSensAdDirs Number of AD directions used for parameter sensitivities
	 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
	 */
	virtual int residualSensFwdWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, active* const adRes, 
		active* const adY, unsigned int numSensAdDirs) = 0;

	/**
	 * @brief Computes consistent initial values (state variables without their time derivatives)
	 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] the initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ have
	 *          to be consistent. This functions updates the initial state \f$ y_0 \f$ and, therefore, provides the
	 *          first step in consistent initialization.
	 *          
	 *          This function is to be used with consistentInitialTimeDerivative(). Do not mix normal and lean
	 *          consistent initialization!
	 * 
	 * @param [in] t Current time point
	 * @param [in] secIdx Index of the current section
	 * @param [in] timeFactor Used for time transformation (pre factor of time derivatives) and to compute parameter derivatives with respect to section length
	 * @param [in,out] vecStateY State vector with initial values that are to be updated for consistency
	 * @param [in,out] adRes Pointer to global residual vector of AD datatypes that can be used for computing the Jacobian (or @c nullptr if AD is disabled)
	 * @param [in,out] adY Pointer to global state vector of AD datatypes that can be used for computing the Jacobian (or @c nullptr if AD is disabled)
	 * @param [in] numSensAdDirs Number of AD directions used for parameter sensitivities
	 * @param [in] errorTol Error tolerance for algebraic equations
	 */
	virtual void consistentInitialState(double t, unsigned int secIdx, double timeFactor, double* const vecStateY, active* const adRes, active* const adY, unsigned int numSensAdDirs, double errorTol) = 0;

	/**
	 * @brief Computes consistent initial time derivatives
	 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] the initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ have
	 *          to be consistent. This functions overwrites the initial time derivatives \f$ \dot{y}_0 \f$ such that
	 *          the residual is zero. Thus, this function provides the final step in consistent initialization.
	 * 
	 *          This function is to be used with consistentInitialState(). Do not mix normal and lean
	 *          consistent initialization!
	 *          
	 * @param [in] t Current time point
	 * @param [in] timeFactor Used for time transformation (pre factor of time derivatives) and to compute parameter derivatives with respect to section length
	 * @param [in,out] vecStateYdot On entry, residual without taking time derivatives into account. On exit, consistent state time derivatives.
	 */
	virtual void consistentInitialTimeDerivative(double t, double timeFactor, double* const vecStateYdot) = 0;

	/**
	 * @brief Computes approximately / partially consistent initial values (state variables without their time derivatives)
	 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] the initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ have
	 *          to be consistent. This functions updates the initial state \f$ y_0 \f$ and, therefore, provides the
	 *          first step in consistent initialization.
	 *          
	 *          This function is possibly faster than consistentInitialState(), but updates only a part of the
	 *          state vector. Hence, the result is not guaranteed to be consistent.
	 * 
	 *          This function is to be used with leanConsistentInitialTimeDerivative(). Do not mix normal and lean
	 *          consistent initialization!
	 *
	 * @param [in] t Current time point
	 * @param [in] secIdx Index of the current section
	 * @param [in] timeFactor Used for time transformation (pre factor of time derivatives) and to compute parameter derivatives with respect to section length
	 * @param [in,out] vecStateY State vector with initial values that are to be updated for consistency
	 * @param [in,out] adRes Pointer to global residual vector of AD datatypes that can be used for computing the Jacobian (or @c nullptr if AD is disabled)
	 * @param [in,out] adY Pointer to global state vector of AD datatypes that can be used for computing the Jacobian (or @c nullptr if AD is disabled)
	 * @param [in] numSensAdDirs Number of AD directions used for parameter sensitivities
	 * @param [in] errorTol Error tolerance for algebraic equations
	 */
	virtual void leanConsistentInitialState(double t, unsigned int secIdx, double timeFactor, double* const vecStateY, active* const adRes, active* const adY, unsigned int numSensAdDirs, double errorTol) = 0;

	/**
	 * @brief Computes approximately / partially consistent initial time derivatives
	 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] the initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ have
	 *          to be consistent. This functions overwrites the initial time derivatives \f$ \dot{y}_0 \f$ such that
	 *          the residual is zero. Thus, this function provides the final step in consistent initialization.
	 * 
	 *          This function is possibly faster than consistentInitialTimeDerivative(), but updates only a part of the
	 *          time derivative vector. Hence, the result is not guaranteed to be consistent.
	 * 
	 *          This function is to be used with leanConsistentInitialState(). Do not mix normal and lean
	 *          consistent initialization!
	 *          
	 * @param [in] t Current time point
	 * @param [in] timeFactor Used for time transformation (pre factor of time derivatives) and to compute parameter derivatives with respect to section length
	 * @param [in,out] vecStateYdot On entry, inconsistent state time derivatives. On exit, partially consistent state time derivatives.
	 * @param [in] res On entry, residual without taking time derivatives into account. The data is overwritten during execution of the function.
	 */
	virtual void leanConsistentInitialTimeDerivative(double t, double timeFactor, double* const vecStateYdot, double* const res) = 0;

	// Explicitly import virtual functions from base class in order to avoid warnings because of possible shadowing / hiding of this functions
	using ISimulatableModel::consistentIntialSensitivity;
	using ISimulatableModel::leanConsistentIntialSensitivity;

	/**
	 * @brief Computes consistent initial conditions for all sensitivity subsystems
	 * @details Given the DAE \f[ F(t, y, \dot{y}, p) = 0, \f] the corresponding (linear) forward sensitivity
	 *          system reads \f[ \frac{\partial F}{\partial y}(t, y, \dot{y}) s + \frac{\partial F}{\partial \dot{y}}(t, y, \dot{y}) \dot{s} + \frac{\partial F}{\partial p}(t, y, \dot{y}) = 0. \f]
	 *          The initial values of \f$ s_0 = \frac{\mathrm{d} y_0}{\mathrm{d}p} \f$ and \f$ \dot{s}_0 = \frac{\mathrm{d} \dot{y}_0}{\mathrm{d}p} \f$
	 *          have to be consistent, that means, they have to satisfy the sensitivity equation. This function computes the correct \f$ s_0 \f$ and \f$ \dot{s}_0 \f$
	 *          given \f$ y_0 \f$ and \f$ s_0 \f$.
	 * 
	 * @param [in] t Current time point
	 * @param [in] secIdx Index of the current section
	 * @param [in] timeFactor Used for time transformation (pre factor of time derivatives) and to compute parameter derivatives with respect to section length
	 * @param [in] vecStateY State vector with consistent initial values of the original system
	 * @param [in] vecStateYdot Time derivative state vector with consistent initial values of the original system
	 * @param [in,out] vecSensY Sensitivity subsystem state vectors
	 * @param [in,out] vecSensYdot Time derivative state vectors of the sensitivity subsystems to be initialized
	 * @param [in] adRes Pointer to global residual vector of AD datatypes with parameter sensitivities
	 */
	virtual void consistentIntialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
		std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes) = 0;

	/**
	 * @brief Computes approximately / partially consistent initial conditions for all sensitivity subsystems
	 * @details Given the DAE \f[ F(t, y, \dot{y}, p) = 0, \f] the corresponding (linear) forward sensitivity
	 *          system reads \f[ \frac{\partial F}{\partial y}(t, y, \dot{y}) s + \frac{\partial F}{\partial \dot{y}}(t, y, \dot{y}) \dot{s} + \frac{\partial F}{\partial p}(t, y, \dot{y}) = 0. \f]
	 *          The initial values of \f$ s_0 = \frac{\mathrm{d} y_0}{\mathrm{d}p} \f$ and \f$ \dot{s}_0 = \frac{\mathrm{d} \dot{y}_0}{\mathrm{d}p} \f$
	 *          have to be consistent, that means, they have to satisfy the sensitivity equation. This function computes the correct \f$ s_0 \f$ and \f$ \dot{s}_0 \f$
	 *          given \f$ y_0 \f$ and \f$ s_0 \f$.
	 *          
	 *          This function is possibly faster than consistentIntialSensitivity(), but updates only a part of the
	 *          vectors. Hence, the result is not guaranteed to be consistent. 
	 * 
	 * @param [in] t Current time point
	 * @param [in] secIdx Index of the current section
	 * @param [in] timeFactor Used for time transformation (pre factor of time derivatives) and to compute parameter derivatives with respect to section length
	 * @param [in] vecStateY State vector with consistent initial values of the original system
	 * @param [in] vecStateYdot Time derivative state vector with consistent initial values of the original system
	 * @param [in,out] vecSensY Sensitivity subsystem state vectors
	 * @param [in,out] vecSensYdot Time derivative state vectors of the sensitivity subsystems to be initialized
	 * @param [in] adRes Pointer to global residual vector of AD datatypes with parameter sensitivities
	 */
	virtual void leanConsistentIntialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
		std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes) = 0;

	/**
	 * @brief Sets external functions for this binding model
	 * @details The external functions are not owned by this IBindingModel.
	 * @param [in] extFuns Pointer to array of IExternalFunction objects of size @p size
	 * @param [in] size Number of elements in the IExternalFunction array @p extFuns
	 */
	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) = 0;
};

} // namespace cadet

#endif  // LIBCADET_IUNITOPERATION_HPP_
