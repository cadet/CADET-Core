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
 * Defines the BindingModel interface.
 */

#ifndef LIBCADET_BINDINGMODELINTERFACE_HPP_
#define LIBCADET_BINDINGMODELINTERFACE_HPP_

#include <unordered_map>

#include "cadet/ParameterProvider.hpp"
#include "cadet/ParameterId.hpp"
#include "linalg/DenseMatrix.hpp"
#include "linalg/BandMatrix.hpp"
#include "AutoDiff.hpp"

namespace cadet
{

class IParameterProvider;
class IExternalFunction;

namespace ad
{
	class IJacobianExtractor;
}

namespace model
{

/**
 * @brief Defines an internal BindingModel interface
 * @details The binding model is responsible for handling bound states and their residuals.
 *          There can be multiple bound states, if the binding model supports it.
 *          Some binding models require a dedicated salt component, which is always given
 *          by component 0.
 *          
 *          The ordering inside the solid phase can be chosen by the binding model to some
 *          extent. In general, the ordering is component-major, that means, all bound
 *          phases for each component are listed one after the other. For a model having
 *          4 components with 1, 0, 3, 2 bound states, respectively, the ordering is
 *          comp0bnd0, comp2bnd0, comp2bnd1, comp2bnd2, comp3bnd0, comp3bnd1.
 *          Note that the second component is not present in the solid phase at all since
 *          it has no bound state. The bound states of the other components are ordered
 *          consecutively.
 *          
 *          Algebraic equations have to appear in an isolated block (i.e., there must not
 *          appear differential equations between two algebraic equations).
 *          
 *          It is assumed that the liquid phase concentrations are located right before the
 *          first bound state of the first component, i.e., the general layout is expected
 *          to be 
 *              comp0, comp1, ..., compN, comp0bnd0, comp0bnd1, ..., comp1bnd0, comp1bnd1, ...
 *          where comp0, comp1, etc. denote liquid phase concentrations of the respective
 *          components. Pointers to state vectors usually point to comp0bnd0 instead of comp0,
 *          but this is made clear in the documentation.
 */
class IBindingModel
{
public:

	virtual ~IBindingModel() CADET_NOEXCEPT { }

	/**
	 * @brief Returns the name of the binding model
	 * @details This name is also used to identify and create the binding model in the factory.
	 * @return Name of the binding model
	 */
	virtual const char* name() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns whether the binding model requires additional parameters supplied by configure()
	 * @details After construction of an IBindingModel object, configureModelDiscretization() is called.
	 *          The binding model may require to read additional parameters from the adsorption group
	 *          of the parameter provider. This opportunity is given by a call to configure().
	 *          However, a binding model may not require this. This function communicates to the outside
	 *          whether a call to configure() is necessary.
	 * @return @c true if configure() has to be called, otherwise @c false
	 */
	virtual bool requiresConfiguration() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Sets the number of components and bound states in the model
	 * @details This function is called prior to configure() by the underlying model.
	 *          It can only be called once. Model parameters are configured by
	 *          configure().
	 * 
	 * @param [in] paramProvider Parameter provider
	 * @param [in] nComp Number of components
	 * @param [in] nBound Array of size @p nComp which contains the number of bound states for each component
	 * @param [in] boundOffset Array of size @p nComp with offsets to the first bound state of each component beginning from the solid phase
	 */
	virtual bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset) = 0;

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
	 * 
	 * @param [in] paramProvider Parameter provider
	 * @param [in] unitOpIdx Index of the unit operation this binding model belongs to
	 * @param [in] parTypeIdx Index of the particle type this binding model belongs to
	 * @return @c true if the configuration was successful, otherwise @c false
	 */
	virtual bool configure(IParameterProvider& paramProvider, unsigned int unitOpIdx, unsigned int parTypeIdx) = 0;

	/**
	 * @brief Returns the ParameterId of all bound phase initial conditions (equations)
	 * @details The array has to be filled in the order of the equations.
	 * @param [out] params Array with ParameterId objects to fill
	 * @param [in] unitOpIdx Index of the unit operation this binding model belongs to
	 * @param [in] parTypeIdx Index of the particle type this binding model belongs to
	 */
	virtual void fillBoundPhaseInitialParameters(ParameterId* params, unsigned int unitOpIdx, unsigned int parTypeIdx) const CADET_NOEXCEPT = 0;

	/**
	 * @brief Sets external functions for this binding model
	 * @details The external functions are not owned by this IBindingModel.
	 * 
	 * @param [in] extFuns Pointer to array of IExternalFunction objects of size @p size
	 * @param [in] size Number of elements in the IExternalFunction array @p extFuns
	 */
	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) = 0;

	/**
	 * @brief Checks whether a given parameter exists
	 * @param [in] pId ParameterId that identifies the parameter uniquely
	 * @return @c true if the parameter exists, otherwise @c false
	 */
	virtual bool hasParameter(const ParameterId& pId) const = 0;

	/**
	 * @brief Returns all parameters with their current values that can be made sensitive
	 * @return Map with all parameters that can be made sensitive along with their current value
	 */
	virtual std::unordered_map<ParameterId, double> getAllParameterValues() const = 0;

	/**
	 * @brief Sets a parameter value
	 * @details The parameter identified by its unique parameter is set to the given value.
	 * 
	 * @param [in] pId ParameterId that identifies the parameter uniquely
	 * @param [in] value Value of the parameter
	 * 
	 * @return @c true if the parameter has been successfully set to the given value,
	 *         otherwise @c false (e.g., if the parameter is not available in this model)
	 */
	virtual bool setParameter(const ParameterId& pId, int value) = 0;
	virtual bool setParameter(const ParameterId& pId, double value) = 0;
	virtual bool setParameter(const ParameterId& pId, bool value) = 0;

	/**
	 * @brief Returns the start index and length of the algebraic equations
	 * @details This function is only called if hasAlgebraicEquations() returns @c true.
	 *          It should return the local index (i.e., relative to the @c res pointer in residual())
	 *          of the first algebraic equation and the number of algebraic equations.
	 *          This is primarily used to consistently initialize the sensitivity state vectors.
	 * @param [out] idxStart Local index of the first algebraic equation
	 * @param [out] len Number of algebraic equations
	 */
	virtual void getAlgebraicBlock(unsigned int& idxStart, unsigned int& len) const = 0;

	/**
	 * @brief Returns a pointer to the parameter identified by the given id
	 * @param [in] pId Parameter Id of the sensitive parameter
	 * @return a pointer to the parameter if the binding model contains the parameter, otherwise @c nullptr
	 */
	virtual active* getParameter(const ParameterId& pId) = 0;

	/**
	 * @brief Returns whether this binding model has a salt component
	 * @details The salt component is given by component 0, if present.
	 * @return @c true if salt is required, otherwise @c false
	 */
	virtual bool hasSalt() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns whether this binding model supports multi-state binding
	 * @return @c true if multi-state binding is supported, otherwise @c false
	 */
	virtual bool supportsMultistate() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns whether this binding model supports non-binding components
	 * @details Non-binding components do not have an entry in the solid phase.
	 * @return @c true if non-binding components are supported, otherwise @c false
	 */
	virtual bool supportsNonBinding() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns whether this binding model has algebraic equations
	 * @return @c true if algebraic equations are present, otherwise @c false
	 */
	virtual bool hasAlgebraicEquations() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns whether this binding model depends on time
	 * @details Binding models may depend on time if external functions are used.
	 * @return @c true if the model is time-dependent, otherwise @c false
	 */
	virtual bool dependsOnTime() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns whether this binding model requires workspace
	 * @details The workspace may be required for consistent initialization and / or evaluation
	 *          of residual and Jacobian. A workspace is a memory buffer whose size is given by
	 *          workspaceSize().
	 * @return @c true if the model requires a workspace, otherwise @c false
	 */
	virtual bool requiresWorkspace() const CADET_NOEXCEPT = 0;

	/**
	 * @brief Returns the size of the required workspace in bytes
	 * @details The memory is required for externally dependent binding models and the 
	 *          nonlinear solver in consistent initialization.
	 * @param [in] nComp Number of components
	 * @param [in] totalNumBoundStates Total number of bound states
	 * @param [in] nBoundStates Array with bound states for each component
	 * @return Size of the workspace in bytes
	 */
	virtual unsigned int workspaceSize(unsigned int nComp, unsigned int totalNumBoundStates, unsigned int const* nBoundStates) const CADET_NOEXCEPT = 0;

	/**
	 * @brief Computes consistent initial state values
	 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] the initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ have
	 *          to be consistent. Finding consistent initial conditions is a two step process.
	 *          This functions performs the first step, that is, it updates the initial state \f$ y_0 \f$ overwriting
	 *          the algebraically determined state variables such that the algebraic equations hold.
	 *          
	 *          This function assumes that the liquid phase concentrations have been determined correctly.
	 *          The bound phase DOFs that are determined by algebraic equations (e.g., quasi-stationary
	 *          binding) are calculated by this functions using the liquid phase concentrations.
	 *          
	 *          This function is called simultaneously from multiple threads.
	 *          
	 *          If both @p adRes and @p adY are valid pointers, the Jacobian is to be computed by AD.
	 * 
	 * @param [in] t Current time point
	 * @param [in] z Axial position in normalized coordinates (column inlet = 0, column outlet = 1)
	 * @param [in] r Radial position in normalized coordinates (outer shell = 1, inner center = 0)
	 * @param [in] secIdx Index of the current section
	 * @param [in,out] vecStateY Pointer to first bound state in state vector with initial values 
	 *                 that are to be updated for consistency
	 * @param [in] errorTol Error tolerance for solving the algebraic equations
	 * @param [in,out] adRes Pointer to residual vector of AD datatypes that can be used to compute the Jacobian
	 * @param [in,out] adY Pointer to state vector of AD datatypes that can be used to compute the Jacobian
	 * @param [in] adEqOffset Offset of @p adRes and @p adY to the current cell
	 * @param [in] adDirOffset Offset to the usable AD directions
	 * @param [in] jacExtractor Jacobian extractor
	 * @param [in,out] workingMemory Working memory for nonlinear equation solvers
	 * @param [in,out] workingMat Working matrix for nonlinear equation solvers with at least as 
	 *                 many rows and columns as number of bound states
	 */
	virtual void consistentInitialState(double t, double z, double r, unsigned int secIdx, double* const vecStateY, double errorTol, active* const adRes, active* const adY,
		unsigned int adEqOffset, unsigned int adDirOffset, const ad::IJacobianExtractor& jacExtractor, double* const workingMemory,
		linalg::detail::DenseMatrixBase& workingMat) const = 0;

	/**
	 * @brief Evaluates the residual for one particle shell
	 * @details The binding model is responsible for implementing the complete bound state equations,
	 *          including adding the time derivatives to the residual. This function populates
	 *          the residuals from @c res[0] to @c res[nBound * nComp].
	 *          
	 *          This function is called simultaneously from multiple threads.
	 * 
	 * @param [in] t Current time point
	 * @param [in] z Axial position in normalized coordinates (column inlet = 0, column outlet = 1)
	 * @param [in] r Radial position in normalized coordinates (outer shell = 1, inner center = 0)
	 * @param [in] secIdx Index of the current section
	 * @param [in] timeFactor Used to compute parameter derivatives with respect to section length,
	 *             originates from time transformation and is premultiplied to time derivatives
	 * @param [in] y Pointer to first bound state of the first component in the current particle shell
	 * @param [in] yDot Pointer to first bound state time derivative of the first component in the current particle shell 
	 *             or @c nullptr if time derivatives shall be left out
	 * @param [out] res Pointer to residual equation of first bound state of the first component in the current particle shell
	 * @param [in,out] workSpace Memory work space
	 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
	 */
	virtual int residual(const active& t, double z, double r, unsigned int secIdx, const active& timeFactor, active const* y,
		double const* yDot, active* res, void* workSpace) const = 0;

	virtual int residual(double t, double z, double r, unsigned int secIdx, double timeFactor, active const* y,
		double const* yDot, active* res, void* workSpace) const = 0;

	virtual int residual(const active& t, double z, double r, unsigned int secIdx, const active& timeFactor, double const* y,
		double const* yDot, active* res, void* workSpace) const = 0;

	virtual int residual(double t, double z, double r, unsigned int secIdx, double timeFactor, double const* y,
		double const* yDot, double* res, void* workSpace) const = 0;

	/**
	 * @brief Evaluates the Jacobian of the bound states for one particle shell analytically
	 * @details The binding model is responsible for implementing the complete bound state equations,
	 *          including adding the time derivatives to the residual.
	 * 
	 *          This function is called simultaneously from multiple threads.
	 *          It can be left out (empty implementation) if AD is used to evaluate the Jacobian.
	 *
	 * @param [in] t Current time point
	 * @param [in] z Axial position in normalized coordinates (column inlet = 0, column outlet = 1)
	 * @param [in] r Radial position in normalized coordinates (outer shell = 1, inner center = 0)
	 * @param [in] secIdx Index of the current section
	 * @param [in] y Pointer to first bound state of the first component in the current particle shell
	 * @param [in,out] jac Row iterator pointing to the first bound states row of the underlying BandMatrix in which the Jacobian is stored
	 * @param [in,out] workSpace Memory work space
	 */
	virtual void analyticJacobian(double t, double z, double r, unsigned int secIdx, double const* y, linalg::BandMatrix::RowIterator jac, void* workSpace) const = 0;
	virtual void analyticJacobian(double t, double z, double r, unsigned int secIdx, double const* y, linalg::DenseBandedRowIterator jac, void* workSpace) const = 0;

	/**
	 * @brief Adds the time-discretized part of the Jacobian to the current Jacobian of the bound phase equations in one particle shell
	 * @details The added time derivatives in jacobian() have to be added to the Jacobian of the original equations in order to get
	 *          the time-discretized equation's Jacobian.
	 * 
	 *          This function is called simultaneously from multiple threads.
	 *          It is always required, even if the Jacobian is computed using AD.
	 *
	 * @param [in] alpha Factor of the time derivatives which comes from a BDF discretization and the time transformation
	 * @param [in,out] jac Row iterator pointing to the first bound states row of the underlying BandMatrix in which the Jacobian is stored
	 */
	virtual void jacobianAddDiscretized(double alpha, linalg::FactorizableBandMatrix::RowIterator jac) const = 0;
	virtual void jacobianAddDiscretized(double alpha, linalg::DenseBandedRowIterator jac) const = 0;

	/**
	 * @brief Multiplies the Jacobian of the model with respect to @f$ \dot{y} @f$ with the given vector @p yDotS
	 * @details The model is a function @f$ f(t, y, \dot{y}, p) @f$. This function computes the matrix-vector product
	 *          @f[ \begin{align*} \frac{\partial f}{\partial \dot{y}}(t, y \dot{y}, p) \dot{y}_s. \end{align*} @f]
	 *          This function is always required, even if the Jacobian is computed using AD.
	 * @param [in] yDotS State vector pointer to first bound state of the first component in the current particle shell
	 * @param [out] res Pointer to vector in which the result of the matrix-vector product is stored
	 * @param [in] timeFactor Factor of the time derivatives that comes from time transformation
	 */
	virtual void multiplyWithDerivativeJacobian(double const* yDotS, double* const res, double timeFactor) const = 0;

	/**
	 * @brief Calculates the time derivative of the algebraic residual equations
	 * @details Calculates @f$ \frac{\partial \text{res}_{\text{alg}}}{\partial t} @f$ for the algebraic equations
	 *          in the residual.
	 *          
	 *          This function is called simultaneously from multiple threads.
	 *          It can be left out (empty implementation) which leads to slightly incorrect initial conditions
	 *          when using externally dependent binding models.
	 *          
	 * @param [in] t Current time point
	 * @param [in] z Axial position in normalized coordinates (column inlet = 0, column outlet = 1)
	 * @param [in] r Radial position in normalized coordinates (outer shell = 1, inner center = 0)
	 * @param [in] secIdx Index of the current section
	 * @param [in] y Pointer to first bound state of the first component in the current particle shell
	 * @param [out] dResDt Pointer to array that stores the time derivative
	 * @param [in,out] workSpace Memory work space
	 */
	virtual void timeDerivativeAlgebraicResidual(double t, double z, double r, unsigned int secIdx, double const* y, double* dResDt, void* workSpace) const = 0;
protected:
};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_BINDINGMODELINTERFACE_HPP_
