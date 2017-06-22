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
 * Defines a BindingModelBase class.
 */

#ifndef LIBCADET_BINDINGMODELBASE_HPP_
#define LIBCADET_BINDINGMODELBASE_HPP_

#include "model/BindingModel.hpp"
#include "model/binding/BindingModelMacros.hpp"
#include "ParamIdUtil.hpp"

#include <vector>
#include <unordered_map>

namespace cadet
{

namespace nonlin
{
	class Solver;
}

namespace model
{

/**
 * @brief Defines a BindingModel base class that can be used to implement other binding models
 * @details This base class can be used as a starting point for new binding models. It assumes
 *          a binding model without salt and multi-state support. Non-binding components, however,
 *          are supported. The model is expected to support kinetic and quasi-stationary binding
 *          modes.
 *          
 *          Some common parameter handling is provided using a hash map (std::unordered_map).
 *          Furthermore, a configurable nonlinear solver is managed which can be used for
 *          consistent initialization.
 */
class BindingModelBase : public IBindingModel
{
public:

	BindingModelBase();

	virtual ~BindingModelBase() CADET_NOEXCEPT;

	virtual bool configure(IParameterProvider& paramProvider, unsigned int unitOpIdx);
	virtual bool reconfigure(IParameterProvider& paramProvider, unsigned int unitOpIdx);
	virtual void configureModelDiscretization(unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset);

	virtual std::unordered_map<ParameterId, double> getAllParameterValues() const;
	virtual bool hasParameter(const ParameterId& pId) const;

	virtual bool setParameter(const ParameterId& pId, int value);
	virtual bool setParameter(const ParameterId& pId, double value);
	virtual bool setParameter(const ParameterId& pId, bool value);

	virtual active* getParameter(const ParameterId& pId);

	virtual bool hasSalt() const CADET_NOEXCEPT { return false; }
	virtual bool supportsMultistate() const CADET_NOEXCEPT { return false; }
	virtual bool supportsNonBinding() const CADET_NOEXCEPT { return true; }

	virtual unsigned int consistentInitializationWorkspaceSize() const;

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { }

	virtual void timeDerivativeAlgebraicResidual(double t, double z, double r, unsigned int secIdx, double const* y, double* dResDt) const { }
protected:
	int _nComp; //!< Number of components
	unsigned int const* _nBoundStates; //!< Array with number of bound states for each component
	bool _kineticBinding; //!< Determines whether binding is kinetic (@c true) or quasi-stationary (@c false)

	std::unordered_map<ParameterId, active*> _parameters; //!< Map used to translate ParameterIds to actual variables
	nonlin::Solver* _nonlinearSolver; //!< Nonlinear equation solver for consistent initialization

	/**
	 * @brief Configures the binding model
	 * @details This function implements the (re-)configuration of a binding model. It is called when
	 *          the binding model is configured or reconfigured. On call the _parameters map will always
	 *          be empty and _kineticBinding is already configured. If @p reconfigure is @c true, then
	 *          the binding model has been configured before and the member variables should have enough
	 *          memory reserved.
	 * @param [in] reconfigure @c true if reconfigure() was called, @c false if configure() was called
	 * @param [in] paramProvider Parameter provider
	 * @param [in] unitOpIdx Unit operation index
	 * @return @c true if the configuration was successful, otherwise @c false
	 */
	virtual bool configureImpl(bool reconfigure, IParameterProvider& paramProvider, unsigned int unitOpIdx) = 0;

	/**
	 * @brief Configures the nonlinear solver which is used for consistent initialization
	 * @param [in] paramProvider Parameter provider
	 * @return @c true if the configuration was successful, otherwise @c false
	 */
	bool configureNonlinearSolver(IParameterProvider& paramProvider);
};

/**
 * @brief Defines a BindingModel base class in which either all or none of the equations are kinetic
 * @details This base adds more functionality to BindingModelBase. It assumes that either all or
 *          none of the equations are kinetic, whereas the normal BindingModelBase does not make
 *          such assumptions and also allows mixed equations.
 *          
 *          As a benefit, this class provides consistent initialization with and without AD (and
 *          Jacobian checking) and handles the Jacobian operations with respect to @f$ \dot{y} @f$.
*/
class PureBindingModelBase : public BindingModelBase
{
public:

	PureBindingModelBase();

	virtual ~PureBindingModelBase() CADET_NOEXCEPT;

	virtual bool hasAlgebraicEquations() const CADET_NOEXCEPT { return !_kineticBinding; }
	virtual void getAlgebraicBlock(unsigned int& idxStart, unsigned int& len) const;

	virtual void consistentInitialState(double t, double z, double r, unsigned int secIdx, double* const vecStateY, double errorTol, 
		active* const adRes, active* const adY, unsigned int adEqOffset, unsigned int adOffset, unsigned int diagDir, 
		unsigned int lowerBandwidth, unsigned int upperBandwidth, double* const workingMemory, linalg::detail::DenseMatrixBase& workingMat) const;

	virtual void analyticJacobian(double t, double z, double r, unsigned int secIdx, double const* y, linalg::BandMatrix::RowIterator jac) const;
	virtual void jacobianAddDiscretized(double alpha, linalg::FactorizableBandMatrix::RowIterator jac) const;
	virtual void multiplyWithDerivativeJacobian(double const* yDotS, double* const res, double timeFactor) const;

protected:

	/**
	 * @brief Computes the residual
	 * @details This function is similar to residual() and should effectively call the
	 *          templatized residualImpl() function. Please refer to IBindingModel::residual()
	 *          for more details. Note that the same assumptions made there apply here.
	 * @param [in] t Current time
	 * @param [in] z Axial position in normalized coordinates (column inlet = 0, column outlet = 1)
	 * @param [in] r Radial position in normalized coordinates (outer shell = 1, inner center = 0)
	 * @param [in] secIdx Index of the current section
	 * @param [in] timeFactor Used to compute parameter derivatives with respect to section length (nominal value should always be 1.0)
	 * @param [in] y Pointer to first bound state of the first component in the current particle shell
	 * @param [in] yCp Pointer to first component in bead liquid phase of the current particle shell
	 * @param [in] yDot Pointer to first bound state time derivative of the first component in the current particle shell 
	 *             or @c nullptr if time derivatives shall be left out
	 * @param [out] res Pointer to residual equation of first bound state of the first component in the current particle shell
	 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
	 */
	virtual int residualCore(double t, double z, double r, unsigned int secIdx, double timeFactor,
		double const* y, double const* yCp, double const* yDot, double* res) const = 0;
	virtual int residualCore(double t, double z, double r, unsigned int secIdx, double timeFactor,
		active const* y, double const* yCp, double const* yDot, active* res) const = 0;

	/**
	 * @brief Evaluates the Jacobian of the bound states for one particle shell analytically
	 * @details This function is similar to analyticJacobian() and should effectively call the
	 *          templatized jacobianImpl() function. Please refer to IBindingModel::analyticJacobian()
	 *          for more details. Note that the same assumptions made there apply here.
	 * @param [in] t Current time point
	 * @param [in] z Axial position in normalized coordinates (column inlet = 0, column outlet = 1)
	 * @param [in] r Radial position in normalized coordinates (outer shell = 1, inner center = 0)
	 * @param [in] secIdx Index of the current section
	 * @param [in] y Pointer to first bound state of the first component in the current particle shell
	 * @param [in] yCp Pointer to first component in bead liquid phase of the current particle shell
	 * @param [in,out] jac Row iterator pointing to the first bound states row of the underlying BandMatrix in which the Jacobian is stored
	 */
	virtual void analyticJacobianCore(double t, double z, double r, unsigned int secIdx, double const* y, 
		double const* yCp, linalg::BandMatrix::RowIterator jac) const = 0;
	virtual void analyticJacobianCore(double t, double z, double r, unsigned int secIdx, double const* y, 
		double const* yCp, linalg::detail::DenseMatrixBase::RowIterator jac) const = 0;
};


/**
 * @brief Inserts implementations for common functions required by PureBindingModelBase forwarding them to templatized functions
 * @details An implementation of PureBindingModelBase has to provide some protected virtual functions.
 *          This macro provides the implementation of those functions by forwarding them to the templatized 
 *          functions residualImpl() and jacobianImpl() which are assumed to be present in the class.
 * 
 * @param CLASSNAME Name of the PureBindingModelBase heir (including template)
 * @param TEMPLATELINE Line before each function that may contain a template<typename TEMPLATENAME> modifier
 */
#define CADET_PUREBINDINGMODELBASE_BOILERPLATE_IMPL_BASE(CLASSNAME, TEMPLATELINE)                            \
	CADET_BINDINGMODEL_RESIDUAL_BOILERPLATE_IMPL_BASE(CLASSNAME, TEMPLATELINE)                               \
	TEMPLATELINE                                                                                             \
	int CLASSNAME::residualCore(double t, double z, double r, unsigned int secIdx, double timeFactor,        \
		double const* y, double const* yCp, double const* yDot, double* res) const                           \
	{                                                                                                        \
		return residualImpl<double, double, double, double>(t, z, r, secIdx, timeFactor, y, yCp, yDot, res); \
	}                                                                                                        \
	                                                                                                         \
	TEMPLATELINE                                                                                             \
	int CLASSNAME::residualCore(double t, double z, double r, unsigned int secIdx, double timeFactor,        \
		active const* y, double const* yCp, double const* yDot, active* res) const                           \
	{                                                                                                        \
		return residualImpl<active, double, active, double>(t, z, r, secIdx, timeFactor, y, yCp, yDot, res); \
	}                                                                                                        \
	                                                                                                         \
	TEMPLATELINE                                                                                             \
	void CLASSNAME::analyticJacobianCore(double t, double z, double r, unsigned int secIdx, double const* y, \
		double const* yCp, linalg::BandMatrix::RowIterator jac) const                                        \
	{                                                                                                        \
		jacobianImpl(t, z, r, secIdx, y, yCp, jac);                                                          \
	}                                                                                                        \
	                                                                                                         \
	TEMPLATELINE                                                                                             \
	void CLASSNAME::analyticJacobianCore(double t, double z, double r, unsigned int secIdx, double const* y, \
		double const* yCp, linalg::detail::DenseMatrixBase::RowIterator jac) const                           \
	{                                                                                                        \
		jacobianImpl(t, z, r, secIdx, y, yCp, jac);                                                          \
	}


/**
 * @brief Inserts implementations for common functions required by PureBindingModelBase forwarding them to templatized functions
 * @details An implementation of PureBindingModelBase has to provide some protected virtual functions.
 *          This macro provides the implementation of those functions by forwarding them to the templatized 
 *          functions residualImpl() and jacobianImpl() which are assumed to be present in the class.
 * 
 * @param CLASSNAME Name of the PureBindingModelBase heir
 */
#define CADET_PUREBINDINGMODELBASE_BOILERPLATE_IMPL(CLASSNAME)       \
	CADET_PUREBINDINGMODELBASE_BOILERPLATE_IMPL_BASE(CLASSNAME,)

/**
 * @brief Inserts implementations for common functions required by PureBindingModelBase forwarding them to templatized functions
 * @details An implementation of PureBindingModelBase has to provide some protected virtual functions.
 *          This macro provides the implementation of those functions by forwarding them to the templatized 
 *          functions residualImpl() and jacobianImpl() which are assumed to be present in the class.
 * 
 * @param CLASSNAME Name of the PureBindingModelBase heir
 * @param TEMPLATENAME Name of the template parameter that handles externally dependent binding models
 */
#define CADET_PUREBINDINGMODELBASE_TEMPLATED_BOILERPLATE_IMPL(CLASSNAME, TEMPLATENAME)                             \
	CADET_PUREBINDINGMODELBASE_BOILERPLATE_IMPL_BASE(CLASSNAME<TEMPLATENAME>, template <typename TEMPLATENAME>)

} // namespace model
} // namespace cadet

#endif  // LIBCADET_BINDINGMODELBASE_HPP_
