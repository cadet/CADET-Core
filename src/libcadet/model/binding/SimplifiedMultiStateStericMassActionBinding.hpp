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
 * Defines the SimplifiedMultiStateStericMassActionBinding class.
 */

#ifndef LIBCADET_SIMPLIFIEDMULTISTATESTERICMASSACTIONBINDING_HPP_
#define LIBCADET_SIMPLIFIEDMULTISTATESTERICMASSACTIONBINDING_HPP_

#include "model/binding/BindingModelBase.hpp"
#include "model/binding/BindingModelMacros.hpp"
#include "SlicedVector.hpp"

#include <vector>

namespace cadet
{

namespace model
{

/**
 * @brief Defines the simplified Multi-State-SMA binding model
 * @details Implements a simplification of the Multi-State-SMA adsorption model defined in cadet::model::MultiStateStericMassActionBinding.
 *          The simplification consists in additional assumptions and a different parameterization which significantly 
 *          reduces the number of parameters. The assumptions are
 *              -# Molecules are only exchanged between two adjacent states, this means no transfer from state @f$ q_i^{(1)} @f$ to state @f$ q_i^{(3)} @f$.
 *              -# Characteristic charge @f$ \nu_i^{(j)} @f$ and shielding factor @f$ \sigma_i^{(j)} @f$ only depend on the index of the state @f$ j @f$.
 *
 *          Thus, the exchange parameters @f$ k^{(i)}_{j\ell} @f$, the characteristic charge @f$ \nu_i^{(j)} @f$, and the shielding @f$ \sigma_i^{(j)} @f$
 *          can be parameterized with few degrees of freedom. For all @f$ i = 1,\dots,N_{\text{comp}} @f$ and @f$ j,\ell = 1,\dots,M_i @f$ let
 *          @f[ \begin{align*}
 *                  k^{(i)}_{j\ell} &= \begin{cases}
 *                      0, & \text{for } \abs{j-\ell} \neq 1 \\
 *                      K^{(i)}_{ws} + j K^{(i)}_{ws,\text{lin}} + K^{(i)}_{ws,\text{quad}} (j-1)(j - M+1), & \text{for } \ell = j+1 \\
 *                      K^{(i)}_{sw} + \ell K^{(i)}_{sw,\text{lin}} + K^{(i)}_{sw,\text{quad}} (\ell-1)(\ell - M+1), & \text{for } \ell = j-1,
 *                  \end{cases}\\
 *                  \nu_i^{(j)} &= \nu_{\text{min},i} + \frac{j-1}{M-1} \left( \nu_{\text{max},i} - \nu_{\text{min},i} \right) + \nu_{\text{quad},i} (j-1) (j-M), \\
 *                  \sigma_i^{(j)} &= \sigma_{\text{min},i} + \frac{j-1}{M-1} \left( \sigma_{\text{max},i} - \sigma_{\text{min},i} \right) + \sigma_{\text{quad},i} (j-1) (j-M).
 *          \end{align*} @f]
 *          The parameterization is orthogonal in the sense that the endpoints of the linear interpolation defined by @f$ \nu_{\text{min},i} @f$
 *          and @f$ \nu_{\text{max},i} @f$ are not changed by @f$ \nu_{\text{quad},i} @f$. Rather the quadratic effect is added on top of
 *          the linear interpolation.
 *          
 *          The internal state vector layout follows the more general MultiStateStericMassActinoBinding, which is component-major.
 *          The state vector contains all components and within each component all bound states are placed.
 */
class SimplifiedMultiStateStericMassActionBinding : public BindingModelBase
{
public:

	SimplifiedMultiStateStericMassActionBinding();
	virtual ~SimplifiedMultiStateStericMassActionBinding() CADET_NOEXCEPT;

	static const char* identifier() { return "SIMPLIFIED_MULTISTATE_STERIC_MASS_ACTION"; }
	virtual const char* name() const CADET_NOEXCEPT { return SimplifiedMultiStateStericMassActionBinding::identifier(); }

	virtual void configureModelDiscretization(unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset);

	virtual void getAlgebraicBlock(unsigned int& idxStart, unsigned int& len) const;

	virtual unsigned int consistentInitializationWorkspaceSize() const;
	virtual void consistentInitialState(double t, double z, double r, unsigned int secIdx, double* const vecStateY, double errorTol, 
		active* const adRes, active* const adY, unsigned int adEqOffset, unsigned int adDirOffset, const ad::IJacobianExtractor& jacExtractor, 
		double* const workingMemory, linalg::detail::DenseMatrixBase& workingMat) const;

	CADET_BINDINGMODEL_RESIDUAL_JACOBIAN_BOILERPLATE

	virtual void multiplyWithDerivativeJacobian(double const* yDotS, double* const res, double timeFactor) const;

	virtual bool hasSalt() const CADET_NOEXCEPT { return true; }
	virtual bool supportsMultistate() const CADET_NOEXCEPT { return true; }
	virtual bool supportsNonBinding() const CADET_NOEXCEPT { return true; }
	virtual bool hasAlgebraicEquations() const CADET_NOEXCEPT { return true; }

	virtual bool dependsOnTime() const CADET_NOEXCEPT { return false; }

protected:
	active _lambda; //! Ionic capacity
	util::SlicedVector<active> _kA; //!< Adsorption rate in component-major ordering
	util::SlicedVector<active> _kD; //!< Desorption rate in component-major ordering
	std::vector<active> _nuMin; //!< Characteristic charge, minimum value
	std::vector<active> _nuMax; //!< Characteristic charge, maximum value
	std::vector<active> _nuQuad; //!< Characteristic charge, quadratic modifier
	std::vector<active> _sigmaMin; //!< Steric factor, minimum value
	std::vector<active> _sigmaMax; //!< Steric factor, maximum value
	std::vector<active> _sigmaQuad; //!< Steric factor, quadratic modifier
	std::vector<active> _kSW; //!< State transition rate from strong to weak state, offset
	std::vector<active> _kSW_lin; //!< State transition rate from strong to weak state, linear factor
	std::vector<active> _kSW_quad; //!< State transition rate from strong to weak state, quadratic factor
	std::vector<active> _kWS; //!< State transition rate from weak to strong state, offset
	std::vector<active> _kWS_lin; //!< State transition rate from weak to strong state, linear factor
	std::vector<active> _kWS_quad; //!< State transition rate from weak to strong state, quadratic factor
	active _refC0; //! Liquid phase reference concentration
	active _refQ; //! Solid phase reference concentration

	virtual bool configureImpl(bool reconfigure, IParameterProvider& paramProvider, unsigned int unitOpIdx);
	
	template <typename ParamType>
	inline ParamType sigma(int comp, double state) const;

	template <typename ParamType>
	inline ParamType nu(int comp, double state) const;

	template <typename ParamType>
	inline ParamType k_sw(int comp, double state) const;

	template <typename ParamType>
	inline ParamType k_ws(int comp, double state) const;

	template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
	int residualImpl(const ParamType& t, double z, double r, unsigned int secIdx, const ParamType& timeFactor,
		StateType const* y, CpStateType const* yCp, double const* yDot, ResidualType* res) const;

	template <typename RowIterator>
	void jacobianImpl(double t, double z, double r, unsigned int secIdx, double const* y, double const* yCp, RowIterator jac) const;

	template <typename RowIterator>
	void jacobianAddDiscretizedImpl(double alpha, RowIterator jac) const;
};


} // namespace model

} // namespace cadet

#endif // LIBCADET_SIMPLIFIEDMULTISTATESTERICMASSACTIONBINDING_HPP_
