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
 * Defines the convection dispersion transport operator.
 */

#ifndef LIBCADET_CONVECTIONDISPERSIONOPERATOR_HPP_
#define LIBCADET_CONVECTIONDISPERSIONOPERATOR_HPP_

#include "ParamIdUtil.hpp"
#include "AutoDiff.hpp"
#include "linalg/BandMatrix.hpp"
#include "MemoryPool.hpp"
#include "Weno.hpp"

#include <unordered_map>
#include <vector>

namespace cadet
{

class IParameterProvider;

namespace model
{

namespace operators
{

/**
 * @brief Convection dispersion transport operator
 * @details Implements the equation
 * 
 * @f[\begin{align}
	\frac{\partial c_i}{\partial t} &= - u \frac{\partial c_i}{\partial z} + D_{\text{ax}} \frac{\partial^2 c_i}{\partial z^2} \\
\end{align} @f]
 * with Danckwerts boundary conditions (see @cite Danckwerts1953)
@f[ \begin{align}
u c_{\text{in},i}(t) &= u c_i(t,0) - D_{\text{ax}} \frac{\partial c_i}{\partial z}(t,0) \\
\frac{\partial c_i}{\partial z}(t,L) &= 0
\end{align} @f]
 * Methods are described in @cite VonLieres2010a (WENO, linear solver), @cite Puttmann2013 @cite Puttmann2016 (forward sensitivities, AD, band compression)
 */
class ConvectionDispersionOperator
{
public:

	ConvectionDispersionOperator();
	~ConvectionDispersionOperator() CADET_NOEXCEPT;

	unsigned int requiredADdirs() const CADET_NOEXCEPT;

	void setFlowRates(const active& in, const active& out, const active& colPorosity) CADET_NOEXCEPT;

	bool configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, unsigned int nComp, unsigned int nCol);
	bool reconfigure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters);
	bool notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active* const adRes, active* const adY, unsigned int adDirOffset);

	int residual(double t, unsigned int secIdx, double timeFactor, double const* y, double const* yDot, double* res, bool wantJac);
	int residual(double t, unsigned int secIdx, double timeFactor, active const* y, double const* yDot, active* res, bool wantJac);
	int residual(const active& t, unsigned int secIdx, const active& timeFactor, double const* y, double const* yDot, active* res, bool wantJac);
	int residual(const active& t, unsigned int secIdx, const active& timeFactor, active const* y, double const* yDot, active* res, bool wantJac);

	void prepareADvectors(active* const adRes, active* const adY, unsigned int adDirOffset) const;
	void extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset);

	bool solveTimeDerivativeSystem(double t, unsigned int secIdx, double timeFactor, double* const rhs);
	void multiplyWithDerivativeJacobian(double t, unsigned int secIdx, double timeFactor, double const* sDot, double* ret) const;

	void addTimeDerivativeToJacobian(double alpha, double timeFactor);
	void assembleDiscretizedJacobian(double alpha, double timeFactor);
	bool assembleAndFactorizeDiscretizedJacobian(double alpha, double timeFactor);
	bool solveDiscretizedJacobian(double* rhs) const;

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	double checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const;
#endif

	const active& columnLength() const CADET_NOEXCEPT { return _colLength; }
	const active& crossSectionArea() const CADET_NOEXCEPT { return _crossSection; }
	const active& currentVelocity() const CADET_NOEXCEPT { return _curVelocity; }

	linalg::BandMatrix& jacobian() CADET_NOEXCEPT { return _jacC; }
	const linalg::BandMatrix& jacobian() const CADET_NOEXCEPT { return _jacC; }

	linalg::FactorizableBandMatrix& jacobianDisc() CADET_NOEXCEPT { return _jacCdisc; }
	const linalg::FactorizableBandMatrix& jacobianDisc() const CADET_NOEXCEPT { return _jacCdisc; }
protected:

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	int residualImpl(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* y, double const* yDot, ResidualType* res);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	int residualForwardsFlow(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* y, double const* yDot, ResidualType* res);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	int residualBackwardsFlow(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* y, double const* yDot, ResidualType* res);

	unsigned int _nComp; //!< Number of components
	unsigned int _nCol; //!< Number of axial cells

	linalg::BandMatrix _jacC; //!< Jacobian
	linalg::FactorizableBandMatrix _jacCdisc; //!< Jacobian with time derivatives from BDF method

	active _colLength; //!< Column length \f$ L \f$
	active _crossSection; //!< Cross section area 

	// Section dependent parameters
	std::vector<active> _colDispersion; //!< Column dispersion (may be section dependent) \f$ D_{\text{ax}} \f$
	std::vector<active> _velocity; //!< Interstitial velocity (may be section dependent) \f$ u \f$
	active _curVelocity; //!< Current interstitial velocity \f$ u \f$ in this time section

	ArrayPool _stencilMemory; //!< Provides memory for the stencil
	double* _wenoDerivatives; //!< Holds derivatives of the WENO scheme
	Weno _weno; //!< The WENO scheme implementation
	double _wenoEpsilon; //!< The @f$ \varepsilon @f$ of the WENO scheme (prevents division by zero)

	// Indexer functionality

	// Strides
	inline const int strideColCell() const CADET_NOEXCEPT { return static_cast<int>(_nComp); }
	inline const int strideColComp() const CADET_NOEXCEPT { return 1; }

	// Offsets
	inline const int offsetC() const CADET_NOEXCEPT { return _nComp; }
	inline const int offsetCp() const CADET_NOEXCEPT { return _nComp * _nCol + offsetC(); }

	// Return pointer to first element of state variable in state vector
	template <typename real_t> inline real_t* c(real_t* const data) const { return data + offsetC(); }
	template <typename real_t> inline real_t const* c(real_t const* const data) const { return data + offsetC(); }

	// Return specific variable in state vector
	template <typename real_t> inline real_t& c(real_t* const data, unsigned int col, unsigned int comp) const { return data[offsetC() + comp + col * strideColCell()]; }
	template <typename real_t> inline const real_t& c(real_t const* const data, unsigned int col, unsigned int comp) const { return data[offsetC() + comp + col * strideColCell()]; }
};

} // namespace operators
} // namespace model
} // namespace cadet

#endif  // LIBCADET_CONVECTIONDISPERSIONOPERATOR_HPP_
