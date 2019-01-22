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
struct AdJacobianParams;
struct SimulationTime;

namespace model
{

namespace parts
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
 * Methods are described in @cite VonLieres2010a (WENO, linear solver), and @cite Puttmann2013, @cite Puttmann2016 (forward sensitivities, AD, band compression)
 * 
 * This class does not store the Jacobian. It only fills existing matrices given to its residual() functions.
 */
class ConvectionDispersionOperatorBase
{
public:

	ConvectionDispersionOperatorBase();
	~ConvectionDispersionOperatorBase() CADET_NOEXCEPT;

	void setFlowRates(const active& in, const active& out, const active& colPorosity) CADET_NOEXCEPT;

	bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int nCol, unsigned int strideCell);
	bool configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters);
	bool notifyDiscontinuousSectionTransition(double t, unsigned int secIdx);

	int residual(double t, unsigned int secIdx, double timeFactor, double const* y, double const* yDot, double* res, linalg::BandMatrix* jac);
	int residual(double t, unsigned int secIdx, double timeFactor, active const* y, double const* yDot, active* res);
	int residual(const active& t, unsigned int secIdx, const active& timeFactor, double const* y, double const* yDot, active* res, linalg::BandMatrix* jac);
	int residual(const active& t, unsigned int secIdx, const active& timeFactor, active const* y, double const* yDot, active* res);

	void multiplyWithDerivativeJacobian(const SimulationTime& simTime, double const* sDot, double* ret) const;
	void addTimeDerivativeToJacobian(double alpha, double timeFactor, linalg::FactorizableBandMatrix& jacDisc);

	inline const active& columnLength() const CADET_NOEXCEPT { return _colLength; }
	inline const active& crossSectionArea() const CADET_NOEXCEPT { return _crossSection; }
	inline const active& currentVelocity() const CADET_NOEXCEPT { return _curVelocity; }
	double previousVelocity(unsigned int secIdx) const CADET_NOEXCEPT;

	inline unsigned int nComp() const CADET_NOEXCEPT { return _nComp; }
	inline unsigned int nCol() const CADET_NOEXCEPT { return _nCol; }
	inline const Weno& weno() const CADET_NOEXCEPT { return _weno; }

	unsigned int jacobianLowerBandwidth() const CADET_NOEXCEPT;
	unsigned int jacobianUpperBandwidth() const CADET_NOEXCEPT;
	unsigned int jacobianDiscretizedBandwidth() const CADET_NOEXCEPT;

protected:

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	int residualImpl(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* y, double const* yDot, ResidualType* res, linalg::BandMatrix* jac);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	int residualForwardsFlow(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* y, double const* yDot, ResidualType* res, linalg::BandMatrix* jacMat);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	int residualBackwardsFlow(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* y, double const* yDot, ResidualType* res, linalg::BandMatrix* jacMat);

	unsigned int _nComp; //!< Number of components
	unsigned int _nCol; //!< Number of axial cells
	unsigned int _strideCell; //!< Number of elements between the same item in two adjacent cells

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
	inline const int strideColCell() const CADET_NOEXCEPT { return static_cast<int>(_strideCell); }
	inline const int strideColComp() const CADET_NOEXCEPT { return 1; }

	// Offsets
	inline const int offsetC() const CADET_NOEXCEPT { return _nComp; }
	inline const int offsetCp() const CADET_NOEXCEPT { return strideColCell() * _nCol + offsetC(); }

	// Return pointer to first element of state variable in state vector
	template <typename real_t> inline real_t* c(real_t* const data) const { return data + offsetC(); }
	template <typename real_t> inline real_t const* c(real_t const* const data) const { return data + offsetC(); }

	// Return specific variable in state vector
	template <typename real_t> inline real_t& c(real_t* const data, unsigned int col, unsigned int comp) const { return data[offsetC() + comp + col * strideColCell()]; }
	template <typename real_t> inline const real_t& c(real_t const* const data, unsigned int col, unsigned int comp) const { return data[offsetC() + comp + col * strideColCell()]; }
};


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
 * Methods are described in @cite VonLieres2010a (WENO, linear solver), and @cite Puttmann2013, @cite Puttmann2016 (forward sensitivities, AD, band compression)
 * 
 * This class wraps ConvectionDispersionOperatorBase and provides all the functionality it does. In addition,
 * the Jacobian is stored and corresponding functions are provided (assembly, factorization, solution, retrieval).
 */
class ConvectionDispersionOperator
{
public:

	ConvectionDispersionOperator();
	~ConvectionDispersionOperator() CADET_NOEXCEPT;

	unsigned int requiredADdirs() const CADET_NOEXCEPT;

	void setFlowRates(const active& in, const active& out, const active& colPorosity) CADET_NOEXCEPT;

	bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int nCol);
	bool configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters);
	bool notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, const AdJacobianParams& adJac);

	int residual(double t, unsigned int secIdx, double timeFactor, double const* y, double const* yDot, double* res, bool wantJac);
	int residual(double t, unsigned int secIdx, double timeFactor, active const* y, double const* yDot, active* res, bool wantJac);
	int residual(const active& t, unsigned int secIdx, const active& timeFactor, double const* y, double const* yDot, active* res, bool wantJac);
	int residual(const active& t, unsigned int secIdx, const active& timeFactor, active const* y, double const* yDot, active* res, bool wantJac);

	void prepareADvectors(const AdJacobianParams& adJac) const;
	void extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset);

	bool solveTimeDerivativeSystem(const SimulationTime& simTime, double* const rhs);
	void multiplyWithDerivativeJacobian(const SimulationTime& simTime, double const* sDot, double* ret) const;

	void addTimeDerivativeToJacobian(double alpha, double timeFactor);
	void assembleDiscretizedJacobian(double alpha, double timeFactor);
	bool assembleAndFactorizeDiscretizedJacobian(double alpha, double timeFactor);
	bool solveDiscretizedJacobian(double* rhs) const;

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	double checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const;
#endif

	inline const active& columnLength() const CADET_NOEXCEPT { return _baseOp.columnLength(); }
	inline const active& crossSectionArea() const CADET_NOEXCEPT { return _baseOp.crossSectionArea(); }
	inline const active& currentVelocity() const CADET_NOEXCEPT { return _baseOp.currentVelocity(); }

	inline linalg::BandMatrix& jacobian() CADET_NOEXCEPT { return _jacC; }
	inline const linalg::BandMatrix& jacobian() const CADET_NOEXCEPT { return _jacC; }

	inline linalg::FactorizableBandMatrix& jacobianDisc() CADET_NOEXCEPT { return _jacCdisc; }
	inline const linalg::FactorizableBandMatrix& jacobianDisc() const CADET_NOEXCEPT { return _jacCdisc; }
protected:

	ConvectionDispersionOperatorBase _baseOp;

	linalg::BandMatrix _jacC; //!< Jacobian
	linalg::FactorizableBandMatrix _jacCdisc; //!< Jacobian with time derivatives from BDF method

	// Indexer functionality

	// Strides
	inline const int strideColCell() const CADET_NOEXCEPT { return static_cast<int>(_baseOp.nComp()); }
	inline const int strideColComp() const CADET_NOEXCEPT { return 1; }

	// Offsets
	inline const int offsetC() const CADET_NOEXCEPT { return _baseOp.nComp(); }
	inline const int offsetCp() const CADET_NOEXCEPT { return _baseOp.nComp() * _baseOp.nCol() + offsetC(); }

	// Return pointer to first element of state variable in state vector
	template <typename real_t> inline real_t* c(real_t* const data) const { return data + offsetC(); }
	template <typename real_t> inline real_t const* c(real_t const* const data) const { return data + offsetC(); }

	// Return specific variable in state vector
	template <typename real_t> inline real_t& c(real_t* const data, unsigned int col, unsigned int comp) const { return data[offsetC() + comp + col * strideColCell()]; }
	template <typename real_t> inline const real_t& c(real_t const* const data, unsigned int col, unsigned int comp) const { return data[offsetC() + comp + col * strideColCell()]; }
};

} // namespace parts
} // namespace model
} // namespace cadet

#endif  // LIBCADET_CONVECTIONDISPERSIONOPERATOR_HPP_
