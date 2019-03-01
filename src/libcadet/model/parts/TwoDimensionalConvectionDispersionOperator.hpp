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
 * Defines the 2D convection dispersion transport operator.
 */

#ifndef LIBCADET_2DCONVECTIONDISPERSIONOPERATOR_HPP_
#define LIBCADET_2DCONVECTIONDISPERSIONOPERATOR_HPP_

#include "ParamIdUtil.hpp"
#include "AutoDiff.hpp"
#include "linalg/CompressedSparseMatrix.hpp"
#include "MemoryPool.hpp"
#include "Weno.hpp"

#include <unordered_map>
#include <vector>

namespace cadet
{

class IParameterProvider;
struct SimulationTime;

namespace model
{

namespace parts
{

/**
 * @brief 2D Convection dispersion transport operator
 * @details Implements the equation
 * 
 * @f[\begin{align}
	\frac{\partial c_i}{\partial t} &= - u \frac{\partial c_i}{\partial z} + D_{\text{ax},i}(\rho) \frac{\partial^2 c_i}{\partial z^2} + \frac{1}{\rho} \frac{\partial}{\partial \rho} \left( \rho D_{\rho} \frac{\partial c_i}{\partial \rho} \right) \\
\end{align} @f]
 * with Danckwerts boundary conditions on the axial boundary (see @cite Danckwerts1953)
@f[ \begin{align}
u c_{\text{in},i}(t) &= u c_i(t,0,\rho) - D_{\text{ax},i}(\rho) \frac{\partial c_i}{\partial z}(t,0,\rho) \\
\frac{\partial c_i}{\partial z}(t,L,\rho) &= 0
\end{align} @f]
 * and Neumann boundary conditions on the radial boundary
@f[ \begin{align}
\frac{\partial c_i}{\partial \rho}(t,z,0) &= 0 \\
\frac{\partial c_i}{\partial z}(t,z,R) &= 0
\end{align} @f]
 * Methods are described in @cite VonLieres2010a (WENO, linear solver), and @cite Puttmann2013, @cite Puttmann2016 (forward sensitivities, AD, band compression)
 */
class TwoDimensionalConvectionDispersionOperator
{
public:

	TwoDimensionalConvectionDispersionOperator();
	~TwoDimensionalConvectionDispersionOperator() CADET_NOEXCEPT;

//	unsigned int requiredADdirs() const CADET_NOEXCEPT;

	void setFlowRates(int compartment, const active& in, const active& out) CADET_NOEXCEPT;

	bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int nCol, unsigned int nRad);
	bool configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters);
	bool notifyDiscontinuousSectionTransition(double t, unsigned int secIdx);

	int residual(double t, unsigned int secIdx, double timeFactor, double const* y, double const* yDot, double* res, bool wantJac);
	int residual(double t, unsigned int secIdx, double timeFactor, active const* y, double const* yDot, active* res, bool wantJac);
	int residual(const active& t, unsigned int secIdx, const active& timeFactor, double const* y, double const* yDot, active* res, bool wantJac);
	int residual(const active& t, unsigned int secIdx, const active& timeFactor, active const* y, double const* yDot, active* res, bool wantJac);

/*
	void prepareADvectors(const AdJacobianParams& adJac) const;
	void extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset);
*/

	bool solveTimeDerivativeSystem(const SimulationTime& simTime, double* const rhs);
	void multiplyWithDerivativeJacobian(const SimulationTime& simTime, double const* sDot, double* ret) const;

	bool assembleAndFactorizeDiscretizedJacobian(double alpha, double timeFactor);
	bool solveDiscretizedJacobian(double* rhs, double const* weight, double const* init, double outerTol) const;

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
//	double checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const;
#endif

	inline const active& columnLength() const CADET_NOEXCEPT { return _colLength; }
	inline const active& columnRadius() const CADET_NOEXCEPT { return _colRadius; }
	inline const active& currentVelocity(int idx) const CADET_NOEXCEPT { return _curVelocity[idx]; }
	inline const active& columnPorosity(int idx) const CADET_NOEXCEPT { return _colPorosities[idx]; }
	inline const active& crossSection(int idx) const CADET_NOEXCEPT { return _crossSections[idx]; }
	inline active const* crossSections() const CADET_NOEXCEPT { return _crossSections.data(); }
	inline active const* radialCenters() const CADET_NOEXCEPT { return _radialCenters.data(); }
//	inline active const* radialCentroids() const CADET_NOEXCEPT { return _radialCentroids.data(); }
	inline active const* radialEdges() const CADET_NOEXCEPT { return _radialEdges.data(); }
	inline bool isCurrentFlowForward(int idx) const CADET_NOEXCEPT { return _curVelocity[idx] >= 0.0; }
	const active& axialDispersion(unsigned int idxSec, int idxRad, int idxComp) const CADET_NOEXCEPT;
	const active& radialDispersion(unsigned int idxSec, int idxRad, int idxComp) const CADET_NOEXCEPT;
	double inletFactor(unsigned int idxSec, int idxRad) const CADET_NOEXCEPT;

	inline linalg::CompressedSparseMatrix& jacobian() CADET_NOEXCEPT { return _jacC; }
	inline const linalg::CompressedSparseMatrix& jacobian() const CADET_NOEXCEPT { return _jacC; }

#ifdef CADET_SPARSE_DIRECT
	inline linalg::FactorizableCompressedSparseMatrix& jacobianDisc() CADET_NOEXCEPT { return _jacCdisc; }
	inline const linalg::FactorizableCompressedSparseMatrix& jacobianDisc() const CADET_NOEXCEPT { return _jacCdisc; }
#endif

protected:

	class LinearSolver;
	class GmresSolver;
	template <typename sparse_t> class SparseDirectSolver;
	class DenseDirectSolver;

	friend int schurComplementMultiplier2DCDO(void* userData, double const* x, double* z);

	void assembleDiscretizedJacobian(double alpha, double timeFactor);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	int residualImpl(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* y, double const* yDot, ResidualType* res);

	void setSparsityPattern();

	void setEquidistantRadialDisc();
	void setEquivolumeRadialDisc();
	void setUserdefinedRadialDisc();
	void updateRadialDisc();

	enum class RadialDiscretizationMode : int
	{
		/**
		 * Equidistant distribution of tubular shell edges
		 */
		Equidistant,

		/**
		 * Volumes of tubular shells are uniform
		 */
		Equivolume,

		/**
		 * Shell edges specified by user
		 */
		UserDefined
	};

	unsigned int _nComp; //!< Number of components
	unsigned int _nCol; //!< Number of axial cells
	unsigned int _nRad; //!< Number of radial cells

	active _colLength; //!< Column length \f$ L \f$
	active _colRadius; //!< Column radius \f$ r_c \f$
	std::vector<active> _radialEdges; //!< Boundaries of the radial compartments
	std::vector<active> _radialCenters; //!< Center of each radial compartment
//	std::vector<active> _radialCentroids; //!< Center of mass of each radial compartment
	std::vector<active> _crossSections; //!< Cross section area of each compartment 
	std::vector<active> _colPorosities; //!< Bulk porosity for each compartment
	RadialDiscretizationMode _radialDiscretizationMode;

	std::vector<active> _axialDispersion; //!< Axial dispersion coefficient \f$ D_{\text{ax}} \f$
	std::vector<active> _radialDispersion; //!< Radial dispersion coefficient \f$ D_{\rho} \f$
	std::vector<active> _velocity; //!< Interstitial velocity parameter
	std::vector<active> _curVelocity; //!< Current interstitial velocity \f$ u \f$

	ArrayPool _stencilMemory; //!< Provides memory for the stencil
	double* _wenoDerivatives; //!< Holds derivatives of the WENO scheme
	Weno _weno; //!< The WENO scheme implementation
	double _wenoEpsilon; //!< The @f$ \varepsilon @f$ of the WENO scheme (prevents division by zero)

	linalg::CompressedSparseMatrix _jacC; //!< Jacobian
	LinearSolver* _linearSolver; //!< Solves linear system with time discretized Jacobian
};

} // namespace parts
} // namespace model
} // namespace cadet

#endif  // LIBCADET_2DCONVECTIONDISPERSIONOPERATOR_HPP_
