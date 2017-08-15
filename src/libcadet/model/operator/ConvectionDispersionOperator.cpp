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

#include "model/operator/ConvectionDispersionOperator.hpp"
#include "cadet/Exceptions.hpp"

#include "Stencil.hpp"
#include "ParamReaderHelper.hpp"
#include "AdUtils.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

namespace cadet
{

namespace model
{

namespace operators
{

/**
 * @brief Creates a ConvectionDispersionOperator
 */
ConvectionDispersionOperator::ConvectionDispersionOperator() : _stencilMemory(sizeof(active) * Weno::maxStencilSize()), 
	_wenoDerivatives(new double[Weno::maxStencilSize()]), _weno()
{
}

ConvectionDispersionOperator::~ConvectionDispersionOperator() CADET_NOEXCEPT
{
	delete[] _wenoDerivatives;
}

/**
 * @brief Returns the number of AD directions required for computing the Jacobian
 * @details Band compression is used to minimize the amount of AD directions.
 * @return Number of required AD directions
 */
unsigned int ConvectionDispersionOperator::requiredADdirs() const CADET_NOEXCEPT
{
	return _jacC.stride();
}

/**
 * @brief Reads parameters and allocates memory
 * @details Has to be called once before the operator is used.
 * @param [in] unitOpIdx Unit operation id of the owning unit operation model
 * @param [in] paramProvider Parameter provider for reading parameters
 * @param [out] parameters Map in which local parameters are inserted
 * @param [in] nComp Number of components
 * @param [in] nCol Number of axial cells
 * @return @c true if configuration went fine, @c false otherwise
 */
bool ConvectionDispersionOperator::configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, unsigned int nComp, unsigned int nCol)
{
	_nComp = nComp;
	_nCol = nCol;

	paramProvider.pushScope("discretization");

	// Read WENO settings and apply them
	paramProvider.pushScope("weno");
	_weno.order(paramProvider.getInt("WENO_ORDER"));
	_weno.boundaryTreatment(paramProvider.getInt("BOUNDARY_MODEL"));
	_wenoEpsilon = paramProvider.getDouble("WENO_EPS");
	paramProvider.popScope();

	paramProvider.popScope();

	// ==== Read model parameters
	reconfigure(unitOpIdx, paramProvider, parameters);

	// Allocate memory

	// Note that we have to increase the lower bandwidth by 1 because the WENO stencil is applied to the
	// right cell face (lower + 1 + upper) and to the left cell face (shift the stencil by -1 because influx of cell i
	// is outflux of cell i-1)
	// We also have to make sure that there's at least one sub and super diagonal for the dispersion term
	const unsigned int lb = std::max(_weno.lowerBandwidth() + 1u, 1u) * strideColCell();
	const unsigned int ub = std::max(_weno.upperBandwidth(), 1u) * strideColCell();
	const unsigned int mb = std::max(lb, ub);

	// Allocate matrices such that bandwidths can be switched (backwards flow support)
	_jacC.resize(_nCol * _nComp, lb, ub);

	_jacCdisc.resize(_nCol * _nComp, mb, mb);
	_jacCdisc.repartition(lb, ub);
	return true;
}

/**
 * @brief Reads model parameters
 * @details Only reads parameters that do not affect model structure (e.g., discretization).
 * @param [in] unitOpIdx Unit operation id of the owning unit operation model
 * @param [in] paramProvider Parameter provider for reading parameters
 * @param [out] parameters Map in which local parameters are inserted
 * @return @c true if configuration went fine, @c false otherwise
 */
bool ConvectionDispersionOperator::reconfigure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters)
{
	// Read geometry parameters
	_colLength = paramProvider.getDouble("COL_LENGTH");

	// Read cross section area or set to -1
	_crossSection = -1.0;
	if (paramProvider.exists("CROSS_SECTION_AREA"))
	{
		_crossSection = paramProvider.getDouble("CROSS_SECTION_AREA");
	}

	// Read section dependent parameters (transport)

	// Read VELOCITY
	_velocity.clear();
	if (paramProvider.exists("VELOCITY"))
	{
		readScalarParameterOrArray(_velocity, paramProvider, "VELOCITY", 1);
	}
	readScalarParameterOrArray(_colDispersion, paramProvider, "COL_DISPERSION", 1);

	if (_velocity.empty() && (_crossSection <= 0.0))
	{
		throw InvalidParameterException("At least one of CROSS_SECTION_AREA and VELOCITY has to be set");
	}

	// Add parameters to map
	parameters[makeParamId(hashString("COL_LENGTH"), unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &_colLength;
	parameters[makeParamId(hashString("CROSS_SECTION_AREA"), unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &_crossSection;
	registerScalarSectionDependentParam(hashString("COL_DISPERSION"), parameters, _colDispersion, unitOpIdx);
	registerScalarSectionDependentParam(hashString("VELOCITY"), parameters, _velocity, unitOpIdx);

	return true;
}

/**
 * @brief Notifies the operator that a discontinuous section transition is in progress
 * @details In addition to changing flow direction internally, if necessary, the function returns whether
 *          the flow direction has changed.
 * @param [in] t Current time point
 * @param [in] secIdx Index of the new section that is about to be integrated
 * @param [in,out] adRes Pointer to unit operation's residual vector of AD datatypes to be set up (or @c nullptr if AD is disabled)
 * @param [in,out] adY Pointer to unit operation's state vector of AD datatypes to be set up (or @c nullptr if AD is disabled)
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 * @return @c true if flow direction has changed, otherwise @c false
 */
bool ConvectionDispersionOperator::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active* const adRes, active* const adY, unsigned int adDirOffset)
{
	// If we don't have cross section area, velocity is given by parameter
	if (_crossSection <= 0.0)
		_curVelocity = getSectionDependentScalar(_velocity, secIdx);
	else if (!_velocity.empty())
	{
		// We have both cross section area and interstitial flow rate
		// _curVelocity has already been set to the network flow rate in setFlowRates()
		// the direction of the flow (i.e., sign of _curVelocity) is given by _velocity
		const double dir = static_cast<double>(getSectionDependentScalar(_velocity, secIdx));
		if (dir < 0.0)
			_curVelocity *= -1.0;
	}

	// Check whether the matrix connecting inlet DOFs to first column cells has to be (re)assembled
	const double u = static_cast<double>(_curVelocity);
	double prevU = -u;

	// Determine previous flow direction
	if ((secIdx != 0) && !_velocity.empty())
		prevU = static_cast<double>(getSectionDependentScalar(_velocity, secIdx - 1));

	// If interstitial velocity is given by network flow rate and _velocity isn't set, assume forward flow
	if (_velocity.empty())
		prevU = u;

	// Exit if we do not need to setup (secIdx == 0) or change (prevU and u differ in sign) matrices
	if ((secIdx != 0) && (prevU * u >= 0.0))
		return false;

	if (u >= 0.0)
	{
		// Forwards flow

		// Repartition column bulk Jacobians
		const unsigned int lb = std::max(_weno.lowerBandwidth() + 1u, 1u) * strideColCell();
		const unsigned int ub = std::max(_weno.upperBandwidth(), 1u) * strideColCell();

		_jacC.repartition(lb, ub);
		_jacCdisc.repartition(lb, ub);
	}
	else
	{
		// Backwards flow

		// Repartition column bulk Jacobians
		const unsigned int lb = std::max(_weno.lowerBandwidth() + 1u, 1u) * strideColCell();
		const unsigned int ub = std::max(_weno.upperBandwidth(), 1u) * strideColCell();

		_jacC.repartition(ub, lb);
		_jacCdisc.repartition(ub, lb);
	}

	// Update AD seed vectors since Jacobian structure has changed (bulk block bandwidths)
	prepareADvectors(adRes, adY, adDirOffset);

	return true;
}

/**
 * @brief Sets the AD seed vectors for the bulk transport variables
 * @param [in,out] adRes Pointer to unit operation's residual vector of AD datatypes to be set up (or @c nullptr if AD is disabled)
 * @param [in,out] adY Pointer to unit operation's state vector of AD datatypes to be set up (or @c nullptr if AD is disabled)
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
void ConvectionDispersionOperator::prepareADvectors(active* const adRes, active* const adY, unsigned int adDirOffset) const
{
	// Early out if AD is disabled
	if (!adY)
		return;

	// Get bandwidths of blocks
	const unsigned int lowerColBandwidth = _jacC.lowerBandwidth();
	const unsigned int upperColBandwidth = _jacC.upperBandwidth();

	// Column block
	ad::prepareAdVectorSeedsForBandMatrix(adY + offsetC(), adDirOffset, _nCol * _nComp, lowerColBandwidth, upperColBandwidth, lowerColBandwidth);
}

/**
 * @brief Sets the flow rates for the current time section
 * @details The flow rates may change due to valve switches.
 * @param [in] in Total volumetric inlet flow rate
 * @param [in] out Total volumetric outlet flow rate
 * @param [in] colPorosity Porosity used for computing interstitial velocity from volumetric flow rate
 */
void ConvectionDispersionOperator::setFlowRates(const active& in, const active& out, const active& colPorosity) CADET_NOEXCEPT
{
	// If we have cross section area, interstitial velocity is given by network flow rates
	if (_crossSection > 0.0)
		_curVelocity = in / (_crossSection * colPorosity);
}

/**
 * @brief Computes the residual of the transport equations
 * @param [in] t Current time point
 * @param [in] secIdx Index of the current section
 * @param [in] timeFactor Used to compute parameter derivatives with respect to section length (nominal value should always be 1.0)
 * @param [in] y Pointer to unit operation's state vector
 * @param [in] yDot Pointer to unit operation's time derivative state vector
 * @param [out] res Pointer to unit operation's residual vector
 * @param [in] wantJac Determines whether analytic Jacobian is computed
 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
 */
int ConvectionDispersionOperator::residual(double t, unsigned int secIdx, double timeFactor, double const* y, double const* yDot, double* res, bool wantJac)
{
	if (wantJac)
		return residualImpl<double, double, double, true>(t, secIdx, timeFactor, y, yDot, res);
	else
		return residualImpl<double, double, double, false>(t, secIdx, timeFactor, y, yDot, res);
}

int ConvectionDispersionOperator::residual(double t, unsigned int secIdx, double timeFactor, active const* y, double const* yDot, active* res, bool wantJac)
{
	if (wantJac)
		return residualImpl<active, active, double, true>(t, secIdx, timeFactor, y, yDot, res);
	else
		return residualImpl<active, active, double, false>(t, secIdx, timeFactor, y, yDot, res);
}

int ConvectionDispersionOperator::residual(const active& t, unsigned int secIdx, const active& timeFactor, double const* y, double const* yDot, active* res, bool wantJac)
{
	if (wantJac)
		return residualImpl<double, active, active, true>(t, secIdx, timeFactor, y, yDot, res);
	else
		return residualImpl<double, active, active, false>(t, secIdx, timeFactor, y, yDot, res);
}

int ConvectionDispersionOperator::residual(const active& t, unsigned int secIdx, const active& timeFactor, active const* y, double const* yDot, active* res, bool wantJac)
{
	if (wantJac)
		return residualImpl<active, active, active, true>(t, secIdx, timeFactor, y, yDot, res);
	else
		return residualImpl<active, active, active, false>(t, secIdx, timeFactor, y, yDot, res);
}

template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
int ConvectionDispersionOperator::residualImpl(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* y, double const* yDot, ResidualType* res)
{
	const ParamType u = static_cast<ParamType>(_curVelocity);
	if (u >= 0.0)
		return residualForwardsFlow<StateType, ResidualType, ParamType, wantJac>(t, secIdx, timeFactor, y, yDot, res);
	else
		return residualBackwardsFlow<StateType, ResidualType, ParamType, wantJac>(t, secIdx, timeFactor, y, yDot, res);
}

template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
int ConvectionDispersionOperator::residualForwardsFlow(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* y, double const* yDot, ResidualType* res)
{
	const ParamType u = static_cast<ParamType>(_curVelocity);
	const ParamType d_c = static_cast<ParamType>(getSectionDependentScalar(_colDispersion, secIdx));
	const ParamType h = static_cast<ParamType>(_colLength) / static_cast<double>(_nCol);
	const ParamType h2 = h * h;

	const int strideCell = strideColCell();

	// The stencil caches parts of the state vector for better spatial coherence
	typedef CachingStencil<StateType, ArrayPool> StencilType;
	StencilType stencil(std::max(_weno.stencilSize(), 3u), _stencilMemory, std::max(_weno.order() - 1, 1));

	// Reset Jacobian
	if (wantJac)
		_jacC.setAll(0.0);

	for (unsigned int comp = 0; comp < _nComp; ++comp)
	{
		// The RowIterator is always centered on the main diagonal.
		// This means that jac[0] is the main diagonal, jac[-1] is the first lower diagonal,
		// and jac[1] is the first upper diagonal. We can also access the rows from left to
		// right beginning with the last lower diagonal moving towards the main diagonal and
		// continuing to the last upper diagonal by using the native() method.
		linalg::BandMatrix::RowIterator jac = _jacC.row(comp);

		// Add time derivative to each cell
		if (yDot)
		{
			for (unsigned int col = 0; col < _nCol; ++col)
				c<ResidualType>(res, col, comp) = timeFactor * c<double>(yDot, col, comp);
		}
		else
		{
			for (unsigned int col = 0; col < _nCol; ++col)
				c<ResidualType>(res, col, comp) = 0.0;
		}

		// Fill stencil (left side with zeros, right side with states)
		for (int i = -std::max(_weno.order(), 2) + 1; i < 0; ++i)
			stencil[i] = 0.0;
		for (int i = 0; i < std::max(_weno.order(), 2); ++i)
			stencil[i] = c<StateType>(y, static_cast<unsigned int>(i), comp);

		// Reset WENO output
		StateType vm(0.0); // reconstructed value
		if (wantJac)
			std::fill(_wenoDerivatives, _wenoDerivatives + _weno.stencilSize(), 0.0);

		int wenoOrder = 0;

		// Iterate over all cells
		for (unsigned int col = 0; col < _nCol; ++col)
		{
			// ------------------- Dispersion -------------------

			// Right side, leave out if we're in the last cell (boundary condition)
			if (cadet_likely(col < _nCol - 1))
			{
				c<ResidualType>(res, col, comp) -= d_c / h2 * (stencil[1] - stencil[0]);
				// Jacobian entries
				if (wantJac)
				{
					jac[0] += static_cast<double>(d_c) / static_cast<double>(h2);
					jac[strideCell] -= static_cast<double>(d_c) / static_cast<double>(h2);
				}
			}

			// Left side, leave out if we're in the first cell (boundary condition)
			if (cadet_likely(col > 0))
			{
				c<ResidualType>(res, col, comp) -= d_c / h2 * (stencil[-1] - stencil[0]);
				// Jacobian entries
				if (wantJac)
				{
					jac[0]  += static_cast<double>(d_c) / static_cast<double>(h2);
					jac[-strideCell] -= static_cast<double>(d_c) / static_cast<double>(h2);
				}
			}

			// ------------------- Convection -------------------

			// Add convection through this cell's left face
			if (cadet_likely(col > 0))
			{
				// Remember that vm still contains the reconstructed value of the previous 
				// cell's *right* face, which is identical to this cell's *left* face!
				c<ResidualType>(res, col, comp) -= u / h * vm;

				// Jacobian entries
				if (wantJac)
				{
					for (int i = 0; i < 2 * wenoOrder - 1; ++i)
						// Note that we have an offset of -1 here (compared to the right cell face below), since
						// the reconstructed value depends on the previous stencil (which has now been moved by one cell)
						jac[(i - wenoOrder) * strideCell] -= static_cast<double>(u) / static_cast<double>(h) * _wenoDerivatives[i];
				}
			}
			else
			{
				// In the first cell we need to apply the boundary condition: inflow concentration
				c<ResidualType>(res, col, comp) -= u / h * y[comp];
			}

			// Reconstruct concentration on this cell's right face
			if (wantJac)
				wenoOrder = _weno.reconstruct<StateType, StencilType>(_wenoEpsilon, col, _nCol, stencil, vm, _wenoDerivatives);
			else
				wenoOrder = _weno.reconstruct<StateType, StencilType>(_wenoEpsilon, col, _nCol, stencil, vm);

			// Right side
			c<ResidualType>(res, col, comp) += u / h * vm;
			// Jacobian entries
			if (wantJac)
			{
				for (int i = 0; i < 2 * wenoOrder - 1; ++i)
					jac[(i - wenoOrder + 1) * strideCell] += static_cast<double>(u) / static_cast<double>(h) * _wenoDerivatives[i];
			}

			// Update stencil
			stencil.advance(c<StateType>(y, col + std::max(_weno.order(), 2), comp));
			jac += strideCell;
		}
	}

	// Film diffusion with flux into beads is added in residualFlux() function

	return 0;
}

template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
int ConvectionDispersionOperator::residualBackwardsFlow(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* y, double const* yDot, ResidualType* res)
{
	const ParamType u = static_cast<ParamType>(_curVelocity);
	const ParamType d_c = static_cast<ParamType>(getSectionDependentScalar(_colDispersion, secIdx));
	const ParamType h = static_cast<ParamType>(_colLength) / static_cast<double>(_nCol);
	const ParamType h2 = h * h;

	const int strideCell = strideColCell();

	// The stencil caches parts of the state vector for better spatial coherence
	typedef CachingStencil<StateType, ArrayPool> StencilType;
	StencilType stencil(std::max(_weno.stencilSize(), 3u), _stencilMemory, std::max(_weno.order() - 1, 1));

	// Reset Jacobian
	if (wantJac)
		_jacC.setAll(0.0);

	for (unsigned int comp = 0; comp < _nComp; ++comp)
	{
		// The RowIterator is always centered on the main diagonal.
		// This means that jac[0] is the main diagonal, jac[-1] is the first lower diagonal,
		// and jac[1] is the first upper diagonal. We can also access the rows from left to
		// right beginning with the last lower diagonal moving towards the main diagonal and
		// continuing to the last upper diagonal by using the native() method.
		linalg::BandMatrix::RowIterator jac = _jacC.row(_nComp * (_nCol - 1) + comp);

		// Add time derivative to each cell
		if (yDot)
		{
			for (unsigned int col = 0; col < _nCol; ++col)
				c<ResidualType>(res, col, comp) = timeFactor * c<double>(yDot, col, comp);
		}
		else
		{
			for (unsigned int col = 0; col < _nCol; ++col)
				c<ResidualType>(res, col, comp) = 0.0;
		}

		// Fill stencil (left side with zeros, right side with states)
		for (int i = -std::max(_weno.order(), 2) + 1; i < 0; ++i)
			stencil[i] = 0.0;
		for (int i = 0; i < std::max(_weno.order(), 2); ++i)
			stencil[i] = c<StateType>(y, _nCol - static_cast<unsigned int>(i) - 1, comp);

		// Reset WENO output
		StateType vm(0.0); // reconstructed value
		if (wantJac)
			std::fill(_wenoDerivatives, _wenoDerivatives + _weno.stencilSize(), 0.0);

		int wenoOrder = 0;

		// Iterate over all cells (backwards)
		// Note that col wraps around to unsigned int's maximum value after 0
		for (unsigned int col = _nCol - 1; col < _nCol; --col)
		{
			// ------------------- Dispersion -------------------

			// Right side, leave out if we're in the first cell (boundary condition)
			if (cadet_likely(col < _nCol - 1))
			{
				c<ResidualType>(res, col, comp) -= d_c / h2 * (stencil[-1] - stencil[0]);
				// Jacobian entries
				if (wantJac)
				{
					jac[0] += static_cast<double>(d_c) / static_cast<double>(h2);
					jac[strideCell] -= static_cast<double>(d_c) / static_cast<double>(h2);
				}
			}

			// Left side, leave out if we're in the last cell (boundary condition)
			if (cadet_likely(col > 0))
			{
				c<ResidualType>(res, col, comp) -= d_c / h2 * (stencil[1] - stencil[0]);
				// Jacobian entries
				if (wantJac)
				{
					jac[0] += static_cast<double>(d_c) / static_cast<double>(h2);
					jac[-strideCell] -= static_cast<double>(d_c) / static_cast<double>(h2);
				}
			}

			// ------------------- Convection -------------------

			// Add convection through this cell's right face
			if (cadet_likely(col < _nCol - 1))
			{
				// Remember that vm still contains the reconstructed value of the previous 
				// cell's *left* face, which is identical to this cell's *right* face!
				c<ResidualType>(res, col, comp) += u / h * vm;

				// Jacobian entries
				if (wantJac)
				{
					for (int i = 0; i < 2 * wenoOrder - 1; ++i)
						// Note that we have an offset of +1 here (compared to the left cell face below), since
						// the reconstructed value depends on the previous stencil (which has now been moved by one cell)
						jac[(wenoOrder - i) * strideCell] += static_cast<double>(u) / static_cast<double>(h) * _wenoDerivatives[i];					
				}
			}
			else
			{
				// In the last cell (z = L) we need to apply the boundary condition: inflow concentration
				c<ResidualType>(res, col, comp) += u / h * y[comp];
			}

			// Reconstruct concentration on this cell's left face
			if (wantJac)
				wenoOrder = _weno.reconstruct<StateType, StencilType>(_wenoEpsilon, col, _nCol, stencil, vm, _wenoDerivatives);
			else
				wenoOrder = _weno.reconstruct<StateType, StencilType>(_wenoEpsilon, col, _nCol, stencil, vm);

			// Left face
			c<ResidualType>(res, col, comp) -= u / h * vm;
			// Jacobian entries
			if (wantJac)
			{
				for (int i = 0; i < 2 * wenoOrder - 1; ++i)
					jac[(wenoOrder - i - 1) * strideCell] -= static_cast<double>(u) / static_cast<double>(h) * _wenoDerivatives[i];				
			}

			// Update stencil (be careful because of wrap-around, might cause reading memory very far away [although never used])
			const unsigned int shift = std::max(_weno.order(), 2);
			if (cadet_likely(col - shift < _nCol))
				stencil.advance(c<StateType>(y, col - shift, comp));
			else
				stencil.advance(0.0);
			jac -= strideCell;
		}
	}

	// Film diffusion with flux into beads is added in residualFlux() function

	return 0;
}

/**
 * @brief Multiplies the time derivative Jacobian @f$ \frac{\partial F}{\partial \dot{y}}\left(t, y, \dot{y}\right) @f$ with a given vector
 * @details The operation @f$ z = \frac{\partial F}{\partial \dot{y}} x @f$ is performed.
 *          The matrix-vector multiplication is transformed matrix-free (i.e., no matrix is explicitly formed).
 *          
 *          Note that this function only performs multiplication with the Jacobian of the (axial) transport equations.
 *          The input vectors are assumed to point to the beginning (including inlet DOFs) of the respective unit operation's arrays.
 * @param [in] t Current time point
 * @param [in] secIdx Index of the current section
 * @param [in] timeFactor Factor which is premultiplied to the time derivatives originating from time transformation
 * @param [in] sDot Vector @f$ x @f$ that is transformed by the Jacobian @f$ \frac{\partial F}{\partial \dot{y}} @f$
 * @param [out] ret Vector @f$ z @f$ which stores the result of the operation
 */
void ConvectionDispersionOperator::multiplyWithDerivativeJacobian(double t, unsigned int secIdx, double timeFactor, double const* sDot, double* ret) const
{
	for (int i = offsetC(); i < offsetCp(); ++i)
		ret[i] = timeFactor * sDot[i];
}

/**
 * @brief Extracts the system Jacobian from band compressed AD seed vectors
 * @details The input vectors are assumed to point to the beginning (including inlet DOFs) of the respective unit operation's arrays.
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
void ConvectionDispersionOperator::extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset)
{
	ad::extractBandedJacobianFromAd(adRes + offsetC(), adDirOffset, _jacC.lowerBandwidth(), _jacC);
}

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN

/**
 * @brief Compares the analytical Jacobian with a Jacobian derived by AD
 * @details The analytical Jacobian is assumed to be stored in the corresponding band matrices.
 *          The input vectors are assumed to point to the beginning (including inlet DOFs) of the respective unit operation's arrays.
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 * @return Maximum elementwise absolute difference between analytic and AD Jacobian
 */
double ConvectionDispersionOperator::checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const
{
	// Column
	const double maxDiffCol = ad::compareBandedJacobianWithAd(adRes + offsetC(), adDirOffset, _jacC.lowerBandwidth(), _jacC);
	LOG(Debug) << "-> Col block diff: " << maxDiffCol;

	return maxDiffCol;
}

#endif

/**
 * @brief Assembles the axial transport Jacobian @f$ J_0 @f$ of the time-discretized equations
 * @details The system \f[ \left( \frac{\partial F}{\partial y} + \alpha \frac{\partial F}{\partial \dot{y}} \right) x = b \f]
 *          has to be solved. The system Jacobian of the original equations,
 *          \f[ \frac{\partial F}{\partial y}, \f]
 *          is already computed (by AD or manually in residualImpl() with @c wantJac = true). This function is responsible
 *          for adding
 *          \f[ \alpha \frac{\partial F}{\partial \dot{y}} \f]
 *          to the system Jacobian, which yields the Jacobian of the time-discretized equations
 *          \f[ F\left(t, y_0, \sum_{k=0}^N \alpha_k y_k \right) = 0 \f]
 *          when a BDF method is used. The time integrator needs to solve this equation for @f$ y_0 @f$, which requires
 *          the solution of the linear system mentioned above (@f$ \alpha_0 = \alpha @f$ given in @p alpha).
 *
 * @param [in] alpha Value of \f$ \alpha \f$ (arises from BDF time discretization)
 * @param [in] timeFactor Factor which is premultiplied to the time derivatives originating from time transformation
 */
void ConvectionDispersionOperator::assembleDiscretizedJacobian(double alpha, double timeFactor)
{
	// Copy normal matrix over to factorizable matrix
	_jacCdisc.copyOver(_jacC);

	// Add time derivatives
	addTimeDerivativeToJacobian(alpha, timeFactor);
}


/**
 * @brief Adds the derivatives with respect to @f$ \dot{y} @f$ of @f$ F(t, y, \dot{y}) @f$ to the Jacobian
 * @details This functions computes 
 *          @f[ \begin{align*} \text{_jacCdisc} = \text{_jacCdisc} + \alpha \frac{\partial F}{\partial \dot{y}}. \end{align*} @f]
 *          The factor @f$ \alpha @f$ is useful when constructing the linear system in the time integration process.
 * @param [in] alpha Factor in front of @f$ \frac{\partial F}{\partial \dot{y}} @f$
 * @param [in] timeFactor Factor which is premultiplied to the time derivatives originating from time transformation
 */
void ConvectionDispersionOperator::addTimeDerivativeToJacobian(double alpha, double timeFactor)
{
	alpha *= timeFactor;

	const int gapCell = strideColCell() - static_cast<int>(_nComp) * strideColComp();
	linalg::FactorizableBandMatrix::RowIterator jac = _jacCdisc.row(0);
	for (unsigned int i = 0; i < _nCol; ++i, jac += gapCell)
	{
		for (unsigned int j = 0; j < _nComp; ++j, ++jac)
		{
			// Add time derivative to main diagonal
			jac[0] += alpha;
		}
	}
}

/**
 * @brief Assembles and factorizes the time discretized Jacobian
 * @details See assembleDiscretizedJacobian() for assembly of the time discretized Jacobian.
 * @param [in] alpha Factor in front of @f$ \frac{\partial F}{\partial \dot{y}} @f$
 * @param [in] timeFactor Factor which is premultiplied to the time derivatives originating from time transformation
 * @return @c true if factorization went fine, otherwise @c false
 */
bool ConvectionDispersionOperator::assembleAndFactorizeDiscretizedJacobian(double alpha, double timeFactor)
{
	assembleDiscretizedJacobian(alpha, timeFactor);
	return _jacCdisc.factorize();
}

/**
 * @brief Solves a (previously factorized) equation system
 * @details The (time discretized) Jacobian matrix has to be factorized before calling this function.
 *          Note that the given right hand side vector @p rhs is not shifted by the inlet DOFs. That
 *          is, it is assumed to point directly to the first axial DOF.
 * 
 * @param [in,out] rhs On entry, right hand side of the equation system. On exit, solution of the system.
 * @return @c true if the system was solved correctly, otherwise @c false
 */
bool ConvectionDispersionOperator::solveDiscretizedJacobian(double* rhs) const
{
	return _jacCdisc.solve(rhs);
}

/**
 * @brief Solves a system with the time derivative Jacobian and given right hand side
 * @details Note that the given right hand side vector @p rhs is not shifted by the inlet DOFs. That
 *          is, it is assumed to point directly to the first axial DOF.
 * @param [in] t Current time point
 * @param [in] secIdx Index of the current section
 * @param [in] timeFactor Factor which is premultiplied to the time derivatives originating from time transformation
 * @param [in,out] rhs On entry, right hand side. On exit, solution of the system.
 * @return @c true if the system was solved correctly, @c false otherwise
 */
bool ConvectionDispersionOperator::solveTimeDerivativeSystem(double t, unsigned int secIdx, double timeFactor, double* const rhs)
{
	// Assemble
	_jacCdisc.setAll(0.0);
	addTimeDerivativeToJacobian(1.0, timeFactor);

	// Factorize
	const bool result = _jacCdisc.factorize();
	if (!result)
	{
		LOG(Error) << "Factorize() failed for bulk block";
		return false;
	}

	// Solve
	const bool result2 = _jacCdisc.solve(rhs);
	if (!result2)
	{
		LOG(Error) << "Solve() failed for bulk block";
		return false;
	}

	return true;
}

}  // namespace operators

}  // namespace model

}  // namespace cadet
