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

#include "model/parts/TwoDimensionalConvectionDispersionOperator.hpp"
#include "cadet/Exceptions.hpp"

#include "Stencil.hpp"
#include "ParamReaderHelper.hpp"
#include "AdUtils.hpp"
#include "SimulationTypes.hpp"
#include "linalg/CompressedSparseMatrix.hpp"
#include "model/parts/ConvectionDispersionKernel.hpp"

#ifdef SUPERLU_FOUND
	#include "linalg/SuperLUSparseMatrix.hpp"
#endif
#ifdef UMFPACK_FOUND
	#include "linalg/UMFPackSparseMatrix.hpp"
#endif

#include "linalg/BandMatrix.hpp"
#include "linalg/Gmres.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <algorithm>

namespace cadet
{

namespace model
{

namespace parts
{

class TwoDimensionalConvectionDispersionOperator::LinearSolver
{
public:

	virtual ~LinearSolver() CADET_NOEXCEPT { }

	virtual bool initialize(IParameterProvider& paramProvider, unsigned int nComp, unsigned int nCol, unsigned int nRad, const Weno& weno) = 0;
	virtual void setSparsityPattern(const cadet::linalg::SparsityPattern& pattern) = 0;
	virtual void assembleDiscretizedJacobian(double alpha, double timeFactor) = 0;
	virtual bool factorize() = 0;
	virtual bool solveDiscretizedJacobian(double* rhs, double const* weight, double const* init, double outerTol) const = 0;
};

int schurComplementMultiplier2DCDO(void* userData, double const* x, double* z);

class TwoDimensionalConvectionDispersionOperator::GmresSolver : public TwoDimensionalConvectionDispersionOperator::LinearSolver
{
public:

	GmresSolver(linalg::CompressedSparseMatrix const* jacC) : _jacC(jacC) { }
	virtual ~GmresSolver() CADET_NOEXCEPT { }

	virtual bool initialize(IParameterProvider& paramProvider, unsigned int nComp, unsigned int nCol, unsigned int nRad, const Weno& weno)
	{
		_gmres.initialize(nCol * nComp * nRad, 0, linalg::toOrthogonalization(1), 0);
		_gmres.matrixVectorMultiplier(&schurComplementMultiplier2DCDO, this);
		_cache.resize(nCol * nComp * nRad, 0.0);

		return true;
	}

	virtual void setSparsityPattern(const linalg::SparsityPattern& pattern) { }

	virtual void assembleDiscretizedJacobian(double alpha, double timeFactor)
	{
		_alpha = alpha * timeFactor;
	}

	virtual bool factorize() { return true; }

	virtual bool solveDiscretizedJacobian(double* rhs, double const* weight, double const* init, double outerTol) const 
	{
		if (init)
			std::copy(init, init + _cache.size(), _cache.begin());

		const int gmresResult = _gmres.solve(0.05 * outerTol * std::sqrt(_jacC->rows()), weight, rhs, _cache.data());
		std::copy(_cache.begin(), _cache.end(), rhs);

		return gmresResult == 0;
	}

protected:
	linalg::CompressedSparseMatrix const* const _jacC;
	double _alpha;
	mutable linalg::Gmres _gmres; //!< GMRES algorithm for the Schur-complement in linearSolve()
	mutable std::vector<double> _cache; //!< GMRES cache for result

	int schurComplementMatrixVector(double const* x, double* z) const
	{
		std::fill(z, z + _jacC->rows(), _alpha);
		_jacC->multiplyVector(x, 1.0, 1.0, z);
		return 0;
	}

	// Wrapper for calling the corresponding function in this class
	friend int schurComplementMultiplier2DCDO(void* userData, double const* x, double* z);
};

int schurComplementMultiplier2DCDO(void* userData, double const* x, double* z)
{
	TwoDimensionalConvectionDispersionOperator::GmresSolver* const cdo = static_cast<TwoDimensionalConvectionDispersionOperator::GmresSolver*>(userData);
	return cdo->schurComplementMatrixVector(x, z);
}

#if defined(UMFPACK_FOUND) || defined(SUPERLU_FOUND) 

	template <typename sparse_t>
	class TwoDimensionalConvectionDispersionOperator::SparseDirectSolver : public TwoDimensionalConvectionDispersionOperator::LinearSolver
	{
	public:

		SparseDirectSolver(linalg::CompressedSparseMatrix const* jacC) : _jacC(jacC) { }
		virtual ~SparseDirectSolver() CADET_NOEXCEPT { }

		virtual bool initialize(IParameterProvider& paramProvider, unsigned int nComp, unsigned int nCol, unsigned int nRad, const Weno& weno)
		{
			return true;
		}

		virtual void setSparsityPattern(const linalg::SparsityPattern& pattern)
		{
			_jacCdisc.assignPattern(pattern);
			_jacCdisc.prepare();
		}

		virtual void assembleDiscretizedJacobian(double alpha, double timeFactor)
		{
			// Copy normal matrix over to factorizable matrix
			_jacCdisc.copyFromSamePattern(*_jacC);

			// Add time derivatives
			alpha *= timeFactor;

			for (unsigned int i = 0; i < _jacC->rows(); ++i)
				_jacCdisc.centered(i, 0) += alpha;
		}

		virtual bool factorize()
		{
			return _jacCdisc.factorize();
		}

		virtual bool solveDiscretizedJacobian(double* rhs, double const* weight, double const* init, double outerTol) const 
		{
			return _jacCdisc.solve(rhs);
		}

	protected:
		linalg::CompressedSparseMatrix const* const _jacC;
		sparse_t _jacCdisc;
	};

#endif

class TwoDimensionalConvectionDispersionOperator::DenseDirectSolver : public TwoDimensionalConvectionDispersionOperator::LinearSolver
{
public:

	DenseDirectSolver(linalg::CompressedSparseMatrix const* jacC) : _jacC(jacC) { }
	virtual ~DenseDirectSolver() CADET_NOEXCEPT { }

	virtual bool initialize(IParameterProvider& paramProvider, unsigned int nComp, unsigned int nCol, unsigned int nRad, const Weno& weno)
	{
		// Note that we have to increase the lower bandwidth by 1 because the WENO stencil is applied to the
		// right cell face (lower + 1 + upper) and to the left cell face (shift the stencil by -1 because influx of cell i
		// is outflux of cell i-1)
		// We also have to make sure that there's at least one sub and super diagonal for the dispersion term
		const unsigned int lb = std::max(weno.lowerBandwidth() + 1u, 1u) * nComp * nRad;

		// We have to make sure that there's at least one sub and super diagonal for the dispersion term
		const unsigned int ub = std::max(weno.upperBandwidth(), 1u) * nComp * nRad;

		// When flow direction is changed, the bandwidths of the Jacobian swap.
		// Hence, we have to reserve memory such that the swapped Jacobian can fit into the matrix.
		const unsigned int mb = std::max(lb, ub);

		// Allocate matrices such that bandwidths can be switched (backwards flow support)
		_jacCdisc.resize(nCol * nComp * nRad, mb, mb);
		return true;
	}

	virtual void setSparsityPattern(const linalg::SparsityPattern& pattern) { }

	virtual void assembleDiscretizedJacobian(double alpha, double timeFactor)
	{
		alpha *= timeFactor;

		// Copy normal matrix over to factorizable matrix
		_jacCdisc.setAll(0.0);

		linalg::FactorizableBandMatrix::RowIterator jac = _jacCdisc.row(0);
		for (int i = 0; i < _jacC->rows(); ++i, ++jac)
		{
			int const* const colIdx = _jacC->columnIndicesOfRow(i);
			double const* const vals = _jacC->valuesOfRow(i);
			const int nnz = _jacC->numNonZerosInRow(i);

			// Copy row from sparse matrix to banded matrix
			for (int c = 0; c < nnz; ++c)
			{
				const int diag = colIdx[c] - i;
				jac[diag] = vals[c];
			}

			// Add time derivative			
			jac[0] += alpha;
		}
	}

	virtual bool factorize()
	{
		return _jacCdisc.factorize();
	}

	virtual bool solveDiscretizedJacobian(double* rhs, double const* weight, double const* init, double outerTol) const 
	{
		return _jacCdisc.solve(rhs);
	}

protected:
	linalg::CompressedSparseMatrix const* const _jacC;
	linalg::FactorizableBandMatrix _jacCdisc;
};


/**
 * @brief Creates a TwoDimensionalConvectionDispersionOperator
 */
TwoDimensionalConvectionDispersionOperator::TwoDimensionalConvectionDispersionOperator() : _colPorosities(0), _stencilMemory(sizeof(active) * Weno::maxStencilSize()), 
	_wenoDerivatives(new double[Weno::maxStencilSize()]), _weno(), _linearSolver(nullptr)
{
}

TwoDimensionalConvectionDispersionOperator::~TwoDimensionalConvectionDispersionOperator() CADET_NOEXCEPT
{
	delete[] _wenoDerivatives;
	delete _linearSolver;
}

/**
 * @brief Reads parameters and allocates memory
 * @details Has to be called once before the operator is used.
 * @param [in] paramProvider Parameter provider for reading parameters
 * @param [in] nComp Number of components
 * @param [in] nCol Number of axial cells
 * @return @c true if configuration went fine, @c false otherwise
 */
bool TwoDimensionalConvectionDispersionOperator::configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int nCol, unsigned int nRad)
{
	_nComp = nComp;
	_nCol = nCol;
	_nRad = nRad;

	paramProvider.pushScope("discretization");

	// Read WENO settings and apply them
	paramProvider.pushScope("weno");
	_weno.order(paramProvider.getInt("WENO_ORDER"));
	_weno.boundaryTreatment(paramProvider.getInt("BOUNDARY_MODEL"));
	_wenoEpsilon = paramProvider.getDouble("WENO_EPS");
	paramProvider.popScope();

	// Read solver settings
	if (paramProvider.exists("LINEAR_SOLVER_BULK"))
	{
		const std::string sol = paramProvider.getString("LINEAR_SOLVER_BULK");
		if (sol == "DENSE")
			_linearSolver = new DenseDirectSolver(&_jacC);
		else if (sol == "GMRES")
			_linearSolver = new GmresSolver(&_jacC);
#ifdef UMFPACK_FOUND
		else if (sol == "UMFPACK")
			_linearSolver = new SparseDirectSolver<linalg::UMFPackSparseMatrix>(&_jacC);
#endif
#ifdef SUPERLU_FOUND
		else if (sol == "SUPERLU")
			_linearSolver = new SparseDirectSolver<linalg::SuperLUSparseMatrix>(&_jacC);
#endif
	}

	// Default to sparse solver if available (preferably UMFPACK), fall back to dense
	if (!_linearSolver)
	{
#ifdef UMFPACK_FOUND
		_linearSolver = new SparseDirectSolver<linalg::UMFPackSparseMatrix>(&_jacC);
#elif SUPERLU_FOUND
		_linearSolver = new SparseDirectSolver<linalg::SuperLUSparseMatrix>(&_jacC);
#else		
		_linearSolver = new DenseDirectSolver(&_jacC);
#endif
//		LOG(Info) << "Default to dense banded linear solver due to invalid or missing LINEAR_SOLVER_BULK setting";
	}

	paramProvider.popScope();

	_radialEdges.resize(nRad + 1);
	_radialCenters.resize(nRad);
//	_radialCentroids.resize(nRad);
	_crossSections.resize(nRad);
	_curVelocity.resize(nRad);

	setSparsityPattern();

	return _linearSolver->initialize(paramProvider, nComp, nCol, nRad, _weno);
}

/**
 * @brief Reads model parameters
 * @details Only reads parameters that do not affect model structure (e.g., discretization).
 * @param [in] unitOpIdx Unit operation id of the owning unit operation model
 * @param [in] paramProvider Parameter provider for reading parameters
 * @param [out] parameters Map in which local parameters are inserted
 * @return @c true if configuration went fine, @c false otherwise
 */
bool TwoDimensionalConvectionDispersionOperator::configure(UnitOpIdx unitOpIdx, IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters)
{
	// Read geometry parameters
	_colLength = paramProvider.getDouble("COL_LENGTH");
	_colRadius = paramProvider.getDouble("COL_RADIUS");
	readScalarParameterOrArray(_colPorosities, paramProvider, "COL_POROSITY", _nRad);

	if (_colPorosities.size() < _nRad)
		throw InvalidParameterException("Number of elements in field COL_POROSITY is less than NRAD (" + std::to_string(_nRad) + ")");

	// Read radial discretization mode and default to "EQUIDISTANT"
	paramProvider.pushScope("discretization");
	const std::string rdt = paramProvider.getString("RADIAL_DISC_TYPE");
	if (rdt == "EQUIVOLUME")
		_radialDiscretizationMode = RadialDiscretizationMode::Equivolume;
	else if (rdt == "USER_DEFINED")
	{
		_radialDiscretizationMode = RadialDiscretizationMode::UserDefined;
		readScalarParameterOrArray(_radialEdges, paramProvider, "RADIAL_COMPARTMENTS", 1);

		if (_radialEdges.size() < _nRad + 1)
			throw InvalidParameterException("Number of elements in field RADIAL_COMPARTMENTS is less than NRAD + 1 (" + std::to_string(_nRad + 1) + ")");

		registerParam1DArray(parameters, _radialEdges, [=](bool multi, unsigned int i) { return makeParamId(hashString("RADIAL_COMPARTMENTS"), unitOpIdx, CompIndep, i, BoundStateIndep, ReactionIndep, SectionIndep); });
	}
	else
		_radialDiscretizationMode = RadialDiscretizationMode::Equidistant;
	paramProvider.popScope();

	updateRadialDisc();

	// Read section dependent parameters (transport)

	// Read VELOCITY
	_velocity.clear();
	if (paramProvider.exists("VELOCITY"))
	{
		readScalarParameterOrArray(_velocity, paramProvider, "VELOCITY", _nRad);

		if ((_velocity.size() < _nRad) || (_velocity.size() % _nRad != 0))
			throw InvalidParameterException("Number of elements in field VELOCITY is not a positive multiple of NRAD (" + std::to_string(_nRad) + ")");

		registerParam2DArray(parameters, _velocity, [=](bool multi, unsigned int sec, unsigned int compartment) { return makeParamId(hashString("VELOCITY"), unitOpIdx, CompIndep, compartment, BoundStateIndep, ReactionIndep, multi ? sec : SectionIndep); }, _nRad);
	}
	readParameterMatrix(_axialDispersion, paramProvider, "COL_DISPERSION", _nComp * _nRad, 1);
	readParameterMatrix(_radialDispersion, paramProvider, "COL_DISPERSION_RADIAL", _nComp * _nRad, 1);

	if ((_axialDispersion.size() < _nComp * _nRad) || (_axialDispersion.size() % (_nComp * _nRad) != 0))
		throw InvalidParameterException("Number of elements in field COL_DISPERSION is not a positive multiple of NCOMP * NRAD (" + std::to_string(_nComp * _nRad) + ")");
	if ((_radialDispersion.size() < _nComp * _nRad) || (_radialDispersion.size() % (_nComp * _nRad) != 0))
		throw InvalidParameterException("Number of elements in field COL_DISPERSION_RADIAL is not a positive multiple of NCOMP * NRAD (" + std::to_string(_nComp * _nRad) + ")");

	// Add parameters to map
	registerParam3DArray(parameters, _axialDispersion, [=](bool multi, unsigned int sec, unsigned int compartment, unsigned int comp) { return makeParamId(hashString("COL_DISPERSION"), unitOpIdx, comp, compartment, BoundStateIndep, ReactionIndep, multi ? sec : SectionIndep); }, _nComp, _nRad);
	registerParam3DArray(parameters, _radialDispersion, [=](bool multi, unsigned int sec, unsigned int compartment, unsigned int comp) { return makeParamId(hashString("COL_DISPERSION_RADIAL"), unitOpIdx, comp, compartment, BoundStateIndep, ReactionIndep, multi ? sec : SectionIndep); }, _nComp, _nRad);
	parameters[makeParamId(hashString("COL_LENGTH"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_colLength;
	parameters[makeParamId(hashString("COL_RADIUS"), unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_colRadius;
	registerParam1DArray(parameters, _colPorosities, [=](bool multi, unsigned int i) { return makeParamId(hashString("COL_POROSITY"), unitOpIdx, CompIndep, i, BoundStateIndep, ReactionIndep, SectionIndep); });

	return true;
}

/**
 * @brief Notifies the operator that a discontinuous section transition is in progress
 * @details In addition to changing flow direction internally, if necessary, the function returns whether
 *          the flow direction has changed.
 * @param [in] t Current time point
 * @param [in] secIdx Index of the new section that is about to be integrated
 * @return @c true if flow direction has changed, otherwise @c false
 */
bool TwoDimensionalConvectionDispersionOperator::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx)
{
	bool hasChanged = false;

	if (!_velocity.empty())
	{
		// _curVelocity has already been set to the network flow rate in setFlowRates()
		// the direction of the flow (i.e., sign of _curVelocity) is given by _velocity

		active const* const dirs = getSectionDependentSlice(_velocity, _nRad, secIdx);
		if (secIdx > 0)
		{
			active const* const dirsPrev = getSectionDependentSlice(_velocity, _nRad, secIdx - 1);
			for (unsigned int i = 0; i < _nRad; ++i)
			{
				if (dirs[i] * dirsPrev[i] < 0.0)
				{
					hasChanged = true;
					break;
				}
			}
		}

		for (unsigned int i = 0; i < _nRad; ++i)
		{
			if (dirs[i] < 0.0)
				_curVelocity[i] *= -1.0;
		}
	}

	// Change the sparsity pattern if necessary
	if ((secIdx == 0) || hasChanged)
		setSparsityPattern();

	return hasChanged || (secIdx == 0);
}

/**
 * @brief Sets the flow rates for the current time section
 * @details The flow rates may change due to valve switches.
 * @param [in] in Total volumetric inlet flow rate
 * @param [in] out Total volumetric outlet flow rate
 * @param [in] colPorosity Porosity used for computing interstitial velocity from volumetric flow rate
 */
void TwoDimensionalConvectionDispersionOperator::setFlowRates(int compartment, const active& in, const active& out) CADET_NOEXCEPT
{
	_curVelocity[compartment] = in / (_crossSections[compartment] * _colPorosities[compartment]);
}

double TwoDimensionalConvectionDispersionOperator::inletFactor(unsigned int idxSec, int idxRad) const CADET_NOEXCEPT
{
	const double h = static_cast<double>(_colLength) / static_cast<double>(_nCol);
	return -std::abs(static_cast<double>(_curVelocity[idxRad])) / h;
}

const active& TwoDimensionalConvectionDispersionOperator::axialDispersion(unsigned int idxSec, int idxRad, int idxComp) const CADET_NOEXCEPT
{
	return *(getSectionDependentSlice(_axialDispersion, _nRad * _nComp, idxSec) + idxRad * _nComp + idxComp);
}

const active& TwoDimensionalConvectionDispersionOperator::radialDispersion(unsigned int idxSec, int idxRad, int idxComp) const CADET_NOEXCEPT
{
	return *(getSectionDependentSlice(_radialDispersion, _nRad * _nComp, idxSec) + idxRad * _nComp + idxComp);
}

/**
 * @brief Computes the residual of the transport equations
 * @param [in] t Current time point
 * @param [in] secIdx Index of the current section
 * @param [in] timeFactor Used to compute parameter derivatives with respect to section length (nominal value should always be 1.0)
 * @param [in] y Pointer to unit operation's state vector
 * @param [in] yDot Pointer to unit operation's time derivative state vector
 * @param [out] res Pointer to unit operation's residual vector
 * @param [in] wantJac Determines whether the Jacobian is computed or not
 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
 */
int TwoDimensionalConvectionDispersionOperator::residual(double t, unsigned int secIdx, double timeFactor, double const* y, double const* yDot, double* res, bool wantJac)
{
	if (wantJac)
		return residualImpl<double, double, double, true>(t, secIdx, timeFactor, y, yDot, res);
	else
		return residualImpl<double, double, double, false>(t, secIdx, timeFactor, y, yDot, res);
}

int TwoDimensionalConvectionDispersionOperator::residual(double t, unsigned int secIdx, double timeFactor, active const* y, double const* yDot, active* res, bool wantJac)
{
	if (wantJac)
		return residualImpl<active, active, double, true>(t, secIdx, timeFactor, y, yDot, res);
	else
		return residualImpl<active, active, double, false>(t, secIdx, timeFactor, y, yDot, res);
}

int TwoDimensionalConvectionDispersionOperator::residual(const active& t, unsigned int secIdx, const active& timeFactor, double const* y, double const* yDot, active* res, bool wantJac)
{
	if (wantJac)
		return residualImpl<double, active, active, true>(t, secIdx, timeFactor, y, yDot, res);
	else
		return residualImpl<double, active, active, false>(t, secIdx, timeFactor, y, yDot, res);
}

int TwoDimensionalConvectionDispersionOperator::residual(const active& t, unsigned int secIdx, const active& timeFactor, active const* y, double const* yDot, active* res, bool wantJac)
{
	if (wantJac)
		return residualImpl<active, active, active, true>(t, secIdx, timeFactor, y, yDot, res);
	else
		return residualImpl<active, active, active, false>(t, secIdx, timeFactor, y, yDot, res);
}

template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
int TwoDimensionalConvectionDispersionOperator::residualImpl(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* y, double const* yDot, ResidualType* res)
{
	if (wantJac)
	{
		// Reset Jacobian
		_jacC.setAll(0.0);
	}

	// Handle convection, axial dispersion (WENO)
	const ParamType h = static_cast<ParamType>(_colLength) / static_cast<double>(_nCol);
	for (unsigned int i = 0; i < _nRad; ++i)
	{
		active const* const d_c = getSectionDependentSlice(_axialDispersion, _nRad * _nComp, secIdx) + i * _nComp;

		convdisp::FlowParameters<ParamType> fp{
			static_cast<ParamType>(_curVelocity[i]),
			d_c,
			h,
			_wenoDerivatives,
			&_weno,
			&_stencilMemory,
			_wenoEpsilon,
			static_cast<int>(_nComp * _nRad),  // Stride between two cells
			_nComp,
			_nCol,
			_nComp * i,                        // Offset to the first component of the inlet DOFs in the local state vector
			_nComp * (_nRad + i)               // Offset to the first component of the first bulk cell in the local state vector
		};

		if (wantJac)
			convdisp::residualKernel<StateType, ResidualType, ParamType, linalg::BandedSparseRowIterator, true>(convdisp::SimulationTime<ParamType>{t, secIdx, timeFactor}, y, yDot, res, _jacC.row(i * _nComp), fp);
		else
			convdisp::residualKernel<StateType, ResidualType, ParamType, linalg::BandedSparseRowIterator, false>(convdisp::SimulationTime<ParamType>{t, secIdx, timeFactor}, y, yDot, res, _jacC.row(i * _nComp), fp);
	}

	// Handle radial dispersion
	if (cadet_unlikely(_nRad <= 1))
		return 0;

	active const* const d_rho = getSectionDependentSlice(_radialDispersion, _nRad * _nComp, secIdx);
	const unsigned int offsetC = _nRad * _nComp;
	for (unsigned int col = 0; col < _nCol; ++col)
	{
		const unsigned int offsetColBlock = col * _nRad * _nComp;
		ResidualType* const resColBlock = res + offsetC + offsetColBlock;
		StateType const* const yColBlock = y + offsetC + offsetColBlock;

		for (unsigned int rad = 0; rad < _nRad; ++rad)
		{
			const unsigned int offsetToRadBlock = rad * _nComp;
			const unsigned int offsetColRadBlock = offsetColBlock + offsetToRadBlock;
			ResidualType* const resColRadBlock = resColBlock + offsetToRadBlock;
			StateType const* const yColRadBlock = yColBlock + offsetToRadBlock;

			active const* const localD_rho = d_rho + rad * _nComp;
			const ParamType invEps = 1.0 / static_cast<ParamType>(_colPorosities[rad]);
			const ParamType deltaRho = static_cast<ParamType>(_radialEdges[rad+1]) - static_cast<ParamType>(_radialEdges[rad]);
			const ParamType rho = static_cast<ParamType>(_radialCenters[rad]);
//			const ParamType rhoCentroid = static_cast<ParamType>(_radialCentroids[rad]);

			for (unsigned int comp = 0; comp < _nComp; ++comp)
			{
				const unsigned int offsetCur = offsetColRadBlock + comp;
				ResidualType* const resCur = resColRadBlock + comp;
				StateType const* const yCur = yColRadBlock + comp;
				const ParamType curFluxCoeff = static_cast<ParamType>(_colPorosities[rad]) * static_cast<ParamType>(localD_rho[comp]);

				if (cadet_unlikely(curFluxCoeff <= 0.0))
					continue;

				if (cadet_likely(rad > 0))
				{
					// Flow from inner cell
					const ParamType centerDiffRho = rho - static_cast<ParamType>(_radialCenters[rad-1]);
//					const ParamType centerDiffRho = rhoCentroid - static_cast<ParamType>(_radialCentroids[rad-1]);
					const ParamType innerFluxCoeff = static_cast<ParamType>(_colPorosities[rad-1]) * static_cast<ParamType>(localD_rho[-static_cast<int>(_nComp) + static_cast<int>(comp)]);

					if (cadet_likely(innerFluxCoeff > 0.0))
					{
//						const ParamType fluxCoeff = (centerDiffRho * curFluxCoeff * innerFluxCoeff) / (curFluxCoeff * (static_cast<ParamType>(_radialEdges[rad]) - static_cast<ParamType>(_radialCentroids[rad-1])) + innerFluxCoeff * (rhoCentroid - static_cast<ParamType>(_radialEdges[rad])));
						const ParamType fluxCoeff = (centerDiffRho * curFluxCoeff * innerFluxCoeff) / (curFluxCoeff * (static_cast<ParamType>(_radialEdges[rad]) - static_cast<ParamType>(_radialCenters[rad-1])) + innerFluxCoeff * (rho - static_cast<ParamType>(_radialEdges[rad])));
//						const ParamType fluxCoeff = static_cast<ParamType>(_colPorosities[0]) * static_cast<ParamType>(d_rho[0]);
						const ParamType finalFactor = invEps * (static_cast<ParamType>(_radialEdges[rad]) * fluxCoeff / centerDiffRho) / (rho * deltaRho);

						*resCur += finalFactor * (yCur[0] - yCur[-static_cast<int>(_nComp)]);
						if (wantJac)
						{
							_jacC.centered(offsetCur, 0) += static_cast<double>(finalFactor);
							_jacC.centered(offsetCur, -static_cast<int>(_nComp)) -= static_cast<double>(finalFactor);
						}
					}
				}

				if (cadet_likely(rad < _nRad - 1))
				{
					// Flow from outer cell
					const ParamType centerDiffRho = static_cast<ParamType>(_radialCenters[rad+1]) - rho;
//					const ParamType centerDiffRho = static_cast<ParamType>(_radialCentroids[rad+1]) - rhoCentroid;
					const ParamType outerFluxCoeff = static_cast<ParamType>(_colPorosities[rad+1]) * static_cast<ParamType>(localD_rho[_nComp + comp]);

					if (cadet_likely(outerFluxCoeff > 0.0))
					{
//						const ParamType fluxCoeff = (centerDiffRho * curFluxCoeff * outerFluxCoeff) / (curFluxCoeff * (static_cast<ParamType>(_radialCentroids[rad+1] - static_cast<ParamType>(_radialEdges[rad+1]))) + outerFluxCoeff * (static_cast<ParamType>(_radialEdges[rad+1]) - rhoCentroid));
						const ParamType fluxCoeff = (centerDiffRho * curFluxCoeff * outerFluxCoeff) / (curFluxCoeff * (static_cast<ParamType>(_radialCenters[rad+1] - static_cast<ParamType>(_radialEdges[rad+1]))) + outerFluxCoeff * (static_cast<ParamType>(_radialEdges[rad+1]) - rho));
//						const ParamType fluxCoeff = static_cast<ParamType>(_colPorosities[0]) * static_cast<ParamType>(d_rho[0]);
						const ParamType finalFactor = invEps * (static_cast<ParamType>(_radialEdges[rad+1]) * fluxCoeff / centerDiffRho) / (rho * deltaRho);

						*resCur += finalFactor * (yCur[0] - yCur[_nComp]);
						if (wantJac)
						{
							_jacC.centered(offsetCur, 0) += static_cast<double>(finalFactor);
							_jacC.centered(offsetCur, _nComp) -= static_cast<double>(finalFactor);
						}
					}
				}
			}
		}
	}

	return 0;
}

void TwoDimensionalConvectionDispersionOperator::setSparsityPattern()
{
	// Note that we have to increase the lower non-zeros by 1 because the WENO stencil is applied to the
	// right cell face (lower + 1 + upper) and to the left cell face (shift the stencil by -1 because influx of cell i
	// is outflux of cell i-1)
	// We also have to make sure that there's at least one sub and super diagonal for the dispersion term
	const unsigned int lowerNonZeros = std::max(_weno.lowerBandwidth() + 1u, 1u);
	const unsigned int upperNonZeros = std::max(_weno.upperBandwidth(), 1u);
	// Total number of non-zeros per row is WENO stencil (lowerNonZeros + 1u + upperNonZeros) + radial dispersion (2)
	cadet::linalg::SparsityPattern pattern(_nComp * _nCol * _nRad, lowerNonZeros + 1u + upperNonZeros + 2u);

	// Handle convection, axial dispersion (WENO)
	for (unsigned int i = 0; i < _nRad; ++i)
		cadet::model::parts::convdisp::sparsityPattern(pattern.row(i * _nComp), _nComp, _nCol, _nComp * _nRad, static_cast<double>(_curVelocity[i]), _weno);

	// Handle radial dispersion
	if (_nRad > 1)
	{
		for (unsigned int col = 0; col < _nCol; ++col)
		{
			const unsigned int idxColBlock = col * _nRad * _nComp;

			// First and last cell have only one term
			for (unsigned int comp = 0; comp < _nComp; ++comp)
			{
				// Radial cell 1
				pattern.add(idxColBlock + comp, idxColBlock + comp + _nComp);

				// Radial cell nRad - 1
				pattern.add(idxColBlock + (_nRad - 1) * _nComp + comp, idxColBlock + comp + (_nRad - 2) * _nComp);
			}

			// Middle cells have two terms
			for (unsigned int rad = 1; rad < _nRad - 1; ++rad)
			{
				const unsigned int idxColRadBlock = idxColBlock + rad * _nComp;
				for (unsigned int comp = 0; comp < _nComp; ++comp)
				{
					const unsigned int idxCur = idxColRadBlock + comp;
					pattern.add(idxCur, idxCur + _nComp);
					pattern.add(idxCur, idxCur - _nComp);
				}
			}
		}
	}

	_jacC.assignPattern(pattern);
	_linearSolver->setSparsityPattern(pattern);
}

/**
 * @brief Multiplies the time derivative Jacobian @f$ \frac{\partial F}{\partial \dot{y}}\left(t, y, \dot{y}\right) @f$ with a given vector
 * @details The operation @f$ z = \frac{\partial F}{\partial \dot{y}} x @f$ is performed.
 *          The matrix-vector multiplication is transformed matrix-free (i.e., no matrix is explicitly formed).
 *          
 *          Note that this function only performs multiplication with the Jacobian of the (axial) transport equations.
 *          The input vectors are assumed to point to the beginning (including inlet DOFs) of the respective unit operation's arrays.
 * @param [in] simTime Simulation time information (time point, section index, pre-factor of time derivatives)
 * @param [in] sDot Vector @f$ x @f$ that is transformed by the Jacobian @f$ \frac{\partial F}{\partial \dot{y}} @f$
 * @param [out] ret Vector @f$ z @f$ which stores the result of the operation
 */
void TwoDimensionalConvectionDispersionOperator::multiplyWithDerivativeJacobian(const SimulationTime& simTime, double const* sDot, double* ret) const
{
	double* localRet = ret + _nComp * _nRad;
	double const* localSdot = sDot + _nComp * _nRad;
	for (unsigned int i = 0; i < _nCol * _nComp * _nRad; ++i)
		localRet[i] = simTime.timeFactor * localSdot[i];
}

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
void TwoDimensionalConvectionDispersionOperator::assembleDiscretizedJacobian(double alpha, double timeFactor)
{
	_linearSolver->assembleDiscretizedJacobian(alpha, timeFactor);
}

/**
 * @brief Assembles and factorizes the time discretized Jacobian
 * @details See assembleDiscretizedJacobian() for assembly of the time discretized Jacobian.
 * @param [in] alpha Factor in front of @f$ \frac{\partial F}{\partial \dot{y}} @f$
 * @param [in] timeFactor Factor which is premultiplied to the time derivatives originating from time transformation
 * @return @c true if factorization went fine, otherwise @c false
 */
bool TwoDimensionalConvectionDispersionOperator::assembleAndFactorizeDiscretizedJacobian(double alpha, double timeFactor)
{
	assembleDiscretizedJacobian(alpha, timeFactor);
	return _linearSolver->factorize();
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
bool TwoDimensionalConvectionDispersionOperator::solveDiscretizedJacobian(double* rhs, double const* weight, double const* init, double outerTol) const
{
	return _linearSolver->solveDiscretizedJacobian(rhs, weight, init, outerTol);
}

/**
 * @brief Solves a system with the time derivative Jacobian and given right hand side
 * @details Note that the given right hand side vector @p rhs is not shifted by the inlet DOFs. That
 *          is, it is assumed to point directly to the first axial DOF.
 * @param [in] simTime Simulation time information (time point, section index, pre-factor of time derivatives)
 * @param [in,out] rhs On entry, right hand side. On exit, solution of the system.
 * @return @c true if the system was solved correctly, @c false otherwise
 */
bool TwoDimensionalConvectionDispersionOperator::solveTimeDerivativeSystem(const SimulationTime& simTime, double* const rhs)
{
	const double invFactor = 1.0 / simTime.timeFactor;
	for (unsigned int i = 0; i < _nComp * _nRad * _nCol; ++i)
		rhs[i] *= invFactor;

	return true;
}

void TwoDimensionalConvectionDispersionOperator::setEquidistantRadialDisc()
{
	const active h = _colRadius / _nRad;
	const double pi = 3.1415926535897932384626434;

	_radialEdges[0] = 0.0;
	for (unsigned int r = 0; r < _nRad; ++r)
	{
		// Set last edge to _colRadius for exact geometry
		if (r == _nRad - 1)
			_radialEdges[r+1] = _colRadius;
		else
			_radialEdges[r+1] = h * (r + 1);

		_crossSections[r] = pi * (_radialEdges[r+1] * _radialEdges[r+1] - _radialEdges[r] * _radialEdges[r]);
		_radialCenters[r] = 0.5 * (_radialEdges[r+1] + _radialEdges[r]);
//		_radialCentroids[r] = 2.0 / 3.0 * ( _radialEdges[r+1] * _radialEdges[r+1] / (_radialEdges[r+1] + _radialEdges[r]) + _radialEdges[r]);
	}
}

void TwoDimensionalConvectionDispersionOperator::setEquivolumeRadialDisc()
{
	const active volPerCompartment = _colRadius * _colRadius / _nRad;
	const double pi = 3.1415926535897932384626434;

	_radialEdges[0] = 0.0;
	for (unsigned int r = 0; r < _nRad; ++r)
	{
		// Set last edge to _colRadius for exact geometry
		if (r == _nRad - 1)
			_radialEdges[r+1] = _colRadius;
		else
			_radialEdges[r+1] = sqrt(volPerCompartment + _radialEdges[r] * _radialEdges[r]);

		_crossSections[r] = pi * volPerCompartment;
		_radialCenters[r] = 0.5 * (_radialEdges[r+1] + _radialEdges[r]);
//		_radialCentroids[r] = 2.0 / 3.0 * ( _radialEdges[r+1] * _radialEdges[r+1] / (_radialEdges[r+1] + _radialEdges[r]) + _radialEdges[r]);
	}
}

void TwoDimensionalConvectionDispersionOperator::setUserdefinedRadialDisc()
{
	const double pi = 3.1415926535897932384626434;
	for (unsigned int r = 0; r < _nRad; ++r)
	{
		_crossSections[r] = pi * (_radialEdges[r+1] * _radialEdges[r+1] - _radialEdges[r] * _radialEdges[r]);
		_radialCenters[r] = 0.5 * (_radialEdges[r+1] + _radialEdges[r]);
//		_radialCentroids[r] = 2.0 / 3.0 * ( _radialEdges[r+1] * _radialEdges[r+1] / (_radialEdges[r+1] + _radialEdges[r]) + _radialEdges[r]);
	}
}

void TwoDimensionalConvectionDispersionOperator::updateRadialDisc()
{
	if (_radialDiscretizationMode == RadialDiscretizationMode::Equidistant)
		setEquidistantRadialDisc();
	else if (_radialDiscretizationMode == RadialDiscretizationMode::Equivolume)
		setEquivolumeRadialDisc();
	else if (_radialDiscretizationMode == RadialDiscretizationMode::UserDefined)
		setUserdefinedRadialDisc();
}


}  // namespace parts

}  // namespace model

}  // namespace cadet
