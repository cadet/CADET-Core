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

#include "model/LumpedRateModelWithPores.hpp"
#include "BindingModelFactory.hpp"
#include "ParamReaderHelper.hpp"
#include "cadet/Exceptions.hpp"
#include "cadet/ExternalFunction.hpp"
#include "cadet/SolutionRecorder.hpp"
#include "ConfigurationHelper.hpp"
#include "model/BindingModel.hpp"
#include "linalg/DenseMatrix.hpp"
#include "linalg/BandMatrix.hpp"
#include "linalg/Norms.hpp"

#include "Stencil.hpp"
#include "Weno.hpp"
#include "AdUtils.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <algorithm>
#include <functional>

#include "ParallelSupport.hpp"
#ifdef CADET_PARALLELIZE
	#include <tbb/tbb.h>
#endif

namespace cadet
{

namespace model
{

int schurComplementMultiplierLRMPores(void* userData, double const* x, double* z)
{
	LumpedRateModelWithPores* const lrm = static_cast<LumpedRateModelWithPores*>(userData);
	return lrm->schurComplementMatrixVector(x, z);
}


LumpedRateModelWithPores::LumpedRateModelWithPores(UnitOpIdx unitOpIdx) : UnitOperationBase(unitOpIdx),
	_jacInlet(), _analyticJac(true), _jacobianAdDirs(0), _factorizeJacobian(false), _tempState(nullptr),
	_initC(0), _initCp(0), _initQ(0), _initState(0), _initStateDot(0)
{
}

LumpedRateModelWithPores::~LumpedRateModelWithPores() CADET_NOEXCEPT
{
	delete[] _tempState;

	delete[] _disc.nBound;
	delete[] _disc.boundOffset;
}

unsigned int LumpedRateModelWithPores::numDofs() const CADET_NOEXCEPT
{
	// Column bulk DOFs: nCol * nComp
	// Particle DOFs: nCol particles each having nComp (liquid phase) + sum boundStates (solid phase) DOFs
	// Flux DOFs: nCol * nComp (as many as column bulk DOFs)
	// Inlet DOFs: nComp
	return _disc.nCol * (3 * _disc.nComp + _disc.strideBound) + _disc.nComp;
}

unsigned int LumpedRateModelWithPores::numPureDofs() const CADET_NOEXCEPT
{
	// Column bulk DOFs: nCol * nComp
	// Particle DOFs: nCol particles each having nComp (liquid phase) + sum boundStates (solid phase) DOFs
	// Flux DOFs: nCol * nComp (as many as column bulk DOFs)
	return _disc.nCol * (3 * _disc.nComp + _disc.strideBound);
}


bool LumpedRateModelWithPores::usesAD() const CADET_NOEXCEPT
{
#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	// We always need AD if we want to check the analytical Jacobian
	return true;
#else
	// We only need AD if we are not computing the Jacobian analytically
	return !_analyticJac;
#endif
}

bool LumpedRateModelWithPores::configureModelDiscretization(IParameterProvider& paramProvider, IConfigHelper& helper)
{
	// ==== Read discretization
	_disc.nComp = paramProvider.getInt("NCOMP");

	paramProvider.pushScope("discretization");

	_disc.nCol = paramProvider.getInt("NCOL");

	const std::vector<int> nBound = paramProvider.getIntArray("NBOUND");
	if (nBound.size() < _disc.nComp)
		throw InvalidParameterException("Field NBOUND contains too few elements (NCOMP = " + std::to_string(_disc.nComp) + " required)");

	_disc.nBound = new unsigned int[_disc.nComp];
	std::copy_n(nBound.begin(), _disc.nComp, _disc.nBound);

	// Precompute offsets and total number of bound states (DOFs in solid phase)
	_disc.boundOffset = new unsigned int[_disc.nComp];
	_disc.boundOffset[0] = 0;
	for (unsigned int i = 1; i < _disc.nComp; ++i)
	{
		_disc.boundOffset[i] = _disc.boundOffset[i-1] + _disc.nBound[i-1];
	}
	_disc.strideBound = _disc.boundOffset[_disc.nComp-1] + _disc.nBound[_disc.nComp - 1];

	// Determine whether analytic Jacobian should be used but don't set it right now.
	// We need to setup Jacobian matrices first.
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	const bool analyticJac = paramProvider.getBool("USE_ANALYTIC_JACOBIAN");
#else
	const bool analyticJac = false;
#endif

	// Initialize and configure GMRES for solving the Schur-complement
	_gmres.initialize(_disc.nCol * _disc.nComp, paramProvider.getInt("MAX_KRYLOV"), linalg::toOrthogonalization(paramProvider.getInt("GS_TYPE")), paramProvider.getInt("MAX_RESTARTS"));
	_gmres.matrixVectorMultiplier(&schurComplementMultiplierLRMPores, this);
	_schurSafety = paramProvider.getDouble("SCHUR_SAFETY");

	// Allocate space for initial conditions
	_initC.resize(_disc.nComp);
	_initCp.resize(_disc.nComp);
	_initQ.resize(_disc.strideBound);

	paramProvider.popScope();

	const bool transportSuccess = _convDispOp.configureModelDiscretization(paramProvider, _disc.nComp, _disc.nCol);

	// Allocate memory
	Indexer idxr(_disc);

	_jacInlet.resize(_disc.nComp);

	_jacPdisc.resize(_disc.nCol * (_disc.nComp + _disc.strideBound), _disc.nComp + _disc.strideBound - 1, _disc.nComp + _disc.strideBound - 1);
	_jacP.resize(_disc.nCol * (_disc.nComp + _disc.strideBound), _disc.nComp + _disc.strideBound - 1, _disc.nComp + _disc.strideBound - 1);

	_jacPF.resize(_disc.nComp * _disc.nCol);
	_jacFP.resize(_disc.nComp * _disc.nCol);

	_jacCF.resize(_disc.nComp * _disc.nCol);
	_jacFC.resize(_disc.nComp * _disc.nCol);

	// Set whether analytic Jacobian is used
	useAnalyticJacobian(analyticJac);

	// ==== Construct and configure binding model
	delete _binding;

	_binding = helper.createBindingModel(paramProvider.getString("ADSORPTION_MODEL"));
	if (!_binding)
		throw InvalidParameterException("Unknown binding model " + paramProvider.getString("ADSORPTION_MODEL"));

	const bool bindingConfSuccess = _binding->configureModelDiscretization(paramProvider, _disc.nComp, _disc.nBound, _disc.boundOffset);

	// setup the memory for tempState based on state vector or memory needed for consistent initialization of isotherms, whichever is larger
	unsigned int size = numDofs();
	if (_binding->requiresWorkspace())
	{
		// Required memory (number of doubles) for nonlinear solvers
		const unsigned int requiredMem = (_binding->workspaceSize() + sizeof(double) - 1) / sizeof(double) * _disc.nCol;
		if (requiredMem > size)
		{
			size = requiredMem;
		}
	}
	_tempState = new double[size];

	return transportSuccess && bindingConfSuccess;
}

bool LumpedRateModelWithPores::configure(IParameterProvider& paramProvider)
{
	_parameters.clear();

	const bool transportSuccess = _convDispOp.configure(_unitOpIdx, paramProvider, _parameters);

	// Read geometry parameters
	_colPorosity = paramProvider.getDouble("COL_POROSITY");
	_parRadius = paramProvider.getDouble("PAR_RADIUS");
	_parPorosity = paramProvider.getDouble("PAR_POROSITY");

	// Read vectorial parameters (which may also be section dependent; transport)
	readParameterMatrix(_filmDiffusion, paramProvider, "FILM_DIFFUSION", _disc.nComp, 1);

	if (paramProvider.exists("PORE_ACCESSIBILITY"))
		readParameterMatrix(_poreAccessFactor, paramProvider, "PORE_ACCESSIBILITY", _disc.nComp, 1);
	else
		_poreAccessFactor = std::vector<cadet::active>(_disc.nComp, 1.0);

	// Add parameters to map
	_parameters[makeParamId(hashString("COL_POROSITY"), _unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &_colPorosity;
	_parameters[makeParamId(hashString("PAR_RADIUS"), _unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &_parRadius;
	_parameters[makeParamId(hashString("PAR_POROSITY"), _unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &_parPorosity;

	registerComponentSectionDependentParam(hashString("FILM_DIFFUSION"), _parameters, _filmDiffusion, _unitOpIdx, _disc.nComp);
	registerComponentSectionDependentParam(hashString("PORE_ACCESSIBILITY"), _parameters, _poreAccessFactor, _unitOpIdx, _disc.nComp);

	// Register initial conditions parameters
	for (unsigned int i = 0; i < _disc.nComp; ++i)
	{
		_parameters[makeParamId(hashString("INIT_C"), _unitOpIdx, i, BoundPhaseIndep, ReactionIndep, SectionIndep)] = _initC.data() + i;
		_parameters[makeParamId(hashString("INIT_CP"), _unitOpIdx, i, BoundPhaseIndep, ReactionIndep, SectionIndep)] = _initCp.data() + i;
	}

	if (_binding)
	{
		std::vector<ParameterId> initParams(_disc.strideBound);
		_binding->fillBoundPhaseInitialParameters(initParams.data(), _unitOpIdx);

		for (unsigned int i = 0; i < _disc.strideBound; ++i)
			_parameters[initParams[i]] = _initQ.data() + i;
	}

	// Reconfigure binding model
	if (_binding && paramProvider.exists("adsorption") && _binding->requiresConfiguration())
	{
		paramProvider.pushScope("adsorption");
		const bool bindingConfSuccess = _binding->configure(paramProvider, _unitOpIdx);
		paramProvider.popScope();

		return transportSuccess && bindingConfSuccess;
	}

	return transportSuccess;
}

void LumpedRateModelWithPores::useAnalyticJacobian(const bool analyticJac)
{
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	_analyticJac = analyticJac;
	if (!_analyticJac)
		// We need as many directions as the highest bandwidth of the diagonal blocks:
		// The bandwidth of the column block depends on the size of the WENO stencil, whereas
		// the bandwidth of the particle blocks are given by the number of components and bound states.
		_jacobianAdDirs = std::max(_convDispOp.requiredADdirs(), _jacP.stride());
	else
		_jacobianAdDirs = 0;
#else
	_analyticJac = false;
	// We need as many directions as the highest bandwidth of the diagonal blocks:
	// The bandwidth of the column block depends on the size of the WENO stencil, whereas
	// the bandwidth of the particle blocks are given by the number of components and bound states.
	_jacobianAdDirs = std::max(_convDispOp.requiredADdirs(), _jacP.stride());
#endif
}

void LumpedRateModelWithPores::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active* const adRes, active* const adY, unsigned int adDirOffset)
{
	// Setup flux Jacobian blocks at the beginning of the simulation or in case of
	// section dependent film or particle diffusion coefficients
	if ((secIdx == 0) || (_filmDiffusion.size() > _disc.nComp))
		assembleOffdiagJac(t, secIdx);

	Indexer idxr(_disc);

	// ConvectionDispersionOperator tells us whether flow direction has changed
	if (!_convDispOp.notifyDiscontinuousSectionTransition(t, secIdx, adRes, adY, adDirOffset))
		return;

	// Setup the matrix connecting inlet DOFs to first column cells
	_jacInlet.clear();
	const double h = static_cast<double>(_convDispOp.columnLength()) / static_cast<double>(_disc.nCol);
	const double u = static_cast<double>(_convDispOp.currentVelocity());

	if (u >= 0.0)
	{
		// Forwards flow

		// Place entries for inlet DOF to first column cell conversion
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			_jacInlet.addElement(comp * idxr.strideColComp(), comp, -u / h);
	}
	else
	{
		// Backwards flow

		// Place entries for inlet DOF to last column cell conversion
		const unsigned int offset = (_disc.nCol - 1) * idxr.strideColCell();
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			_jacInlet.addElement(offset + comp * idxr.strideColComp(), comp, u / h);
	}
}

void LumpedRateModelWithPores::setFlowRates(const active& in, const active& out) CADET_NOEXCEPT
{
	_convDispOp.setFlowRates(in, out, _colPorosity);
}

void LumpedRateModelWithPores::reportSolution(ISolutionRecorder& recorder, double const* const solution) const
{
	Exporter expr(_disc, solution);
	recorder.beginUnitOperation(_unitOpIdx, *this, expr);
	recorder.endUnitOperation();
}

void LumpedRateModelWithPores::reportSolutionStructure(ISolutionRecorder& recorder) const
{
	Exporter expr(_disc, nullptr);
	recorder.unitOperationStructure(_unitOpIdx, *this, expr);
}


unsigned int LumpedRateModelWithPores::requiredADdirs() const CADET_NOEXCEPT
{
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	return _jacobianAdDirs;
#else
	// If CADET_CHECK_ANALYTIC_JACOBIAN is active, we always need the AD directions for the Jacobian
	return std::max(_convDispOp.requiredADdirs(), _jacP.stride());
#endif
}

void LumpedRateModelWithPores::prepareADvectors(active* const adRes, active* const adY, unsigned int adDirOffset) const
{
	// Early out if AD is disabled
	if (!adY)
		return;

	Indexer idxr(_disc);

	// Get bandwidths of blocks
	const unsigned int lowerParBandwidth = _jacP.lowerBandwidth();
	const unsigned int upperParBandwidth = _jacP.upperBandwidth();

	// Column block	
	_convDispOp.prepareADvectors(adRes, adY, adDirOffset);

	// Particle block
	ad::prepareAdVectorSeedsForBandMatrix(adY + idxr.offsetCp(), adDirOffset, idxr.strideParBlock() * _disc.nCol, lowerParBandwidth, upperParBandwidth, lowerParBandwidth);
}

/**
 * @brief Extracts the system Jacobian from band compressed AD seed vectors
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
void LumpedRateModelWithPores::extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset)
{
	Indexer idxr(_disc);

	// Column
	_convDispOp.extractJacobianFromAD(adRes, adDirOffset);

	// Particles
	ad::extractBandedJacobianFromAd(adRes + idxr.offsetCp(), adDirOffset, _jacP.lowerBandwidth(), _jacP);
}

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN

/**
 * @brief Compares the analytical Jacobian with a Jacobian derived by AD
 * @details The analytical Jacobian is assumed to be stored in the corresponding band matrices.
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
void LumpedRateModelWithPores::checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const
{
	Indexer idxr(_disc);

	LOG(Debug) << "AD dir offset: " << adDirOffset << " DiagDirCol: " << _convDispOp.jacobian().lowerBandwidth() << " DiagDirPar: " << _jacP.lowerBandwidth();

	// Column
	const double maxDiffCol = _convDispOp.checkAnalyticJacobianAgainstAd(adRes, adDirOffset);

	// Particles
	const double maxDiffPar = ad::compareBandedJacobianWithAd(adRes + idxr.offsetCp(), adDirOffset, _jacP.lowerBandwidth(), _jacP);
	LOG(Debug) << "-> Par block diff: " << maxDiffPar;
}

#endif

int LumpedRateModelWithPores::residual(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double* const res)
{
	BENCH_SCOPE(_timerResidual);

	// Evaluate residual do not compute Jacobian or parameter sensitivities
	return residualImpl<double, double, double, false>(t, secIdx, timeFactor, y, yDot, res);
}

int LumpedRateModelWithPores::residualWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, double* const res, active* const adRes, active* const adY, unsigned int adDirOffset)
{
	BENCH_SCOPE(_timerResidual);

	// Evaluate residual, use AD for Jacobian if required but do not evaluate parameter derivatives
	return residual(t, secIdx, timeFactor, y, yDot, res, adRes, adY, adDirOffset, true, false);
}

int LumpedRateModelWithPores::residual(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, double* const res, 
	active* const adRes, active* const adY, unsigned int adDirOffset, bool updateJacobian, bool paramSensitivity)
{
	if (updateJacobian)
	{
		_factorizeJacobian = true;

#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
		if (_analyticJac)
		{
			if (paramSensitivity)
			{
				const int retCode = residualImpl<double, active, active, true>(t, secIdx, timeFactor, y, yDot, adRes);

				// Copy AD residuals to original residuals vector
				if (res)
					ad::copyFromAd(adRes, res, numDofs());

				return retCode;
			}
			else
				return residualImpl<double, double, double, true>(static_cast<double>(t), secIdx, static_cast<double>(timeFactor), y, yDot, res);
		}
		else
		{
			// Compute Jacobian via AD

			// Copy over state vector to AD state vector (without changing directional values to keep seed vectors)
			// and initalize residuals with zero (also resetting directional values)
			ad::copyToAd(y, adY, numDofs());
			// @todo Check if this is necessary
			ad::resetAd(adRes, numDofs());

			// Evaluate with AD enabled
			int retCode = 0;
			if (paramSensitivity)
				retCode = residualImpl<active, active, active, false>(t, secIdx, timeFactor, adY, yDot, adRes);
			else
				retCode = residualImpl<active, active, double, false>(static_cast<double>(t), secIdx, static_cast<double>(timeFactor), adY, yDot, adRes);

			// Copy AD residuals to original residuals vector
			if (res)
				ad::copyFromAd(adRes, res, numDofs());

			// Extract Jacobian
			extractJacobianFromAD(adRes, adDirOffset);

			return retCode;
		}
#else
		// Compute Jacobian via AD

		// Copy over state vector to AD state vector (without changing directional values to keep seed vectors)
		// and initalize residuals with zero (also resetting directional values)
		ad::copyToAd(y, adY, numDofs());
		// @todo Check if this is necessary
		ad::resetAd(adRes, numDofs());

		// Evaluate with AD enabled
		int retCode = 0;
		if (paramSensitivity)
			retCode = residualImpl<active, active, active, false>(t, secIdx, timeFactor, adY, yDot, adRes);
		else
			retCode = residualImpl<active, active, double, false>(static_cast<double>(t), secIdx, static_cast<double>(timeFactor), adY, yDot, adRes);

		// Only do comparison if we have a residuals vector (which is not always the case)
		if (res)
		{
			// Evaluate with analytical Jacobian which is stored in the band matrices
			retCode = residualImpl<double, double, double, true>(static_cast<double>(t), secIdx, static_cast<double>(timeFactor), y, yDot, res);

			// Compare AD with anaytic Jacobian
			checkAnalyticJacobianAgainstAd(adRes, adDirOffset);
		}

		// Extract Jacobian
		extractJacobianFromAD(adRes, adDirOffset);

		return retCode;
#endif
	}
	else
	{
		if (paramSensitivity)
		{
			// Initalize residuals with zero
			// @todo Check if this is necessary
			ad::resetAd(adRes, numDofs());

			const int retCode = residualImpl<double, active, active, false>(t, secIdx, timeFactor, y, yDot, adRes);

			// Copy AD residuals to original residuals vector
			if (res)
				ad::copyFromAd(adRes, res, numDofs());

			return retCode;
		}
		else
			return residualImpl<double, double, double, false>(static_cast<double>(t), secIdx, static_cast<double>(timeFactor), y, yDot, res);
	}
}

template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
int LumpedRateModelWithPores::residualImpl(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* const y, double const* const yDot, ResidualType* const res)
{
	// Reset Jacobian
	if (wantJac)
		_jacP.setAll(0.0);

	BENCH_START(_timerResidualPar);

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(size_t(0), size_t(_disc.nCol + 1), [&](size_t pblk)
#else
	for (unsigned int pblk = 0; pblk < _disc.nCol + 1; ++pblk)
#endif
	{
		if (cadet_unlikely(pblk == 0))
			_convDispOp.residual(t, secIdx, timeFactor, y, yDot, res, wantJac);
		else
			residualParticle<StateType, ResidualType, ParamType, wantJac>(t, pblk-1, secIdx, timeFactor, y, yDot, res);
	} CADET_PARFOR_END;

	BENCH_STOP(_timerResidualPar);

	residualFlux<StateType, ResidualType, ParamType>(t, secIdx, y, yDot, res);

	// Handle inlet DOFs, which are simply copied to res
	for (unsigned int i = 0; i < _disc.nComp; ++i)
	{
		res[i] = y[i];
	}

	return 0;
}

template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
int LumpedRateModelWithPores::residualParticle(const ParamType& t, unsigned int colCell, unsigned int secIdx, const ParamType& timeFactor, StateType const* yBase, double const* yDotBase, ResidualType* resBase)
{
	Indexer idxr(_disc);

	// Go to the particle block of the given column cell
	StateType const* y = yBase + idxr.offsetCp(colCell);
	double const* yDot = yDotBase + idxr.offsetCp(colCell);
	ResidualType* res = resBase + idxr.offsetCp(colCell);

	const unsigned int requiredMem = (_binding->workspaceSize() + sizeof(double) - 1) / sizeof(double);
	double* const buffer = _tempState + requiredMem * colCell;

	// Prepare parameters
	const ParamType radius = static_cast<ParamType>(_parRadius);

	// Midpoint of current column cell (z coordinate) - needed in externally dependent adsorption kinetic
	const double z = 1.0 / static_cast<double>(_disc.nCol) * (0.5 + colCell);

	// Add time derivatives
	if (yDotBase)
	{
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++res, ++y, ++yDot)
		{
			*res = 0.0;
			const unsigned int nBound = _disc.nBound[comp];
			const ParamType invBetaP = (1.0 - static_cast<ParamType>(_parPorosity)) / (static_cast<ParamType>(_poreAccessFactor[comp]) * static_cast<ParamType>(_parPorosity));

			// Ultimately, we need dc_{p,comp} / dt + 1 / beta_p * [ sum_i  dq_comp^i / dt ]
			// Compute the sum in the brackets first, then divide by beta_p and add dc_p / dt

			// Sum dq_comp^1 / dt + dq_comp^2 / dt + ... + dq_comp^{N_comp} / dt
			for (unsigned int i = 0; i < nBound; ++i)
				// Index explanation:
				//   -comp -> go back to beginning of liquid phase
				//   + strideParLiquid() skip to solid phase
				//   + offsetBoundComp() jump to component (skips all bound states of previous components)
				//   + i go to current bound state
				// Remember this, you'll see it quite a lot ...
				*res += yDot[idxr.strideParLiquid() - comp + idxr.offsetBoundComp(comp) + i];

			// Divide by beta_p and add dcp_i / dt
			*res = timeFactor * (yDot[0] + invBetaP * res[0]);
		}
	}
	else
	{
		std::fill(res, res + _disc.nComp, 0.0);

		// Advance over liquid phase
		res += _disc.nComp;
		y += _disc.nComp;
		yDot += _disc.nComp;
	}

	// Bound phases
	if (!yDotBase)
		yDot = nullptr;

	_binding->residual(static_cast<ParamType>(t), z, static_cast<double>(radius) * 0.5, secIdx, static_cast<ParamType>(timeFactor), y, yDot, res, buffer);
	if (wantJac)
	{
		if (cadet_likely(_disc.strideBound > 0))
		{
			// The RowIterator is always centered on the main diagonal.
			// This means that jac[0] is the main diagonal, jac[-1] is the first lower diagonal,
			// and jac[1] is the first upper diagonal. We can also access the rows from left to
			// right beginning with the last lower diagonal moving towards the main diagonal and
			// continuing to the last upper diagonal by using the native() method.
			linalg::BandMatrix::RowIterator jac = _jacP.row(colCell * idxr.strideParBlock() + idxr.strideParLiquid());

			// static_cast should be sufficient here, but this statement is also analyzed when wantJac = false
			_binding->analyticJacobian(static_cast<double>(t), z, static_cast<double>(radius) * 0.5, secIdx, reinterpret_cast<double const*>(y), jac, buffer);
		}
	}

	return 0;
}

template <typename StateType, typename ResidualType, typename ParamType>
int LumpedRateModelWithPores::residualFlux(const ParamType& t, unsigned int secIdx, StateType const* yBase, double const* yDotBase, ResidualType* resBase)
{
	Indexer idxr(_disc);

	const ParamType invBetaC = 1.0 / static_cast<ParamType>(_colPorosity) - 1.0;
	const ParamType epsP = static_cast<ParamType>(_parPorosity);
	const ParamType radius = static_cast<ParamType>(_parRadius);

	active const* const filmDiff = getSectionDependentSlice(_filmDiffusion, _disc.nComp, secIdx);

	const ParamType jacCF_val = invBetaC * 3.0 / radius;
	const ParamType jacPF_val = -3.0 / (epsP * radius);

	// Get offsets
	ResidualType* const resCol = resBase + idxr.offsetC();
	ResidualType* const resPar = resBase + idxr.offsetCp();
	ResidualType* const resFlux = resBase + idxr.offsetJf();

	StateType const* const yCol = yBase + idxr.offsetC();
	StateType const* const yPar = yBase + idxr.offsetCp();
	StateType const* const yFlux = yBase + idxr.offsetJf();

	// J_f block (identity matrix), adds flux state to flux equation
	for (unsigned int i = 0; i < _disc.nComp * _disc.nCol; ++i)
		resFlux[i] = yFlux[i];

	// J_{0,f} block, adds flux to column void / bulk volume equations
	for (unsigned int i = 0; i < _disc.nCol * _disc.nComp; ++i)
	{
		const unsigned int comp = i % _disc.nComp;
		resCol[i] += jacCF_val * static_cast<ParamType>(filmDiff[comp]) * yFlux[i];
	}

	// J_{f,0} block, adds bulk volume state c_i to flux equation
	for (unsigned int bnd = 0; bnd < _disc.nCol; ++bnd)
	{
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			const unsigned int eq = bnd * idxr.strideColCell() + comp * idxr.strideColComp();
			resFlux[eq] -= yCol[eq];
		}
	}

	// J_{p,f} block, adds flux to particle / bead volume equations
	for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
	{
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			const unsigned int eq = pblk * idxr.strideColCell() + comp * idxr.strideColComp();
			resPar[pblk * idxr.strideParBlock() + comp] += jacPF_val / static_cast<ParamType>(_poreAccessFactor[comp]) * static_cast<ParamType>(filmDiff[comp]) * yFlux[eq];
		}
	}

	// J_{f,p} block, adds outer bead shell state c_{p,i} to flux equation
	for (unsigned int pblk = 0; pblk < _disc.nCol; pblk++)
	{
		for (unsigned int comp = 0; comp < _disc.nComp; comp++)
		{
			const unsigned int eq = pblk * idxr.strideColCell() + comp * idxr.strideColComp();
			resFlux[eq] += yPar[comp + pblk * idxr.strideParBlock()];
		}
	}

	return 0;
}

/**
 * @brief Assembles off diagonal Jacobian blocks
 * @details Assembles the fixed blocks @f$ J_{0,f}, \dots, J_{N_p,f} @f$ and @f$ J_{f,0}, \dots, J_{f, N_p}. @f$
 *          The blocks are fixed for each section.
 * @param [in] t Current time
 * @param [in] secIdx Index of the current section
 */
void LumpedRateModelWithPores::assembleOffdiagJac(double t, unsigned int secIdx)
{
	// Clear matrices for new assembly
	_jacCF.clear();
	_jacFC.clear();
	_jacPF.clear();
	_jacFP.clear();

	Indexer idxr(_disc);

	const double invBetaC = 1.0 / static_cast<double>(_colPorosity) - 1.0;
	const double epsP = static_cast<double>(_parPorosity);
	const double radius = static_cast<double>(_parRadius);
	const double jacCF_val = invBetaC * 3.0 / radius;
	const double jacPF_val = -3.0 / (radius * epsP);

	active const* const filmDiff = getSectionDependentSlice(_filmDiffusion, _disc.nComp, secIdx);

	// Note that the J_f block, which is the identity matrix, is treated in the linear solver

	// J_{0,f} block, adds flux to column void / bulk volume equations
	for (unsigned int eq = 0; eq < _disc.nCol * _disc.nComp; ++eq)
	{
		const unsigned int comp = eq % _disc.nComp;

		// Main diagonal corresponds to j_{f,i} (flux) state variable
		_jacCF.addElement(eq, eq, jacCF_val * static_cast<double>(filmDiff[comp]));
	}

	// J_{f,0} block, adds bulk volume state c_i to flux equation
	for (unsigned int eq = 0; eq < _disc.nCol * _disc.nComp; ++eq)
	{
		_jacFC.addElement(eq, eq, -1.0);
	}

	// J_{p,f} block, implements bead boundary condition in outer bead shell equation
	for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
	{
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			const unsigned int eq = pblk * idxr.strideColCell() + comp * idxr.strideColComp();
			const unsigned int col = pblk * idxr.strideParBlock() + comp;
			_jacPF.addElement(col, eq, jacPF_val / static_cast<double>(_poreAccessFactor[comp]) * static_cast<double>(filmDiff[comp]));
		}
	}

	// J_{f,p} block, adds outer bead shell state c_{p,i} to flux equation
	for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
	{
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			const unsigned int eq = pblk * idxr.strideColCell() + comp * idxr.strideColComp();
			const unsigned int col = pblk * idxr.strideParBlock() + comp;
			_jacFP.addElement(eq, col, 1.0);
		}
	}
}

int LumpedRateModelWithPores::residualSensFwdWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot,
	active* const adRes, active* const adY, unsigned int adDirOffset)
{
	BENCH_SCOPE(_timerResidualSens);

	// Evaluate residual for all parameters using AD in vector mode and at the same time update the 
	// Jacobian (in one AD run, if analytic Jacobians are disabled)
	return residual(t, secIdx, timeFactor, y, yDot, nullptr, adRes, adY, adDirOffset, true, true);
}

int LumpedRateModelWithPores::residualSensFwdAdOnly(const active& t, unsigned int secIdx, const active& timeFactor,
	double const* const y, double const* const yDot, active* const adRes)
{
	BENCH_SCOPE(_timerResidualSens);

	// Evaluate residual for all parameters using AD in vector mode
	return residualImpl<double, active, active, false>(t, secIdx, timeFactor, y, yDot, adRes); 
}

int LumpedRateModelWithPores::residualSensFwdCombine(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, 
	const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS, active const* adRes, 
	double* const tmp1, double* const tmp2, double* const tmp3)
{
	BENCH_SCOPE(_timerResidualSens);

	// tmp1 stores result of (dF / dy) * s
	// tmp2 stores result of (dF / dyDot) * sDot

	for (unsigned int param = 0; param < yS.size(); param++)
	{
		// Directional derivative (dF / dy) * s
		multiplyWithJacobian(0.0, 0u, 1.0, nullptr, nullptr, yS[param], 1.0, 0.0, tmp1);

		// Directional derivative (dF / dyDot) * sDot
		multiplyWithDerivativeJacobian(0.0, 0u, static_cast<double>(timeFactor), nullptr, nullptr, ySdot[param], tmp2);

		double* const ptrResS = resS[param];

		BENCH_START(_timerResidualSensPar);

		// Complete sens residual is the sum:
		// TODO: Chunk TBB loop
#ifdef CADET_PARALLELIZE
		tbb::parallel_for(size_t(0), size_t(numDofs()), [&](size_t i)
#else
		for (unsigned int i = 0; i < numDofs(); ++i)
#endif
		{
			ptrResS[i] = tmp1[i] + tmp2[i] + adRes[i].getADValue(param);
		} CADET_PARFOR_END;

		BENCH_STOP(_timerResidualSensPar);
	}

	return 0;
}

/**
 * @brief Multiplies the given vector with the system Jacobian (i.e., @f$ \frac{\partial F}{\partial y}\left(t, y, \dot{y}\right) @f$)
 * @details Actually, the operation @f$ z = \alpha \frac{\partial F}{\partial y} x + \beta z @f$ is performed.
 * 
 *          Note that residual() or one of its cousins has to be called with the requested point @f$ (t, y, \dot{y}) @f$ once
 *          before calling multiplyWithJacobian() as this implementation ignores the given @f$ (t, y, \dot{y}) @f$.
 * @param [in] t Current time point
 * @param [in] secIdx Index of the current section
 * @param [in] timeFactor Used for time transformation (pre factor of time derivatives) and to compute parameter derivatives with respect to section length
 * @param [in] y Pointer to local state vector
 * @param [in] yDot Pointer to local time derivative state vector
 * @param [in] yS Vector @f$ x @f$ that is transformed by the Jacobian @f$ \frac{\partial F}{\partial y} @f$
 * @param [in] alpha Factor @f$ \alpha @f$ in front of @f$ \frac{\partial F}{\partial y} @f$
 * @param [in] beta Factor @f$ \beta @f$ in front of @f$ z @f$
 * @param [in,out] ret Vector @f$ z @f$ which stores the result of the operation
 */
void LumpedRateModelWithPores::multiplyWithJacobian(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double const* yS, double alpha, double beta, double* ret)
{
	Indexer idxr(_disc);

	// Handle identity matrix of inlet DOFs
	for (unsigned int i = 0; i < _disc.nComp; ++i)
	{
		ret[i] = alpha * yS[i] + beta * ret[i];
	}

	// Interstitial block
	_convDispOp.jacobian().multiplyVector(yS + idxr.offsetC(), alpha, beta, ret + idxr.offsetC());
	_jacCF.multiplyVector(yS + idxr.offsetJf(), alpha, 1.0, ret + idxr.offsetC());

	// Particle block
	const int localOffset = idxr.offsetCp();
	_jacP.multiplyVector(yS + localOffset, alpha, beta, ret + localOffset);
	_jacPF.multiplyVector(yS + idxr.offsetJf(), alpha, 1.0, ret + localOffset);

	// Handle flux equation

	// Set fluxes(ret) = fluxes(yS)
	// This applies the identity matrix in the bottom right corner of the Jaocbian (flux equation)
	for (unsigned int i = idxr.offsetJf(); i < numDofs(); ++i)
		ret[i] = alpha * yS[i] + beta * ret[i];

	double* const retJf = ret + idxr.offsetJf();
	_jacFC.multiplyVector(yS + idxr.offsetC(), alpha, 1.0, retJf);
	_jacFP.multiplyVector(yS + idxr.offsetCp(), alpha, 1.0, retJf);

	// Map inlet DOFs to the column inlet (first bulk cells)
	_jacInlet.multiplyAdd(yS, ret + idxr.offsetC(), alpha);
}

/**
 * @brief Multiplies the time derivative Jacobian @f$ \frac{\partial F}{\partial \dot{y}}\left(t, y, \dot{y}\right) @f$ with a given vector
 * @details The operation @f$ z = \frac{\partial F}{\partial \dot{y}} x @f$ is performed.
 *          The matrix-vector multiplication is transformed matrix-free (i.e., no matrix is explicitly formed).
 * @param [in] t Current time point
 * @param [in] secIdx Index of the current section
 * @param [in] timeFactor Factor which is premultiplied to the time derivatives originating from time transformation
 * @param [in] y Pointer to local state vector
 * @param [in] yDot Pointer to local time derivative state vector
 * @param [in] sDot Vector @f$ x @f$ that is transformed by the Jacobian @f$ \frac{\partial F}{\partial \dot{y}} @f$
 * @param [out] ret Vector @f$ z @f$ which stores the result of the operation
 */
void LumpedRateModelWithPores::multiplyWithDerivativeJacobian(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double const* sDot, double* ret)
{
	Indexer idxr(_disc);

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(size_t(0), size_t(_disc.nCol + 1), [&](size_t idx)
#else
	for (unsigned int idx = 0; idx < _disc.nCol + 1; ++idx)
#endif
	{
		if (cadet_unlikely(idx == 0))
		{
			_convDispOp.multiplyWithDerivativeJacobian(t, secIdx, timeFactor, sDot, ret);
		}
		else
		{
			const unsigned int pblk = idx - 1;

			// Particle
			double const* const localSdot = sDot + idxr.offsetCp(pblk);
			double* const localRet = ret + idxr.offsetCp(pblk);

			// Mobile phase
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				// Add derivative with respect to dc_p / dt to Jacobian
				localRet[comp] = timeFactor * localSdot[comp];

				const double invBetaP = (1.0 - static_cast<double>(_parPorosity)) / (static_cast<double>(_poreAccessFactor[comp]) * static_cast<double>(_parPorosity)) * timeFactor;

				// Add derivative with respect to dq / dt to Jacobian (normal equations)
				for (unsigned int i = 0; i < _disc.nBound[comp]; ++i)
				{
					// Index explanation:
					//   nComp -> skip mobile phase
					//   + _disc.boundOffset[comp] skip bound states of all previous components
					//   + i go to current bound state
					localRet[comp] += invBetaP * localSdot[_disc.nComp + _disc.boundOffset[comp] + i];
				}
			}

			// Solid phase
			_binding->multiplyWithDerivativeJacobian(localSdot + _disc.nComp, localRet + _disc.nComp, timeFactor);
		}
	} CADET_PARFOR_END;

	// Handle fluxes (all algebraic)
	double* const dFdyDot = ret + idxr.offsetJf();
	std::fill(dFdyDot, dFdyDot + _disc.nCol * _disc.nComp, 0.0);

	// Handle inlet DOFs (all algebraic)
	std::fill_n(ret, _disc.nComp, 0.0);
}

void LumpedRateModelWithPores::setExternalFunctions(IExternalFunction** extFuns, unsigned int size)
{
	if (_binding)
		_binding->setExternalFunctions(extFuns, size);
}

unsigned int LumpedRateModelWithPores::localOutletComponentIndex() const CADET_NOEXCEPT
{
	// Inlets are duplicated so need to be accounted for
	if (static_cast<double>(_convDispOp.currentVelocity()) >= 0.0)
		// Forward Flow: outlet is last cell
		return _disc.nComp + (_disc.nCol - 1) * _disc.nComp;
	else
		// Backward flow: Outlet is first cell
		return _disc.nComp;
}

unsigned int LumpedRateModelWithPores::localInletComponentIndex() const CADET_NOEXCEPT
{
	return 0;
}

unsigned int LumpedRateModelWithPores::localOutletComponentStride() const CADET_NOEXCEPT
{
	return 1;
}

unsigned int LumpedRateModelWithPores::localInletComponentStride() const CADET_NOEXCEPT
{
	return 1;
}

void LumpedRateModelWithPores::expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut)
{
	// @todo Write this function
}

void registerLumpedRateModelWithPores(std::unordered_map<std::string, std::function<IUnitOperation*(UnitOpIdx)>>& models)
{
	models[LumpedRateModelWithPores::identifier()] = [](UnitOpIdx uoId) { return new LumpedRateModelWithPores(uoId); };
	models["LRMP"] = [](UnitOpIdx uoId) { return new LumpedRateModelWithPores(uoId); };
}

}  // namespace model

}  // namespace cadet
