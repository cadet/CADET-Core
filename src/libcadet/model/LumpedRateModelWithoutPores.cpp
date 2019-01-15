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

#include "model/LumpedRateModelWithoutPores.hpp"
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
#include "model/parts/BindingConsistentInit.hpp"

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

namespace
{
	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	class ConvOpResidual
	{
	public:
		static inline void call(cadet::model::parts::ConvectionDispersionOperatorBase& op, const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* const y, double const* const yDot, ResidualType* const res, cadet::linalg::BandMatrix& jac)
		{
			// This should not be reached
			cadet_assert(false);
		}
	};

	template <typename ResidualType, typename ParamType>
	class ConvOpResidual<double, ResidualType, ParamType, true> 
	{
	public:
		static inline void call(cadet::model::parts::ConvectionDispersionOperatorBase& op, const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, double const* const y, double const* const yDot, ResidualType* const res, cadet::linalg::BandMatrix& jac)
		{
			op.residual(t, secIdx, timeFactor, y, yDot, res, &jac);
		}
	};

	template <typename ResidualType, typename ParamType>
	class ConvOpResidual<double, ResidualType, ParamType, false> 
	{
	public:
		static inline void call(cadet::model::parts::ConvectionDispersionOperatorBase& op, const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, double const* const y, double const* const yDot, ResidualType* const res, cadet::linalg::BandMatrix& jac)
		{
			op.residual(t, secIdx, timeFactor, y, yDot, res, nullptr);
		}
	};

	template <typename ResidualType, typename ParamType>
	class ConvOpResidual<cadet::active, ResidualType, ParamType, false> 
	{
	public:
		static inline void call(cadet::model::parts::ConvectionDispersionOperatorBase& op, const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, cadet::active const* const y, double const* const yDot, ResidualType* const res, cadet::linalg::BandMatrix& jac)
		{
			op.residual(t, secIdx, timeFactor, y, yDot, res);
		}
	};
}

namespace cadet
{

namespace model
{

LumpedRateModelWithoutPores::LumpedRateModelWithoutPores(UnitOpIdx unitOpIdx) : UnitOperationBase(unitOpIdx),
	_jacInlet(), _analyticJac(true), _jacobianAdDirs(0), _factorizeJacobian(false), _tempState(nullptr), _initC(0),
	_initQ(0), _initState(0), _initStateDot(0)
{
}

LumpedRateModelWithoutPores::~LumpedRateModelWithoutPores() CADET_NOEXCEPT
{
	delete[] _tempState;

	delete[] _disc.nBound;
	delete[] _disc.boundOffset;
}

unsigned int LumpedRateModelWithoutPores::numDofs() const CADET_NOEXCEPT
{
	// Column bulk DOFs: nCol * nComp mobile phase and nCol * (sum boundStates) solid phase
	// Inlet DOFs: nComp
	return _disc.nCol * (_disc.nComp + _disc.strideBound) + _disc.nComp;
}

unsigned int LumpedRateModelWithoutPores::numPureDofs() const CADET_NOEXCEPT
{
	// Column bulk DOFs: nCol * nComp mobile phase and nCol * (sum boundStates) solid phase
	return _disc.nCol * (_disc.nComp + _disc.strideBound);
}


bool LumpedRateModelWithoutPores::usesAD() const CADET_NOEXCEPT
{
#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	// We always need AD if we want to check the analytical Jacobian
	return true;
#else
	// We only need AD if we are not computing the Jacobian analytically
	return !_analyticJac;
#endif
}

bool LumpedRateModelWithoutPores::configureModelDiscretization(IParameterProvider& paramProvider, IConfigHelper& helper)
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

	// Allocate space for initial conditions
	_initC.resize(_disc.nComp);
	_initQ.resize(_disc.strideBound);

	paramProvider.popScope();

	const unsigned int strideCell = _disc.nComp + _disc.strideBound;
	const bool transportSuccess = _convDispOp.configureModelDiscretization(paramProvider, _disc.nComp, _disc.nCol, strideCell);

	// Allocate memory
	Indexer idxr(_disc);

	_jacInlet.resize(_disc.nComp);

	const unsigned int lb = _convDispOp.jacobianLowerBandwidth();
	const unsigned int ub = _convDispOp.jacobianUpperBandwidth();
	const unsigned int mb = _convDispOp.jacobianDiscretizedBandwidth();

	// Allocate matrices such that bandwidths can be switched (backwards flow support)
	_jac.resize(_disc.nCol * strideCell, lb, ub);

	_jacDisc.resize(_disc.nCol * strideCell, mb, mb);
	_jacDisc.repartition(lb, ub);

	// Set whether analytic Jacobian is used
	useAnalyticJacobian(analyticJac);

	// ==== Construct and configure binding model
	clearBindingModels();
	_binding.push_back(nullptr);

	_binding[0] = helper.createBindingModel(paramProvider.getString("ADSORPTION_MODEL"));
	if (!_binding[0])
		throw InvalidParameterException("Unknown binding model " + paramProvider.getString("ADSORPTION_MODEL"));

	const bool bindingConfSuccess = _binding[0]->configureModelDiscretization(paramProvider, _disc.nComp, _disc.nBound, _disc.boundOffset);

	// setup the memory for tempState based on state vector or memory needed for consistent initialization of isotherms, whichever is larger
	unsigned int size = numDofs();
	if (_binding[0]->requiresWorkspace())
	{
		// Required memory (number of doubles) for nonlinear solvers
		const unsigned int requiredMem = (_binding[0]->workspaceSize(_disc.nComp, _disc.strideBound, _disc.nBound) + sizeof(double) - 1) / sizeof(double) * _disc.nCol;
		if (requiredMem > size)
		{
			size = requiredMem;
		}
	}
	_tempState = new double[size];

	return transportSuccess && bindingConfSuccess;
}

bool LumpedRateModelWithoutPores::configure(IParameterProvider& paramProvider)
{
	_parameters.clear();

	const bool transportSuccess = _convDispOp.configure(_unitOpIdx, paramProvider, _parameters);

	// Read geometry parameters
	_totalPorosity = paramProvider.getDouble("TOTAL_POROSITY");

	// Add parameters to map
	_parameters[makeParamId(hashString("TOTAL_POROSITY"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_totalPorosity;

	// Register initial conditions parameters
	for (unsigned int i = 0; i < _disc.nComp; ++i)
		_parameters[makeParamId(hashString("INIT_C"), _unitOpIdx, i, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = _initC.data() + i;

	if (_binding[0])
	{
		std::vector<ParameterId> initParams(_disc.strideBound);
		_binding[0]->fillBoundPhaseInitialParameters(initParams.data(), _unitOpIdx, 0);

		for (unsigned int i = 0; i < _disc.strideBound; ++i)
			_parameters[initParams[i]] = _initQ.data() + i;
	}

	// Reconfigure binding model
	if (_binding[0] && paramProvider.exists("adsorption") && _binding[0]->requiresConfiguration())
	{
		paramProvider.pushScope("adsorption");
		const bool bindingConfSuccess = _binding[0]->configure(paramProvider, _unitOpIdx, 0);
		paramProvider.popScope();

		return transportSuccess && bindingConfSuccess;
	}

	return transportSuccess;
}

void LumpedRateModelWithoutPores::useAnalyticJacobian(const bool analyticJac)
{
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	_analyticJac = analyticJac;
	if (!_analyticJac)
		_jacobianAdDirs = _jac.stride();
	else
		_jacobianAdDirs = 0;
#else
	_analyticJac = false;
	_jacobianAdDirs = _jac.stride();
#endif
}

void LumpedRateModelWithoutPores::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active* const adRes, active* const adY, unsigned int adDirOffset)
{
	Indexer idxr(_disc);

	// ConvectionDispersionOperator tells us whether flow direction has changed
	if (!_convDispOp.notifyDiscontinuousSectionTransition(t, secIdx))
		return;

	// Setup the matrix connecting inlet DOFs to first column cells
	_jacInlet.clear();
	const double h = static_cast<double>(_convDispOp.columnLength()) / static_cast<double>(_disc.nCol);
	const double u = static_cast<double>(_convDispOp.currentVelocity());

	const unsigned int lb = _convDispOp.jacobianLowerBandwidth();
	const unsigned int ub = _convDispOp.jacobianUpperBandwidth();
	if (u >= 0.0)
	{
		// Forwards flow

		// Place entries for inlet DOF to first column cell conversion
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			_jacInlet.addElement(comp * idxr.strideColComp(), comp, -u / h);

		// Repartition Jacobians
		_jac.repartition(lb, ub);
		_jacDisc.repartition(lb, ub);
	}
	else
	{
		// Backwards flow

		// Place entries for inlet DOF to last column cell conversion
		const unsigned int offset = (_disc.nCol - 1) * idxr.strideColCell();
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			_jacInlet.addElement(offset + comp * idxr.strideColComp(), comp, u / h);

		// Repartition Jacobians
		_jac.repartition(ub, lb);
		_jacDisc.repartition(ub, lb);
	}

	prepareADvectors(adRes, adY, adDirOffset);	
}

void LumpedRateModelWithoutPores::setFlowRates(const active& in, const active& out) CADET_NOEXCEPT
{
	_convDispOp.setFlowRates(in, out, _totalPorosity);
}

void LumpedRateModelWithoutPores::reportSolution(ISolutionRecorder& recorder, double const* const solution) const
{
	Exporter expr(_disc, solution);
	recorder.beginUnitOperation(_unitOpIdx, *this, expr);
	recorder.endUnitOperation();
}

void LumpedRateModelWithoutPores::reportSolutionStructure(ISolutionRecorder& recorder) const
{
	Exporter expr(_disc, nullptr);
	recorder.unitOperationStructure(_unitOpIdx, *this, expr);
}


unsigned int LumpedRateModelWithoutPores::requiredADdirs() const CADET_NOEXCEPT
{
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	return _jacobianAdDirs;
#else
	// If CADET_CHECK_ANALYTIC_JACOBIAN is active, we always need the AD directions for the Jacobian
	return _jac.stride();
#endif
}

void LumpedRateModelWithoutPores::prepareADvectors(active* const adRes, active* const adY, unsigned int adDirOffset) const
{
	// Early out if AD is disabled
	if (!adY)
		return;

	Indexer idxr(_disc);

	// Get bandwidths
	const unsigned int lowerBandwidth = _jac.lowerBandwidth();
	const unsigned int upperBandwidth = _jac.upperBandwidth();

	ad::prepareAdVectorSeedsForBandMatrix(adY + _disc.nComp, adDirOffset, _jac.rows(), lowerBandwidth, upperBandwidth, lowerBandwidth);
}

/**
 * @brief Extracts the system Jacobian from band compressed AD seed vectors
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
void LumpedRateModelWithoutPores::extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset)
{
	Indexer idxr(_disc);
	ad::extractBandedJacobianFromAd(adRes + idxr.offsetC(), adDirOffset, _jac.lowerBandwidth(), _jac);
}

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN

/**
 * @brief Compares the analytical Jacobian with a Jacobian derived by AD
 * @details The analytical Jacobian is assumed to be stored in the corresponding band matrices.
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
void LumpedRateModelWithoutPores::checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const
{
	Indexer idxr(_disc);

	const double maxDiff = ad::compareBandedJacobianWithAd(adRes + idxr.offsetC(), adDirOffset, _jac.lowerBandwidth(), _jac);
	LOG(Debug) << "AD dir offset: " << adDirOffset << " DiagDirCol: " << _jac.lowerBandwidth() << " MaxDiff: " << maxDiff;
}

#endif

int LumpedRateModelWithoutPores::residual(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double* const res)
{
	BENCH_SCOPE(_timerResidual);

	// Evaluate residual do not compute Jacobian or parameter sensitivities
	return residualImpl<double, double, double, false>(t, secIdx, timeFactor, y, yDot, res);
}

int LumpedRateModelWithoutPores::residualWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, double* const res, active* const adRes, active* const adY, unsigned int adDirOffset)
{
	BENCH_SCOPE(_timerResidual);

	// Evaluate residual, use AD for Jacobian if required but do not evaluate parameter derivatives
	return residual(t, secIdx, timeFactor, y, yDot, res, adRes, adY, adDirOffset, true, false);
}

int LumpedRateModelWithoutPores::residual(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, double* const res, 
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
int LumpedRateModelWithoutPores::residualImpl(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* const y, double const* const yDot, ResidualType* const res)
{
	ConvOpResidual<StateType, ResidualType, ParamType, wantJac>::call(_convDispOp, t, secIdx, timeFactor, y, yDot, res, _jac);

	Indexer idxr(_disc);
	const unsigned int requiredMem = (_binding[0]->workspaceSize(_disc.nComp, _disc.strideBound, _disc.nBound) + sizeof(double) - 1) / sizeof(double);

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(size_t(0), size_t(_disc.nCol), [&](size_t col)
#else
	for (unsigned int col = 0; col < _disc.nCol; ++col)
#endif
	{
		StateType const* const localY = y + idxr.offsetC() + idxr.strideColCell() * col;
		ResidualType* const localRes = res + idxr.offsetC() + idxr.strideColCell() * col;
		double* const buffer = _tempState + requiredMem * col;

		if (yDot)
		{
			const double invBeta = 1.0 / static_cast<double>(_totalPorosity) - 1.0;
			double const* const localYdot = yDot + idxr.offsetC() + idxr.strideColCell() * col;

			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				// Sum dq_comp^1 / dt + dq_comp^2 / dt + ... + dq_comp^{N_comp} / dt
				double sumQdot = 0.0;
				for (unsigned int i = 0; i < _disc.nBound[comp]; ++i)
					// Index explanation:
					//   + strideColLiquid() skip to solid phase
					//   + offsetBoundComp() jump to component (skips all bound states of previous components)
					//   + i go to current bound state
					sumQdot += localYdot[idxr.strideColLiquid() + idxr.offsetBoundComp(comp) + i];

				// Divide by beta and add to dc_i / dt
				localRes[comp] += timeFactor * invBeta * sumQdot;
			}
		}

		// Midpoint of current column cell (z coordinate) - needed in externally dependent adsorption kinetic
		const double z = 1.0 / static_cast<double>(_disc.nCol) * (0.5 + col);

		double const* const localQdot = yDot ? yDot + idxr.offsetC() + idxr.strideColCell() * col + idxr.strideColLiquid() : nullptr;
		_binding[0]->residual(t, z, 0.0, secIdx, timeFactor, localY + idxr.strideColLiquid(), localY, localQdot, localRes + idxr.strideColLiquid(), buffer);
		if (wantJac)
		{
			if (cadet_likely(_disc.strideBound > 0))
			{
				// The RowIterator is always centered on the main diagonal.
				// This means that jac[0] is the main diagonal, jac[-1] is the first lower diagonal,
				// and jac[1] is the first upper diagonal. We can also access the rows from left to
				// right beginning with the last lower diagonal moving towards the main diagonal and
				// continuing to the last upper diagonal by using the native() method.
				linalg::BandMatrix::RowIterator jac = _jac.row(col * idxr.strideColCell() + idxr.strideColLiquid());

				// static_cast should be sufficient here, but this statement is also analyzed when wantJac = false
				_binding[0]->analyticJacobian(static_cast<double>(t), z, 0.0, secIdx, reinterpret_cast<double const*>(localY) + idxr.strideColLiquid(), idxr.strideColLiquid(), jac, buffer);
			}
		}

	} CADET_PARFOR_END;

	BENCH_STOP(_timerResidualPar);

	// Handle inlet DOFs, which are simply copied to res
	for (unsigned int i = 0; i < _disc.nComp; ++i)
	{
		res[i] = y[i];
	}

	return 0;
}

int LumpedRateModelWithoutPores::residualSensFwdWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot,
	active* const adRes, active* const adY, unsigned int adDirOffset)
{
	BENCH_SCOPE(_timerResidualSens);

	// Evaluate residual for all parameters using AD in vector mode and at the same time update the 
	// Jacobian (in one AD run, if analytic Jacobians are disabled)
	return residual(t, secIdx, timeFactor, y, yDot, nullptr, adRes, adY, adDirOffset, true, true);
}

int LumpedRateModelWithoutPores::residualSensFwdAdOnly(const active& t, unsigned int secIdx, const active& timeFactor,
	double const* const y, double const* const yDot, active* const adRes)
{
	BENCH_SCOPE(_timerResidualSens);

	// Evaluate residual for all parameters using AD in vector mode
	return residualImpl<double, active, active, false>(t, secIdx, timeFactor, y, yDot, adRes); 
}

int LumpedRateModelWithoutPores::residualSensFwdCombine(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, 
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
void LumpedRateModelWithoutPores::multiplyWithJacobian(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double const* yS, double alpha, double beta, double* ret)
{
	Indexer idxr(_disc);

	// Handle identity matrix of inlet DOFs
	for (unsigned int i = 0; i < _disc.nComp; ++i)
	{
		ret[i] = alpha * yS[i] + beta * ret[i];
	}

	// Main Jacobian
	_jac.multiplyVector(yS + idxr.offsetC(), alpha, beta, ret + idxr.offsetC());

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
void LumpedRateModelWithoutPores::multiplyWithDerivativeJacobian(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double const* sDot, double* ret)
{
	Indexer idxr(_disc);
	const double invBeta = (1.0 / static_cast<double>(_totalPorosity) - 1.0) * timeFactor;

	_convDispOp.multiplyWithDerivativeJacobian(t, secIdx, timeFactor, sDot, ret);
	for (unsigned int col = 0; col < _disc.nCol; ++col)
	{
		const unsigned int localOffset = idxr.offsetC() + col * idxr.strideColCell();
		double const* const localSdot = sDot + localOffset;
		double* const localRet = ret + localOffset;

		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			// Add derivative with respect to dq / dt to Jacobian (normal equations)
			for (unsigned int i = 0; i < _disc.nBound[comp]; ++i)
			{
				// Index explanation:
				//   nComp -> skip mobile phase
				//   + _disc.boundOffset[comp] skip bound states of all previous components
				//   + i go to current bound state
				localRet[comp] += invBeta * localSdot[_disc.nComp + _disc.boundOffset[comp] + i];
			}
		}

		// Solid phase
		_binding[0]->multiplyWithDerivativeJacobian(localSdot + _disc.nComp, localRet + _disc.nComp, timeFactor);
	}

	// Handle inlet DOFs (all algebraic)
	std::fill_n(ret, _disc.nComp, 0.0);
}

void LumpedRateModelWithoutPores::setExternalFunctions(IExternalFunction** extFuns, unsigned int size)
{
	if (_binding[0])
		_binding[0]->setExternalFunctions(extFuns, size);
}

unsigned int LumpedRateModelWithoutPores::localOutletComponentIndex() const CADET_NOEXCEPT
{
	// Inlets are duplicated so need to be accounted for
	if (static_cast<double>(_convDispOp.currentVelocity()) >= 0.0)
		// Forward Flow: outlet is last cell
		return _disc.nComp + (_disc.nCol - 1) * (_disc.nComp + _disc.strideBound);
	else
		// Backward flow: Outlet is first cell
		return _disc.nComp;
}

unsigned int LumpedRateModelWithoutPores::localInletComponentIndex() const CADET_NOEXCEPT
{
	return 0;
}

unsigned int LumpedRateModelWithoutPores::localOutletComponentStride() const CADET_NOEXCEPT
{
	return 1;
}

unsigned int LumpedRateModelWithoutPores::localInletComponentStride() const CADET_NOEXCEPT
{
	return 1;
}

void LumpedRateModelWithoutPores::expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut)
{
	// @todo Write this function
}

/**
 * @brief Computes the solution of the linear system involving the system Jacobian
 * @details The system \f[ \left( \frac{\partial F}{\partial y} + \alpha \frac{\partial F}{\partial \dot{y}} \right) x = b \f]
 *          has to be solved. The right hand side \f$ b \f$ is given by @p rhs, the Jacobians are evaluated at the
 *          point \f$(y, \dot{y})\f$ given by @p y and @p yDot. The residual @p res at this point, \f$ F(t, y, \dot{y}) \f$,
 *          may help with this. Error weights (see IDAS guide) are given in @p weight. The solution is returned in @p rhs.
 *
 * @param [in] t Current time point
 * @param [in] timeFactor Factor which is premultiplied to the time derivatives originating from time transformation
 * @param [in] alpha Value of \f$ \alpha \f$ (arises from BDF time discretization)
 * @param [in] outerTol Error tolerance for the solution of the linear system from outer Newton iteration
 * @param [in,out] rhs On entry the right hand side of the linear equation system, on exit the solution
 * @param [in] weight Vector with error weights
 * @param [in] y Pointer to state vector at which the Jacobian is evaluated
 * @param [in] yDot Pointer to time derivative state vector at which the Jacobian is evaluated
 * @return @c 0 on success, @c -1 on non-recoverable error, and @c +1 on recoverable error
 */
int LumpedRateModelWithoutPores::linearSolve(double t, double timeFactor, double alpha, double outerTol, double* const rhs, double const* const weight,
	double const* const y, double const* const yDot)
{
	BENCH_SCOPE(_timerLinearSolve);

	Indexer idxr(_disc);

	bool success = true;

	// Factorize Jacobian only if required
	if (_factorizeJacobian)
	{
		// Assemble
		assembleDiscretizedJacobian(alpha, idxr, timeFactor);

		// Factorize
		success = _jacDisc.factorize();
		if (cadet_unlikely(!success))
		{
			LOG(Error) << "Factorize() failed for par block";
		}

		// Do not factorize again at next call without changed Jacobians
		_factorizeJacobian = false;
	}

	// Handle inlet DOFs
	_jacInlet.multiplySubtract(rhs, rhs + idxr.offsetC());

	// Solve
	const bool result = _jacDisc.solve(rhs + idxr.offsetC());
	if (cadet_unlikely(!result))
	{
		LOG(Error) << "Solve() failed for bulk block";
	}

	return (success && result) ? 0 : 1;;
}

/**
 * @brief Assembles the Jacobian of the time-discretized equations
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
 * @param [in] idxr Indexer
 * @param [in] timeFactor Factor which is premultiplied to the time derivatives originating from time transformation
 */
void LumpedRateModelWithoutPores::assembleDiscretizedJacobian(double alpha, const Indexer& idxr, double timeFactor)
{
	// Copy normal matrix over to factorizable matrix
	_jacDisc.copyOver(_jac);

	// Handle transport equations (dc_i / dt terms)
	_convDispOp.addTimeDerivativeToJacobian(alpha, timeFactor, _jacDisc);

	// Add time derivatives to cells
	const double invBeta = 1.0 / static_cast<double>(_totalPorosity) - 1.0;
	linalg::FactorizableBandMatrix::RowIterator jac = _jacDisc.row(0);
	for (unsigned int j = 0; j < _disc.nCol; ++j)
	{
		// Mobile phase (advances jac accordingly)
		addMobilePhaseTimeDerivativeToJacobianCell(jac, idxr, alpha, invBeta, timeFactor);

		// Stationary phase
		_binding[0]->jacobianAddDiscretized(alpha * timeFactor, jac);

		// Advance pointers over all bound states
		jac += _disc.strideBound;
	}
}

/**
 * @brief Adds Jacobian @f$ \frac{\partial F}{\partial \dot{y}} @f$ to cell of system Jacobian
 * @details Actually adds @f$ \alpha \frac{\partial F}{\partial \dot{y}} @f$, which is useful
 *          for constructing the linear system in BDF time discretization.
 * @param [in,out] jac On entry, RowIterator pointing to the beginning of a cell;
 *                     on exit, the iterator points to the end of the mobile phase
 * @param [in] idxr Indexer
 * @param [in] alpha Value of \f$ \alpha \f$ (arises from BDF time discretization)
 * @param [in] invBeta Inverse porosity term @f$\frac{1}{\beta}@f$
 * @param [in] timeFactor Factor which is premultiplied to the time derivatives originating from time transformation
 */
void LumpedRateModelWithoutPores::addMobilePhaseTimeDerivativeToJacobianCell(linalg::FactorizableBandMatrix::RowIterator& jac, const Indexer& idxr, double alpha, double invBeta, double timeFactor)
{
	// Compute total factor
	alpha *= timeFactor;

	// Mobile phase
	for (int comp = 0; comp < static_cast<int>(_disc.nComp); ++comp, ++jac)
	{
		// Add derivative with respect to dq / dt to Jacobian
		for (int i = 0; i < static_cast<int>(_disc.nBound[comp]); ++i)
		{
			// Index explanation:
			//   -comp -> go back to beginning of liquid phase
			//   + strideColLiquid() skip to solid phase
			//   + offsetBoundComp() jump to component (skips all bound states of previous components)
			//   + i go to current bound state
			jac[idxr.strideColLiquid() - comp + idxr.offsetBoundComp(comp) + i] += alpha * invBeta;
		}
	}
}

void LumpedRateModelWithoutPores::applyInitialCondition(double* const vecStateY, double* const vecStateYdot) const
{
	Indexer idxr(_disc);

	// Check whether full state vector is available as initial condition
	if (!_initState.empty())
	{
		std::fill(vecStateY, vecStateY + idxr.offsetC(), 0.0);
		std::copy(_initState.data(), _initState.data() + numPureDofs(), vecStateY + idxr.offsetC());

		if (!_initStateDot.empty())
		{
			std::fill(vecStateYdot, vecStateYdot + idxr.offsetC(), 0.0);
			std::copy(_initStateDot.data(), _initStateDot.data() + numPureDofs(), vecStateYdot + idxr.offsetC());
		}
		else
			std::fill(vecStateYdot, vecStateYdot + numDofs(), 0.0);

		return;
	}

	double* const stateYbulk = vecStateY + idxr.offsetC();

	// Loop over column cells
	for (unsigned int col = 0; col < _disc.nCol; ++col)
	{
		const unsigned int localOffset = col * idxr.strideColCell();

		// Loop over components in cell
		for (unsigned comp = 0; comp < _disc.nComp; ++comp)
			stateYbulk[localOffset + comp * idxr.strideColComp()] = static_cast<double>(_initC[comp]);

		// Initialize q
		for (unsigned int bnd = 0; bnd < _disc.strideBound; ++bnd)
			stateYbulk[localOffset + idxr.strideColLiquid() + bnd] = static_cast<double>(_initQ[bnd]);
	}
}

void LumpedRateModelWithoutPores::readInitialCondition(IParameterProvider& paramProvider)
{
	_initState.clear();
	_initStateDot.clear();

	// Check if INIT_STATE is present
	if (paramProvider.exists("INIT_STATE"))
	{
		const std::vector<double> initState = paramProvider.getDoubleArray("INIT_STATE");
		_initState = std::vector<double>(initState.begin(), initState.begin() + numPureDofs());

		// Check if INIT_STATE contains the full state and its time derivative
		if (initState.size() >= 2 * numPureDofs())
			_initStateDot = std::vector<double>(initState.begin() + numPureDofs(), initState.begin() + 2 * numPureDofs());
		return;
	}

	const std::vector<double> initC = paramProvider.getDoubleArray("INIT_C");
	std::vector<double> initQ;

	if (paramProvider.exists("INIT_Q"))
		initQ = paramProvider.getDoubleArray("INIT_Q");

	if (initC.size() < _disc.nComp)
		throw InvalidParameterException("INIT_C does not contain enough values for all components");

	if ((_disc.strideBound > 0) && (initQ.size() < _disc.strideBound))
		throw InvalidParameterException("INIT_Q does not contain enough values for all bound states");

	ad::copyToAd(initC.data(), _initC.data(), _disc.nComp);
	if (!initQ.empty())
		ad::copyToAd(initQ.data(), _initQ.data(), _disc.strideBound);
}

/**
 * @brief Computes consistent initial values (state variables without their time derivatives)
 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] the initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ have
 *          to be consistent. This functions updates the initial state \f$ y_0 \f$ and overwrites the time
 *          derivative \f$ \dot{y}_0 \f$ such that they are consistent.
 *          
 *          The process works in two steps:
 *          <ol>
 *              <li>Solve all algebraic equations in the model (e.g., quasi-stationary isotherms, reaction equilibria).</li>
 *              <li>Compute the time derivatives of the state @f$ \dot{y} @f$ such that the residual is 0.
 *                 However, because of the algebraic equations, we need additional conditions to fully determine
 *                 @f$ \dot{y}@f$. By differentiating the algebraic equations with respect to time, we get the
 *                 missing linear equations (recall that the state vector @f$ y @f$ is fixed).
 * 
 *     The right hand side of the linear system is given by the negative residual without contribution 
 *     of @f$ \dot{y} @f$ for differential equations and 0 for algebraic equations 
 *     (@f$ -\frac{\partial F}{\partial t}@f$, to be more precise).</li>
 *          </ol>
 *     
 *     This function performs step 1. See consistentInitialTimeDerivative() for step 2.
 *     
 * 	   This function is to be used with consistentInitialTimeDerivative(). Do not mix normal and lean
 *     consistent initialization!
 *     
 * @param [in] t Current time point
 * @param [in] secIdx Index of the current section
 * @param [in] timeFactor Used for time transformation (pre factor of time derivatives) and to compute parameter derivatives with respect to section length
 * @param [in,out] vecStateY State vector with initial values that are to be updated for consistency
 * @param [in,out] adRes Pointer to residual vector of AD datatypes that can be used for computing the Jacobian (or @c nullptr if AD is disabled)
 * @param [in,out] adY Pointer to state vector of AD datatypes that can be used for computing the Jacobian (or @c nullptr if AD is disabled)
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 * @param [in] errorTol Error tolerance for algebraic equations
 * @todo Decrease amount of allocated memory by partially using temporary vectors (state and Schur complement)
 */
void LumpedRateModelWithoutPores::consistentInitialState(double t, unsigned int secIdx, double timeFactor, double* const vecStateY, 
	active* const adRes, active* const adY, unsigned int adDirOffset, double errorTol)
{
	BENCH_SCOPE(_timerConsistentInit);

	// TODO: Check memory consumption and offsets
	// Round up
	const unsigned int requiredMem = (_binding[0]->workspaceSize(_disc.nComp, _disc.strideBound, _disc.nBound) + sizeof(double) - 1) / sizeof(double);

	Indexer idxr(_disc);

	// Step 1: Solve algebraic equations
	if (_binding[0]->hasAlgebraicEquations())
	{
		ad::BandedJacobianExtractor jacExtractor(_jac.lowerBandwidth(), _jac.lowerBandwidth(), _jac.upperBandwidth());

		//Problem capturing variables here
#ifdef CADET_PARALLELIZE
		BENCH_SCOPE(_timerConsistentInitPar);
		tbb::parallel_for(size_t(0), size_t(_disc.nCol), [&](size_t col)
#else
		for (unsigned int col = 0; col < _disc.nCol; ++col)
#endif
		{
			// Reuse memory of band matrix for dense matrix
			linalg::DenseMatrixView jacobianMatrix(_jacDisc.data() + col * _disc.strideBound * _disc.strideBound, _jacDisc.pivot() + col * _disc.strideBound, _disc.strideBound, _disc.strideBound);

			// Midpoint of current column cell (z coordinate) - needed in externally dependent adsorption kinetic
			const double z = 1.0 / static_cast<double>(_disc.nCol) * (0.5 + col);

			const int localOffsetToCell = idxr.offsetC() + col * idxr.strideColCell();
			const int localOffsetInCell = idxr.strideColLiquid();

			// Get pointer to q variables in cell
			double* const qShell = vecStateY + localOffsetToCell + localOffsetInCell;
			active* const localAdRes = adRes ? adRes + localOffsetToCell : nullptr;
			active* const localAdY = adY ? adY + localOffsetToCell : nullptr;

			// We are essentially creating a 2d vector of blocks out of a linear strip of memory
			const unsigned int offset = requiredMem * col;

			// Solve algebraic variables
			_binding[0]->consistentInitialState(t, z, 0.0, secIdx, qShell, qShell - localOffsetInCell, errorTol, localAdRes, localAdY,
				localOffsetInCell, adDirOffset, jacExtractor, _tempState + offset, jacobianMatrix);
		} CADET_PARFOR_END;
	}
}

/**
 * @brief Computes consistent initial time derivatives
 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] the initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ have
 *          to be consistent. This functions updates the initial state \f$ y_0 \f$ and overwrites the time
 *          derivative \f$ \dot{y}_0 \f$ such that they are consistent.
 *          
 *          The process works in two steps:
 *          <ol>
 *              <li>Solve all algebraic equations in the model (e.g., quasi-stationary isotherms, reaction equilibria).</li>
 *              <li>Compute the time derivatives of the state @f$ \dot{y} @f$ such that the residual is 0.
 *                 However, because of the algebraic equations, we need additional conditions to fully determine
 *                 @f$ \dot{y}@f$. By differentiating the algebraic equations with respect to time, we get the
 *                 missing linear equations (recall that the state vector @f$ y @f$ is fixed). 
 *
 *     The right hand side of the linear system is given by the negative residual without contribution 
 *     of @f$ \dot{y} @f$ for differential equations and 0 for algebraic equations 
 *     (@f$ -\frac{\partial F}{\partial t}@f$, to be more precise).</li>
 *          </ol>
 *
 *     This function performs step 2. See consistentInitialState() for step 1.
 *     
 * 	   This function is to be used with consistentInitialState(). Do not mix normal and lean
 *     consistent initialization!
 *     
 * @param [in] t Current time point
 * @param [in] secIdx Index of the current section
 * @param [in] timeFactor Used for time transformation (pre factor of time derivatives) and to compute parameter derivatives with respect to section length
 * @param [in] vecStateY Consistently initialized state vector
 * @param [in,out] vecStateYdot On entry, residual without taking time derivatives into account. On exit, consistent state time derivatives.
 */
void LumpedRateModelWithoutPores::consistentInitialTimeDerivative(double t, unsigned int secIdx, double timeFactor, double const* vecStateY, double* const vecStateYdot)
{
	BENCH_SCOPE(_timerConsistentInit);

	Indexer idxr(_disc);

	// Step 2: Compute the correct time derivative of the state vector

	// Note that the residual has not been negated, yet. We will do that now.
	for (unsigned int i = 0; i < numDofs(); ++i)
		vecStateYdot[i] = -vecStateYdot[i];

	_jacDisc.setAll(0.0);

	// Handle transport equations (dc_i / dt terms)
	_convDispOp.addTimeDerivativeToJacobian(1.0, timeFactor, _jacDisc);

	const double invBeta = 1.0 / static_cast<double>(_totalPorosity) - 1.0;
	for (unsigned int col = 0; col < _disc.nCol; ++col)
	{
		// Midpoint of current column cell (z coordinate) - needed in externally dependent adsorption kinetic
		const double z = 1.0 / static_cast<double>(_disc.nCol) * (0.5 + col);

		// Assemble
		linalg::FactorizableBandMatrix::RowIterator jac = _jacDisc.row(idxr.strideColCell() * col);

		// Mobile phase (advances jac accordingly)
		addMobilePhaseTimeDerivativeToJacobianCell(jac, idxr, 1.0, invBeta, timeFactor);

		// Stationary phase
		// Populate matrix with time derivative Jacobian first
		_binding[0]->jacobianAddDiscretized(timeFactor, jac);

		// Overwrite rows corresponding to algebraic equations with the Jacobian and set right hand side to 0
		if (_binding[0]->hasAlgebraicEquations())
		{
			parts::BindingConsistentInitializer::consistentInitialTimeDerivative(_binding[0], timeFactor, jac,
				_jac.row(col * idxr.strideColCell() + idxr.strideColLiquid()),
				vecStateYdot + idxr.offsetC() + col * idxr.strideColCell() + idxr.strideColLiquid(),
				t, z, 0.5, secIdx, _tempState);
		}
	}

	// Precondition
	double* const scaleFactors = _tempState + idxr.offsetC();
	_jacDisc.rowScaleFactors(scaleFactors);
	_jacDisc.scaleRows(scaleFactors);

	// Factorize
	const bool result = _jacDisc.factorize();
	if (!result)
	{
		LOG(Error) << "Factorize() failed for par block";
	}

	const bool result2 = _jacDisc.solve(scaleFactors, vecStateYdot + idxr.offsetC());
	if (!result2)
	{
		LOG(Error) << "Solve() failed for par block";
	}
}

/**
 * @brief Computes approximately / partially consistent initial values (state variables without their time derivatives)
 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] the initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ have
 *          to be consistent. This functions updates the initial state \f$ y_0 \f$ and overwrites the time
 *          derivative \f$ \dot{y}_0 \f$ such that they are consistent.
 *          
 *          This function performs a relaxed consistent initialization: Only parts of the vectors are updated
 *          and, hence, consistency is not guaranteed. Since there is less work to do, it is probably faster than
 *          the standard process represented by consistentInitialState().
 *          
 *          The process works in two steps:
 *          <ol>
 *              <li>Keep state and time derivative vectors as they are (i.e., do not solve algebraic equations).</li>
 *              <li>Compute the time derivatives of the state @f$ \dot{y} @f$ such that the residual is 0 for the
 *              mobile phase variables.</li>
 *          </ol>
 *     This function performs step 1. See leanConsistentInitialTimeDerivative() for step 2.
 *     
 * 	   This function is to be used with leanConsistentInitialTimeDerivative(). Do not mix normal and lean
 *     consistent initialization!
 *     
 * @param [in] t Current time point
 * @param [in] secIdx Index of the current section
 * @param [in] timeFactor Used for time transformation (pre factor of time derivatives) and to compute parameter derivatives with respect to section length
 * @param [in,out] vecStateY State vector with initial values that are to be updated for consistency
 * @param [in,out] adRes Pointer to residual vector of AD datatypes that can be used for computing the Jacobian (or @c nullptr if AD is disabled)
 * @param [in,out] adY Pointer to state vector of AD datatypes that can be used for computing the Jacobian (or @c nullptr if AD is disabled)
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 * @param [in] errorTol Error tolerance for algebraic equations
 */
void LumpedRateModelWithoutPores::leanConsistentInitialState(double t, unsigned int secIdx, double timeFactor, double* const vecStateY, 
	active* const adRes, active* const adY, unsigned int adDirOffset, double errorTol)
{
}

/**
 * @brief Computes approximately / partially consistent initial time derivatives
 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] the initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ have
 *          to be consistent. This functions updates the initial state \f$ y_0 \f$ and overwrites the time
 *          derivative \f$ \dot{y}_0 \f$ such that they are consistent.
 *          
 *          This function performs a relaxed consistent initialization: Only parts of the vectors are updated
 *          and, hence, consistency is not guaranteed. Since there is less work to do, it is probably faster than
 *          the standard process represented by consistentInitialState().
 *          
 *          The process works in two steps:
 *          <ol>
 *              <li>Keep state and time derivative vectors as they are (i.e., do not solve algebraic equations).</li>
 *              <li>Compute the time derivatives of the state @f$ \dot{y} @f$ such that the residual is 0 for the
 *              mobile phase variables.</li>
 *          </ol>
 *     This function performs step 2. See leanConsistentInitialState() for step 1.
 *     
 * 	   This function is to be used with leanConsistentInitialState(). Do not mix normal and lean
 *     consistent initialization!
 *     
 * @param [in] t Current time point
 * @param [in] timeFactor Used for time transformation (pre factor of time derivatives) and to compute parameter derivatives with respect to section length
 * @param [in] vecStateY (Lean) consistently initialized state vector
 * @param [in,out] vecStateYdot On entry, inconsistent state time derivatives. On exit, partially consistent state time derivatives.
 * @param [in] res On entry, residual without taking time derivatives into account. The data is overwritten during execution of the function.
 */
void LumpedRateModelWithoutPores::leanConsistentInitialTimeDerivative(double t, double timeFactor, double const* const vecStateY, double* const vecStateYdot, double* const res)
{
	BENCH_SCOPE(_timerConsistentInit);

	Indexer idxr(_disc);

	// Step 2: Compute the correct time derivative of the state vector (only mobile phase DOFs)

	const double invBeta = (1.0 / static_cast<double>(_totalPorosity) - 1.0) * timeFactor;
	for (unsigned int col = 0; col < _disc.nCol; ++col)
	{
		// Offset to current cell's c and q variables
		const unsigned int localOffset = idxr.offsetC() + col * idxr.strideColCell();
		const unsigned int localOffsetQ = localOffset + idxr.strideColLiquid();

		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			// dq_{i,j} / dt is assumed to be fixed, so bring it on the right hand side
			for (unsigned int i = 0; i < _disc.nBound[comp]; ++i)
			{
				res[localOffset + comp] += invBeta * vecStateYdot[localOffsetQ + _disc.boundOffset[comp] + i];
			}

			vecStateYdot[localOffset + comp] = -res[localOffset + comp] / timeFactor;
		}
	}
}

void LumpedRateModelWithoutPores::initializeSensitivityStates(const std::vector<double*>& vecSensY) const
{
	Indexer idxr(_disc);
	for (unsigned int param = 0; param < vecSensY.size(); ++param)
	{
		double* const stateYbulk = vecSensY[param] + idxr.offsetC();

		// Loop over column cells
		for (unsigned int col = 0; col < _disc.nCol; ++col)
		{
			const unsigned int localOffset = col * idxr.strideColCell();

			// Loop over components in cell
			for (unsigned comp = 0; comp < _disc.nComp; ++comp)
				stateYbulk[localOffset + comp * idxr.strideColComp()] = _initC[comp].getADValue(param);

			// Initialize q
			for (unsigned int bnd = 0; bnd < _disc.strideBound; ++bnd)
				stateYbulk[localOffset + idxr.strideColLiquid() + bnd] = _initQ[bnd].getADValue(param);
		}
	}
}

/**
 * @brief Computes consistent initial values and time derivatives of sensitivity subsystems
 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] and initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$,
 *          the sensitivity system for a parameter @f$ p @f$ reads
 *          \f[ \frac{\partial F}{\partial y}(t, y, \dot{y}) s + \frac{\partial F}{\partial \dot{y}}(t, y, \dot{y}) \dot{s} + \frac{\partial F}{\partial p}(t, y, \dot{y}) = 0. \f]
 *          The initial values of this linear DAE, @f$ s_0 = \frac{\partial y_0}{\partial p} @f$ and @f$ \dot{s}_0 = \frac{\partial \dot{y}_0}{\partial p} @f$
 *          have to be consistent with the sensitivity DAE. This functions updates the initial sensitivity\f$ s_0 \f$ and overwrites the time
 *          derivative \f$ \dot{s}_0 \f$ such that they are consistent.
 *          
 *          The process follows closely the one of consistentInitialConditions() and, in fact, is a linearized version of it.
 *          This is necessary because the initial conditions of the sensitivity system \f$ s_0 \f$ and \f$ \dot{s}_0 \f$ are
 *          related to the initial conditions \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ of the original DAE by differentiating them
 *          with respect to @f$ p @f$: @f$ s_0 = \frac{\partial y_0}{\partial p} @f$ and @f$ \dot{s}_0 = \frac{\partial \dot{y}_0}{\partial p}. @f$
 *          <ol>
 *              <li>Solve all algebraic equations in the model (e.g., quasi-stationary isotherms, reaction equilibria).
 *                 Let @f$ \mathcal{I}_a @f$ be the index set of algebraic equations, then, at this point, we have
 *                 \f[ \left( \frac{\partial F}{\partial y}(t, y_0, \dot{y}_0) s + \frac{\partial F}{\partial p}(t, y_0, \dot{y}_0) \right)_{\mathcal{I}_a} = 0. \f]</li>
 *              <li>Compute the time derivatives of the sensitivity @f$ \dot{s} @f$ such that the differential equations hold.
 *                 However, because of the algebraic equations, we need additional conditions to fully determine
 *                 @f$ \dot{s}@f$. By differentiating the algebraic equations with respect to time, we get the
 *                 missing linear equations (recall that the sensitivity vector @f$ s @f$ is fixed).
 *                
 *     Let @f$ \mathcal{I}_d @f$ denote the index set of differential equations.
 *     The right hand side of the linear system is given by @f[ -\frac{\partial F}{\partial y}(t, y, \dot{y}) s - \frac{\partial F}{\partial p}(t, y, \dot{y}), @f]
 *     which is 0 for algebraic equations (@f$ -\frac{\partial^2 F}{\partial t \partial p}@f$, to be more precise).</li>
 *          </ol>
 *     This function requires the parameter sensitivities to be computed beforehand and up-to-date Jacobians.
 * @param [in] t Current time point
 * @param [in] secIdx Index of the current section
 * @param [in] timeFactor Used for time transformation (pre factor of time derivatives) and to compute parameter derivatives with respect to section length
 * @param [in] vecStateY State vector with consistent initial values of the original system
 * @param [in] vecStateYdot Time derivative state vector with consistent initial values of the original system
 * @param [in,out] vecSensY Sensitivity subsystem state vectors
 * @param [in,out] vecSensYdot Time derivative state vectors of the sensitivity subsystems to be initialized
 * @param [in] adRes Pointer to residual vector of AD datatypes with parameter sensitivities
 * @todo Decrease amount of allocated memory by partially using temporary vectors (state and Schur complement)
 */
void LumpedRateModelWithoutPores::consistentInitialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes)
{
	BENCH_SCOPE(_timerConsistentInit);

	Indexer idxr(_disc);

	for (unsigned int param = 0; param < vecSensY.size(); ++param)
	{
		double* const sensY = vecSensY[param];
		double* const sensYdot = vecSensYdot[param];

		// Copy parameter derivative dF / dp from AD and negate it
		for (unsigned int i = _disc.nComp; i < numDofs(); ++i)
			sensYdot[i] = -adRes[i].getADValue(param);

		// Step 1: Solve algebraic equations

		if (_binding[0]->hasAlgebraicEquations())
		{
#ifdef CADET_PARALLELIZE
			BENCH_SCOPE(_timerConsistentInitPar);
			tbb::parallel_for(size_t(0), size_t(_disc.nCol), [&](size_t col)
#else
			for (unsigned int col = 0; col < _disc.nCol; ++col)
#endif
			{
				// Get algebraic block
				unsigned int algStart = 0;
				unsigned int algLen = 0;
				_binding[0]->getAlgebraicBlock(algStart, algLen);

				// Reuse memory of band matrix for dense matrix
				linalg::DenseMatrixView jacobianMatrix(_jacDisc.data() + col * _disc.strideBound * _disc.strideBound, _jacDisc.pivot() + col * _disc.strideBound, algLen, algLen);

				const unsigned int jacRowOffset = idxr.strideColCell() * col + static_cast<unsigned int>(idxr.strideColLiquid());
				const int localC = idxr.offsetC() + col * idxr.strideColCell();

				parts::BindingConsistentInitializer::consistentInitialSensitivityState(algStart, algLen, jacobianMatrix,
					_jac, localC, jacRowOffset, idxr.strideColLiquid(), _disc.strideBound,
					sensY, sensYdot, _tempState + localC);
			} CADET_PARFOR_END;
		}

		// Step 2: Compute the correct time derivative of the state vector

		// Compute right hand side by adding -dF / dy * s = -J * s to -dF / dp which is already stored in sensYdot
		multiplyWithJacobian(static_cast<double>(t), secIdx, static_cast<double>(timeFactor), vecStateY, vecStateYdot, sensY, -1.0, 1.0, sensYdot);

		// Note that we have correctly negated the right hand side

		_jacDisc.setAll(0.0);

		// Handle transport equations (dc_i / dt terms)
		_convDispOp.addTimeDerivativeToJacobian(1.0, static_cast<double>(timeFactor), _jacDisc);

		const double invBeta = 1.0 / static_cast<double>(_totalPorosity) - 1.0;
		for (unsigned int col = 0; col < _disc.nCol; ++col)
		{
			// Assemble
			linalg::FactorizableBandMatrix::RowIterator jac = _jacDisc.row(idxr.strideColCell() * col);

			// Mobile phase (advances jac accordingly)
			addMobilePhaseTimeDerivativeToJacobianCell(jac, idxr, 1.0, invBeta, static_cast<double>(timeFactor));

			// Stationary phase
			// Populate matrix with time derivative Jacobian first
			_binding[0]->jacobianAddDiscretized(static_cast<double>(timeFactor), jac);

			// Overwrite rows corresponding to algebraic equations with the Jacobian and set right hand side to 0
			if (_binding[0]->hasAlgebraicEquations())
			{
				parts::BindingConsistentInitializer::consistentInitialSensitivityTimeDerivative(_binding[0], jac,
					_jac.row(col * idxr.strideColCell() + idxr.strideColLiquid()),
					sensYdot + idxr.offsetC() + col * idxr.strideColCell() + idxr.strideColLiquid());
			}
		}

		// Precondition
		double* const scaleFactors = _tempState + idxr.offsetC();
		_jacDisc.rowScaleFactors(scaleFactors);
		_jacDisc.scaleRows(scaleFactors);

		// Factorize
		const bool result = _jacDisc.factorize();
		if (!result)
		{
			LOG(Error) << "Factorize() failed for par block";
		}

		const bool result2 = _jacDisc.solve(scaleFactors, sensYdot + idxr.offsetC());
		if (!result2)
		{
			LOG(Error) << "Solve() failed for par block";
		}
	}
}

/**
 * @brief Computes approximately / partially consistent initial values and time derivatives of sensitivity subsystems
 * @details Given the DAE \f[ F(t, y, \dot{y}) = 0, \f] and initial values \f$ y_0 \f$ and \f$ \dot{y}_0 \f$,
 *          the sensitivity system for a parameter @f$ p @f$ reads
 *          \f[ \frac{\partial F}{\partial y}(t, y, \dot{y}) s + \frac{\partial F}{\partial \dot{y}}(t, y, \dot{y}) \dot{s} + \frac{\partial F}{\partial p}(t, y, \dot{y}) = 0. \f]
 *          The initial values of this linear DAE, @f$ s_0 = \frac{\partial y_0}{\partial p} @f$ and @f$ \dot{s}_0 = \frac{\partial \dot{y}_0}{\partial p} @f$
 *          have to be consistent with the sensitivity DAE. This functions updates the initial sensitivity\f$ s_0 \f$ and overwrites the time
 *          derivative \f$ \dot{s}_0 \f$ such that they are consistent.
 *          
 *          The process follows closely the one of leanConsistentInitialConditions() and, in fact, is a linearized version of it.
 *          This is necessary because the initial conditions of the sensitivity system \f$ s_0 \f$ and \f$ \dot{s}_0 \f$ are
 *          related to the initial conditions \f$ y_0 \f$ and \f$ \dot{y}_0 \f$ of the original DAE by differentiating them
 *          with respect to @f$ p @f$: @f$ s_0 = \frac{\partial y_0}{\partial p} @f$ and @f$ \dot{s}_0 = \frac{\partial \dot{y}_0}{\partial p}. @f$
 *          <ol>
 *              <li>Keep state and time derivative vectors as they are (i.e., do not solve algebraic equations).</li>
 *              <li>Compute the time derivatives of the state @f$ \dot{y} @f$ such that the residual is 0 for the
 *              mobile phase variables.</li>
 *          </ol>
 *     This function requires the parameter sensitivities to be computed beforehand and up-to-date Jacobians.
 * @param [in] t Current time point
 * @param [in] secIdx Index of the current section
 * @param [in] timeFactor Used for time transformation (pre factor of time derivatives) and to compute parameter derivatives with respect to section length
 * @param [in] vecStateY State vector with consistent initial values of the original system
 * @param [in] vecStateYdot Time derivative state vector with consistent initial values of the original system
 * @param [in,out] vecSensY Sensitivity subsystem state vectors
 * @param [in,out] vecSensYdot Time derivative state vectors of the sensitivity subsystems to be initialized
 * @param [in] adRes Pointer to residual vector of AD datatypes with parameter sensitivities
 * @todo Decrease amount of allocated memory by partially using temporary vectors (state and Schur complement)
 */
void LumpedRateModelWithoutPores::leanConsistentInitialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes)
{
	BENCH_SCOPE(_timerConsistentInit);

	Indexer idxr(_disc);

	for (unsigned int param = 0; param < vecSensY.size(); ++param)
	{
		double* const sensY = vecSensY[param];
		double* const sensYdot = vecSensYdot[param];

		// Copy parameter derivative from AD to tempState and negate it
		for (unsigned int col = 0; col < _disc.nCol; ++col)
		{
			const unsigned int localOffset = idxr.offsetC() + col * idxr.strideColCell();
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				_tempState[localOffset + comp] = -adRes[localOffset + comp].getADValue(param);
			}
		}

		// Step 2: Compute the correct time derivative of the state vector

		// Compute right hand side by adding -dF / dy * s = -J * s to -dF / dp which is already stored in _tempState
		multiplyWithJacobian(static_cast<double>(t), secIdx, static_cast<double>(timeFactor), vecStateY, vecStateYdot, sensY, -1.0, 1.0, _tempState);

		const double invBeta = (1.0 / static_cast<double>(_totalPorosity) - 1.0) * static_cast<double>(timeFactor);
		for (unsigned int col = 0; col < _disc.nCol; ++col)
		{
			// Offset to current cell's c and q variables
			const unsigned int localOffset = idxr.offsetC() + col * idxr.strideColCell();
			const unsigned int localOffsetQ = localOffset + idxr.strideColLiquid();

			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				// dq_{i,j} / dt is assumed to be fixed, so bring it on the right hand side
				for (unsigned int i = 0; i < _disc.nBound[comp]; ++i)
				{
					_tempState[localOffset + comp] -= invBeta * sensYdot[localOffsetQ + _disc.boundOffset[comp] + i];
				}

				sensYdot[localOffset + comp] = _tempState[localOffset + comp] / static_cast<double>(timeFactor);
			}
		}
	}
}

void registerLumpedRateModelWithoutPores(std::unordered_map<std::string, std::function<IUnitOperation*(UnitOpIdx)>>& models)
{
	models[LumpedRateModelWithoutPores::identifier()] = [](UnitOpIdx uoId) { return new LumpedRateModelWithoutPores(uoId); };
	models["LRM"] = [](UnitOpIdx uoId) { return new LumpedRateModelWithoutPores(uoId); };
	models["DPFR"] = [](UnitOpIdx uoId) { return new LumpedRateModelWithoutPores(uoId); };
}

}  // namespace model

}  // namespace cadet
