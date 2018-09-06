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

#include "model/GeneralRateModel.hpp"
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

int schurComplementMultiplierGRM(void* userData, double const* x, double* z)
{
	GeneralRateModel* const grm = static_cast<GeneralRateModel*>(userData);
	return grm->schurComplementMatrixVector(x, z);
}


GeneralRateModel::GeneralRateModel(UnitOpIdx unitOpIdx) : UnitOperationBase(unitOpIdx),
	_jacP(nullptr), _jacPdisc(nullptr), _jacPF(nullptr), _jacFP(nullptr), _jacInlet(),
	_analyticJac(true), _jacobianAdDirs(0), _factorizeJacobian(false), _tempState(nullptr),
	_initC(0), _initCp(0), _initQ(0), _initState(0), _initStateDot(0)
{
}

GeneralRateModel::~GeneralRateModel() CADET_NOEXCEPT
{
	delete[] _tempState;

	delete[] _jacPF;
	delete[] _jacFP;

	delete[] _jacP;
	delete[] _jacPdisc;

	delete[] _disc.nBound;
	delete[] _disc.boundOffset;
}

unsigned int GeneralRateModel::numDofs() const CADET_NOEXCEPT
{
	// Column bulk DOFs: nCol * nComp
	// Particle DOFs: nCol particles each having nComp (liquid phase) + sum boundStates (solid phase) DOFs
	//                in each shell; there are nPar shells
	// Flux DOFs: nCol * nComp (as many as column bulk DOFs)
	// Inlet DOFs: nComp
	return _disc.nCol * (2 * _disc.nComp + _disc.nPar * (_disc.nComp + _disc.strideBound)) + _disc.nComp;
}

unsigned int GeneralRateModel::numPureDofs() const CADET_NOEXCEPT
{
	// Column bulk DOFs: nCol * nComp
	// Particle DOFs: nCol particles each having nComp (liquid phase) + sum boundStates (solid phase) DOFs
	//                in each shell; there are nPar shells
	// Flux DOFs: nCol * nComp (as many as column bulk DOFs)
	return _disc.nCol * (2 * _disc.nComp + _disc.nPar * (_disc.nComp + _disc.strideBound));
}


bool GeneralRateModel::usesAD() const CADET_NOEXCEPT
{
#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	// We always need AD if we want to check the analytical Jacobian
	return true;
#else
	// We only need AD if we are not computing the Jacobian analytically
	return !_analyticJac;
#endif
}

bool GeneralRateModel::configureModelDiscretization(IParameterProvider& paramProvider, IConfigHelper& helper)
{
	// ==== Read discretization
	_disc.nComp = paramProvider.getInt("NCOMP");

	paramProvider.pushScope("discretization");

	_disc.nCol = paramProvider.getInt("NCOL");
	_disc.nPar = paramProvider.getInt("NPAR");

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

	// Configure particle discretization
	_parCellSize.resize(_disc.nPar);
	_parCenterRadius.resize(_disc.nPar);
	_parOuterSurfAreaPerVolume.resize(_disc.nPar);
	_parInnerSurfAreaPerVolume.resize(_disc.nPar);

	// Read particle discretization mode and default to "EQUIDISTANT_PAR"
	_parDiscType = ParticleDiscretizationMode::Equidistant;
	const std::string pdt = paramProvider.getString("PAR_DISC_TYPE");
	if (pdt == "EQUIVOLUME_PAR")
		_parDiscType = ParticleDiscretizationMode::Equivolume;
	else if (pdt == "USER_DEFINED_PAR")
		_parDiscType = ParticleDiscretizationMode::UserDefined;

	if (paramProvider.exists("PAR_DISC_VECTOR"))
		_parDiscVector = paramProvider.getDoubleArray("PAR_DISC_VECTOR");

	// Determine whether analytic Jacobian should be used but don't set it right now.
	// We need to setup Jacobian matrices first.
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	const bool analyticJac = paramProvider.getBool("USE_ANALYTIC_JACOBIAN");
#else
	const bool analyticJac = false;
#endif

	// Initialize and configure GMRES for solving the Schur-complement
	_gmres.initialize(_disc.nCol * _disc.nComp, paramProvider.getInt("MAX_KRYLOV"), linalg::toOrthogonalization(paramProvider.getInt("GS_TYPE")), paramProvider.getInt("MAX_RESTARTS"));
	_gmres.matrixVectorMultiplier(&schurComplementMultiplierGRM, this);
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

	_jacP = new linalg::BandMatrix[_disc.nCol];
	_jacPdisc = new linalg::FactorizableBandMatrix[_disc.nCol];
	for (unsigned int i = 0; i < _disc.nCol; ++i)
	{
		_jacPdisc[i].resize(_disc.nPar * (_disc.nComp + _disc.strideBound), _disc.nComp + _disc.strideBound, _disc.nComp + 2 * _disc.strideBound);
		_jacP[i].resize(_disc.nPar * (_disc.nComp + _disc.strideBound), _disc.nComp + _disc.strideBound, _disc.nComp + 2 * _disc.strideBound);
	}

	_jacPF = new linalg::DoubleSparseMatrix[_disc.nCol];
	_jacFP = new linalg::DoubleSparseMatrix[_disc.nCol];
	for (unsigned int i = 0; i < _disc.nCol; ++i)
	{
		_jacPF[i].resize(_disc.nComp);
		_jacFP[i].resize(_disc.nComp);
	}

	_jacCF.resize(_disc.nComp * _disc.nCol);
	_jacFC.resize(_disc.nComp * _disc.nCol);

	_discParFlux.resize(sizeof(active) * _disc.nComp);

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
		const unsigned int requiredMem = (_binding->workspaceSize() + sizeof(double) - 1) / sizeof(double) * _disc.nPar * _disc.nCol;
		if (requiredMem > size)
		{
			size = requiredMem;
		}
	}
	_tempState = new double[size];

	return transportSuccess && bindingConfSuccess;
}

bool GeneralRateModel::configure(IParameterProvider& paramProvider)
{
	_parameters.clear();

	const bool transportSuccess = _convDispOp.configure(_unitOpIdx, paramProvider, _parameters);

	// Read geometry parameters
	_colPorosity = paramProvider.getDouble("COL_POROSITY");
	_parRadius = paramProvider.getDouble("PAR_RADIUS");
	_parPorosity = paramProvider.getDouble("PAR_POROSITY");

	// Let _parCoreRadius default to 0.0 for backwards compatibility
	if (paramProvider.exists("PAR_CORERADIUS"))
		_parCoreRadius = paramProvider.getDouble("PAR_CORERADIUS");
	else
		_parCoreRadius = 0.0;
        
	// Read vectorial parameters (which may also be section dependent; transport)
	readParameterMatrix(_filmDiffusion, paramProvider, "FILM_DIFFUSION", _disc.nComp, 1);
	readParameterMatrix(_parDiffusion, paramProvider, "PAR_DIFFUSION", _disc.nComp, 1);
	readParameterMatrix(_parSurfDiffusion, paramProvider, "PAR_SURFDIFFUSION", _disc.nComp * _disc.strideBound, 1);

	if (paramProvider.exists("PORE_ACCESSIBILITY"))
		readParameterMatrix(_poreAccessFactor, paramProvider, "PORE_ACCESSIBILITY", _disc.nComp, 1);
	else
		_poreAccessFactor = std::vector<cadet::active>(_disc.nComp, 1.0);

	// Add parameters to map
	_parameters[makeParamId(hashString("COL_POROSITY"), _unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &_colPorosity;
	_parameters[makeParamId(hashString("PAR_RADIUS"), _unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &_parRadius;
	_parameters[makeParamId(hashString("PAR_CORERADIUS"), _unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &_parCoreRadius;
	_parameters[makeParamId(hashString("PAR_POROSITY"), _unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &_parPorosity;

	// Calculate the particle radial discretization variables (_parCellSize, _parCenterRadius, etc.)
	updateRadialDisc();

	registerComponentSectionDependentParam(hashString("FILM_DIFFUSION"), _parameters, _filmDiffusion, _unitOpIdx, _disc.nComp);
	registerComponentSectionDependentParam(hashString("PAR_DIFFUSION"), _parameters, _parDiffusion, _unitOpIdx, _disc.nComp);
	registerComponentSectionDependentParam(hashString("PORE_ACCESSIBILITY"), _parameters, _poreAccessFactor, _unitOpIdx, _disc.nComp);

	// Register particle surface diffusion in this ordering:
	// sec0bnd0comp0, sec0bnd1comp0, sec0bnd2comp0, sec0bnd0comp1, sec0bnd1comp1
	// sec1bnd0comp0, sec1bnd1comp0, sec1bnd2comp0, sec1bnd0comp1, sec1bnd1comp1, ...
	if (_disc.strideBound > 0)
	{
		if (_parSurfDiffusion.size() > _disc.strideBound)
		{
			const unsigned int numSec = _parSurfDiffusion.size() / _disc.strideBound;
			unsigned int idx = 0;
			for (unsigned int sec = 0; sec < numSec; ++sec)
			{
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				{
					for (unsigned int bnd = 0; bnd < _disc.nBound[comp]; ++bnd, ++idx)
						_parameters[makeParamId(hashString("PAR_SURFDIFFUSION"), _unitOpIdx, comp, bnd, ReactionIndep, sec)] = &_parSurfDiffusion[idx];
				}
			}
		}
		else
		{
			unsigned int idx = 0;
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				for (unsigned int bnd = 0; bnd < _disc.nBound[comp]; ++bnd, ++idx)
					_parameters[makeParamId(hashString("PAR_SURFDIFFUSION"), _unitOpIdx, comp, bnd, ReactionIndep, SectionIndep)] = &_parSurfDiffusion[idx];
			}
		}
	}

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

void GeneralRateModel::useAnalyticJacobian(const bool analyticJac)
{
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	_analyticJac = analyticJac;
	if (!_analyticJac)
		// We need as many directions as the highest bandwidth of the diagonal blocks:
		// The bandwidth of the column block depends on the size of the WENO stencil, whereas
		// the bandwidth of the particle blocks are given by the number of components and bound states.
		_jacobianAdDirs = std::max(_convDispOp.requiredADdirs(), _jacP[0].stride());
	else
		_jacobianAdDirs = 0;
#else
	_analyticJac = false;
	// We need as many directions as the highest bandwidth of the diagonal blocks:
	// The bandwidth of the column block depends on the size of the WENO stencil, whereas
	// the bandwidth of the particle blocks are given by the number of components and bound states.
	_jacobianAdDirs = std::max(_convDispOp.requiredADdirs(), _jacP[0].stride());
#endif
}

void GeneralRateModel::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active* const adRes, active* const adY, unsigned int adDirOffset)
{
	// Setup flux Jacobian blocks at the beginning of the simulation or in case of
	// section dependent film or particle diffusion coefficients
	if ((secIdx == 0) || (_filmDiffusion.size() > _disc.nComp) || (_parDiffusion.size() > _disc.nComp))
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

void GeneralRateModel::setFlowRates(const active& in, const active& out) CADET_NOEXCEPT
{
	_convDispOp.setFlowRates(in, out, _colPorosity);
}

void GeneralRateModel::reportSolution(ISolutionRecorder& recorder, double const* const solution) const
{
	Exporter expr(_disc, solution);
	recorder.beginUnitOperation(_unitOpIdx, *this, expr);
	recorder.endUnitOperation();
}

void GeneralRateModel::reportSolutionStructure(ISolutionRecorder& recorder) const
{
	Exporter expr(_disc, nullptr);
	recorder.unitOperationStructure(_unitOpIdx, *this, expr);
}


unsigned int GeneralRateModel::requiredADdirs() const CADET_NOEXCEPT
{
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	return _jacobianAdDirs;
#else
	// If CADET_CHECK_ANALYTIC_JACOBIAN is active, we always need the AD directions for the Jacobian
	return std::max(_convDispOp.requiredADdirs(), _jacP[0].stride());
#endif
}

void GeneralRateModel::prepareADvectors(active* const adRes, active* const adY, unsigned int adDirOffset) const
{
	// Early out if AD is disabled
	if (!adY)
		return;

	Indexer idxr(_disc);

	// Get bandwidths of blocks
	const unsigned int lowerParBandwidth = _jacP[0].lowerBandwidth();
	const unsigned int upperParBandwidth = _jacP[0].upperBandwidth();

	// Column block	
	_convDispOp.prepareADvectors(adRes, adY, adDirOffset);

	// Particle blocks
	for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
	{
		ad::prepareAdVectorSeedsForBandMatrix(adY + idxr.offsetCp(pblk), adDirOffset, idxr.strideParBlock(), lowerParBandwidth, upperParBandwidth, lowerParBandwidth);
	}
}

/**
 * @brief Extracts the system Jacobian from band compressed AD seed vectors
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
void GeneralRateModel::extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset)
{
	Indexer idxr(_disc);

	// Column
	_convDispOp.extractJacobianFromAD(adRes, adDirOffset);

	// Particles
	for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
		ad::extractBandedJacobianFromAd(adRes + idxr.offsetCp(pblk), adDirOffset, _jacP[pblk].lowerBandwidth(), _jacP[pblk]);
}

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN

/**
 * @brief Compares the analytical Jacobian with a Jacobian derived by AD
 * @details The analytical Jacobian is assumed to be stored in the corresponding band matrices.
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
void GeneralRateModel::checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const
{
	Indexer idxr(_disc);

	LOG(Debug) << "AD dir offset: " << adDirOffset << " DiagDirCol: " << _convDispOp.jacobian().lowerBandwidth() << " DiagDirPar: " << _jacP[0].lowerBandwidth();

	// Column
	const double maxDiffCol = _convDispOp.checkAnalyticJacobianAgainstAd(adRes, adDirOffset);

	// Particles
	double maxDiffPar = 0.0;
	for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
	{
		const double localDiff = ad::compareBandedJacobianWithAd(adRes + idxr.offsetCp(pblk), adDirOffset, _jacP[pblk].lowerBandwidth(), _jacP[pblk]);
		LOG(Debug) << "-> Par block diff " << pblk << ": " << localDiff;
		maxDiffPar = std::max(maxDiffPar, localDiff);
	}
}

#endif

int GeneralRateModel::residual(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double* const res)
{
	BENCH_SCOPE(_timerResidual);

	// Evaluate residual do not compute Jacobian or parameter sensitivities
	return residualImpl<double, double, double, false>(t, secIdx, timeFactor, y, yDot, res);
}

int GeneralRateModel::residualWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, double* const res, active* const adRes, active* const adY, unsigned int adDirOffset)
{
	BENCH_SCOPE(_timerResidual);

	// Evaluate residual, use AD for Jacobian if required but do not evaluate parameter derivatives
	return residual(t, secIdx, timeFactor, y, yDot, res, adRes, adY, adDirOffset, true, false);
}

int GeneralRateModel::residual(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, double* const res, 
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
int GeneralRateModel::residualImpl(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* const y, double const* const yDot, ResidualType* const res)
{
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
int GeneralRateModel::residualParticle(const ParamType& t, unsigned int colCell, unsigned int secIdx, const ParamType& timeFactor, StateType const* yBase, double const* yDotBase, ResidualType* resBase)
{
	Indexer idxr(_disc);

	// Go to the particle block of the given column cell
	StateType const* y = yBase + idxr.offsetCp(colCell);
	double const* yDot = yDotBase + idxr.offsetCp(colCell);
	ResidualType* res = resBase + idxr.offsetCp(colCell);

	const unsigned int requiredMem = (_binding->workspaceSize() + sizeof(double) - 1) / sizeof(double);
	double* const buffer = _tempState + requiredMem * colCell;

	// Prepare parameters
	active const* const parDiff = getSectionDependentSlice(_parDiffusion, _disc.nComp, secIdx);

	// Ordering of particle surface diffusion:
	// bnd0comp0, bnd0comp1, bnd0comp2, bnd1comp0, bnd1comp1, bnd1comp2
	active const* const parSurfDiff = getSectionDependentSlice(_parSurfDiffusion, idxr.strideParBound(), secIdx);

	// Midpoint of current column cell (z coordinate) - needed in externally dependent adsorption kinetic
	const double z = (0.5 + static_cast<double>(colCell)) / static_cast<double>(_disc.nCol);

	// Reset Jacobian
	if (wantJac)
		_jacP[colCell].setAll(0.0);

	// The RowIterator is always centered on the main diagonal.
	// This means that jac[0] is the main diagonal, jac[-1] is the first lower diagonal,
	// and jac[1] is the first upper diagonal. We can also access the rows from left to
	// right beginning with the last lower diagonal moving towards the main diagonal and
	// continuing to the last upper diagonal by using the native() method.
	linalg::BandMatrix::RowIterator jac = _jacP[colCell].row(0);

	// Loop over particle cells
	for (unsigned int par = 0; par < _disc.nPar; ++par)
	{
		// Geometry
		const ParamType outerAreaPerVolume = static_cast<ParamType>(_parOuterSurfAreaPerVolume[par]);
		const ParamType innerAreaPerVolume = static_cast<ParamType>(_parInnerSurfAreaPerVolume[par]);

		// Mobile phase
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++res, ++y, ++yDot, ++jac)
		{
			*res = 0.0;
			const unsigned int nBound = _disc.nBound[comp];
			const ParamType invBetaP = (1.0 - static_cast<ParamType>(_parPorosity)) / (static_cast<ParamType>(_poreAccessFactor[comp]) * static_cast<ParamType>(_parPorosity));

			// Add time derivatives
			if (yDotBase)
			{
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

			const ParamType dp = static_cast<ParamType>(parDiff[comp]);

			// Add flow through outer surface
			// Note that inflow boundary conditions are handled in residualFlux().
			if (cadet_likely(par != 0))
			{
				// Difference between two cell-centers
				const ParamType dr = static_cast<ParamType>(_parCenterRadius[par - 1]) - static_cast<ParamType>(_parCenterRadius[par]);

				// Molecular diffusion contribution
				const ResidualType gradCp = (y[-idxr.strideParShell()] - y[0]) / dr;
				*res -= outerAreaPerVolume * dp * gradCp;

				// Surface diffusion contribution
				for (unsigned int i = 0; i < nBound; ++i)
				{
					// See above for explanation of curIdx value
					const int curIdx = idxr.strideParLiquid() - comp + idxr.offsetBoundComp(comp) + i;
					const ResidualType gradQ = (y[-idxr.strideParShell() + curIdx] - y[curIdx]) / dr;
					*res -= outerAreaPerVolume * static_cast<ParamType>(parSurfDiff[idxr.offsetBoundComp(comp) + i]) * invBetaP * gradQ;
				}

				if (wantJac)
				{
					const double localInvBetaP = static_cast<double>(invBetaP);
					const double ouApV = static_cast<double>(outerAreaPerVolume);
					const double ldr = static_cast<double>(dr);

					// Liquid phase
					jac[0] += ouApV * static_cast<double>(dp) / ldr; // dres / dc_p,i^(p,j)
					jac[-idxr.strideParShell()] = -ouApV * static_cast<double>(dp) / ldr; // dres / dc_p,i^(p,j-1)

					// Solid phase
					for (unsigned int i = 0; i < nBound; ++i)
					{
						// See above for explanation of curIdx value
						const int curIdx = idxr.strideParLiquid() - comp + idxr.offsetBoundComp(comp) + i;
						jac[curIdx] += ouApV * localInvBetaP * static_cast<double>(parSurfDiff[idxr.offsetBoundComp(comp) + i]) / ldr; // dres / dq_i^(p,j)
						jac[-idxr.strideParShell() + curIdx] = -ouApV * localInvBetaP * static_cast<double>(parSurfDiff[idxr.offsetBoundComp(comp) + i]) / ldr; // dres / dq_i^(p,j-1)
					}
				}
			}

			// Add flow through inner surface
			// Note that this term vanishes for the most inner shell due to boundary conditions
			if (cadet_likely(par != _disc.nPar - 1))
			{
				// Difference between two cell-centers
				const ParamType dr = static_cast<ParamType>(_parCenterRadius[par]) - static_cast<ParamType>(_parCenterRadius[par + 1]);

				// Molecular diffusion contribution
				const ResidualType gradCp = (y[0] - y[idxr.strideParShell()]) / dr;
				*res += innerAreaPerVolume * dp * gradCp;

				// Surface diffusion contribution
				for (unsigned int i = 0; i < nBound; ++i)
				{
					// See above for explanation of curIdx value
					const unsigned int curIdx = idxr.strideParLiquid() - comp + idxr.offsetBoundComp(comp) + i;
					const ResidualType gradQ = (y[curIdx] - y[idxr.strideParShell() + curIdx]) / dr;
					*res += innerAreaPerVolume * static_cast<ParamType>(parSurfDiff[idxr.offsetBoundComp(comp) + i]) * invBetaP * gradQ;
				}

				if (wantJac)
				{
					const double localInvBetaP = static_cast<double>(invBetaP);
					const double inApV = static_cast<double>(innerAreaPerVolume);
					const double ldr = static_cast<double>(dr);

					// Liquid phase
					jac[0] += inApV * static_cast<double>(dp) / ldr; // dres / dc_p,i^(p,j)
					jac[idxr.strideParShell()] = -inApV * static_cast<double>(dp) / ldr; // dres / dc_p,i^(p,j+1)

					// Solid phase
					for (unsigned int i = 0; i < nBound; ++i)
					{
						// See above for explanation of curIdx value
						const int curIdx = idxr.strideParLiquid() - comp + idxr.offsetBoundComp(comp) + i;
						jac[curIdx] += inApV * localInvBetaP * static_cast<double>(parSurfDiff[idxr.offsetBoundComp(comp) + i]) / ldr; // dres / dq_i^(p,j)
						jac[idxr.strideParShell() + curIdx] = -inApV * localInvBetaP * static_cast<double>(parSurfDiff[idxr.offsetBoundComp(comp) + i]) / ldr; // dres / dq_i^(p,j-1)
					}
				}
			}
		}

		// Bound phases
		if (!yDotBase)
			yDot = nullptr;

		_binding->residual(t, z, static_cast<double>(_parCenterRadius[par]) / static_cast<double>(_parRadius), secIdx, timeFactor, y, yDot, res, buffer);
		if (wantJac)
		{
			// static_cast should be sufficient here, but this statement is also analyzed when wantJac = false
			_binding->analyticJacobian(static_cast<double>(t), z, static_cast<double>(_parCenterRadius[par]) / static_cast<double>(_parRadius), secIdx, reinterpret_cast<double const*>(y), jac, buffer);
		}

		// Advance pointers over all bound states
		y += idxr.strideParBound();
		yDot += idxr.strideParBound();
		res += idxr.strideParBound();
		jac += idxr.strideParBound();
	}
	return 0;
}

template <typename StateType, typename ResidualType, typename ParamType>
int GeneralRateModel::residualFlux(const ParamType& t, unsigned int secIdx, StateType const* yBase, double const* yDotBase, ResidualType* resBase)
{
	Indexer idxr(_disc);

	const ParamType invBetaC = 1.0 / static_cast<ParamType>(_colPorosity) - 1.0;
	const ParamType epsP = static_cast<ParamType>(_parPorosity);

	active const* const filmDiff = getSectionDependentSlice(_filmDiffusion, _disc.nComp, secIdx);
	// Ordering of particle surface diffusion:
	// bnd0comp0, bnd0comp1, bnd0comp2, bnd1comp0, bnd1comp1, bnd1comp2
	active const* const parDiff = getSectionDependentSlice(_parDiffusion, _disc.nComp, secIdx);

	const ParamType surfaceToVolumeRatio = 3.0 / static_cast<ParamType>(_parRadius);
	const ParamType outerAreaPerVolume = static_cast<ParamType>(_parOuterSurfAreaPerVolume[0]);

	const ParamType jacCF_val = invBetaC * surfaceToVolumeRatio;
	const ParamType jacPF_val = -outerAreaPerVolume / epsP;

	// Discretized film diffusion kf for finite volumes
	ParamType* const kf_FV = _discParFlux.create<ParamType>(_disc.nComp);

	const ParamType absOuterShellHalfRadius = 0.5 * static_cast<ParamType>(_parCellSize[0]);
	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
	{
		kf_FV[comp] = 1.0 / (absOuterShellHalfRadius / epsP / static_cast<ParamType>(_poreAccessFactor[comp]) / static_cast<ParamType>(parDiff[comp]) + 1.0 / static_cast<ParamType>(filmDiff[comp]));
	}

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
		resCol[i] += jacCF_val * yFlux[i];

	// J_{f,0} block, adds bulk volume state c_i to flux equation
	for (unsigned int bnd = 0; bnd < _disc.nCol; ++bnd)
	{
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			const unsigned int eq = bnd * idxr.strideColCell() + comp * idxr.strideColComp();
			resFlux[eq] -= kf_FV[comp] * yCol[eq];
		}
	}

	// J_{p,f} block, implements bead boundary condition in outer bead shell equation
	for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
	{
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			const unsigned int eq = pblk * idxr.strideColCell() + comp * idxr.strideColComp();
			resPar[pblk * idxr.strideParBlock() + comp] += jacPF_val / static_cast<ParamType>(_poreAccessFactor[comp]) * yFlux[eq];
		}
	}

	// J_{f,p} block, adds outer bead shell state c_{p,i} to flux equation
	for (unsigned int pblk = 0; pblk < _disc.nCol; pblk++)
	{
		for (unsigned int comp = 0; comp < _disc.nComp; comp++)
		{
			const unsigned int eq = pblk * idxr.strideColCell() + comp * idxr.strideColComp();
			resFlux[eq] += kf_FV[comp] * yPar[comp + pblk * idxr.strideParBlock()];
		}
	}

	_discParFlux.destroy<ParamType>();
	return 0;
}

/**
 * @brief Assembles off diagonal Jacobian blocks
 * @details Assembles the fixed blocks @f$ J_{0,f}, \dots, J_{N_p,f} @f$ and @f$ J_{f,0}, \dots, J_{f, N_p}. @f$
 *          The blocks are fixed for each section.
 * @param [in] t Current time
 * @param [in] secIdx Index of the current section
 */
void GeneralRateModel::assembleOffdiagJac(double t, unsigned int secIdx)
{
	// Clear matrices for new assembly
	_jacCF.clear();
	_jacFC.clear();
	for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
	{
		_jacPF[pblk].clear();
		_jacFP[pblk].clear();
	}

	Indexer idxr(_disc);

	const double invBetaC = 1.0 / static_cast<double>(_colPorosity) - 1.0;
	const double epsP = static_cast<double>(_parPorosity);

	active const* const filmDiff = getSectionDependentSlice(_filmDiffusion, _disc.nComp, secIdx);
	// Ordering of particle diffusion:
	// sec0comp0, sec0comp1, sec0comp2, sec1comp0, sec1comp1, sec1comp2
	active const* const parDiff = getSectionDependentSlice(_parDiffusion, _disc.nComp, secIdx);

	const double surfaceToVolumeRatio = 3.0 / static_cast<double>(_parRadius);
	const double outerAreaPerVolume   = static_cast<double>(_parOuterSurfAreaPerVolume[0]);

	const double jacCF_val = invBetaC * surfaceToVolumeRatio;
	const double jacPF_val = -outerAreaPerVolume / epsP;
	const double absOuterShellHalfRadius = 0.5 * static_cast<double>(_parCellSize[0]);

	// Discretized film diffusion kf for finite volumes
	double* const kf_FV = _discParFlux.create<double>(_disc.nComp);
	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		kf_FV[comp] = 1.0 / (absOuterShellHalfRadius / epsP / static_cast<double>(_poreAccessFactor[comp]) / static_cast<double>(parDiff[comp]) + 1.0 / static_cast<double>(filmDiff[comp]));

	// Note that the J_f block, which is the identity matrix, is treated in the linear solver

	// J_{0,f} block, adds flux to column void / bulk volume equations
	for (unsigned int eq = 0; eq < _disc.nCol * _disc.nComp; ++eq)
	{
		// Main diagonal corresponds to j_{f,i} (flux) state variable
		_jacCF.addElement(eq, eq, jacCF_val);
	}

	// J_{f,0} block, adds bulk volume state c_i to flux equation
	for (unsigned int bnd = 0; bnd < _disc.nCol; ++bnd)
	{
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			// Main diagonal corresponds to c_i state variable in each column cell
			const unsigned int eq = bnd * idxr.strideColCell() + comp * idxr.strideColComp();
			_jacFC.addElement(eq, eq, -kf_FV[comp]);
		}
	}

	// J_{p,f} block, implements bead boundary condition in outer bead shell equation
	for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
	{
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			const unsigned int eq = pblk * idxr.strideColCell() + comp * idxr.strideColComp();
			_jacPF[pblk].addElement(comp, eq, jacPF_val / static_cast<double>(_poreAccessFactor[comp]));
		}
	}

	// J_{f,p} block, adds outer bead shell state c_{p,i} to flux equation
	for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
	{
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		{
			const unsigned int eq = pblk * idxr.strideColCell() + comp * idxr.strideColComp();
			_jacFP[pblk].addElement(eq, comp, kf_FV[comp]);
		}
	}

	_discParFlux.destroy<double>();
}

int GeneralRateModel::residualSensFwdWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot,
	active* const adRes, active* const adY, unsigned int adDirOffset)
{
	BENCH_SCOPE(_timerResidualSens);

	// Evaluate residual for all parameters using AD in vector mode and at the same time update the 
	// Jacobian (in one AD run, if analytic Jacobians are disabled)
	return residual(t, secIdx, timeFactor, y, yDot, nullptr, adRes, adY, adDirOffset, true, true);
}

int GeneralRateModel::residualSensFwdAdOnly(const active& t, unsigned int secIdx, const active& timeFactor,
	double const* const y, double const* const yDot, active* const adRes)
{
	BENCH_SCOPE(_timerResidualSens);

	// Evaluate residual for all parameters using AD in vector mode
	return residualImpl<double, active, active, false>(t, secIdx, timeFactor, y, yDot, adRes); 
}

int GeneralRateModel::residualSensFwdCombine(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, 
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
void GeneralRateModel::multiplyWithJacobian(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double const* yS, double alpha, double beta, double* ret)
{
	Indexer idxr(_disc);

	// Handle identity matrix of inlet DOFs
	for (unsigned int i = 0; i < _disc.nComp; ++i)
	{
		ret[i] = alpha * yS[i] + beta * ret[i];
	}

#ifdef CADET_PARALLELIZE
	tbb::parallel_for(size_t(0), size_t(_disc.nCol + 1), [&](size_t idx)
#else
	for (unsigned int idx = 0; idx < _disc.nCol + 1; ++idx)
#endif
	{
		if (cadet_unlikely(idx == 0))
		{
			_convDispOp.jacobian().multiplyVector(yS + idxr.offsetC(), alpha, beta, ret + idxr.offsetC());
			_jacCF.multiplyVector(yS + idxr.offsetJf(), alpha, 1.0, ret + idxr.offsetC());
		}
		else
		{
			const unsigned int pblk = idx - 1;
			const int localOffset = idxr.offsetCp(pblk);
			_jacP[pblk].multiplyVector(yS + localOffset, alpha, beta, ret + localOffset);
			_jacPF[pblk].multiplyVector(yS + idxr.offsetJf(), alpha, 1.0, ret + localOffset);
		}
	} CADET_PARFOR_END;

	// Handle flux equation

	// Set fluxes(ret) = fluxes(yS)
	// This applies the identity matrix in the bottom right corner of the Jaocbian (flux equation)
	for (unsigned int i = idxr.offsetJf(); i < numDofs(); ++i)
		ret[i] = alpha * yS[i] + beta * ret[i];

	double* const retJf = ret + idxr.offsetJf();
	_jacFC.multiplyVector(yS + idxr.offsetC(), alpha, 1.0, retJf);

	for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
		_jacFP[pblk].multiplyVector(yS + idxr.offsetCp(pblk), alpha, 1.0, retJf);

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
void GeneralRateModel::multiplyWithDerivativeJacobian(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double const* sDot, double* ret)
{
	Indexer idxr(_disc);
	const double invBetaP = (1.0 / static_cast<double>(_parPorosity) - 1.0) * timeFactor;

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
			for (unsigned int shell = 0; shell < _disc.nPar; ++shell)
			{
				double const* const localSdot = sDot + idxr.offsetCp(pblk) + shell * idxr.strideParShell();
				double* const localRet = ret + idxr.offsetCp(pblk) + shell * idxr.strideParShell();

				// Mobile phase
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				{
					// Add derivative with respect to dc_p / dt to Jacobian
					localRet[comp] = timeFactor * localSdot[comp];

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
		}
	} CADET_PARFOR_END;

	// Handle fluxes (all algebraic)
	double* const dFdyDot = ret + idxr.offsetJf();
	std::fill(dFdyDot, dFdyDot + _disc.nCol * _disc.nComp, 0.0);

	// Handle inlet DOFs (all algebraic)
	std::fill_n(ret, _disc.nComp, 0.0);
}

void GeneralRateModel::setExternalFunctions(IExternalFunction** extFuns, unsigned int size)
{
	if (_binding)
		_binding->setExternalFunctions(extFuns, size);
}

unsigned int GeneralRateModel::localOutletComponentIndex() const CADET_NOEXCEPT
{
	// Inlets are duplicated so need to be accounted for
	if (static_cast<double>(_convDispOp.currentVelocity()) >= 0.0)
		// Forward Flow: outlet is last cell
		return _disc.nComp + (_disc.nCol - 1) * _disc.nComp;
	else
		// Backward flow: Outlet is first cell
		return _disc.nComp;
}

unsigned int GeneralRateModel::localInletComponentIndex() const CADET_NOEXCEPT
{
	return 0;
}

unsigned int GeneralRateModel::localOutletComponentStride() const CADET_NOEXCEPT
{
	return 1;
}

unsigned int GeneralRateModel::localInletComponentStride() const CADET_NOEXCEPT
{
	return 1;
}

void GeneralRateModel::expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut)
{
	// @todo Write this function
}

/**
 * @brief Computes equidistant radial nodes in the beads
 */
void GeneralRateModel::setEquidistantRadialDisc()
{
	const active radius = _parRadius - _parCoreRadius;
	const active dr = radius / static_cast<double>(_disc.nPar);
	_parCellSize.assign(_disc.nPar, dr);

	for (unsigned int cell = 0; cell < _disc.nPar; cell++)
	{
		const active r_out = radius - static_cast<active>(cell) * dr;
		const active r_in = radius - static_cast<active>(cell + 1) * dr;

		_parCenterRadius[cell] = radius - (0.5 + static_cast<active>(cell)) * dr;

		// Compute denominator -> corresponding to cell volume
		const active vol = pow(r_out, 3.0) - pow(r_in, 3.0);

		_parOuterSurfAreaPerVolume[cell] = 3.0 * sqr(r_out) / vol;
		_parInnerSurfAreaPerVolume[cell] = 3.0 * sqr(r_in) / vol;
	}
}

/**
 * @brief Computes the radial nodes in the beads in such a way that all shells have the same volume
 */
void GeneralRateModel::setEquivolumeRadialDisc()
{
	active r_out = _parRadius;
	active r_in = _parCoreRadius;
	const active volumePerShell = (pow(_parRadius, 3.0) - pow(_parCoreRadius, 3.0)) / static_cast<double>(_disc.nPar);

	for (unsigned int cell = 0; cell < _disc.nPar; ++cell)
	{
		if (cell != (_disc.nPar - 1))
			r_in = pow(pow(r_out, 3.0) - volumePerShell, (1.0 / 3.0));
		else
			r_in = _parCoreRadius;

		_parCellSize[cell] = r_out - r_in;
		_parCenterRadius[cell] = (r_out + r_in) * 0.5;

		_parOuterSurfAreaPerVolume[cell] = 3.0 * sqr(r_out) / volumePerShell;
		_parInnerSurfAreaPerVolume[cell] = 3.0 * sqr(r_in) / volumePerShell;

		// For the next cell: r_out == r_in of the current cell
		r_out = r_in;
	}
}

/**
 * @brief Computes all helper quantities for radial bead discretization from given radial cell boundaries
 * @details Calculates surface areas per volume for every shell and the radial shell centers.
 */
void GeneralRateModel::setUserdefinedRadialDisc()
{
	// Care for the right ordering and include 0.0 / 1.0 if not already in the vector.
	std::vector<double> orderedInterfaces = _parDiscVector;

	if (std::find(orderedInterfaces.begin(), orderedInterfaces.end(), 0.0) == orderedInterfaces.end())
		orderedInterfaces.push_back(0.0);
	if (std::find(orderedInterfaces.begin(), orderedInterfaces.end(), 1.0) == orderedInterfaces.end())
		orderedInterfaces.push_back(1.0);

	// Sort in descending order
	std::sort(orderedInterfaces.begin(), orderedInterfaces.end(), std::greater<double>());

	// Map [0, 1] -> [core radius, particle radius] via linear interpolation
	for (unsigned int cell = 0; cell < _disc.nPar; ++cell)
		orderedInterfaces[cell] = orderedInterfaces[cell] * (static_cast<double>(_parRadius) - static_cast<double>(_parCoreRadius)) + static_cast<double>(_parCoreRadius);

	for (unsigned int cell = 0; cell < _disc.nPar; ++cell)
	{
		_parCellSize[cell] = orderedInterfaces[cell] - orderedInterfaces[cell + 1];
		_parCenterRadius[cell] = (orderedInterfaces[cell] + orderedInterfaces[cell + 1]) * 0.5;

		// Compute denominator -> corresponding to cell volume
		const active vol = std::pow(orderedInterfaces[cell], 3.0) - std::pow(orderedInterfaces[cell + 1], 3.0);

		_parOuterSurfAreaPerVolume[cell] = 3.0 * sqr(orderedInterfaces[cell]) / vol;
		_parInnerSurfAreaPerVolume[cell] = 3.0 * sqr(orderedInterfaces[cell + 1]) / vol;
	}
}

void GeneralRateModel::updateRadialDisc()
{
	if (_parDiscType == ParticleDiscretizationMode::Equidistant)
		setEquidistantRadialDisc();
	else if (_parDiscType == ParticleDiscretizationMode::Equivolume)
		setEquivolumeRadialDisc();
	else if (_parDiscType == ParticleDiscretizationMode::UserDefined)
		setUserdefinedRadialDisc();
}

bool GeneralRateModel::setParameter(const ParameterId& pId, double value)
{
	const bool result = UnitOperationBase::setParameter(pId, value);

	// Check whether particle radius or core radius has changed and update radial discretization if necessary
	if (result && ((pId.name == hashString("PAR_RADIUS")) || (pId.name == hashString("PAR_CORERADIUS"))))
		updateRadialDisc();

	return result;
}

void GeneralRateModel::setSensitiveParameterValue(const ParameterId& pId, double value)
{
	UnitOperationBase::setSensitiveParameterValue(pId, value);

	// Check whether particle radius or core radius has changed and update radial discretization if necessary
	if ((pId.name == hashString("PAR_RADIUS")) || (pId.name == hashString("PAR_CORERADIUS")))
		updateRadialDisc();
}

bool GeneralRateModel::setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue)
{
	const bool result = UnitOperationBase::setSensitiveParameter(pId, adDirection, adValue);

	// Check whether particle radius or core radius has been set active and update radial discretization if necessary
	// Note that we need to recompute the radial discretization variables (_parCellSize, _parCenterRadius, _parOuterSurfAreaPerVolume, _parInnerSurfAreaPerVolume)
	// because their gradient has changed (although their nominal value has not changed).
	if ((pId.name == hashString("PAR_RADIUS")) || (pId.name == hashString("PAR_CORERADIUS")))
		updateRadialDisc();

	return result;
}

void registerGeneralRateModel(std::unordered_map<std::string, std::function<IUnitOperation*(UnitOpIdx)>>& models)
{
	models[GeneralRateModel::identifier()] = [](UnitOpIdx uoId) { return new GeneralRateModel(uoId); };
	models["GRM"] = [](UnitOpIdx uoId) { return new GeneralRateModel(uoId); };
}

}  // namespace model

}  // namespace cadet
