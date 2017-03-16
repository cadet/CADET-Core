// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//
//  Copyright © 2008-2016: The CADET Authors
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
#include "linalg/DenseMatrix.hpp"
#include "linalg/BandMatrix.hpp"
#include "linalg/Norms.hpp"

#include "Stencil.hpp"
#include "Weno.hpp"
#include "AdUtils.hpp"
#include "ParamIdUtil.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <algorithm>
#include <functional>

#include "OpenMPSupport.hpp"

namespace
{
	template <class Elem_t>
	inline bool contains(const typename std::vector<Elem_t>& vec, const Elem_t& item)
	{
		const typename std::vector<Elem_t>::const_iterator it = std::find(vec.begin(), vec.end(), item);
		return it != vec.end();
	}

	template <class Elem_t>
	inline bool contains(const typename std::unordered_set<Elem_t>& set, const Elem_t& item)
	{
		return set.find(item) != set.end();
	}
}

namespace cadet
{

namespace model
{

int schurComplementMultiplier(void* userData, double const* x, double* z)
{
	GeneralRateModel* const grm = static_cast<GeneralRateModel*>(userData);
	return grm->schurComplementMatrixVector(x, z);
}


GeneralRateModel::GeneralRateModel(UnitOpIdx unitOpIdx) : _unitOpIdx(unitOpIdx), _binding(nullptr),
	_jacC(nullptr), _jacP(nullptr), _jacPF(nullptr), _jacFP(nullptr), _jacCdisc(nullptr), _jacPdisc(nullptr),
	_analyticJac(true), _stencilMemory(sizeof(active) * Weno::maxStencilSize()), _wenoDerivatives(new double[Weno::maxStencilSize()]),
	_weno(), _jacobianAdDirs(0), _factorizeJacobian(false), _tempState(nullptr)
{

}

GeneralRateModel::~GeneralRateModel() CADET_NOEXCEPT
{
	delete[] _tempState;

	delete[] _wenoDerivatives;

	delete[] _jacPF;
	delete[] _jacFP;

	delete[] _jacC;
	delete[] _jacCdisc;

	delete[] _jacP;
	delete[] _jacPdisc;

	delete _binding;

	delete[] _disc.nBound;
	delete[] _disc.boundOffset;
}

unsigned int GeneralRateModel::numDofs() const CADET_NOEXCEPT
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

bool GeneralRateModel::configure(IParameterProvider& paramProvider, IConfigHelper& helper)
{
	// ==== Read discretization
	_disc.nComp = paramProvider.getInt("NCOMP");

	paramProvider.pushScope("discretization");

	_disc.nCol = paramProvider.getInt("NCOL");
	_disc.nPar = paramProvider.getInt("NPAR");

	const std::vector<int> nBound = paramProvider.getIntArray("NBOUND");
	_disc.nBound = new unsigned int[_disc.nComp];
	std::copy(nBound.begin(), nBound.end(), _disc.nBound);

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

	const std::string parDiscType = paramProvider.getString("PAR_DISC_TYPE");
	if (parDiscType == "EQUIVOLUME_PAR")
		setEquivolumeRadialDisc();
	else if (parDiscType == "USER_DEFINED_PAR")
	{
		const std::vector<double> parInterfaces = paramProvider.getDoubleArray("PAR_DISC_VECTOR");
		setUserdefinedRadialDisc(parInterfaces);
	}
	else // Handle parDiscType == "EQUIDISTANT_PAR" and default
		setEquidistantRadialDisc();

	// Read WENO settings and apply them
	paramProvider.pushScope("weno");
	_weno.order(paramProvider.getInt("WENO_ORDER"));
	_weno.boundaryTreatment(paramProvider.getInt("BOUNDARY_MODEL"));
	_wenoEpsilon = paramProvider.getDouble("WENO_EPS");
	paramProvider.popScope();

	// Determine whether analytic Jacobian should be used but don't set it right now.
	// We need to setup Jacobian matrices first.
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	const bool analyticJac = paramProvider.getInt("USE_ANALYTIC_JACOBIAN");
#else
	const bool analyticJac = false;
#endif

	// Initialize and configure GMRES for solving the Schur-complement
	_gmres.initialize(_disc.nCol * _disc.nComp, paramProvider.getInt("MAX_KRYLOV"), linalg::toOrthogonalization(paramProvider.getInt("GS_TYPE")), paramProvider.getInt("MAX_RESTARTS"));
	_gmres.matrixVectorMultiplier(&schurComplementMultiplier, this);
	_schurSafety = paramProvider.getDouble("SCHUR_SAFETY");

	paramProvider.popScope();

	// ==== Read model parameters
	reconfigure(paramProvider);

	// Allocate memory
	Indexer idxr(_disc);

	_jacC = new linalg::BandMatrix[_disc.nComp];
	_jacCdisc = new linalg::FactorizableBandMatrix[_disc.nComp];
	for (unsigned int i = 0; i < _disc.nComp; ++i)
	{
		// Note that we have to increase the lower bandwidth by 1 because the WENO stencil is applied to the
		// right cell face (lower + 1 + upper) and to the left cell face (shift the stencil by -1 because influx of cell i
		// is outflux of cell i-1)
		_jacC[i].resize(_disc.nCol, _weno.lowerBandwidth() + 1, _weno.upperBandwidth());
		_jacCdisc[i].resize(_disc.nCol, _weno.lowerBandwidth() + 1, _weno.upperBandwidth());
	}

	_jacP = new linalg::BandMatrix[_disc.nCol];
	_jacPdisc = new linalg::FactorizableBandMatrix[_disc.nCol];
	for (unsigned int i = 0; i < _disc.nCol; ++i)
	{
		_jacPdisc[i].resize(_disc.nPar * (_disc.nComp + _disc.strideBound), _disc.nComp + _disc.strideBound, _disc.nComp + 2 * _disc.strideBound);
		_jacP[i].resize(_disc.nPar * (_disc.nComp + _disc.strideBound), _disc.nComp + _disc.strideBound, _disc.nComp + 2 * _disc.strideBound);
	}

	_jacPF = new linalg::SparseMatrix[_disc.nCol];
	_jacFP = new linalg::SparseMatrix[_disc.nCol];
	for (unsigned int i = 0; i < _disc.nCol; ++i)
	{
		_jacPF[i].resize(_disc.nComp);
		_jacFP[i].resize(_disc.nComp);
	}

	_jacCF.resize(_disc.nComp * _disc.nCol);
	_jacFC.resize(_disc.nComp * _disc.nCol);

	_discParFlux.resize(sizeof(active) * _disc.nComp);

	_tempState = new double[numDofs()];

	// Set whether analytic Jacobian is used
	useAnalyticJacobian(analyticJac);

	// ==== Construct and configure binding model
	delete _binding;

	_binding = helper.createBindingModel(paramProvider.getString("ADSORPTION_MODEL"));
	if (!_binding)
		throw InvalidParameterException("Unknown binding model " + paramProvider.getString("ADSORPTION_MODEL"));

	_binding->configureModelDiscretization(_disc.nComp, _disc.nBound, _disc.boundOffset);

	paramProvider.pushScope("adsorption");
	const bool bindingConfSuccess = _binding->configure(paramProvider, _unitOpIdx);
	paramProvider.popScope();

	return bindingConfSuccess;
}

bool GeneralRateModel::reconfigure(IParameterProvider& paramProvider)
{
	// Read geometry parameters
	_colLength = paramProvider.getDouble("COL_LENGTH");
	_colPorosity = paramProvider.getDouble("COL_POROSITY");
	_parRadius = paramProvider.getDouble("PAR_RADIUS");
	_parPorosity = paramProvider.getDouble("PAR_POROSITY");

	// Read section dependent parameters (transport)
	readScalarParameterOrArray(_colDispersion, paramProvider, "COL_DISPERSION", 1);
	readScalarParameterOrArray(_velocity, paramProvider, "VELOCITY", 1);

	// Read vectorial parameters (which may also be section dependent; transport)
	readParameterMatrix(_filmDiffusion, paramProvider, "FILM_DIFFUSION", _disc.nComp, 1);
	readParameterMatrix(_parDiffusion, paramProvider, "PAR_DIFFUSION", _disc.nComp, 1);
	readParameterMatrix(_parSurfDiffusion, paramProvider, "PAR_SURFDIFFUSION", _disc.nComp * _disc.strideBound, 1);

	// Add parameters to map
	_parameters.clear();
	_parameters[makeParamId(hashString("COL_LENGTH"), _unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &_colLength;
	_parameters[makeParamId(hashString("COL_POROSITY"), _unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &_colPorosity;
	_parameters[makeParamId(hashString("PAR_RADIUS"), _unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &_parRadius;
	_parameters[makeParamId(hashString("PAR_POROSITY"), _unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &_parPorosity;

	registerScalarSectionDependentParam(hashString("COL_DISPERSION"), _parameters, _colDispersion, _unitOpIdx);
	registerScalarSectionDependentParam(hashString("VELOCITY"), _parameters, _velocity, _unitOpIdx);

	registerComponentSectionDependentParam(hashString("FILM_DIFFUSION"), _parameters, _filmDiffusion, _unitOpIdx, _disc.nComp);
	registerComponentSectionDependentParam(hashString("PAR_DIFFUSION"), _parameters, _parDiffusion, _unitOpIdx, _disc.nComp);

	// Register particle surface diffusion in this ordering:
	// sec0bnd0comp0, sec0bnd1comp0, sec0bnd2comp0, sec0bnd0comp1, sec0bnd1comp1
	// sec1bnd0comp0, sec1bnd1comp0, sec1bnd2comp0, sec1bnd0comp1, sec1bnd1comp1, ...
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

	// Reconfigure binding model
	if (_binding)
		return _binding->reconfigure(paramProvider, _unitOpIdx);

	return true;
}

std::unordered_map<ParameterId, double> GeneralRateModel::getAllParameterValues() const
{
	std::unordered_map<ParameterId, double> data;
	std::transform(_parameters.begin(), _parameters.end(), std::inserter(data, data.end()),
	               [](const std::pair<const ParameterId, active*>& p) { return std::make_pair(p.first, static_cast<double>(*p.second)); });

	if (!_binding)
		return data;

	const std::unordered_map<ParameterId, double> localData = _binding->getAllParameterValues();
	for (const std::pair<ParameterId, double>& val : localData)
		data[val.first] = val.second;

	return data;
}

bool GeneralRateModel::hasParameter(const ParameterId& pId) const
{
	const bool hasParam = _parameters.find(pId) != _parameters.end();
	if (_binding)
		return hasParam || _binding->hasParameter(pId);
	return hasParam;
}

bool GeneralRateModel::setParameter(const ParameterId& pId, int value)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return false;

	if (_binding)
		return _binding->setParameter(pId, value);
	return false;
}

bool GeneralRateModel::setParameter(const ParameterId& pId, double value)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return false;

	auto paramHandle = _parameters.find(pId);
	if (paramHandle != _parameters.end())
	{
		paramHandle->second->setValue(value);
		return true;
	}
	else if (_binding)
		return _binding->setParameter(pId, value);

	return false;
}

bool GeneralRateModel::setParameter(const ParameterId& pId, bool value)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return false;

	if (_binding)
		return _binding->setParameter(pId, value);
	return false;
}

void GeneralRateModel::setSensitiveParameterValue(const ParameterId& pId, double value)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return;

	// Check our own parameters
	auto paramHandle = _parameters.find(pId);
	if ((paramHandle != _parameters.end()) && contains(_sensParams, paramHandle->second))
	{
		paramHandle->second->setValue(value);
		return;
	}

	// Check binding model parameters
	if (_binding)
	{
		active* const val = _binding->getParameter(pId);
		if (val && contains(_sensParams, val))
		{
			val->setValue(value);
			return;
		}
	}
}

bool GeneralRateModel::setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return false;

	// Check own parameters
	auto paramHandle = _parameters.find(pId);
	if (paramHandle != _parameters.end())
	{
		LOG(Debug) << "Found parameter " << pId << " in GRM: Dir " << adDirection << " is set to " << adValue;

		// Register parameter and set AD seed / direction
		_sensParams.insert(paramHandle->second);
		paramHandle->second->setADValue(adDirection, adValue);
		return true;
	}

	// Check binding model parameters
	if (_binding)
	{
		active* const paramBinding = _binding->getParameter(pId);
		if (paramBinding)
		{
			LOG(Debug) << "Found parameter " << pId << " in AdsorptionModel: Dir " << adDirection << " is set to " << adValue;

			// Register parameter and set AD seed / direction
			_sensParams.insert(paramBinding);
			paramBinding->setADValue(adDirection, adValue);
			return true;
		}
	}

	return false;
}

void GeneralRateModel::clearSensParams()
{
	// Remove AD directions from parameters
	for (auto sp : _sensParams)
		sp->setADValue(0.0);

	_sensParams.clear();
}

void GeneralRateModel::useAnalyticJacobian(const bool analyticJac)
{
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	_analyticJac = analyticJac;
	if (!_analyticJac)
		// We need as many directions as the highest bandwidth of the diagonal blocks:
		// The bandwidth of the column block depends on the size of the WENO stencil, whereas
		// the bandwidth of the particle blocks are given by the number of components and bound states.
		_jacobianAdDirs = std::max(_jacC[0].stride(), _jacP[0].stride());
	else
		_jacobianAdDirs = 0;
#else
	_analyticJac = false;
	// We need as many directions as the highest bandwidth of the diagonal blocks:
	// The bandwidth of the column block depends on the size of the WENO stencil, whereas
	// the bandwidth of the particle blocks are given by the number of components and bound states.
	_jacobianAdDirs = std::max(_jacC[0].stride(), _jacP[0].stride());
#endif
}

void GeneralRateModel::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx)
{
	// Setup flux Jacobian blocks at the beginning of the simulation
	if (secIdx == 0)
		assembleOffdiagJac(t, secIdx);

	// Flux Jacobian blocks only change for section dependent film or particle diffusion coefficients
	if ((_filmDiffusion.size() > _disc.nComp) || (_parDiffusion.size() > _disc.nComp))
		assembleOffdiagJac(t, secIdx);
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
	return std::max(_jacC[0].stride(), _jacP[0].stride());
#endif
}

void GeneralRateModel::prepareADvectors(active* const adRes, active* const adY, unsigned int numSensAdDirs) const
{
	// Early out if AD is disabled
	if (!adY)
		return;

	Indexer idxr(_disc);

	// Get bandwidths of blocks
	const unsigned int lowerColBandwidth = _jacC[0].lowerBandwidth();
	const unsigned int upperColBandwidth = _jacC[0].upperBandwidth();
	const unsigned int lowerParBandwidth = _jacP[0].lowerBandwidth();
	const unsigned int upperParBandwidth = _jacP[0].upperBandwidth();

	// Column block	
	for (int comp = 0; comp < static_cast<int>(_disc.nComp); ++comp)
	{
		ad::prepareAdVectorSeedsForBandMatrix(adY + comp * idxr.strideColComp(), numSensAdDirs, _disc.nCol, lowerColBandwidth, upperColBandwidth, lowerColBandwidth);
	}

	// Particle blocks
	for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
	{
		ad::prepareAdVectorSeedsForBandMatrix(adY + idxr.offsetCp(pblk), numSensAdDirs, idxr.strideParBlock(), lowerParBandwidth, upperParBandwidth, lowerParBandwidth);
	}
}

/**
 * @brief Extracts the system Jacobian from band compressed AD seed vectors
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] numSensAdDirs Number of AD directions used for parameter sensitivities
 */
void GeneralRateModel::extractJacobianFromAD(active const* const adRes, unsigned int numSensAdDirs)
{
	Indexer idxr(_disc);

	// Column
	for (int comp = 0; comp < static_cast<int>(_disc.nComp); ++comp)
		ad::extractBandedJacobianFromAd(adRes + comp * idxr.strideColComp(), numSensAdDirs, _jacC[comp].lowerBandwidth(), _jacC[comp]);

	// Particles
	for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
		ad::extractBandedJacobianFromAd(adRes + idxr.offsetCp(pblk), numSensAdDirs, _jacP[pblk].lowerBandwidth(), _jacP[pblk]);
}

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN

/**
 * @brief Compares the analytical Jacobian with a Jacobian derived by AD
 * @details The analytical Jacobian is assumed to be stored in the corresponding band matrices.
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] numSensAdDirs Number of AD directions used for parameter sensitivities
 */
void GeneralRateModel::checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int numSensAdDirs) const
{
	Indexer idxr(_disc);

	double maxDiffCol = 0.0;
	double maxDiffPar = 0.0;

	LOG(Debug) << "AD dir offset: " << numSensAdDirs << " DiagDirCol: " << _jacC[0].lowerBandwidth() << " DiagDirPar: " << _jacP[0].lowerBandwidth();

	// Column
	for (int comp = 0; comp < static_cast<int>(_disc.nComp); ++comp)
	{
		const double localDiff = ad::compareBandedJacobianWithAd(adRes + comp * idxr.strideColComp(), numSensAdDirs, _jacC[comp].lowerBandwidth(), _jacC[comp]);
		LOG(Debug) << "-> Col block diff " << comp << ": " << localDiff;
		maxDiffCol = std::max(maxDiffCol, localDiff);
	}

	// Particles
	for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
	{
		const double localDiff = ad::compareBandedJacobianWithAd(adRes + idxr.offsetCp(pblk), numSensAdDirs, _jacP[pblk].lowerBandwidth(), _jacP[pblk]);
		LOG(Debug) << "-> Par block diff " << pblk << ": " << localDiff;
		maxDiffPar = std::max(maxDiffPar, localDiff);
	}
}

#endif

int GeneralRateModel::residual(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double* const res)
{
	LOG(Trace) << "======= RESIDUAL ========== t = " << static_cast<double>(t) << " sec = " << secIdx << " dt = " << static_cast<double>(timeFactor);
	BENCH_SCOPE(_timerResidual);

	// Evaluate residual do not compute Jacobian or parameter sensitivities
	return residualImpl<double, double, double, false>(t, secIdx, timeFactor, y, yDot, res);
}

int GeneralRateModel::residualWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, double* const res, active* const adRes, active* const adY, unsigned int numSensAdDirs)
{
	LOG(Trace) << "======= RESIDUAL ========== t = " << static_cast<double>(t) << " sec = " << secIdx << " dt = " << static_cast<double>(timeFactor);
	BENCH_SCOPE(_timerResidual);

	// Evaluate residual, use AD for Jacobian if required but do not evaluate parameter derivatives
	return residual(t, secIdx, timeFactor, y, yDot, res, adRes, adY, numSensAdDirs, true, false);
}

int GeneralRateModel::residual(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, double* const res, 
	active* const adRes, active* const adY, unsigned int numSensAdDirs, bool updateJacobian, bool paramSensitivity)
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
			extractJacobianFromAD(adRes, numSensAdDirs);

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
			checkAnalyticJacobianAgainstAd(adRes, numSensAdDirs);
		}

		// Extract Jacobian
		extractJacobianFromAD(adRes, numSensAdDirs);

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

double GeneralRateModel::residualNorm(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot)
{
	// We use the _tempState vector to store the residual
	residualImpl<double, double, double, false>(t, secIdx, timeFactor, y, yDot, _tempState);

//	printStateVector("Residual", _tempState, _disc, Indexer(_disc));
//	printVector("Consistency residual", _tempState, numDofs());

	return linalg::linfNorm(_tempState, numDofs());
}

template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
int GeneralRateModel::residualImpl(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* const y, double const* const yDot, ResidualType* const res)
{
	LOG(Debug) << "t = " << t << " timeFactor = " << timeFactor;

	BENCH_START(_timerResidualPar);

	#pragma omp parallel for schedule(static)
	for (ompuint_t pblk = 0; pblk <= _disc.nCol; ++pblk)
	{
		if (cadet_unlikely(pblk == 0))
			residualBulk<StateType, ResidualType, ParamType, wantJac>(t, secIdx, timeFactor, y, yDot, res);
		else
			residualParticle<StateType, ResidualType, ParamType, wantJac>(t, pblk-1, secIdx, timeFactor, y, yDot, res);
	}

	BENCH_STOP(_timerResidualPar);

	residualFlux<StateType, ResidualType, ParamType>(t, secIdx, y, yDot, res);

	return 0;
}

template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
int GeneralRateModel::residualBulk(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* y, double const* yDot, ResidualType* res)
{
	const ParamType u = static_cast<ParamType>(getSectionDependentScalar(_velocity, secIdx));
	const ParamType d_c = static_cast<ParamType>(getSectionDependentScalar(_colDispersion, secIdx));
	const ParamType h = static_cast<ParamType>(_colLength) / static_cast<double>(_disc.nCol);
	const ParamType h2 = h * h;

	Indexer idxr(_disc);

	// The stencil caches parts of the state vector for better spatial coherence
	typedef CachingStencil<StateType, ArrayPool> StencilType;
	StencilType stencil(_weno.stencilSize(), _stencilMemory, _weno.order() - 1);

	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
	{
		// Reset Jacobian
		if (wantJac)
			_jacC[comp].setAll(0.0);

		// The RowIterator is always centered on the main diagonal.
		// This means that jac[0] is the main diagonal, jac[-1] is the first lower diagonal,
		// and jac[1] is the first upper diagonal. We can also access the rows from left to
		// right beginning with the last lower diagonal moving towards the main diagonal and
		// continuing to the last upper diagonal by using the native() method.
		linalg::BandMatrix::RowIterator jac = _jacC[comp].row(0);

		// Add time derivative to each cell
		if (yDot)
		{
			for (unsigned int col = 0; col < _disc.nCol; ++col)
				idxr.c<ResidualType>(res, col, comp) = timeFactor * idxr.c<double>(yDot, col, comp);
		}
		else
		{
			for (unsigned int col = 0; col < _disc.nCol; ++col)
				idxr.c<ResidualType>(res, col, comp) = 0.0;
		}

		// Fill stencil (left side with zeros, right side with states)
		for (int i = -_weno.order() + 1; i < 0; ++i)
			stencil[i] = 0.0;
		for (int i = 0; i < _weno.order(); ++i)
			stencil[i] = idxr.c<StateType>(y, static_cast<unsigned int>(i), comp);

		// Reset WENO output
		StateType vm(0.0); // reconstructed value
		for (unsigned int i = 0; i < _weno.stencilSize(); ++i)
			_wenoDerivatives[i] = 0.0;

		int wenoOrder = 0;

		// Iterate over all cells
		for (unsigned int col = 0; col < _disc.nCol; ++col)
		{
			// ------------------- Dispersion -------------------

			// Right side, leave out if we're in the last cell (boundary condition)
			if (cadet_likely(col < _disc.nCol - 1))
			{
				idxr.c<ResidualType>(res, col, comp) -= d_c / h2 * (stencil[1] - stencil[0]);
				// Jacobian entries
				if (wantJac)
				{
					jac[0] += static_cast<double>(d_c) / static_cast<double>(h2);
					jac[1] -= static_cast<double>(d_c) / static_cast<double>(h2);
				}
			}

			// Left side, leave out if we're in the first cell (boundary condition)
			if (cadet_likely(col > 0))
			{
				idxr.c<ResidualType>(res, col, comp) -= d_c / h2 * (stencil[-1] - stencil[0]);
				// Jacobian entries
				if (wantJac)
				{
					jac[0]  += static_cast<double>(d_c) / static_cast<double>(h2);
					jac[-1] -= static_cast<double>(d_c) / static_cast<double>(h2);
				}
			}

			// ------------------- Convection -------------------

			// Add convection through this cell's left face
			if (cadet_likely(col > 0))
			{
				// Remember that vm still contains the reconstructed value of the previous 
				// cell's *right* face, which is identical to this cell's *left* face!
				idxr.c<ResidualType>(res, col, comp) -= u / h * vm;

				// Jacobian entries
				if (wantJac)
				{
					for (int i = 0; i < 2 * wenoOrder - 1; ++i)
						// Note that we have an offset of -1 here (compared to the right cell face below), since
						// the reconstructed value depends on the previous stencil (which has now been moved by one cell)
						jac[i - wenoOrder] -= static_cast<double>(u) / static_cast<double>(h) * static_cast<double>(_wenoDerivatives[i]);
				}
			}
			else
			{
				// In the first cell we need to apply the boundary condition: inflow concentration
//				idxr.c<ResidualType>(res, col, comp) -= u / h * inlet[comp];
			}

			// Reconstruct concentration on this cell's right face
			wenoOrder = _weno.reconstruct<StateType, StencilType, wantJac>(_wenoEpsilon, col, _disc.nCol, stencil, vm, _wenoDerivatives);

			// Right side
			idxr.c<ResidualType>(res, col, comp) += u / h * vm;
			// Jacobian entries
			if (wantJac)
			{
				for (int i = 0; i < 2 * wenoOrder - 1; ++i)
					jac[i - wenoOrder + 1] += static_cast<double>(u) / static_cast<double>(h) * _wenoDerivatives[i];
			}

			// Update stencil
			stencil.advance(idxr.c<StateType>(y, col + _weno.order(), comp));
			++jac;
		}
	}

	// Film diffusion with flux into beads is added in residualFlux() function

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

	// Prepare parameters
	const ParamType radius = static_cast<ParamType>(_parRadius);
	const ParamType invBetaP = 1.0 / static_cast<ParamType>(_parPorosity) - 1.0;

	active const* const parDiff = getSectionDependentSlice(_parDiffusion, _disc.nComp, secIdx);

	// Ordering of particle surface diffusion:
	// bnd0comp0, bnd0comp1, bnd0comp2, bnd1comp0, bnd1comp1, bnd1comp2
	active const* const parSurfDiff = getSectionDependentSlice(_parSurfDiffusion, idxr.strideParBound(), secIdx);

	// Midpoint of current column cell (z coordinate) - needed in externally dependent adsorption kinetic
	const double z = 1.0 / static_cast<double>(_disc.nCol) * (0.5 + colCell);

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
		const ParamType outerAreaPerVolume = _parOuterSurfAreaPerVolume[par] / radius;
		const ParamType innerAreaPerVolume = _parInnerSurfAreaPerVolume[par] / radius;

		// Mobile phase
		for (unsigned int comp = 0; comp < _disc.nComp; ++comp, ++res, ++y, ++yDot, ++jac)
		{
			*res = 0.0;

			// Add time derivatives
			if (yDotBase)
			{
				// Ultimately, we need dc_{p,comp} / dt + 1 / beta_p * [ sum_i  dq_comp^i / dt ]
				// Compute the sum in the brackets first, then divide by beta_p and add dc_p / dt

				// Sum dq_comp^1 / dt + dq_comp^2 / dt + ... + dq_comp^{N_comp} / dt
				for (unsigned int i = 0; i < _disc.nBound[comp]; ++i)
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
				const ParamType dr = (_parCenterRadius[par - 1] - _parCenterRadius[par]) * radius;

				// Molecular diffusion contribution
				const ResidualType gradCp = (y[-idxr.strideParShell()] - y[0]) / dr;
				*res -= outerAreaPerVolume * dp * gradCp;

				// Surface diffusion contribution
				for (unsigned int i = 0; i < _disc.nBound[comp]; ++i)
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
					for (unsigned int i = 0; i < _disc.nBound[comp]; ++i)
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
				const ParamType dr = (_parCenterRadius[par] - _parCenterRadius[par + 1]) * radius;

				// Molecular diffusion contribution
				const ResidualType gradCp = (y[0] - y[idxr.strideParShell()]) / dr;
				*res += innerAreaPerVolume * dp * gradCp;

				// Surface diffusion contribution
				for (unsigned int i = 0; i < _disc.nBound[comp]; ++i)
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
					for (unsigned int i = 0; i < _disc.nBound[comp]; ++i)
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

		_binding->residual(t, z, _parCenterRadius[par], secIdx, timeFactor, y, yDot, res);
		if (wantJac)
		{
			// static_cast should be sufficient here, but this statement is also analyzed when wantJac = false
			_binding->analyticJacobian(static_cast<double>(t), z, _parCenterRadius[par], secIdx, reinterpret_cast<double const*>(y), jac);
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
	const ParamType radius = static_cast<ParamType>(_parRadius);

	active const* const filmDiff = getSectionDependentSlice(_filmDiffusion, _disc.nComp, secIdx);
	// Ordering of particle surface diffusion:
	// bnd0comp0, bnd0comp1, bnd0comp2, bnd1comp0, bnd1comp1, bnd1comp2
	active const* const parDiff = getSectionDependentSlice(_parDiffusion, _disc.nComp, secIdx);

	const ParamType surfaceToVolumeRatio = 3.0 / radius;
	const ParamType outerAreaPerVolume = _parOuterSurfAreaPerVolume[0] / radius;

	const ParamType jacCF_val = invBetaC * surfaceToVolumeRatio;
	const ParamType jacPF_val = -outerAreaPerVolume / epsP;

	// Discretized film diffusion kf for finite volumes
	ParamType* const kf_FV = _discParFlux.create<ParamType>(_disc.nComp);

	const double relOuterShellHalfRadius = 0.5 * _parCellSize[0];
	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
	{
		kf_FV[comp] = 1.0 / (radius * relOuterShellHalfRadius / epsP / static_cast<ParamType>(parDiff[comp]) + 1.0 / static_cast<ParamType>(filmDiff[comp]));
	}

	// Get offsets
	ResidualType* const resCol = resBase;
	ResidualType* const resPar = resBase + idxr.offsetCp();
	ResidualType* const resFlux = resBase + idxr.offsetJf();

	StateType const* const yCol = yBase;
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
			resPar[pblk * idxr.strideParBlock() + comp] += jacPF_val * yFlux[eq];
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
	const double radius = static_cast<double>(_parRadius);

	active const* const filmDiff = getSectionDependentSlice(_filmDiffusion, _disc.nComp, secIdx);
	// Ordering of particle diffusion:
	// sec0comp0, sec0comp1, sec0comp2, sec1comp0, sec1comp1, sec1comp2
	active const* const parDiff = getSectionDependentSlice(_parDiffusion, _disc.nComp, secIdx);

	const double surfaceToVolumeRatio = 3.0 / radius;
	const double outerAreaPerVolume   = _parOuterSurfAreaPerVolume[0] / radius;

	const double jacCF_val = invBetaC * surfaceToVolumeRatio;
	const double jacPF_val = -outerAreaPerVolume / epsP;
	const double relOuterShellHalfRadius = 0.5 * _parCellSize[0];

	// Discretized film diffusion kf for finite volumes
	double* const kf_FV = _discParFlux.create<double>(_disc.nComp);
	for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
		kf_FV[comp] = 1.0 / (radius * relOuterShellHalfRadius / epsP / static_cast<double>(parDiff[comp]) + 1.0 / static_cast<double>(filmDiff[comp]));

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
			_jacPF[pblk].addElement(comp, eq, jacPF_val);
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

void GeneralRateModel::residualSensFwdNorm(unsigned int nSens, const active& t, unsigned int secIdx,
		const active& timeFactor, double const* const y, double const* const yDot,
		const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, double* const norms,
		active* const adRes, double* const tmp)
{
	// Evaluate residual for all parameters using AD in vector mode
	residualImpl<double, active, active, false>(t, secIdx, timeFactor, y, yDot, adRes);

	for (unsigned int param = 0; param < yS.size(); param++)
	{
		// Directional derivative (dF / dy) * s
		multiplyWithJacobian(yS[param], tmp);

		// Directional derivative (dF / dyDot) * sDot
		multiplyWithDerivativeJacobian(ySdot[param], _tempState, static_cast<double>(timeFactor));

		// Complete sens residual is the sum
		norms[param] = 0.0;
		for (unsigned int i = 0; i < numDofs(); i++)
		{
			tmp[i] += _tempState[i] + adRes[i].getADValue(param);
			norms[param] = std::max(std::abs(tmp[i]), norms[param]);
		}
		LOG(Debug) << "sensRes = " << cadet::log::VectorPtr<double>(tmp, numDofs());
	}
}

int GeneralRateModel::residualSensFwdWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot,
	active* const adRes, active* const adY, unsigned int numSensAdDirs)
{
/*
	LOG(Debug) << "======= RESIDUAL SENS ========== t = " << static_cast<double>(t) << std::endl; // << " sec = " << secIdx;
	LOG(Debug) << "y = " << cadet::log::VectorPtr<double>(y, numDofs()) << "\n"
	           << "yDot = " << cadet::log::VectorPtr<double>(yDot, numDofs());
*/

	BENCH_SCOPE(_timerResidualSens);

	// Evaluate residual for all parameters using AD in vector mode and at the same time update the 
	// Jacobian (in one AD run, if analytic Jacobians are disabled)
	return residual(t, secIdx, timeFactor, y, yDot, nullptr, adRes, adY, numSensAdDirs, true, true);
}

int GeneralRateModel::residualSensFwdAdOnly(const active& t, unsigned int secIdx, const active& timeFactor,
	double const* const y, double const* const yDot, active* const adRes)
{
/*
	LOG(Debug) << "======= RESIDUAL SENS ========== t = " << static_cast<double>(t) << std::endl; // << " sec = " << secIdx;
	LOG(Debug) << "y = " << cadet::log::VectorPtr<double>(y, numDofs()) << "\n"
	           << "yDot = " << cadet::log::VectorPtr<double>(yDot, numDofs());
*/

	BENCH_SCOPE(_timerResidualSens);

	// Evaluate residual for all parameters using AD in vector mode
	return residualImpl<double, active, active, false>(t, secIdx, timeFactor, y, yDot, adRes); 
}

int GeneralRateModel::residualSensFwdCombine(const active& timeFactor, const std::vector<const double*>& yS, const std::vector<const double*>& ySdot,
	const std::vector<double*>& resS, active const* adRes, double* const tmp1, double* const tmp2, double* const tmp3)
{
/*
	LOG(Debug) << "s = " << cadet::log::VectorPtr<double>(yS[0], numDofs()) << "\n"
	           << "sDot = " << cadet::log::VectorPtr<double>(ySdot[0], numDofs());
*/

	BENCH_SCOPE(_timerResidualSens);

	// tmp1 stores result of (dF / dy) * s
	// tmp2 stores result of (dF / dyDot) * sDot

	for (unsigned int param = 0; param < yS.size(); param++)
	{
		// Directional derivative (dF / dy) * s
		multiplyWithJacobian(yS[param], tmp1);

		// Directional derivative (dF / dyDot) * sDot
		multiplyWithDerivativeJacobian(ySdot[param], tmp2, static_cast<double>(timeFactor));

		double* const ptrResS = resS[param];

		BENCH_START(_timerResidualSensPar);

		// Complete sens residual is the sum:
		#pragma omp parallel for schedule(static)
		for (ompuint_t i = 0; i < numDofs(); i++)
			ptrResS[i] = tmp1[i] + tmp2[i] + adRes[i].getADValue(param);

		BENCH_STOP(_timerResidualSensPar);

/*
		LOG(Debug) << "tmp1 = " << cadet::log::VectorPtr<double>(tmp1, numDofs()) << "\n"
		           << "tmp2 = " << cadet::log::VectorPtr<double>(tmp2, numDofs()) << "\n"
		           << "adRes = " << cadet::log::VectorPtr<active>(adRes, numDofs()) << "\n"
		           << "sensRes = " << cadet::log::VectorPtr<double>(ptrResS, numDofs());
*/
	}

	return 0;
}

int GeneralRateModel::residualSensFwd(unsigned int nSens, const active& t, unsigned int secIdx,
	const active& timeFactor, double const* const y, double const* const yDot, double const* const res,
	const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS,
	active* const adRes, double* const tmp1, double* const tmp2, double* const tmp3)
{
	BENCH_START(_timerResidualSens);
	
	residualSensFwdAdOnly(t, secIdx, timeFactor, y, yDot, adRes);
	residualSensFwdCombine(timeFactor, yS, ySdot, resS, adRes, tmp1, tmp2, tmp3);
	
	BENCH_STOP(_timerResidualSens);
	
	return 0;
}

/**
 * @brief Multiplies the given vector with the system Jacobian (i.e., @f$ \frac{\partial F}{\partial y} @f$)
 * @details Actually, the operation @f$ z = \alpha \frac{\partial F}{\partial y} x + \beta z @f$ is performed.
 * @param [in] yS Vector @f$ x @f$ that is transformed by the Jacobian @f$ \frac{\partial F}{\partial y} @f$
 * @param [in] alpha Factor @f$ \alpha @f$ in front of @f$ \frac{\partial F}{\partial y} @f$
 * @param [in] beta Factor @f$ \beta @f$ in front of @f$ z @f$
 * @param [in,out] ret Vector @f$ z @f$ which stores the result of the operation
 */
void GeneralRateModel::multiplyWithJacobian(double const* yS, double alpha, double beta, double* ret)
{
	Indexer idxr(_disc);

	// Set fluxes(ret) = fluxes(yS)
	// This applies the identity matrix in the bottom right corner of the Jaocbian (flux equation)
	for (unsigned int i = idxr.offsetJf(); i < numDofs(); ++i)
		ret[i] = alpha * yS[i] + beta * ret[i];

	BENCH_START(_timerResidualSensPar);

	#pragma omp parallel
	{
		// Threads that are done with multiplying with the bulk column blocks can proceed
		// to the particle blocks
		#pragma omp for schedule(static) nowait
		for (int comp = 0; comp < static_cast<int>(_disc.nComp); ++comp)
		{
			_jacC[comp].multiplyVector(yS + comp * idxr.strideColComp(), alpha, beta, ret + comp * idxr.strideColComp());
		}

		#pragma omp for schedule(static)
		for (ompuint_t pblk = 0; pblk < _disc.nCol; ++pblk)
		{
			const int localOffset = idxr.offsetCp(pblk);
			_jacP[pblk].multiplyVector(yS + localOffset, alpha, beta, ret + localOffset);
			_jacPF[pblk].multiplyVector(yS + idxr.offsetJf(), alpha, 1.0, ret + localOffset);
		}
	}

	BENCH_STOP(_timerResidualSensPar);

	// Multiply with the flux block in the column equation
	_jacCF.multiplyVector(yS + idxr.offsetJf(), alpha, 1.0, ret);

	// Handle flux equation
	double* const retJf = ret + idxr.offsetJf();
	_jacFC.multiplyVector(yS, alpha, 1.0, retJf);

	for (unsigned int pblk = 0; pblk < _disc.nCol; ++pblk)
		_jacFP[pblk].multiplyVector(yS + idxr.offsetCp(pblk), alpha, 1.0, retJf);
}

/**
 * @brief Multiplies the time derivative Jacobian @f$ \frac{\partial F}{\partial \dot{y}} @f$ with a given vector
 * @details The operation @f$ z = \frac{\partial F}{\partial \dot{y}} x @f$ is performed.
 *          The matrix-vector multiplication is transformed matrix-free (i.e., no matrix is explicitly formed).
 * @param [in] sDot Vector @f$ x @f$ that is transformed by the Jacobian @f$ \frac{\partial F}{\partial \dot{y}} @f$
 * @param [out] ret Vector @f$ z @f$ which stores the result of the operation
 * @param [in] timeFactor Factor which is premultiplied to the time derivatives originating from time transformation
 */
void GeneralRateModel::multiplyWithDerivativeJacobian(double const* sDot, double* ret, double timeFactor)
{
	Indexer idxr(_disc);
	const double invBetaP = (1.0 / static_cast<double>(_parPorosity) - 1.0) * timeFactor;

	BENCH_START(_timerResidualSensPar);

	#pragma omp parallel for schedule(static)
	for (int pblk = -1; pblk < static_cast<int>(_disc.nCol); ++pblk)
	{
		if (cadet_unlikely(pblk == -1))
		{
			// Column
			for (int i = 0; i < idxr.offsetCp(0); ++i)
				ret[i] = timeFactor * sDot[i];
		}
		else
		{
			// Particle
			for (unsigned int shell = 0; shell < _disc.nPar; ++shell)
			{
				double const* const localSdot = sDot + idxr.offsetCp(pblk) + shell * idxr.strideParShell();
				double* const localRet = ret + idxr.offsetCp(pblk) + shell * idxr.strideParShell();

				// Mobile phase
				for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				{
					// Add derviative with respect to dc_p / dt to Jacobian
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
	}

	BENCH_STOP(_timerResidualSensPar);

	// Handle fluxes (all algebraic)
	double* const dFdyDot = ret + idxr.offsetJf();
	std::fill(dFdyDot, dFdyDot + _disc.nCol * _disc.nComp, 0.0);
}

void GeneralRateModel::setExternalFunctions(IExternalFunction** extFuns, unsigned int size)
{
	if (_binding)
		_binding->setExternalFunctions(extFuns, size);
}

active GeneralRateModel::inletConnectionFactorActive(unsigned int compIdx, unsigned int secIdx) const CADET_NOEXCEPT
{
	return -getSectionDependentScalar(_velocity, secIdx) / _colLength * static_cast<double>(_disc.nCol);
}

double GeneralRateModel::inletConnectionFactor(unsigned int compIdx, unsigned int secIdx) const CADET_NOEXCEPT
{
	const double u = static_cast<double>(getSectionDependentScalar(_velocity, secIdx));
	const double h = static_cast<double>(_colLength) / static_cast<double>(_disc.nCol);
	return -u / h;
}

unsigned int GeneralRateModel::localOutletComponentIndex() const CADET_NOEXCEPT
{
	return _disc.nCol - 1;
}

unsigned int GeneralRateModel::localInletComponentIndex() const CADET_NOEXCEPT
{
	return 0;
}

unsigned int GeneralRateModel::localOutletComponentStride() const CADET_NOEXCEPT
{
	return _disc.nCol;
}

unsigned int GeneralRateModel::localInletComponentStride() const CADET_NOEXCEPT
{
	return _disc.nCol;
}

void GeneralRateModel::expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut)
{
	// @todo Write this function
}

/**
 * @brief Computes equidistant radial nodes in the beads
 * @details Normalized coordinates are used (i.e., outer bead shell has radius @c 1.0). The radial size
 *          of each shell is @f$ \Delta r = \frac{1}{n} @f$, where @f$ n @f$ is the number of shells.
 *          Thus, the @f$i@f$-th shell (from outside to inside of the bead starting from 0) has outer radius
 *          @f$ r_{\text{out}} = \left(1 - \Delta r i \right) @f$ and inner radius
 *          @f$ r_{\text{in}} = \left(1 - \Delta r (i+1) \right). @f$
 */
void GeneralRateModel::setEquidistantRadialDisc()
{
	const double dr = 1.0 / static_cast<double>(_disc.nPar);
	_parCellSize.assign(_disc.nPar, dr);

	for (unsigned int cell = 0; cell < _disc.nPar; cell++)
	{
		_parCenterRadius[cell] = 1.0 - (0.5 + static_cast<double>(cell)) * dr;

		// Compute denominator -> corresponding to cell volume
		const double vol = std::pow(1.0 - static_cast<double>(cell) * dr, 3.0) - std::pow(1.0 - static_cast<double>(cell + 1) * dr, 3.0);

		_parOuterSurfAreaPerVolume[cell] = 3.0 * sqr(1.0 - static_cast<double>(cell) * dr) / vol;
		_parInnerSurfAreaPerVolume[cell] = 3.0 * sqr(1.0 - static_cast<double>(cell + 1) * dr) / vol;
	}
}

/**
 * @brief Computes the radial nodes in the beads in such a way that all shells have the same volume
 * @details Normalized coordinates are used (i.e., outer bead shell has radius @c 1.0). The full
 *          normalized bead volume is @f$ 1^3 = 1.@f$ Thus, each shell should have a volume of @f$ \frac{1}{n}@f$,
 *          where @f$ n @f$ is the number of shells. Inner and outer radius, which bound a shell, have to
 *          satisfy @f[ r_{\text{out}}^3 - r_{\text{in}}^3 = \frac{1}{n}. @f]
 *          Solving for @f$ r_{\text{in}} @f$, we get
 *          @f[ r_{\text{in}} = \left( r_{\text{out}}^3 - \frac{1}{n} \right)^{\frac{1}{3}}. @f]
 */
void GeneralRateModel::setEquivolumeRadialDisc()
{
	double r_out = 1.0;
	double r_in = 0.0;

	for (unsigned int cell = 0; cell < _disc.nPar; ++cell)
	{
		if (cell != (_disc.nPar - 1))
			r_in = std::pow(std::pow(r_out, 3.0) - 1.0 / static_cast<double>(_disc.nPar), (1.0 / 3.0));

		_parCellSize[cell] = r_out - r_in;
		_parCenterRadius[cell] = r_out - 0.5 * _parCellSize[cell];

		const double invVol = static_cast<double>(_disc.nPar);
		_parOuterSurfAreaPerVolume[cell] = 3.0 * sqr(r_out) * invVol;
		_parInnerSurfAreaPerVolume[cell] = 3.0 * sqr(r_in) * invVol;

		// For the next cell: r_out == r_in of the current cell
		r_out = r_in;
	}
}

/**
 * @brief Computes all helper quantities for radial bead discretization from given radial cell boundaries
 * @details Calculates surface areas per volume for every shell and the radial shell centers.
 * @param cellInterfaces Vector with normalized radial cell boundaries (starting with @c 0.0 and ending with @c 1.0)
 */
void GeneralRateModel::setUserdefinedRadialDisc(const std::vector<double>& cellInterfaces)
{
	// Care for the right ordering and include 0.0 / 1.0 if not already in the vector.
	std::vector<double> orderedInterfaces = cellInterfaces;

	if (std::find(orderedInterfaces.begin(), orderedInterfaces.end(), 0.0) == orderedInterfaces.end())
		orderedInterfaces.push_back(0.0);
	if (std::find(orderedInterfaces.begin(), orderedInterfaces.end(), 1.0) == orderedInterfaces.end())
		orderedInterfaces.push_back(1.0);

	// Sort in descending order
	std::sort(orderedInterfaces.begin(), orderedInterfaces.end(), std::greater<double>());

	for (unsigned int cell = 0; cell < _disc.nPar; ++cell)
	{
		_parCellSize[cell] = orderedInterfaces[cell] - orderedInterfaces[cell + 1];
		_parCenterRadius[cell] = orderedInterfaces[cell] - 0.5 * _parCellSize[cell];

		// Compute denominator -> corresponding to cell volume
		const double vol = std::pow(orderedInterfaces[cell], 3.0) - std::pow(orderedInterfaces[cell + 1], 3.0);

		_parOuterSurfAreaPerVolume[cell] = 3.0 * sqr(orderedInterfaces[cell]) / vol;
		_parInnerSurfAreaPerVolume[cell] = 3.0 * sqr(orderedInterfaces[cell + 1]) / vol;
	}
}


}  // namespace model

}  // namespace cadet
