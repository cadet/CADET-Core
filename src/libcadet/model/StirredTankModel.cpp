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

#include "model/StirredTankModel.hpp"
#include "ParamReaderHelper.hpp"
#include "ParamReaderScopes.hpp"
#include "cadet/Exceptions.hpp"
#include "cadet/SolutionRecorder.hpp"
#include "model/BindingModel.hpp"
#include "model/parts/BindingConsistentInit.hpp"
#include "SimulationTypes.hpp"
#include "ParallelSupport.hpp"

#include "ConfigurationHelper.hpp"
#include "linalg/Norms.hpp"
#include "AdUtils.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <algorithm>
#include <functional>

namespace cadet
{

namespace model
{

CSTRModel::CSTRModel(UnitOpIdx unitOpIdx) : UnitOperationBase(unitOpIdx), _nComp(0), _nParType(0), _nBound(nullptr), _boundOffset(nullptr), _strideBound(nullptr), _offsetParType(nullptr), 
	_totalBound(0), _analyticJac(true), _jac(), _jacFact(), _factorizeJac(false), _initConditions(0), _initConditionsDot(0)
{
	// Mutliplexed binding models make no sense in CSTR
	_singleBinding = false;
}

CSTRModel::~CSTRModel() CADET_NOEXCEPT
{
	delete[] _boundOffset;
	delete[] _nBound;
	delete[] _strideBound;
	delete[] _offsetParType;
}

unsigned int CSTRModel::numDofs() const CADET_NOEXCEPT
{
	return 2 *_nComp + _totalBound + 1;
}

unsigned int CSTRModel::numPureDofs() const CADET_NOEXCEPT
{
	return _nComp + _totalBound + 1;
}

bool CSTRModel::usesAD() const CADET_NOEXCEPT
{
#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	// We always need AD if we want to check the analytical Jacobian
	return true;
#else
	// We only need AD if we are not computing the Jacobian analytically
	return !_analyticJac;
#endif
}

void CSTRModel::setFlowRates(active const* in, active const* out) CADET_NOEXCEPT 
{ 
	_flowRateIn = in[0];
	_flowRateOut = out[0];
}

bool CSTRModel::configureModelDiscretization(IParameterProvider& paramProvider, IConfigHelper& helper)
{
	_nComp = paramProvider.getInt("NCOMP");

	if (paramProvider.exists("NBOUND"))
	{
		const std::vector<int> nBound = paramProvider.getIntArray("NBOUND");
		_nParType = nBound.size() / _nComp;
		
		_nBound = new unsigned int[_nComp * _nParType];
		std::copy(nBound.begin(), nBound.begin() + _nComp * _nParType, _nBound);

		// Precompute offsets and total number of bound states (DOFs in solid phase)
		_boundOffset = new unsigned int[_nComp * _nParType];
		_offsetParType = new unsigned int[_nParType];
		_strideBound = new unsigned int[_nParType];
		_totalBound = 0;
		_offsetParType[0] = 0;
		for (unsigned int j = 0; j < _nParType; ++j)
		{
			unsigned int* const bo = _boundOffset + j * _nComp;
			unsigned int* const nb = _nBound + j * _nComp;

			bo[0] = 0;
			for (unsigned int i = 1; i < _nComp; ++i)
			{
				bo[i] = bo[i - 1] + nb[i - 1];
			}
			_strideBound[j] = bo[_nComp-1] + nb[_nComp - 1];
			_totalBound += _strideBound[j];

			if (j < _nParType - 1)
				_offsetParType[j + 1] = _totalBound;
		}
	}
	else
	{
		_nParType = 0;
		_totalBound = 0;
		_nBound = nullptr;
		_boundOffset = nullptr;
		_offsetParType = nullptr;
	}

	if ((_nParType > 1) && (_totalBound == 0))
		throw InvalidParameterException("Multiple particle types but not a single bound state detected");

	// Make sure each particle type has at least one bound state
	for (unsigned int i = 0; i < _nParType; ++i)
	{
		if (_strideBound[i] == 0)
			throw InvalidParameterException("Particle type " + std::to_string(i) + " does not have a bound state");
	}

	// Allocate Jacobian
	const unsigned int nVar = _nComp + _totalBound + 1;
	_jac.resize(nVar, nVar);
	_jacFact.resize(nVar, nVar);

	// Allocate space for initial conditions
	_initConditions.resize(nVar);
	_initConditionsDot.resize(nVar, 0.0);

	// Determine whether analytic Jacobian should be used
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	bool analyticJac = true;
	if (paramProvider.exists("USE_ANALYTIC_JACOBIAN"))
		analyticJac = paramProvider.getBool("USE_ANALYTIC_JACOBIAN");
#else
	// Default to AD Jacobian when analytic Jacobian is to be checked
	const bool analyticJac = false;
#endif
	useAnalyticJacobian(analyticJac);

	// ==== Construct and configure binding models
	clearBindingModels();
	_binding = std::vector<IBindingModel*>(_nParType, nullptr);

	if ((_nParType > 0) && !paramProvider.exists("ADSORPTION_MODEL"))
		throw InvalidParameterException("Binding model required when using at least one particle type");

	bool bindingConfSuccess = true;
	if (_nParType > 0)
	{
		const std::vector<std::string> bindModelNames = paramProvider.getStringArray("ADSORPTION_MODEL");
		if (bindModelNames.size() < _nParType)
			throw InvalidParameterException("Field ADSORPTION_MODEL contains too few elements (" + std::to_string(_nParType) + " required)");

		for (unsigned int i = 0; i < _nParType; ++i)
		{
			_binding[i] = helper.createBindingModel(bindModelNames[i]);
			if (!_binding[i])
				throw InvalidParameterException("Unknown binding model " + bindModelNames[i]);

			MultiplexedScopeSelector scopeGuard(paramProvider, "adsorption", _nParType == 1, i, _nParType == 1, _binding[i]->usesParamProviderInDiscretizationConfig());
			bindingConfSuccess = _binding[i]->configureModelDiscretization(paramProvider, _nComp, _nBound + i * _nComp, _boundOffset + i * _nComp) && bindingConfSuccess;
		}
	}

	return bindingConfSuccess;
}

bool CSTRModel::configure(IParameterProvider& paramProvider)
{
	_curFlowRateFilter = 0.0;
	_flowRateFilter.clear();
	const bool hasFlowrateFilter = paramProvider.exists("FLOWRATE_FILTER");
	if (hasFlowrateFilter)
		readScalarParameterOrArray(_flowRateFilter, paramProvider, "FLOWRATE_FILTER", 1);

	_porosity = 1.0;
	if (paramProvider.exists("POROSITY"))
		_porosity = paramProvider.getDouble("POROSITY");

	if (_totalBound > 0)
	{
		// Let PAR_TYPE_VOLFRAC default to 1.0 for backwards compatibility
		if (paramProvider.exists("PAR_TYPE_VOLFRAC"))
			readScalarParameterOrArray(_parTypeVolFrac, paramProvider, "PAR_TYPE_VOLFRAC", 1);
		else
			_parTypeVolFrac.resize(1, 1.0);

		if (_nParType != _parTypeVolFrac.size())
			throw InvalidParameterException("Number of elements in field PAR_TYPE_VOLFRAC does not match number of particle types");

		// Check that particle volume fractions sum to 1.0
		const double volFracSum = std::accumulate(_parTypeVolFrac.begin(), _parTypeVolFrac.end(), 0.0, 
			[](double a, const active& b) -> double { return a + static_cast<double>(b); });
		if (std::abs(1.0 - volFracSum) > 1e-10)
			throw InvalidParameterException("Sum of field PAR_TYPE_VOLFRAC differs from 1.0 (is " + std::to_string(volFracSum) + ")");
	}

	_parameters.clear();
	if (hasFlowrateFilter)
		registerScalarSectionDependentParam(hashString("FLOWRATE_FILTER"), _parameters, _flowRateFilter, _unitOpIdx, ParTypeIndep);
	_parameters[makeParamId(hashString("POROSITY"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = &_porosity;
	if (_totalBound > 0)
		registerParam1DArray(_parameters, _parTypeVolFrac, [=](bool multi, unsigned int type) { return makeParamId(hashString("PAR_TYPE_VOLFRAC"), _unitOpIdx, CompIndep, type, BoundStateIndep, ReactionIndep, SectionIndep); });

	// Register initial conditions parameters
	for (unsigned int i = 0; i < _nComp; ++i)
		_parameters[makeParamId(hashString("INIT_C"), _unitOpIdx, i, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = _initConditions.data() + i;

	for (unsigned int type = 0; type < _nParType; ++type)
	{
		if (!_binding[type])
			continue;

		std::vector<ParameterId> initParams(_strideBound[type]);
		_binding[type]->fillBoundPhaseInitialParameters(initParams.data(), _unitOpIdx, type);

		active* const ic = _initConditions.data() + _nComp + _offsetParType[type];
		for (unsigned int i = 0; i < _strideBound[type]; ++i)
			_parameters[initParams[i]] = ic + i;
	}

	_parameters[makeParamId(hashString("INIT_VOLUME"), _unitOpIdx, CompIndep, ParTypeIndep, BoundStateIndep, ReactionIndep, SectionIndep)] = _initConditions.data() + _nComp + _totalBound;

	// Reconfigure binding model
	if (!_binding.empty())
	{
		bool bindingConfSuccess = true;
		for (unsigned int type = 0; type < _nParType; ++type)
		{
 			if (!_binding[type] || !_binding[type]->requiresConfiguration())
 				continue;

			MultiplexedScopeSelector scopeGuard(paramProvider, "adsorption", type, _nParType == 1, false);
			if (scopeGuard.isActive())
				continue;
			
			bindingConfSuccess = _binding[type]->configure(paramProvider, _unitOpIdx, type) && bindingConfSuccess;
		}

		return bindingConfSuccess;
	}

	return true;
}

unsigned int CSTRModel::threadLocalMemorySize() const CADET_NOEXCEPT
{
	const unsigned int nVar = _nComp + _totalBound + 1;
	unsigned int bindingRequiredMem = 2 * nVar * sizeof(double);

	if (_nParType == 0)
		return bindingRequiredMem;

	for (unsigned int i = 0; i < _nParType; ++i)
	{
		if (_binding[i]->requiresWorkspace())
			bindingRequiredMem = std::max(bindingRequiredMem, _binding[i]->workspaceSize(_nComp, _strideBound[i], _nBound + i * _nComp));
	}

	return bindingRequiredMem;
}

void CSTRModel::setSectionTimes(double const* secTimes, bool const* secContinuity, unsigned int nSections)
{
}

void CSTRModel::useAnalyticJacobian(const bool analyticJac)
{
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	_analyticJac = analyticJac;
#else
	// Use AD Jacobian if analytic Jacobian is to be checked
	_analyticJac = false;
#endif
}

void CSTRModel::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, const AdJacobianParams& adJac)
{
	if (_flowRateFilter.size() > 1)
	{
		_curFlowRateFilter = _flowRateFilter[secIdx];
	}
	else if (_flowRateFilter.size() == 1)
	{
		_curFlowRateFilter = _flowRateFilter[0];
	}
}

void CSTRModel::reportSolution(ISolutionRecorder& recorder, double const* const solution) const
{
	Exporter expr(_nComp, _nParType, _nBound, _strideBound, _boundOffset, _totalBound, solution);
	recorder.beginUnitOperation(_unitOpIdx, *this, expr);
	recorder.endUnitOperation();
}

void CSTRModel::reportSolutionStructure(ISolutionRecorder& recorder) const
{
	Exporter expr(_nComp, _nParType, _nBound, _strideBound, _boundOffset, _totalBound, nullptr);
	recorder.unitOperationStructure(_unitOpIdx, *this, expr);
}

unsigned int CSTRModel::requiredADdirs() const CADET_NOEXCEPT
{
#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
	// We don't need AD if analytic Jaocbian is enabled
	if (_analyticJac)
		return 0;
#endif

	// We always need AD if CADET_CHECK_ANALYTIC_JACOBIAN and if analytic Jacobian is disabled
	return _nComp + _totalBound + 1;
}

void CSTRModel::prepareADvectors(const AdJacobianParams& adJac) const
{
	// Early out if AD is disabled
	if (!adJac.adY)
		return;

	ad::prepareAdVectorSeedsForDenseMatrix(adJac.adY + _nComp, adJac.adDirOffset, _jac.columns());
}

void CSTRModel::applyInitialCondition(const SimulationState& simState) const
{
	// Inlet DOFs
	std::fill(simState.vecStateY, simState.vecStateY + _nComp, 0.0);
	std::fill(simState.vecStateYdot, simState.vecStateYdot + _nComp, 0.0);

	const unsigned int nDof = numPureDofs();
	ad::copyFromAd(_initConditions.data(), simState.vecStateY + _nComp, nDof);

	std::copy(_initConditionsDot.data(), _initConditionsDot.data() + nDof, simState.vecStateYdot + _nComp);
}

void CSTRModel::readInitialCondition(IParameterProvider& paramProvider)
{
	// Clear time derivative
	std::fill(_initConditionsDot.begin(), _initConditionsDot.end(), 0.0);

	// Check if INIT_STATE is present
	if (paramProvider.exists("INIT_STATE"))
	{
		const unsigned int nDof = numPureDofs();
		const std::vector<double> initState = paramProvider.getDoubleArray("INIT_STATE");

		ad::copyToAd(initState.data(), _initConditions.data(), nDof);

		// Check if INIT_STATE contains the full state and its time derivative
		if (initState.size() >= 2 * nDof)
			std::copy(initState.data() + nDof, initState.data() + 2 * nDof, _initConditionsDot.data());

		return;
	}

	const std::vector<double> initC = paramProvider.getDoubleArray("INIT_C");

	if (initC.size() < _nComp)
		throw InvalidParameterException("INIT_C does not contain enough values for all components");
	
	ad::copyToAd(initC.data(), _initConditions.data(), _nComp);

	if (paramProvider.exists("INIT_Q"))
	{
		const std::vector<double> initQ = paramProvider.getDoubleArray("INIT_Q");
		ad::copyToAd(initQ.data(), _initConditions.data() + _nComp, _totalBound);
	}
	else
		ad::fillAd(_initConditions.data() + _nComp, _totalBound, 0.0);

	if (paramProvider.exists("INIT_VOLUME"))
		_initConditions[_nComp + _totalBound].setValue(paramProvider.getDouble("INIT_VOLUME"));
	else
		_initConditions[_nComp + _totalBound].setValue(0.0);
}

void CSTRModel::consistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, const util::ThreadLocalArray& threadLocalMem)
{
	double * const c = vecStateY + _nComp;
	const double v = c[_nComp + _totalBound];

	// Check if volume is 0
	if (v == 0.0)
	{
		const double flowIn = static_cast<double>(_flowRateIn);
		const double flowOut = static_cast<double>(_flowRateOut);

		// Volume: \dot{V} = F_{in} - F_{out} - F_{filter}
		const double vDot = flowIn - flowOut - static_cast<double>(_curFlowRateFilter);

		// We have the equation
		//    \dot{V} * (c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]) + V * (\dot{c}_i + 1 / beta * [sum_j sum_m d_j \dot{q}_{i,m}]) = c_{in,i} * F_in + c_i * F_out
		// which is now algebraic wrt. c due to V = 0:
		//    \dot{V} * (c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]) = c_{in,i} * F_in + c_i * F_out
		// Separating knowns from unknowns gives
		//    (\dot{V} + F_out) * c + \dot{V} / beta * [sum_j sum_m d_j q_{i,m}] = c_in * F_in

		// We assume that all binding models are fully dynamic (i.e., they do not have
		// algebraic equations). 
		// TODO: Figure out what to do when at least one binding model has algebraic equations.

		// If the binding models do not have algebraic equations, the
		// bound states q_{i,j} are dynamic and, thus, already determined
		//    (\dot{V} + F_out) * c = c_in * F_in - \dot{V} / beta * [sum_j sum_m d_j q_{j,i,m}]
		// This finally leads to
		//    c = (c_in * F_in - \dot{V} / beta * [sum_j sum_m d_j q_{j,i,m}]) / (\dot{V} + F_out)

		// Note that if the denominator (\dot{V} + F_out) were 0, we had
		//    0 = \dot{V} + F_out = F_{in} - F_{filter}
		// which leads to
		//    F_in = F_filter
		// Since F_out >= 0 and \dot{V} = -F_out, we get
		//    \dot{V} <= 0
		// Assuming a valid configuration, we obtain \dot{V} = 0
		// as the tank would get a negative volume otherwise.
		// Concluding, we arrive at \dot{V} = F_out = 0.
		// In this situation, F_in = F_filter = 0 has to hold
		// as otherwise the liquid (solvent) is immediately and
		// fully taken out, leaving only the pure dry components.
		// We, hence, assume that this doesn't happen and simply
		// do nothing leaving the initial conditions in place.

		const double denom = vDot + flowOut;
		if (denom == 0.0)
			return;

		for (unsigned int i = 0; i < _nComp; ++i)
			c[i] = vecStateY[i] * flowIn / denom;;

		const double qFactor = vDot * (1.0 / static_cast<double>(_porosity) - 1.0) / denom;
		for (unsigned int type = 0; type < _nParType; ++type)
		{
			unsigned int const* const bo = _boundOffset + type * _nComp;
			unsigned int const* const nb = _nBound + type * _nComp;
			double const* const typeQ = c + _nComp + _offsetParType[type];
			for (unsigned int i = 0; i < _nComp; ++i)
			{
				double qSum = 0.0;
				double const* const localQ = typeQ + bo[i];
				for (unsigned int j = 0; j < nb[i]; ++j)
					qSum += localQ[j];

				c[i] -= qFactor * static_cast<double>(_parTypeVolFrac[type]) * qSum;
			}
		}
	}
	else
	{
		double* const buffer = threadLocalMem.get();
		for (unsigned int type = 0; type < _nParType; ++type)
		{
			// Compute quasi-stationary binding model state
			if (_binding[type]->hasAlgebraicEquations())
			{
				const unsigned int offset = _nComp + _offsetParType[type];
				active* const localAdRes = adJac.adRes ? adJac.adRes + offset : nullptr;
				active* const localAdY = adJac.adY ? adJac.adY + offset : nullptr;

				linalg::DenseMatrixView jacobianMatrix(_jacFact.data(), _jacFact.pivotData(), _strideBound[type], _strideBound[type]);

				// Solve algebraic variables
				_binding[type]->consistentInitialState(simTime.t, simTime.secIdx, ColumnPosition{0.0, 0.0, 0.0}, c + offset, c, errorTol, localAdRes, localAdY,
					offset, adJac.adDirOffset, ad::DenseJacobianExtractor(), buffer, jacobianMatrix);
			}
		}
	}
}

void CSTRModel::consistentInitialTimeDerivative(const SimulationTime& simTime, double const* vecStateY, double* const vecStateYdot, const util::ThreadLocalArray& threadLocalMem) 
{
	double const* const c = vecStateY + _nComp;
	double* const cDot = vecStateYdot + _nComp;
	const double v = c[_nComp + _totalBound];

	const double flowIn = static_cast<double>(_flowRateIn);
	const double flowOut = static_cast<double>(_flowRateOut);

	// Note that the residual has not been negated, yet. We will do that now.
	for (unsigned int i = 0; i < numDofs(); ++i)
		vecStateYdot[i] = -vecStateYdot[i];

	// Assemble time derivative Jacobian
	_jacFact.setAll(0.0);
	addTimeDerivativeJacobian(simTime.t, simTime.timeFactor, ConstSimulationState{vecStateY, nullptr}, _jacFact);

	// Check if volume is 0
	if (v == 0.0)
	{
		// Volume: \dot{V} = F_{in} - F_{out} - F_{filter}
		const double vDot = flowIn - flowOut - static_cast<double>(_curFlowRateFilter);
		_jacFact.native(_nComp + _totalBound, _nComp + _totalBound) = 1.0;
		cDot[_nComp + _totalBound] = vDot;

		// We have the equation
		//    V * \dot{c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]} + \dot{V} * (c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]) = c_{in,i} * F_in - c_i * F_out
		// which is now algebraic wrt. c due to V = 0:
		//    \dot{V} * (c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]) = c_{in,i} * F_in - c_i * F_out
		// So we take the derivative wrt. to time t on both sides
		//    2 * \dot{V} * \dot{c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]} + V * \ddot{c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]} + \ddot{V} * (c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]) = \dot{c}_{in,i} * F_in - \dot{c} * F_out
		// and use the fact that \ddot{V} = 0 and V = 0 to arrive at
		//    2 * \dot{V} * \dot{c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]} = \dot{c}_{in,i} * F_in - \dot{c} * F_out
		// Separating knowns from unknowns gives
		//    (2 * \dot{V} + F_out) * \dot{c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]} = \dot{c}_{in,i} * F_in
		// which finally yields
		//    \dot{c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]} = \dot{c}_{in,i} * F_in / (2 * \dot{V} + F_out)

		// Note that if the denominator were 0, we had
		//    0 = 2 * \dot{V} + F_out = 2 * F_in - 2 * F_filter - F_out
		// which leads to
		//    F_out = 2 * F_in - 2 * F_filter                       (*)
		// Plugging this back into the \dot{V} equation gives
		//    \dot{V} = F_in - F_out - F_filter = -F_in + F_filter
		// Since V = 0, a valid choice of parameters has to ensure
		// \dot{V} >= 0. In this case, we obtain
		//    F_in <= F_filter
		// On the other hand, we infer from (*) that
		//    0 <= F_out <= 0   =>   F_out = 0   =>   F_in = F_filter
		// This, in turn, concludes \dot{V} = 0. In this situation, 
		//    F_in = F_filter = 0
		// has to hold as otherwise the liquid (solvent) is immediately
		// and fully taken out, leaving only the pure dry components.
		// We, hence, assume that this doesn't happen. Summarizing, we
		// have
		//    \dot{V} = 0, F_in = F_filter = F_out = 0
		// and nothing can happen or change. Therefore, \dot{c} is set
		// to 0.0.

		typename linalg::DenseMatrix::RowIterator itRow = _jacFact.row(0);

		const double denom = 2.0 * vDot + flowOut;
		if (denom == 0.0)
		{
			// Assume F_in = F_filter = 0, set cDot to 0
			for (unsigned int i = 0; i < _nComp; ++i, ++itRow)
			{
				itRow.setAll(0.0);
				itRow[0] = 1.0;
			}
			std::fill(cDot, cDot + _nComp, 0.0);
		}
		else
		{
			const double factor = flowIn / denom;
			const double invBeta = 1.0 / static_cast<double>(_porosity) - 1.0;
			for (unsigned int i = 0; i < _nComp; ++i, ++itRow)
			{
				itRow.setAll(0.0);

				// Jacobian of c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]

				// d/dc_i
				itRow[0] = 1.0;

				// Advance to bound states
				const int bndIdx = _nComp - i;
				for (unsigned int type = 0; type < _nParType; ++type)
				{
					const int innerIdx = bndIdx + _offsetParType[type] + _boundOffset[type * _nComp + i];
					const double innerFactor = static_cast<double>(_parTypeVolFrac[type]) * invBeta;
					for (unsigned int j = 0; j < _nBound[type * _nComp + i]; ++j)
					{
						itRow[innerIdx + j] = innerFactor;
					}
				}

				// TODO: This is wrong as vecStateYdot does not contain \dot{c}_in (on entry)
				// This scenario violates the assumption that every outlet DOF is dynamic
				// which is key to the consistent initialization algorithm. Fixing this problem
				// requires fundamental changes to the consistent initialization concept
				// implemented so far.
				vecStateYdot[i] = 0.0;
				cDot[i] = vecStateYdot[i] * factor;
			}
		}
	}

	double* const buffer = threadLocalMem.get();

	// Overwrite rows corresponding to algebraic equations with the Jacobian and set right hand side to 0
	for (unsigned int type = 0; type < _nParType; ++type)
	{
		if (_binding[type]->hasAlgebraicEquations())
		{
			parts::BindingConsistentInitializer::consistentInitialTimeDerivative(_binding[type], simTime, _jacFact.row(_nComp + _offsetParType[type]),
				_jac.row(_nComp + _offsetParType[type]), vecStateYdot + 2 * _nComp + _offsetParType[type], ColumnPosition{0.0, 0.0, 0.0}, buffer);
		}
	}

	// Factorize
	const bool result = _jacFact.robustFactorize(buffer);
	if (!result)
	{
		LOG(Error) << "Factorize() failed";
	}

	// Solve
	const bool result2 = _jacFact.robustSolve(vecStateYdot + _nComp, buffer);
	if (!result2)
	{
		LOG(Error) << "Solve() failed";
	}
}

void CSTRModel::leanConsistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, const util::ThreadLocalArray& threadLocalMem)
{
	// It is assumed that the bound states are (approximately) correctly initialized.
	// Thus, only the liquid phase has to be initialized

	double * const c = vecStateY + _nComp;
	const double v = c[_nComp];

	// Check if volume is 0
	if (v == 0.0)
	{
		const double flowIn = static_cast<double>(_flowRateIn);
		const double flowOut = static_cast<double>(_flowRateOut);

		// Volume: \dot{V} = F_{in} - F_{out} - F_{filter}
		const double vDot = flowIn - flowOut - static_cast<double>(_curFlowRateFilter);

		// We have the equation
		//    \dot{V} * (c + 1 / beta * [sum_j q_j]) + V * (\dot{c} + 1 / beta * [sum_j \dot{q}_j]) = c_in * F_in + c * F_out
		// which is now algebraic wrt. c due to V = 0:
		//    \dot{V} * (c + 1 / beta * [sum_j q_j]) = c_in * F_in + c * F_out
		// Separating knowns from unknowns gives
		//    (\dot{V} + F_out) * c = c_in * F_in - \dot{V} / beta * [sum_j q_j]
		// Hence, we obtain
		//    c = (c_in * F_in - \dot{V} / beta * [sum_j q_j]) / (\dot{V} + F_out)

		// Note that if the denominator were 0, we had
		//    0 = \dot{V} + F_out = F_in - F_filter
		// which leads to
		//    F_in = F_filter
		// Since F_out >= 0 and \dot{V} = -F_out, we get
		//    \dot{V} <= 0
		// Assuming a valid configuration, we obtain \dot{V} = 0
		// as the tank would get a negative volume otherwise.
		// Concluding, we arrive at \dot{V} = F_out = 0.
		// In this situation, F_in = F_filter = 0 has to hold
		// as otherwise the liquid (solvent) is immediately and
		// fully taken out, leaving only the pure dry components.
		// We, hence, assume that this doesn't happen and simply
		// do nothing leaving the initial conditions in place.

		const double denom = vDot + flowOut;
		if (denom == 0.0)
			return;

		for (unsigned int i = 0; i < _nComp; ++i)
			c[i] = vecStateY[i] * flowIn / denom;;

		const double qFactor = vDot * (1.0 / static_cast<double>(_porosity) - 1.0) / denom;
		for (unsigned int type = 0; type < _nParType; ++type)
		{
			unsigned int const* const bo = _boundOffset + type * _nComp;
			unsigned int const* const nb = _nBound + type * _nComp;
			double const* const typeQ = c + _nComp + _offsetParType[type];
			for (unsigned int i = 0; i < _nComp; ++i)
			{
				double qSum = 0.0;
				double const* const localQ = typeQ + bo[i];
				for (unsigned int j = 0; j < nb[i]; ++j)
					qSum += localQ[j];

				c[i] -= qFactor * static_cast<double>(_parTypeVolFrac[type]) * qSum;
			}
		}
	}
}

void CSTRModel::leanConsistentInitialTimeDerivative(double t, double timeFactor, double const* const vecStateY, double* const vecStateYdot, double* const res, const util::ThreadLocalArray& threadLocalMem)
{
	// It is assumed that the bound states are (approximately) correctly initialized.
	// Thus, only the liquid phase has to be initialized

	double const* const c = vecStateY + _nComp;
	double* const cDot = vecStateYdot + _nComp;
	double* const resC = res + _nComp;
	const double v = c[_nComp];

	const double flowIn = static_cast<double>(_flowRateIn);
	const double flowOut = static_cast<double>(_flowRateOut);

	// Volume: \dot{V} = F_{in} - F_{out} - F_{filter}
	const double vDot = flowIn - flowOut - static_cast<double>(_curFlowRateFilter);
	cDot[_nComp] = vDot;

	// Check if volume is 0
	if (v == 0.0)
	{
		// We have the equation
		//    V * \dot{c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]} + \dot{V} * (c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]) = c_{in,i} * F_in - c_i * F_out
		// which is now algebraic wrt. c due to V = 0:
		//    \dot{V} * (c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]) = c_{in,i} * F_in - c_i * F_out
		// So we take the derivative wrt. to time t on both sides
		//    2 * \dot{V} * \dot{c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]} + V * \ddot{c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]} + \ddot{V} * (c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]) = \dot{c}_{in,i} * F_in - \dot{c} * F_out
		// and use the fact that \ddot{V} = 0 and V = 0 to arrive at
		//    2 * \dot{V} * \dot{c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]} = \dot{c}_{in,i} * F_in - \dot{c} * F_out
		// Separating knowns from unknowns gives
		//    (2 * \dot{V} + F_out) * \dot{c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]} = \dot{c}_{in,i} * F_in
		// which finally yields
		//    \dot{c_i + 1 / beta * [sum_j sum_m d_j q_{i,m}]} = \dot{c}_{in,i} * F_in / (2 * \dot{V} + F_out)

		// Note that if the denominator were 0, we had
		//    0 = 2 * \dot{V} + F_out = 2 * F_in - 2 * F_filter - F_out
		// which leads to
		//    F_out = 2 * F_in - 2 * F_filter                       (*)
		// Plugging this back into the \dot{V} equation gives
		//    \dot{V} = F_in - F_out - F_filter = -F_in + F_filter
		// Since V = 0, a valid choice of parameters has to ensure
		// \dot{V} >= 0. In this case, we obtain
		//    F_in <= F_filter
		// On the other hand, we infer from (*) that
		//    0 <= F_out <= 0   =>   F_out = 0   =>   F_in = F_filter
		// This, in turn, concludes \dot{V} = 0. In this situation, 
		//    F_in = F_filter = 0
		// has to hold as otherwise the liquid (solvent) is immediately
		// and fully taken out, leaving only the pure dry components.
		// We, hence, assume that this doesn't happen. Summarizing, we
		// have
		//    \dot{V} = 0, F_in = F_filter = F_out = 0
		// and nothing can happen or change. Therefore, \dot{c} is set
		// to 0.0.

		const double denom = 2.0 * vDot + flowOut;
		if (denom == 0.0)
		{
			// Assume F_in = F_filter = 0
			std::fill(cDot, cDot + _nComp, 0.0);
		}
		else
		{
			const double qFactor = 2.0 * vDot * (1.0 / static_cast<double>(_porosity) - 1.0) / denom;
			const double factor = flowIn / denom;

			for (unsigned int i = 0; i < _nComp; ++i)
			{
				// TODO: This is wrong as vecStateYdot does not contain \dot{c}_in (on entry)
				// This scenario violates the assumption that every outlet DOF is dynamic
				// which is key to the consistent initialization algorithm. Fixing this problem
				// requires fundamental changes to the consistent initialization concept
				// implemented so far.
				vecStateYdot[i] = 0.0;

				cDot[i] = vecStateYdot[i] * factor;

				for (unsigned int type = 0; type < _nParType; ++type)
				{
					unsigned int const* const bo = _boundOffset + type * _nComp;
					unsigned int const* const nb = _nBound + type * _nComp;
					double const* const localQdot = cDot + _nComp + _offsetParType[type] + bo[i];

					double qDotSum = 0.0;
					for (unsigned int j = 0; j < nb[i]; ++j)
						qDotSum += localQdot[j];

					cDot[i] -= qFactor * static_cast<double>(_parTypeVolFrac[type]) * qDotSum;
				}
			}
		}
	}
	else
	{
		const double invBeta = 1.0 / static_cast<double>(_porosity) - 1.0;

		// Concentrations: V * (\dot{c} + 1 / beta * [sum_j sum_m d_j \dot{q}_{j,m}]) = c_in * F_in + c * F_out - \dot{V} * (c + 1 / beta * [sum_j sum_m d_j q_{j,m}])
		//             <=> V * \dot{c} = c_in * F_in + c * F_out - \dot{V} * (c + 1 / beta * [sum_j sum_m d_j q_{j,m}]) - V / beta * [sum_j sum_m d_j \dot{q}_{j,m}]
		//                             = -res - \dot{V} * (c + 1 / beta * [sum_j sum_m d_j q_{j,m}]) - V / beta * [sum_j sum_m d_j \dot{q}_{j,m}]
		// => \dot{c} = (-res - \dot{V} * (c + 1 / beta * [sum_j sum_m d_j q_{j,m}]) - V / beta * [sum_j sum_m d_j \dot{q}_m]) / V
		for (unsigned int i = 0; i < _nComp; ++i)
		{
			double qSum = 0.0;
			double qDotSum = 0.0;
			for (unsigned int type = 0; type < _nParType; ++type)
			{
				double const* const localQ = c + _nComp + _offsetParType[type] + _boundOffset[type * _nComp + i];
				double const* const localQdot = cDot + _nComp +  _offsetParType[type] +_boundOffset[type * _nComp + i];
				double qSumType = 0.0;
				double qDotSumType = 0.0;
				for (unsigned int j = 0; j < _nBound[type * _nComp + i]; ++j)
				{
					qSumType += localQ[j];
					qDotSumType += localQdot[j];
				}
				qSum += static_cast<double>(_parTypeVolFrac[type]) * qSumType;
				qDotSum += static_cast<double>(_parTypeVolFrac[type]) * qDotSumType;
			}

			cDot[i] = (-resC[i] - vDot * (c[i] + invBeta * qSum) - v * invBeta * qDotSum) / v;
		}
	}
}

void CSTRModel::leanConsistentInitialSensitivity(const ActiveSimulationTime& simTime, const ConstSimulationState& simState,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes, const util::ThreadLocalArray& threadLocalMem)
{
	consistentInitialSensitivity(simTime, simState, vecSensY, vecSensYdot, adRes, threadLocalMem);
}

int CSTRModel::residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const util::ThreadLocalArray& threadLocalMem)
{
	return residualImpl<double, double, double, false>(simTime.t, simTime.secIdx, simTime.timeFactor, simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
}

template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
int CSTRModel::residualImpl(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* const y, double const* const yDot, ResidualType* const res, const util::ThreadLocalArray& threadLocalMem)
{
	StateType const* const cIn = y;
	StateType const* const c = y + _nComp;
	const StateType& v = y[2 * _nComp + _totalBound];

	double const* const cDot = yDot ? yDot + _nComp : nullptr;
	const double vDot = yDot ? yDot[2 * _nComp + _totalBound] : 0.0;

	const ParamType flowIn = static_cast<ParamType>(_flowRateIn);
	const ParamType flowOut = static_cast<ParamType>(_flowRateOut);

	// Inlet DOF
	for (unsigned int i = 0; i < _nComp; ++i)
	{
		res[i] = cIn[i];
	}

	// Concentrations: \dot{V} * (c + 1 / beta * [sum_j q_j]) + V * (\dot{c} + 1 / beta * [sum_j \dot{q}_j]) = c_in * F_in - c * F_out
	const ParamType invBeta = 1.0 / static_cast<ParamType>(_porosity) - 1.0;
	ResidualType* const resC = res + _nComp;
	for (unsigned int i = 0; i < _nComp; ++i)
	{
		resC[i] = 0.0;

		// Add time derivatives
		if (cadet_likely(yDot))
		{
			// Ultimately, we need (dc_{i} / dt + 1 / beta * [ sum_j sum_m d_j dq_{j,i,m} / dt ]) * V
			// and (c_{i} + 1 / beta * [ sum_j sum_m d_j q_{i,m} ]) * dV / dt
			// Compute the sum in the brackets first, then divide by beta and add remaining term

			// Sum q_{i,1} + q_{i,2} + ... + q_{i,N_i}
			// and sum dq_{i,1} / dt + dq_{i,2} / dt + ... + dq_{i,N_i} / dt
			StateType qSum = 0.0;
			double qDotSum = 0.0;

			for (unsigned int type = 0; type < _nParType; ++type)
			{
				StateType const* const q = c + _nComp + _offsetParType[type] + _boundOffset[type * _nComp + i];
				double const* const qDot = cDot + _nComp + _offsetParType[type] + _boundOffset[type * _nComp + i];
				
				StateType qSumType = 0.0;
				double qDotSumType = 0.0;
				for (unsigned int j = 0; j < _nBound[type * _nComp + i]; ++j)
				{
					qSumType += q[j];
					qDotSumType += qDot[j];
				}

				const double pvf = static_cast<double>(_parTypeVolFrac[type]);
				qSum += pvf * qSumType;
				qDotSum += pvf * qDotSumType;
			}

			// Divide by beta and add c_i and dc_i / dt
			resC[i] = timeFactor * ((cDot[i] + invBeta * qDotSum) * v + vDot * (c[i] + invBeta * qSum));
		}

		resC[i] += -flowIn * cIn[i] + flowOut * c[i];
	}

	// Bound states
	double* const buffer = threadLocalMem.get();
	for (unsigned int type = 0; type < _nParType; ++type)
	{
		double const* const qDot = yDot ? yDot + 2 * _nComp + _offsetParType[type] : nullptr;
		_binding[type]->residual(t, secIdx, timeFactor, ColumnPosition{0.0, 0.0, 0.0}, c + _nComp + _offsetParType[type], c, qDot, res + 2 * _nComp + _offsetParType[type], buffer);
	}

	// Volume: \dot{V} = F_{in} - F_{out} - F_{filter}
	res[2 * _nComp + _totalBound] = vDot - flowIn + flowOut + static_cast<ParamType>(_curFlowRateFilter);

	if (wantJac)
	{
		_jac.setAll(0.0);

		// Assemble Jacobian: dRes / dy

		// Concentrations: \dot{V} * (c_i + 1 / beta * [sum_j sum_m d_j q_{j,i,m}]) + V * (\dot{c}_i + 1 / beta * [sum_j sum_m d_j \dot{q}_{j,i,m}]) - c_{in,i} * F_in + c_i * F_out == 0
		const double vDotTimeFactor = static_cast<double>(vDot) * static_cast<double>(timeFactor);
		for (unsigned int i = 0; i < _nComp; ++i)
		{
			_jac.native(i, i) = vDotTimeFactor + static_cast<double>(flowOut);

			if (cadet_likely(yDot))
			{
				double qDotSum = 0.0;
				const double vDotInvBeta = vDotTimeFactor * static_cast<double>(invBeta);
				for (unsigned int type = 0; type < _nParType; ++type)
				{
					double const* const qiDot = cDot + _nComp + _offsetParType[type] + _boundOffset[type * _nComp + i];
					const unsigned int localOffset = _nComp + _offsetParType[type] + _boundOffset[type * _nComp + i];

					const double vDotInvBetaParVolFrac = vDotInvBeta * static_cast<double>(_parTypeVolFrac[type]);
					double qDotSumType = 0.0;
					for (unsigned int j = 0; j < _nBound[type * _nComp + i]; ++j)
					{
						_jac.native(i, localOffset + j) = vDotInvBetaParVolFrac;
						// + _nComp: Moves over liquid phase components
						// + _offsetParType[type]: Moves to particle type
						// + _boundOffset[i]: Moves over bound states of previous components
						// + j: Moves to current bound state j of component i

						qDotSumType += qiDot[j];
					}

					qDotSum += static_cast<double>(_parTypeVolFrac[type]) * qDotSumType;
				}

				_jac.native(i, _nComp + _totalBound) = (cDot[i] + static_cast<double>(invBeta) * qDotSum) * static_cast<double>(timeFactor);
			}
		}

		// Bound states
		for (unsigned int type = 0; type < _nParType; ++type)
			_binding[0]->analyticJacobian(static_cast<double>(t), secIdx, ColumnPosition{0.0, 0.0, 0.0}, reinterpret_cast<double const*>(y) + 2 * _nComp + _offsetParType[type], _nComp + _offsetParType[type], _jac.row(_nComp + _offsetParType[type]), buffer);

		// Volume: \dot{V} - F_{in} + F_{out} + F_{filter} == 0
	}

	return 0;
}

int CSTRModel::residual(const ActiveSimulationTime& simTime, const ConstSimulationState& simState, double* const res,
	const AdJacobianParams& adJac, const util::ThreadLocalArray& threadLocalMem, bool updateJacobian, bool paramSensitivity)
{
	if (updateJacobian)
	{
		_factorizeJac = true;

#ifndef CADET_CHECK_ANALYTIC_JACOBIAN
		if (_analyticJac)
		{
			if (paramSensitivity)
			{
				const int retCode = residualImpl<double, active, active, true>(simTime.t, simTime.secIdx, simTime.timeFactor, simState.vecStateY, simState.vecStateYdot, adJac.adRes, threadLocalMem);

				// Copy AD residuals to original residuals vector
				if (res)
					ad::copyFromAd(adJac.adRes, res, numDofs());

				return retCode;
			}
			else
				return residualImpl<double, double, double, true>(static_cast<double>(simTime.t), simTime.secIdx, static_cast<double>(simTime.timeFactor), simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
		}
		else
		{
			// Compute Jacobian via AD

			// Copy over state vector to AD state vector (without changing directional values to keep seed vectors)
			// and initalize residuals with zero (also resetting directional values)
			ad::copyToAd(simState.vecStateY, adJac.adY, numDofs());
			// @todo Check if this is necessary
			ad::resetAd(adJac.adRes, numDofs());

			// Evaluate with AD enabled
			int retCode = 0;
			if (paramSensitivity)
				retCode = residualImpl<active, active, active, false>(simTime.t, simTime.secIdx, simTime.timeFactor, adJac.adY, simState.vecStateYdot, adJac.adRes, threadLocalMem);
			else
				retCode = residualImpl<active, active, double, false>(static_cast<double>(simTime.t), simTime.secIdx, static_cast<double>(simTime.timeFactor), adJac.adY, simState.vecStateYdot, adJac.adRes, threadLocalMem);

			// Copy AD residuals to original residuals vector
			if (res)
				ad::copyFromAd(adJac.adRes, res, numDofs());

			// Extract Jacobian
			extractJacobianFromAD(adJac.adRes, adJac.adDirOffset);

			return retCode;
		}
#else
		// Compute Jacobian via AD

		// Copy over state vector to AD state vector (without changing directional values to keep seed vectors)
		// and initalize residuals with zero (also resetting directional values)
		ad::copyToAd(simState.vecStateY, adJac.adY, numDofs());
		// @todo Check if this is necessary
		ad::resetAd(adJac.adRes, numDofs());

		// Evaluate with AD enabled
		int retCode = 0;
		if (paramSensitivity)
			retCode = residualImpl<active, active, active, false>(simTime.t, simTime.secIdx, simTime.timeFactor, adJac.adY, simState.vecStateYdot, adJac.adRes, threadLocalMem);
		else
			retCode = residualImpl<active, active, double, false>(static_cast<double>(simTime.t), simTime.secIdx, static_cast<double>(simTime.timeFactor), adJac.adY, simState.vecStateYdot, adJac.adRes, threadLocalMem);

		// Only do comparison if we have a residuals vector (which is not always the case)
		if (res)
		{
			// Evaluate with analytical Jacobian which is stored in the band matrices
			retCode = residualImpl<double, double, double, true>(static_cast<double>(simTime.t), simTime.secIdx, static_cast<double>(simTime.timeFactor), simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);

			// Compare AD with anaytic Jacobian
			checkAnalyticJacobianAgainstAd(adJac.adRes, adJac.adDirOffset);
		}

		// Extract Jacobian
		extractJacobianFromAD(adJac.adRes, adJac.adDirOffset);

		return retCode;
#endif
	}
	else
	{
		if (paramSensitivity)
		{
			// Initalize residuals with zero
			// @todo Check if this is necessary
			ad::resetAd(adJac.adRes, numDofs());

			const int retCode = residualImpl<double, active, active, false>(simTime.t, simTime.secIdx, simTime.timeFactor, simState.vecStateY, simState.vecStateYdot, adJac.adRes, threadLocalMem);

			// Copy AD residuals to original residuals vector
			if (res)
				ad::copyFromAd(adJac.adRes, res, numDofs());

			return retCode;
		}
		else
			return residualImpl<double, double, double, false>(static_cast<double>(simTime.t), simTime.secIdx, static_cast<double>(simTime.timeFactor), simState.vecStateY, simState.vecStateYdot, res, threadLocalMem);
	}
}

int CSTRModel::residualWithJacobian(const ActiveSimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, const util::ThreadLocalArray& threadLocalMem)
{
	return residual(simTime, simState, res, adJac, threadLocalMem, true, false);
}

int CSTRModel::residualSensFwdAdOnly(const ActiveSimulationTime& simTime, const ConstSimulationState& simState, active* const adRes, const util::ThreadLocalArray& threadLocalMem)
{
	// Evaluate residual for all parameters using AD in vector mode
	return residualImpl<double, active, active, false>(simTime.t, simTime.secIdx, simTime.timeFactor, simState.vecStateY, simState.vecStateYdot, adRes, threadLocalMem);
}

int CSTRModel::residualSensFwdWithJacobian(const ActiveSimulationTime& simTime, const ConstSimulationState& simState, const AdJacobianParams& adJac, const cadet::util::ThreadLocalArray& threadLocalMem)
{
	// Evaluate residual for all parameters using AD in vector mode and at the same time update the 
	// Jacobian (in one AD run, if analytic Jacobians are disabled)
	return residual(simTime, simState, nullptr, adJac, threadLocalMem, true, true);
}

void CSTRModel::initializeSensitivityStates(const std::vector<double*>& vecSensY) const
{
	const unsigned int nDof = numPureDofs();
	for (unsigned int param = 0; param < vecSensY.size(); ++param)
	{
		double* const sensY = vecSensY[param];
		ad::copyFromAdDirection(_initConditions.data(), sensY + _nComp, nDof, param);
	}
}

void CSTRModel::consistentInitialSensitivity(const ActiveSimulationTime& simTime, const ConstSimulationState& simState,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes, const util::ThreadLocalArray& threadLocalMem)
{
	// TODO: Handle v = 0

	double* const buffer = threadLocalMem.get();
	for (unsigned int param = 0; param < vecSensY.size(); ++param)
	{
		double* const sensY = vecSensY[param];
		double* const sensYdot = vecSensYdot[param];

		// Copy parameter derivative dF / dp from AD and negate it
		for (unsigned int i = _nComp; i < numDofs(); ++i)
			sensYdot[i] = -adRes[i].getADValue(param);

		// Step 1: Solve algebraic equations

		// Step 1a: Compute quasi-stationary binding model state
		for (unsigned int type = 0; type < _nParType; ++type)
		{
			if (_binding[type]->hasAlgebraicEquations())
			{
				// Get algebraic block
				unsigned int algStart = 0;
				unsigned int algLen = 0;
				_binding[type]->getAlgebraicBlock(algStart, algLen);

				// Reuse memory for dense matrix
				linalg::DenseMatrixView jacobianMatrix(_jacFact.data(), _jacFact.pivotData(), algLen, algLen);

				// We used to apply robustFactorize() and robustSolve() here, but the model part will only use factorize() and solve()
				parts::BindingConsistentInitializer::consistentInitialSensitivityState(algStart, algLen, jacobianMatrix,
					_jac, _nComp, _nComp + _offsetParType[type], _nComp + _offsetParType[type], _strideBound[type], sensY, sensYdot, buffer);
			}
		}

		// Step 2: Compute the correct time derivative of the state vector

		// Step 2a: Assemble, factorize, and solve diagonal blocks of linear system

		// Compute right hand side by adding -dF / dy * s = -J * s to -dF / dp which is already stored in sensYdot
		multiplyWithJacobian(toSimple(simTime), simState, sensY, -1.0, 1.0, sensYdot);

		// Note that we have correctly negated the right hand side

		// Assemble dF / dYdot
		_jacFact.setAll(0.0);
		addTimeDerivativeJacobian(static_cast<double>(simTime.t), static_cast<double>(simTime.timeFactor), simState, _jacFact);

		// Overwrite rows corresponding to algebraic equations with the Jacobian and set right hand side to 0
		for (unsigned int type = 0; type < _nParType; ++type)
		{
			if (_binding[type]->hasAlgebraicEquations())
			{
				parts::BindingConsistentInitializer::consistentInitialSensitivityTimeDerivative(_binding[type], _jacFact.row(_nComp + _offsetParType[type]),
					_jac.row(_nComp + _offsetParType[type]), sensYdot + 2 * _nComp + _offsetParType[type]);
			}
		}

		// Factorize
		const bool result = _jacFact.robustFactorize(buffer);
		if (!result)
		{
			LOG(Error) << "Factorize() failed";
		}

		// Solve
		const bool result2 = _jacFact.robustSolve(sensYdot + _nComp, buffer);
		if (!result2)
		{
			LOG(Error) << "Solve() failed";
		}
	}
}

void CSTRModel::multiplyWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* yS, double alpha, double beta, double* ret)
{
	const double flowIn = static_cast<double>(_flowRateIn);
	double* const resTank = ret + _nComp;

	// Inlet DOFs
	for (unsigned int i = 0; i < _nComp; ++i)
	{
		ret[i] = alpha * yS[i] + beta * ret[i];
	}

	// Multiply with main body Jacobian: dRes / dy
	_jac.multiplyVector(yS + _nComp, alpha, beta, resTank);

	// Map inlet DOFs to the tank (tank cells)
	for (unsigned int i = 0; i < _nComp; ++i)
	{
		resTank[i] -= alpha * flowIn * yS[i];
	}
}

void CSTRModel::multiplyWithDerivativeJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* sDot, double* ret)
{
	// Handle inlet DOFs (all algebraic)
	std::fill_n(ret, _nComp, 0.0);

	// Handle actual ODE DOFs
	double const* const c = simState.vecStateY + _nComp;
	double const* const q = simState.vecStateY + 2 * _nComp;
	const double v = simState.vecStateY[2 * _nComp + _totalBound];
	const double invBeta = 1.0 / static_cast<double>(_porosity) - 1.0;
	const double vInvBeta = simTime.timeFactor * v * invBeta;
	const double timeV = simTime.timeFactor * v;
	double* const r = ret + _nComp;
	double const* const s = sDot + _nComp;

	// Concentrations: \dot{V} * (c_i + 1 / beta * [sum_j sum_m d_j q_{j,i,m}]) + V * (\dot{c}_i + 1 / beta * [sum_j sum_m d_j \dot{q}_{j,i,m}]) - c_{in,i} * F_in + c_i * F_out == 0
	for (unsigned int i = 0; i < _nComp; ++i)
	{
		r[i] = timeV * s[i];

		double qSum = 0.0;
		for (unsigned int type = 0; type < _nParType; ++type)
		{
			double const* const qi = q + _offsetParType[type] + _boundOffset[type * _nComp + i];
			const unsigned int localOffset = _nComp + _offsetParType[type] + _boundOffset[type * _nComp + i];
			const double vInvBetaParVolFrac = vInvBeta * static_cast<double>(_parTypeVolFrac[type]);
			double qSumType = 0.0;
			for (unsigned int j = 0; j < _nBound[type * _nComp + i]; ++j)
			{
				r[i] += vInvBetaParVolFrac * s[localOffset + j];
				// + _nComp: Moves over liquid phase components
				// + _offsetParType[type]: Moves to particle type
				// + _boundOffset[i]: Moves over bound states of previous components
				// + j: Moves to current bound state j of component i

				qSumType += qi[j];
			}

			qSum += static_cast<double>(_parTypeVolFrac[type]) * qSumType;
		}
		r[i] += simTime.timeFactor * (c[i] + invBeta * qSum) * s[_nComp + _totalBound];
	}

	// Bound states
	for (unsigned int type = 0; type < _nParType; ++type)
		_binding[type]->multiplyWithDerivativeJacobian(s + _nComp + _offsetParType[type], r + _nComp + _offsetParType[type], simTime.timeFactor);

	// Volume: \dot{V} - F_{in} + F_{out} + F_{filter} == 0
	r[_nComp + _totalBound] = simTime.timeFactor * s[_nComp + _totalBound];
}

int CSTRModel::linearSolve(double t, double timeFactor, double alpha, double tol, double* const rhs, double const* const weight,
	const ConstSimulationState& simState)
{
	const double flowIn = static_cast<double>(_flowRateIn);

	// Handle inlet equations by backsubstitution
	for (unsigned int i = 0; i < _nComp; ++i)
	{
		rhs[i + _nComp] += flowIn * rhs[i];
	}

	bool success = true;
	if (_factorizeJac)
	{
		// Factorization is necessary
		_factorizeJac = false;
		_jacFact.copyFrom(_jac);

		addTimeDerivativeJacobian(t, timeFactor * alpha, simState, _jacFact);
		success = _jacFact.factorize();
	}
	success = success && _jacFact.solve(rhs + _nComp);

	// Return 0 on success and 1 on failure
	return success ? 0 : 1;
}

template <typename MatrixType>
void CSTRModel::addTimeDerivativeJacobian(double t, double timeFactor, const ConstSimulationState& simState, MatrixType& mat)
{
	double const* const c = simState.vecStateY + _nComp;
	double const* const q = simState.vecStateY + 2 * _nComp;
	const double v = simState.vecStateY[2 * _nComp + _totalBound];
	const double invBeta = 1.0 / static_cast<double>(_porosity) - 1.0;
	const double vInvBeta = timeFactor * v * invBeta;
	const double timeV = timeFactor * v;

	// Assemble Jacobian: dRes / dyDot

	// Concentrations: \dot{V} * (c_i + 1 / beta * [sum_j sum_m d_j q_{j,i,m}]) + V * (\dot{c}_i + 1 / beta * [sum_j sum_m d_j \dot{q}_{j,i,m}]) - c_{in,i} * F_in + c_i * F_out == 0
	for (unsigned int i = 0; i < _nComp; ++i)
	{
		mat.native(i, i) += timeV;

		double qSum = 0.0;
		for (unsigned int type = 0; type < _nParType; ++type)
		{
			double const* const qi = q + _offsetParType[type] + _boundOffset[type * _nComp + i];
			const unsigned int localOffset = _nComp + _offsetParType[type] + _boundOffset[type * _nComp + i];
			const double vInvBetaParVolFrac = vInvBeta * static_cast<double>(_parTypeVolFrac[type]);
			double qSumType = 0.0;
			for (unsigned int j = 0; j < _nBound[type * _nComp + i]; ++j)
			{
				mat.native(i, localOffset + j) += vInvBetaParVolFrac;
				// + _nComp: Moves over liquid phase components
				// + _offsetParType[type]: Moves to particle type
				// + _boundOffset[i]: Moves over bound states of previous components
				// + j: Moves to current bound state j of component i

				qSumType += qi[j];
			}

			qSum += static_cast<double>(_parTypeVolFrac[type]) * qSumType;
		}
		mat.native(i, _nComp + _totalBound) += timeFactor * (c[i] + invBeta * qSum);
	}

	// Bound states
	for (unsigned int type = 0; type < _nParType; ++type)
		_binding[type]->jacobianAddDiscretized(timeFactor, mat.row(_nComp + _offsetParType[type]));

	// Volume: \dot{V} - F_{in} + F_{out} + F_{filter} == 0
	mat.native(_nComp + _totalBound, _nComp + _totalBound) += timeFactor;
}

/**
 * @brief Extracts the system Jacobian from AD seed vectors
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
void CSTRModel::extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset)
{
	ad::extractDenseJacobianFromAd(adRes + _nComp, adDirOffset, _jac);
}

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN

/**
 * @brief Compares the analytical Jacobian with a Jacobian derived by AD
 * @details The analytical Jacobian is assumed to be stored in the dense matrix.
 * @param [in] adRes Residual vector of AD datatypes with band compressed seed vectors
 * @param [in] adDirOffset Number of AD directions used for non-Jacobian purposes (e.g., parameter sensitivities)
 */
void CSTRModel::checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const
{
	const double diff = ad::compareDenseJacobianWithAd(adRes + _nComp, adDirOffset, _jac);
	LOG(Debug) << "AD dir offset: " << adDirOffset << " diff: " << diff;
}

#endif

void registerCSTRModel(std::unordered_map<std::string, std::function<IUnitOperation*(UnitOpIdx)>>& models)
{
	models[CSTRModel::identifier()] = [](UnitOpIdx uoId) { return new CSTRModel(uoId); };
}

}  // namespace model

}  // namespace cadet
