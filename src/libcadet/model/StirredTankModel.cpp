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

#include "model/StirredTankModel.hpp"
#include "ParamReaderHelper.hpp"
#include "cadet/Exceptions.hpp"
#include "cadet/SolutionRecorder.hpp"

#include "ConfigurationHelper.hpp"
#include "ParamIdUtil.hpp"
#include "linalg/Norms.hpp"

#include <algorithm>
#include <functional>

#include "LoggingUtils.hpp"
#include "Logging.hpp"
#include "AdUtils.hpp"
#include "SensParamUtil.hpp"

namespace cadet
{

namespace model
{

CSTRModel::CSTRModel(UnitOpIdx unitOpIdx) : _unitOpIdx(unitOpIdx), _jac()
{
}

CSTRModel::~CSTRModel() CADET_NOEXCEPT
{
}

unsigned int CSTRModel::numDofs() const CADET_NOEXCEPT
{
	return 2 *_nComp + 1;
}

unsigned int CSTRModel::numPureDofs() const CADET_NOEXCEPT
{
	return _nComp + 1;
}

bool CSTRModel::usesAD() const CADET_NOEXCEPT
{
	return false;
}

void CSTRModel::setFlowRates(const active& in, const active& out) CADET_NOEXCEPT 
{ 
	_flowRateIn = in;
	_flowRateOut = out;
}

bool CSTRModel::configure(IParameterProvider& paramProvider, IConfigHelper& helper)
{
	_nComp = paramProvider.getInt("NCOMP");

	const unsigned int nVar = _nComp + 1;
	_jac.resize(nVar, nVar);

	return reconfigure(paramProvider);
}

bool CSTRModel::reconfigure(IParameterProvider& paramProvider)
{
	_curFlowRateFilter = 0.0;
	_flowRateFilter.clear();
	if (paramProvider.exists("FLOWRATE_FILTER"))
	{
		readScalarParameterOrArray(_flowRateFilter, paramProvider, "FLOWRATE_FILTER", 1);
	}

	_parameters.clear();
	registerScalarSectionDependentParam(hashString("FLOWRATE_FILTER"), _parameters, _flowRateFilter, _unitOpIdx);
	return true;
}

void CSTRModel::setSectionTimes(double const* secTimes, bool const* secContinuity, unsigned int nSections)
{
}

std::unordered_map<ParameterId, double> CSTRModel::getAllParameterValues() const
{
	std::unordered_map<ParameterId, double> data;
	std::transform(_parameters.begin(), _parameters.end(), std::inserter(data, data.end()),
	               [](const std::pair<const ParameterId, active*>& p) { return std::make_pair(p.first, static_cast<double>(*p.second)); });
	return data;
}

bool CSTRModel::hasParameter(const ParameterId& pId) const
{
	return _parameters.find(pId) != _parameters.end();
}

bool CSTRModel::setParameter(const ParameterId& pId, int value)
{
	return false;
}

bool CSTRModel::setParameter(const ParameterId& pId, double value)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return false;

	auto paramHandle = _parameters.find(pId);
	if (paramHandle != _parameters.end())
	{
		paramHandle->second->setValue(value);
		return true;
	}

	return false;
}

bool CSTRModel::setParameter(const ParameterId& pId, bool value)
{
	return false;
}

void CSTRModel::setSensitiveParameterValue(const ParameterId& pId, double value)
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
}

bool CSTRModel::setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue)
{
	if ((pId.unitOperation != _unitOpIdx) && (pId.unitOperation != UnitOpIndep))
		return false;

	// Check own parameters
	auto paramHandle = _parameters.find(pId);
	if (paramHandle != _parameters.end())
	{
		LOG(Debug) << "Found parameter " << pId << " in CSTR: Dir " << adDirection << " is set to " << adValue;

		// Register parameter and set AD seed / direction
		_sensParams.insert(paramHandle->second);
		paramHandle->second->setADValue(adDirection, adValue);
		return true;
	}

	return false;
}

void CSTRModel::clearSensParams()
{
	// Remove AD directions from parameters
	for (auto sp : _sensParams)
		sp->setADValue(0.0);

	_sensParams.clear();
}

void CSTRModel::useAnalyticJacobian(const bool analyticJac) { }
void CSTRModel::notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active* const adRes, active* const adY, unsigned int adDirOffset)
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
	Exporter expr(_nComp, solution);
	recorder.beginUnitOperation(_unitOpIdx, *this, expr);
	recorder.endUnitOperation();
}

void CSTRModel::reportSolutionStructure(ISolutionRecorder& recorder) const
{
	Exporter expr(_nComp, nullptr);
	recorder.unitOperationStructure(_unitOpIdx, *this, expr);
}

unsigned int CSTRModel::requiredADdirs() const CADET_NOEXCEPT
{
	return _nComp + 1;
}

void CSTRModel::prepareADvectors(active* const adRes, active* const adY, unsigned int numSensAdDirs) const
{
	//TODO: Assemble seed vectors
}

void CSTRModel::applyInitialCondition(double* const vecStateY, double* const vecStateYdot)
{
	std::fill(vecStateY, vecStateY + numDofs(), 0.0);
	std::fill(vecStateYdot, vecStateYdot + numDofs(), 0.0);
}

void CSTRModel::applyInitialCondition(IParameterProvider& paramProvider, double* const vecStateY, double* const vecStateYdot)
{
	// Check if INIT_STATE is present
	if (paramProvider.exists("INIT_STATE"))
	{
		const std::vector<double> initState = paramProvider.getDoubleArray("INIT_STATE");
		std::copy(initState.data(), initState.data() + numDofs(), vecStateY);

		// Check if INIT_STATE contains the full state and its time derivative
		if (initState.size() >= 2 * numDofs())
		{
			double const* const srcYdot = initState.data() + numDofs();
			std::copy(srcYdot, srcYdot + numDofs(), vecStateYdot);
		}
		return;
	}

	const std::vector<double> initC = paramProvider.getDoubleArray("INIT_C");

	if (initC.size() < _nComp)
		throw InvalidParameterException("INIT_C does not contain enough values for all components");
	
	std::copy_n(initC.begin(), _nComp, vecStateY + _nComp);

	if (paramProvider.exists("INIT_VOLUME"))
		vecStateY[2 * _nComp] = paramProvider.getDouble("INIT_VOLUME");
	else
		vecStateY[2 * _nComp] = 0.0;
}


void CSTRModel::consistentInitialState(double t, unsigned int secIdx, double timeFactor, double* const vecStateY, active* const adRes, active* const adY, unsigned int numSensAdDirs, double errorTol)
{
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
		//    V * \dot{c} + \dot{V} * c = c_in * F_in - c * F_out
		// which is now algebraic wrt. c due to V = 0:
		//    \dot{V} * c = c_in * F_in - c * F_out
		// Separating knowns from unknowns gives
		//    (\dot{V} + F_out) * c = c_in * F_in
		// Hence, we obtain
		//    c = c_in * F_in / (\dot{V} + F_out)

		// Note that if the denominator were 0, we had
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
		if (denom != 0.0)
		{
			const double factor = flowIn / denom;
			for (unsigned int i = 0; i < _nComp; i++)
			{
				c[i] = vecStateY[i] * factor;
			}
		}
	}
}

void CSTRModel::consistentInitialTimeDerivative(double t, unsigned int secIdx, double timeFactor, double const* vecStateY, double* const vecStateYdot) 
{
	double const* const c = vecStateY + _nComp;
	double* const cDot = vecStateYdot + _nComp;
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
		//    V * \dot{c} + \dot{V} * c = c_in * F_in - c * F_out
		// which is now algebraic wrt. c due to V = 0:
		//    \dot{V} * c = c_in * F_in - c * F_out
		// So we take the derivative wrt. to time t on both sides
		//    2 * \dot{V} * \dot{c} + V * \ddot{c} + \ddot{V} * c = \dot{c}_in * F_in - \dot{c} * F_out
		// and use the fact that \ddot{V} = 0 and V = 0 to arrive at
		//    2 * \dot{V} * \dot{c} = \dot{c}_in * F_in - \dot{c} * F_out
		// Separating knowns from unknowns gives
		//    (2 * \dot{V} + F_out) * \dot{c} = \dot{c}_in * F_in
		// which finally yields
		//    \dot{c} = \dot{c}_in * F_in / (2 * \dot{V} + F_out)

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
			const double factor = flowIn / denom;
			for (unsigned int i = 0; i < _nComp; i++)
			{
				// TODO: This is wrong as vecStateYdot does not contain \dot{c}_in (on entry)
				vecStateYdot[i] = 0.0;
				cDot[i] = vecStateYdot[i] * factor;
			}
		}
	}
	else
	{
		// Concentrations: V * \dot{c} = c_in * F_in - c * F_out - \dot{V} * c
		//                             = -vecStateYdot - \dot{V} * c
		// => \dot{c} = (-vecStateYdot - \dot{V} * c) / V
		for (unsigned int i = 0; i < _nComp; i++)
		{
			cDot[i] = (-cDot[i] - vDot * c[i]) / v;
		}
	}


}

void CSTRModel::leanConsistentInitialState(double t, unsigned int secIdx, double timeFactor, double* const vecStateY, active* const adRes, active* const adY, unsigned int numSensAdDirs, double errorTol)
{
	consistentInitialState(t, secIdx, timeFactor, vecStateY, adRes, adY, numSensAdDirs, errorTol);
}

void CSTRModel::leanConsistentInitialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes)
{
	consistentInitialSensitivity(t, secIdx, timeFactor, vecStateY, vecStateYdot, vecSensY, vecSensYdot, adRes);
}

int CSTRModel::residual(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double* const res)
{
	return residualImpl<double, double, double>(t, secIdx, timeFactor, y, yDot, res);
}

template <typename StateType, typename ResidualType, typename ParamType>
int CSTRModel::residualImpl(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* const y, double const* const yDot, ResidualType* const res)
{
	StateType const* const cIn = y;
	StateType const* const c = y + _nComp;
	const StateType& v = y[2 * _nComp];

	double const* const cDot = yDot + _nComp;
	const double vDot = yDot ? yDot[2 * _nComp] : 0.0;

	const ParamType flowIn = static_cast<ParamType>(_flowRateIn);
	const ParamType flowOut = static_cast<ParamType>(_flowRateOut);

	// Inlet DOF
	for (unsigned int i = 0; i < _nComp; ++i)
	{
		res[i] = cIn[i];
	}

	// Concentrations: \dot{V} * c + V * \dot{c} = c_in * F_in - c * F_out
	for (unsigned int i = 0; i < _nComp; ++i)
	{
		if (yDot)
			res[i + _nComp] = v * cDot[i] + vDot * c[i] - flowIn * cIn[i] + flowOut * c[i];
		else
			res[i + _nComp] = - flowIn * cIn[i] + flowOut * c[i];
	}

	// Volume: \dot{V} = F_{in} - F_{out} - F_{filter}
	res[2 * _nComp] = vDot - flowIn + flowOut + static_cast<ParamType>(_curFlowRateFilter);

	return 0;
}

int CSTRModel::residual(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, double* const res,
	active* const adRes, active* const adY, unsigned int numSensAdDirs, bool paramSensitivity)
{
	//This method has the same inteface as the GRM but there is no Jacobian to generate or update
	if (paramSensitivity)
	{
		const int retCode = residualImpl<double, active, active>(t, secIdx, timeFactor, y, yDot, adRes);

		// Copy AD residuals to original residuals vector
		if (res)
		{
			ad::copyFromAd(adRes, res, numDofs());
		}

		return retCode;
	}
	else
		return residualImpl<double, double, double>(static_cast<double>(t), secIdx, static_cast<double>(timeFactor), y, yDot, res);
}

int CSTRModel::residualWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, double* const res, active* const adRes, active* const adY, unsigned int numSensAdDirs)
{
	return residual(t, secIdx, timeFactor, y, yDot, res, adRes, adY, numSensAdDirs, false);
}

int CSTRModel::residualSensFwdAdOnly(const active& t, unsigned int secIdx, const active& timeFactor,
	double const* const y, double const* const yDot, active* const adRes)
{
	// Evaluate residual for all parameters using AD in vector mode
	return residualImpl<double, active, active>(t, secIdx, timeFactor, y, yDot, adRes);
}

int CSTRModel::residualSensFwdCombine(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, 
	const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS, active const* adRes, 
	double* const tmp1, double* const tmp2, double* const tmp3)
{
	for (unsigned int param = 0; param < yS.size(); param++)
	{
		// Directional derivative (dF / dy) * s
		multiplyWithJacobian(static_cast<double>(t), secIdx, static_cast<double>(timeFactor), y, yDot, yS[param], 1.0, 0.0, tmp1);

		// Directional derivative (dF / dyDot) * sDot
		multiplyWithDerivativeJacobian(static_cast<double>(t), secIdx, static_cast<double>(timeFactor), y, yDot, ySdot[param], tmp2);

		// Complete sens residual is the sum:
		double* const ptrResS = resS[param];
		for (unsigned int i = 0; i < numDofs(); ++i)
		{
			ptrResS[i] = tmp1[i] + tmp2[i] + adRes[i].getADValue(param);
		}
	}
	return 0;
}

int CSTRModel::residualSensFwdWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y,
	double const* const yDot, active* const adRes, active* const adY, unsigned int numSensAdDirs)
{
	// Evaluate residual for all parameters using AD in vector mode and at the same time update the 
	// Jacobian (in one AD run, if analytic Jacobians are disabled)
	return residual(t, secIdx, timeFactor, y, yDot, nullptr, adRes, adY, numSensAdDirs, true);
}

void CSTRModel::consistentInitialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
	std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes)
{
	for (unsigned int param = 0; param < vecSensY.size(); ++param)
	{
		double* const sensY = vecSensY[param];
		double* const sensYdot = vecSensYdot[param];

		// Calculate -(dF / dY) * s - (dF / dP)
		multiplyWithJacobian(static_cast<double>(t), secIdx, static_cast<double>(timeFactor), vecStateY, vecStateYdot, sensY, -1.0, 0.0, sensYdot);

		for (unsigned int i = _nComp; i < numDofs(); ++i)
			sensYdot[i] -= adRes[i].getADValue(param);

		// Solve for \dot{s}
		assembleTimeDerivativeJacobian(static_cast<double>(t), secIdx, static_cast<double>(timeFactor), vecStateY, vecStateYdot, _jac);
		_jac.factorize();
		_jac.solve(sensYdot + _nComp);
	}
}

void CSTRModel::multiplyWithJacobian(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double const* yS, double alpha, double beta, double* ret)
{
	const double flowIn = static_cast<double>(_flowRateIn);
	double* const resTank = ret + _nComp;

	// Inlet DOFs
	for (unsigned int i = 0; i < _nComp; ++i)
	{
		ret[i] = alpha * yS[i] + beta * ret[i];
	}

	// Assemble Jacobian: dRes / dy
	assembleJacobian(t, secIdx, 0.0, y, yDot, _jac);
	_jac.multiplyVector(yS + _nComp, alpha, beta, resTank);

	// Map inlet DOFs to the tank (tank cells)
	for (unsigned int i = 0; i < _nComp; ++i)
	{
		resTank[i] -= alpha * flowIn * yS[i];
	}
}

void CSTRModel::multiplyWithDerivativeJacobian(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double const* sDot, double* ret)
{
	assembleTimeDerivativeJacobian(t, secIdx, timeFactor, y, yDot, _jac);

	// Handle inlet DOFs (all algebraic)
	std::fill_n(ret, _nComp, 0.0);
	// Multiply main body
	_jac.multiplyVector(sDot + _nComp, ret + _nComp);
}

int CSTRModel::linearSolve(double t, double timeFactor, double alpha, double tol, double* const rhs, double const* const weight,
	double const* const y, double const* const yDot, double const* const res)
{
	const double flowIn = static_cast<double>(_flowRateIn);
	const double flowOut = static_cast<double>(_flowRateOut);

	// Handle inlet equations by backsubstitution
	for (unsigned int i = 0; i < _nComp; i++)
	{
		rhs[i + _nComp] += flowIn * rhs[i];
	}

	// Assemble Jacobian: dRes / dy + alpha * dRes / dyDot
	assembleJacobian(t, 0u, timeFactor * alpha, y, yDot, _jac);
	const bool success = _jac.factorize() && _jac.solve(rhs + _nComp);

	// Return 0 on success and 1 on failure
	return success ? 0 : 1;
}

template <typename MatrixType>
void CSTRModel::assembleJacobian(double t, unsigned int secIdx, double timeFactor, double const* y, double const* yDot, MatrixType& mat)
{
	mat.setAll(0.0);

	const double flowIn = static_cast<double>(_flowRateIn);
	const double flowOut = static_cast<double>(_flowRateOut);

	double const* const c = y + _nComp;
	double const* const cDot = yDot + _nComp;
	const double v = y[2 * _nComp];
	const double vDot = yDot[2 * _nComp];

	// Assemble Jacobian: dRes / dy + alpha * dRes / dyDot

	// Concentrations: \dot{V} * c + V * \dot{c} - c_in * F_in + c * F_out == 0
	for (unsigned int i = 0; i < _nComp; i++)
	{
		mat.native(i, i) = vDot + timeFactor * v + flowOut;
		mat.native(i, _nComp) = cDot[i] + timeFactor * c[i];
	}

	// Volume: \dot{V} - F_{in} + F_{out} + F_{filter} == 0
	mat.native(_nComp, _nComp) = timeFactor;
}

template <typename MatrixType>
void CSTRModel::assembleTimeDerivativeJacobian(double t, unsigned int secIdx, double timeFactor, double const* y, double const* yDot, MatrixType& mat)
{
	mat.setAll(0.0);

	double const* const c = y + _nComp;
	const double v = y[2 * _nComp];

	// Assemble Jacobian: dRes / dyDot

	// Concentrations: \dot{V} * c + V * \dot{c} - c_in * F_in + c * F_out == 0
	for (unsigned int i = 0; i < _nComp; i++)
	{
		mat.native(i, i) = timeFactor * v;
		mat.native(i, _nComp) = timeFactor * c[i];
	}

	// Volume: \dot{V} - F_{in} + F_{out} + F_{filter} == 0
	mat.native(_nComp, _nComp) = timeFactor;
}

}  // namespace model

}  // namespace cadet
