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

#include "model/binding/BindingModelBase.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "nonlin/Solver.hpp"
#include "ParamReaderHelper.hpp"

#include "AdUtils.hpp"
#include "linalg/Norms.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <iterator>
#include <algorithm>

namespace cadet
{

namespace model
{

BindingModelBase::BindingModelBase() : _nComp(0), _nBoundStates(nullptr), _nonlinearSolver(nullptr) { }
BindingModelBase::~BindingModelBase() CADET_NOEXCEPT
{
	delete _nonlinearSolver;
}

void BindingModelBase::configureModelDiscretization(unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
{
	_nComp = nComp;
	_nBoundStates = nBound;
	if (hasMultipleBoundStates(nBound, nComp) && !supportsMultistate())
		throw InvalidParameterException("Binding model does not support multiple bound states");
}

bool BindingModelBase::configure(IParameterProvider& paramProvider, unsigned int unitOpIdx)
{
	// Read binding dynamics (quasi-stationary, kinetic)
	_kineticBinding = paramProvider.getInt("IS_KINETIC");

	if (!configureNonlinearSolver(paramProvider))
		return false;

	return configureImpl(false, paramProvider, unitOpIdx);
}

bool BindingModelBase::reconfigure(IParameterProvider& paramProvider, unsigned int unitOpIdx)
{
	// Read binding dynamics (quasi-stationary, kinetic)
	_kineticBinding = paramProvider.getInt("IS_KINETIC");

	// Clear all parameters and reconfigure
	_parameters.clear();
	return configureImpl(true, paramProvider, unitOpIdx);
}

bool BindingModelBase::configureNonlinearSolver(IParameterProvider& paramProvider)
{
	if (paramProvider.exists("consistency_solver"))
	{
		paramProvider.pushScope("consistency_solver");

		const std::string name = paramProvider.getString("SOLVER_NAME");
		delete _nonlinearSolver;
		_nonlinearSolver = nonlin::createSolver(name);
		_nonlinearSolver->configure(paramProvider);

		paramProvider.popScope();
	}
	else if (!_nonlinearSolver)
	{
		// Use default solver with default settings
		_nonlinearSolver = nonlin::createSolver("");
	}
	return true;	
}

std::unordered_map<ParameterId, double> BindingModelBase::getAllParameterValues() const
{
	std::unordered_map<ParameterId, double> data;
	std::transform(_parameters.begin(), _parameters.end(), std::inserter(data, data.end()),
	               [](const std::pair<const ParameterId, active*>& p) { return std::make_pair(p.first, static_cast<double>(*p.second)); });
	return data;
}

bool BindingModelBase::hasParameter(const ParameterId& pId) const
{
	return _parameters.find(pId) != _parameters.end();
}

bool BindingModelBase::setParameter(const ParameterId& pId, int value)
{
	return false;
}

bool BindingModelBase::setParameter(const ParameterId& pId, double value)
{
	auto paramHandle = _parameters.find(pId);
	if (paramHandle != _parameters.end())
	{
		paramHandle->second->setValue(value);
		return true;
	}

	return false;
}

bool BindingModelBase::setParameter(const ParameterId& pId, bool value)
{
	return false;
}

active* BindingModelBase::getParameter(const ParameterId& pId)
{
	auto paramHandle = _parameters.find(pId);
	if (paramHandle != _parameters.end())
	{
		return paramHandle->second;
	}

	return nullptr;
}

unsigned int BindingModelBase::consistentInitializationWorkspaceSize() const
{
	// Determine problem size
	const unsigned int eqSize = numBoundStates(_nBoundStates, _nComp);
	// Ask nonlinear solver how much memory it needs for this kind of problem
	return _nonlinearSolver->workspaceSize(eqSize);
}

/*
void BindingModelBase::timeDerivativeAlgebraicResidual(double t, double z, double r, unsigned int secIdx, double* const y, double* dResDt) const
{
	if (!hasAlgebraicEquations())
		return;

	unsigned int start = 0;
	unsigned int len = 0;
	getAlgebraicBlock(start, len);

	// Assumes no external dependence
	std::fill_n(dResDt, len, 0.0);
}
*/


PureBindingModelBase::PureBindingModelBase() { }
PureBindingModelBase::~PureBindingModelBase() CADET_NOEXCEPT { }

void PureBindingModelBase::getAlgebraicBlock(unsigned int& idxStart, unsigned int& len) const
{
	idxStart = 0;
	if (_kineticBinding)
		len = 0;
	else
		len = numBoundStates(_nBoundStates, _nComp);
}

void PureBindingModelBase::jacobianAddDiscretized(double alpha, linalg::FactorizableBandMatrix::RowIterator jac) const
{
	// We only add time derivatives for kinetic binding
	if (!_kineticBinding)
		return;

	// All equations are kinetic
	const unsigned int eqSize = numBoundStates(_nBoundStates, _nComp);
	for (unsigned int i = 0; i < eqSize; ++i, ++jac)
	{
		jac[0] += alpha;
	}
}

void PureBindingModelBase::multiplyWithDerivativeJacobian(double const* yDotS, double* const res, double timeFactor) const
{
	// Multiplier is 0 if quasi-stationary and 1 if kinetic binding
	// However, due to premultiplication of time derivatives with the timeFactor (because of time transformation),
	// we set it to timeFactor instead of 1.0 in order to save some multiplications
	const double multiplier = _kineticBinding ? timeFactor : 0.0;

	const unsigned int eqSize = numBoundStates(_nBoundStates, _nComp);
	for (unsigned int i = 0; i < eqSize; ++i)
	{
		res[i] = multiplier * yDotS[i];
	}
}

void PureBindingModelBase::consistentInitialState(double t, double z, double r, unsigned int secIdx, double* const vecStateY, double errorTol, 
	active* const adRes, active* const adY, unsigned int adEqOffset, unsigned int adOffset, unsigned int diagDir, 
	unsigned int lowerBandwidth, unsigned int upperBandwidth, double* const workingMemory, linalg::detail::DenseMatrixBase& workingMat) const
{
	// If we have kinetic binding, there are no algebraic equations and we are done
	if (_kineticBinding)
		return;

	// All equations are algebraic and (except for salt equation) nonlinear
	// Compute the q_i from their corresponding c_{p,i}

	// Determine problem size
	const unsigned int eqSize = numBoundStates(_nBoundStates, _nComp);

	// Check if workingMat satisfies size requirements
	cadet_assert(workingMat.rows() >= eqSize);
	cadet_assert(workingMat.columns() >= eqSize);

	std::fill(workingMemory, workingMemory + _nonlinearSolver->workspaceSize(eqSize), 0.0);

	// Check if workingMat satisfies size requirements
	cadet_assert(workingMat.rows() >= eqSize);
	cadet_assert(workingMat.columns() >= eqSize);

	// Select between analytic and AD Jacobian
	std::function<bool(double const* const, linalg::detail::DenseMatrixBase& jac)> jacobianFunc;
	if (adRes && adY)
	{
		// AD Jacobian
		jacobianFunc = [&](double const* const x, linalg::detail::DenseMatrixBase& mat) -> bool { 
			// Copy over state vector to AD state vector (without changing directional values to keep seed vectors)
			// and initalize residuals with zero (also resetting directional values)
			ad::copyToAd(x, adY + adEqOffset, eqSize);
			// @todo Check if this is necessary
			ad::resetAd(adRes + adEqOffset, eqSize);

			// Call residual with AD enabled
			residualCore(t, z, r, secIdx, 1.0, adY + adEqOffset, vecStateY - _nComp, nullptr, adRes + adEqOffset);
			
#ifdef CADET_CHECK_ANALYTIC_JACOBIAN			
			// Compute analytic Jacobian
			mat.setAll(0.0);
			analyticJacobianCore(t, z, r, secIdx, x, vecStateY - _nComp, mat.row(0));

			// Compare
			const double diff = ad::compareDenseJacobianWithBandedAd(adRes, adEqOffset, adOffset, diagDir, lowerBandwidth, upperBandwidth, mat);
			LOG(Debug) << "MaxDiff " << adEqOffset << ": " << diff;
#endif
			// Extract Jacobian
			ad::extractDenseJacobianFromBandedAd(adRes, adEqOffset, adOffset, diagDir, lowerBandwidth, upperBandwidth, mat);
			return true;
		};
	}
	else
	{
		// Analytic Jacobian
		jacobianFunc = [&](double const* const x, linalg::detail::DenseMatrixBase& mat) -> bool
		{ 
			mat.setAll(0.0);
			analyticJacobianCore(t, z, r, secIdx, x, vecStateY - _nComp, mat.row(0));
			return true;
		};
	}

	const bool conv = _nonlinearSolver->solve([&](double const* const x, double* const res) -> bool
		{
			residualCore(t, z, r, secIdx, 1.0, x, vecStateY - _nComp, nullptr, res); 
			return true; 
		}, 
		jacobianFunc,
		errorTol, vecStateY, workingMemory, workingMat, eqSize);
}

void PureBindingModelBase::analyticJacobian(double t, double z, double r, unsigned int secIdx, double const* y, linalg::BandMatrix::RowIterator jac) const
{
	analyticJacobianCore(t, z, r, secIdx, y, y - _nComp, jac);
}


}  // namespace model

}  // namespace cadet
