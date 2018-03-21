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

#include "model/binding/BindingModelBase.hpp"
#include "model/ExternalFunctionSupport.hpp"
#include "model/binding/BindingModelMacros.hpp"
#include "model/binding/RefConcentrationSupport.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "nonlin/Solver.hpp"
#include "model/Parameters.hpp"
#include "SlicedVector.hpp"
#include "LocalVector.hpp"

#include "AdUtils.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <vector>
#include <unordered_map>
#include <functional>

namespace cadet
{

namespace model
{

/*<codegen>
{
	"name": "BiSMAParamHandler",
	"externalName": "ExtBiSMAParamHandler",
	"parameters":
		[
			{ "type": "ScalarBoundStateDependentParameter", "varName": "lambda", "confName": "BISMA_LAMBDA"},
			{ "type": "ComponentBoundStateMajorDependentParameter", "varName": "kA", "confName": "BISMA_KA"},
			{ "type": "ComponentBoundStateMajorDependentParameter", "varName": "kD", "confName": "BISMA_KD"},
			{ "type": "ComponentBoundStateMajorDependentParameter", "varName": "nu", "confName": "BISMA_NU"},
			{ "type": "ComponentBoundStateMajorDependentParameter", "varName": "sigma", "confName": "BISMA_SIGMA"}
		],
	"constantParameters":
		[
			{ "type": "ReferenceConcentrationBoundStateDependentParameter", "varName": ["refC0", "refQ"], "objName": "_refConcentration", "confPrefix": "BISMA_"}
		]
}
</codegen>*/

/* Parameter description
 ------------------------
 lambda = Ionic capacity in binding site-major ordering
 kA = Adsorption rate in binding site-major ordering
 kD = Desorption rate in binding site-major ordering
 nu = Characteristic charge in binding site-major ordering
 sigma = Steric factor in binding site-major ordering
 refC0, refQ = Reference concentrations
*/

inline const char* BiSMAParamHandler::identifier() CADET_NOEXCEPT { return "BI_STERIC_MASS_ACTION"; }

inline bool BiSMAParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	const unsigned int numStates = firstNonEmptyBoundStates(nBoundStates, nComp);

	if ((_kA.size() != _kD.size()) || (_kA.size() != _nu.size()) || (_kA.size() != _sigma.size()) || (_kA.size() < nComp * numStates))
		throw InvalidParameterException("BISMA_KA, BISMA_KD, BISMA_NU, and BISMA_SIGMA have to have the same size");
	if (_lambda.size() != numStates)
		throw InvalidParameterException("BISMA_LAMBDA has to have as many elements as there are binding sites");

	return true;
}

inline const char* ExtBiSMAParamHandler::identifier() CADET_NOEXCEPT { return "EXT_BI_STERIC_MASS_ACTION"; }

inline bool ExtBiSMAParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	const unsigned int numStates = firstNonEmptyBoundStates(nBoundStates, nComp);

	if ((_kA.size() != _kD.size()) || (_kA.size() != _nu.size()) || (_kA.size() != _sigma.size()) || (_kA.size() < nComp * numStates))
		throw InvalidParameterException("BISMA_KA, BISMA_KD, BISMA_NU, and BISMA_SIGMA have to have the same size");
	if (_lambda.size() != numStates)
		throw InvalidParameterException("BISMA_LAMBDA has to have as many elements as there are binding sites");

	return true;
}


/**
 * @brief Defines the Bi-SMA binding model
 * @details Implements the Bi-SMA adsorption model, which introduces multiple different binding sites and assumes an SMA model for each site: \f[ \begin{align} 
 *              q_0^A &= \Lambda^A - \sum_{j} \nu_j^A q_j^A \\
 *              q_0^B &= \Lambda^B - \sum_{j} \nu_j^B q_j^B \\
 *              \vodts & \vdots
 *              \frac{\mathrm{d}q_i^A}{\mathrm{d}t} &= k_{a,i}^A c_{p,i} \left( \Lambda^A - \sum_j\left( \nu_j^A + \sigma_j^A \right) q_j^A \right)^{\nu_i^A} - k_{d,i}^A q_i^A c_{p,0}^{\nu_i^A} \\
 *              \frac{\mathrm{d}q_i^B}{\mathrm{d}t} &= k_{a,i}^B c_{p,i} \left( \Lambda^B - \sum_j\left( \nu_j^B + \sigma_j^B \right) q_j^B \right)^{\nu_i^B} - k_{d,i}^B q_i^B c_{p,0}^{\nu_i^B} \\  
 *              \vodts & \vdots
 *          \end{align} \f]
 *          Here, several different types of binding sites @f$ q^A @f$, @f$ q^B @f$, etc. are considered. A molecule can either bind to
 *          site A or B (or C, etc.). A direct exchange between the different binding sites does not occur.
 *          While components without bound state (i.e., non-binding components) are supported, all other components must have
 *          the same number of bound states (i.e., binding sites). Component @c 0 is assumed to be salt.
 *          
 *          Since all algebraic variables have to be in one contiguous block, the order of the bound state variables is as follows: \f[ \begin{align} 
 *              q_0^A, q_0^B, q_0^C, \dots, q_1^A, q_2^A, q_3^A, \dots, q_1^B, q_2^B, q_3^B, \dots
 *          \end{align} \f]
 *          First, all the salt components are collected in one block as they are always algebraic.
 *          Then the other components for each bound state follow.
 * @tparam ParamHandler_t Type that can add support for external function dependence
 */
template <class ParamHandler_t>
class BiStericMassActionBindingBase : public BindingModelBase
{
public:

	BiStericMassActionBindingBase() { }
	virtual ~BiStericMassActionBindingBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }
	virtual const char* name() const CADET_NOEXCEPT { return ParamHandler_t::identifier(); }

	virtual void configureModelDiscretization(unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
	{
		BindingModelBase::configureModelDiscretization(nComp, nBound, boundOffset);

		unsigned int numSlices = 0;
		for (unsigned int i = 0; i < nComp; ++i)
		{
			if (nBound[i] == 0)
				continue;

			if (numSlices == 0)
				numSlices = nBound[i];

			if (nBound[i] != numSlices)
				throw InvalidParameterException("BiSMA binding model requires exactly the same bound states for all components");
		}

		_numBindingComp = numBindingComponents(_nBoundStates, _nComp);

		// Allocate space for parameters
		_paramHandler.reserve(nComp * numSlices, numSlices, nComp, nBound);

		// Guarantee that salt has exactly one bound state per binding site
		if (nBound[0] != numSlices)
			throw InvalidParameterException("BiSMA binding model requires exactly one bound state per binding site for salt component");
	}

	virtual void getAlgebraicBlock(unsigned int& idxStart, unsigned int& len) const
	{
		const unsigned int numStates = firstNonEmptyBoundStates(_nBoundStates, _nComp);

		// First equations are Salt, which are always algebraic
		idxStart = 0;
		if (_kineticBinding)
			len = numStates;
		else
			len = numBoundStates(_nBoundStates, _nComp);
	}


	virtual void consistentInitialState(double t, double z, double r, unsigned int secIdx, double* const vecStateY, double errorTol, 
		active* const adRes, active* const adY, unsigned int adEqOffset, unsigned int adDirOffset, const ad::IJacobianExtractor& jacExtractor, 
		double* const workingMemory, linalg::detail::DenseMatrixBase& workingMat) const
	{
		if (!_kineticBinding)
		{
			// All equations are algebraic and (except for salt equation) nonlinear
			// Compute the q_i from their corresponding c_{p,i}

			// Determine problem size
			const unsigned int eqSize = numBoundStates(_nBoundStates, _nComp);
			double* const workSpace = workingMemory + _nonlinearSolver->workspaceSize(eqSize);
			std::fill(workingMemory, workSpace, 0.0);

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
					residualImpl<active, double, active, double>(t, z, r, secIdx, 1.0, adY + adEqOffset, vecStateY - _nComp, nullptr, adRes + adEqOffset, workSpace);
					
#ifdef CADET_CHECK_ANALYTIC_JACOBIAN			
					// Compute analytic Jacobian
					mat.setAll(0.0);
					jacobianImpl(t, z, r, secIdx, x, vecStateY - _nComp, mat.row(0), workSpace); 

					// Compare
					const double diff = jacExtractor.compareWithJacobian(adRes, adEqOffset, adDirOffset, mat);
					LOG(Debug) << "MaxDiff " << adEqOffset << ": " << diff;
#endif
					// Extract Jacobian
					jacExtractor.extractJacobian(adRes, adEqOffset, adDirOffset, mat);
					return true;
				};
			}
			else
			{
				// Analytic Jacobian
				jacobianFunc = [&](double const* const x, linalg::detail::DenseMatrixBase& mat) -> bool { 
					mat.setAll(0.0);
					jacobianImpl(t, z, r, secIdx, x, vecStateY - _nComp, mat.row(0), workSpace); 
					return true;
				};
			}

			const bool conv = _nonlinearSolver->solve([&](double const* const x, double* const res) -> bool {
					residualImpl<double, double, double, double>(t, z, r, secIdx, 1.0, x, vecStateY - _nComp, nullptr, res, workSpace); 
					return true; 
				}, 
				jacobianFunc,
				errorTol, vecStateY, workingMemory, workingMat, eqSize);
		}

		// Compute salt component from given bound states q_j
		// This also corrects invalid salt values from nonlinear solver
		// in case of rapid equilibrium

		const typename ParamHandler_t::params_t& p = _paramHandler.update(t, z, r, secIdx, _nComp, _nBoundStates, workingMemory);
		const unsigned int numStates = p.lambda.size();

		// Loop over all binding site types
		for (unsigned int bndSite = 0; bndSite < numStates; ++bndSite)
		{
			// Salt equation: q_0 - Lambda + Sum[nu_j * q_j, j] == 0 
			//           <=>  q_0 == Lambda - Sum[nu_j * q_j, j] 
			vecStateY[bndSite] = static_cast<double>(p.lambda[bndSite]);

			// Get nu slice for bound state bndSite
			active const* const curNu = p.nu[bndSite];

			unsigned int bndIdx = 1;
			for (int j = 1; j < _nComp; ++j)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[j] == 0)
					continue;

				vecStateY[bndSite] -= static_cast<double>(curNu[j]) * vecStateY[bndIdx * numStates];

				// Next bound component
				++bndIdx;
			}
		}	
	}

	CADET_BINDINGMODEL_RESIDUAL_JACOBIAN_BOILERPLATE
	
	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { _paramHandler.setExternalFunctions(extFuns, size); }

	virtual void multiplyWithDerivativeJacobian(double const* yDotS, double* const res, double timeFactor) const
	{
		// Multiplier is 0 if quasi-stationary and 1 if kinetic binding
		const double multiplier = _kineticBinding ? timeFactor : 0.0;

		// First states are salt (always algebraic)
		const unsigned int numStates = firstNonEmptyBoundStates(_nBoundStates, _nComp);
		std::fill(res, res + numStates, 0.0);

		const unsigned int eqSize = numBoundStates(_nBoundStates, _nComp);
		for (unsigned int i = numStates; i < eqSize; ++i)
			res[i] = multiplier * yDotS[i];
	}

	virtual bool hasSalt() const CADET_NOEXCEPT { return true; }
	virtual bool supportsMultistate() const CADET_NOEXCEPT { return true; }
	virtual bool supportsNonBinding() const CADET_NOEXCEPT { return true; }
	virtual bool hasAlgebraicEquations() const CADET_NOEXCEPT { return true; }
	virtual bool dependsOnTime() const CADET_NOEXCEPT { return ParamHandler_t::dependsOnTime(); }
	virtual bool requiresWorkspace() const CADET_NOEXCEPT { return true; }

protected:
	ParamHandler_t _paramHandler; //!< Handles parameters and their dependence on external functions
	unsigned int _numBindingComp; //!< Number of binding components

	virtual unsigned int paramCacheSize() const CADET_NOEXCEPT { return _paramHandler.cacheSize(); }

	virtual bool configureImpl(bool reconfigure, IParameterProvider& paramProvider, unsigned int unitOpIdx)
	{
		// Read parameters
		_paramHandler.configure(paramProvider, _nComp, _nBoundStates);

		// Register parameters
		_paramHandler.registerParameters(_parameters, unitOpIdx, _nComp, _nBoundStates);

		return true;
	}

	template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
	int residualImpl(const ParamType& t, double z, double r, unsigned int secIdx, const ParamType& timeFactor,
		StateType const* y, CpStateType const* yCp, double const* yDot, ResidualType* res, void* workSpace) const
	{
		const typename ParamHandler_t::params_t& p = _paramHandler.update(static_cast<double>(t), z, r, secIdx, _nComp, _nBoundStates, workSpace);

		const bool hasYdot = yDot;
		const unsigned int numStates = p.lambda.size();

		// Ordering of the states is (q_{comp,state})
		// q_{0,0}, q{0,1}, q_{0,2}, q_{1,0}, q_{1,1}, q_{1,2}, ...
		// A state corresponds to a type of binding site. It is assumed that all components have either 0
		// or the same number of states. Thus, a component is either non-binding or has nSites bound states.
		//
		// The same ordering is used for the equations. That is, q_{0,0}, q_{1,0} and q_{0,1}, q_{1,1} and ... each
		// form one SMA binding model system.

		// Loop over all binding site types
		for (unsigned int bndSite = 0; bndSite < numStates; ++bndSite, ++y, ++yDot, ++res)
		{
			// y, yDot, and res point to q_{0,site}

			active const* const curNu = p.nu[bndSite];
			active const* const curSigma = p.sigma[bndSite];
			active const* const curKa = p.kA[bndSite];
			active const* const curKd = p.kD[bndSite];

			const ParamType refC0 = static_cast<ParamType>(p.refC0[bndSite]);
			const ParamType refQ = static_cast<ParamType>(p.refQ[bndSite]);

			// Salt equation: q_0 - Lambda + Sum[nu_j * q_j, j] == 0 
			//           <=>  q_0 == Lambda - Sum[nu_j * q_j, j] 
			// Also compute \bar{q}_0 = q_0 - Sum[sigma_j * q_j, j]
			res[0] = y[0] - static_cast<ParamType>(p.lambda[bndSite]);
			ResidualType q0_bar = y[0];

			// bndIdx is used as a counter inside one binding site type
			// Getting from one component to another requires a step size of numStates (stride)
			unsigned int bndIdx = 1;
			for (int j = 1; j < _nComp; ++j)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[j] == 0)
					continue;

				res[0] += static_cast<ParamType>(curNu[j]) * y[bndIdx * numStates];
				q0_bar -= static_cast<ParamType>(curSigma[j]) * y[bndIdx * numStates];

				// Next bound component
				++bndIdx;
			}

			const ResidualType q0_bar_divRef = q0_bar / refQ;
			const ResidualType yCp0_divRef = yCp[0] / refC0;

			// Protein equations: dq_i / dt - ( k_{a,i} * c_{p,i} * \bar{q}_0^{nu_i} - k_{d,i} * q_i * c_{p,0}^{nu_i} ) == 0
			//               <=>  dq_i / dt == k_{a,i} * c_{p,i} * \bar{q}_0^{nu_i} - k_{d,i} * q_i * c_{p,0}^{nu_i}
			bndIdx = 1;
			for (int i = 1; i < _nComp; ++i)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[i] == 0)
					continue;

				const ResidualType c0_pow_nu_divRef = pow(yCp0_divRef, static_cast<ParamType>(curNu[i]));
				const ResidualType q0_bar_pow_nu_divRef = pow(q0_bar_divRef, static_cast<ParamType>(curNu[i]));

				// Residual
				res[bndIdx * numStates] = static_cast<ParamType>(curKd[i]) * y[bndIdx * numStates] * c0_pow_nu_divRef - static_cast<ParamType>(curKa[i]) * yCp[i] * q0_bar_pow_nu_divRef;

				// Add time derivative if necessary
				if (_kineticBinding && hasYdot)
				{
					res[bndIdx * numStates] += timeFactor * yDot[bndIdx * numStates];
				}

				// Next bound component
				++bndIdx;
			}
		}
		return 0;
	}

	template <typename RowIterator>
	void jacobianImpl(double t, double z, double r, unsigned int secIdx, double const* y, double const* yCp, RowIterator jac, void* workSpace) const
	{
		const typename ParamHandler_t::params_t& p = _paramHandler.update(t, z, r, secIdx, _nComp, _nBoundStates, workSpace);

		const unsigned int numStates = p.lambda.size();

		// Ordering of the states is (q_{comp,state}, example uses 2 components, 3 binding sites)
		// q_{0,0}, q{0,1}, q_{0,2}, q_{1,0}, q_{1,1}, q_{1,2}, ...
		// Ordering of the equations is the same, that is, we need nSites steps to jump from one equation of an SMA
		// binding model system to the next.

		// Loop over all binding site types
		for (unsigned int bndSite = 0; bndSite < numStates; ++bndSite, ++y)
		{
			// Jump from first row to current Salt row
			RowIterator curJac = jac + bndSite;

			active const* const curNu = p.nu[bndSite];
			active const* const curSigma = p.sigma[bndSite];
			active const* const curKa = p.kA[bndSite];
			active const* const curKd = p.kD[bndSite];
			const double refC0 = static_cast<double>(p.refC0[bndSite]);
			const double refQ = static_cast<double>(p.refQ[bndSite]);
			double q0_bar = y[0];

			// Salt equation: q_0 - Lambda + Sum[nu_j * q_j, j] == 0
			curJac[0] = 1.0;

			int bndIdx = 1;
			for (int j = 1; j < _nComp; ++j)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[j] == 0)
					continue;

				curJac[bndIdx * numStates] = static_cast<double>(curNu[j]);

				// Calculate \bar{q}_0 = q_0 - Sum[sigma_j * q_j, j]
				q0_bar -= static_cast<double>(curSigma[j]) * y[bndIdx * numStates];

				// Next bound component
				++bndIdx;
			}

			// Advance to protein equations
			curJac += numStates;

			const double q0_bar_divRef = q0_bar / refQ;
			const double yCp0_divRef = yCp[0] / refC0;

			// Protein equations: dq_i / dt - ( k_{a,i} * c_{p,i} * \bar{q}_0^{nu_i} - k_{d,i} * q_i * c_{p,0}^{nu_i} ) == 0
			// We have already computed \bar{q}_0 in the loop above
			bndIdx = 1;
			for (int i = 1; i < _nComp; ++i)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[i] == 0)
					continue;

				// Getting to c_{p,0}: -numStates * bndIdx takes us to q_{0,site}, another -bndSite to q_{0,0}. From there, we
				//                     take a -_nComp to reach c_{p,0}.
				//                     This means jac[-bndSite - _nComp - numStates * bndIdx] corresponds to c_{p,0}.
				// Getting to c_{p,i}: Go to c_{p,0} and add i.
				//                     This means jac[i - bndSite - _nComp - numStates * bndIdx] corresponds to c_{p,i}.

				const double ka = static_cast<double>(curKa[i]);
				const double kd = static_cast<double>(curKd[i]);
				const double nu = static_cast<double>(curNu[i]);

				const double c0_pow_nu     = pow(yCp0_divRef, nu);
				const double q0_bar_pow_nu = pow(q0_bar_divRef, nu);
				const double c0_pow_nu_m1_divRef     = pow(yCp0_divRef, nu - 1.0);
				const double q0_bar_pow_nu_m1_divRef = pow(q0_bar_divRef, nu - 1.0) / refQ;

				// dres_i / dc_{p,0}
				curJac[-bndSite - _nComp - numStates * bndIdx] = kd * y[bndIdx * numStates] * nu * c0_pow_nu_m1_divRef / refC0;
				// dres_i / dc_{p,i}
				curJac[i - bndSite - _nComp - numStates * bndIdx] = -ka * q0_bar_pow_nu;
				// dres_i / dq_{0,bndSite}
				curJac[-bndIdx * numStates] = -ka * yCp[i] * nu * q0_bar_pow_nu_m1_divRef;

				// Fill dres_i / dq_{j,bndSite}
				int bndIdx2 = 1;
				for (int j = 1; j < _nComp; ++j)
				{
					// Skip components without bound states (bound state index bndIdx is not advanced)
					if (_nBoundStates[j] == 0)
						continue;

					// dres_i / dq_{j,bndSite}
					curJac[(bndIdx2 - bndIdx) * numStates] = -ka * yCp[i] * nu * q0_bar_pow_nu_m1_divRef * (-static_cast<double>(curSigma[j]));
					// Getting to q_j: -bndIdx * numStates takes us to q_{0,bndSite}, another +bndIdx2 * numStates to
					// q_{j,bndSite}. This means curJac[(bndIdx2 - bndIdx) * numStates] corresponds to q_{j,bndSite}.

					++bndIdx2;
				}

				// Add to dres_i / dq_{i,bndSite}
				curJac[0] += kd * c0_pow_nu;

				// Advance to next equation
				++bndIdx;
				curJac += numStates;
				// Note that there is a spacing of numStates between the equations inside one binding site type
			}
		}
	}

	template <typename RowIterator>
	void jacobianAddDiscretizedImpl(double alpha, RowIterator jac) const
	{
		// We only add time derivatives for kinetic binding
		if (!_kineticBinding)
			return;

		// Skip salt equations which are always algebraic
		const unsigned int numStates = firstNonEmptyBoundStates(_nBoundStates, _nComp);
		jac += numStates;

		const unsigned int eqSize = numBoundStates(_nBoundStates, _nComp) - numStates;
		for (unsigned int i = 0; i < eqSize; ++i, ++jac)
			jac[0] += alpha;
	}
};


typedef BiStericMassActionBindingBase<BiSMAParamHandler> BiStericMassActionBinding;
typedef BiStericMassActionBindingBase<ExtBiSMAParamHandler> ExternalBiStericMassActionBinding;

namespace binding
{
	void registerBiStericMassActionModel(std::unordered_map<std::string, std::function<model::IBindingModel*()>>& bindings)
	{
		bindings[BiStericMassActionBinding::identifier()] = []() { return new BiStericMassActionBinding(); };
		bindings[ExternalBiStericMassActionBinding::identifier()] = []() { return new ExternalBiStericMassActionBinding(); };
	}
}  // namespace binding

}  // namespace model

}  // namespace cadet
