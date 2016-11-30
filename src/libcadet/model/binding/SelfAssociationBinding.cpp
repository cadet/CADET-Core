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

#include "model/binding/BindingModelBase.hpp"
#include "model/binding/ExternalFunctionSupport.hpp"
#include "model/binding/BindingModelMacros.hpp"
#include "model/binding/RefConcentrationSupport.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "nonlin/Solver.hpp"
#include "ParamReaderHelper.hpp"

#include <functional>
#include <unordered_map>
#include <string>
#include <vector>

#include "AdUtils.hpp"
#include "linalg/Norms.hpp"

#include "LoggingUtils.hpp"
#include "Logging.hpp"

#include <iostream>

namespace cadet
{

namespace model
{

/**
 * @brief Handles self association binding model parameters that do not depend on external functions
 */
struct SelfAssociationParamHandler : public BindingParamHandlerBase
{
	static const char* identifier() { return "SELF_ASSOCIATION"; }

	/**
	 * @brief Reads parameters and verifies them
	 * @details See IBindingModel::configure() for details.
	 * @param [in] paramProvider IParameterProvider used for reading parameters
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 * @return @c true if the parameters were read and validated successfully, otherwise @c false
	 */
	inline bool configure(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBoundStates)
	{
		lambda = paramProvider.getDouble("SAI_LAMBDA");
		readParameterMatrix(kA, paramProvider, "SAI_KA1", nComp, 1);
		readParameterMatrix(kA2, paramProvider, "SAI_KA2", nComp, 1);
		readParameterMatrix(kD, paramProvider, "SAI_KD", nComp, 1);
		readParameterMatrix(nu, paramProvider, "SAI_NU", nComp, 1);
		readParameterMatrix(sigma, paramProvider, "SAI_SIGMA", nComp, 1);
		readReferenceConcentrations(paramProvider, "SAI_", refC0, refQ);

		// Check parameters
		if ((kA.size() != kA2.size()) || (kA.size() != kD.size()) || (kA.size() != nu.size()) 
			|| (kA.size() != sigma.size()) || (kA.size() < nComp))
			throw InvalidParameterException("SAI_KA1, SAI_KA2, SAI_KD, SAI_NU, and SAI_SIGMA have to have the same size");

		return true;
	}

	/**
	 * @brief Registers all local parameters in a map for further use
	 * @param [in,out] parameters Map in which the parameters are stored
	 * @param [in] unitOpIdx Index of the unit operation used for registering the parameters
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 */
	inline void registerParameters(std::unordered_map<ParameterId, active*>& parameters, unsigned int unitOpIdx, unsigned int nComp, unsigned int const* nBoundStates)
	{
		parameters[makeParamId(hashString("SAI_LAMBDA"), unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &lambda;
		registerComponentBoundStateDependentParam(hashString("SAI_KA1"), parameters, kA, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("SAI_KA2"), parameters, kA2, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("SAI_KD"), parameters, kD, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("SAI_NU"), parameters, nu, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("SAI_SIGMA"), parameters, sigma, unitOpIdx);
		parameters[makeParamId(hashString("SAI_REFC0"), unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &refC0;
		parameters[makeParamId(hashString("SAI_REFQ"), unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &refQ;
	}

	active lambda; //! Ionic capacity
	std::vector<active> kA; //!< Adsorption rate
	std::vector<active> kA2; //!< Adsorption rate of dimerization
	std::vector<active> kD; //!< Desorption rate
	std::vector<active> nu; //!< Characteristic charge
	std::vector<active> sigma; //!< Steric factor
	active refC0; //! Liquid phase reference concentration
	active refQ; //! Solid phase reference concentration
};

/**
 * @brief Handles self association binding model parameters that depend on an external function
 */
struct ExtSelfAssociationParamHandler : public ExternalBindingParamHandlerBase
{
	static const char* identifier() { return "EXT_SELF_ASSOCIATION"; }

	/**
	 * @brief Reads parameters and verifies them
	 * @details See IBindingModel::configure() for details.
	 * @param [in] paramProvider IParameterProvider used for reading parameters
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 * @return @c true if the parameters were read and validated successfully, otherwise @c false
	 */
	inline bool configure(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBoundStates)
	{
		CADET_READPAR_SCALAR(lambda, paramProvider, "SAI_LAMBDA");
		CADET_READPAR_MATRIX(kA, paramProvider, "SAI_KA1", nComp, 1);
		CADET_READPAR_MATRIX(kA2, paramProvider, "SAI_KA2", nComp, 1);
		CADET_READPAR_MATRIX(kD, paramProvider, "SAI_KD", nComp, 1);
		CADET_READPAR_MATRIX(nu, paramProvider, "SAI_NU", nComp, 1);
		CADET_READPAR_MATRIX(sigma, paramProvider, "SAI_SIGMA", nComp, 1);
		readReferenceConcentrations(paramProvider, "EXT_SAI_", refC0, refQ);

		return ExternalBindingParamHandlerBase::configure(paramProvider, 6);
	}

	/**
	 * @brief Registers all local parameters in a map for further use
	 * @param [in,out] parameters Map in which the parameters are stored
	 * @param [in] unitOpIdx Index of the unit operation used for registering the parameters
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 */
	inline void registerParameters(std::unordered_map<ParameterId, active*>& parameters, unsigned int unitOpIdx, unsigned int nComp, unsigned int const* nBoundStates)
	{
		CADET_REGPAR_SCALAR("SAI_LAMBDA", parameters, lambda, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("SAI_KA1", parameters, kA, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("SAI_KA2", parameters, kA2, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("SAI_KD", parameters, kD, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("SAI_NU", parameters, nu, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("SAI_SIGMA", parameters, sigma, unitOpIdx);
		parameters[makeParamId(hashString("EXT_SAI_REFC0"), unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &refC0;
		parameters[makeParamId(hashString("EXT_SAI_REFQ"), unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &refQ;
	}

	/**
	 * @brief Updates local parameter cache in order to take the external profile into account
	 * @details This function is declared const since the actual parameters are left unchanged by the method.
	 *         The cache is marked as mutable in order to make it writable.
	 * @param [in] t Current time
	 * @param [in] z Axial coordinate in the column
	 * @param [in] r Radial coordinate in the bead
	 * @param [in] secIdx Index of the current section
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 */
	inline void update(double t, double z, double r, unsigned int secIdx, unsigned int nComp, unsigned int const* nBoundStates) const
	{
		evaluateExternalFunctions(t, z, r, secIdx);
		for (unsigned int i = 0; i < nComp; ++i)
		{
			CADET_UPDATE_EXTDEP_VARIABLE_BRACES(kA, i, _extFunBuffer[0]);
			CADET_UPDATE_EXTDEP_VARIABLE_BRACES(kA2, i, _extFunBuffer[1]);
			CADET_UPDATE_EXTDEP_VARIABLE_BRACES(kD, i, _extFunBuffer[2]);
			CADET_UPDATE_EXTDEP_VARIABLE_BRACES(nu, i, _extFunBuffer[3]);
			CADET_UPDATE_EXTDEP_VARIABLE_BRACES(sigma, i, _extFunBuffer[4]);
		}

		CADET_UPDATE_EXTDEP_VARIABLE(lambda, _extFunBuffer[5]);
	}

	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, kA)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, kA2)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, kD)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, nu)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, sigma)
	CADET_DEFINE_EXTDEP_VARIABLE(active, lambda)
	active refC0; //! Liquid phase reference concentration
	active refQ; //! Solid phase reference concentration
};


/**
 * @brief Defines the self association binding model
 * @details Implements the self association adsorption model: \f[ \begin{align} 
 *              q_0 &= \Lambda - \sum_{j} \nu_j q_j \\
 *              \frac{\mathrm{d}q_i}{\mathrm{d}t} &= \left[k_{a1,i} + k_{a2,i} c_{p,i} \right] c_{p,i} \left( \Lambda - \sum_j\left( \nu_j + \sigma_j \right) q_j \right)^{\nu_i} - k_{d,i} q_i c_{p,0}^{\nu_i} 
 *          \end{align} \f]
 *          Here, the coefficient @f$ k_{a2,i} @f$ denotes the adsorption rate of dimerization.
 *          Component @c 0 is assumed to be salt. Multiple bound states are not supported. 
 *          Components without bound state (i.e., non-binding components) are supported.
 *          
 *          See @cite Mollerup2008 and @cite Westerberg2012.
 * @tparam ParamHandler_t Type that can add support for external function dependence
 */
template <class ParamHandler_t>
class SelfAssociationBindingBase : public BindingModelBase
{
public:

	SelfAssociationBindingBase() { }
	virtual ~SelfAssociationBindingBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }
	virtual const char* name() const CADET_NOEXCEPT { return ParamHandler_t::identifier(); }

	virtual void configureModelDiscretization(unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
	{
		BindingModelBase::configureModelDiscretization(nComp, nBound, boundOffset);

		// Guarantee that salt has exactly one bound state
		if (nBound[0] != 1)
			throw InvalidParameterException("Self association binding model requires exactly one bound state for salt component");
	}

	virtual void getAlgebraicBlock(unsigned int& idxStart, unsigned int& len) const
	{
		// First equation is Salt, which is always algebraic
		idxStart = 0;
		if (_kineticBinding)
			len = 1;
		else
			len = numBoundStates(_nBoundStates, _nComp);
	}

	virtual unsigned int consistentInitializationWorkspaceSize() const
	{
		// Determine problem size
		const unsigned int eqSize = numBoundStates(_nBoundStates, _nComp);
		// Ask nonlinear solver how much memory it needs for this kind of problem
		return _nonlinearSolver->workspaceSize(eqSize);
	}

	virtual void consistentInitialState(double t, double z, double r, unsigned int secIdx, double* const vecStateY, double errorTol, 
		active* const adRes, active* const adY, unsigned int adEqOffset, unsigned int adOffset, unsigned int diagDir, 
		unsigned int lowerBandwidth, unsigned int upperBandwidth, double* const workingMemory, linalg::detail::DenseMatrixBase& workingMat) const
	{
		_p.update(t, z, r, secIdx, _nComp, _nBoundStates);

		// If we have kinetic binding, there is only the algebraic salt equation
		if (_kineticBinding)
		{
			// Compute salt component from given bound states q_j
			// Salt equation: q_0 - Lambda + Sum[nu_j * q_j, j] == 0 
			//           <=>  q_0 == Lambda - Sum[nu_j * q_j, j] 
			vecStateY[0] = 0.0;
			double kahanCompensation = 0.0;

			unsigned int bndIdx = 1;
			for (int j = 1; j < _nComp; ++j)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[j] == 0)
					continue;

				const double diff = -static_cast<double>(_p.nu[j]) * vecStateY[bndIdx] - kahanCompensation;
				const double tempSum = vecStateY[0] + diff;
				kahanCompensation = (tempSum - vecStateY[0]) - diff;
				vecStateY[0] = tempSum;

//				vecStateY[0] -= static_cast<double>(_p.nu[j]) * vecStateY[bndIdx];

				// Next bound component
				++bndIdx;
			}

			const double diff = static_cast<double>(_p.lambda) - kahanCompensation;
			const double tempSum = vecStateY[0] + diff;
			vecStateY[0] = tempSum;

			return;
		}

		// All equations are algebraic and (except for salt equation) nonlinear
		// Compute the q_i from their corresponding c_{p,i}

		// Determine problem size
		const unsigned int eqSize = numBoundStates(_nBoundStates, _nComp);
		std::fill(workingMemory, workingMemory + _nonlinearSolver->workspaceSize(eqSize), 0.0);

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
				residualImpl<active, double, active, double>(t, z, r, secIdx, 1.0, adY + adEqOffset, vecStateY - _nComp, nullptr, adRes + adEqOffset);
				
#ifdef CADET_CHECK_ANALYTIC_JACOBIAN			
				// Compute analytic Jacobian
				mat.setAll(0.0);
				jacobianImpl(t, z, r, secIdx, x, vecStateY - _nComp, mat.row(0)); 

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
			jacobianFunc = [&](double const* const x, linalg::detail::DenseMatrixBase& mat) -> bool { 
				mat.setAll(0.0);
				jacobianImpl(t, z, r, secIdx, x, vecStateY - _nComp, mat.row(0)); 
				return true;
			};
		}

		const bool conv = _nonlinearSolver->solve([&](double const* const x, double* const res) -> bool {
				residualImpl<double, double, double, double>(t, z, r, secIdx, 1.0, x, vecStateY - _nComp, nullptr, res); 
				return true; 
			}, 
			jacobianFunc,
			errorTol, vecStateY, workingMemory, workingMat, eqSize);
	}

	virtual int residual(const active& t, double z, double r, unsigned int secIdx, const active& timeFactor, 
		active const* y, double const* yDot, active* res) const;

	virtual int residual(double t, double z, double r, unsigned int secIdx, double timeFactor, 
		active const* y, double const* yDot, active* res) const;

	virtual int residual(const active& t, double z, double r, unsigned int secIdx, const active& timeFactor, 
		double const* y, double const* yDot, active* res) const;

	virtual int residual(double t, double z, double r, unsigned int secIdx, double timeFactor, 
		double const* y, double const* yDot, double* res) const;

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { _p.setExternalFunctions(extFuns, size); }

	virtual void analyticJacobian(double t, double z, double r, unsigned int secIdx, double const* y, linalg::BandMatrix::RowIterator jac) const
	{
		jacobianImpl(t, z, r, secIdx, y, y - _nComp, jac);
	}

	virtual void jacobianAddDiscretized(double alpha, linalg::FactorizableBandMatrix::RowIterator jac) const
	{
		// We only add time derivatives for kinetic binding
		if (!_kineticBinding)
			return;

		// Skip salt equation which is always algebraic
		++jac;

		const unsigned int eqSize = numBoundStates(_nBoundStates, _nComp) - 1;
		for (unsigned int i = 0; i < eqSize; ++i, ++jac)
			jac[0] += alpha;
	}

	virtual void multiplyWithDerivativeJacobian(double const* yDotS, double* const res, double timeFactor) const
	{
		// Multiplier is 0 if quasi-stationary and 1 if kinetic binding
		// However, due to premultiplication of time derivatives with the timeFactor (because of time transformation),
		// we set it to timeFactor instead of 1.0 in order to save some multiplications
		const double multiplier = _kineticBinding ? timeFactor : 0.0;

		// First state is salt (always algebraic)
		res[0] = 0.0;

		const unsigned int eqSize = numBoundStates(_nBoundStates, _nComp);
		for (unsigned int i = 1; i < eqSize; ++i)
			res[i] = multiplier * yDotS[i];
	}

	virtual bool hasSalt() const CADET_NOEXCEPT { return true; }
	virtual bool supportsMultistate() const CADET_NOEXCEPT { return false; }
	virtual bool supportsNonBinding() const CADET_NOEXCEPT { return true; }
	virtual bool hasAlgebraicEquations() const CADET_NOEXCEPT { return true; }

protected:
	ParamHandler_t _p; //!< Handles parameters and their dependence on external functions

	virtual bool configureImpl(bool reconfigure, IParameterProvider& paramProvider, unsigned int unitOpIdx)
	{
		// Read parameters
		_p.configure(paramProvider, _nComp, _nBoundStates);

		// Register parameters
		_p.registerParameters(_parameters, unitOpIdx, _nComp, _nBoundStates);

		return true;
	}

	template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
	int residualImpl(const ParamType& t, double z, double r, unsigned int secIdx, const ParamType& timeFactor,
		StateType const* y, CpStateType const* yCp, double const* yDot, ResidualType* res) const
	{
		_p.update(static_cast<double>(t), z, r, secIdx, _nComp, _nBoundStates);

		// Salt equation: q_0 - Lambda + Sum[nu_j * q_j, j] == 0 
		//           <=>  q_0 == Lambda - Sum[nu_j * q_j, j] 
		// Also compute \bar{q}_0 = q_0 - Sum[sigma_j * q_j, j]
		res[0] = y[0] - static_cast<ParamType>(_p.lambda);
		ResidualType q0_bar = y[0];

		unsigned int bndIdx = 1;
		for (int j = 1; j < _nComp; ++j)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[j] == 0)
				continue;

			res[0] += static_cast<ParamType>(_p.nu[j]) * y[bndIdx];
			q0_bar -= static_cast<ParamType>(_p.sigma[j]) * y[bndIdx];

			// Next bound component
			++bndIdx;
		}

		const ResidualType q0_bar_divRef = q0_bar / static_cast<ParamType>(_p.refQ);
		const ResidualType yCp0_divRef = yCp[0] / static_cast<ParamType>(_p.refC0);

		// Protein equations: dq_i / dt - ( (k_{a,i} + k_{a2,i} * c_{p,i}) * c_{p,i} * \bar{q}_0^{nu_i} - k_{d,i} * q_i * c_{p,0}^{nu_i} ) == 0
		//               <=>  dq_i / dt == (k_{a,i} + k_{a2,i} * c_{p,i}) * c_{p,i} * \bar{q}_0^{nu_i} - k_{d,i} * q_i * c_{p,0}^{nu_i}
		bndIdx = 1;
		for (int i = 1; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			const ResidualType c0_pow_nu = pow(yCp0_divRef, static_cast<ParamType>(_p.nu[i]));
			const ResidualType q0_bar_pow_nu = pow(q0_bar_divRef, static_cast<ParamType>(_p.nu[i]));

			// Residual
			res[bndIdx] = static_cast<ParamType>(_p.kD[i]) * y[bndIdx] * c0_pow_nu - yCp[i] * (static_cast<ParamType>(_p.kA[i]) + yCp[i] * static_cast<ParamType>(_p.kA2[i])) * q0_bar_pow_nu;

			// Add time derivative if necessary
			if (_kineticBinding && yDot)
			{
				res[bndIdx] += timeFactor * yDot[bndIdx];
			}

			// Next bound component
			++bndIdx;
		}

		return 0;
	}

	template <typename RowIterator>
	void jacobianImpl(double t, double z, double r, unsigned int secIdx, double const* y, double const* yCp, RowIterator jac) const
	{
		_p.update(t, z, r, secIdx, _nComp, _nBoundStates);

		double q0_bar = y[0];

		// Salt equation: q_0 - Lambda + Sum[nu_j * q_j, j] == 0
		jac[0] = 1.0;
		int bndIdx = 1;
		for (int j = 1; j < _nComp; ++j)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[j] == 0)
				continue;

			jac[bndIdx] = static_cast<double>(_p.nu[j]);

			// Calculate \bar{q}_0 = q_0 - Sum[sigma_j * q_j, j]
			q0_bar -= static_cast<double>(_p.sigma[j]) * y[bndIdx];

			// Next bound component
			++bndIdx;
		}

		// Advance to protein equations
		++jac;

		const double refC0 = static_cast<double>(_p.refC0);
		const double refQ = static_cast<double>(_p.refQ);
		const double yCp0_divRef = yCp[0] / refC0;
		const double q0_bar_divRef = q0_bar / refQ;

		// Protein equations: dq_i / dt - ( (k_{a,i} + k_{a2,i} * c_{p,i}) * c_{p,i} * \bar{q}_0^{nu_i} - k_{d,i} * q_i * c_{p,0}^{nu_i} ) == 0
		// We have already computed \bar{q}_0 in the loop above
		bndIdx = 1;
		for (int i = 1; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			// Getting to c_{p,0}: -bndIdx takes us to q_0, another -nComp to c_{p,0}. This means jac[-bndIdx - nComp] corresponds to c_{p,0}.
			// Getting to c_{p,i}: -bndIdx takes us to q_0, another -nComp to c_{p,0} and a +i to c_{p,i}.
			//                     This means jac[i - bndIdx - nComp] corresponds to c_{p,i}.

			const double ka = static_cast<double>(_p.kA[i]);
			const double ka2 = static_cast<double>(_p.kA2[i]);
			const double kd = static_cast<double>(_p.kD[i]);
			const double nu = static_cast<double>(_p.nu[i]);

			const double effKa = yCp[i] * (ka + ka2 * yCp[i]);
			const double c0_pow_nu     = pow(yCp0_divRef, nu);
			const double q0_bar_pow_nu = pow(q0_bar_divRef, nu);
			const double c0_pow_nu_m1_divRef     = pow(yCp0_divRef, nu - 1.0) / refC0;
			const double q0_bar_pow_nu_m1_divRef = pow(q0_bar_divRef, nu - 1.0) / refQ;

			// dres_i / dc_{p,0}
			jac[-bndIdx - _nComp] = kd * y[bndIdx] * nu * c0_pow_nu_m1_divRef;
			// dres_i / dc_{p,i}
			jac[i - bndIdx - _nComp] = -(ka + 2.0 * ka2 * yCp[i]) * q0_bar_pow_nu;
			// dres_i / dq_0
			jac[-bndIdx] = -effKa * yCp[i] * nu * q0_bar_pow_nu_m1_divRef;

			// Fill dres_i / dq_j
			int bndIdx2 = 1;
			for (int j = 1; j < _nComp; ++j)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[j] == 0)
					continue;

				// dres_i / dq_j
				jac[bndIdx2 - bndIdx] = -effKa * yCp[i] * nu * q0_bar_pow_nu_m1_divRef * (-static_cast<double>(_p.sigma[j]));
				// Getting to q_j: -bndIdx takes us to q_0, another +bndIdx2 to q_j. This means jac[bndIdx2 - bndIdx] corresponds to q_j.

				++bndIdx2;
			}

			// Add to dres_i / dq_i
			jac[0] += kd * c0_pow_nu;

			// Advance to next equation and Jacobian row
			++bndIdx;
			++jac;
		}
	}
};

CADET_BINDINGMODEL_RESIDUAL_TEMPLATED_BOILERPLATE_IMPL(SelfAssociationBindingBase, ParamHandler_t)


typedef SelfAssociationBindingBase<SelfAssociationParamHandler> SelfAssociationBinding;
typedef SelfAssociationBindingBase<ExtSelfAssociationParamHandler> ExternalSelfAssociationBinding;

namespace binding
{
	void registerSelfAssociationModel(std::unordered_map<std::string, std::function<model::IBindingModel*()>>& bindings)
	{
		bindings[SelfAssociationBinding::identifier()] = []() { return new SelfAssociationBinding(); };
		bindings[ExternalSelfAssociationBinding::identifier()] = []() { return new ExternalSelfAssociationBinding(); };
	}
}  // namespace binding

}  // namespace model

}  // namespace cadet
