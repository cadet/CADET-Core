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
#include "model/binding/ExternalFunctionSupport.hpp"
#include "model/binding/BindingModelMacros.hpp"
#include "model/binding/RefConcentrationSupport.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "nonlin/Solver.hpp"
#include "ParamReaderHelper.hpp"
#include "SlicedVector.hpp"

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

/**
 * @brief Handles MSSMA binding model parameters that do not depend on external functions
 */
struct MSSMAParamHandler : public BindingParamHandlerBase
{
	static const char* identifier() { return "MULTISTATE_STERIC_MASS_ACTION"; }

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
		const unsigned int totalBoundStates = numBoundStates(nBoundStates, nComp);

		lambda = paramProvider.getDouble("MSSMA_LAMBDA");
		readBoundStateDependentParameter<util::SlicedVector<active>, active>(kA, paramProvider, "MSSMA_KA", nComp, nBoundStates);
		readBoundStateDependentParameter<util::SlicedVector<active>, active>(kD, paramProvider, "MSSMA_KD", nComp, nBoundStates);
		readBoundStateDependentParameter<util::SlicedVector<active>, active>(nu, paramProvider, "MSSMA_NU", nComp, nBoundStates);
		readBoundStateDependentParameter<util::SlicedVector<active>, active>(sigma, paramProvider, "MSSMA_SIGMA", nComp, nBoundStates);
		readMatrixValuedBoundStateDependentParameter<util::SlicedVector<active>, active>(kRate, paramProvider, "MSSMA_RATES", nComp, nBoundStates);
		readReferenceConcentrations(paramProvider, "MSSMA_", refC0, refQ);

		// Check parameters
		if ((kA.size() != kD.size()) || (kA.size() != nu.size()) || (kA.size() != sigma.size()) || (kA.size() != totalBoundStates))
			throw InvalidParameterException("MSSMA_KA, MSSMA_KD, MSSMA_NU, and MSSMA_SIGMA have to have the same size");

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
		parameters[makeParamId(hashString("MSSMA_LAMBDA"), unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &lambda;
		registerComponentBoundStateDependentParamCompMajor(hashString("MSSMA_KA"), parameters, kA, unitOpIdx);
		registerComponentBoundStateDependentParamCompMajor(hashString("MSSMA_KD"), parameters, kD, unitOpIdx);
		registerComponentBoundStateDependentParamCompMajor(hashString("MSSMA_NU"), parameters, nu, unitOpIdx);
		registerComponentBoundStateDependentParamCompMajor(hashString("MSSMA_SIGMA"), parameters, sigma, unitOpIdx);
		registerComponentBoundStateDependentParamCompMajor(hashString("MSSMA_RATES"), parameters, kRate, unitOpIdx);
		parameters[makeParamId(hashString("MSSMA_REFC0"), unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &refC0;
		parameters[makeParamId(hashString("MSSMA_REFQ"), unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &refQ;
	}

	/**
	 * @brief Reserves space in the storage of the parameters
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 * @param [in] totalBoundStates Total number of bound states
	 */
	inline void reserve(unsigned int nComp, unsigned int const* nBoundStates, unsigned int totalBoundStates)
	{
		kA.reserve(totalBoundStates, nComp);
		kD.reserve(totalBoundStates, nComp);
		nu.reserve(totalBoundStates, nComp);
		sigma.reserve(totalBoundStates, nComp);

		unsigned int sumSquared = 0;
		for (unsigned int i = 0; i < nComp; ++i)
			sumSquared += nBoundStates[i] * nBoundStates[i];
		kRate.reserve(sumSquared, nComp);
	}

	active lambda; //! Ionic capacity
	util::SlicedVector<active> kA; //!< Adsorption rate in component-major ordering
	util::SlicedVector<active> kD; //!< Desorption rate in component-major ordering
	util::SlicedVector<active> nu; //!< Characteristic charge in component-major ordering
	util::SlicedVector<active> sigma; //!< Steric factor in component-major ordering
	util::SlicedVector<active> kRate; //!< State transition rates @f$ k_{\ell j}^{(i)} @f$ in component-row-major ordering
	active refC0; //! Liquid phase reference concentration
	active refQ; //! Solid phase reference concentration
};

/**
 * @brief Handles MSSMA binding model parameters that depend on an external function
 */
struct ExtMSSMAParamHandler : public ExternalBindingParamHandlerBase
{
	static const char* identifier() { return "EXT_MULTISTATE_STERIC_MASS_ACTION"; }

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
		CADET_READPAR_SCALAR(lambda, paramProvider, "MSSMA_LAMBDA");
		CADET_READPAR_BOUNDSTATEDEP(util::SlicedVector<active>, active, kA, paramProvider, "MSSMA_KA", nComp, nBoundStates);
		CADET_READPAR_BOUNDSTATEDEP(util::SlicedVector<active>, active, kD, paramProvider, "MSSMA_KD", nComp, nBoundStates);
		CADET_READPAR_BOUNDSTATEDEP(util::SlicedVector<active>, active, nu, paramProvider, "MSSMA_NU", nComp, nBoundStates);
		CADET_READPAR_BOUNDSTATEDEP(util::SlicedVector<active>, active, sigma, paramProvider, "MSSMA_SIGMA", nComp, nBoundStates);
		CADET_READPAR_BOUNDSTATEDEP_MATRIX(util::SlicedVector<active>, active, kRate, paramProvider, "MSSMA_RATES", nComp, nBoundStates);
		readReferenceConcentrations(paramProvider, "EXT_MSSMA_", refC0, refQ);

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
		CADET_REGPAR_SCALAR("MSSMA_LAMBDA", parameters, lambda, unitOpIdx);
		CADET_REGPAR_COMPBND_COMPMAJOR("MSSMA_KA", parameters, kA, unitOpIdx);
		CADET_REGPAR_COMPBND_COMPMAJOR("MSSMA_KD", parameters, kD, unitOpIdx);
		CADET_REGPAR_COMPBND_COMPMAJOR("MSSMA_NU", parameters, nu, unitOpIdx);
		CADET_REGPAR_COMPBND_COMPMAJOR("MSSMA_SIGMA", parameters, sigma, unitOpIdx);
		CADET_REGPAR_COMPBND_COMPMAJOR("MSSMA_RATES", parameters, kRate, unitOpIdx);
		parameters[makeParamId(hashString("EXT_MSSMA_REFC0"), unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &refC0;
		parameters[makeParamId(hashString("EXT_MSSMA_REFQ"), unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &refQ;
	}

	/**
	 * @brief Reserves space in the storage of the parameters
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 * @param [in] totalBoundStates Total number of bound states
	 */
	inline void reserve(unsigned int nComp, unsigned int const* nBoundStates, unsigned int totalBoundStates)
	{
		CADET_RESERVE_SPACE2(kA, totalBoundStates, nComp);
		CADET_RESERVE_SPACE2(kD, totalBoundStates, nComp);
		CADET_RESERVE_SPACE2(nu, totalBoundStates, nComp);
		CADET_RESERVE_SPACE2(sigma, totalBoundStates, nComp);

		unsigned int sumSquared = 0;
		for (unsigned int i = 0; i < nComp; ++i)
			sumSquared += nBoundStates[i] * nBoundStates[i];
		CADET_RESERVE_SPACE2(kRate, sumSquared, nComp);
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
		for (unsigned int i = 0; i < kAT0.size(); ++i)
		{
			CADET_UPDATE_EXTDEP_VARIABLE_NATIVE(kA, i, _extFunBuffer[0]);
			CADET_UPDATE_EXTDEP_VARIABLE_NATIVE(kD, i, _extFunBuffer[1]);
			CADET_UPDATE_EXTDEP_VARIABLE_NATIVE(nu, i, _extFunBuffer[2]);
			CADET_UPDATE_EXTDEP_VARIABLE_NATIVE(sigma, i, _extFunBuffer[3]);
		}

		for (unsigned int i = 0; i < kRateT0.size(); ++i)
			CADET_UPDATE_EXTDEP_VARIABLE_NATIVE(kRate, i, _extFunBuffer[4]);

		CADET_UPDATE_EXTDEP_VARIABLE(lambda, _extFunBuffer[5]);
	}

	CADET_DEFINE_EXTDEP_VARIABLE(active, lambda)
	CADET_DEFINE_EXTDEP_VARIABLE(util::SlicedVector<active>, kA)
	CADET_DEFINE_EXTDEP_VARIABLE(util::SlicedVector<active>, kD)
	CADET_DEFINE_EXTDEP_VARIABLE(util::SlicedVector<active>, nu)
	CADET_DEFINE_EXTDEP_VARIABLE(util::SlicedVector<active>, sigma)
	CADET_DEFINE_EXTDEP_VARIABLE(util::SlicedVector<active>, kRate)
	active refC0; //! Liquid phase reference concentration
	active refQ; //! Solid phase reference concentration
};

/**
 * @brief Defines the Multi-State-SMA binding model
 * @details Implements the Multi-State-SMA adsorption model, which introduces @f$ M_i - 1 @f$ additional bound states @f$ q_i^{(j)} @f$ 
 *          for each component @f$ i @f$ with transitions between them: \f[ \begin{align} 
 *              q_0 &= \Lambda - \sum_{i=1}^{n_{\text{comp}}} \sum_{j=1}^{M_i} \nu_i^{(j)} q_i^{(j)} \\
 *              \frac{\mathrm{d}q_i^{(j)}}{\mathrm{d}t} &= k_{a,i}^{(j)} c_{p,i} \bar{q}_0^{\nu_i^{(j)}} - k_{d,i}^{(j)}\: q_i^{(j)}\: c_{p,0}^{\nu_i^{(j)}} + \sum_{\ell = 1}^{j-1} r^{(i)}_{j\ell} + \sum_{\ell = j+1}^{M} r^{(i)}_{j\ell},
 *          \end{align} \f]
 *          where @f$ i = 1, \dots, N_{\text{comp}} @f$ and @f$ j = 1, \dots, M_i @f$.
 *          The net flux @f$ r^{(i)}_{j\ell} @f$ from @f$ q_i^{(\ell)} @f$ into @f$ q_i^{(j)} @f$ is given by \f[ \begin{align*}
 *               r_{j,\ell}^{(i)} = \begin{cases}
 *	                 k_{\ell j}^{(i)} q_i^{(\ell)} \bar{q}_0^{\nu_i^{(j)} - \nu_i^{(\ell)}} - k_{j \ell}^{(i)} q_i^{(j)} c_{p,0}^{\nu_i^{(j)} - \nu_i^{(\ell)}}, & \ell < j, \\
 *	                 k_{\ell j}^{(i)} q_i^{(\ell)} c_{p,0}^{\nu_i^{(\ell)} - \nu_i^{(j)}} - k_{j \ell}^{(i)} q_i^{(j)} \bar{q}_0^{\nu_i^{(\ell)} - \nu_i^{(j)}}, & \ell > j
 *          \end{cases} \end{align*} \f]
 *          for @f$ 1 \leq \ell, j \leq M_i @f$ and @f$ i = 1, \dots, N_{\text{comp}} @f$.
 *          Here, each component possesses several different binding states with increasing strength: @f[ \nu_i^{(j)} \leq \nu_i^{(j+1)} @f]
 *          A molecule of component @f$ i @f$ can transition from bound state @f$ \ell @f$ to state @f$ j @f$ with given rate @f$ k_{\ell j}^{(i)} @f$.
 *          
 *          Components without bound state (i.e., non-binding components) are supported and all components are allowed to have
 *          varying numbers of bound states (except for salt which has only one bound state @f$ q_0 @f$). Component @c 0 is assumed to be salt.
 *
 *          The internal state vector layout is component-major. The state vector contains all components and within each component
 *          all bound states are placed.
 * @tparam ParamHandler_t Type that can add support for external function dependence
 */
template <class ParamHandler_t>
class MultiStateStericMassActionBindingBase : public BindingModelBase
{
public:

	MultiStateStericMassActionBindingBase() { }
	virtual ~MultiStateStericMassActionBindingBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }
	virtual const char* name() const CADET_NOEXCEPT { return ParamHandler_t::identifier(); }

	virtual void configureModelDiscretization(unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
	{
		BindingModelBase::configureModelDiscretization(nComp, nBound, boundOffset);

		if (nBound[0] != 1)
			throw InvalidParameterException("Multi-State SMA binding model requires exactly one bound state for salt");

		const unsigned int totalBoundStates = numBoundStates(nBound, nComp);

		// Allocate space for parameters
		_p.reserve(nComp, nBound, totalBoundStates);
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
		active* const adRes, active* const adY, unsigned int adEqOffset, unsigned int adDirOffset, const ad::IJacobianExtractor& jacExtractor, 
		double* const workingMemory, linalg::detail::DenseMatrixBase& workingMat) const
	{
		_p.update(t, z, r, secIdx, _nComp, _nBoundStates);

		// If we have kinetic binding, there is only the algebraic salt equation
		if (_kineticBinding)
		{
			// Compute salt component from given bound states q_i^j
			// Salt equation: q_0 - Lambda + Sum[Sum[nu_i^j * q_i^j, j], i] == 0 
			//           <=>  q_0 == Lambda - Sum[Sum[nu_i^j * q_i^j, j], i]
			vecStateY[0] = static_cast<double>(_p.lambda);

			// Loop over all components i
			unsigned int bndIdx = 1;
			for (int i = 1; i < _nComp; ++i)
			{
				// Get nu slice for comp
				active const* const curNu = _p.nu[i];

				// Loop over all bound states j
				for (unsigned int j = 0; j < _nBoundStates[i]; ++j)
				{
					vecStateY[0] -= static_cast<double>(curNu[j]) * vecStateY[bndIdx];

					// Next bound component
					++bndIdx;
				}
			}

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

	CADET_BINDINGMODEL_RESIDUAL_JACOBIAN_BOILERPLATE

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { _p.setExternalFunctions(extFuns, size); }
	virtual bool dependsOnTime() const CADET_NOEXCEPT { return ParamHandler_t::dependsOnTime(); }

	virtual void multiplyWithDerivativeJacobian(double const* yDotS, double* const res, double timeFactor) const
	{
		// Multiplier is 0 if quasi-stationary and 1 if kinetic binding
		const double multiplier = _kineticBinding ? timeFactor : 0.0;

		// First state is salt (always algebraic)
		res[0] = 0.0;

		const unsigned int eqSize = numBoundStates(_nBoundStates, _nComp);
		for (unsigned int i = 1; i < eqSize; ++i)
			res[i] = multiplier * yDotS[i];
	}

	virtual bool hasSalt() const CADET_NOEXCEPT { return true; }
	virtual bool supportsMultistate() const CADET_NOEXCEPT { return true; }
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

		// Salt equation: q_0 - Lambda + Sum[Sum[nu_i^j * q_i^j, j], i] == 0 
		//           <=>  q_0 == Lambda - Sum[Sum[nu_i^j * q_i^j, j], i] 
		// Also compute \bar{q}_0 = q_0 - Sum[Sum[sigma_i^j * q_i^j, j], i]
		res[0] = y[0] - static_cast<ParamType>(_p.lambda);
		ResidualType q0_bar = y[0];

		unsigned int bndIdx = 1;

		// Loop over all components i
		for (int i = 1; i < _nComp; ++i)
		{
			active const* const curNu = _p.nu[i];
			active const* const curSigma = _p.sigma[i];

			// Loop over bound states j of component i
			for (unsigned int j = 0; j < _nBoundStates[i]; ++j)
			{
				res[0] += static_cast<ParamType>(curNu[j]) * y[bndIdx];
				q0_bar -= static_cast<ParamType>(curSigma[j]) * y[bndIdx];

				// Next bound component
				++bndIdx;
			}
		}

		const ParamType refC0 = static_cast<ParamType>(_p.refC0);
		const ParamType refQ = static_cast<ParamType>(_p.refQ);
		const ResidualType yCp0_divRef = yCp[0] / refC0;
		const ResidualType q0_bar_divRef = q0_bar / refQ;

		// Protein equations

		// Loop over all components i
		bndIdx = 1;
		for (int i = 1; i < _nComp; ++i)
		{
			active const* const curKa = _p.kA[i];
			active const* const curKd = _p.kD[i];
			active const* const curNu = _p.nu[i];
			active const* const curRates = _p.kRate[i];
			const unsigned int nBoundStatesI = _nBoundStates[i];

			// Loop over bound states j of component i
			for (unsigned int j = 0; j < nBoundStatesI; ++j)
			{
				const ResidualType c0_pow_nu = pow(yCp0_divRef, static_cast<ParamType>(curNu[j]));
				const ResidualType q0_bar_pow_nu = pow(q0_bar_divRef, static_cast<ParamType>(curNu[j]));

				// Calculate residual
				// Adsorption and desorption
				res[bndIdx] = static_cast<ParamType>(curKd[j]) * y[bndIdx] * c0_pow_nu - static_cast<ParamType>(curKa[j]) * yCp[i] * q0_bar_pow_nu;

				// Conversion to and from weaker states
				for (unsigned int l = 0; l < j; ++l)
				{
					const ParamType nuDiff = static_cast<ParamType>(curNu[j]) - static_cast<ParamType>(curNu[l]);
					// Conversion to weaker state
					res[bndIdx] += static_cast<ParamType>(curRates[nBoundStatesI * j + l]) * y[bndIdx] * pow(yCp0_divRef, nuDiff);
					// Conversion from weaker state
					res[bndIdx] -= static_cast<ParamType>(curRates[nBoundStatesI * l + j]) * y[bndIdx - j + l] * pow(q0_bar_divRef, nuDiff);
					// Getting to q_i^l: bndIdx takes us to q_i^j, subtracting j sets the pointer back
					//                   to q_i^0, and adding l brings us to q_i^l.
					//                   This means y[bndIdx - j + l] corresponds to q_i^l
				}

				// Conversion to and from stronger states
				for (unsigned int l = j+1; l < nBoundStatesI; ++l)
				{
					const ParamType nuDiff = static_cast<ParamType>(curNu[l]) - static_cast<ParamType>(curNu[j]);
					// Conversion to stronger state
					res[bndIdx] += static_cast<ParamType>(curRates[nBoundStatesI * j + l]) * y[bndIdx] * pow(q0_bar_divRef, nuDiff);
					// Conversion from stronger state
					res[bndIdx] -= static_cast<ParamType>(curRates[nBoundStatesI * l + j]) * y[bndIdx - j + l] * pow(yCp0_divRef, nuDiff);
				}

				// Add time derivative if necessary
				if (_kineticBinding && yDot)
				{
					res[bndIdx] += timeFactor * yDot[bndIdx];
				}

				// Next bound component
				++bndIdx;
			}
		}
		return 0;
	}

	template <typename RowIterator>
	void jacobianImpl(double t, double z, double r, unsigned int secIdx, double const* y, double const* yCp, RowIterator jac) const
	{
		_p.update(t, z, r, secIdx, _nComp, _nBoundStates);

		double q0_bar = y[0];

		// Salt equation: q_0 - Lambda + Sum[Sum[nu_j * q_i^j, j], i] == 0
		jac[0] = 1.0;
		int bndIdx = 1;
		for (int i = 1; i < _nComp; ++i)
		{
			active const* const curNu = _p.nu[i];
			active const* const curSigma = _p.sigma[i];

			for (unsigned int j = 0; j < _nBoundStates[i]; ++j)
			{
				jac[bndIdx] = static_cast<double>(curNu[j]);

				// Calculate \bar{q}_0 = q_0 - Sum[Sum[sigma_j * q_i^j, j], i]
				q0_bar -= static_cast<double>(curSigma[j]) * y[bndIdx];

				// Next bound component
				++bndIdx;
			}
		}

		// Advance to protein equations
		++jac;

		const double refC0 = static_cast<double>(_p.refC0);
		const double refQ = static_cast<double>(_p.refQ);
		const double yCp0_divRef = yCp[0] / refC0;
		const double q0_bar_divRef = q0_bar / refQ;

		// Protein equations
		// We have already computed \bar{q}_0 in the loop above
		bndIdx = 1;
		for (int i = 1; i < _nComp; ++i)
		{
			active const* const curKa = _p.kA[i];
			active const* const curKd = _p.kD[i];
			active const* const curNu = _p.nu[i];
			active const* const curRates = _p.kRate[i];
			const double refC0 = static_cast<double>(_p.refC0);
			const double refQ = static_cast<double>(_p.refQ);
			const unsigned int nBoundStatesI = _nBoundStates[i];

			for (int j = 0; j < static_cast<int>(nBoundStatesI); ++j)
			{
				// Getting to c_{p,0}: -bndIdx takes us to q_0, another -nComp to c_{p,0}. This means jac[-bndIdx - nComp] corresponds to c_{p,0}.
				// Getting to c_{p,i}: -bndIdx takes us to q_0, another -nComp to c_{p,0} and a +i to c_{p,i}.
				//                     This means jac[i - bndIdx - nComp] corresponds to c_{p,i}.

				const double ka = static_cast<double>(curKa[j]);
				const double kd = static_cast<double>(curKd[j]);
				const double nu = static_cast<double>(curNu[j]);

				const double c0_pow_nu     = pow(yCp0_divRef, nu);
				const double q0_bar_pow_nu = pow(q0_bar_divRef, nu);
				const double c0_pow_nu_m1_divRef     = pow(yCp0_divRef, nu - 1.0) / refC0;
				const double q0_bar_pow_nu_m1_divRef = nu * pow(q0_bar_divRef, nu - 1.0) / refQ;

				// dres_i / dc_{p,0}
				jac[-bndIdx - _nComp] = kd * y[bndIdx] * nu * c0_pow_nu_m1_divRef;
				// dres_i / dc_{p,i}
				jac[i - bndIdx - _nComp] = -ka * q0_bar_pow_nu;
				// dres_i / dq_0
				jac[-bndIdx] = -ka * yCp[i] * q0_bar_pow_nu_m1_divRef;

				// Fill dres_i / dq_i^j (no flux terms, just handle \bar{q}_0^{nu_i} term)
				int bndIdx2 = 1;
				for (int i2 = 1; i2 < _nComp; ++i2)
				{
					active const* const curSigma = _p.sigma[i2];
					for (unsigned int j2 = 0; j2 < _nBoundStates[i2]; ++j2)
					{
						// dres_i / dq_{i2}^{j2}
						jac[bndIdx2 - bndIdx] = ka * yCp[i] * q0_bar_pow_nu_m1_divRef * static_cast<double>(curSigma[j2]);
						// Getting to q_{i2}^{j2}: -bndIdx takes us to q_1^0, another +bndIdx2 to q_{i2}^{j2}.
						// This means jac[bndIdx2 - bndIdx] corresponds to q_{i2}^{j2}.

						++bndIdx2;
					}
				}

				// Add to dres_i / dq_i
				jac[0] += kd * c0_pow_nu;

				// Handle flux terms
				
				// Conversion to and from weaker states
				for (int l = 0; l < j; ++l)
				{
					const double nuDiff = nu - static_cast<double>(curNu[l]);
					const double q0_bar_pow_nudiff_deriv = static_cast<double>(curRates[nBoundStatesI * l + j]) * y[bndIdx - j + l] * nuDiff * pow(q0_bar_divRef, nuDiff - 1.0) / refQ;

					// dres_i / dc_{p,0}
					jac[-bndIdx - _nComp] += static_cast<double>(curRates[nBoundStatesI * j + l]) * y[bndIdx] * nuDiff * pow(yCp0_divRef, nuDiff - 1.0) / refC0;
					// dres_i / dq_0
					jac[-bndIdx] -= q0_bar_pow_nudiff_deriv;
					// dres_i / dq_i^j (current outer loop element)
					jac[0] += static_cast<double>(curRates[nBoundStatesI * j + l]) * pow(yCp0_divRef, nuDiff);
					// dres_i / dq_i^l (without dq0_bar / dq_i^l)
					jac[l - j] -= static_cast<double>(curRates[nBoundStatesI * l + j]) * pow(q0_bar_divRef, nuDiff);
					// dres_i / dq_{i2}^{j2} (accounts for all dq0_bar / dq_{i2}^{j2} including missing dq0_bar / dq_i^l from above)
					bndIdx2 = 1;
					for (int i2 = 1; i2 < _nComp; ++i2)
					{
						active const* const curSigma = _p.sigma[i2];
						for (int j2 = 0; j2 < static_cast<int>(_nBoundStates[i2]); ++j2, ++bndIdx2)
							jac[bndIdx2 - bndIdx] += q0_bar_pow_nudiff_deriv * static_cast<double>(curSigma[j2]);
					}
				}

				// Conversion to and from stronger states
				for (int l = j+1; l < static_cast<int>(nBoundStatesI); ++l)
				{
					const double nuDiff = static_cast<double>(curNu[l]) - nu;
					const double q0_bar_pow_nudiff_deriv = static_cast<double>(curRates[nBoundStatesI * j + l]) * y[bndIdx] * nuDiff * pow(q0_bar_divRef, nuDiff - 1.0) / refQ;

					// dres_i / dc_{p,0}
					jac[-bndIdx - _nComp] -= static_cast<double>(curRates[nBoundStatesI * l + j]) * y[bndIdx - j + l] * nuDiff * pow(yCp0_divRef, nuDiff - 1.0) / refC0;
					// dres_i / dq_0
					jac[-bndIdx] += q0_bar_pow_nudiff_deriv;
					// dres_i / dq_i^j (current outer loop element)
					jac[0] += static_cast<double>(curRates[nBoundStatesI * j + l]) * pow(q0_bar_divRef, nuDiff);
					// dres_i / dq_i^l
					jac[l - j] -= static_cast<double>(curRates[nBoundStatesI * l + j]) * pow(yCp0_divRef, nuDiff);
					// dres_i / dq_{i2}^{j2} (accounts for all dq0_bar / dq_{i2}^{j2})
					bndIdx2 = 1;
					for (int i2 = 1; i2 < _nComp; ++i2)
					{
						active const* const curSigma = _p.sigma[i2];
						for (int j2 = 0; j2 < static_cast<int>(_nBoundStates[i2]); ++j2, ++bndIdx2)
							jac[bndIdx2 - bndIdx] -= q0_bar_pow_nudiff_deriv * static_cast<double>(curSigma[j2]);
					}
				}

				// Advance to next equation and Jacobian row
				++bndIdx;
				++jac;
			}
		}	
	}

	template <typename RowIterator>
	void jacobianAddDiscretizedImpl(double alpha, RowIterator jac) const
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
};


typedef MultiStateStericMassActionBindingBase<MSSMAParamHandler> MultiStateStericMassActionBinding;
typedef MultiStateStericMassActionBindingBase<ExtMSSMAParamHandler> ExternalMultiStateStericMassActionBinding;

namespace binding
{
	void registerMultiStateStericMassActionModel(std::unordered_map<std::string, std::function<model::IBindingModel*()>>& bindings)
	{
		bindings[MultiStateStericMassActionBinding::identifier()] = []() { return new MultiStateStericMassActionBinding(); };
		bindings[ExternalMultiStateStericMassActionBinding::identifier()] = []() { return new ExternalMultiStateStericMassActionBinding(); };
	}
}  // namespace binding

}  // namespace model

}  // namespace cadet
