// SPDX-License-Identifier: AGPL-3.0-or-later
// =================================================================================
//  CADET
//
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Affero General Public
//  License v3.0 (or, at your option, any later version).
// =================================================================================

#include "model/binding/BindingModelBase.hpp"
#include "model/ExternalFunctionSupport.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "model/Parameters.hpp"
#include "LocalVector.hpp"
#include "SimulationTypes.hpp"
#include "MathUtil.hpp"
#include "Logging.hpp"

#include <functional>
#include <unordered_map>
#include <string>
#include <vector>
#include <cmath>
#include <numbers>

using std::numbers::pi;
/*<codegen>
{
	"name": "ColloidalParticleAdsorptionParamHandler",
	"externalName": "ExtColloidalParticleAdsorptionParamHandler",
	"parameters":
		[
			{ "type": "ScalarParameter", "varName": "temperature", "confName": "CPA_TEMPERATURE"},
			{ "type": "ScalarParameter", "varName": "ionicStrength", "confName": "CPA_IONIC_STRENGTH"},
			{ "type": "ScalarParameter", "varName": "permittivity", "confName": "CPA_PERMITTIVITY"},
			{ "type": "ScalarParameter", "varName": "surfaceDensity", "confName": "CPA_LIGAND_DENSITY"},
			{ "type": "ScalarParameter", "varName": "chargeFullLigand", "confName": "CPA_LIGAND_CHARGE_FULL"},
			{ "type": "ScalarParameter", "varName": "pKLigand", "confName": "CPA_LIGAND_PK"},
			{ "type": "ScalarComponentDependentParameter", "varName": "adSurfaceArea", "confName": "CPA_SPECIFIC_SURFACE_AREA"},
			{ "type": "ScalarComponentDependentParameter", "varName": "compRadius", "confName": "CPA_RADIUS"},
			{ "type": "ScalarComponentDependentParameter", "varName": "latCharge", "confName": "CPA_LAT_CHARGE"},
			{ "type": "ScalarParameter", "varName": "refpH", "confName": "CPA_PH_REF"},
			{ "type": "ScalarComponentDependentParameter", "varName": "refDelta", "confName": "CPA_DELTA_REF"},
			{ "type": "ScalarComponentDependentParameter", "varName": "linDelta", "confName": "CPA_DELTA_LIN"},
			{ "type": "ScalarComponentDependentParameter", "varName": "k_kin", "confName": "CPA_KKIN"}
		]
}
</codegen>*/

/* Parameter description
 ------------------------
 CPA_TEMPERATURE:         	Temperature [K]
 CPA_IONIC_STRENGTH:      	Ionic strength [mol/m^3]
 CPA_PERMITTIVITY:        	Relative permittivity [-]
 CPA_LIGAND_DENSITY:     	Ligand surface density Gamma_L [mol/m^2]
 CPA_LIGAND_CHARGE_FULL:  	Charge of fully protonated ligand zeta_L [-]
 CPA_LIGAND_PK:           	pK of the ligand [-]
 CPA_SPECIFIC_SURFACE_AREA: Specific adsorber surface per skeleton volume A_{s,i} [m^-1] (per component)
 CPA_RADIUS:              	Protein radius a_i [m] (per component)
 CPA_EFFECTIVE_CHARGE_COEF:    Protein charge polynomial coefficients z_{i,k}, polynomial-order-row-major [-]
 CPA_LAT_CHARGE:           	Lateral protein charge Z_{lat,i} [-] (per component)
 CPA_PH_REF:               	Reference pH of the protein charge polynomial [-]
 CPA_DELTA_REF:            	Interaction layer parameter Delta_i at pH_ref [-] (per component)
 CPA_DELTA_LIN:            	Slope Delta_{1,i} of log10(Delta_i) w.r.t. |sigma_{I,i}| [m^2/C] (per component)
 CPA_KKIN:                 Kinetic prefactor k^*_{kin,i} [s^-1] (per component)
 CPA_IONIC_VALENCE:        	Ionic valence z_i [-] (per component, optional; enables I_m from the pore phase)
 CPA_PROTON_IDX:           	0-based index of the proton component [-] (optional, default 0)
 CPA_MAXITER:              	Maximum number of Newton iterations for psi_{0,A} [-] (optional, default 100)
*/

namespace cadet
{

namespace model
{

inline const char* ColloidalParticleAdsorptionParamHandler::identifier() CADET_NOEXCEPT { return "COLLOIDAL_PARTICLE_ADSORPTION"; }

inline bool ColloidalParticleAdsorptionParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_compRadius.size() != _latCharge.size())
		|| (_compRadius.size() != _adSurfaceArea.size())
		|| (_compRadius.size() != _refDelta.size())
		|| (_compRadius.size() != _linDelta.size())
		|| (_compRadius.size() != _k_kin.size())
		|| (_compRadius.size() < nComp))
		throw InvalidParameterException("CPA component-dependent parameters must all have the same size (nComp)");

	return true;
}

inline const char* ExtColloidalParticleAdsorptionParamHandler::identifier() CADET_NOEXCEPT { return "EXT_COLLOIDAL_PARTICLE_ADSORPTION"; }

inline bool ExtColloidalParticleAdsorptionParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_compRadius.size() != _latCharge.size())
		|| (_compRadius.size() != _adSurfaceArea.size())
		|| (_compRadius.size() != _refDelta.size())
		|| (_compRadius.size() != _linDelta.size())
		|| (_compRadius.size() != _k_kin.size())
		|| (_compRadius.size() < nComp))
		throw InvalidParameterException("CPA component-dependent parameters must all have the same size (nComp)");

	return true;
}


/**
 * @brief Defines the Colloidal Particle Adsorption (CPA) binding model
 * @details Implements the CPA model by Briskot et al. (2021) for protein adsorption
 *          on ion exchange resins. The model describes adsorption in the linear and
 *          nonlinear regime using:
 *          - Protein-adsorber
 *          - Protein-protein
 *          - Steric blocking via scaled-particle theory (hard-disc ASF)
 *
 *          Kinetic formulation: dq_{v,i}/dt = k_{kin,i} * (K_{v,i} * c_{p,i} - q_{v,i})
 *          Rapid-equilibrium formulation: 0 = q_{v,i} - K_{v,i} * c_{p,i}
 *
 *          Multiple bound states are not supported.
 * @tparam ParamHandler_t Type that can add support for external function dependence
 */
template <class ParamHandler_t>
class ColloidalParticleAdsorptionBindingBase : public ParamHandlerBindingModelBase<ParamHandler_t>
{
public:

	ColloidalParticleAdsorptionBindingBase(){ }
	virtual ~ColloidalParticleAdsorptionBindingBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }

	virtual bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
	{
		const bool res = BindingModelBase::configureModelDiscretization(paramProvider, nComp, nBound, boundOffset);

		return res;
	}

	virtual void timeDerivativeQuasiStationaryFluxes(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yCp, double const* y, double* dResDt, LinearBufferAllocator workSpace) const
	{
		if (!this->hasQuasiStationaryReactions())
			return;

		if (!ParamHandler_t::dependsOnTime())
			return;
	}

	virtual bool hasSalt() const CADET_NOEXCEPT { return false; }
	virtual bool supportsMultistate() const CADET_NOEXCEPT { return false; }
	virtual bool supportsNonBinding() const CADET_NOEXCEPT { return true; }

	CADET_BINDINGMODELBASE_BOILERPLATE

protected:
	using ParamHandlerBindingModelBase<ParamHandler_t>::_paramHandler;
	using ParamHandlerBindingModelBase<ParamHandler_t>::_reactionQuasistationarity;
	using ParamHandlerBindingModelBase<ParamHandler_t>::_nComp;
	using ParamHandlerBindingModelBase<ParamHandler_t>::_nBoundStates;

	// Physical constants
	static constexpr double _elemCharge  = 1.602176634e-19;   // elementaryCharge e [C]
	static constexpr double _avogadroNum = 6.02214076e23;     // avogadro number N_A [1/mol]
	static constexpr double _boltzmann   = 1.380649e-23;      // boltzmann constant k_b [J/K]
	static constexpr double _vacuumPermi = 8.854187818814e-12;  // vacuumPermittivity eps_0 [F/m]
	static constexpr double _daviesFactor = -0.509;
	int _MAXITER = 100; // for newtoniteration in solvePsiAdsorber()
	int _idxProton = 0;

	std::vector<int> _compCharge;
	std::vector<active> _proteinCharge;
	
	/**
	 * @brief Solve for adsorber surface potential psi_{0,A} using Newton*
	 * @details Solves the neutrality condition sigma_{I,A}(psi) = sigma_D(psi)
	 */
	template <typename CpStateType, typename ParamType, typename KappaType>
	typename DoubleActivePromoter<typename DoubleActivePromoter<CpStateType, ParamType>::type, KappaType>::type
	solvePsiAdsorber(CpStateType pH, KappaType kappa, ParamType GammaL,
		ParamType zetaL, ParamType pKL, ParamType eps, ParamType T) const
	{
		using StateParamType = typename DoubleActivePromoter<typename DoubleActivePromoter<CpStateType, ParamType>::type, KappaType>::type;

		const double e = _elemCharge;
		const double kb = _boltzmann;
		const double eps0 = _vacuumPermi;
		const double NA = _avogadroNum;

		StateParamType psi = -0.01; // Initial guess
		bool _converged_psi = false;
		int _iter_psi = 0;
		for (int iter = 0; iter < _MAXITER; ++iter)
		{
			_iter_psi = iter + 1;
			// pH at surface: pH_0 = pH + (e * psi) / (ln(10) * k_b * T)
			const StateParamType pH0 = pH + (e * psi) / (std::log(10.0) * kb * T);
			const StateParamType expTerm = pow(10.0, pKL - pH0);

			// lhs: sigma_{I,A} = e * N_A * Gamma_L * [zeta_L - 1/(1 + 10^{pK_L - pH_0})]
			const StateParamType lhs = e * NA * GammaL * (zetaL - 1.0 / (1.0 + expTerm));

			// rhs: sigma_D = 2 * eps * eps0 * kappa * (k_b*T/e) * sinh(e*psi/(2*k_b*T))
			const StateParamType sinArg = e  / (2.0 * kb * T);
			const StateParamType rhs = 2.0 * eps * eps0 * kappa * (kb * T / e) * sinh(sinArg * psi);

			const StateParamType F = lhs - rhs;

			// Derivatives for Newton step
			// dlhs/dpsi: let x = 10^{pKL-pH0}, dx/dpsi = -x*e/(kb*T)
			// dlhs/dpsi = e*NA*GammaL * dx/dpsi / (1+x)^2 = -e^2*NA*GammaL*x / (kb*T*(1+x)^2)
			const StateParamType dlhs = -e * (e / (kb * T)) * NA * GammaL * expTerm / ((1.0 + expTerm) * (1.0 + expTerm));

			const StateParamType drhs = 2.0 * eps * eps0 * kappa * (kb * T / e) * cosh(sinArg * psi) * sinArg;
			const StateParamType dF = dlhs - drhs;

			if (abs(F) < 1e-14)
			{
				_converged_psi = true;
				break;
			}

			if (abs(dF) < 1e-14)
				break;

			const StateParamType delta = -F / dF;
			psi += delta;

			if (abs(delta) < 1e-14 * (1.0 + abs(psi)))
			{
				_converged_psi = true;
				break;
			}
		}

		if (!_converged_psi)
			LOG(Warning) << "CPA adsorber surface-potential Newton iteration did not converge after " << _iter_psi << " iteration(s); using the last iterate.";

		return psi;
	}

	virtual bool implementsAnalyticJacobian() const CADET_NOEXCEPT { return true; }

	virtual bool requiresWorkspace() const CADET_NOEXCEPT { return true; }

	virtual unsigned int workspaceSize(unsigned int nComp, unsigned int totalNumBoundStates, unsigned int const* nBoundStates) const CADET_NOEXCEPT
	{
		// jacobianImpl needs four scratch arrays of size totalNumBoundStates (q_j / As_j, a_j, As_j, beta_{i,j})
		return ParamHandlerBindingModelBase<ParamHandler_t>::workspaceSize(nComp, totalNumBoundStates, nBoundStates)
			+ 4u * (totalNumBoundStates * sizeof(double) + alignof(double));
	}

	virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
	{
		const bool valid = ParamHandlerBindingModelBase<ParamHandler_t>::configureImpl(paramProvider, unitOpIdx, parTypeIdx);
		readParameterMatrix(_proteinCharge, paramProvider, "CPA_EFFECTIVE_CHARGE_COEF", 1, 1);

		if ((_nComp == 0) || (_proteinCharge.size() < static_cast<std::size_t>(_nComp)) || (_proteinCharge.size() % _nComp != 0))
			throw InvalidParameterException("CPA_EFFECTIVE_CHARGE_COEF must contain complete polynomial-order rows with NCOMP coefficients each");

		const StringHash proteinChargeHash = hashString("CPA_EFFECTIVE_CHARGE_COEF");
		const std::size_t nChargeCoeffs = _proteinCharge.size() / _nComp;
		for (std::size_t coef = 0; coef < nChargeCoeffs; ++coef)
		{
			for (int comp = 0; comp < _nComp; ++comp)
				this->_parameters[makeParamId(proteinChargeHash, unitOpIdx, comp, parTypeIdx, BoundStateIndep, coef, SectionIndep)] = &_proteinCharge[coef * _nComp + comp];
		}

		if (_nComp <= 1)
			throw InvalidParameterException("CPA model: To use PH as a state at least two components need to present");

		_idxProton = 0; // default index is 0
		if(paramProvider.exists("CPA_PROTON_IDX"))
			_idxProton = paramProvider.getInt("CPA_PROTON_IDX");

		if ((_idxProton < 0) || (_idxProton >= _nComp))
			throw InvalidParameterException("CPA binding: PH index must be smaller than the number of components");

		if (_nBoundStates[_idxProton] != 0)
		 	throw InvalidParameterException("PH component must be non-binding (NBOUND = 0)");

		for (int i = 0; i < _nComp; ++i)
		{
			if (_nBoundStates[i] > 1)
				throw InvalidParameterException("Binding model supports at most one bound state per component");
		}

		_MAXITER = 100; // default iteration limit
		if(paramProvider.exists("CPA_MAXITER"))
			_MAXITER = paramProvider.getInt("CPA_MAXITER");

		_compCharge.clear();
		if(paramProvider.exists("CPA_IONIC_VALENCE"))
		{
			_compCharge = paramProvider.getIntArray("CPA_IONIC_VALENCE");
			if (_compCharge.size() != _nComp)
				throw InvalidParameterException("CPA Binding: For every component a charge needs to be provided");

			LOG(Info) << "The definition of component charges is a temporary implementation and will be replaced by a general pH modul in the future.";

			// CPA_TEMPERATURE is given in Kelvin, the Davies correction below is hard-coded for 25 C = 298.15 K.
			// The externally dependent variant does not provide CPA_TEMPERATURE, so the check is skipped there.
			if (paramProvider.exists("CPA_TEMPERATURE"))
			{
				const double temperature = paramProvider.getDouble("CPA_TEMPERATURE");
				if (std::abs(temperature - 298.15) > 1e-10)
					LOG(Warning) << "CPA Binding: The given temperature of " << temperature << " K deviates from 25 C (298.15 K), which is assumed by the Debye-Hueckel coefficient of 0.509 used in the Davies activity correction";
			}
		}
	
		return valid;
	}

	// Returns ionic strength I = 0.5 * sum_i(z_i^2 * c_i), preserving the AD type of yCp
	template <typename CpStateType>
	CpStateType calcIonicStrength(CpStateType const* yCp) const
	{
		CpStateType sum = 0.0;
		for (int i = 0; i < static_cast<int>(_compCharge.size()); ++i)
			sum += yCp[i] * static_cast<double>(_compCharge[i] * _compCharge[i]);
		return 0.5 * sum;
	}

	template <typename CpStateParamType>
	static CpStateParamType daviesActivityCoeff(CpStateParamType ionicStrength, int charge)
	{
		const CpStateParamType ionicStrengthMolar = ionicStrength * 1e-3;
		const CpStateParamType sqrtI = sqrt(ionicStrengthMolar);
		const CpStateParamType logGamma =  _daviesFactor * charge * charge * (sqrtI / (1.0 + sqrtI) - 0.3 * ionicStrengthMolar);
		return pow(10.0, logGamma);
	}

	// d(log10(gamma_i))/d(I_m), needed for the activity-corrected pH Jacobian
	static double dDaviesLogGammaDI(double ionicStrength, int charge)
	{
		if (ionicStrength < 1e-14) return 0.0;
		const double sqrtI = std::sqrt(ionicStrength * 1e-3);
		// d/dI_m [sqrt(I_M)/(1+sqrt(I_M)) - 0.3*I_M], where I_M = 1e-3*I_m
		return _daviesFactor*1e-3 * charge * charge * (1.0 / (2.0 * sqrtI * (1.0 + sqrtI) * (1.0 + sqrtI)) - 0.3);
	}

	template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
	int fluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, StateType const* y,
		CpStateType const* yCp, ResidualType* res, LinearBufferAllocator workSpace) const
	{
		using std::log;

		using CpStateParamType = typename DoubleActivePromoter<CpStateType, ParamType>::type;
		//using StateParamType = typename DoubleActivePromoter<StateType, ParamType>::type;

		typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Physical constants
		const double e    = _elemCharge;
		const double NA   = _avogadroNum;
		const double kb   = _boltzmann;
		const double eps0 = _vacuumPermi;

		// Scalar parameters
		const ParamType T       = static_cast<ParamType>(p->temperature);
		const CpStateParamType Im = !_compCharge.empty()? static_cast<CpStateParamType>(calcIonicStrength(yCp)) : static_cast<CpStateParamType>(p->ionicStrength);
		const ParamType eps     = static_cast<ParamType>(p->permittivity);
		const ParamType GammaL  = static_cast<ParamType>(p->surfaceDensity);
		const ParamType zetaL   = static_cast<ParamType>(p->chargeFullLigand);
		const ParamType pKL     = static_cast<ParamType>(p->pKLigand);
		const ParamType refpH   = static_cast<ParamType>(p->refpH);

		// kappa = sqrt(2 * e^2 * I_m * N_A / (k_b * T * eps * eps0))
		const CpStateParamType kappa = e * sqrt(2.0 * Im * NA / (kb * T * eps * eps0));

		const ParamType kbT = kb * T;

		// pH = -log10(a_H+ [mol/L]) = -log10(gamma_H+ * c_H+ * 1e-3)  (c in mol/m^3, factor 1e-3 converts to mol/L)
		const CpStateParamType gammaHp = !_compCharge.empty() ? daviesActivityCoeff(Im, _compCharge[_idxProton]) : daviesActivityCoeff(Im, 1);
		const CpStateParamType protonActivity = (_compCharge.empty() ? yCp[_idxProton] : gammaHp * yCp[_idxProton]) * 1e-3;
		// Clamp the proton activity at 1e-14 mol/L, which limits the computed pH to 14
		CpStateParamType pH = 14.0;
		if (protonActivity > 1e-14)
			pH = -log(protonActivity) / log(10.0);

		// Solve adsorber surface potential psi_{0,A}
		const CpStateParamType psiA = solvePsiAdsorber(pH, kappa, GammaL, zetaL, pKL, eps, T);

		// beta_{i,j}: e^2 / (4*pi*eps*eps0)
		const ParamType elecPrefactor = e * e / (4.0 * pi * eps * eps0);

		// Theta = pi * N_A * sum_j(a_j^2 * q_j)
		CpStateParamType Theta = 0.0;
		CpStateParamType sumQSurface = 0.0;
		CpStateParamType sumAjQj = 0.0;
		// beta_{i,j} factorizes into g_i * g_j (up to elecPrefactor) with g_i = Zlat_i * exp(kappa*a_i) / (1 + kappa*a_i),
		// so sum_j(beta_{i,j} * q_j) = elecPrefactor * g_i * sum_j(g_j * q_j) and the pairwise loop collapses
		CpStateParamType gQSum = 0.0;

		int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			if (_nBoundStates[i] == 0)
				continue;

			const ParamType a_i  = static_cast<ParamType>(p->compRadius[i]);
			const ParamType As_i = static_cast<ParamType>(p->adSurfaceArea[i]);
			const CpStateParamType q_i = y[bndIdx] / As_i;

			Theta       += a_i * a_i * q_i;
			sumQSurface += q_i;
			sumAjQj     += a_i * q_i;
			gQSum       += static_cast<ParamType>(p->latCharge[i]) * exp(kappa * a_i) / (1.0 + kappa * a_i) * q_i;

			++bndIdx;
		}
		Theta = Theta * (pi * NA);

		// D_hex^2 = 2*sqrt(3) / (3 * N_A * sum_j(q_j))
		CpStateParamType Dhex = 0.0;
		if (static_cast<double>(sumQSurface) > 1e-14)
			Dhex = sqrt(2.0 * std::sqrt(3.0) / (3.0 * NA) / sumQSurface);

		bndIdx = 0;
		const std::size_t nChargeCoeffs = _proteinCharge.size() / _nComp;
		for (int i = 0; i < _nComp; ++i)
		{
			if (_nBoundStates[i] == 0)
				continue;

			const ParamType a_i      = static_cast<ParamType>(p->compRadius[i]);
			const ParamType Zlat_i   = static_cast<ParamType>(p->latCharge[i]);
			const ParamType refZi    = static_cast<ParamType>(_proteinCharge[i]);
			const CpStateParamType pHDifference = pH - refpH;
			CpStateParamType Zi = refZi;
			CpStateParamType pHPower = 1.0;
			for (std::size_t coef = 1; coef < nChargeCoeffs; ++coef)
			{
				pHPower *= pHDifference;
				Zi += static_cast<ParamType>(_proteinCharge[coef * _nComp + i]) * pHPower;
			}

			// 1. Protein surface potential psi_{0,i}
			//    psi_i = (2*k_b*T/e) * asinh(Z_i*e^2 / (8*pi*a_i^2*eps*eps0*kappa*k_b*T))
			//    asinh is evaluated on |psiArg| to avoid the cancellation of log(x + sqrt(x^2+1)) for x << 0
			const CpStateParamType psiArg = Zi * e * e / (8.0 * pi * a_i * a_i * eps * eps0 * kappa * kbT);
			const CpStateParamType absPsiArg = abs(psiArg);
			const CpStateParamType asinhAbs = log(absPsiArg + sqrt(absPsiArg * absPsiArg + 1.0));
			const CpStateParamType psi_i = (2.0 * kbT / e) * (static_cast<double>(psiArg) < 0.0 ? -asinhAbs : asinhAbs);

			// 2. exp(-kappa*delta_{m,i}) at the minimum of u_{A,i}, obtained from d u_{A,i} / dz = 0.
			//    delta_{m,i} itself is never needed: it enters u_{A,i} only through this exponential,
			//    and A_{s,i}*(d*_i - delta_{m,i}) = delta_i by definition of d*_i
			const CpStateParamType ekz = -2.0 * psiA * psi_i / (psiA * psiA + psi_i * psi_i);
			if (ekz < 1e-14)
				throw InvalidParameterException("CPA Binding: While computing delta_m a log(0) would have been calculated, check your parameter settings ");

			// 3. Compute delta_i
			const ParamType dRef_i = static_cast<ParamType>(p->refDelta[i]);
			const ParamType dLin_i = static_cast<ParamType>(p->linDelta[i]);
			const CpStateParamType sigmaI_i = Zi * e / (4.0 * pi * a_i * a_i);
			const ParamType sigmaRef_I = refZi * e / (4.0 * pi * a_i * a_i);
			const CpStateParamType logDelta = log10(dRef_i) + dLin_i * (abs(sigmaI_i) - abs(sigmaRef_I));
			const CpStateParamType delta_i = exp(std::log(10.0) * logDelta);

			// 4. Protein-adsorber interaction: u_{A,i}(delta_{m,i})
			//    u_{A,i}(z) = pi * a_i * eps * eps0 *
			//      [ 2*psi_A*psi_i * ln((1+exp(-kappa*z))/(1-exp(-kappa*z))) - (psi_A^2 + psi_i^2) * ln(1 - exp(-2*kappa*z)) ]
			const CpStateParamType uA_i = pi * a_i * eps * eps0 * (
				2.0 * psiA * psi_i * log((1.0 + ekz) / (1.0 - ekz))
				- (psiA * psiA + psi_i * psi_i) * log(1.0 - ekz * ekz)
			);

			// 5. K_{H,i}
			//    K_{H,i} = (k_b*T / u_{A,i}) * (1 - exp(-u_{A,i} / (k_b*T)))
			CpStateParamType KH_i = 1.0;
			const CpStateParamType uARatio = uA_i / kbT;
			if (std::abs(static_cast<double>(uARatio)) > 1e-14)
				KH_i = (1.0 - exp(-uARatio)) / uARatio;

			// 6. B_i(Theta)
			//    Hard-disc ASF
			//    B_i = (1 - Theta) * exp(
			//      -(pi*a_i^2 * sum_j(q_j*N_A) + 2*pi*a_i * sum_j(a_j*q_j*N_A)) / (1 - Theta)
			//      - pi^2*a_i^2 * (sum_j(a_j*q_j*N_A))^2 / (1 - Theta)^2 )

			CpStateParamType B_i = 0.0;
			if ((Theta < 1.0))
			{
				const CpStateParamType oneMinusTheta = 1.0 - Theta;
				const CpStateParamType nom1 = pi * a_i * a_i * sumQSurface * NA + 2.0 * pi * a_i * sumAjQj * NA;
				const CpStateParamType nom2 = pi * pi * a_i * a_i * (sumAjQj * NA) * (sumAjQj * NA);

				B_i = oneMinusTheta * exp( - nom1 / oneMinusTheta - nom2 / (oneMinusTheta * oneMinusTheta));
			}
			else
				throw InvalidParameterException("CPA Binding: While computing B_i(Theta) Theta must satisfy Theta < 1 check your parameter settings");

			// 7. u_{lat,i}
			//    u_{lat,i} = 3*sqrt(3)*D_hex*N_A
			//                * exp(-kappa*D_hex) / (1 - exp(-3*sqrt(3)/(2*pi)*kappa*D_hex))
			//                * sum_j(q_j * beta_{i,j})
			//    with beta_{i,j} = Zlat_i * Zlat_j * e^2/(4*pi*eps*eps0) * exp(kappa*(a_i+a_j)) / ((1+kappa*a_i)*(1+kappa*a_j))
			//    (Eq. 22), which factorizes into g_i * g_j as described at gQSum above

			CpStateParamType ulat_i = 0.0;

			// sum_j(q_j * beta_{i,j})
			const CpStateParamType g_i = Zlat_i * exp(kappa * a_i) / (1.0 + kappa * a_i);
			const CpStateParamType betaQSum = elecPrefactor * g_i * gQSum;

			const CpStateParamType ulatDenom = 1.0 - exp(-(3.0 * std::sqrt(3.0) / (2.0 * pi)) * kappa * Dhex);

			if (std::abs(static_cast<double>(ulatDenom)) > 1e-14)
			{
				ulat_i = 3.0 * std::sqrt(3.0) * Dhex * NA
					* exp(-kappa * Dhex) / ulatDenom
					* betaQSum;
			}

			// 8. K_{v,i} = As_i * (dstar_i - dm_i) * K_{H,i} * B_i(Theta) * exp(-u_{lat,i} / (k_b*T))
			//    with As_i * (dstar_i - dm_i) = delta_i
			const CpStateParamType Kv_i = delta_i * KH_i * B_i * exp(-ulat_i / kbT);

			if (_reactionQuasistationarity[bndIdx])
				res[bndIdx] = y[bndIdx] - Kv_i * yCp[i];
			else
			{
				// 9. k_{kin,i} = k^*_{kin,i}/2 * (u_A/(k_b*T))^2 / (cosh(u_A/(k_b*T)) - 1)
				const ParamType kKinScale = static_cast<ParamType>(p->k_kin[i]);

				CpStateParamType kKin_i = kKinScale;
				if (std::abs(static_cast<double>(uARatio)) > 1e-14)
					kKin_i = 0.5 * kKinScale * uARatio * uARatio / (cosh(uARatio) - 1.0);
				
				res[bndIdx] = kKin_i * (y[bndIdx] - Kv_i * yCp[i]);
			}

			++bndIdx;
		}

		return 0;
	}

	/**
	 * @brief Compute dpsiA/dpH via implicit function theorem
	 * @details solvePsiAdsorber solves F(psi, pH) = sigma_{I,A}(psi, pH) - sigma_D(psi) = 0.
	 *          By IFT: dpsi/dpH = -(dF/dpH) / (dF/dpsi)
	 */
	double dPsiA_dpH(double psiA, double pH, double kappa, double GammaL,
		double zetaL, double pKL, double eps, double T) const
	{
		const double e = _elemCharge;
		const double kb = _boltzmann;
		const double eps0 = _vacuumPermi;
		const double NA = _avogadroNum;

		const double pH0 = pH + (e * psiA) / (std::log(10.0) * kb * T);
		const double expTerm = std::pow(10.0, pKL - pH0);

		// dF/dpsi: dlhs/dpsi = -e^2*NA*GammaL*x / (kb*T*(1+x)^2), drhs/dpsi = eps*eps0*kappa*cosh(...)
		const double dlhs_dpsi = -e * (e / (kb * T)) * NA * GammaL * expTerm / ((1.0 + expTerm) * (1.0 + expTerm));
		const double sinArg = e / (2.0 * kb * T);
		const double drhs_dpsi = 2.0 * eps * eps0 * kappa * (kb * T / e) * std::cosh(sinArg * psiA) * sinArg;
		const double dF_dpsi = dlhs_dpsi - drhs_dpsi;

		// dF/dpH: dsigma_I/dpH = -e*NA*GammaL*ln10*x/(1+x)^2
		const double dF_dpH = -e * NA * GammaL * std::log(10.0) * expTerm / ((1.0 + expTerm) * (1.0 + expTerm));

		if (std::abs(dF_dpsi) < 1e-14)
			return 0.0;

		return -dF_dpH / dF_dpsi;
	}

	/**
	 * @brief Compute dpsiA/dkappa
	 * @details F(psi, pH, kappa) = sigma_{I,A}(psi, pH) - sigma_D(psi, kappa) = 0.
	 *          By IFT: dpsi/dkappa = -(dF/dkappa) / (dF/dpsi)
	 */
	double dPsiA_dkappa(double psiA, double pH, double kappa, double GammaL,
		double zetaL, double pKL, double eps, double T) const
	{
		const double e = _elemCharge;
		const double kb = _boltzmann;
		const double eps0 = _vacuumPermi;
		const double NA = _avogadroNum;

		const double pH0 = pH + (e * psiA) / (std::log(10.0) * kb * T);
		const double expTerm = std::pow(10.0, pKL - pH0);

		// dF/dpsi (same as in dPsiA_dpH)
		const double dlhs_dpsi = -e * (e / (kb * T)) * NA * GammaL * expTerm / ((1.0 + expTerm) * (1.0 + expTerm));
		const double sinArg = e / (2.0 * kb * T);
		const double drhs_dpsi = 2.0 * eps * eps0 * kappa * (kb * T / e) * std::cosh(sinArg * psiA) * sinArg;
		const double dF_dpsi = dlhs_dpsi - drhs_dpsi;

		// dF/dkappa: sigma_{I,A} is independent of kappa, so dF/dkappa = -d(sigma_D)/dkappa
		// sigma_D = 2*eps*eps0*kappa*(kbT/e)*sinh(e*psi/(2*kbT))
		const double dF_dkappa = -2.0 * eps * eps0 * (kb * T / e) * std::sinh(sinArg * psiA);

		if (std::abs(dF_dpsi) < 1e-14)
			return 0.0;

		return -dF_dkappa / dF_dpsi;
	}

	template <typename RowIterator>
	void jacobianImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* y, double const* yCp, int offsetCp, RowIterator jac, LinearBufferAllocator workSpace) const
	{
		using std::log;
		using std::exp;
		using std::sqrt;
		using std::cosh;
		using std::sinh;

		typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		const double e    = _elemCharge;
		const double NA   = _avogadroNum;
		const double kb   = _boltzmann;
		const double eps0 = _vacuumPermi;

		const double T       = static_cast<double>(p->temperature);
		const double Im      = !_compCharge.empty() ? calcIonicStrength(yCp) : static_cast<double>(p->ionicStrength);
		const double eps     = static_cast<double>(p->permittivity);
		const double GammaL  = static_cast<double>(p->surfaceDensity);
		const double zetaL   = static_cast<double>(p->chargeFullLigand);
		const double pKL     = static_cast<double>(p->pKLigand);

		// pH = -log10(gamma_H+ * c_H+ * 1e-3)  (c in mol/m^3, factor 1e-3 converts to mol/L)
		const double protonActivity = (_compCharge.empty() ? yCp[_idxProton] : daviesActivityCoeff(Im, _compCharge[_idxProton]) * yCp[_idxProton]) * 1e-3;
		// The pH is constant in the clamped branch, so its concentration derivatives vanish there
		const bool isPHClamped = protonActivity <= 1e-14;
		const double pH_val = isPHClamped ? 14.0 : -std::log10(protonActivity);
		const double kappa = sqrt(2.0 * e * e * Im * NA / (kb * T * eps * eps0));
		const double kbT = kb * T;

		// Solve psiA (double version)
		const double psiA = solvePsiAdsorber(pH_val, kappa, GammaL, zetaL, pKL, eps, T);

		// dpsiA / dpH via implicit function theorem
		const double dpsiA_dpH_val = dPsiA_dpH(psiA, pH_val, kappa, GammaL, zetaL, pKL, eps, T);

		const double elecPrefactor = e * e / (4.0 * pi * eps * eps0);

		// Precompute surface concentrations q_j / As_j and auxiliary sums
		const int nTotalBound = std::accumulate(_nBoundStates, _nBoundStates + _nComp, 0);

		// Store per-bound-component data in the work space, this function runs per particle shell and Newton iteration
		BufferedArray<double> qSurfaceArray = workSpace.array<double>(nTotalBound);  // q_j / As_j
		double* const qSurface = static_cast<double*>(qSurfaceArray);
		BufferedArray<double> aArray = workSpace.array<double>(nTotalBound);         // a_j (radius)
		double* const aVec = static_cast<double*>(aArray);
		BufferedArray<double> AsArray = workSpace.array<double>(nTotalBound);        // As_j
		double* const AsVec = static_cast<double*>(AsArray);
		BufferedArray<double> gArray = workSpace.array<double>(nTotalBound);         // g_j of the beta_{i,j} factorization
		double* const gVec = static_cast<double*>(gArray);

		double Theta = 0.0;
		double sumQSurface = 0.0;
		double sumAjQj = 0.0;
		// beta_{i,j} = elecPrefactor * g_i * g_j with g_i = Zlat_i * exp(kappa*a_i) / (1 + kappa*a_i), so that
		// sum_j(beta_{i,j} * q_j) = elecPrefactor * g_i * gQSum and its kappa derivative follows the product rule
		double gQSum = 0.0;
		double dgQSum_dkappa = 0.0;

		int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			if (_nBoundStates[i] == 0)
				continue;

			const double a_i  = static_cast<double>(p->compRadius[i]);
			const double As_i = static_cast<double>(p->adSurfaceArea[i]);
			const double q_i_surf = y[bndIdx] / As_i;
			const double g_i = static_cast<double>(p->latCharge[i]) * exp(kappa * a_i) / (1.0 + kappa * a_i);

			qSurface[bndIdx] = q_i_surf;
			aVec[bndIdx] = a_i;
			AsVec[bndIdx] = As_i;
			gVec[bndIdx] = g_i;

			Theta       += a_i * a_i * q_i_surf;
			sumQSurface += q_i_surf;
			sumAjQj     += a_i * q_i_surf;
			gQSum       += g_i * q_i_surf;
			// dg_i/dkappa = g_i * (a_i - a_i / (1 + kappa*a_i))
			dgQSum_dkappa += g_i * (a_i - a_i / (1.0 + kappa * a_i)) * q_i_surf;

			++bndIdx;
		}
		Theta *= (pi * NA);

		double Dhex = 0.0;
		if (sumQSurface > 1e-14)
			Dhex = sqrt(2.0 * std::sqrt(3.0) / (3.0 * NA) / sumQSurface);

		// --- Main Jacobian loop: one row per bound state ---
		bndIdx = 0;
		const std::size_t nChargeCoeffs = _proteinCharge.size() / _nComp;
		for (int i = 0; i < _nComp; ++i)
		{
			if (_nBoundStates[i] == 0)
				continue;

			const double a_i      = static_cast<double>(p->compRadius[i]);
			const double refZi    = static_cast<double>(_proteinCharge[i]);
			const double refpH    = static_cast<double>(p->refpH);

			const double pHDifference = pH_val - refpH;
			double Zi = refZi;
			double dZi_dpH = 0.0;
			double pHPower = 1.0;
			for (std::size_t coef = 1; coef < nChargeCoeffs; ++coef)
			{
				dZi_dpH += coef * static_cast<double>(_proteinCharge[coef * _nComp + i]) * pHPower;
				pHPower *= pHDifference;
				Zi += static_cast<double>(_proteinCharge[coef * _nComp + i]) * pHPower;
			}

			// --- Protein surface potential psi_i, asinh evaluated on |psiArg| for stability ---
			const double psiArg = Zi * e * e / (8.0 * pi * a_i * a_i * eps * eps0 * kappa * kbT);
			const double psi_i = (2.0 * kbT / e) * std::asinh(psiArg);

			// --- exp(-kappa*delta_{m,i}) at the minimum of u_{A,i} ---
			const double ekz = -2.0 * psiA * psi_i / (psiA * psiA + psi_i * psi_i);

			if (ekz < 1e-14)
				throw InvalidParameterException("CPA Binding: While computing delta_m a log(0) would have been calculated, check your parameter settings ");

			// --- Compute delta_i via Eq. (39) ---
			const double dRef_i = static_cast<double>(p->refDelta[i]);
			const double dLin_i = static_cast<double>(p->linDelta[i]);
			const double sigmaI_i = Zi * e / (4.0 * pi * a_i * a_i);
			const double sigmaRef_I = refZi * e / (4.0 * pi * a_i * a_i);
			const double delta_i = std::pow(10.0, std::log10(dRef_i) + dLin_i * (std::abs(sigmaI_i) - std::abs(sigmaRef_I)));

			// --- Protein-adsorber interaction u_{A,i} ---
			const double logRatio = log((1.0 + ekz) / (1.0 - ekz));
			const double logTerm = log(1.0 - ekz * ekz);

			const double uA_i = pi * a_i * eps * eps0 * (
				2.0 * psiA * psi_i * logRatio
				- (psiA * psiA + psi_i * psi_i) * logTerm
			);

			// Derivatives of uA w.r.t. psiA and psi_i. Since ekz locates the minimum of u_{A,i}, du_{A,i}/d(ekz)
			// vanishes here, so the partials at fixed ekz already are the total derivatives (envelope theorem)
			// and ekz carries no explicit kappa dependence.
			const double duA_dpsiA = pi * a_i * eps * eps0 * (
				2.0 * psi_i * logRatio - 2.0 * psiA * logTerm
			);
			const double duA_dpsi_i = pi * a_i * eps * eps0 * (
				2.0 * psiA * logRatio - 2.0 * psi_i * logTerm
			);

			// --- K_{H,i} ---
			double KH_i = 1.0;
			double dKH_duA = 0.0;
			const double uAratio = uA_i / kbT;

			if (std::abs(uAratio) > 1e-14)
			{
				KH_i = (1.0 - exp(-uAratio))/uAratio;
				// dKH/duA = d/du [ kbT/u * (1 - e^{-u/kbT}) ]
				//         = -kbT/u^2 * (1 - e^{-u/kbT}) + kbT/u * e^{-u/kbT}/kbT
				//         = -KH_i/u + e^{-u/kbT}/u
				dKH_duA = (-KH_i + exp(-uAratio)) / uA_i;
			}

			// --- Kinetic scaling ---
			const bool isQuasiStationary = _reactionQuasistationarity[bndIdx];
			double kKin_i = 1.0;
			double dkKin_duA = 0.0;

			if (!isQuasiStationary)
			{
				const double kKinScale = static_cast<double>(p->k_kin[i]);
				kKin_i = kKinScale;
				if (std::abs(uAratio) > 1e-14)
				{
					const double coshM1 = cosh(uAratio) - 1.0;
					const double sinhUA = sinh(uAratio);
					kKin_i = 0.5 * kKinScale * uAratio * uAratio / coshM1;
					dkKin_duA = 0.5 * kKinScale / kbT * (2.0 * uAratio * coshM1 - uAratio * uAratio * sinhUA) / (coshM1 * coshM1);
				}
			}

			// --- B_i(Theta) ---
			double B_i = 0.0;
			double dBi_dTheta = 0.0;
			double dBi_dsumQ = 0.0;
			double dBi_dsumAjQj = 0.0;

			if ((Theta < 1.0))
			{
				const double oneMinusTheta = 1.0 - Theta;
				const double nom1 = pi * a_i * a_i * sumQSurface * NA + 2.0 * pi * a_i * sumAjQj * NA;
				const double nom2 = pi * pi * a_i * a_i * (sumAjQj * NA) * (sumAjQj * NA);

				const double expArg = -nom1 / oneMinusTheta - nom2 / (oneMinusTheta * oneMinusTheta);
				B_i = oneMinusTheta * exp(expArg);

				// dB_i/dTheta:
				// B_i = (1-T)*exp(f(T)) where f(T) = -nom1/(1-T) - nom2/(1-T)^2
				// dB_i/dT = -exp(f) + (1-T)*exp(f)*f'(T)
				// f'(T) = -nom1/(1-T)^2 - 2*nom2/(1-T)^3
				const double dfTheta = -nom1 / (oneMinusTheta * oneMinusTheta)
					- 2.0 * nom2 / (oneMinusTheta * oneMinusTheta * oneMinusTheta);
				dBi_dTheta = -exp(expArg) + oneMinusTheta * exp(expArg) * dfTheta;

				// dB_i / dsumQSurface (through nom1)
				// dnom1/dsumQ = pi * a_i^2 * NA
				// d(expArg)/dsumQ = -pi*a_i^2*NA / (1-Theta)
				dBi_dsumQ = B_i * (-pi * a_i * a_i * NA / oneMinusTheta);

				// dB_i / dsumAjQj (through nom1 and nom2)
				// dnom1/dsumAjQj = 2*pi*a_i*NA
				// dnom2/dsumAjQj = 2*pi^2*a_i^2 * sumAjQj * NA^2
				const double dexpArg_dsumAjQj = -2.0 * pi * a_i * NA / oneMinusTheta
					- 2.0 * pi * pi * a_i * a_i * sumAjQj * NA * NA / (oneMinusTheta * oneMinusTheta);
				dBi_dsumAjQj = B_i * dexpArg_dsumAjQj;
			}
			else
				throw InvalidParameterException("CPA Binding: While computing B_i(Theta) Theta must satisfy Theta < 1 check your parameter settings");


			// --- u_{lat,i} and its derivatives ---
			double ulat_i = 0.0;

			// betaQSum = sum_j(beta_{i,j} * q_j_surf) = elecPrefactor * g_i * gQSum
			const double g_i = gVec[bndIdx];
			const double betaQSum = elecPrefactor * g_i * gQSum;

			// Dhex-dependent prefactor for u_lat
			const double sqrt3 = std::sqrt(3.0);
			double ulatPrefactor = 0.0;  // 3*sqrt(3)*Dhex*NA * exp(-kappa*Dhex) / denom
			double dulatPrefactor_dDhex = 0.0;

			if (Dhex > 1e-14)
			{
				const double expKD = exp(-kappa * Dhex);
				const double denomArg = 3.0 * sqrt3 / (2.0 * pi) * kappa * Dhex;
				const double denomExp = exp(-denomArg);
				const double denom = 1.0 - denomExp;

				if (std::abs(denom) > 1e-14)
				{
					ulatPrefactor = 3.0 * sqrt3 * Dhex * NA * expKD / denom;
					ulat_i = ulatPrefactor * betaQSum;

					// d(ulatPrefactor)/dDhex
					// let g(D) = 3*sqrt3*D*NA*exp(-k*D) / (1 - exp(-c*k*D)) where c = 3*sqrt3/(2*pi)
					// g'(D) = 3*sqrt3*NA * [exp(-k*D)*(1 - k*D) * denom - D*exp(-k*D)*c*k*exp(-c*k*D)] / denom^2
					// Simplify: g'(D) = 3*sqrt3*NA*exp(-k*D)/denom * [(1-k*D) - D*c*k*exp(-c*k*D)/denom]
					const double c = 3.0 * sqrt3 / (2.0 * pi);
					dulatPrefactor_dDhex = 3.0 * sqrt3 * NA * expKD / denom
						* ((1.0 - kappa * Dhex) - Dhex * c * kappa * denomExp / denom);
				}
			}

			// dDhex/dsumQSurface: Dhex = sqrt(2*sqrt3 / (3*NA*sumQ))
			// dDhex/dsumQ = -0.5 * Dhex / sumQSurface  (if sumQ > 0)
			double dDhex_dsumQ = 0.0;
			if (sumQSurface > 1e-14)
				dDhex_dsumQ = -0.5 * Dhex / sumQSurface;

			// --- Kv_i = delta_i * KH_i * B_i * exp(-ulat_i / kbT), using As_i * (dstar_i - dm_i) = delta_i ---
			const double expUlat = exp(-ulat_i / kbT);
			const double Kv_i = delta_i * KH_i * B_i * expUlat;

			// res_i = scale_i * (y[bndIdx] - Kv_i * yCp[i]), where scale_i is 1 for rapid equilibrium and kKin_i for kinetic binding
			// We need:
			//   dres_i / dc_{p,i}    (direct)
			//   dres_i / dc_{p,pH}   (through psiA -> uA -> KH, kKin)
			//   dres_i / dq_j        (through Theta, sumQ, sumAjQj, Dhex -> B_i, ulat_i)
			//   dres_i / dq_i        (direct + through above)

			// dres_i / dc_{p,i}
			// dres_i / dc_{p,i} = -scale_i * Kv_i
			jac[i - bndIdx - offsetCp] = -kKin_i * Kv_i;

			// dres_i / dc_{p,proton}
			// pH affects: Zi -> psi_i, sigmaI -> delta_i, and psiA (via dPsiA_dpH)
			// These propagate through: uA_i, KH_i, kKin_i, delta_i, Kv_i
			if (!isPHClamped)
			{
				// dpsi_i/dpH: psi_i = (2*kbT/e)*arcsinh(Zi*C1), so dpsi_i/dpH = (2*kbT/e)*C1*dZi/dpH / sqrt(arg^2+1)
				const double C1_psi = e * e / (8.0 * pi * a_i * a_i * eps * eps0 * kappa * kbT);
				const double dpsi_i_dpH = (2.0 * kbT / e) * dZi_dpH * C1_psi / sqrt(psiArg * psiArg + 1.0);

				// ddelta_i/dpH via sigmaI_i
				const double dsigmaI_dpH = dZi_dpH * e / (4.0 * pi * a_i * a_i);
				double ddelta_dpH = 0.0;
				if (std::abs(sigmaI_i) > 1e-14)
					ddelta_dpH = delta_i * std::log(10.0) * dLin_i * (sigmaI_i > 0.0 ? 1.0 : -1.0) * dsigmaI_dpH;

				// Total duA_i/dpH = duA/dpsiA * dpsiA/dpH + duA/dpsi_i * dpsi_i/dpH
				const double duA_dpH = duA_dpsiA * dpsiA_dpH_val
					+ duA_dpsi_i * dpsi_i_dpH;

				// dKH/dpH
				const double dKH_dpH = dKH_duA * duA_dpH;

				// Kinetic scaling only depends on uA_i and is constant for rapid-equilibrium states
				const double dkKin_dpH = dkKin_duA * duA_dpH;

				// dKv/dpH: B_i and ulat_i are independent of pH
				const double dKv_dpH = B_i * expUlat * (ddelta_dpH * KH_i + delta_i * dKH_dpH);

				// dres/dpH
				const double dres_dpH = -(dkKin_dpH * (Kv_i * yCp[i] - y[bndIdx]) + kKin_i * dKv_dpH * yCp[i]);

				// dpH/dc_proton = -1/(c_proton*ln10) - dlog10(gamma_H+)/dIm * dIm/dc_proton
				// (negative signs from pH = -log10(...))
				double dpH_dc_proton = -1.0 / (yCp[_idxProton] * std::log(10.0));
				if (!_compCharge.empty() && std::abs(Im) > 1e-14)
				{
					const double dpH_dIm = -dDaviesLogGammaDI(Im, _compCharge[_idxProton]);
					for (int jj = 0; jj < _nComp; ++jj)
					{
						const double dIm_dcjj = 0.5 * static_cast<double>(_compCharge[jj] * _compCharge[jj]);
						if (std::abs(dIm_dcjj) < 1e-14) continue;
						if (jj == _idxProton)
							dpH_dc_proton += dpH_dIm * dIm_dcjj;
						else
							jac[jj - bndIdx - offsetCp] += dres_dpH * dpH_dIm * dIm_dcjj;
					}
				}
				jac[_idxProton - bndIdx - offsetCp] += dres_dpH * dpH_dc_proton;
			}

			// Im = 0.5 * sum_j(z_j^2 * c_j), so dIm/dc_j = 0.5 * z_j^2
			// dkappa/dc_j = dkappa/dIm * dIm/dc_j = (kappa/(2*Im)) * 0.5 * z_j^2
			if (!_compCharge.empty() && std::abs(Im) > 1e-14)
			{
				const double dkappa_dIm = kappa / (2.0 * Im);

				// dpsiA/dkappa via IFT
				const double dpsiA_dkappa_val = dPsiA_dkappa(psiA, pH_val, kappa, GammaL, zetaL, pKL, eps, T);

				// dpsi_i/dkappa: psiArg = Zi*e^2/(8*pi*a_i^2*eps*eps0*kappa*kbT)
				// dpsiArg/dkappa = -psiArg / kappa
				const double dpsiArg_dkappa = -psiArg / kappa;
				const double dpsi_i_dkappa = (2.0 * kbT / e) * dpsiArg_dkappa / sqrt(psiArg * psiArg + 1.0);

				// duA/dkappa: only through psiA and psi_i, since ekz locates the minimum of u_{A,i}
				const double duA_dkappa = duA_dpsiA * dpsiA_dkappa_val
					+ duA_dpsi_i * dpsi_i_dkappa;

				// dKH/dkappa
				const double dKH_dkappa = dKH_duA * duA_dkappa;

				// dkKin/dkappa: delta_i doesn't depend on kappa, only uA does
				const double dkKin_dkappa = dkKin_duA * duA_dkappa;

				// dulat_i/dkappa: through beta_ij and ulatPrefactor, with
				// d(elecPrefactor*g_i*gQSum)/dkappa = elecPrefactor * (dg_i/dkappa * gQSum + g_i * dgQSum/dkappa)
				const double dg_i_dkappa = g_i * (a_i - a_i / (1.0 + kappa * a_i));
				const double dbetaQSum_dkappa = elecPrefactor * (dg_i_dkappa * gQSum + g_i * dgQSum_dkappa);

				// dulatPrefactor/dkappa: ulatPrefactor = 3*sqrt3*Dhex*NA*exp(-kappa*Dhex) / (1-exp(-c*kappa*Dhex))
				double dulatPrefactor_dkappa = 0.0;
				if (Dhex > 1e-14)
				{
					const double c_lat = 3.0 * sqrt3 / (2.0 * pi);
					const double hLat = exp(-c_lat * kappa * Dhex);
					const double denomLat = 1.0 - hLat;
					if (std::abs(denomLat) > 1e-14)
						dulatPrefactor_dkappa = ulatPrefactor * (-Dhex) * (1.0 + (c_lat - 1.0) * hLat) / denomLat;
				}

				const double dulat_dkappa = dulatPrefactor_dkappa * betaQSum + ulatPrefactor * dbetaQSum_dkappa;

				// dKv/dkappa: Kv = delta_i*KH*B*exp(-ulat/kbT)
				// B_i and delta_i don't depend on kappa
				const double dKv_dkappa = delta_i * B_i * expUlat
					* (dKH_dkappa - KH_i * dulat_dkappa / kbT);

				// dres/dkappa
				const double dres_dkappa = -(dkKin_dkappa * (Kv_i * yCp[i] - y[bndIdx]) + kKin_i * dKv_dkappa * yCp[i]);

				// Distribute over all charged pore-phase components:
				// dres/dc_j = dres/dkappa * dkappa/dIm * dIm/dc_j, with dIm/dc_j = 0.5 * z_j^2
				for (int j = 0; j < _nComp; ++j)
				{
					const double dIm_dcj = 0.5 * static_cast<double>(_compCharge[j] * _compCharge[j]);
					if (std::abs(dIm_dcj) < 1e-14)
						continue;
					jac[j - bndIdx - offsetCp] += dres_dkappa * dkappa_dIm * dIm_dcj;
				}
			}

			// === dres_i / dq_i (direct term from scale_i * y[bndIdx]) ===
			jac[0] = +kKin_i;

			// === dres_i / dq_j (through Kv_i which depends on B_i, ulat_i via Theta, sumQ, sumAjQj, Dhex) ===
			// Kv_i = delta_i*KH * B_i * exp(-ulat/kbT)
			// dKv/dq_j = delta_i*KH * [dB_i/dq_j * exp(-ulat/kbT) + B_i * exp(-ulat/kbT) * (-1/kbT) * dulat/dq_j]
			//          = Kv_i * [dB_i/dq_j / B_i  -  dulat/dq_j / kbT]   (when B_i != 0)

			// For each q_j (bound state index k):
			//   dTheta/dq_k    = pi * NA * a_k^2 / As_k
			//   dsumQ/dq_k     = 1 / As_k
			//   dsumAjQj/dq_k  = a_k / As_k

			// dB_i/dq_k = dB_i/dTheta * dTheta/dq_k + dB_i/dsumQ * dsumQ/dq_k + dB_i/dsumAjQj * dsumAjQj/dq_k
			// dulat_i/dq_k = ulatPrefactor * beta_{ik} / As_k
			//              + dulatPrefactor/dDhex * dDhex/dsumQ * dsumQ/dq_k * betaQSum

			int bndIdx2 = 0;
			for (int j = 0; j < _nComp; ++j)
			{
				if (_nBoundStates[j] == 0)
					continue;

				const double a_k = aVec[bndIdx2];
				const double As_k = AsVec[bndIdx2];

				const double dTheta_dqk = pi * NA * a_k * a_k / As_k;
				const double dsumQ_dqk = 1.0 / As_k;
				const double dsumAjQj_dqk = a_k / As_k;

				// dB_i/dq_k
				const double dBi_dqk = dBi_dTheta * dTheta_dqk
					+ dBi_dsumQ * dsumQ_dqk
					+ dBi_dsumAjQj * dsumAjQj_dqk;

				// dulat_i/dq_k: direct (beta * q term) + indirect (Dhex depends on sumQ)
				const double dulat_direct = ulatPrefactor * elecPrefactor * g_i * gVec[bndIdx2] / As_k;
				const double dulat_indirect = dulatPrefactor_dDhex * dDhex_dsumQ * dsumQ_dqk * betaQSum;
				const double dulat_dqk = dulat_direct + dulat_indirect;

				// dKv_i / dq_k
				double dKv_dqk = 0.0;
				if (std::abs(B_i) > 1e-14)
					dKv_dqk = Kv_i * (dBi_dqk / B_i - dulat_dqk / kbT);
				else
					dKv_dqk = delta_i * KH_i * expUlat * (dBi_dqk - B_i * dulat_dqk / kbT);

				// dres_i / dq_k = -scale_i * dKv_i/dq_k * c_{p,i}
				const double dres_dqk = -kKin_i * dKv_dqk * yCp[i];

				// jac[bndIdx2 - bndIdx] points to q_{bndIdx2}
				jac[bndIdx2 - bndIdx] += dres_dqk;

				++bndIdx2;
			}

			// Advance to next equation
			++bndIdx;
			++jac;
		}
	}
};

typedef ColloidalParticleAdsorptionBindingBase<ColloidalParticleAdsorptionParamHandler> ColloidalParticleAdsorptionBinding;
typedef ColloidalParticleAdsorptionBindingBase<ExtColloidalParticleAdsorptionParamHandler> ExternalColloidalParticleAdsorptionBinding;

namespace binding
{
	void registerColloidalParticleAdsorptionModel(std::unordered_map<std::string, std::function<model::IBindingModel*()>>& bindings)
	{
		bindings[ColloidalParticleAdsorptionBinding::identifier()] = []() { return new ColloidalParticleAdsorptionBinding(); };
		bindings[ExternalColloidalParticleAdsorptionBinding::identifier()] = []() { return new ExternalColloidalParticleAdsorptionBinding(); };
	}
}  // namespace binding

}  // namespace model

}  // namespace cadet
