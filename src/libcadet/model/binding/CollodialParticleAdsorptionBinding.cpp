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
			{ "type": "ScalarParameter", "varName": "surfaceDensity", "confName": "CPA_SURFACE_DENSITY"},
			{ "type": "ScalarParameter", "varName": "chargeFullLigand", "confName": "CPA_CHARGE_FULL_LIGAND"},
			{ "type": "ScalarParameter", "varName": "pKLigand", "confName": "CPA_PK_LIGAND"},
			{ "type": "ScalarComponentDependentParameter", "varName": "adSurfaceArea", "confName": "CPA_SURFACE_AREA"},
			{ "type": "ScalarComponentDependentParameter", "varName": "compRadius", "confName": "CPA_PROTEIN_RADIUS"},
			{ "type": "ScalarComponentDependentParameter", "varName": "latCharge", "confName": "CPA_COMP_LAT_CHARGE"},
			{ "type": "ScalarComponentDependentParameter", "varName": "refCompCharge", "confName": "CPA_COMP_CHARGE_REF"},
			{ "type": "ScalarComponentDependentParameter", "varName": "linCompCharge", "confName": "CPA_COMP_CHARGE_LIN"},
			{ "type": "ScalarComponentDependentParameter", "varName": "quadCompCharge", "confName": "CPA_COMP_CHARGE_QUAD"},
			{ "type": "ScalarParameter", "varName": "refpH", "confName": "CPA_PH_REF"},
			{ "type": "ScalarComponentDependentParameter", "varName": "refDelta", "confName": "CPA_DELTA_REF"},
			{ "type": "ScalarComponentDependentParameter", "varName": "linDelta", "confName": "CPA_DELTA_LIN"},
			{ "type": "ScalarComponentDependentParameter", "varName": "diffCoeff", "confName": "CPA_DIFFUSION_COEFF"}
		],
	"constantParameters":
		[
			{ "type": "ScalarBoolParameter", "varName": "phIdx", "confName": "CPA_PROTON_IDX"}
		]
}
</codegen>*/

/* Parameter description
 ------------------------
 CPA_TEMPERATURE:         Temperature [K]
 CPA_IONIC_STRENGTH:      Ionic strength [mol/m^3]
 CPA_PERMITTIVITY:        Relative permittivity [-]
 CPA_SURFACE_DENSITY:     Ligand surface density Gamma_L [mol/m^2]
 CPA_CHARGE_FULL_LIGAND:  Charge of fully protonated ligand zeta_L [-]
 CPA_PK_LIGAND:           pK of the ligand [-]
 CPA_PH:                  Bulk pH [-]
 CPA_SURFACE_AREA:        Specific adsorber surface per skeleton volume A_{s,i} [m^-1] (per component)
 CPA_PROTEIN_RADIUS:      Protein radius a_i [m] (per component)
 CPA_COMP_CHARGE_REF:	
 CPA_COMP_CHARGE_LIN:
 CPA_COMP_CHARGE_QUAD:
 CPA_PH_REF:
 CPA_DELTA_REF:
 CPA_DELTA_LIN:
 CPA_COMPONENT_CHARGE: (temporally)
*/

namespace cadet
{

namespace model
{

inline const char* ColloidalParticleAdsorptionParamHandler::identifier() CADET_NOEXCEPT { return "COLLOIDAL_PARTICLE_ADSORPTION"; }

inline bool ColloidalParticleAdsorptionParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_compRadius.size() != _refCompCharge.size())
		|| (_compRadius.size() != _linCompCharge.size())
		|| (_compRadius.size() != _quadCompCharge.size())
		|| (_compRadius.size() != _latCharge.size())
		|| (_compRadius.size() != _adSurfaceArea.size())
		|| (_compRadius.size() != _refDelta.size())
		|| (_compRadius.size() != _linDelta.size())
		|| (_compRadius.size() < nComp))
		throw InvalidParameterException("CPA component-dependent parameters must all have the same size (nComp)");

	return true;
}

inline const char* ExtColloidalParticleAdsorptionParamHandler::identifier() CADET_NOEXCEPT { return "EXT_COLLOIDAL_PARTICLE_ADSORPTION"; }

inline bool ExtColloidalParticleAdsorptionParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_compRadius.size() != _refCompCharge.size())
		|| (_compRadius.size() != _linCompCharge.size())
		|| (_compRadius.size() != _quadCompCharge.size())
		|| (_compRadius.size() != _latCharge.size())
		|| (_compRadius.size() != _adSurfaceArea.size())
		|| (_compRadius.size() != _refDelta.size())
		|| (_compRadius.size() != _linDelta.size())
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
	virtual bool supportsNonBinding() const CADET_NOEXCEPT { return false; }
	virtual bool hasQuasiStationaryReactions() const CADET_NOEXCEPT { return false; }

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
	int _MAXITER = 100; // for newtoniteration in solvePsiAdsorber()
	int _idxProton = 0;

	std::vector<int> _compCharge;
	
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
		for (int iter = 0; iter < _MAXITER; ++iter)
		{
			// pH at surface: pH_0 = pH + (e * psi) / (ln(10) * k_b * T)
			const StateParamType pH0 = pH + (e * psi) / (std::log(10.0) * kb * T);

			// lhs: sigma_{I,A} = e * N_A * Gamma_L * [zeta_L - 1/(1 + 10^{pK_L - pH_0})]
			const StateParamType lhs = e * NA * GammaL * (zetaL - 1.0 / (1.0 + pow(10.0, pKL - pH0)));

			// rhs: sigma_D = 2 * eps * eps0 * kappa * (k_b*T/e) * sinh(e*psi/(2*k_b*T))
			const StateParamType sinArg = e  / (2.0 * kb * T);
			const StateParamType rhs = 2.0 * eps * eps0 * kappa * (kb * T / e) * sinh(sinArg * psi);

			const StateParamType F = lhs - rhs;

			// Derivatives for Newton step
			// dlhs/dpsi: let x = 10^{pKL-pH0}, dx/dpsi = -x*e/(kb*T)
			// dlhs/dpsi = e*NA*GammaL * dx/dpsi / (1+x)^2 = -e^2*NA*GammaL*x / (kb*T*(1+x)^2)
			const StateParamType expTerm = pow(10.0, pKL - pH0);
			const StateParamType dlhs = -e * (e / (kb * T)) * NA * GammaL * expTerm / ((1.0 + expTerm) * (1.0 + expTerm));

			const StateParamType drhs = 2.0 * eps * eps0 * kappa * (kb * T / e) * cosh(sinArg * psi) * sinArg;
			const StateParamType dF = dlhs - drhs;

			if (abs(dF) < 1e-30)
				break;

			const StateParamType delta = -F / dF;
			psi += delta;

			if (abs(delta) < 1e-15 * (1.0 + abs(psi)))
				break;
		}
		return psi;
	}

	virtual bool implementsAnalyticJacobian() const CADET_NOEXCEPT { return true; }

	virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
	{
		const bool valid = ParamHandlerBindingModelBase<ParamHandler_t>::configureImpl(paramProvider, unitOpIdx, parTypeIdx);


		if (_nComp <= 1)
			throw InvalidParameterException("CPA model: To use PH as a state at least two components need to present");

		if(paramProvider.exists("CPA_PROTON_IDX"))// default index is 0
			_idxProton = paramProvider.getInt("CPA_PROTON_IDX");

		if (_nBoundStates[_idxProton] != 0)
		 	throw InvalidParameterException("PH component must be non-binding (NBOUND = 0)");

		for (int i = 0; i < _nComp; ++i)
		{
			if (_nBoundStates[i] > 1)
				throw InvalidParameterException("Binding model supports at most one bound state per component");
		}

		if(paramProvider.exists("CPA_MAXITER"))// default index is 0
			_MAXITER = paramProvider.getInt("CPA_MAXITER");

		if(paramProvider.exists("CPA_COMPONENT_CHARGE"))
		{
			_compCharge = paramProvider.getIntArray("CPA_COMPONENT_CHARGE");
			if (_compCharge.size() != _nComp)
				throw InvalidParameterException("CPA Binding: For every component a charge needs to be provided");

			LOG(Info) << "The definition of component charges is a temporary implementation and will soon be replaced by a general pH modul.";
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

	static double daviesActivityCoeff(double ionicStrength, int charge)
	{
		const double I_M = ionicStrength * 1e-3;
		const double sqrtI = std::sqrt(I_M);
		const double logGamma = -0.509 * charge * charge * (sqrtI / (1.0 + sqrtI) - 0.3 * I_M);
		return std::pow(10.0, logGamma);
	}

	// d(log10(gamma_i))/d(I_m), needed for the activity-corrected pH Jacobian
	static double dDaviesLogGammaDI(double ionicStrength, int charge)
	{
		if (ionicStrength < 1e-30) return 0.0;
		const double I_M = ionicStrength * 1e-3;
		const double sqrtI = std::sqrt(ionicStrength * 1e-3);
		// d/dI [ sqrt(I)/(1+sqrt(I)) - 0.3*I ] = 1/(2*sqrt(I)*(1+sqrt(I))^2) - 0.3
		return -0.509 * charge * charge * (1.0 / (2.0 * sqrtI * (1.0 + sqrtI) * (1.0 + sqrtI)) - 0.3);
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
		const CpStateParamType Im = !_compCharge.empty()
			? static_cast<CpStateParamType>(calcIonicStrength(yCp))
			: static_cast<CpStateParamType>(p->ionicStrength);
		const ParamType eps     = static_cast<ParamType>(p->permittivity);
		const ParamType GammaL  = static_cast<ParamType>(p->surfaceDensity);
		const ParamType zetaL   = static_cast<ParamType>(p->chargeFullLigand);
		const ParamType pKL     = static_cast<ParamType>(p->pKLigand);
		const ParamType refpH   = static_cast<ParamType>(p->refpH);

		// kappa = sqrt(2 * e^2 * I_m * N_A / (k_b * T * eps * eps0))
		const CpStateParamType kappa = e * sqrt(2.0 * Im * NA / (kb * T * eps * eps0));

		const ParamType kbT = kb * T;

		// pH = -log10(a_H+ [mol/L]) = -log10(gamma_H+ * c_H+ * 1e-3)  (c in mol/m^3, factor 1e-3 converts to mol/L)
		const double gammaHp = !_compCharge.empty() ? daviesActivityCoeff(static_cast<double>(Im), _compCharge[_idxProton]) :daviesActivityCoeff(static_cast<double>(Im), 1.0) ;

		const CpStateType pH = _compCharge.empty()
			? -log(yCp[_idxProton] * 1e-3) / log(10.0)
			: -log(gammaHp* yCp[_idxProton] * 1e-3) / log(10.0);

		// Solve adsorber surface potential psi_{0,A}
		const CpStateParamType psiA = solvePsiAdsorber(
			pH, kappa, GammaL, zetaL, pKL, eps, T);

		// beta_{i,j}: e^2 / (4*pi*eps*eps0)
		const ParamType elecPrefactor = e * e / (4.0 * pi * eps * eps0);

		// Theta = pi * N_A * sum_j(a_j^2 * q_j)
		CpStateParamType Theta = 0.0;
		CpStateParamType sumQSurface = 0.0;
		CpStateParamType sumAjQj = 0.0;

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

			++bndIdx;
		}
		Theta = Theta * (pi * NA);

		// D_hex^2 = 2*sqrt(3) / (3 * N_A * sum_j(q_j))
		CpStateParamType Dhex = 0.0;
		if (static_cast<double>(sumQSurface) > 1e-30)
			Dhex = sqrt(2.0 * std::sqrt(3.0) / (3.0 * NA) / sumQSurface);

		bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			if (_nBoundStates[i] == 0)
				continue;

			const ParamType a_i      = static_cast<ParamType>(p->compRadius[i]);
			const ParamType As_i     = static_cast<ParamType>(p->adSurfaceArea[i]);
			const ParamType Zlat_i   = static_cast<ParamType>(p->latCharge[i]);
			const ParamType refZi	 = static_cast<ParamType>(p->refCompCharge[i]);
			const ParamType linZi	 = static_cast<ParamType>(p->linCompCharge[i]);
			const ParamType quadZi   = static_cast<ParamType>(p->quadCompCharge[i]);	

			const CpStateParamType Zi = refZi + linZi * (pH - refpH) + quadZi * (pH - refpH) * (pH - refpH);

			// 1. Protein surface potential psi_{0,i}
			//    psi_i = (2*k_b*T/e) * asinh(Z_i*e^2 / (8*pi*a_i^2*eps*eps0*kappa*k_b*T))
			const CpStateParamType psiArg = Zi * e * e / (8.0 * pi * a_i * a_i * eps * eps0 * kappa * kbT);
			const CpStateParamType psi_i = (2.0 * kbT / e) * log(psiArg + sqrt(psiArg * psiArg + 1.0));

			// 2. Compute delta_m analytically: z^* = -ln( -2*psiA*psi_i / (psiA^2 + psi_i^2))/kappa
			const CpStateParamType dmRatio = -2.0 * psiA * psi_i / (psiA * psiA + psi_i * psi_i);
			const CpStateParamType dm_i = -log(dmRatio) / kappa;

			// 3. Compute delta_i and dstar_i
			const ParamType dRef_i = static_cast<ParamType>(p->refDelta[i]);
			const ParamType dLin_i = static_cast<ParamType>(p->linDelta[i]);
			const CpStateParamType sigmaI_i = Zi * e / (4.0 * pi * a_i * a_i);
			const ParamType sigmaRef_I = refZi * e / (4.0 * pi * a_i * a_i);
			const CpStateParamType logDelta = dRef_i + dLin_i * (abs(sigmaI_i) - abs(sigmaRef_I));
			const CpStateParamType delta_i = exp(std::log(10.0) * logDelta);
			const CpStateParamType dstar_i = dm_i + delta_i / As_i;

			// 4. Protein-adsorber interaction: u_{A,i}(delta_{m,i})
			//    u_{A,i}(z) = pi * a_i * eps * eps0 *
			//      [ 2*psi_A*psi_i * ln((1+exp(-kappa*z))/(1-exp(-kappa*z)))
			//        - (psi_A^2 + psi_i^2) * ln(1 - exp(-2*kappa*z)) ]
			const CpStateParamType ekz = exp(-kappa * dm_i);
			const CpStateParamType uA_i = pi * a_i * eps * eps0 * (
				2.0 * psiA * psi_i * log((1.0 + ekz) / (1.0 - ekz))
				- (psiA * psiA + psi_i * psi_i) * log(1.0 - ekz * ekz)
			);

			// 5. K_{H,i}
			//    K_{H,i} = (k_b*T / u_{A,i}) * (1 - exp(-u_{A,i} / (k_b*T)))
			CpStateParamType KH_i = 0.0;
			if (std::abs(static_cast<double>(uA_i)) > 1e-30)
				KH_i = (kbT / uA_i) * (1.0 - exp(-uA_i / kbT));

			// 6. B_i(Theta)
			//    Hard-disc ASF
			//    B_i = (1 - Theta) * exp(
			//      -(pi*a_i^2 * sum_j(q_j*N_A) + 2*pi*a_i * sum_j(a_j*q_j*N_A)) / (1 - Theta)
			//      - pi*a_i^2 * (sum_j(a_j*q_j*N_A))^2 / (1 - Theta)^2 )

			CpStateParamType B_i = 0.0;
			if (std::abs(static_cast<double>(Theta)) < 1.0)
			{
				const CpStateParamType oneMinusTheta = 1.0 - Theta;
				const CpStateParamType nom1 = pi * a_i * a_i * sumQSurface * NA
					+ 2.0 * pi * a_i * sumAjQj * NA;
				const CpStateParamType nom2 = pi * a_i * a_i
					* (sumAjQj * NA) * (sumAjQj * NA);

				B_i = oneMinusTheta * exp(-nom1 / oneMinusTheta - nom2 / (oneMinusTheta * oneMinusTheta));
			}

			// 7. u_{lat,i} 
			//    u_{lat,i} = 3*sqrt(3)*D_hex*N_A
			//                * exp(-kappa*D_hex) / (1 - exp(-3*sqrt(3)/(2*pi)*kappa*D_hex))
			//                * sum_j(q_j * beta_{i,j})
			//    beta_{i,j} computed inline (Eq. 22)

			CpStateParamType ulat_i = 0.0;

			// sum_j(q_j * beta_{i,j})
			CpStateParamType betaQSum = 0.0;
			int bndIdx2 = 0;
			for (int j = 0; j < _nComp; ++j)
			{
				if (_nBoundStates[j] == 0)
					continue;

				const ParamType Zlat_j = static_cast<ParamType>(p->latCharge[j]);
				const ParamType a_j    = static_cast<ParamType>(p->compRadius[j]);
				const ParamType As_j   = static_cast<ParamType>(p->adSurfaceArea[j]);

				// beta_{i,j} = Zlat_i * Zlat_j * e^2/(4*pi*eps*eps0)
				//              * exp(kappa*(a_i+a_j)) / ((1+kappa*a_i)*(1+kappa*a_j))
				const CpStateParamType beta_ij = Zlat_i * Zlat_j * elecPrefactor
					* exp(kappa * (a_i + a_j))
					/ ((1.0 + kappa * a_i) * (1.0 + kappa * a_j));

				const CpStateParamType q_j = y[bndIdx2] / As_j;
				betaQSum += beta_ij * q_j;

				++bndIdx2;
			}

			const CpStateParamType ulatDenom = 1.0 - exp(-(3.0 * std::sqrt(3.0) / (2.0 * pi)) * kappa * Dhex);

			if (std::abs(static_cast<double>(ulatDenom)) > 1e-30)
			{
				ulat_i = 3.0 * std::sqrt(3.0) * Dhex * NA
					* exp(-kappa * Dhex) / ulatDenom
					* betaQSum;
			}

			// 8. K_{v,i} = As_i * (dstar_i - dm_i) * K_{H,i} * B_i(Theta) * exp(-u_{lat,i} / (k_b*T))
			const CpStateParamType Kv_i = As_i * (dstar_i - dm_i) * KH_i * B_i * exp(-ulat_i / kbT);

			// 9. k_{kin,i} = D_i / (2*Delta^2) * (u_A/(k_bT))^2 / (cosh(u_A/(k_bT)) - 1)
			const ParamType D_i = static_cast<ParamType>(p->diffCoeff[i]); // Typical protein pore diffusion coefficient [m^2/s]
			const CpStateParamType kKin_i_star = D_i / (2.0 * (dstar_i - dm_i) * (dstar_i - dm_i));
			const CpStateParamType uARatio = uA_i / kbT;
			const CpStateParamType kKin_i = kKin_i_star * uARatio * uARatio / (cosh(uARatio) - 1.0);


			res[bndIdx] = kKin_i * (y[bndIdx] - Kv_i * yCp[i]);

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

		if (std::abs(dF_dpsi) < 1e-30)
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

		if (std::abs(dF_dpsi) < 1e-30)
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
		const double pH_val = _compCharge.empty()
			? -std::log10(yCp[_idxProton] * 1e-3)
			: -std::log10(daviesActivityCoeff(Im, _compCharge[_idxProton]) * yCp[_idxProton] * 1e-3);
		const double kappa = sqrt(2.0 * e * e * Im * NA / (kb * T * eps * eps0));
		const double kbT = kb * T;

		// Solve psiA (double version)
		const double psiA = solvePsiAdsorber(pH_val, kappa, GammaL, zetaL, pKL, eps, T);

		// dpsiA / dpH via implicit function theorem
		const double dpsiA_dpH_val = dPsiA_dpH(psiA, pH_val, kappa, GammaL, zetaL, pKL, eps, T);

		const double elecPrefactor = e * e / (4.0 * pi * eps * eps0);

		// Precompute surface concentrations q_j / As_j and auxiliary sums
		const int nTotalBound = std::accumulate(_nBoundStates, _nBoundStates + _nComp, 0);

		// Store per-bound-component data
		std::vector<double> qSurface(nTotalBound, 0.0);  // q_j / As_j
		std::vector<double> aVec(nTotalBound, 0.0);       // a_j (radius)
		std::vector<double> AsVec(nTotalBound, 0.0);       // As_j
		std::vector<int> compIdx(nTotalBound, 0);          // component index for each bound index

		double Theta = 0.0;
		double sumQSurface = 0.0;
		double sumAjQj = 0.0;

		int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			if (_nBoundStates[i] == 0)
				continue;

			const double a_i  = static_cast<double>(p->compRadius[i]);
			const double As_i = static_cast<double>(p->adSurfaceArea[i]);
			const double q_i_surf = y[bndIdx] / As_i;

			qSurface[bndIdx] = q_i_surf;
			aVec[bndIdx] = a_i;
			AsVec[bndIdx] = As_i;
			compIdx[bndIdx] = i;

			Theta       += a_i * a_i * q_i_surf;
			sumQSurface += q_i_surf;
			sumAjQj     += a_i * q_i_surf;

			++bndIdx;
		}
		Theta *= (pi * NA);

		double Dhex = 0.0;
		if (sumQSurface > 1e-30)
			Dhex = sqrt(2.0 * std::sqrt(3.0) / (3.0 * NA) / sumQSurface);

		// --- Main Jacobian loop: one row per bound state ---
		bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			if (_nBoundStates[i] == 0)
				continue;

			const double a_i      = static_cast<double>(p->compRadius[i]);
			const double As_i     = static_cast<double>(p->adSurfaceArea[i]);
			const double refZi    = static_cast<double>(p->refCompCharge[i]);
			const double linZi    = static_cast<double>(p->linCompCharge[i]);
			const double quadZi   = static_cast<double>(p->quadCompCharge[i]);
			const double Zlat_i   = static_cast<double>(p->latCharge[i]);
			const double refpH    = static_cast<double>(p->refpH);

			const double Zi = refZi + linZi * (pH_val - refpH) + quadZi * (pH_val - refpH) * (pH_val - refpH);

			// --- Protein surface potential psi_i ---
			const double psiArg = Zi * e * e / (8.0 * pi * a_i * a_i * eps * eps0 * kappa * kbT);
			const double psi_i = (2.0 * kbT / e) * log(psiArg + sqrt(psiArg * psiArg + 1.0));

			// --- Compute delta_m analytically ---
			const double dmRatio = -2.0 * psiA * psi_i / (psiA * psiA + psi_i * psi_i);
			const double dm_i = -log(dmRatio) / kappa;

			// --- Compute delta_i via Eq. (39) and dstar_i ---
			const double dRef_i = static_cast<double>(p->refDelta[i]);
			const double dLin_i = static_cast<double>(p->linDelta[i]);
			const double sigmaI_i = Zi * e / (4.0 * pi * a_i * a_i);
			const double sigmaRef_I = refZi * e / (4.0 * pi * a_i * a_i);
			const double delta_i = std::pow(10.0, dRef_i + dLin_i * (std::abs(sigmaI_i) - std::abs(sigmaRef_I)));
			const double dstar_i = dm_i + delta_i / As_i;

			// --- Protein-adsorber interaction u_{A,i} ---
			const double ekz = exp(-kappa * dm_i);
			const double logRatio = log((1.0 + ekz) / (1.0 - ekz));
			const double logTerm = log(1.0 - ekz * ekz);

			const double uA_i = pi * a_i * eps * eps0 * (
				2.0 * psiA * psi_i * logRatio
				- (psiA * psiA + psi_i * psi_i) * logTerm
			);

			// Partial derivatives of uA w.r.t. psiA, psi_i, and dm_i
			const double duA_dpsiA = pi * a_i * eps * eps0 * (
				2.0 * psi_i * logRatio - 2.0 * psiA * logTerm
			);
			const double duA_dpsi_i = pi * a_i * eps * eps0 * (
				2.0 * psiA * logRatio - 2.0 * psi_i * logTerm
			);
			const double duA_ddm = pi * a_i * eps * eps0
				* (-2.0 * kappa * ekz / (1.0 - ekz * ekz))
				* (2.0 * psiA * psi_i + (psiA * psiA + psi_i * psi_i) * ekz);

			// --- K_{H,i} ---
			double KH_i = 0.0;
			double dKH_duA = 0.0;
			if (std::abs(uA_i) > 1e-30)
			{
				const double expUA = exp(-uA_i / kbT);
				KH_i = (kbT / uA_i) * (1.0 - expUA);
				// dKH/duA = d/du [ kbT/u * (1 - e^{-u/kbT}) ]
				//         = -kbT/u^2 * (1 - e^{-u/kbT}) + kbT/u * e^{-u/kbT}/kbT
				//         = -KH_i/u + e^{-u/kbT}/u
				dKH_duA = (-KH_i + expUA) / uA_i;
			}

			// --- k_{kin,i} ---
			const double D_i = static_cast<double>(p->diffCoeff[i]);
			const double Delta_i = dstar_i - dm_i;
			const double kKin_star = D_i / (2.0 * Delta_i * Delta_i);
			const double uAratio = uA_i / kbT;
			const double coshUA = cosh(uAratio);
			const double kKin_i = kKin_star * uAratio * uAratio / (coshUA - 1.0);

			// dkKin/duA: let x = uA/kbT, kKin = kKin_star * x^2 / (cosh(x)-1)
			// d/dx [x^2/(cosh(x)-1)] = [2x(cosh(x)-1) - x^2 sinh(x)] / (cosh(x)-1)^2
			double dkKin_duA = 0.0;
			if (std::abs(uA_i) > 1e-30)
			{
				const double sinhUA = sinh(uAratio);
				const double coshM1 = coshUA - 1.0;
				dkKin_duA = kKin_star / kbT * (2.0 * uAratio * coshM1 - uAratio * uAratio * sinhUA) / (coshM1 * coshM1);
			}

			// --- B_i(Theta) ---
			double B_i = 0.0;
			double dBi_dTheta = 0.0;
			double dBi_dsumQ = 0.0;
			double dBi_dsumAjQj = 0.0;

			if (std::abs(Theta) < 1.0)
			{
				const double oneMinusTheta = 1.0 - Theta;
				const double nom1 = pi * a_i * a_i * sumQSurface * NA
					+ 2.0 * pi * a_i * sumAjQj * NA;
				const double nom2 = pi * a_i * a_i * (sumAjQj * NA) * (sumAjQj * NA);

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
				// dnom2/dsumAjQj = 2*pi*a_i^2 * sumAjQj * NA^2
				const double dexpArg_dsumAjQj = -2.0 * pi * a_i * NA / oneMinusTheta
					- 2.0 * pi * a_i * a_i * sumAjQj * NA * NA / (oneMinusTheta * oneMinusTheta);
				dBi_dsumAjQj = B_i * dexpArg_dsumAjQj;
			}

			// --- u_{lat,i} and its derivatives ---
			double ulat_i = 0.0;

			// betaQSum = sum_j(beta_{i,j} * q_j_surf)
			double betaQSum = 0.0;
			std::vector<double> beta_ij_vec(nTotalBound, 0.0);
			int bndIdx2 = 0;
			for (int j = 0; j < _nComp; ++j)
			{
				if (_nBoundStates[j] == 0)
					continue;

				const double Zlat_j = static_cast<double>(p->latCharge[j]);
				const double a_j    = static_cast<double>(p->compRadius[j]);

				const double beta_ij = Zlat_i * Zlat_j * elecPrefactor
					* exp(kappa * (a_i + a_j))
					/ ((1.0 + kappa * a_i) * (1.0 + kappa * a_j));

				beta_ij_vec[bndIdx2] = beta_ij;
				betaQSum += beta_ij * qSurface[bndIdx2];

				++bndIdx2;
			}

			// Dhex-dependent prefactor for u_lat
			const double sqrt3 = std::sqrt(3.0);
			double ulatPrefactor = 0.0;  // 3*sqrt(3)*Dhex*NA * exp(-kappa*Dhex) / denom
			double dulatPrefactor_dDhex = 0.0;

			if (Dhex > 1e-30)
			{
				const double expKD = exp(-kappa * Dhex);
				const double denomArg = 3.0 * sqrt3 / (2.0 * pi) * kappa * Dhex;
				const double denomExp = exp(-denomArg);
				const double denom = 1.0 - denomExp;

				if (std::abs(denom) > 1e-30)
				{
					ulatPrefactor = 3.0 * sqrt3 * Dhex * NA * expKD / denom;
					ulat_i = ulatPrefactor * betaQSum;

					// d(ulatPrefactor)/dDhex
					// let g(D) = 3*sqrt3*D*NA*exp(-k*D) / (1 - exp(-c*k*D)) where c = 3*sqrt3/(2*pi)
					// g'(D) = 3*sqrt3*NA * [exp(-k*D)*(1 - k*D) * denom + D*exp(-k*D)*c*k*exp(-c*k*D)] / denom^2
					// Simplify: g'(D) = 3*sqrt3*NA*exp(-k*D)/denom * [(1-k*D) + D*c*k*exp(-c*k*D)/denom]
					const double c = 3.0 * sqrt3 / (2.0 * pi);
					dulatPrefactor_dDhex = 3.0 * sqrt3 * NA * expKD / denom
						* ((1.0 - kappa * Dhex) + Dhex * c * kappa * denomExp / denom);
				}
			}

			// dDhex/dsumQSurface: Dhex = sqrt(2*sqrt3 / (3*NA*sumQ))
			// dDhex/dsumQ = -0.5 * Dhex / sumQSurface  (if sumQ > 0)
			double dDhex_dsumQ = 0.0;
			if (sumQSurface > 1e-30)
				dDhex_dsumQ = -0.5 * Dhex / sumQSurface;

			// --- Kv_i = As_i * Delta_i * KH_i * B_i * exp(-ulat_i / kbT) ---
			const double expUlat = exp(-ulat_i / kbT);
			const double Kv_i = As_i * Delta_i * KH_i * B_i * expUlat;

			// res_i = kKin_i * (Kv_i * yCp[i] - y[bndIdx])
			// We need:
			//   dres_i / dc_{p,i}    (direct)
			//   dres_i / dc_{p,pH}   (through psiA -> uA -> KH, kKin)
			//   dres_i / dq_j        (through Theta, sumQ, sumAjQj, Dhex -> B_i, ulat_i)
			//   dres_i / dq_i        (direct + through above)

			// dres_i / dc_{p,i}
			// dres_i / dc_{p,i} = -kKin_i * Kv_i
			jac[i - bndIdx - offsetCp] = -kKin_i * Kv_i;

			// dres_i / dc_{p,proton}
			// pH affects: Zi -> psi_i, sigmaI -> delta_i, and psiA (via dPsiA_dpH)
			// These propagate through: dm_i, uA_i, KH_i, kKin_i, Delta_i, Kv_i
			{
				// dZi/dpH
				const double dZi_dpH = linZi + 2.0 * quadZi * (pH_val - refpH);

				// dpsi_i/dpH: psi_i = (2*kbT/e)*arcsinh(Zi*C1), so dpsi_i/dpH = (2*kbT/e)*C1*dZi/dpH / sqrt(arg^2+1)
				const double C1_psi = e * e / (8.0 * pi * a_i * a_i * eps * eps0 * kappa * kbT);
				const double dpsi_i_dpH = (2.0 * kbT / e) * dZi_dpH * C1_psi / sqrt(psiArg * psiArg + 1.0);

				// ddm_i/dpH via R = -2*psiA*psi_i/(psiA^2+psi_i^2)
				const double S_pot = psiA * psiA + psi_i * psi_i;
				const double dR_dpsiA = 2.0 * psi_i * (psiA * psiA - psi_i * psi_i) / (S_pot * S_pot);
				const double dR_dpsi_i = 2.0 * psiA * (psi_i * psi_i - psiA * psiA) / (S_pot * S_pot);
				const double ddm_dpH = (std::abs(dmRatio) > 1e-30)
					? (-1.0 / (dmRatio * kappa)) * (dR_dpsiA * dpsiA_dpH_val + dR_dpsi_i * dpsi_i_dpH)
					: 0.0;

				// ddelta_i/dpH via sigmaI_i
				const double dsigmaI_dpH = dZi_dpH * e / (4.0 * pi * a_i * a_i);
				double ddelta_dpH = 0.0;
				if (std::abs(sigmaI_i) > 1e-30)
					ddelta_dpH = delta_i * std::log(10.0) * dLin_i * (sigmaI_i > 0.0 ? 1.0 : -1.0) * dsigmaI_dpH;
				const double dDelta_dpH = ddelta_dpH / As_i;

				// Total duA_i/dpH = duA/dpsiA * dpsiA/dpH + duA/dpsi_i * dpsi_i/dpH + duA/ddm * ddm/dpH
				const double duA_dpH = duA_dpsiA * dpsiA_dpH_val
					+ duA_dpsi_i * dpsi_i_dpH
					+ duA_ddm * ddm_dpH;

				// dKH/dpH
				const double dKH_dpH = dKH_duA * duA_dpH;

				// dkKin/dpH: kKin depends on Delta_i and uA_i
				const double dkKin_dDelta = (std::abs(Delta_i) > 1e-30) ? -2.0 * kKin_i / Delta_i : 0.0;
				const double dkKin_dpH = dkKin_dDelta * dDelta_dpH + dkKin_duA * duA_dpH;

				// dKv/dpH: B_i and ulat_i are independent of pH
				const double dKv_dpH = As_i * B_i * expUlat * (dDelta_dpH * KH_i + Delta_i * dKH_dpH);

				// dres/dpH
				const double dres_dpH = -(dkKin_dpH * (Kv_i * yCp[i] - y[bndIdx]) + kKin_i * dKv_dpH * yCp[i]);

				// dpH/dc_proton = -1/(c_proton*ln10) - dlog10(gamma_H+)/dIm * dIm/dc_proton
				// (negative signs from pH = -log10(...))
				double dpH_dc_proton = -1.0 / (yCp[_idxProton] * std::log(10.0));
				if (!_compCharge.empty() && std::abs(Im) > 1e-30)
				{
					const double dpH_dIm = -dDaviesLogGammaDI(Im, _compCharge[_idxProton]);
					for (int jj = 0; jj < _nComp; ++jj)
					{
						const double dIm_dcjj = 0.5 * static_cast<double>(_compCharge[jj] * _compCharge[jj]);
						if (std::abs(dIm_dcjj) < 1e-30) continue;
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
			if (!_compCharge.empty() && std::abs(Im) > 1e-30)
			{
				const double dkappa_dIm = kappa / (2.0 * Im);

				// dpsiA/dkappa via IFT
				const double dpsiA_dkappa_val = dPsiA_dkappa(psiA, pH_val, kappa, GammaL, zetaL, pKL, eps, T);

				// dpsi_i/dkappa: psiArg = Zi*e^2/(8*pi*a_i^2*eps*eps0*kappa*kbT)
				// dpsiArg/dkappa = -psiArg / kappa
				const double dpsiArg_dkappa = -psiArg / kappa;
				const double dpsi_i_dkappa = (2.0 * kbT / e) * dpsiArg_dkappa / sqrt(psiArg * psiArg + 1.0);

				// ddm_i/dkappa: dm_i = -log(R)/kappa, R = -2*psiA*psi_i/(psiA^2+psi_i^2)
				// dR/dkappa through psiA and psi_i
				const double S_pot_k = psiA * psiA + psi_i * psi_i;
				const double dR_dpsiA_k = 2.0 * psi_i * (psiA * psiA - psi_i * psi_i) / (S_pot_k * S_pot_k);
				const double dR_dpsi_i_k = 2.0 * psiA * (psi_i * psi_i - psiA * psiA) / (S_pot_k * S_pot_k);
				const double dR_dkappa = dR_dpsiA_k * dpsiA_dkappa_val + dR_dpsi_i_k * dpsi_i_dkappa;
				double ddm_dkappa = 0.0;
				if (std::abs(dmRatio) > 1e-30)
					ddm_dkappa = log(dmRatio) / (kappa * kappa) - dR_dkappa / (dmRatio * kappa);

				// duA/dkappa: through psiA, psi_i, and explicit kappa in ekz
				// ekz = exp(-kappa*dm_i), duA_ddm accounts for d(ekz)/d(dm_i) at constant kappa
				// Extra term for explicit kappa: duA_ddm * dm_i / kappa
				const double duA_dkappa = duA_dpsiA * dpsiA_dkappa_val
					+ duA_dpsi_i * dpsi_i_dkappa
					+ duA_ddm * (ddm_dkappa + dm_i / kappa);

				// dKH/dkappa
				const double dKH_dkappa = dKH_duA * duA_dkappa;

				// dkKin/dkappa: Delta_i = delta_i/As_i doesn't depend on kappa, only uA does
				const double dkKin_dkappa = dkKin_duA * duA_dkappa;

				// dulat_i/dkappa: through beta_ij and ulatPrefactor
				// dbeta_ij/dkappa = beta_ij * [(a_i+a_j) - a_i/(1+kappa*a_i) - a_j/(1+kappa*a_j)]
				double dbetaQSum_dkappa = 0.0;
				int bndIdx3 = 0;
				for (int j = 0; j < _nComp; ++j)
				{
					if (_nBoundStates[j] == 0)
						continue;

					const double a_j = aVec[bndIdx3];
					const double dbeta_dkappa = beta_ij_vec[bndIdx3]
						* ((a_i + a_j) - a_i / (1.0 + kappa * a_i) - a_j / (1.0 + kappa * a_j));
					dbetaQSum_dkappa += dbeta_dkappa * qSurface[bndIdx3];
					++bndIdx3;
				}

				// dulatPrefactor/dkappa: ulatPrefactor = 3*sqrt3*Dhex*NA*exp(-kappa*Dhex) / (1-exp(-c*kappa*Dhex))
				double dulatPrefactor_dkappa = 0.0;
				if (Dhex > 1e-30)
				{
					const double c_lat = 3.0 * sqrt3 / (2.0 * pi);
					const double hLat = exp(-c_lat * kappa * Dhex);
					const double denomLat = 1.0 - hLat;
					if (std::abs(denomLat) > 1e-30)
						dulatPrefactor_dkappa = ulatPrefactor * (-Dhex) * (1.0 + (c_lat - 1.0) * hLat) / denomLat;
				}

				const double dulat_dkappa = dulatPrefactor_dkappa * betaQSum + ulatPrefactor * dbetaQSum_dkappa;

				// dKv/dkappa: Kv = As*Delta*KH*B*exp(-ulat/kbT)
				// B_i and Delta_i don't depend on kappa
				const double dKv_dkappa = As_i * Delta_i * B_i * expUlat
					* (dKH_dkappa - KH_i * dulat_dkappa / kbT);

				// dres/dkappa
				const double dres_dkappa = -(dkKin_dkappa * (Kv_i * yCp[i] - y[bndIdx]) + kKin_i * dKv_dkappa * yCp[i]);

				// Distribute over all charged pore-phase components:
				// dres/dc_j = dres/dkappa * dkappa/dIm * dIm/dc_j, with dIm/dc_j = 0.5 * z_j^2
				for (int j = 0; j < _nComp; ++j)
				{
					const double dIm_dcj = 0.5 * static_cast<double>(_compCharge[j] * _compCharge[j]);
					if (std::abs(dIm_dcj) < 1e-30)
						continue;
					jac[j - bndIdx - offsetCp] += dres_dkappa * dkappa_dIm * dIm_dcj;
				}
			}

			// === dres_i / dq_i (direct term: +kKin_i due to negative sign in res) ===
			jac[0] = +kKin_i;

			// === dres_i / dq_j (through Kv_i which depends on B_i, ulat_i via Theta, sumQ, sumAjQj, Dhex) ===
			// Kv_i = As*Delta*KH * B_i * exp(-ulat/kbT)
			// dKv/dq_j = As*Delta*KH * [dB_i/dq_j * exp(-ulat/kbT) + B_i * exp(-ulat/kbT) * (-1/kbT) * dulat/dq_j]
			//          = Kv_i * [dB_i/dq_j / B_i  -  dulat/dq_j / kbT]   (when B_i != 0)

			// For each q_j (bound state index k):
			//   dTheta/dq_k    = pi * NA * a_k^2 / As_k
			//   dsumQ/dq_k     = 1 / As_k
			//   dsumAjQj/dq_k  = a_k / As_k

			// dB_i/dq_k = dB_i/dTheta * dTheta/dq_k + dB_i/dsumQ * dsumQ/dq_k + dB_i/dsumAjQj * dsumAjQj/dq_k
			// dulat_i/dq_k = ulatPrefactor * beta_{ik} / As_k
			//              + dulatPrefactor/dDhex * dDhex/dsumQ * dsumQ/dq_k * betaQSum

			bndIdx2 = 0;
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
				const double dulat_direct = ulatPrefactor * beta_ij_vec[bndIdx2] / As_k;
				const double dulat_indirect = dulatPrefactor_dDhex * dDhex_dsumQ * dsumQ_dqk * betaQSum;
				const double dulat_dqk = dulat_direct + dulat_indirect;

				// dKv_i / dq_k
				double dKv_dqk = 0.0;
				if (std::abs(B_i) > 1e-30)
					dKv_dqk = Kv_i * (dBi_dqk / B_i - dulat_dqk / kbT);
				else
					dKv_dqk = As_i * Delta_i * KH_i * expUlat * (dBi_dqk - B_i * dulat_dqk / kbT);

				// dres_i / dq_k = -kKin_i * dKv_i/dq_k * c_{p,i}
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
