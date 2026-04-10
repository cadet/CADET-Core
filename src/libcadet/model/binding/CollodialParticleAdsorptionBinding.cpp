// =============================================================================
//  CADET
//
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

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
			{ "type": "ScalarComponentDependentParameter", "varName": "compCharge", "confName": "CPA_COMP_CHARGE"},
			{ "type": "ScalarComponentDependentParameter", "varName": "latCharge", "confName": "CPA_COMP_LAT_CHARGE"},
			{ "type": "ScalarComponentDependentParameter", "varName": "deltaM", "confName": "CPA_DELTA_M"},
			{ "type": "ScalarComponentDependentParameter", "varName": "deltaStar", "confName": "CPA_DELTA_STAR"}
		],
	"constantParameters":
		[
			{ "type": "ScalarBoolParameter", "varName": "phIdx", "confName": "CPA_PROTON_IDX"},
			{ "type": "ScalarBoolParameter", "varName": "isKinetic", "confName": "CPA_IS_KINETIC"}
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
 CPA_COMP_CHARGE:         Net protein charge Z_i [-] (per component)
 CPA_COMP_LAT_CHARGE:     Lateral charge Z_{lat,i} [-] (per component)
 CPA_DELTA_M:             Location of energy minimum delta_{m,i} [m] (per component)
 CPA_DELTA_STAR:          Thickness of interaction boundary layer delta*_i [m] (per component)
 CPA_USE_PH:              Whether pH is a mobile phase component (bool)
*/

namespace cadet
{

namespace model
{

inline const char* ColloidalParticleAdsorptionParamHandler::identifier() CADET_NOEXCEPT { return "COLLOIDAL_PARTICLE_ADSORPTION"; }

inline bool ColloidalParticleAdsorptionParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_compRadius.size() != _compCharge.size())
		|| (_compRadius.size() != _latCharge.size())
		|| (_compRadius.size() != _adSurfaceArea.size())
		|| (_compRadius.size() != _deltaM.size())
		|| (_compRadius.size() != _deltaStar.size())
		|| (_compRadius.size() < nComp))
		throw InvalidParameterException("CPA component-dependent parameters must all have the same size (nComp)");

	return true;
}

inline const char* ExtColloidalParticleAdsorptionParamHandler::identifier() CADET_NOEXCEPT { return "EXT_COLLOIDAL_PARTICLE_ADSORPTION"; }

inline bool ExtColloidalParticleAdsorptionParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_compRadius.size() != _compCharge.size())
		|| (_compRadius.size() != _latCharge.size())
		|| (_compRadius.size() != _adSurfaceArea.size())
		|| (_compRadius.size() != _deltaM.size())
		|| (_compRadius.size() != _deltaStar.size())
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

		// TODO: Compute time derivative of quasi-stationary fluxes
	}

	virtual bool hasSalt() const CADET_NOEXCEPT { return false; }
	virtual bool supportsMultistate() const CADET_NOEXCEPT { return false; }
	virtual bool supportsNonBinding() const CADET_NOEXCEPT { return false; }
	virtual bool hasQuasiStationaryReactions() const CADET_NOEXCEPT { return false; }

	CADET_BINDINGMODEL_RESIDUAL_BOILERPLATE

protected:
	using ParamHandlerBindingModelBase<ParamHandler_t>::_paramHandler;
	using ParamHandlerBindingModelBase<ParamHandler_t>::_reactionQuasistationarity;
	using ParamHandlerBindingModelBase<ParamHandler_t>::_nComp;
	using ParamHandlerBindingModelBase<ParamHandler_t>::_nBoundStates;

	// Physical constants
	//TODO constexpr
	const double _elemCharge 	= 1.602176634e-19;   // elementaryCharge e [C]
	const double _avogadroNum   = 6.02214076e23;     // avogadro number N_A [1/mol]
	const double _boltzmann     = 1.380649e-23;      // boltzmann constant k_b [J/K]
	const double _vacuumPermi 	= 8.8541878128e-12;  // vacuumPermittivity eps_0 [F/m]
	
	const int _MAXITER = 100; // for newtoniteration in solvePsiAdsorber()
	int _idxProd = 0;
	
	/**
	 * @brief Solve for adsorber surface potential psi_{0,A} using Newton
	 * @details Solves the neutrality condition sigma_{I,A}(psi) = sigma_D(psi)
	 */
	template <typename CpStateType, typename ParamType> //TODO psi param type or double
	ParamType solvePsiAdsorber(CpStateType pH, ParamType kappa, ParamType GammaL,
		ParamType zetaL, ParamType pKL, ParamType eps, ParamType T) const
	{
		using StateParamType = typename DoubleActivePromoter<CpStateType, ParamType>::type;

		const double e = _elemCharge;
		const double kb = _boltzmann;
		const double eps0 = _vacuumPermi;
		const double NA = _avogadroNum;

		ParamType psi = -0.01; // Initial guess
		for (int iter = 0; iter < _MAXITER; ++iter)
		{
			// pH at surface: pH_0 = pH + (e * psi) / (ln(10) * k_b * T)
			const StateParamType pH0 = pH + (e * psi) / (std::log(10.0) * kb * T);

			// lhs: sigma_{I,A} = e * N_A * Gamma_L * [zeta_L - 1/(1 + 10^{pK_L - pH_0})]
			const ParamType lhs = e * NA * GammaL * (zetaL - 1.0 / (1.0 + pow(10.0, static_cast<double>(pKL - pH0))));

			// rhs: sigma_D = 2 * eps * eps0 * kappa * (k_b*T/e) * sinh(e*psi/(2*k_b*T))
			const ParamType sinArg = e  / (2.0 * kb * T);
			const ParamType rhs = 2.0 * eps * eps0 * kappa * (kb * T / e) * sinh(sinArg * psi);

			const ParamType F = lhs - rhs;

			// Derivatives for Newton step
			// dlhs/dpsi = -e*NA*GammaL * b * 10^{pKL-pH0} / (1 + 10^{pKL-pH0})^2
			// where b = e/(ln10*kb*T) = d(pH0)/d(psi) * ln10
			const ParamType b = e / (std::log(10.0) * kb * T);
			const double expTerm = std::pow(10.0, static_cast<double>(pKL - pH0));
			const ParamType dlhs = -e * NA * GammaL * b * expTerm / ((1.0 + expTerm) * (1.0 + expTerm));

			const ParamType drhs = 2.0 * eps * eps0 * kappa * (kb * T / e) * cosh(sinArg * psi) * sinArg;
			const ParamType dF = dlhs - drhs;

			if (abs(dF) < 1e-30)
				break;

			const ParamType delta = -F / dF;
			psi += delta;

			if (abs(delta) < 1e-15 * (1.0 + abs(psi)))
				break;
		}
		return psi;
	}

	virtual bool implementsAnalyticJacobian() const CADET_NOEXCEPT { return true; }

	virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
	{
		const bool res = ParamHandlerBindingModelBase<ParamHandler_t>::configureImpl(paramProvider, unitOpIdx, parTypeIdx);

		//TODO ph modelling
		//TODO add is_kinetic


		if (_nComp <= 1)
			throw InvalidParameterException("CPA model: To use PH as a state at least two components need to present");

		if(paramProvider.exists("CPA_PROTON_IDX"))// default index is 0
			_idxProd = paramProvider.getInt("CPA_PROTON_IDX");

		if (_nBoundStates[_idxProd] != 0)
		 	throw InvalidParameterException("PH component must be non-binding (NBOUND = 0)");

		for (int i = 0; i < _nComp; ++i)
		{
			if (_nBoundStates[i] > 1)
				throw InvalidParameterException("Binding model supports at most one bound state per component");
		}

		return res;
	}

	template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
	int fluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, StateType const* y,
		CpStateType const* yCp, ResidualType* res, LinearBufferAllocator workSpace) const
	{
		using std::log;

		using CpStateParamType = typename DoubleActivePromoter<CpStateType, ParamType>::type;
		//using StateParamType = typename DoubleActivePromoter<StateType, ParamType>::type;

		typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Physical constants (double, no AD tracking needed)
		const double e    = _elemCharge;
		const double NA   = _avogadroNum;
		const double kb   = _boltzmann;
		const double eps0 = _vacuumPermi;

		// Scalar parameters (ParamType for AD sensitivity tracking)
		const ParamType T       = static_cast<ParamType>(p->temperature);
		const ParamType Im      = static_cast<ParamType>(p->ionicStrength);
		const ParamType eps     = static_cast<ParamType>(p->permittivity);
		const ParamType GammaL  = static_cast<ParamType>(p->surfaceDensity);
		const ParamType zetaL   = static_cast<ParamType>(p->chargeFullLigand);
		const ParamType pKL     = static_cast<ParamType>(p->pKLigand);

		const CpStateType pH  = log10(yCp[_idxProd]); // pH = log10(c_pH_proxy), c_pH_proxy = 10^pH

		// kappa = sqrt(2 * e^2 * I_m * N_A / (k_b * T * eps * eps0))
		const ParamType kappa = e * sqrt(2.0 * Im * NA / (kb * T * eps * eps0));

		const ParamType kbT = kb * T;

		// Solve adsorber surface potential psi_{0,A}
		const ParamType psiA = solvePsiAdsorber(
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
			const CpStateParamType q_i = y[bndIdx] / As_i;  // BUG FIX: was yCp[bndIdx], must read bound state y

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
			const ParamType Zi       = static_cast<ParamType>(p->compCharge[i]);
			const ParamType Zlat_i   = static_cast<ParamType>(p->latCharge[i]);
			const ParamType dm_i     = static_cast<ParamType>(p->deltaM[i]);
			const ParamType dstar_i  = static_cast<ParamType>(p->deltaStar[i]);

			// 1. Protein surface potential psi_{0,i}
			//    psi_i = (2*k_b*T/e) * asinh(Z_i*e^2 / (8*pi*a_i^2*eps*eps0*kappa*k_b*T))
			const ParamType psiArg = Zi * e * e / (8.0 * pi * a_i * a_i * eps * eps0 * kappa * kbT);
			const ParamType psi_i = (2.0 * kbT / e) * log(psiArg + sqrt(psiArg * psiArg + 1.0));

			// 2. Protein-adsorber interaction: u_{A,i}(delta_{m,i})
			//    u_{A,i}(z) = pi * a_i * eps * eps0 *
			//      [ 2*psi_A*psi_i * ln((1+exp(-kappa*z))/(1-exp(-kappa*z)))
			//        - (psi_A^2 + psi_i^2) * ln(1 - exp(-2*kappa*z)) ]

			const ParamType ekz = exp(-kappa * dm_i);
			const ParamType uA_i = pi * a_i * eps * eps0 * (
				2.0 * psiA * psi_i * log((1.0 + ekz) / (1.0 - ekz))
				- (psiA * psiA + psi_i * psi_i) * log(1.0 - ekz * ekz)
			);

			// 3. K_{H,i} (Eq. 10)
			//    K_{H,i} = (k_b*T / u_{A,i})
			//              * (1 - exp(-u_{A,i} / (k_b*T)))

			ParamType KH_i = 0.0;
			if (std::abs(static_cast<double>(uA_i)) > 1e-30)
				KH_i =(kbT / uA_i) * (1.0 - exp(-uA_i / kbT));

			// 4. B_i(Theta)
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

			// 5. u_{lat,i} 
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
				const ParamType beta_ij = Zlat_i * Zlat_j * elecPrefactor
					* exp(kappa * (a_i + a_j))
					/ ((1.0 + kappa * a_i) * (1.0 + kappa * a_j));

				const CpStateParamType q_j = y[bndIdx2] / As_j;  // BUG FIX: was yCp[bndIdx2], must read bound state y
				betaQSum += beta_ij * q_j;

				++bndIdx2;
			}

			const CpStateParamType denom = 1.0 - exp(-(3.0 * std::sqrt(3.0) / (2.0 * pi)) * kappa * Dhex);

			if (std::abs(static_cast<double>(denom)) > 1e-30)
			{
				ulat_i = 3.0 * std::sqrt(3.0) * Dhex * NA
					* exp(-kappa * Dhex) / denom
					* betaQSum;
			} 

			// 6.  K_{v,i} =  As_i * (dstar_i - dm_i) * K_{H,i} * B_i(Theta) * exp(-u_{lat,i} / (k_b*T))
			const CpStateParamType Kv_i =  As_i * (dstar_i - dm_i) * KH_i * B_i * exp(-ulat_i / kbT);

			// 7. k_{kin,i} =  D_i /Delta^2_i  * 1/2 (u_A,i(detla_m,i)/k_bT)^2 * (cosh(u_A,i(detla_m,i)//k_bT)-1))^{-1}
			const double D_i = 1e-10; // Typical protein pore diffusion coefficient [m^2/s]
			const ParamType kKin_i_star =  D_i / (2 * (dstar_i - dm_i)*(dstar_i - dm_i));
			const ParamType kKin_i = kKin_i_star * (uA_i/kbT)*(uA_i/kbT) * 1/(cosh(uA_i/kbT)-1);
			
			//res: dq/dt - k_{kin,i} * (K_{v,i} * c_{p,i} - q_{v,i})  = 0
			// => flux = k_{kin,i} * (K_{v,i} * c_{p,i} - q_{v,i})
			// For CADET convention: res = -(flux), so add negative sign
			res[bndIdx] = -kKin_i * (Kv_i * yCp[i] - y[bndIdx]);

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
		const double b = e / (std::log(10.0) * kb * T); // dpH0/dpsi

		// dF/dpsi (same as Newton dF)
		const double dlhs_dpsi = -e * NA * GammaL * b * expTerm / ((1.0 + expTerm) * (1.0 + expTerm));
		const double sinArg = e / (2.0 * kb * T);
		const double drhs_dpsi = 2.0 * eps * eps0 * kappa * (kb * T / e) * std::cosh(sinArg * psiA) * sinArg;
		const double dF_dpsi = dlhs_dpsi - drhs_dpsi;

		// dF/dpH: sigma_{I,A} depends on pH through pH0 = pH + ...
		// dsigma_{I,A}/dpH = e * NA * GammaL * ln(10) * 10^{pKL-pH0} / (1 + 10^{pKL-pH0})^2
		// (positive because decreasing pH0 = increasing 10^{pK-pH0} contribution)
		const double dF_dpH = e * NA * GammaL * std::log(10.0) * expTerm / ((1.0 + expTerm) * (1.0 + expTerm));

		if (std::abs(dF_dpsi) < 1e-30)
			return 0.0;

		return -dF_dpH / dF_dpsi;
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
		const double Im      = static_cast<double>(p->ionicStrength);
		const double eps     = static_cast<double>(p->permittivity);
		const double GammaL  = static_cast<double>(p->surfaceDensity);
		const double zetaL   = static_cast<double>(p->chargeFullLigand);
		const double pKL     = static_cast<double>(p->pKLigand);

		const double pH_val = std::log10(yCp[_idxProd]);  // BUG FIX: must apply log10 like in fluxImpl
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
			const double Zi       = static_cast<double>(p->compCharge[i]);
			const double Zlat_i   = static_cast<double>(p->latCharge[i]);
			const double dm_i     = static_cast<double>(p->deltaM[i]);
			const double dstar_i  = static_cast<double>(p->deltaStar[i]);

			// --- Protein surface potential psi_i (depends only on parameters, not state) ---
			const double psiArg = Zi * e * e / (8.0 * pi * a_i * a_i * eps * eps0 * kappa * kbT);
			const double psi_i = (2.0 * kbT / e) * log(psiArg + sqrt(psiArg * psiArg + 1.0));

			// --- Protein-adsorber interaction u_{A,i} ---
			const double ekz = exp(-kappa * dm_i);
			const double logRatio = log((1.0 + ekz) / (1.0 - ekz));
			const double logTerm = log(1.0 - ekz * ekz);

			const double uA_i = pi * a_i * eps * eps0 * (
				2.0 * psiA * psi_i * logRatio
				- (psiA * psiA + psi_i * psi_i) * logTerm
			);

			// du_{A,i}/dpsiA = pi * a_i * eps * eps0 * (2*psi_i*logRatio - 2*psiA*logTerm)
			const double duA_dpsiA = pi * a_i * eps * eps0 * (
				2.0 * psi_i * logRatio - 2.0 * psiA * logTerm
			);

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
			const double D_i = 1e-10; // Typical protein pore diffusion coefficient [m^2/s]
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

			// === dres_i / dc_{p,i} ===
			// dres_i / dc_{p,i} = -kKin_i * Kv_i
			jac[i - bndIdx - offsetCp] = -kKin_i * Kv_i;

			// === dres_i / dc_{p,pH} (pH is at index _idxpH in yCp) ===
			// Chain: res depends on psiA through uA_i, which affects KH_i and kKin_i
			// dres/dpH = (dres/duA * duA/dpsiA + dKv/dpsiA_via_uA) * dpsiA/dpH
			//
			// dKv/duA = As*Delta * dKH/duA * B * exp(-ulat/kbT)
			// dres/duA = -(dkKin/duA * (Kv*cp - q) + kKin * dKv/duA * cp)
			const double dKv_duA = As_i * Delta_i * dKH_duA * B_i * expUlat;
			const double dres_duA = -(dkKin_duA * (Kv_i * yCp[i] - y[bndIdx])
				+ kKin_i * dKv_duA * yCp[i]);
			const double dres_dpH = dres_duA * duA_dpsiA * dpsiA_dpH_val;

			jac[_idxProd - bndIdx - offsetCp] += dres_dpH;

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
