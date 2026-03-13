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
			{ "type": "ScalarParameter", "varName": "pH", "confName": "CPA_PH"},
			{ "type": "ScalarComponentDependentParameter", "varName": "adSurfaceArea", "confName": "CPA_SURFACE_AREA"},
			{ "type": "ScalarComponentDependentParameter", "varName": "compRadius", "confName": "CPA_PROTEIN_RADIUS"},
			{ "type": "ScalarComponentDependentParameter", "varName": "compCharge", "confName": "CPA_COMP_CHARGE"},
			{ "type": "ScalarComponentDependentParameter", "varName": "latCharge", "confName": "CPA_COMP_LAT_CHARGE"},
			{ "type": "ScalarComponentDependentParameter", "varName": "deltaM", "confName": "CPA_DELTA_M"},
			{ "type": "ScalarComponentDependentParameter", "varName": "deltaStar", "confName": "CPA_DELTA_STAR"},
			{ "type": "ScalarComponentDependentParameter", "varName": "kKin", "confName": "CPA_KKIN"}
		],
	"constantParameters":
		[
			{ "type": "ScalarBoolParameter", "varName": "usePh", "confName": "CPA_USE_PH"}
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
 CPA_KKIN:                Kinetic rate constant k_{kin,i} [1/s] (per component)
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
		|| (_compRadius.size() != _kKin.size())
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
		|| (_compRadius.size() != _kKin.size())
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

	ColloidalParticleAdsorptionBindingBase() : _startIdx(1) { }
	virtual ~ColloidalParticleAdsorptionBindingBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }

	virtual void timeDerivativeQuasiStationaryFluxes(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yCp, double const* y, double* dResDt, LinearBufferAllocator workSpace) const
	{
		if (!this->hasQuasiStationaryReactions())
			return;

		if (!ParamHandler_t::dependsOnTime())
			return;

		// TODO: Compute time derivative of quasi-stationary fluxes
	}

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
	
	int _startIdx;
	int _MAXITER = 100;

	
	/**
	 * @brief Solve for adsorber surface potential psi_{0,A} using Newton-Raphson
	 * @details Solves the neutrality condition sigma_{I,A}(psi) = sigma_D(psi)
	 *          
	 */
	inline double solvePsiAdsorber(double pH, double kappa, double GammaL,
		double zetaL, double pKL, double eps, double T) const
	{
		const double e = _elemCharge;
		const double kb = _boltzmann;
		const double eps0 = _vacuumPermi;
		const double NA = _avogadroNum;

		double psi = -0.01; // Initial guess
		for (int iter = 0; iter < _MAXITER; ++iter)
		{
			// pH at surface: pH_0 = pH + (e * psi) / (ln(10) * k_b * T)
			const double pH0 = pH + (e * psi) / (std::log(10.0) * kb * T);

			// lhs: sigma_{I,A} = e * N_A * Gamma_L * [zeta_L - 1/(1 + 10^{pK_L - pH_0})]
			const double lhs = e * NA * GammaL * (zetaL - 1.0 / (1.0 + std::pow(10.0, pKL- pH0)));

			// rhs: sigma_D = 2 * eps * eps0 * kappa * (k_b*T/e) * sinh(e*psi/(2*k_b*T))
			const double sinArg = e  / (2.0 * kb * T);
			const double rhs = 2.0 * eps * eps0 * kappa * (kb * T / e) * std::sinh( sinArg * psi);

			const double F = lhs - rhs;

			// Derivatives for Newton step
			// dlhs/dpsi = -e*NA*GammaL * b * 10^{pKL-pH0} / (1 + 10^{pKL-pH0})^2
			// where b = e/(ln10*kb*T) = d(pH0)/d(psi) * ln10
			const double b = e / (std::log(10.0) * kb * T);
			const double expTerm = std::pow(10.0, pKL - pH0);
			const double dlhs = -e * NA * GammaL * b * expTerm / ((1.0 + expTerm) * (1.0 + expTerm));

			const double drhs = 2.0 * eps * eps0 * kappa * (kb * T / e) * std::cosh(sinArg * psi) * sinArg;
			const double dF = dlhs - drhs;

			if (std::abs(dF) < 1e-30)
				break;

			const double delta = -F / dF;
			psi += delta;

			if (std::abs(delta) < 1e-15 * (1.0 + std::abs(psi)))
				break;
		}
		return psi;
	}

	virtual bool implementsAnalyticJacobian() const CADET_NOEXCEPT { return false; }

	virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
	{
		const bool res = ParamHandlerBindingModelBase<ParamHandler_t>::configureImpl(paramProvider, unitOpIdx, parTypeIdx);

		//TODO ph modelling
		//TODO add is_kinetic
		// if (_nComp <= 1) 
		// 	throw InvalidParameterException("No protein component present");

		// if (_nBoundStates[0] != 0)
		// 	throw InvalidParameterException("Salt component (index 0) must be non-binding (NBOUND = 0)");

		// if (_paramHandler.usePh().get() && (_nComp <= 2))
		// 	throw InvalidParameterException("No protein component present (existing two components are salt and PH)");

		// if (_paramHandler.usePh().get()) 
		// {
		// 	_startIdx = 2;
		// 	if (_nBoundStates[1] != 0)
		// 		throw InvalidParameterException("PH pseudocomponent (index 1) must be non-binding (NBOUND = 0)");
		// }

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

		//using CpStateParamType = typename DoubleActivePromoter<CpStateType, ParamType>::type;
		using StateParamType = typename DoubleActivePromoter<StateType, ParamType>::type;

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
		const ParamType pH_val  = static_cast<ParamType>(p->pH);

		// kappa = sqrt(2 * e^2 * I_m * N_A / (k_b * T * eps * eps0))
		// Eq. (13)
		const ParamType kappa = sqrt(2.0 * e * e * Im * NA
			/ (kb * T * eps * eps0));

		const ParamType kbT = kb * T;

		// Solve adsorber surface potential psi_{0,A} (Eqs. 12, 15, 16, 17)
		const double psiA = solvePsiAdsorber(
			static_cast<double>(pH_val), static_cast<double>(kappa),
			static_cast<double>(GammaL), static_cast<double>(zetaL),
			static_cast<double>(pKL), static_cast<double>(eps),
			static_cast<double>(T));

		// beta_{i,j} (Eq. 22): e^2 / (4*pi*eps*eps0)
		const ParamType elecPrefactor = e * e / (4.0 * pi * eps * eps0);

		// Theta = pi * N_A * sum_j(a_j^2 * q_j)
		StateParamType Theta = 0.0;
		StateParamType sumQSurface = 0.0;
		StateParamType sumAjQj = 0.0;

		int bndIdx = 0;
		for (int i = _startIdx; i < _nComp; ++i)
		{
			if (_nBoundStates[i] == 0)
				continue;

			const ParamType a_i  = static_cast<ParamType>(p->compRadius[i]);
			const ParamType As_i = static_cast<ParamType>(p->adSurfaceArea[i]);
			const StateParamType q_i = y[bndIdx] / As_i;

			Theta       += a_i * a_i * q_i;
			sumQSurface += q_i;
			sumAjQj     += a_i * q_i;

			++bndIdx;
		}
		Theta = Theta * (pi * NA);

		// D_hex^2 = 2*sqrt(3) / (3 * N_A * sum_j(q_j))
		StateParamType Dhex = 0.0;
		if (static_cast<double>(sumQSurface) > 1e-30)
			Dhex = sqrt(2.0 * std::sqrt(3.0) / (3.0 * NA) / sumQSurface);

		bndIdx = 0;
		for (int i = _startIdx; i < _nComp; ++i)
		{
			if (_nBoundStates[i] == 0)
				continue;

			const ParamType a_i      = static_cast<ParamType>(p->compRadius[i]);
			const ParamType As_i     = static_cast<ParamType>(p->adSurfaceArea[i]);
			const ParamType Zi       = static_cast<ParamType>(p->compCharge[i]);
			const ParamType Zlat_i   = static_cast<ParamType>(p->latCharge[i]);
			const ParamType dm_i     = static_cast<ParamType>(p->deltaM[i]);
			const ParamType dstar_i  = static_cast<ParamType>(p->deltaStar[i]);
			const ParamType kKin_i   = static_cast<ParamType>(p->kKin[i]);

			// 1. Protein surface potential psi_{0,i} (Eq. 14)
			//    psi_i = (2*k_b*T/e) * asinh(Z_i*e^2 / (8*pi*a_i^2*eps*eps0*kappa*k_b*T))
			const ParamType psiArg = Zi * e * e / (8.0 * pi * a_i * a_i * eps * eps0 * kappa * kbT);
			const ParamType psi_i = (2.0 * kbT / e) * log(psiArg + sqrt(psiArg * psiArg + 1.0));

			// 2. Protein-adsorber interaction: u_{A,i}(delta_{m,i}) (Eq. 18, 19)
			//    u_{A,i}(z) = pi * a_i * eps * eps0 *
			//      [ 2*psi_A*psi_i * ln((1+exp(-kappa*z))/(1-exp(-kappa*z)))
			//        - (psi_A^2 + psi_i^2) * ln(1 - exp(-2*kappa*z)) ]

			const ParamType ekz = exp(-kappa * dm_i);
			const ParamType uA_i = pi * a_i * eps * eps0 * (
				2.0 * psiA * psi_i * log((1.0 + ekz) / (1.0 - ekz))
				- (psiA * psiA + psi_i * psi_i) * log(1.0 - ekz * ekz)
			);

			// 3. K_{H,i} (Eq. 10)
			//    K_{H,i} = A_{s,i} * (delta*_i - delta_{m,i}) * (k_b*T / u_{A,i})
			//              * (1 - exp(-u_{A,i} / (k_b*T)))

			ParamType KH_i = 0.0;
			if (std::abs(static_cast<double>(uA_i)) > 1e-30)
				KH_i = As_i * (dstar_i - dm_i) * (kbT / uA_i) * (1.0 - exp(-uA_i / kbT));

			// 4. B_i(Theta)
			//    Eq. (28): Hard-disc ASF
			//    B_i = (1 - Theta) * exp(
			//      -(pi*a_i^2 * sum_j(q_j*N_A) + 2*pi*a_i * sum_j(a_j*q_j*N_A)) / (1 - Theta)
			//      - pi*a_i^2 * (sum_j(a_j*q_j*N_A))^2 / (1 - Theta)^2 )

			StateParamType B_i = 0.0;
			if (std::abs(static_cast<double>(Theta)) < 1.0)
			{
				const StateParamType oneMinusTheta = 1.0 - Theta;
				const StateParamType nom1 = pi * a_i * a_i * sumQSurface * NA
					+ 2.0 * pi * a_i * sumAjQj * NA;
				const StateParamType nom2 = pi * a_i * a_i
					* (sumAjQj * NA) * (sumAjQj * NA);

				B_i = oneMinusTheta * exp(-nom1 / oneMinusTheta - nom2 / (oneMinusTheta * oneMinusTheta));
			}

			// 5. u_{lat,i} (Eq. 23)
			//    u_{lat,i} = 3*sqrt(3)*D_hex*N_A
			//                * exp(-kappa*D_hex) / (1 - exp(-3*sqrt(3)/(2*pi)*kappa*D_hex))
			//                * sum_j(q_j * beta_{i,j})
			//    beta_{i,j} computed inline (Eq. 22)

			StateParamType ulat_i = 0.0;

			// sum_j(q_j * beta_{i,j})
			StateParamType betaQSum = 0.0;
			int bndIdx2 = 0;
			for (int j = _startIdx; j < _nComp; ++j)
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

				const StateParamType q_j = y[bndIdx2] / As_j;
				betaQSum += beta_ij * q_j;

				++bndIdx2;
			}

			const StateParamType denom = 1.0 - exp(-(3.0 * std::sqrt(3.0) / (2.0 * pi)) * kappa * Dhex);

			if (std::abs(static_cast<double>(denom)) > 1e-30)
			{
				ulat_i = 3.0 * std::sqrt(3.0) * Dhex * NA
					* exp(-kappa * Dhex) / denom
					* betaQSum;
			} 

			// 6.  K_{v,i} = K_{H,i} * B_i(Theta) * exp(-u_{lat,i} / (k_b*T))
			//    Residual: k_{kin,i} * (K_{v,i} * c_{p,i} - q_{v,i})
			//    Eq. (11)
			const StateParamType Kv_i = KH_i * B_i * exp(-ulat_i / kbT);

			res[bndIdx] = kKin_i * (Kv_i * yCp[i] - y[bndIdx]);

			++bndIdx;
		}

		return 0;
	}

	template <typename RowIterator>
	void jacobianImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* y, double const* yCp, int offsetCp, RowIterator jac, LinearBufferAllocator workSpace) const
	{
		//TODO
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
