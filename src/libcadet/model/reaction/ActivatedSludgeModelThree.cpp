// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "model/reaction/ReactionModelBase.hpp"
#include "model/ExternalFunctionSupport.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "model/Parameters.hpp"
#include "LocalVector.hpp"
#include "SimulationTypes.hpp"
#include "linalg/ActiveDenseMatrix.hpp"
#include "MathUtil.hpp"
#include "Memory.hpp"

#include <algorithm>
#include <unordered_map>
#include <string>
#include <vector>

/*<codegen>
{
	"name": "ActivatedSludgeModelThreeParamHandler",
	"externalName": "ExtActivatedSludgeModelThreeParamHandler",
	"parameters":
		[
			{ "type": "ScalarReactionDependentParameter", "varName": "vMax", "confName": "MM_VMAX"},
			{ "type": "ScalarReactionDependentParameter", "varName": "kMM", "confName": "MM_KMM"},
			{ "type": "ComponentDependentReactionDependentParameter", "varName": "kInhibit", "confName": "MM_KI"}
		]
}
</codegen>*/

/* Parameter description
 ------------------------
*/


namespace cadet
{

namespace model
{

inline const char* ActivatedSludgeModelThreeParamHandler::identifier() CADET_NOEXCEPT { return "ACTIVATED_SLUDGE_MODEL3"; }

inline bool ActivatedSludgeModelThreeParamHandler::validateConfig(unsigned int nReactions, unsigned int nComp, unsigned int const* nBoundStates)
{
	return true;
}

inline const char* ExtActivatedSludgeModelThreeParamHandler::identifier() CADET_NOEXCEPT { return "EXT_ACTIVATED_SLUDGE_MODEL3"; }

inline bool ExtActivatedSludgeModelThreeParamHandler::validateConfig(unsigned int nReactions, unsigned int nComp, unsigned int const* nBoundStates)
{
	return true;
}

namespace
{
	/**
	 * @brief Registers a matrix-valued parameter (row-major storage) with components as rows
	 * @details The matrix-valued parameter has as many rows as there are components in the system.
	 * @param [in,out] parameters Parameter map
	 * @param [in] unitOpIdx Unit operation id
	 * @param [in] parTypeIdx Particle type index
	 * @param [in] paramName Name of the parameter
	 * @param [in] mat Matrix to register
	 */
	inline void registerCompRowMatrix(std::unordered_map<ParameterId, active*>& parameters, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx, const std::string& paramName, cadet::linalg::ActiveDenseMatrix& mat)
	{
		const cadet::StringHash hashName = cadet::hashStringRuntime(paramName);
		cadet::registerParam2DArray(parameters, mat.data(), mat.elements(), [=](bool multi, unsigned int row, unsigned int col)
			{
				return cadet::makeParamId(hashName, unitOpIdx, row, parTypeIdx, cadet::BoundStateIndep, col, cadet::SectionIndep);
			},
			mat.columns()
		);
	}
}


/**
 * @brief Defines a Michaelis-Menten reaction kinetic with simple inhibition
 * @details Implements the Michaelis-Menten kinetics: \f[ \begin{align}
 *              S \nu,
 *          \end{align} \f]
 *          where \f$ S \f$ is the stoichiometric matrix and the fluxes are given by
 *          \f[ \begin{align}
 *              \nu_i = \frac{\mu_{\mathrm{max},i} c_S}{k_{\mathrm{MM},i} + c_S}.
 *          \end{align} \f]
 *          The substrate component \f$ c_S \f$ is identified by the index of the
 *          first negative entry in the stoichiometry of this reaction.
 *			In addition, the reaction might be inhibited by other components. In this
 *          case, the flux has the form
 *          \f[ \begin{align}
 *              \nu_i = \frac{\mu_{\mathrm{max},i} c_S}{k_{\mathrm{MM},i} + c_S} \cdot \frac{1}{1 + \sum_i \frac{1+ k_{\mathrm{I},i,j}}{k_{\mathrm{I},i,j} + c_{\mathrm{I},j}}}.
 *          \end{align} \f]
 *          The value of \f$ k_{\mathrm{I},i,j} \f$ decides whether component \f$ j \f$
 *          inhibits reaction \f$ i \f$. If \f$ k_{\mathrm{I},i,j} \leq 0 \f$, the component
 *          does not inhibit the reaction.
 *          Only reactions in liquid phase are supported (no solid phase or cross-phase reactions).
 * @tparam ParamHandler_t Type that can add support for external function dependence
 */
template <class ParamHandler_t>
class ActivatedSludgeModelThreeBase : public DynamicReactionModelBase
{
public:

	ActivatedSludgeModelThreeBase() : _idxSubstrate(0) { }
	virtual ~ActivatedSludgeModelThreeBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }
	virtual const char* name() const CADET_NOEXCEPT { return ParamHandler_t::identifier(); }

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { _paramHandler.setExternalFunctions(extFuns, size); }
	virtual bool dependsOnTime() const CADET_NOEXCEPT { return ParamHandler_t::dependsOnTime(); }
	virtual bool requiresWorkspace() const CADET_NOEXCEPT { return true; }
	virtual unsigned int workspaceSize(unsigned int nComp, unsigned int totalNumBoundStates, unsigned int const* nBoundStates) const CADET_NOEXCEPT
	{
		return _paramHandler.cacheSize(_stoichiometryBulk.columns(), nComp, totalNumBoundStates) + std::max(_stoichiometryBulk.columns() * sizeof(active), 2 * (_nComp + totalNumBoundStates) * sizeof(double));
	}

	virtual bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
	{
		DynamicReactionModelBase::configureModelDiscretization(paramProvider, nComp, nBound, boundOffset);

		if (paramProvider.exists("MM_STOICHIOMETRY_BULK"))
		{
			const std::size_t numElements = paramProvider.numElements("MM_STOICHIOMETRY_BULK");
			if (numElements % nComp != 0)
				throw InvalidParameterException("Size of field MM_STOICHIOMETRY_BULK must be a positive multiple of NCOMP (" + std::to_string(nComp) + ")");

			const unsigned int nReactions = numElements / nComp;

			_stoichiometryBulk.resize(nComp, nReactions);
			_idxSubstrate = std::vector<int>(nReactions, -1);
		}

		return true;
	}

	virtual unsigned int numReactionsLiquid() const CADET_NOEXCEPT { return _stoichiometryBulk.columns(); }
	virtual unsigned int numReactionsCombined() const CADET_NOEXCEPT { return 0; }

	CADET_DYNAMICREACTIONMODEL_BOILERPLATE

protected:
	ParamHandler_t _paramHandler; //!< Handles parameters and their dependence on external functions

	linalg::ActiveDenseMatrix _stoichiometry;
	


	virtual bool configureStoich(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
	{
		return true;
	}


	virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
	{
		_stoichiometry.resize(_nComp, 13);
		_stoichiometry.setAll(0);

		// SO
		_stoichiometry.native(0, 1) = YSTO_aer - 1;
		_stoichiometry.native(0, 3) = 1 - 1 / YH_aer;
		_stoichiometry.native(0, 5) = -1 * (1 - fXI);
		_stoichiometry.native(0, 7) = -1;
		_stoichiometry.native(0, 9 = -64/14 * 1/YA + 1;
		_stoichiometry.native(0, 10) = -1 * (1 - fXI);
		_stoichiometry.native(0, 12) = 1;

		// SS
		_stoichiometry.native(1, 0) = 1 - fSI;
		_stoichiometry.native(1, 1) = -1;
		_stoichiometry.native(1, 2) = -1;

		// SNH
		_stoichiometry.native(2, 0) = c1n;
		_stoichiometry.native(2, 1) = c2n;
		_stoichiometry.native(2, 2) = c3n;
		_stoichiometry.native(2, 3) = c4n;
		_stoichiometry.native(2, 4) = c5n;
		_stoichiometry.native(2, 5) = c6n;
		_stoichiometry.native(2, 6) = c7n;
		_stoichiometry.native(2, 9) = c10n;
		_stoichiometry.native(2, 10) = c11n;
		_stoichiometry.native(2, 11) = c12n;

		// SNO
		_stoichiometry.native(3, 2) = c3no;
		_stoichiometry.native(3, 4) = c5no;
		_stoichiometry.native(3, 6) = c7no;
		_stoichiometry.native(3, 8) = c9no;
		_stoichiometry.native(3, 9) = c10no;
		_stoichiometry.native(3, 11) = c12no;

		// SN2
		_stoichiometry.native(4, 2) = -c3no;
		_stoichiometry.native(4, 4) = -c5no;
		_stoichiometry.native(4, 6) = -c7no;
		_stoichiometry.native(4, 8) = -c9no;
		_stoichiometry.native(4, 11) = -c12no;

		// SALK
		_stoichiometry.native(5, 0) = c1a;
		_stoichiometry.native(5, 1) = c2a;
		_stoichiometry.native(5, 2) = c3a;
		_stoichiometry.native(5, 3) = c4a;
		_stoichiometry.native(5, 4) = c5a;
		_stoichiometry.native(5, 5) = c6a;
		_stoichiometry.native(5, 6) = c7a;
		_stoichiometry.native(5, 8) = c9a;
		_stoichiometry.native(5, 9) = c10a;
		_stoichiometry.native(5, 10) = c11a;
		_stoichiometry.native(5, 11) = c12a;

		// SI
		_stoichiometry.native(6, 0) = fSI;

		// XI
		_stoichiometry.native(7, 5) = fXI;
		_stoichiometry.native(7, 6) = fXI;
		_stoichiometry.native(7, 10) = fXI;
		_stoichiometry.native(7, 11) = fXI;

		// XS
		_stoichiometry.native(8, 0) = -1;

		// XH
		_stoichiometry.native(9, 3) = 1;
		_stoichiometry.native(9, 4) = 1;
		_stoichiometry.native(9, 5) = -1;
		_stoichiometry.native(9, 6) = -1;

		// XSTO
		_stoichiometry.native(10, 1) = YSO_aer;
		_stoichiometry.native(10, 2) = YSTO_anox;
		_stoichiometry.native(10, 3) = -1 / YH_aer;
		_stoichiometry.native(10, 4) = -1 / YH_anox;
		_stoichiometry.native(10, 7) = -1;
		_stoichiometry.native(10, 8) = -1;

		// XA
		_stoichiometry.native(11, 9) = 1;
		_stoichiometry.native(11, 10) = -1;
		_stoichiometry.native(11, 11) = -1;

		// XMI
		_stoichiometry.native(12, 5) = fXMI_BM;
		_stoichiometry.native(12, 6) = fXMI_BM;
		_stoichiometry.native(12, 10) = fXMI_BM;
		_stoichiometry.native(12, 11) = fXMI_BM;

		//_paramHandler.configure(paramProvider, _stoichiometry.columns(), _nComp, _nBoundStates);
		_paramHandler.registerParameters(_parameters, unitOpIdx, parTypeIdx, _nComp, _nBoundStates);

		// configure stoichiometric matrix from given parameters
		//configureStoich(paramProvider,unitOpIdx, parTypeIdx)

		//registerCompRowMatrix(_parameters, unitOpIdx, parTypeIdx, "MM_STOICHIOMETRY_BULK", _stoichiometryBulk); todo

		return true;
	}

	template <typename StateType, typename ResidualType, typename ParamType, typename FactorType>
	int residualLiquidImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
		StateType const* y, ResidualType* res, const FactorType& factor, LinearBufferAllocator workSpace) const
	{
		typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Calculate fluxes
		typedef typename DoubleActivePromoter<StateType, ParamType>::type flux_t;
		BufferedArray<flux_t> fluxes = workSpace.array<flux_t>(_stoichiometry.columns());
		
		ParamType Kh20 = 10.2;
		ParamType T = 20;
		ParamType iO2 = 0;
		ParamType V = 1;
		ParamType k_sto20 = 13.68;
		ParamType KX = 1;
		ParamType KHO2 = 0.2;
		ParamType KHSS = 3;
		ParamType KHNO3 = 0.5;
		ParamType etaHNO3 = 0.5;
		ParamType KHNH4 = 0.01;
		ParamType KHALK = 0.1;
		ParamType KHSTO = 0.11;
		ParamType muH20 = 3;
		ParamType etaHend = 0.5;
		ParamType bH20 = 0.33;
		ParamType muAUT20 = 1.12;
		ParamType KNO2 = 0.5;
		ParamType KNNH4 = 0.7;
		ParamType KNALK = 0.5;
		ParamType bAUT20 = 0.18;
		ParamType etaNend = 0.5;

		// derived parameters
		double ft04 = exp(-0.04 * (20.0 - static_cast<double>(T)));
		double ft07 = exp(-0.06952 * (20 - static_cast<double>(T)));
		double ft105 = exp(-0.105 * (20 - static_cast<double>(T)));
		double k_sto = static_cast<double>(k_sto20) * ft07;
		double muH = static_cast<double>(muH20) * ft07;
		double bH = static_cast<double>(bH20) * ft07;
		double muAUT = static_cast<double>(muAUT20) * ft105;
		double bAUT = static_cast<double>(bAUT20) * ft105;

		StateType SO = y[0];
		StateType SS = y[1];
		StateType SNH = y[2];
		StateType SNO = y[3];
		StateType SN2 = y[4];
		StateType SALK = y[5];
		// StateType SI = y[6]; // unused
		// StateType XI = y[7]; // unused
		StateType XS = y[8];
		StateType XH = y[9];
		StateType XH_S = XH; // ASM3hC: XH_S = max(XH, 0.1)
		StateType XSTO = y[10];
		StateType XA = y[11];
		// StateType XMI = y[12]; // unused

		
		// p1: Hydrolysis of organic structures
		fluxes[0] = Kh20 * ft04 * XS/XH_S / (XS/XH_S + KX) * XH;
		
		// p2: Aerobic storage of SS
		fluxes[1] = k_sto * SO / (SO + KHO2) * SS / (SS + KHSS) * SNO / (SNO + KHNO3) * XH;

		// p3: Anoxic storage of SS
		fluxes[2] = k_sto * etaHNO3 * KHO2 / (SO + KHO2) * SS / (SS + KHSS) * SNO / (SNO + KHNO3) * XH;

		// p4: Aerobic growth of heterotrophic biomass (XH)
		fluxes[3] = muH * SO / (SO + KHO2) * SNH / (SNH + KHNH4) * SALK / (SALK + KHALK) * XSTO/XH_S / (XSTO/XH_S + KHSTO) * XH;

		// p5: Anoxic growth of heterotrophic biomass (XH, denitrification)
		fluxes[4] = muH * etaHNO3 * KHO2 / (KHO2 + SO) * SNH / (KHNH4 + SNH) * SALK / (KHALK + SALK) * XSTO / XH_S * 1 / (KHSTO + XSTO/XH_S) * SNO / (KHNO3 + SNO) * XH;

		// r6: Aerobic endogenous respiration of heterotroph microorganisms (XH)
		fluxes[5] = bH * SO / (SO + KHO2) * XH;

		// r7: Anoxic endogenous respiration of heterotroph microorganisms (XH)
		fluxes[6] = bH * etaHend * KHO2 / (SO + KHO2) * SNO / (SNO + KHNO3) * XH;

		// r8: Aerobic respiration of internal cell storage products
		fluxes[7] = bH * SO / (SO + KHO2) * XSTO;

		// r9: Anoxic respiration of internal cell storage products 
		fluxes[8] = bH * etaHend * KHO2 / (SO + KHO2) * SNO / (SNO + KHNO3) * XSTO;

		// r10: Aerobic growth of autotrophic biomass (XAUT, nitrification)
		fluxes[9] = muAUT * SO / (SO + KNO2) * SNH / (SNH + KNNH4) * SALK / (SALK + KNALK) * XA;

		// r11: Aerobic endogenous respiration of autotrophic biomass (XAUT)
		fluxes[10] = bAUT * SO / (SO + KHO2) * XA;

		// r12: Anoxic endogenous respiration of autotrophic biomass (XAUT)
		fluxes[11] = bAUT * etaNend * SNO / (SNO + KHNO3) * KHO2 / (SO + KHO2) * XA;

		// r13: Aeration
		fluxes[12] = iO2 / V * 1000;
			
		// Add reaction terms to residual
		_stoichiometry.multiplyVector(static_cast<flux_t*>(fluxes), factor, res);

		return 0;
	}

	template <typename StateType, typename ResidualType, typename ParamType>
	int residualCombinedImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
		StateType const* yLiquid, StateType const* ySolid, ResidualType* resLiquid, ResidualType* resSolid, double factor, LinearBufferAllocator workSpace) const
	{
		std::fill_n(resLiquid, _nComp, 0.0);

		if (_nTotalBoundStates == 0)
			return 0;

		std::fill_n(resSolid, _nTotalBoundStates, 0.0);

		return 0;
	}

	template <typename RowIterator>
	void jacobianLiquidImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* y, double factor, const RowIterator& jac, LinearBufferAllocator workSpace) const
	{
		typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);
		curJac = jac;
			
		// reaction 1: Kh20 * ft04 * XS/XH_S / (XS/XH_S + KX) * XH;
		// dr1/XS = Kh20 * ft04 * XH / (XS/XH_S + KX) - Kh20 * ft04 * XS/XH_S / (XS/XH_S + KX)^2 * XH;
		// dr1/XH_S = -Kh20 * ft04 * XS / (XS/XH_S + KX) / XH_S^2 * XH;
		// dr1/HX = Kh20 * ft04 * XS/XH_S / (XS/XH_S + KX);

		//curJac = jac;
		// curJac[8] =  dr1/XS;
		// curJac[X] =  dr1/XH_S; ??
		// curJac[9] = dr1/XH;

		// curjac++;

		// reaction2: k_sto * SO / (SO + KHO2) * SS / (SS + KHSS) * SNO / (SNO + KHNO3) * XH;
		// dr2/S0 = k_sto * ft07 * SS / (SS + KHSS) * SNO / (SNO + KHNO3) / (SO + KHO2)^2 * XH;
		// dr2/SS = k_sto * ft07 * SO / (SO + KHO2) * SNO / (SNO + KHNO3) / (SS + KHSS)^2 * XH;
		// dr2/SNO = k_sto * ft07 * SO / (SO + KHO2) * SS / (SS + KHSS) / (SNO + KHNO3)^2 * XH;
		// dr2/XH = k_sto * ft07 * SO / (SO + KHO2) * SS / (SS + KHSS) * SNO / (SNO + KHNO3);

		//curJac[-1] = dr2/S0;
		//curJac[0] = dr2/SS;
		//curJac[1] = dr2/SNO;
		//curJac[8] = dr2/XH;

		// curJac++;

		//reaction3: k_sto * etaHNO3 * KHO2 / (SO + KHO2) * SS / (SS + KHSS) * SNO / (SNO + KHNO3) * XH;
		// dr3/S0 = k_sto * ft07 * etaHNO3 * KHO2 / (SO + KHO2)^2 * SS / (SS + KHSS) * SNO / (SNO + KHNO3) * XH;
		// dr3/SS = k_sto * ft07 * etaHNO3 * KHO2 / (SO + KHO2) * SS / (SS + KHSS)^2 * SNO / (SNO + KHNO3) * XH;
		// dr3/SNO = k_sto * ft07 * etaHNO3 * KHO2 / (SO + KHO2) * SS / (SS + KHSS) / (SNO + KHNO3)^2 * XH;
		// dr3/XH = k_sto * ft07 * etaHNO3 * KHO2 / (SO + KHO2) * SS / (SS + KHSS) * SNO / (SNO + KHNO3);

		//reaction4: muH * SO / (SO + KHO2) * SNH / (SNH + KHNH4) * SALK / (SALK + KHALK) * XSTO/XH_S / (XSTO/XH_S + KHSTO) * XH;
		// dr4/S0 = muH * ft07 * SNH / (SNH + KHNH4) * SALK / (SALK + KHALK) * XSTO/XH_S / (XSTO/XH_S + KHSTO) / (SO + KHO2)^2 * XH;
		// dr4/SNH = muH * ft07 * SO / (SO + KHO2) * SALK / (SALK + KHALK) * XSTO/XH_S / (XSTO/XH_S + KHSTO) / (SNH + KHNH4)^2 * XH;
		// dr4/SALK = muH * ft07 * SO / (SO + KHO2) * SNH / (SNH + KHNH4) * XSTO/XH_S / (XSTO/XH_S + KHSTO) / (SALK + KHALK)^2 * XH;
		// dr4/XSTO = muH * ft07 * SO / (SO + KHO2) * SNH / (SNH + KHNH4) * SALK / (SALK + KHALK) / (XSTO/XH_S + KHSTO)^2 * XH;
		// dr4/XH_S = -muH * ft07 * SO / (SO + KHO2) * SNH / (SNH + KHNH4) * SALK / (SALK + KHALK) / XSTO/XH_S^2 * XH;
		// dr4/XH = muH * ft07 * SO / (SO + KHO2) * SNH / (SNH + KHNH4) * SALK / (SALK + KHALK) * XSTO/XH_S / (XSTO/XH_S + KHSTO);

		//reaction5: muH * etaHNO3 * KHO2 / (KHO2 + SO) * SNH / (KHNH4 + SNH) * SALK / (KHALK + SALK) * XSTO / XH_S * 1 / (KHSTO + XSTO/XH_S) * SNO / (KHNO3 + SNO) * XH;
		// dr5/S0 = muH * ft07 * etaHNO3 * KHO2 / (KHO2 + SO)^2 * SNH / (KHNH4 + SNH) * SALK / (KHALK + SALK) * XSTO / XH_S * 1 / (KHSTO + XSTO/XH_S) * SNO / (KHNO3 + SNO) * XH;
		// dr5/SNH = muH * ft07 * etaHNO3 * KHO2 / (KHO2 + SO) * SNH / (KHNH4 + SNH)^2 * SALK / (KHALK + SALK) * XSTO / XH_S * 1 / (KHSTO + XSTO/XH_S) * SNO / (KHNO3 + SNO) * XH;
		// dr5/SALK = muH * ft07 * etaHNO3 * KHO2 / (KHO2 + SO) * SNH / (KHNH4 + SNH) * SALK / (KHALK + SALK)^2 * XSTO / XH_S * 1 / (KHSTO + XSTO/XH_S) * SNO / (KHNO3 + SNO) * XH;
		// dr5/XSTO = muH * ft07 * etaHNO3 * KHO2 / (KHO2 + SO) * SNH / (KHNH4 + SNH) * SALK / (KHALK + SALK) / (XSTO/XH_S + KHSTO)^2 * SNO / (KHNO3 + SNO) * XH;
		// dr5/XH_S = -muH * ft07 * etaHNO3 * KHO2 / (KHO2 + SO) * SNH / (KHNH4 + SNH) * SALK / (KHALK + SALK) / XSTO/XH_S^2 * SNO / (KHNO3 + SNO) * XH;
		// dr5/SNO = muH * ft07 * etaHNO3 * KHO2 / (KHO2 + SO) * SNH / (KHNH4 + SNH) * SALK / (KHALK + SALK) * XSTO / XH_S * 1 / (KHSTO + XSTO/XH_S) / (KHNO3 + SNO)^2 * XH;
		// dr5/XH = muH * ft07 * etaHNO3 * KHO2 / (KHO2 + SO) * SNH / (KHNH4 + SNH) * SALK / (KHALK + SALK) * XSTO / XH_S * 1 / (KHSTO + XSTO/XH_S) * SNO / (KHNO3 + SNO);


		//reaction6: bH * SO / (SO + KHO2) * XH;
		// dr6/S0 = bH * ft07 * XH / (SO + KHO2)^2;
		// dr6/XH = bH * ft07 * SO / (SO + KHO2);
		
		//reaction7: bH * etaHend * KHO2 / (SO + KHO2) * SNO / (SNO + KHNO3) * XH;
		// dr7/S0 = bH * ft07 * etaHend * KHO2 / (SO + KHO2)^2 * SNO / (SNO + KHNO3) * XH;
		// dr7/SNO = bH * ft07 * etaHend * KHO2 / (SO + KHO2) * SNO / (SNO + KHNO3)^2 * XH;
		// dr7/XH = bH * ft07 * etaHend * KHO2 / (SO + KHO2) * SNO / (SNO + KHNO3) * XH;

		//reaction8: bH * SO / (SO + KHO2) * XSTO;
		// dr8/S0 = bH * ft07 * XSTO / (SO + KHO2)^2;
		// dr8/XSTO = bH * ft07 * SO / (SO + KHO2);

		//reaction9: bH * etaHend * KHO2 / (SO + KHO2) * SNO / (SNO + KHNO3) * XSTO;
		// dr9/S0 = bH * ft07 * etaHend * KHO2 / (SO + KHO2)^2 * SNO / (SNO + KHNO3) * XSTO;
		// dr9/SNO = bH * ft07 * etaHend * KHO2 / (SO + KHO2) * SNO / (SNO + KHNO3)^2 * XSTO;

		//reaction10: muAUT * SO / (SO + KNO2) * SNH / (SNH + KNNH4) * SALK / (SALK + KNALK) * XA;
		// dr10/S0 = muAUT * ft105 * SNH / (SNH + KNNH4) * SALK / (SALK + KNALK) * XA / (SO + KNO2)^2;
		// dr10/SNH = muAUT * ft105 * SO / (SO + KNO2) * SALK / (SALK + KNALK) / (SNH + KNNH4)^2 * XA;
		// dr10/SALK = muAUT * ft105 * SO / (SO + KNO2) * SNH / (SNH + KNNH4) / (SALK + KNALK)^2 * XA;

		//reaction11: bAUT * SO / (SO + KHO2) * XA;
		// dr11/S0 = bAUT * ft105 * XA / (SO + KHO2)^2;
		// dr11/XA = bAUT * ft105 * SO / (SO + KHO2);


		//reaction12: bAUT * etaNend * SNO / (SNO + KHNO3) * KHO2 / (SO + KHO2) * XA;
		// dr12/S0 = bAUT * ft105 * etaNend * SNO / (SNO + KHNO3) * KHO2 / (SO + KHO2)^2 * XA;
		// dr12/SNO = bAUT * ft105 * etaNend * KHO2 / (SO + KHO2) * SNO / (SNO + KHNO3)^2 * XA;
		// dr12/XA = bAUT * ft105 * etaNend * SNO / (SNO + KHNO3) * KHO2 / (SO + KHO2) * XA;


	
		// StateType SI = y[6]; // unused
		// StateType XI = y[7]; // unused
		// StateType XMI = y[12]; // unused

	
		
	}

	template <typename RowIteratorLiquid, typename RowIteratorSolid>
	void jacobianCombinedImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yLiquid, double const* ySolid, double factor, const RowIteratorLiquid& jacLiquid, const RowIteratorSolid& jacSolid, LinearBufferAllocator workSpace) const
	{
	}
};

typedef ActivatedSludgeModelThreeBase<ActivatedSludgeModelThreeParamHandler> ActivatedSludgeModelThreeReaction;
typedef ActivatedSludgeModelThreeBase<ExtActivatedSludgeModelThreeParamHandler> ExternalActivatedSludgeModelThreeReaction;

namespace reaction
{
	void registerActivatedSludgeModelThreeReaction(std::unordered_map<std::string, std::function<model::IDynamicReactionModel*()>>& reactions)
	{
		reactions[ActivatedSludgeModelThreeReaction::identifier()] = []() { return new ActivatedSludgeModelThreeReaction(); };
		reactions[ExternalActivatedSludgeModelThreeReaction::identifier()] = []() { return new ExternalActivatedSludgeModelThreeReaction(); };
	}
}  // namespace reaction

}  // namespace model

}  // namespace cadet
