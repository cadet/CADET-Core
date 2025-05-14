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
			{ "type": "ScalarParameter", "varName": "Kh20", "confName": "ASM3_KH_20"},
			{ "type": "ScalarParameter", "varName": "T", "confName": "ASM3_T"},
			{ "type": "ScalarParameter", "varName": "iO2", "confName": "ASM3_I_O2"},
			{ "type": "ScalarParameter", "varName": "V", "confName": "ASM3_V"},
			{ "type": "ScalarParameter", "varName": "k_sto20", "confName": "ASM3_K_STO_20"},
			{ "type": "ScalarParameter", "varName": "KX", "confName": "ASM3_KX"},
			{ "type": "ScalarParameter", "varName": "KHO2", "confName": "ASM3_KHO2"},
			{ "type": "ScalarParameter", "varName": "KHSS", "confName": "ASM3_KHSS"},
			{ "type": "ScalarParameter", "varName": "KHNO3", "confName": "ASM3_KHNO3"},
			{ "type": "ScalarParameter", "varName": "etaHNO3", "confName": "ASM3_ETA_HNO3"},
			{ "type": "ScalarParameter", "varName": "KHNH4", "confName": "ASM3_KHNH4"},
			{ "type": "ScalarParameter", "varName": "KHALK", "confName": "ASM3_KHALK"},
			{ "type": "ScalarParameter", "varName": "KHSTO", "confName": "ASM3_KHSTO"},
			{ "type": "ScalarParameter", "varName": "muH20", "confName": "ASM3_MU_H_20"},
			{ "type": "ScalarParameter", "varName": "etaHend", "confName": "ASM3_ETA_H_END"},
			{ "type": "ScalarParameter", "varName": "bH20", "confName": "ASM3_B_H_20"},
			{ "type": "ScalarParameter", "varName": "muAUT20", "confName": "ASM3_MU_AUT_20"},
			{ "type": "ScalarParameter", "varName": "KNO2", "confName": "ASM3_KNO2"},
			{ "type": "ScalarParameter", "varName": "KNNH4", "confName": "ASM3_KNNH4"},
			{ "type": "ScalarParameter", "varName": "KNALK", "confName": "ASM3_KNALK"},
			{ "type": "ScalarParameter", "varName": "bAUT20", "confName": "ASM3_B_AUT_20"},
			{ "type": "ScalarParameter", "varName": "etaNend", "confName": "ASM3_ETA_N_END"}
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

	ActivatedSludgeModelThreeBase() CADET_NOEXCEPT { }
	virtual ~ActivatedSludgeModelThreeBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }
	virtual const char* name() const CADET_NOEXCEPT { return ParamHandler_t::identifier(); }

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { _paramHandler.setExternalFunctions(extFuns, size); }
	virtual bool dependsOnTime() const CADET_NOEXCEPT { return ParamHandler_t::dependsOnTime(); }
	virtual bool requiresWorkspace() const CADET_NOEXCEPT { return true; }
	virtual unsigned int workspaceSize(unsigned int nComp, unsigned int totalNumBoundStates, unsigned int const* nBoundStates) const CADET_NOEXCEPT
	{
		return _paramHandler.cacheSize(_stoichiometry.columns(), nComp, totalNumBoundStates) + std::max(_stoichiometry.columns() * sizeof(active), 2 * (_nComp + totalNumBoundStates) * sizeof(double));
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

			_stoichiometry.resize(nComp, nReactions);
		}

		return true;
	}

	virtual unsigned int numReactionsLiquid() const CADET_NOEXCEPT { return _stoichiometry.columns(); }
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

		// parameter set ASM3hC
		double iNSI = paramProvider.getDouble("ASM_iNSI");
		double iNSS = paramProvider.getDouble("ASM_iNSS");
		double iNXI = paramProvider.getDouble("ASM_iNXI");
		double iNXS = paramProvider.getDouble("ASM_iNXS");
		double iNBM = paramProvider.getDouble("ASM_iNBM");
		
		double fSI = paramProvider.getDouble("ASM_fSI");
		double fXI = paramProvider.getDouble("ASM_fXI");

		double YH_aer = paramProvider.getDouble("ASM_YH_aer");
		double YH_anox = paramProvider.getDouble("ASM_YH_anox");
		double YSTO_aer = paramProvider.getDouble("ASM_YSTO_aer");
		double YSTO_anox = paramProvider.getDouble("ASM_YSTO_anox");
		double YA = paramProvider.getDouble("ASM_YA");
		
		double fiSS_BM_prod = paramProvider.getDouble("ASM_fiSS_BM_prod");
		double iVSS_BM = paramProvider.getDouble("ASM_iVSS_BM");
		double iTSS_VSS_BM = paramProvider.getDouble("ASM_iTSS_VSS_BM");

		// internal variables
		double fXMI_BM = fiSS_BM_prod * fXI * iVSS_BM * (iTSS_VSS_BM - 1);
		double c1n = iNXS - iNSI * fSI - (1 - fSI) * iNSS;
		double c2n = iNSS;
		double c3n = iNSS;
		double c4n = -iNBM;
		double c5n = c4n;
		double c6n = -fXI * iNXI + iNBM;
		double c7n = c6n;
		double c10n = -1 / YA - iNBM;
		double c11n = -fXI * iNXI + iNBM;
		double c12n = c11n;
		double c3no = (YSTO_anox - 1) / (40 / 14);
		double c5no = (1 - 1/YH_anox) / (40 / 14);
		double c7no = (fXI - 1) / (40 / 14);
		double c9no = -14 / 40;
		double c10no = 1 / YA;
		double c12no = c7no;
		double c1a = c1n / 14;
		double c2a = c2n / 14;
		double c3a = (c3n - c3no) / 14;
		double c4a = c4n / 14;
		double c5a = (c5n - c5no) / 14;
		double c6a = c6n / 14;
		double c7a = (c7n - c7no) / 14;
		double c9a = 1 / 40;
		double c10a = (c10n - c10no) / 14;
		double c11a = c11n / 14;
		double c12a = (c12n - c12no) / 14;

		// SO
		_stoichiometry.native(0, 1) = YSTO_aer - 1;
		_stoichiometry.native(0, 3) = 1 - 1 / YH_aer;
		_stoichiometry.native(0, 5) = -1 * (1 - fXI);
		_stoichiometry.native(0, 7) = -1;
		_stoichiometry.native(0, 9) = -64/14 * 1/YA + 1;
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
		_stoichiometry.native(10, 1) = YSTO_aer;
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
		
		flux_t Kh20		= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->Kh20);
		flux_t T		= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->T);
		flux_t iO2		= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->iO2);
		flux_t V		= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->V);
		flux_t k_sto20	= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->k_sto20);
		flux_t KX		= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->KX);
		flux_t KHO2		= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->KHO2);
		flux_t KHSS		= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->KHSS);
		flux_t KHNO3	= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->KHNO3);
		flux_t etaHNO3	= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->etaHNO3);
		flux_t KHNH4	= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->KHNH4);
		flux_t KHALK	= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->KHALK);
		flux_t KHSTO	= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->KHSTO);
		flux_t muH20	= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->muH20);
		flux_t etaHend	= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->etaHend);
		flux_t bH20		= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->bH20);
		flux_t muAUT20	= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->muAUT20);
		flux_t KNO2		= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->KNO2);
		flux_t KNNH4	= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->KNNH4);
		flux_t KNALK	= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->KNALK);
		flux_t bAUT20	= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->bAUT20);
		flux_t etaNend	= static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->etaNend);

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
		// TODO: is V in litres?
		fluxes[12] = iO2 / V;
			
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
		RowIterator curJac = jac;

		double Kh20		= static_cast<double>(p->Kh20);
		double T		= static_cast<double>(p->T);
		double iO2		= static_cast<double>(p->iO2);
		double V		= static_cast<double>(p->V);
		double k_sto20	= static_cast<double>(p->k_sto20);
		double KX		= static_cast<double>(p->KX);
		double KHO2		= static_cast<double>(p->KHO2);
		double KHSS		= static_cast<double>(p->KHSS);
		double KHNO3	= static_cast<double>(p->KHNO3);
		double etaHNO3	= static_cast<double>(p->etaHNO3);
		double KHNH4	= static_cast<double>(p->KHNH4);
		double KHALK	= static_cast<double>(p->KHALK);
		double KHSTO	= static_cast<double>(p->KHSTO);
		double muH20	= static_cast<double>(p->muH20);
		double etaHend	= static_cast<double>(p->etaHend);
		double bH20		= static_cast<double>(p->bH20);
		double muAUT20	= static_cast<double>(p->muAUT20);
		double KNO2		= static_cast<double>(p->KNO2);
		double KNNH4	= static_cast<double>(p->KNNH4);
		double KNALK	= static_cast<double>(p->KNALK);
		double bAUT20	= static_cast<double>(p->bAUT20);
		double etaNend	= static_cast<double>(p->etaNend);

		// derived parameters
		double ft04 = exp(-0.04 * (20.0 - static_cast<double>(p->T)));
		double ft07 = exp(-0.06952 * (20 - static_cast<double>(p->T)));
		double ft105 = exp(-0.105 * (20 - static_cast<double>(p->T)));
		double k_sto = static_cast<double>(p->k_sto20) * ft07;
		double muH = static_cast<double>(p->muH20) * ft07;
		double bH = static_cast<double>(p->bH20) * ft07;
		double muAUT = static_cast<double>(p->muAUT20) * ft105;
		double bAUT = static_cast<double>(p->bAUT20) * ft105;

		double SO = y[0];
		double SS = y[1];
		double SNH = y[2];
		double SNO = y[3];
		double SN2 = y[4];
		double SALK = y[5];
		double SI = y[6]; // unused
		double XI = y[7]; // unused
		double XS = y[8];
		double XH = y[9];
		double XH_S = XH; // ASM3hC: XH_S = max(XH, 0.1)
		double XSTO = y[10];
		double XA = y[11];
		double XMI = y[12]; // unused

		size_t idxSO = 0;
		size_t idxSS = 1;
		size_t idxSNH = 2;
		size_t idxSNO = 3;
		size_t idxSN2 = 4;
		size_t idxSALK = 5;
		size_t idxSI = 6;
		size_t idxXI = 7;
		size_t idxXS = 8;
		size_t idxXH = 9;
		size_t idxXSTO = 10;
		size_t idxXA = 11;
		size_t idxXMI = 12;

		double d[13][13] = {};
		//memset(d, 0, 13 * 13 * sizeof(double));
		// p1: Kh20 * ft04 * XS/XH_S / (XS/XH_S + KX) * XH;
		d[0][idxXS] = Kh20 * ft04 * XH / (XS/XH_S + KX) - Kh20 * ft04 * XS/XH_S / ((XS/XH_S + KX) * (XS/XH_S + KX)) * XH;
		//double dp1_XH_S = -Kh20 * ft04 * XS / (XS/XH_S + KX) / (XH_S * XH_S) * XH;
		d[0][idxXH] = Kh20 * ft04 * XS/XH_S / (XS/XH_S + KX);

		/*
		curJac = jac;
		size_t curIdxXS = idxXS;
		for (size_t i = 0; i < _stoichiometry.rows(); i++) {
			curJac[idxXS] += _stoichiometry.native(i, 0) * dp1_XS;
			curJac[idxXH] += _stoichiometry.native(i, 0) * dp1_XH;
			curJac++;
			curIdxXS--;
		}
		*/

		// reaction2: k_sto * SO / (SO + KHO2) * SS / (SS + KHSS) * SNO / (SNO + KHNO3) * XH;
		d[1][idxSO] = k_sto * ft07 * SS / (SS + KHSS) * SNO / (SNO + KHNO3) / ((SO + KHO2) * (SO + KHO2)) * XH;
		d[1][idxSS] = k_sto * ft07 * SO / (SO + KHO2) * SNO / (SNO + KHNO3) / ((SS + KHSS) * (SS + KHSS)) * XH;
		d[1][idxSNO] = k_sto * ft07 * SO / (SO + KHO2) * SS / (SS + KHSS) / ((SNO + KHNO3) * (SNO + KHNO3)) * XH;
		d[1][idxXH] = k_sto * ft07 * SO / (SO + KHO2) * SS / (SS + KHSS) * SNO / (SNO + KHNO3);

		//reaction3: k_sto * etaHNO3 * KHO2 / (SO + KHO2) * SS / (SS + KHSS) * SNO / (SNO + KHNO3) * XH;
		d[2][idxSO] = k_sto * ft07 * etaHNO3 * KHO2 / ((SO + KHO2) * (SO + KHO2)) * SS / (SS + KHSS) * SNO / (SNO + KHNO3) * XH;
		d[2][idxSS] = k_sto * ft07 * etaHNO3 * KHO2 / (SO + KHO2) * SS / ((SS + KHSS) * (SS + KHSS)) * SNO / (SNO + KHNO3) * XH;
		d[2][idxSNO] = k_sto * ft07 * etaHNO3 * KHO2 / (SO + KHO2) * SS / (SS + KHSS) / ((SNO + KHNO3) * (SNO + KHNO3)) * XH;
		d[2][idxXH] = k_sto * ft07 * etaHNO3 * KHO2 / (SO + KHO2) * SS / (SS + KHSS) * SNO / (SNO + KHNO3);

		//reaction4: muH * SO / (SO + KHO2) * SNH / (SNH + KHNH4) * SALK / (SALK + KHALK) * XSTO/XH_S / (XSTO/XH_S + KHSTO) * XH;
		d[3][idxSO] = muH * ft07 * SNH / (SNH + KHNH4) * SALK / (SALK + KHALK) * XSTO/XH_S / (XSTO/XH_S + KHSTO) / ((SO + KHO2) * (SO + KHO2)) * XH;
		d[3][idxSNH] = muH * ft07 * SO / (SO + KHO2) * SALK / (SALK + KHALK) * XSTO/XH_S / (XSTO/XH_S + KHSTO) / ((SNH + KHNH4) * (SNH + KHNH4)) * XH;
		d[3][idxSALK] = muH * ft07 * SO / (SO + KHO2) * SNH / (SNH + KHNH4) * XSTO/XH_S / (XSTO/XH_S + KHSTO) / ((SALK + KHALK) * (SALK + KHALK)) * XH;
		d[3][idxXSTO] = muH * ft07 * SO / (SO + KHO2) * SNH / (SNH + KHNH4) * SALK / (SALK + KHALK) / ((XSTO/XH_S + KHSTO) * (XSTO/XH_S + KHSTO)) * XH;
		//double dp4_XH_S = -muH * ft07 * SO / (SO + KHO2) * SNH / (SNH + KHNH4) * SALK / (SALK + KHALK) / (XSTO/XH_S * XSTO/XH_S) * XH;
		d[3][idxXH] = muH * ft07 * SO / (SO + KHO2) * SNH / (SNH + KHNH4) * SALK / (SALK + KHALK) * XSTO/XH_S / (XSTO/XH_S + KHSTO);

		//reaction5: muH * etaHNO3 * KHO2 / (KHO2 + SO) * SNH / (KHNH4 + SNH) * SALK / (KHALK + SALK) * XSTO / XH_S * 1 / (KHSTO + XSTO/XH_S) * SNO / (KHNO3 + SNO) * XH;
		d[4][idxSO] = muH * ft07 * etaHNO3 * KHO2 / ((KHO2 + SO) * (KHO2 + SO)) * SNH / (KHNH4 + SNH) * SALK / (KHALK + SALK) * XSTO / XH_S * 1 / (KHSTO + XSTO/XH_S) * SNO / (KHNO3 + SNO) * XH;
		d[4][idxSNH] = muH * ft07 * etaHNO3 * KHO2 / (KHO2 + SO) * SNH / ((KHNH4 + SNH) * (KHNH4 + SNH)) * SALK / (KHALK + SALK) * XSTO / XH_S * 1 / (KHSTO + XSTO/XH_S) * SNO / (KHNO3 + SNO) * XH;
		d[4][idxSALK] = muH * ft07 * etaHNO3 * KHO2 / (KHO2 + SO) * SNH / (KHNH4 + SNH) * SALK / ((KHALK + SALK) * (KHALK + SALK)) * XSTO / XH_S * 1 / (KHSTO + XSTO/XH_S) * SNO / (KHNO3 + SNO) * XH;
		d[4][idxXSTO] = muH * ft07 * etaHNO3 * KHO2 / (KHO2 + SO) * SNH / (KHNH4 + SNH) * SALK / (KHALK + SALK) / ((XSTO/XH_S + KHSTO) * (XSTO/XH_S + KHSTO)) * SNO / (KHNO3 + SNO) * XH;
		//double dp5_XH_S = -muH * ft07 * etaHNO3 * KHO2 / (KHO2 + SO) * SNH / (KHNH4 + SNH) * SALK / (KHALK + SALK) / (XSTO/XH_S * XSTO/XH_S) * SNO / (KHNO3 + SNO) * XH;
		d[4][idxSNO] = muH * ft07 * etaHNO3 * KHO2 / (KHO2 + SO) * SNH / (KHNH4 + SNH) * SALK / (KHALK + SALK) * XSTO / XH_S * 1 / (KHSTO + XSTO/XH_S) / ((KHNO3 + SNO) * (KHNO3 + SNO)) * XH;
		d[4][idxXH] = muH * ft07 * etaHNO3 * KHO2 / (KHO2 + SO) * SNH / (KHNH4 + SNH) * SALK / (KHALK + SALK) * XSTO / XH_S * 1 / (KHSTO + XSTO/XH_S) * SNO / (KHNO3 + SNO);


		//reaction6: bH * SO / (SO + KHO2) * XH;
		d[5][idxSO] = bH * ft07 * XH / ((SO + KHO2) * (SO + KHO2));
		d[5][idxXH] = bH * ft07 * SO / (SO + KHO2);
		
		//reaction7: bH * etaHend * KHO2 / (SO + KHO2) * SNO / (SNO + KHNO3) * XH;
		d[6][idxSO] = bH * ft07 * etaHend * KHO2 / ((SO + KHO2) * (SO + KHO2)) * SNO / (SNO + KHNO3) * XH;
		d[6][idxSNO] = bH * ft07 * etaHend * KHO2 / (SO + KHO2) * SNO / ((SNO + KHNO3) * (SNO + KHNO3)) * XH;
		d[6][idxXH] = bH * ft07 * etaHend * KHO2 / (SO + KHO2) * SNO / (SNO + KHNO3) * XH;

		//reaction8: bH * SO / (SO + KHO2) * XSTO;
		d[7][idxSO] = bH * ft07 * XSTO / ((SO + KHO2) * (SO + KHO2));
		d[7][idxXSTO] = bH * ft07 * SO / (SO + KHO2);

		//reaction9: bH * etaHend * KHO2 / (SO + KHO2) * SNO / (SNO + KHNO3) * XSTO;
		d[8][idxSO] = bH * ft07 * etaHend * KHO2 / ((SO + KHO2) * (SO + KHO2)) * SNO / (SNO + KHNO3) * XSTO;
		d[8][idxSNO] = bH * ft07 * etaHend * KHO2 / (SO + KHO2) * SNO / ((SNO + KHNO3) * (SNO + KHNO3)) * XSTO;

		//reaction10: muAUT * SO / (SO + KNO2) * SNH / (SNH + KNNH4) * SALK / (SALK + KNALK) * XA;
		d[9][idxSO] = muAUT * ft105 * SNH / (SNH + KNNH4) * SALK / (SALK + KNALK) * XA / ((SO + KNO2) * (SO + KNO2));
		d[9][idxSNH] = muAUT * ft105 * SO / (SO + KNO2) * SALK / (SALK + KNALK) / ((SNH + KNNH4) * (SNH + KNNH4)) * XA;
		d[9][idxSALK] = muAUT * ft105 * SO / (SO + KNO2) * SNH / (SNH + KNNH4) / ((SALK + KNALK) * (SALK + KNALK)) * XA;

		//reaction11: bAUT * SO / (SO + KHO2) * XA;
		d[10][idxSO] = bAUT * ft105 * XA / ((SO + KHO2) * (SO + KHO2));
		d[10][idxXA] = bAUT * ft105 * SO / (SO + KHO2);


		//reaction12: bAUT * etaNend * SNO / (SNO + KHNO3) * KHO2 / (SO + KHO2) * XA;
		d[11][idxSO] = bAUT * ft105 * etaNend * SNO / (SNO + KHNO3) * KHO2 / ((SO + KHO2) * (SO + KHO2)) * XA;
		d[11][idxSNO] = bAUT * ft105 * etaNend * KHO2 / (SO + KHO2) * SNO / ((SNO + KHNO3) * (SNO + KHNO3)) * XA;
		d[11][idxXA] = bAUT * ft105 * etaNend * SNO / (SNO + KHNO3) * KHO2 / (SO + KHO2) * XA;

		for (size_t rIdx = 0; rIdx < 13; rIdx++) {
			for (size_t compIdx = 0; compIdx < 13; compIdx++) {
				// TODO consider configurable component indices
				for (size_t i = 0; i < _stoichiometry.rows(); i++) {
					curJac = &jac[i];
					curJac[compIdx] += _stoichiometry.native(i, rIdx) * d[rIdx][compIdx];
				}
			}
		}
		
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
