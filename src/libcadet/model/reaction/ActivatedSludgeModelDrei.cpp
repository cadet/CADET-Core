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
	"name": "MichaelisMentenParamHandler",
	"externalName": "ExtMichaelisMentenParamHandler",
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

inline const char* MichaelisMentenParamHandler::identifier() CADET_NOEXCEPT { return "MICHAELIS_MENTEN"; }

inline bool MichaelisMentenParamHandler::validateConfig(unsigned int nReactions, unsigned int nComp, unsigned int const* nBoundStates)
{
	return true;
}

inline const char* ExtMichaelisMentenParamHandler::identifier() CADET_NOEXCEPT { return "EXT_MICHAELIS_MENTEN"; }

inline bool ExtMichaelisMentenParamHandler::validateConfig(unsigned int nReactions, unsigned int nComp, unsigned int const* nBoundStates)
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
class MichaelisMentenReactionBase : public DynamicReactionModelBase
{
public:

	MichaelisMentenReactionBase() : _idxSubstrate(0) { }
	virtual ~MichaelisMentenReactionBase() CADET_NOEXCEPT { }

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
		_paramHandler.configure(paramProvider, _stoichiometryBulk.columns(), _nComp, _nBoundStates);
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
		/* | Component | Index |
			| -------- - | ---- - |
			| SO | 0     |
			| SS | 1     |
			| SNH | 2    |
			| SNO | 3    |
			| SN2 | 4    |
			| SALK | 5   |
			| SI | 6     |
			| XI | 7     |
			| XS | 8     |
			| XH | 9     |
			| XSTO | 10  |
			| XA | 11	 |
			| XMI | 12	 | */
		
		// flux 0: kh2o fto(T) HS/XH/(...
		ParamType Kh20 = 0.0;
		ParamType ft04 = 0.0; // eigentlich abhängig von temperatur (später Parameter-Parameter Dependencies) 
		//ParamType KX = 0.0; // TODO: is this a param?

		StateType SO = y[0];
		StateType SS = y[1];
		StateType SNH = y[2];
		StateType SNO = y[3];

		StateType XS = y[8];
		StateType XH = y[9];
		
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





		// for each inhibitor component
		curJac = jac;
		for (int row = 0; row < _nComp; ++row, ++curJac)
		{
			const double colFactor = static_cast<double>(_stoichiometryBulk.native(row, r)) * factor;
			for (int comp = 0; comp < _nComp; ++comp)
			{
				const double kI = static_cast<double>(p->kInhibit[_nComp * r + comp]);
				if (kI <= 0.0)
					continue;

				double dvdi = - (vMax * y[idxSubs] * (kMM + y[idxSubs])) / (denom * denom * kI);
				curJac[comp - static_cast<int>(row)] += colFactor * dvdi;
			}
		}
		
	}

	template <typename RowIteratorLiquid, typename RowIteratorSolid>
	void jacobianCombinedImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yLiquid, double const* ySolid, double factor, const RowIteratorLiquid& jacLiquid, const RowIteratorSolid& jacSolid, LinearBufferAllocator workSpace) const
	{
	}
};

typedef MichaelisMentenReactionBase<MichaelisMentenParamHandler> MichaelisMentenReaction;
typedef MichaelisMentenReactionBase<ExtMichaelisMentenParamHandler> ExternalMichaelisMentenReaction;

namespace reaction
{
	void registerMichaelisMentenReaction(std::unordered_map<std::string, std::function<model::IDynamicReactionModel*()>>& reactions)
	{
		reactions[MichaelisMentenReaction::identifier()] = []() { return new MichaelisMentenReaction(); };
		reactions[ExternalMichaelisMentenReaction::identifier()] = []() { return new ExternalMichaelisMentenReaction(); };
	}
}  // namespace reaction

}  // namespace model

}  // namespace cadet
