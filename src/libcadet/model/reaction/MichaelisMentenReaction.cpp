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
			_idxInhibitor.resize(nReactions);
			_idxSubstrate.resize(nReactions);
			_oldInterface = false;

		}

		return true;
	}

	virtual unsigned int numReactionsLiquid() const CADET_NOEXCEPT { return _stoichiometryBulk.columns(); }
	virtual unsigned int numReactionsCombined() const CADET_NOEXCEPT { return 0; }

	CADET_DYNAMICREACTIONMODEL_BOILERPLATE

protected:
	ParamHandler_t _paramHandler; //!< Handles parameters and their dependence on external functions

	linalg::ActiveDenseMatrix _stoichiometryBulk;
	std::vector<std::vector<int>> _idxSubstrate;
	std::vector<std::vector<int>> _idxInhibitor;
	bool _oldInterface;

	virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
	{
		_paramHandler.configure(paramProvider, _stoichiometryBulk.columns(), _nComp, _nBoundStates);
		_paramHandler.registerParameters(_parameters, unitOpIdx, parTypeIdx, _nComp, _nBoundStates);

		if ((_stoichiometryBulk.columns() > 0) && ((_paramHandler.vMax().size() < _stoichiometryBulk.columns()) || (_paramHandler.kMM().size() < _stoichiometryBulk.columns())))
			throw InvalidParameterException("MM_VMAX and MM_KMM have to have the same size (number of reactions)");
		
		if (_stoichiometryBulk.rows() > 1 && (_paramHandler.kMM().size() == _stoichiometryBulk.columns() * _stoichiometryBulk.rows()))
			_oldInterface = false;
		else
			_oldInterface = true;

		if ((_stoichiometryBulk.columns() > 0) && (_paramHandler.kInhibit().size() < _stoichiometryBulk.columns() * _nComp))
			throw InvalidParameterException("MM_KI have to have the size (number of reactions) x (number of components)");

		if (paramProvider.exists("MM_STOICHIOMETRY_BULK"))
		{
			const std::vector<double> s = paramProvider.getDoubleArray("MM_STOICHIOMETRY_BULK");

			if (s.size() != _stoichiometryBulk.elements())
				throw InvalidParameterException("MM_STOICHIOMETRY_BULK has changed size (number of reactions changed)");

			std::copy(s.begin(), s.end(), _stoichiometryBulk.data());

			// Find substrate index (first educt)
			_idxSubstrate.clear();
			_idxInhibitor.clear();
			const std::vector<double> ki = paramProvider.getDoubleArray("MM_KI");
			for (int r = 0; r < _stoichiometryBulk.columns(); ++r)
			{
				std::vector<int> idxSubstrateReaction_r;
				std::vector<int> idxInhibitorReaction_r;
				for (int c = 0; c < _nComp; ++c)
				{
					if (_stoichiometryBulk.native(c, r) < 0.0)
						idxSubstrateReaction_r.push_back(c);

					double kI = ki[_nComp * r + c];
					if (kI > 0.0)
						idxInhibitorReaction_r.push_back(c);

				}
				if(idxSubstrateReaction_r.empty())
					throw InvalidParameterException("No substrate found in reaction " + std::to_string(r));

				//throw error if inhibitor is substrate
				if(!idxInhibitorReaction_r.empty())
				{
					// check if inhibitor is substrate
					for (int idxSubstrateReaction_r : idxSubstrateReaction_r)
					{
						if(std::find(idxInhibitorReaction_r.begin(), idxInhibitorReaction_r.end(), idxSubstrateReaction_r) != idxInhibitorReaction_r.end())
							throw InvalidParameterException("Inhibitor is also substrate in reaction " + std::to_string(r) + "this is not supportet yet");
					}
				}

				_idxSubstrate.push_back(idxSubstrateReaction_r);
				_idxInhibitor.push_back(idxInhibitorReaction_r);
			}

		}
		registerCompRowMatrix(_parameters, unitOpIdx, parTypeIdx, "MM_STOICHIOMETRY_BULK", _stoichiometryBulk);
		return true;
	}

	template <typename StateType, typename ResidualType, typename ParamType, typename FactorType>
	int residualLiquidImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
		StateType const* y, ResidualType* res, const FactorType& factor, LinearBufferAllocator workSpace) const
	{
		typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Calculate fluxes
		typedef typename DoubleActivePromoter<StateType, ParamType>::type flux_t;
		BufferedArray<flux_t> fluxes = workSpace.array<flux_t>(_stoichiometryBulk.columns());
		for (int r = 0; r < _stoichiometryBulk.columns(); ++r)
		{
			int nSub = _idxSubstrate[r].size();
			int nInh = _idxInhibitor[r].size();
			flux_t vProd = 1.0;
			for (int sIdx = 0; sIdx < nSub; sIdx++)
			{
				int s = _idxSubstrate[r][sIdx];
				flux_t inhSum = 0;
				for (int iIdx = 0; iIdx < nInh; iIdx++)
				{
					int i = _idxInhibitor[r][iIdx];
					flux_t ki_rsi = static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->kInhibit[_nComp * r + i]);
					inhSum += y[i] / ki_rsi;
				}

				flux_t kMM_rs = static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->kMM[r]);
				if (!_oldInterface)
					// kMM_rs = kMM[r][s]
					 kMM_rs = static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->kMM[_nComp * r + s]);

				vProd *= y[s] / ((kMM_rs + y[s]) * (1+inhSum));

			}

			flux_t vmax_r = static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->vMax[r]);
			fluxes[r] = vmax_r * vProd;
		}

		// Add reaction terms to residual
		_stoichiometryBulk.multiplyVector(static_cast<flux_t*>(fluxes), factor, res);

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

		for (int r = 0; r < _stoichiometryBulk.columns(); ++r)
		{
			// Calculate flux and inhibitor sum
			int nSub = _idxSubstrate[r].size();
			int nInh = _idxInhibitor[r].size();

			double flux = 1.0;
			std::vector<double> substratFlux(nSub, 0.0);
			std::vector<double> inhSumOfSub(nSub, 0.0);

			for (int sIdx = 0; sIdx < nSub; sIdx++)
			{
				int s = _idxSubstrate[r][sIdx];
				double inhSum = 0.0;
				for (int iIdx = 0; iIdx < nInh; iIdx++)
				{
					int i = _idxInhibitor[r][iIdx];
					double kI_rsi = static_cast<double>(p->kInhibit[_nComp * r + i]);

					inhSum += y[i] / kI_rsi;
				}

				double kMM_rs = static_cast<double>(p->kMM[r]);
				if (!_oldInterface)
					// kMM_rs = kMM[r][s]
					 kMM_rs = static_cast<double>(p->kMM[_nComp * r + s]);

				double sFlux = y[s] / ((kMM_rs + y[s]) * (1 + inhSum));
				flux *= sFlux;
				// save flux and inhibitor sum for each substrat
				substratFlux[sIdx] = sFlux;
				inhSumOfSub[sIdx] = inhSum;
			}

			double vmax_r = static_cast<double>(p->vMax[r]);
			flux *= vmax_r;

			// calculate derivative for every component jac_flux = dv/dy
			for (int comp = 0; comp < _nComp; ++comp)
			{
				// 1. Check if cIdx is substrate or inhibitor
				bool isSubstrate = false;
				int substrateIdx = -1;
				for (int sIdx = 0; sIdx < nSub; sIdx++)
				{
					if (comp == _idxSubstrate[r][sIdx])
					{
						isSubstrate = true;
						substrateIdx = sIdx;
						break;
					}
				}

				bool isInhibitor = false;
				int inhibitorIdx = -1;
				for (int iIdx = 0; iIdx < nInh; ++iIdx)
				{
					if (comp == _idxInhibitor[r][iIdx])
					{
						isInhibitor = true;
						inhibitorIdx = iIdx;
						break;
					}
				}

				double dvdy = 0.0;

				// case 1: comp is a substrat and not an inhibitor
				// dvds = vmax * prod_i!=s (y[i]/(kMM_rs + y[i])*(1+inhSum)) * df/ds
				// df/ds = ((kMM_rs + y[i])*(1+inhSum)) - y[i](1+inhSum) / ((kMM_rs + y[i])*(1+inhSum))^2
				//       = (kMM_rs) / ((kMM_rs + y[i])*(1+inhSum)^2)
				if (isSubstrate && !isInhibitor)
				{
					int s = comp; // comp is a substrat
					double kMM_rs = static_cast<double>(p->kMM[r]);
					if (!_oldInterface)
						kMM_rs = static_cast<double>(p->kMM[_nComp * r + s]);

					double inhSum = inhSumOfSub[substrateIdx];
					double dnom = (kMM_rs + y[s])* (kMM_rs + y[s])*(1+inhSum);

					double dsubFluxds = kMM_rs/dnom ;
					double factor = flux / substratFlux[substrateIdx];

					dvdy += factor * dsubFluxds;

				}

				//case 2: comp is a inhibitor and not a substrat
				// dvdI = vmax * sum_L(prod_!i=l (y[i]/(kMM_rs + y[i])*(1+inhSum)) * df/dI
				// L is the set of all substrates which are inhibited by I
				// df/dI =  - ((kMM_rs+y[s])/kI_ri)/((kMM_rs + y[s])*(1+inhSum))^2
				//		 =  - (y[s])/((kMM_rs + y[s])*(1+inhSum)^2 kI_ri)
				if (isInhibitor && !isSubstrate)
				{
					int i = comp;
					double kI_ri = static_cast<double>(p->kInhibit[_nComp * r + i]);

					for (int sIdx = 0; sIdx < nSub; sIdx++)
					{
						int s = _idxSubstrate[r][sIdx];
						double kMM_rs = static_cast<double>(p->kMM[r]);
						if (!_oldInterface)
							// kMM_rs = kMM[r][s]
							kMM_rs = static_cast<double>(p->kMM[_nComp * r + s]);
						double inhSum = inhSumOfSub[sIdx];

						double denom = (kMM_rs + y[s])*(1+inhSum)* (1 + inhSum);
						double dinhFluxdi = - (y[s])/ (denom * kI_ri);
						double factor = flux / substratFlux[sIdx];

						dvdy += factor * dinhFluxdi;
					}

				}

				if (isInhibitor && isSubstrate)
				{
					throw(InvalidParameterException("Inhibitor and substrat in the same reaction is not supported yet"));
				}

				if (std::abs(dvdy) < 1e-18)
					continue;

				// Add gradients to Jacobian
				RowIterator curJac = jac;
				for (int row = 0; row < _nComp; ++row, ++curJac)
				{
					const double colFactor = static_cast<double>(_stoichiometryBulk.native(row, r)) * factor;
					curJac[comp - static_cast<int>(row)] += colFactor * dvdy;
				}

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
