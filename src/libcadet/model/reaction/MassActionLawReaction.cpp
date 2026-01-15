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

#include "model/reaction/ReactionModelBase.hpp"
#include "model/ExternalFunctionSupport.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "model/Parameters.hpp"
#include "LocalVector.hpp"
#include "SimulationTypes.hpp"
#include "linalg/ActiveDenseMatrix.hpp"
#include "Memory.hpp"

#include <functional>
#include <algorithm>
#include <unordered_map>
#include <string>
#include <vector>

/*<codegen>
{
	"name": "MassActionLawParamHandler",
	"externalName": "ExtMassActionLawParamHandler",
	"parameters":
		[
			{ "type": "ScalarReactionDependentParameter", "varName": "kFwd", "confName": "MAL_KFWD"},
			{ "type": "ScalarReactionDependentParameter", "varName": "kBwd", "confName": "MAL_KBWD"}
		]
}
</codegen>*/

/* Parameter description
 ------------------------
 kFwd = Forward rate for reactions
 kBwd = Backward rate for reactions 
*/

namespace cadet
{

namespace model
{

inline const char* MassActionLawParamHandler::identifier() CADET_NOEXCEPT { return "MASS_ACTION_LAW"; }

inline bool MassActionLawParamHandler::validateConfig(unsigned int nReactions, unsigned int nComp, unsigned int const* nBoundStates)
{
	return true;
}

inline const char* ExtMassActionLawParamHandler::identifier() CADET_NOEXCEPT { return "EXT_MASS_ACTION_LAW"; }

inline bool ExtMassActionLawParamHandler::validateConfig(unsigned int nReactions, unsigned int nComp, unsigned int const* nBoundStates)
{
	return true;
}


namespace
{
	/**
	 * @brief Registers a matrix-valued parameter (row-major storage) with bound states as rows
	 * @details The matrix-valued parameter has as many rows as there are bound states in the system.
	 * @param [in,out] parameters Parameter map
	 * @param [in] unitOpIdx Unit operation id
	 * @param [in] parTypeIdx Particle type index
	 * @param [in] paramName Name of the parameter
	 * @param [in] mat Matrix to register
	 * @param [in] nComp Number of components
	 * @param [in] boundOffset Array with offsets to bound states of a specific component
	 */
	inline void registerBoundStateRowMatrix(std::unordered_map<ParameterId, active*>& parameters, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx, const std::string& paramName, cadet::linalg::ActiveDenseMatrix& mat,
		unsigned int nComp, unsigned int const* boundOffset)
	{
		const cadet::StringHash hashName = cadet::hashStringRuntime(paramName);
		cadet::registerParam2DArray(parameters, mat.data(), mat.elements(), [=](bool multi, unsigned int row, unsigned int col)
			{
				const unsigned int comp = std::lower_bound(boundOffset, boundOffset + nComp, row) - boundOffset;
				const unsigned int bnd = row - boundOffset[comp];
				return cadet::makeParamId(hashName, unitOpIdx, comp, parTypeIdx, bnd, col, cadet::SectionIndep);
			},
			mat.columns()
		);
	}

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

	/**
	 * @brief Reads reaction rate exponents and registers them
	 * @details Reads a matrix-valued parameter into the given pre-allocated matrix and registers the
	 *          parameters in the map. If @p boundOffset is @c nullptr, the matrix has as many columns
	 *          as there are components. Otherwise, the matrix has as many columns as there are bound states.
	 * @param [in] paramProvider Parameter provider
	 * @param [in,out] parameters Parameter map
	 * @param [in] unitOpIdx Unit operation id
	 * @param [in] parTypeIdx Particle type index
	 * @param [in] paramName Name of the parameter
	 * @param [out] mat Matrix that holds the values (pre-allocated)
	 * @param [in] nComp Number of components
	 * @param [in] boundOffset Array with offsets to bound states of a specific component, or @c nullptr
	 */
	inline void readAndRegisterExponents(cadet::IParameterProvider& paramProvider, std::unordered_map<ParameterId, active*>& parameters, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx,
		const std::string& paramName, cadet::linalg::ActiveDenseMatrix& mat, unsigned int nComp, unsigned int const* boundOffset)
	{
		if (paramProvider.exists(paramName))
		{
			const std::vector<double> s = paramProvider.getDoubleArray(paramName);
			if (static_cast<int>(s.size()) != mat.elements())
				throw InvalidParameterException("Expected " + paramName + " to be of size " + std::to_string(mat.elements()) + "");

			std::copy(s.begin(), s.end(), mat.data());
		}
		else
			mat.setAll(0.0);

		if (boundOffset)
			registerBoundStateRowMatrix(parameters, unitOpIdx, parTypeIdx, paramName, mat, nComp, boundOffset);
		else
			registerCompRowMatrix(parameters, unitOpIdx, parTypeIdx, paramName, mat);
	}

	/**
	 * @brief Calculate gradient of reaction rate (flux) in liquid phase
	 * @param [out] fluxGrad Array that holds the gradient of the reaction rate
	 * @param [in] r Index of the reaction
	 * @param [in] nComp Number of components
	 * @param [in] rate Rate constant
	 * @param [in] exponents Matrix with exponents in the rate law
	 * @param [in] y Array of concentrations
	 */
	inline void fluxGrad(double* fluxGrad, unsigned int r, unsigned int nComp, double rate, const cadet::linalg::ActiveDenseMatrix& exponents, double const* y)
	{
		for (unsigned int c = 0; c < nComp; ++c)
		{
			if (cadet_unlikely(exponents.native(c, r) != 0.0))
				fluxGrad[c] = rate;
			else
				fluxGrad[c] = 0.0;
		}

		// Calculate gradient
		for (unsigned int c = 0; c < nComp; ++c)
		{
			if (cadet_unlikely(exponents.native(c, r) != 0.0))
			{
				const double exponentValue = static_cast<double>(exponents.native(c, r));
				const double v = pow(y[c], exponentValue);
				for (unsigned int j = 0; j < c; ++j)
					fluxGrad[j] *= v;

				fluxGrad[c] *= exponentValue * pow(y[c], exponentValue - 1.0);

				for (unsigned int j = c + 1; j < nComp; ++j)
					fluxGrad[j] *= v;

			}
		}
	}
	inline void fluxGradCombined(double* fluxGrad, unsigned int r, unsigned int nComp, unsigned int nTotalBoundStates, double rate,
		const cadet::linalg::ActiveDenseMatrix& expLiquid, const cadet::linalg::ActiveDenseMatrix& expSolid, double const* yLiquid, double const* ySolid)
	{

	}

	template <typename num_t>
	inline num_t rateConstantOrZero(const num_t& rate, unsigned int idxReaction, const cadet::linalg::ActiveDenseMatrix& exponents, unsigned int nComp)
	{
		for (unsigned int c = 0; c < nComp; ++c)
		{
			if (cadet_unlikely(exponents.native(c, idxReaction) != 0.0))
				return rate;
		}
		return 0.0;
	}

	template <typename num_t>
	inline num_t rateConstantOrZero(const num_t& rate, unsigned int idxReaction, const cadet::linalg::ActiveDenseMatrix& expLiquid, const cadet::linalg::ActiveDenseMatrix& expSolid, unsigned int nComp, unsigned int nTotalBoundStates)
	{
		for (unsigned int c = 0; c < nComp; ++c)
		{
			if (cadet_unlikely(expLiquid.native(c, idxReaction) != 0.0))
				return rate;
		}
		for (unsigned int c = 0; c < nTotalBoundStates; ++c)
		{
			if (cadet_unlikely(expSolid.native(c, idxReaction) != 0.0))
				return rate;
		}
		return 0.0;
	}
}

template <class ParamHandler_t>
class MassActionLawReactionBase : public DynamicReactionModelBase
{
public:

	MassActionLawReactionBase() { }
	virtual ~MassActionLawReactionBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }
	virtual const char* name() const CADET_NOEXCEPT { return ParamHandler_t::identifier(); }

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { _paramHandler.setExternalFunctions(extFuns, size); }
	virtual bool dependsOnTime() const CADET_NOEXCEPT { return ParamHandler_t::dependsOnTime(); }
	virtual bool requiresWorkspace() const CADET_NOEXCEPT { return true; }
	virtual unsigned int workspaceSize(unsigned int nComp, unsigned int totalNumBoundStates, unsigned int const* nBoundStates) const CADET_NOEXCEPT
	{
		return _paramHandler.cacheSize(maxNumReactions(), nComp, totalNumBoundStates) + std::max(maxNumReactions() * sizeof(active), 2 * (_nComp + totalNumBoundStates) * sizeof(double));
	}

	virtual bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
	{
		DynamicReactionModelBase::configureModelDiscretization(paramProvider, nComp, nBound, boundOffset);
	
		if (paramProvider.exists("MAL_STOICHIOMETRY"))
		{
			const std::size_t numElements = paramProvider.numElements("MAL_STOICHIOMETRY");
			if (numElements % nComp != 0)
				throw InvalidParameterException("Size of field MAL_STOICHIOMETRY must be a positive multiple of NCOMP (" + std::to_string(nComp) + ")");

			const unsigned int nReactions = numElements / nComp;

			_stoichiometry.resize(nComp, nReactions);
			_expFwd.resize(nComp, nReactions);
			_expBwd.resize(nComp, nReactions);
		}

		if (!nBound || !boundOffset)
			return true;

		return true;
	}

	virtual unsigned int numReactions() const CADET_NOEXCEPT { return _stoichiometry.columns(); }
	virtual unsigned int numReactionsLiquid() const CADET_NOEXCEPT { return 0; }
	virtual unsigned int numReactionsCombined() const CADET_NOEXCEPT { return 0; }

	CADET_DYNAMICREACTIONMODEL_BOILERPLATE

protected:
	ParamHandler_t _paramHandler; //!< Handles parameters and their dependence on external functions

	linalg::ActiveDenseMatrix _stoichiometry;
	linalg::ActiveDenseMatrix _expFwd;
	linalg::ActiveDenseMatrix _expBwd;
	
	
	inline int maxNumReactions() const CADET_NOEXCEPT
	{
		return _stoichiometry.columns();
	}

	virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
	{
		_paramHandler.configure(paramProvider, maxNumReactions(), _nComp, _nBoundStates);
		_paramHandler.registerParameters(_parameters, unitOpIdx, parTypeIdx, _nComp, _nBoundStates);

		if ((_stoichiometry.columns() > 0) && ((static_cast<int>(_paramHandler.kFwd().size()) < _stoichiometry.columns()) || (static_cast<int>(_paramHandler.kBwd().size()) < _stoichiometry.columns())))
			throw InvalidParameterException("MAL_KFWD and MAL_KBWD have to have the same size (number of reactions)");
		
		if (paramProvider.exists("MAL_STOICHIOMETRY"))
		{
			const std::vector<double> s = paramProvider.getDoubleArray("MAL_STOICHIOMETRY");

			if (static_cast<int>(s.size()) != _stoichiometry.elements())
				throw InvalidParameterException("MAL_STOICHIOMETRY has changed size (number of reactions changed)");

			std::copy(s.begin(), s.end(), _stoichiometry.data());
		}
		registerCompRowMatrix(_parameters, unitOpIdx, parTypeIdx, "MAL_STOICHIOMETRY", _stoichiometry);

		if (paramProvider.exists("MAL_EXPONENTS_FWD"))
		{
			const std::vector<double> s = paramProvider.getDoubleArray("MAL_EXPONENTS_FWD");
			if (static_cast<int>(s.size()) != _stoichiometry.elements())
				throw InvalidParameterException("Expected MAL_EXPONENTS_FWD and MAL_STOICHIOMETRY to be of the same size (" + std::to_string(_stoichiometry.elements()) + ")");

			std::copy(s.begin(), s.end(), _expFwd.data());
		}
		else
		{
			// Obtain default values from stoichiometry
			std::transform(_stoichiometry.data(), _stoichiometry.data() + _stoichiometry.elements(), _expFwd.data(), [](const active& v) { return std::max(0.0, -static_cast<double>(v)); });
		}
		registerCompRowMatrix(_parameters, unitOpIdx, parTypeIdx, "MAL_EXPONENTS_FWD", _expFwd);

		if (paramProvider.exists("MAL_EXPONENTS_BWD"))
		{
			const std::vector<double> s = paramProvider.getDoubleArray("MAL_EXPONENTS_BWD");
			if (static_cast<int>(s.size()) != _stoichiometry.elements())
				throw InvalidParameterException("Expected MAL_EXPONENTS_BWD and MAL_STOICHIOMETRY to be of the same size (" + std::to_string(_stoichiometry.elements()) + ")");

			std::copy(s.begin(), s.end(), _expBwd.data());
		}
		else
		{
			// Obtain default values from stoichiometry
			std::transform(_stoichiometry.data(), _stoichiometry.data() + _stoichiometry.elements(), _expBwd.data(), [](const active& v) { return std::max(0.0, static_cast<double>(v)); });
		}
		registerCompRowMatrix(_parameters, unitOpIdx, parTypeIdx, "MAL_EXPONENTS_BWD", _expBwd);
		

		return true;
	}

	template <typename StateType, typename ResidualType, typename ParamType, typename FactorType>
	int residualFluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
		const unsigned int nStates, StateType const* y, ResidualType* res, const FactorType& factor, LinearBufferAllocator workSpace) const
	{
		typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Calculate fluxes
		typedef typename DoubleActivePromoter<StateType, ParamType>::type flux_t;
		BufferedArray<flux_t> fluxes = workSpace.array<flux_t>(maxNumReactions());

		for (int r = 0; r < _stoichiometry.columns(); ++r)
		{
			flux_t fwd = rateConstantOrZero(static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->kFwd[r]), r, _expFwd, nStates);

			for (int c = 0; c < nStates; ++c)
			{
				if (_expFwd.native(c, r) != 0.0)
				{
					if (static_cast<double>(y[c]) > 0.0)
						fwd *= pow(static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(y[c]),
							static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(_expFwd.native(c, r)));
					else
					{
						fwd *= 0.0;
						break;
					}
				}
			}

			flux_t bwd = rateConstantOrZero(static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(p->kBwd[r]), r, _expBwd, nStates);
			for (int c = 0; c < nStates; ++c)
			{	

				if (_expBwd.native(c, r) != 0.0)
				{
					if (static_cast<double>(y[c]) > 0.0)
						bwd *= pow(static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(y[c]),
							static_cast<typename DoubleActiveDemoter<flux_t, active>::type>(_expBwd.native(c, r)));
					else
					{
						bwd *= 0.0;
						break;
					}
				}
			}

			fluxes[r] = fwd - bwd;
		}

		// Add reaction terms to residual
		_stoichiometry.multiplyVector(static_cast<flux_t*>(fluxes), factor, res);

		return 0;
	}


	template <typename RowIterator>
	void jacobianFluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, const unsigned int nState, double const* y, double factor, const RowIterator& jac, LinearBufferAllocator workSpace) const
	{
		typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		BufferedArray<double> fluxes = workSpace.array<double>(2 * nState);
		double* const fluxGradFwd = static_cast<double*>(fluxes);
		double* const fluxGradBwd = fluxGradFwd + nState;
		
		for (int r = 0; r < _stoichiometry.columns(); ++r)
		{
			// Calculate gradients of forward and backward fluxes
			double kFwd = static_cast<double>(p->kFwd[r]);
			double kBwd = static_cast<double>(p->kBwd[r]);
			

			fluxGrad(fluxGradFwd, r, nState, kFwd, _expFwd, y);
			fluxGrad(fluxGradBwd, r, nState, kBwd, _expBwd, y);

			// Add gradients to Jacobian
			RowIterator curJac = jac;
			for (int row = 0; row < nState; ++row, ++curJac)
			{
				const double colFactor = static_cast<double>(_stoichiometry.native(row, r)) * factor;
				for (int col = 0; col < nState; ++col)
					curJac[col - static_cast<int>(row)] += colFactor * (fluxGradFwd[col] - fluxGradBwd[col]);
			}
		}
	}

		template <typename StateType, typename ResidualType, typename ParamType>
	int residualCombinedImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
		StateType const* yLiquid, StateType const* ySolid, ResidualType* resLiquid, ResidualType* resSolid, double factor, LinearBufferAllocator workSpace) const
	{
		return 0;
	}
	
	template <typename RowIteratorLiquid, typename RowIteratorSolid>
	void jacobianCombinedImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* yLiquid, double const* ySolid, double factor, const RowIteratorLiquid& jacLiquid, const RowIteratorSolid& jacSolid, LinearBufferAllocator workSpace) const
	{

	}


};

typedef MassActionLawReactionBase<MassActionLawParamHandler> MassActionLawReaction;
typedef MassActionLawReactionBase<ExtMassActionLawParamHandler> ExternalMassActionLawReaction;

namespace reaction
{
	void registerMassActionLawReaction(std::unordered_map<std::string, std::function<model::IDynamicReactionModel*()>>& reactions)
	{
		reactions[MassActionLawReaction::identifier()] = []() { return new MassActionLawReaction(); };
		reactions[ExternalMassActionLawReaction::identifier()] = []() { return new ExternalMassActionLawReaction(); };
	}
}  // namespace reaction

}  // namespace model

}  // namespace cadet
