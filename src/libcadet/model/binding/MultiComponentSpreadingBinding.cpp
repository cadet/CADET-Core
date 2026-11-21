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
#include "SlicedVector.hpp"
#include "LocalVector.hpp"
#include "SimulationTypes.hpp"

#include <vector>
#include <unordered_map>
#include <functional>

/*<codegen>
{
	"name": "SpreadingParamHandler",
	"externalName": "ExtSpreadingParamHandler",
	"parameters":
		[
			{ "type": "ComponentBoundStateMajorDependentParameter", "varName": "kA", "confName": "MCSPR_KA"},
			{ "type": "ComponentBoundStateMajorDependentParameter", "varName": "kD", "confName": "MCSPR_KD"},
			{ "type": "ComponentBoundStateMajorDependentParameter", "varName": "qMax", "confName": "MCSPR_QMAX"},
			{ "type": "ScalarComponentDependentParameter", "varName": "k12", "confName": "MCSPR_K12"},
			{ "type": "ScalarComponentDependentParameter", "varName": "k21", "confName": "MCSPR_K21"}
		]
}
</codegen>*/

/* Parameter description
 ------------------------
 kA = Adsorption rate in state-major ordering
 kD = Desorption rate in state-major ordering
 qMax = Capacity in state-major ordering
 k12 = Transition rate from state @f$ q_i^A @f$ to state @f$ q_i^B @f$
 k21 = Transition rate from state @f$ q_i^B @f$ to state @f$ q_i^A @f$
*/

namespace cadet
{

namespace model
{

inline const char* SpreadingParamHandler::identifier() CADET_NOEXCEPT { return "MULTI_COMPONENT_SPREADING"; }

inline bool SpreadingParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_kA.size() != _kD.size()) || (_kA.size() != _qMax.size()) || (_kA.size() < nComp * 2))
		throw InvalidParameterException("MCSPR_KA, MCSPR_KD, and MCSPR_QMAX have to have the same size");
	if ((_k12.size() != _k21.size()) || (_k12.size() < nComp))
		throw InvalidParameterException("MCSPR_K12 and MCSPR_K21 have to have the same size (number of components)");

	return true;
}

inline const char* ExtSpreadingParamHandler::identifier() CADET_NOEXCEPT { return "EXT_MULTI_COMPONENT_SPREADING"; }

inline bool ExtSpreadingParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_kA.size() != _kD.size()) || (_kA.size() != _qMax.size()) || (_kA.size() < nComp * 2))
		throw InvalidParameterException("EXT_MCSPR_KA, EXT_MCSPR_KD, and EXT_MCSPR_QMAX have to have the same size");
	if ((_k12.size() != _k21.size()) || (_k12.size() < nComp))
		throw InvalidParameterException("EXT_MCSPR_K12 and EXT_MCSPR_K21 have to have the same size (number of components)");

	return true;
}


/**
 * @brief Defines the multi component spreading binding model
 * @details Implements the multi component spreading adsorption model: \f[ \begin{align} 
 *                \frac{\mathrm{d} q_i^A}{\mathrm{d} t} &= \left( k_a^A\: c_{p,i} - k_{AB} q_i^A \right) q_{\text{max},i}^A \left( 1 - \sum_{j=1}^{N_{\text{comp}}} \frac{q_j^A}{q_{\text{max},j}^A} - \sum_{j=1}^{N_{\text{comp}}} \frac{q_j^B}{q_{\text{max},j}^B} \right) - k_d^A q_i^A + k_{BA} q_i^B \\
 *                \frac{\mathrm{d} q_i^B}{\mathrm{d} t} &= \left( k_a^B\: c_{p,i} + k_{AB} q_i^A \right) q_{\text{max},i}^A \left( 1 - \sum_{j=1}^{N_{\text{comp}}} \frac{q_j^A}{q_{\text{max},j}^A} - \sum_{j=1}^{N_{\text{comp}}} \frac{q_j^B}{q_{\text{max},j}^B} \right) - \left( k_d^B + k_{BA} \right) q_i^B
 *          \end{align} \f]
 *          Here, a second bound state $q_i^B$ is added to the Langmuir model and the exchange between the two bound 
 *          states $q_i^A$ and $q_i^B$ is allowed. The second bound state may correspond to a different orientation 
 *          on the surface or a different folding state of the molecule.
 *          While components without bound state (i.e., non-binding components) are supported, all other components must have
 *          @c 2 bound states.
 *          
 *          Internal state vector layout is component-major. First, all bound states of component 0 are placed, then all bound states of component 1, etc.
 * @tparam ParamHandler_t Type that can add support for external function dependence
 */
template <class ParamHandler_t>
class MultiComponentSpreadingBindingBase : public ParamHandlerBindingModelBase<ParamHandler_t>
{
public:

	MultiComponentSpreadingBindingBase() { }
	virtual ~MultiComponentSpreadingBindingBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }

	virtual bool configureModelDiscretization(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
	{
		const bool res = BindingModelBase::configureModelDiscretization(paramProvider, nComp, nBound, boundOffset);

		const unsigned int numSlices = 2;
		for (unsigned int i = 0; i < nComp; ++i)
		{
			if (nBound[i] == 0)
				continue;

			if (nBound[i] != numSlices)
				throw InvalidParameterException("Multi component spreading binding model requires exactly two bound states for all (binding) components");
		}

		_numBindingComp = numBindingComponents(_nBoundStates, _nComp);

		// Allocate space for parameters
		_paramHandler.reserve(nComp * numSlices, numSlices, nComp, nBound);

		return res;
	}

	virtual bool hasSalt() const CADET_NOEXCEPT { return false; }
	virtual bool supportsMultistate() const CADET_NOEXCEPT { return true; }
	virtual bool supportsNonBinding() const CADET_NOEXCEPT { return true; }

	CADET_BINDINGMODELBASE_BOILERPLATE

protected:
	using ParamHandlerBindingModelBase<ParamHandler_t>::_paramHandler;
	using ParamHandlerBindingModelBase<ParamHandler_t>::_reactionQuasistationarity;
	using ParamHandlerBindingModelBase<ParamHandler_t>::_nComp;
	using ParamHandlerBindingModelBase<ParamHandler_t>::_nBoundStates;

	unsigned int _numBindingComp; //!< Number of binding components

	virtual bool implementsAnalyticJacobian() const CADET_NOEXCEPT { return true; }

	template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
	int fluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, StateType const* y,
		CpStateType const* yCp, ResidualType* res, LinearBufferAllocator workSpace) const
	{
		typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Ordering of the states is component-major, that is, q_{0,A}, q_{0,B}, q_{1,A}, q_{1,B}, ...
		// Hence, the stride from one binding component to the next is 2 and state B of a component
		// sits right behind its state A. The fluxes use the same ordering.

		ResidualType qSum = 1.0;
		for (int j = 0; j < 2; ++j)
		{
			active const* const localQmax = p->qMax[j];
			int bndIdx = 0;
			for (int i = 0; i < _nComp; ++i)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[i] == 0)
					continue;

				qSum -= y[2 * bndIdx + j] / static_cast<ParamType>(localQmax[i]);

				// Next bound component
				++bndIdx;
			}
		}

		active const* const localKaA = p->kA[0];
		active const* const localKdA = p->kD[0];
		active const* const localKaB = p->kA[1];
		active const* const localKdB = p->kD[1];
		// Still use q_{max}^A for the second bound state (B)
		active const* const localQmax = p->qMax[0];

		int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			// Residual of the first bound state (A)
			res[2 * bndIdx] = static_cast<ParamType>(localKdA[i]) * y[2 * bndIdx] - static_cast<ParamType>(p->k21[i]) * y[2 * bndIdx + 1]
					- (static_cast<ParamType>(localKaA[i]) * yCp[i] - static_cast<ParamType>(p->k12[i]) * y[2 * bndIdx]) * static_cast<ParamType>(localQmax[i]) * qSum;

			// Residual of the second bound state (B)
			res[2 * bndIdx + 1] = (static_cast<ParamType>(localKdB[i]) + static_cast<ParamType>(p->k21[i])) * y[2 * bndIdx + 1]
					- (static_cast<ParamType>(localKaB[i]) * yCp[i] + static_cast<ParamType>(p->k12[i]) * y[2 * bndIdx]) * static_cast<ParamType>(localQmax[i]) * qSum;

			// Next bound component
			++bndIdx;
		}

		return 0;
	}

	template <typename RowIterator>
	void jacobianImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* y, double const* yCp, int offsetCp, RowIterator jac, LinearBufferAllocator workSpace) const
	{
		typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Ordering of the states is component-major, that is, q_{0,A}, q_{0,B}, q_{1,A}, q_{1,B}, ...
		// Ordering of the Jacobian rows (i.e., of the fluxes) is the same. The row of q_{i,A} is the
		// bound state 2 * bndIdx and the row of q_{i,B} is the bound state 2 * bndIdx + 1.

		double qSum = 1.0;
		for (int j = 0; j < 2; ++j)
		{
			active const* const localQmax = p->qMax[j];
			int bndIdx = 0;
			for (int i = 0; i < _nComp; ++i)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[i] == 0)
					continue;

				qSum -= y[2 * bndIdx + j] / static_cast<double>(localQmax[i]);

				// Next bound component
				++bndIdx;
			}
		}

		active const* const localKaA = p->kA[0];
		active const* const localKdA = p->kD[0];
		active const* const localKaB = p->kA[1];
		active const* const localKdB = p->kD[1];
		active const* const qMax1 = p->qMax[0];
		active const* const qMax2 = p->qMax[1];

		int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			const double k12 = static_cast<double>(p->k12[i]);
			const double k21 = static_cast<double>(p->k21[i]);

			// First state (A), which is the row of bound state 2 * bndIdx
			{
				const double ka = static_cast<double>(localKaA[i]);
				const double kd = static_cast<double>(localKdA[i]);

				// dres_i^A / dc_{p,i}
				jac[i - 2 * bndIdx - offsetCp] = -ka * static_cast<double>(qMax1[i]) * qSum;
				// Getting to c_{p,i}: -2 * bndIdx takes us to q_{0,A}, another -offsetCp to c_{p,0} and a +i to c_{p,i}.
				//                     This means jac[i - 2 * bndIdx - offsetCp] corresponds to c_{p,i}.

				// Fill dres_i^A / dq_j^A and dres_i^A / dq_j^B
				const double factor = (ka * yCp[i] - y[2 * bndIdx] * k12) * static_cast<double>(qMax1[i]);
				int bndIdx2 = 0;
				for (int j = 0; j < _nComp; ++j)
				{
					// Skip components without bound states (bound state index bndIdx2 is not advanced)
					if (_nBoundStates[j] == 0)
						continue;

					// dres_i^A / dq_j^A
					jac[2 * bndIdx2 - 2 * bndIdx] = factor / static_cast<double>(qMax1[j]);
					// dres_i^A / dq_j^B
					jac[2 * bndIdx2 + 1 - 2 * bndIdx] = factor / static_cast<double>(qMax2[j]);
					// Getting to q_{j,A}: -2 * bndIdx takes us to q_{0,A}, another +2 * bndIdx2 to q_{j,A}
					//                     and one more step to q_{j,B}.

					++bndIdx2;
				}

				// Add to dres_i^A / dq_i^A
				jac[0] += kd + k12 * static_cast<double>(qMax1[i]) * qSum; // last summand by product rule

				// Add to dres_i^A / dq_i^B
				jac[1] -= k21;

				// Advance to next flux and Jacobian row
				++jac;
			}

			// Second state (B), which is the row of bound state 2 * bndIdx + 1
			{
				const double ka = static_cast<double>(localKaB[i]);
				const double kd = static_cast<double>(localKdB[i]);

				// dres_i^B / dc_{p,i}
				jac[i - 2 * bndIdx - 1 - offsetCp] = -ka * static_cast<double>(qMax1[i]) * qSum;
				// Getting to c_{p,i}: -2 * bndIdx - 1 takes us to q_{0,A}, another -offsetCp to c_{p,0} and a +i to c_{p,i}.

				// Fill dres_i^B / dq_j^A and dres_i^B / dq_j^B
				const double factor = (ka * yCp[i] + y[2 * bndIdx] * k12) * static_cast<double>(qMax1[i]);
				int bndIdx2 = 0;
				for (int j = 0; j < _nComp; ++j)
				{
					// Skip components without bound states (bound state index bndIdx2 is not advanced)
					if (_nBoundStates[j] == 0)
						continue;

					// dres_i^B / dq_j^A
					jac[2 * bndIdx2 - 2 * bndIdx - 1] = factor / static_cast<double>(qMax1[j]);
					// dres_i^B / dq_j^B
					jac[2 * bndIdx2 - 2 * bndIdx] = factor / static_cast<double>(qMax2[j]);

					++bndIdx2;
				}

				// Add to dres_i^B / dq_i^B
				jac[0] += kd + k21;

				// Add to dres_i^B / dq_i^A
				jac[-1] -= k12 * static_cast<double>(qMax1[i]) * qSum;  // last summand by product rule

				// Advance to next flux and Jacobian row
				++jac;
			}

			// Next bound component
			++bndIdx;
		}
	}
};


typedef MultiComponentSpreadingBindingBase<SpreadingParamHandler> MultiComponentSpreadingBinding;
typedef MultiComponentSpreadingBindingBase<ExtSpreadingParamHandler> ExternalMultiComponentSpreadingBinding;

namespace binding
{
	void registerMultiComponentSpreadingModel(std::unordered_map<std::string, std::function<model::IBindingModel*()>>& bindings)
	{
		bindings[MultiComponentSpreadingBinding::identifier()] = []() { return new MultiComponentSpreadingBinding(); };
		bindings[ExternalMultiComponentSpreadingBinding::identifier()] = []() { return new ExternalMultiComponentSpreadingBinding(); };
	}
}  // namespace binding

}  // namespace model

}  // namespace cadet
