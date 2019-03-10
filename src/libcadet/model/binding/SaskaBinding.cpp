// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2019: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
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

#include <functional>
#include <unordered_map>
#include <string>
#include <vector>

namespace cadet
{

namespace model
{

/*<codegen>
{
	"name": "SaskaParamHandler",
	"externalName": "ExtSaskaParamHandler",
	"parameters":
		[
			{ "type": "ScalarComponentDependentParameter", "varName": "h", "confName": "SASKA_H"},
			{ "type": "ComponentDependentComponentVectorParameter", "varName": "k", "confName": "SASKA_K"}
		]
}
</codegen>*/

/* Parameter description
 ------------------------
 h = Henry coefficient
 k = Quadratic factor
*/

inline const char* SaskaParamHandler::identifier() CADET_NOEXCEPT { return "SASKA"; }

inline bool SaskaParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_h.size() != _k.slices()) || (_k.size() != nComp * nComp) || (_h.size() < nComp))
		throw InvalidParameterException("SASKA_K has to have NCOMP^2 and SASKA_H NCOMP many elements");

	return true;
}

inline const char* ExtSaskaParamHandler::identifier() CADET_NOEXCEPT { return "EXT_SASKA"; }

inline bool ExtSaskaParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if ((_h.size() != _k.slices()) || (_k.size() != nComp * nComp) || (_h.size() < nComp))
		throw InvalidParameterException("EXT_SASKA_K has to have NCOMP^2 and EXT_SASKA_H NCOMP many elements");

	return true;
}


/**
 * @brief Defines the Saska binding model
 * @details Implements the Saska adsorption model: \f[ \begin{align} 
 *              \frac{\mathrm{d} q_i}{\mathrm{d} t} = H_i c_{p,i} + \sum_{j=1}^{N_{\text{comp}}} k_{ij} c_{p,i} c_{p,j} - q_i,
 *          \end{align} \f]
 *          where @f$ H_i @f$ denotes the Henry coefficient and @f$ k_{ij} @f$ the quadratic factor.
 *          Multiple bound states are not supported. However, the quadratic factors use the bound phase parameter index for @f$ j @f$. 
 *          Components without bound state (i.e., non-binding components) are supported.
 *          
 *          See @cite Saska1992.
 * @tparam ParamHandler_t Type that can add support for external function dependence
 */
template <class ParamHandler_t>
class SaskaBindingBase : public PureBindingModelBase
{
public:

	SaskaBindingBase() { }
	virtual ~SaskaBindingBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }
	virtual const char* name() const CADET_NOEXCEPT { return ParamHandler_t::identifier(); }

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { _paramHandler.setExternalFunctions(extFuns, size); }
	virtual bool dependsOnTime() const CADET_NOEXCEPT { return ParamHandler_t::dependsOnTime(); }
	virtual bool requiresWorkspace() const CADET_NOEXCEPT { return hasAlgebraicEquations() || ParamHandler_t::requiresWorkspace(); }

	virtual void timeDerivativeAlgebraicResidual(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* y, double* dResDt, void* workSpace) const
	{
		if (!hasAlgebraicEquations())
			return;

		if (!ParamHandler_t::dependsOnTime())
			return;

		// Update external function and compute time derivative of parameters
		const typename ParamHandler_t::params_t& p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);
		const typename ParamHandler_t::params_t dpDt = _paramHandler.updateTimeDerivative(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Pointer to first component in liquid phase
		double const* yCp = y - _nComp;

		// Protein equations: dq_i / dt - ( H_i c_{p,i} + \sum_j k_{ij} c_{p,i} c_{p,j} - q_i ) == 0
		//               <=>  dq_i / dt == H_i c_{p,i} + \sum_j k_{ij} c_{p,i} c_{p,j} - q_i
		unsigned int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			const active h = dpDt.h[i];
			active const* const kSlice = dpDt.k[i];

			// Residual
			dResDt[bndIdx] = -static_cast<double>(h) * yCp[i];
			for (int j = 0; j < _nComp; ++j)
				dResDt[bndIdx] -= static_cast<double>(kSlice[j]) * yCp[i] * yCp[j];

			// Next bound component
			++bndIdx;
		}
	}

	CADET_PUREBINDINGMODELBASE_BOILERPLATE

protected:
	ParamHandler_t _paramHandler; //!< Handles parameters and their dependence on external functions

	virtual unsigned int paramCacheSize(unsigned int nComp, unsigned int totalNumBoundStates, unsigned int const* nBoundStates) const CADET_NOEXCEPT
	{
		return _paramHandler.cacheSize(nComp, totalNumBoundStates, nBoundStates);
	}

	virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
	{
		// Read parameters
		_paramHandler.configure(paramProvider, _nComp, _nBoundStates);

		// Register parameters
		_paramHandler.registerParameters(_parameters, unitOpIdx, parTypeIdx, _nComp, _nBoundStates);

		return true;
	}

	template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
	int residualImpl(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, const ColumnPosition& colPos,
		StateType const* y, CpStateType const* yCp, double const* yDot, ResidualType* res, void* workSpace) const
	{
		const typename ParamHandler_t::params_t& p = _paramHandler.update(static_cast<double>(t), secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Protein equations: dq_i / dt - ( H_i c_{p,i} + \sum_j k_{ij} c_{p,i} c_{p,j} - q_i ) == 0
		//               <=>  dq_i / dt == H_i c_{p,i} + \sum_j k_{ij} c_{p,i} c_{p,j} - q_i
		unsigned int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			const active h = p.h[i];
			active const* const kSlice = p.k[i];

			// Residual
			res[bndIdx] = -static_cast<ParamType>(h) * yCp[i] + y[bndIdx];
			for (int j = 0; j < _nComp; ++j)
				res[bndIdx] -= static_cast<ParamType>(kSlice[j]) * yCp[i] * yCp[j];

			// Add time derivative if necessary
			if (_kineticBinding && yDot)
			{
				res[bndIdx] += timeFactor * yDot[bndIdx];
			}

			// Next bound component
			++bndIdx;
		}

		return 0;
	}

	template <typename RowIterator>
	void jacobianImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* y, double const* yCp, int offsetCp, RowIterator jac, void* workSpace) const
	{
		const typename ParamHandler_t::params_t& p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Protein equations: dq_i / dt - ( H_i c_{p,i} + \sum_j k_{ij} c_{p,i} c_{p,j} - q_i ) == 0
		//               <=>  dq_i / dt - ( H_i c_{p,i} + c_{p,i} * \sum_j k_{ij} c_{p,j} - q_i ) == 0
		int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			const double h = static_cast<double>(p.h[i]);
			active const* const kSlice = p.k[i];

			// dres_i / dc_{p,j}
			for (int j = 0; j < _nComp; ++j)
				jac[j - bndIdx - offsetCp] = -static_cast<double>(kSlice[j]) * yCp[i];
				// Getting to c_{p,i}: -bndIdx takes us to q_0, another -offsetCp to c_{p,0} and a +j to c_{p,j}.
				//                     This means jac[j - bndIdx - offsetCp] corresponds to c_{p,j}.

			// dres_i / dc_{p,i}
			for (int j = 0; j < _nComp; ++j)
				jac[i - bndIdx - offsetCp] -= static_cast<double>(kSlice[j]) * yCp[j];
			
			jac[i - bndIdx - offsetCp] -= h;

			// dres_i / dq_i
			jac[0] = 1.0;

			// Advance to next equation and Jacobian row
			++bndIdx;
			++jac;
		}
	}
};


typedef SaskaBindingBase<SaskaParamHandler> SaskaBinding;
typedef SaskaBindingBase<ExtSaskaParamHandler> ExternalSaskaBinding;

namespace binding
{
	void registerSaskaModel(std::unordered_map<std::string, std::function<model::IBindingModel*()>>& bindings)
	{
		bindings[SaskaBinding::identifier()] = []() { return new SaskaBinding(); };
		bindings[ExternalSaskaBinding::identifier()] = []() { return new ExternalSaskaBinding(); };
	}
}  // namespace binding

}  // namespace model

}  // namespace cadet
