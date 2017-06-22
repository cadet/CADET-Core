// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2017: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "model/binding/BindingModelBase.hpp"
#include "model/binding/ExternalFunctionSupport.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "ParamReaderHelper.hpp"
#include "SlicedVector.hpp"

#include <functional>
#include <unordered_map>
#include <string>
#include <vector>

namespace cadet
{

namespace model
{

/**
 * @brief Handles Saska binding model parameters that do not depend on external functions
 */
struct SaskaParamHandler : public BindingParamHandlerBase
{
	static const char* identifier() { return "SASKA"; }

	/**
	 * @brief Reads parameters and verifies them
	 * @details See IBindingModel::configure() for details.
	 * @param [in] paramProvider IParameterProvider used for reading parameters
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 * @return @c true if the parameters were read and validated successfully, otherwise @c false
	 */
	inline bool configure(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBoundStates)
	{
		// Allocate space for parameters
		h.reserve(nComp);
		k.reserve(nComp * nComp, nComp);

		// Read parameters
		readParameterMatrix(h, paramProvider, "SASKA_H", nComp, 1);
		readBoundStateDependentParameter<util::SlicedVector<active>, active>(k, paramProvider, "SASKA_K", nComp, nComp);

		// Check parameters
		if ((h.size() != k.slices()) || (k.size() != nComp * nComp) || (h.size() < nComp))
			throw InvalidParameterException("SASKA_K has to have NCOMP^2 and SASKA_H NCOMP many elements");

		return true;
	}

	/**
	 * @brief Registers all local parameters in a map for further use
	 * @param [in,out] parameters Map in which the parameters are stored
	 * @param [in] unitOpIdx Index of the unit operation used for registering the parameters
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 */
	inline void registerParameters(std::unordered_map<ParameterId, active*>& parameters, unsigned int unitOpIdx, unsigned int nComp, unsigned int const* nBoundStates)
	{
		registerComponentBoundStateDependentParam(hashString("SASKA_H"), parameters, h, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("SASKA_K"), parameters, k, unitOpIdx);
	}

	std::vector<active> h; //!< Henry coefficient
	util::SlicedVector<active> k; //!< Quadratic factor
};

/**
 * @brief Handles Saska binding model parameters that depend on an external function
 */
struct ExtSaskaParamHandler : public ExternalBindingParamHandlerBase
{
	static const char* identifier() { return "EXT_SASKA"; }

	/**
	 * @brief Reads parameters and verifies them
	 * @details See IBindingModel::configure() for details.
	 * @param [in] paramProvider IParameterProvider used for reading parameters
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 * @return @c true if the parameters were read and validated successfully, otherwise @c false
	 */
	inline bool configure(IParameterProvider& paramProvider, unsigned int nComp, unsigned int const* nBoundStates)
	{
		CADET_RESERVE_SPACE(h, nComp);
		CADET_RESERVE_SPACE2(k, nComp * nComp, nComp);

		CADET_READPAR_MATRIX(h, paramProvider, "SASKA_H", nComp, 1);
		CADET_READPAR_BOUNDSTATEDEP(util::SlicedVector<active>, active, k, paramProvider, "SASKA_K", nComp, nComp);

		return ExternalBindingParamHandlerBase::configure(paramProvider, 2);
	}

	/**
	 * @brief Registers all local parameters in a map for further use
	 * @param [in,out] parameters Map in which the parameters are stored
	 * @param [in] unitOpIdx Index of the unit operation used for registering the parameters
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 */
	inline void registerParameters(std::unordered_map<ParameterId, active*>& parameters, unsigned int unitOpIdx, unsigned int nComp, unsigned int const* nBoundStates)
	{
		CADET_REGPAR_COMPBND_VEC("SASKA_H", parameters, h, unitOpIdx);
		CADET_REGPAR_COMPBND("SASKA_K", parameters, k, unitOpIdx);
	}

	/**
	 * @brief Updates local parameter cache in order to take the external profile into account
	 * @details This function is declared const since the actual parameters are left unchanged by the method.
	 *         The cache is marked as mutable in order to make it writable.
	 * @param [in] t Current time
	 * @param [in] z Axial coordinate in the column
	 * @param [in] r Radial coordinate in the bead
	 * @param [in] secIdx Index of the current section
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 */
	inline void update(double t, double z, double r, unsigned int secIdx, unsigned int nComp, unsigned int const* nBoundStates) const
	{
		evaluateExternalFunctions(t, z, r, secIdx);
		for (unsigned int i = 0; i < nComp; ++i)
			CADET_UPDATE_EXTDEP_VARIABLE_BRACES(h, i, _extFunBuffer[0]);

		for (unsigned int i = 0; i < nComp * nComp; ++i)
			CADET_UPDATE_EXTDEP_VARIABLE_NATIVE(k, i, _extFunBuffer[1]);
	}

	/**
	 * @brief Updates local parameter cache and calculates time derivative in case of external dependence
	 * @details This function is declared const since the actual parameters are left unchanged by the method.
	 *         The cache is marked as mutable in order to make it writable.
	 * @param [in] t Current time
	 * @param [in] z Axial coordinate in the column
	 * @param [in] r Radial coordinate in the bead
	 * @param [in] secIdx Index of the current section
	 * @param [in] nComp Number of components
	 * @param [in] nBoundStates Array with number of bound states for each component
	 */
	inline void updateTimeDerivative(double t, double z, double r, unsigned int secIdx, unsigned int nComp, unsigned int const* nBoundStates) const
	{
		const std::vector<double> extTimeDeriv = evaluateTimeDerivativeExternalFunctions(t, z, r, secIdx);
		for (unsigned int i = 0; i < nComp; ++i)
			CADET_UPDATE_EXTDEP_VARIABLE_TDIFF_BRACES(h, i, _extFunBuffer[0], extTimeDeriv[0]);

		for (unsigned int i = 0; i < nComp * nComp; ++i)
			CADET_UPDATE_EXTDEP_VARIABLE_TDIFF_NATIVE(k, i, _extFunBuffer[1], extTimeDeriv[1]);
	}

	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, h)
	CADET_DEFINE_EXTDEP_VARIABLE(util::SlicedVector<active>, k)
};


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

	virtual int residual(const active& t, double z, double r, unsigned int secIdx, const active& timeFactor,
		active const* y, double const* yDot, active* res) const;

	virtual int residual(double t, double z, double r, unsigned int secIdx, double timeFactor,
		active const* y, double const* yDot, active* res) const;

	virtual int residual(const active& t, double z, double r, unsigned int secIdx, const active& timeFactor,
		double const* y, double const* yDot, active* res) const;

	virtual int residual(double t, double z, double r, unsigned int secIdx, double timeFactor,
		double const* y, double const* yDot, double* res) const;

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { _p.setExternalFunctions(extFuns, size); }
	virtual bool dependsOnTime() const CADET_NOEXCEPT { return ParamHandler_t::dependsOnTime(); }

	virtual void timeDerivativeAlgebraicResidual(double t, double z, double r, unsigned int secIdx, double const* y, double* dResDt) const
	{
		if (!hasAlgebraicEquations())
			return;

		if (!ParamHandler_t::dependsOnTime())
			return;

		// Update external function and compute time derivative of parameters
		_p.update(t, z, r, secIdx, _nComp, _nBoundStates);
		ParamHandler_t dpDt = _p;
		dpDt.updateTimeDerivative(t, z, r, secIdx, _nComp, _nBoundStates);

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

protected:
	ParamHandler_t _p; //!< Handles parameters and their dependence on external functions

	virtual bool configureImpl(bool reconfigure, IParameterProvider& paramProvider, unsigned int unitOpIdx)
	{
		// Read parameters
		_p.configure(paramProvider, _nComp, _nBoundStates);

		// Register parameters
		_p.registerParameters(_parameters, unitOpIdx, _nComp, _nBoundStates);

		return true;
	}

	virtual int residualCore(double t, double z, double r, unsigned int secIdx, double timeFactor,
		double const* y, double const* yCp, double const* yDot, double* res) const;
	virtual int residualCore(double t, double z, double r, unsigned int secIdx, double timeFactor,
		active const* y, double const* yCp, double const* yDot, active* res) const;

	virtual void analyticJacobianCore(double t, double z, double r, unsigned int secIdx, double const* y, 
		double const* yCp, linalg::BandMatrix::RowIterator jac) const;
	virtual void analyticJacobianCore(double t, double z, double r, unsigned int secIdx, double const* y, 
		double const* yCp, linalg::detail::DenseMatrixBase::RowIterator jac) const;

	template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
	int residualImpl(const ParamType& t, double z, double r, unsigned int secIdx, const ParamType& timeFactor,
		StateType const* y, CpStateType const* yCp, double const* yDot, ResidualType* res) const
	{
		_p.update(static_cast<double>(t), z, r, secIdx, _nComp, _nBoundStates);

		// Protein equations: dq_i / dt - ( H_i c_{p,i} + \sum_j k_{ij} c_{p,i} c_{p,j} - q_i ) == 0
		//               <=>  dq_i / dt == H_i c_{p,i} + \sum_j k_{ij} c_{p,i} c_{p,j} - q_i
		unsigned int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			const active h = _p.h[i];
			active const* const kSlice = _p.k[i];

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
	void jacobianImpl(double t, double z, double r, unsigned int secIdx, double const* y, double const* yCp, RowIterator jac) const
	{
		_p.update(t, z, r, secIdx, _nComp, _nBoundStates);

		// Protein equations: dq_i / dt - ( H_i c_{p,i} + \sum_j k_{ij} c_{p,i} c_{p,j} - q_i ) == 0
		//               <=>  dq_i / dt - ( H_i c_{p,i} + c_{p,i} * \sum_j k_{ij} c_{p,j} - q_i ) == 0
		int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			const double h = static_cast<double>(_p.h[i]);
			active const* const kSlice = _p.k[i];

			// dres_i / dc_{p,j}
			for (int j = 0; j < _nComp; ++j)
				jac[j - bndIdx - _nComp] = -static_cast<double>(kSlice[j]) * yCp[i];
				// Getting to c_{p,i}: -bndIdx takes us to q_0, another -nComp to c_{p,0} and a +j to c_{p,j}.
				//                     This means jac[j - bndIdx - nComp] corresponds to c_{p,j}.

			// dres_i / dc_{p,i}
			for (int j = 0; j < _nComp; ++j)
				jac[i - bndIdx - _nComp] -= static_cast<double>(kSlice[j]) * yCp[j];
			
			jac[i - bndIdx - _nComp] -= h;

			// dres_i / dq_i
			jac[0] = 1.0;

			// Advance to next equation and Jacobian row
			++bndIdx;
			++jac;
		}
	}
};

CADET_PUREBINDINGMODELBASE_TEMPLATED_BOILERPLATE_IMPL(SaskaBindingBase, ParamHandler_t)

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
