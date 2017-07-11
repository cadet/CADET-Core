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

#include <vector>
#include <unordered_map>
#include <functional>

namespace cadet
{

namespace model
{

/**
 * @brief Handles multi component spreading binding model parameters that do not depend on external functions
 */
struct SpreadingParamHandler : public BindingParamHandlerBase
{
	static const char* identifier() { return "MULTI_COMPONENT_SPREADING"; }

	/**
	 * @brief Reads parameters and verifies them
	 * @details See IBindingModel::configure() for details.
	 * @param [in] paramProvider IParameterProvider used for reading parameters
	 * @param [in] nComp Number of components
	 * @param [in] numStates Number of binding site types
	 * @return @c true if the parameters were read and validated successfully, otherwise @c false
	 */
	inline bool configure(IParameterProvider& paramProvider, unsigned int nComp, unsigned int numStates)
	{
		readBoundStateDependentParameter<util::SlicedVector<active>, active>(kA, paramProvider, "MCSPR_KA", nComp, 2);
		readBoundStateDependentParameter<util::SlicedVector<active>, active>(kD, paramProvider, "MCSPR_KD", nComp, 2);
		readBoundStateDependentParameter<util::SlicedVector<active>, active>(qMax, paramProvider, "MCSPR_QMAX", nComp, 2);
		readParameterMatrix(k12, paramProvider, "MCSPR_K12", nComp, 1);
		readParameterMatrix(k21, paramProvider, "MCSPR_K21", nComp, 1);

		// Check parameters
		if ((kA.size() != kD.size()) || (kA.size() != qMax.size()) || (kA.size() < nComp * 2))
			throw InvalidParameterException("MCSPR_KA, MCSPR_KD, and MCSPR_QMAX have to have the same size");
		if ((k12.size() != k21.size()) || (k12.size() < nComp))
			throw InvalidParameterException("MCSPR_K12 and MCSPR_K21 have to have the same size (number of components)");

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
		registerComponentBoundStateDependentParam(hashString("MCSPR_KA"), parameters, kA, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("MCSPR_KD"), parameters, kD, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("MCSPR_QMAX"), parameters, qMax, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("MCSPR_K12"), parameters, k12, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("MCSPR_K21"), parameters, k21, unitOpIdx);
	}

	/**
	 * @brief Reserves space in the storage of the parameters
	 * @param [in] nComp Number of components
	 * @param [in] numSlices Number of slices / binding site types
	 */
	inline void reserve(unsigned int nComp, unsigned int numSlices)
	{
		kA.reserve(nComp * numSlices, numSlices);
		kD.reserve(nComp * numSlices, numSlices);
		qMax.reserve(nComp * numSlices, numSlices);
		k21.reserve(nComp);
		k12.reserve(nComp);
	}

	util::SlicedVector<active> kA; //!< Adsorption rate in state-major ordering
	util::SlicedVector<active> kD; //!< Desorption rate in state-major ordering
	util::SlicedVector<active> qMax; //!< Capacity in state-major ordering
	std::vector<active> k12; //<! Transition rate from state @f$ q_i^A @f$ to state @f$ q_i^B @f$
	std::vector<active> k21; //<! Transition rate from state @f$ q_i^B @f$ to state @f$ q_i^A @f$
};

/**
 * @brief Handles multi component spreading binding model parameters that depend on an external function
 */
struct ExtSpreadingParamHandler : public ExternalBindingParamHandlerBase
{
	static const char* identifier() { return "EXT_MULTI_COMPONENT_SPREADING"; }

	/**
	 * @brief Reads parameters and verifies them
	 * @details See IBindingModel::configure() for details.
	 * @param [in] paramProvider IParameterProvider used for reading parameters
	 * @param [in] nComp Number of components
	 * @param [in] numSlices Number of binding site types
	 * @return @c true if the parameters were read and validated successfully, otherwise @c false
	 */
	inline bool configure(IParameterProvider& paramProvider, unsigned int nComp, unsigned int numSlices)
	{
		CADET_READPAR_BOUNDSTATEDEP(util::SlicedVector<active>, active, kA, paramProvider, "MCSPR_KA", nComp, 2);
		CADET_READPAR_BOUNDSTATEDEP(util::SlicedVector<active>, active, kD, paramProvider, "MCSPR_KD", nComp, 2);
		CADET_READPAR_BOUNDSTATEDEP(util::SlicedVector<active>, active, qMax, paramProvider, "MCSPR_QMAX", nComp, 2);
		CADET_READPAR_MATRIX(k12, paramProvider, "MCSPR_K12", nComp, 1);
		CADET_READPAR_MATRIX(k21, paramProvider, "MCSPR_K21", nComp, 1);

		return ExternalBindingParamHandlerBase::configure(paramProvider, 5);
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
		CADET_REGPAR_COMPBND("MCSPR_KA", parameters, kA, unitOpIdx);
		CADET_REGPAR_COMPBND("MCSPR_KD", parameters, kD, unitOpIdx);
		CADET_REGPAR_COMPBND("MCSPR_QMAX", parameters, qMax, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("MCSPR_K12", parameters, k12, unitOpIdx);
		CADET_REGPAR_COMPBND_VEC("MCSPR_K21", parameters, k21, unitOpIdx);
	}

	/**
	 * @brief Reserves space in the storage of the parameters
	 * @param [in] nComp Number of components
	 * @param [in] numSlices Number of slices / binding site types
	 */
	inline void reserve(unsigned int nComp, unsigned int numSlices)
	{
		CADET_RESERVE_SPACE2(kA, nComp * numSlices, numSlices);
		CADET_RESERVE_SPACE2(kD, nComp * numSlices, numSlices);
		CADET_RESERVE_SPACE2(qMax, nComp * numSlices, numSlices);
		CADET_RESERVE_SPACE(k12, nComp);
		CADET_RESERVE_SPACE(k21, nComp);
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
		for (unsigned int i = 0; i < kAT0.size(); ++i)
		{
			CADET_UPDATE_EXTDEP_VARIABLE_NATIVE(kA, i, _extFunBuffer[0]);
			CADET_UPDATE_EXTDEP_VARIABLE_NATIVE(kD, i, _extFunBuffer[1]);
			CADET_UPDATE_EXTDEP_VARIABLE_NATIVE(qMax, i, _extFunBuffer[2]);
		}

		for (unsigned int i = 0; i < k12T0.size(); ++i)
		{
			CADET_UPDATE_EXTDEP_VARIABLE_BRACES(k12, i, _extFunBuffer[3]);
			CADET_UPDATE_EXTDEP_VARIABLE_BRACES(k21, i, _extFunBuffer[4]);
		}
	}

	CADET_DEFINE_EXTDEP_VARIABLE(util::SlicedVector<active>, kA)
	CADET_DEFINE_EXTDEP_VARIABLE(util::SlicedVector<active>, kD)
	CADET_DEFINE_EXTDEP_VARIABLE(util::SlicedVector<active>, qMax)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, k21)
	CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, k12)
};


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
 *          Internal state vector layout is state-major. First, all components of state A are placed, then all components of state B.
 * @tparam ParamHandler_t Type that can add support for external function dependence
 */
template <class ParamHandler_t>
class MultiComponentSpreadingBindingBase : public PureBindingModelBase
{
public:

	MultiComponentSpreadingBindingBase() { }
	virtual ~MultiComponentSpreadingBindingBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }
	virtual const char* name() const CADET_NOEXCEPT { return ParamHandler_t::identifier(); }

	virtual void configureModelDiscretization(unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
	{
		BindingModelBase::configureModelDiscretization(nComp, nBound, boundOffset);

		unsigned int numSlices = 2;
		for (unsigned int i = 0; i < nComp; ++i)
		{
			if (nBound[i] == 0)
				continue;

			if (nBound[i] != numSlices)
				throw InvalidParameterException("Multi component spreading binding model requires exactly two bound states for all (binding) components");
		}

		_numBindingComp = numBindingComponents(_nBoundStates, _nComp);

		// Allocate space for parameters
		_p.reserve(nComp, numSlices);
	}

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

	virtual bool hasSalt() const CADET_NOEXCEPT { return false; }
	virtual bool supportsMultistate() const CADET_NOEXCEPT { return true; }
	virtual bool supportsNonBinding() const CADET_NOEXCEPT { return true; }

protected:
	ParamHandler_t _p; //!< Handles parameters and their dependence on external functions
	unsigned int _numBindingComp; //!< Number of binding components

	virtual bool configureImpl(bool reconfigure, IParameterProvider& paramProvider, unsigned int unitOpIdx)
	{
		const unsigned int numStates = firstNonEmptyBoundStates(_nBoundStates, _nComp);

		// Read parameters
		_p.configure(paramProvider, _nComp, numStates);

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

		ResidualType qSum = 1.0;
		unsigned int bndIdx = 0;
		for (int j = 0; j < 2; ++j)
		{
			active const* const localQmax = _p.qMax[j];
			for (int i = 0; i < _nComp; ++i)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[i] == 0)
					continue;

				qSum -= y[bndIdx] / static_cast<ParamType>(localQmax[i]);

				// Next bound component
				++bndIdx;
			}
		}

		// First bound state (A)
		active const* localKa = _p.kA[0];
		active const* localKd = _p.kD[0];
		active const* localQmax = _p.qMax[0];

		bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			// Residual
			res[bndIdx] = static_cast<ParamType>(localKd[i]) * y[bndIdx] - static_cast<ParamType>(_p.k21[i]) * y[bndIdx + _numBindingComp] 
					- (static_cast<ParamType>(localKa[i]) * yCp[i] - static_cast<ParamType>(_p.k12[i]) * y[bndIdx]) * static_cast<ParamType>(localQmax[i]) * qSum;

			// Add time derivative if necessary
			if (_kineticBinding && yDot)
			{
				res[bndIdx] += timeFactor * yDot[bndIdx];
			}

			// Next bound component
			++bndIdx;
		}

		// Second bound state (B)
		localKa = _p.kA[1];
		localKd = _p.kD[1];
		// Still use q_{max}^A here

		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			// Residual
			res[bndIdx] = (static_cast<ParamType>(localKd[i]) + static_cast<ParamType>(_p.k21[i])) * y[bndIdx] 
					- (static_cast<ParamType>(localKa[i]) * yCp[i] + static_cast<ParamType>(_p.k12[i]) * y[bndIdx - _numBindingComp]) * static_cast<ParamType>(localQmax[i]) * qSum;

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

		double qSum = 1.0;
		int bndIdx = 0;
		for (int j = 0; j < 2; ++j)
		{
			active const* const localQmax = _p.qMax[j];
			for (int i = 0; i < _nComp; ++i)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[i] == 0)
					continue;

				qSum -= y[bndIdx] / static_cast<double>(localQmax[i]);

				// Next bound component
				++bndIdx;
			}
		}

		active const* localKa = _p.kA[0];
		active const* localKd = _p.kD[0];
		active const* const qMax1 = _p.qMax[0];
		active const* const qMax2 = _p.qMax[1];

		// First state (A)
		bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			const double ka = static_cast<double>(localKa[i]);
			const double kd = static_cast<double>(localKd[i]);

			// dres_i / dc_{p,i}
			jac[i - bndIdx - _nComp] = -ka * static_cast<double>(qMax1[i]) * qSum;
			// Getting to c_{p,i}: -bndIdx takes us to q_0, another -nComp to c_{p,0} and a +i to c_{p,i}.
			//                     This means jac[i - bndIdx - nComp] corresponds to c_{p,i}.

			// Fill dres_i / dq_j^A
			const double factor = (ka * yCp[i] - y[bndIdx] * static_cast<double>(_p.k12[i])) * static_cast<double>(qMax1[i]);
			int bndIdx2 = 0;
			for (int j = 0; j < _nComp; ++j)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[j] == 0)
					continue;

				// dres_i / dq_j^A
				jac[bndIdx2 - bndIdx] = factor / static_cast<double>(qMax1[j]);
				// Getting to q_j: -bndIdx takes us to q_0, another +bndIdx2 to q_j. This means jac[bndIdx2 - bndIdx] corresponds to q_j.

				++bndIdx2;
			}

			// Add to dres_i / dq_i^A
			jac[0] += kd + static_cast<double>(_p.k12[i]) * static_cast<double>(qMax1[i]) * qSum; // last summand by product rule

			// Fill dres_i / dq_j^B
			for (int j = 0; j < _nComp; ++j)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[j] == 0)
					continue;

				// dres_i / dq_j^B
				jac[bndIdx2 - bndIdx] = factor / static_cast<double>(qMax2[j]);
				// Getting to q_j: -bndIdx takes us to q_0, another +bndIdx2 to q_j. This means jac[bndIdx2 - bndIdx] corresponds to q_j.

				++bndIdx2;
			}

			// Add to dres_i / dq_i^B
			jac[_numBindingComp] -= static_cast<double>(_p.k21[i]);

			// Advance to next equation and Jacobian row
			++bndIdx;
			++jac;
		}

		// Second state (B)
		localKa = _p.kA[1];
		localKd = _p.kD[1];

		for (int i = 0; i < _nComp; ++i)
		{
			// Skip components without bound states (bound state index bndIdx is not advanced)
			if (_nBoundStates[i] == 0)
				continue;

			const double ka = static_cast<double>(localKa[i]);
			const double kd = static_cast<double>(localKd[i]);

			// dres_i / dc_{p,i}
			jac[i - bndIdx - _nComp] = -ka * static_cast<double>(qMax1[i]) * qSum;
			// Getting to c_{p,i}: -bndIdx takes us to q_0, another -nComp to c_{p,0} and a +i to c_{p,i}.
			//                     This means jac[i - bndIdx - nComp] corresponds to c_{p,i}.

			// Fill dres_i / dq_j^A
			const double factor = (ka * yCp[i] + y[bndIdx - static_cast<int>(_numBindingComp)] * static_cast<double>(_p.k12[i])) * static_cast<double>(qMax1[i]);
			int bndIdx2 = 0;
			for (int j = 0; j < _nComp; ++j)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[j] == 0)
					continue;

				// dres_i / dq_j^A
				jac[bndIdx2 - bndIdx] = factor / static_cast<double>(qMax1[j]);
				// Getting to q_j: -bndIdx takes us to q_0, another +bndIdx2 to q_j. This means jac[bndIdx2 - bndIdx] corresponds to q_j.

				++bndIdx2;
			}

			// Add to dres_i / dq_i^B
			jac[0] += kd + static_cast<double>(_p.k21[i]);

			// Fill dres_i / dq_j^B
			for (int j = 0; j < _nComp; ++j)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[j] == 0)
					continue;

				// dres_i / dq_j^B
				jac[bndIdx2 - bndIdx] = factor / static_cast<double>(qMax2[j]);
				// Getting to q_j: -bndIdx takes us to q_0, another +bndIdx2 to q_j. This means jac[bndIdx2 - bndIdx] corresponds to q_j.

				++bndIdx2;
			}

			// Add to dres_i / dq_i^A
			jac[-static_cast<int>(_numBindingComp)] -= static_cast<double>(_p.k12[i]) * static_cast<double>(qMax1[i]) * qSum;  // last summand by product rule

			// Advance to next equation and Jacobian row
			++bndIdx;
			++jac;
		}
	}
};

CADET_PUREBINDINGMODELBASE_TEMPLATED_BOILERPLATE_IMPL(MultiComponentSpreadingBindingBase, ParamHandler_t)


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
