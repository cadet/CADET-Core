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
 * @brief Handles Bi-Langmuir binding model parameters that do not depend on external functions
 */
struct BiLangmuirParamHandler : public BindingParamHandlerBase
{
	static const char* identifier() { return "MULTI_COMPONENT_BILANGMUIR"; }

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
		readBoundStateDependentParameter<util::SlicedVector<active>, active>(kA, paramProvider, "MCBL_KA", nComp, numSlices);
		readBoundStateDependentParameter<util::SlicedVector<active>, active>(kD, paramProvider, "MCBL_KD", nComp, numSlices);
		readBoundStateDependentParameter<util::SlicedVector<active>, active>(qMax, paramProvider, "MCBL_QMAX", nComp, numSlices);

		// Check parameters
		if ((kA.size() != kD.size()) || (kA.size() != qMax.size()) || (kA.size() < nComp * numSlices))
			throw InvalidParameterException("MCBL_KA, MCBL_KD, and MCBL_QMAX have to have the same size");

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
		registerComponentBoundStateDependentParam(hashString("MCBL_KA"), parameters, kA, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("MCBL_KD"), parameters, kD, unitOpIdx);
		registerComponentBoundStateDependentParam(hashString("MCBL_QMAX"), parameters, qMax, unitOpIdx);
	}

	/**
	 * @brief Reserves space in the storage of the parameters
	 * @param [in] numElem Number of elements (total)
	 * @param [in] numSlices Number of slices / binding site types
	 */
	inline void reserve(unsigned int numElem, unsigned int numSlices)
	{
		kA.reserve(numElem, numSlices);
		kD.reserve(numElem, numSlices);
		qMax.reserve(numElem, numSlices);
	}

	util::SlicedVector<active> kA; //!< Adsorption rate in state-major ordering
	util::SlicedVector<active> kD; //!< Desorption rate in state-major ordering
	util::SlicedVector<active> qMax; //!< Capacity in state-major ordering
};

/**
 * @brief Handles Bi-Langmuir binding model parameters that depend on an external function
 */
struct ExtBiLangmuirParamHandler : public ExternalBindingParamHandlerBase
{
	static const char* identifier() { return "EXT_MULTI_COMPONENT_BILANGMUIR"; }

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
		CADET_READPAR_BOUNDSTATEDEP(util::SlicedVector<active>, active, kA, paramProvider, "MCBL_KA", nComp, numSlices);
		CADET_READPAR_BOUNDSTATEDEP(util::SlicedVector<active>, active, kD, paramProvider, "MCBL_KD", nComp, numSlices);
		CADET_READPAR_BOUNDSTATEDEP(util::SlicedVector<active>, active, qMax, paramProvider, "MCBL_QMAX", nComp, numSlices);

		return ExternalBindingParamHandlerBase::configure(paramProvider, 3);
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
		CADET_REGPAR_COMPBND("MCBL_KA", parameters, kA, unitOpIdx);
		CADET_REGPAR_COMPBND("MCBL_KD", parameters, kD, unitOpIdx);
		CADET_REGPAR_COMPBND("MCBL_QMAX", parameters, qMax, unitOpIdx);
	}

	/**
	 * @brief Reserves space in the storage of the parameters
	 * @param [in] numElem Number of elements (total)
	 * @param [in] numSlices Number of slices / binding site types
	 */
	inline void reserve(unsigned int numElem, unsigned int numSlices)
	{
		CADET_RESERVE_SPACE2(kA, numElem, numSlices);
		CADET_RESERVE_SPACE2(kD, numElem, numSlices);
		CADET_RESERVE_SPACE2(qMax, numElem, numSlices);
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
		{
			CADET_UPDATE_EXTDEP_VARIABLE_TDIFF_NATIVE(kA, i, _extFunBuffer[0], extTimeDeriv[0]);
			CADET_UPDATE_EXTDEP_VARIABLE_TDIFF_NATIVE(kD, i, _extFunBuffer[1], extTimeDeriv[1]);
			CADET_UPDATE_EXTDEP_VARIABLE_TDIFF_NATIVE(qMax, i, _extFunBuffer[2], extTimeDeriv[2]);
		}
	}

	CADET_DEFINE_EXTDEP_VARIABLE(util::SlicedVector<active>, kA)
	CADET_DEFINE_EXTDEP_VARIABLE(util::SlicedVector<active>, kD)
	CADET_DEFINE_EXTDEP_VARIABLE(util::SlicedVector<active>, qMax)
};


/**
 * @brief Defines the multi component Bi-Langmuir binding model
 * @details Implements the Bi-Langmuir adsorption model: \f[ \begin{align} 
 *              \frac{\mathrm{d}q_i^A}{\mathrm{d}t} &= k_{a,i}^A c_{p,i} q_{\text{max},i}^A \left( 1 - \sum_j \frac{q_j^A}{q_{\text{max},j}^A} \right) - k_{d,i}^A q_i^A \\
 *              \frac{\mathrm{d}q_i^B}{\mathrm{d}t} &= k_{a,i}^B c_{p,i} q_{\text{max},i}^B \left( 1 - \sum_j \frac{q_j^B}{q_{\text{max},j}^B} \right) - k_{d,i}^B q_i^B \\
 *              \vodts & \vdots
 *          \end{align} \f]
 *          Here, several different types of binding sites @f$ q^A @f$, @f$ q^B @f$, etc. are considered. A molecule can either bind to
 *          site A or B (or C, etc.). A direct exchange between the different binding sites does not occur.
 *          While components without bound state (i.e., non-binding components) are supported, all other components must have
 *          the same number of bound states (i.e., binding sites).
 *          
 *          Internal state vector order is state-major. The state vector is composed of all bound states and within each bound state
 *          all components are listed.
 *          
 *          See @cite Guiochon2006.
 * @tparam ParamHandler_t Type that can add support for external function dependence
 */
template <class ParamHandler_t>
class BiLangmuirBindingBase : public PureBindingModelBase
{
public:

	BiLangmuirBindingBase() { }
	virtual ~BiLangmuirBindingBase() CADET_NOEXCEPT { }

	static const char* identifier() { return ParamHandler_t::identifier(); }
	virtual const char* name() const CADET_NOEXCEPT { return ParamHandler_t::identifier(); }

	virtual void configureModelDiscretization(unsigned int nComp, unsigned int const* nBound, unsigned int const* boundOffset)
	{
		BindingModelBase::configureModelDiscretization(nComp, nBound, boundOffset);

		unsigned int numSlices = 0;
		for (unsigned int i = 0; i < nComp; ++i)
		{
			if (nBound[i] == 0)
				continue;

			if (numSlices == 0)
				numSlices = nBound[i];

			if (nBound[i] != numSlices)
				throw InvalidParameterException("Bi-Langmuir binding model requires exactly the same bound states for all components");
		}

		_numBindingComp = numBindingComponents(_nBoundStates, _nComp);

		// Allocate space for parameters
		_p.reserve(nComp * numSlices, numSlices);
	}

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { _p.setExternalFunctions(extFuns, size); }

	virtual bool hasSalt() const CADET_NOEXCEPT { return false; }
	virtual bool supportsMultistate() const CADET_NOEXCEPT { return true; }
	virtual bool supportsNonBinding() const CADET_NOEXCEPT { return true; }
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

		// Protein equations: dq_i^j / dt - ( k_{a,i}^j * c_{p,i} * (1 - \sum q_i^j / q_{max,i}^j) - k_{d,i}^j * q_i^j) == 0
		//               <=>  dq_i^j / dt == k_{a,i}^j * c_{p,i} * (1 - \sum q_i^j / q_{max,i}^j) - k_{d,i}^j * q_i^j

		double const* const yCp = y - _nComp;
		const unsigned int nSites = _p.kA.slices();

		// Ordering of the states is (q_{comp,state})
		// q_{0,0}, q{0,1}, q_{0,2}, q_{1,0}, q_{1,1}, q_{1,2}, ...
		// A state corresponds to a type of binding site. It is assumed that all components have either 0
		// or the same number of states. Thus, a component is either non-binding or has nSites bound states.
		//
		// The same ordering is used for the equations. That is, q_{0,0}, q_{1,0} and q_{0,1}, q_{1,1} and ... each
		// form one Langmuir binding model system.

		// Loop over all binding site types
		for (unsigned int site = 0; site < nSites; ++site, ++y, ++dResDt)
		{
			// y, and dResDt point to q_{0,site}

			// Get parameter slice for current binding site type
			active const* const localKa = _p.kA[site];
			active const* const localKd = _p.kD[site];
			active const* const localQmax = _p.qMax[site];
			active const* const localKaT = dpDt.kA[site];
			active const* const localKdT = dpDt.kD[site];
			active const* const localQmaxT = dpDt.qMax[site];

			double qSum = 1.0;
			double qSumT = 0.0;
			unsigned int bndIdx = 0;

			// bndIdx is used as a counter inside one binding site type
			// Getting from one component to another requires a step size of nSites (stride)

			for (int i = 0; i < _nComp; ++i)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[i] == 0)
					continue;

				const double summand = y[bndIdx * nSites] / static_cast<double>(localQmax[i]);
				qSum -= summand;
				qSumT += summand / static_cast<double>(localQmax[i]) * static_cast<double>(localQmaxT[i]);

				// Next bound component
				++bndIdx;
			}

			bndIdx = 0;
			for (int i = 0; i < _nComp; ++i)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[i] == 0)
					continue;

				// Residual
				dResDt[bndIdx * nSites] = static_cast<double>(localKdT[i]) * y[bndIdx * nSites] 
					- yCp[i] * (static_cast<double>(localKaT[i]) * static_cast<double>(localQmax[i]) * qSum 
					           + static_cast<double>(localKa[i]) * static_cast<double>(localQmaxT[i]) * qSum
					           + static_cast<double>(localKa[i]) * static_cast<double>(localQmax[i]) * qSumT);

				// Next bound component
				++bndIdx;
			}
		}
	}

	CADET_PUREBINDINGMODELBASE_BOILERPLATE

protected:
	ParamHandler_t _p; //!< Handles parameters and their dependence on external functions
	unsigned int _numBindingComp; //!< Number of binding components

	virtual bool configureImpl(bool reconfigure, IParameterProvider& paramProvider, unsigned int unitOpIdx)
	{
		const unsigned int numSlices = firstNonEmptyBoundStates(_nBoundStates, _nComp);

		// Read parameters
		_p.configure(paramProvider, _nComp, numSlices);

		// Register parameters
		_p.registerParameters(_parameters, unitOpIdx, _nComp, _nBoundStates);

		return true;
	}

	template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
	int residualImpl(const ParamType& t, double z, double r, unsigned int secIdx, const ParamType& timeFactor,
		StateType const* y, CpStateType const* yCp, double const* yDot, ResidualType* res) const
	{
		_p.update(static_cast<double>(t), z, r, secIdx, _nComp, _nBoundStates);

		const bool hasYdot = yDot;

		// Protein equations: dq_i^j / dt - ( k_{a,i}^j * c_{p,i} * (1 - \sum q_i^j / q_{max,i}^j) - k_{d,i}^j * q_i^j) == 0
		//               <=>  dq_i^j / dt == k_{a,i}^j * c_{p,i} * (1 - \sum q_i^j / q_{max,i}^j) - k_{d,i}^j * q_i^j

		const unsigned int nSites = _p.kA.slices();

		// Ordering of the states is (q_{comp,state})
		// q_{0,0}, q{0,1}, q_{0,2}, q_{1,0}, q_{1,1}, q_{1,2}, ...
		// A state corresponds to a type of binding site. It is assumed that all components have either 0
		// or the same number of states. Thus, a component is either non-binding or has nSites bound states.
		//
		// The same ordering is used for the equations. That is, q_{0,0}, q_{1,0} and q_{0,1}, q_{1,1} and ... each
		// form one Langmuir binding model system.

		// Loop over all binding site types
		for (unsigned int site = 0; site < nSites; ++site, ++y, ++yDot, ++res)
		{
			// y, yDot, and res point to q_{0,site}

			// Get parameter slice for current binding site type
			active const* const localKa = _p.kA[site];
			active const* const localKd = _p.kD[site];
			active const* const localQmax = _p.qMax[site];

			ResidualType qSum = 1.0;
			unsigned int bndIdx = 0;

			// bndIdx is used as a counter inside one binding site type
			// Getting from one component to another requires a step size of nSites (stride)

			for (int i = 0; i < _nComp; ++i)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[i] == 0)
					continue;

				qSum -= y[bndIdx * nSites] / static_cast<ParamType>(localQmax[i]);

				// Next bound component
				++bndIdx;
			}

			bndIdx = 0;
			for (int i = 0; i < _nComp; ++i)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[i] == 0)
					continue;

				// Residual
				res[bndIdx * nSites] = static_cast<ParamType>(localKd[i]) * y[bndIdx * nSites] - static_cast<ParamType>(localKa[i]) * yCp[i] * static_cast<ParamType>(localQmax[i]) * qSum;

				// Add time derivative if necessary
				if (_kineticBinding && hasYdot)
				{
					res[bndIdx * nSites] += timeFactor * yDot[bndIdx * nSites];
				}

				// Next bound component
				++bndIdx;
			}
		}

		return 0;
	}

	template <typename RowIterator>
	void jacobianImpl(double t, double z, double r, unsigned int secIdx, double const* y, double const* yCp, RowIterator jac) const
	{
		_p.update(t, z, r, secIdx, _nComp, _nBoundStates);

		// Protein equations: dq_i^j / dt - ( k_{a,i}^j * c_{p,i} * (1 - \sum q_i^j / q_{max,i}^j) - k_{d,i}^j * q_i^j) == 0

		// Ordering of the states is (q_{comp,state}, example uses 2 components, 3 binding sites)
		// q_{0,0}, q{0,1}, q_{0,2}, q_{1,0}, q_{1,1}, q_{1,2}, ...
		// Ordering of the equations is the same, that is, we need nSites steps to jump from one equation of a Langmuir
		// binding model system to the next.

		const int nSites = static_cast<int>(_p.kA.slices());

		// Loop over all binding site types
		for (unsigned int site = 0; site < nSites; ++site, ++y)
		{
			// Get parameter slice for current binding site type
			active const* const localKa = _p.kA[site];
			active const* const localKd = _p.kD[site];
			active const* const localQmax = _p.qMax[site];

			double qSum = 1.0;
			int bndIdx = 0;
			for (int i = 0; i < _nComp; ++i)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[i] == 0)
					continue;

				qSum -= y[bndIdx * nSites] / static_cast<double>(localQmax[i]);

				// Next bound component
				++bndIdx;
			}

			bndIdx = 0;
			for (int i = 0; i < _nComp; ++i)
			{
				// Skip components without bound states (bound state index bndIdx is not advanced)
				if (_nBoundStates[i] == 0)
					continue;

				const double ka = static_cast<double>(localKa[i]);
				const double kd = static_cast<double>(localKd[i]);

				// dres_i / dc_{p,i}
				jac[i - site - _nComp - nSites * bndIdx] = -ka * static_cast<double>(localQmax[i]) * qSum;
				// Getting to c_{p,i}: -nSites * bndIdx takes us to q_{0,site}, another -site to q_{0,0}. From there, we
				//                     take a -_nComp to reach c_{p,0} and a +i to arrive at c_{p,i}.
				//                     This means jac[i - site - _nComp - nSites * bndIdx] corresponds to c_{p,i}.

				// Fill dres_i / dq_j
				int bndIdx2 = 0;
				for (int j = 0; j < _nComp; ++j)
				{
					// Skip components without bound states (bound state index bndIdx2 is not advanced)
					if (_nBoundStates[j] == 0)
						continue;

					// dres_i / dq_j
					jac[(bndIdx2 - bndIdx) * nSites] = ka * yCp[i] * static_cast<double>(localQmax[i]) / static_cast<double>(localQmax[j]);
					// Getting to q_j: -bndIdx * nSites takes us to q_{0,site}, another +bndIdx2 to q_{j,site}.
					// This means jac[(bndIdx2 - bndIdx) * nSites] corresponds to q_{j,site}.

					++bndIdx2;
				}

				// Add to dres_i / dq_{i,site}
				jac[0] += kd;

				// Advance to next equation and Jacobian row
				++bndIdx;
				jac += nSites;
				// Note that there is a spacing of nSites between the equations inside one binding site type
			}
			// We are at the end of the equations block q_{_numBindingComp,site} and need to jump back to the beginning
			// by using -_numBindingComp * nSites steps. From there, we take one step forward to arrive at q_{0,site+1}.
			jac -= _numBindingComp * nSites - 1;
		}
	}
};

typedef BiLangmuirBindingBase<BiLangmuirParamHandler> BiLangmuirBinding;
typedef BiLangmuirBindingBase<ExtBiLangmuirParamHandler> ExternalBiLangmuirBinding;

namespace binding
{
	void registerBiLangmuirModel(std::unordered_map<std::string, std::function<model::IBindingModel*()>>& bindings)
	{
		bindings[BiLangmuirBinding::identifier()] = []() { return new BiLangmuirBinding(); };
		bindings[ExternalBiLangmuirBinding::identifier()] = []() { return new ExternalBiLangmuirBinding(); };
	}
}  // namespace binding

}  // namespace model

}  // namespace cadet
