// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//
//  Copyright © 2008-2016: The CADET Authors
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

#include <functional>
#include <unordered_map>
#include <string>
#include <vector>

namespace cadet
{

	namespace model
	{

		/**
		* @brief Handles HICWANG binding model parameters that do not depend on external functions
		*/
		struct HICWANGParamHandler : public BindingParamHandlerBase
		{
			static const char* identifier() { return "HICWANG"; }

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
				beta0 = paramProvider.getDouble("BETA0");
				beta1 = paramProvider.getDouble("BETA1");

				readParameterMatrix(kKin, paramProvider, "KKIN", nComp, 1);
				readParameterMatrix(kEQ, paramProvider, "KEQ", nComp, 1);
				readParameterMatrix(nu, paramProvider, "NU", nComp, 1);
				readParameterMatrix(qMax, paramProvider, "QMAX", nComp, 1);

				// Check parameters
				if ((kKin.size() != kEQ.size()) || (kKin.size() != nu.size()) || (kKin.size() != qMax.size()))
					throw InvalidParameterException("NU, KKIN, KEQ, and QMAX must have the same size");

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
				registerComponentBoundStateDependentParam(hashString("KKIN"), parameters, kKin, unitOpIdx);
				registerComponentBoundStateDependentParam(hashString("KEQ"), parameters, kEQ, unitOpIdx);
				registerComponentBoundStateDependentParam(hashString("NU"), parameters, nu, unitOpIdx);
				registerComponentBoundStateDependentParam(hashString("QMAX"), parameters, qMax, unitOpIdx);

				parameters[makeParamId(hashString("BETA0"), unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &beta0;
				parameters[makeParamId(hashString("BETA1"), unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &beta1;
			}

			std::vector<active> nu; //!< NU
			std::vector<active> kEQ; //!< KEQ
			std::vector<active> qMax; //!< qMax
			std::vector<active> kKin; //!< kKIn

			active beta0;
			active beta1;
		};

		/**
		* @brief Handles HICWANG binding model parameters that depend on an external function
		*/
		struct ExtHICWANGParamHandler : public ExternalBindingParamHandlerBase
		{
			static const char* identifier() { return "EXT_HICWANG"; }

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
				CADET_READPAR_MATRIX(kKin, paramProvider, "KKIN", nComp, 1);
				CADET_READPAR_MATRIX(kEQ, paramProvider, "KEQ", nComp, 1);
				CADET_READPAR_MATRIX(nu, paramProvider, "NU", nComp, 1);
				CADET_READPAR_MATRIX(qMax, paramProvider, "QMAX", nComp, 1);

				CADET_READPAR_SCALAR(beta0, paramProvider, "BETA0");
				CADET_READPAR_SCALAR(beta1, paramProvider, "BETA1");

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
				CADET_REGPAR_COMPBND_VEC("KKIN", parameters, kKin, unitOpIdx);
				CADET_REGPAR_COMPBND_VEC("KEQ", parameters, kEQ, unitOpIdx);
				CADET_REGPAR_COMPBND_VEC("NU", parameters, nu, unitOpIdx);
				CADET_REGPAR_COMPBND_VEC("QMAX", parameters, qMax, unitOpIdx);

				parameters[makeParamId(hashString("EXT_BETA0"), unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &beta0;
				parameters[makeParamId(hashString("EXT_BETA1"), unitOpIdx, CompIndep, BoundPhaseIndep, ReactionIndep, SectionIndep)] = &beta1;
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
				{
					CADET_UPDATE_EXTDEP_VARIABLE_BRACES(kKin, i, _extFunBuffer[0]);
					CADET_UPDATE_EXTDEP_VARIABLE_BRACES(kEQ, i, _extFunBuffer[1]);
					CADET_UPDATE_EXTDEP_VARIABLE_BRACES(nu, i, _extFunBuffer[2]);
					CADET_UPDATE_EXTDEP_VARIABLE_BRACES(qMax, i, _extFunBuffer[3]);
				}
				CADET_UPDATE_EXTDEP_VARIABLE(beta0, _extFunBuffer[5]);
				CADET_UPDATE_EXTDEP_VARIABLE(beta1, _extFunBuffer[6]);
			}

			CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, kKin)
				CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, kEQ)
				CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, nu)
				CADET_DEFINE_EXTDEP_VARIABLE(std::vector<active>, qMax)

				CADET_DEFINE_EXTDEP_VARIABLE(active, beta0)
				CADET_DEFINE_EXTDEP_VARIABLE(active, beta1)
		};


		/**
		* @brief Defines the multi component HICWANG binding model
		* @details Implements the HICWANG adsorption model: \f[ \begin{align}
		*              \frac{\mathrm{d}q_i}{\mathrm{d}t} &= k_{a,i} c_{p,i} q_{\text{max},i} \left( 1 - \sum_j \frac{q_j}{q_{\text{max},j}} \right) - k_{d,i} q_i
		*          \end{align} \f]
		*          Multiple bound states are not supported.
		*          Components without bound state (i.e., non-binding components) are supported.
		*
		*          See @cite HICWANG1916.
		* @tparam ParamHandler_t Type that can add support for external function dependence
		*/
		template <class ParamHandler_t>
		class HICWANGBindingBase : public PureBindingModelBase
		{
		public:

			HICWANGBindingBase() { }
			virtual ~HICWANGBindingBase() CADET_NOEXCEPT { }

			static const char* identifier() { return ParamHandler_t::identifier(); }
			virtual const char* name() const CADET_NOEXCEPT { return ParamHandler_t::identifier(); }

			virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { _p.setExternalFunctions(extFuns, size); }

			virtual bool hasSalt() const CADET_NOEXCEPT { return true; }
			virtual bool supportsMultistate() const CADET_NOEXCEPT { return false; }
			virtual bool supportsNonBinding() const CADET_NOEXCEPT { return true; }
			virtual bool hasAlgebraicEquations() const CADET_NOEXCEPT { return !_kineticBinding; }
			virtual bool dependsOnTime() const CADET_NOEXCEPT { return ParamHandler_t::dependsOnTime(); }

		CADET_PUREBINDINGMODELBASE_BOILERPLATE

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

			template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
			int residualImpl(const ParamType& t, double z, double r, unsigned int secIdx, const ParamType& timeFactor,
				StateType const* y, CpStateType const* yCp, double const* yDot, ResidualType* res) const
			{
				_p.update(static_cast<double>(t), z, r, secIdx, _nComp, _nBoundStates);

				auto beta0 = static_cast<ParamType>(_p.beta0);
				auto beta1 = static_cast<ParamType>(_p.beta1);

				auto beta = static_cast<ParamType>(beta0 * exp(beta1 * yCp[0]));
				//LOG(Debug) << "beta\t" << beta;

				ResidualType qSum = 1.0;
				ResidualType qProd = 0.0;
				unsigned int bndIdx = 0;
				for (int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states (bound state index bndIdx is not advanced)
					if (_nBoundStates[i] == 0)
						continue;

					auto q = static_cast<ParamType>(y[bndIdx]);
					qSum -= y[bndIdx] / static_cast<ParamType>(_p.qMax[i]);
					//LOG(Debug) << "y[bndIdx]\t" << y[bndIdx];

					if (q < 0)
					{
						q = 0.0;
					}
					else
					{
						qProd += pow(q, static_cast<ParamType>(_p.nu[i])*beta);
					}
					// Next bound component
					++bndIdx;
				}

				//LOG(Debug) << "qSum\t" << qSum;
				//LOG(Debug) << "qProd\t" << qProd;

				bndIdx = 0;

				for (int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states (bound state index bndIdx is not advanced)
					if (_nBoundStates[i] == 0)
						continue;
					auto qMax = static_cast<ParamType>(_p.qMax[i]);
					auto kKin = static_cast<ParamType>(_p.kKin[i]);
					auto kEQ = static_cast<ParamType>(_p.kEQ[i]);
					auto nu = static_cast<ParamType>(_p.nu[i]);

					if (y[bndIdx] < 0 || yCp[i] < 0)
					{
						res[bndIdx] = kKin * kEQ * yCp[i] * (nu*y[bndIdx]/qMax-1);
					}
					else
					{
						res[bndIdx] = kKin * (y[bndIdx] * static_cast<ParamType>(qProd) - kEQ * static_cast<ParamType>(pow(qSum, nu)) * yCp[i]);
					}

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

			}
		};

		typedef HICWANGBindingBase<HICWANGParamHandler> HICWANGBinding;
		typedef HICWANGBindingBase<ExtHICWANGParamHandler> ExternalHICWANGBinding;

		namespace binding
		{
			void registerHICWANGModel(std::unordered_map<std::string, std::function<model::IBindingModel*()>>& bindings)
			{
				bindings[HICWANGBinding::identifier()] = []() { return new HICWANGBinding(); };
				bindings[ExternalHICWANGBinding::identifier()] = []() { return new ExternalHICWANGBinding(); };
			}
		}  // namespace binding

	}  // namespace model

}  // namespace cadet
