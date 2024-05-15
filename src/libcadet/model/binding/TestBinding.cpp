// =============================================================================
//  CADET
//  
//  Copyright © 2008-2021: The CADET Authors
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
#include "LocalVector.hpp"
#include "SimulationTypes.hpp"

#include <functional>
#include <unordered_map>
#include <string>
#include <vector>

/* In the following lines of code the interface to the binding model is defined.
 In the interface the binding model specific parameters are defined. Please modify configuration name and 
 parameter names according to the naming convection explained in the tutorial.*/

/*<codegen>
{
	"name": "TestParamHandler",
	"externalName": "ExtTestParamHandler",
	"parameters":
		[
			{ "type": "ScalarComponentDependentParameter", "varName": "kA", "confName": "MCL_KA"},
			{ "type": "ScalarComponentDependentParameter", "varName": "kD", "confName": "MCL_KD"},
			{ "type": "ScalarComponentDependentParameter", "varName": "qMax", "confName": "MCL_QMAX"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p1", "confName": "ANN_P1"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p2", "confName": "ANN_P2"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p3", "confName": "ANN_P3"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p4", "confName": "ANN_P4"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p5", "confName": "ANN_P5"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p6", "confName": "ANN_P6"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p7", "confName": "ANN_P7"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p8", "confName": "ANN_P8"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p9", "confName": "ANN_P9"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p10", "confName": "ANN_P10"},
			{ "type": "ScalarComponentDependentParameter", "varName": "p11", "confName": "ANN_P11"},
			{ "type": "ScalarComponentDependentParameter", "varName": "P_norm", "confName": "ANN_P_NORM"},
			{ "type": "ScalarComponentDependentParameter", "varName": "q_norm", "confName": "ANN_Q_NORM"},
			{ "type": "ScalarComponentDependentParameter", "varName": "S_norm", "confName": "ANN_S_NORM"}
		]
}
</codegen>*/

/* Parameter description
 ------------------------
 kA = Adsorption rate
 kD = Desorption rate
 qMax = Capacity
  p1 = ANN param 1
 p2 = ANN param 2
 p3 = ANN param 3
 p4 = ANN param 4
 p5 = ANN param 5
 p6 = ANN param 6
 p7 = ANN param 7
 p8 = ANN param 8
 p9 = ANN param 9
 p10 = ANN param 10
 p11 = ANN param 11
 P_norm = Normative Vaule for Protien pore conc
 q_norm = Normative Vaule for Protien adsorbed conc
 S_norm = Normative Vaule for Salt pore conc
*/

namespace cadet
{

	namespace model
	{

		inline const char* TestParamHandler::identifier() CADET_NOEXCEPT { return "TEST_MULTI_COMPONENT_LANGMUIR"; }

		inline bool TestParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
		{
			if ((_kA.size() != _kD.size()) || (_kA.size() != _qMax.size()) || (_kA.size() != _p1.size()) || (_kA.size() != _p2.size())
				|| (_kA.size() != _p3.size()) || (_kA.size() != _p4.size()) || (_kA.size() != _p5.size()) || (_kA.size() != _p6.size())
				|| (_kA.size() != _p7.size()) || (_kA.size() != _p8.size()) || (_kA.size() != _p9.size()) || (_kA.size() != _p10.size())
				|| (_kA.size() != _p11.size()) || (_kA.size() != _P_norm.size()) || (_kA.size() != _q_norm.size()) || (_kA.size() != _S_norm.size()) || (_kA.size() < nComp))
				throw InvalidParameterException("MCL_KA, MCL_KD, and MCL_QMAX have to have the same size");

			return true;
		}

		inline const char* ExtTestParamHandler::identifier() CADET_NOEXCEPT { return "EXT_TEST_MULTI_COMPONENT_LANGMUIR"; }

		inline bool ExtTestParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
		{
			if ((_kA.size() != _kD.size()) || (_kA.size() != _qMax.size()) || (_kA.size() != _p1.size()) || (_kA.size() != _p2.size())
				|| (_kA.size() != _p3.size()) || (_kA.size() != _p4.size()) || (_kA.size() != _p5.size()) || (_kA.size() != _p6.size())
				|| (_kA.size() != _p7.size()) || (_kA.size() != _p8.size()) || (_kA.size() != _p9.size()) || (_kA.size() != _p10.size())
				|| (_kA.size() != _p11.size()) || (_kA.size() != _P_norm.size()) || (_kA.size() != _q_norm.size()) || (_kA.size() != _S_norm.size()) || (_kA.size() < nComp))
				throw InvalidParameterException("EXT_MCL_KA, EXT_MCL_KD, and EXT_MCL_QMAX have to have the same size");

			return true;
		}


		/**
		 * @brief Defines the multi component Langmuir binding model
		 * @details Implements the Langmuir adsorption model: \f[ \begin{align}
		 *              \frac{\mathrm{d}q_i}{\mathrm{d}t} &= k_{a,i} c_{p,i} q_{\text{max},i} \left( 1 - \sum_j \frac{q_j}{q_{\text{max},j}} \right) - k_{d,i} q_i
		 *          \end{align} \f]
		 *          Multiple bound states are not supported.
		 *          Components without bound state (i.e., non-binding components) are supported.
		 *
		 *          See @cite Langmuir1916.
		 * @tparam ParamHandler_t Type that can add support for external function dependence
		 */
		template <class ParamHandler_t>
		class TestBindingBase : public ParamHandlerBindingModelBase<ParamHandler_t>
		{
		public:

			TestBindingBase() { }
			virtual ~TestBindingBase() CADET_NOEXCEPT { }

			static const char* identifier() { return ParamHandler_t::identifier(); }

			virtual bool implementsAnalyticJacobian() const CADET_NOEXCEPT { return true; }

			CADET_BINDINGMODELBASE_BOILERPLATE

		protected:
			using ParamHandlerBindingModelBase<ParamHandler_t>::_paramHandler;
			using ParamHandlerBindingModelBase<ParamHandler_t>::_reactionQuasistationarity;
			using ParamHandlerBindingModelBase<ParamHandler_t>::_nComp;
			using ParamHandlerBindingModelBase<ParamHandler_t>::_nBoundStates;

			// In the follwing the class method the binding model mass transfer behavior is implemented. 
			// Apart from the mentioned places in the tutorial, do not modify anything in this function. 
			template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
			int fluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, StateType const* y,
				CpStateType const* yCp, ResidualType* res, LinearBufferAllocator workSpace) const
			{
				typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

				// Protein fluxes: -k_{a,i} * c_{p,i} * q_{max,i} * (1 - \sum_j q_j / q_{max,j}) + k_{d,i} * q_i
				ResidualType qSum = 1.0;
				unsigned int bndIdx = 0;
				for (int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states (bound state index bndIdx is not advanced)
					if (_nBoundStates[i] == 0)
						continue;

					qSum -= y[bndIdx] / static_cast<ParamType>(p->qMax[i]);

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
					const double p1 = static_cast<double>(p->p1[i]);
					const double p2 = static_cast<double>(p->p2[i]);
					const double p3 = static_cast<double>(p->p3[i]);
					const double p4 = static_cast<double>(p->p4[i]);
					const double p5 = static_cast<double>(p->p5[i]);
					const double p6 = static_cast<double>(p->p6[i]);
					const double p7 = static_cast<double>(p->p7[i]);
					const double p8 = static_cast<double>(p->p8[i]);
					const double p9 = static_cast<double>(p->p9[i]);
					const double p10 = static_cast<double>(p->p10[i]);
					const double p11 = static_cast<double>(p->p11[i]);
					const double P_norm = static_cast<double>(p->P_norm[i]);
					const double q_norm = static_cast<double>(p->q_norm[i]);
					const double S_norm = static_cast<double>(p->S_norm[i]);

					// const double CF1 = exp(-((p1 * yCp[i] / P_norm) + (p3 * y[bndIdx] / q_norm) + (p5 * yCp[0] / S_norm) + p7));
					// const double CF2 = exp(-((p2 * yCp[i] / P_norm) + (p4 * y[bndIdx] / q_norm) + (p6 * yCp[0] / S_norm) + p8));

					// double ANN = (p11 + (p9 / (1 + CF1)) + (p10 / (1 + CF2)));
					// ANN = (p11 + (p9 / (1 + exp(-((p1 * yCp[i] / P_norm) + (p3 * y[bndIdx] / q_norm) + (p5 * yCp[0] / S_norm) + p7)))) + (p10 / (1 + exp(-((p2 * yCp[i] / P_norm) + (p4 * y[bndIdx] / q_norm) + (p6 * yCp[0] / S_norm) + p8)))));

					res[bndIdx] = static_cast<ParamType>(p->kD[i]) * y[bndIdx] * abs((p11 + (p9 / (1 + exp(-((p1 * yCp[i] / P_norm) + (p3 * y[bndIdx] / q_norm) + (p5 * yCp[0] / S_norm) + p7)))) + (p10 / (1 + exp(-((p2 * yCp[i] / P_norm) + (p4 * y[bndIdx] / q_norm) + (p6 * yCp[0] / S_norm) + p8)))))) - static_cast<ParamType>(p->kA[i]) * yCp[i] * static_cast<ParamType>(p->qMax[i]) * qSum;

					// Next bound component
					++bndIdx;
				}

				return 0;
			}
			// In the follwing the class method the anlytical Jacobian of the binding model should be implemented.
			// Apart from the mentioned places in the tutorial, do not modify anything in this function. 
			template <typename RowIterator>
			void jacobianImpl(double t, unsigned int secIdx, const ColumnPosition& colPos, double const* y, double const* yCp, int offsetCp, RowIterator jac, LinearBufferAllocator workSpace) const
			{
				typename ParamHandler_t::ParamsHandle const p = _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

				// Protein fluxes: -k_{a,i} * c_{p,i} * q_{max,i} * (1 - \sum_j q_j / q_{max,j}) + k_{d,i} * q_i
				double qSum = 1.0;
				int bndIdx = 0;
				for (int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states (bound state index bndIdx is not advanced)
					if (_nBoundStates[i] == 0)
						continue;

					qSum -= y[bndIdx] / static_cast<double>(p->qMax[i]);

					// Next bound component
					++bndIdx;
				}

				bndIdx = 0;
				for (int i = 0; i < _nComp; ++i)
				{
					// Skip components without bound states (bound state index bndIdx is not advanced)
					if (_nBoundStates[i] == 0)
						continue;

					const double ka = static_cast<double>(p->kA[i]);
					const double kd = static_cast<double>(p->kD[i]);
					// ANN Params
					const double p1 = static_cast<double>(p->p1[i]);
					const double p2 = static_cast<double>(p->p2[i]);
					const double p3 = static_cast<double>(p->p3[i]);
					const double p4 = static_cast<double>(p->p4[i]);
					const double p5 = static_cast<double>(p->p5[i]);
					const double p6 = static_cast<double>(p->p6[i]);
					const double p7 = static_cast<double>(p->p7[i]);
					const double p8 = static_cast<double>(p->p8[i]);
					const double p9 = static_cast<double>(p->p9[i]);
					const double p10 = static_cast<double>(p->p10[i]);
					const double p11 = static_cast<double>(p->p11[i]);
					const double P_norm = static_cast<double>(p->P_norm[i]);
					const double q_norm = static_cast<double>(p->q_norm[i]);
					const double S_norm = static_cast<double>(p->S_norm[i]);
					// CF = Common Factor
					// const double CF1 = exp(-((p1 * yCp[i] / P_norm) + (p3 * y[bndIdx] / q_norm) + (p5 * yCp[0] / S_norm) + p7));
					// const double CF2 = exp(-((p2 * yCp[i] / P_norm) + (p4 * y[bndIdx] / q_norm) + (p6 * yCp[0] / S_norm) + p8));
					// const double ANN = (p11 + (p9 / (1 + CF1)) + (p10 / (1 + CF2)));
					// const double ANN_CF = (ANN / abs(ANN));

					// dres_i / dc_{p,i}
					// const double dc_pi_ANN = (((p10 * p2 * CF2) / (P_norm * pow(CF2 + 1, 2))) + ((p1 * p9 * CF1) / (P_norm * pow(CF1 + 1, 2)))) * ANN_CF;
					
					jac[i - bndIdx - offsetCp] = kd * y[bndIdx] * (((p10 * p2 * exp(-((p2 * yCp[i] / P_norm) + (p4 * y[bndIdx] / q_norm) + (p6 * yCp[0] / S_norm) + p8))) / (P_norm * pow(exp(-((p2 * yCp[i] / P_norm) + (p4 * y[bndIdx] / q_norm) + (p6 * yCp[0] / S_norm) + p8)) + 1, 2))) + ((p1 * p9 * exp(-((p1 * yCp[i] / P_norm) + (p3 * y[bndIdx] / q_norm) + (p5 * yCp[0] / S_norm) + p7))) / (P_norm * pow(exp(-((p1 * yCp[i] / P_norm) + (p3 * y[bndIdx] / q_norm) + (p5 * yCp[0] / S_norm) + p7)) + 1, 2)))) * ((p11 + (p9 / (1 + exp(-((p1 * yCp[i] / P_norm) + (p3 * y[bndIdx] / q_norm) + (p5 * yCp[0] / S_norm) + p7)))) + (p10 / (1 + exp(-((p2 * yCp[i] / P_norm) + (p4 * y[bndIdx] / q_norm) + (p6 * yCp[0] / S_norm) + p8))))) / abs((p11 + (p9 / (1 + exp(-((p1 * yCp[i] / P_norm) + (p3 * y[bndIdx] / q_norm) + (p5 * yCp[0] / S_norm) + p7)))) + (p10 / (1 + exp(-((p2 * yCp[i] / P_norm) + (p4 * y[bndIdx] / q_norm) + (p6 * yCp[0] / S_norm) + p8))))))) - ka * static_cast<double>(p->qMax[i]) * qSum;
					// Getting to c_{p,i}: -bndIdx takes us to q_0, another -offsetCp to c_{p,0} and a +i to c_{p,i}.
					//                     This means jac[i - bndIdx - offsetCp] corresponds to c_{p,i}.

					// dres_i / dc_{p,0}
					// const double dc_p0_ANN = (((p10 * p6 * CF2) / (S_norm * pow(CF2 + 1, 2))) + ((p5 * p9 * CF1) / (S_norm * pow(CF1 + 1, 2)))) * ANN_CF;

					jac[-bndIdx - offsetCp] = kd * y[bndIdx] * (((p10 * p6 * exp(-((p2 * yCp[i] / P_norm) + (p4 * y[bndIdx] / q_norm) + (p6 * yCp[0] / S_norm) + p8))) / (S_norm * pow(exp(-((p2 * yCp[i] / P_norm) + (p4 * y[bndIdx] / q_norm) + (p6 * yCp[0] / S_norm) + p8)) + 1, 2))) + ((p5 * p9 * exp(-((p1 * yCp[i] / P_norm) + (p3 * y[bndIdx] / q_norm) + (p5 * yCp[0] / S_norm) + p7))) / (S_norm * pow(exp(-((p1 * yCp[i] / P_norm) + (p3 * y[bndIdx] / q_norm) + (p5 * yCp[0] / S_norm) + p7)) + 1, 2)))) * ((p11 + (p9 / (1 + exp(-((p1 * yCp[i] / P_norm) + (p3 * y[bndIdx] / q_norm) + (p5 * yCp[0] / S_norm) + p7)))) + (p10 / (1 + exp(-((p2 * yCp[i] / P_norm) + (p4 * y[bndIdx] / q_norm) + (p6 * yCp[0] / S_norm) + p8))))) / abs((p11 + (p9 / (1 + exp(-((p1 * yCp[i] / P_norm) + (p3 * y[bndIdx] / q_norm) + (p5 * yCp[0] / S_norm) + p7)))) + (p10 / (1 + exp(-((p2 * yCp[i] / P_norm) + (p4 * y[bndIdx] / q_norm) + (p6 * yCp[0] / S_norm) + p8)))))));
					// Getting to c_{p,0}: -bndIdx takes us to q_0, another -offsetCp to c_{p,0}.
					//                     This means jac[bndIdx - offsetCp] corresponds to c_{p,0}.


					// Fill dres_i / dq_j
					int bndIdx2 = 0;
					for (int j = 0; j < _nComp; ++j)
					{
						// Skip components without bound states (bound state index bndIdx is not advanced)
						if (_nBoundStates[j] == 0)
							continue;

						// dres_i / dq_j
						jac[bndIdx2 - bndIdx] = static_cast<double>(p->kA[i]) * yCp[i] * static_cast<double>(p->qMax[i]) / static_cast<double>(p->qMax[j]);
						// Getting to q_j: -bndIdx takes us to q_0, another +bndIdx2 to q_j. This means jac[bndIdx2 - bndIdx] corresponds to q_j.

						++bndIdx2;
					}

					// Add to dres_i / dq_i
					// const double dc_qi_ANN = (((p10 * p4 * CF2) / (q_norm * pow(CF2 + 1, 2))) + ((p3 * p9 * CF1) / (q_norm * pow(CF1 + 1, 2)))) * ANN_CF;

					jac[0] += kd * abs((p11 + (p9 / (1 + exp(-((p1 * yCp[i] / P_norm) + (p3 * y[bndIdx] / q_norm) + (p5 * yCp[0] / S_norm) + p7)))) + (p10 / (1 + exp(-((p2 * yCp[i] / P_norm) + (p4 * y[bndIdx] / q_norm) + (p6 * yCp[0] / S_norm) + p8)))))) + kd * y[bndIdx] * (((p10 * p4 * exp(-((p2 * yCp[i] / P_norm) + (p4 * y[bndIdx] / q_norm) + (p6 * yCp[0] / S_norm) + p8))) / (q_norm * pow(exp(-((p2 * yCp[i] / P_norm) + (p4 * y[bndIdx] / q_norm) + (p6 * yCp[0] / S_norm) + p8)) + 1, 2))) + ((p3 * p9 * exp(-((p1 * yCp[i] / P_norm) + (p3 * y[bndIdx] / q_norm) + (p5 * yCp[0] / S_norm) + p7))) / (q_norm * pow(exp(-((p1 * yCp[i] / P_norm) + (p3 * y[bndIdx] / q_norm) + (p5 * yCp[0] / S_norm) + p7)) + 1, 2)))) * ((p11 + (p9 / (1 + exp(-((p1 * yCp[i] / P_norm) + (p3 * y[bndIdx] / q_norm) + (p5 * yCp[0] / S_norm) + p7)))) + (p10 / (1 + exp(-((p2 * yCp[i] / P_norm) + (p4 * y[bndIdx] / q_norm) + (p6 * yCp[0] / S_norm) + p8))))) / abs((p11 + (p9 / (1 + exp(-((p1 * yCp[i] / P_norm) + (p3 * y[bndIdx] / q_norm) + (p5 * yCp[0] / S_norm) + p7)))) + (p10 / (1 + exp(-((p2 * yCp[i] / P_norm) + (p4 * y[bndIdx] / q_norm) + (p6 * yCp[0] / S_norm) + p8)))))));

					// Advance to next flux and Jacobian row
					++bndIdx;
					++jac;
				}
			}
		};

		// Do not forget to make changes in the following lines of code according to the guidelines given in the tutorial.
		typedef TestBindingBase<TestParamHandler> TestBinding;
		typedef TestBindingBase<ExtTestParamHandler> ExternalTestBinding;

		namespace binding
		{
			void registerTestModel(std::unordered_map<std::string, std::function<model::IBindingModel* ()>>& bindings)
			{
				bindings[TestBinding::identifier()] = []() { return new TestBinding(); };
				bindings[ExternalTestBinding::identifier()] = []() { return new ExternalTestBinding(); };
			}
		}  // namespace binding

	}  // namespace model

}  // namespace cadet
