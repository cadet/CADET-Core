// SPDX-License-Identifier: AGPL-3.0-or-later
// =================================================================================
//  CADET
//
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Affero General Public
//  License v3.0 (or, at your option, any later version).
// =================================================================================

#include "model/binding/BindingModelBase.hpp"
#include "model/ExternalFunctionSupport.hpp"
#include "model/binding/BindingModelMacros.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "model/Parameters.hpp"
#include "LocalVector.hpp"
#include "SimulationTypes.hpp"
#include "AdUtils.hpp"
#include "model/binding/NeuralNetwork.h"

#include <functional>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <cmath>

/*<codegen>
{
    "name": "NeuralNetworkParamHandler",
    "externalName": "ExtNeuralNetworkParamHandler",
    "parameters":
        [
            { "type": "ScalarComponentDependentParameter", "varName": "kKin", "confName": "NN_KKIN"}
        ]
}
</codegen>*/

namespace cadet
{
    namespace model
    {

        inline const char* NeuralNetworkParamHandler::identifier() CADET_NOEXCEPT
        {
            return "NEURAL_NETWORK";
        }

        inline bool NeuralNetworkParamHandler::validateConfig(
            unsigned int nComp, unsigned int const* nBoundStates)
        {
            if (_kKin.size() < nComp)
                throw InvalidParameterException("NN_KKIN must have NCOMP entries");
            return true;
        }

        inline const char* ExtNeuralNetworkParamHandler::identifier() CADET_NOEXCEPT
        {
            return "EXT_NEURAL_NETWORK";
        }

        inline bool ExtNeuralNetworkParamHandler::validateConfig(
            unsigned int nComp, unsigned int const* nBoundStates)
        {
            if (_kKin.size() < nComp)
                throw InvalidParameterException("NN_KKIN must have NCOMP entries");
            return true;
        }

        template <class ParamHandler_t>
        class NeuralNetworkBindingBase : public ParamHandlerBindingModelBase<ParamHandler_t>
        {
        public:

            NeuralNetworkBindingBase()
                : _bias0(), _kernel0()
                , _bias1(), _kernel1()
                , _bias2(), _kernel2()
                , _normFactor()
                , _porosFactors()
                , _offsets()
                , _nLayers(1u)
                , _nNodes(0u)
                , _nInput(1u)
                , _nOutput(1u)
                , _annModels()
                , _cpCache()
                , _predCache()
                , _gradCache()
            {
            }

            virtual ~NeuralNetworkBindingBase() CADET_NOEXCEPT {}

            static const char* identifier() { return ParamHandler_t::identifier(); }
            virtual bool implementsAnalyticJacobian() const CADET_NOEXCEPT { return true; }

            virtual unsigned int workspaceSize(unsigned int nComp, unsigned int totalNumBoundStates,
                unsigned int const* nBoundStates) const CADET_NOEXCEPT
            {
                return ParamHandlerBindingModelBase<ParamHandler_t>::workspaceSize(
                    nComp, totalNumBoundStates, nBoundStates);
            }

            CADET_BINDINGMODELBASE_BOILERPLATE

        protected:

            using ParamHandlerBindingModelBase<ParamHandler_t>::_paramHandler;
            using ParamHandlerBindingModelBase<ParamHandler_t>::_reactionQuasistationarity;
            using ParamHandlerBindingModelBase<ParamHandler_t>::_nComp;
            using ParamHandlerBindingModelBase<ParamHandler_t>::_nBoundStates;

            // Network weights and biases — one set of vectors per bound state.
            // Inner vector follows col-major storage (matches training export).
            std::vector<std::vector<double>> _bias0;    // b1: (nNodes,)            per bound state
            std::vector<std::vector<double>> _kernel0;  // W1: (nNodes x nInput)    per bound state
            std::vector<std::vector<double>> _bias1;    // b2: (nNodes,) or (1,)    per bound state
            std::vector<std::vector<double>> _kernel1;  // W2: (nNodes x nNodes) or (1 x nNodes) per bound state
            std::vector<std::vector<double>> _bias2;    // b3: (1,)      — 2-layer only, per bound state
            std::vector<std::vector<double>> _kernel2;  // W3: (1 x nNodes) — 2-layer only, per bound state

            std::vector<double> _normFactor;   // shared normalisation factor per input dimension (length _nInput)
            std::vector<double> _porosFactors; // porosity scaling — one per bound state
            std::vector<double> _offsets;      // ANN(x=0) — one per bound state, subtracted for zero intercept

            unsigned int _nLayers;   // number of hidden layers (1 or 2), shared across all bound states
            unsigned int _nNodes;    // nodes per hidden layer,            shared across all bound states
            unsigned int _nInput;    // input dimensionality (== nComp)
            unsigned int _nOutput;   // output dimensionality (always 1 — scalar isotherm per bound state)

            std::vector<std::unique_ptr<ClassANN>> _annModels;  // one model per bound state

            // Pre-allocated per-call buffers — reused across bound states (sequential evaluation)
            mutable std::vector<double> _cpCache;    // normalised c_p input, length _nInput
            mutable std::vector<double> _predCache;  // ANN output, length _nOutput (== 1)
            mutable std::vector<double> _gradCache;  // ANN gradient dy/dx, length _nInput

            virtual bool configureImpl(IParameterProvider& paramProvider,
                UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
            {
                // Base class handles _paramHandler setup (kKin registration etc.)
                const bool result = ParamHandlerBindingModelBase<ParamHandler_t>::configureImpl(
                    paramProvider, unitOpIdx, parTypeIdx);

                // --- Read shared scalar configuration ---
                _nLayers = static_cast<unsigned int>(paramProvider.getInt("NLAYERS"));
                _nNodes = static_cast<unsigned int>(paramProvider.getInt("NNODES"));

                if (_nLayers != 1u && _nLayers != 2u)
                    throw InvalidParameterException(
                        "NEURAL_NETWORK binding: NLAYERS must be 1 or 2. Got: "
                        + std::to_string(_nLayers));

                if (_nNodes == 0u)
                    throw InvalidParameterException("NEURAL_NETWORK binding: NNODES must be > 0.");

                // Shared normalisation: one factor per input component
                _normFactor = paramProvider.getDoubleArray("NORM_FACTOR");

                // Input/output dimensionality: each bound state's ANN takes the full
                // c_p vector and predicts a scalar solid concentration.
                _nInput = _nComp;
                _nOutput = 1u;

                if (_normFactor.size() < _nInput)
                    throw InvalidParameterException(
                        "NEURAL_NETWORK binding: NORM_FACTOR must have at least NCOMP ("
                        + std::to_string(_nInput) + ") entries.");

                // --- Validate bound state count and compute total ---
                unsigned int totalNumBoundStates = 0u;
                for (unsigned int i = 0u; i < _nComp; ++i)
                {
                    if (_nBoundStates[i] != 1u)
                        throw InvalidParameterException(
                            "NEURAL_NETWORK binding requires exactly one bound state per component "
                            "(component " + std::to_string(i) + " has "
                            + std::to_string(_nBoundStates[i]) + ").");
                    ++totalNumBoundStates;
                }

                // --- Clear and reserve per-bound-state containers ---
                _bias0.clear();    _bias0.reserve(totalNumBoundStates);
                _kernel0.clear();  _kernel0.reserve(totalNumBoundStates);
                _bias1.clear();    _bias1.reserve(totalNumBoundStates);
                _kernel1.clear();  _kernel1.reserve(totalNumBoundStates);
                _bias2.clear();    _bias2.reserve(totalNumBoundStates);
                _kernel2.clear();  _kernel2.reserve(totalNumBoundStates);
                _porosFactors.clear(); _porosFactors.reserve(totalNumBoundStates);
                _offsets.clear();  _offsets.reserve(totalNumBoundStates);
                _annModels.clear(); _annModels.reserve(totalNumBoundStates);

                // --- Read per-bound-state weights ---
                // Scope layout:
                //   adsorption/bound_state_N/layer_0: {KERNEL, BIAS}   W1, b1
                //   adsorption/bound_state_N/layer_1: {KERNEL, BIAS}   W2, b2
                //   adsorption/bound_state_N/layer_2: {KERNEL, BIAS}   W3, b3  (2-layer only)
                //   adsorption/bound_state_N/POROSITY_FACTOR
                for (unsigned int bndIdx = 0u; bndIdx < totalNumBoundStates; ++bndIdx)
                {
					const std::string scopeName = std::format("bound_state_{:03}", bndIdx);
                    paramProvider.pushScope(scopeName);

                    const double pf = paramProvider.getDouble("POROSITY_FACTOR");
                    _porosFactors.push_back(pf);

                    paramProvider.pushScope("layer_0");
                    _kernel0.push_back(paramProvider.getDoubleArray("KERNEL"));
                    _bias0.push_back(paramProvider.getDoubleArray("BIAS"));
                    paramProvider.popScope();

                    paramProvider.pushScope("layer_1");
                    _kernel1.push_back(paramProvider.getDoubleArray("KERNEL"));
                    _bias1.push_back(paramProvider.getDoubleArray("BIAS"));
                    paramProvider.popScope();

                    if (_nLayers == 2u)
                    {
                        paramProvider.pushScope("layer_2");
                        _kernel2.push_back(paramProvider.getDoubleArray("KERNEL"));
                        _bias2.push_back(paramProvider.getDoubleArray("BIAS"));
                        paramProvider.popScope();
                    }

                    paramProvider.popScope();  // bound_state_N

                    // --- Validate weight sizes for this bound state ---
                    // W1: (nNodes x nInput) -> nNodes * nInput elements
                    if (_kernel0.back().size() != _nNodes * _nInput)
                        throw InvalidParameterException(
                            "NEURAL_NETWORK binding: " + scopeName
                            + "/layer_0 KERNEL size must be NNODES * NCOMP = "
                            + std::to_string(_nNodes * _nInput));
                    if (_bias0.back().size() != _nNodes)
                        throw InvalidParameterException(
                            "NEURAL_NETWORK binding: " + scopeName
                            + "/layer_0 BIAS size must be NNODES = "
                            + std::to_string(_nNodes));

                    if (_nLayers == 2u)
                    {
                        // W2: (nNodes x nNodes)
                        if (_kernel1.back().size() != _nNodes * _nNodes)
                            throw InvalidParameterException(
                                "NEURAL_NETWORK binding: " + scopeName
                                + "/layer_1 KERNEL size must be NNODES * NNODES = "
                                + std::to_string(_nNodes * _nNodes));
                        if (_bias1.back().size() != _nNodes)
                            throw InvalidParameterException(
                                "NEURAL_NETWORK binding: " + scopeName
                                + "/layer_1 BIAS size must be NNODES = "
                                + std::to_string(_nNodes));
                        // W3: (nOutput x nNodes) == (1 x nNodes)
                        if (_kernel2.back().size() != _nOutput * _nNodes)
                            throw InvalidParameterException(
                                "NEURAL_NETWORK binding: " + scopeName
                                + "/layer_2 KERNEL size must be NOUTPUT * NNODES = "
                                + std::to_string(_nOutput * _nNodes));
                        if (_bias2.back().size() != _nOutput)
                            throw InvalidParameterException(
                                "NEURAL_NETWORK binding: " + scopeName
                                + "/layer_2 BIAS size must be NOUTPUT = "
                                + std::to_string(_nOutput));
                    }
                    else  // 1-layer: layer_1 is the output layer
                    {
                        // W2: (nOutput x nNodes)
                        if (_kernel1.back().size() != _nOutput * _nNodes)
                            throw InvalidParameterException(
                                "NEURAL_NETWORK binding: " + scopeName
                                + "/layer_1 KERNEL size must be NOUTPUT * NNODES = "
                                + std::to_string(_nOutput * _nNodes));
                        if (_bias1.back().size() != _nOutput)
                            throw InvalidParameterException(
                                "NEURAL_NETWORK binding: " + scopeName
                                + "/layer_1 BIAS size must be NOUTPUT = "
                                + std::to_string(_nOutput));
                    }

                    // Build ANN model for this bound state
                    _annModels.push_back(
                        std::make_unique<ClassANN>(_nLayers, _nNodes, _nInput, _nOutput));
                }

                // --- Pre-allocate shared runtime buffers ---
                _cpCache.assign(_nInput, 0.0);
                _predCache.assign(_nOutput, 0.0);
                _gradCache.assign(_nInput, 0.0);

                // --- Compute per-bound-state offsets = ANN_i(x=0) for zero-intercept correction ---
                // Enforces q_i(c_p = 0) = 0. Verify against training data before use.
                std::fill(_cpCache.begin(), _cpCache.end(), 0.0);
                for (unsigned int bndIdx = 0u; bndIdx < totalNumBoundStates; ++bndIdx)
                {
                    _annForward(bndIdx, _cpCache.data(), _predCache.data());
                    _offsets.push_back(_predCache[0]);
                }

                return result;
            }

            // -------------------------------------------------------------------------
            // Flux: res[bndIdx] = kKin[i] * (q[bndIdx] - q_ANN_bndIdx(c_p))
            //
            // Each bound state has its own ANN trained on the full c_p input vector.
            // At steady state (res = 0): q[bndIdx] = q_ANN_bndIdx(c_p).
            // AD types in yCp are cast to double — Jacobian provided analytically.
            // -------------------------------------------------------------------------
            template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
            int fluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
                StateType const* y, CpStateType const* yCp, ResidualType* res,
                LinearBufferAllocator workSpace) const
            {
                typename ParamHandler_t::ParamsHandle const p =
                    _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

                // Normalise the full c_p vector — shared across all bound-state ANNs
                for (unsigned int d = 0u; d < _nInput; ++d)
                    _cpCache[d] = static_cast<double>(yCp[d]) * _normFactor[d];

                unsigned int bndIdx = 0u;
                for (unsigned int i = 0u; i < _nComp; ++i)
                {
                    if (_nBoundStates[i] == 0u)
                        continue;

                    // Forward pass for this bound state's ANN
                    _annForward(bndIdx, _cpCache.data(), _predCache.data());

                    // q_iso = (ANN_bndIdx(c_p) - ANN_bndIdx(0)) * porosity_factor
                    const double qML = (_predCache[0] - _offsets[bndIdx]) * _porosFactors[bndIdx];

                    // res = kKin * (q - q_iso)
                    res[bndIdx] = static_cast<ParamType>(p->kKin[i])
                        * (static_cast<ResidualType>(y[bndIdx])
                            - static_cast<ResidualType>(qML));

                    ++bndIdx;
                }

                return 0;
            }

            // -------------------------------------------------------------------------
            // Jacobian of res[bndIdx] = kKin[i] * (q[bndIdx] - q_ANN_bndIdx(c_p))
            //
            // dres[bndIdx]/dq[bndIdx]  =  kKin[i]
            // dres[bndIdx]/dc_p[d]     = -kKin[i] * d(q_ANN_bndIdx)/d(c_p_norm[d])
            //                                      * normFactor[d] * porosFactors[bndIdx]
            //
            // Each bound state uses its own ANN (and thus its own gradient).
            // -------------------------------------------------------------------------
            template <typename RowIterator>
            void jacobianImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
                double const* y, double const* yCp, int offsetCp, RowIterator jac,
                LinearBufferAllocator workSpace) const
            {
                typename ParamHandler_t::ParamsHandle const p =
                    _paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

                // Normalise c_p
                for (unsigned int d = 0u; d < _nInput; ++d)
                    _cpCache[d] = yCp[d] * _normFactor[d];

                unsigned int bndIdx = 0u;
                for (int i = 0; i < static_cast<int>(_nComp); ++i)
                {
                    if (_nBoundStates[i] == 0u)
                        continue;

                    const double kkin = static_cast<double>(p->kKin[i]);

                    // Compute gradient d(q_ANN_bndIdx)/d(c_p_normalised) for this bound state
                    _annJacobian(bndIdx, _cpCache.data(), _gradCache.data());

                    // dres[bndIdx]/dc_p[d] for each input dimension d
                    for (int d = 0; d < static_cast<int>(_nInput); ++d)
                        jac[d - static_cast<int>(bndIdx) - offsetCp] =
                        -kkin * _gradCache[d] * _normFactor[d] * _porosFactors[bndIdx];

                    // dres[bndIdx]/dq[bndIdx] = kKin
                    jac[0] = kkin;

                    ++bndIdx;
                    ++jac;
                }
            }

        private:

            // -------------------------------------------------------------------------
            // Forward pass dispatcher for bound state bndIdx.
            // Input x must already be normalised.
            // -------------------------------------------------------------------------
            void _annForward(unsigned int bndIdx, const double* x, double* pred) const
            {
                if (_nLayers == 1u)
                    _annModels[bndIdx]->forward_single_layer(
                        x,
                        _kernel0[bndIdx].data(), _bias0[bndIdx].data(),
                        _kernel1[bndIdx].data(), _bias1[bndIdx].data(),
                        pred);
                else
                    _annModels[bndIdx]->forward_two_layers(
                        x,
                        _kernel0[bndIdx].data(), _bias0[bndIdx].data(),
                        _kernel1[bndIdx].data(), _bias1[bndIdx].data(),
                        _kernel2[bndIdx].data(), _bias2[bndIdx].data(),
                        pred);
            }

            // -------------------------------------------------------------------------
            // Jacobian dispatcher for bound state bndIdx.
            // Input x must already be normalised.
            // -------------------------------------------------------------------------
            void _annJacobian(unsigned int bndIdx, const double* x, double* grad) const
            {
                if (_nLayers == 1u)
                    _annModels[bndIdx]->jacobian_single_layer(
                        x,
                        _kernel0[bndIdx].data(), _bias0[bndIdx].data(),
                        _kernel1[bndIdx].data(),
                        grad);
                else
                    _annModels[bndIdx]->jacobian_two_layers(
                        x,
                        _kernel0[bndIdx].data(), _bias0[bndIdx].data(),
                        _kernel1[bndIdx].data(), _bias1[bndIdx].data(),
                        _kernel2[bndIdx].data(),
                        grad);
            }
        };

        typedef NeuralNetworkBindingBase<NeuralNetworkParamHandler>    NeuralNetworkBinding;
        typedef NeuralNetworkBindingBase<ExtNeuralNetworkParamHandler> ExternalNeuralNetworkBinding;

        namespace binding
        {
            void registerNeuralNetworkModel(
                std::unordered_map<std::string, std::function<model::IBindingModel* ()>>& bindings)
            {
                bindings[NeuralNetworkBinding::identifier()] = []() { return new NeuralNetworkBinding(); };
                bindings[ExternalNeuralNetworkBinding::identifier()] = []() { return new ExternalNeuralNetworkBinding(); };
            }
        }  // namespace binding

    }  // namespace model
}  // namespace cadet