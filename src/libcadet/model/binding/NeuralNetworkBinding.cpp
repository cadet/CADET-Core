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
                : _layerKernels()
                , _layerBiases()
                , _kernelPtrs()
                , _biasPtrs()
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

            // Network weights and biases — indexed [bndIdx][layerIdx].
            // layerIdx runs 0.._nLayers (inclusive): indices 0.._nLayers-1 are hidden
            // layers, index _nLayers is the output layer.
            // Inner vector follows col-major storage (matches training export).
            std::vector<std::vector<std::vector<double>>> _layerKernels; // [bndIdx][layerIdx]
            std::vector<std::vector<std::vector<double>>> _layerBiases;  // [bndIdx][layerIdx]

            // Pre-allocated raw-pointer caches used in the hot path — avoids heap
            // allocation inside fluxImpl / jacobianImpl.
            std::vector<std::vector<const double*>> _kernelPtrs; // [bndIdx][layerIdx]
            std::vector<std::vector<const double*>> _biasPtrs;   // [bndIdx][layerIdx]

            std::vector<double> _normFactor;   // normalisation factor per input dimension (length _nInput)
            std::vector<double> _porosFactors; // porosity scaling — one per bound state
            std::vector<double> _offsets;      // ANN(x=0) — one per bound state, subtracted for zero intercept

            unsigned int _nLayers;   // number of hidden layers (>= 1), shared across all bound states
            unsigned int _nNodes;    // nodes per hidden layer,         shared across all bound states
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
                _nNodes  = static_cast<unsigned int>(paramProvider.getInt("NNODES"));

                if (_nLayers < 1u)
                    throw InvalidParameterException(
                        "NEURAL_NETWORK binding: NLAYERS must be >= 1. Got: "
                        + std::to_string(_nLayers));

                if (_nNodes == 0u)
                    throw InvalidParameterException("NEURAL_NETWORK binding: NNODES must be > 0.");

                // Input/output dimensionality: each bound state's ANN takes the full
                // c_p vector and predicts a scalar solid concentration.
                _nInput  = _nComp;
                _nOutput = 1u;

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

                const unsigned int nTotalLayers = _nLayers + 1u; // hidden layers + output layer

                // --- Clear and reserve containers ---
                _layerKernels.clear(); _layerKernels.reserve(totalNumBoundStates);
                _layerBiases.clear();  _layerBiases.reserve(totalNumBoundStates);
                _kernelPtrs.clear();   _kernelPtrs.reserve(totalNumBoundStates);
                _biasPtrs.clear();     _biasPtrs.reserve(totalNumBoundStates);
                _porosFactors.clear(); _porosFactors.reserve(totalNumBoundStates);
                _offsets.clear();      _offsets.reserve(totalNumBoundStates);
                _annModels.clear();    _annModels.reserve(totalNumBoundStates);
                _normFactor.clear();   _normFactor.reserve(totalNumBoundStates);

                // --- Read parameters for each NN / bound state ---
                for (unsigned int bndIdx = 0u; bndIdx < totalNumBoundStates; ++bndIdx)
                {
                    const std::string scopeName = std::format("bound_state_{:03}", bndIdx);
                    paramProvider.pushScope(scopeName);

                    _porosFactors.push_back(paramProvider.getDouble("POROSITY_FACTOR"));
                    _normFactor.push_back(paramProvider.getDouble("NORM_FACTOR"));

                    _layerKernels.emplace_back();
                    _layerBiases.emplace_back();
                    _layerKernels.back().reserve(nTotalLayers);
                    _layerBiases.back().reserve(nTotalLayers);

                    for (unsigned int l = 0u; l < nTotalLayers; ++l)
                    {
                        paramProvider.pushScope(std::format("layer_{}", l));
                        _layerKernels.back().push_back(paramProvider.getDoubleArray("KERNEL"));
                        _layerBiases.back().push_back(paramProvider.getDoubleArray("BIAS"));
                        paramProvider.popScope();
                    }

                    paramProvider.popScope();  // bound_state_N

                    // --- Validate weight sizes for this bound state ---
                    // Hidden layers 0.._nLayers-1
                    for (unsigned int l = 0u; l < _nLayers; ++l)
                    {
                        const unsigned int prevSize = (l == 0u) ? _nInput : _nNodes;
                        if (_layerKernels.back()[l].size() != _nNodes * prevSize)
                            throw InvalidParameterException(
                                "NEURAL_NETWORK binding: " + scopeName
                                + "/layer_" + std::to_string(l) + " KERNEL size must be "
                                + std::to_string(_nNodes) + " * " + std::to_string(prevSize)
                                + " = " + std::to_string(_nNodes * prevSize));
                        if (_layerBiases.back()[l].size() != _nNodes)
                            throw InvalidParameterException(
                                "NEURAL_NETWORK binding: " + scopeName
                                + "/layer_" + std::to_string(l) + " BIAS size must be NNODES = "
                                + std::to_string(_nNodes));
                    }
                    // Output layer _nLayers
                    if (_layerKernels.back()[_nLayers].size() != _nOutput * _nNodes)
                        throw InvalidParameterException(
                            "NEURAL_NETWORK binding: " + scopeName
                            + "/layer_" + std::to_string(_nLayers) + " KERNEL size must be "
                            + std::to_string(_nOutput) + " * " + std::to_string(_nNodes)
                            + " = " + std::to_string(_nOutput * _nNodes));
                    if (_layerBiases.back()[_nLayers].size() != _nOutput)
                        throw InvalidParameterException(
                            "NEURAL_NETWORK binding: " + scopeName
                            + "/layer_" + std::to_string(_nLayers) + " BIAS size must be NOUTPUT = "
                            + std::to_string(_nOutput));

                    _annModels.push_back(
                        std::make_unique<ClassANN>(_nLayers, _nNodes, _nInput, _nOutput));
                }

                // --- Build pre-allocated pointer caches (avoids allocation in hot path) ---
                _kernelPtrs.resize(totalNumBoundStates);
                _biasPtrs.resize(totalNumBoundStates);
                for (unsigned int bndIdx = 0u; bndIdx < totalNumBoundStates; ++bndIdx)
                {
                    _kernelPtrs[bndIdx].resize(nTotalLayers);
                    _biasPtrs[bndIdx].resize(nTotalLayers);
                    for (unsigned int l = 0u; l < nTotalLayers; ++l)
                    {
                        _kernelPtrs[bndIdx][l] = _layerKernels[bndIdx][l].data();
                        _biasPtrs[bndIdx][l]   = _layerBiases[bndIdx][l].data();
                    }
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
            // Forward pass for bound state bndIdx.
            // Uses pre-allocated pointer caches — no heap allocation.
            // Input x must already be normalised.
            // -------------------------------------------------------------------------
            void _annForward(unsigned int bndIdx, const double* x, double* pred) const
            {
                _annModels[bndIdx]->forward(_kernelPtrs[bndIdx], _biasPtrs[bndIdx], x, pred);
            }

            // -------------------------------------------------------------------------
            // Jacobian for bound state bndIdx.
            // Uses pre-allocated pointer caches — no heap allocation.
            // Input x must already be normalised.
            // -------------------------------------------------------------------------
            void _annJacobian(unsigned int bndIdx, const double* x, double* grad) const
            {
                _annModels[bndIdx]->jacobian(_kernelPtrs[bndIdx], _biasPtrs[bndIdx], x, grad);
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