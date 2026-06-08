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
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <cmath>

/*<codegen>
{
	"name": "MachineLearningParamHandler",
	"externalName": "ExtMachineLearningParamHandler",
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

inline const char* MachineLearningParamHandler::identifier() CADET_NOEXCEPT
{
	return "NEURAL_NETWORK";
}

inline bool MachineLearningParamHandler::validateConfig(
	unsigned int nComp, unsigned int const* nBoundStates)
{
	if (_kKin.size() < nComp)
		throw InvalidParameterException("NN_KKIN must have NCOMP entries");
	return true;
}

inline const char* ExtMachineLearningParamHandler::identifier() CADET_NOEXCEPT
{
	return "EXT_NEURAL_NETWORK";
}

inline bool ExtMachineLearningParamHandler::validateConfig(
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
		, _porosFactor(1.0)
		, _offset(0.0)
		, _nLayers(1u)
		, _nNodes(0u)
		, _nInput(1u)
		, _nOutput(1u)
		, _annModel()
		, _cpCache()
		, _predCache()
		, _gradCache()
	{
	}

	virtual ~NeuralNetworkBindingBase() CADET_NOEXCEPT { }

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

	// Network weights and biases (col-major storage, matches training export)
	std::vector<double> _bias0;    // b1: (nNodes,)
	std::vector<double> _kernel0;  // W1: (nNodes x nInput) col-major
	std::vector<double> _bias1;    // b2: (nNodes,) for 2-layer; (nOutput,) for 1-layer output
	std::vector<double> _kernel1;  // W2: (nNodes x nNodes) col-major for 2-layer;
	                               //     (nOutput x nNodes) col-major for 1-layer output
	std::vector<double> _bias2;    // b3: (nOutput,) — 2-layer only
	std::vector<double> _kernel2;  // W3: (nOutput x nNodes) col-major — 2-layer only

	std::vector<double> _normFactor;  // normalisation factor per input dimension
	double              _porosFactor; // porosity scaling applied to prediction
	double              _offset;      // prediction at x=0, subtracted for zero intercept

	unsigned int _nLayers;   // number of hidden layers (1 or 2)
	unsigned int _nNodes;    // nodes per hidden layer
	unsigned int _nInput;    // input dimensionality (== nComp, single-component => 1)
	unsigned int _nOutput;   // output dimensionality (always 1 for single-component)

	std::unique_ptr<ClassANN>  _annModel;

	// Pre-allocated per-call buffers (mutable — used in const flux/jacobian methods)
	mutable std::vector<double> _cpCache;    // normalised c_p input, length _nInput
	mutable std::vector<double> _predCache;  // ANN output, length _nOutput
	mutable std::vector<double> _gradCache;  // ANN gradient dy/dx, length _nInput

	virtual bool configureImpl(IParameterProvider& paramProvider,
		UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
	{
		// Base class handles _paramHandler setup (kKin registration etc.)
		const bool result = ParamHandlerBindingModelBase<ParamHandler_t>::configureImpl(
			paramProvider, unitOpIdx, parTypeIdx);

		// --- Read scalar configuration ---
		_nLayers = static_cast<unsigned int>(paramProvider.getInt("NLAYERS"));
		_nNodes  = static_cast<unsigned int>(paramProvider.getInt("NNODES"));

		if (_nLayers != 1u && _nLayers != 2u)
			throw InvalidParameterException(
				"NEURAL_NETWORK binding: NLAYERS must be 1 or 2. Got: "
				+ std::to_string(_nLayers));

		if (_nNodes == 0u)
			throw InvalidParameterException("NEURAL_NETWORK binding: NNODES must be > 0.");
		// Normalisation: one factor per input component
		_normFactor = paramProvider.getDoubleArray("NORM_FACTOR");

		// Porosity factor: single scalar
		{
			const std::vector<double> pf = paramProvider.getDoubleArray("POROSITY_FACTOR");
			if (pf.empty())
				throw InvalidParameterException("NEURAL_NETWORK binding: POROSITY_FACTOR must not be empty.");
			_porosFactor = pf[0];
		}

		// For now single-component only — each component's c_p is one input dimension.
		// _nInput == _nComp; _nOutput == 1 (scalar isotherm).
		_nInput  = _nComp;
		_nOutput = 1u;

		if (_normFactor.size() < _nInput)
			throw InvalidParameterException(
				"NEURAL_NETWORK binding: NORM_FACTOR must have at least NCOMP entries.");

		// --- Read weights from HDF5 scope ---
		// Scope layout:
		//   adsorption/layer_0: {KERNEL, BIAS}   W1, b1
		//   adsorption/layer_1: {KERNEL, BIAS}   W2, b2  (hidden for 2-layer; output for 1-layer)
		//   adsorption/layer_2: {KERNEL, BIAS}   W3, b3  (output — 2-layer only)

		paramProvider.pushScope("layer_0");
		_kernel0 = paramProvider.getDoubleArray("KERNEL");
		_bias0   = paramProvider.getDoubleArray("BIAS");
		paramProvider.popScope();  // layer_0

		paramProvider.pushScope("layer_1");
		_kernel1 = paramProvider.getDoubleArray("KERNEL");
		_bias1   = paramProvider.getDoubleArray("BIAS");
		paramProvider.popScope();  // layer_1

		if (_nLayers == 2u)
		{
			paramProvider.pushScope("layer_2");
			_kernel2 = paramProvider.getDoubleArray("KERNEL");
			_bias2   = paramProvider.getDoubleArray("BIAS");
			paramProvider.popScope();  // layer_2
		}

		// --- Validate weight sizes ---
		// W1: (nNodes x nInput) -> nNodes * nInput elements
		if (_kernel0.size() != _nNodes * _nInput)
			throw InvalidParameterException(
				"NEURAL_NETWORK binding: layer_0 KERNEL size must be NODES * NCOMP = "
				+ std::to_string(_nNodes * _nInput));
		if (_bias0.size() != _nNodes)
			throw InvalidParameterException(
				"NEURAL_NETWORK binding: layer_0 BIAS size must be NODES = "
				+ std::to_string(_nNodes));

		if (_nLayers == 2u)
		{
			// W2: (nNodes x nNodes)
			if (_kernel1.size() != _nNodes * _nNodes)
				throw InvalidParameterException(
					"NEURAL_NETWORK binding: layer_1 KERNEL size must be NODES * NODES = "
					+ std::to_string(_nNodes * _nNodes));
			if (_bias1.size() != _nNodes)
				throw InvalidParameterException(
					"NEURAL_NETWORK binding: layer_1 BIAS size must be NODES = "
					+ std::to_string(_nNodes));
			// W3: (nOutput x nNodes)
			if (_kernel2.size() != _nOutput * _nNodes)
				throw InvalidParameterException(
					"NEURAL_NETWORK binding: layer_2 KERNEL size must be NOUTPUT * NODES = "
					+ std::to_string(_nOutput * _nNodes));
			if (_bias2.size() != _nOutput)
				throw InvalidParameterException(
					"NEURAL_NETWORK binding: layer_2 BIAS size must be NOUTPUT = "
					+ std::to_string(_nOutput));
		}
		else  // 1-layer: layer_1 is the output layer
		{
			// W2: (nOutput x nNodes)
			if (_kernel1.size() != _nOutput * _nNodes)
				throw InvalidParameterException(
					"NEURAL_NETWORK binding: layer_1 KERNEL size must be NOUTPUT * NODES = "
					+ std::to_string(_nOutput * _nNodes));
			if (_bias1.size() != _nOutput)
				throw InvalidParameterException(
					"NEURAL_NETWORK binding: layer_1 BIAS size must be NOUTPUT = "
					+ std::to_string(_nOutput));
		}

		// --- Build ANN model ---
		_annModel = std::make_unique<ClassANN>(_nLayers, _nNodes, _nInput, _nOutput);

		// --- Pre-allocate runtime buffers ---
		_cpCache.assign(_nInput, 0.0);
		_predCache.assign(_nOutput, 0.0);
		_gradCache.assign(_nInput, 0.0);

		// --- Compute offset = ANN(x=0) for zero-intercept correction ---
		// Enforces q(c_p = 0) = 0. Only valid if your isotherm passes through
		// the origin. Verify against training data before use.
		std::fill(_cpCache.begin(), _cpCache.end(), 0.0);
		_annForward(_cpCache.data(), _predCache.data());
		_offset = _predCache[0];

		return result;
	}

	// -------------------------------------------------------------------------
	// Flux: res[bndIdx] = kKin[i] * (q[bndIdx] - q_ANN(c_p))
	//
	// At steady state (res = 0): q = q_ANN(c_p).
	// AD types in yCp are cast to double — Jacobian provided analytically.
	// -------------------------------------------------------------------------
	template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
	int fluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
		StateType const* y, CpStateType const* yCp, ResidualType* res,
		LinearBufferAllocator workSpace) const
	{
		typename ParamHandler_t::ParamsHandle const p =
			_paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Normalise c_p and run forward pass
		for (unsigned int d = 0; d < _nInput; ++d)
			_cpCache[d] = static_cast<double>(yCp[d]) * _normFactor[d];

		_annForward(_cpCache.data(), _predCache.data());

		// Apply offset and porosity scaling
		// q_iso = (ANN(c_p) - ANN(0)) * porosity_factor
		const double qML = (_predCache[0] - _offset) * _porosFactor;

		unsigned int bndIdx = 0;
		for (unsigned int i = 0; i < _nComp; ++i)
		{
			if (_nBoundStates[i] == 0)
				continue;

			// res = kKin * (q - q_iso)
			// At equilibrium: q = q_iso(c_p)
			res[bndIdx] = static_cast<ParamType>(p->kKin[i])
				* (static_cast<ResidualType>(y[bndIdx])
				   - static_cast<ResidualType>(qML));

			++bndIdx;
		}

		return 0;
	}

	// -------------------------------------------------------------------------
	// Jacobian of res[bndIdx] = kKin[i] * (q[bndIdx] - q_ANN(c_p))
	//
	// dres/dq[bndIdx]  =  kKin[i]                           -> jac[0]
	// dres/dc_p[d]     = -kKin[i] * d(q_ANN)/dc_p[d] * normFactor[d] * porosFactor
	//
	// The chain rule picks up normFactor because the ANN input is
	// c_p_normalised = c_p * normFactor, so:
	//   d(q_ANN)/dc_p[d] = d(q_ANN)/d(c_p_norm[d]) * normFactor[d] * porosFactor
	//
	// Jacobian index (follows CADET banded layout, same as Langmuir):
	//   jac[d - bndIdx - offsetCp]  corresponds to c_p[d]
	//   jac[0]                      corresponds to q[bndIdx]
	// -------------------------------------------------------------------------
	template <typename RowIterator>
	void jacobianImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
		double const* y, double const* yCp, int offsetCp, RowIterator jac,
		LinearBufferAllocator workSpace) const
	{
		typename ParamHandler_t::ParamsHandle const p =
			_paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Normalise c_p
		for (unsigned int d = 0; d < _nInput; ++d)
			_cpCache[d] = yCp[d] * _normFactor[d];

		// Compute gradient d(q_ANN)/d(c_p_normalised) — length _nInput
		_annJacobian(_cpCache.data(), _gradCache.data());

		unsigned int bndIdx = 0;
		for (int i = 0; i < static_cast<int>(_nComp); ++i)
		{
			if (_nBoundStates[i] == 0)
				continue;

			const double kkin = static_cast<double>(p->kKin[i]);

			// dres/dc_p[d] = -kKin * d(q_ANN)/d(c_p_norm[d]) * normFactor[d] * porosFactor
			for (int d = 0; d < static_cast<int>(_nInput); ++d)
				jac[d - static_cast<int>(bndIdx) - offsetCp] =
					-kkin * _gradCache[d] * _normFactor[d] * _porosFactor;

			// dres/dq[bndIdx] = kKin
			jac[0] = kkin;

			++bndIdx;
			++jac;
		}
	}

private:

	// -------------------------------------------------------------------------
	// Forward pass dispatcher — writes ANN output to pred (length _nOutput).
	// Input x must already be normalised.
	// -------------------------------------------------------------------------
	void _annForward(const double* x, double* pred) const
	{
		if (_nLayers == 1u)
			_annModel->forward_single_layer(
				x,
				_kernel0.data(), _bias0.data(),
				_kernel1.data(), _bias1.data(),
				pred);
		else
			_annModel->forward_two_layers(
				x,
				_kernel0.data(), _bias0.data(),
				_kernel1.data(), _bias1.data(),
				_kernel2.data(), _bias2.data(),
				pred);
	}

	// -------------------------------------------------------------------------
	// Jacobian dispatcher — writes d(ANN)/d(x_normalised) to grad (length _nInput).
	// Input x must already be normalised.
	// -------------------------------------------------------------------------
	void _annJacobian(const double* x, double* grad) const
	{
		if (_nLayers == 1u)
			_annModel->jacobian_single_layer(
				x,
				_kernel0.data(), _bias0.data(),
				_kernel1.data(),
				grad);
		else
			_annModel->jacobian_two_layers(
				x,
				_kernel0.data(), _bias0.data(),
				_kernel1.data(), _bias1.data(),
				_kernel2.data(),
				grad);
	}
};

typedef NeuralNetworkBindingBase<MachineLearningParamHandler>    NeuralNetworkBinding;
typedef NeuralNetworkBindingBase<ExtMachineLearningParamHandler> ExternalNeuralNetworkBinding;

namespace binding
{
	void registerMachineLearningModel(
		std::unordered_map<std::string, std::function<model::IBindingModel*()>>& bindings)
	{
		bindings[NeuralNetworkBinding::identifier()]         = []() { return new NeuralNetworkBinding(); };
		bindings[ExternalNeuralNetworkBinding::identifier()] = []() { return new ExternalNeuralNetworkBinding(); };
	}
}  // namespace binding

}  // namespace model
}  // namespace cadet
