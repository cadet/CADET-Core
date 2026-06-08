// =============================================================================
//  CADET
//
//  Copyright © 2008-present: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

#include "model/binding/BindingModelBase.hpp"
#include "model/ExternalFunctionSupport.hpp"
#include "model/binding/BindingModelMacros.hpp"
#include "model/ModelUtils.hpp"
#include "cadet/Exceptions.hpp"
#include "model/Parameters.hpp"
#include "LocalVector.hpp"
#include "SimulationTypes.hpp"
#include "AdUtils.hpp"
#include "model/binding/GPR_Class.h"

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <iomanip>
#include <sstream>

/*<codegen>
{
	"name": "GPRParamHandler",
	"externalName": "ExtGPRParamHandler",
	"parameters":
		[
			{ "type": "ScalarComponentDependentParameter", "varName": "kKin", "confName": "GPR_KKIN"}
		]
}
</codegen>*/

namespace cadet
{

namespace model
{

inline const char* GPRParamHandler::identifier() CADET_NOEXCEPT { return "GAUSSIAN_PROCESS_REGRESSION"; }

inline bool GPRParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if (_kKin.size() < nComp)
		throw InvalidParameterException("GPR_KKIN has to have NTOTALBOUND entries");

	return true;
}

inline const char* ExtGPRParamHandler::identifier() CADET_NOEXCEPT { return "EXT_GAUSSIAN_PROCESS_REGRESSION"; }

inline bool ExtGPRParamHandler::validateConfig(unsigned int nComp, unsigned int const* nBoundStates)
{
	if (_kKin.size() < nComp)
		throw InvalidParameterException("GPR_KKIN has to have NTOTALBOUND entries");

	return true;
}

// Helper struct to hold bound-state-specific GPR hyperparameters
struct BoundStateGPRParams
{
	double mlp_weight_variance;
	double mlp_bias_variance;
	double mlp_variance;
	double linear_variance;
	double rbf_variance;
	double rbf_lengthscale;
	double gaussian_noise_variance;
	std::string kernel;
};

template <class ParamHandler_t>
class GPRBindingBase : public ParamHandlerBindingModelBase<ParamHandler_t>
{
public:

	GPRBindingBase()
		: _poreConc()
		, _solidConc()
		, _boundStateParams()
		, _gpModels()
		, _kernelMats()
		, _alphas()
		, _offsets()
		, _nTrain(0u)
		, _cpVecCache()
	{
	}

	virtual ~GPRBindingBase() CADET_NOEXCEPT { }

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

	// Training data and GPR state — one model per bound state
	std::vector<double>                             _poreConc;           // X_train: (nTrain x nComp), row-major, shared across all bound states
	std::vector<double>                             _solidConc;          // Y_train: (nTrain x totalNumBoundStates), row-major, bound-state-specific outputs
	std::vector<BoundStateGPRParams>                _boundStateParams;   // hyperparameters per bound state
	std::vector<std::unique_ptr<GP::GPR_Class>>     _gpModels;           // one GPR model per bound state
	std::vector<std::vector<double>>                _kernelMats;         // one kernel matrix per bound state: each (nTrain x nTrain)
	std::vector<std::vector<double>>                _alphas;             // one alpha vector per bound state: each (nTrain,)
	std::vector<double>                             _offsets;            // one offset per bound state (prediction at c_p = 0)
	unsigned int                                    _nTrain;             // number of training points
	mutable std::vector<double>                     _cpVecCache;         // reusable buffer for c_p vector (length _nComp)

	// Helper function to create bound-state-specific parameter name with BNDXXX suffix
	std::string makeBoundStateParamName(const std::string& baseName, unsigned int bndIdx) const
	{
		if (bndIdx == static_cast<unsigned int>(-1))
			return baseName;
		std::ostringstream oss;
		oss << baseName << "_BND_" << std::setw(3) << std::setfill('0') << bndIdx;
		return oss.str();
	}

	// Helper function to read bound-state-specific GPR kernel parameters
	BoundStateGPRParams readBoundStateParams(IParameterProvider& paramProvider, unsigned int bndIdx) const
	{
		BoundStateGPRParams params;

		// Read kernel type for this bound state
		params.kernel = paramProvider.getString("KERNEL");

		// Initialize all parameters to 0.0
		params.mlp_weight_variance = 0.0;
		params.mlp_bias_variance = 0.0;
		params.mlp_variance = 0.0;
		params.linear_variance = 0.0;
		params.rbf_variance = 0.0;
		params.rbf_lengthscale = 0.0;
		params.gaussian_noise_variance = 0.0;

		// Read kernel-specific parameters
		if (params.kernel == "MLP")
		{
			if (_nBoundStates[_nComp - 1] == 1 && paramProvider.exists("MLP_WEIGHT_VARIANCE"))
				bndIdx = -1; // no suffix
			params.mlp_weight_variance = paramProvider.getDouble(
				makeBoundStateParamName("MLP_WEIGHT_VARIANCE", bndIdx));
			params.mlp_bias_variance = paramProvider.getDouble(
				makeBoundStateParamName("MLP_BIAS_VARIANCE", bndIdx));
			params.mlp_variance = paramProvider.getDouble(
				makeBoundStateParamName("MLP_VARIANCE", bndIdx));
			params.gaussian_noise_variance = paramProvider.getDouble(
				makeBoundStateParamName("GAUSSIAN_NOISE_VARIANCE", bndIdx));
		}
		else if (params.kernel == "MLP_Linear")
		{
			if (_nBoundStates[_nComp - 1] == 1 && paramProvider.exists("MLP_WEIGHT_VARIANCE"))
				bndIdx = -1; // no suffix
			params.mlp_weight_variance = paramProvider.getDouble(
				makeBoundStateParamName("MLP_WEIGHT_VARIANCE", bndIdx));
			params.mlp_bias_variance = paramProvider.getDouble(
				makeBoundStateParamName("MLP_BIAS_VARIANCE", bndIdx));
			params.mlp_variance = paramProvider.getDouble(
				makeBoundStateParamName("MLP_VARIANCE", bndIdx));
			params.linear_variance = paramProvider.getDouble(
				makeBoundStateParamName("LINEAR_VARIANCE", bndIdx));
			params.gaussian_noise_variance = paramProvider.getDouble(
				makeBoundStateParamName("GAUSSIAN_NOISE_VARIANCE", bndIdx));
		}
		else if (params.kernel == "RBF")
		{
			if (_nBoundStates[_nComp - 1] == 1 && paramProvider.exists("RBF_VARIANCE"))
				bndIdx = -1; // no suffix
			params.rbf_variance = paramProvider.getDouble(
				makeBoundStateParamName("RBF_VARIANCE", bndIdx));
			params.rbf_lengthscale = paramProvider.getDouble(
				makeBoundStateParamName("RBF_LENGTHSCALE", bndIdx));
			params.gaussian_noise_variance = paramProvider.getDouble(
				makeBoundStateParamName("GAUSSIAN_NOISE_VARIANCE", bndIdx));
		}
		else if (params.kernel == "RBF_Linear")
		{
			if (_nBoundStates[_nComp - 1] == 1 && paramProvider.exists("RBF_VARIANCE"))
				bndIdx = -1; // no suffix
			params.linear_variance = paramProvider.getDouble(
				makeBoundStateParamName("LINEAR_VARIANCE", bndIdx));
			params.rbf_variance = paramProvider.getDouble(
				makeBoundStateParamName("RBF_VARIANCE", bndIdx));
			params.rbf_lengthscale = paramProvider.getDouble(
				makeBoundStateParamName("RBF_LENGTHSCALE", bndIdx));
			params.gaussian_noise_variance = paramProvider.getDouble(
				makeBoundStateParamName("GAUSSIAN_NOISE_VARIANCE", bndIdx));
		}
		else
		{
			throw InvalidParameterException(
				"KERNEL must be one of: MLP, RBF, RBF_Linear, MLP_Linear");
		}

		// Validate hyperparameters
		if (params.mlp_weight_variance < 0.0 || params.mlp_bias_variance < 0.0 ||
			params.mlp_variance < 0.0 || params.linear_variance < 0.0 ||
			params.rbf_variance < 0.0 || params.rbf_lengthscale < 0.0 ||
			params.gaussian_noise_variance < 0.0)
		{
			throw InvalidParameterException(
				"All GPR hyperparameters for bound state " + std::to_string(bndIdx) + " must be non-negative.");
		}

		return params;
	}

	virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
	{
		const bool result = ParamHandlerBindingModelBase<ParamHandler_t>::configureImpl(
			paramProvider, unitOpIdx, parTypeIdx);

		// --- Read training data and configuration ---
		_solidConc = paramProvider.getDoubleArray("CS_VALS");
		_poreConc  = paramProvider.getDoubleArray("CP_VALS");

		const int cpNDim = paramProvider.getInt("CP_NDIM");
		const int csNDim = paramProvider.getInt("CS_NDIM");

		// --- Validation ---

		// CP_NDIM must match NCOMP: each GPR input dimension corresponds to one chromatographic component.
		if (static_cast<unsigned int>(cpNDim) != _nComp)
			throw InvalidParameterException(
				"CP_NDIM must equal NCOMP for GPR binding — each GPR input dimension "
				"corresponds to one chromatographic component.");

		// Calculate total number of bound states
		unsigned int totalNumBoundStates = 0;
		for (unsigned int i = 0; i < _nComp; ++i)
		{
			if(_nBoundStates[i] != 1)
				throw InvalidParameterException(
					"GPR binding requires exactly one bound state for each component");
			totalNumBoundStates += _nBoundStates[i];
		}

		// CS_NDIM must match total number of bound states
		if (static_cast<unsigned int>(csNDim) != totalNumBoundStates)
			throw InvalidParameterException(
				"CS_NDIM must equal total number of bound states (sum of all NBOUND entries) = "
				+ std::to_string(totalNumBoundStates) + ", but got " + std::to_string(csNDim));

		if (totalNumBoundStates == 0)
			throw InvalidParameterException(
				"GPR binding requires at least one bound state.");

		if (_solidConc.empty() || _poreConc.empty())
			throw InvalidParameterException(
				"Training data CS_VALS and CP_VALS must not be empty.");

		// CP_VALS has shape (nTrain x nComp) flattened row-major
		// CS_VALS has shape (nTrain x totalNumBoundStates) flattened row-major — one output column per bound state
		if (_poreConc.size() % _nComp != 0)
			throw InvalidParameterException(
				"CP_VALS size must be a multiple of CP_NDIM.");

		_nTrain = static_cast<unsigned int>(_poreConc.size() / _nComp);

		if (_solidConc.size() != static_cast<std::size_t>(_nTrain) * totalNumBoundStates)
			throw InvalidParameterException(
				"CP_VALS size must equal NTRAIN * NCOMP = " + std::to_string(_nTrain * _nComp) + ". "
				"CS_VALS size must equal NTRAIN * NTOTALBOUND = " + std::to_string(_nTrain * totalNumBoundStates) + ". "
				"CS_VALS is expected as a flat (nTrain x totalNumBoundStates) row-major array, "
				"with one column per bound state.");

		// --- Read bound-state-specific hyperparameters and build GPR models ---
		_boundStateParams.clear();
		_gpModels.clear();
		_kernelMats.clear();
		_alphas.clear();
		_offsets.clear();

		_boundStateParams.reserve(totalNumBoundStates);
		_gpModels.reserve(totalNumBoundStates);
		_kernelMats.reserve(totalNumBoundStates);
		_alphas.reserve(totalNumBoundStates);
		_offsets.reserve(totalNumBoundStates);

		// Loop through all bound states across all components
		unsigned int bndIdx = 0;
		for (unsigned int comp = 0; comp < _nComp; ++comp)
		{
			for (unsigned int bnd = 0; bnd < _nBoundStates[comp]; ++bnd)
			{
				// Read bound-state-specific hyperparameters
				BoundStateGPRParams params = readBoundStateParams(paramProvider, bndIdx);
				_boundStateParams.push_back(params);

				// Extract Y_train for this bound state: column bndIdx of CS_VALS
				std::vector<double> yTrain(_nTrain);
				for (unsigned int row = 0; row < _nTrain; ++row)
					yTrain[row] = _solidConc[row * totalNumBoundStates + bndIdx];

				// Create GPR model for this bound state with its specific hyperparameters
				auto gpModel = std::make_unique<GP::GPR_Class>(
					_nTrain,
					1u,
					_nComp,
					params.mlp_weight_variance,
					params.mlp_bias_variance,
					params.mlp_variance,
					params.linear_variance,
					params.rbf_variance,
					params.rbf_lengthscale,
					params.gaussian_noise_variance,
					params.kernel);

				// Compute training kernel matrix K(X_train, X_train) for this bound state
				std::vector<double> kernelMat(static_cast<std::size_t>(_nTrain) * _nTrain, 0.0);
				gpModel->GPR_kernel(_poreConc.data(), _poreConc.data(), kernelMat.data());

				// Solve for alpha = (K + sigma^2 I)^{-1} y
				std::vector<double> alpha(_nTrain, 0.0);
				gpModel->kernel_inv_y(yTrain.data(), kernelMat.data(), alpha.data());

				// Compute offset = GPR prediction at c_p = 0 (all dimensions zero)
				std::vector<double> cpZero(_nComp, 0.0);
				const double offset = gpModel->prediction(_poreConc.data(), cpZero.data(), alpha.data());

				// Store model, kernel matrix, alpha, and offset for this bound state
				_gpModels.push_back(std::move(gpModel));
				_kernelMats.push_back(std::move(kernelMat));
				_alphas.push_back(std::move(alpha));
				_offsets.push_back(offset);

				++bndIdx;
			}
		}

		// Allocate cache for c_p vector
		_cpVecCache.assign(_nComp, 0.0);

		return result;
	}

	// -------------------------------------------------------------------------
	// Flux: res[bndIdx] = kKin[comp] * (q[bndIdx] - q_GPR_bndIdx(c_p))
	//
	// Each bound state uses its own GPR model trained on its corresponding
	// column in CS_VALS with its own hyperparameters.
	//
	// At steady state (res = 0): q[bndIdx] = q_GPR_bndIdx(c_p).
	//
	// Note: GPR prediction always operates in double. AD types in yCp are cast
	// to double here. The Jacobian w.r.t. yCp is provided analytically in jacobianImpl.
	// -------------------------------------------------------------------------
	template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
	int fluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
		StateType const* y, CpStateType const* yCp, ResidualType* res,
		LinearBufferAllocator workSpace) const
	{
		typename ParamHandler_t::ParamsHandle const p =
			_paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Fill c_p vector from all nComp pore concentrations
		for (unsigned int d = 0; d < _nComp; ++d)
			_cpVecCache[d] = static_cast<double>(yCp[d]);

		unsigned int bndIdx = 0;
		for (unsigned int comp = 0; comp < _nComp; ++comp)
		{
			for (unsigned int bnd = 0; bnd < _nBoundStates[comp]; ++bnd)
			{
				// Evaluate GPR model for this specific bound state (with its own hyperparameters)
				const double qML = _gpModels[bndIdx]->prediction(
					_poreConc.data(), _cpVecCache.data(), _alphas[bndIdx].data())
					- _offsets[bndIdx];

				res[bndIdx] = static_cast<ParamType>(p->kKin[comp])
					* (static_cast<ResidualType>(y[bndIdx]) - static_cast<ResidualType>(qML));

				++bndIdx;
			}
		}

		return 0;
	}

	// -------------------------------------------------------------------------
	// Jacobian of res[bndIdx] = kKin[comp] * (q[bndIdx] - q_GPR_bndIdx(c_p))
	//
	// dres[bndIdx]/dq[bndIdx]  =  kKin[comp]                         -> jac[0]
	// dres[bndIdx]/dc_p[d]     = -kKin[comp] * dq_GPR_bndIdx/dc_p[d] -> jac[d - bndIdx - offsetCp]
	//
	// Each bound state uses its own GPR model (with its own hyperparameters) and
	// thus has its own gradient.
	//
	// Jacobian index convention (follows CADET banded layout, same as Langmuir):
	//   Starting at jac (= row for q[bndIdx]):
	//     -bndIdx          -> q[0]
	//     -bndIdx-offsetCp -> c_p[0]
	//     d - bndIdx - offsetCp -> c_p[d]
	// -------------------------------------------------------------------------
	template <typename RowIterator>
	void jacobianImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
		double const* y, double const* yCp, int offsetCp, RowIterator jac,
		LinearBufferAllocator workSpace) const
	{
		typename ParamHandler_t::ParamsHandle const p =
			_paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Build c_p vector
		std::vector<double> cpVec(_nComp);
		for (unsigned int d = 0; d < _nComp; ++d)
			cpVec[d] = yCp[d];

		unsigned int bndIdx = 0;
		for (unsigned int comp = 0; comp < _nComp; ++comp)
		{
			for (unsigned int bnd = 0; bnd < _nBoundStates[comp]; ++bnd)
			{
				// Compute full gradient dq_GPR_bndIdx/dc_p for this bound state's model
				// (using its own hyperparameters)
				const std::vector<double> dqdc = _gpModels[bndIdx]->GPR_derivative(
					_poreConc.data(), cpVec.data(), _alphas[bndIdx].data());

				const double kkin = static_cast<double>(p->kKin[comp]);

				// dres[bndIdx] / dc_p[d] for each input dimension d
				for (int d = 0; d < static_cast<int>(_nComp); ++d)
					jac[d - static_cast<int>(bndIdx) - offsetCp] = -kkin * dqdc[d];

				// dres[bndIdx] / dq[bndIdx]
				jac[0] = kkin;

				++bndIdx;
				++jac;
			}
		}
	}
};

typedef GPRBindingBase<GPRParamHandler>    GPRBinding;
typedef GPRBindingBase<ExtGPRParamHandler> ExternalGPRBinding;

namespace binding
{
	void registerGaussianProcessRegressionModel(
		std::unordered_map<std::string, std::function<model::IBindingModel*()>>& bindings)
	{
		bindings[GPRBinding::identifier()]         = []() { return new GPRBinding(); };
		bindings[ExternalGPRBinding::identifier()] = []() { return new ExternalGPRBinding(); };
	}
}  // namespace binding

}  // namespace model

}  // namespace cadet
