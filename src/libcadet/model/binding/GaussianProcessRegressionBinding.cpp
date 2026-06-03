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

template <class ParamHandler_t>
class GPRBindingBase : public ParamHandlerBindingModelBase<ParamHandler_t>
{
public:

	GPRBindingBase()
		: _poreConc()
		, _solidConc()
		, _gprParams()
		, _kernelMat()
		, _alpha()
		, _offset(0.0)
		, _nDim(1u)
		, _gpModel()
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

	// Training data and GPR state — populated once in configureImpl
	std::vector<double>            _poreConc;    // X_train: (nTrain x nDim), row-major
	std::vector<double>            _solidConc;   // Y_train: (nTrain,)
	std::vector<double>            _gprParams;   // hyperparameters
	std::vector<double>            _kernelMat;   // K(X_train, X_train): (nTrain x nTrain)
	std::vector<double>            _alpha;       // (K + sigma^2 I)^{-1} y: (nTrain,)
	double                         _offset;      // prediction at c_p = 0 for bias removal
	unsigned int                   _nDim;        // GPR input dimensionality (== nComp)
	std::unique_ptr<GP::GPR_Class> _gpModel;
	mutable std::vector<double>    _cpVecCache;  // reusable buffer for c_p vector (length _nDim)

	virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
	{
		const bool result = ParamHandlerBindingModelBase<ParamHandler_t>::configureImpl(
			paramProvider, unitOpIdx, parTypeIdx);

		// --- Read training data and configuration ---
		_solidConc = paramProvider.getDoubleArray("CS_VALS");
		_poreConc  = paramProvider.getDoubleArray("CP_VALS");

		const std::string kernel = paramProvider.getString("KERNEL");

		// Initialize _gprParams with 7 entries (all kernel hyperparameters)
		// [mlp_weight_var, mlp_bias_var, mlp_var, lin_var, rbf_var, rbf_ls, gaussian_var]
		_gprParams.assign(7, 0.0);

		// Read kernel-specific parameters and use dummy values (0.0) for unused ones
		if (kernel == "MLP")
		{
			_gprParams[0] = paramProvider.getDouble("MLP_WEIGHT_VARIANCE");
			_gprParams[1] = paramProvider.getDouble("MLP_BIAS_VARIANCE");
			_gprParams[2] = paramProvider.getDouble("MLP_VARIANCE");
			// linear_var, rbf_var, rbf_ls remain 0.0 (unused)
			_gprParams[6] = paramProvider.getDouble("GAUSSIAN_NOISE_VARIANCE");
		}
		else if (kernel == "RBF")
		{
			// mlp_weight_var, mlp_bias_var, mlp_var remain 0.0 (unused)
			_gprParams[4] = paramProvider.getDouble("RBF_VARIANCE");
			_gprParams[5] = paramProvider.getDouble("RBF_LENGTHSCALE");
			_gprParams[6] = paramProvider.getDouble("GAUSSIAN_NOISE_VARIANCE");
		}
		else if (kernel == "RBF_Linear")
		{
			// mlp_weight_var, mlp_bias_var, mlp_var remain 0.0 (unused)
			_gprParams[3] = paramProvider.getDouble("LINEAR_VARIANCE");
			_gprParams[4] = paramProvider.getDouble("RBF_VARIANCE");
			_gprParams[5] = paramProvider.getDouble("RBF_LENGTHSCALE");
			_gprParams[6] = paramProvider.getDouble("GAUSSIAN_NOISE_VARIANCE");
		}
		else if (kernel == "MLP_Linear")
		{
			_gprParams[0] = paramProvider.getDouble("MLP_WEIGHT_VARIANCE");
			_gprParams[1] = paramProvider.getDouble("MLP_BIAS_VARIANCE");
			_gprParams[2] = paramProvider.getDouble("MLP_VARIANCE");
			_gprParams[3] = paramProvider.getDouble("LINEAR_VARIANCE");
			// rbf_var, rbf_ls remain 0.0 (unused)
			_gprParams[6] = paramProvider.getDouble("GAUSSIAN_NOISE_VARIANCE");
		}
		else
		{
			throw InvalidParameterException(
				"KERNEL must be one of: MLP, RBF, RBF_Linear, MLP_Linear");
		}

		for (int param = 0; param < 7; ++param)
		{
			if (_gprParams[param] < 0.0)
				throw InvalidParameterException(
					"All GPR hyperparameters must be non-negative.");
		}

		const int nDimRaw = paramProvider.getInt("NDIM");
		if (nDimRaw <= 0)
			throw InvalidParameterException(
				"NDIM must be a positive integer specifying the number of GPR input dimensions.");
		_nDim = static_cast<unsigned int>(nDimRaw);

		// --- Validation ---

		// NDIM must match NCOMP: each GPR input dimension corresponds to one chromatographic component.
		// NDIM is retained as an explicit user-supplied sense check against NCOMP.
		if (_nDim != _nComp)
			throw InvalidParameterException(
				"NDIM must equal NCOMP for GPR binding — each GPR input dimension "
				"corresponds to one chromatographic component.");

		if (_gprParams.size() < 7u)
			throw InvalidParameterException(
				"Internal error: _gprParams must contain 7 entries.");

		if (_solidConc.empty() || _poreConc.empty())
			throw InvalidParameterException(
				"Training data CS_VALS and CP_VALS must not be empty.");

		// CP_VALS has shape (nTrain x nDim) flattened row-major; CS_VALS has shape (nTrain,)
		if (_poreConc.size() != _solidConc.size() * _nDim)
			throw InvalidParameterException(
				"CP_VALS size must equal CS_VALS size * NDIM. "
				"CP_VALS is expected as a flat (nTrain x nDim) row-major array.");

		const unsigned int nTrain = static_cast<unsigned int>(_solidConc.size());

		// --- Build GPR model ---
		_gpModel = std::make_unique<GP::GPR_Class>(
			nTrain,
			1u,
			_nDim,
			_gprParams[0],   // mlp_weight_variance
			_gprParams[1],   // mlp_bias_variance
			_gprParams[2],   // mlp_variance
			_gprParams[3],   // linear_variance
			_gprParams[4],   // rbf_variance
			_gprParams[5],   // rbf_lengthscale (= ls^2 as exported from Python)
			_gprParams[6],   // GAUSSIAN_NOISE_VARIANCE (noise)
			kernel);

		// --- Compute training kernel matrix K(X_train, X_train) ---
		_kernelMat.assign(static_cast<std::size_t>(nTrain) * nTrain, 0.0);
		_gpModel->GPR_kernel(_poreConc.data(), _poreConc.data(), _kernelMat.data());

		// --- Solve for alpha = (K + sigma^2 I)^{-1} y ---
		// kernel_inv_y does NOT destroy _kernelMat
		_alpha.assign(nTrain, 0.0);
		_gpModel->kernel_inv_y(_solidConc.data(), _kernelMat.data(), _alpha.data());

		// --- Compute offset = GPR prediction at c_p = 0 (all dimensions zero) ---
		// Enforces q(c_p = 0) = 0. Only physically valid if the isotherm passes through
		// the origin. Verify against training data before enabling.
		_cpVecCache.assign(_nDim, 0.0);
		_offset = _gpModel->prediction(_poreConc.data(), _cpVecCache.data(), _alpha.data());

		return result;
	}

	// -------------------------------------------------------------------------
	// Flux: res[bndIdx] = kKin[i] * (q[bndIdx] - q_GPR(c_p))
	//
	// The GPR isotherm is evaluated at the full c_p vector (all nDim components).
	// At steady state (res = 0): q = q_GPR(c_p).
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

		// Fill c_p vector from all nDim pore concentrations
		for (unsigned int d = 0; d < _nDim; ++d)
			_cpVecCache[d] = static_cast<double>(yCp[d]);

		const double qML = _gpModel->prediction(_poreConc.data(), _cpVecCache.data(), _alpha.data())
			- _offset;

		unsigned int bndIdx = 0;
		for (unsigned int i = 0; i < _nComp; ++i)
		{
			if (_nBoundStates[i] == 0)
				continue;

			res[bndIdx] = static_cast<ParamType>(p->kKin[i])
				* (static_cast<ResidualType>(y[bndIdx]) - static_cast<ResidualType>(qML));

			++bndIdx;
		}

		return 0;
	}

	// -------------------------------------------------------------------------
	// Jacobian of res[bndIdx] = kKin[i] * (q[bndIdx] - q_GPR(c_p))
	//
	// dres[bndIdx]/dq[bndIdx]  =  kKin[i]                     -> jac[0]
	// dres[bndIdx]/dc_p[d]     = -kKin[i] * dq_GPR/dc_p[d]    -> jac[d - bndIdx - offsetCp]
	//
	// The GPR gradient dq_GPR/dc_p is a vector of length nDim (one partial per
	// input dimension / component), computed analytically from the kernel.
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
		std::vector<double> cpVec(_nDim);
		for (unsigned int d = 0; d < _nDim; ++d)
			cpVec[d] = yCp[d];

		// Compute full gradient dq_GPR/dc_p — length nDim
		const std::vector<double> dqdc =
			_gpModel->GPR_derivative(_poreConc.data(), cpVec.data(), _alpha.data());

		unsigned int bndIdx = 0;
		for (int i = 0; i < static_cast<int>(_nComp); ++i)
		{
			if (_nBoundStates[i] == 0)
				continue;

			const double kkin = static_cast<double>(p->kKin[i]);

			// dres[bndIdx] / dc_p[d] for each input dimension d
			for (int d = 0; d < static_cast<int>(_nDim); ++d)
				jac[d - static_cast<int>(bndIdx) - offsetCp] = -kkin * dqdc[d];

			// dres[bndIdx] / dq[bndIdx]
			jac[0] = kkin;

			++bndIdx;
			++jac;
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
