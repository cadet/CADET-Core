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
		, _kernelName()
		, _nDim(1u)
		, _gpModel()
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
	std::vector<double>            _poreConc;   // X_train: (nTrain x nDim), row-major
	std::vector<double>            _solidConc;  // Y_train: (nTrain,)
	std::vector<double>            _gprParams;  // hyperparameters
	std::vector<double>            _kernelMat;  // K(X_train, X_train): (nTrain x nTrain)
	std::vector<double>            _alpha;      // K^{-1} y: (nTrain,)
	double                         _offset;     // prediction at c_p = 0 for bias removal
	std::string                    _kernelName;
	unsigned int                   _nDim;
	std::unique_ptr<GP::GPR_Class> _gpModel;

	virtual bool configureImpl(IParameterProvider& paramProvider, UnitOpIdx unitOpIdx, ParticleTypeIdx parTypeIdx)
	{
		const bool result = ParamHandlerBindingModelBase<ParamHandler_t>::configureImpl(paramProvider, unitOpIdx, parTypeIdx);

		// Input parameters
		_solidConc = paramProvider.getDoubleArray("CS_VALS");
		_poreConc = paramProvider.getDoubleArray("CP_VALS");
		_gprParams = paramProvider.getDoubleArray("TRAINED_PARAMS");
		_kernelName = paramProvider.getString("KERNEL");
		if (paramProvider.getInt("NDIM") > 0)
			_nDim = paramProvider.getInt("NDIM");
		else
			throw InvalidParameterException("NDIM must be a positive integer specifying the number of dimensions for the GPR input data");

		// --- Validation ---
		if (_nComp != 1u)
			throw InvalidParameterException(
				"Current GPR binding supports single-component only.");

		if (_gprParams.size() < 7u)
			throw InvalidParameterException(
				"TRAINED_PARAMS must contain at least 7 entries: "
				"[mlp_weight_var, mlp_bias_var, mlp_var, lin_var, rbf_var, rbf_ls, gaussian_var]");

		if (_solidConc.empty() || _poreConc.empty())
			throw InvalidParameterException(
				"Training data Q_VALS and C_VALS must not be empty.");

		// C_VALS has shape (nTrain x nDim) flattened, Q_VALS has shape (nTrain,)
		if (_poreConc.size() != _solidConc.size() * _nDim)
			throw InvalidParameterException(
				"C_VALS size must equal Q_VALS size * NDIM.");

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
			_gprParams[5],   // rbf_lengthscale (= ls^2)
			_gprParams[6],   // gaussian_variance (noise)
			_kernelName);

		// --- Compute training kernel matrix K(X_train, X_train) ---
		_kernelMat.assign(static_cast<std::size_t>(nTrain) * nTrain, 0.0);

		if (_kernelName == "MLP")
			_gpModel->MLP_kernel(nTrain, nTrain, _nDim,
				_poreConc.data(), _poreConc.data(), _kernelMat.data());
		else if (_kernelName == "RBF")
			_gpModel->RBF_kernel(_poreConc.data(), _poreConc.data(),
				_kernelMat.data(), nTrain);
		else if (_kernelName == "RBF_Linear")
			_gpModel->RBF_Linear_Kernel(_poreConc.data(), _poreConc.data(),
				_kernelMat.data(), nTrain);
		else if (_kernelName == "MLP_Linear")
			_gpModel->MLP_Linear_Kernel(_poreConc.data(), _poreConc.data(),
				_kernelMat.data(), nTrain);
		else
			throw InvalidParameterException("GPR: unsupported kernel: " + _kernelName);

		// --- Solve for alpha = (K + sigma^2 I)^{-1} y ---
		// kernel_inv_y does NOT destroy _kernelMat (unlike original MKL dposv)
		_alpha.assign(nTrain, 0.0);
		_gpModel->kernel_inv_y(_solidConc.data(), _kernelMat.data(), _alpha.data());

		// --- Compute offset = prediction at c_p = 0 for bias correction ---
		// This enforces q(c_p=0) = 0. Only physically meaningful if your
		// isotherm truly passes through the origin. Verify against training data.
		std::vector<double> zeroPt(_nDim, 0.0);
		_offset = _gpModel->prediction(_poreConc.data(), zeroPt.data(), _alpha.data());

		return result;
	}

	// -------------------------------------------------------------------------
	// Flux: res[bndIdx] = kKin[i] * (q[bndIdx] - q_GPR(c_p))
	//
	// At steady state res = 0, so q = q_GPR(c_p), i.e. the GPR isotherm.
	// The kinetic rate kKin drives q toward the GPR equilibrium prediction.
	// -------------------------------------------------------------------------
	template <typename StateType, typename CpStateType, typename ResidualType, typename ParamType>
	int fluxImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
		StateType const* y, CpStateType const* yCp, ResidualType* res,
		LinearBufferAllocator workSpace) const
	{
		typename ParamHandler_t::ParamsHandle const p =
			_paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Evaluate GPR isotherm at current pore concentration
		// Note: GPR always works in double; AD types in yCp are cast to double here.
		// The Jacobian w.r.t. yCp is provided analytically in jacobianImpl.
		std::vector<double> cpVec(_nDim);
		for (unsigned int i = 0; i < _nDim; ++i)
			cpVec[i] = static_cast<double>(yCp[i]);

		const double qML = _gpModel->prediction(_poreConc.data(), cpVec.data(), _alpha.data())
			- _offset;

		unsigned int bndIdx = 0;
		for (unsigned int i = 0; i < _nComp; ++i)
		{
			if (_nBoundStates[i] == 0)
				continue;

			// res = kKin * (q - q_GPR)
			// At equilibrium (res=0): q = q_GPR(c_p)  ✓
			res[bndIdx] = static_cast<ParamType>(p->kKin[i])
				* (static_cast<ResidualType>(y[bndIdx]) - static_cast<ResidualType>(qML));

			++bndIdx;
		}

		return 0;
	}

	// -------------------------------------------------------------------------
	// Jacobian of res[bndIdx] = kKin[i] * (q[bndIdx] - q_GPR(c_p[i]))
	//
	// dres[bndIdx]/dq[bndIdx]  =  kKin[i]          -> jac[0]
	// dres[bndIdx]/dc_p[i]     = -kKin[i] * dqGPR/dc_p[i]  -> jac[i - bndIdx - offsetCp]
	//
	// No cross-component Jacobian terms (single-component, no coupling).
	// -------------------------------------------------------------------------
	template <typename RowIterator>
	void jacobianImpl(double t, unsigned int secIdx, const ColumnPosition& colPos,
		double const* y, double const* yCp, int offsetCp, RowIterator jac,
		LinearBufferAllocator workSpace) const
	{
		typename ParamHandler_t::ParamsHandle const p =
			_paramHandler.update(t, secIdx, colPos, _nComp, _nBoundStates, workSpace);

		// Build c_p vector for derivative evaluation
		std::vector<double> cpVec(_nDim);
		for (unsigned int i = 0; i < _nDim; ++i)
			cpVec[i] = yCp[i];

		// Compute dq_GPR / dc_p using the analytical kernel derivative
		double dqdc = 0.0;
		if (_kernelName == "MLP")
			dqdc = _gpModel->MLP_derivative(_poreConc.data(), cpVec.data(), _alpha.data());
		else if (_kernelName == "RBF")
			dqdc = _gpModel->RBF_derivative(_poreConc.data(), cpVec.data(), _alpha.data());
		else if (_kernelName == "RBF_Linear")
			dqdc = _gpModel->RBF_Linear_derivative(_poreConc.data(), cpVec.data(), _alpha.data());
		else if (_kernelName == "MLP_Linear")
			dqdc = _gpModel->MLP_Linear_derivative(_poreConc.data(), cpVec.data(), _alpha.data());
		// If kernel is unknown we leave dqdc = 0 — configureImpl would have thrown already.

		unsigned int bndIdx = 0;
		for (int i = 0; i < _nComp; ++i)
		{
			if (_nBoundStates[i] == 0)
				continue;

			const double kkin = static_cast<double>(p->kKin[i]);

			// dres[bndIdx] / dc_p[i]
			// Index: from q[bndIdx], step back bndIdx positions to q[0],
			// then -offsetCp to c_p[0], then +i to c_p[i].
			jac[i - static_cast<int>(bndIdx) - offsetCp] = -kkin * dqdc;

			// dres[bndIdx] / dq[bndIdx]
			jac[0] = kkin;

			++bndIdx;
			++jac;
		}
	}
};

typedef GPRBindingBase<GPRParamHandler> GPRBinding;
typedef GPRBindingBase<ExtGPRParamHandler> ExternalGPRBinding;

namespace binding
{
	void registerGaussianProcessRegressionModel(std::unordered_map<std::string, std::function<model::IBindingModel* ()>>& bindings)
	{
		bindings[GPRBinding::identifier()] = []() { return new GPRBinding(); };
		bindings[ExternalGPRBinding::identifier()] = []() { return new ExternalGPRBinding(); };
	}
}  // namespace binding

}  // namespace model

}  // namespace cadet
