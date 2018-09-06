// =============================================================================
//  CADET - The Chromatography Analysis and Design Toolkit
//  
//  Copyright © 2008-2018: The CADET Authors
//            Please see the AUTHORS and CONTRIBUTORS file.
//  
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
* @file
* Defines the CSTR model as a unit operation.
*/

#ifndef LIBCADET_CSTR_HPP_
#define LIBCADET_CSTR_HPP_

#include "model/UnitOperationBase.hpp"
#include "cadet/SolutionExporter.hpp"
#include "AutoDiff.hpp"
#include "linalg/DenseMatrix.hpp"
#include "model/ModelUtils.hpp"

#include <array>
#include <vector>

namespace cadet
{

namespace model
{

/**
 * @brief Continuously stirred tank (reactor) model
 * @details This is a simple CSTR model with variable volume using the ``well mixed assumption''.
 * @f[\begin{align}
	\frac{\mathrm{d}}{\mathrm{d} t}\left( V \left[ c_i + \frac{1}{\beta} \sum_{j=1}^{N_{\text{bnd},i}} q_{i,j} \right] \right) &= F_{\text{in}} c_{\text{in},i} - F_{\text{out}} c_i \\
	a \frac{\partial q_{i,j}}{\partial t} &= f_{\text{iso},i,j}(c, q) \\
	\frac{\partial V}{\partial t} &= F_{\text{in}} - F_{\text{out}} - F_{\text{filter}}
\end{align} @f]
 * The model can be used as a plain stir tank without any binding states.
 */
class CSTRModel : public UnitOperationBase
{
public:

	CSTRModel(UnitOpIdx unitOpIdx);
	virtual ~CSTRModel() CADET_NOEXCEPT;

	virtual unsigned int numDofs() const CADET_NOEXCEPT;
	virtual unsigned int numPureDofs() const CADET_NOEXCEPT;
	virtual bool usesAD() const CADET_NOEXCEPT;
	virtual unsigned int requiredADdirs() const CADET_NOEXCEPT;

	virtual UnitOpIdx unitOperationId() const CADET_NOEXCEPT { return _unitOpIdx; }
	virtual unsigned int numComponents() const CADET_NOEXCEPT { return _nComp; }
	virtual void setFlowRates(const active& in, const active& out) CADET_NOEXCEPT;
	virtual bool canAccumulate() const CADET_NOEXCEPT { return true; }

	static const char* identifier() { return "CSTR"; }
	virtual const char* unitOperationName() const CADET_NOEXCEPT { return "CSTR"; }

	virtual bool configureModelDiscretization(IParameterProvider& paramProvider, IConfigHelper& helper);
	virtual bool configure(IParameterProvider& paramProvider);
	virtual void notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active* const adRes, active* const adY, unsigned int adDirOffset);

	virtual void useAnalyticJacobian(const bool analyticJac);

	virtual void reportSolution(ISolutionRecorder& recorder, double const* const solution) const;
	virtual void reportSolutionStructure(ISolutionRecorder& recorder) const;

	virtual int residual(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double* const res);
	virtual int residual(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, double* const res, active* const adRes, active* const adY, unsigned int adDirOffset, bool updateJacobian, bool paramSensitivity);
	virtual int residualWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, double* const res, active* const adRes, active* const adY, unsigned int adDirOffset);

	virtual int residualSensFwdAdOnly(const active& t, unsigned int secIdx, const active& timeFactor,
		double const* const y, double const* const yDot, active* const adRes);

	virtual int residualSensFwdWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, active* const adRes, active* const adY, unsigned int adDirOffset);

	virtual int linearSolve(double t, double timeFactor, double alpha, double tol, double* const rhs, double const* const weight,
		double const* const y, double const* const yDot);

	virtual void prepareADvectors(active* const adRes, active* const adY, unsigned int adDirOffset) const;

	virtual void applyInitialCondition(double* const vecStateY, double* const vecStateYdot) const;
	virtual void readInitialCondition(IParameterProvider& paramProvider);

	virtual void consistentInitialState(double t, unsigned int secIdx, double timeFactor, double* const vecStateY, active* const adRes, active* const adY, unsigned int adDirOffset, double errorTol);
	virtual void consistentInitialTimeDerivative(double t, unsigned int secIdx, double timeFactor, double const* vecStateY, double* const vecStateYdot);
	
	virtual void consistentInitialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
		std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes);

	virtual void leanConsistentInitialState(double t, unsigned int secIdx, double timeFactor, double* const vecStateY, active* const adRes, active* const adY, unsigned int adDirOffset, double errorTol);
	virtual void leanConsistentInitialTimeDerivative(double t, double timeFactor, double const* const vecStateY, double* const vecStateYdot, double* const res);
	
	virtual void leanConsistentInitialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
		std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes);

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { }

	virtual void multiplyWithJacobian(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double const* yS, double alpha, double beta, double* ret);
	virtual void multiplyWithDerivativeJacobian(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double const* sDot, double* ret);

	virtual bool hasInlet() const CADET_NOEXCEPT { return true; }
	virtual bool hasOutlet() const CADET_NOEXCEPT { return true; }

	virtual unsigned int localOutletComponentIndex() const CADET_NOEXCEPT { return _nComp; }
	virtual unsigned int localOutletComponentStride() const CADET_NOEXCEPT { return 1; }
	virtual unsigned int localInletComponentIndex() const CADET_NOEXCEPT { return 0; }
	virtual unsigned int localInletComponentStride() const CADET_NOEXCEPT { return 1; }

	virtual void setSectionTimes(double const* secTimes, bool const* secContinuity, unsigned int nSections);

	virtual void expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut) { }

#ifdef CADET_BENCHMARK_MODE
	virtual std::vector<double> benchmarkTimings() const { return std::vector<double>(0); }
	virtual char const* const* benchmarkDescriptions() const { return nullptr; }
#endif

	inline const std::vector<active>& flowRateFilter() const { return _flowRateFilter; }
	inline std::vector<active>& flowRateFilter() { return _flowRateFilter; }
	inline void flowRateFilter(const std::vector<active>& frf) { _flowRateFilter = frf; }

protected:

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	int residualImpl(const ParamType& t, unsigned int secIdx, const ParamType& timeFactor, StateType const* const y, double const* const yDot, ResidualType* const res);

	template <typename MatrixType>
	void addTimeDerivativeJacobian(double t, double timeFactor, double const* y, double const* yDot, MatrixType& mat);

	void extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset);
#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	void checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const;
#endif

	unsigned int _nComp; //!< Number of components
	unsigned int* _nBound; //!< Array with number of bound states for each component
	unsigned int* _boundOffset; //!< Array with offset to the first bound state of each component in the solid phase
	unsigned int _strideBound; //!< Total number of bound states

	active _porosity; //!< Porosity \f$ \varepsilon \f$
	active _flowRateIn; //!< Volumetric flow rate of incoming stream
	active _flowRateOut; //!< Volumetric flow rate of drawn outgoing stream
	active _curFlowRateFilter; //!< Current volumetric flow rate of liquid outtake stream for this section
	std::vector<active> _flowRateFilter; //!< Volumetric flow rate of liquid outtake stream

	bool _analyticJac; //!< Flag that determines whether analytic or AD Jacobian is used
	linalg::DenseMatrix _jac; //!< Jacobian
	linalg::DenseMatrix _jacFact; //!< Factorized Jacobian
	bool _factorizeJac; //!< Flag that tracks whether the Jacobian needs to be factorized
	double* _consistentInitBuffer; //!< Memory for consistent initialization (solving nonlinear equations)

	std::vector<active> _initConditions; //!< Initial conditions, ordering: Liquid phase concentration, solid phase concentration, volume
	std::vector<double> _initConditionsDot; //!< Initial conditions for time derivative

	class Exporter : public ISolutionExporter
	{
	public:

		Exporter(unsigned int nComp, unsigned int const* nBound, unsigned int strideBound, unsigned int const* boundOffset, double const* data) : _data(data), _nComp(nComp), _nBound(nBound), _strideBound(strideBound), _boundOffset(boundOffset) { }

		virtual bool hasMultipleBoundStates() const CADET_NOEXCEPT { return cadet::model::hasMultipleBoundStates(_nBound, _nComp); }
		virtual bool hasNonBindingComponents() const CADET_NOEXCEPT { return cadet::model::hasNonBindingComponents(_nBound, _nComp); }
		virtual bool hasParticleFlux() const CADET_NOEXCEPT { return false; }
		virtual bool hasParticleMobilePhase() const CADET_NOEXCEPT { return false; }
		virtual bool hasSolidPhase() const CADET_NOEXCEPT { return _strideBound > 0; }
		virtual bool hasVolume() const CADET_NOEXCEPT { return true; }

		virtual unsigned int numComponents() const CADET_NOEXCEPT { return _nComp; }
		virtual unsigned int numAxialCells() const CADET_NOEXCEPT { return 1; }
		virtual unsigned int numRadialCells() const CADET_NOEXCEPT { return 0; }
		virtual unsigned int numBoundStates() const CADET_NOEXCEPT { return _strideBound; }
		virtual unsigned int const* numBoundStatesPerComponent() const CADET_NOEXCEPT { return _nBound; }
		virtual unsigned int numBoundStates(unsigned int comp) const CADET_NOEXCEPT { return _nBound[comp]; }
		virtual unsigned int numBulkDofs() const CADET_NOEXCEPT { return _nComp; }
		virtual unsigned int numParticleMobilePhaseDofs() const CADET_NOEXCEPT { return 0; }
		virtual unsigned int numSolidPhaseDofs() const CADET_NOEXCEPT { return _strideBound; }
		virtual unsigned int numFluxDofs() const CADET_NOEXCEPT { return 0; }
		virtual unsigned int numVolumeDofs() const CADET_NOEXCEPT { return 1; }

		virtual double concentration(unsigned int component, unsigned int axialCell) const { return _data[_nComp + component]; }
		virtual double flux(unsigned int component, unsigned int axialCell) const { return 0.0; }
		virtual double mobilePhase(unsigned int component, unsigned int axialCell, unsigned int radialCell) const { return 0.0; }
		virtual double solidPhase(unsigned int component, unsigned int axialCell, unsigned int radialCell, unsigned int boundState) const { return _data[2 * _nComp + _boundOffset[component] + boundState]; }
		virtual double volume(unsigned int dof) const { return _data[2 * _nComp + dof]; }

		virtual double const* concentration() const { return _data + _nComp; }
		virtual double const* flux() const { return nullptr; }
		virtual double const* mobilePhase() const { return nullptr; }
		virtual double const* solidPhase() const { return _data + 2 * _nComp; }
		virtual double const* volume() const { return _data + 2 * _nComp + _strideBound; }
		virtual double const* inlet(unsigned int& stride) const
		{
			stride = 1;
			return _data;
		}
		virtual double const* outlet(unsigned int& stride) const
		{
			stride = 1;
			return _data + _nComp;
		}

		virtual StateOrdering const* concentrationOrdering(unsigned int& len) const
		{
			len = _concentrationOrdering.size();
			return _concentrationOrdering.data();
		}

		virtual StateOrdering const* fluxOrdering(unsigned int& len) const
		{
			len = 0;
			return nullptr;
		}

		virtual StateOrdering const* mobilePhaseOrdering(unsigned int& len) const
		{
			len = 0;
			return nullptr;
		}

		virtual StateOrdering const* solidPhaseOrdering(unsigned int& len) const
		{
			len = _solidOrdering.size();
			return _solidOrdering.data();
		}

		virtual unsigned int bulkMobilePhaseStride() const { return _nComp; }
		virtual unsigned int particleMobilePhaseStride() const { return 0; }
		virtual unsigned int solidPhaseStride() const { return _strideBound; }

	protected:
		double const* const _data;
		unsigned int _nComp;
		unsigned int const* _nBound;
		unsigned int _strideBound;
		unsigned int const* _boundOffset;

		const std::array<StateOrdering, 1> _concentrationOrdering = { { StateOrdering::Component } };
		const std::array<StateOrdering, 2> _solidOrdering = { { StateOrdering::Component, StateOrdering::BoundState } };
	};
};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_CSTR_HPP_
