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
 * Defines the outlet model.
 */

#ifndef LIBCADET_OUTLETMODEL_HPP_
#define LIBCADET_OUTLETMODEL_HPP_

#include "cadet/SolutionExporter.hpp"

#include "UnitOperation.hpp"
#include "AutoDiff.hpp"
#include "ParamIdUtil.hpp"
#include "model/ModelUtils.hpp"

#include <array>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace cadet
{

namespace model
{

/**
 * @brief Outlet model
 * @details This unit operation model is solely for buffering concentration profiles
 *          for readout.
 */
class OutletModel : public IUnitOperation
{
public:

	OutletModel(UnitOpIdx unitOpIdx);
	virtual ~OutletModel() CADET_NOEXCEPT;

	virtual unsigned int numDofs() const CADET_NOEXCEPT;
	virtual unsigned int numPureDofs() const CADET_NOEXCEPT;
	virtual bool usesAD() const CADET_NOEXCEPT;
	virtual unsigned int requiredADdirs() const CADET_NOEXCEPT;

	virtual UnitOpIdx unitOperationId() const CADET_NOEXCEPT { return _unitOpIdx; }
	virtual unsigned int numComponents() const CADET_NOEXCEPT { return _nComp; }
	virtual void setFlowRates(const active& in, const active& out) CADET_NOEXCEPT { }
	virtual bool canAccumulate() const CADET_NOEXCEPT { return true; }

	static const char* identifier() { return "OUTLET"; }
	virtual const char* unitOperationName() const CADET_NOEXCEPT { return "OUTLET"; }

	virtual bool configureModelDiscretization(IParameterProvider& paramProvider, IConfigHelper& helper);
	virtual bool configure(IParameterProvider& paramProvider);
	virtual void notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active* const adRes, active* const adY, unsigned int adDirOffset);
	
	virtual std::unordered_map<ParameterId, double> getAllParameterValues() const;
	virtual bool hasParameter(const ParameterId& pId) const;
	virtual double getParameterDouble(const ParameterId& pId) const;

	virtual bool setParameter(const ParameterId& pId, int value);
	virtual bool setParameter(const ParameterId& pId, double value);
	virtual bool setParameter(const ParameterId& pId, bool value);

	virtual bool setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue);
	virtual void setSensitiveParameterValue(const ParameterId& id, double value);

	virtual void clearSensParams();

	virtual void useAnalyticJacobian(const bool analyticJac);

	virtual void reportSolution(ISolutionRecorder& recorder, double const* const solution) const;
	virtual void reportSolutionStructure(ISolutionRecorder& recorder) const;

	virtual int residual(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double* const res);
	virtual int residualWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, double* const res, active* const adRes, active* const adY, unsigned int adDirOffset);

	virtual int residualSensFwdAdOnly(const active& t, unsigned int secIdx, const active& timeFactor,
		double const* const y, double const* const yDot, active* const adRes);

	virtual int residualSensFwdCombine(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, 
		const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS, active const* adRes, 
		double* const tmp1, double* const tmp2, double* const tmp3);

	virtual int residualSensFwdWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, active* const adRes, active* const adY, unsigned int adDirOffset);

	// linearSolve and assembleAndPrepareDAEJacobian are null operations since there are only inlet DOFs, which are treated by ModelSystem
	virtual int linearSolve(double t, double timeFactor, double alpha, double tol, double* const rhs, double const* const weight,
		double const* const y, double const* const yDot) { return 0; }

	virtual void prepareADvectors(active* const adRes, active* const adY, unsigned int adDirOffset) const;

	virtual void applyInitialCondition(double* const vecStateY, double* const vecStateYdot) const;
	virtual void readInitialCondition(IParameterProvider& paramProvider);

	virtual void consistentInitialState(double t, unsigned int secIdx, double timeFactor, double* const vecStateY, active* const adRes, active* const adY, unsigned int adDirOffset, double errorTol) { }
	virtual void consistentInitialTimeDerivative(double t, unsigned int secIdx, double timeFactor, double const* vecStateY, double* const vecStateYdot) { }

	virtual void consistentInitialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
		std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes);

	virtual void leanConsistentInitialState(double t, unsigned int secIdx, double timeFactor, double* const vecStateY, active* const adRes, active* const adY, unsigned int adDirOffset, double errorTol) { }
	virtual void leanConsistentInitialTimeDerivative(double t, double timeFactor, double const* const vecStateY, double* const vecStateYdot, double* const res) { }

	virtual void leanConsistentInitialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
		std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes);

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size) { }

	virtual void multiplyWithJacobian(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double const* yS, double alpha, double beta, double* ret);
	virtual void multiplyWithDerivativeJacobian(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double const* sDot, double* ret);

	virtual bool hasInlet() const CADET_NOEXCEPT { return true; }
	virtual bool hasOutlet() const CADET_NOEXCEPT { return false; }

	virtual unsigned int localOutletComponentIndex() const CADET_NOEXCEPT { return 0; }
	virtual unsigned int localOutletComponentStride() const CADET_NOEXCEPT { return 0; }
	virtual unsigned int localInletComponentIndex() const CADET_NOEXCEPT { return 0; }
	virtual unsigned int localInletComponentStride() const CADET_NOEXCEPT { return 1; }

	virtual void setSectionTimes(double const* secTimes, bool const* secContinuity, unsigned int nSections);

	virtual void expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut) { }

#ifdef CADET_BENCHMARK_MODE
	virtual std::vector<double> benchmarkTimings() const { return std::vector<double>(0); }
	virtual char const* const* benchmarkDescriptions() const { return nullptr; }
#endif

protected:

	UnitOpIdx _unitOpIdx; //!< Unit operation index
	unsigned int _nComp; //!< Number of components

	class Exporter : public ISolutionExporter
	{
	public:

		Exporter(unsigned int nComp, double const* data) : _data(data), _nComp(nComp) { }

		virtual bool hasMultipleBoundStates() const CADET_NOEXCEPT { return false; }
		virtual bool hasNonBindingComponents() const CADET_NOEXCEPT { return true; }
		virtual bool hasParticleFlux() const CADET_NOEXCEPT { return false; }
		virtual bool hasParticleMobilePhase() const CADET_NOEXCEPT { return false; }
		virtual bool hasSolidPhase() const CADET_NOEXCEPT { return false; }
		virtual bool hasVolume() const CADET_NOEXCEPT { return false; }

		virtual unsigned int numComponents() const CADET_NOEXCEPT { return _nComp; }
		virtual unsigned int numAxialCells() const CADET_NOEXCEPT { return 1; }
		virtual unsigned int numRadialCells() const CADET_NOEXCEPT { return 0; }
		virtual unsigned int numBoundStates() const CADET_NOEXCEPT { return 0; }
		virtual unsigned int const* numBoundStatesPerComponent() const CADET_NOEXCEPT { return nullptr; }
		virtual unsigned int numBoundStates(unsigned int comp) const CADET_NOEXCEPT { return 0; }
		virtual unsigned int numBulkDofs() const CADET_NOEXCEPT { return _nComp; }
		virtual unsigned int numParticleMobilePhaseDofs() const CADET_NOEXCEPT { return 0; }
		virtual unsigned int numSolidPhaseDofs() const CADET_NOEXCEPT { return 0; }
		virtual unsigned int numFluxDofs() const CADET_NOEXCEPT { return 0; }
		virtual unsigned int numVolumeDofs() const CADET_NOEXCEPT { return 0; }
		
		virtual double concentration(unsigned int component, unsigned int axialCell) const { return _data[component]; }
		virtual double flux(unsigned int component, unsigned int axialCell) const { return 0.0; }
		virtual double mobilePhase(unsigned int component, unsigned int axialCell, unsigned int radialCell) const { return 0.0; }
		virtual double solidPhase(unsigned int component, unsigned int axialCell, unsigned int radialCell, unsigned int boundState) const { return 0.0; }
		virtual double volume(unsigned int dof) const { return 0.0; }
		
		virtual double const* concentration() const { return _data; }
		virtual double const* flux() const { return nullptr; }
		virtual double const* mobilePhase() const { return nullptr; }
		virtual double const* solidPhase() const { return nullptr; }
		virtual double const* volume() const { return nullptr; }
		virtual double const* inlet(unsigned int& stride) const
		{
			stride = 1;
			return _data;
		}
		virtual double const* outlet(unsigned int& stride) const
		{
			stride = 1;
			return _data;
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
			len = 0;
			return nullptr;
		}

		virtual unsigned int bulkMobilePhaseStride() const { return _nComp; }
		virtual unsigned int particleMobilePhaseStride() const { return 0; }
		virtual unsigned int solidPhaseStride() const { return 0; }

	protected:
		double const* const _data;
		unsigned int _nComp;

		const std::array<StateOrdering, 1> _concentrationOrdering = { { StateOrdering::Component } };
	};
};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_OUTLETMODEL_HPP_
