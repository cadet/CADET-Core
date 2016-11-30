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

/**
 * @file 
 * Defines the general rate model (GRM).
 */

#ifndef LIBCADET_MODELSYSTEM_IMPL_HPP_
#define LIBCADET_MODELSYSTEM_IMPL_HPP_

#include "cadet/ModelSystem.hpp"
#include "SimulatableModel.hpp"
#include "UnitOperation.hpp"
#include "AutoDiff.hpp"
#include "SlicedVector.hpp"

#include <vector>
#include <unordered_map>

#include "Benchmark.hpp"

namespace cadet
{

namespace model
{

/**
 * @brief Defines a system of unit operations models
 * @details
 */
class ModelSystem : public ISimulatableModel, public IModelSystem
{
public:

	ModelSystem();
	virtual ~ModelSystem() CADET_NOEXCEPT;

	virtual UnitOpIdx maxUnitOperationId() const CADET_NOEXCEPT;

	virtual void addModel(IModel* unitOp);
	virtual IModel* getModel(unsigned int index);
	virtual IModel const* getModel(unsigned int index) const;
	virtual IModel* getUnitOperationModel(unsigned int unitOpIdx);
	virtual IModel const* getUnitOperationModel(unsigned int unitOpIdx) const;
	virtual unsigned int numModels() const CADET_NOEXCEPT;
	virtual void removeModel(IModel const* unitOp);
	virtual IModel* removeModel(UnitOpIdx unitOp);

	virtual unsigned int addExternalFunction(IExternalFunction& extFun);
	virtual IExternalFunction* getExternalFunction(unsigned int index);
	virtual IExternalFunction const* getExternalFunction(unsigned int index) const;
	virtual unsigned int numExternalFunctions() const CADET_NOEXCEPT;
	virtual void removeExternalFunction(IExternalFunction const* extFun);

	virtual unsigned int numDofs() const CADET_NOEXCEPT;
	virtual unsigned int numPureDofs() const CADET_NOEXCEPT;
	virtual bool usesAD() const CADET_NOEXCEPT;
	virtual unsigned int requiredADdirs() const CADET_NOEXCEPT;

	virtual bool configure(IParameterProvider& paramProvider, IConfigHelper& helper);
	virtual bool reconfigure(IParameterProvider& paramProvider);
	virtual void notifyDiscontinuousSectionTransition(double t, unsigned int secIdx);
	virtual bool reconfigureModel(IParameterProvider& paramProvider, unsigned int unitOpIdx);

	virtual bool hasParameter(const ParameterId& pId) const;
	virtual std::unordered_map<ParameterId, double> getAllParameterValues() const;

	virtual bool setParameter(const ParameterId& pId, int value);
	virtual bool setParameter(const ParameterId& pId, double value);
	virtual bool setParameter(const ParameterId& pId, bool value);

	virtual bool setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue);
	virtual void setSensitiveParameterValue(const ParameterId& id, double value);

	virtual void clearSensParams();

	virtual void reportSolution(ISolutionRecorder& recorder, double const* const solution) const;
	virtual void reportSolutionStructure(ISolutionRecorder& recorder) const;

	virtual int residual(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double* const res);
	virtual int residualWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, double* const res, active* const adRes, active* const adY, unsigned int numSensAdDirs);
	virtual double residualNorm(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot);

	virtual int residualSensFwd(unsigned int nSens, const active& t, unsigned int secIdx, const active& timeFactor,
		double const* const y, double const* const yDot, double const* const res,
		const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS,
		active* const adRes, double* const tmp1, double* const tmp2, double* const tmp3);

	virtual void residualSensFwdNorm(unsigned int nSens, const active& t, unsigned int secIdx, const active& timeFactor, 
		double const* const y, double const* const yDot,
		const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, double* const norms,
		active* const adRes, double* const tmp);

	virtual int linearSolve(double t, double timeFactor, double alpha, double tol, double* const rhs, double const* const weight,
		double const* const y, double const* const yDot, double const* const res);

	virtual void prepareADvectors(active* const adRes, active* const adY, unsigned int numSensAdDirs) const;

	virtual void applyInitialCondition(double* const vecStateY, double* const vecStateYdot);
	virtual void applyInitialCondition(IParameterProvider& paramProvider, double* const vecStateY, double* const vecStateYdot);

	virtual void consistentInitialConditions(double t, unsigned int secIdx, double timeFactor, double* const vecStateY, double* const vecStateYdot, active* const adRes, active* const adY, unsigned int numSensAdDirs, double errorTol);

	virtual void consistentIntialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
		std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active* const adRes, active* const adY);

	virtual void leanConsistentInitialConditions(double t, unsigned int secIdx, double timeFactor, double* const vecStateY, double* const vecStateYdot, active* const adRes, active* const adY, unsigned int numSensAdDirs, double errorTol);

	virtual void leanConsistentIntialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
		std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active* const adRes, active* const adY);

	virtual void setSectionTimes(double const* secTimes, bool const* secContinuity, unsigned int nSections);

	virtual void expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut);
	virtual std::vector<double> calculateErrorTolsForAdditionalDofs(double const* errorTol, unsigned int errorTolLength);

#ifdef CADET_BENCHMARK_MODE
	virtual std::vector<double> benchmarkTimings() const
	{
		return std::vector<double>({
			_timerResidual.totalElapsedTime(),
			_timerResidualSens.totalElapsedTime(),
			_timerConsistentInit.totalElapsedTime(),
			_timerLinearSolve.totalElapsedTime()
		});
	}

	virtual char const* const* benchmarkDescriptions() const
	{
		static const char* const desc[] = {
			"Residual",
			"ResidualSens",
			"ConsistentInit",
			"LinearSolve"
		};
		return desc;
	}
#endif

protected:

	void configureSwitches(IParameterProvider& paramProvider);

	template <typename StateType, typename ResidualType, typename ParamType>
	void residualConnectUnitOps(unsigned int secIdx, StateType const* const y, double const* const yDot, ResidualType* const res);

	void rebuildInternalDataStructures();

	template <typename ParamType>
	bool setParameterImpl(const ParameterId& pId, const ParamType value);

	void checkConnectionList(std::vector<int>& conn) const;

	std::vector<IUnitOperation*> _models; //!< Unit operation models
	std::vector<IExternalFunction*> _extFunctions; //!< External functions
	std::vector<unsigned int> _dofOffset; //!< Vector with DOF offsets for each unit operation
	util::SlicedVector<int> _connections; //!< Vector of connection lists for each section
	std::vector<unsigned int> _switchSectionIndex; //!< Holds indices of sections where valves are switched

	unsigned int _curSwitchIndex; //!< Current index in _switchSectionIndex list 

	BENCH_TIMER(_timerResidual)
	BENCH_TIMER(_timerResidualSens)
	BENCH_TIMER(_timerConsistentInit)
	BENCH_TIMER(_timerLinearSolve)
};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_MODELSYSTEM_IMPL_HPP_
