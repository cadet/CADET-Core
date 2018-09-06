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
 * Defines the general rate model (GRM).
 */

#ifndef LIBCADET_MODELSYSTEM_IMPL_HPP_
#define LIBCADET_MODELSYSTEM_IMPL_HPP_

#include "SimulatableModel.hpp"
#include "UnitOperation.hpp"
#include "AutoDiff.hpp"
#include "SlicedVector.hpp"
#include "ParamIdUtil.hpp"

#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>

#ifdef CADET_PARALLELIZE
	#include <tbb/spin_mutex.h>
#endif

#include "linalg/SparseMatrix.hpp"
#include "linalg/Gmres.hpp"

#include "Benchmark.hpp"

namespace cadet
{

namespace model
{

/**
 * @brief Defines a system of unit operations models
 * @details
 */
class ModelSystem : public ISimulatableModel
{
public:

	ModelSystem();
	virtual ~ModelSystem() CADET_NOEXCEPT;

	virtual UnitOpIdx maxUnitOperationId() const CADET_NOEXCEPT;

	virtual void addModel(IModel* unitOp);
	virtual IModel* getModel(unsigned int index);
	virtual IModel const* getModel(unsigned int index) const;
	virtual IUnitOperation* getUnitOperationModel(unsigned int unitOpIdx);
	virtual IUnitOperation const* getUnitOperationModel(unsigned int unitOpIdx) const;
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

	virtual bool configureModelDiscretization(IParameterProvider& paramProvider, IConfigHelper& helper);
	virtual bool configure(IParameterProvider& paramProvider);
	virtual void notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, active* const adRes, active* const adY, unsigned int adDirOffset);
	virtual bool configureModel(IParameterProvider& paramProvider, unsigned int unitOpIdx);

	virtual bool hasParameter(const ParameterId& pId) const;
	virtual std::unordered_map<ParameterId, double> getAllParameterValues() const;
	virtual double getParameterDouble(const ParameterId& pId) const;

	virtual bool setParameter(const ParameterId& pId, int value);
	virtual bool setParameter(const ParameterId& pId, double value);
	virtual bool setParameter(const ParameterId& pId, bool value);

	virtual bool setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue);
	virtual void setSensitiveParameterValue(const ParameterId& id, double value);

	virtual void clearSensParams();

	virtual void reportSolution(ISolutionRecorder& recorder, double const* const solution) const;
	virtual void reportSolutionStructure(ISolutionRecorder& recorder) const;

	virtual int residual(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot, double* const res);

	virtual int residualWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, double* const res, active* const adRes, active* const adY, unsigned int adDirOffset);
	virtual double residualNorm(double t, unsigned int secIdx, double timeFactor, double const* const y, double const* const yDot);

	virtual int residualSensFwd(unsigned int nSens, const active& t, unsigned int secIdx, const active& timeFactor,
		double const* const y, double const* const yDot, double const* const res,
		const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS,
		active* const adRes, double* const tmp1, double* const tmp2, double* const tmp3);

	virtual int residualSensFwdWithJacobian(unsigned int nSens, const active& t, unsigned int secIdx, const active& timeFactor,
		double const* const y, double const* const yDot, double const* const res,
		const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS,
		active* const adRes, active* const adY, unsigned int adDirOffset, double* const tmp1, double* const tmp2, double* const tmp3);

	virtual void residualSensFwdNorm(unsigned int nSens, const active& t, unsigned int secIdx, const active& timeFactor, 
		double const* const y, double const* const yDot,
		const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, double* const norms,
		active* const adRes, double* const tmp);

	virtual int linearSolve(double t, double timeFactor, double alpha, double tol, double* const rhs, double const* const weight,
		double const* const y, double const* const yDot);

	virtual void prepareADvectors(active* const adRes, active* const adY, unsigned int adDirOffset) const;

	virtual void applyInitialCondition(double* const vecStateY, double* const vecStateYdot) const;
	virtual void readInitialCondition(IParameterProvider& paramProvider);

	virtual void consistentInitialConditions(double t, unsigned int secIdx, double timeFactor, double* const vecStateY, double* const vecStateYdot, active* const adRes, active* const adY, unsigned int adDirOffset, double errorTol);

	virtual void consistentInitialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
		std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active* const adRes, active* const adY);

	virtual void leanConsistentInitialConditions(double t, unsigned int secIdx, double timeFactor, double* const vecStateY, double* const vecStateYdot, active* const adRes, active* const adY, unsigned int adDirOffset, double errorTol);

	virtual void leanConsistentInitialSensitivity(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
		std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active* const adRes, active* const adY);
	
	virtual void setSectionTimes(double const* secTimes, bool const* secContinuity, unsigned int nSections);

	int schurComplementMatrixVector(double const* x, double* z, double t, double timeFactor, double alpha, double outerTol, double const* const weight,
		double const* const y, double const* const yDot) const;

	virtual void expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut);
	virtual std::vector<double> calculateErrorTolsForAdditionalDofs(double const* errorTol, unsigned int errorTolLength);

#ifdef CADET_BENCHMARK_MODE
	virtual std::vector<double> benchmarkTimings() const
	{
		return std::vector<double>({
			_timerResidual.totalElapsedTime(),
			_timerResidualSens.totalElapsedTime(),
			_timerConsistentInit.totalElapsedTime(),
			_timerLinearAssemble.totalElapsedTime(),
			_timerLinearSolve.totalElapsedTime(),
			_timerMatVec.totalElapsedTime()
		});
	}

	virtual char const* const* benchmarkDescriptions() const
	{
		static const char* const desc[] = {
			"Residual",
			"ResidualSens",
			"ConsistentInit",
			"LinearAssemble",
			"LinearSolve",
			"MatVec"
		};
		return desc;
	}
#endif

	virtual void genJacobian(double t, unsigned int secIdx, double timeFactor, double const * const y, double const * const yDot);

	virtual void genJacobian(unsigned int nSens, const active& t, unsigned int secIdx,
		const active& timeFactor, double const* const y, double const* const yDot, double const* const res,
		const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS,
		active* const adRes, double* const tmp1, double* const tmp2, double* const tmp3);

	virtual void multiplyWithJacobian(double t, unsigned int secIdx, double timeFactor, double const* y, double const* yDot, double const* yS, double alpha, double beta, double* ret);
	virtual void multiplyWithDerivativeJacobian(double t, unsigned int secIdx, double timeFactor, double const* y, double const* yDot, double const* yS, double* ret);
protected:

	void configureSwitches(IParameterProvider& paramProvider);

	template <typename StateType, typename ResidualType, typename ParamType>
	void residualConnectUnitOps(unsigned int secIdx, StateType const* const y, double const* const yDot, ResidualType* const res) CADET_NOEXCEPT;
	void solveCouplingDOF(double * const vec);

	void multiplyWithMacroJacobian(double const* yS, double alpha, double beta, double* ret);
	inline void multiplyWithMacroJacobian(double const* yS, double* ret)
	{
		multiplyWithMacroJacobian(yS, 1.0, 0.0, ret);
	}

	int dResDpFwdWithJacobian(const active& t, unsigned int secIdx, const active& timeFactor, double const* const y, double const* const yDot, active* const adRes, active* const adY, unsigned int adDirOffset);

	void rebuildInternalDataStructures();
	void allocateSuperStructMatrices();
	void assembleSuperStructMatrices(unsigned int secIdx);

	template <typename ParamType>
	bool setParameterImpl(const ParameterId& pId, const ParamType value);

	void checkConnectionList(const std::vector<double>& conn, std::vector<int>& connOnly, std::vector<double>& totalOutflow, unsigned int idxSwitch) const;

	template <typename tag_t>
	void consistentInitialConditionAlgorithm(double t, unsigned int secIdx, double timeFactor, double* const vecStateY, double* const vecStateYdot, active* const adRes, active* const adY, unsigned int adDirOffset, double errorTol);

	template <typename tag_t>
	void consistentInitialSensitivityAlgorithm(const active& t, unsigned int secIdx, const active& timeFactor, double const* vecStateY, double const* vecStateYdot,
		std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active* const adRes, active* const adY);

	template <bool evalJacobian>
	int residualSensFwdWithJacobianAlgorithm(unsigned int nSens, const active& t, unsigned int secIdx, const active& timeFactor,
		double const* const y, double const* const yDot, double const* const res,
		const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS,
		active* const adRes, active* const adY, unsigned int adDirOffset, double* const tmp1, double* const tmp2, double* const tmp3);

	/**
	 * @brief Returns the number of coupling DOFs
	 * @details The number of coupling DOFs is given by the total number of inlet DOFs in the system.
	 * @return Number of coupling DOFs
	 */
	inline unsigned int numCouplingDOF() const CADET_NOEXCEPT { return _couplingIdxMap.size(); }

	std::vector<IUnitOperation*> _models; //!< Unit operation models
	std::vector<IExternalFunction*> _extFunctions; //!< External functions
	linalg::SparseMatrix<double>* _jacNF; //!< Jacobian block connecting coupling DOF to inlets
	linalg::SparseMatrix<double>* _jacFN; //!< Jacobian block connecting outlets to coupling DOF
	linalg::SparseMatrix<active>* _jacActiveFN; //!< Jacobian block connecting outlets to coupling DOF
	std::vector<unsigned int> _dofOffset; //!< Vector with DOF offsets for each unit operation
	std::vector<unsigned int> _dofs; //!< Vector with DOF for each unit
	util::SlicedVector<int> _connections; //!< Vector of connection lists for each section
	util::SlicedVector<active> _flowRates; //!< Vector of connection flow rates for each section
	std::vector<unsigned int> _switchSectionIndex; //!< Holds indices of sections where valves are switched
	unsigned int _curSwitchIndex; //!< Current index in _switchSectionIndex list 

	mutable std::vector<int> _errorIndicator; //!< Storage for return value of unit operation function calls

	double* _tempState; //!< Temporary storage for the state vector
	std::vector<active> _totalInletFlow; //!< Total flow rate into each inlet at the current section

	std::vector<std::vector<const double*>> _yStemp; //!< Needed to store offsets for unit operations
	std::vector<std::vector<const double*>> _yStempDot;  //!< Needed to store offsets for unit operations
	std::vector<std::vector<double*>> _resSTemp;  //!< Needed to store offsets for unit operations

	std::map<std::pair<unsigned int, unsigned int>, unsigned int> _couplingIdxMap; //!< Maps (UnitOpIdx, CompIdx) to local coupling DOF index

	std::unordered_map<ParameterId, active*> _parameters; //!< Provides access to all parameters
	std::unordered_set<active*> _sensParams; //!< Holds all parameters with activated AD directions

	linalg::Gmres _gmres; //!< GMRES algorithm for the Schur-complement in linearSolve()
	double _schurSafety; //!< Safety factor for Schur-complement solution

	std::vector<unsigned int> _inOutModels; //!< Indices of unit operation models in _models that have inlet and outlet

	std::vector<double> _initState; //!< Initial state vector
	std::vector<double> _initStateDot; //!< Initial time derivative state vector

#ifdef CADET_PARALLELIZE
	typedef tbb::spin_mutex SchurComplementMutex;
	mutable SchurComplementMutex _schurMutex;
#endif

	BENCH_TIMER(_timerResidual)
	BENCH_TIMER(_timerResidualSens)
	BENCH_TIMER(_timerConsistentInit)
	BENCH_TIMER(_timerLinearAssemble)
	BENCH_TIMER(_timerLinearSolve)
	BENCH_TIMER(_timerMatVec)
};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_MODELSYSTEM_IMPL_HPP_
