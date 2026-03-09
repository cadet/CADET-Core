// =============================================================================
//  CADET
//
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Public License v3.0 (or, at
//  your option, any later version) which accompanies this distribution, and
//  is available at http://www.gnu.org/licenses/gpl.html
// =============================================================================

/**
 * @file
 * Defines the radial general rate model (rGRM) using Discontinuous Galerkin (DG) discretization.
 */

#ifndef LIBCADET_RADIALGENERALRATEMODELWITHPORESDG_HPP_
#define LIBCADET_RADIALGENERALRATEMODELWITHPORESDG_HPP_

#include "model/UnitOperationBase.hpp"
#include "model/particle/ParticleModel.hpp"
#include "cadet/StrongTypes.hpp"
#include "cadet/SolutionExporter.hpp"
#include "model/parts/ConvectionDispersionOperatorDG.hpp"
#include "AutoDiff.hpp"
#include "linalg/BandedEigenSparseRowIterator.hpp"
#include "linalg/EigenSolverWrapper.hpp"
#include "Memory.hpp"
#include "model/ModelUtils.hpp"
#include "ParameterMultiplexing.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>
#include <vector>

#include "Benchmark.hpp"

using namespace Eigen;

namespace cadet
{

namespace model
{

namespace parts
{
	namespace cell
	{
		struct CellParameters;
	}
}

class IDynamicReactionModel;
class IParameterStateDependence;

/**
 * @brief Radial General Rate Model with DG discretization
 * @details Implements the radial transport equations with full pore discretization:
 *
 * @f[\begin{align}
		\frac{\partial c_i}{\partial t} &= - \frac{u}{\rho} \frac{\partial c_i}{\partial \rho} + \frac{1}{\rho} \frac{\partial}{\partial \rho}\left( D_{\rho,i} \rho \frac{\partial c_i}{\partial \rho} \right) - \frac{1 - \varepsilon_c}{\varepsilon_c} \frac{3}{r_p} j_{f,i} \\
		\frac{\partial c_{p,i}}{\partial t} + \frac{1 - \varepsilon_p}{\varepsilon_p} \frac{\partial q_{i}}{\partial t} &= D_{p,i} \left( \frac{\partial^2 c_{p,i}}{\partial r^2} + \frac{2}{r} \frac{\partial c_{p,i}}{\partial r} \right) + D_{s,i} \frac{1 - \varepsilon_p}{\varepsilon_p} \left( \frac{\partial^2 q_{i}}{\partial r^2} + \frac{2}{r} \frac{\partial q_{i}}{\partial r} \right) \\
		a \frac{\partial q_i}{\partial t} &= f_{\text{iso}}(c_p, q)
	\end{align} @f]
 * with film diffusion flux @f$ j_{f,i} = k_{f,i}(c_i - c_{p,i}|_{r=r_p}) @f$
 */
class RadialGeneralRateModelDG : public UnitOperationBase
{
public:

	RadialGeneralRateModelDG(UnitOpIdx unitOpIdx);
	virtual ~RadialGeneralRateModelDG() CADET_NOEXCEPT;

	virtual unsigned int numDofs() const CADET_NOEXCEPT;
	virtual unsigned int numPureDofs() const CADET_NOEXCEPT;
	virtual bool usesAD() const CADET_NOEXCEPT;
	virtual unsigned int requiredADdirs() const CADET_NOEXCEPT;

	virtual UnitOpIdx unitOperationId() const CADET_NOEXCEPT { return _unitOpIdx; }
	virtual unsigned int numComponents() const CADET_NOEXCEPT { return _disc.nComp; }
	virtual void setFlowRates(active const* in, active const* out) CADET_NOEXCEPT;
	virtual unsigned int numInletPorts() const CADET_NOEXCEPT { return 1; }
	virtual unsigned int numOutletPorts() const CADET_NOEXCEPT { return 1; }
	virtual bool canAccumulate() const CADET_NOEXCEPT { return false; }

	static const char* identifier() { return "RADIAL_GENERAL_RATE_MODEL_DG"; }
	virtual const char* unitOperationName() const CADET_NOEXCEPT { return identifier(); }

	virtual bool configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper);
	virtual bool configure(IParameterProvider& paramProvider);
	virtual void notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, const ConstSimulationState& simState, const AdJacobianParams& adJac);

	virtual void useAnalyticJacobian(const bool analyticJac);

	virtual void reportSolution(ISolutionRecorder& recorder, double const* const solution) const;
	virtual void reportSolutionStructure(ISolutionRecorder& recorder) const;

	virtual int residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, util::ThreadLocalStorage& threadLocalMem);

	virtual int jacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem);
	virtual int residualWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem);
	virtual int residualSensFwdAdOnly(const SimulationTime& simTime, const ConstSimulationState& simState, active* const adRes, util::ThreadLocalStorage& threadLocalMem);
	virtual int residualSensFwdWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem);

	virtual int residualSensFwdCombine(const SimulationTime& simTime, const ConstSimulationState& simState,
		const std::vector<const double*>& yS, const std::vector<const double*>& ySdot, const std::vector<double*>& resS, active const* adRes,
		double* const tmp1, double* const tmp2, double* const tmp3);

	virtual int linearSolve(double t, double alpha, double tol, double* const rhs, double const* const weight,
		const ConstSimulationState& simState);

	virtual void prepareADvectors(const AdJacobianParams& adJac) const;

	virtual void applyInitialCondition(const SimulationState& simState) const;
	virtual void readInitialCondition(IParameterProvider& paramProvider);

	virtual void consistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, util::ThreadLocalStorage& threadLocalMem);
	virtual void consistentInitialTimeDerivative(const SimulationTime& simTime, double const* vecStateY, double* const vecStateYdot, util::ThreadLocalStorage& threadLocalMem);

	virtual void initializeSensitivityStates(const std::vector<double*>& vecSensY) const;
	virtual void consistentInitialSensitivity(const SimulationTime& simTime, const ConstSimulationState& simState,
		std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes, util::ThreadLocalStorage& threadLocalMem);

	virtual void leanConsistentInitialState(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, util::ThreadLocalStorage& threadLocalMem);
	virtual void leanConsistentInitialTimeDerivative(double t, double const* const vecStateY, double* const vecStateYdot, double* const res, util::ThreadLocalStorage& threadLocalMem);

	virtual void leanConsistentInitialSensitivity(const SimulationTime& simTime, const ConstSimulationState& simState,
		std::vector<double*>& vecSensY, std::vector<double*>& vecSensYdot, active const* const adRes, util::ThreadLocalStorage& threadLocalMem);

	virtual bool hasInlet() const CADET_NOEXCEPT { return true; }
	virtual bool hasOutlet() const CADET_NOEXCEPT { return true; }

	virtual unsigned int localOutletComponentIndex(unsigned int port) const CADET_NOEXCEPT;
	virtual unsigned int localOutletComponentStride(unsigned int port) const CADET_NOEXCEPT;
	virtual unsigned int localInletComponentIndex(unsigned int port) const CADET_NOEXCEPT;
	virtual unsigned int localInletComponentStride(unsigned int port) const CADET_NOEXCEPT;

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size);
	virtual void setSectionTimes(double const* secTimes, bool const* secContinuity, unsigned int nSections) { }

	virtual void expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut);

	virtual void multiplyWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* yS, double alpha, double beta, double* ret);
	virtual void multiplyWithDerivativeJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* sDot, double* ret);

	inline void multiplyWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* yS, double* ret)
	{
		multiplyWithJacobian(simTime, simState, yS, 1.0, 0.0, ret);
	}

	virtual bool setParameter(const ParameterId& pId, double value);
	virtual bool setSensitiveParameter(const ParameterId& pId, unsigned int adDirection, double adValue);
	virtual void setSensitiveParameterValue(const ParameterId& id, double value);

	virtual unsigned int threadLocalMemorySize() const CADET_NOEXCEPT;

#ifdef CADET_BENCHMARK_MODE
	virtual std::vector<double> benchmarkTimings() const
	{
		return std::vector<double>({
			static_cast<double>(numDofs()),
			_timerResidual.totalElapsedTime(),
			_timerResidualPar.totalElapsedTime(),
			_timerResidualSens.totalElapsedTime(),
			_timerResidualSensPar.totalElapsedTime(),
			_timerConsistentInit.totalElapsedTime(),
			_timerConsistentInitPar.totalElapsedTime(),
			_timerLinearSolve.totalElapsedTime()
		});
	}

	virtual char const* const* benchmarkDescriptions() const
	{
		static const char* const desc[] = {
			"DOFs",
			"Residual",
			"ResidualPar",
			"ResidualSens",
			"ResidualSensPar",
			"ConsistentInit",
			"ConsistentInitPar",
			"LinearSolve"
		};
		return desc;
	}
#endif

protected:

	class Indexer;

	int residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem, bool updateJacobian, bool paramSensitivity);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes = true>
	int residualImpl(double t, unsigned int secIdx, StateType const* const y, double const* const yDot, ResidualType* const res, util::ThreadLocalStorage& threadLocalMem);

	void extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset);

	void assembleDiscretizedJacobian(double alpha, const Indexer& idxr);

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	void checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const;
#endif

	class Discretization
	{
	public:
		unsigned int nComp; //!< Number of components
		unsigned int nParType; //!< Number of particle types
		unsigned int nElem; //!< Number of radial bulk elements
		unsigned int polyDeg; //!< polynomial degree for bulk
		unsigned int nNodes; //!< Number of nodes per bulk cell
		unsigned int nPoints; //!< Number of discrete bulk points (nElem * nNodes)

		unsigned int* nBound; //!< Array with number of bound states for each component
		unsigned int* boundOffset; //!< Array with offset to the first bound state of each component in the solid phase
		unsigned int strideBound; //!< Total number of bound states

		unsigned int* nParPoints; //!< Number of particle points per particle type
		unsigned int* parTypeOffset; //!< Offset to particle type DOFs
		unsigned int nTotalParPoints; //!< Total particle points across all types

		int curSection; //!< current section index
		bool newStaticJac; //!< determines whether static analytical jacobian needs to be computed

		~Discretization()
		{
			delete[] nBound;
			delete[] boundOffset;
			delete[] nParPoints;
			delete[] parTypeOffset;
		}
	};

	Discretization _disc; //!< Discretization info

	std::vector<IParticleModel*> _particles; //!< Particle models (one per particle type)

	parts::RadialConvectionDispersionOperatorBaseDG _convDispOp; //!< Convection dispersion operator for radial bulk transport
	IDynamicReactionModel* _dynReactionBulk; //!< Dynamic reactions in the bulk volume
	std::vector<IDynamicReactionModel*> _dynReaction; //!< Dynamic reaction models for particles (owned)

	cadet::linalg::EigenSolverBase* _linearSolver; //!< Linear solver

	Eigen::SparseMatrix<double, RowMajor> _jac; //!< Jacobian
	Eigen::SparseMatrix<double, RowMajor> _jacDisc; //!< Jacobian with time derivatives from BDF method
	Eigen::MatrixXd _jacInlet; //!< Jacobian inlet DOF block matrix

	active _colPorosity; //!< Column porosity (external porosity)
	std::vector<active> _parTypeVolFrac; //!< Particle type volume fractions

	bool _analyticJac; //!< Determines whether AD or analytic Jacobians are used
	unsigned int _jacobianAdDirs; //!< Number of AD seed vectors required for Jacobian computation
	bool _factorizeJacobian; //!< Determines whether the Jacobian needs to be factorized
	double* _tempState; //!< Temporary storage

	std::vector<active> _initC; //!< Liquid bulk phase initial conditions
	std::vector<std::vector<active>> _initCp; //!< Liquid pore phase initial conditions (per particle type)
	std::vector<std::vector<active>> _initCs; //!< Solid phase initial conditions (per particle type)
	std::vector<double> _initState; //!< Initial conditions for state vector if given
	std::vector<double> _initStateDot; //!< Initial conditions for time derivative

	BENCH_TIMER(_timerResidual)
	BENCH_TIMER(_timerResidualPar)
	BENCH_TIMER(_timerResidualSens)
	BENCH_TIMER(_timerResidualSensPar)
	BENCH_TIMER(_timerConsistentInit)
	BENCH_TIMER(_timerConsistentInitPar)
	BENCH_TIMER(_timerLinearSolve)

	class Indexer
	{
	public:
		Indexer(const Discretization& disc) : _disc(disc) { }

		// Strides
		inline int strideColNode() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp); }
		inline int strideColCell() const CADET_NOEXCEPT { return static_cast<int>(_disc.nNodes * strideColNode()); }
		inline int strideColComp() const CADET_NOEXCEPT { return 1; }

		// Offsets
		inline int offsetC() const CADET_NOEXCEPT { return _disc.nComp; }
		inline int offsetCp(unsigned int parType) const CADET_NOEXCEPT
		{
			return _disc.nComp + _disc.nPoints * _disc.nComp + _disc.parTypeOffset[parType];
		}
		inline int offsetBoundComp(unsigned int comp) const CADET_NOEXCEPT { return _disc.boundOffset[comp]; }

		// Particle strides
		inline int strideParBlock(unsigned int parType) const CADET_NOEXCEPT
		{
			return _disc.nParPoints[parType] * (_disc.nComp + _disc.strideBound);
		}
		inline int strideParNode(unsigned int parType) const CADET_NOEXCEPT
		{
			return _disc.nComp + _disc.strideBound;
		}

		// Return pointer to first element of state variable in state vector
		template <typename real_t> inline real_t* c(real_t* const data) const { return data + offsetC(); }
		template <typename real_t> inline real_t const* c(real_t const* const data) const { return data + offsetC(); }

		template <typename real_t> inline real_t* cp(real_t* const data, unsigned int parType) const { return data + offsetCp(parType); }
		template <typename real_t> inline real_t const* cp(real_t const* const data, unsigned int parType) const { return data + offsetCp(parType); }

	protected:
		const Discretization& _disc;
	};

	class Exporter : public ISolutionExporter
	{
	public:

		Exporter(const Discretization& disc, const RadialGeneralRateModelDG& model, double const* data)
			: _disc(disc), _idx(disc), _model(model), _data(data) { }
		Exporter(const Discretization&& disc, const RadialGeneralRateModelDG& model, double const* data) = delete;

		virtual bool hasParticleFlux() const CADET_NOEXCEPT { return false; }
		virtual bool hasParticleMobilePhase() const CADET_NOEXCEPT { return true; }
		virtual bool hasSolidPhase() const CADET_NOEXCEPT { return _disc.strideBound > 0; }
		virtual bool hasVolume() const CADET_NOEXCEPT { return false; }
		virtual bool isParticleLumped(unsigned int parType) const CADET_NOEXCEPT { return false; }
		virtual bool hasPrimaryExtent() const CADET_NOEXCEPT { return true; }

		virtual unsigned int numComponents() const CADET_NOEXCEPT { return _disc.nComp; }
		virtual unsigned int numPrimaryCoordinates() const CADET_NOEXCEPT { return _disc.nPoints; }
		virtual unsigned int numSecondaryCoordinates() const CADET_NOEXCEPT { return 0; }
		virtual unsigned int numInletPorts() const CADET_NOEXCEPT { return 1; }
		virtual unsigned int numOutletPorts() const CADET_NOEXCEPT { return 1; }
		virtual unsigned int numParticleTypes() const CADET_NOEXCEPT { return _disc.nParType; }
		virtual unsigned int numParticleShells(unsigned int parType) const CADET_NOEXCEPT
		{
			return _disc.nParPoints[parType];
		}
		virtual unsigned int numBoundStates(unsigned int parType) const CADET_NOEXCEPT { return _disc.strideBound; }
		virtual unsigned int numMobilePhaseDofs() const CADET_NOEXCEPT { return _disc.nComp * _disc.nPoints; }
		virtual unsigned int numParticleMobilePhaseDofs(unsigned int parType) const CADET_NOEXCEPT
		{
			return _disc.nComp * _disc.nPoints * _disc.nParPoints[parType];
		}
		virtual unsigned int numParticleMobilePhaseDofs() const CADET_NOEXCEPT;
		virtual unsigned int numSolidPhaseDofs(unsigned int parType) const CADET_NOEXCEPT
		{
			return _disc.strideBound * _disc.nPoints * _disc.nParPoints[parType];
		}
		virtual unsigned int numSolidPhaseDofs() const CADET_NOEXCEPT;
		virtual unsigned int numParticleFluxDofs() const CADET_NOEXCEPT { return 0u; }
		virtual unsigned int numVolumeDofs() const CADET_NOEXCEPT { return 0; }

		virtual int writeMobilePhase(double* buffer) const;
		virtual int writeSolidPhase(double* buffer) const;
		virtual int writeParticleMobilePhase(double* buffer) const;
		virtual int writeSolidPhase(unsigned int parType, double* buffer) const;
		virtual int writeParticleMobilePhase(unsigned int parType, double* buffer) const;
		virtual int writeParticleFlux(double* buffer) const { return 0; }
		virtual int writeParticleFlux(unsigned int parType, double* buffer) const { return 0; }
		virtual int writeVolume(double* buffer) const { return 0; }
		virtual int writeInlet(unsigned int port, double* buffer) const;
		virtual int writeInlet(double* buffer) const;
		virtual int writeOutlet(unsigned int port, double* buffer) const;
		virtual int writeOutlet(double* buffer) const;

		virtual int writePrimaryCoordinates(double* coords) const;
		virtual int writeSecondaryCoordinates(double* coords) const { return 0; }
		virtual int writeParticleCoordinates(unsigned int parType, double* coords) const;

	protected:
		const Discretization& _disc;
		const Indexer _idx;
		const RadialGeneralRateModelDG& _model;
		double const* const _data;
	};

	void updateSection(int secIdx) {
		if (_disc.curSection != secIdx) {
			_disc.curSection = secIdx;
			_disc.newStaticJac = true;
		}
	}

	// Jacobian pattern functions
	typedef Eigen::Triplet<double> T;

	void setPattern(Eigen::SparseMatrix<double, RowMajor>& mat, bool stateDer);
	void convDispJacPattern(std::vector<T>& tripletList);
	void particleJacPattern(std::vector<T>& tripletList, unsigned int secIdx);
	void stateDerPattern(std::vector<T>& tripletList);

};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_RADIALGENERALRATEMODELWITHPORESDG_HPP_
