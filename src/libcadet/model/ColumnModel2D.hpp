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
* Defines the 2D LRMP model as a unit operation.
 */

#ifndef LIBCADET_COLUMNMODEL2D_HPP_
#define LIBCADET_COLUMNMODEL2D_HPP_

#include "model/UnitOperationBase.hpp"
#include "model/particle/ParticleModel.hpp"
#include "cadet/StrongTypes.hpp"
#include "cadet/SolutionExporter.hpp"
#include "model/parts/TwoDimensionalConvectionDispersionOperatorDG.hpp"
#include "linalg/BandedEigenSparseRowIterator.hpp"
#include "linalg/EigenSolverWrapper.hpp"
#include "AutoDiff.hpp"
#include "Memory.hpp"
#include "model/ModelUtils.hpp"
#include "model/ParameterMultiplexing.hpp"

#include <array>
#include <vector>

#include "Benchmark.hpp"

#include <Eigen/Sparse>
#include <Eigen/Dense>

namespace cadet
{

namespace model
{

class IDynamicReactionModel;

/**
 * @brief Chromatography column with 2D spatial resolution (axial and radial)
 * @details For the 2D flow equations, refer to the 2D General rate model @cite Guiochon2006, @cite Gu1995, @cite Felinger2004
 * 
 * @f[\begin{align}
	\frac{\partial c_i^b}{\partial t} &= - u \frac{\partial c_i^b}{\partial z} + \frac{\partial}{\partial z} \left(D_{\text{ax},i} \frac{\partial c_i^b}{\partial z} \right) + \frac{1}{\rho} \frac{\partial}{\partial \rho} \left( \rho D_{\text{rad},i} \frac{\partial c_i^b}{\partial \rho} \right)
	- \frac{1 - \varepsilon_c}{\varepsilon_c} \frac{3 k_{f,i}}{r_p} j_{f,i} \\
	\frac{\partial c_{p,i}}{\partial t} + \frac{1 - \varepsilon_p}{\varepsilon_p} \frac{\partial q_{i}}{\partial t} &= \frac{3 k_{f,i}}{\varepsilon_p r_p} (c_i - c_{p,i}) \\
	a \frac{\partial q_i}{\partial t} &= f_{\text{iso}}(c_p, q)
	\end{align} @f]
 * Danckwerts boundary conditions (see @cite Danckwerts1953)
@f[ \begin{align}
u c_{\text{in},i}(t) &= u c_i^b(t,0,.) - D_{\text{ax},i} \frac{\partial c_i^b}{\partial z}(t,0,.) \\
\frac{\partial c_i^b}{\partial z}(t,L,.) &= 0 \\
\frac{\partial c_i^b}{\partial \rho}(t,.,0) &= 0 \\
\frac{\partial c_i^b}{\partial \rho}(t,.,R_c) &= 0 \\
\end{align} @f]
 * Methods are described in @cite @TODO (2D DG paper), @cite Puttmann2013 @cite Puttmann2016 (forward sensitivities, AD, band compression)
 */
class ColumnModel2D : public UnitOperationBase
{
public:

	ColumnModel2D(UnitOpIdx unitOpIdx);
	virtual ~ColumnModel2D() CADET_NOEXCEPT;

	virtual unsigned int numDofs() const CADET_NOEXCEPT;
	virtual unsigned int numPureDofs() const CADET_NOEXCEPT;
	virtual bool usesAD() const CADET_NOEXCEPT;
	virtual unsigned int requiredADdirs() const CADET_NOEXCEPT;

	virtual UnitOpIdx unitOperationId() const CADET_NOEXCEPT { return _unitOpIdx; }
	virtual unsigned int numComponents() const CADET_NOEXCEPT { return _disc.nComp; }
	virtual void setFlowRates(active const* in, active const* out) CADET_NOEXCEPT;
	virtual unsigned int numInletPorts() const CADET_NOEXCEPT { return _disc.radNPoints; }
	virtual unsigned int numOutletPorts() const CADET_NOEXCEPT { return _disc.radNPoints; }
	virtual bool canAccumulate() const CADET_NOEXCEPT { return false; }

	static const char* identifier() { return "COLUMN_MODEL_2D"; }
	virtual const char* unitOperationName() const CADET_NOEXCEPT { return identifier(); }

	virtual bool configureModelDiscretization(IParameterProvider& paramProvider, const IConfigHelper& helper);
	virtual bool configure(IParameterProvider& paramProvider);
	virtual void notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, const ConstSimulationState& simState, const AdJacobianParams& adJac);

	virtual void useAnalyticJacobian(const bool analyticJac);

	virtual void reportSolution(ISolutionRecorder& recorder, double const* const solution) const;
	virtual void reportSolutionStructure(ISolutionRecorder& recorder) const;

	virtual int jacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem);

	virtual int residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, util::ThreadLocalStorage& threadLocalMem);

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
	virtual bool setParameter(const ParameterId& pId, int value);
	virtual bool setParameter(const ParameterId& pId, bool value);
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
			_timerJacobianPar.totalElapsedTime(),
			_timerConsistentInit.totalElapsedTime(),
			_timerConsistentInitPar.totalElapsedTime(),
			_timerLinearSolve.totalElapsedTime(),
			_timerFactorize.totalElapsedTime(),
			_timerFactorizePar.totalElapsedTime()
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
			"JacobianPar",
			"ConsistentInit",
			"ConsistentInitPar",
			"LinearSolve",
			"Factorize",
			"FactorizePar",
			"MatVec",
		};
		return desc;
	}
#endif

protected:

	class Indexer;

	int residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem, bool updateJacobian, bool paramSensitivity);
	
	int calcTransportJacobian(unsigned int secIdx)
	{
		Indexer idxr(_disc);
		// inlet and bulk jacobian
		_convDispOp.assembleConvDispJacobian(_globalJac, _jacInlet, true, 0);

		// particle transport diffusion Jacobian (without isotherm, which is handled in residualKernel)
		for (int colNode = 0; colNode < _disc.nBulkPoints; colNode++)
		{
			for (int parType = 0; parType < _disc.nParType; parType++)
			{
				_particles[parType]->calcParticleDiffJacobian(secIdx, colNode, idxr.offsetCp(ParticleTypeIndex{ static_cast<unsigned int>(parType) }, ParticleIndex{ static_cast<unsigned int>(colNode) }) - idxr.offsetC(), _globalJac);
			}
		}

		// film diffusion Jacobian
		for (unsigned int parType = 0; parType < _disc.nParType; parType++)
		{
			for (unsigned int z = 0; z < _disc.axNPoints; z++)
			{
				for (unsigned int r = 0; r < _disc.radNElem; r++)
				{
					const unsigned int colNode = z * _disc.radNPoints + r * _disc.radNNodes;
					_particles[parType]->calcFilmDiffJacobian(secIdx, idxr.offsetCp(ParticleTypeIndex{ static_cast<unsigned int>(parType) }, ParticleIndex{ static_cast<unsigned int>(colNode) }) - idxr.offsetC(), colNode * idxr.strideColRadialNode(), _disc.radNNodes, _disc.nParType, static_cast<double>(_convDispOp.columnPorosity(r)), &_parTypeVolFrac[colNode * _disc.nParType], _globalJac);
				}
			}
		}

		return _globalJac.isCompressed(); // check if the jacobian estimation fits the pattern
	}

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes = true>
	int residualImpl(double t, unsigned int secIdx, StateType const* const y, double const* const yDot, ResidualType* const res, util::ThreadLocalStorage& threadLocalMem);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
	int residualBulk(double t, unsigned int secIdx, StateType const* y, double const* yDot, ResidualType* res, util::ThreadLocalStorage& threadLocalMem);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
	int residualParticle(double t, unsigned int parType, unsigned int colCell, unsigned int secIdx, StateType const* y, double const* yDot, ResidualType* res, util::ThreadLocalStorage& threadLocalMem);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes>
	int residualFlux(double t, unsigned int secIdx, StateType const* y, double const* yDot, ResidualType* res);

	void extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset);

	void assembleDiscretizedGlobalJacobian(double alpha, Indexer idxr);

	void addTimeDerivativeToJacobianParticleShell(linalg::BandedEigenSparseRowIterator& jac, const Indexer& idxr, double alpha, unsigned int parType);

	unsigned int numAdDirsForJacobian() const CADET_NOEXCEPT;

	int multiplexInitialConditions(const cadet::ParameterId& pId, unsigned int adDirection, double adValue);
	int multiplexInitialConditions(const cadet::ParameterId& pId, double val, bool checkSens);

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	void checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const;
#endif

	struct Discretization
	{
		unsigned int nComp; //!< Number of components
		unsigned int axNPoints; //!< Number of axial discrete points
		unsigned int radNPoints; //!< Number of radial discrete points
		unsigned int radNElem; //!< Number of radial DG elements, which can be used as radial zones for radial parameter inhomogeneities
		unsigned int radNNodes; //!< Number of radial DG nodes per element
		unsigned int nBulkPoints; //!< Number of total bulk discrete points
		unsigned int nParType; //!< Number of particle types
		unsigned int* nParPoints; //!< Array with number of discrete points for each particle type
		unsigned int* parTypeOffset; //!< Array with offsets (in particle block) to particle type, additional last element contains total number of particle DOFs
		unsigned int* nBound; //!< Array with number of bound states for each component and particle type (particle type major ordering)
		unsigned int* boundOffset; //!< Array with offset to the first bound state of each component in the solid phase (particle type major ordering)
		unsigned int* strideBound; //!< Total number of bound states for each particle type, additional last element contains total number of bound states for all types
		unsigned int* nBoundBeforeType; //!< Array with number of bound states before a particle type (cumulative sum of strideBound)
	
		bool newStaticJac; //!< determines whether or not a new static transport Jacobian is desired, ie on section switches
	};

	enum class ParticleDiscretizationMode : int
	{
		/**
		 * Equidistant distribution of shell edges
		 */
		Equidistant,

		/**
		 * Volumes of shells are uniform
		 */
		Equivolume,

		/**
		 * Shell edges specified by user
		 */
		UserDefined
	};

	Discretization _disc; //!< Discretization info
//	IExternalFunction* _extFun; //!< External function (owned by library user)

	std::vector<IParticleModel*> _particles; //!< Particle dispersion operator

	parts::TwoDimensionalConvectionDispersionOperatorDG _convDispOp; //!< Convection dispersion operator for interstitial volume transport
	IDynamicReactionModel* _dynReactionBulk; //!< Dynamic reactions in the bulk volume
	ReactionSystem _reaction;
	std::vector<const ReactionSystem*> _reacParticle;

	cadet::linalg::EigenSolverBase* _linearSolver; //!< Linear solver

	Eigen::MatrixXd _jacInlet; //!< Jacobian inlet DOF block matrix connects inlet DOFs to first bulk cells

	// for FV the bulk jacobians are defined in the ConvDisp operator.
	Eigen::SparseMatrix<double, Eigen::RowMajor> _globalJac; //!< global Jacobian
	Eigen::SparseMatrix<double, Eigen::RowMajor> _globalJacDisc; //!< global Jacobian with time derivative from BDF method

	std::vector<active> _parTypeVolFrac; //!< Volume fraction of each particle type
	MultiplexMode _parTypeVolFracMode; //!< Multiplexing mode of particle volume fractions

	bool _analyticJac; //!< Determines whether AD or analytic Jacobians are used
	unsigned int _jacobianAdDirs; //!< Number of AD seed vectors required for Jacobian computation

	bool _factorizeJacobian; //!< Determines whether the Jacobian needs to be factorized
	double* _tempState; //!< Temporary storage with the size of the state vector or larger if binding models require it

	std::vector<active> _initC; //!< Liquid bulk phase initial conditions
	bool _singleRadiusInitC;
	std::vector<active> _initCp; //!< Liquid particle phase initial conditions
	bool _singleRadiusInitCp;
	std::vector<active> _initCs; //!< Solid phase initial conditions
	bool _singleRadiusInitCs;
	std::vector<double> _initState; //!< Initial conditions for state vector if given
	std::vector<double> _initStateDot; //!< Initial conditions for time derivative

	BENCH_TIMER(_timerResidual)
	BENCH_TIMER(_timerResidualPar)
	BENCH_TIMER(_timerResidualSens)
	BENCH_TIMER(_timerResidualSensPar)
	BENCH_TIMER(_timerJacobianPar)
	BENCH_TIMER(_timerConsistentInit)
	BENCH_TIMER(_timerConsistentInitPar)
	BENCH_TIMER(_timerLinearSolve)
	BENCH_TIMER(_timerFactorize)
	BENCH_TIMER(_timerFactorizePar)
	BENCH_TIMER(_timerMatVec)

	class Indexer
	{
	public:
		Indexer(const Discretization& disc) : _disc(disc) { }

		// Strides
		inline int strideColAxialNode() const CADET_NOEXCEPT { return strideColRadialNode() * static_cast<int>(_disc.radNPoints); }
		inline int strideColRadialNode() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp); }
		inline int strideColRadialZone() const CADET_NOEXCEPT { return static_cast<int>(_disc.radNNodes * strideColRadialNode()); }
		inline int strideColComp() const CADET_NOEXCEPT { return 1; }

		inline int strideParComp() const CADET_NOEXCEPT { return 1; }
		inline int strideParLiquid() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp); }
		inline int strideParBound(int parType) const CADET_NOEXCEPT { return static_cast<int>(_disc.strideBound[parType]); }
		inline int strideParNode(int parType) const CADET_NOEXCEPT { return strideParLiquid() + strideParBound(parType); }
		inline int strideParBlock(int parType) const CADET_NOEXCEPT { return static_cast<int>(_disc.nParPoints[parType]) * strideParNode(parType); }

		// Offsets
		inline int offsetC() const CADET_NOEXCEPT { return _disc.nComp * _disc.radNPoints; }
		inline int offsetCp() const CADET_NOEXCEPT { return offsetC() + _disc.nComp * _disc.nBulkPoints; }
		inline int offsetCp(ParticleTypeIndex pti) const CADET_NOEXCEPT { return offsetCp() + _disc.parTypeOffset[pti.value]; }
		inline int offsetCp(ParticleTypeIndex pti, ParticleIndex pi) const CADET_NOEXCEPT { return offsetCp() + _disc.parTypeOffset[pti.value] + strideParBlock(pti.value) * pi.value; }
		inline int offsetBoundComp(ParticleTypeIndex pti, ComponentIndex comp) const CADET_NOEXCEPT { return _disc.boundOffset[pti.value * _disc.nComp + comp.value]; }

		// Return pointer to first element of state variable in state vector
		template <typename real_t> inline real_t* c(real_t* const data) const { return data + offsetC(); }
		template <typename real_t> inline real_t const* c(real_t const* const data) const { return data + offsetC(); }

		template <typename real_t> inline real_t* cp(real_t* const data) const { return data + offsetCp(); }
		template <typename real_t> inline real_t const* cp(real_t const* const data) const { return data + offsetCp(); }

		template <typename real_t> inline real_t* q(real_t* const data) const { return data + offsetCp() + strideParLiquid(); }
		template <typename real_t> inline real_t const* q(real_t const* const data) const { return data + offsetCp() + strideParLiquid(); }

		// Return specific variable in state vector
		template <typename real_t> inline real_t& c(real_t* const data, unsigned int col, unsigned int rad, unsigned int comp) const { return data[offsetC() + comp + col * strideColAxialNode() + rad * strideColRadialNode()]; }
		template <typename real_t> inline const real_t& c(real_t const* const data, unsigned int col, unsigned int rad, unsigned int comp) const { return data[offsetC() + comp + col * strideColAxialNode() + rad * strideColRadialNode()]; }

	protected:
		const Discretization& _disc;
	};

	class Exporter : public ISolutionExporter
	{
	public:

		Exporter(const Discretization& disc, const ColumnModel2D& model, double const* data) : _disc(disc), _idx(disc), _model(model), _data(data) { }
		Exporter(const Discretization&& disc, const ColumnModel2D& model, double const* data) = delete;

		virtual bool hasParticleFlux() const CADET_NOEXCEPT { return true; }
		virtual bool hasParticleMobilePhase() const CADET_NOEXCEPT { return true; }
		virtual bool hasSolidPhase() const CADET_NOEXCEPT { return _disc.strideBound[_disc.nParType] > 0; }
		virtual bool hasVolume() const CADET_NOEXCEPT { return false; }
		virtual bool isParticleLumped(unsigned int parType) const CADET_NOEXCEPT { return _model._particles[parType]->isParticleLumped(); }
		virtual bool hasPrimaryExtent() const CADET_NOEXCEPT { return true; }
		virtual bool hasSmoothnessIndicator() const CADET_NOEXCEPT { return false; }

		virtual unsigned int numComponents() const CADET_NOEXCEPT { return _disc.nComp; }
		virtual unsigned int numPrimaryCoordinates() const CADET_NOEXCEPT { return _disc.axNPoints; }
		virtual unsigned int numSecondaryCoordinates() const CADET_NOEXCEPT { return _disc.radNPoints; }
		virtual unsigned int numInletPorts() const CADET_NOEXCEPT { return _disc.radNPoints; }
		virtual unsigned int numOutletPorts() const CADET_NOEXCEPT { return _disc.radNPoints; }
		virtual unsigned int numParticleTypes() const CADET_NOEXCEPT { return _disc.nParType; }
		virtual unsigned int numParticleShells(unsigned int parType) const CADET_NOEXCEPT { return _disc.nParPoints[parType]; }
		virtual unsigned int numBoundStates(unsigned int parType) const CADET_NOEXCEPT { return _disc.strideBound[parType]; }
		virtual unsigned int numMobilePhaseDofs() const CADET_NOEXCEPT { return _disc.nComp * _disc.nBulkPoints; }
		virtual unsigned int numParticleMobilePhaseDofs(unsigned int parType) const CADET_NOEXCEPT { return _disc.nComp * _disc.nBulkPoints * _disc.nParPoints[parType]; }
		virtual unsigned int numParticleMobilePhaseDofs() const CADET_NOEXCEPT { return _disc.nComp * _disc.nBulkPoints * std::accumulate(&_disc.nParPoints[0], &_disc.nParPoints[_disc.nParType], 0); }
		virtual unsigned int numSolidPhaseDofs(unsigned int parType) const CADET_NOEXCEPT { return _disc.strideBound[parType] * _disc.nParPoints[parType] * _disc.nBulkPoints; }
		virtual unsigned int numSolidPhaseDofs() const CADET_NOEXCEPT {
			unsigned int nDofPerParType = 0;
			for (unsigned int i = 0; i < _disc.nParType; ++i)
				nDofPerParType += _disc.nParPoints[i] * _disc.strideBound[i];
			return _disc.nBulkPoints * nDofPerParType;
		}
		virtual unsigned int numParticleFluxDofs() const CADET_NOEXCEPT { return _disc.nComp * _disc.axNPoints * _disc.radNPoints * _disc.nParType; }
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

		virtual int writeSmoothnessIndicator(double* indicator) const { return 0; }

		virtual int writePrimaryCoordinates(double* coords) const
		{
			_model._convDispOp.writeAxialCoordinates(coords);
			return _disc.axNPoints;
		}
		virtual int writeSecondaryCoordinates(double* coords) const
		{
			_model._convDispOp.writeRadialCoordinates(coords);
			return _disc.radNPoints;
		}
		virtual int writeParticleCoordinates(unsigned int parType, double* coords) const
		{
			return _model._particles[parType]->writeParticleCoordinates(coords);
		}

	protected:
		const Discretization& _disc;
		const Indexer _idx;
		const ColumnModel2D& _model;
		double const* const _data;
	};

	typedef Eigen::Triplet<double> T;

	/**
	* @brief sets the sparsity pattern of the convection dispersion Jacobian
	*/
	void setJacobianPattern(Eigen::SparseMatrix<double, Eigen::RowMajor>& mat, const int secIdx, bool hasBulkReaction)
	{
		Indexer idxr(_disc);
		std::vector<T> tripletList;

		int nBulkEntries = _convDispOp.nJacEntries();
		if (hasBulkReaction)
			nBulkEntries += _disc.nBulkPoints * _disc.nComp * _disc.nComp; // add nComp entries for every component at each discrete bulk point

		// particle
		int nParEntries = 0;
		for (int type = 0; type < _disc.nParType; type++)
		{
			nParEntries += _disc.nBulkPoints * _particles[type]->jacobianNNZperParticle();
		}

		tripletList.reserve(nBulkEntries + nParEntries);

		_convDispOp.convDispJacPattern(tripletList);

		// bulk reaction jacobian
		if (hasBulkReaction)
		{
			for (unsigned int colNode = 0; colNode < _disc.nBulkPoints; colNode++)
			{
				for (unsigned int comp = 0; comp < _disc.nComp; comp++)
				{
					for (unsigned int toComp = 0; toComp < _disc.nComp; toComp++)
					{
						tripletList.push_back(T(colNode * idxr.strideColRadialNode() + comp * idxr.strideColComp(),
							colNode * idxr.strideColRadialNode() + toComp * idxr.strideColComp(),
							0.0));
					}
				}
			}
		}

		for (int parType = 0; parType < _disc.nParType; parType++)
		{
			for (int colNode = 0; colNode < _disc.nBulkPoints; colNode++)
				_particles[parType]->setParJacPattern(tripletList, idxr.offsetCp(ParticleTypeIndex{ static_cast<unsigned int>(parType) }, ParticleIndex{ static_cast<unsigned int>(colNode) }) - idxr.offsetC(), colNode * idxr.strideColRadialNode(), colNode, secIdx);
		}

		mat.setFromTriplets(tripletList.begin(), tripletList.end());
	}
};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_COLUMNMODEL2D_HPP_
