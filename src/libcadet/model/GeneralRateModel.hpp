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

/**
 * @file 
 * Defines the general rate model (GRM).
 */

#ifndef LIBCADET_GENERALRATEMODEL_HPP_
#define LIBCADET_GENERALRATEMODEL_HPP_

#include "model/UnitOperationBase.hpp"
#include "cadet/StrongTypes.hpp"
#include "cadet/SolutionExporter.hpp"
#include "model/parts/ConvectionDispersionOperator.hpp"
#include "AutoDiff.hpp"
#include "linalg/SparseMatrix.hpp"
#include "linalg/BandMatrix.hpp"
#include "linalg/Gmres.hpp"
#include "Memory.hpp"
#include "model/ModelUtils.hpp"
#include "ParameterMultiplexing.hpp"

#include <array>
#include <vector>

#include "Benchmark.hpp"

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
class IParameterDependence;

/**
 * @brief General rate model of liquid column chromatography
 * @details See @cite Guiochon2006, @cite Gu1995, @cite Felinger2004
 * 
 * @f[\begin{align}
	\frac{\partial c_i}{\partial t} &= - u \frac{\partial c_i}{\partial z} + D_{\text{ax},i} \frac{\partial^2 c_i}{\partial z^2} - \frac{1 - \varepsilon_c}{\varepsilon_c} \frac{3}{r_p} j_{f,i} \\
	\frac{\partial c_{p,i}}{\partial t} + \frac{1 - \varepsilon_p}{\varepsilon_p} \frac{\partial q_{i}}{\partial t} &= D_{p,i} \left( \frac{\partial^2 c_{p,i}}{\partial r^2} + \frac{2}{r} \frac{\partial c_{p,i}}{\partial r} \right) + D_{s,i} \frac{1 - \varepsilon_p}{\varepsilon_p} \left( \frac{\partial^2 q_{i}}{\partial r^2} + \frac{2}{r} \frac{\partial q_{i}}{\partial r} \right) \\
	a \frac{\partial q_i}{\partial t} &= f_{\text{iso}}(c_p, q)
\end{align} @f]
@f[ \begin{align}
	j_{f,i} = k_{f,i} \left( c_i - c_{p,i} \left(\cdot, \cdot, r_p\right)\right)
\end{align} @f]
 * Danckwerts boundary conditions (see @cite Danckwerts1953)
@f[ \begin{align}
u c_{\text{in},i}(t) &= u c_i(t,0) - D_{\text{ax},i} \frac{\partial c_i}{\partial z}(t,0) \\
\frac{\partial c_i}{\partial z}(t,L) &= 0 \\
\varepsilon_p D_{p,i} \frac{\partial c_{p,i}}{\partial r}(\cdot, \cdot, r_p) + (1-\varepsilon_p) D_{s,i} \frac{\partial q_{i}}{\partial r}(\cdot, \cdot, r_p) &= j_{f,i} \\
\frac{\partial c_{p,i}}{\partial r}(\cdot, \cdot, 0) &= 0
\end{align} @f]
 * Methods are described in @cite VonLieres2010a (WENO, linear solver), @cite Puttmann2013 @cite Puttmann2016 (forward sensitivities, AD, band compression)
 */
class GeneralRateModel : public UnitOperationBase
{
public:

	GeneralRateModel(UnitOpIdx unitOpIdx);
	virtual ~GeneralRateModel() CADET_NOEXCEPT;

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

	static const char* identifier() { return "GENERAL_RATE_MODEL"; }
	virtual const char* unitOperationName() const CADET_NOEXCEPT { return "GENERAL_RATE_MODEL"; }

	virtual bool configureModelDiscretization(IParameterProvider& paramProvider, IConfigHelper& helper);
	virtual bool configure(IParameterProvider& paramProvider);
	virtual void notifyDiscontinuousSectionTransition(double t, unsigned int secIdx, const ConstSimulationState& simState, const AdJacobianParams& adJac);

	virtual void useAnalyticJacobian(const bool analyticJac);

	virtual void reportSolution(ISolutionRecorder& recorder, double const* const solution) const;
	virtual void reportSolutionStructure(ISolutionRecorder& recorder) const;

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

	virtual std::unordered_map<ParameterId, double> getAllParameterValues() const;
	virtual double getParameterDouble(const ParameterId& pId) const;
	virtual bool hasParameter(const ParameterId& pId) const;

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
			_timerFactorizePar.totalElapsedTime(),
			_timerMatVec.totalElapsedTime(),
			_timerGmres.totalElapsedTime(),
			static_cast<double>(_gmres.numIterations())
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
			"Gmres",
			"NumGMRESIter"
		};
		return desc;
	}
#endif

protected:

	class Indexer;

	int residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem, bool updateJacobian, bool paramSensitivity);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	int residualImpl(double t, unsigned int secIdx, StateType const* const y, double const* const yDot, ResidualType* const res, util::ThreadLocalStorage& threadLocalMem);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	int residualBulk(double t, unsigned int secIdx, StateType const* y, double const* yDot, ResidualType* res, util::ThreadLocalStorage& threadLocalMem);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	int residualParticle(double t, unsigned int parType, unsigned int colCell, unsigned int secIdx, StateType const* y, double const* yDot, ResidualType* res, util::ThreadLocalStorage& threadLocalMem);

	template <typename StateType, typename ResidualType, typename ParamType>
	int residualFlux(double t, unsigned int secIdx, StateType const* y, double const* yDot, ResidualType* res);

	void assembleOffdiagJac(double t, unsigned int secIdx, double const* vecStateY);
	void assembleOffdiagJacFluxParticle(double t, unsigned int secIdx, double const* vecStateY);
	void extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset);

	int schurComplementMatrixVector(double const* x, double* z) const;
	void assembleDiscretizedJacobianParticleBlock(unsigned int parType, unsigned int pblk, double alpha, const Indexer& idxr);
	
	void setEquidistantRadialDisc(unsigned int parType);
	void setEquivolumeRadialDisc(unsigned int parType);
	void setUserdefinedRadialDisc(unsigned int parType);
	void updateRadialDisc();

	void addTimeDerivativeToJacobianParticleShell(linalg::FactorizableBandMatrix::RowIterator& jac, const Indexer& idxr, double alpha, unsigned int parType);
	void solveForFluxes(double* const vecState, const Indexer& idxr) const;
	
	unsigned int numAdDirsForJacobian() const CADET_NOEXCEPT;

	int multiplexInitialConditions(const cadet::ParameterId& pId, unsigned int adDirection, double adValue);
	int multiplexInitialConditions(const cadet::ParameterId& pId, double val, bool checkSens);

	void clearParDepSurfDiffusion();

	parts::cell::CellParameters makeCellResidualParams(unsigned int parType, int const* qsReaction) const;

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	void checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const;
#endif

	struct Discretization
	{
		unsigned int nComp; //!< Number of components
		unsigned int nCol; //!< Number of column cells
		unsigned int nParType; //!< Number of particle types
		unsigned int* nParCell; //!< Array with number of radial cells in each particle type
		unsigned int* nParCellsBeforeType; //!< Array with total number of radial cells before a particle type (cumulative sum of nParCell), additional last element contains total number of particle shells
		unsigned int* parTypeOffset; //!< Array with offsets (in particle block) to particle type, additional last element contains total number of particle DOFs
		unsigned int* nBound; //!< Array with number of bound states for each component and particle type (particle type major ordering)
		unsigned int* boundOffset; //!< Array with offset to the first bound state of each component in the solid phase (particle type major ordering)
		unsigned int* strideBound; //!< Total number of bound states for each particle type, additional last element contains total number of bound states for all types
		unsigned int* nBoundBeforeType; //!< Array with number of bound states before a particle type (cumulative sum of strideBound)
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
	std::vector<bool> _hasSurfaceDiffusion; //!< Determines whether surface diffusion is present in each particle type
//	IExternalFunction* _extFun; //!< External function (owned by library user)

	parts::ConvectionDispersionOperator _convDispOp; //!< Convection dispersion operator for interstitial volume transport
	IDynamicReactionModel* _dynReactionBulk; //!< Dynamic reactions in the bulk volume

	linalg::BandMatrix* _jacP; //!< Particle jacobian diagonal blocks (all of them)
	linalg::FactorizableBandMatrix* _jacPdisc; //!< Particle jacobian diagonal blocks (all of them) with time derivatives from BDF method

	linalg::DoubleSparseMatrix _jacCF; //!< Jacobian block connecting interstitial states and fluxes (interstitial transport equation)
	linalg::DoubleSparseMatrix _jacFC; //!< Jacobian block connecting fluxes and interstitial states (flux equation)
	linalg::DoubleSparseMatrix* _jacPF; //!< Jacobian blocks connecting particle states and fluxes (particle transport boundary condition)
	linalg::DoubleSparseMatrix* _jacFP; //!< Jacobian blocks connecting fluxes and particle states (flux equation)

	linalg::DoubleSparseMatrix _jacInlet; //!< Jacobian inlet DOF block matrix connects inlet DOFs to first bulk cells

	active _colPorosity; //!< Column porosity (external porosity) \f$ \varepsilon_c \f$
	std::vector<active> _parRadius; //!< Particle radius \f$ r_p \f$
	bool _singleParRadius;
	std::vector<active> _parCoreRadius; //!< Particle core radius \f$ r_c \f$
	bool _singleParCoreRadius;
	std::vector<active> _parPorosity; //!< Particle porosity (internal porosity) \f$ \varepsilon_p \f$
	bool _singleParPorosity;
	std::vector<active> _parTypeVolFrac; //!< Volume fraction of each particle type
	std::vector<ParticleDiscretizationMode> _parDiscType; //!< Particle discretization mode
	std::vector<double> _parDiscVector; //!< Particle discretization shell edges
	std::vector<double> _parGeomSurfToVol; //!< Particle surface to volume ratio factor (i.e., 3.0 for spherical, 2.0 for cylindrical, 1.0 for hexahedral)

	// Vectorial parameters
	std::vector<active> _filmDiffusion; //!< Film diffusion coefficient \f$ k_f \f$
	MultiplexMode _filmDiffusionMode;
	std::vector<active> _parDiffusion; //!< Particle diffusion coefficient \f$ D_p \f$
	MultiplexMode _parDiffusionMode;
	std::vector<active> _parSurfDiffusion; //!< Particle surface diffusion coefficient \f$ D_s \f$
	MultiplexMode _parSurfDiffusionMode;
	std::vector<active> _poreAccessFactor; //!< Pore accessibility factor \f$ F_{\text{acc}} \f$
	MultiplexMode _poreAccessFactorMode;
	std::vector<IParameterDependence*> _parDepSurfDiffusion; //!< Parameter dependencies for particle surface diffusion
	bool _singleParDepSurfDiffusion; //!< Determines whether a single parameter dependence for particle surface diffusion is used
	bool _hasParDepSurfDiffusion; //!< Determines whether particle surface diffusion parameter dependencies are present

	bool _axiallyConstantParTypeVolFrac; //!< Determines whether particle type volume fraction is homogeneous across axial coordinate
	bool _analyticJac; //!< Determines whether AD or analytic Jacobians are used
	unsigned int _jacobianAdDirs; //!< Number of AD seed vectors required for Jacobian computation

	std::vector<active> _parCellSize; //!< Particle cell / shell size
	std::vector<active> _parCenterRadius; //!< Particle cell-centered position for each particle cell
	std::vector<active> _parOuterSurfAreaPerVolume; //!< Particle shell outer sphere surface to volume ratio
	std::vector<active> _parInnerSurfAreaPerVolume; //!< Particle shell inner sphere surface to volume ratio

	ArrayPool _discParFlux; //!< Storage for discretized @f$ k_f @f$ value

	bool _factorizeJacobian; //!< Determines whether the Jacobian needs to be factorized
	double* _tempState; //!< Temporary storage with the size of the state vector or larger if binding models require it
	linalg::Gmres _gmres; //!< GMRES algorithm for the Schur-complement in linearSolve()
	double _schurSafety; //!< Safety factor for Schur-complement solution
	int _colParBoundaryOrder; //!< Order of the bulk-particle boundary discretization

	std::vector<active> _initC; //!< Liquid bulk phase initial conditions
	std::vector<active> _initCp; //!< Liquid particle phase initial conditions
	std::vector<active> _initQ; //!< Solid phase initial conditions
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
	BENCH_TIMER(_timerGmres)

	// Wrapper for calling the corresponding function in GeneralRateModel class
	friend int schurComplementMultiplierGRM(void* userData, double const* x, double* z);

	class Indexer
	{
	public:
		Indexer(const Discretization& disc) : _disc(disc) { }

		// Strides
		inline const int strideColCell() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp); }
		inline const int strideColComp() const CADET_NOEXCEPT { return 1; }

		inline const int strideParComp() const CADET_NOEXCEPT { return 1; }
		inline const int strideParLiquid() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp); }
		inline const int strideParBound(int parType) const CADET_NOEXCEPT { return static_cast<int>(_disc.strideBound[parType]); }
		inline const int strideParShell(int parType) const CADET_NOEXCEPT { return strideParLiquid() + strideParBound(parType); }
		inline const int strideParBlock(int parType) const CADET_NOEXCEPT { return static_cast<int>(_disc.nParCell[parType]) * strideParShell(parType); }

		inline const int strideFluxCell() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp) * static_cast<int>(_disc.nParType); }
		inline const int strideFluxParType() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp); }
		inline const int strideFluxComp() const CADET_NOEXCEPT { return 1; }

		// Offsets
		inline const int offsetC() const CADET_NOEXCEPT { return _disc.nComp; }
		inline const int offsetCp() const CADET_NOEXCEPT { return _disc.nComp * _disc.nCol + offsetC(); }
		inline const int offsetCp(ParticleTypeIndex pti) const CADET_NOEXCEPT { return offsetCp() + _disc.parTypeOffset[pti.value]; }
		inline const int offsetCp(ParticleTypeIndex pti, ParticleIndex pi) const CADET_NOEXCEPT { return offsetCp() + _disc.parTypeOffset[pti.value] + strideParBlock(pti.value) * pi.value; }
		inline const int offsetJf() const CADET_NOEXCEPT { return offsetCp() + _disc.parTypeOffset[_disc.nParType]; }
		inline const int offsetJf(ParticleTypeIndex pti) const CADET_NOEXCEPT { return offsetJf() + pti.value * _disc.nCol * _disc.nComp; }
		inline const int offsetBoundComp(ParticleTypeIndex pti, ComponentIndex comp) const CADET_NOEXCEPT { return _disc.boundOffset[pti.value * _disc.nComp + comp.value]; }

		// Return pointer to first element of state variable in state vector
		template <typename real_t> inline real_t* c(real_t* const data) const { return data + offsetC(); }
		template <typename real_t> inline real_t const* c(real_t const* const data) const { return data + offsetC(); }

		template <typename real_t> inline real_t* cp(real_t* const data) const { return data + offsetCp(); }
		template <typename real_t> inline real_t const* cp(real_t const* const data) const { return data + offsetCp(); }

		template <typename real_t> inline real_t* q(real_t* const data) const { return data + offsetCp() + strideParLiquid(); }
		template <typename real_t> inline real_t const* q(real_t const* const data) const { return data + offsetCp() + strideParLiquid(); }

		template <typename real_t> inline real_t* jf(real_t* const data) const { return data + offsetJf(); }
		template <typename real_t> inline real_t const* jf(real_t const* const data) const { return data + offsetJf(); }

		// Return specific variable in state vector
		template <typename real_t> inline real_t& c(real_t* const data, unsigned int col, unsigned int comp) const { return data[offsetC() + comp + col * strideColCell()]; }
		template <typename real_t> inline const real_t& c(real_t const* const data, unsigned int col, unsigned int comp) const { return data[offsetC() + comp + col * strideColCell()]; }

	protected:
		const Discretization& _disc;
	};

	class Exporter : public ISolutionExporter
	{
	public:

		Exporter(const Discretization& disc, const GeneralRateModel& model, double const* data) : _disc(disc), _idx(disc), _model(model), _data(data) { }
		Exporter(const Discretization&& disc, const GeneralRateModel& model, double const* data) = delete;

		virtual bool hasParticleFlux() const CADET_NOEXCEPT { return true; }
		virtual bool hasParticleMobilePhase() const CADET_NOEXCEPT { return true; }
		virtual bool hasSolidPhase() const CADET_NOEXCEPT { return _disc.strideBound[_disc.nParType] > 0; }
		virtual bool hasVolume() const CADET_NOEXCEPT { return false; }

		virtual unsigned int numComponents() const CADET_NOEXCEPT { return _disc.nComp; }
		virtual unsigned int numAxialCells() const CADET_NOEXCEPT { return _disc.nCol; }
		virtual unsigned int numRadialCells() const CADET_NOEXCEPT { return 0; }
		virtual unsigned int numInletPorts() const CADET_NOEXCEPT { return 1; }
		virtual unsigned int numOutletPorts() const CADET_NOEXCEPT { return 1; }
		virtual unsigned int numParticleTypes() const CADET_NOEXCEPT { return _disc.nParType; }
		virtual unsigned int numParticleShells(unsigned int parType) const CADET_NOEXCEPT { return _disc.nParCell[parType]; }
		virtual unsigned int numBoundStates(unsigned int parType) const CADET_NOEXCEPT { return _disc.strideBound[parType]; }
		virtual unsigned int numBulkDofs() const CADET_NOEXCEPT { return _disc.nComp * _disc.nCol; }
		virtual unsigned int numParticleMobilePhaseDofs(unsigned int parType) const CADET_NOEXCEPT { return _disc.nCol * _disc.nParCell[parType] * _disc.nComp; }
		virtual unsigned int numSolidPhaseDofs(unsigned int parType) const CADET_NOEXCEPT { return _disc.nCol * _disc.nParCell[parType] * _disc.strideBound[parType]; }
		virtual unsigned int numFluxDofs() const CADET_NOEXCEPT { return _disc.nComp * _disc.nCol * _disc.nParType; }
		virtual unsigned int numVolumeDofs() const CADET_NOEXCEPT { return 0; }

		virtual double const* concentration() const { return _idx.c(_data); }
		virtual double const* flux() const { return _idx.jf(_data); }
		virtual double const* particleMobilePhase(unsigned int parType) const { return _data + _idx.offsetCp(ParticleTypeIndex{parType}); }
		virtual double const* solidPhase(unsigned int parType) const { return _data + _idx.offsetCp(ParticleTypeIndex{parType}) + _idx.strideParLiquid(); }
		virtual double const* volume() const { return nullptr; }
		virtual double const* inlet(unsigned int port, unsigned int& stride) const
		{
			stride = _idx.strideColComp();
			return _data;
		}
		virtual double const* outlet(unsigned int port, unsigned int& stride) const
		{
			stride = _idx.strideColComp();
			if (_model._convDispOp.currentVelocity() >= 0)
				return &_idx.c(_data, _disc.nCol - 1, 0);
			else
				return &_idx.c(_data, 0, 0);
		}

		virtual StateOrdering const* concentrationOrdering(unsigned int& len) const
		{
			len = _concentrationOrdering.size();
			return _concentrationOrdering.data();
		}

		virtual StateOrdering const* fluxOrdering(unsigned int& len) const
		{
			len = _fluxOrdering.size();
			return _fluxOrdering.data();
		}

		virtual StateOrdering const* mobilePhaseOrdering(unsigned int& len) const
		{
			len = _particleOrdering.size();
			return _particleOrdering.data();
		}

		virtual StateOrdering const* solidPhaseOrdering(unsigned int& len) const
		{
			len = _solidOrdering.size();
			return _solidOrdering.data();
		}

		virtual unsigned int bulkMobilePhaseStride() const { return _idx.strideColCell(); }
		virtual unsigned int particleMobilePhaseStride(unsigned int parType) const { return _idx.strideParShell(parType); }
		virtual unsigned int solidPhaseStride(unsigned int parType) const { return _idx.strideParShell(parType); }

		virtual void axialCoordinates(double* coords) const
		{
			const double h = static_cast<double>(_model._convDispOp.columnLength()) / static_cast<double>(_disc.nCol);
			for (unsigned int i = 0; i < _disc.nCol; ++i)
				coords[i] = (i + 0.5) * h;
		}
		virtual void radialCoordinates(double* coords) const { }
		virtual void particleCoordinates(unsigned int parType, double* coords) const
		{
			active const* const pcr = _model._parCenterRadius.data() + _disc.nParCellsBeforeType[parType];
			for (unsigned int i = 0; i < _disc.nParCell[parType]; ++i)
				coords[i] = static_cast<double>(pcr[i]);
		}

	protected:
		const Discretization& _disc;
		const Indexer _idx;
		const GeneralRateModel& _model;
		double const* const _data;

		const std::array<StateOrdering, 2> _concentrationOrdering = { { StateOrdering::AxialCell, StateOrdering::Component } };
		const std::array<StateOrdering, 4> _particleOrdering = { { StateOrdering::ParticleType, StateOrdering::AxialCell, StateOrdering::ParticleShell, StateOrdering::Component } };
		const std::array<StateOrdering, 5> _solidOrdering = { { StateOrdering::ParticleType, StateOrdering::AxialCell, StateOrdering::ParticleShell, StateOrdering::Component, StateOrdering::BoundState } };
		const std::array<StateOrdering, 3> _fluxOrdering = { { StateOrdering::ParticleType, StateOrdering::AxialCell, StateOrdering::Component } };
	};
};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_GENERALRATEMODEL_HPP_
