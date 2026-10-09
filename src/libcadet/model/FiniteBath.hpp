// SPDX-License-Identifier: AGPL-3.0-or-later
// =================================================================================
//  CADET
//
//  Copyright © 2008-present: The CADET-Core Authors
//            Please see the AUTHORS.md file.
//
//  All rights reserved. This program and the accompanying materials
//  are made available under the terms of the GNU Affero General Public
//  License v3.0 (or, at your option, any later version).
// =================================================================================

/**
 * @file
 * Defines the finite bath unit operation
 */

#ifndef LIBCADET_FINITEBATH_HPP_
#define LIBCADET_FINITEBATH_HPP_

#include "model/UnitOperationBase.hpp"
#include "model/particle/ParticleModel.hpp"
#include "model/reaction/ReactionSystem.hpp"
#include "cadet/StrongTypes.hpp"
#include "cadet/SolutionExporter.hpp"
#include "AutoDiff.hpp"
#include "linalg/BandedEigenSparseRowIterator.hpp"
#include "linalg/EigenSolverWrapper.hpp"
#include "Memory.hpp"
#include "model/ModelUtils.hpp"
#include "ParameterMultiplexing.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
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
class IParameterStateDependence;

/**
 * @brief Finite bath model
 * @details A well mixed vessel of constant volume that holds a distribution of particle types,
 *          that is, the 0D bulk model combined with the 0D or 1D particle models. It describes
 *          batch uptake experiments and, more generally, any stirred tank in which the transport
 *          into the particles is limited by film diffusion.
 *
 *          Writing the bulk mass balance per unit bulk liquid volume @f$ V^\ell @f$ gives
@f[ \begin{align}
	\frac{\mathrm{d} c^\ell_i}{\mathrm{d} t} &= \frac{1}{V^\ell} \left( F_\text{in} c^\ell_{\text{in},i} - F_\text{out} c^\ell_i \right)
	- \frac{1 - \varepsilon_b}{\varepsilon_b} \sum_{j=1}^{N_\text{par}} d_j a_j k_{f,j,i} \left( c^\ell_i - c^p_{j,i}\left(\cdot, r_{p,j}\right) \right)
	+ f_{\text{react},i}^\ell\left( c^\ell \right),
\end{align} @f]
 *          where @f$ \varepsilon_b = V^\ell / (V^\ell + V^p) @f$ is the bulk porosity, @f$ d_j @f$ the volume
 *          fraction of particle type @f$ j @f$ and @f$ a_j @f$ its surface to volume ratio. The particle
 *          equations are provided by the particle models, see @c IParticleModel.
 *
 *          Since the volume is constant, the flow rates have to cancel, @f$ F_\text{in} = F_\text{out} @f$.
 *          Particles in rapid equilibrium with the bulk liquid are not supported, those are described by
 *          the CSTR, see @c CSTRModel.
 */
class FiniteBath : public UnitOperationBase
{
public:

	FiniteBath(UnitOpIdx unitOpIdx);
	virtual ~FiniteBath() CADET_NOEXCEPT;

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

	static const char* identifier() { return "FINITE_BATH"; }
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

	virtual unsigned int localOutletComponentIndex(unsigned int port) const CADET_NOEXCEPT { return _disc.nComp; }
	virtual unsigned int localOutletComponentStride(unsigned int port) const CADET_NOEXCEPT { return 1; }
	virtual unsigned int localInletComponentIndex(unsigned int port) const CADET_NOEXCEPT { return 0; }
	virtual unsigned int localInletComponentStride(unsigned int port) const CADET_NOEXCEPT { return 1; }

	virtual void setExternalFunctions(IExternalFunction** extFuns, unsigned int size);
	virtual void setSectionTimes(double const* secTimes, bool const* secContinuity, unsigned int nSections) { }

	virtual void expandErrorTol(double const* errorSpec, unsigned int errorSpecSize, double* expandOut) { }

	virtual void multiplyWithJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* yS, double alpha, double beta, double* ret);
	virtual void multiplyWithDerivativeJacobian(const SimulationTime& simTime, const ConstSimulationState& simState, double const* sDot, double* ret);

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
			_timerResidualSens.totalElapsedTime(),
			_timerConsistentInit.totalElapsedTime(),
			_timerLinearSolve.totalElapsedTime(),
			_timerFactorize.totalElapsedTime()
		});
	}

	virtual char const* const* benchmarkDescriptions() const
	{
		static const char* const desc[] = {
			"DOFs",
			"Residual",
			"ResidualSens",
			"ConsistentInit",
			"LinearSolve",
			"Factorize"
		};
		return desc;
	}
#endif

protected:

	class Indexer;

	int residual(const SimulationTime& simTime, const ConstSimulationState& simState, double* const res, const AdJacobianParams& adJac, util::ThreadLocalStorage& threadLocalMem, bool updateJacobian, bool paramSensitivity);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes = true>
	int residualImpl(double t, unsigned int secIdx, StateType const* const y, double const* const yDot, ResidualType* const res, util::ThreadLocalStorage& threadLocalMem);

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes = true>
	int residualBulk(double t, unsigned int secIdx, StateType const* y, double const* yDot, ResidualType* res, util::ThreadLocalStorage& threadLocalMem);

	void extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset);

	void assembleDiscretizedGlobalJacobian(double alpha, Indexer idxr);

	void addTimeDerivativeToJacobianParticleShell(linalg::BandedEigenSparseRowIterator& jac, const Indexer& idxr, double alpha, unsigned int parType);

	unsigned int numAdDirsForJacobian() const CADET_NOEXCEPT;

	int multiplexInitialConditions(const cadet::ParameterId& pId, unsigned int adDirection, double adValue);
	int multiplexInitialConditions(const cadet::ParameterId& pId, double val, bool checkSens);
	void consistentInitialBulkLiquidEquilibrium(const SimulationTime& simTime, double* const vecStateY, double errorTol, util::ThreadLocalStorage& threadLocalMem);
	void consistentInitialBindingEquilibrium(const SimulationTime& simTime, double* const vecStateY, const AdJacobianParams& adJac, double errorTol, util::ThreadLocalStorage& threadLocalMem, unsigned int type);
	void consistentInitialBulkTimeDerivative(const SimulationTime& simTime, double const* vecStateY, double* const vecStateYdot, util::ThreadLocalStorage& threadLocalMem);
	void consistentInitialBindingTimeDerivative(const SimulationTime& simTime, double* const vecStateYdot, util::ThreadLocalStorage& threadLocalMem, unsigned int type);
	void initializeSensitivityBulkStates(double* const sensY, std::size_t param) const;
	void initializeSensitivityParticleStates(double* const sensY, std::size_t param, unsigned int type) const;
	void consistentInitialSensitivityBindingEquilibrium(double* const sensY, double const* const sensYdot, util::ThreadLocalStorage& threadLocalMem, unsigned int type);
	void consistentInitialSensitivityBulkTimeDerivative();
	void consistentInitialSensitivityBindingTimeDerivative(double* const sensYdot, unsigned int type);
	void leanConsistentInitialBindingEquilibrium(const SimulationTime& simTime, double* const vecStateY, util::ThreadLocalStorage& threadLocalMem, unsigned int type);
	void solveBulkTimeDerivativeSystem(double* const rhs);

	parts::cell::CellParameters makeCellResidualParams(unsigned int parType, int const* qsReaction) const;

	void resizeConservedMoietyJacobianBuffer();

	/**
	 * @brief Returns the derivative of a bulk residual with respect to the corresponding inlet concentration
	 */
	inline double inletJacobian() const CADET_NOEXCEPT { return -static_cast<double>(_flowRateIn) / static_cast<double>(_liquidVolume); }

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	void checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const;
#endif

	struct Discretization
	{
		unsigned int nComp; //!< Number of components
		unsigned int nParType; //!< Number of particle types
		unsigned int* nParPoints; //!< Array with number of discrete points for each particle type
		unsigned int* parTypeOffset; //!< Array with offsets (in particle block) to particle type, additional last element contains total number of particle DOFs
		unsigned int* nBound; //!< Array with number of bound states for each component and particle type (particle type major ordering)
		unsigned int* boundOffset; //!< Array with offset to the first bound state of each component in the solid phase (particle type major ordering)
		unsigned int* strideBound; //!< Total number of bound states for each particle type, additional last element contains total number of bound states for all types
		unsigned int* nBoundBeforeType; //!< Array with number of bound states before a particle type (cumulative sum of strideBound)

		int curSection; //!< current time section index

		Discretization() : nComp(0), nParType(0), nParPoints(nullptr), parTypeOffset(nullptr), nBound(nullptr),
			boundOffset(nullptr), strideBound(nullptr), nBoundBeforeType(nullptr), curSection(-1)
		{
		}

		~Discretization() // make sure this memory is freed correctly
		{
			delete[] nParPoints;
			delete[] parTypeOffset;
			delete[] nBound;
			delete[] boundOffset;
			delete[] strideBound;
			delete[] nBoundBeforeType;
		}
	};

	Discretization _disc; //!< Discretization info

	std::vector<IParticleModel*> _particles; //!< Particle models of the particle type distribution

	ReactionSystem _reaction; //!< Reaction system for bulk phase

	std::vector<Eigen::Triplet<double>> _cMJacobianEntries;

	cadet::linalg::EigenSolverBase* _linearSolver; //!< Linear solver

	Eigen::SparseMatrix<double, Eigen::RowMajor> _globalJac; //!< static part of global Jacobian
	Eigen::SparseMatrix<double, Eigen::RowMajor> _globalJacDisc; //!< global Jacobian with time derivative from BDF method

	active _liquidVolume; //!< Constant volume of the bulk liquid \f$ V^\ell \f$
	active _bulkPorosity; //!< Bulk porosity \f$ \varepsilon_b \f$, i.e. the ratio of bulk liquid volume to total volume
	std::vector<active> _parTypeVolFrac; //!< Volume fraction of each particle type

	active _flowRateIn; //!< Volumetric flow rate of incoming stream
	active _flowRateOut; //!< Volumetric flow rate of drawn outgoing stream

	bool _analyticJac; //!< Determines whether AD or analytic Jacobians are used
	unsigned int _jacobianAdDirs; //!< Number of AD seed vectors required for Jacobian computation

	bool _factorizeJacobian; //!< Determines whether the Jacobian needs to be factorized
	double* _tempState; //!< Temporary storage with the size of the state vector or larger if binding models require it

	std::vector<active> _initC; //!< Liquid bulk phase initial conditions
	std::vector<active> _initCp; //!< Liquid particle phase initial conditions
	std::vector<active> _initCs; //!< Solid phase initial conditions
	std::vector<double> _initState; //!< Initial conditions for state vector if given
	std::vector<double> _initStateDot; //!< Initial conditions for time derivative

	BENCH_TIMER(_timerResidual)
	BENCH_TIMER(_timerResidualSens)
	BENCH_TIMER(_timerConsistentInit)
	BENCH_TIMER(_timerLinearSolve)
	BENCH_TIMER(_timerFactorize)

	class Indexer
	{
	public:
		Indexer(const Discretization& disc) : _disc(disc) { }

		// Strides
		inline int strideColComp() const CADET_NOEXCEPT { return 1; }

		inline int strideParComp() const CADET_NOEXCEPT { return 1; }
		inline int strideParLiquid() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp); }
		inline int strideParBound(int parType) const CADET_NOEXCEPT { return static_cast<int>(_disc.strideBound[parType]); }
		inline int strideParNode(int parType) const CADET_NOEXCEPT { return strideParLiquid() + strideParBound(parType); }
		inline int strideParBlock(int parType) const CADET_NOEXCEPT { return static_cast<int>(_disc.nParPoints[parType]) * strideParNode(parType); }

		// Offsets
		inline int offsetC() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp); }
		inline int offsetCp() const CADET_NOEXCEPT { return offsetC() + static_cast<int>(_disc.nComp); }
		inline int offsetCp(ParticleTypeIndex pti) const CADET_NOEXCEPT { return offsetCp() + _disc.parTypeOffset[pti.value]; }
		inline int offsetBoundComp(ParticleTypeIndex pti, ComponentIndex comp) const CADET_NOEXCEPT { return _disc.boundOffset[pti.value * _disc.nComp + comp.value]; }

		// Return pointer to first element of state variable in state vector
		template <typename real_t> inline real_t* c(real_t* const data) const { return data + offsetC(); }
		template <typename real_t> inline real_t const* c(real_t const* const data) const { return data + offsetC(); }

		template <typename real_t> inline real_t* cp(real_t* const data) const { return data + offsetCp(); }
		template <typename real_t> inline real_t const* cp(real_t const* const data) const { return data + offsetCp(); }

	protected:
		const Discretization& _disc;
	};

	class Exporter : public ISolutionExporter
	{
	public:

		Exporter(const Discretization& disc, const FiniteBath& model, double const* data) : _disc(disc), _idx(disc), _model(model), _data(data) { }
		Exporter(const Discretization&& disc, const FiniteBath& model, double const* data) = delete;

		virtual bool hasParticleFlux() const CADET_NOEXCEPT { return false; }
		virtual bool hasParticleMobilePhase() const CADET_NOEXCEPT { return _disc.nParType > 0; }
		virtual bool hasSolidPhase() const CADET_NOEXCEPT { return _disc.strideBound[_disc.nParType] > 0; }
		virtual bool hasVolume() const CADET_NOEXCEPT { return false; }
		virtual bool isParticleLumped(unsigned int parType) const CADET_NOEXCEPT { return _model._particles[parType]->isParticleLumped(); }
		virtual bool hasPrimaryExtent() const CADET_NOEXCEPT { return false; }
		virtual bool discHasSmoothnessIndicator() const CADET_NOEXCEPT { return false; }

		virtual unsigned int numComponents() const CADET_NOEXCEPT { return _disc.nComp; }
		virtual unsigned int numPrimaryCoordinates() const CADET_NOEXCEPT { return 1; }
		virtual unsigned int numSecondaryCoordinates() const CADET_NOEXCEPT { return 0; }
		virtual unsigned int numInletPorts() const CADET_NOEXCEPT { return 1; }
		virtual unsigned int numOutletPorts() const CADET_NOEXCEPT { return 1; }
		virtual unsigned int numParticleTypes() const CADET_NOEXCEPT { return _disc.nParType; }
		virtual unsigned int numParticleShells(unsigned int parType) const CADET_NOEXCEPT { return _disc.nParPoints[parType]; }
		virtual unsigned int numBoundStates(unsigned int parType) const CADET_NOEXCEPT { return _disc.strideBound[parType]; }
		virtual unsigned int numBoundStates(unsigned int parType, unsigned int comp) const CADET_NOEXCEPT { return _disc.nBound[parType * _disc.nComp + comp]; }
		virtual unsigned int numMobilePhaseDofs() const CADET_NOEXCEPT { return _disc.nComp; }
		virtual unsigned int numParticleMobilePhaseDofs() const CADET_NOEXCEPT
		{
			unsigned int nDof = 0;
			for (unsigned int i = 0; i < _disc.nParType; ++i)
				nDof += _disc.nParPoints[i] * _disc.nComp;
			return nDof;
		}
		virtual unsigned int numParticleMobilePhaseDofs(unsigned int parType) const CADET_NOEXCEPT { return _disc.nParPoints[parType] * _disc.nComp; }
		virtual unsigned int numSolidPhaseDofs() const CADET_NOEXCEPT
		{
			unsigned int nDof = 0;
			for (unsigned int i = 0; i < _disc.nParType; ++i)
				nDof += _disc.nParPoints[i] * _disc.strideBound[i];
			return nDof;
		}
		virtual unsigned int numSolidPhaseDofs(unsigned int parType) const CADET_NOEXCEPT { return _disc.nParPoints[parType] * _disc.strideBound[parType]; }
		virtual unsigned int numParticleFluxDofs() const CADET_NOEXCEPT { return 0; }
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
		virtual int writeSmoothnessIndicator(double* buffer) const { return 0; }

		virtual int writePrimaryCoordinates(double* coords) const
		{
			coords[0] = 0.0;
			return 1;
		}

		virtual int writeSecondaryCoordinates(double* coords) const { return 0; }

		virtual int writeParticleCoordinates(unsigned int parType, double* coords) const
		{
			return _model._particles[parType]->writeParticleCoordinates(coords);
		}

	protected:
		const Discretization& _disc;
		const Indexer _idx;
		const FiniteBath& _model;
		double const* const _data;
	};

	typedef Eigen::Triplet<double> T;

	void setJacobianPattern(Eigen::SparseMatrix<double, Eigen::RowMajor>& globalJ, unsigned int secIdx, bool hasBulkReaction)
	{
		Indexer idxr(_disc);
		std::vector<T> tripletList;

		// The bulk block is dense if liquid reactions couple the components, otherwise the film diffusion
		// only contributes the main diagonal
		const int bulkEntries = hasBulkReaction ? _disc.nComp * _disc.nComp : _disc.nComp;

		int particleEntries = 0;
		for (unsigned int type = 0; type < _disc.nParType; ++type)
		{
			// every bound state might depend on every bound and liquid state
			const int isothermNNZ = idxr.strideParNode(type) * _disc.nParPoints[type] * _disc.strideBound[type];
			particleEntries += _particles[type]->jacobianNNZperParticle() + isothermNNZ;
		}

		tripletList.reserve(bulkEntries + particleEntries);

		if (hasBulkReaction)
		{
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
			{
				for (unsigned int toComp = 0; toComp < _disc.nComp; ++toComp)
					tripletList.push_back(T(idxr.offsetC() + comp, idxr.offsetC() + toComp, 0.0));
			}
		}
		else
		{
			for (unsigned int comp = 0; comp < _disc.nComp; ++comp)
				tripletList.push_back(T(idxr.offsetC() + comp, idxr.offsetC() + comp, 0.0));
		}

		// particle jacobian (including film diffusion, isotherm, reaction and time derivative)
		for (unsigned int type = 0; type < _disc.nParType; ++type)
			_particles[type]->setParJacPattern(tripletList, idxr.offsetCp(ParticleTypeIndex{ type }), idxr.offsetC(), 0u, secIdx);

		const auto& cm = _reaction.conservedMoieties("liquid");
		if (cm.isEnabled() && cm.numEquilibriumReactions() > 0)
			cm.addPatternToBlocks(tripletList, _disc.nComp, idxr.offsetC(), 1u, _disc.nComp, 0, globalJ.cols());

		globalJ.setFromTriplets(tripletList.begin(), tripletList.end());
	}

	/**
	 * @brief Computes the Jacobian via finite differences (testing purpose)
	 */
	Eigen::MatrixXd calcFDJacobian(const double* y_, const double* yDot_, const SimulationTime simTime, util::ThreadLocalStorage& threadLocalMem, double alpha)
	{
		Eigen::Map<const Eigen::VectorXd> hmpf(y_, numDofs());
		Eigen::VectorXd y = hmpf;
		Eigen::VectorXd yDot;
		if (yDot_)
		{
			Eigen::Map<const Eigen::VectorXd> hmpf2(yDot_, numDofs());
			yDot = hmpf2;
		}
		else
			return Eigen::MatrixXd::Zero(numDofs(), numDofs());

		Eigen::VectorXd res = Eigen::VectorXd::Zero(numDofs());
		const double* yPtr = &y[0];
		const double* yDotPtr = &yDot[0];
		double* resPtr = &res[0];

		Eigen::MatrixXd jacobian = Eigen::MatrixXd::Zero(numDofs(), numDofs());
		const double epsilon = 0.01;

		residualImpl<double, double, double, false>(simTime.t, simTime.secIdx, yPtr, yDotPtr, resPtr, threadLocalMem);

		for (int col = 0; col < jacobian.cols(); ++col)
			jacobian.col(col) = -(1.0 + alpha) * res;

		for (int dof = 0; dof < jacobian.cols(); ++dof)
		{
			y[dof] += epsilon;
			residualImpl<double, double, double, false>(simTime.t, simTime.secIdx, yPtr, yDotPtr, resPtr, threadLocalMem);
			y[dof] -= epsilon;
			jacobian.col(dof) += res;
		}

		for (int dof = 0; dof < jacobian.cols(); ++dof)
		{
			yDot[dof] += epsilon;
			residualImpl<double, double, double, false>(simTime.t, simTime.secIdx, yPtr, yDotPtr, resPtr, threadLocalMem);
			yDot[dof] -= epsilon;
			jacobian.col(dof) += alpha * res;
		}

		jacobian /= epsilon;

		return jacobian;
	}
};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_FINITEBATH_HPP_
