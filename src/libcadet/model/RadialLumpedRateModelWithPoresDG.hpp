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
 * Defines the radial lumped rate model with pores (rLRMP) using a Discontinous Galerkin (DG) scheme.
 */

#ifndef LIBCADET_RADIALLUMPEDRATEMODELWITHPORESDG_HPP_
#define LIBCADET_RADIALLUMPEDRATEMODELWITHPORESDG_HPP_

#include "BindingModel.hpp"
#include "ParallelSupport.hpp"

#include "UnitOperationBase.hpp"
#include "cadet/SolutionExporter.hpp"
#include "model/parts/ConvectionDispersionOperatorDG.hpp"
#include "model/particle/ParticleModel.hpp"
#include "AutoDiff.hpp"
#include "linalg/SparseMatrix.hpp"
#include "linalg/BandMatrix.hpp"
#include "linalg/BandedEigenSparseRowIterator.hpp"
#include "linalg/EigenSolverWrapper.hpp"
#include "linalg/Gmres.hpp"
#include "Memory.hpp"
#include "model/ModelUtils.hpp"
#include "ParameterMultiplexing.hpp"

#include <array>
#include <vector>

#include "Benchmark.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;

namespace cadet
{

	namespace model
	{

		/**
		* @brief Radial lumped rate model of liquid column chromatography with pores
		* @details Implements the radial transport equations with DG discretization:
		*
		* @f[\begin{align}
				\frac{\partial c_i}{\partial t} &= - \frac{u}{\rho} \frac{\partial c_i}{\partial \rho} + \frac{1}{\rho} \frac{\partial}{\partial \rho}\left( D_{\rho,i} \rho \frac{\partial c_i}{\partial \rho} \right) - \frac{1 - \varepsilon_c}{\varepsilon_c} \frac{3 k_{f,i}}{r_p} j_{f,i} \\
				\frac{\partial c_{p,i}}{\partial t} + \frac{1 - \varepsilon_p}{\varepsilon_p} \frac{\partial q_i}{\partial t} &= \frac{3 k_{f,i}}{\varepsilon_p r_p} j_{f,i} \\
				a \frac{\partial q_i}{\partial t} &= f_{\text{iso}}(c_p, q)
			\end{align} @f]
		* with film diffusion flux @f$ j_{f,i} = c_i - c_{p,i} @f$ and Danckwerts boundary conditions
		* @f[ \begin{align}
			u c_{\text{in},i}(t) &= u c_i(t,\rho_{in}) - D_{\rho,i} \frac{\partial c_i}{\partial \rho}(t,\rho_{in}) \\
			\frac{\partial c_i}{\partial \rho}(t,\rho_{out}) &= 0
			\end{align} @f]
		*/
		class RadialLumpedRateModelWithPoresDG : public UnitOperationBase
		{
		public:

			RadialLumpedRateModelWithPoresDG(UnitOpIdx unitOpIdx);
			virtual ~RadialLumpedRateModelWithPoresDG() CADET_NOEXCEPT;

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

			static const char* identifier() { return "RADIAL_LUMPED_RATE_MODEL_WITH_PORES_DG"; }
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

			template <typename StateType, typename ResidualType, typename ParamType, bool wantJac, bool wantRes = true>
			int residualBulk(double t, unsigned int secIdx, StateType const* yBase, double const* yDotBase, ResidualType* resBase, util::ThreadLocalStorage& threadLocalMem);

			template <typename StateType, typename ResidualType, typename ParamType>
			int residualFlux(double t, unsigned int secIdx, StateType const* y, double const* yDot, ResidualType* res);

			void extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset);

			void assembleDiscretizedJacobian(double alpha, const Indexer& idxr);
			void addTimeDerivativeToJacobianNode(linalg::BandedEigenSparseRowIterator& jac, const Indexer& idxr, double alpha, double invBetaP, unsigned int node) const;

			/**
			 * @brief Computes film diffusion matrices M_K for each cell
			 * @detail M_K[cell][i,j] = ∫ L_i * L_j * ρ * k_f(ρ) dξ
			 *         For constant k_f: M_K = k_f * M_ρ
			 *         For variable k_f: Uses quadrature integration
			 */
			void computeFilmDiffusionMatrices();

			/**
			 * @brief Updates film diffusion coefficients at interfaces
			 * @detail Called when section-dependent or spatially-dependent k_f is used
			 */
			void updateFilmDiffusionValues(unsigned int secIdx, unsigned int comp);

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
			void checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const;
#endif

			class Discretization
			{
			public:
				unsigned int nComp; //!< Number of components
				unsigned int nElem; //!< Number of radial elements
				unsigned int polyDeg; //!< polynomial degree
				unsigned int nNodes; //!< Number of nodes per cell
				unsigned int nPoints; //!< Number of discrete Points

				unsigned int* nBound; //!< Array with number of bound states for each component
				unsigned int* boundOffset; //!< Array with offset to the first bound state of each component in the solid phase
				unsigned int strideBound; //!< Total number of bound states

				int curSection; //!< current section index
				bool newStaticJac; //!< determines wether static analytical jacobian needs to be computed (every section)

				~Discretization() // make sure this memory is freed correctly
				{
					delete[] nBound;
					delete[] boundOffset;
				}
			};

			Discretization _disc; //!< Discretization info

			// Radial convection dispersion operator for interstitial volume transport
			parts::RadialConvectionDispersionOperatorBaseDG _convDispOp;

			cadet::linalg::EigenSolverBase* _linearSolver; //!< Linear solver

			Eigen::SparseMatrix<double, RowMajor> _jac; //!< Jacobian
			Eigen::SparseMatrix<double, RowMajor> _jacDisc; //!< Jacobian with time derivatives from BDF method
			Eigen::MatrixXd _jacInlet; //!< Jacobian inlet DOF block matrix connects inlet DOFs to first bulk cells

			// Porosity parameters
			active _colPorosity; //!< Column porosity (external porosity) \f$ \varepsilon_c \f$
			active _parPorosity; //!< Particle porosity (internal porosity) \f$ \varepsilon_p \f$

			// Particle parameters
			active _parRadius; //!< Particle radius \f$ r_p \f$
			double _parGeomSurfToVol; //!< Particle surface to volume ratio (3.0 for spheres)

			// Film diffusion parameters
			std::vector<active> _filmDiffusion; //!< Film diffusion coefficient \f$ k_f \f$
			MultiplexMode _filmDiffusionMode; //!< Multiplexing mode for film diffusion
			IParameterParameterDependence* _filmDiffDep; //!< Film diffusion dependency on position

			// Film diffusion DG matrices (per cell)
			std::vector<Eigen::MatrixXd> _M_K; //!< Film diffusion mass matrices M_K for each cell
			std::vector<Eigen::MatrixXd> _filmDiffCoupling; //!< Precomputed M_ρ^{-1} * M_K for each cell
			bool _variableFilmDiff; //!< Flag for spatially-varying film diffusion
			std::vector<Eigen::VectorXd> _filmDiffAtNodes; //!< Per-cell film diffusion values at DG nodes

			// Pore accessibility
			std::vector<active> _poreAccessFactor; //!< Pore accessibility factor \f$ F_{\text{acc}} \f$
			MultiplexMode _poreAccessFactorMode; //!< Multiplexing mode for pore access factor

			std::vector<IDynamicReactionModel*> _dynReaction; //!< Dynamic reaction models (owned)
			ReactionSystem _reaction; //!< Reaction system wrapper for CellParameters (non-owning refs)

			bool _analyticJac; //!< Determines whether AD or analytic Jacobians are used
			unsigned int _jacobianAdDirs; //!< Number of AD seed vectors required for Jacobian computation
			bool _factorizeJacobian; //!< Determines whether the Jacobian needs to be factorized
			double* _tempState; //!< Temporary storage with the size of the state vector or larger if binding models require it

			std::vector<active> _initC; //!< Liquid bulk phase initial conditions
			std::vector<active> _initCp; //!< Liquid pore phase initial conditions
			std::vector<active> _initCs; //!< Solid phase initial conditions
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

				inline int strideParNode() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp + _disc.strideBound); }
				inline int strideParCell() const CADET_NOEXCEPT { return static_cast<int>(_disc.nNodes * strideParNode()); }
				inline int strideParLiquid() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp); }
				inline int strideParBound() const CADET_NOEXCEPT { return static_cast<int>(_disc.strideBound); }

				// Offsets
				inline int offsetC() const CADET_NOEXCEPT { return _disc.nComp; }
				inline int offsetCp() const CADET_NOEXCEPT { return _disc.nComp + _disc.nPoints * _disc.nComp; }
				inline int offsetBoundComp(unsigned int comp) const CADET_NOEXCEPT { return _disc.boundOffset[comp]; }

				// Return pointer to first element of state variable in state vector
				template <typename real_t> inline real_t* c(real_t* const data) const { return data + offsetC(); }
				template <typename real_t> inline real_t const* c(real_t const* const data) const { return data + offsetC(); }

				template <typename real_t> inline real_t* cp(real_t* const data) const { return data + offsetCp(); }
				template <typename real_t> inline real_t const* cp(real_t const* const data) const { return data + offsetCp(); }

				template <typename real_t> inline real_t* q(real_t* const data) const { return data + offsetCp() + strideParLiquid(); }
				template <typename real_t> inline real_t const* q(real_t const* const data) const { return data + offsetCp() + strideParLiquid(); }

				// Return specific variable in state vector
				template <typename real_t> inline real_t& c(real_t* const data, unsigned int node, unsigned int comp) const { return data[offsetC() + comp * strideColComp() + node * strideColNode()]; }
				template <typename real_t> inline const real_t& c(real_t const* const data, unsigned int node, unsigned int comp) const { return data[offsetC() + comp * strideColComp() + node * strideColNode()]; }

				template <typename real_t> inline real_t& cp(real_t* const data, unsigned int node, unsigned int comp) const { return data[offsetCp() + comp + node * strideParNode()]; }
				template <typename real_t> inline const real_t& cp(real_t const* const data, unsigned int node, unsigned int comp) const { return data[offsetCp() + comp + node * strideParNode()]; }

			protected:
				const Discretization& _disc;
			};

			class Exporter : public ISolutionExporter
			{
			public:

				Exporter(const Discretization& disc, const RadialLumpedRateModelWithPoresDG& model, double const* data) : _disc(disc), _idx(disc), _model(model), _data(data) { }
				Exporter(const Discretization&& disc, const RadialLumpedRateModelWithPoresDG& model, double const* data) = delete;

				virtual bool hasParticleFlux() const CADET_NOEXCEPT { return false; }
				virtual bool hasParticleMobilePhase() const CADET_NOEXCEPT { return true; }
				virtual bool hasSolidPhase() const CADET_NOEXCEPT { return _disc.strideBound > 0; }
				virtual bool hasVolume() const CADET_NOEXCEPT { return false; }
				virtual bool isParticleLumped(unsigned int parType) const CADET_NOEXCEPT { return true; }
				virtual bool hasPrimaryExtent() const CADET_NOEXCEPT { return true; }

				virtual unsigned int numComponents() const CADET_NOEXCEPT { return _disc.nComp; }
				virtual unsigned int numPrimaryCoordinates() const CADET_NOEXCEPT { return _disc.nPoints; }
				virtual unsigned int numSecondaryCoordinates() const CADET_NOEXCEPT { return 0; }
				virtual unsigned int numInletPorts() const CADET_NOEXCEPT { return 1; }
				virtual unsigned int numOutletPorts() const CADET_NOEXCEPT { return 1; }
				virtual unsigned int numParticleTypes() const CADET_NOEXCEPT { return 1; }
				virtual unsigned int numParticleShells(unsigned int parType) const CADET_NOEXCEPT { return 1; }
				virtual unsigned int numBoundStates(unsigned int parType) const CADET_NOEXCEPT { return _disc.strideBound; }
				virtual unsigned int numMobilePhaseDofs() const CADET_NOEXCEPT { return _disc.nComp * _disc.nPoints; }
				virtual unsigned int numParticleMobilePhaseDofs(unsigned int parType) const CADET_NOEXCEPT { return _disc.nComp * _disc.nPoints; }
				virtual unsigned int numParticleMobilePhaseDofs() const CADET_NOEXCEPT { return _disc.nComp * _disc.nPoints; }
				virtual unsigned int numSolidPhaseDofs(unsigned int parType) const CADET_NOEXCEPT { return _disc.strideBound * _disc.nPoints; }
				virtual unsigned int numSolidPhaseDofs() const CADET_NOEXCEPT { return _disc.strideBound * _disc.nPoints; }
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
				/**
				 * @brief calculates the physical node coordinates of the DG discretization
				 */
				virtual int writePrimaryCoordinates(double* coords) const
				{
					for (unsigned int i = 0; i < _disc.nElem; i++) {
						for (unsigned int j = 0; j < _disc.nNodes; j++) {
							// mapping for radial coordinates
							coords[i * _disc.nNodes + j] = _model._convDispOp.elemLeftBound(i) +
								0.5 * (static_cast<double>(_model._convDispOp.columnLength()) / static_cast<double>(_disc.nElem)) *
								(1.0 + _model._convDispOp.LGLnodes()[j]);
						}
					}
					return _disc.nPoints;
				}

				virtual int writeSecondaryCoordinates(double* coords) const { return 0; }
				virtual int writeParticleCoordinates(unsigned int parType, double* coords) const { return 0; }

			protected:
				const Discretization& _disc;
				const Indexer _idx;
				const RadialLumpedRateModelWithPoresDG& _model;
				double const* const _data;
			};

			/**
			 * @brief sets the current section index and section dependend velocity, dispersion
			 */
			void updateSection(int secIdx) {

				if (_disc.curSection != secIdx) {
					_disc.curSection = secIdx;
					_disc.newStaticJac = true;
				}
			}

			// ==========================================================================================================================================================  //
			// ========================================						DG Jacobian							=========================================================  //
			// ==========================================================================================================================================================  //

			typedef Eigen::Triplet<double> T;

			/**
			 * @brief sets the sparsity pattern of the static Jacobian
			 * @detail DG ConvDisp pattern, film diffusion coupling, and isotherm pattern.
			 * @param [in] stateDer bool if state derivative pattern should be added (for _jacDisc)
			 */
			void setPattern(Eigen::SparseMatrix<double, RowMajor>& mat, bool stateDer, bool has_reaction);

			/**
			 * @brief computes the convection dispersion part of the state derivative pattern of the jacobian.
			 */
			void stateDerPattern(std::vector<T>& tripletList);

			/**
			* @brief sets the sparsity pattern of the isotherm and reaction Jacobian
			*/
			void bindingAndReactionPattern(std::vector<T>& tripletList, bool has_reaction);

			/**
			 * @brief sets the sparsity pattern of the film diffusion coupling Jacobian
			 */
			void filmDiffusionPattern(std::vector<T>& tripletList);

		};

	} // namespace model
} // namespace cadet

#endif  // LIBCADET_RADIALLUMPEDRATEMODELWITHPORESDG_HPP_
