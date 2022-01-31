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
 * Defines the lumped rate model without pores (LRM) for DG approach.
 */

#ifndef LIBCADET_LUMPEDRATEMODELWITHOUTPORESDG_HPP_
#define LIBCADET_LUMPEDRATEMODELWITHOUTPORESDG_HPP_

#include "UnitOperationBase.hpp"
#include "cadet/SolutionExporter.hpp"
#include "model/parts/ConvectionDispersionOperator.hpp"
#include "AutoDiff.hpp"
#include "linalg/SparseMatrix.hpp"
#include "linalg/BandMatrix.hpp"
#include "linalg/Gmres.hpp"
#include "Memory.hpp"
#include "model/ModelUtils.hpp"

#include <array>
#include <vector>

#include "Benchmark.hpp"
#include <numbers>
//#include "DGspecific.hpp" TODO: source out DG specifics
#include "C:\Users\jmbr\Cadet\libs\eigen-3.4.0\Eigen\Dense.hpp"	// use LA lib Eigen for Matrix operations
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#include <math.h>
#endif
using namespace Eigen;

namespace cadet
{

namespace model
{

/**
 * @brief Lumped rate model of liquid column chromatography without pores
 * @details See @cite Guiochon2006, @cite Gu1995, @cite Felinger2004
 * 
 * @f[\begin{align}
	\frac{\partial c_i}{\partial t} + \frac{1 - \varepsilon_t}{\varepsilon_t} \frac{\partial q_{i}}{\partial t} &= - u \frac{\partial c_i}{\partial z} + D_{\text{ax},i} \frac{\partial^2 c_i}{\partial z^2} \\
	a \frac{\partial q_i}{\partial t} &= f_{\text{iso}}(c, q)
\end{align} @f]
 * Danckwerts boundary conditions (see @cite Danckwerts1953)
@f[ \begin{align}
u c_{\text{in},i}(t) &= u c_i(t,0) - D_{\text{ax},i} \frac{\partial c_i}{\partial z}(t,0) \\
\frac{\partial c_i}{\partial z}(t,L) &= 0
\end{align} @f]
 * Methods are described in @cite VonLieres2010a (WENO, linear solver), @cite Puttmann2013 @cite Puttmann2016 (forward sensitivities, AD, band compression)
 */
class LumpedRateModelWithoutPoresDG : public UnitOperationBase
{
public:

	LumpedRateModelWithoutPoresDG(UnitOpIdx unitOpIdx);
	virtual ~LumpedRateModelWithoutPoresDG() CADET_NOEXCEPT;

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

	static const char* identifier() { return "LUMPED_RATE_MODEL_WITHOUT_PORES_DG"; }
	virtual const char* unitOperationName() const CADET_NOEXCEPT { return identifier(); }

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

	template <typename StateType, typename ResidualType, typename ParamType, bool wantJac>
	int residualImpl(double t, unsigned int secIdx, StateType const* const y, double const* const yDot, ResidualType* const res, util::ThreadLocalStorage& threadLocalMem);

	void extractJacobianFromAD(active const* const adRes, unsigned int adDirOffset);

	void assembleDiscretizedJacobian(double alpha, const Indexer& idxr);
	void addTimeDerivativeToJacobianCell(linalg::FactorizableBandMatrix::RowIterator& jac, const Indexer& idxr, double alpha, double invBetaP) const;

#ifdef CADET_CHECK_ANALYTIC_JACOBIAN
	void checkAnalyticJacobianAgainstAd(active const* const adRes, unsigned int adDirOffset) const;
#endif

	struct Discretization
	{
		unsigned int nComp; //!< Number of components
		unsigned int nCol; //!< Number of column cells
		unsigned int polyDeg; //!< polynomial degree
		unsigned int nNodes; //!< Number of nodes in each cell = polynomial degree - 1
		Eigen::VectorXd nodes;		//!< Array with positions of nodes in reference element
		Eigen::VectorXd weights; //!< Array with weights for numerical quadrature of size nNodes
		Eigen::MatrixXd polyDerM; //!< Array with polynomial derivative Matrix of size nNodes^2
		Eigen::MatrixXd polyDerMtranspose; //!< Array with D^T
		unsigned int* nBound; //!< Array with number of bound states for each component
		unsigned int* boundOffset; //!< Array with offset to the first bound state of each component in the solid phase
		unsigned int strideBound; //!< Total number of bound states

		/**
		 * @brief sets DG specific coefficients
		 */
		void initializeDG(cadet::model::LumpedRateModelWithoutPoresDG::Discretization _disc) {
			 lglNodesWeights(_disc.polyDeg, _disc.nodes, _disc.weights);
			 polynomialDerivativeMatrix(_disc.polyDeg, _disc.nodes, _disc.polyDerM, _disc.polyDerMTrans);
		 }

		/**
		 * @brief computes the Legendre polynomial L_N and q = L_N+1 - L_N-2 and q' at point x
		 * @param [in] polyDeg polynomial degree of spatial Discretization
		 * @param [in] x evaluation point
		 * @param [in] L pointer for allocated memory of L(x)
		 * @param [in] q pointer for allocated memory of q(x)
		 * @param [in] qder pointer for allocated memory of q'(x)
		 */
		void qAndL(int polyDeg, double x, double* L, double* q, double* qder) {
			// auxiliary variables (Legendre polynomials)
			double L_2 = 1.0;
			double L_1 = x;
			double Lder_2 = 0.0;
			double Lder_1 = 1.0;
			double Lder = 0.0;
			for (double k = 2; k <= polyDeg; k++) {
				*L = ((2 * k - 1) * x * L_1 - (k - 1) * L_2) / k;
				Lder = Lder_2 + (2 * k - 1) * L_1;
				L_2 = L_1;
				L_1 = *L;
				Lder_2 = Lder_1;
				Lder_1 = Lder;
			}
			*q = ((2.0 * polyDeg + 1) * x * *L - polyDeg * L_2) / (polyDeg + 1.0) - L_2;
			*qder = Lder_1 + (2.0 * polyDeg + 1) * L_1 - Lder_2;
		}

		/**
		 * @brief computes and assigns the Legendre-Gauss nodes and weights
		 * @param [in] polyDeg polynomial degree of spatial Discretization
		 * @param [in] pre-initialized nodes array
		 * @param [in] pre-initialized weights array
		 */
		void lglNodesWeights(int polyDeg, VectorXd& nodes, VectorXd& weights) {
			// tolerance and max #iterations for Newton iteration
			int nIterations = 10;
			double tolerance = 1e-15;
			// Legendre polynomial and derivative
			double L = 0;
			double q = 0;
			double qder = 0;
			switch (polyDeg) {
			case 0:
				throw std::invalid_argument("Polynomial degree must be at least 1 !");
				break;
			case 1:
				nodes[0] = -1;
				weights[0] = 1;
				nodes[1] = 1;
				weights[1] = 1;
				break;
			default:
				nodes[0] = -1;
				nodes[polyDeg] = 1;
				weights[0] = 2.0 / (polyDeg * (polyDeg + 1.0));
				weights[polyDeg] = weights[0];
				// use symmetrie, only compute half of points and weights
				for (int j = 1; j <= floor((polyDeg + 1) / 2) - 1; j++) {
					//  first guess for Newton iteration
					nodes[j] = -cos(M_PI * (j + 0.25) / polyDeg - 3 / (8.0 * polyDeg * M_PI * (j + 0.25)));
					// Newton iteration to find zero points of Legendre Polynomial
					for (int k = 0; k <= nIterations; k++) {
						qAndL(polyDeg, nodes[j], &L, &q, &qder);
						nodes[j] = nodes[j] - q / qder;
						if (abs(q / qder) <= tolerance * abs(nodes[j])) {
							break;
						}
					}
					// calculate weights
					qAndL(polyDeg, nodes[j], &L, &q, &qder);
					weights[j] = 2 / (polyDeg * (polyDeg + 1.0) * pow(L, 2));
					nodes[polyDeg - j] = -nodes[j]; // copy to second half of points and weights
					weights[polyDeg - j] = weights[j];
				}
			}
			if (polyDeg % 2 == 0) { // for even polyDeg we have an odd number of points which include 0.0
				qAndL(polyDeg, 0.0, &L, &q, &qder);
				nodes[polyDeg / 2] = 0;
				weights[polyDeg / 2] = 2 / (polyDeg * (polyDeg + 1.0) * pow(L, 2));
			}
		}

		/**
		 * @brief computation of barycentric weights for fast polynomial evaluation
		 * @param [in] polyDeg polynomial degree of spatial Discretization
		 * @param [in] nodes node vector
		 * @param [in] weights pre-allocated vector for barycentric weights. Must be set to ones!
		 */
		void barycentricWeights(int polyDeg, VectorXd& nodes, VectorXd& baryWeights) {
			for (int j = 1; j <= polyDeg; j++) {
				for (int k = 0; k <= j - 1; k++) {
					baryWeights[k] = baryWeights[k] * (nodes[k] - nodes[j]) * 1.0;
					baryWeights[j] = baryWeights[j] * (nodes[j] - nodes[k]) * 1.0;
				}
			}
			for (int j = 0; j <= polyDeg; j++) {
				baryWeights[j] = 1 / baryWeights[j];
			}
		}

		/**
		 * @brief initialization of polynomial derivative matrix D and D^T
		 * @param [in] polyDeg polynomial degree of spatial Discretization
		 * @param [in] nodes LGL node vector
		 * @param [in] D pre-allocated derivative matrix
		 * @param [in] DT pre-allocated transposed derivative matrix
		 */
		void polynomialDerivativeMatrix(int polyDeg, VectorXd& nodes, MatrixXd& D, MatrixXd& DT) {
			VectorXd baryWeights = VectorXd::Ones(polyDeg + 1);
			barycentricWeights(polyDeg, nodes, baryWeights);
			for (int i = 0; i <= polyDeg; i++) {
				for (int j = 0; j <= polyDeg; j++) {
					if (i != j) {
						D(i, j) = baryWeights[j] / (baryWeights[i] * (nodes[i] - nodes[j]));
						D(i, i) += -D(i, j);
					}
				}
			}
			DT = D.transpose();
		}

		/** TODO add polynomial evaluation function with baryw?
		* @brief fast polynomial evaluation
		*/
	};

	Discretization _disc; //!< Discretization info
//	IExternalFunction* _extFun; //!< External function (owned by library user)

	parts::ConvectionDispersionOperatorBase _convDispOp; //!< Convection dispersion operator for interstitial volume transport

	linalg::BandMatrix _jac; //!< Jacobian
	linalg::FactorizableBandMatrix _jacDisc; //!< Jacobian with time derivatives from BDF method

	linalg::DoubleSparseMatrix _jacInlet; //!< Jacobian inlet DOF block matrix connects inlet DOFs to first bulk cells

	active _totalPorosity; //!< Total porosity \f$ \varepsilon_t \f$

	bool _analyticJac; //!< Determines whether AD or analytic Jacobians are used
	unsigned int _jacobianAdDirs; //!< Number of AD seed vectors required for Jacobian computation

	bool _factorizeJacobian; //!< Determines whether the Jacobian needs to be factorized
	double* _tempState; //!< Temporary storage with the size of the state vector or larger if binding models require it

	std::vector<active> _initC; //!< Liquid phase initial conditions
	std::vector<active> _initQ; //!< Solid phase initial conditions
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
		inline int strideColCell() const CADET_NOEXCEPT { return static_cast<int>((_disc.nComp + _disc.strideBound) * _disc.nNodes); }
		inline int strideColNode() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp + _disc.strideBound); }
		inline int strideColComp() const CADET_NOEXCEPT { return 1; }

		inline int strideColLiquid() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp); }
		inline int strideColBound() const CADET_NOEXCEPT { return static_cast<int>(_disc.strideBound); }

		// Offsets
		inline int offsetC() const CADET_NOEXCEPT { return _disc.nComp; }
		inline int offsetBoundComp(unsigned int comp) const CADET_NOEXCEPT { return _disc.boundOffset[comp]; }

		// Return pointer to first element of state variable in state vector
		template <typename real_t> inline real_t* c(real_t* const data) const { return data + offsetC(); }
		template <typename real_t> inline real_t const* c(real_t const* const data) const { return data + offsetC(); }

		template <typename real_t> inline real_t* q(real_t* const data) const { return data + offsetC() + strideColLiquid(); }
		template <typename real_t> inline real_t const* q(real_t const* const data) const { return data + offsetC() + strideColLiquid(); }

		// Return specific variable in state vector
		template <typename real_t> inline real_t& c(real_t* const data, unsigned int col, unsigned int comp) const { return data[offsetC() + comp + col * strideColCell()]; }
		template <typename real_t> inline const real_t& c(real_t const* const data, unsigned int col, unsigned int comp) const { return data[offsetC() + comp + col * strideColCell()]; }

	protected:
		const Discretization& _disc;
	};

	class Exporter : public ISolutionExporter
	{
	public:

		Exporter(const Discretization& disc, const LumpedRateModelWithoutPoresDG& model, double const* data) : _disc(disc), _idx(disc), _model(model), _data(data) { }
		Exporter(const Discretization&& disc, const LumpedRateModelWithoutPoresDG& model, double const* data) = delete;

		virtual bool hasParticleFlux() const CADET_NOEXCEPT { return false; }
		virtual bool hasParticleMobilePhase() const CADET_NOEXCEPT { return false; }
		virtual bool hasSolidPhase() const CADET_NOEXCEPT { return _disc.strideBound > 0; }
		virtual bool hasVolume() const CADET_NOEXCEPT { return false; }

		virtual unsigned int numComponents() const CADET_NOEXCEPT { return _disc.nComp; }
		virtual unsigned int numAxialCells() const CADET_NOEXCEPT { return _disc.nCol; }
		virtual unsigned int numRadialCells() const CADET_NOEXCEPT { return 0; }
		virtual unsigned int numInletPorts() const CADET_NOEXCEPT { return 1; }
		virtual unsigned int numOutletPorts() const CADET_NOEXCEPT { return 1; }
		virtual unsigned int numParticleTypes() const CADET_NOEXCEPT { return 1u; }
		virtual unsigned int numParticleShells(unsigned int parType) const CADET_NOEXCEPT { return 0u; }
		virtual unsigned int numBoundStates(unsigned int parType) const CADET_NOEXCEPT { return _disc.strideBound; }
		virtual unsigned int numBulkDofs() const CADET_NOEXCEPT { return _disc.nComp * _disc.nCol; }
		virtual unsigned int numParticleMobilePhaseDofs(unsigned int parType) const CADET_NOEXCEPT { return 0u; }
		virtual unsigned int numSolidPhaseDofs(unsigned int parType) const CADET_NOEXCEPT { return _disc.strideBound * _disc.nCol; }
		virtual unsigned int numFluxDofs() const CADET_NOEXCEPT { return 0u; }
		virtual unsigned int numVolumeDofs() const CADET_NOEXCEPT { return 0u; }

		virtual double const* concentration() const { return _idx.c(_data); }
		virtual double const* flux() const { return nullptr; }
		virtual double const* particleMobilePhase(unsigned int parType) const { return nullptr; }
		virtual double const* solidPhase(unsigned int parType) const { return _idx.q(_data); }
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

		virtual unsigned int bulkMobilePhaseStride() const { return _idx.strideColCell(); }
		virtual unsigned int particleMobilePhaseStride(unsigned int parType) const { return 0; }
		virtual unsigned int solidPhaseStride(unsigned int parType) const { return _idx.strideColCell(); }

		virtual void axialCoordinates(double* coords) const
		{
			const double h = static_cast<double>(_model._convDispOp.columnLength()) / static_cast<double>(_disc.nCol);
			for (unsigned int i = 0; i < _disc.nCol; ++i)
				coords[i] = (i + 0.5) * h;
		}
		virtual void radialCoordinates(double* coords) const { }
		virtual void particleCoordinates(unsigned int parType, double* coords) const { }

	protected:
		const Discretization& _disc;
		const Indexer _idx;
		const LumpedRateModelWithoutPoresDG& _model;
		double const* const _data;

		const std::array<StateOrdering, 2> _concentrationOrdering = { { StateOrdering::AxialCell, StateOrdering::Component } };
		const std::array<StateOrdering, 3> _solidOrdering = { { StateOrdering::AxialCell, StateOrdering::Component, StateOrdering::BoundState } };
	};
};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_LUMPEDRATEMODELWITHOUTPORESDG_HPP_
