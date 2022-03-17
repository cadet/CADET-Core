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
 * Defines the lumped rate model without pores (LRM) using a Discontinous Galerkin (DG) scheme.
 */

#ifndef LIBCADET_LUMPEDRATEMODELWITHOUTPORESDG_HPP_
#define LIBCADET_LUMPEDRATEMODELWITHOUTPORESDG_HPP_

#include "BindingModel.hpp"
#include "ParallelSupport.hpp"

#include "UnitOperationBase.hpp"
#include "cadet/SolutionExporter.hpp"
#include "model/parts/ConvectionDispersionOperator.hpp"
#include "AutoDiff.hpp"
#include "linalg/SparseMatrix.hpp"
#include "linalg/BandMatrix.hpp"
#include "linalg/BandedEigenSparseRowIterator.hpp"
#include "linalg/Gmres.hpp"
#include "Memory.hpp"
#include "model/ModelUtils.hpp"

#include <array>
#include <vector>

#include "Benchmark.hpp"
#include <numbers>

#include "C:\Users\jmbr\Cadet\libs\eigen-3.4.0\Eigen\Dense.hpp"	// use LA lib Eigen for Matrix operations
#include "C:\Users\jmbr\Cadet\libs\eigen-3.4.0\Eigen\Sparse.hpp"

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

	class Discretization // @TODO: separate convDisp DGoperator class
	{ 	// NOTE: no different Riemann solvers or boundary conditions
	public:
		unsigned int nComp; //!< Number of components
		unsigned int nCol; //!< Number of column cells
		unsigned int polyDeg; //!< polynomial degree
		unsigned int nNodes; //!< Number of nodes per cell
		unsigned int nPoints; //!< Number of discrete Points
		double deltaZ;
		bool modal;	//!< bool switch: 1 for modal basis, 0 for nodal basis
		Eigen::VectorXd nodes; //!< Array with positions of nodes in reference element
		Eigen::MatrixXd polyDerM; //!< Array with polynomial derivative Matrix
		Eigen::VectorXd invWeights; //!< Array with weights for numerical quadrature of size nNodes
		Eigen::MatrixXd invMM; //!< dense !INVERSE! mass matrix for modal (exact) integration

		unsigned int* nBound; //!< Array with number of bound states for each component
		unsigned int* boundOffset; //!< Array with offset to the first bound state of each component in the solid phase
		unsigned int strideBound; //!< Total number of bound states

		// vgl. convDispOperator
		Eigen::VectorXd dispersion; //!< Column dispersion (may be section and component dependent)
		bool _dispersionCompIndep; //!< Determines whether dispersion is component independent
		double velocity; //!< Interstitial velocity (may be section dependent) \f$ u \f$
		double curVelocity;
		int curSection; //!< current section index

		double length_;
		double porosity;
		// Section dependent parameters

		Eigen::VectorXd g;
		Eigen::VectorXd h;
		Eigen::VectorXd surfaceFlux; //!< stores the surface flux values
		Eigen::Vector4d boundary; //!< stores the boundary values from Danckwert boundary conditions

		std::vector<bool> isKinetic;

		bool newStaticJac; //!< determines wether static analytical jacobian needs to be computed (every section)

		/**
		* @brief computes LGL nodes, integration weights, polynomial derivative matrix
		*/
		void initializeDG() {

			nNodes = polyDeg + 1;
			nPoints = nNodes * nCol;
			// Allocate space for DG discretization
			nodes.resize(nNodes);
			nodes.setZero();
			invWeights.resize(nNodes);
			invWeights.setZero();
			polyDerM.resize(nNodes, nNodes);
			polyDerM.setZero();
			invMM.resize(nNodes, nNodes);
			invMM.setZero();

			g.resize(nPoints);
			g.setZero();
			h.resize(nPoints);
			h.setZero();
			boundary.setZero();
			surfaceFlux.resize(nCol + 1);
			surfaceFlux.setZero();

			newStaticJac = true;

			lglNodesWeights();
			// @TODO: make modal/nodal switch during calculation possible?
			if (modal) {
				derivativeMatrix_JACOBI();
				invMMatrix_JACOBI();
			}
			else {
				derivativeMatrix_LAGRANGE();
				invMMatrix_JACOBI(); // modal/nodal switch
			}
		}

	private:

		/* ===================================================================================
		*   Lagrange Basis operators and auxiliary functions
		* =================================================================================== */

		/**
		* @brief computes the Legendre polynomial L_N and q = L_N+1 - L_N-2 and q' at point x
		* @param [in] polyDeg polynomial degree of spatial Discretization
		* @param [in] x evaluation point
		* @param [in] L pre-allocated L(x)
		* @param [in] q pre-allocated q(x)
		* @param [in] qder pre-allocated q'(x)
		*/
		void qAndL(const double x, double& L, double& q, double& qder) {
			// auxiliary variables (Legendre polynomials)
			double L_2 = 1.0;
			double L_1 = x;
			double Lder_2 = 0.0;
			double Lder_1 = 1.0;
			double Lder = 0.0;
			for (double k = 2; k <= polyDeg; k++) {
				L = ((2 * k - 1) * x * L_1 - (k - 1) * L_2) / k;
				Lder = Lder_2 + (2 * k - 1) * L_1;
				L_2 = L_1;
				L_1 = L;
				Lder_2 = Lder_1;
				Lder_1 = Lder;
			}
			q = ((2.0 * polyDeg + 1) * x * L - polyDeg * L_2) / (polyDeg + 1.0) - L_2;
			qder = Lder_1 + (2.0 * polyDeg + 1) * L_1 - Lder_2;
		}

		/**
		 * @brief computes and assigns the Legendre-Gauss nodes and (inverse) weights
		 */
		void lglNodesWeights() {
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
				invWeights[0] = 1;
				nodes[1] = 1;
				invWeights[1] = 1;
				break;
			default:
				nodes[0] = -1;
				nodes[polyDeg] = 1;
				invWeights[0] = 2.0 / (polyDeg * (polyDeg + 1.0));
				invWeights[polyDeg] = invWeights[0];
				// use symmetrie, only compute half of points and weights
				for (unsigned int j = 1; j <= floor((polyDeg + 1) / 2) - 1; j++) {
					//  first guess for Newton iteration
					nodes[j] = -cos(M_PI * (j + 0.25) / polyDeg - 3 / (8.0 * polyDeg * M_PI * (j + 0.25)));
					// Newton iteration to find roots of Legendre Polynomial
					for (unsigned int k = 0; k <= nIterations; k++) {
						qAndL(nodes[j], L, q, qder);
						nodes[j] = nodes[j] - q / qder;
						if (abs(q / qder) <= tolerance * abs(nodes[j])) {
							break;
						}
					}
					// calculate weights
					qAndL(nodes[j], L, q, qder);
					invWeights[j] = 2.0 / (polyDeg * (polyDeg + 1.0) * pow(L, 2.0));
					nodes[polyDeg - j] = -nodes[j]; // copy to second half of points and weights
					invWeights[polyDeg - j] = invWeights[j];
				}
			}
			if (polyDeg % 2 == 0) { // for even polyDeg we have an odd number of points which include 0.0
				qAndL(0.0, L, q, qder);
				nodes[polyDeg / 2] = 0;
				invWeights[polyDeg / 2] = 2.0 / (polyDeg * (polyDeg + 1.0) * pow(L, 2.0));
			}
			// inverse the weights
			invWeights = invWeights.cwiseInverse();
		}

		/**
		 * @brief computation of barycentric weights for fast polynomial evaluation
		 * @param [in] weights pre-allocated vector for barycentric weights. Must be set to ones!
		 */
		void barycentricWeights(VectorXd& baryWeights) {
			for (unsigned int j = 1; j <= polyDeg; j++) {
				for (unsigned int k = 0; k <= j - 1; k++) {
					baryWeights[k] = baryWeights[k] * (nodes[k] - nodes[j]) * 1.0;
					baryWeights[j] = baryWeights[j] * (nodes[j] - nodes[k]) * 1.0;
				}
			}
			for (unsigned int j = 0; j <= polyDeg; j++) {
				baryWeights[j] = 1 / baryWeights[j];
			}
		}

		/**
		 * @brief computation of LAGRANGE polynomial derivative matrix
		 */
		void derivativeMatrix_LAGRANGE() {
			VectorXd baryWeights = VectorXd::Ones(polyDeg + 1u);
			barycentricWeights(baryWeights);
			for (unsigned int i = 0; i <= polyDeg; i++) {
				for (unsigned int j = 0; j <= polyDeg; j++) {
					if (i != j) {
						polyDerM(i, j) = baryWeights[j] / (baryWeights[i] * (nodes[i] - nodes[j]));
						polyDerM(i, i) += -polyDerM(i, j);
					}
				}
			}
		}

		/* ===================================================================================
		*   Jacobi Basis operators and auxiliary functions
		* =================================================================================== */

		/*
		* @brief computation of normalized Jacobi polynomial P of order N at nodes x
		*/
		void jacobiPolynomial(const double alpha, const double beta, const int N, MatrixXd& P, int index) {
			// factor needed to normalize the Jacobi polynomial using the gamma function
			double gamma0 = std::pow(2.0, alpha + beta + 1.0) / (alpha + beta + 1.0) * std::tgamma(alpha + 1.0) * std::tgamma(beta + 1) / std::tgamma(alpha + beta + 1);
			MatrixXd PL(N + 1, nodes.size());
			for (unsigned int i = 0; i < nodes.size(); ++i) {
				PL(0, i) = 1.0 / std::sqrt(gamma0);
			}
			if (N == 0) {
				for (unsigned int i = 0; i < nodes.size(); ++i) {
					P(i, index) = PL(0, i);
				}
				return;
			}
			double gamma1 = (alpha + 1) * (beta + 1) / (alpha + beta + 3) * gamma0;
			for (unsigned int i = 0; i < nodes.size(); ++i) {
				PL(1, i) = ((alpha + beta + 2) * nodes(i) / 2 + (alpha - beta) / 2) / std::sqrt(gamma1);
			}
			if (N == 1) {
				for (unsigned int i = 0; i < nodes.size(); ++i) {
					P(i, index) = PL(1, i);
				}
				return;
			}
			double a = 2.0 / (2.0 + alpha + beta) * std::sqrt((alpha + 1) * (beta + 1) / (alpha + beta + 3));
			for (unsigned int i = 0; i < N - 1; ++i) {
				double j = i + 1.0;
				double h1 = 2.0 * (i + 1.0) + alpha + beta;
				double a_aux = 2.0 / (h1 + 2.0) * std::sqrt((j + 1.0) * (j + 1.0 + alpha + beta) * (j + 1.0 + alpha) * (j + 1.0 + beta) / (h1 + 1) / (h1 + 3));
				double b = -(std::pow(alpha, 2) - std::pow(beta, 2)) / h1 / (h1 + 2);

				for (unsigned int k = 0; k < nodes.size(); ++k) {
					PL(i + 2, k) = 1 / a_aux * (-a * PL(i, k) + (nodes(k) - b) * PL(i + 1, k));
				}
				a = a_aux;
			}
			for (unsigned int i = 0; i < nodes.size(); ++i) {
				P(i, index) = PL(N, i);
			}
		}

		/**
		* @brief returns generalized Jacobi Vandermonde matrix
		* @detail normalized Legendre Vandermonde matrix
		*/
		MatrixXd getVandermonde_JACOBI() {

			MatrixXd V(nodes.size(), nodes.size());

			for (unsigned int j = 0; j < V.cols(); ++j) {
				// legendre polynomial: alpha = beta = 0.0
				jacobiPolynomial(0.0, 0.0, j, V, j);
			}
			return V;
		}

		/*
		* @brief computes the gradient vandermonde matrix of orthonormal Legendre polynomials
		* @V_x [in] pre-allocated gradient Vandermonde matrix
		*/
		void GradVandermonde(MatrixXd& V_x) {
			// Legendre polynomial
			double alpha = 0.0;
			double beta = 0.0;

			for (unsigned int order = 1; order < nodes.size(); order++) {
				jacobiPolynomial(alpha + 1, beta + 1, order - 1, V_x, order);
				V_x.block(0, order, V_x.rows(), 1) *= std::sqrt(order * (order + alpha + beta + 1));
			}
		}

		//@TODO?: Not needed? as the D matrix is (approx) the same as for the nodal approach ! (and computed without matrix inversion or linear solve)
		//        D_nod is approx D_mod to 1e-14 but not up to std::numeric_limits<double>::epsilon()... relevant ?
		/**
		 * @brief computes the Jacobi polynomial derivative matrix D
		 * @detail computes the normalized Legendre polynomial derivative matrix D
		 */
		void derivativeMatrix_JACOBI() {

			// Compute the gradient Vandermonde matrix and transpose
			const int Np = nodes.size();
			GradVandermonde(polyDerM);
			MatrixXd V_xT = polyDerM.transpose();

			// Instead of using matrix inversion, solve the linear systems to obtain D using SVD decomposition (slow but accurate)
			for (unsigned int i = 0; i < Np; ++i) {
				polyDerM.block(i, 0, 1, Np) = (getVandermonde_JACOBI().transpose().jacobiSvd(Eigen::ComputeFullU |
					Eigen::ComputeFullV).solve(V_xT.block(0, i, Np, 1))).transpose();
			}
		}

		/**
		* @brief returns Jacobi polynomial induced mass matrix
		* @detail returns normalized Legendre polynomial induced mass matrix
		* @param [in] nodes (LGL)
		*/
		void invMMatrix_JACOBI() {
			invMM = (getVandermonde_JACOBI() * (getVandermonde_JACOBI().transpose()));
		}

	};

	Discretization _disc; //!< Discretization info

//	IExternalFunction* _extFun; //!< External function (owned by library user)

	// used as auxiliary supplier 
	parts::ConvectionDispersionOperatorBase _convDispOp; //!< Convection dispersion operator for interstitial volume transport

	Eigen::SparseMatrix<double, RowMajor> _jac; //!< Jacobian
	Eigen::SparseMatrix<double, RowMajor> _jacDisc; //!< Jacobian with time derivatives from BDF method
	Eigen::SparseMatrix<double, RowMajor> _jacInlet; //!< Jacobian inlet DOF block matrix connects inlet DOFs to first bulk cells

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
		inline int strideColNode() const CADET_NOEXCEPT { return static_cast<int>(_disc.nComp + _disc.strideBound); }
		inline int strideColCell() const CADET_NOEXCEPT { return static_cast<int>(_disc.nNodes * strideColNode()); }
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
		template <typename real_t> inline real_t& c(real_t* const data, unsigned int node, unsigned int comp) const { return data[offsetC() + comp * strideColComp() + node * strideColNode()]; }
		template <typename real_t> inline const real_t& c(real_t const* const data, unsigned int node, unsigned int comp) const { return data[offsetC() + comp * strideColComp() + node * strideColNode()]; }

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
		// @TODO?? actually we need number of axial discrete points here, not number of axial cells !
		virtual unsigned int numAxialCells() const CADET_NOEXCEPT { return _disc.nCol * _disc.nNodes; }
		virtual unsigned int numRadialCells() const CADET_NOEXCEPT { return 0; }
		virtual unsigned int numInletPorts() const CADET_NOEXCEPT { return 1; }
		virtual unsigned int numOutletPorts() const CADET_NOEXCEPT { return 1; }
		virtual unsigned int numParticleTypes() const CADET_NOEXCEPT { return 1u; }
		virtual unsigned int numParticleShells(unsigned int parType) const CADET_NOEXCEPT { return 0u; }
		virtual unsigned int numBoundStates(unsigned int parType) const CADET_NOEXCEPT { return _disc.strideBound; }
		virtual unsigned int numBulkDofs() const CADET_NOEXCEPT { return _disc.nComp * _disc.nPoints; }
		virtual unsigned int numParticleMobilePhaseDofs(unsigned int parType) const CADET_NOEXCEPT { return 0u; }
		virtual unsigned int numSolidPhaseDofs(unsigned int parType) const CADET_NOEXCEPT { return _disc.strideBound * _disc.nPoints; }
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
				return &_idx.c(_data, _disc.nPoints - 1, 0);
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

		virtual unsigned int bulkMobilePhaseStride() const { return _idx.strideColNode(); }
		virtual unsigned int particleMobilePhaseStride(unsigned int parType) const { return 0; }
		virtual unsigned int solidPhaseStride(unsigned int parType) const { return _idx.strideColNode(); }

		/**
		* @brief calculates the physical node coordinates of the DG discretization with double! interface nodes
		*/
		virtual void axialCoordinates(double* coords) const {
			VectorXd x_l = VectorXd::LinSpaced(_disc.nCol + 1u, 0.0, _disc.length_);
			for (unsigned int i = 0; i < _disc.nCol; i++) {
				for (unsigned int j = 0; j < _disc.nNodes; j++) {
					// mapping 
					coords[i * _disc.nNodes + j] = x_l[i] + 0.5 * (_disc.length_ / static_cast<double>(_disc.nCol)) * (1.0 + _disc.nodes[j]);
				}
			}
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

	/**
	* @brief sets the current section index and section dependend velocity, dispersion
	*/
	void updateSection(int secIdx) {

		if (_disc.curSection != secIdx) {

			_disc.curSection = secIdx;
			_disc.newStaticJac = true;

			// update velocity and dispersion
			_disc.velocity = static_cast<double>(_convDispOp.currentVelocity());
			if (_convDispOp.dispersionCompIndep())
				for (int comp = 0; comp < _disc.nComp; comp++) {
					_disc.dispersion[comp] = static_cast<double>(_convDispOp.currentDispersion(secIdx)[0]);
				}
			else {
				for (int comp = 0; comp < _disc.nComp; comp++) {
					_disc.dispersion[comp] = static_cast<double>(_convDispOp.currentDispersion(secIdx)[comp]);
				}
			}

			//std::cout << "NEW SECTION: " << _disc.curSection << std::endl;
			//std::cout << "v: " << _disc.velocity << std::endl;
			//std::cout << "D_ax: " << _disc.dispersion << std::endl;
		}
	}

	// ==========================================================================================================================================================  //
	// ========================================						DG Residual, RHS						======================================================  //
	// ==========================================================================================================================================================  //

	/**
	* @brief calculates the volume Integral of the auxiliary equation
	* @param [in] current state vector
	* @param [in] stateDer vector to be changed
	* @param [in] aux true if auxiliary, else main equation
	*/
	void volumeIntegral(Eigen::Map<const VectorXd, 0, InnerStride<Dynamic>>& state, Eigen::Map<VectorXd, 0, InnerStride<Dynamic>>& stateDer) {
		// comp-cell-node state vector: use of Eigen lib performance
		for (unsigned int Cell = 0; Cell < _disc.nCol; Cell++) {
			stateDer.segment(Cell * _disc.nNodes, _disc.nNodes)
				-= _disc.polyDerM * state.segment(Cell * _disc.nNodes, _disc.nNodes);
		}
	}

	/*
	* @brief calculates the interface fluxes h* of Convection Dispersion equation
	*/
	void InterfaceFlux(Eigen::Map<const VectorXd, 0, InnerStride<Dynamic>>& C, const VectorXd& g, unsigned int comp) {

		// component-wise strides
		unsigned int strideCell = _disc.nNodes;
		unsigned int strideNode = 1u;

		// Conv.Disp. flux: h* = h*_conv + h*_disp = numFlux(v c_l, v c_r) + 0.5 sqrt(D_ax) (S_l + S_r)

		// calculate inner interface fluxes
		for (unsigned int Cell = 1; Cell < _disc.nCol; Cell++) {
			// h* = h*_conv + h*_disp
			_disc.surfaceFlux[Cell] // inner interfaces
				= _disc.velocity * (C[Cell * strideCell - strideNode])
				- 0.5 * std::sqrt(_disc.dispersion[comp]) * (g[Cell * strideCell - strideNode] // left cell
														       + g[Cell * strideCell]);
		}

		// boundary fluxes
			// left boundary interface
		_disc.surfaceFlux[0]
			= _disc.velocity * _disc.boundary[0];

		// right boundary interface
		_disc.surfaceFlux[_disc.nCol]
			= _disc.velocity * (C[_disc.nCol * strideCell - strideNode])
			- std::sqrt(_disc.dispersion[comp]) * 0.5 * (g[_disc.nCol * strideCell - strideNode] // last cell last node
														   + _disc.boundary[3]); // right boundary value S
	}

	/**
	* @brief calculates and fills the surface flux values for auxiliary equation
	*/
	void InterfaceFluxAuxiliary(Eigen::Map<const VectorXd, 0, InnerStride<Dynamic>>& C) {

		// component-wise strides
		unsigned int strideCell = _disc.nNodes;
		unsigned int strideNode = 1u;

		// Auxiliary flux: c* = 0.5 (c_l + c_r)

		// calculate inner interface fluxes
		for (unsigned int Cell = 1; Cell < _disc.nCol; Cell++) {
			_disc.surfaceFlux[Cell] // left interfaces
				= 0.5 * (C[Cell * strideCell - strideNode] + // left node
						 C[Cell * strideCell]); // right node
		}
		// calculate boundary interface fluxes

		_disc.surfaceFlux[0] // left boundary interface
			= 0.5 * (C[0] + // boundary value
					 C[0]); // first cell first node

		_disc.surfaceFlux[(_disc.nCol)] // right boundary interface
			= 0.5 * (C[_disc.nCol * strideCell - strideNode] + // last cell last node
					 C[_disc.nCol * strideCell - strideNode]);// // boundary value
	}

	/**
	* @brief calculates the surface Integral, depending on the approach (modal/nodal)
	* @param [in] state relevant state vector
	* @param [in] stateDer state derivative vector the solution is added to
	* @param [in] aux true for auxiliary equation, false for main equation
		surfaceIntegral(cPtr, &(disc.g[0]), disc,&(disc.h[0]), resPtrC, 0, secIdx);
	*/
	void surfaceIntegral(Eigen::Map<const VectorXd, 0, InnerStride<Dynamic>>& C, Eigen::Map<const VectorXd, 0, InnerStride<Dynamic>>& state, Eigen::Map<VectorXd, 0, InnerStride<Dynamic>>& stateDer, bool aux, unsigned int Comp) {

		// component-wise strides
		unsigned int strideCell = _disc.nNodes;
		unsigned int strideNode = 1u;

		// calc numerical flux values c* or h* depending on equation switch aux
		(aux == 1) ? InterfaceFluxAuxiliary(C) : InterfaceFlux(C, _disc.g, Comp);
		if (_disc.modal) { // modal approach -> dense mass matrix
			for (unsigned int Cell = 0; Cell < _disc.nCol; Cell++) {
				// strong surface integral -> M^-1 B [state - state*]
				for (unsigned int Node = 0; Node < _disc.nNodes; Node++) {
					stateDer[Cell * strideCell + Node * strideNode]
						-= _disc.invMM(Node, 0) * (state[Cell * strideCell]
												   - _disc.surfaceFlux[Cell])
						- _disc.invMM(Node, _disc.polyDeg) * (state[Cell * strideCell + _disc.polyDeg * strideNode]
															  - _disc.surfaceFlux[(Cell + 1)]);
				}
			}
		}
		else { // nodal approach -> diagonal mass matrix
			for (unsigned int Cell = 0; Cell < _disc.nCol; Cell++) {
				// strong surface integral -> M^-1 B [state - state*]
				stateDer[Cell * strideCell] // first cell node
					-= _disc.invWeights[0] * (state[Cell * strideCell] // first node
											  - _disc.surfaceFlux(Cell));
				stateDer[Cell * strideCell + _disc.polyDeg * strideNode] // last cell node
					+= _disc.invWeights[_disc.polyDeg] * (state[Cell * strideCell + _disc.polyDeg * strideNode]
														  - _disc.surfaceFlux(Cell + 1));
			}
		}
	}

	/**
	* @brief calculates the substitute h = vc - sqrt(D_ax) g(c)
	*/
	void calcH(Eigen::Map<const VectorXd, 0, InnerStride<>>& C, unsigned int Comp) {
		_disc.h = _disc.velocity * C - std::sqrt(_disc.dispersion[Comp]) * _disc.g;
	}

	/**
	* @brief calculates the isotherm right hand side
	*/
	void calcRHSq_DG(double t, unsigned int secIdx, const double* yPtr, double* const resPtr, util::ThreadLocalStorage& threadLocalMem) {

		Indexer idx(_disc);

		const double* localC = yPtr + idx.offsetC();
		const double* localQ = yPtr + idx.offsetC() + idx.strideColLiquid() + idx.offsetBoundComp(0);
		double* localQRes = resPtr + idx.offsetC() + idx.strideColLiquid() + idx.offsetBoundComp(0);

		for (unsigned int point = 0; point < _disc.nPoints; point++) {

			double z = _disc.deltaZ * std::floor(point / _disc.nNodes)
				+ 0.5 * _disc.deltaZ * (1 + _disc.nodes[point % _disc.nNodes]);

			_binding[0]->flux(t, secIdx, ColumnPosition{ z, 0.0, 0.0 }, localQ, localC, localQRes, threadLocalMem.get());

			localC += idx.strideColNode(); // next solid concentration
			localQ += idx.strideColNode(); // next liquid concentration
			localQRes += idx.strideColNode(); // next liquid concentration

		}
	}

	/**
	* @brief applies the inverse Jacobian of the mapping
	*/
	void applyMapping(Eigen::Map<VectorXd, 0, InnerStride<>>& state) {
		state *= (2.0 / _disc.deltaZ);
	}
	/**
	* @brief applies the inverse Jacobian of the mapping and auxiliary factor -1
	*/
	void applyMapping_Aux(Eigen::Map<VectorXd, 0, InnerStride<>>& state, unsigned int Comp) {
		state *= (-2.0 / _disc.deltaZ) * ((_disc.dispersion[Comp] == 0.0) ? 1.0 : std::sqrt(_disc.dispersion[Comp]));
	}

	void ConvDisp_DG(Eigen::Map<const VectorXd, 0, InnerStride<Dynamic>>& C, Eigen::Map<VectorXd, 0, InnerStride<Dynamic>>& resC, double t, unsigned int Comp) {

		// ===================================//
		// reset cache                        //
		// ===================================//

		resC.setZero();
		_disc.h.setZero();
		_disc.g.setZero();
		_disc.surfaceFlux.setZero();
		// get Map objects of auxiliary variable memory
		Eigen::Map<VectorXd, 0, InnerStride<>> g(&_disc.g[0], _disc.nPoints, InnerStride<>(1));
		Eigen::Map<const VectorXd, 0, InnerStride<>> h(&_disc.h[0], _disc.nPoints, InnerStride<>(1));

		// ======================================//
		// solve auxiliary system g = d c / d x  //
		// ======================================//

		volumeIntegral(C, g); // DG volumne integral in strong form

		surfaceIntegral(C, C, g, 1, Comp); // surface integral in strong form

		applyMapping_Aux(g, Comp); // inverse mapping from reference space and auxiliary factor

		_disc.surfaceFlux.setZero(); // reset surface flux storage as it is used twice

		// ======================================//
		// solve main equation w_t = d h / d x   //
		// ======================================//

		calcH(C, Comp); // calculate the substitute h(S(c), c) = sqrt(D_ax) g(c) - v c

		volumeIntegral(h, resC); // DG volumne integral in strong form

		calcBoundaryValues(C);// update boundary values including auxiliary variable g

		surfaceIntegral(C, h, resC, 0, Comp); // DG surface integral in strong form

		applyMapping(resC); // inverse mapping to reference space

	}

	/**
	* @brief computes ghost nodes used to implement Danckwerts boundary conditions
	*/
	void calcBoundaryValues(Eigen::Map<const VectorXd, 0, InnerStride<>>& C) {

		//cache.boundary[0] = c_in -> inlet DOF idas suggestion
		_disc.boundary[1] = C[_disc.nPoints - 1]; // c_r outlet
		_disc.boundary[2] = -_disc.g[0]; // S_l inlet
		_disc.boundary[3] = -_disc.g[_disc.nPoints - 1]; // g_r outlet
	}

	// ==========================================================================================================================================================  //
	// ========================================						DG Jacobian							=========================================================  //
	// ==========================================================================================================================================================  //

	typedef Eigen::Triplet<double> T;

	/**
	* @brief sets the sparsity pattern of the static Jacobian
	* @detail inlet DOFs plus ConvDisp pattern DG and isotherm pattern. Independent of the isotherm,
		all liquid and solid entries at a discrete point are set.
	* @param [in] stateDer bool if state derivative pattern should be added (for _jacDisc)
	*/
	void setPattern(Eigen::SparseMatrix<double, RowMajor>& mat, bool stateDer) {

		std::vector<T> tripletList;
		// TODO?: convDisp NNZ times two for now, but Convection NNZ < Dispersion NNZ
		// inlet + ConvDispDG + isotherm
		tripletList.reserve(_disc.nComp + 2u * calcConvDispNNZ(_disc) + (_disc.strideBound + _disc.nComp) * _disc.nPoints);

		// inlet DOFs Jacobian ( forward flow! ) //_jacInlet;
		for (unsigned int i = 0; i < _disc.nComp; i++) {
			tripletList.push_back(T(i, i, 0.0));
		}

		if (_disc.modal)
			ConvDispModalPattern(tripletList);
		else
			ConvDispNodalPattern(tripletList);

		isothermPattern(tripletList);

		if (stateDer)
			stateDerPattern(tripletList); // only adds [ d convDisp / d q_t ] because main diagonal is already included !

		mat.setFromTriplets(tripletList.begin(), tripletList.end());

	}

	/**
	* @brief computes the convection dispersion part of the state derivative pattern of the jacobian.
	* @detail the main also diagonal belongs to the state derivative pattern but is not set here, so it should be set previously.
	*/
	void stateDerPattern(std::vector<T>& tripletList) {

		Indexer idxr(_disc);
		unsigned int offC = idxr.offsetC();

		for (unsigned int point = 0; point < _disc.nPoints; point++) {
			for (unsigned int comp = 0; comp < _disc.nComp; comp++) {
				if (_disc.nBound[comp]) // either one or null
					// row:    jump over inlet, go node and comp strides to find all liquid states
					// column: jump over inlet and liquid states, add component offset and go node strides to find corresponding solid states
					tripletList.push_back(T(offC + point * idxr.strideColNode() + comp * idxr.strideColComp(),
											offC + idxr.strideColLiquid() + idxr.offsetBoundComp(comp) + point * idxr.strideColNode(),
											0.0));
			}
		}

	}

	/**
	* @brief sets the sparsity pattern of the convection dispersion Jacobian for the nodal DG scheme
	*/
	int ConvDispNodalPattern(std::vector<T>& tripletList) {

		Indexer idx(_disc);

		int sNode = idx.strideColNode();
		int sCell = idx.strideColCell();
		int sComp = idx.strideColComp();
		int offC = idx.offsetC();

		unsigned int nNodes = _disc.nNodes;
		unsigned int polyDeg = _disc.polyDeg;
		unsigned int nCells = _disc.nCol;
		unsigned int nComp = _disc.nComp;

		// @TODO: special cases?
		if (nCells < 3)
			throw std::invalid_argument("Nodal Jacobian special case for nCells < 3 not implemented (yet?)");

		/*======================================================*/
		/*			Define Convection Jacobian Block			*/
		/*======================================================*/

		// Convection block [ d RHS_conv / d c ], also depends on first entry of previous cell

		// special inlet DOF treatment for first cell
		for (unsigned int comp = 0; comp < nComp; comp++) {
			for (unsigned int i = 0; i < nNodes; i++) {
				tripletList.push_back(T(offC + comp * sComp + i * sNode, comp * sComp, 0.0));
				for (unsigned int j = 1; j < nNodes + 1; j++) {
					tripletList.push_back(T(offC + comp * sComp + i * sNode,
											offC + comp * sComp + (j - 1) * sNode,
											0.0));
				}
			}
		}
		for (unsigned int cell = 1; cell < nCells; cell++) {
			for (unsigned int comp = 0; comp < nComp; comp++) {
				for (unsigned int i = 0; i < nNodes; i++) {
					for (unsigned int j = 0; j < nNodes + 1; j++) {
						// row: jump over inlet DOFs and previous cells, add component offset and go node strides from there for each convection block entry
						// col: jump over inlet DOFs and previous cells, go back one node, add component offset and go node strides from there for each convection block entry
						tripletList.push_back(T(offC + cell * sCell + comp * sComp + i * sNode,
												offC + cell * sCell - sNode + comp * sComp + j * sNode,
												0.0));
					}
				}
			}
		}

		/*======================================================*/
		/*			Define Dispersion Jacobian Block			*/
		/*======================================================*/

		/*		Inner cell dispersion blocks		*/


		// Dispersion block [ d RHS_disp / d c ], depends on whole previous and subsequent cell

		// insert Blocks to Jacobian inner cells
		for (unsigned int cell = 1; cell < nCells - 1; cell++) {
			for (unsigned int comp = 0; comp < nComp; comp++) {
				for (unsigned int i = 0; i < nNodes; i++) {
					for (unsigned int j = 0; j < 3 * nNodes; j++) {
						// row: jump over inlet DOFs and previous cells, add component offset and go node strides from there for each dispersion block entry
						// col: jump over inlet DOFs and previous cells, go back one cell, add component offset and go node strides from there for each dispersion block entry
						tripletList.push_back(T(offC + cell * sCell + comp * sComp + i * sNode,
												offC + (cell - 1) * sCell + comp * sComp + j * sNode,
												0.0));
					}
				}
			}
		}

		/*				Boundary cell Dispersion blocks			*/

		for (unsigned int comp = 0; comp < nComp; comp++) {
			for (unsigned int i = 0; i < nNodes; i++) {
				for (unsigned int j = nNodes; j < 3 * nNodes; j++) {
					tripletList.push_back(T(offC + comp * sComp + i * sNode,
											offC + comp * sComp + (j - nNodes) * sNode,
											0.0));
				}
			}
		}

		for (unsigned int comp = 0; comp < nComp; comp++) {
			for (unsigned int i = 0; i < nNodes; i++) {
				for (unsigned int j = 0; j < 2 * nNodes; j++) {
					tripletList.push_back(T(offC + (nCells - 1) * sCell + comp * sComp + i * sNode,
											offC + (nCells - 1 - 1) * sCell + comp * sComp + j * sNode,
											0.0));
				}
			}
		}

		return 0;
	}
			
	/**
	* @brief sets the sparsity pattern of the convection dispersion Jacobian for the modal DG scheme
	*/
	int ConvDispModalPattern(std::vector<T>& tripletList) {

		Indexer idx(_disc);

		int sNode = idx.strideColNode();
		int sCell = idx.strideColCell();
		int sComp = idx.strideColComp();
		int offC = idx.offsetC();

		unsigned int nNodes = _disc.nNodes;
		unsigned int nCells = _disc.nCol;
		unsigned int nComp = _disc.nComp;

		// @TODO: special cases?
		if (nCells < 5)
			throw std::invalid_argument("Modal Jacobian special case for nCells < 5 not implemented (yet?)");

		/*======================================================*/
		/*			Define Convection Jacobian Block			*/
		/*======================================================*/

		// Convection block [ d RHS_conv / d c ], additionally depends on first entry of previous cell
		//MatrixXd convBlock = MatrixXd::Zero(nNodes, nNodes + 1);
		// special inlet DOF treatment for first cell
		for (unsigned int comp = 0; comp < nComp; comp++) {
			for (unsigned int i = 0; i < nNodes; i++) {
				tripletList.push_back(T(offC + comp * sComp + i * sNode, comp * sComp, 0.0));
				for (unsigned int j = 1; j < nNodes + 1; j++) {
					tripletList.push_back(T(offC + comp * sComp + i * sNode,
											offC + comp * sComp + (j - 1) * sNode,
											0.0));
				}
			}
		}
		for (unsigned int cell = 1; cell < nCells; cell++) {
			for (unsigned int comp = 0; comp < nComp; comp++) {
				for (unsigned int i = 0; i < nNodes; i++) {
					for (unsigned int j = 0; j < nNodes + 1; j++) {
						// row: jump over inlet DOFs and previous cells, add component offset and go node strides from there for each convection block entry
						// col: jump over inlet DOFs and previous cells, go back one node, add component offset and go node strides from there for each convection block entry
						tripletList.push_back(T(offC + cell * sCell + comp * sComp + i * sNode,
												offC + cell * sCell - sNode + comp * sComp + j * sNode,
												0.0));
					}
				}
			}
		}

		/*======================================================*/
		/*			Define Dispersion Jacobian Block			*/
		/*======================================================*/

		/* Inner cells */

		// Inner dispersion block [ d RHS_disp / d c ], depends on whole previous and subsequent cell plus first entries of subsubsequent cells
		//MatrixXd dispBlock = MatrixXd::Zero(nNodes, 3 * nNodes + 2);
		for (unsigned int cell = 2; cell < nCells - 2; cell++) {
			for (unsigned int comp = 0; comp < nComp; comp++) {
				for (unsigned int i = 0; i < nNodes; i++) {
					for (unsigned int j = 0; j < 3 * nNodes + 2; j++) {
						// row: jump over inlet DOFs and previous cells, add component offset and go node strides from there for each dispersion block entry
						// col: jump over inlet DOFs and previous cells, go back one cell and one node, add component offset and go node strides from there for each dispersion block entry
						tripletList.push_back(T(offC + cell * sCell + comp * sComp + i * sNode,
												offC + cell * sCell - (nNodes + 1) * sNode + comp * sComp + j * sNode,
												0.0));
					}
				}
			}
		}

		/*		boundary cell neighbours		*/

		// left boundary cell neighbour
		for (unsigned int comp = 0; comp < nComp; comp++) {
			for (unsigned int i = 0; i < nNodes; i++) {
				for (unsigned int j = 1; j < 3 * nNodes + 2; j++) {
					// row: jump over inlet DOFs and previous cell, add component offset and go node strides from there for each dispersion block entry
					// col: jump over inlet DOFs, add component offset and go node strides from there for each dispersion block entry. Also adjust for iterator j (-1)
					tripletList.push_back(T(offC + nNodes * sNode + comp * sComp + i * sNode,
											offC + comp * sComp + (j - 1) * sNode,
											0.0));
				}
			}
		}
		// right boundary cell neighbour
		for (unsigned int comp = 0; comp < nComp; comp++) {
			for (unsigned int i = 0; i < nNodes; i++) {
				for (unsigned int j = 0; j < 3 * nNodes + 2 - 1; j++) {
					// row: jump over inlet DOFs and previous cells, add component offset and go node strides from there for each dispersion block entry
					// col: jump over inlet DOFs and previous cells, go back one cell and one node, add component offset and go node strides from there for each dispersion block entry.
					tripletList.push_back(T(offC + (nCells - 2) * sCell + comp * sComp + i * sNode,
											offC + (nCells - 2) * sCell - (nNodes + 1) * sNode + comp * sComp + j * sNode,
											0.0));
				}
			}
		}

		/*			boundary cells			*/

		// left boundary cell
		for (unsigned int comp = 0; comp < nComp; comp++) {
			for (unsigned int i = 0; i < nNodes; i++) {
				for (unsigned int j = nNodes + 1; j < 3 * nNodes + 2; j++) {
					// row: jump over inlet DOFs, add component offset and go node strides from there for each dispersion block entry
					// col: jump over inlet DOFs and previous cells, add component offset, adjust for iterator j (-Nnodes-1) and go node strides from there for each dispersion block entry.
					tripletList.push_back(T(offC + comp * sComp + i * sNode,
											offC + comp * sComp + (j - (nNodes + 1)) * sNode,
											0.0));
				}
			}
		}
		// right boundary cell
		for (unsigned int comp = 0; comp < nComp; comp++) {
			for (unsigned int i = 0; i < nNodes; i++) {
				for (unsigned int j = 0; j < 2 * nNodes + 1; j++) {
					// row: jump over inlet DOFs and previous cells, add component offset and go node strides from there for each dispersion block entry
					// col: jump over inlet DOFs and previous cells, go back one cell and one node, add component offset and go node strides from there for each dispersion block entry.
					tripletList.push_back(T(offC + (nCells - 1) * sCell + comp * sComp + i * sNode,
											offC + (nCells - 1) * sCell - (nNodes + 1) * sNode + comp * sComp + j * sNode,
											0.0));
				}
			}
		}

		return 0;
	}

	/**
	* @brief sets the sparsity pattern of the isotherm Jacobian
	* @detail Independent of the isotherm, all liquid and solid entries (so all entries, the isotherm could theoretically depend on) at a discrete point are set.
	*/
	void isothermPattern(std::vector<T>& tripletList) {
		
		Indexer idxr(_disc);
		// loop over all discrete points and solid states and add all liquid plus solid entries at that solid state at that discrete point
		for (unsigned int point = 0; point < _disc.nPoints; point++) {
			for (unsigned int solid = 0; solid < _disc.strideBound; solid++) {
				for (unsigned int conc = 0; conc < _disc.nComp + _disc.strideBound; conc++) {
					// column: jump over inlet, previous discrete points, liquid concentration and add the offset of the current bound state
					// row:    jump over inlet and previous discrete points. add entries for all liquid and bound concentrations (conc)
					tripletList.push_back(T(idxr.offsetC() + idxr.strideColNode() * point + idxr.strideColLiquid() + solid * idxr.offsetBoundComp(solid),
											idxr.offsetC() + idxr.strideColNode() * point + conc,
											0.0));
				}
			}
		}

	}

	/**
	* @brief analytically calculates the (static) state jacobian
	* @return 1 if jacobain estimation fits the predefined pattern of the jacobian, 0 if not.
	*/
	int calcStaticAnaJacobian(double t, unsigned int secIdx, const double* const y, util::ThreadLocalStorage& threadLocalMem) {

		// not necessary: set _jac to zero but keep pattern ! (.setZero() deletes pattern)
		//double* vPtr = _jac.valuePtr();
		//for (int k = 0; k < _jac.nonZeros(); k++) {
		//	vPtr[k] = 0.0;
		//}

			// DG convection dispersion Jacobian
		if (_disc.modal)
			calcConvDispModalJacobian();
		else
			calcConvDispNodalJacobian();

		// inlet Jacobian (forward flow)
		double* valPtr = _jac.valuePtr();
		for (unsigned int i = 0; i < _disc.nComp; i++) {
			valPtr[i] = 1.0;
		}


		if (!_jac.isCompressed()) // if matrix lost its compressed storage, the pattern did not fit.
			return 0;

		//MatrixXd mat = _jac.toDense()*10;
		//for (int i = 0; i < mat.rows(); i++) {
		//	for (int j = 0; j < mat.cols(); j++) {
		//		if (mat(i, j) != 0) {
		//			mat(i, j) = 1.0;
		//		}
		//	}
		//}
		//std::cout << std::fixed << std::setprecision(0) << "JAC\n" << mat << std::endl; // #include <iomanip>

		return 1;
	}

	/**
	* @brief calculates the number of non zeros for the DG convection dispersion jacobian
	* @detail only dispersion entries are relevant as the convection entries are a subset of these
	*/
	unsigned int calcConvDispNNZ(Discretization disc) {

		if (disc.modal) {
			return disc.nComp * ((3u * disc.nCol - 2u) * disc.nNodes * disc.nNodes + (2u * disc.nCol - 3u) * disc.nNodes);
		}
		else {
			return disc.nComp * (disc.nCol * disc.nNodes * disc.nNodes + 8u * disc.nNodes);
		}
	}

	/**
	* @brief analytically calculates the convection dispersion jacobian for the nodal DG scheme
	*/
	int calcConvDispNodalJacobian() {

		Indexer idx(_disc);

		int sNode = idx.strideColNode();
		int sCell = idx.strideColCell();
		int sComp = idx.strideColComp();
		int offC = idx.offsetC();

		unsigned int nNodes = _disc.nNodes;
		unsigned int polyDeg = _disc.polyDeg;
		unsigned int nCells = _disc.nCol;
		unsigned int nComp = _disc.nComp;

		// @TODO: special cases?
		if (nCells < 3)
			throw std::invalid_argument("Nodal Jacobian special case for nCells < 3 not implemented (yet?)");

		/*======================================================*/
		/*			Compute Dispersion Jacobian Block			*/
		/*======================================================*/

		/*		Inner cell dispersion blocks		*/

		// auxiliary Block for [ d g(c) / d c ], needed in Dispersion block
		MatrixXd GBlock = MatrixXd::Zero(nNodes, nNodes + 2);
		GBlock.block(0, 1, nNodes, nNodes) = _disc.polyDerM;
		GBlock(0, 0) -= 0.5 * _disc.invWeights[0];
		GBlock(0, 1) += 0.5 * _disc.invWeights[0];
		GBlock(nNodes - 1, nNodes) -= 0.5 * _disc.invWeights[polyDeg];
		GBlock(nNodes - 1, nNodes + 1) += 0.5 * _disc.invWeights[polyDeg];
		GBlock *= 2 / _disc.deltaZ;

		// Dispersion block [ d RHS_disp / d c ], depends on whole previous and subsequent cell
		MatrixXd dispBlock = MatrixXd::Zero(nNodes, 3 * nNodes); //
		// NOTE: N = polyDeg
		//cell indices : 0	 , ..., nNodes - 1;	nNodes, ..., 2 * nNodes - 1;	2 * nNodes, ..., 3 * nNodes - 1
		//			j  : -N-1, ..., -1		  ; 0     , ..., N			   ;	N + 1, ..., 2N + 1
		dispBlock.block(0, nNodes - 1, nNodes, nNodes + 2) = _disc.polyDerM * GBlock;
		dispBlock(0, nNodes - 1) += -_disc.invWeights[0] * (-0.5 * GBlock(0, 0) + 0.5 * GBlock(nNodes - 1, nNodes)); // G_N,N		i=0, j=-1
		dispBlock(0, nNodes) += -_disc.invWeights[0] * (-0.5 * GBlock(0, 1) + 0.5 * GBlock(nNodes - 1, nNodes + 1)); // G_N,N+1	i=0, j=0
		dispBlock.block(0, nNodes + 1, 1, nNodes) += -_disc.invWeights[0] * (-0.5 * GBlock.block(0, 2, 1, nNodes)); // G_i,j		i=0, j=1,...,N+1
		dispBlock.block(0, 0, 1, nNodes - 1) += -_disc.invWeights[0] * (0.5 * GBlock.block(nNodes - 1, 1, 1, nNodes - 1)); // G_N,j+N+1		i=0, j=-N-1,...,-2
		dispBlock.block(nNodes - 1, nNodes - 1, 1, nNodes) += _disc.invWeights[nNodes - 1] * (-0.5 * GBlock.block(nNodes - 1, 0, 1, nNodes)); // G_i,j+N+1		i=N, j=--1,...,N-1
		dispBlock(nNodes - 1, 2 * nNodes - 1) += _disc.invWeights[nNodes - 1] * (-0.5 * GBlock(nNodes - 1, nNodes) + 0.5 * GBlock(0, 0)); // G_i,j		i=N, j=N
		dispBlock(nNodes - 1, 2 * nNodes) += _disc.invWeights[nNodes - 1] * (-0.5 * GBlock(nNodes - 1, nNodes + 1) + 0.5 * GBlock(0, 1)); // G_i,j		i=N, j=N+1
		dispBlock.block(nNodes - 1, 2 * nNodes + 1, 1, nNodes - 1) += _disc.invWeights[nNodes - 1] * (0.5 * GBlock.block(0, 2, 1, nNodes - 1)); // G_0,j-N-1		i=N, j=N+2,...,2N+1
		dispBlock *= 2 / _disc.deltaZ;

		// insert Blocks to Jacobian inner cells
		for (unsigned int cell = 1; cell < nCells - 1; cell++) {
			for (unsigned int comp = 0; comp < nComp; comp++) {
				for (unsigned int i = 0; i < dispBlock.rows(); i++) {
					for (unsigned int j = 0; j < dispBlock.cols(); j++) {
						// row: jump over inlet DOFs and previous cells, add component offset and go node strides from there for each dispersion block entry
						// col: jump over inlet DOFs and previous cells, go back one cell, add component offset and go node strides from there for each dispersion block entry
						_jac.coeffRef(offC + cell * sCell + comp * sComp + i * sNode,
									  offC + (cell - 1) * sCell + comp * sComp + j * sNode)
									  = -dispBlock(i, j) * _disc.dispersion[comp];
					}
				}
			}
		}

		/*				Boundary cell Dispersion blocks			*/

		/* left cell */
		// adjust auxiliary Block [ d g(c) / d c ] for left boundary cell
		MatrixXd GBlockBound = GBlock;
		GBlockBound(0, 1) -= 0.5 * _disc.invWeights[0] * 2 / _disc.deltaZ;

		// estimate dispersion block ( j < 0 not needed)
		dispBlock.setZero();
		dispBlock.block(0, nNodes - 1, nNodes, nNodes + 2) = _disc.polyDerM * GBlockBound;
		dispBlock.block(0, nNodes - 1, 1, nNodes + 2) += -_disc.invWeights[0] * (-GBlockBound.block(0, 0, 1, nNodes + 2)); // G_N,N		i=0, j=-1,...,N+1
		dispBlock.block(nNodes - 1, nNodes - 1, 1, nNodes) += _disc.invWeights[nNodes - 1] * (-0.5 * GBlockBound.block(nNodes - 1, 0, 1, nNodes)); // G_i,j+N+1		i=N, j=--1,...,N-1
		dispBlock(nNodes - 1, 2 * nNodes - 1) += _disc.invWeights[nNodes - 1] * (-0.5 * GBlockBound(nNodes - 1, nNodes) + 0.5 * GBlockBound(0, 0)); // G_i,j		i=N, j=N
		dispBlock(nNodes - 1, 2 * nNodes) += _disc.invWeights[nNodes - 1] * (-0.5 * GBlockBound(nNodes - 1, nNodes + 1) + 0.5 * GBlock(0, 1)); // G_i,j		i=N, j=N+1
		dispBlock.block(nNodes - 1, 2 * nNodes + 1, 1, nNodes - 1) += _disc.invWeights[nNodes - 1] * (0.5 * GBlock.block(0, 2, 1, nNodes - 1)); // G_0,j-N-1		i=N, j=N+2,...,2N+1
		dispBlock *= 2 / _disc.deltaZ;
		// copy *-1 to Jacobian
		for (unsigned int comp = 0; comp < nComp; comp++) {
			for (unsigned int i = 0; i < dispBlock.rows(); i++) {
				for (unsigned int j = nNodes; j < dispBlock.cols(); j++) {
					_jac.coeffRef(offC + comp * sComp + i * sNode,
								  offC + comp * sComp + (j - nNodes) * sNode)
								  = -dispBlock(i, j) * _disc.dispersion[comp];
				}
			}
		}

		/* right cell */
	   // adjust auxiliary Block [ d g(c) / d c ] for left boundary cell
		GBlockBound(0, 1) += 0.5 * _disc.invWeights[0] * 2 / _disc.deltaZ; 	// reverse change from left boundary
		GBlockBound(nNodes - 1, nNodes) += 0.5 * _disc.invWeights[polyDeg] * 2 / _disc.deltaZ;

		// estimate dispersion block (only estimation differences to inner cell at N = 0 and j > N not needed)
		dispBlock.block(0, nNodes - 1, nNodes, nNodes + 2) = _disc.polyDerM * GBlockBound;
		dispBlock(0, nNodes - 1) += -_disc.invWeights[0] * (-0.5 * GBlockBound(0, 0) + 0.5 * GBlock(nNodes - 1, nNodes)); // G_N,N		i=0, j=-1
		dispBlock(0, nNodes) += -_disc.invWeights[0] * (-0.5 * GBlockBound(0, 1) + 0.5 * GBlock(nNodes - 1, nNodes + 1)); // G_N,N+1	i=0, j=0
		dispBlock.block(0, nNodes + 1, 1, nNodes) += -_disc.invWeights[0] * (-0.5 * GBlockBound.block(0, 2, 1, nNodes)); // G_i,j		i=0, j=1,...,N+1
		dispBlock.block(0, 0, 1, nNodes - 1) += -_disc.invWeights[0] * (0.5 * GBlock.block(nNodes - 1, 1, 1, nNodes - 1)); // G_N,j+N+1		i=0, j=-N-1,...,-2
		dispBlock.block(nNodes - 1, nNodes - 1, 1, nNodes + 2) += _disc.invWeights[nNodes - 1] * (-GBlockBound.block(nNodes - 1, 0, 1, nNodes + 2)); // G_i,j+N+1		i=N, j=--1,...,N+1
		dispBlock *= 2 / _disc.deltaZ;
		// copy *-1 to Jacobian
		for (unsigned int comp = 0; comp < nComp; comp++) {
			for (unsigned int i = 0; i < dispBlock.rows(); i++) {
				for (unsigned int j = 0; j < 2 * nNodes; j++) {
					_jac.coeffRef(offC + (nCells - 1) * sCell + comp * sComp + i * sNode,
								  offC + (nCells - 1 - 1) * sCell + comp * sComp + j * sNode)
								  = -dispBlock(i, j) * _disc.dispersion[comp];
				}
			}
		}

		/*======================================================*/
		/*			Compute Convection Jacobian Block			*/
		/*======================================================*/

		// Convection block [ d RHS_conv / d c ], also depends on first entry of previous cell
		MatrixXd convBlock = MatrixXd::Zero(nNodes, nNodes + 1);
		convBlock.block(0, 1, nNodes, nNodes) -= _disc.polyDerM;
		convBlock(0, 0) += _disc.invWeights[0];
		convBlock(0, 1) -= _disc.invWeights[0];
		convBlock *= 2 * _disc.velocity / _disc.deltaZ;

		// special inlet DOF treatment for first cell
		for (unsigned int comp = 0; comp < nComp; comp++) {
			for (unsigned int i = 0; i < convBlock.rows(); i++) {
				_jac.coeffRef(offC + comp * sComp + i * sNode, comp * sComp) = -convBlock(i, 0);
				for (unsigned int j = 1; j < convBlock.cols(); j++) {
					_jac.coeffRef(offC + comp * sComp + i * sNode,
								  offC + comp * sComp + (j - 1) * sNode)
								  += -convBlock(i, j);
				}
			}
		}
		for (unsigned int cell = 1; cell < nCells; cell++) {
			for (unsigned int comp = 0; comp < nComp; comp++) {
				for (unsigned int i = 0; i < convBlock.rows(); i++) {
					//Eigen::SparseMatrix<double, RowMajor>::InnerIterator it(_jac, offC + cell * sCell + comp * sComp + i * sNode);
					//it += _disc.polyDeg; // jump over dispersion block entries
					for (unsigned int j = 0; j < convBlock.cols(); j++) {
						// row: jump over inlet DOFs and previous cells, add component offset and go node strides from there for each convection block entry
						// col: jump over inlet DOFs and previous cells, go back one node, add component offset and go node strides from there for each convection block entry
						_jac.coeffRef(offC + cell * sCell + comp * sComp + i * sNode,
									  offC + cell * sCell - sNode + comp * sComp + j * sNode)
									  += -convBlock(i, j);
						//it.valueRef() += -convBlock(i, j);
						//++it;
					}
				}
			}
		}

		return 0;
	}
	/**
	* @brief analytically calculates the convection dispersion jacobian for the modal DG scheme
	*/
	int calcConvDispModalJacobian() {

		Indexer idx(_disc);

		int sNode = idx.strideColNode();
		int sCell = idx.strideColCell();
		int sComp = idx.strideColComp();
		int offC = idx.offsetC();

		unsigned int nNodes = _disc.nNodes;
		unsigned int nCells = _disc.nCol;
		unsigned int nComp = _disc.nComp;

		// @TODO: special cases?
		if (nCells < 5)
			throw std::invalid_argument("Modal Jacobian special case for nCells < 5 not implemented (yet?)");

		/*======================================================*/
		/*			Compute Dispersion Jacobian Block			*/
		/*======================================================*/

		/* Inner cells */

		// Auxiliary Block [ d g(c) / d c ], additionally depends on boundary entries of neighbouring cells
		MatrixXd gBlock = MatrixXd::Zero(nNodes, nNodes + 2);
		gBlock.block(0, 1, nNodes, nNodes) = _disc.polyDerM;
		gBlock.block(0, 0, nNodes, 1) -= 0.5 * _disc.invMM.block(0, 0, nNodes, 1);
		gBlock.block(0, 1, nNodes, 1) += 0.5 * _disc.invMM.block(0, 0, nNodes, 1);
		gBlock.block(0, nNodes, nNodes, 1) -= 0.5 * _disc.invMM.block(0, nNodes - 1, nNodes, 1);
		gBlock.block(0, nNodes + 1, nNodes, 1) += 0.5 * _disc.invMM.block(0, nNodes - 1, nNodes, 1);
		gBlock *= 2 / _disc.deltaZ;

		// B matrix from DG scheme
		MatrixXd B = MatrixXd::Zero(nNodes, nNodes);
		B(0, 0) = -1.0;
		B(nNodes - 1, nNodes - 1) = 1.0;

		// Inner dispersion block [ d RHS_disp / d c ], depends on whole previous and subsequent cell plus first entries of subsubsequent cells
		MatrixXd dispBlock = MatrixXd::Zero(nNodes, 3 * nNodes + 2); //
		// auxiliary block [ d g^* / d c ], depends on whole previous and subsequent cell plus first entries of subsubsequent cells
		MatrixXd gStarDC = MatrixXd::Zero(nNodes, 3 * nNodes + 2);
		// NOTE: N = polyDeg
		//  indices  gStarDC    :     0   ,   1   , ..., nNodes; nNodes+1, ..., 2 * nNodes;	2*nNodes+1, ..., 3 * nNodes; 3*nNodes+1
		//	derivative index j  : -(N+1)-1, -(N+1),... ,  -1   ;   0     , ...,		N	 ;	  N + 1	  , ..., 2N + 2    ; 2(N+1) +1
		// auxiliary block [d g^* / d c]
		gStarDC.block(0, nNodes, 1, nNodes + 2) += gBlock.block(0, 0, 1, nNodes + 2);
		gStarDC.block(0, 0, 1, nNodes + 2) += gBlock.block(nNodes - 1, 0, 1, nNodes + 2);
		gStarDC.block(nNodes - 1, nNodes, 1, nNodes + 2) += gBlock.block(nNodes - 1, 0, 1, nNodes + 2);
		gStarDC.block(nNodes - 1, 2 * nNodes, 1, nNodes + 2) += gBlock.block(0, 0, 1, nNodes + 2);
		gStarDC *= 0.5;

		//  indices  dispBlock :   0	 ,   1   , ..., nNodes;	nNodes+1, ..., 2 * nNodes;	2*nNodes+1, ..., 3 * nNodes; 3*nNodes+1
		//	derivative index j  : -(N+1)-1, -(N+1),...,	 -1	  ;   0     , ...,		N	 ;	  N + 1	  , ..., 2N + 2    ; 2(N+1) +1
		dispBlock.block(0, nNodes, nNodes, nNodes + 2) += _disc.polyDerM * gBlock - _disc.invMM * B * gBlock;
		dispBlock += _disc.invMM * B * gStarDC;
		dispBlock *= 2 / _disc.deltaZ;

		for (unsigned int cell = 2; cell < nCells - 2; cell++) {
			for (unsigned int comp = 0; comp < nComp; comp++) {
				for (unsigned int i = 0; i < dispBlock.rows(); i++) {
					for (unsigned int j = 0; j < dispBlock.cols(); j++) {
						// row: jump over inlet DOFs and previous cells, add component offset and go node strides from there for each dispersion block entry
						// col: jump over inlet DOFs and previous cells, go back one cell and one node, add component offset and go node strides from there for each dispersion block entry
						_jac.coeffRef(offC + cell * sCell + comp * sComp + i * sNode,
									  offC + cell * sCell - (nNodes + 1) * sNode + comp * sComp + j * sNode)
									  = -dispBlock(i, j) * _disc.dispersion[comp];
					}
				}
			}
		}

		/*		boundary cell neighbours		*/

		// left boundary cell neighbour

		// boundary auxiliary block [ d g(c) / d c ]
		MatrixXd GBlockBound_l = MatrixXd::Zero(nNodes, nNodes + 2);
		GBlockBound_l.block(0, 1, nNodes, nNodes) += _disc.polyDerM;
		GBlockBound_l.block(0, nNodes, nNodes, 1) -= 0.5 * _disc.invMM.block(0, nNodes - 1, nNodes, 1);
		GBlockBound_l.block(0, nNodes + 1, nNodes, 1) += 0.5 * _disc.invMM.block(0, nNodes - 1, nNodes, 1);
		GBlockBound_l *= 2 / _disc.deltaZ;
		// auxiliary block [d g^* / d c]
		gStarDC.setZero();
		gStarDC.block(0, nNodes, 1, nNodes + 2) += gBlock.block(0, 0, 1, nNodes + 2);
		gStarDC.block(0, 0, 1, nNodes + 2) += GBlockBound_l.block(nNodes - 1, 0, 1, nNodes + 2);
		gStarDC.block(nNodes - 1, nNodes, 1, nNodes + 2) += gBlock.block(nNodes - 1, 0, 1, nNodes + 2);
		gStarDC.block(nNodes - 1, 2 * nNodes, 1, nNodes + 2) += gBlock.block(0, 0, 1, nNodes + 2);
		gStarDC *= 0.5;

		dispBlock.setZero();
		dispBlock.block(0, nNodes, nNodes, nNodes + 2) += _disc.polyDerM * gBlock - _disc.invMM * B * gBlock;
		dispBlock += _disc.invMM * B * gStarDC;
		dispBlock *= 2 / _disc.deltaZ;

		for (unsigned int comp = 0; comp < nComp; comp++) {
			for (unsigned int i = 0; i < dispBlock.rows(); i++) {
				for (unsigned int j = 1; j < dispBlock.cols(); j++) {
					// row: jump over inlet DOFs and previous cell, add component offset and go node strides from there for each dispersion block entry
					// col: jump over inlet DOFs, add component offset and go node strides from there for each dispersion block entry. Also adjust for iterator j (-1)
					_jac.coeffRef(offC + nNodes * sNode + comp * sComp + i * sNode,
								  offC + comp * sComp + (j - 1) * sNode)
								  = -dispBlock(i, j) * _disc.dispersion[comp];
				}
			}
		}

		// right boundary cell neighbour
		// boundary auxiliary block [ d g(c) / d c ]
		MatrixXd GBlockBound_r = MatrixXd::Zero(nNodes, nNodes + 2);
		GBlockBound_r.block(0, 1, nNodes, nNodes) += _disc.polyDerM;
		GBlockBound_r.block(0, 0, nNodes, 1) -= 0.5 * _disc.invMM.block(0, 0, nNodes, 1);
		GBlockBound_r.block(0, 1, nNodes, 1) += 0.5 * _disc.invMM.block(0, 0, nNodes, 1);
		GBlockBound_r *= 2 / _disc.deltaZ;
		// auxiliary block [d g^* / d c]
		gStarDC.setZero();
		gStarDC.block(0, nNodes, 1, nNodes + 2) += gBlock.block(0, 0, 1, nNodes + 2);
		gStarDC.block(0, 0, 1, nNodes + 2) += gBlock.block(nNodes - 1, 0, 1, nNodes + 2);
		gStarDC.block(nNodes - 1, nNodes, 1, nNodes + 2) += gBlock.block(nNodes - 1, 0, 1, nNodes + 2);
		gStarDC.block(nNodes - 1, 2 * nNodes, 1, nNodes + 2) += GBlockBound_r.block(0, 0, 1, nNodes + 2);
		gStarDC *= 0.5;

		dispBlock.setZero();
		dispBlock.block(0, nNodes, nNodes, nNodes + 2) += _disc.polyDerM * gBlock - _disc.invMM * B * gBlock;
		dispBlock += _disc.invMM * B * gStarDC;
		dispBlock *= 2 / _disc.deltaZ;

		for (unsigned int comp = 0; comp < nComp; comp++) {
			for (unsigned int i = 0; i < dispBlock.rows(); i++) {
				for (unsigned int j = 0; j < dispBlock.cols() - 1; j++) {
					// row: jump over inlet DOFs and previous cells, add component offset and go node strides from there for each dispersion block entry
					// col: jump over inlet DOFs and previous cells, go back one cell and one node, add component offset and go node strides from there for each dispersion block entry.
					_jac.coeffRef(offC + (nCells - 2) * sCell + comp * sComp + i * sNode,
								  offC + (nCells - 2) * sCell - (nNodes + 1) * sNode + comp * sComp + j * sNode)
								  = -dispBlock(i, j) * _disc.dispersion[comp];
				}
			}
		}

		/*			boundary cells			*/

		// left boundary cell
		dispBlock.setZero();
		gStarDC.setZero();
		gStarDC.block(nNodes - 1, nNodes, 1, nNodes + 2) += GBlockBound_l.block(nNodes - 1, 0, 1, nNodes + 2);
		gStarDC.block(nNodes - 1, 2 * nNodes, 1, nNodes + 2) += gBlock.block(0, 0, 1, nNodes + 2);
		gStarDC *= 0.5;
		dispBlock.block(0, nNodes, nNodes, nNodes + 2) += _disc.polyDerM * GBlockBound_l - _disc.invMM * B * GBlockBound_l;
		dispBlock.block(0, nNodes + 1, nNodes, 2 * nNodes + 1) += _disc.invMM * B * gStarDC.block(0, nNodes + 1, nNodes, 2 * nNodes + 1);
		dispBlock *= 2 / _disc.deltaZ;

		for (int comp = 0; comp < nComp; comp++) {
			for (int i = 0; i < dispBlock.rows(); i++) {
				for (int j = nNodes + 1; j < dispBlock.cols(); j++) {
					// row: jump over inlet DOFs, add component offset and go node strides from there for each dispersion block entry
					// col: jump over inlet DOFs and previous cells, add component offset, adjust for iterator j (-Nnodes-1) and go node strides from there for each dispersion block entry.
					_jac.coeffRef(offC + comp * sComp + i * sNode,
								  offC + comp * sComp + (j - (nNodes + 1)) * sNode)
								  = -dispBlock(i, j) * _disc.dispersion[comp];
				}
			}
		}

		// right boundary cell
		dispBlock.setZero();
		gStarDC.setZero();
		gStarDC.block(0, nNodes, 1, nNodes + 2) += GBlockBound_r.block(0, 0, 1, nNodes + 2);
		gStarDC.block(0, 0, 1, nNodes + 2) += gBlock.block(nNodes - 1, 0, 1, nNodes + 2);
		gStarDC *= 0.5;
		dispBlock.block(0, nNodes, nNodes, nNodes + 2) += _disc.polyDerM * GBlockBound_r - _disc.invMM * B * GBlockBound_r;
		dispBlock += _disc.invMM * B * gStarDC;
		dispBlock *= 2 / _disc.deltaZ;

		for (int comp = 0; comp < nComp; comp++) {
			for (int i = 0; i < dispBlock.rows(); i++) {
				for (int j = 0; j < 2 * nNodes + 1; j++) {
					// row: jump over inlet DOFs and previous cells, add component offset and go node strides from there for each dispersion block entry
					// col: jump over inlet DOFs and previous cells, go back one cell and one node, add component offset and go node strides from there for each dispersion block entry.
					_jac.coeffRef(offC + (nCells - 1) * sCell + comp * sComp + i * sNode,
								  offC + (nCells - 1) * sCell - (nNodes + 1) * sNode + comp * sComp + j * sNode)
								  = -dispBlock(i, j) * _disc.dispersion[comp];
				}
			}
		}

		/*======================================================*/
		/*			Compute Convection Jacobian Block			*/
		/*======================================================*/

		// Convection block [ d RHS_conv / d c ], additionally depends on first entry of previous cell
		MatrixXd convBlock = MatrixXd::Zero(nNodes, nNodes + 1);
		convBlock.block(0, 0, nNodes, 1) += _disc.invMM.block(0, 0, nNodes, 1);
		convBlock.block(0, 1, nNodes, nNodes) -= _disc.polyDerM;
		convBlock.block(0, 1, nNodes, 1) -= _disc.invMM.block(0, 0, nNodes, 1);
		convBlock *= 2 * _disc.velocity / _disc.deltaZ;

		// special inlet DOF treatment for first cell
		for (unsigned int comp = 0; comp < nComp; comp++) {
			for (unsigned int i = 0; i < convBlock.rows(); i++) {
				_jac.coeffRef(offC + comp * sComp + i * sNode, comp * sComp) = -convBlock(i, 0);
				for (unsigned int j = 1; j < convBlock.cols(); j++) {
					_jac.coeffRef(offC + comp * sComp + i * sNode,
								  offC + comp * sComp + (j - 1) * sNode)
								  -= convBlock(i, j);
				}
			}
		}
		for (unsigned int cell = 1; cell < nCells; cell++) {
			for (unsigned int comp = 0; comp < nComp; comp++) {
				for (unsigned int i = 0; i < convBlock.rows(); i++) {
					for (unsigned int j = 0; j < convBlock.cols(); j++) {
						// row: jump over inlet DOFs and previous cells, add component offset and go node strides from there for each convection block entry
						// col: jump over inlet DOFs and previous cells, go back one node, add component offset and go node strides from there for each convection block entry
						_jac.coeffRef(offC + cell * sCell + comp * sComp + i * sNode,
									  offC + cell * sCell - sNode + comp * sComp + j * sNode)
									  -= convBlock(i, j);
					}
				}
			}
		}

		return 0;
	}

	/**
	* @brief analytically calculates the isotherm jacobian
	* @return 1 if jacobain estimation fits the predefined pattern of the jacobian, 0 if not.
	*/
	int calcIsothermJacobian(double t, unsigned int secIdx, const double* const y, util::ThreadLocalStorage& threadLocalMem) {

		Indexer idxr(_disc);

		// set rowIterator and local state pointer to first solid concentration
		linalg::BandedEigenSparseRowIterator rowIterator(_jac, idxr.offsetC() + idxr.strideColLiquid());
		const double* yLocal = y + idxr.offsetC() + idxr.strideColLiquid();

		// Offset from the first component of the mobile phase to the first bound state
		unsigned int offSetCp = _disc.nComp;

		for (unsigned int point = 0; point < _disc.nPoints; point++) {

			double z = _disc.deltaZ * std::floor(point / _disc.nNodes)
				+ 0.5 * _disc.deltaZ * (1 + _disc.nodes[point % _disc.nNodes]);

			_binding[0]->analyticJacobian(t, secIdx, ColumnPosition{ z, 0.0, 0.0 }, yLocal, offSetCp, rowIterator, threadLocalMem.get());

			// set rowIterator and y to first bound concentration of next discrete point
			rowIterator += idxr.strideColNode();
			yLocal += idxr.strideColNode();

		}

		if (!_jac.isCompressed()) // if matrix lost its compressed storage, the pattern did not fit.
			return 0;

		return 1;
	}

	/**
	* @brief adds time derivative to the jacobian
	*/
	void addTimederJacobian(double alpha) {

		Indexer idxr(_disc);

		int sNode = idxr.strideColNode();
		int offC = idxr.offsetC();

		// =================================================================================================== //
		//	 Time derivative Jacobian: d Residual / d y_t   												   //
		// =================================================================================================== //

		linalg::BandedEigenSparseRowIterator jac(_jacDisc, idxr.offsetC());

		double Beta = (1.0 - _disc.porosity) / _disc.porosity;

		for (unsigned int point = 0; point < _disc.nPoints; point++) {
			for (unsigned int comp = 0; comp < _disc.nComp; comp++) {
				// d convDispRHS / d c_t
				jac[0] += alpha;

				if (_disc.nBound[comp]) { // either one or null; no loop necessary
					// d convDispRHS / d q_t
					jac[idxr.strideColLiquid() - comp + idxr.offsetBoundComp(comp)] += alpha * Beta;
				}
				++jac;
			}
			for (unsigned int comp = 0; comp < _disc.nComp; comp++) {
				if (_disc.nBound[comp]) { // either one or null; no loop over bound states necessary
					if (_disc.isKinetic[idxr.offsetBoundComp(comp)]) {
						// d isotherm / d q_t
						jac[0] += alpha;
					}
					++jac;
				}
			}
		}

	}

};

} // namespace model
} // namespace cadet

#endif  // LIBCADET_LUMPEDRATEMODELWITHOUTPORESDG_HPP_
